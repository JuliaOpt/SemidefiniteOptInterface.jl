module SemidefiniteOptInterface

using MathOptInterface
const MOI = MathOptInterface

using MathOptInterfaceUtilities
const MOIU = MathOptInterfaceUtilities
MOIU.@instance SDInstance () (EqualTo, GreaterThan, LessThan) (Zeros, Nonnegatives, Nonpositives, SecondOrderCone, PositiveSemidefiniteConeTriangle) () (SingleVariable,) (ScalarAffineFunction,) (VectorOfVariables,) (VectorAffineFunction,)

abstract type AbstractSDSolver <: MOI.AbstractSolver end

MOI.getattribute(m::AbstractSDSolver, ::MOI.SupportsDuals) = true

abstract type AbstractSDSolverInstance <: MOI.AbstractSolverInstance end

include("interface.jl")

const SVF = MOI.SingleVariable
const VVF = MOI.VectorOfVariables
const VF  = Union{SVF, VVF}
const SAF{T} = MOI.ScalarAffineFunction{T}
const VAF{T} = MOI.VectorAffineFunction{T}
const AF{T}  = Union{SAF{T}, VAF{T}}
const ASF{T} = Union{SVF, SAF{T}}
const AVF{T} = Union{VVF, VAF{T}}

const ZS = Union{MOI.EqualTo, MOI.Zeros}
const NS = Union{MOI.GreaterThan, MOI.Nonnegatives}
const PS = Union{MOI.LessThan, MOI.Nonpositives}
const SO = MOI.SecondOrderCone
const DS = MOI.PositiveSemidefiniteConeTriangle
const SupportedSets = Union{ZS, NS, PS, DS}
const BridgedSets = Union{SO}

const VR = MOI.VariableReference
const CR{FT, ST} = MOI.ConstraintReference{FT, ST}

include("setbridges.jl")

mutable struct SOItoMOIBridge{ST <: AbstractSDSolver, SIT <: AbstractSDSolverInstance} <: MOI.AbstractSolverInstance
    solver::ST
    sdinstance::SDInstance{Float64}
    sdsolver::SIT
    objshift::Float64
    nconstrs::Int
    constr::Int
    nblocks::Int
    blockdims::Vector{Int}
    free::IntSet
    varmap::Vector{Vector{Tuple{Int, Int, Int, Float64, Float64}}} # Variable Reference value vi -> blk, i, j, coef, shift # x = sum coef * X[blk][i, j] + shift
    constrmap::Vector{UnitRange{Int}} # Constraint Reference value ci -> cs
    slackmap::Vector{Tuple{Int, Int, Int, Float64}} # c -> blk, i, j, coef
    # Bridges
    bridgemap::Vector{Int}
    soc::Vector{SOCtoPSDCBridge{Float64}}
    function SOItoMOIBridge(solver::ST, sdsolver::SIT) where {ST, SIT}
        new{ST, SIT}(solver, SDInstance{Float64}(), sdsolver,
            0.0, 0, 0, 0,
            Int[],
            IntSet(),
            Vector{Tuple{Int,Int,Int,Float64}}[],
            UnitRange{Int}[],
            Tuple{Int, Int, Int, Float64}[],
            Int[],
            SOCtoPSDCBridge{Float64}[])
    end
    function SOItoMOIBridge(solver::ST) where ST
        SOItoMOIBridge(solver, SDSolverInstance(solver))
    end
end

MOI.SolverInstance(s::AbstractSDSolver) = SOItoMOIBridge(s)

MOI.supportsproblem(m::AbstractSDSolver, ::Type{<:MOI.AbstractFunction}, constrtypes) = false

_supports(x::Type, y::Type) = false
function _supports(::Type{<:Union{VF, AF}}, ::Type{<:Union{SupportedSets, BridgedSets}})
    true
end
MOI.supportsproblem(m::AbstractSDSolver, ::Type{<:SAF}, constrtypes) = reduce((b, ct) -> b & _supports(ct...), true, constrtypes)

include("load.jl")

# Variables

MOI.addvariable!(m::SOItoMOIBridge) = MOI.addvariable!(m.sdinstance)
MOI.addvariables!(m::SOItoMOIBridge, n::Integer) = MOI.addvariables!(m.sdinstance, n)

# Constraints

MOI.delete!(m::SOItoMOIBridge, r::Union{MOI.VariableReference, CR}) = MOI.delete!(m.sdinstance, r)
MOI.addconstraint!(m::SOItoMOIBridge, f, s) = MOI.addconstraint!(m.sdinstance, f, s)

const InstanceAttributeRef = Union{MOI.ConstraintFunction, MOI.ConstraintSet}
MOI.cangetattribute(m::SOItoMOIBridge, a::InstanceAttributeRef, ref) = MOI.cangetattribute(m.sdinstance, a, ref)
MOI.getattribute(m::SOItoMOIBridge, a::InstanceAttributeRef, ref) = MOI.getattribute(m.sdinstance, a, ref)

const InstanceAttribute = Union{MOI.NumberOfVariables,
                                MOI.NumberOfConstraints,
                                MOI.ListOfConstraints,
                                MOI.ObjectiveFunction,
                                MOI.Sense}
MOI.cangetattribute(m::SOItoMOIBridge, a::InstanceAttribute) = MOI.cangetattribute(m.sdinstance, a)

function MOI.getattribute(m::SOItoMOIBridge, a::InstanceAttribute)
    MOI.getattribute(m.sdinstance, a)
end

function MOI.optimize!(m::SOItoMOIBridge)
    loadprimal!(m)
    MOI.optimize!(m.sdsolver)
end

function MOI.modifyconstraint!(m::SOItoMOIBridge, cr::CR, change)
    MOI.modifyconstraint!(m.sdinstance, cr, change)
end

# Objective

MOI.setobjective!(m::SOItoMOIBridge, sense::MOI.OptimizationSense, f) = MOI.setobjective!(m.sdinstance, sense, f)
MOI.modifyobjective!(m::SOItoMOIBridge, change::MOI.AbstractFunctionModification) = MOI.modifyobjective!(m.sdinstance, change)

_objsgn(m) = m.sdinstance.sense == MOI.MinSense ? -1 : 1
MOI.cangetattribute(m::SOItoMOIBridge, ::MOI.ObjectiveValue) = true
function MOI.getattribute(m::SOItoMOIBridge, ::MOI.ObjectiveValue)
    m.objshift + _objsgn(m) * getprimalobjectivevalue(m.sdsolver) + m.sdinstance.objective.constant
end

# Attributes

MOI.cangetattribute(m::AbstractSDSolverInstance, ::MOI.TerminationStatus) = true
const SolverStatus = Union{MOI.TerminationStatus, MOI.PrimalStatus, MOI.DualStatus}
MOI.cangetattribute(m::SOItoMOIBridge, s::SolverStatus) = MOI.cangetattribute(m.sdsolver, s)
MOI.getattribute(m::SOItoMOIBridge, s::SolverStatus) = MOI.getattribute(m.sdsolver, s)


MOI.cangetattribute(m::SOItoMOIBridge, ::Union{MOI.VariablePrimal,
                                               MOI.ConstraintPrimal,
                                               MOI.ConstraintDual}, ref) = true

function MOI.getattribute(m::SOItoMOIBridge, ::MOI.VariablePrimal, vr::MOI.VariableReference)
    X = getX(m.sdsolver)
    x = 0.0
    for (blk, i, j, coef, shift) in m.varmap[vr.value]
        x += shift
        if blk != 0
            x += X[blk][i, j] * sign(coef)
        end
    end
    x
end
function MOI.getattribute(m::SOItoMOIBridge, vp::MOI.VariablePrimal, vr::Vector{MOI.VariableReference})
    MOI.getattribute.(m, vp, vr)
end

function _getattribute(m::SOItoMOIBridge, cr::CR{<:ASF}, f)
    cs = m.constrmap[cr.value]
    @assert length(cs) == 1
    f(m, first(cs))
end
function _getattribute(m::SOItoMOIBridge, cr::CR{<:AVF}, f)
    f.(m, m.constrmap[cr.value])
end

function getslack(m::SOItoMOIBridge, c::Int)
    X = getX(m.sdsolver)
    blk, i, j, coef = m.slackmap[c]
    if iszero(blk)
        0.0
    else
        X[blk][i, j]
    end
end

function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR)
    _getattribute(m, cr, getslack) + _getconstant(m, MOI.getattribute(m, MOI.ConstraintSet(), cr))
end

function getvardual(m::SOItoMOIBridge, vi::UInt64)
    Z = getZ(m.sdsolver)
    z = 0.
    for (blk, i, j, coef) in m.varmap[vi]
        if blk != 0
            z += Z[blk][i, j] * sign(coef)
        end
    end
    z
end
getvardual(m::SOItoMOIBridge, f::SVF) = getvardual(m, f.variable.value)
getvardual(m::SOItoMOIBridge, f::VVF) = map(vr -> getvardual(m, vr.value), f.variables)
function MOI.getattribute(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR{<:VF, <:ZS})
    _getattribute(m, cr, getdual) + getvardual(m, MOI.getattribute(m, MOI.ConstraintFunction(), cr))
end
function MOI.getattribute(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR{<:VF, <:SupportedSets})
    getvardual(m, MOI.getattribute(m, MOI.ConstraintFunction(), cr))
end

function getdual(m::SOItoMOIBridge, c::Int)
    if c == 0
        0.
    else
        -gety(m.sdsolver)[c]
    end
end
function MOI.getattribute(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR)
    _getattribute(m, cr, getdual)
end

function MOI.getattribute{F}(m::SOItoMOIBridge, a::MOI.ConstraintDual, cr::CR{F, MOI.SecondOrderCone})
    MOI.getattribute(m, a, m.soc[m.bridgemap[cr.value]])
end
function MOI.getattribute{F}(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR{F, MOI.SecondOrderCone})
    MOI.getattribute(m, a, m.soc[m.bridgemap[cr.value]])
end

end # module
