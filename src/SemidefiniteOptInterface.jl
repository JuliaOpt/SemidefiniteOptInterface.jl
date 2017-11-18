module SemidefiniteOptInterface

using MathOptInterface
const MOI = MathOptInterface

using MathOptInterfaceUtilities
const MOIU = MathOptInterfaceUtilities
MOIU.@instance SDInstance () (EqualTo, GreaterThan, LessThan, Interval) (Zeros, Nonnegatives, Nonpositives, SecondOrderCone, RotatedSecondOrderCone, PositiveSemidefiniteConeTriangle, PositiveSemidefiniteConeScaled) () (SingleVariable,) (ScalarAffineFunction,) (VectorOfVariables,) (VectorAffineFunction,)

abstract type AbstractSDSolver <: MOI.AbstractSolver end

MOI.get(m::AbstractSDSolver, ::Union{MOI.SupportsDuals,
                                     MOI.SupportsAddVariableAfterSolve,
                                     MOI.SupportsAddConstraintAfterSolve,
                                     MOI.SupportsDeleteVariable,
                                     MOI.SupportsDeleteConstraint}) = true

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
const IL = MOI.Interval
const CS = MOI.PositiveSemidefiniteConeScaled
const SO = MOI.SecondOrderCone
const RS = MOI.RotatedSecondOrderCone
const DS = MOI.PositiveSemidefiniteConeTriangle
const SupportedSets = Union{ZS, NS, PS, DS}
const BridgedSets = Union{IL, CS, SO, RS}

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
    int::Vector{SplitIntervalBridge{Float64}}
    psdcs::Vector{PSDCScaledBridge{Float64}}
    soc::Vector{SOCtoPSDCBridge{Float64}}
    rsoc::Vector{RSOCtoPSDCBridge{Float64}}
    double::Vector{CR} # created when there are two cones for same variable
    function SOItoMOIBridge(solver::ST, sdsolver::SIT) where {ST, SIT}
        new{ST, SIT}(solver, SDInstance{Float64}(), sdsolver,
            0.0, 0, 0, 0,
            Int[],
            IntSet(),
            Vector{Tuple{Int,Int,Int,Float64}}[],
            UnitRange{Int}[],
            Tuple{Int, Int, Int, Float64}[],
            Int[],
            SplitIntervalBridge{Float64}[],
            PSDCScaledBridge{Float64}[],
            SOCtoPSDCBridge{Float64}[],
            RSOCtoPSDCBridge{Float64}[],
            CR[])
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

MOI.isvalid(m::SOItoMOIBridge, r::Union{VR, CR}) = MOI.isvalid(m.sdinstance, r)
MOI.candelete(m::SOItoMOIBridge, r::Union{VR, CR}) = MOI.candelete(m.sdinstance, r)
MOI.delete!(m::SOItoMOIBridge, r::Union{VR, CR}) = MOI.delete!(m.sdinstance, r)
MOI.addconstraint!(m::SOItoMOIBridge, f::Union{ASF, AVF}, s) = MOI.addconstraint!(m.sdinstance, f, s)

const InstanceAttributeRef = Union{MOI.ConstraintFunction, MOI.ConstraintSet}
MOI.canget(m::SOItoMOIBridge, a::InstanceAttributeRef, ref::CR) = MOI.canget(m.sdinstance, a, ref)
MOI.get(m::SOItoMOIBridge, a::InstanceAttributeRef, ref::CR) = MOI.get(m.sdinstance, a, ref)

const InstanceAttribute = Union{MOI.NumberOfVariables,
                                MOI.NumberOfConstraints,
                                MOI.ListOfConstraints,
                                MOI.ObjectiveFunction,
                                MOI.ObjectiveSense}
MOI.canget(m::SOItoMOIBridge, a::InstanceAttribute) = MOI.canget(m.sdinstance, a)

function MOI.get(m::SOItoMOIBridge, a::InstanceAttribute)
    MOI.get(m.sdinstance, a)
end

function MOI.optimize!(m::SOItoMOIBridge)
    loadprimal!(m)
    MOI.optimize!(m.sdsolver)
end

MOI.canmodifyconstraint(m::SOItoMOIBridge, cr::CR, change) = MOI.canmodifyconstraint(m.sdinstance, cr, change)
function MOI.modifyconstraint!(m::SOItoMOIBridge, cr::CR, change)
    MOI.modifyconstraint!(m.sdinstance, cr, change)
end

# Objective

MOI.set!(m::SOItoMOIBridge, att::Union{MOI.ObjectiveSense, MOI.ObjectiveFunction}, arg) = MOI.set!(m.sdinstance, att, arg)
MOI.canmodifyobjective(m::SOItoMOIBridge, change::MOI.AbstractFunctionModification) = MOI.canmodifyobjective(m.sdinstance, change)
MOI.modifyobjective!(m::SOItoMOIBridge, change::MOI.AbstractFunctionModification) = MOI.modifyobjective!(m.sdinstance, change)

_objsgn(m) = m.sdinstance.sense == MOI.MinSense ? -1 : 1
MOI.canget(m::SOItoMOIBridge, ::MOI.ObjectiveValue) = true
function MOI.get(m::SOItoMOIBridge, ::MOI.ObjectiveValue)
    m.objshift + _objsgn(m) * getprimalobjectivevalue(m.sdsolver) + m.sdinstance.objective.constant
end

# Attributes

MOI.canget(m::AbstractSDSolverInstance, ::MOI.TerminationStatus) = true
const SolverStatus = Union{MOI.TerminationStatus, MOI.PrimalStatus, MOI.DualStatus}
MOI.canget(m::SOItoMOIBridge, s::SolverStatus) = MOI.canget(m.sdsolver, s)
MOI.get(m::SOItoMOIBridge, s::SolverStatus) = MOI.get(m.sdsolver, s)


MOI.canget(m::SOItoMOIBridge, ::MOI.ResultCount) = true
MOI.get(m::SOItoMOIBridge, ::MOI.ResultCount) = 1

MOI.canget(m::SOItoMOIBridge, ::Union{MOI.VariablePrimal,
                                      MOI.ConstraintPrimal,
                                      MOI.ConstraintDual}, ref::Union{CR, VR}) = true

MOI.canget(m::SOItoMOIBridge, ::Union{MOI.VariablePrimal,
                                      MOI.ConstraintPrimal,
                                      MOI.ConstraintDual}, ref::Vector{R}) where R <: Union{CR, VR} = true


function MOI.get(m::SOItoMOIBridge, ::MOI.VariablePrimal, vr::VR)
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
function MOI.get(m::SOItoMOIBridge, vp::MOI.VariablePrimal, vr::Vector{VR})
    MOI.get.(m, vp, vr)
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
        if i != j
            coef *= 2 # We should take X[blk][i, j] + X[blk][j, i] but they are equal
        end
        coef * X[blk][i, j]
    end
end

function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR)
    _getattribute(m, cr, getslack) + _getconstant(m, MOI.get(m, MOI.ConstraintSet(), cr))
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
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR{<:VF, <:ZS})
    _getattribute(m, cr, getdual) + getvardual(m, MOI.get(m, MOI.ConstraintFunction(), cr))
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR{<:VF, <:SupportedSets})
    getvardual(m, MOI.get(m, MOI.ConstraintFunction(), cr))
end

function getdual(m::SOItoMOIBridge, c::Int)
    if c == 0
        0.
    else
        -gety(m.sdsolver)[c]
    end
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR)
    _getattribute(m, cr, getdual)
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR{F, DS}) where F
    scalevec!(_getattribute(m, cr, getdual), 1/2)
end

# Bridges

function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintDual, cr::CR{F, IL{Float64}}) where F
    MOI.get(m, a, m.int[m.bridgemap[cr.value]])
end
function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR{F, IL{Float64}}) where F
    MOI.get(m, a, m.int[m.bridgemap[cr.value]])
end

function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintDual, cr::CR{F, CS}) where F
    MOI.get(m, a, m.psdcs[m.bridgemap[cr.value]])
end
function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR{F, CS}) where F
    MOI.get(m, a, m.psdcs[m.bridgemap[cr.value]])
end

function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintDual, cr::CR{F, SO}) where F
    MOI.get(m, a, m.soc[m.bridgemap[cr.value]])
end
function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR{F, SO}) where F
    MOI.get(m, a, m.soc[m.bridgemap[cr.value]])
end

function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintDual, cr::CR{F, RS}) where F
    MOI.get(m, a, m.rsoc[m.bridgemap[cr.value]])
end
function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR{F, RS}) where F
    MOI.get(m, a, m.rsoc[m.bridgemap[cr.value]])
end


end # module
