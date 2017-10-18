module SemidefiniteOptInterface

using MathOptInterface
const MOI = MathOptInterface

using MathOptInterfaceUtilities
const MOIU = MathOptInterfaceUtilities
MOIU.@instance SDInstance () (EqualTo, GreaterThan, LessThan, Interval) (Zeros, Nonnegatives, Nonpositives, SecondOrderCone, RotatedSecondOrderCone, PositiveSemidefiniteConeTriangle, PositiveSemidefiniteConeScaled) () (SingleVariable,) (ScalarAffineFunction,) (VectorOfVariables,) (VectorAffineFunction,)

abstract type AbstractSDSolver <: MOI.AbstractSolver end

MOI.getattribute(m::AbstractSDSolver, ::Union{MOI.SupportsDuals,
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

mutable struct SOItoMOIBridge{T, ST <: AbstractSDSolver, SIT <: AbstractSDSolverInstance} <: MOI.AbstractSolverInstance
    solver::ST
    sdinstance::SDInstance{T}
    sdsolver::SIT
    objshift::T
    nconstrs::Int
    constr::Int
    nblocks::Int
    blockdims::Vector{Int}
    free::IntSet
    varmap::Vector{Vector{Tuple{Int, Int, Int, T, T}}} # Variable Reference value vi -> blk, i, j, coef, shift # x = sum coef * X[blk][i, j] + shift
    constrmap::Vector{UnitRange{Int}} # Constraint Reference value ci -> cs
    slackmap::Vector{Tuple{Int, Int, Int, T}} # c -> blk, i, j, coef
    # Bridges
    bridgemap::Vector{Int}
    int::Vector{SplitIntervalBridge{T}}
    psdcs::Vector{PSDCScaledBridge{T}}
    soc::Vector{SOCtoPSDCBridge{T}}
    rsoc::Vector{RSOCtoPSDCBridge{T}}
    double::Vector{CR} # created when there are two cones for same variable
    function SOItoMOIBridge(T, solver::ST, sdsolver::SIT) where {ST, SIT}
        new{T, ST, SIT}(solver, SDInstance{T}(), sdsolver,
            zero(T), 0, 0, 0,
            Int[],
            IntSet(),
            Vector{Tuple{Int,Int,Int,T}}[],
            UnitRange{Int}[],
            Tuple{Int, Int, Int, T}[],
            Int[],
            SplitIntervalBridge{T}[],
            PSDCScaledBridge{T}[],
            SOCtoPSDCBridge{T}[],
            RSOCtoPSDCBridge{T}[],
            CR[])
    end
    function SOItoMOIBridge(solver::ST) where ST
        #T = method_exists(MOI.coefficienttype, Tuple{ST}) ? MOI.coefficienttype(solver) : Float64
        SOItoMOIBridge(MOI.coefficienttype(solver), solver, SDSolverInstance(solver))
    end
end


MOI.SolverInstance(s::AbstractSDSolver) = SOItoMOIBridge(s)

MOI.coefficienttype(::SOItoMOIBridge{T}) where T = T

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
MOI.cangetattribute(m::SOItoMOIBridge, a::InstanceAttributeRef, ref::CR) = MOI.cangetattribute(m.sdinstance, a, ref)
MOI.getattribute(m::SOItoMOIBridge, a::InstanceAttributeRef, ref::CR) = MOI.getattribute(m.sdinstance, a, ref)

const InstanceAttribute = Union{MOI.NumberOfVariables,
                                MOI.NumberOfConstraints,
                                MOI.ListOfConstraints,
                                MOI.ObjectiveFunction,
                                MOI.ObjectiveSense}
MOI.cangetattribute(m::SOItoMOIBridge, a::InstanceAttribute) = MOI.cangetattribute(m.sdinstance, a)

function MOI.getattribute(m::SOItoMOIBridge, a::InstanceAttribute)
    MOI.getattribute(m.sdinstance, a)
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

MOI.setattribute!(m::SOItoMOIBridge, att::Union{MOI.ObjectiveSense, MOI.ObjectiveFunction}, arg) = MOI.setattribute!(m.sdinstance, att, arg)
MOI.canmodifyobjective(m::SOItoMOIBridge, change::MOI.AbstractFunctionModification) = MOI.canmodifyobjective(m.sdinstance, change)
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


MOI.cangetattribute(m::SOItoMOIBridge, ::MOI.ResultCount) = true
MOI.getattribute(m::SOItoMOIBridge, ::MOI.ResultCount) = 1

MOI.cangetattribute(m::SOItoMOIBridge, ::Union{MOI.VariablePrimal,
                                               MOI.ConstraintPrimal,
                                               MOI.ConstraintDual}, ref::Union{CR, VR}) = true

MOI.cangetattribute(m::SOItoMOIBridge, ::Union{MOI.VariablePrimal,
                                               MOI.ConstraintPrimal,
                                               MOI.ConstraintDual}, ref::Vector{R}) where R <: Union{CR, VR} = true


function MOI.getattribute(m::SOItoMOIBridge{T}, ::MOI.VariablePrimal, vr::VR) where T
    X = getX(m.sdsolver)
    x = zero(T)
    for (blk, i, j, coef, shift) in m.varmap[vr.value]
        x += shift
        if blk != 0
            x += X[blk][i, j] * sign(coef)
        end
    end
    x
end
function MOI.getattribute(m::SOItoMOIBridge, vp::MOI.VariablePrimal, vr::Vector{VR})
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

function getslack(m::SOItoMOIBridge{T}, c::Int) where T
    X = getX(m.sdsolver)
    blk, i, j, coef = m.slackmap[c]
    if iszero(blk)
        zero(T)
    else
        if i != j
            coef *= 2 # We should take X[blk][i, j] + X[blk][j, i] but they are equal
        end
        coef * X[blk][i, j]
    end
end

function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR)
    _getattribute(m, cr, getslack) + _getconstant(m, MOI.getattribute(m, MOI.ConstraintSet(), cr))
end

function getvardual(m::SOItoMOIBridge{T}, vi::UInt64) where T
    Z = getZ(m.sdsolver)
    z = zero(T)
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

function getdual(m::SOItoMOIBridge{T}, c::Int) where T
    if c == 0
        zero(T)
    else
        -gety(m.sdsolver)[c]
    end
end
function MOI.getattribute(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR)
    _getattribute(m, cr, getdual)
end
function MOI.getattribute(m::SOItoMOIBridge{T}, ::MOI.ConstraintDual, cr::CR{F, DS}) where {T, F}
    scalevec!(_getattribute(m, cr, getdual), one(T)/2)
end

# Bridges

function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintDual, cr::CR{F, IL{Float64}}) where F
    MOI.getattribute(m, a, m.int[m.bridgemap[cr.value]])
end
function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR{F, IL{Float64}}) where F
    MOI.getattribute(m, a, m.int[m.bridgemap[cr.value]])
end

function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintDual, cr::CR{F, CS}) where F
    MOI.getattribute(m, a, m.psdcs[m.bridgemap[cr.value]])
end
function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR{F, CS}) where F
    MOI.getattribute(m, a, m.psdcs[m.bridgemap[cr.value]])
end

function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintDual, cr::CR{F, SO}) where F
    MOI.getattribute(m, a, m.soc[m.bridgemap[cr.value]])
end
function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR{F, SO}) where F
    MOI.getattribute(m, a, m.soc[m.bridgemap[cr.value]])
end

function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintDual, cr::CR{F, RS}) where F
    MOI.getattribute(m, a, m.rsoc[m.bridgemap[cr.value]])
end
function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR{F, RS}) where F
    MOI.getattribute(m, a, m.rsoc[m.bridgemap[cr.value]])
end


end # module
