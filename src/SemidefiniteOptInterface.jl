module SemidefiniteOptInterface

using MathOptInterface
const MOI = MathOptInterface

using MathOptInterfaceUtilities
const MOIU = MathOptInterfaceUtilities
MOIU.@instance SDInstance () (EqualTo, GreaterThan, LessThan) (Zeros, Nonnegatives, Nonpositives, PositiveSemidefiniteConeTriangle) () (SingleVariable,) (ScalarAffineFunction,) (VectorOfVariables,) (VectorAffineFunction,)

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
const DS = MOI.PositiveSemidefiniteConeTriangle
const SupportedSets = Union{ZS, NS, PS, DS}
const BridgedSets = Union{MOI.Interval,
                          MOI.SecondOrderCone,
                          MOI.RotatedSecondOrderCone,
                          MOI.PositiveSemidefiniteConeScaled}

const VR = MOI.VariableReference
const CR{FT, ST} = MOI.ConstraintReference{FT, ST}

mutable struct SOItoMOIBridge{T, SIT <: AbstractSDSolverInstance} <: MOI.AbstractSolverInstance
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
    double::Vector{CR} # created when there are two cones for same variable
    function SOItoMOIBridge{T}(sdsolver::SIT) where {T, SIT}
        new{T, SIT}(SDInstance{T}(), sdsolver,
            zero(T), 0, 0, 0,
            Int[],
            IntSet(),
            Vector{Tuple{Int, Int, Int, T}}[],
            UnitRange{Int}[],
            Tuple{Int, Int, Int, T}[],
            Float64[])
    end
end

SDOIInstance(sdsolver::AbstractSDSolverInstance, T=Float64) = PSDCScaled{T}(RSOCtoPSDC{T}(SOCtoPSDC{T}(SplitInterval{T}(SOItoMOIBridge{T}(sdsolver)))))

include("setbridges.jl")
@bridge SplitInterval MOIU.SplitIntervalBridge () (Interval,) () () () (ScalarAffineFunction,) () ()
@bridge PSDCScaled PSDCScaledBridge () () (PositiveSemidefiniteConeScaled,) () () () (VectorOfVariables,) (VectorAffineFunction,)
@bridge SOCtoPSDC SOCtoPSDCBridge () () (SecondOrderCone,) () () () (VectorOfVariables,) (VectorAffineFunction,)
@bridge RSOCtoPSDC RSOCtoPSDCBridge () () (RotatedSecondOrderCone,) () () () (VectorOfVariables,) (VectorAffineFunction,)

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


function MOI.get(m::SOItoMOIBridge{T}, ::MOI.VariablePrimal, vr::VR) where T
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

function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR)
    _getattribute(m, cr, getslack) + _getconstant(m, MOI.get(m, MOI.ConstraintSet(), cr))
end
# These constraints do not create any constraint or slack, it is just variables
_getvarprimal(m, sv::SVF) = MOI.get(m, MOI.VariablePrimal(), sv.variable)
_getvarprimal(m, vv::VVF) = MOI.get(m, MOI.VariablePrimal(), vv.variables)
function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR{<:VF, <:Union{NS, PS, DS}})
    _getvarprimal(m, MOI.get(m, MOI.ConstraintFunction(), cr))
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
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR{<:VF, <:ZS})
    _getattribute(m, cr, getdual) + getvardual(m, MOI.get(m, MOI.ConstraintFunction(), cr))
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR{<:VF, <:SupportedSets})
    getvardual(m, MOI.get(m, MOI.ConstraintFunction(), cr))
end

function getdual(m::SOItoMOIBridge{T}, c::Int) where T
    if c == 0
        zero(T)
    else
        -gety(m.sdsolver)[c]
    end
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR)
    _getattribute(m, cr, getdual)
end
function MOI.get(m::SOItoMOIBridge{T}, ::MOI.ConstraintDual, cr::CR{F, DS}) where {T,F}
    scalevec!(_getattribute(m, cr, getdual), one(T)/2)
end

end # module
