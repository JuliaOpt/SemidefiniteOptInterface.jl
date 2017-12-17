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
                          MOI.RotatedSecondOrderCone}

const VI = MOI.VariableIndex
const CI{FT, ST} = MOI.ConstraintIndex{FT, ST}

mutable struct SOItoMOIBridge{T, SIT <: AbstractSDSolverInstance} <: MOI.AbstractSolverInstance
    sdinstance::SDInstance{T}
    sdsolver::SIT
    objshift::T
    nconstrs::Int
    nblocks::Int
    blockdims::Vector{Int}
    free::Set{VI}
    varmap::Dict{VI, Vector{Tuple{Int, Int, Int, T, T}}} # Variable Index vi -> blk, i, j, coef, shift # x = sum coef * X[blk][i, j] + shift
    constrmap::Dict{CI, UnitRange{Int}} # Constraint Index ci -> cs
    slackmap::Vector{Tuple{Int, Int, Int, T}} # c -> blk, i, j, coef
    double::Vector{CI} # created when there are two cones for same variable
    function SOItoMOIBridge{T}(sdsolver::SIT) where {T, SIT}
        new{T, SIT}(SDInstance{T}(), sdsolver,
            zero(T), 0, 0,
            Int[],
            Set{VI}(),
            Dict{VI, Tuple{Int, Int, Int, T}}(),
            Dict{CI, UnitRange{Int}}(),
            Tuple{Int, Int, Int, T}[],
            Float64[])
    end
end

SDOIInstance(sdsolver::AbstractSDSolverInstance, T=Float64) = RootDet{T}(GeoMean{T}(RSOCtoPSDC{T}(SOCtoPSDC{T}(SplitInterval{T}(SOItoMOIBridge{T}(sdsolver))))))

include("data.jl")

include("setbridges.jl")
@bridge SplitInterval MOIU.SplitIntervalBridge () (Interval,) () () () (ScalarAffineFunction,) () ()
@bridge SOCtoPSDC SOCtoPSDCBridge () () (SecondOrderCone,) () () () (VectorOfVariables,) (VectorAffineFunction,)
@bridge RSOCtoPSDC RSOCtoPSDCBridge () () (RotatedSecondOrderCone,) () () () (VectorOfVariables,) (VectorAffineFunction,)
@bridge GeoMean MOIU.GeoMeanBridge () () (GeometricMeanCone,) () () () (VectorOfVariables,) (VectorAffineFunction,)
@bridge RootDet MOIU.RootDetBridge () () (RootDetConeTriangle,) () () () (VectorOfVariables,) (VectorAffineFunction,)

include("load.jl")

# Constraints

function MOI.optimize!(m::SOItoMOIBridge)
    loadprimal!(m)
    MOI.optimize!(m.sdsolver)
end

# Objective

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
                                      MOI.ConstraintDual}, ref::Union{CI, VI}) = true

MOI.canget(m::SOItoMOIBridge, ::Union{MOI.VariablePrimal,
                                      MOI.ConstraintPrimal,
                                      MOI.ConstraintDual}, ref::Vector{R}) where R <: Union{CI, VI} = true


function MOI.get(m::SOItoMOIBridge{T}, ::MOI.VariablePrimal, vi::VI) where T
    X = getX(m.sdsolver)
    x = zero(T)
    for (blk, i, j, coef, shift) in m.varmap[vi]
        x += shift
        if blk != 0
            x += X[blk][i, j] * sign(coef)
        end
    end
    x
end
function MOI.get(m::SOItoMOIBridge, vp::MOI.VariablePrimal, vi::Vector{VI})
    MOI.get.(m, vp, vi)
end

function _getattribute(m::SOItoMOIBridge, ci::CI{<:ASF}, f)
    cs = m.constrmap[ci]
    @assert length(cs) == 1
    f(m, first(cs))
end
function _getattribute(m::SOItoMOIBridge, ci::CI{<:AVF}, f)
    f.(m, m.constrmap[ci])
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

function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CI)
    _getattribute(m, cr, getslack) + _getconstant(m, MOI.get(m, MOI.ConstraintSet(), cr))
end
# These constraints do not create any constraint or slack, it is just variables
_getvarprimal(m, sv::SVF) = MOI.get(m, MOI.VariablePrimal(), sv.variable)
_getvarprimal(m, vv::VVF) = MOI.get(m, MOI.VariablePrimal(), vv.variables)
function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CI{<:VF, <:Union{NS, PS, DS}})
    _getvarprimal(m, MOI.get(m, MOI.ConstraintFunction(), cr))
end

function getvardual(m::SOItoMOIBridge{T}, vi::VI) where T
    Z = getZ(m.sdsolver)
    z = zero(T)
    for (blk, i, j, coef) in m.varmap[vi]
        if blk != 0
            z += Z[blk][i, j] * sign(coef)
        end
    end
    z
end
getvardual(m::SOItoMOIBridge, f::SVF) = getvardual(m, f.variable)
getvardual(m::SOItoMOIBridge, f::VVF) = map(vi -> getvardual(m, vi), f.variables)
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CI{<:VF, <:ZS})
    _getattribute(m, cr, getdual) + getvardual(m, MOI.get(m, MOI.ConstraintFunction(), cr))
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CI{<:VF, <:SupportedSets})
    getvardual(m, MOI.get(m, MOI.ConstraintFunction(), cr))
end

function getdual(m::SOItoMOIBridge{T}, c::Int) where T
    if c == 0
        zero(T)
    else
        -gety(m.sdsolver)[c]
    end
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CI)
    _getattribute(m, cr, getdual)
end
function scalevec!(v, c)
    d = div(isqrt(1+8length(v))-1, 2)
    @assert div(d*(d+1), 2) == length(v)
    i = 1
    for j in 1:d
        for k in i:(i+j-2)
            v[k] *= c
        end
        i += j
    end
    v
end
function MOI.get(m::SOItoMOIBridge{T}, ::MOI.ConstraintDual, cr::CI{F, DS}) where {T,F}
    scalevec!(_getattribute(m, cr, getdual), one(T)/2)
end

end # module
