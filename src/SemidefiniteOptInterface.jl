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

const VI = MOI.VariableIndex
const CI{FT, ST} = MOI.ConstraintIndex{FT, ST}

mutable struct SOItoMOIBridge{T, SIT <: AbstractSDSolverInstance} <: MOI.AbstractSolverInstance
    sdinstance::SDInstance{T} # Will be removed when
    idxmap::MOIU.IndexMap     # InstanceManager is ready
    sdsolver::SIT
    objsign::Int
    objshift::T
    nconstrs::Int
    nblocks::Int
    blockdims::Vector{Int}
    free::IntSet
    varmap::Vector{Vector{Tuple{Int, Int, Int, T, T}}} # Variable Index vi -> blk, i, j, coef, shift # x = sum coef * X[blk][i, j] + shift
    zeroblock::Dict{CI, Int}
    constrmap::Dict{CI, UnitRange{Int}} # Constraint Index ci -> cs
    slackmap::Vector{Tuple{Int, Int, Int, T}} # c -> blk, i, j, coef
    double::Vector{CI} # created when there are two cones for same variable
    function SOItoMOIBridge{T}(sdsolver::SIT) where {T, SIT}
        new{T, SIT}(SDInstance{T}(), MOIU.IndexMap(), sdsolver,
            1, zero(T), 0, 0,
            Int[],
            IntSet(),
            Vector{Tuple{Int, Int, Int, T}}[],
            Dict{CI, Int}(),
            Dict{CI, UnitRange{Int}}(),
            Tuple{Int, Int, Int, T}[],
            CI[])
    end
end
varmap(instance::SOItoMOIBridge, vi::VI) = instance.varmap[vi.value]
function setvarmap!(instance::SOItoMOIBridge{T}, vi::VI, v::Tuple{Int, Int, Int, T, T}) where T
    setvarmap!(instance, vi, [v])
end
function setvarmap!(instance::SOItoMOIBridge{T}, vi::VI, vs::Vector{Tuple{Int, Int, Int, T, T}}) where T
    instance.varmap[vi.value] = vs
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

function MOI.empty!(instance::SOItoMOIBridge{T}) where T
    for s in instance.double
        MOI.delete!(m, s)
    end
    instance.double = CI[]
    instance.objsign = 1
    instance.objshift = zero(T)
    instance.nconstrs = 0
    instance.nblocks = 0
    instance.blockdims = Int[]
    instance.free = IntSet()
    instance.varmap = Vector{Tuple{Int, Int, Int, T}}[]
    instance.zeroblock = Dict{CI, Int}()
    instance.constrmap = Dict{CI, UnitRange{Int}}()
    instance.slackmap = Tuple{Int, Int, Int, T}[]
end

MOI.copy!(dest::SOItoMOIBridge, src::MOI.AbstractInstance) = MOIU.allocateload!(dest, src)

# Constraints

function MOI.optimize!(m::SOItoMOIBridge)
    res = MOI.copy!(m, m.sdinstance)
    @assert res.status == MOI.CopySuccess
    m.idxmap = res.indexmap
    MOI.optimize!(m.sdsolver)
end

# Objective

MOI.canget(m::SOItoMOIBridge, ::MOI.ObjectiveValue) = true
function MOI.get(m::SOItoMOIBridge, ::MOI.ObjectiveValue)
    m.objshift + m.objsign * getprimalobjectivevalue(m.sdsolver) + m.sdinstance.objective.constant
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
                                      MOI.ConstraintDual}, ::Type{<:MOI.Index}) = true

function _getblock(M, blk::Int, s::Type{<:Union{NS, ZS}})
    diag(M[blk])
end
function _getblock(M, blk::Int, s::Type{<:PS})
    -diag(M[blk])
end
# Vectorized length for matrix dimension d
sympackedlen(d::Integer) = (d*(d+1)) >> 1
function _getblock(M::AbstractMatrix{T}, blk::Int, s::Type{<:DS}) where T
    B = M[blk]
    d = Base.LinAlg.checksquare(B)
    n = sympackedlen(d)
    v = Vector{T}(n)
    k = 0
    for j in 1:d
        for i in 1:j
            k += 1
            v[k] = B[i, j]
        end
    end
    @assert k == n
    v
end
function getblock(M, blk::Int, s::Type{<:MOI.AbstractScalarSet})
    vd = _getblock(M, blk, s)
    @assert length(vd) == 1
    vd[1]
end
function getblock(M, blk::Int, s::Type{<:MOI.AbstractVectorSet})
    _getblock(M, blk, s)
end

getvarprimal(m::SOItoMOIBridge, blk::Int, S) = getblock(getX(m.sdsolver), blk, S)
getvardual(m::SOItoMOIBridge, blk::Int, S) = getblock(getZ(m.sdsolver), blk, S)

function MOI.get(m::SOItoMOIBridge{T}, ::MOI.VariablePrimal, vi::VI) where T
    vi = m.idxmap[vi]
    X = getX(m.sdsolver)
    x = zero(T)
    for (blk, i, j, coef, shift) in varmap(m, vi)
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

function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, ci::CI)
    constant = _getconstant(m, MOI.get(m, MOI.ConstraintSet(), ci))
    ci = m.idxmap[ci]
    _getattribute(m, ci, getslack) + constant
end
# These constraints do not create any constraint or slack, it is just variables
#_getvarprimal(m, sv::SVF) = MOI.get(m, MOI.VariablePrimal(), sv.variable)
#_getvarprimal(m, vv::VVF) = MOI.get(m, MOI.VariablePrimal(), vv.variables)
function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, ci::CI{<:VF, S}) where S <: Union{NS, PS, DS}
    ci = m.idxmap[ci]
    if ci.value >= 0
        blk = m.zeroblock[ci] # TODO double
    else
        blk = -ci.value
    end
    getvarprimal(m, blk, S)
    #_getvarprimal(m, MOI.get(m, MOI.ConstraintFunction(), ci))
end

function getvardual(m::SOItoMOIBridge{T}, vi::VI) where T
    Z = getZ(m.sdsolver)
    z = zero(T)
    for (blk, i, j, coef) in varmap(m, vi)
        if blk != 0
            z += Z[blk][i, j] * sign(coef)
        end
    end
    z
end
getvardual(m::SOItoMOIBridge, f::SVF) = getvardual(m, f.variable)
getvardual(m::SOItoMOIBridge, f::VVF) = map(vi -> getvardual(m, vi), f.variables)
#function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, ci::CI{<:VF, S})
#    ci = m.idxmap[ci]
#    _getattribute(m, ci, getdual) + getvardual(m, MOI.get(m, MOI.ConstraintFunction(), ci))
#end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, ci::CI{<:VF, S}) where S<:SupportedSets
    ci = m.idxmap[ci]
    if ci.value < 0
        getvardual(m, -ci.value, S)
    else
        dual = _getattribute(m, ci, getdual)
        if haskey(m.zeroblock, ci) # ZS
            dual + getvardual(m, m.zeroblock[ci], S)
        else # var constraint on unfree constraint
            dual
        end
    end
end

function getdual(m::SOItoMOIBridge{T}, c::Int) where T
    if c == 0
        zero(T)
    else
        -gety(m.sdsolver)[c]
    end
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, ci::CI)
    ci = m.idxmap[ci]
    _getattribute(m, ci, getdual)
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
function MOI.get(m::SOItoMOIBridge{T}, ::MOI.ConstraintDual, ci::CI{F, DS}) where {T,F}
    ci = m.idxmap[ci]
    scalevec!(_getattribute(m, ci, getdual), one(T)/2)
end

end # module
