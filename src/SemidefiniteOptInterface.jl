module SemidefiniteOptInterface

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

using Compat
using Compat.LinearAlgebra # for diag

abstract type AbstractSDOptimizer <: MOI.AbstractOptimizer end

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
const CI{F, S} = MOI.ConstraintIndex{F, S}

mutable struct SOItoMOIBridge{T, SIT <: AbstractSDOptimizer} <: MOI.AbstractOptimizer
    sdoptimizer::SIT
    setconstant::Dict{Int64, T}
    blkconstant::Dict{Int, T}
    objconstant::T
    objsign::Int
    objshift::T
    nconstrs::Int
    nblocks::Int
    blockdims::Vector{Int}
    free::BitSet
    varmap::Vector{Vector{Tuple{Int, Int, Int, T, T}}} # Variable Index vi -> blk, i, j, coef, shift # x = sum coef * block(X, blk)[i, j] + shift
    zeroblock::Dict{CI, Int}
    constrmap::Dict{CI, UnitRange{Int}} # Constraint Index ci -> cs
    slackmap::Vector{Tuple{Int, Int, Int, T}} # c -> blk, i, j, coef
    double::Vector{CI} # created when there are two cones for same variable
    function SOItoMOIBridge{T}(sdoptimizer::SIT) where {T, SIT}
        new{T, SIT}(sdoptimizer, Dict{Int64, T}(), Dict{Int, T}(),
            zero(T), 1, zero(T), 0, 0,
            Int[],
            BitSet(),
            Vector{Tuple{Int, Int, Int, T}}[],
            Dict{CI, Int}(),
            Dict{CI, UnitRange{Int}}(),
            Tuple{Int, Int, Int, T}[],
            CI[])
    end
end
varmap(optimizer::SOItoMOIBridge, vi::VI) = optimizer.varmap[vi.value]
function setvarmap!(optimizer::SOItoMOIBridge{T}, vi::VI, v::Tuple{Int, Int, Int, T, T}) where T
    setvarmap!(optimizer, vi, [v])
end
function setvarmap!(optimizer::SOItoMOIBridge{T}, vi::VI, vs::Vector{Tuple{Int, Int, Int, T, T}}) where T
    optimizer.varmap[vi.value] = vs
end

SDOIOptimizer(sdoptimizer::AbstractSDOptimizer, T=Float64) = SOItoMOIBridge{T}(sdoptimizer)

include("load.jl")

function MOI.is_empty(optimizer::SOItoMOIBridge)
    isempty(optimizer.double) &&
    isempty(optimizer.setconstant) &&
    isempty(optimizer.blkconstant) &&
    iszero(optimizer.objconstant) &&
    optimizer.objsign == 1 &&
    iszero(optimizer.objshift) &&
    iszero(optimizer.nconstrs) &&
    iszero(optimizer.nblocks) &&
    isempty(optimizer.blockdims) &&
    isempty(optimizer.free) &&
    isempty(optimizer.varmap) &&
    isempty(optimizer.zeroblock) &&
    isempty(optimizer.constrmap) &&
    isempty(optimizer.slackmap)
end
function MOI.empty!(optimizer::SOItoMOIBridge{T}) where T
    for s in optimizer.double
        MOI.delete(m, s)
    end
    MOI.empty!(optimizer.sdoptimizer)
    optimizer.double = CI[]
    optimizer.setconstant = Dict{Int64, T}()
    optimizer.blkconstant = Dict{Int, T}()
    optimizer.objconstant = zero(T)
    optimizer.objsign = 1
    optimizer.objshift = zero(T)
    optimizer.nconstrs = 0
    optimizer.nblocks = 0
    optimizer.blockdims = Int[]
    optimizer.free = BitSet()
    optimizer.varmap = Vector{Tuple{Int, Int, Int, T}}[]
    optimizer.zeroblock = Dict{CI, Int}()
    optimizer.constrmap = Dict{CI, UnitRange{Int}}()
    optimizer.slackmap = Tuple{Int, Int, Int, T}[]
end

function setconstant!(optimizer::SOItoMOIBridge, ci::CI, s) end
function setconstant!(optimizer::SOItoMOIBridge, ci::CI, s::MOI.AbstractScalarSet)
    optimizer.setconstant[ci.value] = MOIU.getconstant(s)
end
function addsetconstant(optimizer::SOItoMOIBridge, ci::CI{<:Any, <:MOI.AbstractScalarSet}, x)
    x + optimizer.setconstant[ci.value]
end
function addsetconstant(optimizer::SOItoMOIBridge, ci::CI, x)
    x
end
function addblkconstant(optimizer::SOItoMOIBridge, ci::CI{<:Any, <:Union{NS, PS}}, x)
    blk = -ci.value
    x .+ optimizer.blkconstant[blk]
end
function addblkconstant(optimizer::SOItoMOIBridge, ci::CI, x)
    x
end

function MOI.supports(optimizer::SOItoMOIBridge{T},
                      ::Union{MOI.ObjectiveSense,
                              MOI.ObjectiveFunction{<:Union{MOI.SingleVariable,
                                                            MOI.ScalarAffineFunction{T}}}}) where T
    return true
end

function MOI.supports_constraint(::SOItoMOIBridge{T},
                                 ::Type{<:Union{VF, AF{T}}},
                                 ::Type{<:SupportedSets}) where T
    return true
end

function MOI.copy_to(dest::SOItoMOIBridge, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end
MOIU.supports_allocate_load(::SOItoMOIBridge, copy_names::Bool) = !copy_names

MOI.optimize!(m::SOItoMOIBridge) = MOI.optimize!(m.sdoptimizer)

# Objective

function MOI.get(m::SOItoMOIBridge, ::MOI.ObjectiveValue)
    m.objshift + m.objsign * getprimalobjectivevalue(m.sdoptimizer) + m.objconstant
end

# Attributes

const SolverStatus = Union{MOI.TerminationStatus, MOI.PrimalStatus, MOI.DualStatus}
MOI.get(m::SOItoMOIBridge, s::SolverStatus) = MOI.get(m.sdoptimizer, s)


MOI.get(m::SOItoMOIBridge, ::MOI.ResultCount) = 1

function _getblock(M, blk::Int, s::Type{<:Union{NS, ZS}})
    diag(block(M, blk))
end
function _getblock(M, blk::Int, s::Type{<:PS})
    -diag(block(M, blk))
end
# Vectorized length for matrix dimension d
sympackedlen(d::Integer) = (d*(d+1)) >> 1
function _getblock(M::AbstractMatrix{T}, blk::Int, s::Type{<:DS}) where T
    B = block(M, blk)
    d = Compat.LinearAlgebra.checksquare(B)
    n = sympackedlen(d)
    v = Vector{T}(undef, n)
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

getvarprimal(m::SOItoMOIBridge, blk::Int, S) = getblock(getX(m.sdoptimizer), blk, S)
function getvardual(m::SOItoMOIBridge, blk::Int, S)
    getblock(getZ(m.sdoptimizer), blk, S)
end

function MOI.get(m::SOItoMOIBridge{T}, ::MOI.VariablePrimal, vi::VI) where T
    X = getX(m.sdoptimizer)
    x = zero(T)
    for (blk, i, j, coef, shift) in varmap(m, vi)
        x += shift
        if blk != 0
            x += block(X, blk)[i, j] * sign(coef)
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
    X = getX(m.sdoptimizer)
    blk, i, j, coef = m.slackmap[c]
    if iszero(blk)
        zero(T)
    else
        if i != j
            coef *= 2 # We should take block(X, blk)[i, j] + block(X, blk)[j, i] but they are equal
        end
        coef * block(X, blk)[i, j]
    end
end

function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, ci::CI{F, S}) where {F, S}
    if ci.value >= 0
        addsetconstant(m, ci, _getattribute(m, ci, getslack))
    else
        # Variable Function-in-S with S different from Zeros and EqualTo and not a double variable constraint
        blk = -ci.value
        addblkconstant(m, ci, getvarprimal(m, blk, S))
    end
end

function getvardual(m::SOItoMOIBridge{T}, vi::VI) where T
    Z = getZ(m.sdoptimizer)
    z = zero(T)
    for (blk, i, j, coef) in varmap(m, vi)
        if blk != 0
            z += block(Z, blk)[i, j] * sign(coef)
        end
    end
    z
end
getvardual(m::SOItoMOIBridge, f::SVF) = getvardual(m, f.variable)
getvardual(m::SOItoMOIBridge, f::VVF) = map(vi -> getvardual(m, vi), f.variables)
#function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, ci::CI{<:VF, S})
#    _getattribute(m, ci, getdual) + getvardual(m, MOI.get(m, MOI.ConstraintFunction(), ci))
#end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, ci::CI{<:VF, S}) where S<:SupportedSets
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
        -gety(m.sdoptimizer)[c]
    end
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, ci::CI)
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
function MOI.get(m::SOItoMOIBridge{T}, ::MOI.ConstraintDual,
                 ci::CI{<:AF{T}, DS}) where T
    scalevec!(_getattribute(m, ci, getdual), one(T)/2)
end

include("sdpa.jl")
include("mock.jl")

end # module
