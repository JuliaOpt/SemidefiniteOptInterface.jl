module SemidefiniteOptInterface

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

using LinearAlgebra # for diag

abstract type AbstractSDOptimizer <: MOI.AbstractOptimizer end

include("interface.jl")

const SAF{T} = MOI.ScalarAffineFunction{T}

const SupportedSets = Union{MOI.Nonnegatives, MOI.PositiveSemidefiniteConeTriangle}

const CI{F, S} = MOI.ConstraintIndex{F, S}
const AFFEQ{T} = MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.EqualTo{T}}

mutable struct SOItoMOIBridge{T, SIT <: AbstractSDOptimizer} <: MOI.AbstractOptimizer
    sdoptimizer::SIT
    b::Vector{T}
    objconstant::T
    objsign::Int
    blockdims::Vector{Int}
    varmap::Vector{Tuple{Int, Int, Int}} # Variable Index vi -> blk, i, j
    function SOItoMOIBridge{T}(sdoptimizer::SIT) where {T, SIT}
        new{T, SIT}(sdoptimizer, T[], zero(T), 1,
                    Int[], Vector{Tuple{Int, Int, Int, T}}[])
    end
end
varmap(optimizer::SOItoMOIBridge, vi::MOI.VariableIndex) = optimizer.varmap[vi.value]
function setvarmap!(optimizer::SOItoMOIBridge{T}, vi::MOI.VariableIndex, v::Tuple{Int, Int, Int}) where T
    setvarmap!(optimizer, vi, v)
end

SDOIOptimizer(sdoptimizer::AbstractSDOptimizer, T=Float64) = SOItoMOIBridge{T}(sdoptimizer)

function MOIU.allocate(optimizer::SOItoMOIBridge, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    # To be sure that it is done before load(optimizer, ::ObjectiveFunction, ...), we do it in allocate
    optimizer.objsign = sense == MOI.MIN_SENSE ? -1 : 1
end
function MOIU.allocate(::SOItoMOIBridge, ::MOI.ObjectiveFunction, ::Union{MOI.SingleVariable, MOI.ScalarAffineFunction}) end

function MOIU.load(::SOItoMOIBridge, ::MOI.ObjectiveSense, ::MOI.OptimizationSense) end
# Loads objective coefficient α * vi
function load_objective_term!(optimizer::SOItoMOIBridge, α, vi::MOI.VariableIndex)
    blk, i, j = varmap(optimizer, vi)
    coef = optimizer.objsign * α
    if i != j
        coef /= 2
    end
    # in SDP format, it is max and in MPB Conic format it is min
    setobjectivecoefficient!(optimizer.sdoptimizer, coef, blk, i, j)
end
function MOIU.load(optimizer::SOItoMOIBridge, ::MOI.ObjectiveFunction, f::MOI.ScalarAffineFunction)
    obj = MOIU.canonical(f)
    optimizer.objconstant = f.constant
    for t in obj.terms
        if !iszero(t.coefficient)
            load_objective_term!(optimizer, t.coefficient, t.variable_index)
        end
    end
end
function MOIU.load(optimizer::SOItoMOIBridge{T}, ::MOI.ObjectiveFunction, f::MOI.SingleVariable) where T
    load_objective_term!(optimizer, one(T), f.variable)
end

function new_block(optimizer::SOItoMOIBridge, set::MOI.Nonnegatives)
    push!(optimizer.blockdims, -MOI.dimension(set))
    blk = length(optimizer.blockdims)
    for i in 1:MOI.dimension(set)
        push!(optimizer.varmap, (blk, i, i))
    end
end

function new_block(optimizer::SOItoMOIBridge, set::MOI.PositiveSemidefiniteConeTriangle)
    push!(optimizer.blockdims, set.side_dimension)
    blk = length(optimizer.blockdims)
    for i in 1:set.side_dimension
        for j in 1:i
            push!(optimizer.varmap, (blk, i, j))
        end
    end
end

function MOIU.allocate_constrained_variables(
    optimizer::SOItoMOIBridge,
    set::Union{MOI.Nonnegatives, MOI.PositiveSemidefiniteConeTriangle}
)
    offset = length(optimizer.varmap)
    new_block(optimizer, set)
    ci = MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(set)}(offset + 1)
    return [MOI.VariableIndex(i) for i in offset .+ (1:MOI.dimension(set))], ci
end

function MOIU.load_constrained_variables(
    optimizer::SOItoMOIBridge, vis::Vector{MOI.VariableIndex},
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables},
    set::Union{MOI.Nonnegatives, MOI.PositiveSemidefiniteConeTriangle})
end

function MOIU.load_variables(optimizer::SOItoMOIBridge, nvars)
    @assert nvars == length(optimizer.varmap)
    init!(optimizer.sdoptimizer, optimizer.blockdims, length(optimizer.b))
end

function MOIU.allocate_constraint(optimizer::SOItoMOIBridge{T},
                                  func::MOI.ScalarAffineFunction{T},
                                  set::MOI.EqualTo{T}) where T
    push!(optimizer.b, MOI.constant(set))
    return AFFEQ{T}(length(optimizer.b))
end

function MOIU.load_constraint(m::SOItoMOIBridge, ci::AFFEQ,
                              f::MOI.ScalarAffineFunction, s::MOI.EqualTo)
    f = MOIU.canonical(f) # sum terms with same variables and same outputindex
    for t in f.terms
        if !iszero(t.coefficient)
            blk, i, j = varmap(m, t.variable_index)
            coef = t.coefficient
            if i != j
                coef /= 2
            end
            setconstraintcoefficient!(m.sdoptimizer, coef, ci.value, blk, i, j)
        end
    end
    setconstraintconstant!(m.sdoptimizer, MOI.constant(s) - MOI.constant(f), ci.value)
end

function MOI.supports(optimizer::SOItoMOIBridge, attr::MOI.AbstractOptimizerAttribute)
    return MOI.supports(optimizer.sdoptimizer, attr)
end
function MOI.get(optimizer::SOItoMOIBridge, attr::MOI.AbstractOptimizerAttribute)
    return MOI.get(optimizer.sdoptimizer, attr)
end
function MOI.set(optimizer::SOItoMOIBridge,
                 attr::MOI.AbstractOptimizerAttribute, value)
    return MOI.set(optimizer.sdoptimizer, attr, value)
end

function MOI.is_empty(optimizer::SOItoMOIBridge)
    isempty(optimizer.b) &&
    iszero(optimizer.objconstant) &&
    optimizer.objsign == 1 &&
    isempty(optimizer.blockdims) &&
    isempty(optimizer.varmap)
end
function MOI.empty!(optimizer::SOItoMOIBridge{T}) where T
    MOI.empty!(optimizer.sdoptimizer)
    optimizer.b = T[]
    optimizer.objconstant = zero(T)
    optimizer.objsign = 1
    optimizer.blockdims = Int[]
    optimizer.varmap = Tuple{Int, Int, Int}[]
end

function block(optimizer::SOItoMOIBridge, ci::MOI.ConstraintIndex{MOI.VectorOfVariables})
    return optimizer.varmap[ci.value][1]
end
function dimension(optimizer::SOItoMOIBridge, ci::MOI.ConstraintIndex{MOI.VectorOfVariables})
    blockdim = optimizer.blockdims[block(optimizer, ci)]
    if blockdim < 0
        return -blockdim
    else
        return MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(blockdim))
    end
end

function MOI.supports(
    optimizer::SOItoMOIBridge{T},
    ::Union{MOI.ObjectiveSense,
            MOI.ObjectiveFunction{<:Union{MOI.SingleVariable,
                                          MOI.ScalarAffineFunction{T}}}}) where T
    return true
end

function MOI.supports_constraint(
    ::SOItoMOIBridge, ::Type{MOI.VectorOfVariables},
    ::Type{<:Union{MOI.Nonnegatives,
                   MOI.PositiveSemidefiniteConeTriangle}})
    return true
end
function MOI.supports_constraint(
    ::SOItoMOIBridge{T}, ::Type{MOI.ScalarAffineFunction{T}},
    ::Type{MOI.EqualTo{T}}) where T
    return true
end

function MOI.copy_to(dest::SOItoMOIBridge, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end
MOIU.supports_allocate_load(::SOItoMOIBridge, copy_names::Bool) = !copy_names

MOI.optimize!(m::SOItoMOIBridge) = MOI.optimize!(m.sdoptimizer)

# Objective

function MOI.get(m::SOItoMOIBridge, attr::Union{MOI.ObjectiveValue, MOI.DualObjectiveValue})
    return m.objsign * MOI.get(m.sdoptimizer, attr) + m.objconstant
end

# Attributes
function MOI.get(m::SOItoMOIBridge,
                 attr::Union{MOI.RawStatusString, MOI.TerminationStatus,
                             MOI.PrimalStatus, MOI.DualStatus, MOI.SolveTime})
    return MOI.get(m.sdoptimizer, attr)
end

MOI.get(m::SOItoMOIBridge, ::MOI.ResultCount) = 1

function vectorize_block(M, blk::Integer, s::Type{MOI.Nonnegatives})
    return diag(block(M, blk))
end
function vectorize_block(M::AbstractMatrix{T}, blk::Integer, s::Type{MOI.PositiveSemidefiniteConeTriangle}) where T
    B = block(M, blk)
    d = LinearAlgebra.checksquare(B)
    n = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(d))
    v = Vector{T}(undef, n)
    k = 0
    for j in 1:d
        for i in 1:j
            k += 1
            v[k] = B[i, j]
        end
    end
    @assert k == n
    return v
end

function MOI.get(m::SOItoMOIBridge{T}, ::MOI.VariablePrimal, vi::MOI.VariableIndex) where T
    X = getX(m.sdoptimizer)
    blk, i, j = varmap(m, vi)
    return block(X, blk)[i, j]
end

function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.VectorOfVariables, S}) where S<:SupportedSets
    return vectorize_block(getX(m.sdoptimizer), block(m, ci), S)
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintPrimal, ci::AFFEQ)
    return m.b[ci.value]
end

function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual,
                 ci::CI{MOI.VectorOfVariables, S}) where S<:SupportedSets
    return vectorize_block(getZ(m.sdoptimizer), block(m, ci), S)
end
function MOI.get(optimizer::SOItoMOIBridge, ::MOI.ConstraintDual, ci::AFFEQ)
    return -gety(optimizer.sdoptimizer)[ci.value]
end

include("sdpa.jl")
include("mock.jl")

end # module
