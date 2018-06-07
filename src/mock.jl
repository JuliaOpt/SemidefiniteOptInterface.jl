mutable struct MockSDOptimizer{T} <: AbstractSDOptimizer
    nconstrs::Int
    blkdims::Vector{Int}
    constraint_constants::Vector{T}
    objective_coefficients::Vector{Tuple{T, Int, Int, Int}}
    constraint_coefficients::Vector{Vector{Tuple{T, Int, Int, Int}}}
end
MockSDOptimizer{T}() where T = MockSDOptimizer{T}(0, Int[], T[], Tuple{T, Int, Int, Int}[], Vector{Tuple{T, Int, Int, Int}}[])
mockSDoptimizer(T::Type) = SDOIOptimizer(MockSDOptimizer{T}(), T)
coefficienttype(::MockSDOptimizer{T}) where T = T

getnumberofconstraints(optimizer::MockSDOptimizer) = optimizer.nconstrs
getnumberofblocks(optimizer::MockSDOptimizer) = length(optimizer.blkdims)
getblockdimension(optimizer::MockSDOptimizer, blk) = optimizer.blkdims[blk]
function init!(optimizer::MockSDOptimizer{T}, blkdims::Vector{Int}, nconstrs::Integer) where T
    optimizer.nconstrs = nconstrs
    optimizer.blkdims = blkdims
    optimizer.constraint_constants = zeros(T, nconstrs)
    optimizer.objective_coefficients = Tuple{T, Int, Int, Int}[]
    optimizer.constraint_coefficients = map(i -> Tuple{T, Int, Int, Int}[], 1:nconstrs)
end

getconstraintconstant(optimizer::MockSDOptimizer, c) = optimizer.constraint_constants[c]
function setconstraintconstant!(optimizer::MockSDOptimizer, val, c::Integer)
    optimizer.constraint_constants[c] = val
end

getobjectivecoefficients(optimizer::MockSDOptimizer) = optimizer.objective_coefficients
function setobjectivecoefficient!(optimizer::MockSDOptimizer, val, blk::Integer, i::Integer, j::Integer)
    push!(optimizer.objective_coefficients, (val, blk, i, j))
end

getconstraintcoefficients(optimizer::MockSDOptimizer, c) = optimizer.constraint_coefficients[c]
function setconstraintcoefficient!(optimizer::MockSDOptimizer, val, c::Integer, blk::Integer, i::Integer, j::Integer)
    push!(optimizer.constraint_coefficients[c], (val, blk, i, j))
end
