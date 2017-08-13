# l ≤ ⟨a, x⟩ + α ≤ u
# is transformed into
# ⟨a, x⟩ + α ≥ l
# ⟨a, x⟩ + α ≤ u
struct SplitIntervalBridge{T}
    lower::CR{MOI.ScalarAffineFunction{T}, MOI.GreaterThan{T}}
    upper::CR{MOI.ScalarAffineFunction{T}, MOI.LessThan{T}}
end
function SplitIntervalBridge(m, f, s::MOI.Interval)
    lower = MOI.addconstraint!(m, f, MOI.GreaterThan(s.lower))
    upper = MOI.addconstraint!(m, f, MOI.LessThan(s.upper))
    SplitIntervalBridge(lower, upper)
end
function MOI.getattribute(m, a::MOI.ConstraintPrimal, c::SplitIntervalBridge)
    # lower and upper should give the same value
    MOI.getattribute(m, MOI.ConstraintPrimal(), c.lower)
end
function MOI.getattribute(m, a::MOI.ConstraintDual, c::SplitIntervalBridge)
    lowd = MOI.getattribute(m, MOI.ConstraintDual(), c.lower) # Should be nonnegative
    uppd = MOI.getattribute(m, MOI.ConstraintDual(), c.upper) # Should be nonpositive
    if lowd > -uppd
        lowd
    else
        uppd
    end
end
function MOI.delete!(m, c::SplitIntervalBridge)
    MOI.delete!(m, c.lower)
    MOI.delete!(m, c.upper)
end

function tofun(f::MOI.VectorOfVariables)
    nv = length(f.variables)
    MOI.VectorAffineFunction(collect(1:nv), f.variables, ones(nv), zeros(nv))
end

# [x x x]
# [x x x] in PSD cone scaled
# [x x x]
# is transformed into
# [x    x/√2 x/√2]
# [x/√2 x    x/√2] in PSD cone triangle
# [x/√2 x/√2 x   ]
struct PSDCScaledBridge{T}
    dim::Int
    diagidx::IntSet
    cr::CR{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}
end
unscalefunction(f::MOI.VectorOfVariables, diagidx) = unscalefunction(tofun(f), diagidx)
function unscalefunction(f::MOI.VectorAffineFunction, diagidx)
    outputindex = f.outputindex
    variables = f.variables
    coefficients = copy(f.coefficients)
    constant = copy(f.constant)
    s2 = sqrt(2)
    scalevec!(constant, 1/s2)
    for i in eachindex(outputindex)
        if !(outputindex[i] in diagidx)
            coefficients[i] /= s2
        end
    end
    g = MOI.VectorAffineFunction(outputindex, variables, coefficients, constant)
end
function PSDCScaledBridge(m, f, s::MOI.PositiveSemidefiniteConeScaled)
    dim = MOI.dimension(s)
    diagidx = IntSet()
    i = 1
    j = dim
    while j > 0
        push!(diagidx, i)
        i += j
        j -= 1
    end
    cr = MOI.addconstraint!(m, unscalefunction(f, diagidx), MOI.PositiveSemidefiniteConeTriangle(dim))
    PSDCScaledBridge(dim, diagidx, cr)
end
function scalevec!(v, c)
    d = div(isqrt(1+8length(v))-1, 2)
    @assert div(d*(d+1), 2) == length(v)
    i = 1
    for j in d:-1:1
        for k in (i+1):(i+j-1)
            v[k] *= c
        end
        i += j
    end
    v
end
function MOI.getattribute(m, a::MOI.ConstraintPrimal, c::PSDCScaledBridge)
    scalevec!(MOI.getattribute(m, MOI.ConstraintPrimal(), c.cr), sqrt(2))
end
function MOI.getattribute(m, a::MOI.ConstraintDual, c::PSDCScaledBridge)
    scalevec!(MOI.getattribute(m, MOI.ConstraintDual(), c.cr), sqrt(2))
end

"""
    _SOCtoPSDCaff{T}(f::MOI.VectorAffineFunction{T}, g::MOI.ScalarAffineFunction{T})

Builds a VectorAffineFunction representing the upper (or lower) triangular part of the matrix
[ f[1]     f[2:end]^T ]
[ f[2:end] g * I      ]
"""
function _SOCtoPSDCaff{T}(f::MOI.VectorAffineFunction{T}, g::MOI.ScalarAffineFunction{T})
    dim = length(f.constant)
    n = div(dim * (dim+1), 2)
    # Needs to add t*I
    N0 = length(f.variables)
    Ni = length(g.variables)
    N = N0 + (dim-1) * Ni
    outputindex  = Vector{Int}(N); outputindex[1:N0]  = f.outputindex
    variables    = Vector{VR}(N);  variables[1:N0]    = f.variables
    coefficients = Vector{T}(N);   coefficients[1:N0] = f.coefficients
    constant = [f.constant; zeros(T, n - length(f.constant))]
    cur = N0
    j = dim
    i = dim + 1
    while i <= n
        outputindex[cur+(1:Ni)]  = i
        variables[cur+(1:Ni)]    = g.variables
        coefficients[cur+(1:Ni)] = g.coefficients
        constant[i] = g.constant
        j -= 1
        i += j
        cur += Ni
    end
    MOI.VectorAffineFunction(outputindex, variables, coefficients, constant)
end

function Base.getindex(f::MOI.VectorAffineFunction, i)
    I = find(oi -> oi == i, f.outputindex)
    MOI.ScalarAffineFunction(f.variables[I], f.coefficients[I], f.constant[i])
end

# (t, x) is transformed into the matrix
# [t x^T]
# [x t*I]
# Indeed by the Schur Complement, it is positive definite iff
# tI ≻ 0
# t - x^T * (t*I)^(-1) * x ≻ 0
# which is equivalent to
# t > 0
# t^2 > x^T * x
struct SOCtoPSDCBridge{T}
    dim::Int
    cr::CR{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}
end
function SOCtoPSDCBridge(m, f, s::MOI.SecondOrderCone)
    d = MOI.dimension(s)
    cr = MOI.addconstraint!(m, _SOCtoPSDCaff(f), MOI.PositiveSemidefiniteConeTriangle(d))
    SOCtoPSDCBridge(d, cr)
end

_SOCtoPSDCaff(f::MOI.VectorOfVariables) = _SOCtoPSDCaff(tofun(f))
_SOCtoPSDCaff(f::MOI.VectorAffineFunction) = _SOCtoPSDCaff(f, f[1])

function MOI.getattribute(m, a::MOI.ConstraintPrimal, c::SOCtoPSDCBridge)
    MOI.getattribute(m, MOI.ConstraintPrimal(), c.cr)[1:c.dim]
end
function MOI.getattribute(m, a::MOI.ConstraintDual, c::SOCtoPSDCBridge)
    dual = MOI.getattribute(m, MOI.ConstraintDual(), c.cr)
    tdual = 0.0
    j = c.dim
    i = 1
    n = div(c.dim * (c.dim+1), 2)
    while i <= n
        tdual += dual[i]
        i += j
        j -= 1
    end
    [tdual; dual[2:c.dim]]
end

# (t, u, x) is transformed into the matrix
# [t  x^T ]
# [x u*I/2]
# Indeed by the Schur Complement, it is positive definite iff
# uI ≻ 0
# t - x^T * (u*I/2)^(-1) * x ≻ 0
# which is equivalent to
# u > 0
# 2t*u > x^T * x
struct RSOCtoPSDCBridge{T}
    dim::Int
    cr::CR{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}
end
function RSOCtoPSDCBridge(m, f, s::MOI.RotatedSecondOrderCone)
    d = MOI.dimension(s)-1
    cr = MOI.addconstraint!(m, _RSOCtoPSDCaff(f), MOI.PositiveSemidefiniteConeTriangle(d))
    RSOCtoPSDCBridge(d, cr)
end

_RSOCtoPSDCaff(f::MOI.VectorOfVariables) = _RSOCtoPSDCaff(tofun(f))
function _RSOCtoPSDCaff(f::MOI.VectorAffineFunction)
    n = length(f.constant)
    g = f[2]
    g = MOI.ScalarAffineFunction(g.variables, g.coefficients*2, g.constant*2)
    _SOCtoPSDCaff(f[[1; 3:n]], g)
end

function MOI.getattribute(m, a::MOI.ConstraintPrimal, c::RSOCtoPSDCBridge)
    MOI.getattribute(m, MOI.ConstraintPrimal(), c.cr)[1:c.dim]
end
function MOI.getattribute(m, a::MOI.ConstraintDual, c::RSOCtoPSDCBridge)
    dual = MOI.getattribute(m, MOI.ConstraintDual(), c.cr)
    udual = 0.0
    j = c.dim
    i = c.dim + 1
    n = div(c.dim * (c.dim+1), 2)
    while i <= n
        udual += dual[i]
        j -= 1
        i += j
    end
    [dual[1]; 2udual; dual[2:c.dim]]
end

function MOI.delete!(m, c::Union{PSDCScaledBridge, SOCtoPSDCBridge, RSOCtoPSDCBridge})
    MOI.delete!(m, c.cr)
end
