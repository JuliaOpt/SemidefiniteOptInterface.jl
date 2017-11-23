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
function PSDCScaledBridge{T}(instance, f, s::MOI.PositiveSemidefiniteConeScaled) where T
    dim = MOI.dimension(s)
    diagidx = IntSet()
    i = 0
    for j in 1:dim
        i += j
        push!(diagidx, i)
    end
    cr = MOI.addconstraint!(instance, unscalefunction(f, diagidx), MOI.PositiveSemidefiniteConeTriangle(dim))
    PSDCScaledBridge(dim, diagidx, cr)
end
MOI.get(::PSDCScaledBridge, ::MOI.NumberOfConstraints) = 0
MOI.get(::PSDCScaledBridge{T}, ::MOI.NumberOfConstraints{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}) where T = 1
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
MOI.canget(instance::MOI.AbstractInstance, a::MOI.ConstraintPrimal, c::PSDCScaledBridge) = true
function MOI.get(instance::MOI.AbstractInstance, a::MOI.ConstraintPrimal, c::PSDCScaledBridge)
    scalevec!(MOI.get(instance, MOI.ConstraintPrimal(), c.cr), sqrt(2))
end
MOI.canget(instance::MOI.AbstractInstance, a::MOI.ConstraintDual, c::PSDCScaledBridge) = true
function MOI.get(instance::MOI.AbstractInstance, a::MOI.ConstraintDual, c::PSDCScaledBridge)
    scalevec!(MOI.get(instance, MOI.ConstraintDual(), c.cr), sqrt(2))
end

function trimap(i, j)
    @assert j <= i
    div((i-1)*i, 2) + j
end

"""
    _SOCtoPSDCaff{T}(f::MOI.VectorAffineFunction{T}, g::MOI.ScalarAffineFunction{T})

Builds a VectorAffineFunction representing the upper (or lower) triangular part of the matrix
[ f[1]     f[2:end]' ]
[ f[2:end] g * I     ]
"""
function _SOCtoPSDCaff{T}(f::MOI.VectorAffineFunction{T}, g::MOI.ScalarAffineFunction{T})
    dim = length(f.constant)
    n = div(dim * (dim+1), 2)
    # Needs to add t*I
    N0 = length(f.variables)
    Ni = length(g.variables)
    N = N0 + (dim-1) * Ni
    outputindex  = Vector{Int}(N); outputindex[1:N0]  = trimap.(f.outputindex, 1)
    variables    = Vector{VR}(N);  variables[1:N0]    = f.variables
    coefficients = Vector{T}(N);   coefficients[1:N0] = f.coefficients
    constant = [f.constant; zeros(T, n - length(f.constant))]
    cur = N0
    for i in 2:dim
        k = trimap(i, i)
        outputindex[cur+(1:Ni)]  = k
        variables[cur+(1:Ni)]    = g.variables
        coefficients[cur+(1:Ni)] = g.coefficients
        constant[k] = g.constant
        cur += Ni
    end
    MOI.VectorAffineFunction(outputindex, variables, coefficients, constant)
end

# (t, x) is transformed into the matrix
# [t  x']
# [x t*I]
# Indeed by the Schur Complement, it is positive definite iff
# tI ≻ 0
# t - x' * (t*I)^(-1) * x ≻ 0
# which is equivalent to
# t > 0
# t^2 > x' * x
struct SOCtoPSDCBridge{T}
    dim::Int
    cr::CR{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}
end
function SOCtoPSDCBridge{T}(instance, f, s::MOI.SecondOrderCone) where T
    d = MOI.dimension(s)
    cr = MOI.addconstraint!(instance, _SOCtoPSDCaff(f), MOI.PositiveSemidefiniteConeTriangle(d))
    SOCtoPSDCBridge(d, cr)
end
MOI.get(::SOCtoPSDCBridge, ::MOI.NumberOfConstraints) = 0
MOI.get(::SOCtoPSDCBridge{T}, ::MOI.NumberOfConstraints{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}) where T = 1

_SOCtoPSDCaff(f::MOI.VectorOfVariables) = _SOCtoPSDCaff(tofun(f))
_SOCtoPSDCaff(f::MOI.VectorAffineFunction) = _SOCtoPSDCaff(f, MOIU.eachscalar(f)[1])

MOI.canget(instance::MOI.AbstractInstance, a::MOI.ConstraintPrimal, c::SOCtoPSDCBridge) = true
function MOI.get(instance::MOI.AbstractInstance, a::MOI.ConstraintPrimal, c::SOCtoPSDCBridge)
    MOI.get(instance, MOI.ConstraintPrimal(), c.cr)[trimap.(1:c.dim, 1)]
end
MOI.canget(instance::MOI.AbstractInstance, a::MOI.ConstraintDual, c::SOCtoPSDCBridge) = true
function MOI.get(instance::MOI.AbstractInstance, a::MOI.ConstraintDual, c::SOCtoPSDCBridge)
    dual = MOI.get(instance, MOI.ConstraintDual(), c.cr)
    tdual = sum(i -> dual[trimap(i, i)], 1:c.dim)
    [tdual; dual[trimap.(2:c.dim, 1)]*2]
end

# (t, u, x) is transformed into the matrix
# [t   x']
# [x 2u*I]
# Indeed by the Schur Complement, it is positive definite iff
# uI ≻ 0
# t - x' * (2u*I)^(-1) * x ≻ 0
# which is equivalent to
# u > 0
# 2t*u > x' * x
struct RSOCtoPSDCBridge{T}
    dim::Int
    cr::CR{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}
end
function RSOCtoPSDCBridge{T}(instance, f, s::MOI.RotatedSecondOrderCone) where T
    d = MOI.dimension(s)-1
    cr = MOI.addconstraint!(instance, _RSOCtoPSDCaff(f), MOI.PositiveSemidefiniteConeTriangle(d))
    RSOCtoPSDCBridge(d, cr)
end
MOI.get(::RSOCtoPSDCBridge, ::MOI.NumberOfConstraints) = 0
MOI.get(::RSOCtoPSDCBridge{T}, ::MOI.NumberOfConstraints{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}) where T = 1

_RSOCtoPSDCaff(f::MOI.VectorOfVariables) = _RSOCtoPSDCaff(tofun(f))
function _RSOCtoPSDCaff(f::MOI.VectorAffineFunction)
    n = length(f.constant)
    g = MOIU.eachscalar(f)[2]
    g = MOI.ScalarAffineFunction(g.variables, g.coefficients*2, g.constant*2)
    _SOCtoPSDCaff(MOIU.eachscalar(f)[[1; 3:n]], g)
end

MOI.canget(instance::MOI.AbstractInstance, a::MOI.ConstraintPrimal, c::RSOCtoPSDCBridge) = true
function MOI.get(instance::MOI.AbstractInstance, a::MOI.ConstraintPrimal, c::RSOCtoPSDCBridge)
    x = MOI.get(instance, MOI.ConstraintPrimal(), c.cr)[[trimap(1, 1); trimap(2, 2); trimap.(2:c.dim, 1)]]
    x[2] /= 2 # It is (2u*I)[1,1] so it needs to be divided by 2 to get u
    x
end
MOI.canget(instance::MOI.AbstractInstance, a::MOI.ConstraintDual, c::RSOCtoPSDCBridge) = true
function MOI.get(instance::MOI.AbstractInstance, a::MOI.ConstraintDual, c::RSOCtoPSDCBridge)
    dual = MOI.get(instance, MOI.ConstraintDual(), c.cr)
    udual = sum(i -> dual[trimap(i, i)], 2:c.dim)
    [dual[1]; 2udual; dual[trimap.(2:c.dim, 1)]*2]
end

function MOI.delete!(instance::MOI.AbstractInstance, c::Union{PSDCScaledBridge, SOCtoPSDCBridge, RSOCtoPSDCBridge})
    MOI.delete!(instance, c.cr)
end
