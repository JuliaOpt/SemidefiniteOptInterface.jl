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
    n = div(d * (d+1), 2)
    cr = MOI.addconstraint!(m, _SOCtoPSDCaff(f, d, n), MOI.PositiveSemidefiniteConeTriangle(d))
    SOCtoPSDCBridge(d, cr)
end
function _SOCtoPSDCaff(f::MOI.VectorOfVariables, dim, n)
    nv = length(f.variables)
    _SOCtoPSDCaff(MOI.VectorAffineFunction(collect(1:nv), f.variables, ones(nv), zeros(nv)), dim, n)
end
function _SOCtoPSDCaff{T}(f::MOI.VectorAffineFunction{T}, dim, n)
    constant = [f.constant; zeros(T, n - length(f.constant))]
    # Needs to add t*I
    tI = find(oi -> oi == 1, f.outputindex)
    tvars = f.variables[tI]
    tcoefs = f.coefficients[tI]
    N0 = length(f.outputindex)
    Ni = length(tI)
    N = N0 + (dim-1) * Ni
    outputindex  = Vector{Int}(N); outputindex[1:N0]  = f.outputindex
    variables    = Vector{VR}(N);  variables[1:N0]    = f.variables
    coefficients = Vector{T}(N);   coefficients[1:N0] = f.coefficients
    cur = N0
    j = dim-1
    i = 1+dim
    while i <= n
        outputindex[cur+(1:Ni)]  = i
        variables[cur+(1:Ni)]    = tvars
        coefficients[cur+(1:Ni)] = tcoefs
        i += j
        j -= 1
        cur += Ni
    end
    MOI.VectorAffineFunction(outputindex, variables, coefficients, constant)
end
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
function MOI.delete!(m, c::SOCtoPSDCBridge)
    MOI.delete!(m, c.cr)
end
