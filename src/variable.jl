const VIS = Union{VI, Vector{VI}}

function newblock(m::SOItoMOIBridge, n)
    push!(m.blockdims, n)
    m.nblocks += 1
end

isfree(m, v::VI) = v in m.free
isfree(m, v::Vector{VI}) = all(isfree.(m, v))
function unfree(m, v)
    @assert isfree(m, v)
    delete!(m.free, v)
end

function loadvariable!(m::SOItoMOIBridge{T}, vs::VIS, s::ZS) where T
    blk = newblock(m, -_length(vs))
    for (i, v) in _enumerate(vs)
        m.varmap[v] = [(blk, i, i, one(T), _getconstant(m, s))]
        unfree(m, v)
    end
end
vscaling(::Type{<:NS}) = 1
vscaling(::Type{<:PS}) = -1
_length(vi::VI) = 1
_length(vi::Vector{VI}) = length(vi)
_enumerate(vi::VI) = enumerate((vi,))
_enumerate(vi::Vector{VI}) = enumerate(vi)
function loadvariable!{S<:Union{NS, PS}}(m::SOItoMOIBridge, vs::VIS, s::S)
    blk = newblock(m, -_length(vs))
    for (i, v) in _enumerate(vs)
        m.varmap[v] = [(blk, i, i, vscaling(S), _getconstant(m, s))]
        unfree(m, v)
    end
end
function getmatdim(k::Integer)
    # n*(n+1)/2 = k
    # n^2+n-2k = 0
    # (-1 + sqrt(1 + 8k))/2
    n = div(isqrt(1 + 8k) - 1, 2)
    if n * (n+1) != 2*k
        error("sd dim not consistent")
    end
    n
end
function loadvariable!(m::SOItoMOIBridge{T}, vs::VIS, ::DS) where T
    d = getmatdim(length(vs))
    k = 0
    blk = newblock(m, d)
    for i in 1:d
        for j in 1:i
            k += 1
            m.varmap[vs[k]] = [(blk, i, j, i == j ? one(T) : one(T)/2, zero(T))]
            unfree(m, vs[k])
        end
    end
end
function loadvariable!(m::SOItoMOIBridge{T}, cr, constr::SVF, s) where T
    vi = constr.variable
    if isfree(m, vi)
        loadvariable!(m, vi, s)
    else
        push!(m.double, MOI.addconstraint!(m, MOI.ScalarAffineFunction([constr.variable], [one(T)], zero(T)), s))
    end
end
function loadvariable!(m::SOItoMOIBridge, cr, constr::VVF, s)
    vis = constr.variables
    if isfree(m, vis)
        loadvariable!(m, vis, s)
    else
        push!(m.double, MOI.addconstraint!(m, MOI.VectorAffineFunction{Float64}(constr), s))
    end
end

function loadvariable!(m::SOItoMOIBridge, cr, constr::AF, s) end

function loadfreevariables!(m::SOItoMOIBridge{T}) where T
    for vi in m.free
        blk = newblock(m, -2)
        # x free transformed into x = y - z with y, z >= 0
        m.varmap[vi] = [(blk, 1, 1, one(T), zero(T)), (blk, 2, 2, -one(T), zero(T))]
    end
end
