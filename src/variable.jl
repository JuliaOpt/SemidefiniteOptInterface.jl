const VIS = Union{UInt64, Vector{UInt64}}

function newblock(m::SOItoMOIBridge, n)
    push!(m.blockdims, n)
    m.nblocks += 1
end

isfree(m, v::UInt64) = v in m.free
isfree(m, v::Vector{UInt64}) = all(isfree.(m, v))
function unfree(m, v)
    @assert isfree(m, v)
    delete!(m.free, v)
end

function loadvariable!(m::SOItoMOIBridge{T}, vs::VIS, s::ZS) where T
    blk = newblock(m, -length(vs))
    for (i, v) in enumerate(vs)
        m.varmap[v] = [(blk, i, i, one(T), _getconstant(m, s))]
        unfree(m, v)
    end
end
vscaling(::Type{<:NS}) = 1
vscaling(::Type{<:PS}) = -1
function loadvariable!{S<:Union{NS, PS}}(m::SOItoMOIBridge, vs::VIS, s::S)
    blk = newblock(m, -length(vs))
    for (i, v) in enumerate(vs)
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
    vs = constr.variable.value
    if isfree(m, vs)
        loadvariable!(m, vs, s)
    else
        push!(m.double, MOI.addconstraint!(m, MOI.ScalarAffineFunction([constr.variable], [one(T)], zero(T)), s))
    end
end
function loadvariable!(m::SOItoMOIBridge, cr, constr::VVF, s)
    vs = map(v -> v.value, constr.variables)
    if isfree(m, vs)
        loadvariable!(m, vs, s)
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
