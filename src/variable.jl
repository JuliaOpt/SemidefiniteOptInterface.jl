const VIS = Union{UInt64, Vector{UInt64}}

function newblock(m::SOItoMOIBridge, n)
    push!(m.blockdims, n)
    m.nblocks += 1
end

function loadconstraint!(m::SOItoMOIBridge, vs::VIS, ::Type{<:ZS})
    for vi in vs
        m.varmap[vi] = []
        pop!(m.free, vi)
    end
end
vscaling(::Type{<:NS}) = 1.
vscaling(::Type{<:PS}) = -1.
function loadconstraint!(m::SOItoMOIBridge, vs::VIS, ::Type{T}) where T <: Union{NS, PS}
    blk = newblock(m, -length(vs))
    for (i, v) in enumerate(vs)
        m.varmap[v] = [(blk, i, i, vscaling(T))]
        pop!(m.free, v)
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
end
function loadconstraint!(m::SOItoMOIBridge, vs::VIS, ::Type{DS})
    d = getmatdim(length(vs))
    k = 0
    blk = newblock(m, d)
    for i in 1:d
        for j in i:d
            k += 1
            m.varmap[vs[k]] = [(blk, i, j, 1.0)]
            pop!(m.free, vs[k])
        end
    end
end
function loadconstraint!{S}(m::SOItoMOIBridge, ci, constr::SVF, ::Type{S})
    loadconstraint!(m, constr.variable.value, S)
end
function loadconstraint!{S}(m::SOItoMOIBridge, ci, constr::VVF, ::Type{S})
    loadconstraint!(m, map(v -> v.value, constr.variables), S)
end

function loadfreevariables!(m::SOItoMOIBridge)
    for vi in m.free
        blk = newblock(m, -2)
        # x free transformed into x = y - z with y, z >= 0
        m.varmap[vi] = [(blk, 1, 1, 1.), (blk, 2, 2, -1.)]
    end
end
