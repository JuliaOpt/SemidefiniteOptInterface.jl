const VIS = Union{UInt64, Vector{UInt64}}

function newblock(m::SOItoMOIBridge, n)
    push!(m.blockdims, n)
    m.nblocks += 1
end

function loadvariable!(m::SOItoMOIBridge, vs::VIS, s::ZS)
    blk = newblock(m, -length(vs))
    for (i, v) in enumerate(vs)
        m.varmap[v] = [(blk, i, i, 1.0, _getconstant(m, s))]
        pop!(m.free, v)
    end
end
vscaling(::Type{<:NS}) = 1.
vscaling(::Type{<:PS}) = -1.
function loadvariable!{S<:Union{NS, PS}}(m::SOItoMOIBridge, vs::VIS, s::S)
    blk = newblock(m, -length(vs))
    for (i, v) in enumerate(vs)
        m.varmap[v] = [(blk, i, i, vscaling(S), _getconstant(m, s))]
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
    n
end
function loadvariable!(m::SOItoMOIBridge, vs::VIS, ::DS)
    d = getmatdim(length(vs))
    k = 0
    blk = newblock(m, d)
    for i in 1:d
        for j in i:d
            k += 1
            m.varmap[vs[k]] = [(blk, i, j, 1.0, 0.0)]
            pop!(m.free, vs[k])
        end
    end
end
function loadvariable!(m::SOItoMOIBridge, cr, constr::SVF, s)
    loadvariable!(m, constr.variable.value, s)
end
function loadvariable!(m::SOItoMOIBridge, cr, constr::VVF, s)
    loadvariable!(m, map(v -> v.value, constr.variables), s)
end

function loadvariable!(m::SOItoMOIBridge, cr, constr::AF, s) end

function loadfreevariables!(m::SOItoMOIBridge)
    for vi in m.free
        blk = newblock(m, -2)
        # x free transformed into x = y - z with y, z >= 0
        m.varmap[vi] = [(blk, 1, 1, 1., 0.0), (blk, 2, 2, -1., 0.0)]
    end
end
