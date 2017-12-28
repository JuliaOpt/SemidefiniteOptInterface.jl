const VIS = Union{VI, Vector{VI}}

function newblock(m::SOItoMOIBridge, n)
    push!(m.blockdims, n)
    m.nblocks += 1
end

isfree(m, v::VI) = v.value in m.free
isfree(m, v::Vector{VI}) = all(isfree.(m, v))
function unfree(m, v)
    @assert isfree(m, v)
    delete!(m.free, v.value)
end

function _constraintvariable!(m::SOItoMOIBridge{T}, vs::VIS, s::ZS) where T
    blk = newblock(m, -_length(vs))
    for (i, v) in _enumerate(vs)
        setvarmap!(m, v, (blk, i, i, one(T), _getconstant(m, s)))
        unfree(m, v)
    end
    blk
end
vscaling(::Type{<:NS}) = 1.
vscaling(::Type{<:PS}) = -1.
_length(vi::VI) = 1
_length(vi::Vector{VI}) = length(vi)
_enumerate(vi::VI) = enumerate((vi,))
_enumerate(vi::Vector{VI}) = enumerate(vi)
function _constraintvariable!(m::SOItoMOIBridge, vs::VIS, s::S) where S<:Union{NS, PS}
    blk = newblock(m, -_length(vs))
    for (i, v) in _enumerate(vs)
        setvarmap!(m, v, (blk, i, i, vscaling(S), _getconstant(m, s)))
        unfree(m, v)
    end
    blk
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
function _constraintvariable!(m::SOItoMOIBridge{T}, vs::VIS, ::DS) where T
    d = getmatdim(length(vs))
    k = 0
    blk = newblock(m, d)
    for i in 1:d
        for j in 1:i
            k += 1
            setvarmap!(m, vs[k], (blk, i, j, i == j ? one(T) : one(T)/2, zero(T)))
            unfree(m, vs[k])
        end
    end
    blk
end
_var(f::SVF) = f.variable
_var(f::VVF) = f.variables
function MOIU.allocateconstraint!(m::SOItoMOIBridge{T}, f::VF, s) where T
    vis = _var(f)
    fr = isfree(m, vis)
    if fr
        blk = _constraintvariable!(m, vis, s)
        if isa(s, ZS)
            ci = _allocateconstraint!(m, f, s)
            m.zeroblock[ci] = blk
            ci
        else
            CI{typeof(f), typeof(s)}(-blk)
        end
    else
        _allocateconstraint!(m, f, s)
    end
end


function loadfreevariables!(m::SOItoMOIBridge{T}) where T
    for vi in m.free
        blk = newblock(m, -2)
        # x free transformed into x = y - z with y, z >= 0
        setvarmap!(m, VI(vi), [(blk, 1, 1, one(T), zero(T)), (blk, 2, 2, -one(T), zero(T))])
    end
end
