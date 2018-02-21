function createslack!(m::SOItoMOIBridge{T}, cs, ::ZS) where T
    m.slackmap[cs] = (0, 0, 0, zero(T))
end
function createslack!(m::SOItoMOIBridge, cs, ::S) where S <: Union{NS, PS}
    blk = newblock(m, -length(cs))
    for (i, c) in enumerate(cs)
        m.slackmap[c] = (blk, i, i, vscaling(S))
    end
end
function createslack!(m::SOItoMOIBridge{T}, cs, ::DS) where T
    d = getmatdim(length(cs))
    k = 0
    blk = newblock(m, d)
    for i in 1:d
        for j in 1:i
            k += 1
            m.slackmap[cs[k]] = (blk, i, j, i == j ? one(T) : one(T)/2)
        end
    end
end

function createslack!(m::SOItoMOIBridge, ci::CI, f, s)
    cs = m.constrmap[ci]
    createslack!(m, cs, s)
end

nconstraints(f::Union{SVF, SAF}, s) = 1
nconstraints(f::VVF, s) = length(f.variables)
nconstraints(f::VAF, s) = length(f.constant)

MOIU.canallocateconstraint(::SOItoMOIBridge{T}, f::Type{<:Union{VF, AF{T}}}, ::Type{<:SupportedSets}) where T = true
function _allocateconstraint!(m::SOItoMOIBridge, f, s)
    ci = CI{typeof(f), typeof(s)}(m.nconstrs)
    n = nconstraints(f, s)
    m.constrmap[ci] = m.nconstrs + (1:n)
    m.nconstrs += n
    resize!(m.slackmap, m.nconstrs)
    createslack!(m, ci, f, s)
    ci
end
function MOIU.allocateconstraint!(m::SOItoMOIBridge, f::AF, s::SupportedSets)
    _allocateconstraint!(m::SOItoMOIBridge, f, s)
end

function loadslack!(m::SOItoMOIBridge, c::Integer)
    blk, i, j, coef = m.slackmap[c]
    if blk != 0
        @assert !iszero(coef)
        setconstraintcoefficient!(m.sdoptimizer, -coef, c, blk, i, j)
    end
end
function loadslacks!(m::SOItoMOIBridge, cs)
    for c in cs
        loadslack!(m, c)
    end
end

_row(f::SAF, i) = 1
_row(f::VAF, i) = f.outputindex[i]

_getconstant(m::SOItoMOIBridge, s::MOI.AbstractScalarSet) = MOIU.getconstant(s)
_getconstant(m::SOItoMOIBridge{T}, s::MOI.AbstractSet) where T = zero(T)

function loadcoefficients!(m::SOItoMOIBridge, cs::UnitRange, f::AF, s)
    f = MOIU.canonical(f) # sum terms with same variables and same outputindex
    if !isempty(cs)
        rhs = _getconstant(m, s) - f.constant
        for i in 1:length(f.variables)
            val = f.coefficients[i]
            if !iszero(val)
                row = _row(f, i)
                c = cs[row]
                for (blk, i, j, coef, shift) in varmap(m, f.variables[i])
                    if !iszero(blk)
                        @assert !iszero(coef)
                        setconstraintcoefficient!(m.sdoptimizer, val*coef, c, blk, i, j)
                    end
                    if isa(rhs, Vector)
                        rhs[row] -= val * shift
                    else
                        rhs -= val * shift
                    end
                end
            end
        end
        for j in 1:length(f.constant)
            c = cs[j]
            setconstraintconstant!(m.sdoptimizer, rhs[j], c)
        end
    end
end

MOIU.canloadconstraint(::SOItoMOIBridge{T}, f::Type{<:Union{VF, AF{T}}}, ::Type{<:SupportedSets}) where T = true
function MOIU.loadconstraint!(m::SOItoMOIBridge, ci::CI, f::AF, s::SupportedSets)
    setconstant!(m, ci, s)
    cs = m.constrmap[ci]
    @assert !isempty(cs)
    loadslacks!(m, cs)
    loadcoefficients!(m, cs, f, s)
end
