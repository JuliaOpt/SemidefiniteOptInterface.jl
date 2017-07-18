function createslack!(m::SOItoMOIBridge, cs, ::Type{<:ZS})
    m.slackmap[cs] = (0, 0, 0, 0.)
end
function createslack!(m::SOItoMOIBridge, cs, ::Type{T}) where T <: Union{NS, PS}
    blk = newblock(m, -length(cs))
    for (i, c) in enumerate(cs)
        m.slackmap[c] = (blk, i, i, vscaling(T))
    end
end
function createslack!(m::SOItoMOIBridge, cs, ::Type{DS})
    d = getmatdim(length(cs))
    k = 0
    blk = newblock(m, d)
    for i in 1:d
        for j in i:d
            k += 1
            m.slackmap[cs[k]] = (blk, i, j, i == j ? 1. : .5)
        end
    end
end

function createslack!{S}(m::SOItoMOIBridge, ci, constr::AF, ::Type{S})
    cs = m.constrmap[ci]
    createslack!(m, cs, S)
end

function loadslack!(m::SOItoMOIBridge, c::Integer)
    blk, i, j, coef = m.slackmap[c]
    if blk != 0
        @assert !iszero(coef)
        setconstraintcoefficient!(m.sdsolver, -coef, c, blk, i, j)
    end
end
function loadslacks!(m::SOItoMOIBridge, cs)
    for c in cs
        loadslack!(m, c)
    end
end

_row(f::SAF, i) = 1
_row(f::VAF, i) = f.outputindex[i]

function loadcoefficients!(m::SOItoMOIBridge, cs::UnitRange, f::AF)
    if !isempty(cs)
        for j in 1:nconstraints(f)
            c = cs[j]
            setconstraintconstant!(m.sdsolver, -f.constant[j], c)
        end
        for i in 1:length(f.variables)
            val = f.coefficients[i]
            if !iszero(val)
                c = cs[_row(f, i)]
                for (blk, i, j, coef) in m.varmap[f.variables[i].value]
                    @assert coef != 0
                    setconstraintcoefficient!(m.sdsolver, val*coef, c, blk, i, j)
                end
            end
        end
    end
end

function loadconstraint!{S}(m::SOItoMOIBridge, ci::UInt64, af::AF, ::Type{S})
    cs = m.constrmap[ci]
    loadslacks!(m, cs)
    loadcoefficients!(m, cs, af)
end
