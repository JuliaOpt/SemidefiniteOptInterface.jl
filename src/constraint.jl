function createslack!(m::SOItoMOIBridge, cs, ::ZS)
    m.slackmap[cs] = (0, 0, 0, 0.)
end
function createslack!(m::SOItoMOIBridge, cs, ::S) where S <: Union{NS, PS}
    blk = newblock(m, -length(cs))
    for (i, c) in enumerate(cs)
        m.slackmap[c] = (blk, i, i, vscaling(S))
    end
end
function createslack!(m::SOItoMOIBridge, cs, ::DS)
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

function createslack!(m::SOItoMOIBridge, cr::CR, constr::VF, s) end

function createslack!(m::SOItoMOIBridge, cr::CR, constr::AF, s)
    cs = m.constrmap[cr.value]
    createslack!(m, cs, s)
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

_getconstant(m::SOItoMOIBridge, s::MOI.AbstractScalarSet) = MOIU.getconstant(s)
_getconstant(m::SOItoMOIBridge, s::MOI.AbstractSet) = 0.0

function loadcoefficients!(m::SOItoMOIBridge, cs::UnitRange, f::AF, s)
    if !isempty(cs)
        rhs = _getconstant(m, s) - f.constant
        for i in 1:length(f.variables)
            val = f.coefficients[i]
            if !iszero(val)
                row = _row(f, i)
                c = cs[row]
                for (blk, i, j, coef, shift) in m.varmap[f.variables[i].value]
                    if !iszero(blk)
                        @assert !iszero(coef)
                        setconstraintcoefficient!(m.sdsolver, val*coef, c, blk, i, j)
                    end
                    if isa(rhs, Vector)
                        rhs[row] -= coef * shift
                    else
                        rhs -= coef * shift
                    end
                end
            end
        end
        for j in 1:length(f.constant)
            c = cs[j]
            setconstraintconstant!(m.sdsolver, rhs[j], c)
        end
    end
end

function loadconstraint!(m::SOItoMOIBridge, cr::CR, af::VF, s) end
function loadconstraint!(m::SOItoMOIBridge, cr::CR, af::AF, s)
    cs = m.constrmap[cr.value]
    loadslacks!(m, cs)
    loadcoefficients!(m, cs, af, s)
end
