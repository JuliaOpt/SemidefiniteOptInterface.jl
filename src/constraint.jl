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
                for (blk, i, j, coef, shift) in m.varmap[f.variables[i].value]
                    if !iszero(blk)
                        @assert !iszero(coef)
                        setconstraintcoefficient!(m.sdsolver, val*coef, c, blk, i, j)
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
            setconstraintconstant!(m.sdsolver, rhs[j], c)
        end
    end
end

_var(f::SVF, j) = f.variable
_var(f::VVF, j) = f.variables[j]
function loadconstraint!(m::SOItoMOIBridge{T}, cr::CR, f::VF, s::ZS) where T
    cs = m.constrmap[cr.value]
    for j in length(cs)
        vm = m.varmap[_var(f, j).value]
        @assert length(vm) == 1
        (blk, i, j, coef, shift) = first(vm)
        @assert coef == one(T)
        @assert s isa MOI.Zeros || shift == s.value
        c = cs[j]
        setconstraintcoefficient!(m.sdsolver, one(T), c, blk, i, j)
        setconstraintconstant!(m.sdsolver, zero(T), c)
    end
end
function loadconstraint!(m::SOItoMOIBridge, cr::CR, f::VF, s) end
function loadconstraint!(m::SOItoMOIBridge, cr::CR, af::AF, s::SupportedSets)
    cs = m.constrmap[cr.value]
    loadslacks!(m, cs)
    loadcoefficients!(m, cs, af, s)
end
