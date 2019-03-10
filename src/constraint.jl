nconstraints(f::Union{SVF, SAF}, s) = 1
nconstraints(f::VVF, s) = length(f.variables)

function _allocate_constraint(m::SOItoMOIBridge, f, s)
    ci = CI{typeof(f), typeof(s)}(m.nconstrs)
    n = nconstraints(f, s)
    # Fails on Julia v0.6
    #m.constrmap[ci] = m.nconstrs .+ (1:n)
    m.constrmap[ci] = (m.nconstrs + 1):(m.nconstrs + n)
    m.nconstrs += n
    return ci
end
function MOIU.allocate_constraint(m::SOItoMOIBridge, f::SAF, s::SupportedSets)
    _allocate_constraint(m::SOItoMOIBridge, f, s)
end

output_index(::MOI.ScalarAffineTerm)  = 1
output_index(t::MOI.VectorAffineTerm) = t.output_index
scalar_term(t::MOI.ScalarAffineTerm) = t
scalar_term(t::MOI.VectorAffineTerm) = t.scalar_term

_getconstant(m::SOItoMOIBridge, s::MOI.AbstractScalarSet) = MOIU.getconstant(s)
_getconstant(m::SOItoMOIBridge{T}, s::MOI.AbstractSet) where T = zero(T)

function loadcoefficients!(m::SOItoMOIBridge, cs::UnitRange, f::SAF, s)
    f = MOIU.canonical(f) # sum terms with same variables and same outputindex
    if !isempty(cs)
        rhs = _getconstant(m, s) .- MOI._constant(f)
        for t in f.terms
            st = scalar_term(t)
            if !iszero(st.coefficient)
                c = cs[output_index(t)]
                for (blk, i, j, coef, shift) in varmap(m, st.variable_index)
                    if !iszero(blk)
                        @assert !iszero(coef)
                        setconstraintcoefficient!(m.sdoptimizer, st.coefficient*coef, c, blk, i, j)
                    end
                    if isa(rhs, Vector)
                        rhs[output_index(t)] -= st.coefficient * shift
                    else
                        rhs -= st.coefficient * shift
                    end
                end
            end
        end
        for j in 1:MOI.output_dimension(f)
            c = cs[j]
            setconstraintconstant!(m.sdoptimizer, rhs[j], c)
        end
    end
end

function MOIU.load_constraint(m::SOItoMOIBridge, ci::CI, f::SAF, s::SupportedSets)
    setconstant!(m, ci, s)
    cs = m.constrmap[ci]
    @assert !isempty(cs)
    loadcoefficients!(m, cs, f, s)
end
