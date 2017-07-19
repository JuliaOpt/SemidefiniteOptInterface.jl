include("variable.jl")
include("constraint.jl")

function numberconstraint!(m::SOItoMOIBridge, ci, constr, ::Type)
    n = nconstraints(constr)
    m.constrmap[ci] = m.constr + (1:n)
    m.constr += n
end

for f in (:loadconstraint, :createslack, :numberconstraint)
    funs = Symbol(string(f) * "s!")
    fun = Symbol(string(f) * "!")
    @eval begin
        function $funs{FT, ST}(m::SOItoMOIBridge, constrs::Vector{Tuple{UInt64, FT}}, ::Type{FT}, ::Type{ST})
            for (ci, constr) in constrs
                $fun(m, ci, constr, ST)
            end
        end

        function $funs{FT}(m::SOItoMOIBridge, cs::SDConstraints{FT})
            $funs(m, cs.zeros, FT, MOI.Zeros)
            $funs(m, cs.nnegs, FT, MOI.Nonnegatives)
            $funs(m, cs.nposs, FT, MOI.Nonpositives)
            $funs(m, cs.psdcs, FT, MOI.PositiveSemidefiniteConeTriangle)
        end
    end
end

function loadobjective!(m::SOItoMOIBridge)
    obj = m.sdinstance.objective
    sgn = _objsgn(m)
    for (vr, val) in zip(obj.variables, obj.coefficients)
        vi = vr.value
        if !iszero(val)
            for (blk, i, j, coef) in m.varmap[vi]
                # in SDP format, it is max and in MPB Conic format it is min
                setobjectivecoefficient!(m.sdsolver, sgn * coef * val, blk, i, j)
            end
        end
    end
end

nconstraints(f::VAF) = length(f.constant)
nconstraints(f::SAF) = 1

nconstraints(ci::UInt64, f::MOI.AbstractFunction) = nconstraints(f)
nconstraints(c::Tuple) = nconstraints(c...)

nconstraints{FT<:VAF}(cs::Vector{Tuple{UInt64, FT}}) = sum(nconstraints.(cs))

nconstraints{FT<:SAF}(cs::Vector{Tuple{UInt64, FT}}) = length(cs)

function nconstraints(cs::SDConstraints)
    nconstraints(cs.zeros) + nconstraints(cs.nnegs) + nconstraints(cs.nposs) + nconstraints(cs.psdcs)
end

function init!(m::SOItoMOIBridge)
    m.nconstrs = nconstraints(m.sdinstance.sa) + nconstraints(m.sdinstance.va)
    m.varmap = Vector{Vector{Tuple{Int,Int,Int,Float64}}}(m.sdinstance.nvars)
    m.constrmap = Vector{UnitRange{Int}}(m.sdinstance.nconstrs)
    m.slackmap = Vector{Tuple{Int, Int, Int, Float64}}(m.nconstrs)
    m.constr = 0
    m.free = IntSet(1:m.sdinstance.nvars)
end

function loadprimal!(m::SOItoMOIBridge)
    init!(m)
    loadconstraints!(m, m.sdinstance.sv)
    loadconstraints!(m, m.sdinstance.vv)
    loadfreevariables!(m)
    numberconstraints!(m, m.sdinstance.sa)
    numberconstraints!(m, m.sdinstance.va)
    createslacks!(m, m.sdinstance.sa)
    createslacks!(m, m.sdinstance.va)
    initinstance!(m.sdsolver, m.blockdims, m.nconstrs)
    loadconstraints!(m, m.sdinstance.sa)
    loadconstraints!(m, m.sdinstance.va)
    loadobjective!(m)
end
