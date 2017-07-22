include("variable.jl")
include("constraint.jl")

function numberconstraint!(m::SOItoMOIBridge, cr, f, s)
    n = nconstraints(f)
    m.constrmap[cr.value] = m.constr + (1:n)
    m.constr += n
end

for f in (:loadvariable, :loadconstraint, :createslack, :numberconstraint)
    funs = Symbol(string(f) * "s!")
    fun = Symbol(string(f) * "!")
    @eval begin
        function $funs(m::SOItoMOIBridge, constrs::Vector)
            for constr in constrs
                $fun(m, constr...)
            end
        end

        function $funs(m::SOItoMOIBridge, cs::SDInstanceScalarConstraints)
            $funs(m, cs.equalto)
            $funs(m, cs.greaterthan)
            $funs(m, cs.lessthan)
        end
        function $funs(m::SOItoMOIBridge, cs::SDInstanceVectorConstraints)
            $funs(m, cs.zeros)
            $funs(m, cs.nonnegatives)
            $funs(m, cs.nonpositives)
            $funs(m, cs.positivesemidefiniteconetriangle)
        end
    end
end

function loadobjective!(m::SOItoMOIBridge)
    obj = m.sdinstance.objective
    sgn = _objsgn(m)
    for (vr, val) in zip(obj.variables, obj.coefficients)
        vi = vr.value
        if !iszero(val)
            for (blk, i, j, coef, shift) in m.varmap[vi]
                if !iszero(blk)
                    # in SDP format, it is max and in MPB Conic format it is min
                    setobjectivecoefficient!(m.sdsolver, sgn * coef * val, blk, i, j)
                end
                m.objshift += coef * val * shift
            end
        end
    end
end

nconstraints(f::VAF) = length(f.constant)
nconstraints(f::SAF) = 1

nconstraints(cr::CR, f::MOI.AbstractFunction, s::MOI.AbstractSet) = nconstraints(f)
nconstraints(c::Tuple) = nconstraints(c...)

#nconstraints(cs::Vector) = sum(nconstraints.(cs))
function nconstraints(cs::Vector)
    s = 0
    for c in cs
        s += nconstraints(c)
    end
    s
end

nconstraints{F<:SAF, S}(cs::Vector{MOIU.C{F, S}}) = length(cs)
nconstraints{F<:VVF, S}(cs::Vector{MOIU.C{F, S}}) = 0

function nconstraints(cs::SDInstanceScalarConstraints)
    nconstraints(cs.equalto) + nconstraints(cs.greaterthan) + nconstraints(cs.lessthan)
end
function nconstraints(cs::SDInstanceVectorConstraints)
    nconstraints(cs.zeros) + nconstraints(cs.nonnegatives) + nconstraints(cs.nonpositives) + nconstraints(cs.positivesemidefiniteconetriangle)
end


function init!(m::SOItoMOIBridge)
    m.nconstrs = nconstraints(m.sdinstance.sa) + nconstraints(m.sdinstance.va)
    m.objshift = 0.0
    m.constr = 0
    m.nblocks = 0
    m.blockdims = Int[]
    m.free = IntSet(1:m.sdinstance.nvars)
    m.varmap = Vector{Vector{Tuple{Int,Int,Int,Float64,Float64}}}(m.sdinstance.nvars)
    m.constrmap = Vector{UnitRange{Int}}(m.sdinstance.nconstrs)
    m.slackmap = Vector{Tuple{Int, Int, Int, Float64}}(m.nconstrs)
end

function loadprimal!(m::SOItoMOIBridge)
    init!(m)
    loadvariables!(m, m.sdinstance.sv)
    loadvariables!(m, m.sdinstance.vv)
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
