include("variable.jl")
include("constraint.jl")

function numberconstraint!(m::SOItoMOIBridge, ci, f, s)
    n = nconstraints(f, s)
    m.constrmap[ci] = m.constr + (1:n)
    m.constr += n
end

for (f, SS) in ((:loadvariable, :SupportedSets), (:loadconstraint, :SupportedSets), (:createslack, :SupportedSets), (:numberconstraint, :SupportedSets))
    funs = Symbol(string(f) * "s!")
    fun = Symbol(string(f) * "!")
    @eval begin
        function $funs(m::SOItoMOIBridge, constrs::Vector) end
        function $funs{F, S<:$SS}(m::SOItoMOIBridge, constrs::Vector{MOIU.C{F, S}})
            for constr in constrs
                $fun(m, constr...)
            end
        end
    end
end

function loadobjective!(m::SOItoMOIBridge)
    obj = MOIU.canonical(m.sdinstance.objective)
    sgn = _objsgn(m)
    for (vi, val) in zip(obj.variables, obj.coefficients)
        if !iszero(val)
            for (blk, i, j, coef, shift) in m.varmap[vi]
                if !iszero(blk)
                    # in SDP format, it is max and in MPB Conic format it is min
                    setobjectivecoefficient!(m.sdsolver, sgn * coef * val, blk, i, j)
                end
                m.objshift += val * shift
            end
        end
    end
end

nconstraints(f::SVF, s::MOI.EqualTo) = 1
nconstraints(f::VVF, s::MOI.Zeros) = length(f.variables)
nconstraints(f::VF, s) = 0
nconstraints(f::VAF, s) = length(f.constant)
nconstraints(f::SAF, s) = 1

nconstraints(cr::CI, f::MOI.AbstractFunction, s::MOI.AbstractSet) = nconstraints(f, s)
nconstraints(c::Tuple) = nconstraints(c...)

nconstraints(cs::Vector) = sum(nconstraints.(cs))

nconstraints{F<:SAF, S<:SupportedSets}(cs::Vector{MOIU.C{F, S}}) = length(cs)
nconstraints{F<:SVF, S<:ZS}(cs::Vector{MOIU.C{F, S}}) = length(cs)
nconstraints{F<:SVF, S<:SupportedSets}(cs::Vector{MOIU.C{F, S}}) = 0

function initconstraints!(m::SOItoMOIBridge)
    m.nconstrs = sum(MOIU.broadcastvcat(nconstraints, m.sdinstance))
    m.slackmap = Vector{Tuple{Int, Int, Int, Float64}}(m.nconstrs)
end

_broadcastcall(f, m) = MOIU.broadcastcall(constrs -> f(m, constrs), m.sdinstance)

function _empty!(m::SOItoMOIBridge{T}) where T
    for s in m.double
        MOI.delete!(m, s)
    end
    m.double = CI[]
    m.objshift = zero(T)
    m.constr = 0
    m.nblocks = 0
    m.blockdims = Int[]
    m.free = Set{VI}()
    m.varmap = Dict{VI, Tuple{Int,Int,Int,T,T}}()
    m.constrmap = Dict{CI, UnitRange{Int}}()
end

function _allocatevariables!(m::SOItoMOIBridge, vis::Vector{VI})
    for vi in vis
        push!(m.free, vi)
    end
end

function loadprimal!(m::SOItoMOIBridge)
    _empty!(m)
    _allocatevariables!(m, MOI.get(m, MOI.ListOfVariableIndices()))
    _broadcastcall(loadvariables!, m)
    loadfreevariables!(m)
    initconstraints!(m)
    _broadcastcall(numberconstraints!, m)
    _broadcastcall(createslacks!, m)
    initinstance!(m.sdsolver, m.blockdims, m.nconstrs)
    _broadcastcall(loadconstraints!, m)
    loadobjective!(m)
end
