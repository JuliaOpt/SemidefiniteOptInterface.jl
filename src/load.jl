include("variable.jl")
include("constraint.jl")

function numberconstraint!(m::SOItoMOIBridge, cr, f, s)
    n = nconstraints(f, s)
    m.constrmap[cr.value] = m.constr + (1:n)
    m.constr += n
end

function bridgeconstraint!(m, cr, f, s::MOI.Interval)
    push!(m.int, SplitIntervalBridge(m, f, s))
    m.bridgemap[cr.value] = length(m.int)
end

function bridgeconstraint!(m, cr, f, s::MOI.PositiveSemidefiniteConeScaled)
    push!(m.psdcs, PSDCScaledBridge(m, f, s))
    m.bridgemap[cr.value] = length(m.psdcs)
end

function bridgeconstraint!(m, cr, f, s::MOI.SecondOrderCone)
    push!(m.soc, SOCtoPSDCBridge(m, f, s))
    m.bridgemap[cr.value] = length(m.soc)
end

function bridgeconstraint!(m, cr, f, s::MOI.RotatedSecondOrderCone)
    push!(m.rsoc, RSOCtoPSDCBridge(m, f, s))
    m.bridgemap[cr.value] = length(m.rsoc)
end

for (f, SS) in ((:loadvariable, :SupportedSets), (:loadconstraint, :SupportedSets), (:createslack, :SupportedSets), (:numberconstraint, :SupportedSets), (:bridgeconstraint, :BridgedSets))
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
    for (vr, val) in zip(obj.variables, obj.coefficients)
        vi = vr.value
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

nconstraints(cr::CR, f::MOI.AbstractFunction, s::MOI.AbstractSet) = nconstraints(f, s)
nconstraints(c::Tuple) = nconstraints(c...)

nconstraints(cs::Vector) = sum(nconstraints.(cs))

nconstraints{F<:SAF, S<:SupportedSets}(cs::Vector{MOIU.C{F, S}}) = length(cs)
nconstraints{F<:SVF, S<:ZS}(cs::Vector{MOIU.C{F, S}}) = length(cs)
nconstraints{F<:SVF, S<:SupportedSets}(cs::Vector{MOIU.C{F, S}}) = 0
nconstraints{F, S<:BridgedSets}(cs::Vector{MOIU.C{F, S}}) = 0

function resetbridges!(m)
    for s in m.int
        MOI.delete!(m, s)
    end
    for s in m.psdcs
        MOI.delete!(m, s)
    end
    for s in m.soc
        MOI.delete!(m, s)
    end
    for s in m.rsoc
        MOI.delete!(m, s)
    end
    for s in m.double
        MOI.delete!(m, s)
    end
    m.bridgemap = Vector{Int}(m.sdinstance.nextconstraintid)
    m.int = SplitIntervalBridge{Float64}[]
    m.psdcs = PSDCScaledBridge{Float64}[]
    m.soc = SOCtoPSDCBridge{Float64}[]
    m.rsoc = RSOCtoPSDCBridge{Float64}[]
    m.double = CR[]
end

function initvariables!(m::SOItoMOIBridge)
    m.objshift = 0.0
    m.constr = 0
    m.nblocks = 0
    m.blockdims = Int[]
    m.free = IntSet(1:m.sdinstance.nextvariableid)
    m.varmap = Vector{Vector{Tuple{Int,Int,Int,Float64,Float64}}}(m.sdinstance.nextvariableid)
end

function initconstraints!(m::SOItoMOIBridge)
    m.nconstrs = sum(MOIU.broadcastvcat(nconstraints, m.sdinstance))
    m.constrmap = Vector{UnitRange{Int}}(m.sdinstance.nextconstraintid)
    m.slackmap = Vector{Tuple{Int, Int, Int, Float64}}(m.nconstrs)
end

_broadcastcall(f, m) = MOIU.broadcastcall(constrs -> f(m, constrs), m.sdinstance)

function loadprimal!(m::SOItoMOIBridge)
    resetbridges!(m)
    _broadcastcall(bridgeconstraints!, m)
    initvariables!(m)
    _broadcastcall(loadvariables!, m)
    loadfreevariables!(m)
    initconstraints!(m)
    _broadcastcall(numberconstraints!, m)
    _broadcastcall(createslacks!, m)
    initinstance!(m.sdsolver, m.blockdims, m.nconstrs)
    _broadcastcall(loadconstraints!, m)
    loadobjective!(m)
end
