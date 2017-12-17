include("variable.jl")
include("constraint.jl")

for f in (:loadvariable, :loadconstraint, :allocateconstraint)
    funs = Symbol(string(f) * "s!")
    fun = Symbol(string(f) * "!")
    @eval begin
        function $funs(m::SOItoMOIBridge, constrs::Vector) end
        function $funs{F, S<:SupportedSets}(m::SOItoMOIBridge, constrs::Vector{MOIU.C{F, S}})
            for constr in constrs
                $fun(m, constr...)
            end
        end
    end
end

_broadcastcall(f, m) = MOIU.broadcastcall(constrs -> f(m, constrs), m.sdinstance)

function _loadobjective!(m::SOItoMOIBridge, sense::MOI.OptimizationSense, f::MOI.ScalarAffineFunction)
    obj = MOIU.canonical(f)
    m.objsign = sense == MOI.MinSense ? -1 : 1
    for (vi, val) in zip(obj.variables, obj.coefficients)
        if !iszero(val)
            for (blk, i, j, coef, shift) in m.varmap[vi]
                if !iszero(blk)
                    # in SDP format, it is max and in MPB Conic format it is min
                    setobjectivecoefficient!(m.sdsolver, m.objsign * coef * val, blk, i, j)
                end
                m.objshift += val * shift
            end
        end
    end
end

function _empty!(m::SOItoMOIBridge{T}) where T
    for s in m.double
        MOI.delete!(m, s)
    end
    m.double = CI[]
    m.objsign = 1
    m.objshift = zero(T)
    m.nconstrs = 0
    m.nblocks = 0
    m.blockdims = Int[]
    m.free = Set{VI}()
    m.varmap = Dict{VI, Tuple{Int,Int,Int,T,T}}()
    m.constrmap = Dict{CI, UnitRange{Int}}()
    m.slackmap = Tuple{Int, Int, Int, T}[]
end

function _allocatevariables!(m::SOItoMOIBridge, vis::Vector{VI})
    for vi in vis
        push!(m.free, vi)
    end
end

function _setup!(m::SOItoMOIBridge)
    loadfreevariables!(m)
    initinstance!(m.sdsolver, m.blockdims, m.nconstrs)
end

function loadprimal!(m::SOItoMOIBridge)
    _empty!(m)
    _allocatevariables!(m, MOI.get(m, MOI.ListOfVariableIndices()))
    _broadcastcall(loadvariables!, m)
    _broadcastcall(allocateconstraints!, m)
    _setup!(m)
    _broadcastcall(loadconstraints!, m)
    _loadobjective!(m, MOI.get(m, MOI.ObjectiveSense()), MOI.get(m, MOI.ObjectiveFunction()))
end
