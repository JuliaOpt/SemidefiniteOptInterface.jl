include("variable.jl")
include("constraint.jl")

MOIU.canallocate(::SOItoMOIBridge, ::MOI.ObjectiveSense) = true
function MOIU.allocate!(instance::SOItoMOIBridge, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    # To be sure that it is done before load!(instance, ::ObjectiveFunction, ...), we do it in allocate!
    instance.objsign = sense == MOI.MinSense ? -1 : 1
end
MOIU.canallocate(::SOItoMOIBridge, ::MOI.ObjectiveFunction) = true
function MOIU.allocate!(::SOItoMOIBridge, ::MOI.ObjectiveFunction, ::MOI.ScalarAffineFunction) end

MOIU.canload(m::SOItoMOIBridge, ::MOI.ObjectiveSense) = true
function MOIU.load!(::SOItoMOIBridge, ::MOI.ObjectiveSense, ::MOI.OptimizationSense) end
MOIU.canload(m::SOItoMOIBridge, ::MOI.ObjectiveFunction) = true
function MOIU.load!(instance::SOItoMOIBridge, ::MOI.ObjectiveFunction, f::MOI.ScalarAffineFunction)
    obj = MOIU.canonical(f)
    for (vi, val) in zip(obj.variables, obj.coefficients)
        if !iszero(val)
            for (blk, i, j, coef, shift) in varmap(instance, vi)
                if !iszero(blk)
                    # in SDP format, it is max and in MPB Conic format it is min
                    setobjectivecoefficient!(instance.sdsolver, instance.objsign * coef * val, blk, i, j)
                end
                instance.objshift += val * shift
            end
        end
    end
end

function MOIU.allocatevariables!(instance::SOItoMOIBridge{T}, nvars) where T
    instance.free = IntSet(1:nvars)
    instance.varmap = Vector{Vector{Tuple{Int, Int, Int, T, T}}}(nvars)
    VI.(1:nvars)
end

function MOIU.loadvariables!(instance::SOItoMOIBridge, nvars)
    @assert nvars == length(instance.varmap)
    loadfreevariables!(instance)
    initinstance!(instance.sdsolver, instance.blockdims, instance.nconstrs)
end
