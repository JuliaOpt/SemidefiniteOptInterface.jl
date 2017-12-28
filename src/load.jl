include("variable.jl")
include("constraint.jl")

function MOIU.allocateobjective!(m::SOItoMOIBridge, sense::MOI.OptimizationSense, f::MOI.ScalarAffineFunction) end

function MOIU.loadobjective!(m::SOItoMOIBridge, sense::MOI.OptimizationSense, f::MOI.ScalarAffineFunction)
    obj = MOIU.canonical(f)
    m.objsign = sense == MOI.MinSense ? -1 : 1
    for (vi, val) in zip(obj.variables, obj.coefficients)
        if !iszero(val)
            for (blk, i, j, coef, shift) in varmap(m, vi)
                if !iszero(blk)
                    # in SDP format, it is max and in MPB Conic format it is min
                    setobjectivecoefficient!(m.sdsolver, m.objsign * coef * val, blk, i, j)
                end
                m.objshift += val * shift
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
