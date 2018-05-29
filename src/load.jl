include("variable.jl")
include("constraint.jl")

MOIU.canallocate(::SOItoMOIBridge, ::MOI.ObjectiveSense) = true
function MOIU.allocate!(optimizer::SOItoMOIBridge, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    # To be sure that it is done before load!(optimizer, ::ObjectiveFunction, ...), we do it in allocate!
    optimizer.objsign = sense == MOI.MinSense ? -1 : 1
end
MOIU.canallocate(::SOItoMOIBridge{T}, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}) where T = true
function MOIU.allocate!(::SOItoMOIBridge, ::MOI.ObjectiveFunction, ::MOI.ScalarAffineFunction) end

MOIU.canload(m::SOItoMOIBridge, ::MOI.ObjectiveSense) = true
function MOIU.load!(::SOItoMOIBridge, ::MOI.ObjectiveSense, ::MOI.OptimizationSense) end
MOIU.canload(m::SOItoMOIBridge{T}, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}) where T = true
function MOIU.load!(optimizer::SOItoMOIBridge, ::MOI.ObjectiveFunction, f::MOI.ScalarAffineFunction)
    obj = MOIU.canonical(f)
    optimizer.objconstant = f.constant
    for t in obj.terms
        if !iszero(t.coefficient)
            for (blk, i, j, coef, shift) in varmap(optimizer, t.variable_index)
                if !iszero(blk)
                    # in SDP format, it is max and in MPB Conic format it is min
                    setobjectivecoefficient!(optimizer.sdoptimizer, optimizer.objsign * coef * t.coefficient, blk, i, j)
                end
                optimizer.objshift += t.coefficient * shift
            end
        end
    end
end

function MOIU.allocatevariables!(optimizer::SOItoMOIBridge{T}, nvars) where T
    optimizer.free = IntSet(1:nvars)
    optimizer.varmap = Vector{Vector{Tuple{Int, Int, Int, T, T}}}(nvars)
    VI.(1:nvars)
end

function MOIU.loadvariables!(optimizer::SOItoMOIBridge, nvars)
    @assert nvars == length(optimizer.varmap)
    loadfreevariables!(optimizer)
    init!(optimizer.sdoptimizer, optimizer.blockdims, optimizer.nconstrs)
end
