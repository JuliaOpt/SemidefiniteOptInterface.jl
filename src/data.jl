# Index
MOI.candelete(instance::SOItoMOIBridge, r::MOI.Index) = MOI.candelete(instance.sdinstance, r)
MOI.isvalid(instance::SOItoMOIBridge, r::MOI.Index) = MOI.isvalid(instance.sdinstance, r)
MOI.delete!(instance::SOItoMOIBridge, r::MOI.Index) = MOI.delete!(instance.sdinstance, r)

# Attributes
for f in (:canget, :canset, :set!, :get, :get!)
    @eval begin
        MOI.$f(instance::SOItoMOIBridge, attr::MOI.AnyAttribute) = MOI.$f(instance.sdinstance, attr)
        # Objective function
        MOI.$f(instance::SOItoMOIBridge, attr::MOI.AnyAttribute, arg::Union{MOI.OptimizationSense, MOI.AbstractScalarFunction}) = MOI.$f(instance.sdinstance, attr, arg)
    end
end
for f in (:canget, :canset)
    @eval begin
        MOI.$f(instance::SOItoMOIBridge, attr::MOI.AnyAttribute, index::Type{<:MOI.Index}) = MOI.$f(instance.sdinstance, attr, index)
    end
end
for f in (:set!, :get, :get!)
    @eval begin
        MOI.$f(instance::SOItoMOIBridge, attr::MOI.AnyAttribute, index::MOI.Index) = MOI.$f(instance.sdinstance, attr, index)
        MOI.$f(instance::SOItoMOIBridge, attr::MOI.AnyAttribute, indices::Vector{<:MOI.Index}) = MOI.$f(instance.sdinstance, attr, indices)
    end
end

# Constraints
MOI.canaddconstraint(instance::SOItoMOIBridge, f::MOI.AbstractFunction, s::MOI.AbstractSet) = MOI.canaddconstraint(instance.sdinstance, f, s)
MOI.addconstraint!(instance::SOItoMOIBridge, f::MOI.AbstractFunction, s::MOI.AbstractSet) = MOI.addconstraint!(instance.sdinstance, f, s)
MOI.canmodifyconstraint(instance::SOItoMOIBridge, cr::CI, change) = MOI.canmodifyconstraint(instance.sdinstance, cr, change)
MOI.modifyconstraint!(instance::SOItoMOIBridge, cr::CI, change) = MOI.modifyconstraint!(instance.sdinstance, cr, change)

# Objective
MOI.canmodifyobjective(instance::SOItoMOIBridge, change::MOI.AbstractFunctionModification) = MOI.canmodifyobjective(instance.sdinstance, change)
MOI.modifyobjective!(instance::SOItoMOIBridge, change::MOI.AbstractFunctionModification) = MOI.modifyobjective!(instance.sdinstance, change)

# Variables
MOI.addvariable!(instance::SOItoMOIBridge) = MOI.addvariable!(instance.sdinstance)
MOI.addvariables!(instance::SOItoMOIBridge, n) = MOI.addvariables!(instance.sdinstance, n)
