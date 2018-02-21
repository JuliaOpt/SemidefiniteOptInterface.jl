# Index
MOI.candelete(instance::SOItoMOIBridge, r::MOI.Index) = MOI.candelete(instance.sdmodel, r)
MOI.isvalid(instance::SOItoMOIBridge, r::MOI.Index) = MOI.isvalid(instance.sdmodel, r)
MOI.delete!(instance::SOItoMOIBridge, r::MOI.Index) = MOI.delete!(instance.sdmodel, r)

# Attributes
for f in (:canget, :canset, :set!, :get, :get!)
    @eval begin
        MOI.$f(instance::SOItoMOIBridge, attr::MOI.AnyAttribute) = MOI.$f(instance.sdmodel, attr)
        # Objective function
        MOI.$f(instance::SOItoMOIBridge, attr::MOI.AnyAttribute, arg::Union{MOI.OptimizationSense, MOI.AbstractScalarFunction}) = MOI.$f(instance.sdmodel, attr, arg)
    end
end
for f in (:canget, :canset)
    @eval begin
        MOI.$f(instance::SOItoMOIBridge, attr::MOI.AnyAttribute, index::Type{<:MOI.Index}) = MOI.$f(instance.sdmodel, attr, index)
    end
end
for f in (:set!, :get, :get!)
    @eval begin
        MOI.$f(instance::SOItoMOIBridge, attr::MOI.AnyAttribute, index::MOI.Index) = MOI.$f(instance.sdmodel, attr, index)
        MOI.$f(instance::SOItoMOIBridge, attr::MOI.AnyAttribute, indices::Vector{<:MOI.Index}) = MOI.$f(instance.sdmodel, attr, indices)
    end
end

# Constraints
MOI.canaddconstraint(instance::SOItoMOIBridge, f::MOI.AbstractFunction, s::MOI.AbstractSet) = MOI.canaddconstraint(instance.sdmodel, f, s)
MOI.addconstraint!(instance::SOItoMOIBridge, f::MOI.AbstractFunction, s::MOI.AbstractSet) = MOI.addconstraint!(instance.sdmodel, f, s)
MOI.canmodifyconstraint(instance::SOItoMOIBridge, cr::CI, change) = MOI.canmodifyconstraint(instance.sdmodel, cr, change)
MOI.modifyconstraint!(instance::SOItoMOIBridge, cr::CI, change) = MOI.modifyconstraint!(instance.sdmodel, cr, change)

# Objective
MOI.canmodifyobjective(instance::SOItoMOIBridge, change::MOI.AbstractFunctionModification) = MOI.canmodifyobjective(instance.sdmodel, change)
MOI.modifyobjective!(instance::SOItoMOIBridge, change::MOI.AbstractFunctionModification) = MOI.modifyobjective!(instance.sdmodel, change)

# Variables
MOI.addvariable!(instance::SOItoMOIBridge) = MOI.addvariable!(instance.sdmodel)
MOI.addvariables!(instance::SOItoMOIBridge, n) = MOI.addvariables!(instance.sdmodel, n)
