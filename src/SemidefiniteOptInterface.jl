module SemidefiniteOptInterface

using MathOptInterface
const MOI = MathOptInterface

abstract type AbstractSDSolver <: MOI.AbstractSolver end

MOI.getattribute(m::AbstractSDSolver, ::MOI.SupportsDuals) = true

abstract type AbstractSDSolverInstance <: MOI.AbstractSolverInstance end

include("interface.jl")

const SVF = MOI.ScalarVariablewiseFunction
const VVF = MOI.VectorVariablewiseFunction
const VF  = Union{SVF, VVF}
const SAF{T} = MOI.ScalarAffineFunction{T}
const VAF{T} = MOI.VectorAffineFunction{T}
const AF{T}  = Union{SAF{T}, VAF{T}}
const ASF{T} = Union{SVF, SAF{T}}
const AVF{T} = Union{VVF, VAF{T}}

const CR{FT, ST} = MOI.ConstraintReference{FT, ST}

include("sdinstance.jl")

mutable struct SOItoMOIBridge{ST <: AbstractSDSolver, SIT <: AbstractSDSolverInstance} <: MOI.AbstractSolverInstance
    solver::ST
    sdinstance::SDInstance{Float64}
    sdsolver::SIT
    nconstrs::Int
    constr::Int
    nblocks::Int
    blockdims::Vector{Int}
    free::IntSet
    varmap::Vector{Vector{Tuple{Int, Int, Int, Float64}}} # Variable Reference value vi -> blk, i, j, coef
    constrmap::Vector{UnitRange{Int}} # Constraint Reference value ci -> cs
    slackmap::Vector{Tuple{Int, Int, Int, Float64}} # c -> blk, i, j, coef
    function SOItoMOIBridge(solver::ST, sdsolver::SIT) where {ST, SIT}
        new{ST, SIT}(solver, SDInstance{Float64}(), sdsolver,
            0, 0, 0,
            Int[],
            IntSet(),
            Vector{Tuple{Int,Int,Int,Float64}}[],
            UnitRange{Int}[],
            Tuple{Int, Int, Int, Float64}[])
    end
    function SOItoMOIBridge(solver::ST) where ST
        SOItoMOIBridge(solver, SDSolverInstance(solver))
    end
end

MOI.SolverInstance(s::AbstractSDSolver) = SOItoMOIBridge(s)

MOI.supportsproblem(m::AbstractSDSolver, ::Type{<:MOI.AbstractFunction}, constrtypes) = false

function _supports(x::DataType, y::DataType)
    false
end
function _supports{FT<:Union{VF, AF}, ST<:Union{ZS, NS, PS, DS}}(::Type{FT}, ::Type{ST})
    true
end
MOI.supportsproblem(m::AbstractSDSolver, ::Type{<:SAF}, constrtypes) = reduce((b, ct) -> b & _supports(ct...), true, constrtypes)

include("load.jl")

# Variables

MOI.addvariable!(m::SOItoMOIBridge) = MOI.addvariable!(m.sdinstance)
MOI.addvariables!(m::SOItoMOIBridge, n::Integer) = MOI.addvariables!(m.sdinstance, n)

# Constraints

MOI.delete!(m::SOItoMOIBridge, cr::CR) = MOI.delete!(m.sdinstance, cr)
MOI.addconstraint!(m::SOItoMOIBridge, f, s) = MOI.addconstraint!(m.sdinstance, f, s)

MOI.cangetattribute(m::SOItoMOIBridge, a::Union{MOI.ConstraintFunction, MOI.ConstraintSet}, ref) = MOI.cangetattribute(m.sdinstance, a, ref)
MOI.getattribute(m::SOItoMOIBridge, a::Union{MOI.ConstraintFunction, MOI.ConstraintSet}, ref) = MOI.getattribute(m.sdinstance, a, ref)

function MOI.getattribute(m::SOItoMOIBridge, a::Union{MOI.ListOfConstraints,
                                                      MOI.NumberOfConstraints,
                                                      MOI.NumberOfVariables,
                                                      MOI.Sense})
    MOI.getattribute(m.sdinstance, a)
end

function MOI.optimize!(m::SOItoMOIBridge)
    loadprimal!(m)
    MOI.optimize!(m.sdsolver)
end

function MOI.modifyconstraint!(m::SOItoMOIBridge, cr::CR, change)
    MOI.modifyconstraint!(m.sdinstance, cr, change)
end

# Objective

function MOI.setobjective!(m::SOItoMOIBridge, sense::MOI.OptimizationSense, f)
    m.sdinstance.sense = sense
    m.sdinstance.objective = f
end

_objsgn(m) = m.sdinstance.sense == MOI.MinSense ? -1 : 1
function MOI.getattribute(m, ::MOI.ObjectiveValue)
    _objsgn(m) * getprimalobjectivevalue(m.sdsolver) + m.sdinstance.objective.constant
end
function MOI.modifyobjective!(m, change::MOI.AbstractFunctionModification)
    m.sdinstance.objective = MOI.modifyfunction(m.sdinstance.objective, change)
end

# Attributes

MOI.cangetattribute(m::SOItoMOIBridge, s::Union{MOI.PrimalStatus,
                                                MOI.DualStatus}) = MOI.cangetattribute(m.sdsolver, s)

MOI.cangetattribute(m::SOItoMOIBridge, ::Union{MOI.TerminationStatus,
                                               MOI.ObjectiveValue,
                                               MOI.NumberOfVariables,
                                               MOI.NumberOfConstraints,
                                               MOI.ListOfConstraints}) = true

MOI.cangetattribute(m::SOItoMOIBridge, ::Union{MOI.VariablePrimal,
                                               MOI.ConstraintPrimal,
                                               MOI.ConstraintDual}, ref) = true

MOI.getattribute(m::SOItoMOIBridge, s::Union{MOI.TerminationStatus,
                                             MOI.PrimalStatus,
                                             MOI.DualStatus}) = MOI.getattribute(m.sdsolver, s)

function MOI.getattribute(m::SOItoMOIBridge, ::MOI.VariablePrimal, vr::MOI.VariableReference)
    X = getX(m.sdsolver)
    x = 0.0
    for (blk, i, j, coef) in m.varmap[vr.value]
        x += X[blk][i, j] * coef
    end
    x
end
function MOI.getattribute(m::SOItoMOIBridge, vp::MOI.VariablePrimal, vr::Vector{MOI.VariableReference})
    MOI.getattribute.(m, vp, vr)
end

function _getattribute(m::SOItoMOIBridge, cr::CR{<:ASF}, f)
    cs = m.constrmap[cr.value]
    @assert length(cs) == 1
    f(m, first(cs))
end
function _getattribute(m::SOItoMOIBridge, cr::CR{<:AVF}, f)
    f.(m, m.constrmap[cr.value])
end

function getslack(m::SOItoMOIBridge, c::Int)
    X = getX(m.sdsolver)
    blk, i, j, coef = m.slackmap[c]
    if iszero(blk)
        0.0
    else
        X[blk][i, j] * coef
    end
end

function MOI.getattribute(m::SOItoMOIBridge, a::MOI.ConstraintPrimal, cr::CR)
    _getattribute(m, cr, getslack)
end

function getvardual(m::SOItoMOIBridge, vi::UInt64)
    Z = getZ(m.sdsolver)
    z = 0.
    for (blk, i, j, coef) in m.varmap[vi]
        z += Z[blk][i, j] / coef
    end
    z
end
getvardual(m::SOItoMOIBridge, f::SVF) = getvardual(m, f.variable.value)
getvardual(m::SOItoMOIBridge, f::VVF) = map(vr -> getvardual(m, vr.value), f.variables)
function MOI.getattribute(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR{<:VF})
    getvardual(m, MOI.getattribute(m, MOI.ConstraintFunction(), cr))
end

function getdual(m::SOItoMOIBridge, c::Int)
    if c == 0
        0.
    else
        gety(m.sdsolver)[c]
    end
end
function MOI.getattribute(m::SOItoMOIBridge, ::MOI.ConstraintDual, cr::CR)
    _getattribute(m, cr, getdual)
end

end # module
