struct SDConstraints{FT}
    zeros::Vector{Tuple{UInt64, FT}}
    nnegs::Vector{Tuple{UInt64, FT}}
    nposs::Vector{Tuple{UInt64, FT}}
    psdcs::Vector{Tuple{UInt64, FT}}
    function SDConstraints{FT}() where FT
        new{FT}(Tuple{UInt64, FT}[], Tuple{UInt64, FT}[], Tuple{UInt64, FT}[], Tuple{UInt64, FT}[])
    end
end

const ZS = Union{MOI.EqualTo, MOI.Zeros}
const NS = Union{MOI.GreaterThan, MOI.Nonnegatives}
const PS = Union{MOI.LessThan, MOI.Nonpositives}
const DS = MOI.PositiveSemidefiniteConeTriangle

const SS = Union{MOI.EqualTo, MOI.LessThan, MOI.GreaterThan}

_getconstant(s::MOI.EqualTo) = s.value
_getconstant(s::MOI.LessThan) = s.upper
_getconstant(s::MOI.GreaterThan) = s.lower

_cst(f, s) = f
function _cst{T}(f::SAF{T}, s::SS)
    SAF{T}(f.variables, f.coefficients, f.constant - _cst(s))
end

function _addconstraint!{FT}(constrs::Vector{Tuple{UInt64, FT}}, ci::UInt64, f::FT, s)
    push!(constrs, (ci, _cst(f, s)))
    length(constrs)
end
_addconstraint!(m::SDConstraints, ci::UInt64, f, s::ZS) = _addconstraint!(m.zeros, ci, f, s)
_addconstraint!(m::SDConstraints, ci::UInt64, f, s::NS) = _addconstraint!(m.nnegs, ci, f, s)
_addconstraint!(m::SDConstraints, ci::UInt64, f, s::PS) = _addconstraint!(m.nposs, ci, f, s)
_addconstraint!(m::SDConstraints, ci::UInt64, f, s::DS) = _addconstraint!(m.psdcs, ci, f, s)

function _modifyconstraint!{FT}(constrs::Vector{Tuple{UInt64, FT}}, i::Int, cr::CR{FT}, change::FT)
    constrs[i] = change
end
function _modifyconstraint!{FT}(constrs::Vector{Tuple{UInt64, FT}}, i::Int, cr::CR{FT}, change::MOI.AbstractFunctionModification)
    constrs[i] = MOI.modifyfunction(constrs[i], change)
end
function _modifyconstraint!{FT}(constrs::Vector{Tuple{UInt64, FT}}, i::Int, cr::CR{FT}, change::Function)
    constrs[i] = change(constrs[i])
end
_modifyconstraint!{FT}(m::SDConstraints, i::Int, cr::CR{FT, <:ZS}, change) = _modifyconstraint!(m.zeros, i, cr, change)
_modifyconstraint!{FT}(m::SDConstraints, i::Int, cr::CR{FT, <:NS}, change) = _modifyconstraint!(m.nnegs, i, cr, change)
_modifyconstraint!{FT}(m::SDConstraints, i::Int, cr::CR{FT, <:PS}, change) = _modifyconstraint!(m.nposs, i, cr, change)
_modifyconstraint!{FT}(m::SDConstraints, i::Int, cr::CR{FT, <:DS}, change) = _modifyconstraint!(m.psdcs, i, cr, change)

_zstype{T}(::Type{<:ASF{T}}) = MOI.EqualTo{T}
_zstype(::Type{<:AVF}) = MOI.Zeros
_nstype{T}(::Type{<:ASF{T}}) = MOI.GreaterThan{T}
_nstype(::Type{<:AVF}) = MOI.Nonnegatives
_pstype{T}(::Type{<:ASF{T}}) = MOI.LessThan{T}
_pstype(::Type{<:AVF}) = MOI.Nonpositives
_dstype(::Type{<:MOI.AbstractFunction}) = MOI.PositiveSemidefiniteConeTriangle
_zlist{FT}(m::SDConstraints{FT}) = isempty(m.zeros) ? [] : [(FT, _zstype(FT))]
_nlist{FT}(m::SDConstraints{FT}) = isempty(m.nnegs) ? [] : [(FT, _nstype(FT))]
_plist{FT}(m::SDConstraints{FT}) = isempty(m.nposs) ? [] : [(FT, _pstype(FT))]
_dlist{FT}(m::SDConstraints{FT}) = isempty(m.psdcs) ? [] : [(FT, _dstype(FT))]
function MOI.getattribute{FT}(m::SDConstraints{FT}, loc::MOI.ListOfConstraints)
    [_zlist(m); _nlist(m); _plist(m); _dlist(m)]
end


MOI.getattribute(m::SDConstraints, noc::MOI.NumberOfConstraints{<:Any, <:ZS}) = length(m.zeros)
MOI.getattribute(m::SDConstraints, noc::MOI.NumberOfConstraints{<:Any, <:NS}) = length(m.nnegs)
MOI.getattribute(m::SDConstraints, noc::MOI.NumberOfConstraints{<:Any, <:PS}) = length(m.nposs)
MOI.getattribute(m::SDConstraints, noc::MOI.NumberOfConstraints{<:Any, <:DS}) = length(m.psdcs)

mutable struct SDInstance{T}
    sense::MOI.OptimizationSense
    objective::SAF{T}
    sv::SDConstraints{SVF}
    sa::SDConstraints{SAF{T}}
    vv::SDConstraints{VVF}
    va::SDConstraints{VAF{T}}
    rhs::Vector{T}
    nvars::UInt64
    nconstrs::UInt64
    constrmap::Vector{Int} # Constraint Reference value ci -> index in array in sdinstance
    function SDInstance{T}() where T
        new{T}(true, SAF{T}(MOI.VariableReference[], T[], zero(T)),
               SDConstraints{SVF}(), SDConstraints{SAF{T}}(), SDConstraints{VVF}(), SDConstraints{VAF{T}}(),
               T[], 0, 0, Int[])
    end
end

# Variables
MOI.getattribute(m::SDInstance, ::MOI.NumberOfVariables) = m.nvars
MOI.addvariable!(m::SDInstance) = MOI.VariableReference(m.nvars += 1)
function MOI.addvariables!(m::SDInstance, n::Integer)
    [MOI.addvariable!(m) for i in 1:n]
end

# Constraints
_addconstraint!(m::SDInstance, ci::UInt64, f::SVF, s) = _addconstraint!(m.sv, ci, f, s)
_addconstraint!(m::SDInstance, ci::UInt64, f::SAF, s) = _addconstraint!(m.sa, ci, f, s)
_addconstraint!(m::SDInstance, ci::UInt64, f::VVF, s) = _addconstraint!(m.vv, ci, f, s)
_addconstraint!(m::SDInstance, ci::UInt64, f::VAF, s) = _addconstraint!(m.va, ci, f, s)

function _setrhs(m::SDInstance, ci::UInt64, s) end
function _setrhs(m::SDInstance, ci::UInt64, s::SS)
    m.rhs[ci] = _getconstant(s)
end

function MOI.addconstraint!{FT, ST}(m::SDInstance, f::FT, s::ST)
    push!(m.constrmap, _addconstraint!(m, m.nconstrs += 1, f, s))
    _setrhs(m, m.nconstrs, s)
    CR{FT, ST}(m.nconstrs)
end

_modifyconstraint!(m::SDInstance, i::Int, cr::CR{<:SVF}, change) = _modifyconstraint!(m.sv, i, cr, change)
_modifyconstraint!(m::SDInstance, i::Int, cr::CR{<:SAF}, change) = _modifyconstraint!(m.sa, i, cr, change)
_modifyconstraint!(m::SDInstance, i::Int, cr::CR{<:VVF}, change) = _modifyconstraint!(m.vv, i, cr, change)
_modifyconstraint!(m::SDInstance, i::Int, cr::CR{<:VAF}, change) = _modifyconstraint!(m.va, i, cr, change)

function _modifyconstraint!{FT, ST}(m::SDInstance, i::Int, cr::CR{FT, ST}, change::ST) end
function _modifyconstraint!(m::SDInstance, i::Int, cr::CR{<:SVF}, change::SS)
    m.rhs[cr.value] = _getconstant(change)
end
function _modifyconstraint!(m::SDInstance, i::Int, cr::CR{<:SAF}, change::SS)
    x = _getconstant(change)
    Δ = x - m.rhs[cr.value]
    m.rhs[cr.value] = x
    _modifyconstant(m, i, cr, f -> MOI.modifyfunction(f, MOI.ScalarConstantChange(f.constant + Δ)))
end

function MOI.modifyconstraint!(m::SDInstance, cr::CR, change)
    _modifyconstraint!(m, m.constrmap[cr.value], cr, change)
end

function MOI.getattribute(m::SDInstance, loc::MOI.ListOfConstraints)
    [MOI.getattribute(m.sv, loc); MOI.getattribute(m.sa, loc); MOI.getattribute(m.vv, loc); MOI.getattribute(m.va, loc)]
end

MOI.getattribute(m::SDInstance, noc::MOI.NumberOfConstraints{<:SVF}) = MOI.getattribute(m.sv, noc)
MOI.getattribute(m::SDInstance, noc::MOI.NumberOfConstraints{<:SAF}) = MOI.getattribute(m.sa, noc)
MOI.getattribute(m::SDInstance, noc::MOI.NumberOfConstraints{<:VVF}) = MOI.getattribute(m.vv, noc)
MOI.getattribute(m::SDInstance, noc::MOI.NumberOfConstraints{<:VAF}) = MOI.getattribute(m.va, noc)
