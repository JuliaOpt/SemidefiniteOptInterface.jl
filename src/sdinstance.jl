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

_cst(s::MOI.EqualTo) = s.value
_cst(s::MOI.LessThan) = s.upper
_cst(s::MOI.GreaterThan) = s.lower
_cst(f, s) = f
function _cst{T}(f::SAF{T}, s)
    SAF{T}(f.variables, f.coefficients, f.constant - _cst(s))
end

function _addconstraint!{FT, ST}(constrs::Vector{Tuple{UInt64, FT}}, ci::UInt64, f::FT, s::ST)
    push!(constrs, (ci, _cst(f, s)))
    MOI.ConstraintReference{FT, ST}(ci)
end
_addconstraint!(m::SDConstraints, ci::UInt64, f, s::ZS) = _addconstraint!(m.zeros, ci, f, s)
_addconstraint!(m::SDConstraints, ci::UInt64, f, s::NS) = _addconstraint!(m.nnegs, ci, f, s)
_addconstraint!(m::SDConstraints, ci::UInt64, f, s::PS) = _addconstraint!(m.nposs, ci, f, s)
_addconstraint!(m::SDConstraints, ci::UInt64, f, s::DS) = _addconstraint!(m.psdcs, ci, f, s)

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
    function SDInstance{T}() where T
        new{T}(true, SAF{T}(MOI.VariableReference[], T[], zero(T)), SDConstraints{SVF}(), SDConstraints{SAF{T}}(), SDConstraints{VVF}(), SDConstraints{VAF{T}}())
    end
end

_addconstraint!(m::SDInstance, ci::UInt64, f::SVF, s) = _addconstraint!(m.sv, ci, f, s)
_addconstraint!(m::SDInstance, ci::UInt64, f::SAF, s) = _addconstraint!(m.sa, ci, f, s)
_addconstraint!(m::SDInstance, ci::UInt64, f::VVF, s) = _addconstraint!(m.vv, ci, f, s)
_addconstraint!(m::SDInstance, ci::UInt64, f::VAF, s) = _addconstraint!(m.va, ci, f, s)

function MOI.getattribute(m::SDInstance, loc::MOI.ListOfConstraints)
    [MOI.getattribute(m.sv, loc); MOI.getattribute(m.sa, loc); MOI.getattribute(m.vv, loc); MOI.getattribute(m.va, loc)]
end

MOI.getattribute(m::SDInstance, noc::MOI.NumberOfConstraints{<:SVF}) = MOI.getattribute(m.sv, noc)
MOI.getattribute(m::SDInstance, noc::MOI.NumberOfConstraints{<:SAF}) = MOI.getattribute(m.sa, noc)
MOI.getattribute(m::SDInstance, noc::MOI.NumberOfConstraints{<:VVF}) = MOI.getattribute(m.vv, noc)
MOI.getattribute(m::SDInstance, noc::MOI.NumberOfConstraints{<:VAF}) = MOI.getattribute(m.va, noc)
