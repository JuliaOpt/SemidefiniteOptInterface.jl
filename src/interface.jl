"""
    SDSolverInstance(s::AbstractSDSolver)

Creates an solver instance.
"""
function SDSolverInstance end

"""
    initinstance!(s::AbstractSDSolverInstance, blkdims::Vector{Int}, nconstrs::Integer)

Initialize the instance with nconstrs constraints and blkdims blocks.
"""
function initinstance! end

"""
    setconstraintconstant!(s::AbstractSDSolverInstance, val, constr::Integer)

Sets the entry `constr` of `b` to `val`.
"""
function setconstraintconstant! end

"""
    setconstraintcoefficient!(m::AbstractSDSolverInstance, val, constr::Integer, blk::Integer, i::Integer, j::Integer)

Sets the entry `i`, `j` of the block `blk` of the matrix of the constraint `constr` to `val`.
"""
function setconstraintcoefficient! end

"""
    setobjectivecoefficient!(m::AbstractSDSolverInstance, val, blk::Integer, i::Integer, j::Integer)

Sets the entry `i`, `j` of the block `blk` of the objective matrix to `val`.
"""
function setobjectivecoefficient! end

"""
    getX(m::AbstractSDSolverInstance)

Returns the solution X as a block matrix.
"""
function getX end

"""
    gety(m::AbstractSDSolverInstance)

Returns the solution y.
"""
function gety end

"""
    getZ(m::AbstractSDSolverInstance)

Returns the solution Z.
"""
function getZ end

"""
    getprimalobjectivevalue(m::AbstractSDSolverInstance)

Returns the primal objective value.
"""
function getprimalobjectivevalue end

"""
    getdualobjectivevalue(m::AbstractSDSolverInstance)

Returns the dual objective value.
"""
function getdualobjectivevalue end
