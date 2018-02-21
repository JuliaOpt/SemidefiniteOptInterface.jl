"""
    init!(optimizer::AbstractSDOptimizer, blkdims::Vector{Int}, nconstrs::Integer)

Initialize the optimizer with nconstrs constraints and blkdims blocks.
"""
function init! end

"""
    setconstraintconstant!(optimizer::AbstractSDOptimizer, val, constr::Integer)

Sets the entry `constr` of `b` to `val`.
"""
function setconstraintconstant! end

"""
    setconstraintcoefficient!(optimizer::AbstractSDOptimizer, val, constr::Integer, blk::Integer, i::Integer, j::Integer)

Sets the entry `i`, `j` of the block `blk` of the matrix of the constraint `constr` to `val`.
"""
function setconstraintcoefficient! end

"""
    setobjectivecoefficient!(optimizer::AbstractSDOptimizer, val, blk::Integer, i::Integer, j::Integer)

Sets the entry `i`, `j` of the block `blk` of the objective matrix to `val`.
"""
function setobjectivecoefficient! end

"""
    getX(optimizer::AbstractSDOptimizer)

Returns the solution X as a block matrix.
"""
function getX end

"""
    gety(optimizer::AbstractSDOptimizer)

Returns the solution y.
"""
function gety end

"""
    getZ(optimizer::AbstractSDOptimizer)

Returns the solution Z.
"""
function getZ end

"""
    getprimalobjectivevalue(optimizer::AbstractSDOptimizer)

Returns the primal objective value.
"""
function getprimalobjectivevalue end

"""
    getdualobjectivevalue(optimizer::AbstractSDOptimizer)

Returns the dual objective value.
"""
function getdualobjectivevalue end
