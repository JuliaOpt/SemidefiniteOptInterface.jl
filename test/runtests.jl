using SemidefiniteOptInterface
const SDOI = SemidefiniteOptInterface

using Compat
using Compat.Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIB = MOI.Bridges

include("sdpa.jl")

const MOIU = MOI.Utilities
MOIU.@model SDModelData () (EqualTo, GreaterThan, LessThan) (Zeros, Nonnegatives, Nonpositives, PositiveSemidefiniteConeTriangle) () (SingleVariable,) (ScalarAffineFunction,) (VectorOfVariables,) (VectorAffineFunction,)

mock = SDOI.MockSDOptimizer{Float64}()
mock_optimizer = SDOI.SDOIOptimizer(mock, Float64)
cached_mock_optimizer = MOIU.CachingOptimizer(SDModelData{Float64}(), mock_optimizer)
config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

include("unit.jl")
include("contlinear.jl")
include("contconic.jl")
