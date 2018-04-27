using SemidefiniteOptInterface
using Base.Test

import CSDP
optimizers = [CSDP.CSDPOptimizer(printlevel=0)]

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIB = MOI.Bridges

const MOIU = MOI.Utilities
MOIU.@model SDModelData () (EqualTo, GreaterThan, LessThan) (Zeros, Nonnegatives, Nonpositives, PositiveSemidefiniteConeTriangle) () (SingleVariable,) (ScalarAffineFunction,) (VectorOfVariables,) (VectorAffineFunction,)

const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

# name too long with $optimizer, waiting for https://github.com/JuliaOpt/MathOptInterface.jl/pull/201
@testset "Linear tests with optimizer" for optimizer in optimizers
    @test MOI.isempty(optimizer)
    MOI.empty!(optimizer)
    @test MOI.isempty(optimizer)
    MOIT.contlineartest(MOIB.SplitInterval{Float64}(MOIU.CachingOptimizer(SDModelData{Float64}(), optimizer)), config)
end
@testset "Conic tests with optimizer" for optimizer in optimizers
    MOI.empty!(optimizer)
    @test MOI.isempty(optimizer)
    MOIT.contconictest(MOIB.RootDet{Float64}(MOIB.GeoMean{Float64}(MOIB.RSOCtoPSD{Float64}(MOIB.SOCtoPSD{Float64}(MOIU.CachingOptimizer(SDModelData{Float64}(), optimizer))))), config, ["psds", "rootdets", "logdet", "exp"])
end
