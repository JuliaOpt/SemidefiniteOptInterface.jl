using SemidefiniteOptInterface
using Base.Test

import CSDP
optimizers = [CSDP.CSDPOptimizer(printlevel=0)]

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

MOIU.@model SDModelData () (EqualTo, GreaterThan, LessThan) (Zeros, Nonnegatives, Nonpositives, PositiveSemidefiniteConeTriangle) () (SingleVariable,) (ScalarAffineFunction,) (VectorOfVariables,) (VectorAffineFunction,)

using MathOptInterfaceBridges
const MOIB = MathOptInterfaceBridges

MOIB.@bridge SplitInterval MOIB.SplitIntervalBridge () (Interval,) () () () (ScalarAffineFunction,) () ()
MOIB.@bridge SOCtoPSDC MOIB.SOCtoPSDCBridge () () (SecondOrderCone,) () () () (VectorOfVariables,) (VectorAffineFunction,)
MOIB.@bridge RSOCtoPSDC MOIB.RSOCtoPSDCBridge () () (RotatedSecondOrderCone,) () () () (VectorOfVariables,) (VectorAffineFunction,)
MOIB.@bridge GeoMean MOIB.GeoMeanBridge () () (GeometricMeanCone,) () () () (VectorOfVariables,) (VectorAffineFunction,)
MOIB.@bridge RootDet MOIB.RootDetBridge () () (RootDetConeTriangle,) () () () (VectorOfVariables,) (VectorAffineFunction,)

const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

@testset "Linear tests with $optimizer" for optimizer in optimizers
    MOIT.contlineartest(SplitInterval{Float64}(MOIU.CachingOptimizer(SDModelData{Float64}(), optimizer)), config)
end
@testset "Conic tests with $optimizer" for optimizer in optimizers
    MOIT.contconictest(RootDet{Float64}(GeoMean{Float64}(RSOCtoPSDC{Float64}(SOCtoPSDC{Float64}(MOIU.CachingOptimizer(SDModelData{Float64}(), optimizer))))), config, ["logdet", "exp"])
end
