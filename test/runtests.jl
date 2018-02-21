using SemidefiniteOptInterface
using Base.Test

import CSDP
optimizers = [CSDP.CSDPOptimizer(printlevel=0)]

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test

const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

@testset "Linear tests with $optimizer" for optimizer in optimizers
    MOIT.contlineartest(optimizer, config)
end
@testset "Conic tests with $optimizer" for optimizer in optimizers
    MOIT.contconictest(optimizer, config, ["logdet", "exp"])
end
