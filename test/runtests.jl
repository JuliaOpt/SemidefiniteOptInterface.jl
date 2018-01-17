using SemidefiniteOptInterface
using Base.Test

import CSDP
solvers = [() -> CSDP.CSDPInstance(printlevel=0)]

using MathOptInterfaceTests
const MOIT = MathOptInterfaceTests

const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

@testset "Linear tests with $solver" for solver in solvers
    MOIT.contlineartest(solver, config)
end
@testset "Conic tests with $solver" for solver in solvers
    MOIT.contconictest(solver, config, ["logdet", "exp"])
end
