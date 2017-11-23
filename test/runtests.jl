using SemidefiniteOptInterface
using Base.Test

import CSDP
solvers = [() -> CSDP.CSDPInstance(printlevel=0)]

using MathOptInterfaceTests
const MOIT = MathOptInterfaceTests

# Need 1e-5 for geomean test
const config = MOIT.TestConfig(1e-5, 1e-5, true, true, true)

@testset "Linear tests with $solver" for solver in solvers
    MOIT.contlineartest(solver, config)
end
@testset "Conic tests with $solver" for solver in solvers
    MOIT.contconictest(solver, config, ["logdet", "rootdet1qv", "rootdet1qf"])
end
