using SemidefiniteOptInterface
using Base.Test

using CSDP
solvers = [CSDP.CSDPSolver(printlevel=0)]

using MathOptInterfaceTests
const MOIT = MathOptInterfaceTests

const config = MOIT.TestConfig(1e-7, 1e-7, true, true, true)

@testset "Linear tests with $solver" for solver in solvers
    MOIT.contlineartest(solver, config)
end
@testset "Conic tests with $solver" for solver in solvers
    MOIT.contconictest(solver, config)
end
