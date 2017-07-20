using SemidefiniteOptInterface
using Base.Test

using CSDP
solvers = [CSDP.CSDPSolver()]

@testset "Linear tests with $solver" for solver in solvers
    include(joinpath(Pkg.dir("MathOptInterface"), "test", "contlinear.jl"))
    contlineartest(solver, 1e-7)
end
#@testset "Conic tests with $solver" for solver in solvers
#    include(joinpath(Pkg.dir("MathOptInterface"), "test", "contconic.jl"))
#    contconictest(solver, 1e-7)
#end
