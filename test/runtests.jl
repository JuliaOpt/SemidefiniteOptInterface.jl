using SemidefiniteOptInterface
using Base.Test

using CSDP
@testset "Conic tests" begin
    include(joinpath(Pkg.dir("MathOptInterface"), "test", "contconic.jl"))
    contconictest(CSDP.CSDPSolver(), 1e-7)
end
