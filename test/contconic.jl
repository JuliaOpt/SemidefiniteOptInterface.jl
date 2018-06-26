@testset "MOI Continuous Conic" begin
    @testset "Linear" begin
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 2.0], [2.8924 0.0; 0.0 1.8924], [2.3501 0.0; 0.0 2.3501], [3.6263 0.0; 0.0 1.6263]],
                                                                    [3.0, 1.0, -0.0, -2.0, -0.0]))
        MOIT.lin1ftest(cached_mock_optimizer, config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 2.0]],
                                                                    [3.0, 1.0]))
        MOIT.lin1vtest(cached_mock_optimizer, config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[16.0], [3.0], [13.9279 0.0; 0.0 17.9279], [14.3814 0.0; 0.0 17.3814], [25.3959 0.0; 0.0 9.3959], [15.8211 0.0; 0.0 15.8211]],
                                                                    [-7.0, -2.0, 4.0, -7.0, -0.0, 0.0]))
        MOIT.lin2ftest(cached_mock_optimizer, config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[0.0], [16.0], [3.0], [75.9147 0.0; 0.0 79.9147]],
                                                                    [70.2006, -7.0, -2.0, 4.0]))
        MOIT.lin2vtest(cached_mock_optimizer, config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    tuple(),
                                                                    [-0.5, 0.5]))
        MOIT.lin3test(cached_mock_optimizer, config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    tuple(),
                                                                    [-1.0]))
        MOIT.lin4test(cached_mock_optimizer, config)
    end
    @testset "SOC" begin
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[1.0 0.7071 0.7071; 0.7071 1.0 -0.0; 0.7071 -0.0 1.0], [2.2375 0.0; 0.0 1.2375], [2.0662 0.0; 0.0 1.3591], [2.0662 0.0; 0.0 1.3591]],
                                                                    [1.4142, -0.7071, 1.0, -0.3536, 1.0, -0.7071, -0.3536]))
        MOIT.soc1ftest(MOIB.SOCtoPSD{Float64}(cached_mock_optimizer), config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[1.0 0.7071 0.7071; 0.7071 1.0 -0.0; 0.7071 -0.0 1.0], [2.2375 0.0; 0.0 1.2375], [2.0662 0.0; 0.0 1.3591], [2.0662 0.0; 0.0 1.3591]],
                                                                    [1.4142, -0.7071, 1.0, -0.3536, 1.0, -0.7071, -0.3536]))
        MOIT.soc1vtest(MOIB.SOCtoPSD{Float64}(cached_mock_optimizer), config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[0.0], [1.0 -0.7071 0.7071; -0.7071 1.0 0.0; 0.7071 0.0 1.0], [0.4046 0.0; 0.0 1.1117], [1.0895 0.0; 0.0 0.3824], [1.3079 0.0; 0.0 0.3079]],
                                                                    [-1.4142, -1.0, -0.7071, -1.0, -0.3536, 1.0, 0.7071, -0.3536]))
        MOIT.soc2ntest(MOIB.SOCtoPSD{Float64}(cached_mock_optimizer), config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[0.0], [1.0 -0.7071 0.7071; -0.7071 1.0 0.0; 0.7071 0.0 1.0], [0.4046 0.0; 0.0 1.1117], [1.0895 0.0; 0.0 0.3824], [1.3079 0.0; 0.0 0.3079]],
                                                                    [-1.4142, 1.0, -0.7071, -1.0, -0.3536, 1.0, 0.7071, -0.3536]))
        MOIT.soc2ptest(MOIB.SOCtoPSD{Float64}(cached_mock_optimizer), config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    tuple(),
                                                                    [-1.0, 1.0, -0.5, 1.0, -0.5]))
        MOIT.soc3test(MOIB.SOCtoPSD{Float64}(cached_mock_optimizer), config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[1.0 0.8944 0.4472; 0.8944 1.0 0.0; 0.4472 0.0 1.0], [1.9332 0.0; 0.0 0.9332], [1.868 0.0; 0.0 0.9736], [1.5795 0.0; 0.0 1.1323], [1.868 0.0; 0.0 0.9736], [1.5795 0.0; 0.0 1.1323]],
                                                                    [2.2361, 2.0, 1.0, -1.118, 2.0, -0.8944, 1.0, -0.8944, -0.2236]))
        MOIT.soc4test(MOIB.SOCtoPSD{Float64}(cached_mock_optimizer), config)
    end
    @testset "RSOC" begin
            MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                        [[0.5 0.7071 0.7071; 0.7071 2.0 0.0; 0.7071 0.0 2.0], [0.8667 0.0; 0.0 0.1596], [0.8667 0.0; 0.0 0.1596]],
                                                                        [-1.4142, 1.0, -0.1768, 1.0, -0.3536, -0.1768]))
            MOIT.rotatedsoc1ftest(MOIB.RSOCtoPSD{Float64}(cached_mock_optimizer), config)
            MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                        [[0.0], [0.0], [0.5 0.7071 0.7071; 0.7071 2.0 0.0; 0.7071 0.0 2.0], [2.6164 0.0; 0.0 1.9093], [2.6164 0.0; 0.0 1.9093]],
                                                                        [5183.15, 5182.44, -1.4142, 1.0, -0.1768, 1.0, -0.3536, -0.1768]))
            MOIT.rotatedsoc1vtest(MOIB.RSOCtoPSD{Float64}(cached_mock_optimizer), config)
            MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                        tuple(),
                                                                        [141.088, -47.8864, 47.5533, -46.2201]))
            MOIT.rotatedsoc2test(MOIB.RSOCtoPSD{Float64}(cached_mock_optimizer), config)
            MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                        [[0.0], [0.0], [3.0], [0.0], [1.0 0.0; 0.0 0.0], [0.7071 1.0 0.0001; 1.0 1.4142 -0.0; 0.0001 -0.0 1.4142], [0.7071 1.7321; 1.7321 4.2426], [4.5908 0.0; 0.0 2.8588]],
                                                                        [0.0, 0.0, 0.2887, -0.6124, 0.866, -0.3062, 0.0, -0.0001, -0.0, -1.2247, 1.0, -0.2041]))
            MOIT.rotatedsoc3test(MOIB.RSOCtoPSD{Float64}(cached_mock_optimizer), config)
    end
    @testset "GeoMean" begin
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[0.0], [0.0], [1.4142 2.0; 2.0 2.8284], [1.0 1.4142; 1.4142 2.0], [1.0 1.4142; 1.4142 2.0], [1.8324 0.0; 0.0 0.8324], [1.7758 0.0; 0.0 0.7758], [1.7772 0.0; 0.0 0.7772], [1.7728 0.0; 0.0 0.7728], [2.5017 0.0; 0.0 0.5017], [2.0744 0.0; 0.0 0.6602], [2.0671 0.0; 0.0 0.6529]],
                                                                    [1.0, 0.3333, -0.4714, 0.6667, -0.2357, -0.3333, 0.4714, -0.1667, -0.3333, 0.4714, -0.1667]))
        MOIT.geomean1ftest(MOIB.GeoMean{Float64}(MOIB.RSOCtoPSD{Float64}(cached_mock_optimizer)), config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[0.0], [0.0], [1.4142 2.0; 2.0 2.8284], [1.0 1.4142; 1.4142 2.0], [1.0 1.4142; 1.4142 2.0], [1.8324 0.0; 0.0 0.8324], [1.7758 0.0; 0.0 0.7758], [1.7772 0.0; 0.0 0.7772], [1.7728 0.0; 0.0 0.7728], [2.5017 0.0; 0.0 0.5017], [2.0744 0.0; 0.0 0.6602], [2.0671 0.0; 0.0 0.6529]],
                                                                    [1.0, 0.3333, -0.4714, 0.6667, -0.2357, -0.3333, 0.4714, -0.1667, -0.3333, 0.4714, -0.1667]))
        MOIT.geomean1vtest(MOIB.GeoMean{Float64}(MOIB.RSOCtoPSD{Float64}(cached_mock_optimizer)), config)
    end
    @testset "PSD" begin
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[1.0 1.0; 1.0 1.0], [2.2633 0.0; 0.0 1.2633], [2.2669 0.0; 0.0 1.2669], [2.2633 0.0; 0.0 1.2633]],
                                                                    [-2.0, -1.0, 2.0, -1.0]))
        MOIT.psdt0ftest(cached_mock_optimizer, config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[1.0 1.0; 1.0 1.0]],
                                                                    [-2.0]))
        MOIT.psdt0vtest(cached_mock_optimizer, config)
        # PSD1: see comments in MOI/src/Test/contconic.jl to see where these constants come from
        δ = √(1 + (3*√2+2)*√(-116*√2+166) / 14) / 2
        ε = √((1 - 2*(√2-1)*δ^2) / (2-√2))
        y2 = 1 - ε*δ
        y1 = 1 - √2*y2
        obj = y1 + y2/2
        k = -2*δ/ε
        x2 = ((3-2obj)*(2+k^2)-4) / (4*(2+k^2)-4*√2)
        α = √(3-2obj-4x2)/2
        β = k*α
        Xv = [α^2, α*β, β^2, α^2, α*β, α^2]
        xv = [√2*x2, x2, x2]
        cX0 = 1+(√2-1)*y2
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[Xv[1] Xv[2] Xv[1]; Xv[2] Xv[3] Xv[2]; Xv[1] Xv[2] Xv[1]],
                                                                    [0.254409 0.179895 0.179895; 0.179895 0.254409 0.0; 0.179895 0.0 0.254409],
                                                                    [0.349485 0; 0 0.132235],
                                                                    [0.114062 0; 0 0.374032],
                                                                    [0.419535 0; 0 0.108444],
                                                                    [0.313002 0; 0 0.0957514],
                                                                    [0.114062 0; 0 0.374032],
                                                                    [0.349485 0; 0 0.132235],
                                                                    [0.354357 0; 0 0.0999473],
                                                                    [0.332125 0; 0 0.152231],
                                                                    [0.332125 0; 0 0.152231]],
                                                                    [y2*√2-1, -y2, -cX0, -1.35619, -cX0, 2y2, -1.35619, -cX0, -y2/√2, y2, -y2/√2/2, y2, -y2/√2, -y2/√2/2]))
        MOIT.psdt1ftest(MOIB.SOCtoPSD{Float64}(cached_mock_optimizer), config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[Xv[1] Xv[2] Xv[1]; Xv[2] Xv[3] Xv[2]; Xv[1] Xv[2] Xv[1]],
                                                                     [xv[1] xv[2] xv[2]; xv[2] xv[1] 0.0; xv[2] 0.0 xv[1]],
                                                                     [0.486623 0; 0 0.486623-xv[1]],
                                                                     [0.4475 0; 0 0.4475-xv[2]],
                                                                     [0.4475 0; 0 0.4475-xv[2]]],
                                                                    [y2*√2-1, -y2, -y2/√2, y2, -y2/√2/2, y2, -y2/√2, -y2/√2/2]))
        MOIT.psdt1vtest(MOIB.SOCtoPSD{Float64}(cached_mock_optimizer), config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[0.0], [20/3 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 10/3 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0], [4.0981 -2.1213; -2.1213 1.0981], [7.3275 0.0; 0.0 0.6608], [1.7945 0.0; 0.0 1.7945], [4.1881 0.0; 0.0 0.8547], [1.7945 0.0; 0.0 1.7945], [1.7945 0.0; 0.0 1.7945], [1.7945 0.0; 0.0 1.7945], [2.8504 0.0; 0.0 0.9485]],
                                                                    [0.0, -0.1902, -0.0, 0.126, -0.0, 0.1426, 0.1426, 0.0127, -0.2113, -0.8165, -0.7887]))
        MOIT.psdt2test(cached_mock_optimizer, config)
    end
end
