@testset "Linear" begin
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[1.0 0 0; 0 0 0; 0 0 2.0], [2.89241 0; 0 1.89241], [2.35008 0; 0 2.35008], [3.62629 0; 0 1.62629]],
        [0.0, -2.0, 0.0, 3.0, 1.0]))
    MOIT.lin1ftest(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 2.0]],
        [3.0, 1.0]))
    MOIT.lin1vtest(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[16.0], [3.0], [13.9279 0.0; 0.0 17.9279], [14.3814 0.0; 0.0 17.3814], [25.3959 0.0; 0.0 9.3959], [15.8211 0.0; 0.0 15.8211]],
        [-7.0, -2.0, 4.0, 0.0, 0.0, -7.0]))
    MOIT.lin2ftest(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[0.0], [16.0], [3.0], [75.9147 0.0; 0.0 79.9147]],
        [-7.0, -2.0, 4.0, 70.2306]))
    MOIT.lin2vtest(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        MOI.INFEASIBLE, MOI.INFEASIBLE_POINT, [-0.5, 0.5]))
    MOIT.lin3test(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        MOI.INFEASIBLE, MOI.INFEASIBLE_POINT, [-1.0]))
    MOIT.lin4test(bridged, config)
end
@testset "SOC" begin
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[1.0 0.7071 0.7071; 0.7071 1.0 -0.0; 0.7071 -0.0 1.0], [2.2375 0.0; 0.0 1.2375], [2.0662 0.0; 0.0 1.3591], [2.0662 0.0; 0.0 1.3591]],
        [1.4142, -0.7071, 1.0, -0.3536, 1.0, -0.7071, -0.3536]))
    MOIT.soc1ftest(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[1.0 0.7071 0.7071; 0.7071 1.0 -0.0; 0.7071 -0.0 1.0], [2.2375 0.0; 0.0 1.2375], [2.0662 0.0; 0.0 1.3591], [2.0662 0.0; 0.0 1.3591]],
        [1.4142, -0.7071, 1.0, -0.3536, 1.0, -0.7071, -0.3536]))
    MOIT.soc1vtest(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[0.0], [1.0 -0.7071 0.7071; -0.7071 1.0 0.0; 0.7071 0.0 1.0], [0.4046 0.0; 0.0 1.1117], [1.0895 0.0; 0.0 0.3824], [1.3079 0.0; 0.0 0.3079]],
        [-1, -√2, -√2/2, -1, -√2/4, 1, √2/2, -√2/4]))
    MOIT.soc2ntest(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[0.0], [1.0 -0.7071 0.7071; -0.7071 1.0 0.0; 0.7071 0.0 1.0], [0.4046 0.0; 0.0 1.1117], [1.0895 0.0; 0.0 0.3824], [1.3079 0.0; 0.0 0.3079]],
        [1, -√2, -√2/2, -1, -√2/4, 1, √2/2, -√2/4]))
    MOIT.soc2ptest(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                MOI.INFEASIBLE,
                                                                MOI.INFEASIBLE_POINT,
                                                                [-1.0, 1.0, -0.5, 1.0, -0.5]))
    MOIT.soc3test(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.0 0.8944 0.4472; 0.8944 1.0 0.0; 0.4472 0.0 1.0], [1.9332 0.0; 0.0 0.9332], [1.868 0.0; 0.0 0.9736], [1.5795 0.0; 0.0 1.1323], [1.868 0.0; 0.0 0.9736], [1.5795 0.0; 0.0 1.1323]],
                                                                [2.2361, 2.0, 1.0, -1.118, 2.0, -0.8944, 1.0, -0.8944, -0.2236]))
    MOIT.soc4test(bridged, config)
end
@testset "RSOC" begin
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    [[0.5 0.7071 0.7071; 0.7071 2.0 0.0; 0.7071 0.0 2.0], [0.8667 0.0; 0.0 0.1596], [0.8667 0.0; 0.0 0.1596]],
                                                                    [-1.4142, 1.0, -0.1768, 1.0, -0.3536, -0.1768]))
        MOIT.rotatedsoc1ftest(bridged, config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                            [[0.5 1/√2 1/√2; 1/√2 2.0 0.0; 1/√2 0.0 2.0], [0.0], [0.0], [3.65568 0; 0 2.94857], [3.65568 0; 0 2.94857]],
                                                            [-√2, 1.0, -0.176777, 1.0, -0.353555, -0.176777, 2189.72, 2189.01]))
        MOIT.rotatedsoc1vtest(bridged, config)
        MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                    MOI.INFEASIBLE,
                                                                    tuple(),
                                                                    [-47.8864, 47.5533, -46.2201, 141.088]))
        MOIT.rotatedsoc2test(bridged, config)
        # FIXME u >= 0.0 followed by u <= ub. We need to drop support for
        #       SingleVariable-in-LessThan but we need variable bridge otherwise,
        #       it creates a slack
#            MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
#                                                                        [[0.0], [0.0], [3.0], [0.0], [1.0 0.0; 0.0 0.0], [0.7071 1.0 0.0001; 1.0 1.4142 -0.0; 0.0001 -0.0 1.4142], [0.7071 1.7321; 1.7321 4.2426], [4.5908 0.0; 0.0 2.8588]],
#                                                                        [0.0, 0.0, 0.2887, -0.6124, 0.866, -0.3062, 0.0, -0.0001, -0.0, -1.2247, 1.0, -0.2041]))
#            MOIT.rotatedsoc3test(bridged, config)
end
@testset "GeoMean" begin
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[√2 2.0; 2.0 2*√2], [1.0 √2; √2 2.0], [1.0 √2; √2 2.0], [0.0], [0.0],
         [1.83238 0; 0 0.832376], [1.77579 0; 0 0.775793],
         [1.77719 0; 0 0.777194], [1.77283 0; 0 0.772829],
         [2.50173 0; 0 0.501735], [2.07437 0; 0 0.660159],
         [2.06708 0; 0 0.652869]],
        [1.0, -0.471405, 2/3, -0.235702, -1/3, 0.471405, -1/6, -1/3, 0.471405,
         -1/6, 1/3]))
    MOIT.geomean1ftest(bridged, config)
    MOIT.geomean1vtest(bridged, config)
end
@testset "PSD" begin
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[1.0 1.0; 1.0 1.0], [2.2633 0.0; 0.0 1.2633], [2.2669 0.0; 0.0 1.2669], [2.2633 0.0; 0.0 1.2633]],
        [-1.0, 2.0, -1.0, -2.0]))
    MOIT.psdt0ftest(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[1.0 1.0; 1.0 1.0]],
        [-2.0]))
    MOIT.psdt0vtest(bridged, config)
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
        [-cX0, -1.35619, -cX0, 2y2, -1.35619, -cX0, -y2/√2, y2, -y2/√2/2, y2, -y2/√2, -y2/√2/2, y2*√2-1, -y2]))
    MOIT.psdt1ftest(bridged, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[Xv[1] Xv[2] Xv[1]; Xv[2] Xv[3] Xv[2]; Xv[1] Xv[2] Xv[1]],
         [xv[1] xv[2] xv[2]; xv[2] xv[1] 0.0; xv[2] 0.0 xv[1]],
         [0.486623 0; 0 0.486623-xv[1]],
         [0.4475 0; 0 0.4475-xv[2]],
         [0.4475 0; 0 0.4475-xv[2]]],
        [-y2/√2, y2, -y2/√2/2, y2, -y2/√2, -y2/√2/2, y2*√2-1, -y2]))
    MOIT.psdt1vtest(bridged, config)
	MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
        [[20/3 0 0 0 0 0; 0 0 0 0 0 0; 0 0 10/3 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0],
         [4.09808 -2.12132; -2.12132 1.09808], [0.0], [7.32746 0; 0 0.660794],
         [1.79451 0; 0 1.79451], [4.18806 0; 0 0.854733], [1.79451 0; 0 1.79451],
         [1.79451 0; 0 1.79451], [1.79451 0; 0 1.79451], [2.85043 0; 0 0.94851]],
         [-0.190192, 0.0, 0.125977, 0.0, 0.142644, 0.142644, 0.0127405, -0.211325, -0.816497, -0.788675, 0.0]))
    MOIT.psdt2test(bridged, config)
end
