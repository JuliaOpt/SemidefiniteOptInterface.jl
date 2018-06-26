@testset "MOI Continuous Linear" begin
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.0], [0.0], [0.0]],
                                                                [1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.0], [0.0], [0.0]],
                                                                [1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0], [1.0], [0.0]],
                                                                [2.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0], [2.0], [0.0]],
                                                                [2.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [1.0], [0.0], [0.0]],
                                                                [54.4045, 1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [2.0], [0.0]],
                                                                [61.2652, 1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0], [2.0]],
                                                                [79.6537, 2.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [1.0], [1.0], [0.0]],
                                                                [30.9277, 1.5, -0.5]))
    MOIT.linear1test(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.0], [0.0], [0.0]],
                                                                [1.0]))
    MOIT.linear2test(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[3.0], [0.0]],
                                                                [-1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [3.0]],
                                                                [-0.0]))
    MOIT.linear3test(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0]],
                                                                [-0.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0]],
                                                                [-0.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0]],
                                                                [-0.0]))
    MOIT.linear4test(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.3333], [1.3333], [0.0], [0.0]],
                                                                [0.3333, 0.3333]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[2.0], [0.0], [0.0], [2.0]],
                                                                [0.5, -0.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[4.0], [0.0], [0.0]],
                                                                [1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[2.0], [0.0]],
                                                                [0.5]))
    MOIT.linear5test(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0], [0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0]],
                                                                [-1.0, 1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0], [168.953 0.0; 0.0 68.9528], [107.78 0.0; 0.0 107.78]],
                                                                [-1.0, 1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0], [250.152 0.0; 0.0 150.152], [150.152 0.0; 0.0 250.152]],
                                                                [-1.0, 1.0]))
    MOIT.linear6test(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0], [0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0]],
                                                                [-1.0, 1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0], [168.953 0.0; 0.0 68.9528], [107.78 0.0; 0.0 107.78]],
                                                                [-1.0, 1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0], [250.152 0.0; 0.0 150.152], [150.152 0.0; 0.0 250.152]],
                                                                [-1.0, 1.0]))
    MOIT.linear7test(cached_mock_optimizer, config)

    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                tuple(),
                                                                [1.0]))
    MOIT.linear8atest(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                (MOI.InfeasibilityCertificate, [[0.7709], [0.2291], [0.3126]])))
    MOIT.linear8btest(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                (MOI.InfeasibilityCertificate, [[0.5], [0.5]])))
    MOIT.linear8ctest(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[29.0909], [36.3636], [4.5455], [0.0], [0.0001]],
                                                                [-0.0, 11.3636, 0.8636]))
    MOIT.linear9test(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[5.0], [5.0], [5.0], [0.0]],
                                                                [-0.0, 1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[2.5], [2.5], [0.0], [5.0]],
                                                                [-1.0, 0.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.0], [1.0], [0.0], [10.0]],
                                                                [-1.0, -0.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[6.0], [6.0], [10.0], [0.0]],
                                                                [-0.0, 1.0]))
    MOIT.linear10test(MOIB.SplitInterval{Float64}(cached_mock_optimizer), config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.0], [0.0], [3.2439 0.0; 0.0 2.2439], [3.2439 0.0; 0.0 2.2439]],
                                                                [0.0, -1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [1.0], [1.3765 0.0; 0.0 0.8765], [1.3765 0.0; 0.0 0.8765]],
                                                                [-1.0, 0.0]))
    MOIT.linear11test(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                tuple(),
                                                                [1.0, 3.0]))
    MOIT.linear12test(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[24.0266], [27.6758 0.0; 0.0 22.6705], [27.6758 0.0; 0.0 22.6705]],
                                                                [0.0, 0.0]))
    MOIT.linear13test(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.5], [1.0], [0.0], [0.0]],
                                                                [2.0, 1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.0], [0.0]],
                                                                [1.0]))
    MOIT.linear14test(cached_mock_optimizer, config)
end
