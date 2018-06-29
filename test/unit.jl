@testset "Unit" begin
    MOIT.add_variable(cached_mock_optimizer, config)
    MOIT.add_variables(cached_mock_optimizer, config)
    MOIT.delete_variable(cached_mock_optimizer, config)
    MOIT.delete_variable(cached_mock_optimizer, config)
    MOIT.feasibility_sense(cached_mock_optimizer, config)
    MOIT.getconstraint(cached_mock_optimizer, config)
    MOIT.getvariable(cached_mock_optimizer, config)
    MOIT.max_sense(cached_mock_optimizer, config)
    MOIT.min_sense(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [0.0 0; 0 0.0]],
                                                                [1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.0], [0.0], [0.624034 0; 0 0.624034]],
                                                                [0.0, 1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [2.73683 0; 0 1.73683]],
                                                                [1.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.88804e-9], [1.0], [3.25375 0; 0 2.25375]],
                                                                [1.0, 0.0]),
                                  (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [4.47626 0; 0 2.47626]],
                                                                [1.0]))
    MOIT.solve_affine_deletion_edge_cases(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.72037 0; 0 1.22037]],
                                                                [-0.5]))
    MOIT.solve_affine_equalto(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [2.10102 0; 0 1.60102]],
                                                                [-0.5]))
    MOIT.solve_affine_greaterthan(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[3.0], [0.0], [4.88931 0; 0 2.88931]],
                                                                [0.0, 1.5]))
    MOIT.solve_affine_interval(MOIB.SplitInterval{Float64}(cached_mock_optimizer), config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [2.06719 0; 0 1.56719]],
                                                                [0.5]))
    MOIT.solve_affine_lessthan(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0]],
                                                                [0.0]))
    MOIT.solve_constant_obj(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0]],
                                                                [0.0]))
    MOIT.solve_singlevariable_obj(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0], [1.0]],
                                                                [0.0]))
    MOIT.solve_with_lowerbound(cached_mock_optimizer, config)
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[1.0], [0.0]],
                                                                [2.0]))
    MOIT.solve_with_upperbound(cached_mock_optimizer, config)
    MOIT.variablenames(cached_mock_optimizer, config)
end
