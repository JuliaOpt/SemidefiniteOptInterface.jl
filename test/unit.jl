@testset "Unit" begin
    MOIU.set_mock_optimize!(mock, (mock) -> MOIU.mock_optimize!(mock,
                                                                [[0.0]],
                                                                [0.0]))
    MOIT.solve_singlevariable_obj(cached_mock_optimizer, config)
end
