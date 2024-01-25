@testset "Interface AubChar" begin
    dt = 1e-3
    solver = AubChar(dt)
    @test solver.tableau.A[5,4] == -1.0025739247772099238420
    @test solver.tableau.b[5] == 0.130002156105755533849
    @test solver.tableau.c[5] == 1
    print(solver.info)
end 