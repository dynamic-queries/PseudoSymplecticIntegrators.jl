@testset "Interface AubChar" begin
    dt = 1e-3
    solver = AubChar(dt)
    @test solver.tableau.A[5,4] == -1.0025739247772099238420
    @test solver.tableau.b[5] == 0.130002156105755533849
    @test solver.tableau.c[5] == 1
    print(solver.info)
end

@testset "ODEProblem Interface scalar case" begin
    function du(u,p,t)
        du = -p[1]*u[1]
        return [du]
    end

    tspan = (0,10)
    params = [0.1]
    t = LinRange(tspan[1],tspan[2],100)
    dt = 0.5
    u0 = [1.0]
    func = SimpleFunction(du)
    prob = ODEProblem(func, u0, tspan, params, t)
    @test prob.f  == func
    @test prob.u0 == u0
    @test prob.tspan == tspan
    @test prob.params == params
    # @show size(du(u0,params,t))
end

@testset "ODEProblem Interface vector case" begin
    A = rand(3,3)
    function du(u,p,t)
        du = p[1] .* A*u
        return du
    end 
    tspan = (0,10)
    params = [0.1]
    t = LinRange(tspan[1],tspan[2],100)
    dt = 0.5
    u0 = [1.0,1.0,1.0]
    func = SimpleFunction(du)
    prob = ODEProblem(func, u0, tspan, params, t)
    @test prob.f  == func
    @test prob.u0 == u0
    @test prob.tspan == tspan
    @test prob.params == params
    # @show size(du(u0,params,t))
end

@testset "ODEProblem Interface matrix case" begin
    A = rand(3,3)
    function du(u,p,t)
        du = p[1] .* A*u
        return du
    end
    tspan = (0,10)
    params = [0.1]
    t = LinRange(tspan[1],tspan[2],100)
    dt = 0.5
    u0 = rand(3,3)
    func = SimpleFunction(du)
    prob = ODEProblem(func, u0, tspan, params, t)
    @test prob.f  == func
    @test prob.u0 == u0
    @test prob.tspan == tspan
    @test prob.params == params
    # @show size(du(u0,params,t))
end