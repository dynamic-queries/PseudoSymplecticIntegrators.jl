struct ODEProblem{StateType}
    f::Union{SimpleFunction, PartionedFunction}
    u0::StateType
    tspan::Tuple
    t::Union{StepRange, LinRange, StepRangeLen, Vector}
    params::Union{Tuple,Vector,Float64}
    
    function ODEProblem(f, u0::Union{Vector,Matrix}, tspan, params, t)
        new{typeof(u0)}(f,u0,tspan,t,params)
    end
end

struct ODESolution{StateType}
    prob::ODEProblem
    ts::Vector{Float64}
    us::Vector{StateType}
    tquery::Union{StepRange, LinRange, StepRangeLen, Vector}
    sol::Vector{StateType}
end 

function solve(prob::ODEProblem, solver::PseudoSymIntegrator)
    @unpack u0, tspan, param, t, f = prob
    @unpack dt = solver
    us = [u0]
    ts = [tspan[1]]
    t_iter = tspan[1]
    while t_iter < tspan[2]
        u_iter = perform_step(last(us), last(ts), params, f, solver)
        t_iter += dt
        push!(us, u_iter)
        push!(ts, t_iter) 
    end
    if last(ts) != tspan[2]
        u_iter = perform_step(last(us), last(ts), params, f, solver)
        t_iter = last(ts) + dt
    end
    u_interp = interpolate(us,ts,f,p,t)
    return ODESolution(prob, ts, us, t, u_interp)
end

# Use Hermite interpolation
function interpolate(u::Vector{Any}, t::Vector{Float64}, f::Union{SimpleFunction, PartionedFunction}, params::Tuple,tquery::Union{StepRange, StepRangeLen, Vector})
    
end