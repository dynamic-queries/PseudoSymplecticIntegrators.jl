abstract type AbstractFunction end

struct SimpleFunction <: AbstractFunction 
    f::Function
end

struct PartionedFunction <: AbstractFunction
    f1::Function
    f2::Function
end

struct ODEProblem{StateType}
    f::Union{SimpleFunction, PartionedFunction}
    u0::StateType
    tspan::Tuple
    t::Union{StepRange, StepRangeLen, Vector}
    params::Tuple
    
    function ODEProblem(f, tspan, params, t)
        new(f,tspan,t,params)
    end
end

struct ODESolution{StateType}
    prob::ODEProblem
    ts::Vector{Float64}
    us::Vector{StateType}
    tquery::Union{StepRange, StepRangeLen, Vector}
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
    u_interp = interpolate(us,ts,t)
    return ODESolution(prob, ts, us, t, u_interp)
end

# Use Hermite interpolation
function interpolate(u::Vector{Any}, t::Vector{Float64}, tquery::Union{StepRange, StepRangeLen, Vector})
    
end