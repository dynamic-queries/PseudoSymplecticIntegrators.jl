abstract type PseudoSymIntegrator end
abstract type Tableau end

struct AubChar <: PseudoSymIntegrator
    info::String
    nstages::Int
    tableau::Tableau
    dt::Float64
    
    function AubChar(dt)
        info = string(
        "---------------------------------------------\n",
        "AubChar: Explicit RK pseudosymplectic integrator \nIntegration Order 3\nSymplecticity order 6\n",
        "---------------------------------------------")
        nstages = 5
        tab = AubCharTableau(Float64)
        new(info, nstages, tab, dt)
    end 
end

struct AubCharTableau <: Tableau
    A::SparseMatrixCSC
    b::SparseVector
    c::SparseVector
end 

function AubCharTableau(type::Type)
    A = spzeros(type,5,5)
    B = spzeros(type,5)
    C = spzeros(type,5)
    
    A[2,1] = 0.13502027922908531468
    A[3,1] = -0.47268213605236986919
    A[3,2] = 1.05980250415418968199
    A[4,1] = -1.21650460595688538935
    A[4,2] = 2.16217630216752533012
    A[4,3] = -0.372345924265360030384
    A[5,1] = 0.3327444303638736757818
    A[5,2] = -0.2088266829658723128357
    A[5,3] = 1.8786561773792085608959
    A[5,4] = -1.0025739247772099238420
    
    B[1] = 0.04113894457091769183
    B[2] = 0.26732123194413937348
    B[3] = 0.86700906289954518480
    B[4] = -0.30547139552035758861
    B[5] = 0.130002156105755533849
    
    C[1] = 0.0
    C[2] = 0.13502027922908531468
    C[3] = 0.58712036810181981280
    C[4] = 0.57332577194527991038
    C[5] = 1
    return AubCharTableau(A,B,C)
end

function perform_step(ui, ti, p, f::SimpleFunction, solver::AubChar)
    @unpack dt, tableau, nstages = solver
    @unpack A,b,c = tableau
    k1 = f.f(ui,p,t)
    k2 = f.f(ui + dt*A[2,1]*k1,p,t+c[2]*dt)
    k3 = f.f(ui + dt*(A[3,1]*k1 + A[3,2]*k2),p,t+c[3]*dt)
    k4 = f.f(ui + dt*(A[4,1]*k1 + A[4,2]*k2 + A[4,3]*k3),p,t+c[4]*dt)
    k5 = f.f(ui + dt*(A[5,1]*k1 + A[5,2]*k2 + A[5,3]*k3 + A[5,4]*k4),p,t+c[5]*dt)
    ui1 = b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4 + b[5]*k5
    return ui1
end

function split(v::Vector)
    n = length(v)
    l = floor(Int,n/2)
    v[1:l], v[l+1:end]
end

function split(m::Matrix)
    n = size(m,1)
    l = floor(Int, n/2)
    m[1:l,:], m[l+1:end,:]
end 

function perform_step(ui, ti, p, f::PartionedFunction, solver::AubChar)
    @unpack dt, tableau, nstages = solver
    @unpack A,b,c = tableau
    
    function g(x,t)
        y,z = split(x)
        vcat(f.f1(y,p,t), f.f2(z,p,t))
    end 
    
    k1 = g(ui,t+c[1]*dt)
    k2 = g(ui + dt*(A[2,1]*k1),t+c[2]*dt)
    k3 = f.g(ui + dt*(A[3,1]*k1 + A[3,2]*k2),t+c[3]*dt)
    k4 = f.g(ui + dt*(A[4,1]*k1 + A[4,2]*k2 + A[4,3]*k3),t+c[4]*dt)
    k5 = f.g(ui + dt*(A[5,1]*k1 + A[5,2]*k2 + A[5,3]*k3 + A[5,4]*k4),t+c[5]*dt)
    ui1 = b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4 + b[5]*k5
    return ui1
end