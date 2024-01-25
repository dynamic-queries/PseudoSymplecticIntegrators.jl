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

function perform_step(ui, ti, p, f, solver::AubChar)
    @unpack dt, tableau, nstages = solver
    k1 = f()
end 