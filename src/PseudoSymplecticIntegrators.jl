module PseudoSymplecticIntegrators
    using SparseArrays
    using UnPack
    
    include("types.jl")
    export SimpleFunction, PartionedFunction

    include("algorithms.jl")
    export AubChar
    
    include("iterator.jl")
    export ODEProblem
end
