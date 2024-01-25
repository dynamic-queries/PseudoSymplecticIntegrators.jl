module PseudoSymplecticIntegrators
    using SparseArrays
    using UnPack
    
    include("algorithms.jl")
    export AubChar
    
    include("iterator.jl")
end
