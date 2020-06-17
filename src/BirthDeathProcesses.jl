module BirthDeathProcesses
    using InverseLaplace
    using Lazy
    export GeneralBDP, LinearBDP, tp

    include("gcf.jl")
    include("bdp.jl")

end # module
