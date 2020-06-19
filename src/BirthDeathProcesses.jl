module BirthDeathProcesses  # AKA BDP
    using InverseLaplace
    using Lazy
    using Distributions
    using Parameters
    import Distributions: pdf, logpdf
    export GeneralBDP, LinearBDP, Transient, tp, pgf, pdf, logpdf

    abstract type BDP{T} end
    abstract type Transient <: DiscreteUnivariateDistribution end

    include("gcf.jl")     # continued fractions
    include("gbdp.jl")    # general BDP
    include("linear.jl")  # linear BDP
end # module
