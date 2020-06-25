# [![Build Status](https://travis-ci.com/arzwa/BirthDeathProcesses.jl.svg?branch=master)](https://travis-ci.com/arzwa/BirthDeathProcesses.jl)

# # (General) birth death processes

# This is a module with some functionalities related to (general) birth-death processes. A birth-death process `X(t)` is a continuous-time Markov process on the non-negative integers where the only allowed transitions are the birth of a new 'particle' and 'death' of an extant particle. The transition probabilities are defined as `P(X(t+Δt) = i+1|X(t)=i) = λᵢΔt + o(Δt)` and `P(X(t+Δt)=i-1|X(t)=i) = μᵢΔt + o(Δt)`. An important special case is where `λᵢ = iλ` and `μᵢ = iμ`, in which case the process is called a linear birth-death process.

# Currently the module contains functions for working with linear and general birth-death processes. For the latter, the approach of Crawford & Suchard (2012) is used to compute the Laplace transform of the transition probabilities, using continued fractions, which is numerically inverted using `InverseLaplace.jl`.

# The methods for general BDPs are not yet very efficiently implemented.

# ## Example
using BirthDeathProcesses

# ### Linear BDP
lbdp = LinearBDP(0.3, 0.2)
tp(lbdp, 0.2, 6, 9)  # transition probability

# Transient distribution for the linear BDP for `t=2` and `X₀=2`
tbdp = lbdp(2., 2)

# This acts like a distribution from `Distributions.jl`
pdf.(tbdp, 0:10)  # probabilities to end in states 0,1,...,10

# The probability generating function and a sampler are implemented for the linear BDP
pgf(tbdp, 0.)  # the extinction probability == pgf evaluated at 0
rand(tbdp)

# ### General BDP

# There are two implementations, `GeneralBDP` is a fast implementation which uses a bound on the continued fraction (default at 100 terms) and uses Arrays as data structures:
λ = (k)->0.3k^2*exp(-1.2*k)
μ = (k)->0.3k
gbdp = GeneralBDP(λ, μ)
tp(gbdp, 1.0, 19, 27)

# If one exceeds the depth of the BDP, an error is raised
try tp(gbdp, 1.0, 101, 101); catch ex @show ex; end

# So we should increase the BDP depth
gbdp = GeneralBDP(λ, μ, 200)
@time tp(gbdp, 1.0, 101, 102)

# There is another implementation using Lazy lists, which is more elegant but vastly less efficient
import BirthDeathProcesses: LazyGeneralBDP
gbdp = LazyGeneralBDP(λ, μ)
@time tp(gbdp, 1.0, 101, 102)


using Literate #src
Literate.markdown(@__FILE__, joinpath(@__DIR__, ".."), execute=true, documenter=false) #src
