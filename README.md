[![Build Status](https://travis-ci.com/arzwa/BirthDeathProcesses.jl.svg?branch=master)](https://travis-ci.com/arzwa/BirthDeathProcesses.jl)

# (General) birth death processes

This is a module with some functionalities related to (general) birth-death processes. A birth-death process `X(t)` is a continuous-time Markov process on the non-negative integers where the only allowed transitions are the birth of a new 'particle' and 'death' of an extant particle. The transition probabilities are defined as `P(X(t+Δt) = i+1|X(t)=i) = λᵢΔt + o(Δt)` and `P(X(t+Δt)=i-1|X(t)=i) = μᵢΔt + o(Δt)`. An important special case is where `λᵢ = iλ` and `μᵢ = iμ`, in which case the process is called a linear birth-death process.

Currently the module contains functions for working with linear and general birth-death processes. For the latter, the approach of Crawford & Suchard (2012) is used to compute the Laplace transform of the transition probabilities, using continued fractions, which is numerically inverted using `InverseLaplace.jl`.

The methods for general BDPs are not yet very efficiently implemented.

## Example

```julia
using BirthDeathProcesses
```

### Linear BDP

```julia
lbdp = LinearBDP(0.3, 0.2)
tp(lbdp, 0.2, 6, 9)  # transition probability
```

```
0.0059196166724168565
```

Transient distribution for the linear BDP for `t=2` and `X₀=2`

```julia
tbdp = lbdp(2., 2)
```

```
BirthDeathProcesses.TransientLinearBDP{Float64}(λ=0.3, μ=0.2, a=2, t=2.0, α=0.26607578096471346, β=0.3991136714470701)
```

This acts like a distribution from `Distributions.jl`

```julia
pdf.(tbdp, 0:10)  # probabilities to end in states 0,1,...,10
```

```
11-element Array{Float64,1}:
 0.07079632121598217
 0.23468151522042963
 0.2881500371272334
 0.1926264156371718
 0.1078597561048541
 0.05541281292948958
 0.02705085606317386
 0.012765930525135248
 0.005881137338864456
 0.00266097756552092
 0.0011872485532721843
```

The probability generating function and a sampler are implemented for the linear BDP

```julia
pgf(tbdp, 0.)  # the extinction probability == pgf evaluated at 0
rand(tbdp)
```

```
2
```

### General BDP

There are two implementations, `GeneralBDP` is a fast implementation which uses a bound on the continued fraction (default at 100 terms) and uses Arrays as data structures:

```julia
λ = (k)->0.3k^2*exp(-1.2*k)
μ = (k)->0.3k
gbdp = GeneralBDP(λ, μ)
tp(gbdp, 1.0, 19, 27)
```

```
1.0442021641887911e-84
```

If one exceeds the depth of the BDP, an error is raised

```julia
try tp(gbdp, 1.0, 101, 101); catch ex @show ex; end
```

```
AssertionError("(101, 101) exceeds BDP depth (100)")
```

So we should increase the BDP depth

```julia
gbdp = GeneralBDP(λ, μ, 200)
@time tp(gbdp, 1.0, 101, 102)
```

```
6.611263593795198e-63
```

There is another implementation using Lazy lists, which is more elegant but vastly less efficient

```julia
import BirthDeathProcesses: LazyGeneralBDP
gbdp = LazyGeneralBDP(λ, μ)
@time tp(gbdp, 1.0, 101, 102)
```

```
6.1535547441041e-63
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

