# (General) birth death processes

This is a module with some functionalities related to (general) birth-death processes. A birth-death process `X(t)` is a continuous-time Markov process on the non-negative integers where the only allowed transitions are the birth of a new 'particle' and 'death' of an extant particle. The transition probabilities are defined as `P(X(t+Δt) = i+1|X(t)=i) = λᵢΔt + o(Δt)` and `P(X(t+Δt)=i-1|X(t)=i) = μᵢΔt + o(Δt)`. An important special case is where `λᵢ = iλ` and `μᵢ = iμ`, in which case the process is called a linear birth-death process.

Currently the module contains functions for working with linear and general birth-death processes. For the latter, the approach of Crawford & Suchard (2012) is used to compute the Laplace transform of the transition probabilities, using continued fractions, which is numerically inverted using `InverseLaplace.jl`.

The methods for general BDPs are not yet very efficiently implemented.

## Example

```julia
using BirthDeathProcesses
```

Linear BDP

```julia
lbdp = LinearBDP(0.3, 0.2)
tp(lbdp, 0.2, 6, 9)  # transition probability
```

```
0.0059196166724168565
```

Transient distribution for the linear BDP for `t=2` and `X₀=2`

```julia
tr = Transient(LinearBDP(0.3, 0.2), 1.2, 2)
pdf.(tr, 0:10)  # probabilities to end in states 0,1,...,10
pgf(tr, 0.)  # the extinction probability == pgf evaluated at 0
rand(tr)
```

```
7
```

General BDP

```julia
λ = (k)->0.3k^2*exp(-1.2*k)
μ = (k)->0.3k
gbdp = GeneralBDP(λ, μ)
@time tp(gbdp, 19, 27, 1.0)
```

```
1.044202164188786e-84
```

as noted above, the general BDPs need some work to get more efficient (and robust). The papers from Crawford & Suchard (2012); Crawford, Minin & Suchard (2014) and Ho et al. (2018) should be helpful.

```julia
using Literate
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

