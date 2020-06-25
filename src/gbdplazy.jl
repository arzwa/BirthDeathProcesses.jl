# General BDPs using Lazy lists for continued fractions
# Lazy stuff is elegant IMHO, but inefficient
# doing it.
"""
    LazyGeneralBDP(λ::Function, μ::Function)

A general birth-death process. `λ` and `μ` should be provided as
functions that take an integer `k` and return `λₖ` (or `μₖ`) i.e.
the birth (death) rate when the total population is of size `k`.
"""
struct LazyGeneralBDP{T} <: BDP{T}
    λ::Function
    μ::Function
    LazyGeneralBDP(f, g) = new{typeof(f(1))}(f, g)
end

LazyGeneralBDP(λ::T, μ::T) where T<:Real = LazyGeneralBDP(k->k*λ, k->k*μ)

# NOTE μ(0) should be defined as 0.
_a(p::LazyGeneralBDP) = 1.:@>> Lazy.range(2) map(k->-p.λ(k-2)*p.μ(k-1))
_b(p::LazyGeneralBDP, s) = @>> Lazy.range(1) map(k->s + p.λ(k-1) + p.μ(k-1))

gcf_f00(p::LazyGeneralBDP, s::T) where T = LazyContinuedFraction(0., _a(p), _b(p, s))

# Laplace transform of transition probability i-> j
function tf(p::LazyGeneralBDP{T}, i, j, s) where T
    if j <= i
        factor1 = prod([p.μ(k) for k=j+1:i])
    else
        factor1 = prod([p.λ(k) for k=i:j-1])
        j_ = j; j = i; i = j_
    end
    f00 = gcf_f00(p, s)
    i == j == 0 && return lentz(f00)
    AB = convergents(f00, i+1)
    as = (j == 0 ? 1. : AB[j,2]) : AB[i,2]*f00.a[i+2] : drop(i+2, f00.a)
    bs = AB[i+1,2] : drop(i+1, f00.b)
    gcf = LazyContinuedFraction(zero(T), as, bs)
    return factor1 * lentz(gcf)
end

function tp(p::LazyGeneralBDP, i, j)
    ilt = ILT(s->tf(p, i, j, s))
    (t)->real(ilt(t))
end

tp(p::LazyGeneralBDP, t, i, j) = real(tp(p, i, j)(t))

# implement a Transient struct for the General BDP, which should implement more
# efficient computation of transition probabilities, because the laplace
# transforms are related recursively? I mean, for a single `s` we should be able
# to more economically calculate the different `fₘₙ` values?

# Also, in applications one generally needs a transition probability matrix on some subspace? (note that in phylogenetic applications we need to use a pruning algorithm I think, since the algorithm of Csuros & Milkos is tied to the linear BDP, i.e. independent evolution of lineages is required for the recursions to hold)
