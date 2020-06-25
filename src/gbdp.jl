# General BDPs using arrays for continued fractions
"""
    GeneralBDP(λ::Function, μ::Function)

A general birth-death process. `λ` and `μ` should be provided as
functions that take an integer `k` and return `λₖ` (or `μₖ`) i.e.
the birth (death) rate when the total population is of size `k`.
"""
struct GeneralBDP{T} <: BDP{T}
    fλ ::Function
    fμ ::Function
    λ  ::Vector{T}
    μ  ::Vector{T}
    a  ::Vector{T}
    b  ::Vector{T}
end

function GeneralBDP(l, m, depth=100)
    λ = l.(0:depth)
    μ = m.(0:depth)
    a, b = getab(λ, μ)
    GeneralBDP(l, m, λ, μ, a, b)
end

depth(p::GeneralBDP) = length(p.λ) -1

GeneralBDP(λ::T, μ::T) where T<:Real = GeneralBDP(k->k*λ, k->k*μ)

# NOTE μ(0) should be defined as 0.
function getab(λ::Vector{T}, μ::Vector{T}) where T
    a = vcat(one(T), map(k->-λ[k-1]*μ[k], 2:length(μ)))
    b = λ .+ μ
    a, b
end

# _a(p::GeneralBDP) = 1.:@>> Lazy.range(2) map(k->-p.λ(k-2)*p.μ(k-1))
# _b(p::GeneralBDP, s) = @>> Lazy.range(1) map(k->s + p.λ(k-1) + p.μ(k-1))

# gcf_f00(p::GeneralBDP, s::T) where T = LazyContinuedFraction(0., _a(p), _b(p, s))

# Laplace transform of transition probability i-> j
function tf(p::GeneralBDP{T}, i, j, s::V) where {T,V}
    d = depth(p)
    @assert (i <= d || j <= d) "($i, $j) exceeds BDP depth ($d)"
    if j <= i
        factor1 = prod([p.μ[k+1] for k=j+1:i])
    else
        factor1 = prod([p.λ[k+1] for k=i:j-1])
        j_ = j; j = i; i = j_
    end
    f00 = ContinuedFraction(zero(V), p.a, p.b .+ s)
    i == j == 0 && return lentz(f00)
    AB = convergents(f00, i+1)
    a1 = j == 0 ? Complex(one(T)) : AB[j,2]
    as = vcat(a1, AB[i,2]*f00.a[i+2], f00.a[i+3:end])
    bs = vcat(AB[i+1,2], f00.b[i+2:end])
    gcf = ContinuedFraction(Complex(zero(T)), as, bs)
    return factor1 * lentz(gcf)
end

function tp(p::GeneralBDP, i, j)
    ilt = ILT(s->tf(p, i, j, s))
    (t)->real(ilt(t))
end

tp(p::GeneralBDP, t, i, j) = real(InverseLaplace.ilt(s->tf(p, i, j, s), t))

# implement a Transient struct for the General BDP, which should implement more
# efficient computation of transition probabilities, because the laplace
# transforms are related recursively? I mean, for a single `s` we should be able
# to more economically calculate the different `fₘₙ` values?

# Also, in applications one generally needs a transition probability matrix on some subspace? (note that in phylogenetic applications we need to use a pruning algorithm I think, since the algorithm of Csuros & Milkos is tied to the linear BDP, i.e. independent evolution of lineages is required for the recursions to hold)
