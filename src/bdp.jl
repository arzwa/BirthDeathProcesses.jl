# General BDPs using Lazy lists ad continued fractions
# Very elegant, but not sure about efficiency of the Lazy list way of doing it
"""
    GeneralBDP(λ::Function, μ::Function)
"""
struct GeneralBDP{T}
    λ::Function
    μ::Function
    GeneralBDP(f, g) = new{typeof(f(1))}(f, g)
end

GeneralBDP(λ::T, μ::T) where T<:Real = GeneralBDP(k->k*λ, k->k*μ)

# NOTE μ(0) should be defined as 0.
_a(p::GeneralBDP) = 1.:@>> Lazy.range(2) map(k->-p.λ(k-2)*p.μ(k-1))
_b(p::GeneralBDP, s) = @>> Lazy.range(1) map(k->s + p.λ(k-1) + p.μ(k-1))

gcf_f00(p::GeneralBDP, s::T) where T = ContinuedFraction(0., _a(p), _b(p, s))

# Laplace transform of transition probability i-> j
function tf(p::GeneralBDP{T}, i, j, s) where T
    if j <= i
        factor1 = prod([p.μ(k) for k=j+1:i])
    else
        factor1 = prod([p.λ(k) for k=i:j-1])
        j_ = j; j = i; i = j_
    end
    f00 = gcf_f00(p, s)
    AB = convergents(f00, i+1)
    as = AB[j,2] : AB[i,2]*f00.a[i+2] : drop(i+2, f00.a)
    bs = AB[i+1,2] : drop(i+1, f00.b)
    gcf = ContinuedFraction(zero(T), as, bs)
    return factor1 * lentz(gcf)
end

function tp(p::GeneralBDP, i, j)
    ilt = ILT(s->tf(p, i, j, s))
    (t)->real(ilt(t))
end

tp(p::GeneralBDP, i, j, t) = tp(p, i, j)(t)

"""
    LinearBDP(λ, μ)

Linear (Kendall) Birth-Death process (λₖ = kλ, μₖ = kμ).
"""
struct LinearBDP{T}
    λ::T
    μ::T
end

LinearBDP(λ, μ) = LinearBDP(promote(λ, μ)...)

tp(p::LinearBDP, i, j, t) = (i == j == 0) ? 1.0 : plinbdp(i, j, p.λ, p.μ, t)

# transition probabilities for linear bdp
ϕ(t, λ, μ) = λ ≈ μ ? λ*t/(1 + λ*t) : μ*(exp(t*(λ-μ)) - 1)/(λ*exp(t*(λ-μ)) - μ)
ψ(t, λ, μ) = λ ≈ μ ? λ*t/(1 + λ*t) : (λ/μ)*ϕ(t, λ, μ)
ξ(i, j, k, t, λ, μ) = binomial(i, k)*binomial(i+j-k-1,i-1)*
    ϕ(t, λ, μ)^(i-k)*ψ(t, λ, μ)^(j-k)*(1-ϕ(t, λ, μ)-ψ(t, λ, μ))^k
plinbdp(i, j, λ, μ, t) = sum([ξ(i, j, k, t, λ, μ) for k = 0:min(i,j)])
