# Linear birth-death process.
# The linear BDP has a nice closed form pgf and a shifted geometric transient
# distribution.
const ΛMATOL = 1e-6

struct LinearBDP{T} <: BDP{T}
    λ::T
    μ::T
end

getαβ(p::LinearBDP, t) = getαβ(p.λ, p.μ, t)
function getαβ(λ, μ, t)
    α = getα(λ, μ, t)
    α, (λ/μ)*α
end

function getα(λ, μ, t)
    x = t*(λ - μ)
    isapprox(x, zero(x), atol=ΛMATOL) ?
        λ*t/(one(t)+λ*t) :
        μ*(exp(x)-one(t))/(λ*exp(x)-μ)
end

# # transition probability for the linear BDP
function tp(p::LinearBDP, t, a, b)
    a == b == zero(a) && return one(t)
    α, β = getαβ(p, t)
    tp(a, b, α, β)
end

tp(a::Int, b::Int, α, β) = b == zero(b) ? α^a :
    mapreduce(k->_ξ(a, b, k, α, β), +, 0:min(a,b))

# logarithm of the factorial using Stirlings approximation
_stir(n) = n*log(n) - n + log(2π*n)/2

# non-overflowing binomial
_bin(n, k) =
    k == 0 ? 1. :
    n <= 60 ? float(binomial(n, k)) :
    exp(_stir(n)-_stir(k)-_stir(n - k))

# helper function for BDP transition probability
_ξ(i, j, k, α, β) = _bin(i, k)*_bin(i+j-k-1,i-1)*α^(i-k)*β^(j-k)*(one(α)-α-β)^k

"""
    TransientLinearBDP(process, t, [a=1])

Construct the transient distribution for the process at time `t`
with initial value `a`.
"""
struct TransientLinearBDP{T} <: Transient
    λ::T
    μ::T
    a::Int64
    t::T
    α::T
    β::T
    TransientLinearBDP(λ::T, μ::T, t::T, a=1) where T =
        new{T}(λ, μ, a, t, getαβ(λ, μ, t)...)
end

(p::LinearBDP)(t, a=1) = TransientLinearBDP(p.λ, p.μ, t, a)

# sum of `a` shifted geometric rv's
function Base.rand(d::TransientLinearBDP{T}) where T
    mapreduce((_)->rand() < d.α ? 0 : rand(Geometric(d.β)) + 1, +, 1:d.a)
end

function Distributions.pdf(d::TransientLinearBDP{T}, x::Int) where T
    x == d.a == zero(x) && return one(t)
    tp(d.a, x, d.α, d.β)
end

Distributions.logpdf(d::TransientLinearBDP, x) = log(pdf(d))

"""
    pgf(::Transient{LinearBDP}, s)

The probability generating function for the linear BDP.
"""
function pgf(d::TransientLinearBDP{T}, s) where T
    @unpack λ, μ, t, a = d
    x = t*(μ - λ)
    ρ = exp(x)
    f = isapprox(x, zero(x), atol=ΛMATOL) ?
        (one(T) - (λ*t - one(T))*(s - one(T)))/(one(T) - λ*t*(s - one(T))) :
        (ρ*(λ*s - μ) - μ*(s-one(T)))/(ρ*(λ*s - μ) - λ*(s-one(T)))
    f^a
end
