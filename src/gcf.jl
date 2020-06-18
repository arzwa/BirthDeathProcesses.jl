# Generalized continued fractions
# We use lazy lists
const TINY = 10^-30
ifzerotiny(x::T) where T = iszero(x) ? T(TINY) : x

"""
    ContinuedFraction(b₀, a, b)

A generalized contiued fraction represented by Lazy lists.
    b₀ +       a₁
        ----------------
        b₁ +     a₂
            ------------
            b₂ +  a₃/...
"""
struct ContinuedFraction{T}
    b₀::T
    a ::List
    b ::List
end

Base.show(io::IO, gcf::ContinuedFraction{T}) where T =
    write(io, "GCF{$T}: $(terms_string(gcf, 2))")

function terms_string(gcf::ContinuedFraction, n)
    str = "($(gcf.b₀) + "
    for i=1:n
        str *= "$(gcf.a[i]) / ($(gcf.b[i]) + "
    end
    str * " … " * ")"^(n+1)
end

"""
    lentz(c::GeneralizedContinuedFraction)

Modified Lentz method, as in section 5.2 of 'Numerical recipes in C'.
"""
function lentz(gcf::ContinuedFraction{T}, ϵ=10^-7, nmax=1000) where T
    C = f = ifzerotiny(gcf.b₀)
    D = zero(C)
    Δ = Inf
    for i=1:nmax
        D = one(D)/ifzerotiny(gcf.b[i] + gcf.a[i]*D)
        C = ifzerotiny(gcf.b[i] + gcf.a[i]/C)
        Δ = C*D
        f *= Δ
        abs(Δ-one(Δ)) < ϵ && break
    end
    return f
end

# should be able to do it fully lazily
# FIXME
function _lentz(gcf::ContinuedFraction{T}, ϵ=10^-7, nmax=1000) where T
    C = f = ifzerotiny(gcf.b₀)
    Δ = Inf
    D = @lazy zero(T):(@>> (takelast(D, 1) * gcf.a + gcf.b) map(x->one(T)/ifzerotiny(x)))
    # D = @lazy zero(T):takelast(D, 1)
    take(D, 20)
end
#= benchmarks:
# (1) indexing lazy list naively:
# BenchmarkTools.Trial:
# memory estimate:  31.23 KiB
# allocs estimate:  1743
# --------------
# minimum time:     443.350 μs (0.00% GC)
# median time:      452.875 μs (0.00% GC)
# mean time:        494.379 μs (0.74% GC)
# maximum time:     7.698 ms (90.59% GC)
# --------------
# samples:          10000
# evals/sample:     1

# (2) take first nmax entries
# BenchmarkTools.Trial:
# memory estimate:  38.43 KiB
# allocs estimate:  2083
# --------------
# minimum time:     516.871 μs (0.00% GC)
# median time:      522.330 μs (0.00% GC)
# mean time:        544.990 μs (0.66% GC)
# maximum time:     6.233 ms (89.03% GC)
# --------------
# samples:          9137
# evals/sample:     1 =#

"""
    convergents(gcf::ContinuedFraction, k::Int)

Compute the convergents up to index `k`. Returns a (k × 2) matrix with
the first and second columns the numerator and denominator respectively.
"""
function convergents(gcf::ContinuedFraction{T}, k::Int) where T
    AB = zeros(Complex, k, 2)
    A = C = ifzerotiny(gcf.b₀)
    B = one(T)
    D = zero(T)
    for j=1:k
        D = one(D)/ifzerotiny(gcf.b[j] + gcf.a[j]*D)
        C = ifzerotiny(gcf.b[j] + gcf.a[j]/C)
        AB[j, 1] = A = A*C
        AB[j, 2] = B = B/D
    end
    return AB
end
