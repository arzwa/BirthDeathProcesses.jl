
struct ContinuedFraction{T,V}
    b₀::T
    a ::Vector{V}
    b ::Vector{T}
end

depth(cf::ContinuedFraction) = length(cf.a)

Base.show(io::IO, gcf::ContinuedFraction{T}) where T =
    write(io, "GCF{$T} truncated at k=$(depth(gcf)): $(terms_string(gcf, 2))")

function terms_string(gcf::ContinuedFraction, n)
    str = "($(gcf.b₀) + "
    for i=1:n
        str *= "$(gcf.a[i]) / ($(gcf.b[i]) + "
    end
    str * " … " * ")"^(n+1)
end

ContinuedFraction(cf::LazyContinuedFraction{T}, depth=100) where T =
    ContinuedFraction{T}(cf.b₀,
        collect(take(cf.a, depth)),
        collect(take(cf.b, depth)))

"""
    lentz(c::ContinuedFraction)

Modified Lentz method, as in section 5.2 of 'Numerical recipes in C'.
"""
function lentz(gcf::ContinuedFraction{T}, ϵ=10^-7, nmax=depth(gcf)) where T
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
    abs(Δ-one(Δ)) > ϵ &&
        @warn "Reached `nmax=$nmax` in `lentz` ($(abs(Δ-one(Δ))) > ϵ)"
    return f
end

"""
    convergents(gcf::ContinuedFraction, k::Int)

Compute the convergents up to index `k`. Returns a (k × 2) matrix with
the first and second columns the numerator and denominator respectively.
"""
function convergents(gcf::ContinuedFraction{T}, k::Int) where T
    if k > depth(gcf)
        @warn "Convergent $i > depth $(depth(gcf))"
        k = depth(gcf)
    end
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
