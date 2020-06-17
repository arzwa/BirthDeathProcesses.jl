using Test, Random, BirthDeathProcesses, Lazy
import BirthDeathProcesses: ContinuedFraction, lentz, convergents
Random.seed!(19081994)

@testset begin "Continued fraction for eˣ"
    afun(x) = @>> 1:-Lazy.range(1) map(y->y*x)
    bfun(x) = 1.:@>> Lazy.range(2) map(y->y+x)
    gcf_exp(x) = ContinuedFraction(1., afun(x), bfun(x));
    for x in randn(10)
        @test isapprox(lentz(gcf_exp(x)), exp(x), atol=1e-6)
    end
end

@testset begin "Continued fraction for tan(x)"
    afun(x) = x:@>> constantly(-1.) map(y->y*x^2)
    bfun(x) = @>> Lazy.range(1) map(y->y*2-1)
    gcf_tan(x) = ContinuedFraction(0., afun(x), bfun(x))
    for x in randn(10)
        @test isapprox(lentz(gcf_tan(x)), tan(x))
        A, B = convergents(gcf_tan(x), 10)[end,:]
        @test isapprox(tan(x), A/B, atol=1e-6)
    end
end

@testset begin "Transition probabilities"
    gbdp = GeneralBDP(0.2, 0.3)
    lbdp = LinearBDP(0.2, 0.3)
    for (i, j) in zip(rand(1:20, 10),rand(1:20, 10))
        t = exp(randn())
        ilt = tp(gbdp, i, j)
        @test isapprox(real(ilt(t)), tp(lbdp, i, j, t), atol=1e-4)
    end
end
