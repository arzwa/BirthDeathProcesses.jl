using Test, Random, BirthDeathProcesses, Lazy
import BirthDeathProcesses: ContinuedFraction, lentz, convergents
Random.seed!(19081994)

@testset "BDP.jl" begin
    @testset "Continued fraction for eˣ" begin
        afun(x) = @>> 1:-Lazy.range(1) map(y->y*x)
        bfun(x) = 1.:@>> Lazy.range(2) map(y->y+x)
        gcf_exp(x) = ContinuedFraction(1., afun(x), bfun(x));
        for x in randn(10)
            @test isapprox(lentz(gcf_exp(x)), exp(x), atol=1e-6)
        end
    end

    @testset "Continued fraction for tan(x)" begin
        afun(x) = x:@>> constantly(-1.) map(y->y*x^2)
        bfun(x) = @>> Lazy.range(1) map(y->y*2-1)
        gcf_tan(x) = ContinuedFraction(0., afun(x), bfun(x))
        for x in randn(10)
            @test isapprox(lentz(gcf_tan(x)), tan(x))
            A, B = convergents(gcf_tan(x), 10)[end,:]
            @test isapprox(tan(x), A/B, atol=1e-6)
        end
    end

    @testset "Transition probabilities" begin
        for i=1:20
            λ = exp(randn())
            μ = exp(log(λ) + randn()/2)
            gbdp = GeneralBDP(λ, μ)
            lbdp = LinearBDP(λ, μ)
            for (i, j) in zip(rand(0:20, 10),rand(0:20, 10))
                t = exp(randn())
                ilt = tp(gbdp, i, j)
                @test isapprox(real(ilt(t)), tp(lbdp, t, i, j), atol=1e-3)
            end
        end
    end

    @testset "LinearBDP transient distribution" begin
        for i=1:100
            λ = exp(randn())
            μ = exp(log(λ) + randn()/2)
            d = Transient(LinearBDP(λ, μ), rand(), 1)
            @test isapprox(sum(pdf.(d, 0:100)), one(λ), atol=1e-3)
        end
    end
end
