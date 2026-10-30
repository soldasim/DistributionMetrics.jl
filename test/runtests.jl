using DistributionMetrics
using Test
using Aqua
using Distributions
using KernelFunctions
using OptimizationPRIMA
using Random
using Statistics

@testset "DistributionMetrics.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(DistributionMetrics)
    end

    @testset "TVMetric" begin
        xs = collect(range(-6, 6; length=2001))'
        log_ws = fill(log(12.), size(xs, 2))
        logpdf_(d) = x -> logpdf(d, only(x))

        p = Normal(0, 1)
        q = Normal(1, 1)
        tv_true = 2 * cdf(Normal(), 0.5) - 1

        tv = TVMetric(; grid=xs, log_ws, true_logpost=logpdf_(p))
        @test calculate_metric(tv, logpdf_(p), logpdf_(q)) ≈ tv_true atol=1e-3
        @test calculate_metric(tv, logpdf_(p), logpdf_(p)) ≈ 0 atol=1e-12

        @testset "symmetry" begin
            tv_nocache = TVMetric(; grid=xs, log_ws)
            @test calculate_metric(tv_nocache, logpdf_(p), logpdf_(q)) ≈ calculate_metric(tv_nocache, logpdf_(q), logpdf_(p)) atol=1e-10
        end

        @testset "unnormalized densities" begin
            shifted(f, c) = x -> f(x) + c
            @test calculate_metric(tv, logpdf_(p), shifted(logpdf_(q), 7.)) ≈ calculate_metric(tv, logpdf_(p), logpdf_(q)) atol=1e-10
        end

        @testset "construction" begin
            tv_ws = TVMetric(; grid=xs, ws=fill(12., size(xs, 2)), true_logpost=logpdf_(p))
            @test tv_ws.log_ws ≈ tv.log_ws
            @test tv.true_logvals ≈ logpdf_(p).(eachcol(xs))
            @test_throws AssertionError TVMetric(; grid=xs)
            @test_throws AssertionError TVMetric(; grid=xs, ws=ones(2001), log_ws=zeros(2001))
        end

        @testset "cached true values" begin
            fail(_) = error("true log-density must not be evaluated")
            @test calculate_metric(tv, fail, logpdf_(q)) ≈ tv_true atol=1e-3

            tv_nocache = TVMetric(; grid=xs, log_ws)
            @test tv_nocache.true_logvals === nothing
            @test calculate_metric(tv_nocache, logpdf_(p), logpdf_(q)) ≈ tv_true atol=1e-3
        end

        @testset "non-normalizable approximation" begin
            @test calculate_metric(tv, logpdf_(p), x -> -Inf) ≈ 0.5 atol=1e-10
        end
    end

    @testset "MMDMetric" begin
        rng = MersenneTwister(1)
        kernel = with_lengthscale(SqExponentialKernel(), [1., 1.])
        mmd = MMDMetric(; kernel)

        X = randn(rng, 2, 200)
        scores = [calculate_metric(mmd, X, X .+ s) for s in [0.5, 1., 2.]]

        @test calculate_metric(mmd, X, X) ≈ 0 atol=1e-12
        @test all(>(0), scores)
        @test issorted(scores)
        @test calculate_metric(mmd, X, X .+ 1.) ≈ calculate_metric(mmd, X .+ 1., X)

        @testset "explicit formula" begin
            k(a, b) = exp(-sum(abs2, a .- b) / 2)
            A = [0. 1.; 0. 0.]
            B = reshape([2., 1.], 2, 1)
            expected = mean(k(a, a′) for a in eachcol(A), a′ in eachcol(A)) +
                mean(k(b, b′) for b in eachcol(B), b′ in eachcol(B)) -
                2 * mean(k(a, b) for a in eachcol(A), b in eachcol(B))
            @test calculate_metric(mmd, A, B) ≈ expected
        end

        @test_logs calculate_metric(mmd, X, X .+ 1.)
    end

    @testset "OptMMDMetric" begin
        rng = MersenneTwister(2)
        X = randn(rng, 2, 100)
        Y = randn(rng, 2, 100) .+ 0.7
        bounds = ([-3., -3.], [3., 3.])

        opt = OptMMDMetric(; kernel=SqExponentialKernel(), bounds, algorithm=BOBYQA())
        fixed = MMDMetric(; kernel=with_lengthscale(SqExponentialKernel(), (bounds[2] .- bounds[1]) ./ 3))

        @test calculate_metric(opt, X, Y) >= calculate_metric(fixed, X, Y) - 1e-8
        @test calculate_metric(opt, X, X) ≈ 0 atol=1e-12

        opt_kw = OptMMDMetric(; kernel=SqExponentialKernel(), bounds, algorithm=BOBYQA(), maxiters=10)
        @test haskey(opt_kw.kwargs, :maxiters)
        @test calculate_metric(opt_kw, X, Y) >= 0

        @test_logs calculate_metric(opt, X, Y)
        @test_logs (:info, r"lengthscales") calculate_metric(opt, X, Y; info=true)
    end

    @testset "interface" begin
        @test TVMetric <: PDFMetric <: DistributionMetric
        @test MMDMetric <: SampleMetric <: DistributionMetric
        @test OptMMDMetric <: SampleMetric
    end
end
