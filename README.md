# DistributionMetrics

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://soldasim.github.io/DistributionMetrics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://soldasim.github.io/DistributionMetrics.jl/dev/)
[![Build Status](https://github.com/soldasim/DistributionMetrics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/soldasim/DistributionMetrics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/soldasim/DistributionMetrics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/soldasim/DistributionMetrics.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Metrics for measuring how well an approximate probability distribution matches a reference distribution,
either from samples or from (possibly unnormalized) log-densities.
Originally part of [BOSIP.jl](https://github.com/soldasim/BOSIP.jl); usable on its own with any method.

| Metric | Type | Input | Notes |
|---|---|---|---|
| `TVMetric` | `PDFMetric` | log-densities | Total Variation distance approximated on a precomputed weighted grid. Densities may be unnormalized. |
| `MMDMetric` | `SampleMetric` | samples | Maximum Mean Discrepancy with a given kernel. |
| `OptMMDMetric` | `SampleMetric` | samples | MMD with kernel lengthscales optimized at each evaluation. |

All metrics are evaluated by `calculate_metric` (lower is better). Samples are matrices with one sample per column.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/soldasim/DistributionMetrics.jl")
```

## Usage

```julia
using DistributionMetrics, KernelFunctions

# Samples: Maximum Mean Discrepancy
true_samples = randn(2, 500)
approx_samples = randn(2, 500) .+ 0.5
mmd = MMDMetric(; kernel=with_lengthscale(SqExponentialKernel(), [1., 1.]))
calculate_metric(mmd, true_samples, approx_samples)

# Log-densities: Total Variation distance on a grid
grid = rand(2, 10_000) .* 12 .- 6                    # uniform on [-6, 6]^2
log_ws = fill(log(12^2), size(grid, 2))              # -log q(x) for the uniform grid density q
true_logpdf(x) = -sum(abs2, x) / 2
approx_logpdf(x) = -sum(abs2, x .- 0.5) / 2
tv = TVMetric(; grid, log_ws, true_logpost=true_logpdf)
calculate_metric(tv, true_logpdf, approx_logpdf)
```

The `true_logpost` keyword of `TVMetric` caches the true log-density values on the grid,
which makes repeated evaluations cheap.
