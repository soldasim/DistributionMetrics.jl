const AbstractBounds = Tuple{<:AbstractVector{<:Real}, <:AbstractVector{<:Real}}

"""
    OptMMDMetric(; kernel, bounds, algorithm, kwargs...)

Like [`MMDMetric`](@ref), but the kernel lengthscales are optimized
to maximize the MMD during each evaluation of the metric.

## Keywords
- `kernel::Kernel`: The kernel *without* lengthscales.
- `bounds`: Tuple `(lb, ub)` of the bounds of the sampled domain.
- `algorithm`: The optimization algorithm (an `Optimization.jl` solver) used to optimize the lengthscales.
- `kwargs...`: Passed to the optimization algorithm.
"""
struct OptMMDMetric <: SampleMetric
    kernel::Kernel
    bounds::AbstractBounds
    algorithm::Any
    kwargs::Base.Pairs{Symbol, <:Any}
end
function OptMMDMetric(;
    kernel,
    bounds,
    algorithm,
    kwargs...
)
    return OptMMDMetric(kernel, bounds, algorithm, kwargs)
end

function calculate_metric(mmd::OptMMDMetric, true_samples::AbstractMatrix{<:Real}, approx_samples::AbstractMatrix{<:Real}; info::Bool=false)
    info && @info "Optimizing kernel lengthscales for the MMD metric ..."

    function _calc_mmd(λ::AbstractVector{<:Real})
        kernel = with_lengthscale(mmd.kernel, λ)
        return calc_mmd(kernel, eachcol(true_samples), eachcol(approx_samples))
    end

    f(λ, _) = (-1) * _calc_mmd(λ)
    x0 = (mmd.bounds[2] .- mmd.bounds[1]) ./ 3

    prob = OptimizationProblem(f, x0, nothing;
        lb = zero(mmd.bounds[1]),
        ub = mmd.bounds[2] .- mmd.bounds[1],
    )
    sol = solve(prob, mmd.algorithm; mmd.kwargs...)

    return _calc_mmd(sol.u)
end
