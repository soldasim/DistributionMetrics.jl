"""
    TVMetric(; grid, ws | log_ws, true_logpost | true_logvals)

Approximates the Total Variation (TV) distance between the true and the approximate distribution
by a weighted sum over a precomputed grid of points.

## Keywords
- `grid::Matrix{Float64}`: The points used to approximate the TV integral (one point per column).
- `log_ws::Vector{Float64}`: The log-weights of the grid points. Should be `-log q(x)`,
    where `q` is the density used to sample the grid
    (`-log(area)` is appropriate for an evenly distributed grid).
    The non-logarithmic weights `ws` can be given instead.
- `true_logpost::Function`: The (possibly unnormalized) true log-density.
    If given, its values on the grid are computed once and cached.
- `true_logvals::Vector{Float64}`: The (possibly unnormalized) true log-density values on the grid.
    Used directly if the values are precomputed or can be computed more efficiently in a batch.
"""
struct TVMetric <: PDFMetric
    grid::Matrix{Float64}
    log_ws::Vector{Float64}
    true_logvals::Union{Nothing, Vector{Float64}}
end
function TVMetric(;
    grid,
    ws = nothing,
    log_ws = nothing,
    true_logpost = nothing,
    true_logvals = nothing,
)
    @assert xor(isnothing(ws), isnothing(log_ws))
    isnothing(log_ws) && (log_ws = log.(ws))

    if isnothing(true_logvals) && !isnothing(true_logpost)
        true_logvals = true_logpost.(eachcol(grid))
    end

    return TVMetric(grid, log_ws, true_logvals)
end

function calculate_metric(tv::TVMetric, true_logpost::Function, approx_logpost::Function; info::Bool=false)
    if isnothing(tv.true_logvals)
        true_logvals = true_logpost.(eachcol(tv.grid))
    else
        true_logvals = tv.true_logvals
    end
    approx_logvals = approx_logpost.(eachcol(tv.grid))

    return calc_tv(tv.log_ws, approx_logvals, true_logvals)
end

function calc_tv(log_ws, approx_logvals, true_logvals)
    ev_approx = logmeanexp(log_ws .+ approx_logvals)
    ev_true = logmeanexp(log_ws .+ true_logvals)

    # Non-normalizable distributions (e.g. due to numerical issues) get zero density everywhere.
    approx_logvals_norm = isinf(ev_approx) ? fill(-Inf, length(approx_logvals)) : approx_logvals .- ev_approx
    true_logvals_norm = isinf(ev_true) ? fill(-Inf, length(true_logvals)) : true_logvals .- ev_true

    diffs = @. abs( exp(approx_logvals_norm) - exp(true_logvals_norm) )
    tv = (1/2) * exp(logmeanexp(log_ws .+ log.(diffs)))

    return tv
end
