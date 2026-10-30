"""
Subtypes of `DistributionMetric` evaluate the quality of an approximate distribution
with respect to a reference (true) distribution.

The metrics are grouped into two categories; [`SampleMetric`](@ref) and [`PDFMetric`](@ref).
"""
abstract type DistributionMetric end

"""
A [`DistributionMetric`](@ref) computed from samples of the true and the approximate distribution.

Each subtype *should* implement:
- `calculate_metric(::SampleMetric, true_samples::AbstractMatrix{<:Real}, approx_samples::AbstractMatrix{<:Real}; info=false) -> ::Real`
"""
abstract type SampleMetric <: DistributionMetric end

"""
A [`DistributionMetric`](@ref) computed from the (possibly unnormalized) log-densities
of the true and the approximate distribution.

Each subtype *should* implement:
- `calculate_metric(::PDFMetric, true_logpdf::Function, approx_logpdf::Function; info=false) -> ::Real`
"""
abstract type PDFMetric <: DistributionMetric end

"""
    score = calculate_metric(::SampleMetric, true_samples, approx_samples; info=false)
    score = calculate_metric(::PDFMetric, true_logpdf, approx_logpdf; info=false)

Evaluate the accuracy of the approximate distribution. Lower scores are better.

Samples are matrices with one sample per column.
Log-densities are functions mapping a vector to a real number.

## Keywords
- `info::Bool=false`: Print progress messages.
"""
function calculate_metric end
