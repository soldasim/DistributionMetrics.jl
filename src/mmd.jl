"""
    MMDMetric(; kernel)

Approximates the Maximum Mean Discrepancy (MMD) between samples of the true and the approximate distribution.

## Keywords
- `kernel::Kernel`: The kernel used to calculate the MMD.
    The kernel lengthscales have to be chosen appropriately (see [`OptMMDMetric`](@ref) for automatic tuning).
"""
@kwdef struct MMDMetric <: SampleMetric
    kernel::Kernel
end

function calculate_metric(mmd::MMDMetric, true_samples::AbstractMatrix{<:Real}, approx_samples::AbstractMatrix{<:Real}; info::Bool=false)
    return calc_mmd(mmd.kernel, eachcol(true_samples), eachcol(approx_samples))
end

function calc_mmd(kernel::Kernel, X::AbstractVector, Y::AbstractVector)
    val_X = mean(kernelmatrix(kernel, X))
    val_Y = mean(kernelmatrix(kernel, Y))
    val_XY = mean(kernelmatrix(kernel, X, Y))
    return val_X + val_Y - 2*val_XY
end
