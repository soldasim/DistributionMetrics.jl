module DistributionMetrics

using KernelFunctions
using Optimization
using Statistics

export DistributionMetric, SampleMetric, PDFMetric
export calculate_metric
export TVMetric
export MMDMetric, OptMMDMetric

include("metric.jl")
include("utils.jl")
include("tv.jl")
include("mmd.jl")
include("opt_mmd.jl")

end
