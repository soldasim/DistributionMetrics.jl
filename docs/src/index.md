```@meta
CurrentModule = DistributionMetrics
```

# DistributionMetrics

Metrics for measuring how well an approximate probability distribution matches a reference distribution,
either from samples or from (possibly unnormalized) log-densities.

All metrics are evaluated by [`calculate_metric`](@ref) (lower is better).
Samples are matrices with one sample per column.

## Metric types

```@docs
DistributionMetric
SampleMetric
PDFMetric
calculate_metric
```

## Total Variation

```@docs
TVMetric
```

## Maximum Mean Discrepancy

```@docs
MMDMetric
OptMMDMetric
```
