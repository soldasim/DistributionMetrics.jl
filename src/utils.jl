function logsumexp(x::AbstractVector{<:Real})
    x_max = maximum(x)
    isinf(x_max) && return x_max
    return x_max + log(sum(exp.(x .- x_max)))
end

logmeanexp(x::AbstractVector{<:Real}) = logsumexp(x) - log(length(x))
