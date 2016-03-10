immutable NoiseModelGaussian{T<:AbstractFloat,V<:AbstractVector} <: AbstractNoiseModel
    variance::V
    distributions::Vector{Distributions.Normal}
    NoiseModelGaussian(v::AbstractVector{T},d::Vector{Distributions.Normal}) = new(v,d)
end

@inline _noise(::Type{Val{:gaussian}},v::AbstractVector) = NoiseModelGaussian{eltype(v),typeof(v)}(v,map(Distributions.Normal,zeros(v),sqrt(v)))

function loglikelihood{T<:AbstractFloat}(n::NoiseModelGaussian,measurements::AbstractArray{T},modeldata::AbstractArray{T})
    r = zero(T)
    for j=1:length(n.distributions)
        @simd for i=1:size(measurements,1)
            @inbounds r+= convert(T,logpdf(n.distributions[j],modeldata[i,j]-measurements[i,j]))
        end
    end
    r
end

function applynoise!{T<:AbstractFloat}(n::NoiseModelGaussian,d::AbstractArray{T})
    for j=1:length(n.distributions)
        @simd for i=1:size(d,1)
            @inbounds d[i,j] += convert(T,rand(n.distributions[j]))
        end
    end
    d
end

noisemodelname(::NoiseModelGaussian) = "Gaussian Noise Model"
