immutable IndicatorStationary{T<:AbstractFloat} <: AbstractIndicatorMatrix
    stationary::Vector
    samples::Vector
    IndicatorStationary(s::Vector{T},i::Vector{Int}) = new(s,i)
end

_indicator{T<:AbstractFloat}(::Type{Val{:stationary}},numproposals::Int,numsamples::Int,::Type{T}) = IndicatorStationary{T}(zeros(T,numproposals+1),zeros(Int,numsamples+1))

numproposals(indicator::IndicatorStationary) = length(indicator.stationary) - 1

function transitionprobability!{T<:AbstractFloat}(indicator::IndicatorStationary{T},indicatoracceptance::AbstractVector{T},auxiliaryacceptances::AbstractVector)
    icounter = 0
    for k=1:length(auxiliaryacceptances)
        acceptanceratio = auxiliaryacceptances[k]
        @simd for i=1:length(acceptanceratio)
            @inbounds indicator.stationary[icounter+=1] = acceptanceratio[i] + indicatoracceptance[1]
        end
    end
    indicator.stationary[end] = zero(T)
    _calculatetransition(indicator)
end

function transitionprobability!{T<:AbstractFloat}(indicator::IndicatorStationary{T},acceptanceratio::AbstractVector{T})
    @simd for i=1:length(acceptanceratio)
        @inbounds indicator.stationary[i] = acceptanceratio[i]
    end
    indicator.stationary[end] = zero(T)
    _calculatetransition(indicator)
end

function _calculatetransition(indicator::IndicatorStationary)
    maxacc = maximum(indicator.stationary)
    @simd for i=1:length(indicator.stationary)
        @inbounds indicator.stationary[i] = exp(indicator.stationary[i] - maxacc)
    end
    sumacc = sum(indicator.stationary)
    @simd for i=1:length(indicator.stationary)
        @inbounds indicator.stationary[i] /= sumacc
    end
    indicator.stationary
end

function sampleindicator!(indicator::IndicatorStationary)
    c = Distributions.Categorical(indicator.stationary)
    indicator.samples[1] = length(indicator.stationary)
    @simd for i=2:length(indicator.samples)
        @inbounds indicator.samples[i] = rand(c)
    end
    indicator.samples
end

function show(io::IO,indicator::IndicatorStationary)
    np = numproposals(indicator)
    ns = numsamples(indicator)
    println(io,"IndicatorStationary with $(np) proposal",np>1?"s":""," and $(ns) sample",ns>1?"s":"")
    println(io," Fields: :stationary, :samples")
end
