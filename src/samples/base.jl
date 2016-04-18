### Generic BaseSample type used by most samplers
immutable BaseSample{N<:Number,T<:AbstractFloat,V<:AbstractVector,A<:AbstractArray} <: AbstractSample{ZeroOrder}
    values::A
    loglikelihood::V
    logprior::V
    BaseSample(v::AbstractArray{N},ll::AbstractVector{T},lp::AbstractVector{T}) = new(v,ll,lp)
end

### Constructors
function BaseSample{N<:Number,T<:AbstractFloat}(nparas::Integer,nsamples::Integer,::Type{N},::Type{T})
    v = zeros(N,_valuestuple(nparas,nsamples))
    ll = Vector{T}(nsamples)
    lp = Vector{T}(nsamples)
    BaseSample{N,T,typeof(ll),typeof(v)}(v,ll,lp)
end

### Factory function
_samples{N<:Number,T<:AbstractFloat}(::Type{Val{:base}},nparas::Integer,nsamples::Integer,::Type{N},::Type{T}) = BaseSample(nparas,nsamples,N,T)

@inline sampletypename(s::BaseSample) = :base




