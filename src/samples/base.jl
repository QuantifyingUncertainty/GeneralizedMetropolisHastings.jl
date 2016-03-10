### Generic BaseSample type used by most samplers
immutable BaseSample{N<:Number,T<:AbstractFloat,V<:AbstractVector,A<:AbstractArray} <: AbstractSample{ZeroOrder}
    values::A
    loglikelihood::V
    logprior::V
    BaseSample(v::AbstractArray{N},ll::AbstractVector{T},lp::AbstractVector{T}) = new(v,ll,lp)
end

### Factory function
@inline _samples{N<:Number,T<:AbstractFloat}(::Type{Val{:base}},nparas::Integer,nsamples::Integer,::Type{N},::Type{T}) =
    BaseSample{N,T,Vector,Array}(zeros(N,_valuestuple(nparas,nsamples)),Vector{T}(nsamples),Vector{T}(nsamples))


