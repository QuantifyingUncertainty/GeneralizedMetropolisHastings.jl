### GradientSample is used by samplers that compute (up to) the gradient of the log-target (for ex HMC)
immutable GradientSample{N<:Number,T<:AbstractFloat,V<:AbstractVector,A<:AbstractArray} <: AbstractSample{FirstOrder}
    values::A
    loglikelihood::V
    logprior::V
    gradloglikelihood::A
    gradlogprior::A
    GradientSample(v::AbstractArray{N},ll::AbstractVector{T},lp::AbstractVector{T},gl::AbstractArray{T},gp::AbstractArray{T}) = new(v,ll,lp,gl,gp)
end

function GradientSample{N<:Number,T<:AbstractFloat}(nparas::Integer,nsamples::Integer,::Type{N},::Type{T})
    v = zeros(N,_valuestuple(nparas,nsamples))
    GradientSample{N,T,Vector,Array}(v,Vector{T}(nsamples),Vector{T}(nsamples),similar(v,T),similar(v,T))
end

### Factory function
_samples{N<:Number,T<:AbstractFloat}(::Type{Val{:gradient}},nparas::Integer,nsamples::Integer,::Type{N},::Type{T}) = GradientSample(nparas,nsamples,N,T)

@inline sampletypename(s::GradientSample) = :gradient







