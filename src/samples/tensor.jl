
#############################################################################################################################################################
### TensorSample is used by samplers that compute (up to) the tensor of the log-target (for ex SmMALA)

immutable TensorSample{N<:Number,T<:AbstractFloat,V<:AbstractVector,A<:AbstractArray,D<:AbstractArray} <: AbstractSample{SecondOrder}
    values::A
    loglikelihood::V
    logprior::V
    gradloglikelihood::A
    gradlogprior::A
    tensorloglikelihood::D
    tensorlogprior::D
    TensorSample(v::AbstractArray{N},ll::AbstractVector{T},lp::AbstractVector{T},gl::AbstractArray{T},gp::AbstractArray{T},
                 tl::AbstractArray{T},tp::AbstractArray{T}) = new(v,ll,lp,gl,gp,tl,tp)
end

### Factory function
@inline function _samples{N<:Number,T<:AbstractFloat}(::Type{Val{:tensor}},nparas::Integer,nsamples::Integer,::Type{N},::Type{T})
    v = zeros(N,_valuestuple(nparas,nsamples))
    u = _tensortuple(Val{ndims(v)},nparas,nparas,nsamples)
    TensorSample{N,T,Vector,Array,Array}(v,Vector{T}(nsamples),Vector{T}(nsamples),similar(v,T),similar(v,T),similar(v,T,u),similar(v,T,u))
end

################################################################################################################################################################
### TangentTensorSample is used by samplers that approximate the calculation of the tensor of the log-target (for ex TrSmMALARandom)

immutable TangentTensorSample{N<:Number,T<:AbstractFloat,V<:AbstractVector,A<:AbstractArray,D<:AbstractArray} <: AbstractSample{SecondOrder}
    values::A
    loglikelihood::V
    logprior::V
    gradloglikelihood::A
    gradlogprior::A
    tensorloglikelihood::D
    tensorlogprior::D
    tangentvectors::D
    TangentTensorSample(v::AbstractArray{N},ll::AbstractVector{T},lp::AbstractVector{T},gl::AbstractArray{T},gp::AbstractArray{T},
                        tl::AbstractArray{T},tp::AbstractArray{T},tv::AbstractArray{T}) = new(v,ll,lp,gl,gp,tl,tp,tv)
end

@inline function _samples{N<:Number,T<:AbstractFloat}(::Type{Val{:tangent}},nparas::Integer,nsamples::Integer,ntangents::Integer,::Type{N},::Type{T})
    v = zeros(N,_valuestuple(nparas,nsamples))
    u1 = _tensortuple(Val{ndims(v)},nparas,nparas,nsamples)
    u2 = _tensortuple(Val{ndims(v)},nparas,ntangents,nsamples)
    TangentTensorSample{N,T,Vector,Array,Array}(v,Vector{T}(nsamples),Vector{T}(nsamples),similar(v,T),similar(v,T),similar(v,T,u1),similar(v,T,u1),similar(v,T,u2))
end

### Additional functionality for tangent tensor samples
numtangents(s::TangentTensorSample) = size(s.tangentvectors,2)
