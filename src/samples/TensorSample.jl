### TensorSample is used by samplers that compute (up to) the tensor of the log-target (for ex SmMALA)
type TensorSample{T<:Number,V<:AbstractVector,M<:AbstractMatrix,D<:AbstractArray{3}} <: Sample
  values::M
  loglikelihood::V
  logprior::V
  gradloglikelihood::M
  gradlogprior::M
  tensorloglikelihood::D
  tensorlogprior::D
  TensorSample(v::AbstractArray{T,2},ll::AbstractArray{T,1},lp::AbstractArray{T,1},gl::AbstractArray{T,2},gp::AbstractArray{T,2},tl::AbstractArray{T,3},tp::AbstractArray{T,3}) = new(v,ll,lp,gl,gp,tl,tp)
end

function TensorSample(v::AbstractArray)
  e = eltype(v)
  t1 = typeof(similar(v,1))
  t2 = typeof(v)
  t3 = typeof(similar(v,1,1,1))
  s1 = size(v,1)
  s2 = size(v,2)
  TensorSample{e,t2,t1,t3}(v,similar(v,s2),similar(v,s2),similar(v),similar(v),similar(v,s1,s1,s2),similar(v,s1,s1,s2))
end

@inline _samples{T<:Number}(::Type{Val{:tensor}},nparas::Integer,nsamples::Integer,::Type{T}) = TensorSample(Array{T}(nparas,nsamples))
