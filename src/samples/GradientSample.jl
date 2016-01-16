### GradientSample is used by samplers that compute (up to) the gradient of the log-target (for ex HMC)
type GradientSample{T<:Number,V<:AbstractVector,M<:AbstractMatrix} <: Sample
  values::M
  loglikelihood::V
  logprior::V
  gradloglikelihood::M
  gradlogprior::M
  GradientSample(v::AbstractArray{T,2},ll::AbstractArray{T,1},lp::AbstractArray{T,1},gl::AbstractArray{T,2},gp::AbstractArray{T,2}) = new(v,ll,lp,gl,gp)
end

### External constructor
function GradientSample(v::AbstractMatrix)
  e = eltype(v)
  t1 = typeof(similar(v,1))
  t2 = typeof(v)
  s2 = size(v,2)
  GradientSample{e,t1,t2}(v,similar(v,s2),similar(v,s2),similar(v),similar(v))
end

### Factory function
@inline _samples{T<:Number}(::Type{Val{:gradient}},nparas::Integer,nsamples::Integer,::Type{T}) = GradientSample(Array{T}(nparas,nsamples))




