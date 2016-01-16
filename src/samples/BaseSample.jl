### Generic BaseSample type used by most samplers
type BaseSample{T<:Number,V<:AbstractVector,M<:AbstractMatrix} <: Sample{}
  values::M
  loglikelihood::V
  logprior::V
  BaseSample(v::AbstractArray{T,2},ll::AbstractArray{T,1},lp::AbstractArray{T,1}) = new(v,ll,lp)
end

### External constructor
function BaseSample(v::AbstractMatrix)
  e = eltype(v)
  t1 = typeof(similar(v,1))
  t2 = typeof(v)
  s2 = size(v,2)
  BaseSample{e,t1,t2}(v,similar(v,s2),similar(v,s2))
end

### Factory function
@inline _samples{T<:Number}(::Type{Val{:base}},nparas::Integer,nsamples::Integer,::Type{T}) = BaseSample(Array{T}(nparas,nsamples))
