### TangentTensorSample is used by samplers that approximate the calculation of the tensor of the log-target (for ex TrSmMALARandom)
type TangentTensorSample{T<:Number,V<:AbstractVector,M<:AbstractMatrix,D<:AbstractArray{3}} <: Sample
  values::M
  loglikelihood::V
  logprior::V
  gradloglikelihood::M
  gradlogprior::M
  tensorloglikelihood::D
  tensorlogprior::D
  tangentvectors::D
  TangentTensorSample(v::AbstractArray{T,2},ll::AbstractArray{T,1},lp::AbstractArray{T,1},gl::AbstractArray{T,2},gp::AbstractArray{T,2},tl::AbstractArray{T,3},tp::AbstractArray{T,3},tv::AbstractArray{T,3}) = new(v,ll,lp,gl,gp,tl,tp,tv)
end

function TangentTensorSample(v::AbstractArray,ntangents::Integer)
  e = eltype(v)
  t1 = typeof(similar(v,1))
  t2 = typeof(v)
  t3 = typeof(similar(v,1,1,1))
  s1 = size(v,1)
  s2 = size(v,2)
  TangentTensorSample{e,t2,t1,t3}(v,similar(v,s2),similar(v,s2),similar(v),similar(v),similar(v,s1,s1,s2),similar(v,s1,s1,s2),similar(v,s1,ntangents,s2))
end

### Additional functionality for tangent tensor samples
numtangents(s::TangentTensorSample) = size(s.tangentvectors,2)
@inline gettangent(s::TangentTensorSample,para::Integer,tangent::Integer,sample::Integer) = s.tangentvectors[para,tangent,sample]
@inline settangent!(s::TangentTensorSample,para::Integer,tangent::Integer,sample::Integer,val::Number) = (s.tangentvectors[para,tangent,sample] = val)



