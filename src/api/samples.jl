### Generic BaseSample type used by most samplers

type BaseSample <: MCSample{NullOrder}
  values::Vector{Float64}
  loglikelihood::Float64
  logprior::Float64
end

BaseSample(s::Vector{Float64}) = BaseSample(s,NaN,NaN)
BaseSample(n::Int) = BaseSample(zeros(n))

### GradientSample is used by samplers that compute (up to) the gradient of the log-target (for ex HMC)

type GradientSample <: MCSample{FirstOrder}
  values::Vector{Float64}
  loglikelihood::Float64
  logprior::Float64
  gradloglikelihood::Vector{Float64}
  gradlogprior::Vector{Float64}
end

GradientSample(s::Vector{Float64}) = GradientSample(s,NaN,NaN,fill(NaN,length(s)),fill(NaN,length(s)))
GradientSample(n::Int) = GradientSample(zeros(n))

### TensorSample is used by samplers that compute (up to) the tensor of the log-target (for ex SmMALA)

type TensorSample <: MCSample{SecondOrder}
  values::Vector{Float64}
  loglikelihood::Float64
  logprior::Float64
  gradloglikelihood::Vector{Float64}
  gradlogprior::Vector{Float64}
  tensorloglikelihood::Matrix{Float64}
  tensorlogprior::Matrix{Float64}
end

TensorSample(s::Vector{Float64}) = TensorSample(s,NaN,NaN,fill(NaN,length(s)),fill(NaN,length(s)),fill(NaN,length(s),length(s)),fill(NaN,length(s),length(s)))
TensorSample(n::Int) = TensorSample(zeros(n))

### ApproximateTensorSample is used by samplers that approximate the calculation of the tensor of the log-target (for ex TrSmMALARandom)

type ApproximateTensorSample <: MCSample{SecondOrder}
  values::Vector{Float64}
  loglikelihood::Float64
  logprior::Float64
  gradloglikelihood::Vector{Float64}
  gradlogprior::Vector{Float64}
  tensorloglikelihood::Matrix{Float64}
  tensorlogprior::Matrix{Float64}
  tangentvectors::Matrix{Float64}
end

ApproximateTensorSample(s::Vector{Float64},ntangent::Int) = ApproximateTensorSample(s,NaN,NaN,fill(NaN,length(s)),fill(NaN,length(s)),fill(NaN,length(s),length(s)),fill(NaN,length(s),length(s)),fill(NaN,length(s),ntangent))
ApproximateTensorSample(n::Int,ntangent::Int) = ApproximateTensorSample(zeros(n),ntangent)

numparas(s::MCSample) = length(s.values)

### Functionality from Base package
import Base.==
==(s1::BaseSample,s2::BaseSample) = (isequal(s1.values,s2.values) && isequal(s1.loglikelihood,s2.loglikelihood) && isequal(s1.logprior,s2.logprior))
==(s1::GradientSample,s2::GradientSample) = (isequal(s1.values,s2.values) && isequal(s1.loglikelihood,s2.loglikelihood) && isequal(s1.logprior,s2.logprior) &&
                                             isequal(s1.gradloglikelihood,s2.gradloglikelihood) && isequal(s1.gradlogprior,s2.gradlogprior))
==(s1::TensorSample,s2::TensorSample) = (isequal(s1.values,s2.values) && isequal(s1.loglikelihood,s2.loglikelihood) && isequal(s1.logprior,s2.logprior) &&
                                         isequal(s1.gradloglikelihood,s2.gradloglikelihood) && isequal(s1.gradlogprior,s2.gradlogprior) &&
                                         isequal(s1.tensorloglikelihood,s2.tensorloglikelihood) && isequal(s1.tensorlogprior,s2.tensorlogprior))
==(s1::ApproximateTensorSample,s2::ApproximateTensorSample) = (isequal(s1.values,s2.values) && isequal(s1.loglikelihood,s2.loglikelihood) && isequal(s1.logprior,s2.logprior) &&
                                                              isequal(s1.gradloglikelihood,s2.gradloglikelihood) && isequal(s1.gradlogprior,s2.gradlogprior) &&
                                                              isequal(s1.tensorloglikelihood,s2.tensorloglikelihood) && isequal(s1.tensorlogprior,s2.tensorlogprior) &&
                                                              isequal(s1.tangentvectors,s2.tangentvectors))

function Base.show(io::IO,s::BaseSample,p1::AbstractString = "",p2::AbstractString = "")
  println(io,p1,typeof(s)," with")
  println(io,p2," values: ",s.values)
  println(io,p2," loglikelihood: ",s.loglikelihood)
  println(io,p2," logprior: ",s.logprior)
  println(io)
  nothing
end

function Base.show(io::IO,s::GradientSample,p1::AbstractString = "",p2::AbstractString = "")
  println(io,p1,typeof(s)," with")
  println(io,p2," values: ",s.values)
  println(io,p2," loglikelihood: ",s.loglikelihood)
  println(io,p2," logprior: ",s.logprior)
  println(io,p2," gradloglikelihood: ",s.gradloglikelihood)
  println(io,p2," gradlogprior: ",s.gradlogprior)
  println(io)
  nothing
end

function Base.show(io::IO,s::TensorSample,p1::AbstractString = "",p2::AbstractString = "")
  println(io,p1,typeof(s)," with")
  println(io,p2," values: ",s.values)
  println(io,p2," loglikelihood: ",s.loglikelihood)
  println(io,p2," logprior: ",s.logprior)
  println(io,p2," gradloglikelihood: ",s.gradloglikelihood)
  println(io,p2," gradlogprior: ",s.gradlogprior)
  println(io,p2," tensorloglikelihood: ")
  show(io,s.tensorloglikelihood)
  println(io)
  println(io,p2," tensorlogprior: ")
  show(io,s.tensorlogprior)
  println(io)
  println(io)
  nothing
end

function Base.show(io::IO,s::ApproximateTensorSample,p1::AbstractString = "",p2::AbstractString = "")
  println(io,p1,typeof(s)," with")
  println(io,p2," values: ",s.values)
  println(io,p2," loglikelihood: ",s.loglikelihood)
  println(io,p2," logprior: ",s.logprior)
  println(io,p2," gradloglikelihood: ",s.gradloglikelihood)
  println(io,p2," gradlogprior: ",s.gradlogprior)
  println(io,p2," tensorloglikelihood: ")
  show(io,s.tensorloglikelihood)
  println(io)
  println(io,p2," tensorlogprior: ")
  show(io,s.tensorlogprior)
  println(io)
  println(io,p2," tangentvectors: ")
  display(s.tangentvectors)
  println(io)
  println(io)
  nothing
end

function Base.show{S<:MCSample}(io::IO,v::Array{S})
  for i=1:length(v)
    show(io,v[i],"[$i]"," ")
  end
end



