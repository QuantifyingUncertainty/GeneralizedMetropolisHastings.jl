### Generic BaseSample type used by most samplers

type BaseSample <: MCSample{NullOrder}
  values::Vector{Float64}
  loglikelihood::Float64
  logprior::Float64
end

BaseSample(s::Vector{Float64}) = BaseSample(s,NaN,NaN)
BaseSample(n::Int) = BaseSample(zeros(n))

==(s1::BaseSample,s2::BaseSample) = (isequal(s1.values,s2.values) && isequal(s1.loglikelihood,s2.loglikelihood) && isequal(s1.logprior,s2.logprior))

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

==(s1::GradientSample,s2::GradientSample) = (isequal(s1.values,s2.values) && isequal(s1.loglikelihood,s2.loglikelihood) && isequal(s1.logprior,s2.logprior) &&
                                               isequal(s1.gradloglikelihood,s2.gradloglikelihood) && isequal(s1.gradlogprior,s2.gradlogprior))

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

==(s1::TensorSample,s2::TensorSample) = (isequal(s1.values,s2.values) && isequal(s1.loglikelihood,s2.loglikelihood) && isequal(s1.logprior,s2.logprior) &&
                                               isequal(s1.gradloglikelihood,s2.gradloglikelihood) && isequal(s1.gradlogprior,s2.gradlogprior) &&
                                           isequal(s1.tensorloglikelihood,s2.tensorloglikelihood) && isequal(s1.tensorlogprior,s2.tensorlogprior))

function Base.show(io::IO,s::BaseSample)
  println(io,"BaseSample with")
  println(io,"values: ",s.values)
  println(io,"loglikelihood: ",s.loglikelihood)
  println(io,"logprior: ",s.logprior)
  println(io)
  nothing
end

function Base.show(io::IO,s::GradientSample)
  println(io,"GradientSample with")
  println(io,"values: ",s.values)
  println(io,"loglikelihood: ",s.loglikelihood)
  println(io,"logprior: ",s.logprior)
  println(io,"gradloglikelihood: ",s.gradloglikelihood)
  println(io,"gradlogprior: ",s.gradlogprior)
  println(io)
  nothing
end

function Base.show(io::IO,s::TensorSample)
  println(io,"TensorSample with")
  println(io,"values: ",s.values)
  println(io,"loglikelihood: ",s.loglikelihood)
  println(io,"logprior: ",s.logprior)
  println(io,"gradloglikelihood: ",s.gradloglikelihood)
  println(io,"gradlogprior: ",s.gradlogprior)
  println(io,"tensorloglikelihood: ")
  show(io,s.tensorloglikelihood)
  println(io)
  println(io,"tensorlogprior: ")
  show(io,s.tensorlogprior)
  println(io)
  println(io)
  nothing
end

function Base.show{S<:MCSample}(io::IO,v::Array{S})
  println(typeof(v)," with ",length(v)," samples")
  for i=1:length(v)
    show(io,v[i])
  end
end



