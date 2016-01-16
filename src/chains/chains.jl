### MCChain stores the output of a Monte Carlo iteration
type MCChainOrdinary <: MCChain
  values::Matrix{Float64}
  logprior::Vector{Float64}
  loglikelihood::Vector{Float64}
  acceptance::Int
  nsample::Int
  runtime::Float64
end

type MCChainWithGradient <: MCChain
  values::Matrix{Float64}
  logprior::Vector{Float64}
  loglikelihood::Vector{Float64}
  gradloglikelihood::Matrix{Float64}
  acceptance::Int
  nsample::Int
  runtime::Float64
end

function create_chain(nparas::Int,nsamples::Int;storegradient::Bool=false,runtime::Float64=0.0)
  c::MCChain
  if !storegradient
    c = MCChainOrdinary(zeros(nparas,nsamples),zeros(nsamples),zeros(nsamples),0,0,runtime)
  else
    c = MCChainWithGradient(zeros(nparas,nsamples),zeros(nsamples),zeros(nsamples),zeros(nparas,nsamples),0,0,runtime)
  end
  return c;
end

function store_results!(c::MCChainOrdinary,h::MCHeap,indicators::Vector{Int})
  for i = 2:length(indicators)
    c.nsample += 1
    store_common!(c,h.samples[indicators[i]])
  end
  add_acceptance!(c,indicators)
  return c;
end

function store_results!(c::MCChainWithGradient,h::MCHeap,indicators::Vector{Int})
  for i = 2:length(indicators)
    c.nsample += 1
    store_common!(c,h.samples[indicators[i]])
    store_gradient!(c,h.samples[indicators[i]])
  end
  add_acceptance!(c,indicators)
  return c;
end

function store_common!(c::MCChain,s::MCSample)
  c.values[:,c.nsample] = s.values
  c.logprior[c.nsample] = s.logprior
  c.loglikelihood[c.nsample] = s.loglikelihood
end

store_gradient!(c::MCChain,s::MCSample) = (c.gradloglikelihood[:,c.nsample] = s.gradloglikelihood)

add_acceptance!(c::MCChain,indicators::Vector{Int}) = @simd for i=2:length(indicators) @inbounds indicators[i]!=indicators[i-1]?c.acceptance+=1:nothing end

import Base.==
==(c1::MCChainOrdinary,c2::MCChainOrdinary) = (isequal(c1.values,c2.values) && isequal(c1.logprior,c2.logprior) && isequal(c1.loglikelihood,c2.loglikelihood))
==(c1::MCChainWithGradient,c2::MCChainWithGradient) = (isequal(c1.values,c2.values) && isequal(c1.logprior,c2.logprior) && isequal(c1.loglikelihood,c2.loglikelihood) && isequal(c1.gradloglikelihood,c2.gradloglikelihood))

function Base.show(io::IO,c::MCChain)
  println("MCChain with nparas: ",size(c.values,1)," and nsamples: ",size(c.values,2))
  println("Samples performed: ",c.nsample,", total accepted: ",c.acceptance,", acceptance rate: ",c.acceptance/c.nsample)
  println("Runtime: ",c.runtime)
end
