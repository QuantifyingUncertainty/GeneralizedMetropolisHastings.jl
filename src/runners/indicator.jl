function sample_indicator(h::MCHeap,indicator::Int,indicatefunction::IndicatorMatrixFunction)

  nsamples = numel(h)
  v = zeros(nsamples)

  for i = 1:nsamples
    v[i] = h.samples[i].loglikelihood + h.samples[i].logprior
  end

  a = (exp(v - maximum(v)) / sum(exp(v - maximum(v))))
  A = create_indicator_matrix(indicatefunction,a)

  # Sample indicator variable N times and return the result
  sample_indicator_matrix(A,nsamples,indicator)
end

create_indicator_matrix(::IndicatorMatrixStationary,v::Vector{Float64}) = repmat(v',length(v))
create_indicator_matrix(::IndicatorMatrixOptimal,v::Vector{Float64}) = eye(length(v))

function sample_indicator_matrix(a::Matrix{Float64},nsamples::Int,indicator::Int)

  #initialise variables
  indicators = zeros(Int,nsamples)

  #sample the indicator transition matrix nsamples times
  indicators[1] = indicator
  for j = 2:nsamples
    indicators[j] = rand(Distributions.Categorical(vec(a[indicators[j-1],:])))
  end

  #return the indicator samples
  indicators
end
