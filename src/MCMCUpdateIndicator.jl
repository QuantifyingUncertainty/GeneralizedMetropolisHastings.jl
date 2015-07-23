function sample_indicator(chain::MarkovChain; matrixcreator::IndicatorMatrixCreateFunction = IndicatorMatrixStandard())
  # Input: chain
  # Optional: matrixcreator
  # Output: Vector of samples of indicator variable

  # Initialise
  n = chain.NumOfProposals
  v = zeros(n+1)

  for i = 1:n+1
    v[i] = chain.Geometry[i].LL + chain.Geometry[i].LogPrior
  end

  v = (exp(v - maximum(v)) / sum(exp(v - maximum(v))))
  a = create_indicator_matrix(v,n,matrixcreator)

  # Sample indicator variable N times
  result = sample_indicator_matrix(a,n,chain.SampleIndicator)

  # Update the sample indicator
  chain.SampleIndicator = result[end]

  # Return the n results
  result

end

create_indicator_matrix(::IndicatorMatrixStationary,v::Vector{Float64},n::Int64) = repmat(v',n+1)
create_indicator_matrix(::IndicatorMatrixOptimal,v::Vector{Float64},n::Int64) = eye(n+1)

function sample_indicator_matrix(a::Array{Float64,2},numsamples::Int64,startindex::Int64)

  #initialise variables
  r = zeros(Int64,numsamples)
  i = startindex

  #sample the indicator transition matrix numprop times
  for j = 1:numsamples
    i = r[j] = rand(Categorical(vec(a[i,:])))
  end

  #return the indicator samples
  result

end
