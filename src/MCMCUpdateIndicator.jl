function SampleIndicator(Chain::MarkovChain)
# Input: Chain
# Output: Vector of samples of indicator variable

# Initialise the acceptance probabilities
A   = zeros(1,Chain.NumOfProposals + 1)

for i = 1 : Chain.NumOfProposals + 1
    A[i] = Chain.Geometry[i].LL + Chain.Geometry[i].LogPrior
end
Acc = (exp(A - maximum(A)) / sum(exp(A - maximum(A))))

# Sample indicator variable N times
result = sample_indicator_matrix(repmat(Acc,Chain.NumOfProposals+1,1),Chain.NumOfProposals,Chain.SampleIndicator)

# Update the sample indicator
Chain.SampleIndicator = result[end]

result

end

function sample_indicator_matrix(A::Array{Float64,2},numproposals::Int64,startindex::Int64)

  #initialise variables
  result = zeros(Int64,numproposals + 1)
  result[1] = startindex

  #sample the indicator transition matrix numprop times
  for i = 1:numproposals
    result[i+1] = rand(Categorical(vec(A[result[i],:])))
  end

  #return the indicator samples
  result

end
