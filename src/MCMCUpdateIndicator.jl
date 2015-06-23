function SampleIndicator(Chain::MarkovChain)
# Input: Chain
# Output: Vector of samples of indicator variable

# Initialise the acceptance probabilities
A   = zeros(Chain.NumOfProposals + 1)
Acc = zeros(Chain.NumOfProposals + 1)

for i = 1 : Chain.NumOfProposals + 1
    A[i] = Chain.Geometry[i].LL + Chain.Geometry[i].LogPrior
end
Acc = vec(exp(A - maximum(A)) / sum(exp(A - maximum(A))))

# Sample indicator variable N times
IndicatorSamples = rand(Categorical(Acc), Chain.NumOfProposals) # Note Acc must sum to 1

# Update the sample indicator
Chain.SampleIndicator = IndicatorSamples[end]

return IndicatorSamples

end
