function SampleIndicator(Chain::MarkovChain)
# Input: Chain
# Output: Vector of samples of indicator variable

# Initialise the acceptance probabilities
A   = zeros(Chain.NumOfProposals + 1)
Acc = zeros(Chain.NumOfProposals + 1)

Posteriors = zeros(Chain.NumOfProposals + 1)
for i = 1 : Chain.NumOfProposals + 1
    Posteriors[i] = Chain.Geometry[i].LL + Chain.Geometry[i].LogPrior
end
Acc = vec(exp(Posteriors - maximum(Posteriors)) / sum(exp(Posteriors - maximum(Posteriors))))

# Sample indicator variable N times
IndicatorSamples = rand(Categorical(Acc), Chain.NumOfProposals) # Note Acc must sum to 1

# Update the sample indicator
Chain.SampleIndicator = IndicatorSamples[end]

return IndicatorSamples

end
