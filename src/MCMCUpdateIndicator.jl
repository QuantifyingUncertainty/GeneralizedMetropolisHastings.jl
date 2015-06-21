function SampleIndicator(Model, Chain::MarkovChain)
# Input: Model, Chain
# Output: Vector of samples of indicator variable


if Chain.Sampler == "MH" ||
   Chain.Sampler == "BivariateStructuredMarginal" ||
   Chain.Sampler == "BivariateStructuredMetropolis" ||
   Chain.Sampler == "SmMALA" ||
   Chain.Sampler == "TrSmMALA" ||
   Chain.Sampler == "TrSmMALA_Random" ||
   Chain.Sampler == "AdaptiveMH"
    ### Standard Metropolis-Hastings sampler ###

    # Initialise the acceptance probabilities
    A   = zeros(Chain.NumOfProposals + 1)
    Acc = zeros(Chain.NumOfProposals + 1)

    if Chain.NumOfProposals == 1

        # Only 1 proposal so use optimal Metropolis-Hastings ratio
        A[1] = (Chain.Geometry[1].LL + Chain.Geometry[1].LogPrior + Chain.Geometry[1].ProposalProbability[2])
        A[2] = (Chain.Geometry[2].LL + Chain.Geometry[2].LogPrior + Chain.Geometry[2].ProposalProbability[1])

        if Chain.SampleIndicator == 1
            # Proposed point is 2
            Acc[2] = min(1, exp(A[2] - A[1]))
            Acc[1] = 1 - Acc[2]
        else
            # Proposed point is 1
            Acc[1] = min(1, exp(A[1] - A[2]))
            Acc[2] = 1 - Acc[1]
        end
    else

        # Multiple proposals, so use heatbath ratios

        # NEED TO UPDATE!!!

    end

    # Sample indicator variable N times
    IndicatorSamples = rand(Categorical(Acc), Chain.NumOfProposals) # Note Acc must sum to 1


elseif Chain.Sampler == "Gibbs"
    ### Gibbs sampler ###

    # Everything accepted since it's a Gibbs sampler
    IndicatorSamples = Chain.SampleIndicator

else
    # Code for other samplers

end



# Update the sample indicator
Chain.SampleIndicator = IndicatorSamples[end]

return IndicatorSamples



end
