module GeneralizedMetropolisHastings

using Distributions
using Sundials

#imports

export
	###types
	ODEModel,TargetOnlyModel,GaussianBivariate,MarkovChain,MarkovChainGeometry,MarkovChainProposal,MCMCSimulation,
	ProposalDistributionMH, ProposalDistributionSmMALA, ProposalDistributionSmMALARandom, ProposalDistributionAdaptiveMH,
	###functions
	MCMCRun

include("MCMCTypes.jl")
include("MCMCUpdateParameters.jl")
include("MCMCUpdateIndicator.jl")
include("MCMCRun.jl")

end # module