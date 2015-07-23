abstract AdaptiveMetropolisHastings <: MCSampler #derive from this for non-symmetric proposal densities
abstract AdaptiveMetropolis <: MetropolisHastings #derive from this for symmetric proposal densities

### Type holding the parameters for a random-walk Metrolopolis sampler
immutable AdaptiveMHNormal <: AdaptiveMetropolis #mormal distribution is symmetric Metropolis sampler
  initialscaling::Float64
end

### Construct MHNormal based on number of parameters and initial scaling
AdaptiveMHNormal(nparas::Int,s::Float64) = MHNormal(eye(nparas),s)
nparas(s::MHNormal) = size(s.covariance,1)

### Type holding the state of the Markov Chain for a Generalized M-H sampler
type MHHeap <: MCHeap{BaseSample}
  samples::Vector{BaseSample}
  sampledensity::ProposalDensity
  fromvalues::Vector{Float64}
  runningmean::Vector{Float64}
  runningcovariance::Matrix{Float64}
  scaling::Float64
end


        # Update the proposal mean and covariance
        #println(CurIter)
        if CurIter == 1
            Chain.ProposalDistribution.RunningCov  = (1/CurIter)*(Chain.Geometry[Chain.SampleIndicator].Parameters-Chain.ProposalDistribution.RunningMean)*(Chain.Geometry[Chain.SampleIndicator].Parameters-Chain.ProposalDistribution.RunningMean)';
            Chain.ProposalDistribution.RunningMean = ((CurIter-1)*Chain.ProposalDistribution.RunningMean + Chain.Geometry[Chain.SampleIndicator].Parameters)/CurIter;
            Chain.ProposalDistribution.Density     = MvNormal(zeros(Model.NumOfParas),
                                                              eye(Model.NumOfParas)*(Chain.ProposalDistribution.StepSize^2))
        else
            Chain.ProposalDistribution.RunningCov  = ((CurIter-2)/(CurIter-1))*Chain.ProposalDistribution.RunningCov + (1/CurIter)*(Chain.Geometry[Chain.SampleIndicator].Parameters-Chain.ProposalDistribution.RunningMean)*(Chain.Geometry[Chain.SampleIndicator].Parameters-Chain.ProposalDistribution.RunningMean)';
            Chain.ProposalDistribution.RunningMean = ((CurIter-1)*Chain.ProposalDistribution.RunningMean + Chain.Geometry[Chain.SampleIndicator].Parameters)/CurIter;
            Chain.ProposalDistribution.Density     = MvNormal(zeros(Model.NumOfParas),
                                                              Chain.ProposalDistribution.RunningCov + eye(Model.NumOfParas)*(Chain.ProposalDistribution.StepSize^2))
        end


        # Propose the new parameters values
        for PropNum = 1:Chain.NumOfProposals+1
            # Update for all except the indicated sample
            if PropNum != Chain.SampleIndicator

                # Propose new point
                #println("Proposing new point...")
                Chain.Geometry[PropNum].Parameters    = Chain.Geometry[Chain.SampleIndicator].Parameters + rand(Chain.ProposalDistribution.Density) # using same proposal
                end

