function MCMCRun( MySimulation::MCMCSimulation )

#####################
# Initialise chains #
#####################

Sampler         = MySimulation.Sampler
NumOfProposals  = MySimulation.NumOfProposals

# Indicate which sample is the current, and which is the proposal
SampleIndex = 1    # Initialise at 1

AttemptedProposal = 0
AcceptedProposal  = 0
AttemptedExchange = 0
AcceptedExchange  = 0

# Initialise these to zero - we'll calculate once object exists
LL                  = 0.0
LogPrior            = 0.0

# Set up matrix for storing probabilities of proposals
ProposalProbability = zeros(MySimulation.NumOfProposals+1)

GradLL              = zeros(MySimulation.Model.NumOfParas,1)
GradLogPrior        = zeros(MySimulation.Model.NumOfParas,1)
HessianLL           = zeros(MySimulation.Model.NumOfParas, MySimulation.Model.NumOfParas)
HessianGradLogPrior = zeros(MySimulation.Model.NumOfParas, MySimulation.Model.NumOfParas)

# Define array of chain geometries for each point in the state space
ChainGeometry = Array(MarkovChainGeometry, MySimulation.NumOfProposals + 1)
for i = 1:MySimulation.NumOfProposals + 1
    # Initialise parameters either from prior or from default paras
    if MySimulation.InitialiseFromPrior
        # NEED TO UPDATE!!!
        Parameters = copy(MySimulation.Model.DefaultParas)
    else
        Parameters = copy(MySimulation.Model.DefaultParas)
    end

    ChainGeometry[i] = MarkovChainGeometry(Parameters, LL, LogPrior, ProposalProbability, GradLL, GradLogPrior, HessianLL, HessianGradLogPrior)
end



# Set up Markov chain proposal variables
println("Setting up Markov chain proposal...")

StepSize = MySimulation.InitialStepSize

if Sampler == "MH"

    # Set the proposal distribution for each of these samplers
    if MySimulation.Model.NumOfParas == 1
        # Univariate
        Density = Normal(0, MySimulation.ProposalCovariance[1])
    else
        # Set the multivariate proposal distribution for each of these samplers
        Density = MvNormal(zeros(MySimulation.Model.NumOfParas), StepSize*MySimulation.ProposalCovariance)
    end

    # Create proposal distribution object for MH
    ProposalDistribution = ProposalDistributionMH(Density, StepSize, MySimulation.ProposalCovariance)

elseif Sampler == "SmMALA"  || Sampler == "TrSmMALA"

    # Set the proposal distribution for each of these samplers
    if MySimulation.Model.NumOfParas == 1
        # Univariate
        Density = Normal(0, MySimulation.ProposalCovariance[1])
    else
        # Set the multivariate proposal distribution for each of these samplers
        Density = MvNormal(zeros(MySimulation.Model.NumOfParas), MySimulation.ProposalCovariance)
    end

    # Create proposal distribution object for SmMALA and TrSmMALA
    ProposalDistribution = Array(ProposalDistributionSmMALA, MySimulation.NumOfProposals + 1)
    for i = 1:(MySimulation.NumOfProposals + 1)
      ProposalDistribution[i] = ProposalDistributionSmMALA(Density, StepSize)
    end

elseif Sampler == "TrSmMALA_Random"

    # Set the proposal distribution for each of these samplers
    if MySimulation.Model.NumOfParas == 1
        # Univariate
        Density = Normal(0, MySimulation.ProposalCovariance[1])
    else
        # Set the multivariate proposal distribution for each of these samplers
        Density = MvNormal(zeros(MySimulation.Model.NumOfParas), MySimulation.ProposalCovariance)
    end

    NumberOfVectors   = MySimulation.Model.AuxiliaryVars
    TangentVectors = zeros(MySimulation.Model.NumOfParas, NumberOfVectors)

    # Create proposal distribution object for SmMALA and TrSmMALA
    ProposalDistribution = Array(ProposalDistributionSmMALARandom, MySimulation.NumOfProposals + 1)
    for i = 1:(MySimulation.NumOfProposals + 1)
      ProposalDistribution[i] = ProposalDistributionSmMALARandom(Density, StepSize, NumberOfVectors, TangentVectors)
    end

elseif Sampler == "AdaptiveMH"
    # Set the proposal distribution for each of these samplers
    if MySimulation.Model.NumOfParas == 1
        # Univariate
        Density = Normal(0, MySimulation.ProposalCovariance[1])
    else
        # Set the multivariate proposal distribution for each of these samplers
        Density = MvNormal(zeros(MySimulation.Model.NumOfParas), MySimulation.ProposalCovariance)
    end

    RunningMean = zeros(MySimulation.Model.NumOfParas)
    RunningCov  = zeros(MySimulation.Model.NumOfParas)

    # Create proposal distribution object for adaptive MH
    ProposalDistribution = ProposalDistributionAdaptiveMH(Density, StepSize, RunningMean, RunningCov)

elseif Sampler == "Gibbs" || Sampler == "BivariateStructuredMarginal" || Sampler == "BivariateStructuredMetropolis"
    # Custom proposal so leave empty
    ProposalDistribution = []
end

CurrentIteration = 1;

# Create the chain object
println("Creating chain object...")
Chain = MarkovChain( Sampler,
                     NumOfProposals,
                     SampleIndex,
                     CurrentIteration,
                     AttemptedProposal,
                     AcceptedProposal,
                     AttemptedExchange,
                     AcceptedExchange,
                     ChainGeometry,
                     ProposalDistribution)


# Initialise the LL and any other quantities needed for the chosen sampler
Initialise(MySimulation.Model, Chain)

# Initialise variable for storing results
#SavedSamples = Array(Float64, MySimulation.Model.NumOfParas, (MySimulation.NumOfIterations*Chain.NumOfProposals))
TotalNumOfSamples = MySimulation.NumOfIterations*Chain.NumOfProposals;
SaveIteration = min(100000, TotalNumOfSamples)
SavedSamples = Array(Float64, MySimulation.Model.NumOfParas, SaveIteration)

################
# Main routine #
################
println("Starting main MCMC routine...")
for IterationNum = 1:MySimulation.NumOfIterations

    # Keep track of the iteration number within the chain
    Chain.CurrentIteration = IterationNum;

    #println(IterationNum)

    if mod(IterationNum, 1000) == 0
        println(IterationNum)

        # Calculate acceptance rate
        Counter = 0;
        for i = 2:Idx
            if SavedSamples[:, i] != SavedSamples[:, i-1]
              Counter = Counter + 1;
            end
        end

        # Print acceptance rate
        println(Counter/(Idx-1))
    end

    # Propose parameters and caluclate new geometry - Gibbs step 1
    UpdateParameters(MySimulation.Model, Chain)

    # Sample the indicator variable to select samples - Gibbs step 2
    IndicatorSamples = SampleIndicator(MySimulation.Model, Chain)

    #println(IndicatorSamples)
    #println(Chain.SampleIndicator)

    # Save the samples
    # ONLY VALID FOR SINGLE PROPOSAL!
    if mod(IterationNum,SaveIteration) == 0
        PartNum = IterationNum/SaveIteration;
        Idx = SaveIteration
        SavedSamples[:, Idx] = Chain.Geometry[ IndicatorSamples[1] ].Parameters

        # Output the results to a file
        writedlm(string(MySimulation.OutputID, "_Part", PartNum), SavedSamples, ',')
        Idx = 1;
    else
        Idx = mod(IterationNum, SaveIteration);
        SavedSamples[:, Idx] = Chain.Geometry[ IndicatorSamples[1] ].Parameters
    end



end

########################
# Post-simulation work #
########################

# Print summary statistics of run
for i = 1:MySimulation.Model.NumOfParas
    println("")
    println( string("Parameter ", i, " Mean: ", mean(SavedSamples[i,:])) )
    println( string("Parameter ", i, " Var: ", var(SavedSamples[i,:])) )
    println( string("Parameter ", i, " Std: ", sqrt(var(SavedSamples[i,:]))) )
end


  # Output the remaining results to a file
  writedlm(string(MySimulation.OutputID, "_Part_Final"), SavedSamples, ',')

end
