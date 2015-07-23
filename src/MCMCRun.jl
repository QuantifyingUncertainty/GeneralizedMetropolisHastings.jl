function MCMCRun( MySimulation::MCMCSimulation )

  stepsize = MySimulation.InitialStepSize
  numparas = MySimulation.Model.NumOfParas
  numproposals = MySimulation.NumOfProposals

#####################
# Initialise chains #
#####################

  # Define array of chain geometries for each point in the state space
  ChainGeometry = Array(MarkovChainGeometry,numproposals + 1)

  for i = 1:numproposals + 1

    # Initialise parameters either from prior or from default paras
    if MySimulation.InitialiseFromPrior
        # NEED TO UPDATE!!!
        Parameters = copy(MySimulation.Model.DefaultParas)
    else
        Parameters = copy(MySimulation.Model.DefaultParas)
    end

    ChainGeometry[i] = MarkovChainGeometry(numproposals,numparas)
  end

#######################
# Initialise samplers #
#######################

  # Set up Markov chain proposal variables
  println("Setting up Markov chain proposal...")

  Sampler = Array{AbstractSampler,numproposals + 1}
  for i = 1:numproposals + 1
    Sampler[i] = create_sampler
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
    IndicatorSamples = SampleIndicator(Chain)

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
