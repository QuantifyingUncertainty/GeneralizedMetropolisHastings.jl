
abstract AbstractModel





# Define Markov chain object
type MarkovChain

    NumOfProposals::Int64
    SampleIndicator::Int64
    CurrentIteration::Int64


    Samples::Vector{S<:Sample}
    Sampler::Vector{S<:Sampler}

end

# Define the MCMC Simulation object
type MCMCSimulation

    OutputID
    Sampler # e.g. MH, MALA, SmMALA, TrSmMALA etc
    NumOfProposals::Int64
    NumOfIterations::Int64
    InitialStepSize::Float64
    InitialiseFromPrior::Bool

    # This can be any model, as long as corresponding UpdateParameters function is defined
    Model::AbstractModel

end
