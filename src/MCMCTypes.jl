# Define the statistical model object
type ODEModel

    ModelType
    ModelName
    NumOfParas::Int64
    ParaNames::Array
    DefaultParas::Array{Float64}
    Priors::Array{Any}

    # Data
    DataMeasurements::Array{Float64}
    DataTimePoints::Array{Float64}
    DataNoiseVariance::Array{Float64}
    NumOfTimePoints::Int64

    # ODE specific
    ODEFunction::Function
    DefaultInitialConditions::Array{Float64}
    NumOfSpecies::Int64
    ObservedSpecies::Array{Int64}
    UnobservedSpecies::Array{Int64}
    abstol::Float64
    reltol::Float64
    InferInitialConditions::Bool

    AuxiliaryVars::Any

end

# Define the statistical model object
type TargetOnlyModel

    ModelType
    ModelName
    NumOfParas::Int64
    ParaNames::Array
    DefaultParas::Array{Float64}
    Priors::Array{Any}
    LLEval::Function
end


# Define Markov chain geometry object
type MarkovChainGeometry

    Parameters::Array{Float64}

    LL::Float64
    LogPrior::Float64
    ProposalProbability::Array{Float64}

    GradLL::Array{Float64}              # vector
    GradLogPrior::Array{Float64}        # vector
    HessianLL::Array{Float64}           # matrix
    HessianGradLogPrior::Array{Float64} # matrix

end


# Define proposal object for M-H MCMC
type ProposalDistributionMH

    Density::Any                        # This will be a distribution type
    StepSize::Float64
    ProposalCov::Array{Float64}

end

# Define proposal object for SmMALA MCMC
type ProposalDistributionSmMALA

    Density::Any                        # This will be a distribution type
    StepSize::Float64

end

# Define proposal object for SmMALA MCMC
type ProposalDistributionSmMALARandom

    Density::Any                        # This will be a distribution type
    StepSize::Float64
    NumberOfVectors::Int64
    TangentVectors::Array{Float64}

end

# Define proposal object for Adaptive MCMC
type ProposalDistributionAdaptiveMH

    Density::Any                        # This will be a distribution type
    StepSize::Float64
    RunningMean::Array{Float64}
    RunningCov::Array{Float64}

end


# Define Markov chain object
type MarkovChain

    Sampler
    NumOfProposals::Int64
    SampleIndicator::Int64
    CurrentIteration::Int64

    AttemptedProposal::Int64
    AcceptedProposal::Int64
    AttemptedExchange::Int64
    AcceptedExchange::Int64

    Geometry::Array{MarkovChainGeometry}

    ProposalDistribution::Any           # Distribution for the proposal

end


# Define the MCMC Simulation object
type MCMCSimulation

    OutputID
    Sampler # e.g. MH, MALA, SmMALA, TrSmMALA etc
    NumOfProposals::Int64
    NumOfIterations::Int64
    InitialStepSize::Float64
    ProposalCovariance::Array{Float64}
    InitialiseFromPrior::Bool

    # This can be any model, as long as corresponding UpdateParameters function is defined
    Model
end
