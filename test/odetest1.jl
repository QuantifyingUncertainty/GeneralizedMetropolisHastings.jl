using GeneralizedMetropolisHastings
using Distributions

# Define settings for simulation
OutputID            = "TestOutputFile"
Sampler             = "MH"
NumOfProposals      = 2
NumOfIterations     = 1000
InitialStepSize     = 0.01
ProposalCovariance  = eye(3);
InitialiseFromPrior = false # Sample starting parameters from prior

#######################
# Create an ODE model #
#######################

# Define the ODEModel object
ModelType    = "ODE"
ModelName    = "FHN"
NumOfParas   = 3
ParaNames    = ["a", "b", "c"]
DefaultParas = [0.2, 0.2, 3.0]

Prior        = Array(Distribution, NumOfParas)
Prior[1]     = Uniform(0, 10)
Prior[2]     = Uniform(0, 10)
Prior[3]     = Uniform(0, 10)

# Data
DataMeasurements  = [1.0 1.0;
                     1.0 1.0;
                     1.0 1.0]
DataTimePoints    = [0.2 0.2;
                     0.5 0.5;
                     0.8 0.8]
DataNoiseVariance = ones(3,2) # Same size as the dataset
DataNoiseVariance[:,1] = 0.0215
DataNoiseVariance[:,2] = 0.0048

NumOfTimePoints = 3

# ODE specific
# Specify differential equation function
function FHN(t,y,ydot,p)
    V = y[1]
    R = y[2]

    a = p[1]
    b = p[2]
    c = p[3]

    dV = c*(V-(V^3)/3+R)
    dR = -(V-a+b*R)/c

    ydot[1] = dV
    ydot[2] = dR
end

FHNFun                   = FHN
DefaultInitialConditions = [-1.0, 1.0]
NumOfStates             = 2
ObservedStates          = [1, 2]
UnobservedStates        = []
abstol                   = 1e-8
reltol                   = 1e-8
InferInitialConditions   = false

AuxiliaryVars = 20; # Number of tangent samples (We'll store this in the model for now)

MyODEModel = ODEModel( ModelType,
                       ModelName,
                       NumOfParas,
                       ParaNames,
                       DefaultParas,
                       Prior,
                       DataMeasurements,
                       DataTimePoints,
                       DataNoiseVariance,
                       NumOfTimePoints,
                       FHNFun,
                       DefaultInitialConditions,
                       NumOfStates,
                       ObservedStates,
                       UnobservedStates,
                       abstol,
                       reltol,
                       InferInitialConditions,
                       AuxiliaryVars)

############################
# Create simulation object #
############################

MySimulation = MCMCSimulation( OutputID,
                               Sampler,
                               NumOfProposals,
                               NumOfIterations,
                               InitialStepSize,
                               ProposalCovariance,
                               InitialiseFromPrior,
                               MyODEModel )

# Run the MCMC code, passing in the MCMCSimulation object
MCMCRun( MySimulation )
