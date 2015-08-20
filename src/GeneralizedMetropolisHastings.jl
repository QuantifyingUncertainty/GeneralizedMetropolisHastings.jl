module GeneralizedMetropolisHastings

using Distributions
using Sundials

if VERSION < v"0.4.0-dev"
    using Docile
end

@docstrings

#imports

export
  ###Policy types
  ValuesFromDefault,ValuesFromPrior,
  ProposalFromIndicator,ProposalFromAuxiliary,
  IndicatorMatrixStationary,IndicatorMatrixOptimal,
  GenericPolicy,
  ###Types
  ModelParameters,
  BaseSample,GradientSample,TensorSample,ApproximateTensorSample,
  NormalDensity,
  MHNormal,SmMALANormal,TrSmMALANormal,TrSmMALARandomNormal,
  MHHeap,SmMALAHeap,
  TargetModel,ODEModel,
  MCChain,
  GMHRunner,
  ###Functions
  values,named,anonymous,numel, #from parameters.jl
  numparas, #from samplers.jl
  logprobability, #from densities.jl
  evaluate,loglikelihood,gradloglikelihood,tensorvalue, #specify for specific models
  logprior!,gradlogprior!,tensorlogprior!, #from models.jl
  loglikelihood!,gradient!,tensor!,update_geometry!,
  initialize!,iterate!,run!,sample_indicator, #form runners.jl
  print_gmh_module_loaded

include("api/api.jl")
include("api/policies.jl")
include("api/parameters.jl")
include("api/samples.jl")
include("api/chains.jl")
include("densities/densities.jl")
include("densities/NormalDensity.jl")
include("samplers/samplers.jl")
include("samplers/MetropolisHastings.jl")
include("samplers/SmMALA.jl")
include("models/models.jl")
include("models/TargetModel.jl")
include("models/ODEModel.jl")
include("runners/indicator.jl")
include("runners/GMHRunner.jl")

end # module
