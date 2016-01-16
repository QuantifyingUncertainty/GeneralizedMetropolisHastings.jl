module GeneralizedMetropolisHastings

using Compat

import Distributions
import Sundials

import Base:
  ==

export
  ###Policy types
  InitializeFrom,InitializeFromDefault,InitializeFromPrior,
  ProposeFrom,ProposeFromIndicator,ProposeFromAuxiliary,
  GenerateIndicator,IndicatorStationary,IndicatorCyclical,
  RuntimePolicy,GenericPolicy,
  ###Types
  Parameter,ParameterDefault,ParameterPrior,
  Sample,BaseSample,GradientSample,TensorSample,TangentTensorSample,
#   NormalDensity,LogNormalDensity,
#   MHNormal,SmMALANormal,TrSmMALANormal,TrSmMALARandomNormal,
#   MHHeap,SmMALAHeap,
#   TargetModel,ODEModel,
#   MCChain,
#   GMHRunner,
  ###Functions
  policy, #form policies.jl
  parameter,parameters,initvalues!,initvalues,logprior, #from parameters.jl
  sample,samples,numparas,numsamples,numtangents, #from samplers.jl
  getvalue,getloglikelihood,getlogprior,getgradloglikelihood,getgradlogprior,gettensorloglikelihood,gettensorlogprior,gettangent,
  setvalue!,setloglikelihood!,setlogprior!,setgradloglikelihood!,setgradlogprior!,settensorloglikelihood!,settensorlogprior!,settangent!,
#   logprobability, #from densities.jl
#   evaluate,loglikelihood,gradloglikelihood,tensorvalue, #specify for specific models
#   logprior!,gradlogprior!,tensorlogprior!, #from models.jl
#   loglikelihood!,gradient!,tensor!,update_geometry!,
#   initialize!,iterate!,run!,sample_indicator, #form runners.jl
  print_gmh_module_loaded

include("policies/policies.jl")
include("policies/GenericPolicy.jl")
include("parameters/parameters.jl")
include("samples/samples.jl")
include("samples/BaseSample.jl")
include("samples/GradientSample.jl")
include("samples/TensorSample.jl")
include("samples/TangentTensorSample.jl")

# include("api/chains.jl")
# include("distributions/Normal.jl")
# include("distributions/LogNormal.jl")
# include("proposals/proposals.jl")
# include("proposals/SymmetricProposal.jl")
# include("proposals/AsymmetricProposal.jl")
# include("samplers/samplers.jl")
# include("samplers/MetropolisHastings.jl")
# include("samplers/SmMALA.jl")
# include("models/models.jl")
# include("models/TargetModel.jl")
# include("models/ODEModel.jl")
# include("runners/indicator.jl")
# include("runners/GMHRunner.jl")

end # module
