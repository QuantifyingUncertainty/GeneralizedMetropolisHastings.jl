module GeneralizedMetropolisHastings

using Compat

import StatsFuns
import Distributions
import Sundials

import Base: ==
import Base: size, length, eltype, show, display, time, copy!

import Distributions: MvNormal,MvLogNormal
import Distributions: rand!,logpdf!,logpdf,location,location!

export
    ###Policy trait types
    AbstractPolicyTrait,InitializeFrom,ProposeFrom,IndicatorType,SamplerStates,
    AbstractPolicy,GMHPolicy,
    ###Types
    AbstractParameter,ParameterDefault,ParameterUnivariate,
    AbstractData,DataArray,DataFunction,
    AbstractNoiseModel,NoiseModelGaussian,
    AbstractSample,BaseSample,GradientSample,TensorSample,TangentTensorSample,
    AbstractProposalDensity,SymmetricDensity,ASymmetricDensity,DistributionWrapper,NormalDensity,LogNormalDensity,
    AbstractSampler,MHNormal,MHLogNormal,
    AbstractSamplerState,
#   MHNormal,SmMALANormal,TrSmMALANormal,TrSmMALARandomNormal,
    AbstractTuner,MonitorTuner,ScaleTuner,
    AbstractTunerState,
    AbstractModel,TargetModel,ODEModel,
    AbstractChain,ChainStandard,ChainGradient,
    AbstractIndicatorMatrix,IndicatorStationary,
    AbstractRunner,GMHRunner,
  ###Functions
    trait,traitvalue,traittype,policy, #form policies.jl
    parameter,parameters,initvalues!,initvalues,logprior!,logprior, #from parameters.jl
    data,numvalues,numvars,generate!,dataindex,datavalues, #from data.jl
    noise,loglikelihood,applynoise!, #from noise.jl
    sample,samples,numparas,numsamples,numtangents,copy!, #from samples.jl
    density,condition!,propose!,logprobability,logprobability!,issymmetric, #from densities.jl
    sampler,samplerstate,setfrom!,propose!,acceptanceratio!,acceptanceratio,tune!,from,proposals, #from samplers
    tuner,tunerstate,rate,accepted,proposed,totalproposed,tune,resetburnin!,period,verbose, #from tuners.jl
    model,geometry!,evaluate!,loglikelihood,measurements,#from models.jl
    chain,store!,accepted!,logposterior, #from chains.jl
    indicator,numproposals,transitionprobability!,sampleindicator!,indicatorsamples,accepted, #from indicators.jl
    runner,samplerstates,initialize!,iterate!,updatefrom!,iterateandstore!,iterateandtune!,run!,
    print_gmh_module_loaded

include("policies/traits.jl")
include("policies/policies.jl")
include("policies/GMHPolicy.jl")
include("parameters/parameters.jl")
include("data/data.jl")
include("data/DataArray.jl")
include("data/DataFunction.jl")
include("noise/noise.jl")
include("noise/NoiseModelGaussian.jl")
include("samples/samples.jl")
include("samples/base.jl")
include("samples/gradient.jl")
include("samples/tensor.jl")
include("proposals/proposals.jl")
include("proposals/DistributionWrapper.jl")
include("proposals/NormalDensity.jl")
include("proposals/LogNormalDensity.jl")
include("samplers/samplers.jl")
include("samplers/MetropolisHastings.jl")
# include("samplers/SmMALA.jl")
include("indicators/indicators.jl")
include("indicators/IndicatorStationary.jl")
include("models/models.jl")
include("models/TargetModel.jl")
include("models/ODEModel.jl")
include("tuners/tunefunctions.jl")
include("tuners/tuners.jl")
include("tuners/MonitorTuner.jl")
include("tuners/ScaleTuner.jl")
include("chains/chains.jl")
include("chains/ChainStandard.jl")
include("chains/ChainGradient.jl")
include("runners/runners.jl")
include("runners/GMHRunner.jl")

function print_gmh_module_loaded()
  println("$module_name(current_module()) module loaded successfully")
end

end # module
