module GeneralizedMetropolisHastings

import Compat

import PDMats
import Distributions
import Sundials

import Base: ==, size, length, eltype, show, display, time, similar, copy!, copy
import Base: mean, var, rand, rand!, scale!

import Distributions: pdf,logpdf!,logpdf,location,location!,scale

export
    ###Abstract Types
    AbstractPolicyTrait,AbstractPolicy,
    AbstractParameter,
    AbstractData,
    AbstractNoiseModel,
    AbstractSample,
    AbstractProposalDensity,SymmetricDensity,ASymmetricDensity,
    AbstractSampler,AbstractSamplerState,
    AbstractTuner,AbstractTunerState,
    AbstractModel,
    AbstractChain,
    AbstractIndicatorMatrix,
    AbstractJobSegment,AbstractRemoteSegments,
    AbstractRunner,AbstractMHRunner,
    ###Functions
    trait,traitvalue,traittype,
    policy,
    parameters,
    data,dataindex,datavalues,
    distribution,distributions,
    density,
    noise,
    sampler,
    tuner,
    model,evaluate!,applynoise!,measurements,numparas,
    runner,
    run!,
    samples,numsamples,loglikelihood,logprior,logposterior,
    print_gmh_module_loaded

include("policies/traits.jl")
include("policies/policies.jl")
include("policies/MHRuntimePolicy.jl")
include("parameters/parameters.jl")
include("data/data.jl")
include("data/DataArray.jl")
include("data/DataFunction.jl")
include("distributions/distributions.jl")
include("distributions/Bactrian.jl")
include("proposals/proposals.jl")
include("proposals/DistributionWrapper.jl")
include("noise/noise.jl")
include("noise/NoiseModelGaussian.jl")
include("samples/samples.jl")
include("samples/base.jl")
include("samples/gradient.jl")
include("samples/tensor.jl")
include("samplers/samplers.jl")
include("samplers/MetropolisHastings.jl")
include("samplers/AdaptiveMetropolis.jl")
include("samplers/SmMALA.jl")
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
include("jobsegments/jobsegments.jl")
include("jobsegments/GMHSegment.jl")
include("jobsegments/remotesegments.jl")
include("jobsegments/MHRemoteSegments.jl")
include("runners/runners.jl")
include("runners/MHRunner.jl")
include("runners/SMHRunner.jl")
include("runners/GMHRunner.jl")
include("testutil/util.jl")
include("testutil/sincos.jl")
include("testutil/springmass.jl")

function print_gmh_module_loaded()
  println("$module_name(current_module()) module loaded successfully")
end

end # module
