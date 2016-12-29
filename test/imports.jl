import Base.Test: @test,@test_approx_eq,@test_approx_eq_eps,@test_throws

@everywhere import GeneralizedMetropolisHastings:
    ###Types
    MHRunnerType,InitializeFrom,ProposeFrom,IndicatorType,JobSegments,ChainType,StoreDuring,
    MHRuntimePolicy,
    ###Types
    AbstractParameter,ParameterDefault,ParameterUnivariate,
    DataArray,DataFunction,
    AbstractBactrian,Bactrian,
    AbstractNoiseModel,NoiseModelGaussian,
    BaseSample,GradientSample,TensorSample,TangentTensorSample,
    SymmetricDensity,ASymmetricDensity,DistributionWrapper,NormalDensity,LogNormalDensity,
    MonitorTuner,ScaleTuner,
    AbstractModel,TargetModel,ODEModel,
    ChainStandard,ChainGradient,
    IndicatorStationary,
    GMHSegment,MHRemoteSegments,
    AbstractMHRunner,SMHRunner,GMHRunner,
    ###Functions
    showsamplerstatevars, #from util.jl
    trait,traitvalue,traittype,policy, #form policies.jl
    parameter,parameters,initvalues!,initvalues,logprior!,logprior, #from parameters.jl
    data,numvalues,numvars,generate!,dataindex,datavalues, #from data.jl
    noise,loglikelihood,applynoise!, #from noise.jl
    samples,numparas,numsamples,numtangents,sampletype,calculationtype,similar,offset!, #from samples.jl
    distribution,distributions,recenter,rescale,update, #from distributions.jl
    density,condition!,scale!,update!,propose!,logprobability,logprobability!, #from densities.jl
    sampler,samplerstate,setfrom!,propose!,acceptance!,acceptance,tune!,from,proposals, #from samplers.jl
    prepare!,prepareindicator!,prepareauxiliary!,getsamplerstatevars,#from samplers.jl
    tuner,tunerstate,rate,accepted,proposed,total,index,numtunesteps,current,tune,nextindex!,period,verbose,needstuning,accepted!,showstep, #from tuners.jl
    model,geometry!,evaluate!,loglikelihood,measurements,noisemodel,#from models.jl
    chain,store!,accepted!,logposterior, #from chains.jl
    indicator,numproposals,transitionprobability!,sampleindicator!,indicatorsamples,accepted, #from indicators.jl
    segment,numproposals,iterate!,getsamples, #from jobsegments.jl
    remotesegments,numsegments,numproposalspersegment,numtotalproposals,iterate!,retrievesamples!,
    runner,run!,initialize!,preparenext!,auxiliary!
