### Standard Metropolis-Hastings Monte Carlo runner
immutable SMHRunner <: AbstractMHRunner
    numburnin::Int
    numiterations::Int
    policy::MHRuntimePolicy #policies determining runtime behaviour of the MCMC

    function SMHRunner(nburnin::Int,niterations::Int,policy::MHRuntimePolicy)
        @assert 0 <= nburnin && 0 < niterations
        new(nburnin,niterations,policy)
    end
end

###Factory functions
_runner(::Type{Val{:standard}},policy::MHRuntimePolicy,niterations::Int;numburnin =0) = SMHRunner(numburnin,niterations,policy)

function run!(runner_::SMHRunner,model_::AbstractModel,sampler_::AbstractSampler,tuner_::AbstractTuner)
    samplerstate_ = samplerstate(sampler_,1,runner_.policy.sampletype,runner_.policy.calculationtype)
    indicator_,tunerstate_,chain_ = createcommon(runner_,tuner_,numparas(model_),1,1)
    tic()
    burnin!(runner_,model_,sampler_,samplerstate_,tuner_,tunerstate_,indicator_,chain_)
    main!(runner_,model_,sampler_,samplerstate_,tuner_,tunerstate_,indicator_,chain_)
    chain_.runtime = toq()
    chain_
end

function burnin!(runner_::SMHRunner,model_::AbstractModel,sampler_::AbstractSampler,samplerstate_::AbstractSamplerState,
                 tuner_::AbstractTuner,tunerstate_::AbstractTunerState,indicator_::AbstractIndicatorMatrix,chain_::AbstractChain)
    storeduring = _storeduring(:burnin,runner_.policy)
    initialize!(runner_,model_,samplerstate_,chain_,storeduring)
    for i=1:runner_.numburnin
        iterate!(runner_,model_,samplerstate_,indicator_,chain_,storeduring)
        accepted!(tunerstate_,indicator_)
        needstuning(tuner_,i)?tune!(runner_,samplerstate_,tuner_,tunerstate_):nothing
        preparenext!(runner_,samplerstate_,indicator_) #prepare for the next iteration
    end
end

function main!(runner_::SMHRunner,model_::AbstractModel,sampler_::AbstractSampler,samplerstate_::AbstractSamplerState,
               tuner_::AbstractTuner,tunerstate_::AbstractTunerState,indicator_::AbstractIndicatorMatrix,chain_::AbstractChain)
    storeduring = _storeduring(:main,runner_.policy)
    for i=1:runner_.numiterations
        iterate!(runner_,model_,samplerstate_,indicator_,chain_,storeduring)
        needstuning(tuner_,i)?println("Iteration $i/$(runner_.numiterations)"):nothing
        preparenext!(runner_,samplerstate_,indicator_)
    end
end

function iterate!(runner_::SMHRunner,model_::AbstractModel,samplerstate_::AbstractSamplerState,
                  indicator_::AbstractIndicatorMatrix,chain_::AbstractChain,storeduring::Bool)
    propose!(samplerstate_)
    geometry!(model_,proposals(samplerstate_))
    a = acceptance!(samplerstate_)
    transitionprobability!(indicator_,a)
    sampleindicator!(indicator_)
    storeduring?store!(runner_,samplerstate_,indicator_,chain_):nothing
end

function store!(runner_::SMHRunner,samplerstate_::AbstractSamplerState,
                indicator_::AbstractIndicatorMatrix,chain_::AbstractChain)
    if indicatorsamples(indicator_)[2] == 2
        store!(chain_,from(samplerstate_),1)
    else
        store!(chain_,proposals(samplerstate_),1)
    end
    accepted!(chain_,indicator_)
end

preparenext!(runner_::SMHRunner,samplerstate_::AbstractSamplerState,indicator_::AbstractIndicatorMatrix) = prepare!(samplerstate_,indicatorsamples(indicator_)[2]!=2)

function tune!(runner_::SMHRunner,samplerstate_::AbstractSamplerState,tuner_::AbstractTuner,tunerstate_::AbstractTunerState)
    tvals = tune(tuner_,tunerstate_)
    showstep(tuner_,tunerstate_)
    tune!(samplerstate_,tvals...)
    nextindex!(tunerstate_)
end

function show(io::IO,r::SMHRunner)
    println(io,"Standard Metropolis-Hastings runner with:")
    println(io," numburnin: ",r.numburnin)
    println(io," numiterations: ",r.numiterations)
    print(io," policy: ")
    show(io,r.policy)
    println(io)
    nothing
end

