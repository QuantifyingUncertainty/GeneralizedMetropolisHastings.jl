### Generalized Metropolis-Hastings Monte Carlo runner
immutable GMHRunner <: AbstractMHRunner
    numburnin::Int
    numiterations::Int
    numproposals::Int
    numindicatorsamples::Int
    policy::MHRuntimePolicy #policies determining runtime behaviour of the MCMC

    function GMHRunner(nburnin::Int,niterations::Int,nproposals::Int,nindsamples::Int,policy::MHRuntimePolicy)
        @assert 0 <= nburnin && 0 < niterations && 0 < nproposals && 0 < nindsamples
        new(nburnin,niterations,nproposals,nindsamples,policy)
    end
end

###Factory functions
_runner(::Type{Val{:generalized}},policy::MHRuntimePolicy,niterations::Int,nproposals::Int;numburnin =0,numindicatorsamples =nproposals) = GMHRunner(numburnin,niterations,nproposals,numindicatorsamples,policy)

function run!(runner_::GMHRunner,model_::AbstractModel,sampler_::AbstractSampler,tuner_::AbstractTuner)
    indicatorstate_ = samplerstate(sampler_,1,runner_.policy.sampletype,runner_.policy.calculationtype)
    segments_ = remotesegments(runner_.policy,model_,sampler_,runner_.numproposals)
    show(segments_)
    indicator_,tunerstate_,chain_ = createcommon(runner_,tuner_,numparas(model_),numtotalproposals(segments_),runner_.numindicatorsamples)
    tic()
    burnin!(runner_,model_,sampler_,indicatorstate_,segments_,tuner_,tunerstate_,indicator_,chain_)
    main!(runner_,model_,sampler_,indicatorstate_,segments_,tuner_,tunerstate_,indicator_,chain_)
    chain_.runtime = toq()
    chain_
end


function burnin!(runner_::GMHRunner,model_::AbstractModel,sampler_::AbstractSampler,indicatorstate_::AbstractSamplerState,segments_::AbstractRemoteSegments,
                 tuner_::AbstractTuner,tunerstate_::AbstractTunerState,indicator_::AbstractIndicatorMatrix,chain_::AbstractChain)
    storeduring = _storeduring(:burnin,runner_.policy)
    initialize!(runner_,model_,from(indicatorstate_),chain_,storeduring)
    for i=1:runner_.numburnin
        iterate!(runner_,model_,indicatorstate_,segments_,indicator_,chain_,storeduring)
        accepted!(tunerstate_,indicator_)
        needstuning(tuner_,i)?tune!(runner_,sampler_,indicatorstate_,segments_,tuner_,tunerstate_):nothing
    end
end

function main!(runner_::GMHRunner,model_::AbstractModel,sampler_::AbstractSampler,indicatorstate_::AbstractSamplerState,segments_::AbstractRemoteSegments,
               tuner_::AbstractTuner,tunerstate_::AbstractTunerState,indicator_::AbstractIndicatorMatrix,chain_::AbstractChain)
    storeduring = _storeduring(:main,runner_.policy)
    for i=1:runner_.numiterations
        iterate!(runner_,model_,indicatorstate_,segments_,indicator_,chain_,storeduring)
        needstuning(tuner_,i)?println("Iteration $i/$(runner_.numiterations)"):nothing
    end
end

function iterate!(runner_::GMHRunner,model_::AbstractModel,indicatorstate_::AbstractSamplerState,segments_::AbstractRemoteSegments,
                  indicator_::AbstractIndicatorMatrix,chain_::AbstractChain,storeduring::Bool)
    auxiliary!(runner_,model_,indicatorstate_)
    segmentacceptances = iterate!(segments_,proposals(indicatorstate_))
    transitionprobability!(indicator_,acceptanceratio!(indicatorstate_),segmentacceptances)
    try
        sampleindicator!(indicator_)
    catch e
        dump(indicatorstate_)
        dump(segments_)
        throw(e)
    end
    storeduring?store!(runner_,indicatorstate_,segments_,indicator_,chain_):nothing
    updatefrom!(runner_,indicatorstate_,segments_,indicator_)
end

function auxiliary!(runner_::GMHRunner,model_::AbstractModel,indicatorstate_::AbstractSamplerState)
    auxcounter = 0
    while ~isfinite((propose!(indicatorstate_) ; geometry!(model_,proposals(indicatorstate_))).logprior[1])
        if (auxcounter+=1) > maxinitialcounter
            error("Problems generating auxiliary proposal, giving up. Please check the parameter priors.")
        end
    end
    indicatorstate_
end

function store!(runner_::GMHRunner,indicatorstate_::AbstractSamplerState,segments_::AbstractRemoteSegments,
                indicator_::AbstractIndicatorMatrix,chain_::AbstractChain)
    nproposals = numproposals(indicator_)
    nsamples = numsamples(indicator_)
    indicatorsamples_ = indicatorsamples(indicator_)
    retrievesamples!(segments_,indicatorsamples_)
    from_ = from(indicatorstate_)
    for i=2:nsamples+1
        @inbounds if indicatorsamples_[i] == nproposals + 1
            store!(chain_,from_,1)
        else
            @inbounds store!(segments_,chain_,indicatorsamples_[i])
        end
    end
    accepted!(chain_,indicator_)
end

function updatefrom!(runner_::GMHRunner,indicatorstate_::AbstractSamplerState,
                     segments_::AbstractRemoteSegments,indicator_::AbstractIndicatorMatrix)
    indicatorend = indicatorsamples(indicator_)[end]
    if indicatorend != numproposals(indicator_) + 1
        setfrom!(indicatorstate_,getsamples(segments_,indicatorend))
    end
    indicatorstate_
end

function tune!(runner_::AbstractRunner,sampler_::AbstractSampler,indicatorstate_::AbstractSamplerState,
               segments_::AbstractRemoteSegments,tuner_::AbstractTuner,tunerstate_::AbstractTunerState)
    tvals = tune(tuner_,tunerstate_)
    showstep(tuner_,tunerstate_)
    tune!(sampler_,indicatorstate_,tvals...)
    tune!(segments_,sampler_,tvals...)
    nextindex!(tunerstate_)
end

function show(io::IO,r::GMHRunner)
    println(io,"Generalized Metropolis-Hastings runner with:")
    println(io," numburnin: ",r.numburnin)
    println(io," numiterations: ",r.numiterations)
    println(io," numproposals: ",r.numproposals)
    println(io," numindicatorsamples: ",r.numindicatorsamples)
    print(io," policy: ")
    show(io,r.policy)
    println(io)
    nothing
end
