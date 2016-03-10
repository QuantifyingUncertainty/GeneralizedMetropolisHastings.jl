### Generalized Metropolis-Hastings Monte Carlo runner
immutable GMHRunner <: AbstractRunner
    numiterations::Int
    numproposals::Int
    numindicatorsamples::Int
    numtotalsamples::Int
    numburnin::Int
    policy::GMHPolicy #policies determining runtime behaviour of the MCMC

    function GMHRunner(numit::Int,numprop::Int,numindsamp::Int,numtotsamp::Int,numburn::Int,pol::GMHPolicy)
        @assert 0 < numit && 0 < numprop && 0 < numindsamp && 0 < numtotsamp && 0 <= numburn && numburn < numit
        @assert numit*numindsamp == numtotsamp
        new(numit,numprop,numindsamp,numtotsamp,numburn,pol)
    end
end

###Factory functions
_runner(::Type{Val{:gmh}},numit::Int,numprop::Int,pol::GMHPolicy;numburnin =0,numindicatorsamples =numprop) = GMHRunner(numit,numprop,numindicatorsamples,numit*numindicatorsamples,numburnin,pol)

###Functions used in all cases

###Generic entry points of the run function
run!(runner_::GMHRunner,model_::AbstractModel,sampler_::AbstractSampler;tuner::AbstractTuner =_tuner(Val{:monitor},100)) = run!(runner_,model_,sampler_,tuner)

function run!(runner_::GMHRunner,model_::AbstractModel,sampler_::AbstractSampler,tuner::AbstractTuner)
    chain_ = chain(:standard,numparas(model_),runner_.numtotalsamples,runner_.policy.sampletype,runner_.policy.calculationtype)
    tunerstate_ = tunerstate(tuner)
    samplerstates_ = samplerstates(runner_,sampler_)
    indicator_ = _indicator(traittype(runner_.policy.indicator),numsamples(samplerstates_[1].proposals)*length(samplerstates_),runner_.numproposals,runner_.policy.calculationtype)
    tic()
    run!(traittype(runner_.policy.propose),runner_,model_,sampler_,samplerstates_,tuner,tunerstate_,indicator_,chain_)
    chain_.runtime = toq()
    chain_
end

function run!(::Type{Val{:indicator}},runner_::GMHRunner,model_::AbstractModel,sampler_::AbstractSampler,samplerstates_::Vector,
              tuner_::AbstractTuner,tunerstate_::AbstractTunerState,indicator_::AbstractIndicatorMatrix,chain_::AbstractChain)
    @assert length(samplerstates_) == 1 && numsamples(samplerstates_[1].proposals) == 1
    initialproposalcounter = 0
    while ~isfinite(initialize!(runner_,model_,from(samplerstates_[1])).logprior[1])
        warn("Initial proposal falls outside parameter priors. Retrying...")
        if (initialproposalcounter+=1) > 100
            error("Difficulties generating initial proposal, giving up. Please check the default and prior values of the parameters.")
        end
    end
    for i=1:runner_.numburnin
        iterateandstore!(runner_,model_,samplerstates_,indicator_,chain_,false)
        accepted!(tunerstate_,indicator_)
        mod(i,period(tuner_))==0?tune!(runner_,sampler_,samplerstates_[1],[],tuner_,tunerstate_):nothing
    end
    for i=1:runner_.numiterations
        iterateandstore!(runner_,model_,samplerstates_,indicator_,chain_)
    end
    chain_
end

function run!(::Type{Val{:auxiliary}},runner_::GMHRunner,model_::AbstractModel,sampler_::AbstractSampler,samplerstates_::Vector,
              tuner_::AbstractTuner,tunerstate_::AbstractTunerState,indicator_::AbstractIndicatorMatrix,chain_::AbstractChain)
    @assert length(samplerstates_) > 0
    indicatorstate_ = samplerstate(sampler_,1,runner_.policy.sampletype,runner_.policy.calculationtype)
    initialproposalcounter = 0
    while ~isfinite(initialize!(runner_,model_,from(indicatorstate_)).logprior[1])
        warn("Initial proposal falls outside parameter priors. Retrying...")
        if (initialproposalcounter+=1) > 100
            error("Difficulties generating initial proposal, giving up. Please check the default and prior values of the parameters.")
        end
    end
    for i=1:runner_.numburnin
        iterateandstore!(runner_,model_,indicatorstate_,samplerstates_,indicator_,chain_,false)
        accepted!(tunerstate_,indicator_)
        mod(i,period(tuner_))==0?tune!(runner_,sampler_,indicatorstate_,samplerstates_,tuner_,tunerstate_):nothing
    end
    for i=1:runner_.numiterations
        iterateandstore!(runner_,model_,indicatorstate_,samplerstates_,indicator_,chain_)
    end
    chain_
end

function tune!(runner_::AbstractRunner,sampler_::AbstractSampler,indicatorstate_::AbstractSamplerState,samplerstates_::Vector,tuner_::AbstractTuner,tunerstate_::AbstractTunerState)
    tvals = tune(tuner_,tunerstate_)
    if verbose(tuner_)
        show(tunerstate_)
        println("   scalefactor=$tvals")
    end
    tune!(sampler_,indicatorstate_,tvals...)
    for i=1:length(samplerstates_)
        tune!(sampler_,samplerstates_[i],tvals...)
    end
    resetburnin!(tunerstate_)
end

function initialize!(runner_::GMHRunner,model_::AbstractModel,sample_::AbstractSample)
    initialize!(runner_.policy.initialize,model_,sample_)
    geometry!(model_,sample_)
    sample_
end

function iterateandstore!(runner_::GMHRunner,model_::AbstractModel,indicatorstate_::AbstractSamplerState,samplerstates_::Vector,indicator_::AbstractIndicatorMatrix,chain_::AbstractChain,store::Bool =true)
    outsidepriorcounter = 0
    while ~isfinite(proposals(iterate!(runner_,model_,indicatorstate_)).logprior[1])
        warn("Auxiliary proposal falls outside parameter priors. Retrying...")
        if (outsidepriorcounter+=1) > 100
            error("Difficulties generating auxiliary proposal, giving up. Please check the default and prior values of the parameters.")
        end
    end
    map!((s)->setfrom!(s,proposals(indicatorstate_)),samplerstates_)
    samplerstates_ = iterate!(runner_,model_,samplerstates_)
    transitionprobability!(indicator_,acceptanceratio(indicatorstate_),map(acceptanceratio,samplerstates_))
    sampleindicator!(indicator_)
    store?store!(runner_,indicatorstate_,samplerstates_,indicator_,chain_):nothing
    updatefrom!(runner_,indicatorstate_,samplerstates_,indicator_)
end

function iterateandstore!(runner_::GMHRunner,model_::AbstractModel,samplerstates_::Vector,indicator_::AbstractIndicatorMatrix,chain_::AbstractChain,store::Bool =true)
    iterate!(runner_,model_,samplerstates_[1])
    transitionprobability!(indicator_,acceptanceratio(samplerstates_[1]))
    sampleindicator!(indicator_)
    store?store!(runner_,samplerstates_[1],samplerstates_,indicator_,chain_):nothing
    updatefrom!(runner_,samplerstates_[1],samplerstates_,indicator_)
end

function iterate!(runner_::GMHRunner,model_::AbstractModel,samplerstates_::Vector)
    pmap((s)->iterate!(runner_,model_,s),samplerstates_)
end

function iterate!(runner_::GMHRunner,model_::AbstractModel,samplerstate_::AbstractSamplerState)
    propose!(samplerstate_)
    geometry!(model_,proposals(samplerstate_))
    acceptanceratio!(samplerstate_)
    samplerstate_
end

function store!(runner_::GMHRunner,indicatorstate_::AbstractSamplerState,samplerstates_::Vector,indicator_::AbstractIndicatorMatrix,chain_::AbstractChain)
    nprops = numproposals(indicator_)
    nsamp = numsamples(indicator_)
    indicatorsamples_ = indicatorsamples(indicator_)
    from_ = from(indicatorstate_)
    i2stuple = _ind2subtuple(samplerstates_)
    for i=2:nsamp+1
        if indicatorsamples_[i] == nprops + 1
            store!(chain_,from_,1)
        else
            s,t = ind2sub(i2stuple,indicatorsamples_[i])
            store!(chain_,proposals(samplerstates_[t]),s)
        end
    end
    accepted!(chain_,indicator_)
end

function updatefrom!(runner_::GMHRunner,indicatorstate_::AbstractSamplerState,samplerstates_::Vector,indicator_::AbstractIndicatorMatrix)
    indicatorend = indicatorsamples(indicator_)[end]
    i2stuple = _ind2subtuple(samplerstates_)
    if indicatorend != numproposals(indicator_) + 1
        s,t = ind2sub(i2stuple,indicatorend)
        setfrom!(indicatorstate_,proposals(samplerstates_[t]),s)
    end
    indicatorstate_
end

@inline _numsamplerstates(::Type{Val{:nprocs}}) = nprocs()
@inline _numsamplerstates(::Type{Val{:nworkers}}) = nworkers()
@inline _numsamplerstates(::Type{Val{:test}}) = 3

@inline _ind2subtuple(samplerstates_::Vector) = (numsamples(samplerstates_[1].proposals),length(samplerstates_))

function samplerstates(runner_::GMHRunner,sampler_::AbstractSampler)
    nsamplerstates = min(runner_.numproposals,_numsamplerstates(traittype(runner_.policy.samplerstates)))
    nproposalsperstate = ceil(Int,runner_.numproposals/nsamplerstates)
    AbstractSamplerState[samplerstate(sampler_,nproposalsperstate,runner_.policy.sampletype,runner_.policy.calculationtype) for i=1:nsamplerstates]
end

function show(io::IO,r::GMHRunner)
    println(io,"Generalized Metropolis-Hastings runner with:")
    println(io," numiterations: ",r.numiterations)
    println(io," numproposals: ",r.numproposals)
    println(io," numindicatorsamples: ",r.numindicatorsamples)
    println(io," numtotalsamples: ",r.numtotalsamples)
    println(io," numburnin: ",r.numburnin)
    print(io," policy: ")
    show(io,r.policy)
    println(io)
    nothing
end
