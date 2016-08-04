abstract AbstractMHRunner <: AbstractRunner #an abstract super type for both Standard and Generalized Metropolis-Hastings

const maxinitialcounter = 100::Int

_runner(p::MHRuntimePolicy,args...;keyargs...) = _runner(traittype(p.runner),p,args...;keyargs...)

_numsamplesinchain(::Type{Val{:burnin}},nburnin::Int,niterations::Int,nindicatorsamples::Int) = nburnin*nindicatorsamples,nburnin*nindicatorsamples
_numsamplesinchain(::Type{Val{:main}},nburnin::Int,niterations::Int,nindicatorsamples::Int) = niterations*nindicatorsamples,0
_numsamplesinchain(::Type{Val{:all}},nburnin::Int,niterations::Int,nindicatorsamples::Int) = (nburnin+niterations)*nindicatorsamples+1,nburnin*nindicatorsamples+1

_storeduring(phase::Symbol,p::MHRuntimePolicy) = in(traitvalue(p.store),[phase,:all])

function createcommon(runner_::AbstractMHRunner,tuner_::AbstractTuner,nparas::Int,nproposals::Int,nindicatorsamples::Int)
    p = runner_.policy
    nchainsamples,nburninsamples = _numsamplesinchain(traittype(p.store),runner_.numburnin,runner_.numiterations,nindicatorsamples)
    i = indicator(traitvalue(p.indicator),nproposals,nindicatorsamples,p.calculationtype)
    a = tunerstate(tuner_,runner_.numburnin,p.calculationtype)
    c = chain(traitvalue(p.chain),nparas,nchainsamples,p.sampletype,p.calculationtype;numburnin=nburninsamples)
    i,a,c
end

function initialize!(runner_::AbstractMHRunner,model_::AbstractModel,state_::AbstractSamplerState,chain_::AbstractChain,store::Bool)
    initialcounter = 0
    sample_ = proposals(state_)
    #initialize the proposals field with values generated from the defaults or from the priors
    while ~isfinite((initialize!(runner_.policy.initialize,model_,sample_) ; geometry!(model_,sample_)).logprior[1])
        if (initialcounter+=1) > maxinitialcounter
            error("Problems generating initial proposal, giving up. Please check the parameter priors.")
        end
    end
    #prepare the distribution for the first proposal (this should result in the proposals field to be copied to the from field, plus other sampler-dependent setup)
    prepare!(state_,true)
    store?store!(chain_,sample_,1):nothing
    state_
end
