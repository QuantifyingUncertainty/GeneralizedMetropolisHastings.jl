### Abstract parent types for Metropolis-Hastings samplers
abstract AbstractMetropolisHastingsSampler <: AbstractSampler
abstract AbstractMetropolisHastingsSamplerState <: AbstractSamplerState{BaseSample}

################################################################################################################################

### Type holding the parameters for a scalable Metropolis-Hastings sampler (i.e., any proposal density with a scale parameter)
immutable ScalableMetropolisHastingsSampler{T<:AbstractFloat,M<:AbstractArray} <: AbstractMetropolisHastingsSampler
    density::Symbol
    initialscalefactor::T
    scaleparameters::M
    extraargs::Tuple
    ScalableMetropolisHastingsSampler(d::Symbol,i::T,p::AbstractArray{T},args...) = new(d,i,p,args)
end

###Factory functions
@inline _sampler(::Type{Val{:mh}},den::Symbol,i::AbstractFloat,m::AbstractArray,args...) = ScalableMetropolisHastingsSampler{typeof(i),typeof(m)}(den,i,m,args...)
@inline _sampler(::Type{Val{:mh}},den::Symbol,i::AbstractFloat,d::Integer,args...) = _sampler(Val{:mh},den,i,ones(typeof(i),d),args...)

###Size functions
@inline numparas(sampler::ScalableMetropolisHastingsSampler) = size(sampler.scaleparameters,1)

###Name functions (used in println for generic AbstractSampler)
@inline samplername(sampler::ScalableMetropolisHastingsSampler) = "Scalable Metropolis-Hastings"

################################################################################################################################

### Type holding the state of the Markov Chain for a Generalized Metropolis-Hastings sampler
type MetropolisHastingsSamplerState{D<:AbstractProposalDensity} <: AbstractMetropolisHastingsSamplerState
    density::D
    from::BaseSample
    proposals::BaseSample
    acceptance::Vector
    MetropolisHastingsSamplerState(d::D,f::BaseSample,p::BaseSample,a::Vector) = new(d,f,p,a)
end

###Factory functions
@inline function _samplerstate{N<:Number,T<:AbstractFloat}(s::ScalableMetropolisHastingsSampler,nsamples::Integer,::Type{N},::Type{T},::Bool)
    nparas = numparas(s)
    f = _samples(Val{:base},nparas,1,N,T)
    p = _samples(Val{:base},nparas,nsamples,N,T)
    d = _density(Val{s.density},zeros(N,nparas),deepcopy(s.scaleparameters),deepcopy(s.extraargs)...)
    s.initialscalefactor!=one(T)?scale!(d,s.initialscalefactor):nothing #rescale the distribution for a scale factor different from 1.0
    MetropolisHastingsSamplerState{typeof(d)}(d,f,p,zeros(T,nsamples))
end

###Name functions (used in println for AbstractSamplerState)
@inline samplerstatename{D<:SymmetricDensity}(state::MetropolisHastingsSamplerState{D}) = "Metropolis"
@inline samplerstatename{D<:ASymmetricDensity}(state::MetropolisHastingsSamplerState{D}) = "Metropolis-Hastings"

###update the sampler state, called from SMHRunner
@inline function prepare!(state::MetropolisHastingsSamplerState,updatefrom::Bool=false)
    if updatefrom
        #set the current point around which we are sampling
        setfrom!(state,state.proposals)
        #condition the density on the current point
        condition!(state.density,state.from.values)
    end
    #return the state
    state
end

###update the auxiliary samplerstate with information from the indicator samplerstate
@inline function prepareauxiliary!(indicatorstate::MetropolisHastingsSamplerState,auxiliarystate::MetropolisHastingsSamplerState)
    #set the auxiliary point as the from point in the auxiliarystate
    setfrom!(auxiliarystate,indicatorstate.proposals)
    #condition the density on the current point
    condition!(auxiliarystate.density,auxiliarystate.from.values)
    #return the auxiliary state
    auxiliarystate
end

###when the "from" field does not need updating
@inline prepareindicator!(indicatorstate::MetropolisHastingsSamplerState) = indicatorstate

###when the "from" field has been updated
@inline function prepareindicator!(indicatorstate::MetropolisHastingsSamplerState,auxiliarystate::MetropolisHastingsSamplerState,i::Int)
    #set the current point around which we are sampling
    setfrom!(indicatorstate,auxiliarystate.proposals,i)
    #condition the density on the current point
    condition!(indicatorstate.density,indicatorstate.from.values)
    #return the indicator state
    indicatorstate
end

###propose new samples from the current point
@inline function propose!(state::MetropolisHastingsSamplerState)
    #propose n values from the current point
    propose!(state.density,state.proposals.values)
    #return the proposals
    state
end

###for a symmetric proposal density, no adjustment of the acceptance rate is needed
function acceptance!{D<:SymmetricDensity}(state::MetropolisHastingsSamplerState{D})
    fl = state.from.loglikelihood[1] + state.from.logprior[1]
    for i=1:numsamples(state.proposals)
        @inbounds state.acceptance[i] = state.proposals.loglikelihood[i] + state.proposals.logprior[i] - fl
    end
    state.acceptance
end

###for an asymmetric proposal density, an adjustment needs to be performed
###based on the difference between K(from,val)/K(val,from)
function acceptance!{D<:ASymmetricDensity}(state::MetropolisHastingsSamplerState{D})
    #calculate the log-probability for each of the proposed points K(from,val)
    logprobability!(state.acceptance,state.density,state.proposals.values)
    #now calculate the inverse probability K(val,from)
    temp1 = similar(state.acceptance,1)
    fl = state.from.loglikelihood[1] + state.from.logprior[1]
    for i = 1:numsamples(state.proposals)
        #condition the proposal density on each proposed point
        @inbounds condition!(state.density,state.proposals.values[:,i])
        #calculate the logprobability K(val,from)
        logprobability!(temp1,state.density,state.from.values)
        #adjust the result accordingly
        @inbounds state.acceptance[i] += (state.proposals.loglikelihood[i] + state.proposals.logprior[i] - temp1[1] - fl)
    end
    #restore the previous state of the proposal density
    condition!(state.density,state.from.values)
    #return the result
    state.acceptance
end

@inline function tune!(state::MetropolisHastingsSamplerState,scalefactor::AbstractFloat)
    scale!(state.density,scalefactor)
    state
end
