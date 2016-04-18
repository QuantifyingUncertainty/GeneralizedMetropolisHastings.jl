### Abstract parent types for Metropolis-Hastings samplers
abstract MHSampler <: AbstractSampler
abstract MHSamplerState <: AbstractSamplerState{BaseSample}

################################################################################################################################

### Type holding the parameters for a Normal-based Metrolopolis sampler
immutable MHNormal{T<:AbstractFloat,M<:AbstractArray} <: MHSampler
    initialscalefactor::T
    covariance::M
    MHNormal(s::T,c::AbstractArray{T}) = new(s,c)
end

###External constructors
@inline _sampler(::Type{Val{:mh}},::Type{Val{:normal}},scalefactor::AbstractFloat,m::AbstractArray) = MHNormal{typeof(scalefactor),typeof(m)}(scalefactor,m)
@inline _sampler(::Type{Val{:mh}},::Type{Val{:normal}},scalefactor::AbstractFloat,d::Integer) = _sampler(Val{:mh},Val{:normal},scalefactor,ones(typeof(scalefactor),d))

@inline numparas(sampler::MHNormal) = size(sampler.covariance,1)
@inline samplername(sampler::MHNormal) = "MHNormal"

################################################################################################################################

###Type holding the parameters for a LogNormal-based Metropolis-Hastings sampler
immutable MHLogNormal{T<:AbstractFloat,M<:AbstractArray} <: MHSampler
    initialscalefactor::T
    scale::M
    MHLogNormal(s::T,c::AbstractArray{T}) = new(s,c)
end

###External constructors
@inline _sampler(::Type{Val{:mh}},::Type{Val{:lognormal}},scalefactor::AbstractFloat,m::AbstractArray) = MHLogNormal{typeof(scalefactor),typeof(m)}(scalefactor,m)
@inline _sampler(::Type{Val{:mh}},::Type{Val{:lognormal}},scalefactor::AbstractFloat,d::Integer) = _sampler(Val{:mh},Val{:lognormal},scalefactor,ones(typeof(scalefactor),d))

@inline numparas(sampler::MHLogNormal) = size(sampler.scale,1)
@inline samplername(sampler::MHLogNormal) = "MHLogNormal"

################################################################################################################################

### Type holding the state of the Markov Chain for a Generalized Metropolis-Hastings sampler
type MHSymmetricSamplerState{T<:AbstractFloat,D<:SymmetricDensity} <: MHSamplerState
    density::D
    scalefactor::T
    from::BaseSample
    proposals::BaseSample
    acceptance::Vector
    MHSymmetricSamplerState(density::D,scalefactor::T,from::BaseSample{T},proposals::BaseSample{T},ac::Vector{T}) = new(density,scalefactor,from,proposals,ac)
end

type MHASymmetricSamplerState{T<:AbstractFloat,D<:ASymmetricDensity} <: MHSamplerState
    density::D
    scalefactor::T
    from::BaseSample
    proposals::BaseSample
    acceptance::Vector
    MHASymmetricSamplerState(density::D,scalefactor::T,from::BaseSample{T},proposals::BaseSample{T},ac::Vector{T}) = new(density,scalefactor,from,proposals,ac)
end

@inline samplerstatename(state::MHSymmetricSamplerState) = "MHSymmetricSamplerState"
@inline samplerstatename(state::MHASymmetricSamplerState) = "MHASymmetricSamplerState"

@inline _scale(scalefactor,densityparam::AbstractMatrix) = Base.scale(scalefactor,densityparam) #for covariance matrices of the proposal densities
@inline _scale(scalefactor,densityparam::AbstractVector) = sqrt(Base.scale(scalefactor,densityparam)) #for variance vectors of the proposal densities

@inline function _samplerstate{N<:Number}(sampler::MHNormal,from::BaseSample{N},proposals::BaseSample{N})
    d = _density(Val{:normal},zeros(N,numparas(from)),_scale(sampler.initialscalefactor,sampler.covariance))
    MHSymmetricSamplerState{N,typeof(d)}(d,sampler.initialscalefactor,from,proposals,zeros(N,numsamples(proposals)))
end

@inline function _samplerstate{N<:Number}(sampler::MHLogNormal,from::BaseSample{N},proposals::BaseSample{N})
    d = _density(Val{:lognormal},zeros(N,numparas(from)),_scale(sampler.initialscalefactor,sampler.scale))
    MHASymmetricSamplerState{N,typeof(d)}(d,sampler.initialscalefactor,from,proposals,zeros(N,numsamples(proposals)))
end

@inline function _samplerstate{N<:Number,T<:AbstractFloat}(sampler::MHSampler,nsamples::Integer,::Type{N},::Type{T})
    from = _samples(Val{:base},numparas(sampler),1,N,T)
    proposals = _samples(Val{:base},numparas(sampler),nsamples,N,T)
    _samplerstate(sampler,from,proposals)
end

@inline function setfrom!(state::MHSamplerState,from::BaseSample)
    copy!(state.from,from)
    state
end

@inline function setfrom!(state::MHSamplerState,from::BaseSample,i::Integer)
    copy!(state.from,1,from,i)
    state
end

###propose new samples from the current point
@inline function propose!(state::MHSamplerState)
    #set the current from point to sample from
    condition!(state.density,state.from.values)
    #propose n values from this point
    propose!(state.density,state.proposals)
    #return the proposals
    state
end

acceptanceratio(state::MHSamplerState) = state.acceptance

###for a symmetric proposal density, no adjustment of the acceptance rate is needed
function acceptanceratio!(state::MHSymmetricSamplerState)
    fl = state.from.loglikelihood[1] + state.from.logprior[1]
    for i=1:numsamples(state.proposals)
        @inbounds state.acceptance[i] = state.proposals.loglikelihood[i] + state.proposals.logprior[i] - fl
    end
    state.acceptance
end

###for an asymmetric proposal density, an adjustment needs to be performed
###based on the difference between K(from,val)/K(val,from)
function acceptanceratio!(state::MHASymmetricSamplerState)
    #condition the proposal density on the current proposal center
    condition!(state.density,state.from.values)
    #calculate the log-probability for each of the proposed points K(from,val)
    logprobability!(state.acceptance,state.density,state.proposals.values)
    #now calculate the inverse probability K(val,from)
    t = similar(state.acceptance,1)
    fl = state.from.loglikelihood[1] + state.from.logprior[1]
    for i = 1:numsamples(state.proposals)
        #condition the proposal density on each proposed point
        @inbounds condition!(state.density,state.proposals.values[:,i])
        #calculate the logprobability K(val,from)
        logprobability!(t,state.density,state.from.values)
        #adjust the result accordingly
        @inbounds state.acceptance[i] += (state.proposals.loglikelihood[i] + state.proposals.logprior[i] - t[1] - fl)
    end
    #restore the previous state of the proposal density
    condition!(state.density,state.from.values)
    #return the result
    state.acceptance
end

@inline _tune!(sampler::MHNormal,state::MHSymmetricSamplerState) = (condition!(state.density,state.from.values,_scale(state.scalefactor,sampler.covariance)) ; state)
@inline _tune!(sampler::MHLogNormal,state::MHASymmetricSamplerState) = (condition!(state.density,state.from.values,_scale(state.scalefactor,sampler.scale)) ; state)
@inline tune!(sampler::MHSampler,state::MHSamplerState,scalefactor::AbstractFloat) = (state.scalefactor *= scalefactor ; _tune!(sampler,state))
