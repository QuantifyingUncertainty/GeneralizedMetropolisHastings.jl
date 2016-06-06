### Abstract parent types for Simplified Manifold MALA samplers
abstract AbstractSmMALASampler <: AbstractSampler
abstract AbstractSmMALASamplerState{S<:AbstractTensorSample} <: AbstractSamplerState{S}

abstract AbstractSmMALATrait <: AbstractPolicyTrait

immutable SmMALAType <: AbstractSmMALATrait
    trait::Symbol
    SmMALAType(s::Symbol) = (@assert in(s,[:standard,:tr]) ; new(s))
end

immutable TensorType <: AbstractSmMALATrait
    trait::Symbol
    TensorType(s::Symbol) = (@assert in(s,[:full,:tangent]) ; new(s))
end

@inline _traitname(::Type{SmMALAType}) = :smmala
@inline _traitname(::Type{TensorType}) = :tensor

@inline _bool2smmala(tr::Bool) = trait(:smmala,tr?(:tr):(:standard))

_trait(::Type{_traitnametype(SmMALAType)},t::Symbol) = SmMALAType(t)
_trait(::Type{_traitnametype(TensorType)},t::Symbol) = TensorType(t)

immutable SmMALAPolicy <: AbstractPolicy
    smmala::SmMALAType
    tensor::TensorType
    SmMALAPolicy(s,t) = new(s,t)
end

SmMALAPolicy(tr::Bool,t::Symbol) = SmMALAPolicy(_bool2smmala(tr),trait(:tensor,t))

################################################################################################################################

### Type holding the parameters for scalable SmMALA samplers
immutable SmMALATensorSampler{T<:AbstractFloat} <: AbstractSmMALASampler
    initialscalefactor::T
    nparas::Int
    tr::Bool #trust region sampler?
    SmMALATensorSampler(s::T,np::Int,tr::Bool) = new(s,np,tr)
end

immutable SmMALATangentSampler{T<:AbstractFloat} <: AbstractSmMALASampler
    initialscalefactor::T
    nparas::Int
    ntangents::Int
    SmMALATangentSampler(s::T,np::Int,nt::Int) = new(s,np,nt)
end

###Factory functions
@inline _sampler(::Type{Val{:smmala}},i::AbstractFloat,np::Integer,tr::Bool) = SmMALATensorSampler{typeof(i)}(i,np,tr)
@inline _sampler(::Type{Val{:smmala}},i::AbstractFloat,np::Integer,nt::Integer) = SmMALATangentSampler{typeof(i)}(i,np,nt)

###Size functions
@inline numparas(sampler::AbstractSmMALASampler) = sampler.nparas

###Name functions (used in println for generic AbstractSampler)
@inline samplername(sampler::SmMALATensorSampler) = string(sampler.tr?"Trust Region ":"","Simplified mMALA")
@inline samplername(sampler::SmMALATangentSampler) = "Trust Region Simplified mMALA with Tangent Tensor"

################################################################################################################################

### Type holding the state of the Markov Chain for a Generalized SmMALA sampler
type SmMALASamplerState{T<:AbstractFloat,S<:AbstractTensorSample} <: AbstractSmMALASamplerState{S}
    policy::SmMALAPolicy
    density::NormalDensity
    from::S
    proposals::S
    acceptance::Vector
    scale::T
    SmMALASamplerState(pol::SmMALAPolicy,d::NormalDensity,f::S,p::S,a::Vector,s::T) = new(pol,d,f,p,a,s)
end

### Helper factory functions to create the samples for different SmMALA samplers
@inline function _samples{N<:Number,T<:AbstractFloat}(s::SmMALATensorSampler,nsamples::Int,::Type{N},::Type{T})
    SmMALAPolicy(s.tr,:full),_samples(Val{:tensor},s.nparas,1,N,T),_samples(Val{:tensor},s.nparas,nsamples,N,T)
end

@inline function _samples{N<:Number,T<:AbstractFloat}(s::SmMALATangentSampler,nsamples::Int,::Type{N},::Type{T})
    SmMALAPolicy(true,:tangent),_samples(Val{:tangent},s.nparas,1,N,T,s.ntangents),_samples(Val{:tangent},s.nparas,nsamples,N,T,s.ntangents)
end

### Factory functions
@inline function _samplerstate{N<:Number,T<:AbstractFloat}(s::AbstractSmMALASampler,nsamples::Integer,::Type{N},::Type{T},::Bool)
    pol,f,p = _samples(s,nsamples,N,T)
    d = _density(Val{:normal},zeros(N,s.nparas),eye(N,s.nparas))
    SmMALASamplerState{N,typeof(f)}(pol,d,f,p,zeros(T,nsamples),s.initialscalefactor)
end

###Name functions (used in println for generic AbstractSampler)
@inline _samplerstatename(::Type{Val{:standard}},::Type{Val{:full}}) = "Simplified mMALA"
@inline _samplerstatename(::Type{Val{:tr}},::Type{Val{:full}}) = "Trust Region Simplified mMALA"
@inline _samplerstatename(::Type{Val{:tr}},::Type{Val{:tangent}}) = "Trust Region Simplified mMALA with Tangent Tensor"
@inline samplerstatename(sampler::SmMALASamplerState) = _samplerstatename(traittype(sampler.policy.smmala),traittype(sampler.policy.tensor))


###update the sampler state, called from SMHRunner
@inline function prepare!(state::SmMALASamplerState,updatefrom::Bool=false)
    if updatefrom
        #set the current point around which we are sampling
        setfrom!(state,state.proposals)
        #update the mean and covariance for the proposal distribution
        _updatedensity!(state)
    end
    #return the state
    state
end

###update the auxiliary samplerstate with information from the indicator samplerstate
@inline function prepareauxiliary!(indicatorstate::SmMALASamplerState,auxiliarystate::SmMALASamplerState)
    #set the auxiliary point as the from point in the auxiliarystate
    setfrom!(auxiliarystate,indicatorstate.proposals)
    #update the mean and covariance for the proposal distribution
    _updatedensity!(auxiliarystate)
    #return the auxiliary state
    auxiliarystate
end

###when the "from" field does not need updating
@inline prepareindicator!(indicatorstate::SmMALASamplerState) = indicatorstate

###when the "from" field has been updated
@inline function prepareindicator!(indicatorstate::SmMALASamplerState,auxiliarystate::SmMALASamplerState,i::Int)
    #set the current point around which we are sampling
    setfrom!(indicatorstate,auxiliarystate.proposals,i)
    #update the mean and covariance for the proposal distribution
    _updatedensity!(indicatorstate)
    #return the indicator state
    indicatorstate
end

###propose new samples from the current point
@inline function propose!(state::SmMALASamplerState)
    #propose n values from the current point
    propose!(state.density,state.proposals.values)
    #return the proposals
    state
end

###for an asymmetric proposal density, an adjustment needs to be performed
###based on the difference between K(from,val)/K(val,from)
function acceptance!(state::SmMALASamplerState)
    #calculate the log-probability for each of the proposed points K(from,val)
    logprobability!(state.acceptance,state.density,state.proposals.values)
    #now calculate the inverse probability K(val,from)
    temp1 = similar(state.acceptance,1)
    fl = state.from.loglikelihood[1] + state.from.logprior[1]
    for i = 1:numsamples(state.proposals)
        _updatedensity!(state,i)
        #calculate the logprobability K(val,from)
        logprobability!(temp1,state.density,state.from.values)
        #adjust the result accordingly
        @inbounds state.acceptance[i] += (state.proposals.loglikelihood[i] + state.proposals.logprior[i] - temp1[1] - fl)
    end
    #restore the density to its original state
    _updatedensity!(state)
    #return the result
    state.acceptance
end

#############LOCAL FUNCTIONS##########################

### Helper function for trust region calculations of mean and covariance
@inline _adddiag!{T<:AbstractFloat}(t::AbstractMatrix{T},v::T) = (@simd for i=1:size(t,1) @inbounds t[i,i]+=v end ; t)
@inline _addfrom!(m::Vector,f::Vector) = (@simd for i=1:length(m) @inbounds m[i]+=f[i] end ; m)

### Calculate the SmMALA mean and covariance
@inline function _meancov{T<:AbstractFloat}(::Type{Val{:standard}},f::Vector{T},g::AbstractVector{T},t::AbstractMatrix{T},s::T)
    _addfrom!(scale!(t\g,s*s/2),f),scale!(inv(t),s*s)
end

### Calculate the Trust Region SmMALA mean and covariance
@inline function _meancov{T<:AbstractFloat}(::Type{Val{:tr}},f::Vector{T},g::AbstractVector{T},t::AbstractMatrix{T},s::T)
    _adddiag!(t,1/s)
    m = _addfrom!(t\g,f)
    c = inv(t)
    _adddiag!(t,-1/s)
    m,c
end

@inline function _meancov(s::SmMALASamplerState)
    _meancov(traittype(s.policy.smmala),s.from.values,s.from.gradloglikelihood,s.from.tensorloglikelihood,s.scale)
end

@inline function _meancov(s::SmMALASamplerState,i::Int)
    _meancov(traittype(s.policy.smmala),s.proposals.values[:,i],s.proposals.gradloglikelihood[:,i],s.proposals.tensorloglikelihood[:,:,i],s.scale)
end

@inline _updatedensity!(s::SmMALASamplerState) = update!(s.density,_meancov(s)...)
@inline _updatedensity!(s::SmMALASamplerState,i::Int) = update!(s.density,_meancov(s,i)...)

