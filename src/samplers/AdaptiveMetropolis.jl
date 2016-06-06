### Abstract parent types for Adaptive samplers
abstract AbstractAdaptiveSampler <: AbstractSampler
abstract AbstractAdaptiveSamplerState <: AbstractSamplerState{BaseSample}

################################################################################################################################

### Type holding the parameters for an adaptive Metropolis sampler (i.e with a normal density the covariance of which is adapted to the landscape)
immutable AdaptiveNormalSampler{T<:AbstractFloat} <: AbstractAdaptiveSampler
    initialscalefactor::T
    nparas::Int
    AdaptiveNormalSampler(s::T,n::Int) = new(s,n)
end

###Factory functions
@inline _sampler(::Type{Val{:adaptive}},s::AbstractFloat,n::Integer) = AdaptiveNormalSampler{typeof(s)}(s,n)

###Size functions
@inline numparas(sampler::AdaptiveNormalSampler) = sampler.nparas

###Name functions (used in println for generic AbstractSampler)
@inline samplername(sampler::AdaptiveNormalSampler) = "Adaptive Normal Metropolis"

################################################################################################################################

### Type holding the running mean and covariance for the adaptive sampler with normal density
type NormalRunningState{T<:AbstractFloat}
    iteration::Int
    runningmean::Vector{T}
    runningcov::Matrix{T}
    scale::T
    tempdiff::Vector{T}
    densitycov::Matrix{T}
    NormalRunningState(i::Int,m::Vector{T},c::Matrix{T},s::T) = new(i,m,c,s,similar(m),similar(c))
end

### Factory function
@inline function _runningstate{T<:AbstractFloat}(::Type{Val{:normal}},n::Int,s::T)
    NormalRunningState{T}(0,zeros(T,n),zeros(T,n,n),s)
end

@inline _covupdatefactor{T<:AbstractFloat}(i::Int,::Type{T}) = i==1?zero(T):T((i-2)/(i-1))

@inline function _updatetempdiff!{T<:AbstractFloat}(d::Vector{T},m::Vector{T},v::Vector{T})
    @simd for i=1:length(v)
        @inbounds d[i] = v[i] - m[i]
    end
    d
end

@inline function _updaterunningcov!{T<:AbstractFloat}(c::Matrix{T},d::Vector{T},ci::Int)
    f = _covupdatefactor(ci,T)
    nrows = length(d)
    for j=1:size(c,2)
        @simd for i=1:nrows
            @inbounds c[i,j] = f*c[i,j] + d[i]*d[j]/ci
        end
    end
    c
end

@inline function _updaterunningmean!{T<:AbstractFloat}(m::Vector{T},v::Vector{T},ci::Int)
    f::T = T(ci-1)
    @simd for i=1:length(m)
        @inbounds m[i] = (f*m[i] + v[i])/ci
    end
    m
end

@inline function _updatedensitycov!{T<:AbstractFloat}(d::Matrix{T},c::Matrix{T},s::T)
    f = s*s
    copy!(d,c)
    @simd for i=1:size(d,1)
        @inbounds d[i,i] += f
    end
    d
end

### Update running mean and covariance
### The formula for running covariance can be derived from the recurrence formulas
### for running covariance, derived by Knuth (see www.johndcook.com/blog/standard_deviation)
@inline function updaterunning!(r::NormalRunningState,v::Vector)
    r.iteration += 1
    _updatetempdiff!(r.tempdiff,r.runningmean,v)
    _updaterunningcov!(r.runningcov,r.tempdiff,r.iteration)
    _updaterunningmean!(r.runningmean,v,r.iteration)
    _updatedensitycov!(r.densitycov,r.runningcov,r.scale)
    r
end

### Type holding the state of the Markov Chain for an Adaptive Metropolis Sampler (with a normal density)
type AdaptiveNormalSamplerState{T<:AbstractFloat} <: AbstractAdaptiveSamplerState
    density::NormalDensity
    from::BaseSample
    proposals::BaseSample
    acceptance::Vector
    #element to keep track of running mean and covariance
    #this is only needed for the main state, not the auxiliary states, hence wrapping it into a "nullable"
    runningstate::Nullable{NormalRunningState{T}}
    AdaptiveNormalSamplerState(n::NormalDensity,f::BaseSample,p::BaseSample,a::Vector,r::Nullable{NormalRunningState{T}}) = new(n,f,p,a,r)
end

###Factory functions
@inline function _samplerstate{N<:Number,T<:AbstractFloat}(s::AdaptiveNormalSampler,nsamples::Integer,::Type{N},::Type{T},auxiliary::Bool)
    nparas = numparas(s)
    f = _samples(Val{:base},nparas,1,N,T)
    p = _samples(Val{:base},nparas,nsamples,N,T)
    n = _density(Val{:normal},s.initialscalefactor*s.initialscalefactor*eye(N,nparas))
    r = auxiliary?Nullable{NormalRunningState{N}}():Nullable(_runningstate(Val{:normal},nparas,s.initialscalefactor))
    AdaptiveNormalSamplerState{N}(n,f,p,zeros(T,nsamples),r)
end

###Name functions (used in println for AbstractSamplerState)
@inline samplerstatename(state::AdaptiveNormalSamplerState) = "Adaptive Normal Metropolis"

@inline function updaterunning!(state::AdaptiveNormalSamplerState)
    if !isnull(state.runningstate)
        r = get(state.runningstate)
        updaterunning!(r,state.from.values)
        state.density = update!(state.density,r.densitycov)
    end
    state
end

###update the sampler state, called from SMHRunner
@inline function prepare!(state::AdaptiveNormalSamplerState,updatefrom::Bool=false)
    if updatefrom
        #set the current point around which we are sampling
        setfrom!(state,state.proposals)
    end
    updaterunning!(state)
end

###update the auxiliary samplerstate with information from the indicator samplerstate
@inline function prepareauxiliary!(indicatorstate::AdaptiveNormalSamplerState,auxiliarystate::AdaptiveNormalSamplerState)
    #set the auxiliary point
    setfrom!(auxiliarystate,indicatorstate.proposals)
    #set the proposal density of the auxiliary state to that of the indicatorstate
    auxiliarystate.density = indicatorstate.density
    #return the auxiliary state
    auxiliarystate
end

@inline prepareindicator!(indicatorstate::AdaptiveNormalSamplerState) = updaterunning!(indicatorstate)

###update the indicator state when starting an iteration of GMHRunner
@inline function prepareindicator!(indicatorstate::AdaptiveNormalSamplerState,auxiliarystate::AdaptiveNormalSamplerState,i::Int)
    #set the current point around which we are sampling
    setfrom!(indicatorstate,auxiliarystate.proposals,i)
    #update the running mean and coveriance
    updaterunning!(indicatorstate)
end

###propose new samples from the current point
@inline function propose!(state::AdaptiveNormalSamplerState)
    #propose n values using the covariance around mean 0.0
    propose!(state.density,state.proposals.values)
    #add the current from to the proposals
    offset!(state.proposals,state.from.values)
    #return the proposals
    state
end

###for a symmetric proposal density, no adjustment of the acceptance rate is needed
function acceptance!(state::AdaptiveNormalSamplerState)
    fl = state.from.loglikelihood[1] + state.from.logprior[1]
    for i=1:numsamples(state.proposals)
        @inbounds state.acceptance[i] = state.proposals.loglikelihood[i] + state.proposals.logprior[i] - fl
    end
    state.acceptance
end

@inline function tune!(state::AdaptiveNormalSamplerState,scalefactor::AbstractFloat)
    if !isnull(state.runningstate) #only the main state needs updating
        get(state.runningstate).scale *= scalefactor
    end
    state
end


