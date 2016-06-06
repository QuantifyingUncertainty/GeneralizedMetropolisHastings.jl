### Sampler types hold the parameters and components that fully specify a Monte Carlo sampler
abstract AbstractSampler

### SamplerState types hold the temporary components used by a Monte Carlo sampler during its run
### This means that SamplerState types represent the internal state ("local variables") of a Monte Carlo sampler
abstract AbstractSamplerState{S<:AbstractSample}

### Factory functions for samplers
### Currently implemented samplers: :mh (Metropolis-Hastings), :adaptive (Adaptive Metropolis)
sampler(name::Symbol,args...;keyargs...) = _sampler(Val{name},args...;keyargs...)

### Factory functions for SamplerStates
samplerstate{N<:Number,T<:AbstractFloat}(s::AbstractSampler,nsamples::Integer,::Type{N},::Type{T},auxiliary::Bool=false) = _samplerstate(s,nsamples,N,T,auxiliary)

@inline samplername(s::AbstractSampler) = "Unknown Sampler"
@inline samplerstatename(s::AbstractSamplerState) = "Unknown Sampler State"

### Common access functions defined for all samplers
@inline numparas(s::AbstractSamplerState) = numparas(s.proposals)
@inline numsamples(s::AbstractSamplerState) = numsamples(s.proposals)
@inline from(s::AbstractSamplerState) = s.from
@inline proposals(s::AbstractSamplerState) = s.proposals
@inline acceptance(s::AbstractSamplerState) = s.acceptance
@inline setfrom!(state::AbstractSamplerState,sample::AbstractSample) = (copy!(state.from,sample) ; state)
@inline setfrom!(state::AbstractSamplerState,sample::AbstractSample,i::Integer) = (copy!(state.from,1,sample,i) ; state)

###Functions that need to be implemented for specific samplers
#This function should be implemented for normal Metropolis-Hastings
prepare!(samplerstate::AbstractSamplerState,updatefrom::Bool =false) = throw(MethodError(prepare!, (samplerstate,updatefrom)))
#These functions need to be implemented for generalized Metropolis-Hastings
#prepare the auxiliary state (called just before propose! is called on the auxiliary state)
prepareauxiliary!(indicator::AbstractSamplerState,auxiliary::AbstractSamplerState) = throw(MethodError(prepareauxiliary!, (indicator,auxiliary)))
#prepare the indicator state
prepareindicator!(indicator::AbstractSamplerState) = throw(MethodError(prepareindicator!, (indicatorstate,)))
prepareindicator!(indicator::AbstractSamplerState,auxiliary::AbstractSamplerState,i::Int) = throw(MethodError(prepareindicator!, (indicator,auxiliary,i)))

propose!(state::AbstractSamplerState) = throw(MethodError(propose!, (sampler,state)))
acceptance!(state::AbstractSamplerState) = throw(MethodError(acceptanceratio!, (state)))
tune!(state::AbstractSamplerState,args...) = nothing

###Generic show function for any sampler
function show(io::IO,s::AbstractSampler)
    println(io,"$(samplername(s)) Sampler with fields:")
    println(io,"  $(fieldnames(s))")
    nothing
end

###Generic show for samepler states
function show(io::IO,s::AbstractSamplerState)
    println(io,"$(samplerstatename(s)) SamplerState with fields: ")
    println(io,"  $(fieldnames(s))")
    nothing
end

