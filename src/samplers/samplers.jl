### Sampler types hold the components that fully specify a Monte Carlo sampler
abstract AbstractSampler

### SamplerState types hold the temporary components used by a Monte Carlo sampler during its run
### This means that SamplerState types represent the internal state ("local variables") of a Monte Carlo sampler
abstract AbstractSamplerState{S<:AbstractSample}

### Factory functions for samplers
### Currently implemented samplers: :mh (Metropolis-Hastings)
### Currently implemented proposal densities: :normal, :lognormal
sampler(name::Symbol,density::Symbol,args...) = _sampler(Val{name},Val{density},args...)

### Factory functions for SamplerStates
samplerstate{N<:Number,T<:AbstractFloat}(s::AbstractSampler,nsamples::Integer,::Type{N},::Type{T}) = _samplerstate(s,nsamples,N,T)

###Get the from field of the samplerstate
from(s::AbstractSamplerState) = s.from
proposals(s::AbstractSamplerState) = s.proposals

###Functions that need to be implemented for specific samplers
setfrom!(state::AbstractSamplerState,sample::AbstractSample) = throw(MethodError(setfrom!, (state,sample)))
setfrom!(state::AbstractSamplerState,sample::AbstractSample,i::Integer) = throw(MethodError(setfrom!, (state,sample,i)))
propose!(state::AbstractSamplerState) = throw(MethodError(propose!, (sampler,state)))
acceptanceratio!(state::AbstractSamplerState) = throw(MethodError(acceptanceratio!, (state)))
acceptanceratio(state::AbstractSamplerState) = throw(MethodError(acceptanceratio, (state)))
tune!(sampler::AbstractSampler,state::AbstractSamplerState,args...) = nothing

###Generic show function for any sampler
function show(io::IO,s::AbstractSampler)
    println(io,"AbstractSampler $(samplername(s)) with fields:")
    println(io,"  $(fieldnames(s))")
    nothing
end

###Generic show for heaps
function show(io::IO,s::AbstractSamplerState)
    println(io,"AbstractSamplerState $(samplerstatename(s)) with fields: ")
    println(io,"  $(fieldnames(s))")
    nothing
end

