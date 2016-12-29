### A segment of a Generalized Metropolis-Hastings run, to be executed on a single worker
immutable GMHSegment <: AbstractJobSegment
    model::AbstractModel
    samplerstate::AbstractSamplerState
end

### Factory function
function _segment(p::MHRuntimePolicy,m::AbstractModel,s::AbstractSampler,nproposals::Integer)
    @assert traitvalue(p.runner) == :generalized
    state = samplerstate(s,nproposals,p.sampletype,p.calculationtype,true)
    GMHSegment(m,state)
end

### Perform one iteration
function iterate!(seg::GMHSegment,indicator::AbstractSamplerState)
    prepareauxiliary!(indicator,seg.samplerstate)
    propose!(seg.samplerstate)
    geometry!(seg.model,proposals(seg.samplerstate))
    acceptance!(seg.samplerstate)
end

### Update the indicatorstate with the sample of this auxiliary state
prepare!(seg::GMHSegment,indicator::AbstractSamplerState,i::Int) = prepareindicator!(indicator,seg.samplerstate,i)

### Get the samples
getsamples(seg::GMHSegment,sampleindex) = copy(proposals(seg.samplerstate),sampleindex)

### Tune the samplerstate
tune!(seg::GMHSegment,tvals...) = tune!(seg.samplerstate,tvals...)

### Get the state variables of the encapsulated samplerstate
getsamplerstatevars(seg::GMHSegment) = getsamplerstatevars(seg.samplerstate)

###
function show(io::IO,s::GMHSegment)
    println(io,"GMHSegment with: ")
    println(io,"model:")
    show(io,s.model)
    println(io)
    println(io,"samplerstate:")
    show(io,s.samplerstate)
    println(io)
    nothing
end
