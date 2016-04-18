### A segment of a Generalized Metropolis-Hastings run, to be executed on a single worker
immutable GMHSegment <: AbstractJobSegment
    model::AbstractModel
    samplerstate::AbstractSamplerState
end

### Factory function
function _segment(p::MHRuntimePolicy,m::AbstractModel,s::AbstractSampler,nproposals::Integer)
    @assert traitvalue(p.runner) == :generalized
    state = samplerstate(s,nproposals,p.sampletype,p.calculationtype)
    GMHSegment(m,state)
end

### Perform one iteration
function iterate!(seg::GMHSegment,f::AbstractSample)
    setfrom!(seg.samplerstate,f)
    propose!(seg.samplerstate)
    geometry!(seg.model,proposals(seg.samplerstate))
    a = acceptanceratio!(seg.samplerstate)
    a
end

### Get the samples
getsamples(seg::GMHSegment,sampleindex) = copy(proposals(seg.samplerstate),sampleindex)

### Tune the samplerstate
tune!(seg::GMHSegment,s::AbstractSampler,tvals...) = tune!(s,seg.samplerstate,tvals...)

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
