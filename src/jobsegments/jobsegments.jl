abstract AbstractJobSegment

### Factory function
segment(p::AbstractPolicy,args...;keyargs...) = _segment(p,args...;keyargs...)

numproposals(seg::AbstractJobSegment) = numsamples(proposals(seg.samplerstate))

### Interface
iterate!(seg::AbstractJobSegment,indicatorstate::AbstractSamplerState) = throw(MethodError(iterate!, (seg,indicatorstate)))
prepare!(seg::AbstractJobSegment,indicatorstate::AbstractSamplerState,i::Int) = throw(MethodError(prepare!, (seg,indicatorstate,i)))
tune!(seg::AbstractJobSegment,s::AbstractSampler,tvals...) = throw(MethodError(tune!, (seg,s.tvals...)))
getsamples(seg::AbstractJobSegment,sampleindex) = throw(MethodError(getsamples, (seg,sampleindex)))

### Interface with RemoteRefs, mapping to underlying AbstractJobSegments
iterate!(remote::RemoteRef,indicatorstate::AbstractSamplerState) = iterate!(fetch(remote),indicatorstate)
prepare!(remote::RemoteRef,indicatorstate::AbstractSamplerState,i::Int) = prepare!(fetch(remote),indicatorstate,i)
tune!(remote::RemoteRef,s::AbstractSampler,tvals...) = tune!(fetch(remote),s,tvals...)
getsamples(remote::RemoteRef,sampleindex) = getsamples(fetch(remote),sampleindex)



