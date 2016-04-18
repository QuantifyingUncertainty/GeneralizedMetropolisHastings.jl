abstract AbstractJobSegment

### Factory function
segment(p::AbstractPolicy,args...;keyargs...) = _segment(p,args...;keyargs...)

numproposals(seg::AbstractJobSegment) = numsamples(proposals(seg.samplerstate))

### Interface
iterate!(seg::AbstractJobSegment,from::AbstractSample) = throw(MethodError(iterate!, (seg,from)))
tune!(seg::AbstractJobSegment,s::AbstractSampler,tvals...) = throw(MethodError(tune!, (seg,s.tvals...)))
getsamples(seg::AbstractJobSegment,sampleindex) = throw(MethodError(getsamples, (seg,sampleindex)))

### Interface with RemoteRefs, mapping to underlying AbstractJobSegments
iterate!(remote::RemoteRef,from::AbstractSample) = iterate!(fetch(remote),from)
tune!(remote::RemoteRef,s::AbstractSampler,tvals...) = tune!(fetch(remote),s,tvals...)
getsamples(remote::RemoteRef,sampleindex) = getsamples(fetch(remote),sampleindex)



