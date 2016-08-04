abstract AbstractJobSegment

### Factory function
segment(p::AbstractPolicy,args...;keyargs...) = _segment(p,args...;keyargs...)

numproposals(seg::AbstractJobSegment) = numsamples(proposals(seg.samplerstate))

### Interface for AbstractJobSegments
iterate!(seg::AbstractJobSegment,indicatorstate::AbstractSamplerState) = throw(MethodError(iterate!, (seg,indicatorstate)))
prepare!(seg::AbstractJobSegment,indicatorstate::AbstractSamplerState,i::Int) = throw(MethodError(prepare!, (seg,indicatorstate,i)))
tune!(seg::AbstractJobSegment,tvals...) = throw(MethodError(tune!, (seg,tvals...)))
getsamples(seg::AbstractJobSegment,sampleindex) = throw(MethodError(getsamples, (seg,sampleindex)))

### Interface which maps from RemoteRefs to AbstractJobSegments
iterate!(remote::RemoteRef,indicatorstate::AbstractSamplerState) = iterate!(fetch(remote),indicatorstate)
prepare!(remote::RemoteRef,indicatorstate::AbstractSamplerState,i::Int) = prepare!(fetch(remote),indicatorstate,i)
tune!(remote::RemoteRef,tvals...) = tune!(fetch(remote),tvals...)
getsamples(remote::RemoteRef,sampleindex) = getsamples(fetch(remote),sampleindex)



