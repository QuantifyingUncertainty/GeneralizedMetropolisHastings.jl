abstract AbstractJobSegment

### Factory function
segment(p::AbstractPolicy,args...;keyargs...) = _segment(p,args...;keyargs...)

numproposals(seg::AbstractJobSegment) = numsamples(proposals(seg.samplerstate))

### Interface for AbstractJobSegments
iterate!(seg::AbstractJobSegment,indicatorstate::AbstractSamplerState) = throw(MethodError(iterate!, (seg,indicatorstate)))
prepare!(seg::AbstractJobSegment,indicatorstate::AbstractSamplerState,i::Int) = throw(MethodError(prepare!, (seg,indicatorstate,i)))
tune!(seg::AbstractJobSegment,tvals...) = throw(MethodError(tune!, (seg,tvals...)))
getsamples(seg::AbstractJobSegment,sampleindex) = throw(MethodError(getsamples, (seg,sampleindex)))
getsamplerstatevars(seg::AbstractJobSegment) = throw(MethodError(getsamplerstatevars, (seg)))

### Interface which maps from Future remote references to AbstractJobSegments
@compat function iterate!(remote::Future,indicatorstate::AbstractSamplerState)
    r = fetch(remote)
    if ~isa(r,RemoteException)
        iterate!(r,indicatorstate)
    else
        throw(r)
    end
end

@compat function prepare!(remote::Future,indicatorstate::AbstractSamplerState,i::Int)
    r = fetch(remote)
    if ~isa(r,RemoteException)
        prepare!(r,indicatorstate,i)
    else
        throw(r)
    end
end

@compat function tune!(remote::Future,tvals...)
    r = fetch(remote)
    if ~isa(r,RemoteException)
        tune!(r,tvals...)
    else
        throw(r)
    end
end

@compat function getsamples(remote::Future,sampleindex)
    r = fetch(remote)
    if ~isa(r,RemoteException)
        getsamples(r,sampleindex)
    else
        throw(r)
    end
end

@compat function getsamplerstatevars(remote::Future)
    r = fetch(remote)
    if ~isa(r,RemoteException)
        getsamplerstatevars(r)
    else
        throw(r)
    end
end
