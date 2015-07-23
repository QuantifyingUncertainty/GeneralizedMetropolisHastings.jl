###Vectorized versions of the in-place proposal function for arbitrary samplers and heaps
propose!(s::MCSampler,h::MCHeap,indicator::Int,from::MCSample) = (set_from!(s,h,from); @simd for j=1:length(h.samples) @inbounds j!=indicator?propose!(s,h,h.samples[j]):nothing end)
propose!(s::MCSampler,h::MCHeap,indicator::Int) = propose!(s,h,indicator,h.samples[indicator]) #propose from the current iterator

###Vectorized versions of the in-place update function for arbitrary samplers and heaps
update_proposals!(s::MCSampler,h::MCHeap,indicator::Int) = @simd for j=1:length(h.samples) @inbounds j!=indicator?update_proposal!(s,h,j):nothing end

###Per-sample update function for arbitrary samplers and heaps is empty function
update_proposal!(::MCSampler,::MCHeap,::Int) = ()

###Vectorized version of show
function Base.show(io::IO,h::MCHeap)
  println(io,typeof(h)," with ",length(h.samples)," samples")
  for i=1:length(h.samples)
    show(io,h.samples[i])
  end
end
