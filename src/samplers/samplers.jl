###Get the number of parameters
numparas(s::MCSampler) = s.nparas

###Get the number of elements in the heap
numel(h::MCHeap) = length(h.samples)

###Vectorized versions of the in-place proposal function for arbitrary samplers and heaps
propose!(s::MCSampler,h::MCHeap,indicator::Int,from::MCSample) = (set_from!(s,h,from); @simd for j=1:numel(h) @inbounds j!=indicator?propose!(s,h,h.samples[j]):nothing end)
propose!(s::MCSampler,h::MCHeap,indicator::Int) = propose!(s,h,indicator,h.samples[indicator]) #propose from the current iterator

###Vectorized versions of the in-place update function for arbitrary samplers and heaps
update_proposals!(s::MCSampler,h::MCHeap,indicator::Int) = @simd for j=1:numel(h) @inbounds j!=indicator?update_proposal!(s,h,j):nothing end

###Generic show function for any sampler
Base.show(io::IO,s::MCSampler) = dump(io,s)

###Generic show for heaps
function Base.show(io::IO,h::MCHeap)
  println(io,typeof(h)," with ",length(h.samples)," samples")
  show(io,h.samples)
end
