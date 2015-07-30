###Vectorized versions of the in-place proposal function for arbitrary densities and samples
propose!{D<:ProposalDensity,S<:MCSample}(d::D,v::Vector{S}) = @simd for j=1:length(v) @inbounds propose!(d,v[j]) end
propose!{D<:ProposalDensity,S<:MCSample}(d::Vector{D},v::Vector{S}) = @simd for j=1:length(v) @inbounds propose!(d[j],v[j]) end

###Vectorized versions of the calculation of the logprobability of a point given a density
logprobability{D<:ProposalDensity,S<:MCSample}(d::Vector{D},s::S) = map((dj)->logprobability(dj,s),d)

###Vectorized version of show
function Base.show{P<:ProposalDensity}(io::IO,p::Array{P})
   println(typeof(p)," with ",length(p)," proposal densities")
  for i=1:length(p)
    show(io,p[i])
  end
end
