###Vectorized versions of the in-place proposal function for arbitrary densities and samples
propose!{D<:ProposalDensity,S<:MCSample}(d::D,v::Vector{S}) = @simd for j=1:length(v) @inbounds propose!(d,v[j]) end
propose!{D<:ProposalDensity,S<:MCSample}(d::Vector{D},v::Vector{S}) = @simd for j=1:length(v) @inbounds propose!(d[j],v[j]) end

###Vectorized version of show
function Base.show{P<:ProposalDensity}(io::IO,p::Array{P})
   println(typeof(p)," with ",length(p)," proposal densities")
  for i=1:length(p)
    show(io,p[i])
  end
end
