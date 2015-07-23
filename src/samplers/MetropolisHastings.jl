abstract MetropolisHastings <: MCSampler #derive from this for non-symmetric proposal densities
abstract Metropolis <: MetropolisHastings #derive from this for symmetric proposal densities

### Type holding the parameters for a random-walk Metrolopolis sampler
immutable MHNormal <: Metropolis #mormal distribution is symmetric Metropolis sampler
  covariance::Matrix{Float64}
  initialscaling::Float64
end
### Construct MHNormal based on number of parameters and initial scaling
MHNormal(nparas::Int,s::Float64) = MHNormal(eye(nparas),s)
nparas(s::MHNormal) = size(s.covariance,1)

### Type holding the state of the Markov Chain for a Generalized M-H sampler
type MHHeap <: MCHeap{BaseSample}
  samples::Vector{BaseSample}
  sampledensity::ProposalDensity
  fromvalues::Vector{Float64}
  scaling::Float64
end

### Construct an MHHeap from a MHNormal sampler and the number of proposals per iteration
MHHeap(s::MHNormal,nprops::Int) = MHHeap([BaseSample(nparas(s)) for i=1:nprops],NormalDensity(s.initialscaling*s.covariance),zeros(nparas(s)),s.initialscaling)

### Functions to set the current point to propose from
set_from!(s::MetropolisHastings,heap::MHHeap,from::BaseSample) = (heap.fromvalues = copy(from.values))

### Propose a new set of values for this BaseSample
function propose!(s::MetropolisHastings,h::MHHeap,b::BaseSample)
  propose!(h.sampledensity,b)
  @simd for i=1:length(h.fromvalues) #explicit for loop to avoid memory allocation of += operator
    @inbounds b.values[i] += h.fromvalues[i]
  end
end

### Base functionality
==(h1::MHHeap,h2::MHHeap) = (h1.samples == h2.samples && h1.sampledensity == h2.sampledensity && h1.fromvalues == h2.fromvalues && h1.scaling == h2.scaling)
