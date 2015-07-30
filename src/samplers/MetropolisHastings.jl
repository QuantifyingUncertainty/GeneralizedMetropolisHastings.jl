abstract MHSampler <: MCSampler #common base class for MH Samplers
abstract Metropolis <: MHSampler #derive from this for symmetric proposal densities
abstract MetropolisHastings <: MHSampler

### Type holding the parameters for a random-walk Metrolopolis sampler
immutable MHNormal <: Metropolis #mormal distribution is a symmetric Metropolis sampler
  nparas::Int
  covariance::Matrix{Float64}
  initialscaling::Float64
  function MHNormal(n::Int,c::Matrix{Float64},s::Float64)
    @assert n == size(c,1) == size(c,2) "Number of parameters should equal size of (square) covariance matrix"
    new(n,c,s)
  end
end

### Construct MHNormal based on number of parameters and initial scaling
MHNormal(Σ::Matrix{Float64},s::Float64) = MHNormal(size(Σ,1),Σ,s)
MHNormal(nparas::Int,s::Float64) = MHNormal(nparas,eye(nparas),s)

### Type holding the state of the Markov Chain for a Generalized M-H sampler
type MHHeap{P<:ProposalDensity} <: MCHeap{BaseSample}
  samples::Vector{BaseSample}
  sampledensities::Vector{P} #will remain empty for symmetric proposal densities
  fromdensity::P
  scaling::Float64
end

###Construct an MHHeap from a MHNormal sampler and the number of proposals per iteration
MHHeap(s::MHNormal,nprops::Int) = MHHeap{NormalDensity}([BaseSample(numparas(s)) for i=1:nprops],[NormalDensity(numparas(s)) for i=1:nprops],NormalDensity(s.initialscaling*s.initialscaling*s.covariance),s.initialscaling)

###Functions to set the current point to propose from
set_from!(s::MHSampler,h::MHHeap,from::BaseSample) = update_density!(h.fromdensity,from.values)

###Propose a new set of values for this BaseSample
propose!(s::MHSampler,h::MHHeap,b::BaseSample) = propose!(h.fromdensity,b)

###Per-sample update function for metropolis samplers and heaps is an empty function
update_proposal!(::Metropolis,::MCHeap,::Int) = ()

###Per-sample update function for MetropolisHastings samplers updates the
update_proposal!(s::MetropolisHastings,h::MCHeap,i::Int) = update_density!(h.sampledensities[i],h.samples[i])

### Base functionality
==(h1::MHHeap,h2::MHHeap) = (h1.samples == h2.samples && h1.sampledensities == h2.sampledensities && h1.fromdensity == h2.fromdensity && h1.scaling == h2.scaling)
