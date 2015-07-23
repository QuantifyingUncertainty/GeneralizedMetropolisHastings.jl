abstract SimpleManifoldMALA <: MCSampler
abstract SmMALA <: SimpleManifoldMALA
abstract TrSmMALA <: SimpleManifoldMALA
abstract TrSmMALARandom <: SimpleManifoldMALA

### Type holding the parameters for normal-density SmMALA-family samplers
### Most SmMALA samplers are used to simulate physical processes with Brownian motion
type SmMALANormal <: SmMALA
  nparas::Int
  initialscaling::Float64
end
type TrSmMALANormal <: TrSmMALA
  nparas::Int
  initialscaling::Float64
end
type TrSmMALARandomNormal <: TrSmMALARandom
  nparas::Int
  initialscaling::Float64
end
typealias SmMALANormalFamily Union(SmMALANormal,TrSmMALANormal,TrSmMALARandomNormal)
typealias TrSmMALANormalFamily Union(TrSmMALANormal,TrSmMALARandomNormal)

### Number of a parameters
nparas(s::SmMALANormalFamily) = s.nparas

### Type holding the state of the Markov Chain for a Generalized SmMALA sampler
type SmMALAHeap <: MCHeap{TensorSample}
  samples::Vector{TensorSample}
  sampledensities::Vector{ProposalDensity}
  fromdensity::ProposalDensity
  scaling::Float64
end

### Construct an SmMALAHeap with normal proposal densities and the number of proposals per iteration
SmMALAHeap(s::SmMALANormalFamily,nprops::Int) = SmMALAHeap([TensorSample(nparas(s)) for i=1:nprops],[NormalDensity(nparas(s)) for i=1:nprops],NormalDensity(nparas(s)),s.initialscaling)

### Set the point to sample from
set_from!(s::SimpleManifoldMALA,h::SmMALAHeap,from::TensorSample) = update_proposal!(s,h.fromdensity,from,h.scaling)

### Propose a new set of values for sample j
propose!(s::SimpleManifoldMALA,h::SmMALAHeap,t::TensorSample) = propose!(h.fromdensity,t)

### Update the sampling density from the values of the TensorSample
update_proposal!(s::SimpleManifoldMALA,h::SmMALAHeap,j::Int) = update_proposal!(s,h.sampledensities[j],h.samples[j],h.scaling)

### Base functionality
==(h1::SmMALAHeap,h2::SmMALAHeap) = (h1.samples == h2.samples && h1.sampledensities == h2.sampledensities && h1.fromdensity == h2.fromdensity && h1.scaling == h2.scaling)

#############LOCAL FUNCTIONS##########################

### Update function for normal densities
update_proposal!(s::SmMALANormalFamily,d::NormalDensity,t::TensorSample,scaling::Float64) = update_density!(d,smmalamean(s,t,scaling),smmalacov(s,t,scaling))

### Helper functions to calculate SmMALA mean and covariance for normal densities
smmalamean(s::SmMALANormal,t::TensorSample,scaling::Float64) = t.values + (t.tensorloglikelihood\t.gradloglikelihood)*scaling*scaling/2
smmalacov(s::SmMALANormal,t::TensorSample,scaling::Float64) = inv(t.tensorloglikelihood)*scaling*scaling

### Helper functions to calculate the trust region SmMALA mean and covariance for normal densities
smmalamean(s::TrSmMALANormalFamily,t::TensorSample,scaling::Float64) = t.values + (t.tensorloglikelihood+eye(s.nparas)/scaling)\t.gradloglikelihood
smmalacov(s::TrSmMALANormalFamily,t::TensorSample,scaling::Float64) = inv(t.tensorloglikelihood + eye(s.nparas)/scaling)
