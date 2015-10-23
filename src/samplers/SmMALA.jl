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
  ntangent::Int
  initialscaling::Float64
end
@compat typealias SmMALANormalFamily Union{SmMALANormal,TrSmMALANormal,TrSmMALARandomNormal}
@compat typealias TrSmMALANormalFamily Union{TrSmMALANormal,TrSmMALARandomNormal}
@compat typealias SmMALAFullTensorNormalFamily Union{SmMALANormal,TrSmMALANormal}

### Number of a parameters
nparas(s::SmMALANormalFamily) = s.nparas

### Type holding the state of the Markov Chain for a Generalized SmMALA sampler
type SmMALAHeap{T<:MCSample{SecondOrder},P<:ProposalDensity} <: MCHeap{MCSample{SecondOrder}}
  samples::Vector{T}
  sampledensities::Vector{P}
  fromdensity::P
  scaling::Float64
end

### Construct an SmMALAHeap with normal proposal densities and the number of proposals per iteration
SmMALAHeap(s::SmMALAFullTensorNormalFamily,nprops::Int) =
  SmMALAHeap{TensorSample,NormalDensity}([TensorSample(nparas(s)) for i=1:nprops],[NormalDensity(nparas(s)) for i=1:nprops],NormalDensity(nparas(s)),s.initialscaling)
SmMALAHeap(s::TrSmMALARandomNormal,nprops::Int) =
  SmMALAHeap{ApproximateTensorSample,NormalDensity}([ApproximateTensorSample(nparas(s),s.ntangent) for i=1:nprops],[NormalDensity(nparas(s)) for i=1:nprops],NormalDensity(nparas(s)),s.initialscaling)

###Define the create_heap factory function for SmMALAHeap
create_heap(s::SmMALAFullTensorNormalFamily,nprops::Int) = SmMALAHeap(s,nprops)
create_heap(s::TrSmMALARandomNormal,nprops::Int) = SmMALAHeap(s,nprops)

### Set the point to sample from
set_from!(s::SimpleManifoldMALA,h::SmMALAHeap,from::TensorSample) = update_proposal!(s,h.fromdensity,from,h.scaling)

### Propose a new set of values for sample j
propose!{T<:MCSample{SecondOrder}}(s::SimpleManifoldMALA,h::SmMALAHeap,t::T) = propose!(h.fromdensity,t)

### Calculate the logprobability of a proposal

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
