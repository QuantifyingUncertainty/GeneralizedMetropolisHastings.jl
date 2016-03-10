############################################################################
### ProposalDensities wrapping a distribution from the Distributions package
############################################################################
type SymmetricDistributionWrapper{D<:Distributions.Distribution} <: SymmetricDensity
    distribution::D
end

type ASymmetricDistributionWrapper{D<:Distributions.Distribution} <: ASymmetricDensity
    distribution::D
end

typealias DistributionWrapper Union{SymmetricDistributionWrapper, ASymmetricDistributionWrapper}

### Internal factory methods
@inline _density(::Type{Val{:true}},d::Distributions.Distribution) = SymmetricDistributionWrapper{typeof(d)}(d)
@inline _density(::Type{Val{:false}},d::Distributions.Distribution) = ASymmetricDistributionWrapper{typeof(d)}(d)

###Propose points by calling rand! on the corresponding distribution
@inline propose!(d::DistributionWrapper,s::AbstractSample) = (rand!(d.distribution,s.values) ; s)

###Calculate the logprobability by calling logpdf! on the corresponding distribution
@inline logprobability!(r::AbstractVector,d::DistributionWrapper,v::AbstractVector) = (@assert length(r) == 1 ; @inbounds r[1] = logpdf(d.distribution,v) ; r)
@inline logprobability!(r::AbstractVector,d::DistributionWrapper,v::AbstractMatrix) = (logpdf!(r,d.distribution,v) ; r)
@inline logprobability(d::DistributionWrapper,v::AbstractVector) = collect(logpdf(d.distribution,v))
@inline logprobability(d::DistributionWrapper,v::AbstractMatrix) = logprobability!(zeros(eltype(v),size(v,2)),d,v)

densityname(d::DistributionWrapper) = "DistributionWrapper"

###Base.show for DistributionWrappers
function show(io::IO,d::DistributionWrapper)
    println(io,"$(densityname(d)) with distribution:")
    show(io,d.distribution)
    println(io)
    nothing
end
