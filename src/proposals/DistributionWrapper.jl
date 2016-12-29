###########################################################################
### ProposalDensities wrapping distributions from the Distributions package
###########################################################################
type SingleSymmetricDistributionWrapper{D<:Distributions.Distribution} <: SymmetricDensity
    distribution::D
end

type SingleASymmetricDistributionWrapper{D<:Distributions.Distribution} <: ASymmetricDensity
    distribution::D
end

type CompoundSymmetricDistributionWrapper{D<:Distributions.UnivariateDistribution} <: SymmetricDensity
    distributions::Vector{D}
end

type CompoundASymmetricDistributionWrapper{D<:Distributions.UnivariateDistribution} <: ASymmetricDensity
    distributions::Vector{D}
end

typealias SingleDistributionWrapper Union{SingleSymmetricDistributionWrapper, SingleASymmetricDistributionWrapper}
typealias CompoundDistributionWrapper Union{CompoundSymmetricDistributionWrapper, CompoundASymmetricDistributionWrapper}
typealias DistributionWrapper Union{SingleDistributionWrapper,CompoundDistributionWrapper}

densityname(d::DistributionWrapper) = "DistributionWrapper"

typealias NormalDensity SingleSymmetricDistributionWrapper{Distributions.MvNormal}
typealias LogNormalDensity SingleASymmetricDistributionWrapper{Distributions.MvLogNormal}
typealias UniformDensity CompoundSymmetricDistributionWrapper{Distributions.Uniform}
typealias LaplaceDensity CompoundSymmetricDistributionWrapper{Distributions.Laplace}
typealias TriangularDensity CompoundSymmetricDistributionWrapper{Distributions.SymTriangularDist}
typealias BactrianDensity CompoundSymmetricDistributionWrapper{Bactrian}

_density(::Type{Val{:normal}},x::AbstractVector,s::AbstractArray) = NormalDensity(_distribution(Val{:normal},x,s))
_density(::Type{Val{:normal}},s::AbstractArray) = NormalDensity(_distribution(Val{:normal},s))
_density(::Type{Val{:lognormal}},x::AbstractVector,s::AbstractArray) = LogNormalDensity(_distribution(Val{:lognormal},x,s))
_density(::Type{Val{:uniform}},x::AbstractVector,s::AbstractVector) = UniformDensity(_distributions(Val{:uniform},x,s))
_density(::Type{Val{:laplace}},x::AbstractVector,s::AbstractVector) = LaplaceDensity(_distributions(Val{:laplace},x,s))
_density(::Type{Val{:triangular}},x::AbstractVector,s::AbstractVector) = TriangularDensity(_distributions(Val{:triangular},x,s))
_density(::Type{Val{:bactrian}},x::AbstractVector,s::AbstractVector,subtype::Symbol,mixing::AbstractFloat) = BactrianDensity(_distributions(Val{:bactrian},x,s,subtype,mixing))

###Propose points by calling rand! on the corresponding distribution
@inline propose!(d::SingleDistributionWrapper,v::AbstractArray) = (rand!(d.distribution,v) ; v)

@inline function propose!(d::CompoundDistributionWrapper,v::AbstractArray)
    nrows = size(v,1)
    for j=1:size(v,2)
        @simd for i=1:nrows
            @inbounds v[i,j] = rand(d.distributions[i])
        end
    end
    v
end

###Calculate the logprobability by calling logpdf! on the corresponding distribution
@inline logprobability!(r::AbstractVector,d::SingleDistributionWrapper,v::AbstractVector) = (@inbounds r[1] = logpdf(d.distribution,v) ; r)
@inline logprobability!(r::AbstractVector,d::SingleDistributionWrapper,v::AbstractMatrix) = (logpdf!(r,d.distribution,v) ; r)
@inline logprobability(d::SingleDistributionWrapper,v::AbstractArray) = logprobability!(zeros(eltype(v),size(v,2)),d,v)

@inline function logprobability!(r::AbstractVector,d::CompoundDistributionWrapper,v::AbstractArray)
    nrows = size(v,1)
    for j=1:size(v,2)
        @inbounds r[j] = zero(eltype(v))
        @simd for i=1:nrows
            @inbounds r[j] += logpdf(d.distributions[i],v[i,j])
        end
    end
    r
end

@inline logprobability(d::CompoundDistributionWrapper,v::AbstractArray) = logprobability!(zeros(eltype(v),size(v,2)),d,v)

#update functions for single distributions
condition!(d::SingleDistributionWrapper,x::AbstractVector) = (d.distribution = recenter(d.distribution,x) ; d)
scale!(d::SingleDistributionWrapper,s::AbstractFloat) = (d.distribution = rescale(d.distribution,s) ; d)
update!(d::SingleDistributionWrapper,x::AbstractVector,s::AbstractArray) = (d.distribution = update(d.distribution,x,s) ; d)
update!(d::NormalDensity,s::AbstractArray) = (d.distribution = update(d.distribution,s) ; d)

#update functions for compound distributions
condition!(d::CompoundDistributionWrapper,x::AbstractVector) = (map!(recenter,d.distributions,d.distributions,x) ; d)
scale!(d::CompoundDistributionWrapper,s::AbstractFloat) = (map!((d1)->rescale(d1,s),d.distributions,d.distributions) ; d)
update!(d::CompoundDistributionWrapper,x::AbstractVector,s::AbstractVector) = (map!(update,d.distributions,d.distributions,x,s) ; d)

###Base.show for DistributionWrappers
function show(io::IO,d::DistributionWrapper)
    println(io,"$(densityname(d)) with fields:")
    for f in fieldnames(d)
        println("  ",f,": ",getfield(d,f))
    end
    nothing
end
