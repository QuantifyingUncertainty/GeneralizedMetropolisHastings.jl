typealias NormalDensity SymmetricDistributionWrapper{MvNormal}

###Internal factory function for a Normal proposal density function
@inline _density(::Type{Val{:normal}},x::AbstractVector,S::AbstractArray) = NormalDensity(MvNormal(x,S))

###The Normal distribution can be conditioned on a point by updating its mean parameter
@inline condition!(d::NormalDensity,x::AbstractArray) = copy!(d.distribution.μ,x)

###If we also need to update the covariance, then we need to recreate the distribution
@inline condition!(d::NormalDensity,x::AbstractArray,S::AbstractArray) = (d.distribution = MvNormal(x,S) ; d)

densityname(d::NormalDensity) = "NormalDensity"

###Nothing else is needed for a Normal proposal density:
###    - propose! is implemented for the super type DistributionWrapper
###    - logprobability! is implemented for the super type DistributionWrapper
###    - show is implemented for the super type DistributionWrapper
