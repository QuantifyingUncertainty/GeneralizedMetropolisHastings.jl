typealias LogNormalDensity ASymmetricDistributionWrapper{MvLogNormal}

###Internal factory function for a LogNormal proposal density
@inline _density(::Type{Val{:lognormal}},x::AbstractVector,S::AbstractArray) = LogNormalDensity(MvLogNormal(log(x),S))

###Conditioning the LogNormal on a point: conditioned point is the median, scale parameter is kept unchanged
@inline condition!(d::LogNormalDensity,x::AbstractArray) = map!(log,d.distribution.normal.μ,x)

###If we also need to update the scale, then we need to recreate the distribution
@inline condition!(d::LogNormalDensity,x::AbstractArray,S::AbstractArray) = (d.distribution = MvLogNormal(log(x),S) ; d)

densityname(d::LogNormalDensity) = "LogNormalDensity"

###Nothing else is needed for a LogNormal proposal density:
###    - propose! is implemented for the super type DistributionWrapper
###    - logprobability! is implemented for the super type DistributionWrapper
###    - show is implemented for the super type DistributionWrapper
