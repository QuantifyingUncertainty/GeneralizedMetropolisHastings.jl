abstract AbstractBactrian{F<:Distributions.VariateForm} <: Distributions.ContinuousDistribution{F}

type Bactrian{T<:AbstractFloat} <: AbstractBactrian{Distributions.Univariate}
    mixture::Distributions.UnivariateMixture
    μ::T
    σ::T
    m::T
end

@inline _bactrianscale(::Type{Val{:normal}},s::AbstractFloat,m::AbstractFloat) = s*sqrt(1-m*m)
@inline _bactrianscale(::Type{Val{:laplace}},s::AbstractFloat,m::AbstractFloat) = s*sqrt(sqrt(2)*(1-m*m))
@inline _bactrianscale(::Type{Val{:triangular}},s::AbstractFloat,m::AbstractFloat) = s*sqrt(6*(1-m*m))

@inline _bactrianscale(d::Distributions.Normal,s::AbstractFloat,m::AbstractFloat) = _bactrianscale(Val{:normal},s,m)
@inline _bactrianscale(d::Distributions.Laplace,s::AbstractFloat,m::AbstractFloat) = _bactrianscale(Val{:laplace},s,m)
@inline _bactrianscale(d::Distributions.SymTriangularDist,s::AbstractFloat,m::AbstractFloat) = _bactrianscale(Val{:triangular},s,m)


@inline function _mixturemodel(subtype::Symbol,x::AbstractFloat,s::AbstractFloat,m::AbstractFloat)
    bs = _bactrianscale(Val{subtype},s,m)
    Distributions.MixtureModel([_distribution(Val{subtype},x-m*s,bs),_distribution(Val{subtype},x+m*s,bs)],[1/2,1/2])
end

@inline function _distribution(::Type{Val{:bactrian}},x::AbstractFloat,s::AbstractFloat,subtype::Symbol,mixing::AbstractFloat)
    Bactrian{typeof(x)}(_mixturemodel(subtype,x,s,mixing),x,s,mixing)
end

function recenter(b::Bactrian,x)
    b.μ = x
    b.mixture.components[1] = recenter(b.mixture.components[1],x-b.m*b.σ)
    b.mixture.components[2] = recenter(b.mixture.components[2],x+b.m*b.σ)
    b
end

function rescale(b::Bactrian,s)
    b.σ = s*b.σ
    bs = _bactrianscale(b.mixture.components[1],b.σ,b.m)
    b.mixture.components[1] = update(b.mixture.components[1],b.μ-b.m*b.σ,bs)
    b.mixture.components[2] = update(b.mixture.components[2],b.μ+b.m*b.σ,bs)
    b
end

function update(b::Bactrian,x,s)
    b.μ = x
    b.σ = s
    bs = _bactrianscale(b.mixture.components[1],s,b.m)
    b.mixture.components[1] = update(b.mixture.components[1],x-b.m*s,bs)
    b.mixture.components[2] = update(b.mixture.components[2],x+b.m*s,bs)
    b
end

rand(b::Bactrian) = rand(b.mixture)
rand!(b::Bactrian,x::AbstractArray) = rand!(b.mixture,x)
pdf(b::Bactrian,x::Real) = pdf(b.mixture,x)

mean(b::Bactrian) = mean(b.mixture)
var(b::Bactrian) = var(b.mixture)
scale(b::Bactrian) = b.σ


