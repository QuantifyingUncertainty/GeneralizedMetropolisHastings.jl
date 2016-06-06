### Generic factory function for distributions
distribution(s::Symbol,args...;keyargs...) = _distribution(Val{s},args...;keyargs...)
distributions(s::Symbol,args...;keyargs...) = _distributions(Val{s},args...;keyargs...) #vectorized version

### Specific factory functions for known distributions
### Normal distribution (univariate and multivariate)
@inline _distribution(::Type{Val{:normal}},x::Real,s::Real) = Distributions.Normal(x,s)
@inline _distribution(::Type{Val{:normal}},x::AbstractVector,s::AbstractArray) = Distributions.MvNormal(copy(x),copy(s))
@inline _distribution(::Type{Val{:normal}},s::AbstractArray) = Distributions.MvNormal(copy(s))

### Lognormal distribution (univariate and multivariate)
@inline _distribution(::Type{Val{:lognormal}},x::Real,s::Real) = Distributions.LogNormal(log(x),s)
@inline _distribution(::Type{Val{:lognormal}},x::AbstractVector,s::AbstractArray) = Distributions.MvLogNormal(log(x),copy(s))

### Uniform distribution (univariate)
@inline _distribution(::Type{Val{:uniform}},μ,s) = Distributions.Uniform(μ-s/2,μ+s/2)
@inline _distribution(::Type{Val{:laplace}},μ,s) = Distributions.Laplace(μ,s)
@inline _distribution(::Type{Val{:triangular}},μ,s) = Distributions.SymTriangularDist(μ,s)

### Vectorized distribution factory function for location/scale univariate distributions
@inline function _distributions(distributiontype::DataType,location::AbstractVector,scale::AbstractVector,args...)
    @assert length(location) == length(scale)
    map((μ,σ)->_distribution(distributiontype,μ,σ,args...),location,scale)
end

#Update the center point for distributions
recenter(d::Distributions.Normal,x) = Distributions.Normal(x,Distributions.std(d))
recenter(d::Distributions.LogNormal,x) = Distributions.LogNormal(log(x),Distributions.params(d)[2])
recenter(d::Distributions.Uniform,x) = _distribution(Val{:uniform},x,Distributions.scale(d))
recenter(d::Distributions.Laplace,x) = _distribution(Val{:laplace},x,Distributions.scale(d))
recenter(d::Distributions.SymTriangularDist,x) = _distribution(Val{:triangular},x,Distributions.scale(d))
recenter(d::Distributions.MvNormal,x) = (copy!(d.μ,x) ; d)
recenter(d::Distributions.MvLogNormal,x) = (map!(log,d.normal.μ,x) ; d)

#Scale distributions
@inline _rescale(::PDMats.PDMat,d::Distributions.MvNormal,s) = (scale!(d.Σ.mat,s*s) ; scale!(d.Σ.chol.factors,s*s) ; d)
@inline _rescale(::PDMats.PDiagMat,d::Distributions.MvNormal,s) = (scale!(d.Σ.diag,s*s) ; scale!(d.Σ.inv_diag,1/s/s) ; d)
@inline _rescale(::Union{PDMats.PDMat,PDMats.PDiagMat},d::Distributions.MvLogNormal,s) = (rescale(d.normal,s) ; d)

rescale(d::Distributions.Normal,s) = Distributions.Normal(d.μ,d.σ*s)
rescale(d::Distributions.MvNormal,s) = _rescale(d.Σ,d,s)
rescale(d::Distributions.LogNormal,s) = Distributions.LogNormal(d.μ,d.σ*s)
rescale(d::Distributions.MvLogNormal,s) = _rescale(d.normal.Σ,d,s)
rescale(d::Distributions.Uniform,s) = (sold = scale(d) ; x = (d.a + d.b)/2 ; Distributions.Uniform(x-sold*s/2,x+sold*s/2))
rescale(d::Distributions.Laplace,s) = Distributions.Laplace(d.μ,s*d.θ)
rescale(d::Distributions.SymTriangularDist,s) = Distributions.SymTriangularDist(d.μ,s*d.σ)

#Update distributions (generally means recreating them)
update(d::Union{Distributions.Normal,Distributions.MvNormal},x,s) = _distribution(Val{:normal},x,s)
update(d::Distributions.MvNormal,s) = _distribution(Val{:normal},s)
update(d::Union{Distributions.LogNormal,Distributions.MvLogNormal},x,s) = _distribution(Val{:lognormal},x,s)
update(d::Distributions.Uniform,x,s) = _distribution(Val{:uniform},x,s)
update(d::Distributions.Laplace,x,s) = _distribution(Val{:laplace},x,s)
update(d::Distributions.SymTriangularDist,x,s) = _distribution(Val{:triangular},x,s)
