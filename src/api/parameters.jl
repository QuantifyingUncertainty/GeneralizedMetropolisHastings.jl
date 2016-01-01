###MCMC Parameter definition which is a scalar
###Has a default value and a univariate prior, which can be discrete or continuous
immutable MCParameterScalar{S<:Distributions.ValueSupport} <: MCParameter
  name::AbstractString
  default::eltype(S)
  prior::Distributions.UnivariateDistribution{S}
end

###Local, non-exported helper types defining default behaviour
typealias ContinuousType eltype(Distributions.Continuous)
typealias DiscreteType eltype(Distributions.Discrete)
typealias ContinuousDefaultPrior Distributions.Uniform
typealias DiscreteDefaultPrior Distributions.DiscreteUniform

###Local, non-exported helper functions defining default behaviour
const defaultparametername::AbstractString = "Parameter"
defaultparametervalue{S<:Distributions.ValueSupport}(::Type{S}) = zero(eltype(S))
defaultparameterprior(l::ContinuousType,h::ContinuousType) = ContinuousDefaultPrior(l,h)
defaultparameterprior(l::DiscreteType,h::DiscreteType) = DiscreteDefaultPrior(l,h)
defaultparameterprior(::Type{Distributions.Continuous}) = ContinuousDefaultPrior(-1e43,1e43) #a continuous uniform covering quasi all 64-bit floating-point numbers
defaultparameterprior(::Type{Distributions.Discrete}) = DiscreteDefaultPrior(-2^61,2^61) #a continuous uniform covering quasi all 64-bit integers

###Constructor functions of scalar parameter definitions
parameter{S<:Distributions.ValueSupport}(n::AbstractString =defaultparametername,d::eltype{S} =defaultparametervalue(S),p::Distributions.UnivariateDistribution{S} =defaultparameterprior(S)) = (@assert insupport(p,d) ; MCParameterScalar(n,d,p))
parameter{S<:Distributions.ValueSupport}(n::AbstractString,d::eltype{S},l::eltype{S},h::eltype{S}) = parameter(n,d,defaultparameterprior(l,h))
parameter{S<:Distributions.ValueSupport}(d::eltype{S},p::Distributions.UnivariateDistribution{S} =defaultparameterprior(S)) = parameter(defaultparametername,d,p)
parameter{S<:Distributions.ValueSupport}(d::eltype{S},l::eltype{S},h::eltype{S}) = parameter(defaultparametername,d,defaultparameterprior(l,h))

###Vectorized constructor functions of scalar parameter definitions
parameters{S<:Distributions.ValueSupport}(n::Vector{AbstractString},d::Vector{elval(S)},p::Vector{Distributions.UnivariateDistribution{S}}) = (@assert length(n) == length(d) == length(p) ; map(parameter,n,d,p))
parameters{S<:Distributions.ValueSupport}(n::Vector{AbstractString},d::Vector{elval(S)},l::Vector{elval(S)},h::Vector{elval(S)}) = (@assert length(n) == length(d) == length(l) == length(h) ; map(parameter,n,d,l,h))
parameters{S<:Distributions.ValueSupport}(n::Vector{AbstractString},d::Vector{elval(S)}) = (@assert length(n) == length(d) ; map(parameter,n,d))
parameters{S<:Distributions.ValueSupport}(d::Vector{elval(S)},p::Vector{Distributions.UnivariateDistribution{S}}) = (@assert length(d) == length(p) ; map(parameter,d,p))
parameters{S<:Distributions.ValueSupport}(d::Vector{elval(S)},l::Vector{elval(S)},h::Vector{elval(S)}) = (@assert length(d) == length(l) == length(h) ; map(parameter,d,l,h))
parameters{S<:Distributions.ValueSupport}(n::Vector{AbstractString}) = map(parameter,n)
parameters{S<:Distributions.ValueSupport}(d::Vector{elval(S)}) = map(parameter,d)
parameters{S<:Distributions.ValueSupport}(::Type{S},nparas::Int) = [parameter{S}() for i=1:nparas]

###Functionality to initialize values from parameter definitions as floating point values
_initvalue{T<:AbstractFloat}(::Type{ValuesFromDefault},p::MCParameter) = T(p.default)
_initvalue{T<:AbstractFloat}(::Type{ValueFromPrior},p::MCParameter) = T(rand(p.prior))
_initvalues!{V<:ValuesFrom,T<:AbstractFloat}(::Type{V},p::Vector{MCParameter},v::Vector{T}) = (@simd for i=1:length(p) @inbounds v[i] = _initvalue(V,p[i]) end ; v)
initvalues!{V<:ValuesFrom,T<:AbstractFloat}(::Type{V},p::Vector{MCParameter},v::Vector{T}) = (@assert length(p) == length(v) ; _initvalues(V,p,v))
initvalues{V<:ValuesFrom,T<:AbstractFloat}(::Type{V},::Type{T},p::Vector{MCParameter}) = _initvalues(V,p,Vector{T}(length(p)))
initvalues{V<:ValuesFrom}(::Type{V},p::Vector{MCParameter}) = initvalues(V,Float64,p)


###Overloaded functions and operators from Base package
function Base.show(io::IO,p::MCParameterScalar)
  println(io,"Scalar Parameter with name: \"",p.name,"\", default: ",p.default," prior: ",p.prior)
  nothing
end





