###Local, non-exported helper types defining default behaviour
typealias ContinuousType eltype(Distributions.Continuous)
typealias DiscreteType eltype(Distributions.Discrete)
typealias ValueType Union{ContinuousType,DiscreteType}
typealias UnivariatePriorType Distributions.UnivariateDistribution

###Local, non-exported helper functions defining default behaviour
@inline _defaultparametervalue{T<:AbstractFloat}(::Type{T}) = zero(T)
@inline _defaultparametervalue(p::UnivariatePriorType,::Type{ContinuousType}) = median(p)
@inline _defaultparametervalue(p::UnivariatePriorType,::Type{DiscreteType}) = trunc(DiscreteType,median(p))
@inline _defaultparameterprior(l::ContinuousType,h::ContinuousType) = Distributions.Uniform(l,h)
@inline _defaultparameterprior(l::DiscreteType,h::DiscreteType) = Distributions.DiscreteUniform(l,h)

###Unnamed parameter which has a default value but no prior
immutable UnnamedParameterDefault{T<:AbstractFloat} <: MCParameter
  default::T
end

###Named MCMC parameter which has a default value but no prior
immutable NamedParameterDefault{N<:AbstractString,T<:AbstractFloat} <: MCParameter
  name::AbstractString
  default::T
  NamedParameterDefault(n::N,d::T) = new(n,d)
end

###Unnamed parameter which has a prior distribution and/or a default value
immutable UnnamedParameterUnivariate{P<:UnivariatePriorType,V<:ValueType} <: MCParameter
  prior::P
  default::V
  UnnamedParameterUnivariate(p::P,d::V) = (@assert Distributions.insupport(p,d) ; new(p,d))
end

###Named MCMC parameter which has a prior distribution and/or a default value
immutable NamedParameterUnivariate{N<:AbstractString,P<:UnivariatePriorType,V<:ValueType} <: MCParameter
  name::N
  prior::P
  default::V
  NamedParameterUnivariate(n::N,p::P,d::V) = (@assert Distributions.insupport(p,d) ; new(n,p,d))
end

typealias ParameterDefault Union{UnnamedParameterDefault,NamedParameterDefault}
typealias ParameterUnivariate Union{UnnamedParameterUnivariate,NamedParameterUnivariate}

###Factory functions for unnamed parameters with only a default value
parameter(d::AbstractFloat) = UnnamedParameterDefault{typeof(d)}(d)
parameter{T<:AbstractFloat}(::Type{T} =DefaultFloatType) = parameter(_defaultparametervalue(T))

###Factory functions for named parameters with only a default value
parameter(n::AbstractString,d::AbstractFloat) = NamedParameterDefault{typeof(n),typeof(d)}(n,d)
parameter{T<:AbstractFloat}(n::AbstractString,::Type{T} =DefaultFloatType) = parameter(n,_defaultparametervalue(T))

###Factory functions for unnamed parameters with a prior distribution
parameter(p::UnivariatePriorType,d::ValueType) = UnnamedParameterUnivariate{typeof(p),eltype(p)}(p,d)
parameter(p::UnivariatePriorType) = parameter(p,_defaultparametervalue(p,eltype(p)))
parameter{V<:ValueType}(l::V,h::V,d::V) = parameter(_defaultparameterprior(l,h),d)
parameter{V<:ValueType}(l::V,h::V) = parameter(_defaultparameterprior(l,h))

###Factory functions for named parameters with a prior distribution
parameter(n::AbstractString,p::UnivariatePriorType,d::ValueType) = NamedParameterUnivariate{typeof(n),typeof(p),eltype(p)}(n,p,d)
parameter(n::AbstractString,p::UnivariatePriorType) = parameter(n,p,_defaultparametervalue(p,eltype(p)))
parameter{V<:ValueType}(n::AbstractString,l::V,h::V,d::V) = parameter(n,_defaultparameterprior(l,h),d)
parameter{V<:ValueType}(n::AbstractString,l::V,h::V) = parameter(n,_defaultparameterprior(l,h))

###Vectorised factory functions
parameters{T<:AbstractFloat}(v::AbstractVector{T}) = map(parameter,v)
parameters{T<:AbstractFloat}(n::Integer,t::Type{T} =DefaultFloatType) = [parameter(t) for i=1:n]
parameters{N<:AbstractString,T<:AbstractFloat}(n::AbstractVector{N},v::AbstractVector{T}) = (@assert length(n) == length(v) ; map(parameter,n,v))
parameters{N<:AbstractString,T<:AbstractFloat}(n::AbstractVector{N},t::Type{T} =DefaultFloatType) = [parameter(ni,t) for ni in n]

parameters{P<:UnivariatePriorType,V<:ValueType}(p::AbstractVector{P},v::AbstractVector{V}) = (@assert length(p) == length(v) ; map(parameter,p,v))
parameters{P<:UnivariatePriorType}(p::AbstractVector{P}) = map(parameter,p)
parameters{V<:ValueType}(l::Vector{V},h::Vector{V},d::Vector{V}) = (@assert length(l) == length(h) == length(d) ; map(parameter,l,h,d))
parameters{V<:ValueType}(l::Vector{V},h::Vector{V}) = (@assert length(l) == length(h) ; map(parameter,l,h))

parameters{N<:AbstractString,P<:UnivariatePriorType,V<:ValueType}(n::AbstractVector{N},p::AbstractVector{P},v::AbstractVector{V}) = (@assert length(n) == length(p) == length(v) ; map(parameter,n,p,v))
parameters{N<:AbstractString,P<:UnivariatePriorType}(n::AbstractVector{N},p::AbstractVector{P}) = (@assert length(n) == length(p) ; map(parameter,n,p))
parameters{N<:AbstractString,V<:ValueType}(n::AbstractVector{N},l::Vector{V},h::Vector{V},d::Vector{V}) = (@assert length(n) == length(l) == length(h) == length(d) ; map(parameter,n,l,h,d))
parameters{N<:AbstractString,V<:ValueType}(n::AbstractVector{N},l::Vector{V},h::Vector{V}) = (@assert length(n) == length(l) == length(h) ; map(parameter,n,l,h))

###Functionality to initialize MCMC from parameter definitions
#for individual parameters
@inline _initvalue{P<:ParameterUnivariate,T<:AbstractFloat}(::Type{ValuesFromPrior},p::P,::Type{T} =DefaultFloatType) = T(rand(p.prior)) #special case for initializing from prior
@inline _initvalue{V<:ValuesFrom,T<:AbstractFloat}(::Type{V},p::MCParameter,::Type{T} =DefaultFloatType) = T(p.default) #general case initializes to default value

#vectorized internal function
@inline _initvalues!{V<:ValuesFrom,P<:MCParameter,T<:AbstractFloat}(::Type{V},p::AbstractVector{P},v::AbstractVector{T}) = (@simd for i=1:length(p) @inbounds v[i] = _initvalue(V,p[i],T) end ; v)

#exported initvalue functions
initvalues!{V<:ValuesFrom,P<:MCParameter,T<:AbstractFloat}(::Type{V},p::AbstractVector{P},v::AbstractVector{T}) = (@assert length(p) == length(v) ; _initvalues!(V,p,v))
initvalues{V<:ValuesFrom,P<:MCParameter,T<:AbstractFloat}(::Type{V},p::AbstractVector{P},::Type{T} =DefaultFloatType) = _initvalues!(V,p,Vector{T}(length(p)))

###Functionality to calculate the logprior
@inline _logprior{P<:ParameterUnivariate,T<:AbstractFloat}(p::P,v::T) = T(Distributions.logpdf(p.prior,v))
@inline _logprior{P<:ParameterDefault,T<:AbstractFloat}(p::P,v::T) = zero(T)

###vectorized internal logprior function
@inline _logprior{P<:MCParameter,T<:AbstractFloat}(p::AbstractVector{P},v::AbstractVector{T}) = (s = zero(T) ; @simd for i=1:length(p) @inbounds s += _logprior(p[i],v[i]) end ; s)
@inline _logprior{P<:ParameterDefault,T<:AbstractFloat}(p::AbstractVector{P},v::AbstractVector{T}) = zero(T)

###exported logprior function
logprior{P<:MCParameter,T<:AbstractFloat}(p::P,v::T) = _logprior(p,v)
logprior{P<:MCParameter,T<:AbstractFloat}(p::AbstractVector{P},v::AbstractVector{T}) = (@assert length(p) == length(v) ; _logprior(p,v))

###Overloaded functions and operators from Base package
import Base.==
function =={P<:MCParameter}(p1::P,p2::P)
  for f in fieldnames(p1)
    isequal(getfield(p1,f),getfield(p2,f))?continue:(return false)
  end
  return true
end

function Base.show{P<:MCParameter}(io::IO,p::P,n::AbstractString ="")
  println(io,n,typeof(p).name.name)
  for f in fieldnames(p)
    println(io," ",f,": ",getfield(p,f))
  end
  nothing
end

function Base.show{P<:MCParameter}(io::IO,v::AbstractVector{P})
  e = isa(eltype(v),Union)?"MCParameter":eltype(v).name.name
  println(io)
  println(io,"Array{$e} with")
  for i=1:length(v)
    show(io,v[i],"[$i] ")
  end
  println(io)
end

Base.display{P<:MCParameter}(io::IO,v::AbstractVector{P}) = Base.show{P}(io,v)







