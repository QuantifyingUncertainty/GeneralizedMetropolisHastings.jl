###Local, non-exported helper types defining default behaviour
typealias ValueType Union{Distributions.Continuous,Distributions.Discrete}

### Paramater base type
abstract Parameter

###Named MCMC parameter which has a default value but no prior
immutable ParameterDefault{N<:Number} <: Parameter
  key::AbstractString
  default::N
  ParameterDefault(key::AbstractString,def::N) = new(key,def)
end

###Named MCMC parameter which has a prior distribution and/or a default value
immutable ParameterUnivariate{U<:Distributions.UnivariateDistribution,V<:ValueType} <: Parameter
  key::AbstractString
  prior::U
  default::V
  ParameterUnivariate(key::AbstractString,pri::U,def::V) = (@assert Distributions.insupport(pri,def) ; new(key,pri,def))
end

###Local, non-exported helper functions defining default behaviour
@inline _defaultkey() = ""
@inline _defaultvalue{N<:Number}(::Type{N}) = zero(N)
@inline _defaultvalue(::Type{Distributions.Continuous}) = mean(p)
@inline _defaultvalue(::Type{Distributions.Discrete}) = trunc(DiscreteType,mean(p))
@inline _defaultprior(l::Distributions.Continuous,h::Distributions.Continuous) = Distributions.Uniform(l,h)
@inline _defaultprior(l::Distributions.Discrete,h::Distributions.Discrete) = Distributions.DiscreteUniform(l,h)

###Factory functions for named parameters with only a default value
parameter(k::AbstractString,d::Number) = ParameterDefault{typeof(d)}(k,d)
parameter{N<:Number}(k::AbstractString,::Type{N}) = parameter(k,_defaultvalue(N))

###Factory functions for unnamed parameters with only a default value
parameter(d::Number) = parameter(_defaultkey(),d)
parameter{N<:Number}(::Type{N}) = parameter(_defaultvalue(N))

###Factory functions for named parameters with a prior distribution
parameter(k::AbstractString,p::Distributions.UnivariateDistribution,d::ValueType) = ParameterUnivariate{typeof(p),eltype(p)}(k,p,d)
parameter(k::AbstractString,p::Distributions.UnivariateDistribution) = parameter(k,p,_defaultvalue(p,eltype(p)))
parameter{V<:ValueType}(k::AbstractString,l::V,h::V,d::V) = parameter(k,_defaultprior(l,h),d)
parameter{V<:ValueType}(k::AbstractString,l::V,h::V) = parameter(k,_defaultprior(l,h))

###Factory functions for unnamed parameters with a prior distribution
parameter(p::Distributions.UnivariateDistribution,d::ValueType) = parameter(_defaultkey(),p,d)
parameter(p::Distributions.UnivariateDistribution) = parameter(p,_defaultvalue(p,eltype(p)))
parameter{V<:ValueType}(l::V,h::V,d::V) = parameter(_defaultprior(l,h),d)
parameter{V<:ValueType}(l::V,h::V) = parameter(_defaultprior(l,h))

###Vectorised factory functions
parameters{N<:Number}(k::AbstractVector,v::AbstractVector{N}) = (@assert length(k) == length(v) ; map(parameter,k,v))
parameters{N<:Number}(k::AbstractVector,::Type{N}) = [parameter(ki,N) for ki in k]
parameters{N<:Number}(v::AbstractVector{N}) = map(parameter,v)
parameters{N<:Number}(n::Integer,::Type{N}) = [parameter(N) for i=1:n]

parameters{P<:Distributions.Distribution,V<:ValueType}(k::AbstractVector,p::AbstractVector{P},v::AbstractVector{V}) = (@assert length(k) == length(p) == length(v) ; map(parameter,k,p,v))
parameters{P<:Distributions.Distribution}(k::AbstractVector,p::AbstractVector{P}) = (@assert length(k) == length(p) ; map(parameter,k,p))
parameters{V<:ValueType}(k::AbstractVector,l::Vector{V},h::Vector{V},d::Vector{V}) = (@assert length(k) == length(l) == length(h) == length(d) ; map(parameter,k,l,h,d))
parameters{V<:ValueType}(k::AbstractVector,l::Vector{V},h::Vector{V}) = (@assert length(k) == length(l) == length(h) ; map(parameter,k,l,h))
parameters{P<:Distributions.Distribution,V<:ValueType}(p::AbstractVector{P},v::AbstractVector{V}) = (@assert length(p) == length(v) ; map(parameter,p,v))
parameters{P<:Distributions.Distribution}(p::AbstractVector{P}) = map(parameter,p)
parameters{V<:ValueType}(l::Vector{V},h::Vector{V},d::Vector{V}) = (@assert length(l) == length(h) == length(d) ; map(parameter,l,h,d))
parameters{V<:ValueType}(l::Vector{V},h::Vector{V}) = (@assert length(l) == length(h) ; map(parameter,l,h))

###Functionality to initialize MCMC from parameter definitions
#for individual parameters
@inline _initvalue{I<:InitializeFrom,N<:Number}(::Type{I},p::Parameter,::Type{N}) = N(p.default) #general case initializes to default value
@inline _initvalue{P<:ParameterUnivariate,N<:Number}(::Type{InitializeFromPrior},p::P,::Type{N}) = N(rand(p.prior)) #special case for initializing from prior

#vectorized internal function
@inline _initvalues!{I<:InitializeFrom,P<:Parameter,N<:Number}(::Type{I},p::AbstractVector{P},v::AbstractVector{N}) = (@simd for i=1:length(p) @inbounds v[i] = _initvalue(V,p[i],N) end ; v)

#exported initvalue functions
initvalues!{I<:InitializeFrom,P<:Parameter,N<:Number}(::Type{I},p::AbstractVector{P},v::AbstractVector{N}) = (@assert length(p) == length(v) ; _initvalues!(I,p,v))
initvalues{I<:InitializeFrom,P<:Parameter,N<:Number}(::Type{I},p::AbstractVector{P},::Type{N}) = _initvalues!(I,p,Vector{N}(length(p)))

###Functionality to calculate the logprior
@inline _logprior{N<:Number}(p::ParameterUnivariate,v::N) = N(Distributions.logpdf(p.prior,v))
@inline _logprior{N<:Number}(p::ParameterDefault,v::N) = zero(N)

###vectorized internal logprior function
@inline _logprior{P<:Parameter,N<:Number}(p::AbstractVector{P},v::AbstractVector{N}) = (s = zero(N) ; @simd for i=1:length(p) @inbounds s += _logprior(p[i],v[i]) end ; s)
@inline _logprior{P<:ParameterDefault,N<:Number}(p::AbstractVector{P},v::AbstractVector{N}) = zero(N)

###exported logprior function
logprior{P<:Parameter,N<:Number}(p::P,v::N) = _logprior(p,v)
logprior{P<:Parameter,N<:Number}(p::AbstractVector{P},v::AbstractVector{N}) = (@assert length(p) == length(v) ; _logprior(p,v))

###Overloaded functions and operators from Base package
import Base.==
function =={P<:Parameter}(p1::P,p2::P)
  for f in fieldnames(p1)
    isequal(getfield(p1,f),getfield(p2,f))?continue:(return false)
  end
  return true
end

function Base.show(io::IO,p::Parameter,n::AbstractString ="")
  println(io,n,typeof(p).name.name)
  for f in fieldnames(p)
    println(io," ",f,": ",getfield(p,f))
  end
  nothing
end

function Base.show{P<:Parameter}(io::IO,v::AbstractVector{P})
  e = isa(eltype(v),Union)?"Parameter":eltype(v).name.name
  println(io)
  println(io,"Array{$e} with")
  for i=1:length(v)
    show(io,v[i],"[$i] ")
  end
  println(io)
end

Base.display{P<:Parameter}(io::IO,v::AbstractVector{P}) = Base.show{P}(io,v)







