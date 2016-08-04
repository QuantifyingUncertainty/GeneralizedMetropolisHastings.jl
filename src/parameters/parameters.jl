### Paramater base type
abstract AbstractParameter

parameters(s::Symbol,args...;keyargs...) = _parameters(Val{s},args...;keyargs...)
_parameters(::Type{Val{:nothing}},args...;keyargs...) = nothing #define at least one _parameters function as hook

###Named MCMC parameter which has a default value but no prior
immutable ParameterDefault{N<:Number} <: AbstractParameter
  key::Symbol
  default::N
  ParameterDefault(key::Symbol,def::N) = new(key,def)
end

###Named MCMC parameter which has a prior distribution and/or a default value
immutable ParameterUnivariate{U<:Distributions.UnivariateDistribution,N<:Number} <: AbstractParameter
  key::Symbol
  prior::U
  default::N
  ParameterUnivariate(key::Symbol,pri::U,def::N) = (@assert Distributions.insupport(pri,def) ; new(key,pri,def))
end

###Local, non-exported helper functions defining default behaviour
@inline _defaultkey() = :param
@inline _defaultvalue{N<:Number}(::Type{N}) = zero(N)
@inline _defaultvalue{P<:Distributions.ContinuousDistribution}(p::P) = mean(p)
@inline _defaultvalue{P<:Distributions.DiscreteDistribution}(p::P) = trunc(eltype(Distributions.Discrete),mean(p))
@inline _defaultunivariate(l::eltype(Distributions.Continuous),h::eltype(Distributions.Continuous)) = Distributions.Uniform(l,h)
@inline _defaultunivariate(l::eltype(Distributions.Discrete),h::eltype(Distributions.Discrete)) = Distributions.DiscreteUniform(l,h)

###Factory functions for named parameters with only a default value
parameter(key::Symbol,default::Number) = ParameterDefault{typeof(default)}(key,default)
parameter{N<:Number}(key::Symbol,::Type{N}) = parameter(key,_defaultvalue(N))

###Factory functions for unnamed parameters with only a default value
parameter(default::Number) = parameter(_defaultkey(),default)
parameter{N<:Number}(::Type{N}) = parameter(_defaultvalue(N))

###Factory functions for named parameters with a prior distribution
parameter(key::Symbol,prior::Distributions.UnivariateDistribution,default::Number) = ParameterUnivariate{typeof(prior),eltype(prior)}(key,prior,default)
parameter(key::Symbol,prior::Distributions.UnivariateDistribution) = parameter(key,prior,_defaultvalue(prior))
parameter{N<:Number}(key::Symbol,low::N,high::N,default::N) = parameter(key,_defaultunivariate(low,high),default)
parameter{N<:Number}(key::Symbol,low::N,high::N) = parameter(key,_defaultunivariate(low,high))

###Factory functions for unnamed parameters with a prior distribution
parameter(prior::Distributions.UnivariateDistribution,default::Number) = parameter(_defaultkey(),prior,default)
parameter(prior::Distributions.UnivariateDistribution) = parameter(prior,_defaultvalue(prior))
parameter{N<:Number}(l::N,h::N,d::N) = parameter(_defaultunivariate(l,h),d)
parameter{N<:Number}(l::N,h::N) = parameter(_defaultunivariate(l,h))

###Vectorised factory functions
parameters{N<:Number}(keys::Vector{Symbol},defaults::Vector{N}) = (@assert length(keys) == length(defaults) ; map(parameter,keys,defaults))
parameters{N<:Number}(keys::Vector{Symbol},::Type{N}) = [parameter(key,N) for key in keys]
parameters{N<:Number}(defaults::Vector{N}) = map(parameter,defaults)
parameters{N<:Number}(n::Integer,::Type{N}) = [parameter(N) for i=1:n]

parameters{P<:Distributions.Distribution,N<:Number}(keys::Vector{Symbol},priors::Vector{P},defaults::Vector{N}) =
  (@assert length(keys) == length(priors) == length(defaults) ; map(parameter,keys,priors,defaults))
parameters{P<:Distributions.Distribution}(keys::Vector{Symbol},priors::Vector{P}) =
  (@assert length(keys) == length(priors) ; map(parameter,keys,priors))
parameters{N<:Number}(keys::Vector{Symbol},lows::Vector{N},highs::Vector{N},defaults::Vector{N}) =
  (@assert length(keys) == length(lows) == length(highs) == length(defaults) ; map(parameter,keys,lows,highs,defaults))
parameters{N<:Number}(keys::Vector{Symbol},lows::Vector{N},highs::Vector{N}) =
  (@assert length(keys) == length(lows) == length(highs) ; map(parameter,keys,lows,highs))
parameters{P<:Distributions.Distribution,N<:Number}(priors::Vector{P},defaults::Vector{N}) =
  (@assert length(priors) == length(defaults) ; map(parameter,priors,defaults))
parameters{P<:Distributions.Distribution}(p::Vector{P}) = map(parameter,p)
parameters{N<:Number}(lows::Vector{N},highs::Vector{N},defaults::Vector{N}) =
  (@assert length(lows) == length(highs) == length(defaults) ; map(parameter,lows,highs,defaults))
parameters{N<:Number}(lows::Vector{N},highs::Vector{N}) =
  (@assert length(lows) == length(highs) ; map(parameter,lows,highs))

###Main functionality

#to initialize values from individual parameter definitions
@inline _initializetodefault{N<:Number}(p::AbstractParameter,::Type{N}) = N(p.default)
@inline _initvalue{N<:Number}(::Type{Val{:default}},p::AbstractParameter,::Type{N}) = _initializetodefault(p,N)
@inline _initvalue{N<:Number}(::Type{Val{:prior}},p::ParameterDefault,::Type{N}) = _initializetodefault(p,N) #can't initialize a ParameterDefault to a prior value
@inline _initvalue{N<:Number}(::Type{Val{:prior}},p::ParameterUnivariate,::Type{N}) = N(rand(p.prior))

#vectorized internal function
@inline _initvalues!{I<:Union{Val{:default},Val{:prior}},P<:AbstractParameter,N<:Number}(::Type{I},p::Vector{P},v::AbstractVector{N}) = (@simd for i=1:length(p) @inbounds v[i] = _initvalue(I,p[i],N) end ; v)

#exported initvalues functions
initvalues!{P<:AbstractParameter,N<:Number}(i::InitializeFrom,p::Vector{P},v::AbstractVector{N}) = (@assert length(p) == length(v) ; _initvalues!(traittype(i),p,v))
initvalues{P<:AbstractParameter,N<:Number}(i::InitializeFrom,p::Vector{P},::Type{N}) = _initvalues!(traittype(i),p,Vector{N}(length(p)))

###Functionality to calculate the logprior
@inline _logprior{T<:AbstractFloat}(p::ParameterUnivariate,v::Number,::Type{T}) = T(Distributions.logpdf(p.prior,v))
@inline _logprior{T<:AbstractFloat}(p::ParameterDefault,v::Number,::Type{T}) = zero(T)

###vectorized internal logprior function
@inline function _logprior{P<:AbstractParameter,T<:AbstractFloat}(p::Vector{P},v::AbstractVector,::Type{T})
    s = zero(T)
    for i=1:length(p)
        @inbounds s += _logprior(p[i],v[i],T)
        isfinite(s)?continue:break
    end
    s
end

@inline function _logprior!{P<:AbstractParameter,T<:AbstractFloat}(r::AbstractVector{T},p::AbstractVector{P},m::AbstractArray)
    for j=1:length(r)
        s = zero(T)
        for i=1:length(p)
            @inbounds s += _logprior(p[i],m[i,j],T)
            isfinite(s)?continue:break
        end
        @inbounds r[j] = s
    end
    r
end

@inline _logprior{P<:ParameterDefault,T<:AbstractFloat}(p::AbstractVector{P},v::AbstractVector,::Type{T}) = zero(T)
@inline _logprior!{P<:ParameterDefault,T<:AbstractFloat}(r::AbstractVector{T},p::AbstractVector{P},v::AbstractArray) = (@simd for i=1:length(r) @inbounds r[i] = zero(T) end ; r)

###exported logprior function
logprior!{P<:AbstractParameter,T<:AbstractFloat}(r::AbstractVector{T},p::AbstractVector{P},m::AbstractArray) = _logprior!(r,p,m)
logprior{P<:AbstractParameter,T<:AbstractFloat}(p::P,v::Number,::Type{T}) = _logprior(p,v,T)
logprior{P<:AbstractParameter,T<:AbstractFloat}(p::AbstractVector{P},v::AbstractVector,::Type{T}) = _logprior(p,v,T)
logprior{P<:AbstractParameter,T<:AbstractFloat}(p::AbstractVector{P},m::AbstractMatrix,::Type{T}) = _logprior!(zeros(T,size(m,2)),p,m)

###Overloaded functions and operators from Base package
function =={P<:AbstractParameter}(p1::P,p2::P)
  for f in fieldnames(p1)
    isequal(getfield(p1,f),getfield(p2,f))?continue:(return false)
  end
  return true
end

function show(io::IO,p::AbstractParameter,n::AbstractString ="")
  println(io,n,typeof(p).name.name)
  for f in fieldnames(p)
    println(io," ",f,": ",getfield(p,f))
  end
  nothing
end

function show{P<:AbstractParameter}(io::IO,v::AbstractVector{P})
  e = isa(eltype(v),Union)?"AbstractParameter":eltype(v).name.name
  println(io,"Array{$e} with")
  for i=1:length(v)
    show(io,v[i],"[$i] ")
  end
end

display{P<:AbstractParameter}(io::IO,v::AbstractVector{P}) = Base.show{P}(io,v)







