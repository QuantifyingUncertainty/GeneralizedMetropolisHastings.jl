typealias Prior UnivariateDistribution

immutable ModelParameters{P<:Prior,S<:String}
  names::Vector{S}
  priors::Vector{P}
  defaults::Vector{Float64}
  function ModelParameters(n::Vector{S},p::Vector{P},d::Vector{Float64})
    @assert length(n) == length(p) == length(d) "Number of names, priors and default values not equal"
    new(n,p,d)
  end
end

#Constructors
ModelParameters{P<:Prior,S<:String}(d::Vector{Float64},p::Vector{P}; names::Vector{S} =String[],index::Vector{Int} =Int[]) =
  ModelParameters{P,S}(expand_names(length(d),names,index),p,d)

ModelParameters{P<:Prior,S<:String}(d::Vector{Float64},p::P =Uniform(-1e43,1e43); names::Vector{S} =String[],index::Vector{Int} =Int[]) =
  ModelParameters{P,S}(expand_names(length(d),names,index),fill(p,length(d)),d)

ModelParameters{P<:Prior,S<:String}(p::Vector{P}; names::Vector{S} =String[],index::Vector{Int} =Int[]) =
  ModelParameters{P,S}(expand_names(length(p),names,index),p,fill(0.0,length(p)))

ModelParameters{P<:Prior,S<:String}(numparas::Int =0,p::P =Uniform(-1e43,1e43); names::Vector{S} =String[],index::Vector{Int} =Int[]) =
  ModelParameters{P,S}(expand_names(numparas,names,index),fill(p,numparas),fill(0.0,numparas))

#utility functions
values(::ValuesFromDefault,p::ModelParameters) = copy(p.defaults)
values(::ValuesFromPrior,p::ModelParameters) = map((x)->rand(x),p.priors)
named(p::ModelParameters) = find(p.names .!= "")
anonymous(p::ModelParameters) = find(p.names .== "")
numel(p::ModelParameters) = length(p.names)

#overloaded functions and operators from Base package
==(x::ModelParameters,y::ModelParameters) = isequal(x.names,y.names) && isequal(x.defaults,y.defaults) && isequal(x.priors,y.priors)

function Base.show(io::IO,p::ModelParameters)
  println(io,"ModelParameters with $(numel(p)) elements")
  for i = 1:numel(p)
    println(io,"  name: \"",p.names[i],"\", default: ",p.defaults[i]," prior: ",p.priors[i])
  end
  nothing
end

#Helper functions that are not needed outside the module
function expand_names{S<:String}(numparas::Int,names::Vector{S},index::Vector{Int} =Int[])
  @assert length(names) == length(index) || isempty(index) && numparas == length(names) "Number of names must equal number of index elements or equal total number of variables"
  @assert isempty(index) || maximum(index) <= numparas "Maximum index of named parameters is larger than total number of parameters"
  if isempty(index) && length(names) == numparas
    index = [1:numparas]
  end
  nc = 0
  S[in(j,index)?names[nc+=1]:"" for j=1:numparas]
end


