abstract DerivateOrder
type ZeroOrder <: DerivateOrder end
type FirstOrder <: DerivateOrder end
type SecondOrder <: DerivateOrder end
type ThirdOrder <: DerivateOrder end
typealias GradientOrder Union{FirstOrder,SecondOrder,ThirdOrder}
typealias TensorOrder Union{SecondOrder,ThirdOrder}
typealias DTensorOrder ThirdOrder

abstract Sample{O<:DerivativeOrder}

#Factory functions for samples of different derivative orders.
#Currently implemented as symbols are: :base, :gradient, :tensor, :tangent
samples{T<:Number}(s::Symbol,nparas::Integer,nsamples::Integer,::Type{T}) = _samples(Val{s},nparas,nsamples,T)
samples{T<:Number}(s::Symbol,nparas::Integer,nsamples::Integer,ntangents::Integer,::Type{T}) = _samples(Val{s},nparas,nsamples,ntangents,T)

### Functionality
numparas(s::Sample) = size(s.values,1)
numsamples(s::Sample) = size(s.values,2)

### Access functions defined for all sample types
@inline getvalue{S<:Sample}(s::S,para::Integer,sample::Integer) = s.values[para,sample]
@inline setvalue!{S<:Sample}(s::S,para::Integer,sample::Integer,val::Number) = (s.values[para,sample] = val)
@inline getloglikelihood{S<:Sample}(s::S,sample::Integer) = s.loglikelihood[sample]
@inline setloglikelihood!{S<:Sample}(s::S,sample::Integer,val::Number) = (s.loglikelihood[sample] = val)
@inline getlogprior{S<:Sample}(s::S,sample::Integer) = s.logprior[sample]
@inline setlogprior!{S<:Sample}(s::S,sample::Integer,val::Number) = (s.logprior[sample] = val)

### Access functions defined for first order sample types and above
@inline getgradloglikelihood{O<:GradientOrder}(s::Sample{O},para::Integer,sample::Integer) = s.gradloglikelihood[para,sample]
@inline setgradloglikelihood!{O<:GradientOrder}(s::Sample{O},para::Integer,sample::Integer,val::Number) = (s.gradloglikelihood[para,sample] = val)
@inline getgradlogprior{O<:GradientOrder}(s::Sample{O},para::Integer,sample::Integer) = s.gradlogprior[para,sample]
@inline setgradlogprior!{O<:GradientOrder}(s::Sample{O},para::Integer,sample::Integer,val::Number) = (s.gradlogprior[para,sample] = val)

### Access functiond defined for second order sample types and above
@inline gettensorloglikelihood{O<:TensorOrder}(s::Sample{O},row::Integer,col::Integer,sample::Integer) = s.tensorloglikelihood[row,col,sample]
@inline settensorloglikelihood!{O<:TensorOrder}(s::Sample{O},row::Integer,col::Integer,sample::Integer,val::Number) = (s.tensorloglikelihood[row,col,sample] = val)
@inline gettensorlogprior{O<:TensorOrder}(s::Sample{O},row::Integer,col::Integer,sample::Integer) = s.tensorlogprior[row,col,sample]
@inline settensorlogprior!{O<:TensorOrder}(s::Sample{O},row::Integer,col::Integer,sample::Integer,val::Number) = (s.tensorlogprior[row,col,sample] = val)

### Functionality from Base package
function =={S<:Sample}(s1::S,s2::S)
  for f in fieldnames(s1)
    isequal(getfield(s1,f),getfield(s2,f))?continue:(return false)
  end
  return true
end

function Base.show{S<:Sample}(io::IO,s::S,n::AbstractString ="")
  println(io,n,typeof(s).name.name," with $numparas(s) parameters and $numsamples(s) samples")
  for f in fieldnames(s)
    print(io," ",f,": ")
    ndims(getfield(s,f))>1?(println(io);show(io,getfield(s,f))):print(io,getfield(s,f))
    println(io)
  end
  nothing
end

function Base.show{S<:Sample}(io::IO,v::AbstractVector{S})
  e = isa(eltype(v),Union)?"Sample":eltype(v).name.name
  println(io)
  println(io,"Array{$e} with")
  for i=1:length(v)
    show(io,v[i],"[$i] ")
  end
end

Base.display{S<:Sample}(io::IO,v::AbstractVector{S}) = Base.show{S}(io,v)



