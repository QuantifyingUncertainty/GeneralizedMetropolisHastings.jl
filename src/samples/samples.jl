#Abstract type definitions
abstract DerivativeOrder
abstract AbstractSample{O<:DerivativeOrder}

#Concrete derivated order types
type ZeroOrder <: DerivativeOrder end
type FirstOrder <: DerivativeOrder end
type SecondOrder <: DerivativeOrder end
type ThirdOrder <: DerivativeOrder end
typealias GradientOrder Union{FirstOrder,SecondOrder,ThirdOrder}
typealias TensorOrder Union{SecondOrder,ThirdOrder}
typealias DTensorOrder ThirdOrder

#Factory functions for all different types of samples
#Currently implemented as symbols are: :base, :gradient, :tensor, :tangent
samples{N<:Number,T<:AbstractFloat}(s::Symbol,nparams::Integer,nsamples::Integer,::Type{N},::Type{T}) = _samples(Val{s},nparams,nsamples,N,T)
samples{N<:Number,T<:AbstractFloat}(s::Symbol,nparams::Integer,nsamples::Integer,ntangents::Integer,::Type{N},::Type{T}) = _samples(Val{s},nparams,nsamples,ntangents,N,T)

### Functionality
numparas(s::AbstractSample) = size(s.values,1)
numsamples(s::AbstractSample) = size(s.values,2)
eltype(s::AbstractSample) = eltype(s.values)

function copy!{S<:AbstractSample}(dest::S,src::S)
    for f in fieldnames(dest)
        copy!(getfield(dest,f),getfield(src,f))
    end
end

function copy!{S<:AbstractSample}(dest::S,src::S,sindex::Integer)
    @assert numsamples(dest) == 1
    for f in fieldnames(dest)
        destvals = getfield(dest,f)
        l = length(destvals)
        copy!(destvals,1,getfield(src,f),(sindex-1)*l+1,l)
    end
end

function =={S<:AbstractSample}(s1::S,s2::S)
    for f in fieldnames(s1)
        isequal(getfield(s1,f),getfield(s2,f))?continue:return false
    end
    return true
end

function show{S<:AbstractSample}(io::IO,s::S,n::AbstractString ="")
    nump = numparas(s)
    nums = numsamples(s)
    println(io,n,typeof(s).name.name," with ",nump," parameter",nump>1?"s":""," and ",nums," sample",nums>1?"s":"")
    for f in fieldnames(s)
        print(io," ",f,": ",typeof(getfield(s,f)))
        println(io)
    end
    nothing
end

function show{S<:AbstractSample}(io::IO,v::AbstractVector{S})
    e = isa(eltype(v),Union)?"AbstractSample":eltype(v).name.name
    println(io)
    println(io,"Array{$e} with")
    for i=1:length(v)
        show(io,v[i],"[$i] ")
    end
end

display{S<:AbstractSample}(io::IO,v::AbstractVector{S}) = Base.show{S}(io,v)

######################################################
@inline _valuestuple(s1::Integer,s2::Integer) = (s2>1?(s1,s2):(s1,)) #if values is one-dimensional, then tensor can be two-dimensional
@inline _tensortuple(::Type{Val{1}},s1::Integer,s2::Integer,s3::Integer) = (s1,s2) #if values is one-dimensional, then tensor can be two-dimensional
@inline _tensortuple(::Type{Val{2}},s1::Integer,s2::Integer,s3::Integer) = (s1,s2,s3) #if values is two-dimensional, then tensor should be three-dimensional



