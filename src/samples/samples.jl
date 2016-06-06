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
samples{N<:Number,T<:AbstractFloat}(s::Symbol,nparams::Integer,nsamples::Integer,::Type{N},::Type{T},args...) = _samples(Val{s},nparams,nsamples,N,T,args...)

### Functionality
numparas(s::AbstractSample) = size(s.values,1)
numsamples(s::AbstractSample) = size(s.values,2)
sampletype(s::AbstractSample) = eltype(s.values)
calculationtype(s::AbstractSample) = eltype(s.loglikelihood)
similar(s::AbstractSample,nsamples::Integer =numsamples(s),tosamplename::Symbol =sampletypename(s)) = _samples(Val{tosamplename},numparas(s),nsamples,sampletype(s),calculationtype(s))

function _copy!(dest::AbstractSample,src::AbstractSample)
    for f in fieldnames(dest)
        copy!(getfield(dest,f),getfield(src,f))
    end
    dest
end

function _copy!(dest::AbstractSample,destindex::Integer,src::AbstractSample,srcindex::Integer)
    ndestsamples = numsamples(dest)
    for f in fieldnames(dest)
        destvals = getfield(dest,f)
        copylength = div(length(destvals),ndestsamples)
        copy!(destvals,(destindex-1)*copylength+1,getfield(src,f),(srcindex-1)*copylength+1,copylength)
    end
    dest
end

function _copy!(dest::AbstractSample,destindex::AbstractVector,src::AbstractSample,srcindex::AbstractVector)
    ndestsamples = numsamples(dest)
    ncopysamples = length(destindex)
    for f in fieldnames(dest)
        destvals = getfield(dest,f)
        copylength = ncopysamples>0?div(length(destvals),ndestsamples):0
        for i=1:ncopysamples
            copy!(destvals,(destindex[i]-1)*copylength+1,getfield(src,f),(srcindex[i]-1)*copylength+1,copylength)
        end
    end
    dest
end

function copy!(dest::AbstractSample,src::AbstractSample)
    @assert numparas(dest) == numparas(src) && numsamples(dest) == numsamples(src)
    _copy!(dest,src)
end

function copy!(dest::AbstractSample,destindex::Integer,src::AbstractSample,srcindex::Integer)
    @assert numparas(dest) == numparas(src) && destindex <= numsamples(dest) && srcindex <= numsamples(src)
    _copy!(dest,destindex,src,srcindex)
end

function copy!(dest::AbstractSample,destindex::AbstractVector,src::AbstractSample,srcindex::AbstractVector)
    @assert numparas(dest) == numparas(src) && maximum(destindex) <= numsamples(dest) && maximum(srcindex) <= numsamples(src)
    _copy!(dest,destindex,src,srcindex)
end

function copy!(dest::AbstractSample,src::AbstractSample,srcindex::AbstractVector)
    @assert numparas(dest) == numparas(src) && numsamples(dest) == length(srcindex) && maximum(srcindex) <= numsamples(src)
    _copy!(dest,1:length(srcindex),src,srcindex)
end

copy(src::AbstractSample) = _copy!(similar(src),src)
copy(src::AbstractSample,i::Integer) = _copy!(similar(src,1),1,src,i)
copy(src::AbstractSample,srcindex::AbstractVector) = _copy!(similar(src,length(srcindex)),1:length(srcindex),src,srcindex)
copy(src::AbstractSample,srcindex::AbstractVector,to::Symbol) = _copy!(similar(src,length(srcindex),to),1:length(srcindex),src,srcindex)

function offset!(s::AbstractSample,v::Vector)
    nparas = numparas(s)
    for j=1:numsamples(s)
        @simd for i=1:nparas
            @inbounds s.values[i,j] += v[i]
        end
    end
    s
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
@inline _valuestuple(s1::Integer,s2::Integer) = s2==1?(s1,):(s1,s2) #if s2 == 1, then values should be one-dimensional
@inline _tensortuple(s1::Integer,s2::Integer,s3::Integer) = s3==1?(s1,s2):(s1,s2,s3) #if s3 == 1, then tensor should be two-dimensional



