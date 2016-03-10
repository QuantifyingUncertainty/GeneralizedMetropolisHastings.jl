immutable DataFunction{T<:Number,D<:Number,A<:AbstractArray} <: AbstractData
    index::Vector
    values::A
    modifyvalues::Bool
    generator::Function
    args::Tuple
    DataFunction(i::AbstractVector{T},d::AbstractArray{D},w::Bool,g::Function,args...) = new(collect(i),d,w,g,args)
end

function _data(::Type{Val{:function!}},i::AbstractVector,d::AbstractArray,g::Function,args...)
    @assert length(i) == size(d,1)
    DataFunction{eltype(i),eltype(d),typeof(d)}(i,d,true,g,args...)
end

function _data(::Type{Val{:function}},i::AbstractVector,g::Function,args...)
    d = g(args...)::AbstractArray
    @assert length(i) == size(d,1)
    DataFunction{eltype(i),eltype(d),typeof(d)}(i,d,false,g,args...)
end

@inline _generate(::Type{Val{true}},d::DataFunction) = d.generator(d.values,d.args...)
@inline _generate(::Type{Val{false}},d::DataFunction) = copy!(d.values,d.generator(d.args...))

numvalues(d::DataFunction) = size(d.values,1)
numvars(d::DataFunction) = size(d.values,2)
eltype(d::DataFunction) = eltype(d.values)
generate!(d::DataFunction) = (_generate(Val{d.modifyvalues},d) ; d.values)
dataindex(d::DataFunction) = d.index
datavalues(d::DataFunction) = d.values
datatypename(d::DataFunction) = "DataFunction"
