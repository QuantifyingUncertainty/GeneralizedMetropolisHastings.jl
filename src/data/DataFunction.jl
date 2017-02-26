"""
    DataFunction

A type wrapping a data generating function.
"""
immutable DataFunction{T<:Number,D<:Number,A<:AbstractArray} <: AbstractData
    index::Vector
    values::A
    modifyvalues::Bool
    generator::Function
    args::Tuple
    DataFunction(i::AbstractVector{T},d::AbstractArray{D},w::Bool,g::Function,args...) = new(collect(i),d,w,g,args)
end

"""
    data(:function!,index::AbstractVector,out::AbstractArray,f!::Function,[,...args])

Create a `DataFunction` object which encapsulates an in-place data generating function.

The `:function!` variant expects a pre-allocated output array `out`,
and a data generating function which takes `out` as its first argument.

# Arguments
* `index::AbstractVector`: an index for the rows of the result (e.g., a time index)
* `out::AbstractArray`: array with one row per data point and one column per variable
* `f!::Function`: a data generating function taking as first argument `out`
* `args...`: additional arguments of the data generating function
---
"""
function data(::Type{Val{:function!}},i::AbstractVector,d::AbstractArray,g::Function,args...)
    @assert length(i) == size(d,1)
    g(d,args...)
    DataFunction{eltype(i),eltype(d),typeof(d)}(i,d,true,g,args...)
end

"""
    data(:function,index::AbstractVector,f::Function[,args...])

Create a `DataFunction` object which encapsulates a data generating function.

The `:function` variant can be used if the data generating function
cannot generate its result in a pre-allocated output array.
The result should contain one row per data point and one column per variable.

# Arguments
* `index::AbstractVector`: an index for the rows of the result (e.g., a time index)
* `f::Function`: a function returning an array of values as a result
* `args...`: additional arguments of the data generating function
---
"""
function data(::Type{Val{:function}},i::AbstractVector,g::Function,args...)
    d = g(args...)
    @assert length(i) == size(d,1)
    DataFunction{eltype(i),eltype(d),typeof(d)}(i,d,false,g,args...)
end

@inline _generate(::Type{Val{true}},d::DataFunction) = d.generator(d.values,d.args...)
@inline _generate(::Type{Val{false}},d::DataFunction) = copy!(d.values,d.generator(d.args...))

numvalues(d::DataFunction) = size(d.values,1)
numvars(d::DataFunction) = size(d.values,2)
eltype(d::DataFunction) = eltype(d.values)
generate!(d::DataFunction) = (_generate(Val{d.modifyvalues},d) ; d)
dataindex(d::DataFunction) = d.index
datavalues(d::DataFunction) = d.values
datatypename(d::DataFunction) = "DataFunction"
