"""
    DataArray

A type wrapping a predefined array of measurement data.
"""
immutable DataArray{T<:Number,D<:Number,A<:AbstractArray} <: AbstractData
    index::Vector
    values::A
    DataArray(i::AbstractVector{T},d::AbstractArray{D}) = new(collect(i),d)
end

"""
    data(:array,index::AbstractVector,values::AbstractArray)

Create a data object wrapping an array of fixed, predefined data values.

# Arguments
* `index::AbstractVector`: an index for each row of the values array (e.g., a time index)
* `values::AbstractArray`: each row is a data point and each column a variable
---
"""
function data(::Type{Val{:array}},i::AbstractVector,d::AbstractArray)
    @assert length(i) == size(d,1)
    DataArray{eltype(i),eltype(d),typeof(d)}(i,d)
end

numvalues(d::DataArray) = size(d.values,1)
numvars(d::DataArray) = size(d.values,2)
eltype(d::DataArray) = eltype(d.values)
generate!(d::DataArray) = d
dataindex(d::DataArray) = d.index
datavalues(d::DataArray) = d.values
datatypename(d::DataArray) = "DataArray"
