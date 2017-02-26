"""
Base type for objects wrapping measurement data.
"""
abstract AbstractData

"""
    data(s::Symbol,args...)

Create data objects.

Currently implemented variants:
---

"""
data(s::Symbol,args...) = data(Val{s},args...)

"""
    datatypename(d)
Return the name of the data type (used by `show()`).
"""
function datatypename(d::AbstractData) end

"""
    numvalues(d)
Return the number of data values.
"""
function numvalues(d::AbstractData) end

"""
    numvars(d)
Return the number of variables.
"""
function numvars(d::AbstractData) end

"""
    eltype(d)
Return the type of the data values.
"""
function eltype(d::AbstractData) end

"""
    generate!(d)
Generate a set of data values. The function should return its argument d.
"""
function generate!(d::AbstractData) end

"""
    dataindex(d)
Retrieve the data index.
"""
function dataindex(d::AbstractData) end

"""
    datavalues(d)
Retrieve the data values.
"""
function datavalues(d::AbstractData) end

function show(io::IO,d::AbstractData)
    nvals = (n = numvalues(d) ; n>1?string(n," values"):string(n," value"))
    nvars = (n = numvars(d) ; n>1?string(n," variables"):string(n," variable"))
    println(io,"$(datatypename(d)) with $(nvars) and $(nvals)")
    println(io," index: ",dataindex(d))
    println(io," values: ",datavalues(d))
    nothing
end
