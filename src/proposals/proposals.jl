abstract AbstractProposalDensity

abstract SymmetricDensity <: AbstractProposalDensity
abstract ASymmetricDensity <: AbstractProposalDensity

issymmetric{D<:SymmetricDensity}(::Union{D,Type{D}}) = true
issymmetric{D<:ASymmetricDensity}(::Union{D,Type{D}}) = false

###Factory function for known distributions (currently :normal and :lognormal)
density(s::Symbol,args...) = _density(Val{s},args...)

###Factory function to wrap arbitrary distributions from the distributions package
density(d::Distributions.Distribution; symmetric::Bool =false) = _density(Val{symmetric},d)

###For some densities, conditioning on a point can be faster than creating the density
###But in the general case, we don't know how to condition a density so throw an error
condition!(d::AbstractProposalDensity,x::AbstractArray,args...) = throw(MethodError(condition!, (d,x,args...)))

###For generic proposal densities, we don't know how to propose points, so throw an error
propose!(d::AbstractProposalDensity,s::AbstractSample) = throw(MethodError(propose!, (d,s)))

###For generic proposal densities, we don't know how to calculate the logprobability of a sample, so throw an error
logprobability!(r::AbstractVector,d::AbstractProposalDensity,v::AbstractArray) = throw(MethodError(logprobability!, (r,d,v)))
logprobability(d::AbstractProposalDensity,v::AbstractArray) = throw(MethodError(logprobability, (d,v)))

function show(io::IO,d::AbstractProposalDensity)
    println(io,"$(typeof(d)) with fields: $(fieldnames(d))")
    nothing
end
