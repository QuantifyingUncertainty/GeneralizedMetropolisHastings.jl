abstract AbstractModel

### Factory function
model(s::Symbol,args...;keyargs...) = _model(Val{s},args...;keyargs...)

### Initialize is the same for all models
initialize!(i::InitializeFrom,m::AbstractModel,s::AbstractSample) = (initvalues!(i,parameters(m),s.values) ; s)

### For each model, an evaluation function needs to be defined
evaluate!(m::AbstractModel,v::AbstractVector) = throw(MethodError(evaluate!, (m,v)))

### Common utility functions
numparas(m::AbstractModel) = length(m.parameters)
parameters(m::AbstractModel) = m.parameters

### Utility functions that need to be implemented for each model
@inline dataindex(m::AbstractModel) = throw(MethodError(dataindex, (m,)))
@inline measurements(m::AbstractModel) = throw(MethodError(measurements, (m,)))
@inline noisemodel(m::AbstractModel) = throw(MethodError(noisemodel, (m,)))








