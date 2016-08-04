immutable ODEModel{P<:AbstractParameter,D<:AbstractData,N<:AbstractNoiseModel} <: AbstractModel
    #Generic model specs
    name::AbstractString
    parameters::Vector

    #Data + noise
    measurements::D
    noisemodel::N

    #ODE specific
    ode::Function
    initial::Vector
    numstates::Int
    observed::Vector{Int}
    unobserved::Vector{Int}
    abstol::Real
    reltol::Real

    ODEModel(name::AbstractString,parameters::Vector{P},measurements::D,noisemodel::N,
             ode::Function,initial::Vector,numstates::Int,observed::Vector{Int},unobserved::Vector{Int},abstol::Real,reltol::Real) =
        new(name,parameters,measurements,noisemodel,ode,initial,numstates,observed,unobserved,abstol,reltol)
end

###Factory function with default values for tolerance
function _model(::Type{Val{:ode}},parameters::Vector,data::AbstractData,noise::AbstractNoiseModel,
       ode::Function,initial::Vector,numstates::Int,observed::Vector{Int};name ="ODE",abstol =1e-6,reltol =1e-6)
    P = eltype(parameters)
    D = typeof(data)
    N = typeof(noise)
    ODEModel{P,D,N}(name,parameters,data,noise,ode,initial,numstates,observed,setdiff(1:numstates,observed),abstol,reltol)
end

###Evaluate the model for the given parameter values
function evaluate!(m::ODEModel,vals::AbstractVector)
    o(t,y,ydot) = m.ode(t,y,ydot,vals)
    sub(Sundials.cvode(o,m.initial,dataindex(m.measurements);reltol=m.reltol,abstol=m.abstol),:,m.observed)
end

###Utility functions used in generic implementations in AbstractModel
@inline dataindex(m::ODEModel) = dataindex(m.measurements)
@inline measurements(m::ODEModel) = datavalues(m.measurements)
@inline noisemodel(m::ODEModel) = m.noisemodel

###Base functionality
function show(io::IO,m::ODEModel)
    println(io,"ODEModel ",m.name)
    print(io,"parameters: ") ; show(io,m.parameters)
    print(io,"measurements: ") ; show(io,m.measurements)
    print(io,"noisemodel: ") ; show(io,m.noisemodel)
    println(io,"ode: ",m.ode)
    println(io,"initial: ",m.initial)
    println(io,"numstates: ",m.numstates)
    println(io,"observed: ",m.observed)
    if ~isempty(m.unobserved)
        println(io,"unobserved: ",m.unobserved)
    end
    println(io,"abstol: ",m.abstol)
    println(io,"reltol: ",m.reltol)
    nothing
end
