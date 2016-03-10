immutable ODEModel{T<:AbstractFloat,P<:AbstractParameter,D<:AbstractData,N<:AbstractNoiseModel} <: AbstractModel
    #Generic model specs
    name::AbstractString
    parameters::Vector
    gradientepsilon::T

    #Data + noise
    measurements::D
    noisemodel::N

    #ODE specific
    ode::Function
    initial::Vector{T}
    numstates::Int
    observed::Vector{Int}
    unobserved::Vector{Int}
    abstol::T
    reltol::T

    #temp location to store model data
    modeldata::Array

    ODEModel(name::AbstractString,parameters::Vector{P},gradientepsilon::T,measurements::D,noisemodel::N,
             ode::Function,initial::Vector{T},numstates::Int,observed::Vector{Int},unobserved::Vector{Int},abstol::T,reltol::T) =
        new(name,parameters,gradientepsilon,measurements,noisemodel,ode,initial,numstates,observed,unobserved,abstol,reltol,similar(datavalues(measurements)))
end

###Factory function with default values for tolerance
function _model(::Type{Val{:ode}},parameters::Vector,data::AbstractData,noise::AbstractNoiseModel,
       ode::Function,initial::Vector,numstates::Int,observed::Vector{Int};name ="ODE",abstol =1e-6,reltol =1e-6,gradientepsilon =1e-4)
    T = eltype(data)
    P = eltype(parameters)
    D = typeof(data)
    N = typeof(noise)
    ODEModel{T,P,D,N}(name,parameters,convert(T,gradientepsilon),data,noise,ode,initial,numstates,observed,setdiff(1:numstates,observed),convert(T,abstol),convert(T,reltol))
end

###Evaluate the model for the given parameter values
function evaluate!(m::ODEModel,vals::AbstractVector)
    o(t,y,ydot) = m.ode(t,y,ydot,vals)
    copy!(m.modeldata,sub(Sundials.cvode(o,m.initial,dataindex(m.measurements);reltol=m.reltol,abstol=m.abstol),:,m.observed))
    m.modeldata
end

###Utility functions used in generic implementations in AbstractModel
@inline dataindex(m::ODEModel) = dataindex(m.measurements)
@inline measurements(m::ODEModel) = datavalues(m.measurements)
@inline noisemodel(m::ODEModel) = m.noisemodel

#gradlogprior!(m::ODEModel,s::Sample{}) = (s.gradlogprior = map((p,v)->((logpdf(p,v+m.gradientepsilon)-logpdf(p,v))/m.gradientepsilon),m.parameters.priors,s.values))

###Calculate the loglikelihood function

###Calculate the gradient of the loglikelihood
#gradienthelper(m::ODEModel,data::Matrix{Float64},sol::Matrix{Float64},grad::Array{Float64,3}) = vec(sum((grad.*(data-sol)./m.variance),(1,2))) #.* does automatic broadcast(), summing result along first 2 dimensions
#gradloglikelihood(m::ODEModel,sol::Matrix{Float64},grad::Array{Float64,3}) = gradienthelper(m,m.measurements,sol,grad)

###Calculate the metric tensor
#tensorvalue(m::ODEModel,grad::Array{Float64,3},i::Int,j::Int) = sum(grad[:,:,i].*grad[:,:,j]./m.variance)

###Generate pseudodata for the approximate metric tensor calculations
pseudodata(m::ODEModel,d::AbstractArray) = applynoise!(m.noise,d)
#tangentvector(m::ODEModel,sol::Matrix{Float64},grad::Array{Float64,3}) = gradienthelper(m,pseudodata(m,sol),sol,grad)

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
    println(io,"gradientepsilon: ",m.gradientepsilon)
    println(io,"abstol: ",m.abstol)
    println(io,"reltol: ",m.reltol)
    nothing
end
