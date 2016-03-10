immutable TargetModel{T<:Number,P<:AbstractParameter,D<:AbstractData,N<:AbstractNoiseModel} <: AbstractModel

    #Generic model specs
    name::AbstractString
    parameters::Vector

    measurements::D
    noisemodel::N

    target::Function #target function
    args::Tuple #additional arguments for the target function

    ###temp location to store model data
    modifydata::Bool
    modeldata::AbstractArray

    TargetModel(name::AbstractString,parameters::Vector{P},measurements::D,noisemodel::N,
                modify::Bool,modeldata::AbstractArray{T},target::Function,args...) = new(name,parameters,measurements,noisemodel,target,args,modify,modeldata)
end

function TargetModel(name::AbstractString,parameters::Vector,measurements::AbstractData,noisemodel::AbstractNoiseModel,modify::Bool,target::Function,args...)
    T = eltype(measurements)
    P = eltype(parameters)
    D = typeof(measurements)
    N = typeof(noisemodel)
    TargetModel{T,P,D,N}(name,parameters,measurements,noisemodel,modify,similar(datavalues(measurements)),target,args...)
end

###Factory functions
function _model(::Type{Val{:target!}},parameters::Vector,measurements::AbstractData,noisemodel::AbstractNoiseModel,target::Function,args...;name="TARGET!")
    TargetModel(name,parameters,measurements,noisemodel,true,target,args...)
end

function _model(::Type{Val{:target}},parameters::Vector,measurements::AbstractData,noisemodel::AbstractNoiseModel,target::Function,args...;name="TARGET")
    TargetModel(name,parameters,measurements,noisemodel,false,target,args...)
end

@inline _evaluate(::Type{Val{true}},m::TargetModel,vals::AbstractVector) = m.target(m.modeldata,m.args...,vals)
@inline _evaluate(::Type{Val{false}},m::TargetModel,vals::AbstractVector) = copy!(m.modeldata,m.target(m.args...,vals))

###Evaluate the model for the given parameter values
evaluate!(m::TargetModel,vals::AbstractVector) = (_evaluate(Val{m.modifydata},m,vals) ; m.modeldata)

###Utility functions used in generic implementations in AbstractModel
@inline dataindex(m::TargetModel) = dataindex(m.measurements)
@inline measurements(m::TargetModel) = datavalues(m.measurements)
@inline noisemodel(m::TargetModel) = m.noisemodel

###Base.show function
function show(io::IO,m::TargetModel)
    println(io,"TargetModel ",m.name)
    print(io,"parameters: ") ; show(io,m.parameters)
    print(io,"measurements: ") ; show(io,m.measurements)
    print(io,"noisemodel: ") ; show(io,m.noisemodel)
    println(io,"target: ",m.target)
    if ~isempty(m.args)
        println(io,"Additional target function arguments: ")
        show(io,m.args)
        println(io)
    end
    nothing
end
