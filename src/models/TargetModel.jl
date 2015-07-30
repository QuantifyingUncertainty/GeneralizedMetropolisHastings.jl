immutable TargetModel <: MCModel

  name::String
  parameters::ModelParameters

  target::Function #target function
  args::Dict #additional arguments for the target function
end

#TargetModel(t::Function;name::String =String[],parameters::ModelParameters =ModelParameters(0),args...) = TargetModel(name,parameters,t,args)

function Base.show(io::IO,t::TargetModel)
  println(io,"TargetModel ",t.name)
  show(io,t.parameters)
  println(io,"target: ",t.target)
  println(io,"args: ",t.args)
  nothing
end
