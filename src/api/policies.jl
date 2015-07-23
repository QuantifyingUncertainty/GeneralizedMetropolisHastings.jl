#Policy types to generate initial parameter values
abstract ValuesFrom
immutable ValuesFromDefault <: ValuesFrom end #generate initial values from the defaults
immutable ValuesFromPrior <: ValuesFrom end #generate initial values from the priors

#Policy types to propose new samples
abstract ProposalFunction
immutable ProposalFromIndicator <: ProposalFunction end #used when proposing when number of proposals == 1 (standard Metropolis-Hastings)
immutable ProposalFromAuxiliary <: ProposalFunction end #auxiliary variable used when number of proposals > 1 (generalized Metropolis-Hastings)

#Policy types for generating the indicator matrix
abstract IndicatorMatrixFunction
immutable IndicatorMatrixStationary <: IndicatorMatrixFunction end #repmat the stationary distribution
immutable IndicatorMatrixOptimal <: IndicatorMatrixFunction end #construct an optimal MC matrix from the stationary distribution

#Abstract runtime policy
abstract RuntimePolicy

#Generic runtime policy
immutable GenericPolicy <: RuntimePolicy
  initialize::ValuesFrom
  indicate::IndicatorMatrixFunction
  propose::ProposalFunction
  function GenericPolicy(v::ValuesFrom,i::IndicatorMatrixFunction,np::Int)
    @assert np > 0 "Number of proposals should be > 0"
    new(v,i,np==1?ProposalFromIndicator():ProposalFromAuxiliary())
  end
end

GenericPolicy(v::ValuesFrom,np::Int) = GenericPolicy(v,IndicatorMatrixStationary(),np)
GenericPolicy(np::Int) = GenericPolicy(ValuesFromDefault(),IndicatorMatrixStationary(),np)

function Base.show(io::IO,p::GenericPolicy)
  println(io,"GenericPolicy with following policy types:")
  println(io,"initialize = ",typeof(p.initialize))
  println(io,"indicate = ",typeof(p.indicate))
  println(io,"propose = ",typeof(p.propose))
  println(io)
  nothing
end
