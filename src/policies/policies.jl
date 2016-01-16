#Policy types to generate initial parameter values
abstract InitializeFrom
immutable InitializeFromDefault <: InitializeFrom end #generate initial values from the defaults
immutable InitializeFromPrior <: InitializeFrom end #generate initial values from the priors

#Policy types to propose new samples
abstract ProposeFrom
immutable ProposeFromIndicator <: ProposeFrom end #used when proposing when number of proposals == 1 (standard Metropolis-Hastings)
immutable ProposeFromAuxiliary <: ProposeFrom end #auxiliary variable used when number of proposals > 1 (generalized Metropolis-Hastings)

#Policy types for generating the indicator matrix
abstract GenerateIndicator
immutable IndicatorStationary <: GenerateIndicator end #repmat the stationary distribution
immutable IndicatorCyclical <: GenerateIndicator end #construct a cyclical indicator matrix from the stationary distribution

#The default number type to run MCMC with
typealias DefaultNumberType Float64

#Abstract runtime policy
abstract RuntimePolicy{N<:Number}

#Local function mapping number of proposals to the ProposeFrom type
_num2proposefrom(n::Integer) = (@assert n > 0 ; n>1?ProposeFromAuxiliary:ProposeFromIndicator)

#Factory function, currently supported type is :generic
policy(s::Symbol,i::DataType,g::DataType,nproposals::Integer) = _policy(Val{s},i,_num2proposefrom(nproposals),g,DefaultNumberType)
policy(s::Symbol,i::DataType,g::DataType,nproposals::Integer,n::DataType) = _policy(Val{s},i,_num2proposefrom(nproposals),g,n)
