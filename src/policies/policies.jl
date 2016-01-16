#Policy types to generate initial parameter values
abstract InitializeFrom
immutable InitializeFromDefault <: IntializeFrom end #generate initial values from the defaults
immutable InitializeFromPrior <: InitializeFrom end #generate initial values from the priors

#Policy types to propose new samples
abstract ProposeFrom
immutable ProposeFromIndicator <: ProposeFrom end #used when proposing when number of proposals == 1 (standard Metropolis-Hastings)
immutable ProposeFromAuxiliary <: ProposeFrom end #auxiliary variable used when number of proposals > 1 (generalized Metropolis-Hastings)

#Policy types for generating the indicator matrix
abstract GenerateIndicator
immutable IndicatorStationary <: GenerateIndicator end #repmat the stationary distribution
immutable IndicatorCyclical <: GenerateIndicator end #construct a cyclical indicator matrix from the stationary distribution

#Abstract runtime policy
abstract RuntimePolicy{T<:Number}

#Factory function, currently supported type is :generic
policy(s::Symbol,i::InitializeFrom,g::GenerateIndicator,nproposals::Integer) = _policy(Val{s},i,g,nproposals,DefaultNumberType)
policy{T<:Number}(s::Symbol,i::InitializeFrom,g::GenerateIndicator,nproposals::Integer,::Type{T}) = _policy(Val{s},i,g,nproposals,T)
