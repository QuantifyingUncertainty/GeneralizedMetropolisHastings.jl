abstract AbstractPolicyTrait

immutable MHRunnerType <: AbstractPolicyTrait
    trait::Symbol
    MHRunnerType(s::Symbol) = (@assert in(s,[:standard,:generalized]) ; new(s))
end

immutable ModelType <: AbstractPolicyTrait
    trait::Symbol
    ModelType(s::Symbol) = (@assert in(s,[:deterministic,:stochastic]) ; new(s))
end

immutable InitializeFrom <: AbstractPolicyTrait
    trait::Symbol
    InitializeFrom(s::Symbol) = (@assert in(s,[:default,:prior]) ; new(s))
end

immutable ProposeFrom <: AbstractPolicyTrait
    trait::Symbol
    ProposeFrom(s::Symbol) = (@assert in(s,[:indicator,:auxiliary]) ; new(s))
end

immutable IndicatorType <: AbstractPolicyTrait
    trait::Symbol
    IndicatorType(s::Symbol) = (@assert in(s,[:stationary,:cyclical]) ; new(s))
end

immutable JobSegments <: AbstractPolicyTrait
    trait::Symbol
    JobSegments(s::Symbol) = (@assert in(s,[:procs,:workers,:none,:one,:two]) ; new(s))
end

immutable ChainType <: AbstractPolicyTrait
    trait::Symbol
    ChainType(s::Symbol) = (@assert in(s,[:standard,:gradient]) ; new(s))
end

immutable StoreDuring <: AbstractPolicyTrait
    trait::Symbol
    StoreDuring(s::Symbol) = (@assert in(s,[:burnin,:main,:all]) ; new(s))
end

@inline _traitname(::Type{MHRunnerType}) = :mhrunner
@inline _traitname(::Type{ModelType}) = :model
@inline _traitname(::Type{InitializeFrom}) = :initialize
@inline _traitname(::Type{ProposeFrom}) = :propose
@inline _traitname(::Type{IndicatorType}) = :indicator
@inline _traitname(::Type{JobSegments}) = :jobsegments
@inline _traitname(::Type{ChainType}) = :chain
@inline _traitname(::Type{StoreDuring}) = :store
@inline _traitnametype(t) = Val{_traitname(t)}

trait(n::Symbol,t::Symbol) = _trait(Val{n},t)
_trait(::Type{_traitnametype(MHRunnerType)},t::Symbol)  = MHRunnerType(t)
_trait(::Type{_traitnametype(ModelType)},t::Symbol) = ModelType(t)
_trait(::Type{_traitnametype(InitializeFrom)},t::Symbol) = InitializeFrom(t)
_trait(::Type{_traitnametype(ProposeFrom)},t::Symbol) = ProposeFrom(t)
_trait(::Type{_traitnametype(IndicatorType)},t::Symbol) = IndicatorType(t)
_trait(::Type{_traitnametype(JobSegments)},t::Symbol) = JobSegments(t)
_trait(::Type{_traitnametype(ChainType)},t::Symbol) = ChainType(t)
_trait(::Type{_traitnametype(StoreDuring)},t::Symbol) = StoreDuring(t)

traitvalue(a::AbstractPolicyTrait) = a.trait
traittype(a::AbstractPolicyTrait) = Val{a.trait}
