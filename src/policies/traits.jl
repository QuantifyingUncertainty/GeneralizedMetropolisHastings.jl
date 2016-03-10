abstract AbstractPolicyTrait

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

immutable SamplerStates <: AbstractPolicyTrait
    trait::Symbol
    SamplerStates(s::Symbol) = (@assert in(s,[:nprocs,:nworkers,:test]) ; new(s))
end

trait(n::Symbol,t::Symbol) = _trait(Val{n},t)
_trait(::Type{Val{:initialize}},t::Symbol) = InitializeFrom(t)
_trait(::Type{Val{:propose}},t::Symbol) = ProposeFrom(t)
_trait(::Type{Val{:indicator}},t::Symbol) = IndicatorType(t)
_trait(::Type{Val{:samplerstates}},t::Symbol) = SamplerStates(t)

traitvalue(a::AbstractPolicyTrait) = a.trait
traittype(a::AbstractPolicyTrait) = Val{a.trait}
