#Generic runtime policy
immutable GMHPolicy{N<:Number,T<:AbstractFloat} <: AbstractPolicy
    initialize::InitializeFrom
    propose::ProposeFrom
    indicator::IndicatorType
    samplerstates::SamplerStates
    sampletype::Type{N}
    calculationtype::Type{T}
    GMHPolicy(i,p,d,s) = new(i,p,d,s,N,T)
end

#Local function mapping number of proposals to the ProposeFrom type
_num2propose(n::Integer) = (@assert n > 0 ; n>1?trait(:propose,:auxiliary):trait(:propose,:indicator))

_policy(::Type{Val{:gmh}},nprops::Int;
        initialize =:prior,
        indicator =:stationary,
        samplerstates =:nworkers,
        sampletype =Float64,
        calculationtype =Float64) =
    GMHPolicy{sampletype,calculationtype}(trait(:initialize,initialize),_num2propose(nprops),trait(:indicator,indicator),trait(:samplerstates,samplerstates))

function show(io::IO,p::GMHPolicy)
    println(io,"GMHPolicy with traits:")
    println(io,"  initialize = ",traitvalue(p.initialize))
    println(io,"  propose = ",traitvalue(p.propose))
    println(io,"  indicator = ",traitvalue(p.indicator))
    println(io,"  samplerstates = ",traitvalue(p.samplerstates))
    println(io,"  sampletype = ",p.sampletype)
    println(io,"  calculationtype = ",p.calculationtype)
    println(io)
    nothing
end
