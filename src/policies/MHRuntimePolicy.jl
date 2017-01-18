#Generic runtime policy for Metropolis-Hastings runner
immutable MHRuntimePolicy{N<:Number,T<:AbstractFloat} <: AbstractPolicy
    runner::MHRunnerType
    model::ModelType
    initialize::InitializeFrom
    propose::ProposeFrom
    indicator::IndicatorType
    jobsegments::JobSegments
    chain::ChainType
    store::StoreDuring
    sampletype::Type{N}
    calculationtype::Type{T}
    MHRuntimePolicy(t,m,i,p,d,j,c,s) = new(t,m,i,p,d,j,c,s,N,T)
end

#Local function mapping number of proposals to the ProposeFrom type
_num2runner(n::Integer) =  (@assert n > 0 ; n>1?trait(:mhrunner,:generalized):trait(:mhrunner,:standard))
_num2propose(n::Integer) = (@assert n > 0 ; n>1?trait(:propose,:auxiliary):trait(:propose,:indicator))
_num2segments(n::Integer,s::Symbol) = (@assert n > 0 ; n>1?trait(:jobsegments,s):trait(:jobsegments,:none))

_policy(::Type{Val{:mh}},nprops::Int;model = :deterministic,initialize = :prior,indicator = :stationary,jobsegments = :workers,
        chain = :standard,store = :main,sampletype =Float64,calculationtype =Float64) =
    MHRuntimePolicy{sampletype,calculationtype}(_num2runner(nprops),trait(:model,model),trait(:initialize,initialize),_num2propose(nprops),
                                                trait(:indicator,indicator),_num2segments(nprops,jobsegments),trait(:chain,chain),trait(:store,store))

function show(io::IO,p::MHRuntimePolicy)
    println(io,"MHRuntimePolicy with traits:")
    println(io,"  runner = ",traitvalue(p.runner))
    println(io,"  model = ",traitvalue(p.model))
    println(io,"  initialize = ",traitvalue(p.initialize))
    println(io,"  propose = ",traitvalue(p.propose))
    println(io,"  indicator = ",traitvalue(p.indicator))
    println(io,"  jobsegments = ",traitvalue(p.jobsegments))
    println(io,"  chain = ",traitvalue(p.chain))
    println(io,"  store = ",traitvalue(p.store))
    println(io,"  sampletype = ",p.sampletype)
    println(io,"  calculationtype = ",p.calculationtype)
    println(io)
    nothing
end
