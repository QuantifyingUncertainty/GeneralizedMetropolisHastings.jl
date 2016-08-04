### ScaleTuner

# ScaleTuner tunes the scaling or stepsize parameters of a sampler's proposal step

immutable ScaleTuner <: AbstractTuner
    period::Int
    verbose::Bool
    targetrate::Real
    score::Function
    scoreargs::Tuple
    function ScaleTuner(period::Int,verbose::Bool,targetrate::Real,score::Function,scoreargs...)
        @assert 0 < period "Tuning period should be positive"
        @assert 0 <= targetrate <= 1 "Target acceptance rate should be between 0 and 1"
        new(period,verbose,targetrate,score,scoreargs)
    end
end

_tunerfunc(::Type{Val{:logistic}}) = (tunelogistic,7.0)
_tunerfunc(::Type{Val{:erf}}) = (tuneerf,3.0)

_tuner(::Type{Val{:scale}},period::Int,targetrate::Real,score::Symbol;verbose::Bool =true) = ScaleTuner(period,verbose,targetrate,_tunerfunc(Val{score})...)
_tuner(::Type{Val{:scale}},period::Int,targetrate::Real,score::Function,scoreargs...;verbose::Bool =true) = ScaleTuner(period,verbose,targetrate,score,scoreargs...)

##########################################################################################################

type ScaleTunerState{T<:AbstractFloat} <: AbstractTunerState
    accepted::Vector{Int}
    proposed::Vector{Int}
    scalefactor::Vector{T}
    index::Int
    function ScaleTunerState(accepted::Vector{Int},proposed::Vector{Int},scalefactor::Vector{T},index::Int)
        @assert index > 0 "Total number of proposed MCMC samples should be non-negative"
        new(accepted,proposed,scalefactor,index)
    end
end

function _tunerstate{T<:AbstractFloat}(tuner::ScaleTuner,nburnin::Int,::Type{T})
    nsteps = ceil(Int,nburnin/tuner.period)
    ScaleTunerState{T}(zeros(Int,nsteps),zeros(Int,nsteps),zeros(T,nsteps),1)
end

function tune(tuner::ScaleTuner,state::ScaleTunerState)
    r = state.accepted[state.index]/state.proposed[state.index]
    s = tuner.score(r-tuner.targetrate,tuner.scoreargs...)
    state.scalefactor[state.index] = convert(eltype(state.scalefactor),s)
end

function show(io::IO,tuner::ScaleTuner)
    println(io,"ScaleTuner: period = $(tuner.period), verbose = $(tuner.verbose) ,targetrate = $(tuner.targetrate)")
    println(io,"            score = $(tuner.score)",isempty(tuner.scoreargs)?"":", scoreargs = $(tuner.scoreargs)")
end

function show(io::IO,state::ScaleTunerState)
    println(io,"ScaleTunerState: ")
    println(io," accepted = $(accepted(state))")
    println(io," proposed = $(proposed(state))")
    println(io," scalefactor = $(round(state.scalefactor,3))")
    println(io," acceptance rate = $(round(rate(state),3))")
    println(io," total proposed = $(total(state))")
end

function showstep(t::ScaleTuner,state::ScaleTunerState)
    if verbose(t) && state.index <= numtunesteps(state)
        a,p = current(state)
        println("  accepted/proposed = $a/$p")
        println("  acceptance rate = $(round(a/p,3))")
        println("  scalefactor = $(round(state.scalefactor[state.index],3))")
        println("  cummulative scaling = $(round((r = one(eltype(state.scalefactor)) ; for i=1:state.index r*=state.scalefactor[i] end ; r),3))")
    end
end
