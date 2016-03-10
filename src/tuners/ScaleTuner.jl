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

type ScaleTunerState <: AbstractTunerState
    accepted::Int
    proposed::Int
    totalproposed::Int
    function ScaleTunerState(accepted::Int,proposed::Int,totalproposed::Int)
        @assert accepted >= 0 "Number of accepted MCMC samples should be non-negative"
        @assert proposed >= 0 "Number of proposed MCMC samples should be non-negative"
        @assert totalproposed >= 0 "Total number of proposed MCMC samples should be non-negative"
        new(accepted,proposed,totalproposed)
    end
end

_tunerstate(tuner::ScaleTuner,a::Int,p::Int,t::Int) = ScaleTunerState(a,p,t)
_tunerstate(tuner::ScaleTuner) = ScaleTunerState(0,0,0)

tune(tuner::ScaleTuner,state::ScaleTunerState) = max(1/2,tuner.score(state.accepted/state.proposed-tuner.targetrate,tuner.scoreargs...))

function show(io::IO,tuner::ScaleTuner)
    println(io,"ScaleTuner: period = $(tuner.period), verbose = $(tuner.verbose) ,targetrate = $(tuner.targetrate)")
    println(io,"            score = $(tuner.score)",isempty(tuner.scoreargs)?"":", scoreargs = $(tuner.scoreargs)")
end

function show(io::IO,state::ScaleTunerState)
    println(io,"ScaleTunerState: accepted = $(state.accepted), proposed = $(state.proposed), totalproposed = $(state.totalproposed)")
end
