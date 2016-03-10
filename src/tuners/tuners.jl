abstract AbstractTuner

###Factory function, currently implemented are :monitor and :scale
tuner(s::Symbol,p::Int,args...;keyargs...) = _tuner(Val{s},p,args...;keyargs...)

period(tuner::AbstractTuner) = tuner.period
verbose(tuner::AbstractTuner) = tuner.verbose

############################################################################ver

abstract AbstractTunerState

###Factory function for tunerstates
tunerstate(t::AbstractTuner,args...) = _tunerstate(t,args...)

tune(tuner::AbstractTuner,state::AbstractTunerState) = ()

accepted(state::AbstractTunerState) = state.accepted
proposed(state::AbstractTunerState) = state.proposed
totalproposed(state::AbstractTunerState) = state.totalproposed
rate(state::AbstractTunerState) = state.accepted/state.proposed

@inline function accepted!(state::AbstractTunerState,indicator::AbstractIndicatorMatrix)
    state.accepted += accepted(indicator)
    state.proposed += numsamples(indicator)
    state.totalproposed += numsamples(indicator)
end

@inline function resetburnin!(state::AbstractTunerState)
    (state.accepted, state.proposed) = (0, 0)
    state
end



