abstract AbstractTuner

###Factory function, currently implemented are :monitor and :scale
tuner(s::Symbol,p::Int,args...;keyargs...) = _tuner(Val{s},p,args...;keyargs...)

period(tuner::AbstractTuner) = tuner.period
verbose(tuner::AbstractTuner) = tuner.verbose
needstuning(tuner::AbstractTuner,i::Int) = mod(i,period(tuner))==0

############################################################################

abstract AbstractTunerState

###Factory function for tunerstates
tunerstate{T<:AbstractFloat}(t::AbstractTuner,nburnin::Int,::Type{T}) = _tunerstate(t,nburnin,T)

accepted(state::AbstractTunerState) = state.accepted
proposed(state::AbstractTunerState) = state.proposed
rate(state::AbstractTunerState) = accepted(state)./proposed(state)
total(state::AbstractTunerState) = sum(proposed(state))
index(state::AbstractTunerState) = state.index
numtunesteps(state::AbstractTunerState) = length(state.accepted)

current(state::AbstractTunerState) = (a = state.accepted[state.index] ; p = state.proposed[state.index] ; tuple(a,p))

tune(tuner::AbstractTuner,state::AbstractTunerState) = ()

@inline function accepted!(state::AbstractTunerState,indicator::AbstractIndicatorMatrix)
    if state.index <= numtunesteps(state)
        state.accepted[state.index] += accepted(indicator)
        state.proposed[state.index] += numsamples(indicator)
    else
        warn("TunerState is full. Additional tuning steps cannot be stored")
    end
end

@inline nextindex!(state::AbstractTunerState) = (state.index+=1)







