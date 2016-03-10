### MonitorTuner

# MonitorTuner is the most elemental tuner type
# It only monitors the acceptance rate of the current period

immutable MonitorTuner <: AbstractTuner
    period::Int # Tuning period over which acceptance rate is computed
    verbose::Bool

    function MonitorTuner(period::Int,verbose::Bool)
      @assert period > 0 "Tuning period should be positive"
      new(period,verbose)
    end
end

_tuner(::Type{Val{:monitor}},period::Int;verbose::Bool =true) = MonitorTuner(period,verbose)

##########################################################################################################

type MonitorTunerState <: AbstractTunerState
    accepted::Int # Number of accepted MCMC samples during current tuning period
    proposed::Int # Number of proposed MCMC samples during current tuning period
    totalproposed::Int # Total number of proposed MCMC samples during burnin

    function MonitorTunerState(accepted::Int,proposed::Int,totalproposed::Int)
        @assert accepted >= 0 "Number of accepted MCMC samples should be non-negative"
        @assert proposed >= 0 "Number of proposed MCMC samples should be non-negative"
        @assert totalproposed >= 0 "Total number of proposed MCMC samples should be non-negative"
        new(accepted, proposed, totalproposed)
    end
end

_tunerstate(tuner::MonitorTuner,accepted::Int,proposed::Int,totalproposed::Int) = MonitorTunerState(accepted,proposed,totalproposed)
_tunerstate(tuner::MonitorTuner) = MonitorTunerState(0,0,0)

show(io::IO,tuner::MonitorTuner) = println(io,"MonitorTuner: period = $(tuner.period), verbose = $(tuner.verbose)")

function show(io::IO,state::MonitorTunerState)
    println(io,"MonitorTunerState: accepted = $(state.accepted), proposed = $(state.proposed), totalproposed = $(state.totalproposed)")
end


