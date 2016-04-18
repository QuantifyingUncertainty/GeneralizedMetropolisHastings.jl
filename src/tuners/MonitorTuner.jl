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
    accepted::Vector{Int} # Number of accepted MCMC samples during current tuning period
    proposed::Vector{Int} # Number of proposed MCMC samples during current tuning period
    index::Int

    function MonitorTunerState(accepted::Vector{Int},proposed::Vector{Int},index::Int)
        @assert index > 0 "Total number of proposed MCMC samples should be non-negative"
        new(accepted,proposed,index)
    end
end

function _tunerstate{T<:AbstractFloat}(tuner::MonitorTuner,nburnin::Int,::Type{T})
    nsteps = ceil(Int,nburnin/tuner.period)
    MonitorTunerState(zeros(Int,nsteps),zeros(Int,nsteps),1)
end

show(io::IO,tuner::MonitorTuner) = println(io,"MonitorTuner: period = $(tuner.period), verbose = $(tuner.verbose)")

function show(io::IO,state::MonitorTunerState)
    println(io,"MonitorTunerState: ")
    println(io," accepted = $(accepted(state))")
    println(io," proposed = $(proposed(state))")
    println(io," acceptance rate = $(round(rate(state),3))")
    println(io," total proposed = $(total(state))")
end

function showstep(t::MonitorTuner,state::MonitorTunerState)
    if verbose(t) && state.index <= numtunesteps(state)
        a,p = current(state)
        println("  accepted/proposed = $(a)/$(p)")
        println("  acceptance rate = $(round(a/p,3))")
    end
end




