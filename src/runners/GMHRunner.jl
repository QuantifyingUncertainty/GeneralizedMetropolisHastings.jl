### Generalized Metropolis-Hastings Monte Carlo runner
type GMHRunner <: MCRunner
  nburnin::Int #number of burnin iterations (tuning only during burnin)
  nsamples::Int #total number of iterations
  nproposals::Int #number of proposals calculated in parallel
  niterations::Int #number of parallel iterations executed (derived from numsamples and numproposals)
  policy::RuntimePolicy #policies determining runtime behaviour of the MCMC
  indicator::Int #the indicator variable

  function GMHRunner(nb::Int,ns::Int,np::Int,pol::RuntimePolicy,i::Int)
    @assert nb >= 0 "Number of burn-in iterations should be non-negative."
    @assert ns > nb "Total number of MCMC iterations should be greater than number of burn-in iterations."
    @assert np >= 0 && np < ns "Number of proposals should be larger >= 0 and smaller than total of MCMC iterations"
    ni::Int = ns%np?(ns÷np)+1:(ns÷np) #number of iterations is given by integer division of ns/np (rounded up)
    new(nb,ns,np,ni,pol,i)
  end
end

GMHRunner(ns::Int,np::Int;burnin::Int = 0,initialize::ValuesFrom = ValuesFromDefault(),indicate::IndicatorMatrixFunction = IndicatorMatrixStationary(),indicator::Int =1)
  = GMHRunner(burnin,ns,np,GenericPolicy(initialize,indicate,np),indicator)

function run!(r::GMHRunner,m::MCModel,s::MCSampler,h::MCHeap)
  tic()
  for i in 1:r.niterations
    iterate!(r,m,s,h)
    #TODO: store results in mcchain
  end

  mcchain.diagnostics, mcchain.runtime = ds, toq()
  mcchain
end

function iterate!(r::GMHRunner,m::MCModel,s::MCSampler,h::MCHeap)
  propose!(r.policy.indicate,m,s,h,r.indicator)
  update_geometry!(m,s,h,r.indicator)
  update_proposals!(s,h,r.indicator)
  sample_indicator!(m,s,h,r.indicator)
end

#Propose from indicator, will only be called if numproposals == 1 (standard M-H)
function propose!(::ProposalFromIndicator,m::MCModel,s::MCSampler,h::MCHeap,indicator::Int)
  propose!(s,h,indicator)
end

#Propose from auxiliary, will be called if numproposals > 1 (generalized M-H)
function propose!(::ProposalFromAuxiliary,m::MCModel,s::MCSampler,h::MCHeap,indicator::Int)
  aux = copy(h.samples[indicator])
  set_from!(s,h,h.samples[indicator])
  propose!(s,h,aux) #calls a single propose, storing the result in auxiliary
  update_geometry!(m,s,h,aux)
  update_proposal!(m,s,aux)
  set_from!(s,h,aux) #now set the location to propose from
  propose!(s,h,indicator) #calls the vectorized propose function (see samplers.jl)
end

###Main entry point in geometry calculations, parallel function using pmap because calculating the geometry can be costly
update_geometry!(m::MCModel,s::MCSampler,h::MCHeap,indicator::Int) = pmap((j)->(j!=indicator?update_geometry!(m,s,h,h.samples[j]):nothing,1:length(h.samples)))

function resume!(m::MCModel, s::MCSampler, r::SerialMC, c::MCChain, t::MCTuner=VanillaMCTuner(), j::Symbol=:task;
  nsteps::Int=100)
  m.init = vec(c.samples[end, :])
  mcrunner::SerialMC = SerialMC(burnin=0, thinning=r.thinning, nsteps=nsteps, storegradlogtarget=r.storegradlogtarget)
  run(m, s, mcrunner, t, j)
end

resume(m::MCModel, s::MCSampler, r::SerialMC, c::MCChain, t::MCTuner=VanillaMCTuner(), j::Symbol=:task;
  nsteps::Int=100) =
  resume!(deepcopy(m), s, r, c, t, j; nsteps=nsteps)

#tuning, resumiing, noise models on the data, uniform data interface different samplers
