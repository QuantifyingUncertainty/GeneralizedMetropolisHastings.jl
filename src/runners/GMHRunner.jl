### Generalized Metropolis-Hastings Monte Carlo runner
immutable GMHRunner <: MCRunner
  nburnin::Int #number of burnin iterations (tuning only during burnin)
  nsamples::Int #total number of iterations
  nproposals::Int #number of proposals calculated in parallel
  niterations::Int #number of parallel iterations executed (derived from numsamples and numproposals)
  policy::RuntimePolicy #policies determining runtime behaviour of the MCMC

  function GMHRunner(nb::Int,ns::Int,np::Int,pol::RuntimePolicy)
    @assert nb >= 0 "Number of burn-in iterations should be non-negative."
    @assert ns > nb "Total number of MCMC iterations should be greater than number of burn-in iterations."
    @assert np >= 0 && np < ns "Number of proposals should be larger >= 0 and smaller than total of MCMC iterations"
    ni::Int = ns%np!=0?(ns÷np)+1:(ns÷np) #number of iterations is given by integer division of ns/np (rounded up)
    new(nb,np*ni,np,ni,pol)
  end
end

GMHRunner(ns::Int,np::Int; burnin =0, initialize =ValuesFromPrior(), indicate =IndicatorMatrixStationary()) = GMHRunner(burnin,ns,np,GenericPolicy(initialize,indicate,np))

function run!(r::GMHRunner,m::MCModel,s::MCSampler)
  h = create_heap(s,r.nproposals + 1)
  c = create_chain(numparas(s),r.nsamples)
  tic()
  indicator::Int = 1
  initialize!(r,m,s,h,indicator)
  for i = 1:r.niterations
    indicators = iterate!(r,m,s,h,indicator)
    store_results!(c,h,indicators)
    indicator = indicators[end]
  end
  c.runtime = toq()
  c
end

###Initialize the sample of the indicator variable
function initialize!(r::GMHRunner,m::MCModel,s::MCSampler,h::MCHeap,indicator)
  h.samples[indicator].values = initvalues(r.policy.initialize,m.parameters) #initialize the indicator sample with values
  update_geometry!(m,h.samples[indicator]) #update the geometry for the indicator sample
end

function iterate!(r::GMHRunner,m::MCModel,s::MCSampler,h::MCHeap,indicator::Int)
  propose!(r.policy.propose,m,s,h,indicator)
  update_geometry!(r,m,h,indicator)
  update_proposals!(s,h,indicator)
  sample_indicator(h,indicator,r.policy.indicate)
end

#Propose from indicator, will only be called if numproposals == 1 (standard M-H)
function propose!(::ProposalFromIndicator,m::MCModel,s::MCSampler,h::MCHeap,indicator::Int)
  propose!(s,h,indicator)
end

#Propose from auxiliary, will be called if numproposals > 1 (generalized M-H)
function propose!(::ProposalFromAuxiliary,m::MCModel,s::MCSampler,h::MCHeap,indicator::Int)
  aux = deepcopy(h.samples[indicator])
  set_from!(s,h,h.samples[indicator])
  propose!(s,h,aux) #calls a single propose, storing the result in auxiliary
  update_geometry!(m,aux)
  #update_proposal!(s,h,aux) TODO
  set_from!(s,h,aux) #now set the location to propose from
  propose!(s,h,indicator) #calls the vectorized propose function (see samplers.jl)
end

###Main entry point in geometry calculations, parallel function using pmap because calculating the geometry can be costly
function update_geometry!(r::GMHRunner,m::MCModel,h::MCHeap,indicator::Int)
  res = pmap((j)->(j!=indicator?update_geometry!(m,h.samples[j]):nothing),1:length(h.samples))
  for i=1:length(res)
    if (i != indicator)
      h.samples[i] = res[i]
    end
  end
end

function Base.show(io::IO,r::GMHRunner)
  println(io,"Generalized Metropolis-Hastings runner with:")
  println(io,"  nburnin: ",r.nburnin)
  println(io,"  nsamples: ",r.nsamples)
  println(io,"  nproposals: ",r.nproposals)
  println(io,"  niterations: ",r.niterations)
  println(io,"  policy: ")
  show(io,r.policy,"   ")
  println(io)
  nothing
end
#tuning, resumiing, noise models on the data, uniform data interface different samplers
