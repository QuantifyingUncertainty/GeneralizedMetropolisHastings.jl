import GeneralizedMetropolisHastings.MCModel,
  GeneralizedMetropolisHastings.MCSampler,
  GeneralizedMetropolisHastings.MCHeap

function test_individual_steps!(r::GMHRunner,m::MCModel,s::MCSampler,h::MCHeap)

  ###Test the initialize! function
  for i=1:numel(h) @test h.samples[i].values == zeros(numparas(s)) end
  initialize!(r,m,s,h,1)
  @test h.samples[1].values != zeros(numparas(s)) && isfinite(h.samples[1].loglikelihood) && isfinite(h.samples[1].logprior)
  for i=2:numel(h) @test h.samples[i].values == zeros(numparas(s)) end

  ###Test the propose! function
  @test h.samples[1].values != zeros(numparas(s)) && h.samples[2].values == zeros(numparas(s))
  propose!(r.policy.propose,m,s,h,1)
  @test h.samples[1].values != zeros(numparas(s)) && h.samples[2].values != zeros(numparas(s))

  ###Test the update_geometry! function
  @test ~isfinite(h.samples[2].loglikelihood) && ~isfinite(h.samples[2].logprior)
  update_geometry!(r,m,h,1)
  @test isfinite(h.samples[2].loglikelihood) && isfinite(h.samples[2].logprior)

  ###Test the sample_indicator function
  sample_indicator(h,1,r.policy.indicate)

end
###Test the runner for an ODEModel using a spring/mass dynamic system
#Initialize data and variables
nsamples = 20
nprops1 = 1
nprops2 = 2
nprops3 = 4
t1 = [0:0.1:10]
i1 = [-1.0,1.0]
p1 = [100.0,10.0]
n1 = [1e-1 1e-1]
m1 = spring_mass_model(t1,i1,p1,n1)
s1 = MHNormal(eye(2),0.01)
h1 = MHHeap(s1,nprops1+1) #to test Standard Metropolis-Hastings
h2 = MHHeap(s1,nprops2+1) #to test Generalized Metropolis-Hastings
h3 = MHHeap(s1,nprops3+1)

###Test the GMHRunner constructors
r1 = GMHRunner(nsamples,nprops1)
@test r1 == GMHRunner(0,nsamples,nprops1,GenericPolicy(ValuesFromPrior(),IndicatorMatrixStationary(),nprops1)) && r1.policy.propose == ProposalFromIndicator()
r2 = GMHRunner(nsamples,nprops2)
@test r2 == GMHRunner(0,nsamples,nprops2,GenericPolicy(ValuesFromPrior(),IndicatorMatrixStationary(),nprops2)) && r2.policy.propose == ProposalFromAuxiliary()
r3 = GMHRunner(nsamples,nprops3)
@test r3 == GMHRunner(0,nsamples,nprops3,GenericPolicy(ValuesFromPrior(),IndicatorMatrixStationary(),nprops3))
@test r3.niterations == div(nsamples,nprops3) && r3.nsamples == nsamples
r4 = GMHRunner(100,6)
@test r4.niterations == 17 && r4.nsamples == 102
@test GMHRunner(nsamples,nprops1) == GMHRunner(nsamples,nprops1; burnin = 0, initialize = ValuesFromPrior(), indicate = IndicatorMatrixStationary())

println("====================")
println("Test show() function")
show(r1)
show(r2)
show(r3)
println("End  show() function")
println("====================")
println()

###Test the Standard MH with a MHNormal sampler
srand(0)
i1a = test_individual_steps!(r1,m1,s1,h1)
h1a = deepcopy(h1)
srand(0)
initialize!(r1,m1,s1,h1,1)
i1b = iterate!(r1,m1,s1,h1,1)
@test h1 == h1a && i1a == i1b == [1,2]
i1c = iterate!(r1,m1,s1,h1,i1a[end])
r1 = GMHRunner(nsamples,nprops1)
srand(0)
c1 = run!(r1,m1,s1,h1)
show(c1)

###Test Generalized MH with a MHNormal sampler
srand(0)
i2a = test_individual_steps!(r2,m1,s1,h2)
h2a = deepcopy(h2)
srand(0)
initialize!(r2,m1,s1,h2,1)
i2b = iterate!(r2,m1,s1,h2,1)
@test h2 == h2a && i2a == i2b == [1,3,3]
i2c = iterate!(r2,m1,s1,h2,i2a[end])
r2 = GMHRunner(nsamples,nprops2)
srand(0)
c2 = run!(r2,m1,s1,h2)
show(c2)

###Test Generalized MH with a MHNormal sampler
srand(0)
i3a = test_individual_steps!(r3,m1,s1,h3)
h3a = deepcopy(h3)
srand(0)
initialize!(r3,m1,s1,h3,1)
i3b = iterate!(r3,m1,s1,h3,1)
@test h3 == h3a && i3a == i3b == [1,5,5,5,5]
i3c = iterate!(r3,m1,s1,h3,i3a[end])
r3 = GMHRunner(nsamples,nprops3)
srand(0)
c3 = run!(r3,m1,s1,h3)
show(c3)

nothing


