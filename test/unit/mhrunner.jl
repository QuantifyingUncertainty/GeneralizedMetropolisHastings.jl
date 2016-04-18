@test GeneralizedMetropolisHastings._numsamplesinchain(Val{:burnin},10,20,5) == (50,50)
@test GeneralizedMetropolisHastings._numsamplesinchain(Val{:main},10,20,5) == (100,0)
@test GeneralizedMetropolisHastings._numsamplesinchain(Val{:all},10,20,5) == (151,51)

mhpolicy1 = policy(:mh,1)
mhpolicy2 = policy(:mh,1,store=:all)
mhpolicy3 = policy(:mh,1,store=:burnin,initialize=:default)

@test GeneralizedMetropolisHastings._storeduring(:burnin,mhpolicy1) == false
@test GeneralizedMetropolisHastings._storeduring(:burnin,mhpolicy2) == true
@test GeneralizedMetropolisHastings._storeduring(:burnin,mhpolicy3) == true
@test GeneralizedMetropolisHastings._storeduring(:main,mhpolicy1) == true
@test GeneralizedMetropolisHastings._storeduring(:main,mhpolicy2) == true
@test GeneralizedMetropolisHastings._storeduring(:main,mhpolicy3) == false

type TestMHRunner <: AbstractMHRunner
    numburnin::Int
    numiterations::Int
    policy::MHRuntimePolicy
end

sprops1 = 1
gprops1 = 2
nburnin1 = 2
niter1 = 4
sampler1 = sampler(:mh,:normal,0.5,[1.0 0.0;0.0 0.1])
samplerstate1 = samplerstate(sampler1,sprops1,Float64,Float64)
samplerstate2 = samplerstate(sampler1,gprops1,Float64,Float64)
tuner1 = tuner(:monitor,1)
tuner2 = tuner(:scale,2,0.5,:logistic)

#create all required objects
mhrunner1 = TestMHRunner(nburnin1,niter1,mhpolicy1)
mhrunner2 = TestMHRunner(nburnin1,niter1,mhpolicy2)
mhrunner3 = TestMHRunner(nburnin1,niter1,mhpolicy3)
indicator1,tunerstate1,chain1 = GeneralizedMetropolisHastings.createcommon(mhrunner1,tuner1,nparas1,sprops1,sprops1)
indicator2,tunerstate2,chain2 = GeneralizedMetropolisHastings.createcommon(mhrunner2,tuner2,nparas1,gprops1,gprops1)

@test numproposals(indicator1) == sprops1 && numsamples(indicator1) == sprops1
@test numproposals(indicator2) == gprops1 && numsamples(indicator2) == gprops1
@test numtunesteps(tunerstate1) == 2
@test numtunesteps(tunerstate2) == 1
@test numparas(chain1) == nparas1 && numsamples(chain1) == niter1
@test numparas(chain2) == nparas1 && numsamples(chain2) == (nburnin1 + niter1)*gprops1 + 1

#initialize the from field
initialize!(mhrunner1,model1,from(samplerstate1),chain1,true)
loglikelihood1 = [loglikelihood(model1,evaluate!(model1,defaults1))]
testfrom(samplerstate1,defaults1,loglikelihood1)

initialize!(mhrunner1,model1,from(samplerstate2),chain2,true)
loglikelihood2 = [loglikelihood(model1,evaluate!(model1,defaults1))]
testfrom(samplerstate2,defaults1,loglikelihood2)
