#create all required objects
smhpolicy1 = policy(:mh,1)
smhpolicy2 = policy(:mh,1,store=:all)
smhrunner1 = runner(smhpolicy1,niter1;numburnin = nburnin1)
smhrunner2 = runner(smhpolicy2,niter1;numburnin = nburnin1)
indicator1,tunerstate1,chain1 = GeneralizedMetropolisHastings.createcommon(smhrunner1,tuner1,nparas1,sprops1,sprops1)
samplerstate1 = samplerstate(sampler1,sprops1,Float64,Float64)
@test numproposals(indicator1) == 1 && numsamples(indicator1) == 1
@test numtunesteps(tunerstate1) == 2
@test numparas(chain1) == 2 && numsamples(chain1) == 4

#test updating the from field
s1 = [1.0,2.0]
ll1 = [-100.0]
s2 = [3.0,4.0]
ll2 = [-50.0]
setfrom(samplerstate1,s1,ll1) ; setproposals(samplerstate1,s2,ll2) ; setindicator(indicator1,[1.0,0.0],[2,1])
updatefrom!(smhrunner1,samplerstate1,indicator1)
testfrom(samplerstate1,s2,ll2)
setfrom(samplerstate1,s1,ll1) ; setproposals(samplerstate1,s2,ll2) ; setindicator(indicator1,[0.0,1.0],[2,2])
updatefrom!(smhrunner1,samplerstate1,indicator1)
testfrom(samplerstate1,s1,ll1)

#test storing a value into the chain
setfrom(samplerstate1,s1,ll1) ; setproposals(samplerstate1,s2,ll2) ; setindicator(indicator1,[1.0,0.0],[2,1])
store!(smhrunner1,samplerstate1,indicator1,chain1)
teststore(chain1,s2,ll2[1],1,1,1)
setfrom(samplerstate1,s1,ll1) ; setproposals(samplerstate1,s2,ll2) ; setindicator(indicator1,[0.0,1.0],[2,2])
store!(smhrunner1,samplerstate1,indicator1,chain1)
teststore(chain1,s1,ll1[1],2,1,2)

#iterate once
println("============================")
println("Testing of iterate! function")
println("============================")
initialize!(smhrunner1,model1,from(samplerstate1),chain1,false)
loglikelihood1 = [loglikelihood(model1,evaluate!(model1,defaults1))]
testfrom(samplerstate1,defaults1,loglikelihood1)
srand(435) ; iterate!(smhrunner1,model1,samplerstate1,indicator1,chain1,true)
srand(435) ; sample1 = propose!(samplerstate1.density,samples(:base,nparas1,sprops1,Float64,Float64)) ; geometry!(model1,sample1)
testproposals(samplerstate1,sample1.values,sample1.loglikelihood)
if indicator1.samples[2] == 2
    testfrom(samplerstate1,defaults1,loglikelihood1)
    teststore(chain1,default1,loglikelihood1[1],3,1,3)
else
    testfrom(samplerstate1,sample1.values,sample1.loglikelihood)
    teststore(chain1,sample1.values,sample1.loglikelihood[1],3,2,3)
end

#test the tuning function
println("=========================")
println("Testing of tune! function")
println("=========================")
accepted!(tunerstate1,indicator1)
accepted!(tunerstate2,indicator1)
tvals1 = tune(tuner1,tunerstate1)
tvals2 = tune(tuner2,tunerstate2)
cov1 = cov(samplerstate1.density.distribution)

@test index(tunerstate1) == 1
tune!(smhrunner1,sampler1,samplerstate1,tuner1,tunerstate1)
@test index(tunerstate1) == 2 && cov1 == cov(samplerstate1.density.distribution)

@test index(tunerstate2) == 1
tune!(smhrunner1,sampler1,samplerstate1,tuner2,tunerstate2)
@test index(tunerstate2) == 2 && cov1*tvals2 == cov(samplerstate1.density.distribution)

#run the whole mcmc
println("========================")
println("Testing of run! function")
println("========================")
srand(435) ; smhc1 = run!(smhrunner1,model1,sampler1,tuner1)
srand(435) ; smhc2 = run!(smhrunner2,model1,sampler1,tuner2)
@test numsamples(smhc1) == 4
@test numsamples(smhc2) == 7
@test samples(smhc2)[:,1] == defaults1
@test samples(smhc2)[:,2] == sample1.values
@test samples(smhc1) == samples(smhc2)[:,4:end]
@test loglikelihood(smhc1) == loglikelihood(smhc2)[4:end]
@test logprior(smhc1,model1) == logprior(smhc2,model1)[4:end]
@test logposterior(smhc1,model1) == logposterior(smhc2,model1)[4:end]

println("==============================")
println("Testing of SMHRunner completed")
println("==============================")
