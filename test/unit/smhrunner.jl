#create all required objects
smhpolicy1 = policy(:mh,1)
smhpolicy2 = policy(:mh,1,store=:all)
smhrunner1 = runner(smhpolicy1,rniter1;numburnin = rnburnin1)
smhrunner2 = runner(smhpolicy2,rniter1,1;numburnin = rnburnin1)
smhindicator1,smhtunerstate1,smhchain1 = GeneralizedMetropolisHastings.createcommon(smhrunner1,rtuner1,rnparas1,smhnprops1,smhnprops1)
smhindicator2,smhtunerstate2,smhchain2 = GeneralizedMetropolisHastings.createcommon(smhrunner1,rtuner2,rnparas1,smhnprops1,smhnprops1)
smhsamplerstate1 = samplerstate(rsampler1,smhnprops1,Float64,Float64)
smhsamplerstate2 = samplerstate(rsampler1,smhnprops1,Float64,Float64)
@test numproposals(smhindicator1) == 1 && numsamples(smhindicator1) == 1
@test indicatorsamples(smhindicator1) == [0,0]
@test indicatorsamples(smhindicator2) == [0,0]
@test numtunesteps(smhtunerstate1) == 2
@test numparas(smhchain1) == rnparas1 && numsamples(smhchain1) == 4

#test the preparenext! function
sv1 = [1.0]
sll1 = [-100.0]
sv2 = [3.0]
sll2 = [-50.0]
setfrom(smhsamplerstate1,sv1,sll1) ; setproposals(smhsamplerstate1,sv2,sll2) ; setindicator(smhindicator1,[1.0,0.0],[2,1])
preparenext!(smhrunner1,smhsamplerstate1,smhindicator1)
testfrom(smhsamplerstate1,sv2,sll2) #in this condition, value should have been copied
@test_approx_eq mean(smhsamplerstate1.density.distribution) sv2 #and the normal distribution of the mh sampler should also have been conditioned
setfrom(smhsamplerstate1,sv1,sll1) ; setproposals(smhsamplerstate1,sv2,sll2) ; setindicator(smhindicator1,[0.0,1.0],[2,2])
preparenext!(smhrunner1,smhsamplerstate1,smhindicator1)
testfrom(smhsamplerstate1,sv1,sll1) #in this condition, value should have stayed the same
@test_approx_eq mean(smhsamplerstate1.density.distribution) sv2 #and the normal distribution should have been left unchanged from previous test

#test storing a value into the chain
setfrom(smhsamplerstate1,sv1,sll1) ; setproposals(smhsamplerstate1,sv2,sll2) ; setindicator(smhindicator1,[1.0,0.0],[2,1])
store!(smhrunner1,smhsamplerstate1,smhindicator1,smhchain1)
teststore(smhchain1,sv2,sll2[1],1,1,1)
setfrom(smhsamplerstate1,sv1,sll1) ; setproposals(smhsamplerstate1,sv2,sll2) ; setindicator(smhindicator1,[0.0,1.0],[2,2])
store!(smhrunner1,smhsamplerstate1,smhindicator1,smhchain1)
teststore(smhchain1,sv1,sll1[1],2,1,2)
@test smhchain1.values == [sv2[1] sv1[1] 0.0 0.0] #test the exact values

#iterate once
println("============================")
println("Testing of iterate! function")
println("============================")
#reset the samplerstates used below
smhsamplerstate1 = samplerstate(rsampler1,smhnprops1,Float64,Float64)
smhsamplerstate2 = samplerstate(rsampler1,smhnprops1,Float64,Float64)
#carry out all the steps of one iteration
srand(rinitseed1) ; initialize!(smhrunner1,rmodel1,smhsamplerstate1,smhchain1,false)
@test smhsamplerstate1.from.values == smhsamplerstate1.proposals.values
propose!(smhsamplerstate1)
@test smhsamplerstate1.from.values != smhsamplerstate1.proposals.values
geometry!(rmodel1,proposals(smhsamplerstate1))
transitionprobability!(smhindicator1,acceptance!(smhsamplerstate1))
sampleindicator!(smhindicator1)
store!(smhrunner1,smhsamplerstate1,smhindicator1,smhchain1)
#perform an iteration
srand(rinitseed1) ; initialize!(smhrunner1,rmodel1,smhsamplerstate2,smhchain2,false)
iterate!(smhrunner2,rmodel1,smhsamplerstate2,smhindicator2,smhchain2,true)
#compare the values from both steps
testfrom(smhsamplerstate2,smhsamplerstate1.from.values,smhsamplerstate1.from.loglikelihood)
testproposals(smhsamplerstate2,smhsamplerstate1.proposals.values,smhsamplerstate1.proposals.loglikelihood)


#test the tuning function
println("=========================")
println("Testing of tune! function")
println("=========================")
accepted!(smhtunerstate1,smhindicator1)
accepted!(smhtunerstate2,smhindicator1)
tvals1 = tune(rtuner1,smhtunerstate1)
tvals2 = tune(rtuner2,smhtunerstate2)
cov1 = cov(smhsamplerstate1.density.distribution)

@test index(smhtunerstate1) == 1
tune!(smhrunner1,smhsamplerstate1,rtuner1,smhtunerstate1)
@test index(smhtunerstate1) == 2
@test cov1 == cov(smhsamplerstate1.density.distribution)

@test index(smhtunerstate2) == 1
tune!(smhrunner1,smhsamplerstate1,rtuner2,smhtunerstate2)
@test index(smhtunerstate2) == 2
@test_approx_eq cov1*tvals2*tvals2 cov(smhsamplerstate1.density.distribution)

#run the whole mcmc
println("========================")
println("Testing of run! function")
println("========================")
srand(rinitseed1) ; smhc1 = run!(smhrunner1,rmodel1,rsampler1,rtuner1)
srand(rinitseed1) ; smhc2 = run!(smhrunner2,rmodel1,rsampler1,rtuner2)
@test numsamples(smhc1) == 4
@test numsamples(smhc2) == 7
@test samples(smhc2)[:,1] == rpriorvals1
@test samples(smhc1) == samples(smhc2)[:,4:end]
@test loglikelihood(smhc1) == loglikelihood(smhc2)[4:end]
@test logprior(smhc1,rmodel1) == logprior(smhc2,rmodel1)[4:end]
@test logposterior(smhc1,rmodel1) == logposterior(smhc2,rmodel1)[4:end]

println("==============================")
println("Testing of SMHRunner completed")
println("==============================")
