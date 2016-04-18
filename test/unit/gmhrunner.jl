println("++++++++++++++++++++++++++++++++")
println("Testing GMHRunner with $(nworkers()) workers")
println("++++++++++++++++++++++++++++++++")

s1 = [1.0,2.0]
ll1 = [-100.0]
s2 = [3.0,4.0]
ll2 = [-50.0]
s3 = [5.0,6.0]
ll3 = [-10.0]

#create all required objects
gmhpolicies = [policy(:mh,gprops1),policy(:mh,gprops1,store=:all,jobsegments=:test)]
gmhrunners = [runner(p,niter1,gprops1;numburnin = nburnin1) for p in gmhpolicies]
gmhchainsamples = [8,13]

for i=1:length(gmhpolicies)
    indicator1,tunerstate1,chain1 = GeneralizedMetropolisHastings.createcommon(gmhrunners[i],tuner1,nparas1,gprops1,gprops1)
    indicator2,tunerstate2,chain2 = GeneralizedMetropolisHastings.createcommon(gmhrunners[i],tuner2,nparas1,gprops1,gprops1)
    indicatorstate1 = samplerstate(sampler1,1,Float64,Float64)
    gmhsegment1 = segment(gmhpolicies[i],model1,sampler1,gprops1)
    remotesegments1 = remotesegments(gmhpolicies[i],model1,sampler1,gprops1)
    @test numproposals(indicator1) == numsamples(indicator1) == gprops1
    @test numsteps(tunerstate1) == 2
    @test numparas(chain1) == nparas1 && numsamples(chain1) == gmhchainsamples[i]
    @test remotesegments1.numsegments == GeneralizedMetropolisHastings._numjobsegments(gmhpolicies[i],gprops1)
    @test remotesegments1.numproposalspersegment == GeneralizedMetropolisHastings._numproposalspersegment(gprops1,remotesegments1.numsegments)
    for j=1:remotesegments1.numsegments
        @test remotesegments1.remote[j].where == collect(GeneralizedMetropolisHastings._processnumbers(gmhpolicies[i],remotesegments1.numsegments))[j]
    end
    println("Testing GMHRunner with $(gprops1) proposals and $(remotesegments1.numsegments) segments")
    println("Testing of updatefrom! function")
    setfrom(indicatorstate1,s1,ll1) ; setproposals(indicatorstate1,s2,ll2) ; setsegment(remotesegments1,s2,ll2,s3,ll3) ; setindicator(indicator1,[1.0,0.0,0.0],[3,1,1])
    updatefrom!(gmhrunners[i],indicatorstate1,remotesegments1,indicator1)
    testfrom(indicatorstate1,s3,ll3)
    setfrom(indicatorstate1,s1,ll1) ; setproposals(indicatorstate1,s2,ll2) ; setsegment(remotesegments1,s2,ll2,s3,ll3) ; setindicator(indicator1,[0.0,0.0,1.0],[3,3,3])
    updatefrom!(gmhrunners[i],indicatorstate1,remotesegments1,indicator1)
    testfrom(indicatorstate1,s1,ll1)

    println("Testing of store! function")
    #test storing a value into the chain
    setfrom(indicatorstate1,s1,ll1) ; setproposals(indicatorstate1,s2,ll2) ; setsegment(remotesegments1,s2,ll2,s3,ll3) ; setindicator(indicator1,[1.0,0.0,0.0],[3,1,1])
    store!(gmhrunners[i],indicatorstate1,remotesegments1,indicator1,chain1)
    teststore(chain1,s3,ll3[1],2,1,2)
    setfrom(indicatorstate1,s1,ll1) ; setproposals(indicatorstate1,s2,ll2) ; setsegment(remotesegments1,s2,ll2,s3,ll3) ; setindicator(indicator1,[0.0,0.0,1.0],[3,3,3])
    store!(gmhrunners[i],indicatorstate1,remotesegments1,indicator1,chain1)
    teststore(chain1,s1,ll1[1],4,1,4)

    println("Testing of iterate! function")
    initialize!(gmhrunners[i],model1,from(indicatorstate1),chain1,false)
    loglikelihood1 = [loglikelihood(model1,evaluate!(model1,defaults1))]
    testfrom(indicatorstate1,defaults1,loglikelihood1)
    srand(435) ; iterate!(gmhrunners[i],model1,indicatorstate1,remotesegments1,indicator1,chain1,true)
    srand(435) ; sample1 = propose!(indicatorstate1.density,samples(:base,nparas1,1,Float64,Float64)) ; geometry!(model1,sample1) ; iterate!(gmhsegment1,sample1)
    testproposals(indicatorstate1,sample1.values,sample1.loglikelihood)
    if indicator1.samples[gprops1+1] == gprops1+1
        testfrom(indicatorstate1,defaults1,loglikelihood1)
        teststore(chain1,defaults1,loglikelihood1[1],6,1,6)
    else
        newfrom = getsamples(remotesegments1,indicator1.samples[gprops1+1])
        testfrom(indicatorstate1,newfrom.values,newfrom.loglikelihood)
        teststore(chain1,newfrom.values,newfrom.loglikelihood[1],6,2,6)
    end

    println("Testing of tune! function")
    accepted!(tunerstate1,indicator1)
    accepted!(tunerstate2,indicator2)
    tvals1 = tune(tuner1,tunerstate1)
    tvals2 = tune(tuner2,tunerstate2)
    cov1 = cov(indicatorstate1.density.distribution)

    @test index(tunerstate1) == 1
    tune!(gmhrunners[i],sampler1,indicatorstate1,remotesegments1,tuner1,tunerstate1)
    @test index(tunerstate1) == 2 && cov1 == cov(indicatorstate1.density.distribution)

    @test index(tunerstate2) == 1
    tune!(gmhrunners[i],sampler1,indicatorstate1,remotesegments1,tuner2,tunerstate2)
    @test index(tunerstate2) == 2 && cov1*tvals2 == cov(indicatorstate1.density.distribution)

    for j=1:remotesegments1.numsegments
        @test cov1*tvals2 == cov(fetch(remotesegments1.remote[j]).samplerstate.density.distribution)
    end

    #run the whole mcmc
    println("Testing of run! function")
    srand(435) ; gmhc = run!(gmhrunners[i],model1,sampler1,tuner1)
    @test numsamples(gmhc) == gmhchainsamples[i]
    println("Results of the run!")
    show(round(samples(gmhc),3)) ; println()
    show(round(loglikelihood(gmhc),3)) ; println()
end

gmhmodel1 = sincosmodel!(0:0.1:10,[0.5,0.5],[0.01,0.01],[0.4,0.4],[0.6,0.6])
gmhsamplerstate1 = samplerstate(sampler1,1,Float64,Float64)
initialize!(gmhrunners[1],gmhmodel1,from(gmhsamplerstate1),chain(:standard,2,10,Float64,Float64),false)
auxiliary!(gmhrunners[1],gmhmodel1,gmhsamplerstate1)

println("====================")
println("Test show() function")
show(gmhrunners[1])
show(gmhrunners[2])
println("End  show() function")
println("====================")
println()

println("==============================")
println("Testing of GMHRunner completed")
println("==============================")

# nothing
