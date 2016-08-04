gv1 = [1.0]
gll1 = [-100.0]
gv2 = [3.0]
gll2 = [-50.0]
gv3 = [5.0]
gll3 = [-10.0]
nprocesses = [1,2,3]

for n in nprocesses
    try
        if nprocs() != n
            addprocs(n-nprocs())
            include(joinpath(Pkg.dir("GeneralizedMetropolisHastings"),"test","imports.jl"))
        end

        println("++++++++++++++++++++++++++++++++")
        println("Testing GMHRunner with $(nprocs()) processes")
        println("++++++++++++++++++++++++++++++++")

        #create all required objects
        gmhpolicies = [policy(:mh,gmhnprops1),policy(:mh,gmhnprops1,store=:all,jobsegments=:test)]
        gmhrunners = [runner(p,rniter1,gmhnprops1;numburnin = rnburnin1) for p in gmhpolicies]
        gmhchainsamples = [8,13]

        for i=1:length(gmhpolicies)
            indicator1,tunerstate1,chain1 = GeneralizedMetropolisHastings.createcommon(gmhrunners[i],rtuner1,rnparas1,gmhnprops1,gmhnprops1)
            indicator2,tunerstate2,chain2 = GeneralizedMetropolisHastings.createcommon(gmhrunners[i],rtuner2,rnparas1,gmhnprops1,gmhnprops1)
            indicatorstate1 = samplerstate(rsampler1,1,Float64,Float64)
            remotesegments1 = remotesegments(gmhpolicies[i],rmodel1,rsampler1,gmhnprops1)

            println("Testing GMHRunner with $(gmhnprops1) proposals and $(remotesegments1.numsegments) segments")
            @test numproposals(indicator1) == numsamples(indicator1) == gmhnprops1
            @test numtunesteps(tunerstate1) == 2
            @test numparas(chain1) == rnparas1 && numsamples(chain1) == gmhchainsamples[i]
            @test remotesegments1.numsegments == GeneralizedMetropolisHastings._numjobsegments(gmhpolicies[i],gmhnprops1)
            @test remotesegments1.numproposalspersegment == GeneralizedMetropolisHastings._numproposalspersegment(gmhnprops1,remotesegments1.numsegments)
            for j=1:remotesegments1.numsegments
                @test remotesegments1.remote[j].where == collect(GeneralizedMetropolisHastings._processnumbers(gmhpolicies[i],remotesegments1.numsegments))[j]
            end

            println("Testing of preparenext! function")
            setfrom(indicatorstate1,gv1,gll1) ; setproposals(indicatorstate1,gv2,gll2) ; setsegment(remotesegments1,gv2,gll2,gv3,gll3,1) ; setindicator(indicator1,[1.0,0.0,0.0],[3,1,1])
            indicatorstate1 = preparenext!(gmhrunners[i],indicatorstate1,remotesegments1,indicator1)
            testfrom(indicatorstate1,gv3,gll3)
            setfrom(indicatorstate1,gv1,gll1) ; setproposals(indicatorstate1,gv2,gll2) ; setsegment(remotesegments1,gv2,gll2,gv3,gll3,1) ; setindicator(indicator1,[0.0,0.0,1.0],[3,3,3])
            indicatorstate1 = preparenext!(gmhrunners[i],indicatorstate1,remotesegments1,indicator1)
            testfrom(indicatorstate1,gv1,gll1)

            println("Testing of store! function")
            #test storing a value into the chain
            setfrom(indicatorstate1,gv1,gll1) ; setproposals(indicatorstate1,gv2,gll2) ; setsegment(remotesegments1,gv2,gll2,gv3,gll3,2) ; setindicator(indicator1,[1.0,0.0,0.0],[3,1,2])
            store!(gmhrunners[i],indicatorstate1,remotesegments1,indicator1,chain1)
            teststore(chain1,gv3,gll3[1],2,2,2)
            setfrom(indicatorstate1,gv1,gll1) ; setproposals(indicatorstate1,gv2,gll2) ; setsegment(remotesegments1,gv2,gll2,gv3,gll3,2) ; setindicator(indicator1,[0.0,0.0,1.0],[3,3,3])
            store!(gmhrunners[i],indicatorstate1,remotesegments1,indicator1,chain1)
            teststore(chain1,gv1,gll1[1],4,2,4)

            #reset the samplerstates used below
            indicatorstate1 = samplerstate(rsampler1,1,Float64,Float64)
            indicatorstate2 = samplerstate(rsampler1,1,Float64,Float64)
            println("Testing of auxiliary! function")
            srand(rinitseed1) ; initialize!(gmhrunners[i],rmodel1,indicatorstate1,chain1,false) ; propose!(indicatorstate1) ; geometry!(rmodel1,proposals(indicatorstate1))
            srand(rinitseed1) ; initialize!(gmhrunners[i],rmodel1,indicatorstate2,chain1,false) ; auxiliary!(gmhrunners[i],rmodel1,indicatorstate2)
            testfrom(indicatorstate2,indicatorstate1.from.values,indicatorstate1.from.loglikelihood)
            testproposals(indicatorstate2,indicatorstate1.proposals.values,indicatorstate1.proposals.loglikelihood)

            println("Testing of iterate! function")
            #carry out the individual steps of an interation
            srand(rinitseed1)
            initialize!(gmhrunners[i],rmodel1,indicatorstate1,chain1,false)
            auxiliary!(gmhrunners[i],rmodel1,indicatorstate1)
            segmentacceptances = iterate!(remotesegments1,indicatorstate1)
            transitionprobability!(indicator1,acceptance!(indicatorstate1),segmentacceptances)
            sampleindicator!(indicator1)
            srand(rinitseed1)
            initialize!(gmhrunners[i],rmodel1,indicatorstate2,chain1,false)
            iterate!(gmhrunners[i],rmodel1,indicatorstate2,remotesegments1,indicator1,chain1,false)
            #compare the values from both steps
            testfrom(indicatorstate2,indicatorstate1.from.values,indicatorstate1.from.loglikelihood)
            testproposals(indicatorstate2,indicatorstate1.proposals.values,indicatorstate1.proposals.loglikelihood)

            println("Testing of tune! function")
            accepted!(tunerstate1,indicator1)
            accepted!(tunerstate2,indicator2)
            tvals1 = tune(rtuner1,tunerstate1)
            tvals2 = tune(rtuner2,tunerstate2)
            cov1 = cov(indicatorstate1.density.distribution)

            @test index(tunerstate1) == 1
            tune!(gmhrunners[i],rsampler1,indicatorstate1,remotesegments1,rtuner1,tunerstate1)
            @test index(tunerstate1) == 2 && cov1 == cov(indicatorstate1.density.distribution)

            @test index(tunerstate2) == 1
            tune!(gmhrunners[i],rsampler1,indicatorstate1,remotesegments1,rtuner2,tunerstate2)
            @test index(tunerstate2) == 2
            @test_approx_eq cov1*tvals2*tvals2 cov(indicatorstate1.density.distribution)

            for j=1:remotesegments1.numsegments
                @test_approx_eq cov1*tvals2*tvals2 cov(fetch(remotesegments1.remote[j]).samplerstate.density.distribution)
            end

            #run the whole mcmc
            println("Testing of run! function")
            srand(435) ; gmhc = run!(gmhrunners[i],rmodel1,rsampler1,rtuner1)
            @test numsamples(gmhc) == gmhchainsamples[i]
            println("Results of the run!")
            show(round(samples(gmhc),3)) ; println()
            show(round(loglikelihood(gmhc),3)) ; println()
        end
    catch e
        rmprocs(workers())
        throw(e)
    end
end

rmprocs(workers())
gmhpolicy = policy(:mh,gmhnprops1)
gmhrunner = runner(gmhpolicy,rniter1,gmhnprops1;numburnin = rnburnin1)
println("====================")
println("Test show() function")
show(gmhrunner)
println("End  show() function")
println("====================")
println()

println("==============================")
println("Testing of GMHRunner completed")
println("==============================")

# nothing
