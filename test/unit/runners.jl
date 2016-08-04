println("Number of processes running: ",nprocs())

#Define some quantities and a very simple TargetModel for all tests below
rinitseed1 = 879
rnburnin1 = 2
rniter1 = 4
smhnprops1 = 1 #number of test proposals for standard mh
gmhnprops1 = 2 #number of test proposals for generalized mh
rtime1 = linspace(0.0,10.0,100)
rlower1 = [2.5]
rupper1 = [3.5]
rdefault1 = [3.0]
rcov1 = [0.01]
rparas1 = parameters([:a],rlower1,rupper1,rdefault1)
rdata1 = data(:array,rtime1,sin(3*rtime1))
rnoise1 = noise(:gaussian,rcov1)
rmodel1 = model(:target,rparas1,rdata1,rnoise1,(t,p)->sin(p[1]*t);name="RunnersTestModel")
rnparas1 = numparas(rmodel1)
rsampler1 = sampler(:mh,:normal,0.1,eye(1))
rtuner1 = tuner(:monitor,1)
rtuner2 = tuner(:scale,2,0.5,:logistic)

function testfrom(state,vals,ll)
    @test state.from.values == vals
    @test state.from.loglikelihood == ll
end

function testproposals(state,vals,ll)
    @test state.proposals.values == vals
    @test state.proposals.loglikelihood == ll
end

function teststore(chain,vals,ll,i,a,p)
    @test chain.values[:,i] == vals
    @test chain.loglikelihood[i] == ll
    @test chain.accepted == a
    @test chain.proposed == p
end

function setfrom(state,vals,loglikelihood)
    copy!(state.from.values,vals)
    copy!(state.from.loglikelihood,loglikelihood)
end

function setproposals(state,vals,loglikelihood)
    copy!(state.proposals.values,vals)
    copy!(state.proposals.loglikelihood,loglikelihood)
end

function setindicator(indicator,probs,samples)
    copy!(indicator.stationary,probs)
    copy!(indicator.samples,samples)
end

function setsegment(rsegs,from,ll1,prop,ll2,j)
    l = length(prop)
    p,s = GeneralizedMetropolisHastings._prop2seg(rsegs,j)
    r = rsegs.remote[s]
    @sync begin
        @spawnat r.where copy!(fetch(r).samplerstate.from.values,from)
        @spawnat r.where copy!(fetch(r).samplerstate.from.loglikelihood,ll1)
        @spawnat r.where copy!(fetch(r).samplerstate.proposals.values,l*(p-1)+1,prop,1,l)
        @spawnat r.where copy!(fetch(r).samplerstate.proposals.loglikelihood,p,ll2,1,1)
    end
    rsegs.prop2collected[j] = (p,s)
end

function testtransitionprobability(states,indicator)
    acc = Float64[]
    for x in states
        acc = vcat(acc,acceptanceratio(x))
    end
    acc = vcat(acc,[0.0])
    acc = exp(acc - maximum(acc))
    acc /= sum(acc)
    @test_approx_eq indicator.stationary acc
end



function testupdatefrom(indicatorstate,states,indicator)
    indend = indicatorsamples(indicator)[end]
    if indend != numproposals(indicator) + 1
        s,t = ind2sub(GeneralizedMetropolisHastings._ind2subtuple(states),indend)
        @test from(indicatorstate).values == proposals(states[t]).values[:,s]
        @test from(indicatorstate).loglikelihood[1] == proposals(states[t]).loglikelihood[s]
        @test from(indicatorstate).logprior[1] == proposals(states[t]).logprior[s]
    end
end

function dumpiterate1(samplerstates,indicator)
    dump(samplerstates[1])
    println("===============================================")
    dump(indicator1)
    println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
end

function runiterateandstore1!(runner,model,samplerstates,indicator,chain)
    prevacc = chain.accepted
    prevprop = chain.proposed
    iterateandstore!(runner,model,samplerstates,indicator,chain)
    teststore(samplerstates[1],samplerstates,indicator,chain,prevacc,prevprop)
    testupdatefrom(samplerstates1[1],samplerstates1,indicator1)
    dumpiterate1(samplerstates,indicator)
end

function dumpiterate2(indicatorstate,samplerstates,indicator)
    dump(indicatorstate)
    println("000000000000000000000000000000000000000000000000000000000000000000000000000")
    dump(samplerstates[1])
    println("111111111111111111111111111111111111111111111111111111111111111111111111111")
    dump(samplerstates[2])
    println("222222222222222222222222222222222222222222222222222222222222222222222222222")
    dump(samplerstates[3])
    println("333333333333333333333333333333333333333333333333333333333333333333333333333")
    dump(indicator)
    println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
end

function runiterateandstore2!(runner,model,indicatorstate,samplerstates,indicator,chain)
    prevacc = chain.accepted
    prevprop = chain.proposed
    iterateandstore!(runner,model,indicatorstate,samplerstates,indicator,chain)
    teststore(indicatorstate,samplerstates,indicator,chain,prevacc,prevprop)
    testupdatefrom(indicatorstate,samplerstates,indicator)
    dumpiterate2(indicatorstate,samplerstates,indicator)
end

