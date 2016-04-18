println("Number of processes running: ",nprocs())

srand(897)
defaults1 = [110.0,8.0]
init1 = [-1.0,1.0]
paras1 = [100.0,10.0]
noise1 = [1e-2,9e-2]
model1 = springmassmodel(linspace(0.0,10.0,100),init1,paras1,noise1,defaults1)
nparas1 = numparas(model1)

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

function setsegment(segment,from,ll1,prop,ll2)
    r = segment.remote[1]
    @sync begin
        @spawnat r.where copy!(fetch(r).samplerstate.from.values,from)
        @spawnat r.where copy!(fetch(r).samplerstate.from.loglikelihood,ll1)
        @spawnat r.where copy!(fetch(r).samplerstate.proposals.values,1,prop,1,length(prop))
        @spawnat r.where copy!(fetch(r).samplerstate.proposals.loglikelihood,1,ll2,1,1)
    end
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

