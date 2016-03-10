include("../util.jl")
println("Number of processes running: ",nprocs())

function testnumsamplerstates(states,pol,nprops)
    @test length(states) == min(GeneralizedMetropolisHastings._numsamplerstates(traittype(pol.samplerstates)),nprops) && numsamples(states[1].proposals) == ceil(Int,nprops/length(states))
    i2stuple = GeneralizedMetropolisHastings._ind2subtuple(states)
    @test i2stuple == (numsamples(states[1].proposals),length(states))
end

function testinitialize(model,state,vals,loglikelihood)
    @test state.from.values == vals
    @test state.from.loglikelihood == loglikelihood
end

function testiterate(model,state,vals,loglikelihood)
    @test state.proposals.values == vals
    @test state.proposals.loglikelihood == loglikelihood
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

function teststore(indicatorstate,states,indicator,chain,prevacc,prevprop)
    ind = indicatorsamples(indicator)
    acc = sum(ind[2:end] - ind[1:end-1] .!= 0)
    @test chain.accepted == acc+prevacc && chain.proposed == numsamples(indicator)+prevprop
    for i=2:numsamples(indicator) + 1
        if ind[i] == numproposals(indicator)+1
            @test chain.values[:,i+prevprop-1] == indicatorstate.density.distribution.μ
        else
            s,t = ind2sub(GeneralizedMetropolisHastings._ind2subtuple(states),ind[i])
            @test chain.values[:,i+prevprop-1] == proposals(states[t]).values[:,s]
            @test chain.loglikelihood[i+prevprop-1] == proposals(states[t]).loglikelihood[s]
        end
    end
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

###The same model is used with all test cases below, a two-parameter spring-mass model
defaults1 = [110.0,8.0]
srand(897) ; model1 = springmassmodel(linspace(0.0,10.0,100),[-1.0,1.0],[100.0,10.0],[1e-2,9e-2],defaults1)
###################Test a runner with 1 proposal per iteration (normal Metropolis-Hastings)
nprops1 = 1
nparas1 = 2
niterations1 = 5
policy1 = policy(:gmh,nprops1)
runner1 = runner(:gmh,niterations1,nprops1,policy1)
sampler1 = sampler(:mh,:normal,0.5,[1.0 0.0;0.0 0.1])
chain1 = chain(:standard,nparas1,nprops1*niterations1,policy1.sampletype,policy1.calculationtype)

samplerstates1 = samplerstates(runner1,sampler1)
testnumsamplerstates(samplerstates1,policy1,nprops1)

indicator1 = indicator(traitvalue(runner1.policy.indicator),length(samplerstates1)*numsamples(samplerstates1[1].proposals),nprops1,Float64)
@test numproposals(indicator1) == nprops1 && numsamples(indicator1) == nprops1

srand(345) ; initialize!(runner1,model1,from(samplerstates1[1]))
srand(345) ; loglikelihood1 = [loglikelihood(model1,evaluate!(model1,defaults1))]
testinitialize(model1,samplerstates1[1],defaults1,loglikelihood1)

srand(456) ; samplerstates1[1] = iterate!(runner1,model1,samplerstates1[1])
srand(456) ; proposal1 = rand!(samplerstates1[1].density.distribution,zeros(2)) ; loglikelihood1 = [loglikelihood(model1,evaluate!(model1,proposal1))]
testiterate(model1,samplerstates1[1],proposal1,loglikelihood1)

transitionprobability!(indicator1,acceptanceratio(samplerstates1[1]))
testtransitionprobability(samplerstates1,indicator1)
srand(567) ; sampleindicator!(indicator1)
srand(567) ; indicatorsamples1 = vcat([nprops1+1],[rand(Distributions.Categorical(indicator1.stationary)) for k=1:nprops1])
@test indicatorsamples(indicator1) == indicatorsamples1
store!(runner1,samplerstates1[1],samplerstates1,indicator1,chain1)
teststore(samplerstates1[1],samplerstates1,indicator1,chain1,0,0)
updatefrom!(runner1,samplerstates1[1],samplerstates1,indicator1)
testupdatefrom(samplerstates1[1],samplerstates1,indicator1)
prevacc1 = chain1.accepted
prevprop1 = chain1.proposed

dumpiterate1(samplerstates1,indicator1)

srand(789) ; iterateandstore!(runner1,model1,samplerstates1,indicator1,chain1)
srand(789) ; proposal1 = rand!(samplerstates1[1].density.distribution,zeros(2)) ; loglikelihood1 = [loglikelihood(model1,evaluate!(model1,proposal1))] ; indicatorsamples1 = vcat([nprops1+1],[rand(Distributions.Categorical(indicator1.stationary)) for i=1:nprops1])
testiterate(model1,samplerstates1[1],proposal1,loglikelihood1)
@test indicatorsamples(indicator1) == indicatorsamples1
teststore(samplerstates1[1],samplerstates1,indicator1,chain1,prevacc1,prevprop1)
testupdatefrom(samplerstates1[1],samplerstates1,indicator1)

dumpiterate1(samplerstates1,indicator1)

runiterateandstore1!(runner1,model1,samplerstates1,indicator1,chain1)
runiterateandstore1!(runner1,model1,samplerstates1,indicator1,chain1)
runiterateandstore1!(runner1,model1,samplerstates1,indicator1,chain1)


###################Test a runner with 2 proposals per iteration (Generalized Metropolis-Hastings)
nprops2 = 6
nparas2 = 2
niterations2 = 5
policy2 = policy(:gmh,nprops2;samplerstates=:test)
runner2 = runner(:gmh,niterations2,nprops2,policy2)
chain2 = chain(:standard,nparas2,nprops2*niterations2,policy2.sampletype,policy2.calculationtype)

samplerstates2 = samplerstates(runner2,sampler1)
testnumsamplerstates(samplerstates2,policy2,nprops2)

indicatorstate2 = samplerstate(sampler1,1,Float64,Float64)

indicator2 = indicator(traitvalue(runner2.policy.indicator),length(samplerstates2)*numsamples(samplerstates2[1].proposals),nprops2,Float64)
@test numproposals(indicator2) == nprops2 && numsamples(indicator2) == nprops2

srand(345) ; initialize!(runner2,model1,from(indicatorstate2))
srand(345) ; loglikelihood2 = [loglikelihood(model1,evaluate!(model1,defaults1))]
testinitialize(model1,indicatorstate2,defaults1,loglikelihood2)

srand(8493) ; indicatorstate2 = iterate!(runner2,model1,indicatorstate2)
srand(8493) ; proposal2 = rand!(indicatorstate2.density.distribution,zeros(2)) ; loglikelihood2 = [loglikelihood(model1,evaluate!(model1,proposal2))]
testiterate(model1,indicatorstate2,proposal2,loglikelihood2)

map!((s)->setfrom!(s,proposals(indicatorstate2)),samplerstates2)
for i=1:length(samplerstates2)
    @test proposals(indicatorstate2) == from(samplerstates2[i])
end
srand(7463) ; samplerstates2 = iterate!(runner2,model1,samplerstates2)
transitionprobability!(indicator2,acceptanceratio(indicatorstate2),map(acceptanceratio,samplerstates2))
testtransitionprobability(samplerstates2,indicator2)
srand(5673) ; sampleindicator!(indicator2)
srand(5673) ; indicatorsamples2 = vcat([nprops2+1],[rand(Distributions.Categorical(indicator2.stationary)) for i=1:nprops2])
store!(runner2,indicatorstate2,samplerstates2,indicator2,chain2)
teststore(indicatorstate2,samplerstates2,indicator2,chain2,0,0)
updatefrom!(runner2,indicatorstate2,samplerstates2,indicator2)
testupdatefrom(indicatorstate2,samplerstates2,indicator2)

dumpiterate2(indicatorstate2,samplerstates2,indicator2)

runiterateandstore2!(runner2,model1,indicatorstate2,samplerstates2,indicator2,chain2)
runiterateandstore2!(runner2,model1,indicatorstate2,samplerstates2,indicator2,chain2)
runiterateandstore2!(runner2,model1,indicatorstate2,samplerstates2,indicator2,chain2)
runiterateandstore2!(runner2,model1,indicatorstate2,samplerstates2,indicator2,chain2)

println("#################################################################################################")
println("#################################################################################################")

runner3 = runner(:gmh,niterations1,nprops1,policy1;numburnin=4)
runner4 = runner(:gmh,niterations2,nprops2,policy2;numburnin=4)
tuner3 = tuner(:scale,2,0.5,:logistic)
tuner4 = tuner(:scale,2,0.5,:erf)

chain3 = run!(runner3,model1,sampler1;tuner=tuner3)
chain4 = run!(runner4,model1,sampler1;tuner=tuner4)

println("====================")
println("Test show() function")
show(runner1)
show(runner2)
show(runner3)
show(runner4)
show(chain1)
show(chain2)
show(chain3)
show(chain4)
println("End  show() function")
println("====================")
println()

nothing
