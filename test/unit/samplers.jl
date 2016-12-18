#############Test Metropolis-Hastings samplers####################

fv1 = [0.1,0.2]
fv2 = [0.3,0.4]
fr1 = samples(:base,2,1,Float64,Float64) ; copy!(fr1.values,fv1)
fr2 = samples(:base,2,2,Float64,Float64) ; copy!(fr2.values,3,fv2,1,2)

for args in [((:mh,:normal,0.1,eye(2)),4,"Scalable Metropolis-Hastings","Metropolis"),
             ((:mh,:normal,0.2,ones(2)),2,"Scalable Metropolis-Hastings","Metropolis"),
             ((:mh,:lognormal,0.1,2),1,"Scalable Metropolis-Hastings","Metropolis-Hastings"),
             ((:mh,:bactrian,0.1,0.5*ones(2),:normal,0.95),3,"Scalable Metropolis-Hastings","Metropolis")]
    s = sampler(args[1]...)
    t = samplerstate(s,args[2],Float64,Float64)
    @test numparas(t.from) == 2
    @test numsamples(t.from) == 1
    @test numparas(t.proposals) == 2
    @test numsamples(t.proposals) == args[2]
    @test numparas(t) == 2
    @test numsamples(t) == args[2]
    @test from(t) == t.from
    @test proposals(t) == t.proposals
    @test acceptance(t) == zeros(args[2])
    @test GeneralizedMetropolisHastings.samplername(s) == args[3]
    @test GeneralizedMetropolisHastings.samplerstatename(t) == args[4]
    setfrom!(t,fr1)
    @test from(t).values == fv1
    setfrom!(t,fr2,1)
    @test from(t).values == zeros(2)
    setfrom!(t,fr2,2)
    @test from(t).values == fv2
end

function testprepare(t,f1,f2,x1,x2,s)
    @test t.from.values == x1
    @test_approx_eq f1(t.density.distribution) x2
    @test_approx_eq f2(t.density.distribution) s
end

function testsym(t,a)
    f = t.from.loglikelihood[1] + t.from.logprior[1]
    @test a == t.proposals.loglikelihood + t.proposals.logprior - f
end

function testasym(t,a)
    f = t.from.loglikelihood[1] + t.from.logprior[1]
    for i=1:length(a)
        d21 = deepcopy(t.density.distribution)
        d21 = recenter(d21,t.proposals.values[:,i])
        k12 = Distributions.logpdf(t.density.distribution,t.proposals.values[:,i])
        k21 = Distributions.logpdf(d21,t.from.values)
        @test_approx_eq a[i] (t.proposals.loglikelihood[i] +t.proposals.logprior[i] -f + k12 - k21)
    end
end

function testsmmala(t,a)
    f = t.from.loglikelihood[1] + t.from.logprior[1]
    for i=1:length(a)
        d21 = deepcopy(t.density.distribution)
        d21 = update(d21,GeneralizedMetropolisHastings._meancov(t,i)...)
        k12 = Distributions.logpdf(t.density.distribution,t.proposals.values[:,i])
        k21 = Distributions.logpdf(d21,t.from.values)
        @test_approx_eq a[i] (t.proposals.loglikelihood[i] +t.proposals.logprior[i] -f + k12 - k21)
    end
end


#test prepare! functions
for p in [(:normal,mean,cov,identity,testsym),
          (:lognormal,Distributions.location,Distributions.scale,log,testasym)]
    s = sampler(:mh,p[1],0.1,eye(2))
    tind = samplerstate(s,1,Float64,Float64)
    taux = samplerstate(s,3,Float64,Float64)
    copy!(tind.proposals.values,fv1)
    copy!(taux.proposals.values,3,2fv1,1,2)
    #test the prepare! functions
    prepare!(tind)
    testprepare(tind,p[2],p[3],zeros(2),p[4](zeros(2)),0.01eye(2))
    prepare!(tind,true)
    testprepare(tind,p[2],p[3],fv1,p[4](fv1),0.01eye(2))
    prepareauxiliary!(tind,taux)
    testprepare(taux,p[2],p[3],fv1,p[4](fv1),0.01eye(2))
    prepareindicator!(tind)
    testprepare(tind,p[2],p[3],fv1,p[4](fv1),0.01eye(2))
    prepareindicator!(tind,taux,2)
    testprepare(tind,p[2],p[3],2fv1,p[4](2fv1),0.01eye(2))
    prepare!(tind,true)
    testprepare(tind,p[2],p[3],fv1,p[4](fv1),0.01eye(2))
    #test the propose! function
    d = distribution(p[1],fv1,eye(2)) ; d = rescale(d,0.1)
    srand(46732) ; propose!(tind)
    srand(46732) ; r1 = rand!(d,similar(tind.proposals.values))
    srand(46733) ; propose!(taux)
    srand(46733) ; r2 = rand!(d,similar(taux.proposals.values))
    @test tind.proposals.values == r1
    @test taux.proposals.values == r2
    #test the acceptance function
    copy!(tind.from.loglikelihood,[-1.0])
    copy!(tind.from.logprior,[-5.0])
    copy!(tind.proposals.loglikelihood,[-0.5])
    copy!(tind.proposals.logprior,[-4.0])
    copy!(taux.from.loglikelihood,[-1.0])
    copy!(taux.from.logprior,[-5.0])
    copy!(taux.proposals.loglikelihood,[-0.1,-0.2,-0.3])
    copy!(taux.proposals.logprior,[-3.0,-3.0,-3.0])
    a1 = acceptance!(tind)
    a2 = acceptance!(taux)
    @test length(a1) == 1
    @test length(a2) == 3
    p[5](tind,a1)
    p[5](taux,a2)
    tune!(tind,2.0)
    tune!(taux,0.5)
    @test_approx_eq p[3](tind.density.distribution) 0.04*eye(2)
    @test_approx_eq p[3](taux.density.distribution) 0.0025*eye(2)
    println("Test show() function for ",p[1])
    show(s)
    show(tind)
    show(taux)
    println("======================================")
end

#############Test Adaptive Metropolis samplers####################

ru1 = GeneralizedMetropolisHastings._runningstate(Val{:normal},2,0.1)
v = [1.0,2.0]
@test ru1.iteration == 0
@test ru1.runningmean == zeros(2)
@test ru1.runningcov == zeros(2,2)
@test ru1.scale == 0.1
@test size(ru1.tempdiff) == (2,)
@test size(ru1.densitycov) == (2,2)

@test GeneralizedMetropolisHastings._covupdatefactor(1,Float64) == 0.0
@test GeneralizedMetropolisHastings._covupdatefactor(2,Float64) == 0.0
@test GeneralizedMetropolisHastings._covupdatefactor(3,Float64) == 0.5

@test_approx_eq GeneralizedMetropolisHastings._updatetempdiff!(zeros(2),v/3,v) 2v/3
@test_approx_eq GeneralizedMetropolisHastings._updaterunningcov!(eye(2),2v/3,3) eye(2)/2+2v/3*2v'/3/3
@test_approx_eq GeneralizedMetropolisHastings._updaterunningmean!(v/3,v,3) (2*v/3+v)/3
@test_approx_eq GeneralizedMetropolisHastings._updatedensitycov!(zeros(2,2),ones(2,2),0.1) ones(2,2)+0.1*0.1*eye(2)

niterations = 10
rv = rand!(zeros(2,niterations))
for i=1:niterations
    GeneralizedMetropolisHastings.updaterunning!(ru1,rv[:,i])
end

#test showing that the covariance is indeed the correct sample covariance
@test ru1.iteration == niterations
@test_approx_eq ru1.runningmean mean(rv,2)
@test_approx_eq ru1.runningcov cov(rv')'

as1 = sampler(:adaptive,0.1,2)
atind = samplerstate(as1,1,Float64,Float64)
ataux = samplerstate(as1,3,Float64,Float64,true)

@test numparas(as1) == 2
@test GeneralizedMetropolisHastings.samplername(as1) == "Adaptive Normal Metropolis"
@test GeneralizedMetropolisHastings.samplerstatename(atind) == "Adaptive Normal Metropolis"
@test isnull(atind.runningstate) == false
@test isnull(ataux.runningstate) == true

copy!(atind.from.values,v)
copy!(ataux.from.values,v)
GeneralizedMetropolisHastings.updaterunning!(atind)
GeneralizedMetropolisHastings.updaterunning!(ataux)
@test mean(atind.density.distribution) == zeros(2)
@test_approx_eq cov(atind.density.distribution) get(atind.runningstate).densitycov
@test mean(ataux.density.distribution) == zeros(2)
@test_approx_eq cov(ataux.density.distribution) 0.01*eye(2)

copy!(atind.proposals.values,2.5v)
copy!(ataux.proposals.values,3,3.5v,1,2)
#test the prepare! functions
prepare!(atind)
@test get(atind.runningstate).iteration == 2
testprepare(atind,mean,cov,v,zeros(2),0.01eye(2))
prepare!(atind,true)
@test get(atind.runningstate).iteration == 3
testprepare(atind,mean,cov,2.5v,zeros(2),get(atind.runningstate).densitycov)
prepareauxiliary!(atind,ataux)
testprepare(ataux,mean,cov,2.5v,zeros(2),get(atind.runningstate).densitycov)
prepareindicator!(atind)
@test get(atind.runningstate).iteration == 4
testprepare(atind,mean,cov,2.5v,zeros(2),get(atind.runningstate).densitycov)
prepareindicator!(atind,ataux,2)
@test get(atind.runningstate).iteration == 5
testprepare(atind,mean,cov,3.5v,zeros(2),get(atind.runningstate).densitycov)
#test the propose! function
dind = distribution(:normal,atind.from.values,cov(atind.density.distribution))
daux = distribution(:normal,ataux.from.values,cov(ataux.density.distribution))
srand(46732) ; propose!(atind)
srand(46732) ; r1 = rand!(dind,similar(atind.proposals.values))
srand(46733) ; propose!(ataux)
srand(46733) ; r2 = rand!(daux,similar(ataux.proposals.values))
@test atind.proposals.values == r1
@test ataux.proposals.values == r2
#test the acceptance function
copy!(atind.from.loglikelihood,[-1.0])
copy!(atind.from.logprior,[-5.0])
copy!(atind.proposals.loglikelihood,[-0.5])
copy!(atind.proposals.logprior,[-4.0])
copy!(ataux.from.loglikelihood,[-1.0])
copy!(ataux.from.logprior,[-5.0])
copy!(ataux.proposals.loglikelihood,[-0.1,-0.2,-0.3])
copy!(ataux.proposals.logprior,[-3.0,-3.0,-3.0])
a1 = acceptance!(atind)
a2 = acceptance!(ataux)
@test length(a1) == 1
@test length(a2) == 3
testsym(atind,a1)
testsym(ataux,a2)
tune!(atind,2.0)
tune!(ataux,0.5)
@test get(atind.runningstate).scale == 0.2
println("Test show() function for ",:adaptive)
show(as1)
show(atind)
show(ataux)
println("======================================")

########################
### Tests of SmMALA
for p in [(false,:full,:standard,:full),
          (true,:full,:tr,:full),
          (true,:tangent,:tr,:tangent)]
    smpol = GeneralizedMetropolisHastings.SmMALAPolicy(p[1],p[2])
    @test traitvalue(smpol.smmala) == p[3]
    @test traitvalue(smpol.tensor) == p[4]
end

@test_approx_eq GeneralizedMetropolisHastings._adddiag!(eye(2),0.5) 1.5*eye(2)
@test_approx_eq GeneralizedMetropolisHastings._addfrom!(ones(2),0.5ones(2)) 1.5*ones(2)

from1 = 0.5ones(2)
grad1 = [-0.5,0.3]
ten1 = [1.4 1.2;1.2 1.4]
m1,c1 = GeneralizedMetropolisHastings._meancov(Val{:standard},from1,grad1,ten1,0.1)
m2,c2 = GeneralizedMetropolisHastings._meancov(Val{:tr},from1,grad1,ten1,0.1)
@test_approx_eq m1 from1+ten1\grad1*0.1*0.1/2
@test_approx_eq m2 from1+(ten1+eye(2)/0.1)\grad1
@test_approx_eq c1 inv(ten1)*0.01
@test_approx_eq c2 inv(ten1+eye(2)/0.1)

for p in [(false,"Simplified mMALA",:standard,:full),
          (true,"Trust Region Simplified mMALA",:tr,:full),
          (10,"Trust Region Simplified mMALA with Tangent Tensor",:tr,:tangent)]
    d = density(:normal,zeros(2),eye(2))
    sm = sampler(:smmala,0.1,2,p[1])
    smind = samplerstate(sm,1,Float64,Float64)
    smaux = samplerstate(sm,3,Float64,Float64)

    @test numparas(sm) == 2
    @test GeneralizedMetropolisHastings.samplername(sm) == p[2]
    @test GeneralizedMetropolisHastings.samplerstatename(smind) == p[2]

    @test traitvalue(smind.policy.smmala) == p[3]
    @test traitvalue(smind.policy.tensor) == p[4]

    copy!(smind.from.values,from1)
    copy!(smind.from.gradloglikelihood,grad1)
    copy!(smind.from.tensorloglikelihood,ten1)
    copy!(smind.proposals.values,2from1)
    copy!(smind.proposals.gradloglikelihood,grad1)
    copy!(smind.proposals.tensorloglikelihood,ten1)
    copy!(smaux.from.values,3from1)
    copy!(smaux.from.gradloglikelihood,grad1)
    copy!(smaux.from.tensorloglikelihood,ten1)
    copy!(smaux.proposals.values,[from1 2from1 3from1])
    copy!(smaux.proposals.gradloglikelihood,[grad1 0.5grad1 0.3grad1])
    copy!(smaux.proposals.tensorloglikelihood,cat(3,3ten1,1.5ten1,ten1))

    @test GeneralizedMetropolisHastings._meancov(smind) == GeneralizedMetropolisHastings._meancov(Val{p[3]},from1,grad1,ten1,0.1)
    @test GeneralizedMetropolisHastings._meancov(smind,1) == GeneralizedMetropolisHastings._meancov(Val{p[3]},2from1,grad1,ten1,0.1)
    @test GeneralizedMetropolisHastings._meancov(smaux) == GeneralizedMetropolisHastings._meancov(Val{p[3]},3from1,grad1,ten1,0.1)
    @test GeneralizedMetropolisHastings._meancov(smaux,2) == GeneralizedMetropolisHastings._meancov(Val{p[3]},2from1,0.5grad1,1.5ten1,0.1)

    GeneralizedMetropolisHastings.update!(d,GeneralizedMetropolisHastings._meancov(Val{p[3]},from1,grad1,ten1,0.1)...)
    GeneralizedMetropolisHastings._updatedensity!(smind)
    @test_approx_eq mean(smind.density.distribution) mean(d.distribution)
    @test_approx_eq cov(smind.density.distribution) cov(d.distribution)
    GeneralizedMetropolisHastings.update!(d,GeneralizedMetropolisHastings._meancov(Val{p[3]},2from1,0.5grad1,1.5ten1,0.1)...)
    GeneralizedMetropolisHastings._updatedensity!(smaux,2)
    @test_approx_eq mean(smaux.density.distribution) mean(d.distribution)
    @test_approx_eq cov(smaux.density.distribution) cov(d.distribution)

    #test the prepare! functions
    prepare!(smind)
    testprepare(smind,mean,cov,from1,GeneralizedMetropolisHastings._meancov(Val{p[3]},from1,grad1,ten1,0.1)...)
    prepare!(smind,true)
    testprepare(smind,mean,cov,2from1,GeneralizedMetropolisHastings._meancov(Val{p[3]},2from1,grad1,ten1,0.1)...)
    prepareauxiliary!(smind,smaux)
    testprepare(smaux,mean,cov,2from1,GeneralizedMetropolisHastings._meancov(Val{p[3]},2from1,grad1,ten1,0.1)...)
    prepareindicator!(smind)
    testprepare(smind,mean,cov,2from1,GeneralizedMetropolisHastings._meancov(Val{p[3]},2from1,grad1,ten1,0.1)...)
    prepareindicator!(smind,smaux,2)
    testprepare(smind,mean,cov,2from1,GeneralizedMetropolisHastings._meancov(Val{p[3]},2from1,0.5grad1,1.5ten1,0.1)...)
    prepare!(smind,true)
    testprepare(smind,mean,cov,2from1,GeneralizedMetropolisHastings._meancov(Val{p[3]},2from1,grad1,ten1,0.1)...)

    #test the propose! function
    dind = distribution(:normal,GeneralizedMetropolisHastings._meancov(Val{p[3]},2from1,grad1,ten1,0.1)...)
    daux = distribution(:normal,GeneralizedMetropolisHastings._meancov(Val{p[3]},2from1,grad1,ten1,0.1)...)
    srand(46732) ; propose!(smind)
    srand(46732) ; r1 = rand!(dind,similar(smind.proposals.values))
    srand(46733) ; propose!(smaux)
    srand(46733) ; r2 = rand!(daux,similar(smaux.proposals.values))
    @test smind.proposals.values == r1
    @test smaux.proposals.values == r2

    #test the acceptance function
    copy!(smind.from.loglikelihood,[-1.0])
    copy!(smind.from.logprior,[-5.0])
    copy!(smind.proposals.loglikelihood,[-0.5])
    copy!(smind.proposals.logprior,[-4.0])
    copy!(smaux.from.loglikelihood,[-1.0])
    copy!(smaux.from.logprior,[-5.0])
    copy!(smaux.proposals.loglikelihood,[-0.1,-0.2,-0.3])
    copy!(smaux.proposals.logprior,[-3.0,-3.0,-3.0])
    a1 = acceptance!(smind)
    a2 = acceptance!(smaux)
    @test length(a1) == 1
    @test length(a2) == 3
    testsmmala(smind,a1)
    testsmmala(smaux,a2)

    println("Test show() function for smmala $(p[3]) $(p[4])")
    show(sm)
    show(smind)
    show(smaux)
    println("======================================")
end






