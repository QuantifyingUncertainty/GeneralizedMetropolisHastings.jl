#############Test Metropolis-Hastings samplers####################

for args in [((:mh,:normal,0.1,eye(2)),4),
            ((:mh,:normal,0.1,ones(2)),1),
            ((:mh,:normal,0.1,2),10)]
    s = sampler(args[1]...)
    t = samplerstate(s,args[2],Float64,Float64)
    @test numparas(t.from) == 2
    @test numsamples(t.from) == 1
    @test numparas(t.proposals) == 2
    @test numsamples(t.proposals) == args[2]
    @test cov(t.density.distribution) == 0.1*eye(2)
end

for args in [((:mh,:lognormal,0.1,eye(2)),4),
            ((:mh,:lognormal,0.1,ones(2)),1),
            ((:mh,:lognormal,0.1,2),10)]
    s = sampler(args[1]...)
    t = samplerstate(s,args[2],Float64,Float64)
    @test numparas(t.from) == 2
    @test numsamples(t.from) == 1
    @test numparas(t.proposals) == 2
    @test numsamples(t.proposals) == args[2]
    @test scale(t.density.distribution) == 0.1*eye(2)
end

s1 = sampler(:mh,:normal,0.1,eye(2))
s2 = sampler(:mh,:lognormal,0.1,2)
t1 = samplerstate(s1,1,Float64,Float64)
t2 = samplerstate(s2,5,Float64,Float64)

setfrom!(t1,BaseSample{Float64,Float64,Vector,Array}(0.5*ones(2),zeros(1),zeros(1)))
setfrom!(t2,BaseSample{Float64,Float64,Vector,Array}([0.1,0.2],zeros(1),zeros(1)))
@test t1.from.values == [0.5,0.5]
@test t2.from.values == [0.1,0.2]

srand(46732) ; propose!(t1) ; srand(46732) ; r1 = rand!(Distributions.MvNormal([0.5,0.5],0.1*eye(2)),similar(t1.proposals.values))
srand(46732) ; propose!(t2) ; srand(46732) ; r2 = rand!(Distributions.MvLogNormal(log([0.1,0.2]),0.1*eye(2)),similar(t2.proposals.values))
@test t1.proposals.values == r1
@test t2.proposals.values == r2

copy!(t1.from.loglikelihood,[-1.0])
copy!(t1.from.logprior,[-5.0])
copy!(t1.proposals.loglikelihood,[-0.5])
copy!(t1.proposals.logprior,[-4.0])

copy!(t2.from.loglikelihood,[-1.0])
copy!(t2.from.logprior,[-5.0])
copy!(t2.proposals.loglikelihood,[-0.1,-0.2,-0.3,-0.4,-0.5])
copy!(t2.proposals.logprior,[-4.0,-4.0,-4.0,-4.0,-4.0])

a1 = acceptanceratio!(t1)
a2 = acceptanceratio!(t2)

@test a1 == [t1.proposals.loglikelihood[1] + t1.proposals.logprior[1] - t1.from.loglikelihood[1] - t1.from.logprior[1]]
for i=1:5
    k12 = Distributions.logpdf(t2.density.distribution,t2.proposals.values[:,i])
    k21 = Distributions.logpdf(Distributions.MvLogNormal(log(t2.proposals.values[:,i]),0.1*eye(2)),t2.from.values)
    @test_approx_eq a2[i] (t2.proposals.loglikelihood[i] +t2.proposals.logprior[i] - t2.from.loglikelihood[1] - t2.from.logprior[1] + k12 - k21)
end

@test_approx_eq t1.scalefactor 0.1
@test_approx_eq t2.scalefactor 0.1
tune!(s1,t1,0.4)
tune!(s2,t2,0.9)
@test_approx_eq t1.scalefactor 0.04
@test_approx_eq t2.scalefactor 0.09
@test_approx_eq cov(t1.density.distribution) 0.04*eye(2)
@test_approx_eq scale(t2.density.distribution) 0.09*eye(2)

###test the show functions
println()
println("====================")
println("Test show() function")
println("====================")
show(s1)
show(s2)
show(t1)
show(t2)
println("====================")
println("End  show() function")
println("====================")
println()


