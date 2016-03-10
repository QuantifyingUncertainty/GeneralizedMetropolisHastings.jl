println("======================")
println("Tests of NormalDensity")
println("======================")

@test NormalDensity <: SymmetricDensity

###Test construction
d1 = density(:normal,ones(2),0.1*eye(2))
@test isa(d1,NormalDensity)
@test issymmetric(d1)
@test mean(d1.distribution) == ones(2)
@test cov(d1.distribution) == 0.1*eye(2)

###Test conditioning of the distribution
condition!(d1,2*ones(2))
@test mean(d1.distribution) == 2*ones(2)
@test cov(d1.distribution) == 0.1*eye(2)

condition!(d1,3*ones(2),0.2*eye(2))
@test mean(d1.distribution) == 3*ones(2)
@test cov(d1.distribution) == 0.2*eye(2)

@test_throws MethodError condition!(d1,ones(2),d1)

###Test proposing and calculatinng the logprobability
for i in [1,3]
    s1 = samples(:base,2,i,Float64,Float64)
    z1 = zeros(Float64,GeneralizedMetropolisHastings._valuestuple(2,i))
    @test (srand(i+23) ; propose!(d1,s1) ; s1.values) == (srand(i+23) ; rand!(d1.distribution,z1))
    p1 = Distributions.logpdf(d1.distribution,s1.values) ; isa(p1,Number)?collect(p1):p1
    @test_approx_eq logprobability(d1,s1.values) p1
end


#################################################################################

println("=========================")
println("Tests of LogNormalDensity")
println("=========================")

@test LogNormalDensity <: ASymmetricDensity

###Test construction
d2 = density(:lognormal,exp(ones(2)),0.1*eye(2))
@test isa(d2,LogNormalDensity)
@test !issymmetric(d2)
@test_approx_eq Distributions.median(d2.distribution) exp(ones(2))
@test_approx_eq Distributions.location(d2.distribution) ones(2)
@test_approx_eq Distributions.scale(d2.distribution) 0.1*eye(2)

###Test conditioning of the distribution
condition!(d2,exp(2*ones(2)))
@test_approx_eq Distributions.median(d2.distribution) exp(2*ones(2))
@test_approx_eq Distributions.location(d2.distribution) 2*ones(2)
@test_approx_eq Distributions.scale(d2.distribution) 0.1*eye(2)

condition!(d2,exp(ones(2)),0.01*eye(2))
@test_approx_eq Distributions.location(d2.distribution) ones(2)
@test_approx_eq Distributions.scale(d2.distribution) 0.01*eye(2)

###Test proposing and calculatinng the logprobability
for i in [1,3]
    s2 = samples(:base,2,i,Float64,Float64)
    z2 = zeros(Float64,GeneralizedMetropolisHastings._valuestuple(2,i))
    @test (srand(i+23) ; propose!(d2,s2) ; s2.values) == (srand(i+23) ; rand!(d2.distribution,z2))
    p2 = Distributions.logpdf(d2.distribution,s2.values) ; isa(p2,Number)?collect(p2):p2
    @test_approx_eq logprobability(d2,s2.values) p2
end

println("============================")
println("Tests of DistributionWrapper")
println("============================")

###Test construction
d3 = density(Distributions.Uniform(0.0,10.0);symmetric = true)
d4 = density(Distributions.LogNormal(1.0,0.1);symmetric = false)
d5 = density(Distributions.Geometric(0.3))
@test isa(d3,DistributionWrapper)
@test isa(d4,DistributionWrapper)
@test isa(d5,DistributionWrapper)
@test issymmetric(d3)
@test !issymmetric(d4)
@test !issymmetric(d5)

###Test proposing
srand(467) ; s3 = propose!(d3,samples(:base,1,1,Float64,Float64)) ; s4 = propose!(d4,samples(:base,1,2,Float64,Float64)) ; s5 = propose!(d5,samples(:base,1,3,Float64,Float64))
srand(467) ; r3 = rand!(d3.distribution,zeros(1)) ; r4 = rand!(d4.distribution,zeros(1,2)) ; r5 = rand!(d5.distribution,zeros(1,3))
@test s3.values == r3
@test s4.values == r4
@test s5.values == r5

###Test calculating the logpdf
@test logprobability(d3,s3.values) == Distributions.logpdf!(zeros(1),d3.distribution,s3.values)
@test logprobability(d4,s4.values) == Distributions.logpdf!(zeros(2),d4.distribution,s4.values)
@test logprobability(d5,s5.values) == Distributions.logpdf!(zeros(3),d5.distribution,s5.values)

###test the show functions
println()
println("====================")
println("Test show() function")
println("====================")
show(d1)
show(d2)
show(d3)
show(d4)
show(d5)
println("====================")
println("End  show() function")
println("====================")
println()

nothing
