println("======================")
println("Tests of NormalDensity")
println("======================")

@test NormalDensity <: SymmetricDensity

###Test construction
for T in [Float64,Float32]
    println("Testing NormalDensity with ",T)
    n1 = distribution(:normal,ones(T,2),0.1*eye(T,2))
    d1 = density(:normal,ones(T,2),0.1*eye(T,2))
    @test isa(d1,NormalDensity)
    @test issymmetric(d1)
    @test mean(d1.distribution) == mean(n1) && cov(d1.distribution) == cov(n1)

    ###Test conditioning of the distribution
    n1 = recenter(n1,2*ones(T,2))
    condition!(d1,2*ones(T,2))
    @test mean(d1.distribution) == mean(n1) && cov(d1.distribution) == cov(n1)

    n1 = rescale(n1,T(1/2))
    scale!(d1,T(1/2))
    @test mean(d1.distribution) == mean(n1) && cov(d1.distribution) == cov(n1)

    n1 = update(n1,3*ones(T,2),3*eye(T,2))
    update!(d1,3*ones(T,2),3*eye(T,2))
    @test mean(d1.distribution) == mean(n1) && cov(d1.distribution) == cov(n1)

    ###Test proposing and calculatinng the logprobability
    for i in [1,3]
        s1 = samples(:base,2,i,T,T)
        z1 = zeros(T,GeneralizedMetropolisHastings._valuestuple(2,i))
        @test (srand(i+23) ; propose!(d1,s1.values)) == (srand(i+23) ; z1 = rand!(n1,z1))
        p1 = Distributions.logpdf(d1.distribution,s1.values) ; isa(p1,Number)?collect(p1):p1
        @test_approx_eq logprobability(d1,s1.values) p1
    end

    ###Test the normal distribution with zero mean
    n2 = distribution(:normal,0.1*eye(T,2))
    d2 = density(:normal,0.1*eye(T,2))
    @test isa(d2,NormalDensity)
    @test issymmetric(d2)
    @test mean(d2.distribution) == mean(n2) && cov(d2.distribution) == cov(n2)
    @test_throws MethodError condition!(d2,zeros(T,2))
    n2 = rescale(n2,T(2.0))
    scale!(d2,T(2.0))
    @test mean(d2.distribution) == mean(n2) && cov(d2.distribution) == cov(n2)
    n2 = update(n2,eye(T,2))
    update!(d2,eye(T,2))
    @test mean(d2.distribution) == mean(n2) && cov(d2.distribution) == cov(n2)
end

#################################################################################

println("=========================")
println("Tests of LogNormalDensity")
println("=========================")

@test LogNormalDensity <: ASymmetricDensity

###Test construction
l1 = distribution(:lognormal,ones(2),0.1*eye(2))
d2 = density(:lognormal,ones(2),0.1*eye(2))
@test isa(d2,LogNormalDensity)
@test !issymmetric(d2)
@test Distributions.location(d2.distribution) == Distributions.location(l1)
@test Distributions.scale(d2.distribution) == Distributions.scale(l1)

###Test conditioning of the distribution
l1 = recenter(l1,2*ones(2))
condition!(d2,2*ones(2))
@test Distributions.location(d2.distribution) == Distributions.location(l1)
@test Distributions.scale(d2.distribution) == Distributions.scale(l1)

l1 = rescale(l1,1/2)
scale!(d2,1/2)
@test Distributions.location(d2.distribution) == Distributions.location(l1)
@test Distributions.scale(d2.distribution) == Distributions.scale(l1)

l1 = update(l1,3*ones(2),3*eye(2))
update!(d2,3*ones(2),3*eye(2))
@test Distributions.location(d2.distribution) == Distributions.location(l1)
@test Distributions.scale(d2.distribution) == Distributions.scale(l1)

###Test proposing and calculatinng the logprobability
for i in [1,3]
    s2 = samples(:base,2,i,Float64,Float64)
    z2 = zeros(Float64,GeneralizedMetropolisHastings._valuestuple(2,i))
    @test (srand(i+23) ; propose!(d2,s2.values)) == (srand(i+23) ; rand!(d2.distribution,z2))
    p2 = Distributions.logpdf(d2.distribution,s2.values) ; isa(p2,Number)?collect(p2):p2
    @test_approx_eq logprobability(d2,s2.values) p2
end

println("====================================")
println("Tests of CompoundDistributionWrapper")
println("====================================")

for p in [(:uniform,()),
          (:laplace,()),
          (:triangular,()),
          (:bactrian,(:normal,0.95)),
          (:bactrian,(:laplace,0.90))]
    println(p)
    dis1 = distributions(p[1],[1.0,2.0],[1.0,0.5],p[2]...)
    den1 = density(p[1],[1.0,2.0],[1.0,0.5],p[2]...)
    @test issymmetric(den1)
    for k=1:length(dis1)
        @test mean(dis1[k]) == mean(den1.distributions[k]) && scale(dis1[k]) == scale(den1.distributions[k])
    end
    condition!(den1,2*ones(2))
    for k=1:length(dis1)
        dis1[k] = recenter(dis1[k],2.0)
        @test mean(dis1[k]) == mean(den1.distributions[k]) && scale(dis1[k]) == scale(den1.distributions[k])
    end
    scale!(den1,2.0)
    for k=1:length(dis1)
        dis1[k] = rescale(dis1[k],2.0)
        @test mean(dis1[k]) == mean(den1.distributions[k]) && scale(dis1[k]) == scale(den1.distributions[k])
    end
    update!(den1,3*ones(2),3.0*ones(2))
    for k=1:length(dis1)
        dis1[k] = update(dis1[k],3.0,3.0)
        @test mean(dis1[k]) == mean(den1.distributions[k]) && scale(dis1[k]) == scale(den1.distributions[k])
    end
    for i in [1,3]
        s3 = samples(:base,2,i,Float64,Float64)
        z3 = zeros(Float64,GeneralizedMetropolisHastings._valuestuple(2,i))
        p3 = zeros(i)
        @test (srand(i+23) ; propose!(den1,s3.values)) == (srand(i+23) ; for k=1:i for j=1:length(dis1) z3[j,k] = rand(dis1[j]) end end ; z3)
        for k=1:i for j=1:length(dis1) p3[k] += Distributions.logpdf(dis1[j],s3.values[j,k]) end end
        @test_approx_eq logprobability(den1,s3.values) p3
    end
end



###test the show functions
println()
println("====================")
println("Test show() function")
println("====================")
show(d1)
show(d2)
println("====================")
println("End  show() function")
println("====================")
println()

nothing
