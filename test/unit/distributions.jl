_getlocation(d::Distributions.LogNormal) = d.μ
_getlocation(d::Distributions.MvLogNormal) = Distributions.location(d)
_getscale(d::Distributions.Normal) = var(d)
_getscale(d::Distributions.MvNormal) = cov(d)
_getscale(d::Distributions.LogNormal) = d.σ^2
_getscale(d::Distributions.MvLogNormal) = Distributions.scale(d)

#testing of the normal distribution
for p in [(1.0,2.0,4.0,2.0),
          (ones(2),4.0*eye(2),4.0*eye(2),2.0*ones(2)),
          (ones(2),2.0*ones(2),4.0*eye(2),2.0*ones(2))]
    cp2 = copy(p[2])
    n = distribution(:normal,p[1],p[2])
    @test_approx_eq mean(n) p[1]
    n = recenter(n,p[4])
    @test_approx_eq mean(n) p[4]
    @test_approx_eq _getscale(n) p[3]
    n = rescale(n,2.0)
    @test_approx_eq _getscale(n) 2.0*2.0*p[3]
    n = update(n,p[1],cp2)
    @test_approx_eq mean(n) p[1]
    @test_approx_eq _getscale(n) p[3]
end

for p in [(eye(2),4.0*eye(2)),
          (ones(2),4.0*eye(2))]
    n = distribution(:normal,eye(2))
    @test mean(n) == zeros(2) && cov(n) == eye(2)
    @test_throws MethodError recenter(n,ones(2))
    n = rescale(n,2.0)
    @test_approx_eq cov(n) 4.0*eye(2)
    n = update(n,p[1])
    @test_approx_eq cov(n) eye(2)
end

#testing of the lognormal distribution
for p in [(1.0,2.0,4.0,2.0),
          (ones(2),4.0*eye(2),4.0*eye(2),2.0*ones(2)),
          (ones(2),2.0*ones(2),4.0*eye(2),2.0*ones(2))]
    cp2 = copy(p[2])
    l1 = distribution(:lognormal,p[1],p[2])
    @test_approx_eq _getlocation(l1) log(p[1])
    l1 = recenter(l1,p[4])
    @test_approx_eq _getlocation(l1) log(p[4])
    @test_approx_eq _getscale(l1) p[3]
    l1 = rescale(l1,2.0)
    @test_approx_eq _getscale(l1) 4.0*p[3]
    l1 = update(l1,p[1],cp2)
    @test_approx_eq _getlocation(l1) log(p[1])
    @test_approx_eq _getscale(l1) p[3]
end

#testing various uniform distributions
for s in [:uniform,:laplace,:triangular]
    l1 = distribution(s,1.0,0.5)
    @test mean(l1) == 1.0 && scale(l1) == 0.5
    l1 = recenter(l1,2.0)
    @test mean(l1) == 2.0 && scale(l1) == 0.5
    l1 = rescale(l1,3.0)
    @test mean(l1) == 2.0 && scale(l1) == 1.5
    l1 = update(l1,1.0,1.0)
    @test mean(l1) == 1.0 && scale(l1) == 1.0
end

#testing the bactrian distribution
for p in [(:normal,1.0),
          (:laplace,sqrt(2.0)),
          (:triangular,6.0)]
    bs1 = GeneralizedMetropolisHastings._bactrianscale(Val{p[1]},0.5,0.95)
    bs2 = GeneralizedMetropolisHastings._bactrianscale(Val{p[1]},1.0,0.95)
    @test_approx_eq bs1 0.5*sqrt(p[2]*(1-0.95*0.95))
    @test_approx_eq bs2 1.0*sqrt(p[2]*(1-0.95*0.95))
    m1 = GeneralizedMetropolisHastings._mixturemodel(p[1],2.0,0.5,0.95)
    @test m1.components[1] == distribution(p[1],2.0-0.5*0.95,bs1)
    @test m1.components[2] == distribution(p[1],2.0+0.5*0.95,bs1)
    b1 = distribution(:bactrian,2.0,0.5,p[1],0.95)
    @test mean(b1.mixture) == 2.0 && b1.μ == 2.0 && b1.σ == 0.5 && b1.m == 0.95
    @test_approx_eq Distributions.params(b1.mixture.components[1])[1] 2.0-0.5*0.95
    @test_approx_eq Distributions.params(b1.mixture.components[2])[1] 2.0+0.5*0.95
    @test_approx_eq Distributions.params(b1.mixture.components[1])[2] bs1
    @test_approx_eq Distributions.params(b1.mixture.components[2])[2] bs1
    b1 = recenter(b1,1.0)
    @test mean(b1.mixture) == 1.0 && b1.μ == 1.0 && b1.σ == 0.5 && b1.m == 0.95
    @test_approx_eq Distributions.params(b1.mixture.components[1])[1] 1.0-0.5*0.95
    @test_approx_eq Distributions.params(b1.mixture.components[2])[1] 1.0+0.5*0.95
    @test_approx_eq Distributions.params(b1.mixture.components[1])[2] bs1
    @test_approx_eq Distributions.params(b1.mixture.components[2])[2] bs1
    b1 = rescale(b1,2.0)
    @test mean(b1.mixture) == 1.0 && b1.μ == 1.0 && b1.σ == 1.0 && b1.m == 0.95
    @test_approx_eq Distributions.params(b1.mixture.components[1])[1] 1.0-1.0*0.95
    @test_approx_eq Distributions.params(b1.mixture.components[2])[1] 1.0+1.0*0.95
    @test_approx_eq Distributions.params(b1.mixture.components[1])[2] bs2
    @test_approx_eq Distributions.params(b1.mixture.components[2])[2] bs2
    b1 = update(b1,2.0,0.5)
    @test mean(b1.mixture) == 2.0 && b1.μ == 2.0 && b1.σ == 0.5 && b1.m == 0.95
    @test_approx_eq Distributions.params(b1.mixture.components[1])[1] 2.0-0.5*0.95
    @test_approx_eq Distributions.params(b1.mixture.components[2])[1] 2.0+0.5*0.95
    @test_approx_eq Distributions.params(b1.mixture.components[1])[2] bs1
    @test_approx_eq Distributions.params(b1.mixture.components[2])[2] bs1
    @test_approx_eq Distributions.logpdf(b1,1.0) log(Distributions.pdf(b1.mixture,1.0))
    @test mean(b1) == mean(b1.mixture)
    @test var(b1) == var(b1.mixture)
    @test (srand(324) ; Distributions.rand(b1)) == (srand(324) ; Distributions.rand(b1.mixture))
    @test (srand(324) ; Distributions.rand!(b1,zeros(2))) == (srand(324) ; Distributions.rand!(b1.mixture,zeros(2)))
end

for p in [(:normal,Distributions.Normal{Float64}),
          (:lognormal,Distributions.LogNormal{Float64}),
          (:uniform,Distributions.Uniform{Float64}),
          (:laplace,Distributions.Laplace{Float64}),
          (:triangular,Distributions.SymTriangularDist{Float64})]
    d = distributions(p[1],[0.0,1.0],[1.0,2.0])
    @test eltype(d) == p[2]
end

#a test to ensure that all multivariate distributions copy their matrix arguments upon creation
#so that creating two from the same matrix does not introduce side-effects when recentering or rescaling
for p in [(:normal,ones(2),eye(2),2ones(2),mean,2ones(2),ones(2),cov,0.01*eye(2),eye(2)),
          (:lognormal,ones(2),eye(2),2ones(2),Distributions.location,log(2ones(2)),zeros(2),Distributions.scale,0.01*eye(2),eye(2)),
          (:normal,ones(Float32,2),eye(Float32,2),2ones(Float32,2),mean,2ones(Float32,2),ones(Float32,2),cov,0.01f0*eye(Float32,2),eye(Float32,2))]
    println("Testing multivariate ",p[1]," copying behaviour")
    d1 = distribution(p[1],p[2],p[3])
    d2 = distribution(p[1],p[2],p[3])
    d1 = recenter(d1,p[4])
    @test_approx_eq p[5](d1) p[6]
    @test_approx_eq p[5](d2) p[7]
    d1 = rescale(d1,0.1)
    @test_approx_eq p[8](d1) p[9]
    @test_approx_eq p[8](d2) p[10]
end


for p in [(:normal,eye(2),cov,0.01eye(2),eye(2))]
    println("Testing multivariate ",p[1]," (zero mean) copying behaviour")
    d1 = distribution(p[1],p[2])
    d2 = distribution(p[1],p[2])
    d1 = rescale(d1,0.1)
    @test_approx_eq p[3](d1) p[4]
    @test_approx_eq p[3](d2) p[5]
end

#test of vectorized distribution factory function
d = distributions(:bactrian,[0.0,1.0],[1.0,2.0],:normal,0.95)
@test eltype(d) == GeneralizedMetropolisHastings.Bactrian{Float64}

