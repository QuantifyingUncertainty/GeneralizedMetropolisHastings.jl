import GeneralizedMetropolisHastings.update_density!
import GeneralizedMetropolisHastings.propose!

m1 = zeros(2)
m2 = ones(2)
m3 = [4.0,2.0,0.0]
Σ1 = eye(2)
Σ2 = [1.0 0.0;0.0 4.0]
Σ3 = 4.0*eye(2)

println("======================")
println("Tests of NormalDensity")
println("======================")

###test the constructors
n1 = NormalDensity(Σ1)
n2 = NormalDensity(m2,Σ2)
n3 = NormalDensity(m1,Σ3)
n4 = NormalDensity(m3,eye(3))

@test length(n1.normal) == 2 && mean(n1.normal) == m1 && cov(n1.normal) == Σ1
@test length(n2.normal) == 2 && mean(n2.normal) == m2 && cov(n2.normal) == Σ2
@test length(n3.normal) == 2 && mean(n3.normal) == m1 && cov(n3.normal) == Σ3

###test the equality operators
@test NormalDensity(eye(3)) == NormalDensity(zeros(3),eye(3))
@test NormalDensity(eye(3)) != NormalDensity(ones(3),eye(3))
@test NormalDensity(eye(3)) != NormalDensity(eye(2))

###test the update functions
@test n1 != n2 != n3
@test update_density!(n1,m2) == update_density!(n2,Σ1) == update_density!(n3,m2,Σ1)
@test mean(update_density!(n1,m1).normal) == m1

###test the propose function
s1 = BaseSample(2)
s2 = BaseSample(2)
srand(0) ; v1 = m1 + randn(2) ; v2 = m1 + randn(2)
srand(0) ; propose!(n1,s1) ; propose!(n2,s2)
@test v1 == s1.values

###test the vectorized propose functions
vs1 = [BaseSample(2) for i=1:2]
@test vs1[1].values == vs1[2].values == zeros(2)
srand(0) ; propose!(n1,vs1)
@test vs1[1].values == v1 && vs1[2].values == v2
srand(0) ; propose!([n1,n2],vs1)
@test vs1[1].values == s1.values && vs1[2].values == s2.values

###testing the sampling behaviour
imax = 100000
vs2 = [BaseSample(3) for i=1:imax]
propose!(n4,vs2)
@test_approx_eq_eps mapreduce((x)->(x.values),+,vs2)/imax m3 0.01

###test the calculcation of proposal probabilities
n10 = NormalDensity([1.0,1.0],eye(2))
n11 = NormalDensity([2.0,2.0],eye(2))
s10 = BaseSample([1.0,1.0])
s11 = BaseSample([2.0,2.0])
@test logprobability(n10,s11) == logprobability(n11,s10) #symmetric distribution
@test logprobability(n10,s10) != logprobability(n10,s11)
@test logprobability([n10,n11],s11) == logprobability([n11,n10],s10)
@test_approx_eq_eps logprobability(n10,[s10,s11]) [-1.8379,-2.8379] 1e-4

###test the show functions
println()
println("====================")
println("Test show() function")
show(n1)
show(n4)
println("End  show() function")
println("====================")
println()

#################################################################################

println("=========================")
println("Tests of LogNormalDensity")
println("=========================")

###test the constructors
l1 = LogNormalDensity(Σ1)
l2 = LogNormalDensity(m2,Σ2)
l3 = LogNormalDensity(m1,Σ3)
l4 = LogNormalDensity(m3,eye(3))

@test length(l1.μ) == 2 && l1.μ == m1 && l1.Σ == Σ1
@test length(l2.μ) == 2 && l2.μ == m2 && l2.Σ == Σ2
@test length(l3.μ) == 2 && l3.μ == m1 && l3.Σ == Σ3

###test the equality operators
@test LogNormalDensity(eye(3)) == LogNormalDensity(zeros(3),eye(3))
@test LogNormalDensity(eye(3)) != LogNormalDensity(ones(3),eye(3))
@test LogNormalDensity(eye(3)) != LogNormalDensity(eye(2))

###test the update functions
@test l1 != l2 != l3
@test update_density!(l1,m2) == update_density!(l2,Σ1) == update_density!(l3,m2,Σ1)
@test update_density!(l1,m1).μ == m1

###test the propose function
s1 = BaseSample(2)
s2 = BaseSample(2)
srand(0) ; v1 = m1 + randn(2) ; v2 = m1 + randn(2)
srand(0) ; propose!(l1,s1) ; propose!(l2,s2)
@test_approx_eq v1 log(s1.values)

###test the vectorized propose functions
vs1 = [BaseSample(2) for i=1:2]
@test vs1[1].values == vs1[2].values == zeros(2)
srand(0) ; propose!(l1,vs1)
@test_approx_eq v1 log(vs1[1].values)
@test_approx_eq v2 log(vs1[2].values)
srand(0) ; propose!([l1,l2],vs1)
@test vs1[1].values == s1.values && vs1[2].values == s2.values

###testing the sampling behaviour
imax = 200000
vs2 = [BaseSample(3) for i=1:imax]
propose!(l4,vs2)
@test_approx_eq_eps log(mapreduce((x)->(x.values),+,vs2)/imax) m3+diag(eye(3))/2 0.01

###test the calculcation of proposal probabilities
l10 = LogNormalDensity([1.0,1.0],eye(2))
l11 = LogNormalDensity([2.0,2.0],eye(2))
s10 = BaseSample([1.0,1.0])
s11 = BaseSample([2.0,2.0])
@test logprobability(l10,s11) != logprobability(l11,s10) #asymmetric distribution
@test_approx_eq logprobability(LogNormalDensity([1.0],eye(1)),[BaseSample([1.0]),BaseSample([2.0])]) [-1.4189,-1.6592]

###test the show functions
println()
println("====================")
println("Test show() function")
show(l1)
show(l4)
println("End  show() function")
println("====================")
println()

nothing
