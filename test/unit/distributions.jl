import GeneralizedMetropolisHastings.Normal
import GeneralizedMetropolisHastings.LogNormal

println("==================")
println("Tests of LogNormal")
println("==================")

###test the constructors
l1 = LogNormal(zeros(1),eye(1))
l2 = LogNormal(eye(1))
l3 = LogNormal(1,1.0)
l4 = LogNormal(ones(1),4.0*eye(1))
l5 = LogNormal([1.0,2.0],[0.04 0.01;0.01 0.01])
l6 = LogNormal([1.0,1.0],[0.0225 -0.005;-0.005 0.01])

@test mean(l1.normal) == mean(l2.normal) == mean(l3.normal) == zeros(1)
@test cov(l1.normal) == cov(l2.normal) == cov(l3.normal) == eye(1)
@test mean(l4.normal) == ones(1) && cov(l4.normal) == 4.0*eye(1)
@test mean(l5.normal) == [1.0,2.0] && cov(l5.normal) == [0.04 0.01;0.01 0.01]
@test mean(l6.normal) == [1.0,1.0] && cov(l6.normal) == [0.0225 -0.005;-0.005 0.01]
@test length(l1) == length(l2) == length(l3) == length(4) == 1
@test length(l5) == 2

###test the support of the distribution (positive numbers)
@test Distributions.insupport(l1,[-1.0]) == false
@test Distributions.insupport(l1,[0.0]) == false
@test Distributions.insupport(l1,[0.1]) == true
@test Distributions.insupport(l1,[Inf]) == false
@test Distributions.insupport(l5,[1.0,2.0]) == true
@test Distributions.insupport(l5,[0.0,1.0]) == false

###test the random generator via the _rand()! function
srand(0) ; x1 = Distributions.rand(GeneralizedMetropolisHastings.Normal(zeros(1),eye(1)))
srand(0) ; x2 = zeros(1) ; Distributions._rand!(l1,x2)
@test_approx_eq exp(x1) x2[1]

###test the batch generator rand() function
srand(0) ; x1 = Distributions.rand(GeneralizedMetropolisHastings.Normal([1.0,2.0],[0.04 0.01;0.01 0.01]))
srand(0) ; x2 = zeros(2) ; Distributions.rand!(l5,x2)
@test_approx_eq exp(x1) x2

###test the probability density function
@test Distributions._pdf(l1,[-1.0]) == zeros(1)
@test Distributions._pdf(l1,[0.0]) == zeros(1)
@test Distributions._pdf(l1,[Inf]) == zeros(1)
@test_approx_eq_eps Distributions._pdf(l1,[1.0]) [0.3989] 1e-4
@test_approx_eq_eps Distributions._pdf(l4,[1.0]) [0.1760] 1e-4
@test_approx_eq_eps Distributions._pdf(l5,[exp(1.0),exp(2.0)]) [0.4575] 1e-4
@test_approx_eq_eps Distributions._pdf(l6,[exp(1.0),exp(1.0)]) [1.5231] 1e-4

###test the log probability density function
@test Distributions._logpdf(l1,[-1.0]) == [-Inf]
@test Distributions._logpdf(l1,[0.0]) == [-Inf]
@test Distributions._logpdf(l1,[Inf]) == [-Inf]
@test_approx_eq_eps Distributions._logpdf(l1,[1.0]) [-0.9189] 1e-4
@test_approx_eq_eps Distributions._logpdf(l4,[1.0]) [-1.7371] 1e-4
@test_approx_eq_eps Distributions._logpdf(l5,[exp(1.0),exp(2.0)]) [-0.7820] 1e-4
@test_approx_eq_eps Distributions._logpdf(l6,[exp(1.0),exp(1.0)]) [0.4207] 1e-4

@test_approx_eq mean(l4) [exp(3.0)]
@test_approx_eq cov(l4) [exp(6.0)*(exp(4.0)-1.0)]
@test_approx_eq var(l4) [exp(6.0)*(exp(4.0)-1.0)]

###test the show functions
println()
println("====================")
println("Test show() function")
show(l1)
show(l5)
println("End  show() function")
println("====================")
println()

