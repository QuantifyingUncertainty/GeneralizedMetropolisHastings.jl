b1 = GeneralizedMetropolisHastings.BaseSample([1.0,2.0])
b2 = GeneralizedMetropolisHastings.BaseSample(zeros(Float32,2,3))
b3 = samples(:base,3,2,Float64)
b4 = samples(:base,2,1,Float32)

@test typeof(b1) <: GeneralizedMetropolisHastings.BaseSample && eltype(b1.values) == Float64 && isequal(b1.values,[1.0,2.0]) && size(b1.loglikelihood) == (1,) && size(b1.logprior) == (1,)
@test typeof(b2) <: GeneralizedMetropolisHastings.BaseSample && eltype(b2.values) == Float32 && isequal(b2.values,zeros(Float32,2,3)) && size(b2.loglikelihood) == (3,) && size(b2.logprior) == (3,)
@test typeof(b3) <: GeneralizedMetropolisHastings.BaseSample && eltype(b3.values) == Float64 && size(b3.values) == (3,2) && size(b3.loglikelihood) == (2,) && size(b3.logprior) == (2,)
@test typeof(b4) <: GeneralizedMetropolisHastings.BaseSample && eltype(b4.values) == Float32 && size(b4.values) == (2,1) && size(b4.loglikelihood) == (1,) && size(b4.logprior) == (1,)

g1 = GeneralizedMetropolisHastings.GradientSample([1.0,2.0,3.0])
g2 = GeneralizedMetropolisHastings.GradientSample([1.0,2.0;3.0,2.0])
g3 = samples(:gradient,2,4)

@test typeof(g1) <: GeneralizedMetropolisHastings.GradientSample && isequal(g1.values,[1.0,2.0,3.0]) && size(g1.loglikelihood) == (1,) && size(g1.logprior) == (1,) && size(g1.gradloglikelihood) == (3,1) && size(g1.gradlogprior) == (3,1)
@test typeof(g2) <: GeneralizedMetropolisHastings.GradientSample && isequal(g2.values,[1.0,2.0;3.0,2.0]) && size(g2.loglikelihood) == (2,) && size(g2.logprior) == (2,) && size(g2.gradloglikelihood,1) == (2,2) && size(g2.gradlogprior,1) == (2,2)
@test typeof(g3) <: GeneralizedMetropolisHastings.GradientSample && size(g3.values) == (2,4) && size(g3.loglikelihood) == (4,) && size(g3.logprior) == (4,) && size(g3.gradloglikelihood) == (2,4) && size(g3.gradlogprior) == (2,4)

t1 = GeneralizedMetropolisHastings.TensorSample([1.0,3.0])
t2 = GeneralizedMetropolisHastings.TensorSample(zeros(5,2))
t3 = samples(:tensor,2.3)

@test typeof(t1) <: GeneralizedMetropolisHastings.TensorSample && isequal(t1.values,[1.0,3.0]) && size(t1.loglikelihood) == (1,) && size(t1.logprior) == (1,) && size(t1.gradloglikelihood) == (2,) && size(t1.gradlogprior) == (2,) && size(t1.tensorloglikelihood) == (2,2,1) && size(t1.tensorlogprior) == (2,2,1)
@test typeof(t2) <: GeneralizedMetropolisHastings.TensorSample && isequal(t2.values,zeros(5,2)) && size(t2.loglikelihood) == (2,) && size(t2.logprior) == (2,) && size(t2.gradloglikelihood) == (5,2) && size(t2.gradlogprior) == (5,2) && size(t2.tensorloglikelihood) == (5,5,2) && size(t2.tensorlogprior) == (5,5,2)
@test typeof(t3) <: GeneralizedMetropolisHastings.TensorSample && size(t3.values) == (2,3) && size(t3.loglikelihood) == (3,) && size(t3.logprior) == (3,) && size(t3.gradloglikelihood) == (2,3) && size(t3.gradlogprior) == (2,3) && size(t3.tensorloglikelihood) == (2,2,3) && size(t3.tensorlogprior) == (2,2,3)

a1 = GeneralizedMetropolisHastings.TangentTensorSample([1.0,3.0],2)
a2 = GeneralizedMetropolisHastings.TangentTensorSample(zeros(5,2),3)
a3 = samples(:tangent,2,3,4)

@test typeof(a1) <: GeneralizedMetropolisHastings.TangentTensorSample && isequal(a1.values,[1.0,3.0]) && size(a1.gradloglikelihood) == (1,) && size(a1.gradlogprior) == (1,) && size(a1.tensorloglikelihood) == (2,2,1) && size(a1.tensorlogprior) == (2,2,1) && size(a1.tangentvectors) == (2,2,1)
@test typeof(a2) <: GeneralizedMetropolisHastings.TangentTensorSample && isequal(a2.values,zeros(5,2)) && size(a2.gradloglikelihood) == (2,) && size(a2.gradlogprior) == (2,) && size(a2.tensorloglikelihood) == (5,5,2) && size(a2.tensorlogprior) == (5,5,2) && size(a2.tangentvectors) == (5,3,2)
@test typeof(a3) <: GeneralizedMetropolisHastings.TangentTensorSample && size(a3.values) == (2,3) && size(a3.gradloglikelihood) == (3,) && size(a3.gradlogprior) == (3,) && size(a3.tensorloglikelihood) == (2,2,3) && size(a3.tensorlogprior) == (2,2,3) && size(a3.tangentvectors) == (2,4,3)

@test numparas(b1) == 2 && numsamples(b1) == 1
@test numparas(b2) == 2 && numsamples(b2) == 3
@test numparas(g3) == 2 && numsamples(g3) == 4
@test numparas(t3) == 2 && numsamples(t3) == 3
@test numparas(a3) == 2 && numsamples(a3) == 3 && numtangents(a3) == 4

println()
println("====================")
println("Test show() function")
show(b1)
show(g1)
show(t1)
show(a1)
println("End  show() function")
println("====================")
println()

nothing
