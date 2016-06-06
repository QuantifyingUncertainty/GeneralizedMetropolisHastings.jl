@test GeneralizedMetropolisHastings._valuestuple(1,0) == (1,0)
@test GeneralizedMetropolisHastings._valuestuple(1,1) == (1,)
@test GeneralizedMetropolisHastings._valuestuple(1,2) == (1,2)
@test GeneralizedMetropolisHastings._valuestuple(2,0) == (2,0)
@test GeneralizedMetropolisHastings._valuestuple(2,1) == (2,)
@test GeneralizedMetropolisHastings._valuestuple(2,2) == (2,2)

@test ndims(zeros(GeneralizedMetropolisHastings._valuestuple(1,0))) == 2
@test ndims(zeros(GeneralizedMetropolisHastings._valuestuple(1,1))) == 1
@test ndims(zeros(GeneralizedMetropolisHastings._valuestuple(1,2))) == 2
@test ndims(zeros(GeneralizedMetropolisHastings._valuestuple(2,0))) == 2
@test ndims(zeros(GeneralizedMetropolisHastings._valuestuple(2,1))) == 1
@test ndims(zeros(GeneralizedMetropolisHastings._valuestuple(2,2))) == 2

@test GeneralizedMetropolisHastings._tensortuple(1,1,0) == (1,1,0)
@test GeneralizedMetropolisHastings._tensortuple(1,1,1) == (1,1)
@test GeneralizedMetropolisHastings._tensortuple(1,1,2) == (1,1,2)
@test GeneralizedMetropolisHastings._tensortuple(2,1,0) == (2,1,0)
@test GeneralizedMetropolisHastings._tensortuple(2,1,1) == (2,1)
@test GeneralizedMetropolisHastings._tensortuple(2,1,2) == (2,1,2)

@test ndims(zeros(GeneralizedMetropolisHastings._tensortuple(1,1,0))) == 3
@test ndims(zeros(GeneralizedMetropolisHastings._tensortuple(1,1,1))) == 2
@test ndims(zeros(GeneralizedMetropolisHastings._tensortuple(1,1,2))) == 3
@test ndims(zeros(GeneralizedMetropolisHastings._tensortuple(1,2,0))) == 3
@test ndims(zeros(GeneralizedMetropolisHastings._tensortuple(1,2,1))) == 2
@test ndims(zeros(GeneralizedMetropolisHastings._tensortuple(1,2,2))) == 3

v1 = [1.0 2.0 3.0;4.0 5.0 6.0]
ll1 = [1.0,2.0,3.0]
lp1 = [0.1,0.2,0.3]
gll1 = 2*v1
glp1 = 3*v1
tll1 = similar(eye(2),2,2,3) ; tll1[:,:,1] = eye(2) ; tll1[:,:,2] = eye(2) ; tll1[:,:,3] = eye(2) ; tll1
tlp1 = similar(eye(2),2,2,3) ; tlp1[:,:,1] = ones(2,2) ; tlp1[:,:,2] = ones(2,2) ; tlp1[:,:,3] = ones(2,2) ; tlp1
tt1 = similar(eye(2,3),2,4,3) ; tt1[:,:,1] = ones(2,4) ; tt1[:,:,2] = 2*ones(2,4) ; tt1[:,:,3] = 3*ones(2,4) ; tt1

b1 = samples(:base,2,3,Float64,Float64)
b2 = samples(:base,2,1,Int,Float64)
b3 = samples(:base,2,3,Float32,Float32)
b4 = samples(:base,2,1,Int16,Float32)

@test typeof(b1) <: BaseSample && eltype(b1.values) == Float64 && size(b1.values) == (2,3) && size(b1.loglikelihood) == (3,) && size(b1.logprior) == (3,)
@test typeof(b2) <: BaseSample && eltype(b2.values) == Int && size(b2.values) == (2,) && size(b2.loglikelihood) == (1,) && size(b2.logprior) == (1,)
@test typeof(b3) <: BaseSample && eltype(b3.values) == Float32 && size(b3.values) == (2,3) && size(b3.loglikelihood) == (3,) && size(b3.logprior) == (3,)
@test typeof(b4) <: BaseSample && eltype(b4.values) == Int16 && size(b4.values) == (2,) && size(b4.loglikelihood) == (1,) && size(b4.logprior) == (1,)
@test ndims(b1.values) == 2 && ndims(b2.values) == 1

copy!(b1.values,v1)
copy!(b1.loglikelihood,ll1)
copy!(b1.logprior,lp1)
bc = samples(:base,2,3,Float64,Float64)
br = samples(:base,2,1,Float64,Float64)
@test b1 != bc
copy!(bc,b1)
@test b1 == bc
copy!(br,1,b1,3)
@test br.values == b1.values[:,3] && br.loglikelihood[1] == b1.loglikelihood[3] && br.logprior[1] == b1.logprior[3]
bo = samples(:base,2,1,Float64,Float64)
copy!(bo,1,br,1)
@test bo.values == br.values && bo.loglikelihood == br.loglikelihood && bo.logprior == br.logprior

bs1 = similar(b1)
bs2 = similar(b1,2)

@test isa(bs1,BaseSample) && numparas(bs1) == numparas(b1) && numsamples(bs1) == numsamples(b1) && sampletype(bs1) == sampletype(b1) && calculationtype(bs1) == calculationtype(b1)
@test isa(bs2,BaseSample) && numparas(bs2) == numparas(b1) && numsamples(bs2) == 2 && sampletype(bs2) == sampletype(b1) && calculationtype(bs2) == calculationtype(b1)

bc1 = copy(b1)
bc2 = copy(b1,1)
bc3 = copy(b1,[3,2])
bc4 = similar(b1) ; copy!(bc4,1,bc2,1) ; copy!(bc4,[2,3],bc3,[2,1])

@test bc1 == b1
@test bc4 == b1
@test numsamples(copy(b1,[])) == 0

offset!(bc1,[1.0,2.0])
@test bc1.values == (b1.values + [1.0 1.0 1.0;2.0 2.0 2.0])

###############################################################################################################################################################
g1 = samples(:gradient,2,3,Float64,Float64)
g2 = samples(:gradient,2,1,Float32,Float32)

@test typeof(g1) <: GradientSample && size(g1.values) == (2,3) && size(g1.loglikelihood) == (3,) && size(g1.logprior) == (3,) && size(g1.gradloglikelihood) == (2,3) && size(g1.gradlogprior) == (2,3)
@test typeof(g2) <: GradientSample && size(g2.values) == (2,) && size(g2.loglikelihood) == (1,) && size(g2.logprior) == (1,) && size(g2.gradloglikelihood) == (2,) && size(g2.gradlogprior) == (2,)
@test ndims(g1.values) == 2 && ndims(g2.values) == 1

copy!(g1.values,v1)
copy!(g1.loglikelihood,ll1)
copy!(g1.logprior,lp1)
copy!(g1.gradloglikelihood,gll1)
copy!(g1.gradlogprior,glp1)
gc = samples(:gradient,2,3,Float64,Float64)
gr = samples(:gradient,2,1,Float64,Float64)
@test g1 != gc
copy!(gc,g1)
@test g1 == gc
copy!(gr,1,g1,2)
@test gr.values == g1.values[:,2] && gr.loglikelihood[1] == g1.loglikelihood[2] && gr.logprior[1] == g1.logprior[2] && gr.gradloglikelihood == g1.gradloglikelihood[:,2] && gr.gradlogprior == g1.gradlogprior[:,2]

gs1 = similar(g1)
gs2 = similar(g1,2)
gs3 = similar(g1,2,:base)

@test isa(gs1,GradientSample) && numparas(gs1) == numparas(g1) && numsamples(gs1) == numsamples(g1) && sampletype(gs1) == sampletype(g1) && calculationtype(gs1) == calculationtype(g1)
@test isa(gs2,GradientSample) && numparas(gs2) == numparas(g1) && numsamples(gs2) == 2 && sampletype(gs2) == sampletype(g1) && calculationtype(gs2) == calculationtype(g1)
@test isa(gs3,BaseSample) && numparas(gs3) == numparas(g1) && numsamples(gs3) == 2 && sampletype(gs3) == sampletype(g1) && calculationtype(gs3) == calculationtype(g1)

gc1 = copy(g1)
gc2 = copy(g1,1)
gc3 = copy(g1,[3,2])
gc4 = similar(g1) ; copy!(gc4,1,gc2,1) ; copy!(gc4,[2,3],gc3,[2,1])
gc5 = copy(g1,[3,2],:base)

@test gc1 == g1
@test gc4 == g1
@test isa(gc5,BaseSample) && g1.values[:,[3,2]] == gc5.values
@test numsamples(copy(g1,[])) == 0

#################################################################################################################################################################

t1 = samples(:tensor,2,3,Float32,Float32)
t2 = samples(:tensor,2,1,Float64,Float64)

@test typeof(t1) <: TensorSample && size(t1.values) == (2,3) && size(t1.loglikelihood) == (3,) && size(t1.logprior) == (3,) && size(t1.gradloglikelihood) == (2,3) && size(t1.gradlogprior) == (2,3) && size(t1.tensorloglikelihood) == (2,2,3) && size(t1.tensorlogprior) == (2,2,3)
@test typeof(t2) <: TensorSample && size(t2.values) == (2,) && size(t2.loglikelihood) == (1,) && size(t2.logprior) == (1,) && size(t2.gradloglikelihood) == (2,) && size(t2.gradlogprior) == (2,) && size(t2.tensorloglikelihood) == (2,2) && size(t2.tensorlogprior) == (2,2)
@test ndims(t1.values) == 2 && ndims(t2.values) == 1

copy!(t1.values,v1)
copy!(t1.loglikelihood,ll1)
copy!(t1.logprior,lp1)
copy!(t1.gradloglikelihood,gll1)
copy!(t1.gradlogprior,glp1)
copy!(t1.tensorloglikelihood,tll1)
copy!(t1.tensorlogprior,tlp1)
tc = samples(:tensor,2,3,Float32,Float32)
tr = samples(:tensor,2,1,Float32,Float32)
@test t1 != tc
copy!(tc,t1)
@test t1 == tc
copy!(tr,1,t1,2)
@test tr.values == t1.values[:,2] && tr.loglikelihood[1] == t1.loglikelihood[2] && tr.logprior[1] == t1.logprior[2]
@test tr.gradloglikelihood == t1.gradloglikelihood[:,2] && tr.gradlogprior == t1.gradlogprior[:,2] && tr.tensorloglikelihood == t1.tensorloglikelihood[:,:,2] && tr.tensorlogprior == t1.tensorlogprior[:,:,2]

ts1 = similar(t1)
ts2 = similar(t1,2)
ts3 = similar(t1,2,:base)

@test isa(ts1,TensorSample) && numparas(ts1) == numparas(t1) && numsamples(ts1) == numsamples(t1) && sampletype(ts1) == sampletype(t1) && calculationtype(ts1) == calculationtype(t1)
@test isa(ts2,TensorSample) && numparas(ts2) == numparas(t1) && numsamples(ts2) == 2 && sampletype(ts2) == sampletype(t1) && calculationtype(ts2) == calculationtype(t1)
@test isa(ts3,BaseSample) && numparas(ts3) == numparas(t1) && numsamples(ts3) == 2 && sampletype(ts3) == sampletype(t1) && calculationtype(ts3) == calculationtype(t1)

@test numsamples(copy(t1,[])) == 0

#####################################################################################################################################################################

a1 = samples(:tangent,2,3,Float64,Float64,4)
a2 = samples(:tangent,2,1,Float64,Float64,4)

@test typeof(a1) <: TangentTensorSample && size(a1.values) == (2,3) && size(a1.loglikelihood) == (3,) && size(a1.logprior) == (3,) && size(a1.gradloglikelihood) == (2,3) && size(a1.gradlogprior) == (2,3) && size(a1.tensorloglikelihood) == (2,2,3) && size(a1.tensorlogprior) == (2,2,3) && size(a1.tangentvectors) == (2,4,3)
@test typeof(a2) <: TangentTensorSample && size(a2.values) == (2,) && size(a2.loglikelihood) == (1,) && size(a2.logprior) == (1,) && size(a2.gradloglikelihood) == (2,) && size(a2.gradlogprior) == (2,) && size(a2.tensorloglikelihood) == (2,2) && size(a2.tensorlogprior) == (2,2) && size(a2.tangentvectors) == (2,4)
@test ndims(a1.values) == 2 && ndims(a2.values) == 1

copy!(a1.values,v1)
copy!(a1.loglikelihood,ll1)
copy!(a1.logprior,lp1)
copy!(a1.gradloglikelihood,gll1)
copy!(a1.gradlogprior,glp1)
copy!(a1.tensorloglikelihood,tll1)
copy!(a1.tensorlogprior,tlp1)
copy!(a1.tangentvectors,tt1)
ac = samples(:tangent,2,3,Float64,Float64,4)
ar = samples(:tangent,2,1,Float64,Float64,4)
@test a1 != ac
copy!(ac,a1)
@test a1 == ac
copy!(ar,1,a1,3)
@test ar.values == a1.values[:,3] && ar.loglikelihood[1] == a1.loglikelihood[3] && ar.logprior[1] == a1.logprior[3] && ar.tangentvectors == a1.tangentvectors[:,:,3]
@test ar.gradloglikelihood == a1.gradloglikelihood[:,3] && ar.gradlogprior == a1.gradlogprior[:,3] && ar.tensorloglikelihood == a1.tensorloglikelihood[:,:,3] && ar.tensorlogprior == a1.tensorlogprior[:,:,3]

as1 = similar(a1)
as2 = similar(a1,2)
as3 = similar(a1,2,:base)

@test isa(as1,TangentTensorSample) && numparas(as1) == numparas(a1) && numsamples(as1) == numsamples(a1) && sampletype(as1) == sampletype(a1) && calculationtype(as1) == calculationtype(a1)
@test isa(as2,TangentTensorSample) && numparas(as2) == numparas(a1) && numsamples(as2) == 2 && sampletype(as2) == sampletype(a1) && calculationtype(as2) == calculationtype(a1)
@test isa(as3,BaseSample) && numparas(as3) == numparas(a1) && numsamples(as3) == 2 && sampletype(as3) == sampletype(a1) && calculationtype(as3) == calculationtype(a1)

@test numparas(b1) == 2 && numsamples(b1) == 3
@test numparas(b2) == 2 && numsamples(b2) == 1
@test numparas(g1) == 2 && numsamples(g1) == 3
@test numparas(t1) == 2 && numsamples(t1) == 3
@test numparas(a1) == 2 && numsamples(a1) == 3 && numtangents(a1) == 4

@test numsamples(copy(a1,[])) == 0

println()
println("====================")
println("Test show() function")
show(b1)
show(g1)
show(t1)
show(a1)
show(b2)
show(g2)
show(t2)
show(a2)
println("End  show() function")
println("====================")
println()

nothing
