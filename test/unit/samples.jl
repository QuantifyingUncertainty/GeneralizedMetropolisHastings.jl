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
copy!(br,b1,3)
@test br.values == b1.values[:,3] && br.loglikelihood[1] == b1.loglikelihood[3] && br.logprior[1] == b1.logprior[3]
bo = samples(:base,2,1,Float64,Float64)
copy!(bo,br,1)
@test bo.values == br.values && bo.loglikelihood == br.loglikelihood && bo.logprior == br.logprior

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
copy!(gr,g1,2)
@test gr.values == g1.values[:,2] && gr.loglikelihood[1] == g1.loglikelihood[2] && gr.logprior[1] == g1.logprior[2] && gr.gradloglikelihood == g1.gradloglikelihood[:,2] && gr.gradlogprior == g1.gradlogprior[:,2]

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
copy!(tr,t1,2)
@test tr.values == t1.values[:,2] && tr.loglikelihood[1] == t1.loglikelihood[2] && tr.logprior[1] == t1.logprior[2]
@test tr.gradloglikelihood == t1.gradloglikelihood[:,2] && tr.gradlogprior == t1.gradlogprior[:,2] && tr.tensorloglikelihood == t1.tensorloglikelihood[:,:,2] && tr.tensorlogprior == t1.tensorlogprior[:,:,2]

a1 = samples(:tangent,2,3,4,Float64,Float64)
a2 = samples(:tangent,2,1,4,Float64,Float64)

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
ac = samples(:tangent,2,3,4,Float64,Float64)
ar = samples(:tangent,2,1,4,Float64,Float64)
@test a1 != ac
copy!(ac,a1)
@test a1 == ac
copy!(ar,a1,3)
@test ar.values == a1.values[:,3] && ar.loglikelihood[1] == a1.loglikelihood[3] && ar.logprior[1] == a1.logprior[3] && ar.tangentvectors == a1.tangentvectors[:,:,3]
@test ar.gradloglikelihood == a1.gradloglikelihood[:,3] && ar.gradlogprior == a1.gradlogprior[:,3] && ar.tensorloglikelihood == a1.tensorloglikelihood[:,:,3] && ar.tensorlogprior == a1.tensorlogprior[:,:,3]

@test numparas(b1) == 2 && numsamples(b1) == 3
@test numparas(b2) == 2 && numsamples(b2) == 1
@test numparas(g1) == 2 && numsamples(g1) == 3
@test numparas(t1) == 2 && numsamples(t1) == 3
@test numparas(a1) == 2 && numsamples(a1) == 3 && numtangents(a1) == 4

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
