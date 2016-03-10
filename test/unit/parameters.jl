###test local default functions
@test GeneralizedMetropolisHastings._defaultkey() == :param
@test GeneralizedMetropolisHastings._defaultvalue(Float64) == 0.0
@test GeneralizedMetropolisHastings._defaultvalue(Float32) == 0.0f0
@test GeneralizedMetropolisHastings._defaultvalue(Distributions.Uniform(0.0,11.0)) == 5.5
@test GeneralizedMetropolisHastings._defaultvalue(Distributions.DiscreteUniform(0,11)) == 5
@test GeneralizedMetropolisHastings._defaultunivariate(0.0,1.0) == Distributions.Uniform(0.0,1.0)
@test GeneralizedMetropolisHastings._defaultunivariate(0,1) == Distributions.DiscreteUniform(0,1)
@test_throws MethodError GeneralizedMetropolisHastings._defaultunivariate(0.0,1)

###test the parameter factory functions for ParameterDefault
d1 = parameter(:d,3.0)
d2 = parameter(:d,Float64)
d3 = parameter(:d,Float32)
d4 = parameter(3.0f0)
d5 = parameter(Float64)
d6 = parameter(Float32)
@test typeof(d1) <: ParameterDefault{Float64} && d1.key == :d && d1.default == 3.0
@test typeof(d2) <: ParameterDefault{Float64} && d2.key == :d && d2.default == 0.0
@test typeof(d3) <: ParameterDefault{Float32} && d3.key == :d && d3.default == 0.0f0
@test typeof(d4) <: ParameterDefault{Float32} && d4.default == 3.0f0
@test typeof(d5) <: ParameterDefault{Float64} && d5.default == 0.0
@test typeof(d6) <: ParameterDefault{Float32} && d6.default == 0.0
@test d2 == parameter(:d,Float64) #equality operator
@test d6 == parameter(Float32)
@test d1 != d2
@test d1 != d5

###test parameter factory functions for ParameterUnivariate
p1 = parameter(:p,Distributions.Normal(),3.0)
p2 = parameter(:p,Distributions.Binomial(),0)
p3 = parameter(:p,Distributions.Uniform(0.0,11.0))
p4 = parameter(:p,Distributions.DiscreteUniform(0,9))
p5 = parameter(:p,0.0,11.0,3.0)
p6 = parameter(:p,0.0,11.0)
p7 = parameter(Distributions.Uniform(0.0,11.0))
p8 = parameter(Distributions.DiscreteUniform(0,9))
p9 = parameter(0.0,11.0,3.0)
p10 = parameter(0.0,11.0)
@test typeof(p1) <: ParameterUnivariate && p1.key == :p && p1.default == 3.0 && p1.prior == Distributions.Normal()
@test typeof(p2) <: ParameterUnivariate && p2.key == :p && p2.default == 0 && p2.prior == Distributions.Binomial()
@test typeof(p3) <: ParameterUnivariate && p3.key == :p && p3.default == 5.5 && p3.prior == Distributions.Uniform(0.0,11.0)
@test typeof(p4) <: ParameterUnivariate && p4.key == :p && p4.default == 4 && p4.prior == Distributions.DiscreteUniform(0,9)
@test typeof(p5) <: ParameterUnivariate && p5.key == :p && p5.default == 3.0 && p5.prior == Distributions.Uniform(0.0,11.0)
@test typeof(p6) <: ParameterUnivariate && p6.key == :p && p6.default == 5.5 && p6.prior == Distributions.Uniform(0.0,11.0)
@test typeof(p7) <: ParameterUnivariate && p7.default == 5.5 && p7.prior == Distributions.Uniform(0.0,11.0)
@test typeof(p8) <: ParameterUnivariate && p8.default == 4 && p8.prior == Distributions.DiscreteUniform(0,9)
@test typeof(p9) <: ParameterUnivariate && p9.default == 3.0 && p9.prior == Distributions.Uniform(0.0,11.0)
@test typeof(p10) <: ParameterUnivariate && p10.default == 5.5 && p10.prior == Distributions.Uniform(0.0,11.0)
@test p3 == p6
@test p7 == p10
@test p1 != p2
@test p1 != p7
@test_throws MethodError parameter(:u,Distributions.Normal(),3)
@test_throws MethodError parameter(:u,Distributions.Binomial(),1.0)
@test_throws MethodError parameter(:u,1.0,2.0,3)
@test_throws AssertionError parameter(:u,Distributions.Uniform(1.0,2.0),0.0)
@test_throws MethodError parameter(Distributions.Normal(),3)
@test_throws MethodError parameter(Distributions.Binomial(),1.0)
@test_throws MethodError parameter(1.0,3)
@test_throws AssertionError parameter(Distributions.Uniform(1.0,2.0),0.0)

@test length(parameters([:d1,:d2],[3.0,1.0])) == 2
@test length(parameters([:d1,:d2],Float32)) == 2
@test length(parameters([1.0,2.0,3.0])) == 3
@test length(parameters(2,Float32)) == 2
vec1 = parameters([:d1,:d2],[Distributions.Normal(),Distributions.Uniform(0.0,11.0)],[0.0,2.0])
vec2 = parameters([:d1,:d2],[0.0,5.0],[6.0,10.0],[3.0,7.5])
vec3 = parameters([:d1,:d2],[Distributions.Normal(),Distributions.Uniform(0.0,11.0)])
vec4 = parameters([:d1,:d2],[0.0,5.0],[6.0,10.0])
vec5 = parameters([Distributions.Normal(),Distributions.Uniform(0.0,11.0)],[0.0,2.0])
vec6 = parameters([0.0,5.0],[6.0,10.0],[3.0,7.5])
vec7 = parameters([Distributions.Normal(),Distributions.Uniform(0.0,11.0)])
vec8 = parameters([0.0,5.0],[6.0,10.0])
vec9 = parameters([3.0,7.5])
@test vec2 == vec4
@test vec6 == vec8

###test the initvalues functions
v1 = GeneralizedMetropolisHastings._initvalue(Val{:default},d1,Float64)
v3 = GeneralizedMetropolisHastings._initvalue(Val{:default},d3,Float32)
v4 = GeneralizedMetropolisHastings._initvalue(Val{:default},d4,Float64)
v6 = GeneralizedMetropolisHastings._initvalue(Val{:default},d6,Float32)
@test typeof(v1) == Float64 && v1 == 3.0
@test typeof(v3) == Float32 && v3 == 0.0f0
@test typeof(v4) == Float64 && v4 == 3.0
@test typeof(v6) == Float32 && v6 == 0.0f0
@test GeneralizedMetropolisHastings._initvalue(Val{:default},p1,Float64) == 3.0
@test GeneralizedMetropolisHastings._initvalue(Val{:default},p3,Float64) == 5.5
@test GeneralizedMetropolisHastings._initvalue(Val{:default},p5,Float64) == 3.0
@test GeneralizedMetropolisHastings._initvalue(Val{:default},p8,Float64) == 4.0
@test (srand(0) ; GeneralizedMetropolisHastings._initvalue(Val{:prior},p3,Float64)) == (srand(0) ; rand(Distributions.Uniform(0.0,11.0)))
@test (srand(0) ; GeneralizedMetropolisHastings._initvalue(Val{:prior},p4,Float32)) == (srand(0) ; Float32(rand(Distributions.DiscreteUniform(0,9))))
@test (srand(0) ; GeneralizedMetropolisHastings._initvalue(Val{:prior},p4,Int)) == (srand(0) ; Int(rand(Distributions.DiscreteUniform(0,9))))
@test GeneralizedMetropolisHastings._initvalues!(Val{:default},AbstractParameter[d1,d3,p1,p4],Array{Float32}(4)) == [3.0f0,0.0f0,3.0f0,4.0f0]
@test initvalues!(trait(:initialize,:default),AbstractParameter[d1,d3,p1,p4],Array{Float64}(4)) == [3.0,0.0,3.0,4.0]
@test initvalues(trait(:initialize,:default),AbstractParameter[d1,d3,p1,p4],Float64) == [3.0,0.0,3.0,4.0]
@test (srand(0) ; initvalues(trait(:initialize,:prior),AbstractParameter[d1,d3,p1,p4],Float64)) == (srand(0) ; [3.0,0.0,rand(Distributions.Normal()),Float64(rand(Distributions.DiscreteUniform(0,9)))])
@test_throws AssertionError initvalues!(trait(:initialize,:default),AbstractParameter[p1],Array{Float64}(2))

###test the logprior functions
l1 = parameter(Distributions.Uniform(1.0,3.0))
l2 = parameter(2.0)
l3 = parameter(Distributions.DiscreteUniform(1.0,3.0))
@test GeneralizedMetropolisHastings._logprior(l1,2.0,Float64) == log(1.0/2.0)
@test GeneralizedMetropolisHastings._logprior(l1,4.0,Float64) == -Inf
@test GeneralizedMetropolisHastings._logprior(l2,2.0,Float64) == 0.0
@test GeneralizedMetropolisHastings._logprior(l1,2.0,Float32) == -log(2.0f0)
@test GeneralizedMetropolisHastings._logprior(l1,4.0,Float32) == -Inf32
@test GeneralizedMetropolisHastings._logprior(l2,2.0,Float32) == 0.0f0
@test GeneralizedMetropolisHastings._logprior(l3,2.0,Float64) == log(1.0/3.0)
@test GeneralizedMetropolisHastings._logprior([l1,l2,l3],[2.0,2.0,2.0],Float64) == log(1.0/2.0) + log(1.0/3.0)
@test logprior(l1,2.0,Float64) == log(1.0/2.0)
@test logprior(AbstractParameter[l1,l2,l3],[2.0,2.0,2.0],Float64) == log(1.0/2.0) + log(1.0/3.0)
@test logprior(AbstractParameter[l1,l1],[2.0,4.0],Float64) == -Inf
@test logprior(AbstractParameter[],Float64[],Float64) == 0.0
@test logprior(ParameterDefault[l2,l2],[2.0,4.0],Float64) == 0.0
@test logprior(ParameterUnivariate[l1,l3],[2.0,2.0],Float64) == log(1.0/2.0) + log(1.0/3.0)
@test logprior!(zeros(2),[l1,l2,l3],[2.0 2.0;2.0 2.0;2.0 4.0]) == [log(1.0/2.0)+log(1.0/3.0),-Inf]
@test logprior([l1,l2,l3],[2.0 2.0;2.0 2.0;2.0 4.0],Float64) == [log(1.0/2.0)+log(1.0/3.0),-Inf]
@test logprior(ParameterDefault[l2,l2],[2.0 2.0;2.0 2.0],Float64) == zeros(2)

println("====================")
println("Test show() function")
println()
show(d1)
show(d4)
show(p1)
show(p3)
show(p4)
show(AbstractParameter[l1,l2,l3])
show(ParameterUnivariate[l1,l3])
show(vec1)
show(vec4)
show([vec1;vec4])
show([vec1;vec9])
println()
println("End  show() function")
println("====================")
println()

nothing
