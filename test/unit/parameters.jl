###test local type aliases
@test GeneralizedMetropolisHastings.ContinuousType == Float64
@test GeneralizedMetropolisHastings.DiscreteType == Int
@test GeneralizedMetropolisHastings.ContinuousType <: GeneralizedMetropolisHastings.ValueType
@test GeneralizedMetropolisHastings.DiscreteType <: GeneralizedMetropolisHastings.ValueType

###test local default functions
@test GeneralizedMetropolisHastings._defaultparametervalue(Float64) == 0.0
@test GeneralizedMetropolisHastings._defaultparametervalue(Float32) == 0.0f0
@test_throws MethodError GeneralizedMetropolisHastings._defaultparametervalue(Int)
@test GeneralizedMetropolisHastings._defaultparametervalue(Distributions.Uniform(0.0,11.0),GeneralizedMetropolisHastings.ContinuousType) == 5.5
@test GeneralizedMetropolisHastings._defaultparametervalue(Distributions.DiscreteUniform(0,11),GeneralizedMetropolisHastings.DiscreteType) == 5
@test GeneralizedMetropolisHastings._defaultparameterprior(0.0,1.0) == Distributions.Uniform(0.0,1.0)
@test GeneralizedMetropolisHastings._defaultparameterprior(0,1) == Distributions.DiscreteUniform(0,1)
@test_throws MethodError GeneralizedMetropolisHastings._defaultparameterprior(0.0,1)

###test the parameter factory functions for ParameterDefault
d1 = parameter("d",3.0)
d2 = parameter("d",Float64)
d3 = parameter("d",Float32)
d4 = parameter("d")
d5 = parameter(3.0f0)
d6 = parameter(Float64)
d7 = parameter(Float32)
d8 = parameter()
@test typeof(d1) <: GeneralizedMetropolisHastings.NamedParameterDefault{ASCIIString,Float64} && d1.name == "d" && d1.default == 3.0
@test typeof(d2) <: GeneralizedMetropolisHastings.NamedParameterDefault{ASCIIString,Float64} && d2.name == "d" && d2.default == 0.0
@test typeof(d3) <: GeneralizedMetropolisHastings.NamedParameterDefault{ASCIIString,Float32} && d3.name == "d" && d3.default == 0.0f0
@test typeof(d4) <: GeneralizedMetropolisHastings.NamedParameterDefault{ASCIIString,Float64} && d4.name == "d" && d4.default == 0.0
@test typeof(d5) <: GeneralizedMetropolisHastings.UnnamedParameterDefault{Float32} && d5.default == 3.0f0
@test typeof(d6) <: GeneralizedMetropolisHastings.UnnamedParameterDefault{Float64} && d6.default == 0.0
@test typeof(d7) <: GeneralizedMetropolisHastings.UnnamedParameterDefault{Float32} && d7.default == 0.0
@test typeof(d8) <: GeneralizedMetropolisHastings.UnnamedParameterDefault{Float64} && d8.default == 0.0
@test d2 == d4 #equality operator
@test d6 == d8
@test d1 != d2
@test d1 != d5

###test parameter factory functions for ParameterUnivariate
p1 = parameter("p",Distributions.Normal(),3.0)
p2 = parameter("p",Distributions.Binomial(),0)
p3 = parameter("p",Distributions.Uniform(0.0,11.0))
p4 = parameter("p",Distributions.DiscreteUniform(0,9))
p5 = parameter("p",0.0,11.0,3.0)
p6 = parameter("p",0.0,11.0)
p7 = parameter(Distributions.Uniform(0.0,11.0))
p8 = parameter(Distributions.DiscreteUniform(0,9))
p9 = parameter(0.0,11.0,3.0)
p10 = parameter(0.0,11.0)
@test typeof(p1) <: GeneralizedMetropolisHastings.NamedParameterUnivariate && p1.name == "p" && p1.default == 3.0 && p1.prior == Distributions.Normal()
@test typeof(p2) <: GeneralizedMetropolisHastings.NamedParameterUnivariate && p2.name == "p" && p2.default == 0 && p2.prior == Distributions.Binomial()
@test typeof(p3) <: GeneralizedMetropolisHastings.NamedParameterUnivariate && p3.name == "p" && p3.default == 5.5 && p3.prior == Distributions.Uniform(0.0,11.0)
@test typeof(p4) <: GeneralizedMetropolisHastings.NamedParameterUnivariate && p4.name == "p" && p4.default == 4 && p4.prior == Distributions.DiscreteUniform(0,9)
@test typeof(p5) <: GeneralizedMetropolisHastings.NamedParameterUnivariate && p5.name == "p" && p5.default == 3.0 && p5.prior == Distributions.Uniform(0.0,11.0)
@test typeof(p6) <: GeneralizedMetropolisHastings.NamedParameterUnivariate && p6.name == "p" && p6.default == 5.5 && p6.prior == Distributions.Uniform(0.0,11.0)
@test typeof(p7) <: GeneralizedMetropolisHastings.UnnamedParameterUnivariate && p7.default == 5.5 && p7.prior == Distributions.Uniform(0.0,11.0)
@test typeof(p8) <: GeneralizedMetropolisHastings.UnnamedParameterUnivariate && p8.default == 4 && p8.prior == Distributions.DiscreteUniform(0,9)
@test typeof(p9) <: GeneralizedMetropolisHastings.UnnamedParameterUnivariate && p9.default == 3.0 && p9.prior == Distributions.Uniform(0.0,11.0)
@test typeof(p10) <: GeneralizedMetropolisHastings.UnnamedParameterUnivariate && p10.default == 5.5 && p10.prior == Distributions.Uniform(0.0,11.0)
@test p3 == p6
@test p7 == p10
@test p1 != p2
@test p1 != p7
@test_throws MethodError parameter("u",Distributions.Normal(),3)
@test_throws MethodError parameter("u",Distributions.Binomial(),1.0)
@test_throws MethodError parameter("u",1.0,2.0,3)
@test_throws AssertionError parameter("u",Distributions.Uniform(1.0,2.0),0.0)
@test_throws MethodError parameter(Distributions.Normal(),3)
@test_throws MethodError parameter(Distributions.Binomial(),1.0)
@test_throws MethodError parameter(1.0,3)
@test_throws AssertionError parameter(Distributions.Uniform(1.0,2.0),0.0)

@test length(parameters(["d1","d2"],[3.0,1.0])) == 2
@test length(parameters(["d1","d2"],Float32)) == 2
@test length(parameters([1.0,2.0,3.0])) == 3
@test length(parameters(2,Float32)) == 2
@test length(parameters(5)) == 5
v1 = parameters(["d1","d2"],[Distributions.Normal(),Distributions.Uniform(0.0,11.0)],[0.0,2.0])
v2 = parameters(["d1","d2"],[0.0,5.0],[6.0,10.0],[3.0,7.5])
v3 = parameters(["d1","d2"],[Distributions.Normal(),Distributions.Uniform(0.0,11.0)])
v4 = parameters(["d1","d2"],[0.0,5.0],[6.0,10.0])
v5 = parameters([Distributions.Normal(),Distributions.Uniform(0.0,11.0)],[0.0,2.0])
v6 = parameters([0.0,5.0],[6.0,10.0],[3.0,7.5])
v7 = parameters([Distributions.Normal(),Distributions.Uniform(0.0,11.0)])
v8 = parameters([0.0,5.0],[6.0,10.0])
@test v2 == v4
@test v6 == v8

###test the initvalues functions
@test typeof(GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,d1,Float64)) == Float64 &&  GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,d1,Float64) == 3.0
@test typeof(GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,d3,Float32)) == Float32 && GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,d3,Float32) == 0.0f0
@test typeof(GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,d5,Float32)) == Float32 && GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,d5,Float64) == 3.0
@test typeof(GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,d6,Float32)) == Float32 && GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,d6,Float32) == 0.0f0
@test GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,p1) == 3.0
@test GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,p3) == 5.5
@test GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,p5) == 3.0
@test GeneralizedMetropolisHastings._initvalue(ValuesFromDefault,p8) == 4.0
@test (srand(0) ; GeneralizedMetropolisHastings._initvalue(ValuesFromPrior,p3,Float64)) == (srand(0) ; rand(Distributions.Uniform(0.0,11.0)))
@test (srand(0) ; GeneralizedMetropolisHastings._initvalue(ValuesFromPrior,p4,Float32)) == (srand(0) ; Float32(rand(Distributions.DiscreteUniform(0,9))))
@test (srand(0) ; GeneralizedMetropolisHastings._initvalue(ValuesFromPrior,p3)) == (srand(0) ; rand(Distributions.Uniform(0.0,11.0)))
@test GeneralizedMetropolisHastings._initvalues!(ValuesFromDefault,MCParameter[d1,d4,p1,p4],Array{Float32}(4)) == [3.0f0,0.0f0,3.0f0,4.0f0]
@test initvalues!(ValuesFromDefault,MCParameter[d1,d4,p1,p4],Array{Float64}(4)) == [3.0,0.0,3.0,4.0]
@test initvalues(ValuesFromDefault,MCParameter[d1,d4,p1,p4]) == [3.0,0.0,3.0,4.0]
@test (srand(0) ; initvalues(ValuesFromPrior,MCParameter[d1,d4,p1,p4])) == (srand(0) ; [3.0,0.0,rand(Distributions.Normal()),Float64(rand(Distributions.DiscreteUniform(0,9)))])
@test_throws AssertionError initvalues!(ValuesFromDefault,MCParameter[p1],Array{Float64}(2))

###test the logprior functions
l1 = parameter(Distributions.Uniform(1.0,3.0))
l2 = parameter(2.0)
l3 = parameter(Distributions.DiscreteUniform(1.0,3.0))
@test GeneralizedMetropolisHastings._logprior(l1,2.0) == log(1.0/2.0)
@test GeneralizedMetropolisHastings._logprior(l1,4.0) == -Inf
@test GeneralizedMetropolisHastings._logprior(l2,2.0) == 0.0
@test GeneralizedMetropolisHastings._logprior(l1,2.0f0) == -log(2.0f0)
@test GeneralizedMetropolisHastings._logprior(l1,4.0f0) == -Inf32
@test GeneralizedMetropolisHastings._logprior(l2,2.0f0) == 0.0f0
@test GeneralizedMetropolisHastings._logprior(l3,2.0) == log(1.0/3.0)
@test GeneralizedMetropolisHastings._logprior([l1,l2,l3],[2.0,2.0,2.0]) == log(1.0/2.0) + log(1.0/3.0)
@test logprior(l1,2.0) == log(1.0/2.0)
@test logprior(MCParameter[l1,l2,l3],[2.0,2.0,2.0]) == log(1.0/2.0) + log(1.0/3.0)
@test logprior(MCParameter[l1,l1],[2.0,4.0]) == -Inf
@test logprior(MCParameter[],Float64[]) == 0.0
@test logprior(GeneralizedMetropolisHastings.ParameterDefault[l2,l2],[2.0,4.0]) == 0.0
@test logprior(GeneralizedMetropolisHastings.ParameterUnivariate[l1,l3],[2.0,2.0]) == log(1.0/2.0) + log(1.0/3.0)

println("====================")
println("Test show() function")
println()
show(d1)
show(d4)
show(p1)
show(p3)
show(p4)
show(MCParameter[l1,l2,l3])
show(GeneralizedMetropolisHastings.ParameterUnivariate[l1,l3])
show(GeneralizedMetropolisHastings.UnnamedParameterUnivariate[l1,l3])
show(v1)
show(v4)
show([v1;v4])
show([v1;v5])
println()
println("End  show() function")
println("====================")
println()

nothing
