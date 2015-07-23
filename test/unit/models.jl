g1(x::Float64,;μ =0.0,σ =1.0) = logpdf(Normal(μ,σ),x)
g2(x::Array{Float64,1};μ =zeros(2),Σ = eye(2)) = logpdf(MvNormal(μ,Σ),x)

t1 = TargetModel("SomeName",ModelParameters(),()->(),[])
show(t1)
