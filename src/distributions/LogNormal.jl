#Implement the interface of the Distributions package
#for a multi-variate log-normal distribution
immutable LogNormal <: Distributions.ContinuousMultivariateDistribution
  normal::Normal
end

###Constructors matching the Distributions.MvNormal constructors
LogNormal(μ::Vector{Float64},Σ::Matrix{Float64}) = LogNormal(Normal(μ,Σ))
LogNormal(Σ::Matrix{Float64}) = LogNormal(Normal(Σ))
LogNormal(n::Int,s::Float64) = LogNormal(Normal(n,s))

function Distributions.insupport{T<:Real}(l::LogNormal,x::AbstractVector{T})
  for i=1:length(x)
    0.0<x[i]<Inf?continue:(return false)
  end
  return true
end

function Distributions.length(l::LogNormal)
  Distributions.length(l.normal)
end

function Distributions._rand!{T<:Real}(l::LogNormal,x::AbstractVector{T})
  Distributions._rand!(l.normal,x)
  Distributions.exp!(x)
end

function Distributions._pdf{T<:Real}(l::LogNormal,x::AbstractVector{T})
  if Distributions.insupport(l,x)
    Distributions._pdf(l.normal,log(x))/prod(x)
  else
    zeros(x)
  end
end

function Distributions._logpdf{T<:Real}(l::LogNormal,x::AbstractVector{T})
  if Distributions.insupport(l,x)
    Distributions._logpdf(l.normal,log(x)) - sum(log(x))
  else
    fill(-Inf,length(x))
  end
end

Distributions.mean(l::LogNormal) = exp(Distributions.mean(l.normal) + 0.5*Distributions.var(l.normal))
Distributions.cov(l::LogNormal) = Distributions.mean(l)*Distributions.mean(l)'.*(exp(cov(l.normal))-1)
Distributions.var(l::LogNormal) = diag(Distributions.cov(l))
Distributions.entropy(l::LogNormal) = 0.5 + 0.5*log(2*π*Distributions.var(l.normal)) + Distributions.mean(l.normal)

function Base.show(io::IO,l::LogNormal)
  println(io,"LogNormal(")
  println(io,"dim: ",length(l))
  println(io,"μ: ",mean(l.normal))
  println(io,"Σ: ")
  show(io,cov(l.normal))
  println("")
  println(")")
end





