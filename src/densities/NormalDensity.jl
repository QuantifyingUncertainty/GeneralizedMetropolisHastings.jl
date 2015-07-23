type NormalDensity <: ProposalDensity
  normal::MvNormal
end

###Constructor with mean and covariance matrix
NormalDensity(μ::Vector{Float64},Σ::Matrix{Float64}) = NormalDensity(MvNormal(μ,Σ))
###Constructor with zero mean and covariance matrix
NormalDensity(Σ::Matrix{Float64}) = NormalDensity(MvNormal(Σ))

###Constructor with mean and std vector
NormalDensity(μ::Vector{Float64},σ::Vector{Float64}) = (length(σ)==1?display_info():nothing ; NormalDensity(MvNormal(μ,σ)))
###Constructor with zero mean and std vector
NormalDensity(σ::Vector{Float64}) = (length(σ)==1?display_info():nothing ; NormalDensity(MvNormal(σ)))

###Constructor with mean and isotropic std value
NormalDensity(μ::Vector{Float64},s::Float64) = NormalDensity(MvNormal(μ,s))
###Constructor with zero mean and isotropic std value
NormalDensity(n::Int,s::Float64 =1.0) = NormalDensity(MvNormal(n,s))

###update the covariance of the density
update_density!(d::NormalDensity,Σ::Matrix{Float64}) = (d.normal = MvNormal(Σ))
update_density!(d::NormalDensity,σ::Vector{Float64}) = (d.normal = MvNormal(σ))
update_density!(d::NormalDensity,s::Float64) = (d.normal = MvNormal(length(d.normal),s))

###update the mean and the covariance of the density
update_density!(d::NormalDensity,μ::Vector{Float64},Σ::Matrix{Float64}) = (d.normal = MvNormal(μ,Σ))
update_density!(d::NormalDensity,μ::Vector{Float64},σ::Vector{Float64}) = (d.normal = MvNormal(μ,σ))
update_density!(d::NormalDensity,μ::Vector{Float64},s::Float64) = (d.normal = MvNormal(μ,s))

###draw a proposal and store it inside the values field of the sample
propose!(d::NormalDensity,s::MCSample) = rand!(d.normal,s.values)

###draw n proposals from the density
propose(d::NormalDensity) = rand(d.normal)
propose(d::NormalDensity,n::Int) = rand(d.normal,n)

###base functionality
==(d1::NormalDensity,d2::NormalDensity) = (mean(d1.normal) == mean(d2.normal) && cov(d1.normal) == cov(d2.normal))

function Base.show(io::IO,d::NormalDensity)
  println(io,"NormalDensity of length ",length(d.normal), " with")
  println(io,"mean: ",mean(d.normal))
  println(io,"covariance: ")
  println(io,cov(d.normal))
  println(io)
  nothing
end

###helper functions (not to be exported)
display_info() = info("for size(σ) == (1,), σ is interpreted as a 1D std vector, not a 2D covariance matrix. Use NormalDensity(μ,s*eye(1)) for the covariance constructor")

