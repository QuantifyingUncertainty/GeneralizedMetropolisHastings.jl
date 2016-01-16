#Generic runtime policy
immutable GenericPolicy{T<:Number} <: RuntimePolicy{T}
  initialize::InitializeFrom
  propose::ProposeFrom
  indicate::GenerateIndicator
  numbertype::Type{T}
end

_policy{T<:Number}(::Type{Val{:generic}},i::InitializeFrom)

GenericPolicy() = GenericPolicy(v,i,nproposals>1?ProposalFromAuxiliary():ProposalFromIndicator())
GenericPolicy(v::ValuesFrom,nproposals::Integer) = GenericPolicy(v,IndicatorStationary(),nproposals)
GenericPolicy(nproposals::Integer) = GenericPolicy(ValuesFromPrior(),nproposals)

function Base.show(io::IO,p::GenericPolicy,s::AbstractString ="")
  println(io,s,"GenericPolicy with following policy types:")
  println(io,s,"  initialize = ",typeof(p.initialize))
  println(io,s,"  indicate = ",typeof(p.indicate))
  println(io,s,"  propose = ",typeof(p.propose))
  println(io)
  nothing
end
