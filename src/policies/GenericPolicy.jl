#Generic runtime policy
immutable GenericPolicy{I<:InitializeFrom,P<:ProposeFrom,G<:GenerateIndicator,N<:Number} <: RuntimePolicy{N}
  initialize::Type{I}
  propose::Type{P}
  indicate::Type{G}
  numbertype::Type{N}
  GenericPolicy() = new(I,P,G,N)
end

_policy(::Type{Val{:generic}},I::DataType,P::DataType,G::DataType,N::DataType) = GenericPolicy{I,P,G,N}()

function Base.show(io::IO,p::GenericPolicy)
  println(io,"GenericPolicy with policy traits:")
  println(io,"  initialize = ",p.initialize.name.name)
  println(io,"  propose = ",p.propose.name.name)
  println(io,"  indicate = ",p.indicate.name.name)
  println(io,"  number type = ",p.numbertype)
  println(io)
  nothing
end
