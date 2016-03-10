### ChainStandard

# Stores the output of a Monte Carlo run
type ChainStandard{N<:Number,T<:AbstractFloat,V<:AbstractVector,A<:AbstractArray} <: AbstractChain
    values::A
    loglikelihood::V
    accepted::Int
    proposed::Int
    runtime::Real
    ChainStandard(v::AbstractArray{N},ll::AbstractVector{T},a::Int,p::Int,r::Real) = new(v,ll,a,p,r)
end

###Factory function
@inline function _chain{N<:Number,T<:AbstractFloat}(::Type{Val{:standard}},nparas::Int,nsamples::Int,::Type{N},::Type{T})
    v = zeros(N,_valuestuple(nparas,nsamples))
    ChainStandard{N,T,Vector,Array}(v,Vector{T}(nsamples),0,0,0.0)
end

function store!(c::ChainStandard,s::AbstractSample,j::Int)
    @assert numparas(c) == numparas(s) && j <= numsamples(s)
    if c.proposed < numsamples(c)
        c.proposed += 1
        @simd for i=1:numparas(s)
            @inbounds c.values[i,c.proposed] = s.values[i,j]
        end
        c.loglikelihood[c.proposed] = s.loglikelihood[j]
    else
        warn("Chain is full. Additional results cannot be stored")
    end
end

function show(io::IO,c::ChainStandard)
  println("ChainStandard with numparas = $(numparas(c)) and numsamples = $(numsamples(c))")
  println("Samples proposed = $(c.proposed), samples accepted = $(c.accepted), acceptance rate = $(c.accepted/c.proposed)")
  println("Total runtime = $(c.runtime)")
  println("Additional fields: :values, :loglikelihood")
end
