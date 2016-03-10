### ChainGradient

# Stores the output of a Monte Carlo run
type ChainGradient{N<:Number,T<:AbstractFloat,V<:AbstractVector,A<:AbstractArray} <: AbstractChain
    values::A
    loglikelihood::V
    gradloglikelihood::A
    accepted::Int
    proposed::Int
    runtime::Real
    ChainGradient(v::AbstractArray{N},ll::AbstractVector{T},gl::AbstractArray{T},a::Int,p::Int,r::Real) = new(v,ll,gl,a,p,r)
end

@inline function _chain{N<:Number,T<:AbstractFloat}(::Type{Val{:gradient}},nparas::Int,nsamples::Int,::Type{N},::Type{T})
    v = zeros(N,_valuestuple(nparas,nsamples))
    ChainGradient{N,T,Vector,Array}(v,Vector{T}(nsamples),similar(v,T),0,0,0.0)
end

function store!{O<:GradientOrder}(c::ChainGradient,s::AbstractSample{O},j::Int)
    @assert numparas(c) == numparas(s) && j <= numsamples(s)
    if c.proposed < numsamples(c)
        c.proposed += 1
        @simd for i=1:numparas(s)
            @inbounds c.values[i,c.proposed] = s.values[i,j]
        end
        @simd for i=1:numparas(s)
            @inbounds c.gradloglikelihood[i,c.proposed] = s.gradloglikelihood[i,j]
        end
        c.loglikelihood[c.proposed] = s.loglikelihood[j]
    else
        warn("Chain is full. Additional results cannot be stored")
    end
end

function show(io::IO,c::ChainGradient)
  println("ChainGradient with numparas = $(numparas(c)) and numsamples = $(numsamples(c))")
  println("Samples proposed = $(c.proposed), samples accepted = $(c.accepted), acceptance rate = $(c.accepted/c.proposed)")
  println("Total runtime = $(c.runtime)")
  println("Additional fields: :values, :loglikelihood, :gradloglikelihood")
end
