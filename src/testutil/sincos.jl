###Target functions
function sincos!{T<:AbstractFloat}(r::Matrix{T},t::Vector{T},paras::Vector{T})
    p1 = 2*pi*paras[1]
    p2 = 2*pi*paras[2]
    @simd for i=1:length(t)
        @inbounds r[i,1] = sin(p1*t[i])
    end
    @simd for i=1:length(t)
        @inbounds r[i,2] = cos(p2*t[i])
    end
    r
end

sincos(t::Vector{T},paras::Vector{T}) = sincos!(zeros(T,length(t),2),t,paras)

@inline function _model(s::Symbol,f::Function,n::AbstractString,t::AbstractVector,pv::AbstractVector,variance::AbstractVector,paraminit...)
    p = parameters([:a,:b],paraminit...)
    n = noise(:gaussian,variance)
    d = data(:array,t,applynoise!(n,sincos(t,pv)))
    _model(Val{s},p,d,n,f;name=n)
end

###Create a target function model
function _model(::Type{Val{:sincos}},t::AbstractVector,pv::AbstractVector,variance::AbstractVector,paraminit...)
    _model(:target,sincos,"Sine-Cosine target function",t,pv,variance,paraminit...)
end

###Create a target function model
function _model(::Type{Val{:sincos!}},t::AbstractVector,pv::AbstractVector,variance::AbstractVector,paraminit...)
    _model(:target!,sincos!,"Sine-Cosine in-place target function",t,pv,variance,paraminit...)
end

