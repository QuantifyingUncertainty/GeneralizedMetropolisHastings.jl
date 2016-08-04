###Types to hold intermediate results of geometry calculations
abstract AbstractGeometry{O<:DerivativeOrder}

type BaseGeometry <: AbstractGeometry{ZeroOrder}
end

type GradientGeometry <: AbstractGeometry{FirstOrder}
end

type TensorGeometry <: AbstractGeometry{SecondOrder}
end

### Functions that, given a model and associated evaluation function, calculate the geometry of the parameter space

### Calculate the geometry over all samples
function geometry!(m::AbstractModel,s::AbstractSample)
    ###Calculation of logprior is the same for all types of models and samples
    logprior!(s.logprior,parameters(m),s.values)
    ###For all samaples falling inside the prior, calculate the geometry
    for i=1:numsamples(s)
        if isfinite(s.logprior[i])
            geometry!(m,s,i)
        else
            s.loglikelihood[i] = -Inf
        end
    end
    return s
end



function geometry!(m::AbstractModel,s::AbstractSample{ZeroOrder},i::Integer)
    r = evaluate!(m,s.values[:,i])
    s.loglikelihood[i] = loglikelihood(m,r)
end

# function geometry!(m::AbstractModel,s::AbstractSample{FirstOrder},i::Integer)
#     r = evaluate!(m,s.values[:,i])
#     s.loglikelihood[i] = loglikelihood(m,r)
#     gradlogprior!(m,s,r,i)
#     gradloglikelihood!(m,s,r,i)
# end

# function geometry!(m::AbstractModel,s::AbstractSample{SecondOrder},i::Integer)
#     r::AbstractArray = evaluate!(m,s.values[:,i])
#     s.loglikelihood[i] = loglikelihood(m,r)
#     gradlogprior!(m,s,r,i)
#     gradloglikelihood!(m,s,r,i)j
#     tensorlogprior!(m,s,r,i)
#     tensorloglikelihood!(m,s,r,i)
# end

### Calculate the likelihood of a set of data values given the measurement data
@inline loglikelihood(m::AbstractModel,data::AbstractArray) = loglikelihood(noisemodel(m),measurements(m),data)


# ###Calculation of the gradient from finite differences
# function gradient(m::AbstractModel,s::AbstractSample{GradientOrder},g::AbstractGeometry{GradientOrder})
#     #pre-allocate the array for finite difference results
#     #fds  = similar(evaldata,size(evaldata)...,nparas(s))
# end

#gradlogprior!(m::ODEModel,s::Sample{}) = (s.gradlogprior = map((p,v)->((logpdf(p,v+m.gradientepsilon)-logpdf(p,v))/m.gradientepsilon),m.parameters.priors,s.values))

###Calculate the loglikelihood function

###Calculate the gradient of the loglikelihood
#gradienthelper(m::ODEModel,data::Matrix{Float64},sol::Matrix{Float64},grad::Array{Float64,3}) = vec(sum((grad.*(data-sol)./m.variance),(1,2))) #.* does automatic broadcast(), summing result along first 2 dimensions
#gradloglikelihood(m::ODEModel,sol::Matrix{Float64},grad::Array{Float64,3}) = gradienthelper(m,m.measurements,sol,grad)

###Calculate the metric tensor
#tensorvalue(m::ODEModel,grad::Array{Float64,3},i::Int,j::Int) = sum(grad[:,:,i].*grad[:,:,j]./m.variance)

###Generate pseudodata for the approximate metric tensor calculations
#pseudodata(m::ODEModel,d::AbstractArray) = applynoise!(m.noise,d)
#tangentvector(m::ODEModel,sol::Matrix{Float64},grad::Array{Float64,3}) = gradienthelper(m,pseudodata(m,sol),sol,grad)

###Calculate the gradient of the loglikelihood
#gradienthelper(m::ODEModel,data::Matrix{Float64},sol::Matrix{Float64},grad::Array{Float64,3}) = vec(sum((grad.*(data-sol)./m.variance),(1,2))) #.* does automatic broadcast(), summing result along first 2 dimensions
#gradloglikelihood(m::ODEModel,sol::Matrix{Float64},grad::Array{Float64,3}) = gradienthelper(m,m.measurements,sol,grad)


#   #calculate finite differences
#   for i = 1:length(s.values)
#     np = copy(s.values)
#     np[i]  = np[i] + m.gradientepsilon
#     realepsilon = np[i] - s.values[i]
#     fd[:,:,i] = (evaluate(m,np)-sol)/realepsilon
#   end
#   s.gradloglikelihood = gradloglikelihood(m,sol,fd)

#   #return the finite differences solution
#   fd

# function tensor!(m::MCModel,t::TensorSample,sol::Matrix{Float64},grad::Array{Float64,3})
#   np = length(t.values)
#   for i=1:np
#     for j=i:np
#       t.tensorloglikelihood[i,j] = tensorvalue(m,grad,i,j)
#       t.tensorloglikelihood[j,i] = t.tensorloglikelihood[i,j]
#     end
#   end
# end

# function tensor!(m::MCModel,t::ApproximateTensorSample,sol::Matrix{Float64},grad::Array{Float64,3})
#   for i=1:size(t.tangentvectors,2)
#     t.tangentvectors[:,i] = tangentvector(m,sol,grad)
#   end
#   t.tensorloglikelihood = (t.tangentvectors*t.tangentvectors')/length(t.tangentvectors)
# end
