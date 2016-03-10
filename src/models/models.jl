abstract AbstractModel

### Factory function
model(s::Symbol,args...;keyargs...) = _model(Val{s},args...;keyargs...)

numparas(m::AbstractModel) = length(m.parameters)
parameters(m::AbstractModel) = m.parameters

initialize!(i::InitializeFrom,m::AbstractModel,s::AbstractSample) = (initvalues!(i,parameters(m),s.values) ; s)

###Calculation of logprior is the same for all models and samples
function geometry!(m::AbstractModel,s::AbstractSample)
    logprior!(s.logprior,parameters(m),s.values)
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
    r::AbstractArray = evaluate!(m,s.values[:,i])
    s.loglikelihood[i] = loglikelihood(m,r)
end

function geometry!(m::AbstractModel,s::AbstractSample{FirstOrder},i::Integer)
    r::AbstractArray = evaluate!(m,s.values[:,i])
    s.loglikelihood[i] = loglikelihood(m,r)
    #gradlogprior!(m,s,r,i)
    #gradloglikelihood!(m,s,r,i)
end

function geometry!(m::AbstractModel,s::AbstractSample{SecondOrder},i::Integer)
    r::AbstractArray = evaluate!(m,s.values[:,i])
    s.loglikelihood[i] = loglikelihood(m,r)
    #gradlogprior!(m,s,r,i)
    #gradloglikelihood!(m,s,r,i)
    #tensorlogprior!(m,s,r,i)
    #tensorloglikelihood!(m,s,r,i)
end

###in general, we do not know how to evaluate models to get model data
evaluate!(m::AbstractModel,v::AbstractVector) = throw(MethodError(evaluate!, (m,v)))

###we do have a generic way to calculate the loglikelihood
loglikelihood(m::AbstractModel,modeldata::AbstractArray) = loglikelihood(noisemodel(m),measurements(m),modeldata)

# ###Calculation of the gradient from finite differences
# function gradient!(m::MCModel,s::MCSample,sol::Matrix{Float64})

#   #pre-allocate the array for finite difference results
#   fd  = Array(Float64,size(sol,1),size(sol,2),length(s.values))
#   np = Array(Float64,length(s.values))

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
# end

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




