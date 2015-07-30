###Geometry calculations per sample type (null, first and second order)
function update_geometry!(m::MCModel,b::MCSample{NullOrder})
  if isfinite(logprior!(m,b))
    loglikelihood!(m,b)
  else
    b.loglikelihood = -Inf
  end
  b
end

function update_geometry!(m::MCModel,g::MCSample{FirstOrder})
  if isfinite(logprior!(m,g)) && all(isfinite(gradlogprior!(m,g)))
    sol::Matrix{Float64} = loglikelihood!(m,g)
    gradient!(m,g,sol)
  else
    b.loglikelihool = -Inf
  end
  g
end

function update_geometry!(m::MCModel,t::MCSample{SecondOrder})
  if isfinite(logprior!(m,t)) && all(isfinite(gradlogprior!(m,t))) && all(isfinite(tensorlogprior!(m,t)))
    sol::Matrix{Float64} = loglikelihood!(m,t);
    grad::Array{Float64,3} = gradient!(m,t,sol);
    tensor!(m,t,sol,grad)
  else
    t.loglikelihood = -Inf
  end
  t
end

###Calculation of the logprior for arbitrary models
logprior!(m::MCModel,s::MCSample) = (s.logprior = sum(map((p,v)->logpdf(p,v),m.parameters.priors,s.values)))
gradlogprior!(m::MCModel,s::MCSample) = (s.gradlogprior = map((p,v)->((logpdf(p,v+m.gradientepsilon)-logpdf(p,v))/m.gradientepsilon),m.parameters.priors,s.values))
tensorlogprior!(m::MCModel,s::MCSample) = (s.tensorlogprior = zeros(s.tensorlogprior)) #does this need to be implemented?

###Calculation of the loglikelihood for arbitrary models and samples
function loglikelihood!(m::MCModel,s::MCSample)

  #call the evaluate function for specific models
  sol::Matrix{Float64} = evaluate(m,s.values)

  #calculate the log-likelihood TODO: separate out the noise model
  s.loglikelihood = loglikelihood(m,sol)

  #return the solution
  sol
end

###Calculation of the gradient from finite differences
function gradient!(m::MCModel,s::MCSample,sol::Matrix{Float64})

  #pre-allocate the array for finite difference results
  fd  = Array(Float64,size(sol,1),size(sol,2),length(s.values))
  np = Array(Float64,length(s.values))

  #calculate finite differences
  for i = 1:length(s.values)
    np = copy(s.values)
    np[i]  = np[i] + m.gradientepsilon
    realepsilon = np[i] - s.values[i]
    fd[:,:,i] = (evaluate(m,np)-sol)/realepsilon
  end
  s.gradloglikelihood = gradloglikelihood(m,sol,fd)

  #return the finite differences solution
  fd
end

function tensor!(m::MCModel,t::TensorSample,sol::Matrix{Float64},grad::Array{Float64,3})
  np = length(t.values)
  for i=1:np
    for j=i:np
      t.tensorloglikelihood[i,j] = tensorvalue(m,grad,i,j)
      t.tensorloglikelihood[j,i] = t.tensorloglikelihood[i,j]
    end
  end
end

function tensor!(m::MCModel,t::ApproximateTensorSample,sol::Matrix{Float64},grad::Array{Float64,3})
  for i=1:size(t.tangentvectors,2)
    t.tangentvectors[:,i] = tangentvector(m,sol,grad)
  end
  t.tensorloglikelihood = (t.tangentvectors*t.tangentvectors')/length(t.tangentvectors)
end




