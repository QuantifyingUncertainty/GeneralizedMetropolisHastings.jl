###Geometry calculations per sample type (null, first and second order)
function update_geometry!(m::MCModel,b::MCSample{NullOrder})
  if isfinite(logprior!(m,b))
    loglikelihood!(m,b)
  end
end

function update_geometry!(m::MCModel,g::MCSample{FirstOrder})
  if isfinite(logprior!(m,g))
    sol::Matrix{Float64} = loglikelihood!(m,g)
    gradient!(m,g,res)
  end
end

function update_geometry!(m::MCModel,t::MCSample{SecondOrder})
  if isfinite(logprior!(m,t))
    sol::Matrix{Float64} = loglikelihood!(m,t);
    grad::Array{Float64,3} = gradient!(m,t,sol);
    tensor!(m,t,grad)
  end
end

###Calculation of the logprior for arbitrary models
logprior!(m::MCModel,s::MCSample) = (s.logprior = sum(map((p,v)->logpdf(p,v),m.parameters.priors,s.values)))

###Calculation of the loglikelihood for arbitrary models and samples
function loglikelihood!(m::MCModel,s::MCSample)

  #you should implement the evaluate function for specific models
  sol::Matrix{Float64} = evaluate(m,s.values)

  #calculate the log-likelihood TODO: separate out the noise model
  s.loglikelihood = loglikelihood(m,sol)

  #return the solution
  sol
end

function gradient!(m::MCModel,s::MCSample,sol::Matrix{Float64})

  #pre-allocate the result array
  fd  = Array(Float64,size(sol,1),size(sol,2),length(s.values))
  paras = Array(Float64,length(s.values))

  e = m.abstol*100.0

  for i = 1:length(s.values)
    newparas = copy(s.values)
    newparas[i]  = newparas[i] + e
    true_e = newparas[i] - s.values[i]
    grad[:,:,i] = (evaluate(m,newparas) - sol)/true_e
    s.gradloglikelihood = sum(((m.measurements-sol)./m.variance).*grad(:,:,i))

    ###TODO FROM HERE
    # Add gradient of prior
    GradLogPrior = (logpdf(model.Priors[i],chain.Geometry[propnum].Parameters[i]+Epsilon) - logpdf(model.Priors[i],chain.Geometry[propnum].Parameters[i]) )./Epsilon
    chain.Geometry[propnum].GradLL[i] = chain.Geometry[propnum].GradLL[i] + GradLogPrior
  end

  #return the finite differences solution
  fd
end




################################################################################
###Calculation of the log likelihood for TargetModels and for arbitrary samplers
################################################################################
function loglikelihood!(m::TargetModel,s::MCSample)

  #calculate the log-likelihood using the target function
  s.loglikelihood = m.target(s.values,m.custom)
end

#############################################################################################
###Calculation of the gradient of the log likelihood for ODEModels and for arbitrary samplers
#############################################################################################

function calculate_tensor!(model::ODEModel,sampler::AbstractSampler,chain::MarkovChain,propnum::Int64,sol_fd)

  # Calculation of the hessian of the log likelihood for ODEModels and for arbitrary samplers
  chain.Geometry[propnum].HessianLL = zeros(model.NumOfParas,model.NumOfParas)

  for i = 1:model.NumOfParas
    for j = i:model.NumOfParas
      for StatesNum in Model.ObservedStates
        for t = 1:model.NumOfTimePoints
          chain.Geometry[propnum].HessianLL[i,j] = chain.Geometry[PropNum].HessianLL[i,j] + (sol_fd[t,StatesNum,i]/model.DataNoiseVariance[t,StatesNum]*sol_fd[t,StatesNum,j])
        end
      end
      chain.Geometry[propnum].HessianLL[j,i] = chain.Geometry[propnum].HessianLL[i,j];
    end
  end

end

function calculate_tensor!(model::ODEModel,sampler::ProposalDistributionSmMALARandom,chain::MarkovChain,propnum::Int64,sol_fd)

  # Now sample the some pseudo data and calculate LL
    Chain.ProposalDistribution[PropNum].TangentVectors = zeros(Model.NumOfParas, Chain.ProposalDistribution[PropNum].NumberOfVectors);
    Pseudodata  = zeros(Model.NumOfTimePoints, length(Model.ObservedStates));

    # For each auxiliary tangent vector
    for TangNum = 1:Chain.ProposalDistribution[PropNum].NumberOfVectors

        # Sample pseudodata and calculate log-likelihood of sampling it
        for s in Model.ObservedStates
            for t = 1:Model.NumOfTimePoints
                Pseudodata[t,s] = rand( Normal(Sol[t,s], sqrt(Model.DataNoiseVariance[t,s])) )
            end
        end

        # Now calculate the actual tangent vectors, i.e. derivative of LL using pseudodata
        Temp = 0.0;

        for i = 1:Model.NumOfParas
            for s in Model.ObservedStates
                Temp = ((Sol[:,s]-Pseudodata[:,s])./Model.DataNoiseVariance[:,s])'*Sol_fd[:,s,i]
                Chain.ProposalDistribution[PropNum].TangentVectors[i,TangNum] = Chain.ProposalDistribution[PropNum].TangentVectors[i,TangNum] - Temp[1]
            end
        end

    end

    # Create approximate metric tensor
    Chain.Geometry[PropNum].HessianLL = (1/Chain.ProposalDistribution[PropNum].NumberOfVectors)*(Chain.ProposalDistribution[PropNum].TangentVectors*Chain.ProposalDistribution[PropNum].TangentVectors');

end




