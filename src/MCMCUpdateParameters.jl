
abstract ProposalFunction
immutable ProposalFunctionNormal <: ProposalFunction end
immutable ProposalFunctionAuxiliary <: ProposalFunction end

function propose(sampler::MCSampler,numproposals::Int64; proposalfunction::ProposalFunction = ProposalFunctionAuxiliary())
  propose(sampler,numparas,)
end

function UpdateParameters(Model::TargetOnlyModel, Chain::MarkovChain)

if Chain.Sampler == "MH"
    ### Standard Metropolis-Hastings sampler ###

    #println("Entering UpdateParameters function...")
    ProposalMean       = Chain.Geometry[Chain.SampleIndicator].Parameters
    ProposalCovariance = Chain.ProposalDistribution.StepSize*Chain.ProposalDistribution.ProposalCov

    # Propose the new parameters values
    for j = 1:Chain.NumOfProposals+1
        # Update for all except the indicated sample
        if j!= Chain.SampleIndicator
            # Propose new point
            #println("Proposing new point...")
            Chain.Geometry[j].Parameters    = rand(MvNormal(ProposalMean, ProposalCovariance))
            # Then update the geometry for all proposed points
            #println("Updating geometry for proposed point...")
            Chain.Geometry[j].LL = Model.LLEval( Chain.Geometry[j].Parameters)f
        end
    end

    # Now calculate the transition probabilities
    for j = 1:Chain.NumOfProposals+1
        # Mean centred on current j'th set of parameters
        ProposalMean = Chain.Geometry[j].Parameters
        for i = 1:Chain.NumOfProposals+1
            if i!=j
                # Calculate the probability of proposal
                Chain.Geometry[j].ProposalProbability[i] = logpdf(MvNormal(ProposalMean, ProposalCovariance), Chain.Geometry[i].Parameters)
            end
        end

    end
else
    # Code for other samplers
end

end

function Initialise(Model::TargetOnlyModel, Chain::MarkovChain)

# Initialise the chain by updating geometry for main point
if Chain.Sampler == "MH"
    println("Initialising...")
    println("The chosen sampler is Metropolis-Hastings.")
    Chain.Geometry[Chain.SampleIndicator].LL = Model.LLEval( Chain.Geometry[Chain.SampleIndicator].Parameters)
else
    error("Sampler specified is not a valid option.")
end

end

function UpdateParameters(Model::ODEModel, Chain::MarkovChain)

    if Chain.Sampler == "MH"

        ############################################
        ### Standard Metropolis-Hastings sampler ###
        ############################################

        # Propose the new parameters values
        for j = 1:Chain.NumOfProposals+1
            # Update for all except the indicated sample
            if j!= Chain.SampleIndicator

                # Propose new point
                #println("Proposing new point...")
                Chain.Geometry[j].Parameters    = Chain.Geometry[Chain.SampleIndicator].Parameters + rand(Chain.ProposalDistribution.Density)

                # Calculate the log-likelihood for the ODE model
                UpdateODEGeometry!(Model, Chain, j, false, false) # Update the geometry for current proposal

            end
        end

        # Symmetric proposal so no need to update probabilities


    elseif Chain.Sampler == "SmMALA"

        ################################
        ### Simplified mMALA sampler ###
        ################################

        # Propose the new parameters values
        for PropNum = 1:Chain.NumOfProposals+1
            # Update for all except the indicated sample
            if PropNum != Chain.SampleIndicator

                # Propose new point
                #println("Proposing new point...")
                Chain.Geometry[PropNum].Parameters    = rand(Chain.ProposalDistribution[Chain.SampleIndicator].Density)

                # Calculate the log-likelihood, gradient and FI for the ODE model
                UpdateODEGeometry!(Model, Chain, PropNum, true, true)

                # Update the proposal mean to current parameters
                Chain.ProposalDistribution[PropNum].Density = MvNormal( vec( Chain.Geometry[PropNum].Parameters + (Chain.ProposalDistribution[PropNum].StepSize^2/2)*(Chain.Geometry[PropNum].HessianLL\Chain.Geometry[PropNum].GradLL) ),
                                                                        inv(Chain.Geometry[PropNum].HessianLL).*Chain.ProposalDistribution[PropNum].StepSize^2 )
                # Use canonical form of multivariate Normal proposal, with J = inv(Cov) = FI/StepSize and h = inv(Cov)*mu
                #Chain.Geometry[PropNum].ProposalDistribution = MvNormalCanon( vec( Chain.Geometry[PropNum].HessianLL*(Chain.Geometry[PropNum].Parameters./Chain.Geometry[PropNum].StepSize^2) + (Chain.Geometry[PropNum].GradLL./2) ),
                #                                                              (Chain.Geometry[PropNum].HessianLL./Chain.Geometry[PropNum].StepSize^2) )

            end
        end

        # Now calculate the transition probabilities
        for j = 1:Chain.NumOfProposals+1

            for i = 1:Chain.NumOfProposals+1
                if i!=j
                    Chain.Geometry[j].ProposalProbability[i] = logpdf(Chain.ProposalDistribution[j].Density, Chain.Geometry[i].Parameters)
                end
            end

        end

    elseif Chain.Sampler == "TrSmMALA"

        #############################################
        ### Trust Region Simplified mMALA sampler ###
        #############################################

        # Propose the new parameters values
        for PropNum = 1:Chain.NumOfProposals+1
            # Update for all except the indicated sample
            if PropNum != Chain.SampleIndicator

                # Propose new point
                #println("Proposing new point...")
                Chain.Geometry[PropNum].Parameters    = rand(Chain.ProposalDistribution[Chain.SampleIndicator].Density)

                # Calculate the log-likelihood, gradient and FI for the ODE model
                UpdateODEGeometry!(Model, Chain, PropNum, true, true)

                # Update the proposal based on current parameters
                Chain.ProposalDistribution[PropNum].Density = MvNormal( vec( Chain.Geometry[PropNum].Parameters + ((Chain.Geometry[PropNum].HessianLL + eye(Model.NumOfParas)*(1/Chain.ProposalDistribution[PropNum].StepSize))\Chain.Geometry[PropNum].GradLL) ),
                                                                         inv(Chain.Geometry[PropNum].HessianLL + eye(Model.NumOfParas)*(1/Chain.ProposalDistribution[PropNum].StepSize)) )

                #println("Current Paras:")
                #println(Chain.Geometry[PropNum].Parameters)
                #println("Proposal:")
                #println(Chain.Geometry[PropNum].ProposalDistribution )
            end
        end

        # Now calculate the transition probabilities
        for j = 1:Chain.NumOfProposals+1

            for i = 1:Chain.NumOfProposals+1
                if i!=j
                    Chain.Geometry[j].ProposalProbability[i] = logpdf(Chain.ProposalDistribution[j].Density, Chain.Geometry[i].Parameters)
                end
            end

        end

    elseif Chain.Sampler == "TrSmMALA_Random"

        ####################################################################
        ### Trust Region Simplified mMALA sampler with Randomised Metric ###
        ####################################################################

        Success = UpdateODEGeometryRandom!(Model, Chain, Chain.SampleIndicator, true) # Update the first proposal

        # Propose the new parameters values
        for PropNum = 1:Chain.NumOfProposals+1
            # Update for all except the indicated sample
            if PropNum != Chain.SampleIndicator

                # Propose new point
                #println("Proposing new point...")
                Chain.Geometry[PropNum].Parameters    = rand(Chain.ProposalDistribution[Chain.SampleIndicator].Density)

                # Calculate the log-likelihood, gradient and FI for the ODE model
                UpdateODEGeometryRandom!(Model, Chain, PropNum, true)

                # Update the proposal based on proposed parameters
                Chain.ProposalDistribution[PropNum].Density = MvNormal( vec( Chain.Geometry[PropNum].Parameters + ((Chain.Geometry[PropNum].HessianLL + eye(Model.NumOfParas)*(1/Chain.ProposalDistribution[PropNum].StepSize))\Chain.Geometry[PropNum].GradLL) ),
                                                                         inv(Chain.Geometry[PropNum].HessianLL + eye(Model.NumOfParas)*(1/Chain.ProposalDistribution[PropNum].StepSize)) )

                #println("Current Paras:")
                #println(Chain.Geometry[PropNum].Parameters)
                #println("Proposal:")
                #println(Chain.Geometry[PropNum].ProposalDistribution )
            end
        end

        # Now calculate the transition probabilities
        for j = 1:Chain.NumOfProposals+1

            for i = 1:Chain.NumOfProposals+1
                if i!=j
                    Chain.Geometry[j].ProposalProbability[i] = logpdf(Chain.ProposalDistribution[j].Density, Chain.Geometry[i].Parameters)
                end
            end

        end

    elseif Chain.Sampler == "AdaptiveMH"

        ############################################
        ### Adaptive Metropolis-Hastings sampler ###
        ############################################

        CurIter = Chain.CurrentIteration;

        # Update the proposal mean and covariance
        #println(CurIter)
        if CurIter == 1
            Chain.ProposalDistribution.RunningCov  = (1/CurIter)*(Chain.Geometry[Chain.SampleIndicator].Parameters-Chain.ProposalDistribution.RunningMean)*(Chain.Geometry[Chain.SampleIndicator].Parameters-Chain.ProposalDistribution.RunningMean)';
            Chain.ProposalDistribution.RunningMean = ((CurIter-1)*Chain.ProposalDistribution.RunningMean + Chain.Geometry[Chain.SampleIndicator].Parameters)/CurIter;
            Chain.ProposalDistribution.Density     = MvNormal(zeros(Model.NumOfParas),
                                                              eye(Model.NumOfParas)*(Chain.ProposalDistribution.StepSize^2))
        else
            Chain.ProposalDistribution.RunningCov  = ((CurIter-2)/(CurIter-1))*Chain.ProposalDistribution.RunningCov + (1/CurIter)*(Chain.Geometry[Chain.SampleIndicator].Parameters-Chain.ProposalDistribution.RunningMean)*(Chain.Geometry[Chain.SampleIndicator].Parameters-Chain.ProposalDistribution.RunningMean)';
            Chain.ProposalDistribution.RunningMean = ((CurIter-1)*Chain.ProposalDistribution.RunningMean + Chain.Geometry[Chain.SampleIndicator].Parameters)/CurIter;
            Chain.ProposalDistribution.Density     = MvNormal(zeros(Model.NumOfParas),
                                                              Chain.ProposalDistribution.RunningCov + eye(Model.NumOfParas)*(Chain.ProposalDistribution.StepSize^2))
        end


        # Propose the new parameters values
        for PropNum = 1:Chain.NumOfProposals+1
            # Update for all except the indicated sample
            if PropNum != Chain.SampleIndicator

                # Propose new point
                #println("Proposing new point...")
                Chain.Geometry[PropNum].Parameters    = Chain.Geometry[Chain.SampleIndicator].Parameters + rand(Chain.ProposalDistribution.Density) # using same proposal

                # Calculate the log-likelihood, gradient and FI for the ODE model
                UpdateODEGeometry!(Model, Chain, PropNum, false, false)

                #println("Current Paras:")
                #println(Chain.Geometry[PropNum].Parameters)
                #println("Proposal:")
                #println(Chain.Geometry[PropNum].ProposalDistribution )
            end
        end

        # Symmetric proposal so no need to calculate the transition probabilities!

    else
      # Code for other samplers
    end
end

function Initialise(Model::ODEModel, Chain::MarkovChain)
    # Initialise the chain by updating geometry for main point
    if Chain.Sampler == "MH"

        println("Initialising...")
        println("The chosen sampler is Metropolis-Hastings.")

        # Calculate the log-likelihood for the ODE model
        UpdateODEGeometry!(Model, Chain, 1, false, false) # Update the first proposal

        # Update the proposal mean to current parameters
        Chain.ProposalDistribution.Density = MvNormal(zeros(Model.NumOfParas), (Chain.ProposalDistribution.StepSize^2)*Chain.ProposalDistribution.ProposalCov)

    elseif Chain.Sampler == "MALA"


    elseif Chain.Sampler == "SmMALA"

        println("Initialising...")
        println("The chosen sampler is Simplified mMALA.")

        # Calculate the log-likelihood, gradient and FI for the ODE model
        UpdateODEGeometry!(Model, Chain, 1, true, true) # Update the first proposal

        # Update the proposal based on current parameters
        Chain.ProposalDistribution[1].Density = MvNormal( vec( Chain.Geometry[1].Parameters + (Chain.ProposalDistribution[1].StepSize^2/2)*(Chain.Geometry[1].HessianLL\Chain.Geometry[1].GradLL) ),
                                                          inv(Chain.Geometry[1].HessianLL).*Chain.ProposalDistribution[1].StepSize^2 )

        # Use canonical form of multivariate Normal proposal, with J = inv(Cov) = FI/StepSize and h = inv(Cov)*mu
        #Chain.Geometry[1].ProposalDistribution = MvNormalCanon( vec( Chain.Geometry[1].HessianLL*(Chain.Geometry[1].Parameters./Chain.Geometry[1].StepSize^2) + (Chain.Geometry[1].GradLL./2) ),
        #                                                        (Chain.Geometry[1].HessianLL./Chain.Geometry[1].StepSize^2) )

    elseif Chain.Sampler == "TrSmMALA"

        println("Initialising...")
        println("The chosen sampler is Trust Region Simplified mMALA.")

        # Calculate the log-likelihood, gradient and FI for the ODE model
        Success = UpdateODEGeometry!(Model, Chain, 1, true, true) # Update the first proposal

        # Update the proposal based on current parameters
        Chain.ProposalDistribution[1].Density = MvNormal( vec( Chain.Geometry[1].Parameters + ((Chain.Geometry[1].HessianLL + eye(Model.NumOfParas)*(1/Chain.ProposalDistribution[1].StepSize))\Chain.Geometry[1].GradLL) ),
                                                          inv(Chain.Geometry[1].HessianLL + eye(Model.NumOfParas)*(1/Chain.ProposalDistribution[1].StepSize)) )

        #println("Current Paras:")
        #println(Chain.Geometry[1].Parameters)
        #println("Proposal:")
        #println(Chain.Geometry[1].ProposalDistribution)

    elseif Chain.Sampler == "TrSmMALA_Random"

        println("Initialising...")
        println("The chosen sampler is Trust Region Simplified mMALA with Randomised Metric.")

        # Calculate the log-likelihood, gradient and FI for the ODE model
        Success = UpdateODEGeometryRandom!(Model, Chain, 1, true) # Update the first proposal

        # Update the proposal based on current parameters
        Chain.ProposalDistribution[1].Density = MvNormal( vec( Chain.Geometry[1].Parameters + ((Chain.Geometry[1].HessianLL + eye(Model.NumOfParas)*(1/Chain.ProposalDistribution[1].StepSize))\Chain.Geometry[1].GradLL) ),
                                                          inv(Chain.Geometry[1].HessianLL + eye(Model.NumOfParas)*(1/Chain.ProposalDistribution[1].StepSize)) )

        #println("Current Paras:")
        #println(Chain.Geometry[1].Parameters)
        #println("Proposal:")
        #println(Chain.Geometry[1].ProposalDistribution)

    elseif Chain.Sampler == "AdaptiveMH"

        println("Initialising...")
        println("The chosen sampler is Adaptive MH.")

        # Calculate the log-likelihood for the ODE model
        Success = UpdateODEGeometry!(Model, Chain, 1, false, false) # Update the first proposal
        # Update the proposal based on current parameters
        Chain.ProposalDistribution.Density = MvNormal( zeros(Model.NumOfParas), eye(Model.NumOfParas)*(Chain.ProposalDistribution.StepSize^2))

    else
        error("Sampler specified is not a valid option.")
    end
end

function UpdateODEGeometry!(Model::ODEModel,Sampler::AbstractSampler,Chain::MarkovChain,PropNum::Int64,CalculateGradLL::Bool,CalculateHessianLL::Bool)

  if !calculat_priors!(Model,Sampler,Chain,PropNum)
    return false
  end

  Sol = calculate_loglikelihood!(Model,Sampler,Chain,PropNum)

  # Only if flag is "true" for calculating gradients
  if CalculateGradLL
    Sol_fd = calculate_gradient!(Model,Sampler,Chain.ProposalDis)

    if CalculateHessianLL
      calculate_tensor!(Model,Chain,PropNum,Sol_fd)
    end
  end
end

function propose(sampler::AbstractSampler,chain::MarkovChain,propnum::Int64,pf::ProposalFunctionNormal)

end

function propose(sampler::AbstractSampler,chain::MarkovChain,propnum::Int64,pf::ProposalFunctionAuxiliary)

end


