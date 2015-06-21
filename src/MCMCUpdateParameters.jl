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
            if Model.NumOfParas == 1
                # 1 dimensional case
                Chain.Geometry[j].Parameters[1] = rand(Normal(ProposalMean[1], ProposalCovariance[1]))
            else
                # Multi-dimensional case
                Chain.Geometry[j].Parameters    = rand(MvNormal(ProposalMean, ProposalCovariance))
            end
            # Then update the geometry for all proposed points
            #println("Updating geometry for proposed point...")
            Chain.Geometry[j].LL = Model.LLEval( Chain.Geometry[j].Parameters)
        end
    end

    # Now calculate the transition probabilities
    for j = 1:Chain.NumOfProposals+1
        # Mean centred on current j'th set of parameters
        ProposalMean = Chain.Geometry[j].Parameters
        for i = 1:Chain.NumOfProposals+1
            if i!=j
                if Model.NumOfParas == 1
                    # 1 dimensional case
                    # Calculate the probability of proposing i from j
                    Chain.Geometry[j].ProposalProbability[i] = logpdf(Normal(ProposalMean[1], ProposalCovariance[1]), Chain.Geometry[i].Parameters[1])
                else
                    # Calculate the probability of proposal
                    Chain.Geometry[j].ProposalProbability[i] = logpdf(MvNormal(ProposalMean, ProposalCovariance), Chain.Geometry[i].Parameters)
                end
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
                if Model.NumOfParas == 1
                    # 1 dimensional case
                    Chain.Geometry[j].Parameters[1] = Chain.Geometry[Chain.SampleIndicator].Parameters + rand(Chain.ProposalDistribution.Density)
                else
                    # Multi-dimensional case
                    Chain.Geometry[j].Parameters    = Chain.Geometry[Chain.SampleIndicator].Parameters + rand(Chain.ProposalDistribution.Density)
                end

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
                if Model.NumOfParas == 1
                    # 1 dimensional case
                    Chain.Geometry[PropNum].Parameters[1] = rand(Chain.ProposalDistribution[Chain.SampleIndicator].Density) # NEEDS UPDATED!!!
                else
                    # Multi-dimensional case
                    Chain.Geometry[PropNum].Parameters    = rand(Chain.ProposalDistribution[Chain.SampleIndicator].Density)
                end

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
                    if Model.NumOfParas == 1
                        # 1 dimensional case
                        # Calculate the probability of proposing i from j
                        # NEEDS UPDATED!!!
                        #Chain.ProposalDistribution[j].ProposalProbability[i] = logpdf(Normal(ProposalMean[1], ProposalCovariance[1]), Chain.Geometry[i].Parameters[1])
                    else
                        # Calculate the probability of proposal from j to i
                        Chain.Geometry[j].ProposalProbability[i] = logpdf(Chain.ProposalDistribution[j].Density, Chain.Geometry[i].Parameters)
                    end
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
                if Model.NumOfParas == 1
                    # 1 dimensional case
                    Chain.Geometry[PropNum].Parameters[1] = rand(Normal(ProposalMean[1], ProposalCovariance[1])) # NEEDS UPDATED!!!
                else
                    # Multi-dimensional case
                    Chain.Geometry[PropNum].Parameters    = rand(Chain.ProposalDistribution[Chain.SampleIndicator].Density)
                end

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
                    if Model.NumOfParas == 1
                        # 1 dimensional case
                        # Calculate the probability of proposing i from j
                        # NEEDS UPDATED!!!
                        Chain.Geometry[j].ProposalProbability[i] = logpdf(Normal(ProposalMean[1], ProposalCovariance[1]), Chain.Geometry[i].Parameters[1])
                    else
                        # Calculate the probability of proposal from j to i
                        Chain.Geometry[j].ProposalProbability[i] = logpdf(Chain.ProposalDistribution[j].Density, Chain.Geometry[i].Parameters)
                    end
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
                if Model.NumOfParas == 1
                    # 1 dimensional case
                    Chain.Geometry[PropNum].Parameters[1] = rand(Normal(ProposalMean[1], ProposalCovariance[1])) # NEEDS UPDATED!!!
                else
                    # Multi-dimensional case
                    Chain.Geometry[PropNum].Parameters    = rand(Chain.ProposalDistribution[Chain.SampleIndicator].Density)
                end

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
                    if Model.NumOfParas == 1
                        # 1 dimensional case
                        # Calculate the probability of proposing i from j
                        # NEEDS UPDATED!!!
                        Chain.Geometry[j].ProposalProbability[i] = logpdf(Normal(ProposalMean[1], ProposalCovariance[1]), Chain.Geometry[i].Parameters[1])
                    else
                        # Calculate the probability of proposal from j to i
                        Chain.Geometry[j].ProposalProbability[i] = logpdf(Chain.ProposalDistribution[j].Density, Chain.Geometry[i].Parameters)
                    end
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
                if Model.NumOfParas == 1
                    # 1 dimensional case
                    Chain.Geometry[PropNum].Parameters[1] = Chain.Geometry[Chain.SampleIndicator].Parameters[1] + rand(Chain.ProposalDistribution.Density)
                else
                    # Multi-dimensional case
                    Chain.Geometry[PropNum].Parameters    = Chain.Geometry[Chain.SampleIndicator].Parameters + rand(Chain.ProposalDistribution.Density) # using same proposal
                end

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





function UpdateODEGeometryRandom!(Model::ODEModel, Chain::MarkovChain, PropNum::Int64, CalculateGradLL::Bool)

for k = 1:Model.NumOfParas
    Chain.Geometry[PropNum].LogPrior = logpdf(Model.Priors[k], Chain.Geometry[PropNum].Parameters[k])

    # If the proposed parameter has zero prior probability mass then exit
    if isinf(Chain.Geometry[PropNum].LogPrior)
        println("Proposed parameter outside prior")
        Chain.Geometry[PropNum].LL = -1e200
        return false
    end
end



# Calculate the log-likelihood for the ODE model
#try
    Sol = Sundials.cvode((t,y,ydot)->Model.ODEFunction(t,y,ydot,Chain.Geometry[PropNum].Parameters), Model.DefaultInitialConditions, Model.DataTimePoints[:,1], reltol=Model.reltol, abstol=Model.abstol)
    #Sol = ode23((t,x)->Model.ODEFunction(t,x,zeros(Model.NumOfSpecies),Chain.Geometry[PropNum].Parameters), Model.DefaultInitialConditions, Model.DataTimePoints[:,1])
    #println(typeof(Sol))
  #catch
#    Chain.Geometry[PropNum].LL = -1e200
#    return false
#end

Chain.Geometry[PropNum].LL = 0

for s in Model.ObservedSpecies
    for t = 1:Model.NumOfTimePoints
        Chain.Geometry[PropNum].LL = Chain.Geometry[PropNum].LL + logpdf(Normal(Model.DataMeasurements[t,s], sqrt(Model.DataNoiseVariance[t,s])), Sol[t,s])
    end
end


# Only if flag is "true" for calculating gradients
if CalculateGradLL
    # Calculate the sensitivities of the ODE model wrt each parameter
    Sol_fd  = Array(Float64, Model.NumOfTimePoints, Model.NumOfSpecies, Model.NumOfParas)
    Epsilon = Model.abstol*100

    # Calculate the gradient of the log-likeilhood
    Chain.Geometry[PropNum].GradLL = zeros(Model.NumOfParas,1)

    for i = 1:Model.NumOfParas
        NewParas     = copy(Chain.Geometry[PropNum].Parameters)
        NewParas[i]  = NewParas[i] + Epsilon
        True_Epsilon = NewParas[i] - Chain.Geometry[PropNum].Parameters[i]

        Sol_fd[:,:,i] = Sundials.cvode((t,y,ydot)->Model.ODEFunction(t,y,ydot,NewParas), Model.DefaultInitialConditions, Model.DataTimePoints[:,1], reltol=Model.reltol, abstol=Model.abstol)
        #Sol_fd[:,:,i] = ode23((t,x)->Model.ODEFunction(t,x,zeros(Model.NumOfSpecies),NewParas), Model.DefaultInitialConditions, Model.DataTimePoints[:,1])
        Sol_fd[:,:,i] = (Sol_fd[:,:,i] - Sol)./True_Epsilon

        #for s in Model.ObservedSpecies
        #    for t = 1:Model.NumOfTimePoints
        #        Chain.Geometry[PropNum].GradLL[i] = Chain.Geometry[PropNum].GradLL[i] -
        #                                            ((Sol[t,s]-Model.DataMeasurements[t,s])/Model.DataNoiseVariance[t,s])*Sol_fd[t,s,i]
        #    end
        #end

        for s in Model.ObservedSpecies
            Temp = ((Sol[:,s]-Model.DataMeasurements[:,s])./Model.DataNoiseVariance[:,s])'*Sol_fd[:,s,i]
            Chain.Geometry[PropNum].GradLL[i] = Chain.Geometry[PropNum].GradLL[i] - Temp[1];
        end

        # Add gradient of prior
        GradLogPrior = (logpdf(Model.Priors[i], Chain.Geometry[PropNum].Parameters[i]+Epsilon) - logpdf(Model.Priors[i], Chain.Geometry[PropNum].Parameters[i]) )./Epsilon
        Chain.Geometry[PropNum].GradLL[i] = Chain.Geometry[PropNum].GradLL[i] + GradLogPrior

        # Truncate the log gradient for better sampling behaviour?
        #Chain.Geometry[PropNum].GradLL = sign(Chain.Geometry[PropNum].GradLL).*min(abs(Chain.Geometry[PropNum].GradLL), 1000)
    end


    # Now sample the some pseudo data and calculate LL
    Chain.ProposalDistribution[PropNum].TangentVectors = zeros(Model.NumOfParas, Chain.ProposalDistribution[PropNum].NumberOfVectors);
    Pseudodata  = zeros(Model.NumOfTimePoints, length(Model.ObservedSpecies));

    # For each auxiliary tangent vector
    for TangNum = 1:Chain.ProposalDistribution[PropNum].NumberOfVectors

        # Sample pseudodata and calculate log-likelihood of sampling it
        for s in Model.ObservedSpecies
            for t = 1:Model.NumOfTimePoints
                Pseudodata[t,s] = rand( Normal(Sol[t,s], sqrt(Model.DataNoiseVariance[t,s])) )
            end
        end

        # Now calculate the actual tangent vectors, i.e. derivative of LL using pseudodata
        Temp = 0.0;

        for i = 1:Model.NumOfParas
            for s in Model.ObservedSpecies
                Temp = ((Sol[:,s]-Pseudodata[:,s])./Model.DataNoiseVariance[:,s])'*Sol_fd[:,s,i]
                Chain.ProposalDistribution[PropNum].TangentVectors[i,TangNum] = Chain.ProposalDistribution[PropNum].TangentVectors[i,TangNum] - Temp[1]
            end
        end

    end


    # Create approximate metric tensor
    Chain.Geometry[PropNum].HessianLL = (1/Chain.ProposalDistribution[PropNum].NumberOfVectors)*(Chain.ProposalDistribution[PropNum].TangentVectors*Chain.ProposalDistribution[PropNum].TangentVectors');

    ### FOR DEBUG ONLY
    # Only if flag is "true" for calculating metric tensors
    if 1 == 0
        # Calculate the metric tensor
        Chain.Geometry[PropNum].HessianLL = zeros(Model.NumOfParas, Model.NumOfParas)

        for i = 1:Model.NumOfParas
            for j = i:Model.NumOfParas

                for SpeciesNum in Model.ObservedSpecies
                    for t = 1:Model.NumOfTimePoints
                        Chain.Geometry[PropNum].HessianLL[i,j] = Chain.Geometry[PropNum].HessianLL[i,j] +
                                                                 (Sol_fd[t,SpeciesNum,i]/Model.DataNoiseVariance[t,SpeciesNum]*Sol_fd[t,SpeciesNum,j])
                    end
                end

                Chain.Geometry[PropNum].HessianLL[j,i] = Chain.Geometry[PropNum].HessianLL[i,j];
            end
        end

        #println(PropNum)
        #println("Likelihood...")
        #println(Chain.Geometry[PropNum].LL)
        println("Hessian...")
        println(Chain.Geometry[PropNum].HessianLL)
    end
    ##############

end # end calculate GradientLL

end







function UpdateODEGeometry!(Model::ODEModel, Chain::MarkovChain, PropNum::Int64, CalculateGradLL::Bool, CalculateHessianLL::Bool)


for k = 1:Model.NumOfParas
    Chain.Geometry[PropNum].LogPrior = logpdf(Model.Priors[k], Chain.Geometry[PropNum].Parameters[k])

    # If the proposed parameter has zero prior probability mass then exit
    if isinf(Chain.Geometry[PropNum].LogPrior)
        #println("Proposed parameter outside prior")
        Chain.Geometry[PropNum].LL = -1e200
        return false
    end
end



# Calculate the log-likelihood for the ODE model
#try
    Sol = Sundials.cvode((t,y,ydot)->Model.ODEFunction(t,y,ydot,Chain.Geometry[PropNum].Parameters), Model.DefaultInitialConditions, Model.DataTimePoints[:,1], reltol=Model.reltol, abstol=Model.abstol)
    #println(length(Model.DataTimePoints[:,1]))
    #Sol = ode23((t,x)->Model.ODEFunction(t,x,zeros(Model.NumOfSpecies),Chain.Geometry[PropNum].Parameters), Model.DataTimePoints[:,1], Model.DefaultInitialConditions)
    #println(length(Sol[2][:]))
    #println((Sol[2][:]))
  #catch
#    Chain.Geometry[PropNum].LL = -1e200
#    return false
#end

Chain.Geometry[PropNum].LL = 0

for s in Model.ObservedSpecies
    for t = 1:Model.NumOfTimePoints
        Chain.Geometry[PropNum].LL = Chain.Geometry[PropNum].LL + logpdf(Normal(Model.DataMeasurements[t,s], sqrt(Model.DataNoiseVariance[t,s])), Sol[t,s])
        #Chain.Geometry[PropNum].LL = Chain.Geometry[PropNum].LL + logpdf(Normal(Model.DataMeasurements[t,s], sqrt(Model.DataNoiseVariance[t,s])), Sol[2][t][s])
    end
end


# Only if flag is "true" for calculating gradients
if CalculateGradLL
    # Calculate the sensitivities of the ODE model wrt each parameter
    Sol_fd  = Array(Float64, Model.NumOfTimePoints, Model.NumOfSpecies, Model.NumOfParas)
    Epsilon = Model.abstol*100

    # Calculate the gradient of the log-likeilhood
    Chain.Geometry[PropNum].GradLL = zeros(Model.NumOfParas,1)

    if CalculateHessianLL
        # Setup variable for the metric tensor
        Chain.Geometry[PropNum].HessianLL = zeros(Model.NumOfParas, Model.NumOfParas)
    end

    for i = 1:Model.NumOfParas
        NewParas     = copy(Chain.Geometry[PropNum].Parameters)
        NewParas[i]  = NewParas[i] + Epsilon
        True_Epsilon = NewParas[i] - Chain.Geometry[PropNum].Parameters[i]

        Sol_fd[:,:,i] = Sundials.cvode((t,y,ydot)->Model.ODEFunction(t,y,ydot,NewParas), Model.DefaultInitialConditions, Model.DataTimePoints[:,1], reltol=Model.reltol, abstol=Model.abstol)
        Sol_fd[:,:,i] = (Sol_fd[:,:,i] - Sol)./True_Epsilon
        #Sol_fd = ode23((t,x)->Model.ODEFunction(t,x,zeros(Model.NumOfSpecies),NewParas), Model.DefaultInitialConditions, Model.DataTimePoints[:,1])
        #Sol_fd[2] = (Sol_fd[2] - Sol[2])./True_Epsilon

        for s in Model.ObservedSpecies
            for t = 1:Model.NumOfTimePoints
                Chain.Geometry[PropNum].GradLL[i] = Chain.Geometry[PropNum].GradLL[i] -
                                                    ((Sol[t,s]-Model.DataMeasurements[t,s])/Model.DataNoiseVariance[t,s])*Sol_fd[t,s,i]
                #Chain.Geometry[PropNum].GradLL[i] = Chain.Geometry[PropNum].GradLL[i] -
                #                                    ((Sol[2][t][s]-Model.DataMeasurements[t,s])/Model.DataNoiseVariance[t,s])*Sol_fd[2][t][s]
            end
        end

        # Add gradient of prior
        GradLogPrior = (logpdf(Model.Priors[i], Chain.Geometry[PropNum].Parameters[i]+Epsilon) - logpdf(Model.Priors[i], Chain.Geometry[PropNum].Parameters[i]) )./Epsilon
        Chain.Geometry[PropNum].GradLL[i] = Chain.Geometry[PropNum].GradLL[i] + GradLogPrior

        # Truncate the log gradient for better sampling behaviour?
        #Chain.Geometry[PropNum].GradLL = sign(Chain.Geometry[PropNum].GradLL).*min(abs(Chain.Geometry[PropNum].GradLL), 1000)

        # Only if flag is "true" for calculating metric tensors
        if CalculateHessianLL
            # Calculate the metric tensor
            Chain.Geometry[PropNum].HessianLL = zeros(Model.NumOfParas, Model.NumOfParas)

            for i = 1:Model.NumOfParas
                for j = i:Model.NumOfParas

                    for SpeciesNum in Model.ObservedSpecies
                        for t = 1:Model.NumOfTimePoints
                            Chain.Geometry[PropNum].HessianLL[i,j] = Chain.Geometry[PropNum].HessianLL[i,j] +
                                                                     (Sol_fd[t,SpeciesNum,i]/Model.DataNoiseVariance[t,SpeciesNum]*Sol_fd[t,SpeciesNum,j])
                            #Chain.Geometry[PropNum].HessianLL[i,j] = Chain.Geometry[PropNum].HessianLL[i,j] +
                            #                                         (Sol_fd[2][t][SpeciesNum]/Model.DataNoiseVariance[t,SpeciesNum]*Sol_fd[2][t][SpeciesNum])
                        end
                    end

                    Chain.Geometry[PropNum].HessianLL[j,i] = Chain.Geometry[PropNum].HessianLL[i,j];
                end
            end

        end # end calculate HessianLL


    end



end # end calculate GradientLL

end

