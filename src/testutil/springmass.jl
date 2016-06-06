###ODE for spring-mass dynamic system
function springmassode(t,y,ydot,paras)
    ydot[1] = y[2]
    ydot[2] = -paras[1]/paras[2]*y[1] #-K/M*X
end

###Function to generate data for the spring-mass dynamic system
function springmassdata(t::AbstractVector,y0::AbstractVector,paras::AbstractVector)
    a = sqrt(paras[1]/paras[2])
    y = [y0[1]*cos(a*t)+y0[2]/a*sin(a*t) -a*y0[1]*sin(a*t)+y0[2]*cos(a*t)]
end

springmassnoisy(t::AbstractVector,y0::AbstractVector,pv::AbstractVector,n::AbstractNoiseModel) = applynoise!(n,springmassdata(t,y0,pv))

###Helper function to create a spring-mass MCModel
### time to evaluate the ODE
### initial of the dynamic system
### modelvals the values for the model parameters
### variance the variance of measurement noise
### paraminit the specifiers for the MCMC parameters
function _model(::Type{Val{:springmass}},t::AbstractVector,y0::AbstractVector,pv::AbstractVector,variance::AbstractVector,paraminit...)
    p = parameters([:K,:M],paraminit...)
    n = noise(:gaussian,variance)
    d = data(:function,t,springmassnoisy,t,y0,pv,n)
    _model(Val{:ode},p,d,n,springmassode,initial,2,[1,2];name = "Spring-Mass ODE")
end
