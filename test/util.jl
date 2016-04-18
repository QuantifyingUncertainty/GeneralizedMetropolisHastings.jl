function objectsizetostr(a::Any)
  bytes = Base.summarysize(a)
  result = ""
  if bytes < 10*1024
    result = @sprintf("%6d bytes  ",bytes)
  elseif bytes < 1024*1024
    result = @sprintf("%5.3f kB",bytes/(1024))
  else
    result = @sprintf("%5.3f MB",bytes/1024/1024)
  end
  result
end

##################################################

###ODE for spring-mass dynamic system
@everywhere function springmassode(t,y,ydot,paras)
    ydot[1] = y[2]
    ydot[2] = -paras[1]/paras[2]*y[1] #-K/M*X
end

###Function that can generate date for the spring-mass dynamic system
function springmassdata(t,y0,paras)
    a = sqrt(paras[1]/paras[2])
    y = [y0[1]*cos(a*t)+y0[2]/a*sin(a*t) -a*y0[1]*sin(a*t)+y0[2]*cos(a*t)]
end

springmassnoisy(n::AbstractNoiseModel,t,y0,paras) = applynoise!(n,springmassdata(t,y0,paras))

###Helper function to create a spring-mass MCModel
### time to evaluate the ODE
### initial of the dynamic system
### modelvals the values for the model parameters
### variance the variance of measurement noise
### paraminit the specifiers for the MCMC parameters
function springmassmodel(time,initial,modelvals,variance,paraminit...)
    p = parameters([:K,:M],paraminit...)
    n = noise(:gaussian,variance)
    d = data(:function,time,springmassnoisy,n,time,initial,modelvals)
    model(:ode,p,d,n,springmassode,initial,2,[1,2];name = "Spring-Mass")
end

#####################################################

###Target function
@everywhere function sincos(t,paras)
    hcat(sin(2*pi*paras[1]*t),cos(2*pi*paras[2]*t))
end

@everywhere function sincos!(r,t,paras)
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

sincosdata(t,paras,n) = data(:array,t,applynoise!(n,sincos(t,paras)))

function sincosmodel(t,modelvals,variance,paraminit...)
    p = parameters([:a,:b],paraminit...)
    n = noise(:gaussian,variance)
    d = sincosdata(t,modelvals,n)
    model(:target,p,d,n,sincos,t;name="Sine-Cosine")
end

function sincosmodel!(t,modelvals,variance,paraminit...)
    p = parameters([:a,:b],paraminit...)
    n = noise(:gaussian,variance)
    d = sincosdata(t,modelvals,n)
    model(:target!,p,d,n,sincos!,t;name="Sine-Cosine!")
end


######################################################

nothing

