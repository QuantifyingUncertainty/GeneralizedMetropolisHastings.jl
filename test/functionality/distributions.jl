for m in [0.8,0.95,0.98]
    b1 = distribution(:bactrian,:normal,0.0,2.0,m)
    b2 = distribution(:bactrian,:laplace,0.0,2.0,m)
    b3 = distribution(:bactrian,:triangular,0.0,2.0,m)

    x = -5.0:0.01:5.0
    PyPlot.plot(x,Distributions.pdf(b1.mixture,x))
    PyPlot.plot(x,Distributions.pdf(b2.mixture,x))
    PyPlot.plot(x,Distributions.pdf(b3.mixture,x))
    PyPlot.xlim(x[1],x[end])
end
