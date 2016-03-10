# tunelogistic uses the logistic function to scale the acceptance rate by a factor ranging from 0 to 2
# k gives the curve's steepness. For larger k, the curve becomes more steep
tunelogistic(x::Real,k::Real) = 2*one(x)/(one(x)+exp(-k*x))

# tuneerf uses the error function (erf) to scale the acceptance rate by a factor ranging from 0 to 2
# k gives the curve's steepness. For larger k, the curve becomes more steep
tuneerf(x::Real,k::Real) = erf(k*x) + one(x)

