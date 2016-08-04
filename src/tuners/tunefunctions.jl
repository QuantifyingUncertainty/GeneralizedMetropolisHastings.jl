# tunelogistic uses the logistic function to scale the acceptance rate by a factor ranging from 0 to 2
# k gives the curve's steepness. For larger k, the curve becomes more steep
function tunelogistic(x::Real,k::Real)
    r = zero(x)
    if x < 0
        r = one(x)/(one(x)+exp(-k*x)) + one(x)/2
    else
        r = 2one(x)/(one(x)+exp(-k*x))
    end
    sqrt(r)
end

# tuneerf uses the error function (erf) to scale the acceptance rate by a factor ranging from 0 to 2
# k gives the curve's steepness. For larger k, the curve becomes more steep
function tuneerf(x::Real,k::Real)
    r = zero(x)
    if x < 0
        r = erf(k*x)/2 + one(x)
    else
        r = erf(k*x) + one(x)
    end
    sqrt(r)
end
