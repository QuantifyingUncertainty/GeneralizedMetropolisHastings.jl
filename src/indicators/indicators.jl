abstract AbstractIndicatorMatrix

###Factory function. Currently implemented is :stationary
indicator(s::Symbol,args...) = _indicator(Val{s},args...)

indicatorsamples(indicator::AbstractIndicatorMatrix) = indicator.samples

function accepted(indicator::AbstractIndicatorMatrix)
    acc = 0
    @simd for i=2:length(indicator.samples)
        @inbounds indicator.samples[i]!=indicator.samples[i-1]?acc+=1:nothing
    end
    acc
end

numsamples(indicator::AbstractIndicatorMatrix) = length(indicator.samples) - 1

