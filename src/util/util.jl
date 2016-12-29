###Function to create a string that shows the size of the object in memory
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

function showsamplerstatevars(toshow::Dict,lstr=nothing)
    println("Sampler state variables",lstr==nothing?":":string(" for ",lstr,":"))
    for k in keys(toshow)
        print(" [\"",k,"\"]: ")
        show(toshow[k])
    end
    println()
end

function showsamplerstatevars(toshow::AbstractArray,lstr=nothing)
    println("Array of sampler state variables",lstr==nothing?":":string(" for ",lstr,":"))
    for i=1:length(toshow)
        showsamplerstatevars(toshow[i],string("[",i,"]"))
    end
    println()
end

showsamplerstatevars(toshow,lstr=nothing) = showsamplerstatevars(getsamplerstatevars(toshow),lstr)
