abstract AbstractData

data(s::Symbol,args...) = _data(Val{s},args...)
datatypename(d::AbstractData) = "AbstractData"

function show(io::IO,d::AbstractData)
    nvals = numvalues(d)
    nvars = numvars(d)
    println(io,datatypename(d)," with ",nvars," variable",nvars>1?"s":""," and ",nvals," value",nvals>1?"s":"")
    for f in fieldnames(d)
        fc = getfield(d,f)
        if typeof(fc) <: Tuple && ~isempty(fc)
            print(io," ",f,":")
            show(io,fc)
            println(io)
        else
            println(io," ",f,": ",typeof(getfield(d,f)))
        end
    end
    nothing
end
