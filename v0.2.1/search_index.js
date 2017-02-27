var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#GeneralizedMetropolisHastings.jl-1",
    "page": "Home",
    "title": "GeneralizedMetropolisHastings.jl",
    "category": "section",
    "text": "A parallel Monte-Carlo Markov Chain algorithmCode base for Generalized Metropolis-Hastings (GMH), an intrinsically-parallel Monte-Carlo Markov Chain algorithm (Calderhead, 2014)."
},

{
    "location": "index.html#User-Manual-1",
    "page": "Home",
    "title": "User Manual",
    "category": "section",
    "text": "Follow the steps in the Package Guide to:Install on various platforms\nRun existing MCMC experiments\nWrite new MCMC experiments in existing problem domains\nDevelop new MCMC experiments in your own problem domain\nExtend the GMH package"
},

{
    "location": "index.html#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "Calderhead B. (2014), A general construction for parallelizing Metropolis-Hastings algorithms, PNAS, Vol: 111, Pages: 17408-17413 doi: 10.1073/pnas.1408184111"
},

{
    "location": "man/guide.html#",
    "page": "Guide",
    "title": "Guide",
    "category": "page",
    "text": "CurrentModule = GeneralizedMetropolisHastings"
},

{
    "location": "man/guide.html#Package-Guide-1",
    "page": "Guide",
    "title": "Package Guide",
    "category": "section",
    "text": ""
},

{
    "location": "man/guide.html#Installation-Guidelines-1",
    "page": "Guide",
    "title": "Installation Guidelines",
    "category": "section",
    "text": ""
},

{
    "location": "man/guide.html#Installing-the-GMH-package-to-replicate-published-experiments-1",
    "page": "Guide",
    "title": "Installing the GMH package to replicate published experiments",
    "category": "section",
    "text": "To run the GMH package with published experiments and data sets on Amazon Web Services or JuliaBox, see the detailed installation instructions"
},

{
    "location": "man/guide.html#Installing-the-GMH-package-locally-1",
    "page": "Guide",
    "title": "Installing the GMH package locally",
    "category": "section",
    "text": "The package is registered in METADATA.jl and can be installed using Pkg.add.Pkg.add(\"GeneralizedMetropolisHastings\")You can test if the package is working correctly by runningPkg.test(\"GeneralizedMetropolisHastings\")"
},

{
    "location": "man/guide.html#Running-existing-experiments-1",
    "page": "Guide",
    "title": "Running existing experiments",
    "category": "section",
    "text": "A set of ready-to-run, elementary examples is available in GMHExamples.jl.Information on how to replicate published experiments in specific problem domains is available from Quantifying Uncertainty."
},

{
    "location": "man/guide.html#Writing-new-experiments-in-existing-problem-domains-1",
    "page": "Guide",
    "title": "Writing new experiments in existing problem domains",
    "category": "section",
    "text": "TBD"
},

{
    "location": "man/guide.html#Developing-new-problem-domains-1",
    "page": "Guide",
    "title": "Developing new problem domains",
    "category": "section",
    "text": "The following steps need to be performed to specify a new problem domain.Set up a folder structure\nAdd measurement data\nSpecify measurement noise\nSpecify model functionsIf the above steps have been completed, runnable MCMCM experiments can be defined in julia scripts or IJulia notebooks as outlined in Writing new experiments in existing problem domains."
},

{
    "location": "man/folderstructure.html#",
    "page": "Set up a folder structure",
    "title": "Set up a folder structure",
    "category": "page",
    "text": ""
},

{
    "location": "man/folderstructure.html#Set-up-a-folder-structure-1",
    "page": "Set up a folder structure",
    "title": "Set up a folder structure",
    "category": "section",
    "text": ""
},

{
    "location": "man/folderstructure.html#For-a-simple-project-1",
    "page": "Set up a folder structure",
    "title": "For a simple project",
    "category": "section",
    "text": "Prepare a separate julia package with the recommended folder structure:MyGMHPackage.jl/\n    MyGMHPackage.jl\n    data/\n    models/\n    notebooks/\n    scripts/\n    test/Measurement data sets are kept in data, model functions in model, and runnable MCMC experiments in scripts and notebooks. Keep unit and other tests for your package in test. Define other folders (e.g., for specific analysis code or to store the results of the MCMC runs) as required.MyGMHPackage.jl is the top level module file. It specifies imported packages, included module files, and package exports. An example is available here."
},

{
    "location": "man/folderstructure.html#For-a-larger-problem-domain-1",
    "page": "Set up a folder structure",
    "title": "For a larger problem domain",
    "category": "section",
    "text": "When the problem domain of interest is likely to contain more than one type of estimation experiment, repeat the above folder structure for different experiment types, and group them hierarchically into problems and domains.MyGMHPackage.jl/\n    MyGMHPackage.jl\n    domain1/\n        problem1/\n            data/\n            models/\n            notebooks/\n            scripts/\n        problem2/\n            ...\n        problem3/\n            ...\n    domain2/\n        problem1/\n            ...\n        problem2/\n            ...\n    ...\n    test/Keep the MyGMHPackage.jl module file and test folder at the top level of the package.An example of this folder structure can be seen in GMHExamples.jl."
},

{
    "location": "man/folderstructure.html#Add-to-Julia-search-path-1",
    "page": "Set up a folder structure",
    "title": "Add to Julia search path",
    "category": "section",
    "text": "Add the package to the Julia search path by adding the following line to your .juliarc.jl file.push!(LOAD_PATH,\"/local/path/to/MyGMHPackage.jl\")"
},

{
    "location": "man/data.html#",
    "page": "Add measurement data",
    "title": "Add measurement data",
    "category": "page",
    "text": ""
},

{
    "location": "man/data.html#Add-measurement-data-1",
    "page": "Add measurement data",
    "title": "Add measurement data",
    "category": "section",
    "text": "CurrentModule = GeneralizedMetropolisHastings"
},

{
    "location": "man/data.html#Creating-measurement-data-1",
    "page": "Add measurement data",
    "title": "Creating measurement data",
    "category": "section",
    "text": "The GMH package requires measurement data to be wrapped in an AbstractData object.Two AbstractData subtypes are available. These are likely to cover most if not all requirements for specifying measurement data. The first is useful for predefined datasets, the second if data is generated by a function. This may be the model function itself or a different data generating function.Create objects of these subtypes with the data factory function."
},

{
    "location": "man/data.html#Predefined-datasets-1",
    "page": "Add measurement data",
    "title": "Predefined datasets",
    "category": "section",
    "text": "Create a data object wrapping an array with the data(:array,...) method. The wrapped array values should have 1 row per data point and 1 column per variable. Vector index contains an index value for each data point.using GeneralizedMetropolisHastings\nindex = 0:10;\nvalues = 0:0.5:5.0;\nd = data(:array,index,values)"
},

{
    "location": "man/data.html#Data-generating-functions-1",
    "page": "Add measurement data",
    "title": "Data generating functions",
    "category": "section",
    "text": "Create a data object wrapping a function with data(:function!,...) or data(:function,...).The first variant requires a function which creates data values in a preallocated vector out, taking it as its first argument. In the example below, t is a time index for each data point, and the last argument of the data function is the x argument of data function sin!.using GeneralizedMetropolisHastings\nsin!(r::Vector,x::AbstractVector) = (for i=1:length(x) r[i] = sin(x[i]) end ; r)\nt = 0.0:0.1:1.0;\nout = zeros(t);\nd = data(:function!,t,out,sin!,2π*t)The second variant requires a data generating function returning a array with 1 row per data point and 1 column per variable. t is a time index for each data point. In the example below, the final 3 arguments of the data function are the arguments of the sincos function.using GeneralizedMetropolisHastings\nsincos(x::AbstractVector,a::Real,b::Real) =  cat(2,sin(a*x),cos(b*x));\nt = 0.0:0.1:1.0;\nd = data(:function,t,sincos,2π*t,2.0,3.0)"
},

{
    "location": "man/data.html#Storing-data-inside-the-package-1",
    "page": "Add measurement data",
    "title": "Storing data inside the package",
    "category": "section",
    "text": "Keep predefined data sets in the data/ subfolder, as explained in Set up a folder structure."
},

{
    "location": "man/data.html#Small-data-sets-1",
    "page": "Add measurement data",
    "title": "Small data sets",
    "category": "section",
    "text": "Small data sets can be stored as an array in a .jl file, wrapped in a function returning an (index,values) tuple.function mysmalldataset()\n    index = 0:10;\n    values = cat(2,0:0.5:5.0,10.0:1.0:20.0); #an 11x2 data array\n    return index,values\nendThe data wrapping function can then be used to create an :array data object:julia> include(\"mysmalldataset.jl\")\njulia> d = data(:array,mysmalldataset()...) #... extracts the returned (index,value) Tuple"
},

{
    "location": "man/data.html#Large-data-sets-1",
    "page": "Add measurement data",
    "title": "Large data sets",
    "category": "section",
    "text": "For larger data sets, it is recommended to store them in data format files such as .jld or .mat.For more information see the julia packages JLD and MAT.Write data functions that load and extract the data from the data files. Some examples are provided here."
},

{
    "location": "man/data.html#Using-AbstractData-objects-1",
    "page": "Add measurement data",
    "title": "Using AbstractData objects",
    "category": "section",
    "text": "Retrieve stored data values with datavalues\nRetrieve the index values with dataindex\nGet the number of variables and values with numvars and numvalues\nGenerate a new set of values with internal function generate!using PyPlot\nusing GeneralizedMetropolisHastings\n\n#define the data object\nt = -1.0:0.05:1.0;\nd = data(:function,t,sin,2π*t);\n\n#generate a new set of data values\nGeneralizedMetropolisHastings.generate!(d);\n\n#retrieve data index and values and plot\nplot(dataindex(d),datavalues(d))\ntitle(\"Data with $(numvars(d)) variable and $(numvalues(d)) values\")\nxlabel(\"Time\")\nylabel(\"Values\")\nxlim(-1.0,1.0)\n\nsavefig(\"dataplot.svg\")(Image: )"
},

{
    "location": "man/data.html#Extend-the-AbstractData-subtype-1",
    "page": "Add measurement data",
    "title": "Extend the AbstractData subtype",
    "category": "section",
    "text": "It is possible to define your own data type by extending AbstractData.Import AbstractData and data into your package module file\nDefine a subtype of AbstractData\nPick a unique symbol (e.g. :mydata) for your datatype and define a factory method dataimmutable MyDataType <: AbstractData\n    field1\n    field2\nend\n\nfunction data(::Type{Val{:mydata}},arg1,arg2,arg3,arg4)\n    #steps to defined field1 and field2 out of arg1 to arg4\n    ...\n    MyDataType(field1,field2) #call the constructor of MyDataType\nendWith that in place, create objects of MyDataType as follows:d = data(:mydata,arg1,arg2,arg3,arg4)Import and implement methods for all functions listed in Data - Public Functions and Data - Internal Functions"
},

{
    "location": "man/noise.html#",
    "page": "Specify measurement noise",
    "title": "Specify measurement noise",
    "category": "page",
    "text": "CurrentModule = GeneralizedMetropolisHastings"
},

{
    "location": "man/noise.html#Specify-measurement-noise-1",
    "page": "Specify measurement noise",
    "title": "Specify measurement noise",
    "category": "section",
    "text": "TBD"
},

{
    "location": "man/models.html#",
    "page": "Specify model functions",
    "title": "Specify model functions",
    "category": "page",
    "text": "CurrentModule = GeneralizedMetropolisHastings"
},

{
    "location": "man/models.html#Specify-model-functions-1",
    "page": "Specify model functions",
    "title": "Specify model functions",
    "category": "section",
    "text": "TBD"
},

{
    "location": "lib/public.html#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "lib/public.html#Public-Documentation-1",
    "page": "Public",
    "title": "Public Documentation",
    "category": "section",
    "text": "CurrentModule = GeneralizedMetropolisHastingsDocumentation for GeneralizedMetropolisHastings.jl's public interface.See Internal Documentation for internal package docs."
},

{
    "location": "lib/public.html#Contents-1",
    "page": "Public",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"public.md\"]\nDepth = 3"
},

{
    "location": "lib/public.html#Index-1",
    "page": "Public",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "lib/public.html#Public-Interface-1",
    "page": "Public",
    "title": "Public Interface",
    "category": "section",
    "text": ""
},

{
    "location": "lib/public.html#Data-Interface-1",
    "page": "Public",
    "title": "Data - Interface",
    "category": "section",
    "text": ""
},

{
    "location": "lib/public.html#GeneralizedMetropolisHastings.data",
    "page": "Public",
    "title": "GeneralizedMetropolisHastings.data",
    "category": "Function",
    "text": "data(s::Symbol,args...)\n\nCreate data objects.\n\nCurrently implemented variants:\n\n\n\ndata(:array,index::AbstractVector,values::AbstractArray)\n\nCreate a data object wrapping an array of fixed, predefined data values.\n\nArguments\n\nindex::AbstractVector: an index for each row of the values array (e.g., a time index)\nvalues::AbstractArray: each row is a data point and each column a variable\n\n\n\n\n\ndata(:function!,index::AbstractVector,out::AbstractArray,f!::Function,[,...args])\n\nCreate a DataFunction object which encapsulates an in-place data generating function.\n\nThe :function! variant expects a pre-allocated output array out, and a data generating function which takes out as its first argument.\n\nArguments\n\nindex::AbstractVector: an index for the rows of the result (e.g., a time index)\nout::AbstractArray: array with one row per data point and one column per variable\nf!::Function: a data generating function taking as first argument out\nargs...: additional arguments of the data generating function\n\n\n\n\n\ndata(:function,index::AbstractVector,f::Function[,args...])\n\nCreate a DataFunction object which encapsulates a data generating function.\n\nThe :function variant can be used if the data generating function cannot generate its result in a pre-allocated output array. The result should contain one row per data point and one column per variable.\n\nArguments\n\nindex::AbstractVector: an index for the rows of the result (e.g., a time index)\nf::Function: a function returning an array of values as a result\nargs...: additional arguments of the data generating function\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#GeneralizedMetropolisHastings.dataindex-Tuple{GeneralizedMetropolisHastings.AbstractData}",
    "page": "Public",
    "title": "GeneralizedMetropolisHastings.dataindex",
    "category": "Method",
    "text": "dataindex(d)\n\nRetrieve the data index.\n\n\n\n"
},

{
    "location": "lib/public.html#GeneralizedMetropolisHastings.datavalues-Tuple{GeneralizedMetropolisHastings.AbstractData}",
    "page": "Public",
    "title": "GeneralizedMetropolisHastings.datavalues",
    "category": "Method",
    "text": "datavalues(d)\n\nRetrieve the data values.\n\n\n\n"
},

{
    "location": "lib/public.html#GeneralizedMetropolisHastings.numvalues-Tuple{GeneralizedMetropolisHastings.AbstractData}",
    "page": "Public",
    "title": "GeneralizedMetropolisHastings.numvalues",
    "category": "Method",
    "text": "numvalues(d)\n\nReturn the number of data values.\n\n\n\n"
},

{
    "location": "lib/public.html#GeneralizedMetropolisHastings.numvars-Tuple{GeneralizedMetropolisHastings.AbstractData}",
    "page": "Public",
    "title": "GeneralizedMetropolisHastings.numvars",
    "category": "Method",
    "text": "numvars(d)\n\nReturn the number of variables.\n\n\n\n"
},

{
    "location": "lib/public.html#Data-Public-Functions-1",
    "page": "Public",
    "title": "Data - Public Functions",
    "category": "section",
    "text": "data\ndataindex(::AbstractData)\ndatavalues(::AbstractData)\nnumvalues(::AbstractData)\nnumvars(::AbstractData)"
},

{
    "location": "lib/public.html#GeneralizedMetropolisHastings.AbstractData",
    "page": "Public",
    "title": "GeneralizedMetropolisHastings.AbstractData",
    "category": "Type",
    "text": "Base type for objects wrapping measurement data.\n\n\n\n"
},

{
    "location": "lib/public.html#Data-Public-Types-1",
    "page": "Public",
    "title": "Data - Public Types",
    "category": "section",
    "text": "AbstractData"
},

{
    "location": "lib/internal.html#",
    "page": "Internal",
    "title": "Internal",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internal.html#Internal-Documentation-1",
    "page": "Internal",
    "title": "Internal Documentation",
    "category": "section",
    "text": "CurrentModule = GeneralizedMetropolisHastingsDocumentation for GeneralizedMetropolisHastings.jl's internal functions.You may need to extend some of these types or functions when developing new problem domains or extending the GMH package.See Public Documentation for its public interface."
},

{
    "location": "lib/internal.html#Contents-1",
    "page": "Internal",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"internal.md\"]\nDepth = 3"
},

{
    "location": "lib/internal.html#Index-1",
    "page": "Internal",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"internal.md\"]"
},

{
    "location": "lib/internal.html#Internal-functions-and-types-1",
    "page": "Internal",
    "title": "Internal functions and types",
    "category": "section",
    "text": ""
},

{
    "location": "lib/internal.html#Data-Internal-1",
    "page": "Internal",
    "title": "Data - Internal",
    "category": "section",
    "text": ""
},

{
    "location": "lib/internal.html#GeneralizedMetropolisHastings.datatypename-Tuple{GeneralizedMetropolisHastings.AbstractData}",
    "page": "Internal",
    "title": "GeneralizedMetropolisHastings.datatypename",
    "category": "Method",
    "text": "datatypename(d)\n\nReturn the name of the data type (used by show()).\n\n\n\n"
},

{
    "location": "lib/internal.html#Base.eltype-Tuple{GeneralizedMetropolisHastings.AbstractData}",
    "page": "Internal",
    "title": "Base.eltype",
    "category": "Method",
    "text": "eltype(d)\n\nReturn the type of the data values.\n\n\n\n"
},

{
    "location": "lib/internal.html#GeneralizedMetropolisHastings.generate!-Tuple{GeneralizedMetropolisHastings.AbstractData}",
    "page": "Internal",
    "title": "GeneralizedMetropolisHastings.generate!",
    "category": "Method",
    "text": "generate!(d)\n\nGenerate a set of data values. The function should return its argument d.\n\n\n\n"
},

{
    "location": "lib/internal.html#Data-Internal-Functions-1",
    "page": "Internal",
    "title": "Data - Internal Functions",
    "category": "section",
    "text": "datatypename(::AbstractData)\neltype(::AbstractData)\ngenerate!(::AbstractData)"
},

{
    "location": "lib/internal.html#GeneralizedMetropolisHastings.DataArray",
    "page": "Internal",
    "title": "GeneralizedMetropolisHastings.DataArray",
    "category": "Type",
    "text": "DataArray\n\nA type wrapping a predefined array of measurement data.\n\n\n\n"
},

{
    "location": "lib/internal.html#GeneralizedMetropolisHastings.DataFunction",
    "page": "Internal",
    "title": "GeneralizedMetropolisHastings.DataFunction",
    "category": "Type",
    "text": "DataFunction\n\nA type wrapping a data generating function.\n\n\n\n"
},

{
    "location": "lib/internal.html#Data-Internal-Types-1",
    "page": "Internal",
    "title": "Data - Internal Types",
    "category": "section",
    "text": "DataArray\nDataFunction"
},

]}
