### Load the plot packages
using Gadfly
using PyPlot

### Include the data file
include("../data/fitzhughnagumo.jl")

### Print a help message
println("Plot the FitzHugh Nagumo data using:")
println(" Gadfly package: gadflyplot_fhn()")
println(" PyPlot package: pyplot_fhn()")

### Plot the FitzHugh Nagumo data using the Gadfly package
function gadflyplot_fhn()

  #load the data
  timepoints,measurements = fitzhugh_nagumo_data()

  #plot the data
  Gadfly.plot(x=timepoints,y=measurements)

end

### Plot the FitzHugh Nagumo data using the PyPlot package
function pyplot_fhn()

  #load the data
  timepoints,measurements = fitzhugh_nagumo_data()

  #plot the data
  fig = PyPlot.figure()
  PyPlot.plot(timepoints,measurements[:,1];label="Membrane Potential")
  PyPlot.plot(timepoints,measurements[:,2];label="Recovery Variable")
  PyPlot.xlabel("Time")
  PyPlot.ylabel("Voltage")
  PyPlot.title("FitzHugh-Nagumo Measurement Data")
  PyPlot.grid("on")
  PyPlot.legend(loc="upper right",fancybox="true")

  #return the figure object
  fig
end

# File return statement
nothing
