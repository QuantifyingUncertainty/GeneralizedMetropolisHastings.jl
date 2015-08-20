### Load the plot packages
using Gadfly
using PyPlot

### Include the data file
include("../data/springmass.jl")

### Print a help message
println("Plot spring-mass measurement data using:")
println(" Gadfly package: gadflyplot_springmass(y0,paras,noisevar;timepoints)")
println(" PyPlot package: pyplot_springmass(y0,paras,noisevar;timepoints)")


### Plot the spring-mass data using the Gadfly package
function gadflyplot_springmass(y0,paras,noisevar;timepoints =nothing)

  #get the data
  timepoints,measurements = spring_mass_data(y0,paras,noisevar;timepoints =timepoints)

  #plot the data
  Gadfly.plot(x=timepoints,y=measurements)

end

### Plot the spring-mass data using the PyPlot package
function pyplot_springmass(y0,paras,noisevar;timepoints =nothing)

  #get the data
  timepoints,measurements = spring_mass_data(y0,paras,noisevar;timepoints =timepoints)

  #plot the data
  fig = PyPlot.figure()
  PyPlot.plot(timepoints,measurements[:,1];label="Position")
  PyPlot.plot(timepoints,measurements[:,2];label="Velocity")
  PyPlot.xlabel("Time")
  PyPlot.ylabel("Amplitude")
  PyPlot.title("Spring-Mass Measurement Data")
  PyPlot.grid("on")
  PyPlot.legend(loc="upper right",fancybox="true")

  #return the figure object
  fig
end

# File return statement
nothing
