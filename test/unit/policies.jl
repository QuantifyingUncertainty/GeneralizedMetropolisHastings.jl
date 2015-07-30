#test the basic constructors for generic policies
p1 = GenericPolicy(ValuesFromPrior(),IndicatorMatrixStationary(),2)
p2 = GenericPolicy(ValuesFromDefault(),IndicatorMatrixStationary(),1)

@test typeof(p1.initialize) == ValuesFromPrior && typeof(p1.propose) == ProposalFromAuxiliary && typeof(p1.indicate) == IndicatorMatrixStationary
@test typeof(p2.initialize) == ValuesFromDefault && typeof(p2.propose) == ProposalFromIndicator && typeof(p2.indicate) == IndicatorMatrixStationary

#test default constructors for most often used policies
p3 = GenericPolicy(ValuesFromDefault(),1)
p4 = GenericPolicy(2)

@test p1 == p4
@test p2 == p3

#test equality operators
@test GenericPolicy(ValuesFromPrior(),3) == GenericPolicy(2) #should be equal
@test GenericPolicy(2) != GenericPolicy(1) #the ProposalFunction will be different

println("====================")
println("Test show() function")
show(p1)
show(p2)
println("End  show() function")
println("====================")
println()

nothing
