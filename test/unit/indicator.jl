import GeneralizedMetropolisHastings.create_indicator_matrix
import GeneralizedMetropolisHastings.sample_indicator_matrix

#test the matrix creation functions
@test create_indicator_matrix(IndicatorMatrixStationary(),[0.1,0.4,0.5]) == [0.1 0.4 0.5;0.1 0.4 0.5;0.1 0.4 0.5]
@test create_indicator_matrix(IndicatorMatrixOptimal(),[0.1,0.4,0.5]) == eye(3)

#sample arbitrary transition matrices
@test sample_indicator_matrix([0.0 1.0 0.0;0.0 0.0 1.0;1.0 0.0 0.0],5,2) == [2,3,1,2,3]
@test sample_indicator_matrix(eye(2),3,1) == [1,1,1]
@test sample_indicator_matrix(eye(2),4,2) == [2,2,2,2]

#sample arbitrary transition matrices
#this part of the test depends on the proporties of the default random number generator in Julia
#if the test fails, it could be because the random generator has changed
srand(0)
@test rand() == 0.8236475079774124
@test sample_indicator_matrix(repmat([0.5 0.5],2),4,1) == [1,2,1,1]
@test sample_indicator_matrix(repmat([0.5 0.25 0.125 0.125],4),6,3) == [3,1,1,1,1,1]
@test sample_indicator_matrix(repmat([0.5 0.25 0.125 0.125],4),6,1) == [1,4,2,2,1,4]

