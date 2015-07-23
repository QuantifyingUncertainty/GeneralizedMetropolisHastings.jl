#this part of the test does not depend on the seed or random generator method
@test GeneralizedMetropolisHastings.sample_indicator_matrix([0.0 1.0 0.0;0.0 0.0 1.0;1.0 0.0 0.0],4,2) == [2,3,1,2,3]
@test GeneralizedMetropolisHastings.sample_indicator_matrix(eye(2),2,1) == [1,1,1]
@test GeneralizedMetropolisHastings.sample_indicator_matrix(eye(2),3,2) == [2,2,2,2]

@test GeneralizedMetropolisHastings.create_indicator_matrix([0.2,0.3,0.5],2,GeneralizedMetropolisHastings.IndicatorMatrixSt)

#this part of the test depends on the proporties of the default random number generator in Julia
#at the time of writing this is the MersenneTwister; if this changes it is likely that tests below fill fail
srand(0)
@test rand() == 0.8236475079774124
@test GeneralizedMetropolisHastings.sample_indicator_matrix(repmat([0.5 0.5],2),3,1) == [1,2,1,1]
@test GeneralizedMetropolisHastings.sample_indicator_matrix(repmat([0.5 0.25 0.125 0.125],4),5,3) == [3,1,1,1,1,1]
@test GeneralizedMetropolisHastings.sample_indicator_matrix(repmat([0.5 0.25 0.125 0.125],4),5,1) == [1,4,2,2,1,4]
