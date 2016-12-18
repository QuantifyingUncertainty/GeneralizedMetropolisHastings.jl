#Define some quantities and a very simple TargetModel for all tests below
rinitseed1 = 879
rnburnin1 = 2
rniter1 = 4
smhnprops1 = 1 #number of test proposals for standard mh
gmhnprops1 = 2 #number of test proposals for generalized mh
rtime1 = linspace(0.0,10.0,100)
rlower1 = [2.5]
rupper1 = [3.5]
rdefault1 = [3.0]
rcov1 = [0.01]
rparas1 = parameters([:a],rlower1,rupper1,rdefault1)
rdata1 = data(:array,rtime1,sin(3*rtime1))
rnoise1 = noise(:gaussian,rcov1)
rmodel1 = model(:target,rparas1,rdata1,rnoise1,(t,p)->sin(p[1]*t);name="RunnersTestModel")
rnparas1 = numparas(rmodel1)
rsampler1 = sampler(:mh,:normal,0.1,eye(1))
rtuner1 = tuner(:monitor,1)
rtuner2 = tuner(:scale,2,0.5,:logistic)

