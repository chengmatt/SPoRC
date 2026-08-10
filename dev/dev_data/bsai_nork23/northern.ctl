# Age specific selectivity parameters estimated in phase:
 -1
#Logistic selectivity used in phase:
 1
#Fishery selectivity parameters estimated in phase:
 3
#Survey selectivity parameters estimated in phase:
 3  
 #Sr phase
 -1
# F40 phase
 6
# Projection phase
 7
# Phase in which the historic F is estimated
 -1
# flag to use harmonic mean for McAllister-Ianelli method (1=yes, otherwise use arithmetic mean)
 1 
# flag for constraining the survey sel (0=do not constrain, 1=constrain a50 to fishery, 2= constrain slope to fishery, 3=constain sel at some age to be close to 1)
3
# age for survey selectivity constraint option 3 (ignored if not used)
30
# CV for survey selectivity constraint
0.003
#####	fishery age comp sample sizes 
# number of hauls w/ otoliths 
49	58	60	110	69	44	71	90	125	103	200 251 294 308 553 434 485
#7	7.615773106	7.745966692	10.48808848	8.306623863	6.633249581	8.426149773	9.486832981	11.18033989	10.14889157	14.14213562	15.84297952	17.1464282	17.54992877
#####	fishery lengths sample sizes 
# number of hauls w/ lengths 
16	11	68	14	42	148	219 224 311 559 780
#4	3.31662479	8.246211251	3.741657387	6.480740698	12.16552506	14.79864859	14.96662955	17.63519209	23.64318084
####	unbiased AI survey age comp sample sizes 
# number of hauls w/ otoliths 
14	19	68	92	69	65	57	72	67 60 146 147 135
#3.741657387	4.358898944	8.246211251	9.591663047	8.306623863	8.062257748	7.549834435	8.485281374	8.185352772	7.745966692	12.08304597	12.12435565
###	survey lengths sample sizes (number of hauls w/ lengths)
#number of hauls with lengths  
#240
# Lambda 1 recruitment regularity
 1.0
# Lambda 2 trawl survey abundance trend
 1.0
# Lambda 3 Catch biomass
200.
# Lambda 4 Selectivity curvature
 1.
# Lambda 5 ssq survey abundance
 1.
# Lambda 6 Avg selectivity
 1.
# Lambda 7
 1.
# Lambda 8 penalty on on selectivity dome-shape (controls descending parts across ages) 
 100. #10.
# Lambda 9 penalty on selectivity smoothness [second difference] across ages
 100.0
# Lambda 10 penalty on inter-annual variation across years [first difference]
 800.0  #400.0.
# Lambda 11 penalty on selectivity inter-annual smoothness [second difference]
 800. #400.0
 # Lambda 12-15 (weights for fac, flc, sac, slc comps)
 1 1 1 1





