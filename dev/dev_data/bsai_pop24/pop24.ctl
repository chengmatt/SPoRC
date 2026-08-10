# Age specific selectivity parameters estimated in phase:
 -1
#Logistic selectivity used in phase:
 3 3
#Fishery selectivity parameters estimated in phase:
 3
# slope Survey selectivity parameters estimated in phase (if lambda >0):
#3 3  
# a50 Survey selectivity parameters estimated in phase (if lambda >0): 
#3 3
# Survey selectivity double normal mean parameters estimated in phase (if used):
#3 3  
# Survey selectivity double normal distance parameters estimated in phase (if used):
#3 3  
# Survey selectivity double normal sig1 parameters estimated in phase (if used): 
#3 3 
# Survey selectivity double normal sig2 parameters estimated in phase (if used): 
#3 3 
# survey catchability estimated in phase (if lambda >0):
4 4
# Projection phase
 7
# Phase in which the historic F is estimated
 -1
# Phase for log_avg_M parameter:
4
# Phase for M_devs parameters:
4
# Phase for log_avg_fmort parameter:
1
# Phase for fmort_dev parameter:
2
# Phase for rec_dev parameter:
2
# Phase for log_rinit parameter:
2

# Switch for composition likelihood function (1 = multinomial, 2 = Dirichlet-Multinomial)
1
# phase of Dirichlet theta parameter, survey age comps (if used)
1 1
# phase of Dirichlet theta parameter, survey length comps (if used)
1 1
# phase of Dirichlet theta parameter, fishery age and length comps (if used)
1 1
# flag to use harmonic mean for McAllister-Ianelli method (1=yes, otherwise use arithmetic mean)
 1 
# flags to fit the biomass estimates from the surveys 
1 1
# flags to fit the abundance estimates from the surveys 
0 0
# flags to fit the age comp data from the surveys
1 1
# flags to fit the length comp data from the surveys
1 1
# flags to fit the catch biomass  
1
# flags to fit the proportional catch age comp data 
1
# flags to fit the proportional catch length comp data
1
#####	unbiased fishery age comp sample sizes
## sqrt of read otoliths (#stg1_fish_unbiased_ac_samp) 
23	15	18	29	22	23	22	20	28	24	21	26	24	27	31	32	32 33 33 34 33 33
## number of otoliths read (#nsamples_fish_unbiased_ac)
510 222 328 823 487 524 466 397 784 581 449 675 572 719 932 1043 1016 1080 1093 1150 1057 1102  

#####	fishery lengths sample sizes
## sqrt of fished lengthed (#stg1_fish_lc_samp)
282	307	304	271	313	256	196	204	173	102	123	109	63	26	49	20	41	29	128	195	187	119	108	127	103	64	89	117 127 126 157 157
## number of fishery lengths (#nsamples_fish_lc)
79749 94232 92578 73529 97946 65508 38505 41553 29924 10377 15209 11966 3995 673 2404 395 1672  822 16468 38009 34812 14200 11724 16113 10545 4128 7968    13582   16169   15916  24502 24793

##### survey age comps
## sqrt of read otoliths (#stg1_srv_ac_samp). Ragged array, rows are different surveys
32 29	35	35	18	32	21	31	34 33 33 30 35  ##(ai numbers)
17 21 20 20 22 20 ##(ebs slope numbers)
## number of otoliths read (#nsamples_srv_ac). Ragged array, rows are different surveys
1015 849 1224  1238  337 1031 462 951 1140  1078  1062  918 1204
299 425 413 415 472 400

##### survey length comps
## sqrt of fish lengthed (#stg1_srv_lc_samp). Ragged array, rows are different surveys
128  ##(ai numbers)
#58 ##(ebs slope numbers)
## number of fishery lengths (#nsamples_srv_lc). Ragged array, rows are different surveys
16448 ##(ai numbers)
# 3398   ##(ebs slope numbers)

# Lambda 1 recruitment regularity
 1.0
# Lambda 2 Catch biomass
 500.
#  Lambda 3 penalty on on selectivity dome-shape (controls descending parts across ages) 
  30 ##10.0
# Lambda 4 penalty on selectivity smoothness [second difference] across ages
  10 ## 10.0
# Lambda 5 penalty on inter-annual variation across years [first difference]
  300 ##300.0.
# Lambda 6 penalty on selectivity inter-annual smoothness [second difference]
  300 ##300.0
#


