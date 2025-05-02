# Bayesian adaptive trial designs for data-driven quality improvement in health services: a simulation study

The R script (trial simulation.R) for generating and analysing dataset



1. specify trial design features 
rd_prior: choice of prior (default: neutral prior (precision 0.33))
nth_by  : frequency of analysis
eff     : stopping rule for efficacy threshold
fut     : stopping rule for equivalence threshold
fut_diff: stopping rule for equivalence indifference zone
* For the choice of design features, please refer to the manuscript of the study.



2. output dataframe (final_df) lists the results of all simulations
a_pop_p: true success probability in arm A
b_pop_p: true success probability in arm B
a_data_p: observed success probability in arm A
b_data_p: observed success probability in arm B 
a_n: : number of participants assigned to arm A
b_n: : number of participants assigned to arm B 
total_n: total sample size
a_n_pos: number of participants assigned to arm A who had a successful outcome
b_n_pos: number of participants assigned to arm B who had a successful outcome
a_med_p: median of the posterior distribution of the success probability in arm A
b_med_p: median of the posterior distribution of the success probability in arm B
a_superior: probability that arm A is superior to arm B 
b_superior: probability that arm B is superior to arm A 
rd_mean: mean of the posterior distribution of the effect difference 
rd_ci_l: 95% CI (2.5th percentile) of the posterior distribution of the effect difference 
rd_ci_h: 95% CI (97.5th percentile) of the posterior distribution of the effect difference 
rd_inc_true: does the posterior distribution of the effect difference include the true value?
end1: trial stopped for efficacy? (y/n)
end2: trial stopped for equivalence? (y/n)
end3: trial stopped by reaching the sample size limit? (y/n)



3. output dataframe (final_df_sum) summarizes final_df by the true effect difference between the two arms
a_pop_p: true success probability in arm A
b_pop_p: true success probability in arm B
a_data_p: mean of all simulations observed success probability in arm A
b_data_p: mean of all simulations observed success probability in arm B
a_med_p: mean of all simulations median of the posterior distribution of the success probability in arm A
b_med_p: mean of all simulations median of the posterior distribution of the success probability in arm B
total: number of simulations
end1: total number of simulations stopped for efficacy
end2: total number of simulations stopped for equivalence
end3: total number of simulations reached sample size limit 
rd_mean: mean of all simulations mean of the posterior distribution of the effect difference 
rd_inc_true: total number of simulations where the posterior distribution of the effect difference include the true value
result_a: total number of simulations where arm A is declared superior
result_b: total number of simulations where arm B is declared superior
result_nodiff: total number of simulations where equivalence is declared
result_incoc: total number of simulations where results were inconclusive (i.e. reached sample size limit) 
samplesize_med: median of the sample sizes of all simulations
samplesize_25q: IQR (25th percentile) of the sample sizes of all simulations
samplesize_75q: IQR (75h percentile) of the sample sizes of all simulations
samplesize_min: minimum value of the sample sizes of all simulations
samplesize_max: maximum value of the sample sizes of all simulations
