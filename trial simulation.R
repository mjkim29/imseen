---
# Creation date: 13 Oct 2023
# Purpose: Conduct two-arm adaptive trial simulation
---
  
### library 
library(gtools)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)
library(rjags)
library(coda)
library(progress)



### create data frame

df = data.frame()
final_df = data.frame()



### Create loop that assigns p1 and p2 in each arm 
plist <- c(100:105)
for (p in plist) {set <- p
                  p_range <- c(0:100)
                  p_range_5 <-seq(5,95,by=5)
                  p_combination <- as.data.frame(permutations(n=100,r=2,v=p_range,repeats.allowed=T))
                  p_combination <- p_combination %>%
                                   subset(V2>=V1 & V2-V1<=10 & V1 %in% p_range_5) %>%
                                   mutate(p_a=V1/100,
                                          p_b=V2/100,
                                          diff = V2-V1,
                                          p=c(1:202))
a_pop_p <- p_combination$p_a[p]
b_pop_p <- p_combination$p_b[p]



### specify trial design features 
# rd_prior: choice of prior (default: neutral prior (precision 0.33))
# nth_by  : frequency of analysis
# eff     : stopping rule for efficacy threshold
# fut     : stopping rule for equivalence threshold
# fut_diff: stopping rule for equivalence indifference zone
rd_prior <- 0.33 
nth_by   <- 100
eff <- 0.95
fut <-0.80 
fut_diff <- 0.01



### number of simulations 
n_iter <- 1000  



### random number generator  
set.seed(1004)
rnlist <- runif(500)     



### Create loop that repeats 1000 simulations for the same set of a_pop_p and b_pop_p  
for (t in 1:n_iter) {# set seed
                     set.seed(rnlist[p])
                     seedlist <- sample(1:2^20, 2000)
                     mainseed <- rnlist[p]*seedlist[t]
                     set.seed(mainseed)
  
                     # data frame for random number generation 
                     order <- c(1:1e4)
                     r <- sort(floor(runif(1e4,1,1e9)),decreasing=FALSE)
                     dataselect <- data.frame(order,r)

                     
                     
### initialize data frames
dataline = list()
basedata <- data.frame(id=numeric(),
                       group=numeric(),
                       r=numeric())
  

### create loop that calculates posterior probability at interim analysis, following the order in "dataselect" 
for (i in 1:10000) { # change seed to the random number generated above 
                    set.seed(dataselect$r[i])
    
                    # create a data frame that corresponds to a fixed sample size in "dataselect"
                     id       <- c(1:nth_by) 
                     group    <- rep(1:2,nth_by/2)
                     r        <- runif(nth_by,0,1)
                     adddata  <- data.frame(group,id,r)
                     basedata <- rbind(basedata,adddata)
    
                     # overall data
                     # in Group 1, assign success(1) if random prop is Pop A%, assign failure(0) if else
                     # in Group 2, assign success(1) if random prop is Pop B%, assign failure(0) if else
                     basedata_o <- basedata %>%
                     mutate(outcome_gr1 = ifelse(group==1 & r<=a_pop_p, 1, 0),
                            outcome_gr2 = ifelse(group==2 & r<=b_pop_p, 1, 0),
                            outcome     = ifelse(group==1,outcome_gr1,outcome_gr2))
    
                     # summarize data at Nth observation 
                     basedata_agg <- basedata_o %>%
                                     group_by(group) %>%
                                     summarise(n_enrolled = n(),
                                               n_outcome_y = sum(outcome),
                                               n_outcome_n = n_enrolled - sum(n_outcome_y),
                                               outcome_rate = n_outcome_y / n_enrolled)
    
    

### define the model to estimate posterior distribution 
model1 <- "model{Y0 ~ dbin(r0, X0) # Data model for unexposed
                 logit(r0) <- lr0  # define the logit of r0
 
                 Y1 ~ dbin(r1, X1) # Data model for exposed
                 logit(r1) <- lr1  # define the logit of r1
 
                 lr1 <- lr0 + lor  # define lr1 as unexposed logit(risk) plus log(OR)
       
                 lr0 ~ dnorm(mean.lr0, pr.lr0)  # Priors for logit of risk of unexposed 
                 lor ~ dnorm(mean.lor, pr.lor)  # Priors for log(OR) 
       
                 rd <- r1 - r0  # Computation of comparison statistics: Risk Difference  
                 rr <- r1/r0    # Computation of comparison statistics: Risk Ratio           
          
                 }"
    


### compile the model    
model1_jags <- jags.model(textConnection(model1), 
                           data =  list(X0=basedata_agg$n_enrolled[1], X1=basedata_agg$n_enrolled[2], 
                                        Y0=basedata_agg$n_outcome_y[1], Y1=basedata_agg$n_outcome_y[2], 
                                           mean.lr0=0, pr.lr0=0.3, mean.lor=0, pr.lor=rd_prior),
                          inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 1234))
    


### estimate posterior
invisible(capture.output( sim_r0 <- coda.samples(model = model1_jags, variable.names = c("r0"), n.iter = 10000) ))
invisible(capture.output( sim_r1 <- coda.samples(model = model1_jags, variable.names = c("r1"), n.iter = 10000) ))
invisible(capture.output( sim_rd <- coda.samples(model = model1_jags, variable.names = c("rd"), n.iter = 10000) ))
invisible(capture.output( sim_rr <- coda.samples(model = model1_jags, variable.names = c("rr"), n.iter = 10000) ))
sim_r0_vector <- unlist(sim_r0)
sim_r1_vector <- unlist(sim_r1)
sim_rd_vector <- unlist(sim_rd)
    
    
### calculate probability of being the best arm 
a_superior <- sum(sim_rd_vector<0, na.rm=TRUE) / 10000
b_superior <- sum(sim_rd_vector>0, na.rm=TRUE) / 10000
arm_prob <- data.frame(a_superior,b_superior)
    
    
    
### stop for efficacy 
diff_05  <- quantile(sim_rd_vector, probs = c(1-eff)) #default = 0.05; 5th percentile
diff_95  <- quantile(sim_rd_vector, probs = c(eff)) #default = 0.95; 95th percentile
trial_end <- arm_prob %>%
             rowwise() %>%
              mutate(end1=diff_05>=0 | diff_95<=(0) )
    
    
### stop for equivalence 
diff_025  <- quantile(sim_rd_vector, probs = c((1-fut)/2)) #default = 0.025
diff_975  <- quantile(sim_rd_vector, probs = c(fut+(1-fut)/2)) #default = 0.975
trial_end <- trial_end %>%
             rowwise() %>%
             mutate(end2=diff_025>(-fut_diff) & diff_975<fut_diff)   #default = -0.01 and 0.01
    
    
    
### stop if maximum sample size is reached
a_n <- basedata_agg$n_enrolled[1]
b_n <- basedata_agg$n_enrolled[2]
a_n_pos <- basedata_agg$n_outcome_y[1]
b_n_pos <- basedata_agg$n_outcome_y[2]
total_n <- a_n+b_n     
trial_end <- trial_end %>%
             rowwise() %>%
             mutate(end3=total_n>=100000)     
    
    
    
### summarize
end1 <- as.vector(trial_end$end1) 
end2 <- as.vector(trial_end$end2)   
end3 <- as.vector(trial_end$end3)   
    
med_a    <- median(sim_r0_vector)
med_b    <- median(sim_r1_vector)
    
rd_mean <- mean(sim_rd_vector)                   
rd_ci   <- quantile(sim_rd_vector, probs = c(0.025, 0.975))   # 95% credible interval 
rd_ci_l  <- rd_ci[1]
rd_ci_h  <- rd_ci[2]
rd_inc_true <- (rd_ci_l<=(b_pop_p-a_pop_p) & (b_pop_p-a_pop_p)<=rd_ci_h)
    
a_data_p=basedata_agg$outcome_rate[1]
b_data_p=basedata_agg$outcome_rate[2]
    
dataline[[i]] <- data.frame(a_pop_p,b_pop_p,a_data_p,b_data_p,med_a,med_b,a_n,b_n,total_n,a_n_pos,b_n_pos,a_superior,b_superior,rd_mean,rd_ci_l,rd_ci_h,rd_inc_true,end1,end2,end3) 
    


### summarise information into a single data frame
prob_by_n = do.call(rbind, dataline)
final_n = tail(prob_by_n,n=1)
    


###stop the loop if either of the stopping rules was triggered
i[final_n$end1==TRUE | final_n$end2==TRUE | final_n$end3==TRUE] <- 0
if (i==0) {
           break 
           } 
    
}    
  
  
### summarize one simulation 
a_data_p <- as.vector(final_n$a_data_p)     
b_data_p <- as.vector(final_n$b_data_p)     
  
a_med_p <- as.vector(final_n$med_a)     
b_med_p <- as.vector(final_n$med_b)     

a_n <- as.vector(final_n$a_n)     
b_n <- as.vector(final_n$b_n)   
total_n <- as.vector(final_n$total_n)     
a_n_pos <- as.vector(final_n$a_n_pos)     
b_n_pos <- as.vector(final_n$b_n_pos)        
  
a_superior <- as.vector(final_n$a_superior)     
b_superior <- as.vector(final_n$b_superior)     
  
rd_mean <- as.vector(final_n$rd_mean)    
rd_ci_l  <- as.vector(final_n$rd_ci_l)   
rd_ci_h <- as.vector(final_n$rd_ci_h)          
  
end1 <- as.vector(final_n$end1) 
end2 <- as.vector(final_n$end2)   
end3 <- as.vector(final_n$end3)   

final_n <- data.frame(set,mainseed,rd_prior,nth_by,eff,fut,fut_diff,a_pop_p,b_pop_p,a_data_p,b_data_p,a_n,b_n,total_n,a_n_pos,b_n_pos,a_med_p,b_med_p,a_superior,b_superior,rd_mean,rd_ci_l,rd_ci_h,rd_inc_true,end1,end2,end3)
  
df = rbind(df, final_n)
  
} 

final_df = rbind(final_df, df)
final_df <- distinct(final_df)

}   


### summary of all simulations 
final_df <- final_df %>%
            mutate(result_a=(end1==TRUE & a_med_p>b_med_p),
                   result_b=(end1==TRUE & a_med_p<b_med_p),
                   result_nodiff=(end2==TRUE),
                   result_inconc=(end3==TRUE))



### summary of simulations by effect size 
final_df_sum <- final_df %>%
                group_by(set) %>%
                summarise(rd_prior=mean(rd_prior),
                          nth_by=mean(nth_by),
                          eff=mean(eff),
                          fut=mean(fut),
                          fut_diff=mean(fut_diff),
                          a_pop_p=mean(a_pop_p),
                          b_pop_p=mean(b_pop_p),
                          a_data_p=mean(a_data_p),
                          b_data_p=mean(b_data_p),
                          a_med_p=mean(a_med_p),
                          b_med_p=mean(b_med_p),
                          total=n(),
                          end1=sum(end1),
                          end2=sum(end2),
                          end3=sum(end3),
                          rd_mean=mean(rd_mean),
                          rd_inc_true=sum(rd_inc_true, na.rm = TRUE),
                          result_a=sum(result_a, na.rm = TRUE),
                          result_b=sum(result_b, na.rm = TRUE),
                          result_nodiff=sum(result_nodiff, na.rm = TRUE),
                          result_inconc=sum(result_inconc, na.rm = TRUE),      
                          samplesize_med=median(total_n),
                          samplesize_25q=quantile(total_n, probs = 0.25),
                          samplesize_75q=quantile(total_n, probs = 0.75),
                          samplesize_min=min(total_n),
                          samplesize_max=max(total_n))





