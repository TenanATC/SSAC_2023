##This is the cleaned-up code for Submission to SSAC 2023
##Using the ShotTracker data for Causal Inference


library(mgcv)
library(gratia)
library(lme4)
library(dplyr)
library(WeightIt)
library(cobalt)
library(survey)
library(interactions)
library(RColorBrewer)
library(ggplot2)

#Place the file from Github in the same folder as the R file
men_dat <- read.csv('./sloan_shottracker_pseudodat.csv')

#making sure that team id, opponent id, and player_id are factors
men_dat$tid <- as.factor(men_dat$tid)
men_dat$opp_1 <- as.factor(men_dat$opp_1)
men_dat$pid <- as.factor(men_dat$pid)

#Calculate Propensity score
#Calculate it two different ways as described here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5157938/
#Assess both the marginal stability weight and the cluster-mean stability weight
num.mod <-gam((dist_1) ~ s(pid, bs= 're'), data = men_dat, method = 'REML')  # numerator model for option 1
den.mod <- gam((dist_1) ~ s(pid, bs= 're')+  s(gametime, k=30) + s(fg2_rate, k=30) + s(fg3_rate, k=30) + home_team_name + tid + opp_1,
               data= men_dat, method = 'REML') # denominator model for both options
num.p_op1 <- dnorm((men_dat$dist_1), mean = fitted(num.mod), sd = last(variance_comp(num.mod)$std_dev)) # numerator calculation for option 1 
num.p_op2 <- dnorm((men_dat$dist_1), mean = mean((men_dat$dist_1)), sd = sd((men_dat$dist_1))) # numerator calculation for option 2
den.p <- dnorm((men_dat$dist_1), mean = fitted(den.mod), sd = last(variance_comp(den.mod)$std_dev)) 
Tps_MLM_op1 <- num.p_op1/den.p
Tps_MLM_op2 <- num.p_op2/den.p

men_dat$ps1_opt1 <- Tps_MLM_op1
men_dat$ps1_opt2 <- Tps_MLM_op2



#balance diagnostics from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6484705/
mbb_covs <- men_dat %>% select(gametime, fg2_rate, fg3_rate, home_team_name, tid, opp_1) %>% data.frame() 
ps1_weigh <- as.weightit(men_dat$ps1_opt1, treat= men_dat$dist_1, estimand='ATE', covs = mbb_covs)
ps2_weigh <- as.weightit(men_dat$ps1_opt2, treat= men_dat$dist_1, estimand='ATE', covs = mbb_covs)

bal.tab(ps1_weigh, stats = c("c", "m"), un = TRUE, thresholds = c(cor = .1), poly = 3)
bal.tab(ps2_weigh, stats = c("c", "m"), un = TRUE, thresholds = c(cor = .1), poly = 3)

##Trim at 95% tails and re-check balance!  
ps1_opt1_95 <- quantile(men_dat$ps1_opt1, probs = 0.95)
ps1_opt2_95 <- quantile(men_dat$ps1_opt2, probs = 0.95)

men_dat$ps1_opt1_trim <- ifelse(men_dat$ps1_opt1 > ps1_opt1_95, ps1_opt1_95, men_dat$ps1_opt1)
men_dat$ps1_opt2_trim <- ifelse(men_dat$ps1_opt2 > ps1_opt2_95, ps1_opt2_95, men_dat$ps1_opt2)

ps1_weigh2 <- as.weightit(men_dat$ps1_opt1_trim, treat= (men_dat$dist_1), estimand='ATE', covs = mbb_covs)
ps2_weigh2 <- as.weightit(men_dat$ps1_opt2_trim, treat= (men_dat$dist_1), estimand='ATE', covs = mbb_covs)

bal.tab(ps1_weigh2, stats = c("c", "m"), un = TRUE, thresholds = c(cor = .1), poly = 3)
bal.tab(ps2_weigh2, stats = c("c", "m"), un = TRUE, thresholds = c(cor = .1), poly = 3)

#Both with trimming and not trimming, the cluster-mean stabilized Propensity score was not balanced.
#Therefore, it was decided not to use this stabilization method

#The un-trimmed Marginal Weight was ultimately used as this facilitates the use of the raw propensity score and
#actually resulted in a slightly lower treatment correlation

#Love plot for graphical depiction that would need to be cleaned up if used for full paper
love.plot(ps2_weigh, stats = c("c"), thresholds = c(cor = .1), 
          abs = TRUE, wrap = 20,
          var.order = "unadjusted", line = TRUE)


#Let's do the Causal analysis!
#set up survey design
psw.design <- svydesign(ids=~pid,#cluster ID
                        nest=T,#indicating that data is actually nested
                        weights =~ps1_opt2,#the weight variable
                        data=men_dat)#the dataset

#fit outcome model with Taylor series to estimate standard errors
psw.outcome.model <- svyglm(formula = is_make~ shot_distance*dist_1,#formula
                            family=quasibinomial(),
                            design=psw.design) #Forces package to use the design described above
summary(psw.outcome.model)

#Use 'interactions' package to explore effects
y_hyps <- seq(2, 30, by= 3)

probe_interaction(psw.outcome.model, pred = dist_1, modx = shot_distance, modx.values = y_hyps,
                  cond.int = TRUE, interval = TRUE,  jnplot = TRUE, control.fdr = TRUE, johnson_neyman = TRUE)

#According to the Johnson-Neyman plot, as long as the shot is greater than 30 feet from hoop
#there is a significant effect of defender distance on making the shot.
#My guess is once you get outside 30 feet, there are too many other factors to making the shot (including random chance)

#Something like this is what we want for the Sloan paper
x <- interact_plot(psw.outcome.model, pred = dist_1, modx = shot_distance, modx.values = y_hyps,
                   interval = F, vary.lty = F, colors = RColorBrewer::brewer.pal(10, name = 'RdYlBu'), legend.main = 'Shot Distance\nfrom Hoop (ft)') + 
  xlab('Opponent Proximity to Shot (ft)') +  ylab('Probability of Making Shot') 


ggsave('SSAC_Defense_Plot1.tiff', x, width = 5, height = 5, units = 'in')
