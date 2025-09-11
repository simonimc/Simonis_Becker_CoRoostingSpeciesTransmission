#---
#title: CST_20250819
#author: Molly Simonis
#date: 2025-08-19
#---

#turn on packages
library(deSolve)
library(ggplot2)
library(patchwork)
library(ggtern)
library(metR)


#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

#make general theme for plots later
th=theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(plot.title = element_text(hjust = 0.5))

#remove scientific notation
options(scipen=999)

#read in a tables of parameter combinations 2sp
params_combo_2sp<- read.csv('ParamCombos_v19_2sp.csv', header = T, sep = ',')
#make sure column names are correct3
names(params_combo_2sp)<- c('beta', 'psi', 'gamma', 'epsilon', 'lambda', 'period')

sub_params_2sp<- as.data.frame(rbind(params_combo_2sp[c(158, 221), ]))


################2 species Models#######################
#Create function for two species SIRS model
SIRS_model_2sp<- function(time, current_states, params){
  
  ## state the state variables
  Sa<- current_states[1] 
  Ia<- current_states[2] 
  Ra<- current_states[3] 
  
  Sb<- current_states[4] 
  Ib<- current_states[5] 
  Rb<- current_states[6] 
  
  
  ## write equations with parameters as a list
  with(as.list(params),{
    
    ## ODEs
    
    #rates of change for vampire bats
    dS_a<- ((b0 - (b1*(Sa + Sb + Ia + Ib + Ra + Rb)))*(Sa + Ia + Ra)) - ((beta*Sa*Ia) + (theta*Sa*Ib)) - (mu*Sa) + (epsilon*Ra)
    
    dI_a<- ((beta*Sa*Ia) + (theta*Sa*Ib)) - (Ia*(mu + gamma)) 
    
    dR_a<- (gamma*Ia) - (Ra*(mu + epsilon))
    
    #rates of change for other species 1
    dS_b<-((b0 - (b1*(Sa + Sb + Ia + Ib + Ra + Rb)))*(Sb + Ib + Rb)) - ((beta*Sb*Ib) + (theta*Sb*Ia)) - (mu*Sb) + (epsilon*Rb)
    
    dI_b<-((beta*Sb*Ib) + (theta*Sb*Ia)) - (Ib*(mu + gamma)) 
    
    dR_b<- (gamma*Ib) - (Rb*(mu + epsilon))
    
    
    ## make a vector of state variables, as list
    dy_SIRS<- c(dS_a, dI_a, dR_a, dS_b, dI_b, dR_b)
    list(dy_SIRS)})
}



#Disease free models first
#SIRS 2sp

## starting population values 
Sa0<- 1000 
Ia0<- 0
Ra0<- 0

Sb0<- 1000
Ib0<- 0
Rb0<- 0

Y0<- c(Sa0, Ia0, Ra0, Sb0, Ib0, Rb0)

##lifespan
mu<- (1/(12*365)) #1 death every 12 years typically 

k<- 3000 #carrying capacity 
b0<- 2/365 #2 bats per year birth rate
b1<- (b0-mu)/k

##define simulation run time in days (running for 85 years)
tmax<- 300*365
by<- 30 
time<- seq(0, tmax, by = by)




#make empty lists to put output dataframes into
SIRSeq_datalist<- list() #for equilibria data
SIRS_out_tail<- list() #for checking data reached equilibria

#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params_2sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp[i, 1]
  lambda<- params_combo_2sp[i, 5]
  psi<- params_combo_2sp[i, 2]
  theta<- beta * exp(lambda*(1 - psi))
  gamma<- params_combo_2sp[i, 3]
  epsilon<- params_combo_2sp[i, 4]
  
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, theta = theta, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model_2sp, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", 
                      "Sb", "Ib", "Rb")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib) / 
    (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
       SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb)
  
  SIRS_out$N<- SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
    SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb
  
  SIRS_out$Na<- SIRS_out$Sa + SIRS_out$Ia + SIRS_out$Ra
  
  SIRS_out$Nb<- SIRS_out$Sb + SIRS_out$Ib + SIRS_out$Rb
  
  #get all time series data
  SIRSeq_datalist[[i]]<- SIRS_out
  
  ##check that final equillibrium timesteps gives a difference of <0.001 individuals within the total pop
  SIRS_out_tail[[i]]<- SIRS_out$N[SIRS_out$time == 31025] - SIRS_out$N[SIRS_out$time == 31020]
  
  param_cat<- sub_params_2sp[i, 6]
  
  print(i)
}

#check all parameter combos made it to equilibrium
SIRS_out_tail_df<- do.call(rbind, SIRS_out_tail)


#make time series data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

SIRSeq_df$period<- NA
SIRSeq_df[1:3651, 12]<- 'short'
SIRSeq_df[3652:7302, 12]<- 'long'


SIRSeq_df$period<- factor(SIRSeq_df$period, levels = c('short', 'long'))


#create variable for period facet
SIRSeq_df$periodFacet[SIRSeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SIRSeq_df$periodFacet[SIRSeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SIRSeq_df$periodFacet<- factor(SIRSeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


#SIRS plot
SIRS_2sp_DiseaseFree<- ggplot(data = SIRSeq_df, aes(x = time, y = N)) +
  facet_grid(~periodFacet, labeller = label_parsed) +
  geom_path(aes(color = 'N'), lwd = 1, linetype = 'solid') +
  geom_path(aes(x = time, y = Na, color = 'Na'), linetype = 'solid') +
  geom_path(aes(x = time, y = Nb, color = 'Nb'), linetype = 'dashed') +
  scale_color_manual(name = "Species Population", 
                     values = c("N" = "black", "Na" = '#C7E020FF', 'Nb' = '#287C8CFF')) +
  th +
  labs(y = 'Total Roost Population (N)', x = 'Time (days)', color = "Species Population", 
       linetype = "Species Population", title = "SIRS 2 species")







#Example time series with infected individuals
## starting population values 
Sa0<- 999 
Ia0<- 1
Ra0<- 0

Sb0<- 999
Ib0<- 1
Rb0<- 0

Y0<- c(Sa0, Ia0, Ra0, Sb0, Ib0, Rb0)


#make empty lists to put output dataframes into
SIRSeq_datalist<- list() #for equilibria data
SIRS_out_tail<- list() #for checking data reached equilibria

#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params_2sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp[i, 1]
  lambda<- params_combo_2sp[i, 5]
  psi<- params_combo_2sp[i, 2]
  theta<- beta * exp(lambda*(1 - psi))
  gamma<- params_combo_2sp[i, 3]
  epsilon<- params_combo_2sp[i, 4]
  
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, theta = theta, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model_2sp, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", 
                      "Sb", "Ib", "Rb")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib) / 
    (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
       SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb)
  
  SIRS_out$N<- SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
    SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb
  
  SIRS_out$Na<- SIRS_out$Sa + SIRS_out$Ia + SIRS_out$Ra
  
  SIRS_out$Nb<- SIRS_out$Sb + SIRS_out$Ib + SIRS_out$Rb
  
  SIRS_out$I<- SIRS_out$Ia + SIRS_out$Ib 
  
  #get all time series data
  SIRSeq_datalist[[i]]<- SIRS_out
  
  ##check that final equillibrium timesteps gives a difference of <0.001 individuals within the total pop
  SIRS_out_tail[[i]]<- SIRS_out$N[SIRS_out$time == 31025] - SIRS_out$N[SIRS_out$time == 31020]
  
  param_cat<- sub_params_2sp[i, 6]
  
  print(i)
}

#check all parameter combos made it to equilibrium
SIRS_out_tail_df<- do.call(rbind, SIRS_out_tail)


#make time series data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

SIRSeq_df$period<- NA
SIRSeq_df[1:3651, 13]<- 'short'
SIRSeq_df[3652:7302, 13]<- 'long'


SIRSeq_df$period<- factor(SIRSeq_df$period, levels = c('short', 'long'))


#create variable for period facet
SIRSeq_df$periodFacet[SIRSeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SIRSeq_df$periodFacet[SIRSeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SIRSeq_df$periodFacet<- factor(SIRSeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


#SIRS plot
SIRS_2sp_InfTime<- ggplot(data = SIRSeq_df, aes(x = time, y = I)) +
  facet_grid(~periodFacet, labeller = label_parsed) +
  geom_path(aes(color = 'I'), lwd = 1, linetype = 'solid') +
  geom_path(aes(x = time, y = Ia, color = 'Ia'), linetype = 'solid') +
  geom_path(aes(x = time, y = Ib, color = 'Ib'), linetype = 'dashed') +
  scale_color_manual(name = "Infected Species Population", 
                     values = c("I" = "black", "Ia" = '#C7E020FF', 'Ib' = '#287C8CFF')) +
  th +
  labs(y = 'Total Infected Population (I)', x = 'Time (days)', 
       color = 'Infected Species Population', linetype = 'Infected Species Population', 
       title = 'SIRS 2 species')






#Full parameter space model


#make empty lists to put output dataframes into
SIRSeq_datalist<- list() #for equilibria data
SIRS_out_tail<- list() #for checking data reached equilibria

#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_2sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp[i, 1]
  lambda<- params_combo_2sp[i, 5]
  psi<- params_combo_2sp[i, 2]
  theta<- beta * exp(lambda*(1 - psi))
  gamma<- params_combo_2sp[i, 3]
  epsilon<- params_combo_2sp[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, theta = theta, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model_2sp, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", 
                      "Sb", "Ib", "Rb")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib) / 
    (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
       SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb)
  
  SIRS_out$N<- SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
    SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb
  
  
  #get all equilibrium only data
  SIRSeq_datalist[[i]]<- SIRS_out[SIRS_out$time == 109500, ]
  
  ##check that final equillibrium timesteps gives a difference of <0.001 individuals within the total pop
  SIRS_out_tail[[i]]<- SIRS_out$N[SIRS_out$time == 109500] - SIRS_out$N[SIRS_out$time == 109470]
  
  print(i)
}

#check all parameter combos made it to equilibrium
SIRS_out_tail_df<- do.call(rbind, SIRS_out_tail)


#make equilibrium data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset 
SIRSeq_params_df_2sp<- data.frame(SIRSeq_df, params_combo_2sp)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SIRSeq_params_df_2sp$ips<- 1/SIRSeq_params_df_2sp$gamma

SIRSeq_params_df_2sp$ips_cat<- round(SIRSeq_params_df_2sp$ips)

SIRSeq_params_df_2sp$ips_cat[SIRSeq_params_df_2sp$ips_cat == 5]<- '5 days infectious'
SIRSeq_params_df_2sp$ips_cat[SIRSeq_params_df_2sp$ips_cat == 548]<- '548 days infectious'


SIRSeq_params_df_2sp$imm_dur<- 1/SIRSeq_params_df_2sp$epsilon

SIRSeq_params_df_2sp$imm_dur_cat<- round(SIRSeq_params_df_2sp$imm_dur)

SIRSeq_params_df_2sp$imm_dur_cat[SIRSeq_params_df_2sp$imm_dur_cat == 60]<- '60 days immune'
SIRSeq_params_df_2sp$imm_dur_cat[SIRSeq_params_df_2sp$imm_dur_cat == 1643]<- '1643 days immune'


SIRSeq_params_df_2sp$psi_perc<- SIRSeq_params_df_2sp$psi*100


SIRSeq_params_df_2sp$beta<- as.character(SIRSeq_params_df_2sp$beta)
SIRSeq_params_df_2sp$beta<- factor(SIRSeq_params_df_2sp$beta, 
                                   levels = c('0.0001', '0.0005', '0.001'))

SIRSeq_params_df_2sp$lambda<- factor(SIRSeq_params_df_2sp$lambda, levels = c('-1', '-5', '-10'))

#SIRSeq_params_df_2sp$period<- factor(SIRSeq_params_df_2sp$period, levels = c('short', 'long'))

#SIRSeq_params_df_2sp$lambda<- factor(SIRSeq_params_df_2sp$lambda, levels = c('-1', '-5', '-10'))

#create variable for lambda facet
SIRSeq_params_df_2sp$lambdaFacet[SIRSeq_params_df_2sp$lambda == -1]<- 'lambda*" = -1"'
SIRSeq_params_df_2sp$lambdaFacet[SIRSeq_params_df_2sp$lambda == -5]<- 'lambda*" = -5"'
SIRSeq_params_df_2sp$lambdaFacet[SIRSeq_params_df_2sp$lambda == -10]<- 'lambda*" = -10"'

SIRSeq_params_df_2sp$lambdaFacet<- factor(SIRSeq_params_df_2sp$lambdaFacet, 
                                          levels = c('lambda*" = -10"', 
                                                     'lambda*" = -5"', 'lambda*" = -1"'))
#create variable for period facet
SIRSeq_params_df_2sp$periodFacet[SIRSeq_params_df_2sp$period == 'short']<- 
  '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SIRSeq_params_df_2sp$periodFacet[SIRSeq_params_df_2sp$period == 'long']<- 
  '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SIRSeq_params_df_2sp$periodFacet<- factor(SIRSeq_params_df_2sp$periodFacet, 
                                          levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                                     '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


#SIRS plot
SIRS_2sp<- ggplot(data = SIRSeq_params_df_2sp, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(periodFacet ~ lambdaFacet, scales = 'free', labeller = label_parsed) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  geom_segment(data = SIRSeq_params_df_2sp[SIRSeq_params_df_2sp$lambda == -10 & 
                                            SIRSeq_params_df_2sp$period == 'short',], 
              aes(x = 87.5, xend = 87.5, y = 0, yend = 0.045), 
              color = 'gray45', lwd = 0.75, linetype = 'dashed') +
  geom_segment(data = SIRSeq_params_df_2sp[SIRSeq_params_df_2sp$lambda == -5 & 
                                             SIRSeq_params_df_2sp$period == 'short',], 
               aes(x = 77.5, xend = 77.5, y = 0, yend = 0.045), 
               color = 'gray45', lwd = 0.75, linetype = 'dashed') +
  geom_segment(data = SIRSeq_params_df_2sp[SIRSeq_params_df_2sp$lambda == -1 & 
                                             SIRSeq_params_df_2sp$period == 'short',], 
               aes(x = 0.001, xend = 0.001, y = 0, yend = 0.045), 
               color = 'gray45', lwd = 0.75, linetype = 'dashed') +
  scale_color_viridis_d(begin = 0.2, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = 'SIRS 2 species') 


SIRS_2sp


#revert beta and lambda to a continuous value
SIRSeq_params_df_2sp$beta<- as.character(SIRSeq_params_df_2sp$beta)
SIRSeq_params_df_2sp$beta<- as.numeric(SIRSeq_params_df_2sp$beta)

SIRSeq_params_df_2sp$lambda<- as.character(SIRSeq_params_df_2sp$lambda)
SIRSeq_params_df_2sp$lambda<- as.numeric(SIRSeq_params_df_2sp$lambda)

#redfine theta
SIRSeq_params_df_2sp$theta<- SIRSeq_params_df_2sp$beta * 
                              exp(SIRSeq_params_df_2sp$lambda*(1 - SIRSeq_params_df_2sp$psi))

#quantify next-gen R0 for SIRS 2 spp
SIRSeq_params_df_2sp$R0<- (k*(SIRSeq_params_df_2sp$beta + SIRSeq_params_df_2sp$theta)) / 
                          (2*(mu + SIRSeq_params_df_2sp$gamma))







#expand SIRS grid to make a heatmap
heatmap_df_2sp<- read.csv('heatmap_params.csv', sep = ',', header = T)

#run expanded param space through SIRS 2sp model
SIRSeq_datalist<- list() #for equilibria data
SIRS_out_tail<- list() #for checking data reached equilibria

#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(heatmap_df_2sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- heatmap_df_2sp[i, 1]
  lambda<- heatmap_df_2sp[i, 5]
  psi<- heatmap_df_2sp[i, 2]
  theta<- beta * exp(lambda*(1 - psi))
  gamma<- heatmap_df_2sp[i, 3]
  epsilon<- heatmap_df_2sp[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, theta = theta, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model_2sp, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", 
                      "Sb", "Ib", "Rb")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib) / 
    (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
       SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb)
  
  SIRS_out$N<- SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
    SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb
  
  
  #get all equilibrium only data
  SIRSeq_datalist[[i]]<- SIRS_out[SIRS_out$time == 109500, ]
  
  ##check that final equillibrium timesteps gives a difference of <0.001 individuals within the total pop
  SIRS_out_tail[[i]]<- SIRS_out$N[SIRS_out$time == 109500] - SIRS_out$N[SIRS_out$time == 109470]
  
  print(i)
}

#check all parameter combos made it to equilibrium
SIRS_out_tail_df<- do.call(rbind, SIRS_out_tail)


#make equilibrium data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset 
SIRS_heatmap_2sp<- data.frame(SIRSeq_df, heatmap_df_2sp)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SIRS_heatmap_2sp$ips<- 1/SIRS_heatmap_2sp$gamma

SIRS_heatmap_2sp$ips_cat<- round(SIRS_heatmap_2sp$ips)

SIRS_heatmap_2sp$ips_cat[SIRS_heatmap_2sp$ips_cat == 5]<- '5 days infectious'
SIRS_heatmap_2sp$ips_cat[SIRS_heatmap_2sp$ips_cat == 548]<- '548 days infectious'


SIRS_heatmap_2sp$imm_dur<- 1/SIRS_heatmap_2sp$epsilon

SIRS_heatmap_2sp$imm_dur_cat<- round(SIRS_heatmap_2sp$imm_dur)

SIRS_heatmap_2sp$imm_dur_cat[SIRS_heatmap_2sp$imm_dur_cat == 60]<- '60 days immune'
SIRS_heatmap_2sp$imm_dur_cat[SIRS_heatmap_2sp$imm_dur_cat == 1643]<- '1643 days immune'


SIRS_heatmap_2sp$psi_perc<- SIRS_heatmap_2sp$psi*100

#SIRS_heatmap_2sp$lambda<- factor(SIRS_heatmap_2sp$lambda, levels = c('-1', '-5', '-10'))

#SIRS_heatmap_2sp$period<- factor(SIRS_heatmap_2sp$period, levels = c('short', 'long'))

#SIRS_heatmap_2sp$lambda<- factor(SIRS_heatmap_2sp$lambda, levels = c('-1', '-5', '-10'))

#create variable for lambda facet
SIRS_heatmap_2sp$lambdaFacet[SIRS_heatmap_2sp$lambda == -1]<- 'lambda*" = -1"'
SIRS_heatmap_2sp$lambdaFacet[SIRS_heatmap_2sp$lambda == -5]<- 'lambda*" = -5"'
SIRS_heatmap_2sp$lambdaFacet[SIRS_heatmap_2sp$lambda == -10]<- 'lambda*" = -10"'

SIRS_heatmap_2sp$lambdaFacet<- factor(SIRS_heatmap_2sp$lambdaFacet, 
                                          levels = c('lambda*" = -10"', 
                                                     'lambda*" = -5"', 'lambda*" = -1"'))
#create variable for period facet
SIRS_heatmap_2sp$periodFacet[SIRS_heatmap_2sp$periods == 'short']<- 
  '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SIRS_heatmap_2sp$periodFacet[SIRS_heatmap_2sp$periods == 'long']<- 
  '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SIRS_heatmap_2sp$periodFacet<- factor(SIRS_heatmap_2sp$periodFacet, 
                                          levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                                     '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


#redfine theta
SIRS_heatmap_2sp$theta<- SIRS_heatmap_2sp$beta * 
  exp(SIRS_heatmap_2sp$lambda*(1 - SIRS_heatmap_2sp$psi))

#quantify next-gen R0 for SIRS 2 spp
SIRS_heatmap_2sp$R0<- (k*(SIRS_heatmap_2sp$beta + SIRS_heatmap_2sp$theta)) / 
  (2*(mu + SIRS_heatmap_2sp$gamma))





SIRS_heatmap_2sp<- as.data.frame(SIRS_heatmap_2sp)


#need to make a df so ggplot doesn't read entries as duplicates
SIRS_2sp_contour_df <- SIRS_heatmap_2sp %>%
  group_by(periodFacet, lambdaFacet, periods, lambda, psi_perc, beta) %>%
  summarise(
    R0 = mean(R0),              # safe: one value per grid cell per facet
    prev_roost = mean(prev_roost),
    .groups = "drop")

#rework lambda and period facet levels
SIRS_2sp_contour_df$lambdaFacet[SIRS_2sp_contour_df$lambda == -1]<- 'lambda*" = -1"'
SIRS_2sp_contour_df$lambdaFacet[SIRS_2sp_contour_df$lambda == -5]<- 'lambda*" = -5"'
SIRS_2sp_contour_df$lambdaFacet[SIRS_2sp_contour_df$lambda == -10]<- 'lambda*" = -10"'

SIRS_2sp_contour_df$lambdaFacet<- factor(SIRS_2sp_contour_df$lambdaFacet, 
                                      levels = c('lambda*" = -10"', 
                                                 'lambda*" = -5"', 'lambda*" = -1"'))
#create variable for period facet
SIRS_2sp_contour_df$periodFacet[SIRS_2sp_contour_df$periods == 'short']<- 
  '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SIRS_2sp_contour_df$periodFacet[SIRS_2sp_contour_df$periods == 'long']<- 
  '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SIRS_2sp_contour_df$periodFacet<- factor(SIRS_2sp_contour_df$periodFacet, 
                                      levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                                 '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


SIRS_heatmap_short<- ggplot(data = SIRS_2sp_contour_df[SIRS_2sp_contour_df$periods == 'short', ],
       aes(x = psi_perc, y = beta)) +
  facet_grid(periodFacet ~ lambdaFacet, labeller = label_parsed) +
  geom_tile(aes(fill = prev_roost)) +
  geom_contour(aes(z = R0),color = "gray80", breaks = 1) +
  scale_fill_viridis_c(begin = 0.3, end = 0.9) +
  th +
  labs(x = "% Relatedness", y = expression(beta), fill = "Roost Prevalence", title = "SIRS 2 species")


SIRS_heatmap_long<- ggplot(data = SIRS_2sp_contour_df[SIRS_2sp_contour_df$periods == 'long', ],
                           aes(x = psi_perc, y = beta)) +
  facet_grid(periodFacet ~ lambdaFacet, labeller = label_parsed) +
  geom_tile(aes(fill = prev_roost)) +
  geom_contour(aes(z = R0),color = "gray80", breaks = c(0, 1)) +
  scale_fill_viridis_c(begin = 0.3, end = 0.9) +
  th +
  labs(x = "% Relatedness", y = expression(beta), fill = "Roost Prevalence")


#Figure 3A
SIRS_2sp_heat<- SIRS_heatmap_short / SIRS_heatmap_long + plot_layout(axes = 'collect')



















#Create function for two species SILI model
SILI_model_2sp<- function(time, current_states, params){
  
  ## state the state variables
  Sa<- current_states[1] 
  Ia<- current_states[2] 
  La<- current_states[3] 
  
  Sb<- current_states[4] 
  Ib<- current_states[5] 
  Lb<- current_states[6] 
  
  
  ## write equations with parameters as a list
  with(as.list(params),{
    
    ## ODEs
    
    #rates of change for vampire bats
    dS_a<- ((b0 - b1*(Sa + Sb + Ia + Ib + La + Lb))*(Sa + Ia + La)) - ((beta*Sa*Ia) + (theta*Sa*Ib)) - (mu*Sa) 
    
    dI_a<- ((beta*Sa*Ia) + (theta*Sa*Ib)) - (Ia*(mu + gamma)) + (omega*La)
    
    dL_a<- (gamma*Ia) - (La*(mu + omega))
    
    #rates of change for other species 1
    dS_b<-((b0 - b1*(Sa + Sb + Ia + Ib + La + Lb))*(Sb + Ib + Lb)) - ((beta*Sb*Ib) + (theta*Sb*Ia)) - (mu*Sb) 
    
    dI_b<-((beta*Sb*Ib) + (theta*Sb*Ia)) - (Ib*(mu + gamma)) + (omega*Lb)
    
    dL_b<- (gamma*Ib) - (Lb*(mu + omega))
    
    
    ## make a vector of state variables, as list
    dy_SILI<- c(dS_a, dI_a, dL_a, dS_b, dI_b, dL_b)
    list(dy_SILI)})
}


#disease-free SILI 2 sp
## starting population values 
Sa0<- 1000 
Ia0<- 0
La0<- 0

Sb0<- 1000
Ib0<- 0
Lb0<- 0

Y0<- c(Sa0, Ia0, La0, Sb0, Ib0, Lb0)


#make empty list to put output dataframes into
SILIeq_datalist<- list()
SILI_out_tail<- list()

#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params_2sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp[i, 1]
  lambda<- params_combo_2sp[i, 5]
  psi<- params_combo_2sp[i, 2]
  theta<- beta * exp(lambda*(1 -psi))
  gamma<- params_combo_2sp[i, 3]
  omega<- params_combo_2sp[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, omega = omega, 
           mu = mu, lambda = lambda, k = k, theta = theta, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model_2sp, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", 
                      "Sb", "Ib", "Lb")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib) / 
    (SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
       SILI_out$Ib + SILI_out$La + SILI_out$Lb)
  
  SILI_out$N<- SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
    SILI_out$Ib + SILI_out$La + SILI_out$Lb
  
  SILI_out$Na<- SILI_out$Sa + SILI_out$Ia + SILI_out$La 
  
  SILI_out$Nb<- SILI_out$Sb + SILI_out$Ib + SILI_out$Lb
  
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- SILI_out
  
  ##check that final equillibrium timesteps gives a difference of <0.001 individuals within the total pop
  SILI_out_tail[[i]]<- SILI_out$N[SILI_out$time == 109500] - SILI_out$N[SILI_out$time == 109470]
  
  print(i)
}

#check all parameter combos made it to equilibrium
SILI_out_tail_df<- do.call(rbind, SILI_out_tail)

#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)



SILIeq_df$period<- NA
SILIeq_df[1:3651, 12]<- 'short'
SILIeq_df[3652:7302, 12]<- 'long'


SILIeq_df$period<- factor(SILIeq_df$period, levels = c('short', 'long'))

#create variable for period facet
SILIeq_df$periodFacet[SILIeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~omega*" = 1/60"'
SILIeq_df$periodFacet[SILIeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'

SILIeq_df$periodFacet<- factor(SILIeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~omega*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'))


#SILI plot
SILI_2sp_DiseaseFree<- ggplot(data = SILIeq_df, aes(x = time, y = N)) +
  facet_grid(~periodFacet, labeller = label_parsed) +
  geom_path(aes(color = 'N'), lwd = 1, linetype = 'solid') +
  geom_path(aes(x = time, y = Na, color = 'Na'), linetype = 'solid') +
  geom_path(aes(x = time, y = Nb, color = 'Nb'), linetype = 'dashed') +
  scale_color_manual(name = "Species Population", 
                     values = c("N" = "black", "Na" = '#C7E020FF', 'Nb' = '#287C8CFF')) +
  th +
  labs(y = 'Total Roost Population (N)', x = 'Time (days)', color = "Species Population", 
       linetype = "Species Population", title = 'SILI 2 species')




#Example with infected indidviduals
## starting population values 
Sa0<- 999 
Ia0<- 1
La0<- 0

Sb0<- 999
Ib0<- 1
Lb0<- 0

Y0<- c(Sa0, Ia0, La0, Sb0, Ib0, Lb0)


#make empty list to put output dataframes into
SILIeq_datalist<- list()
SILI_out_tail<- list()

#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params_2sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp[i, 1]
  lambda<- params_combo_2sp[i, 5]
  psi<- params_combo_2sp[i, 2]
  theta<- beta * exp(lambda*(1 - psi))
  gamma<- params_combo_2sp[i, 3]
  omega<- params_combo_2sp[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, omega = omega, 
           mu = mu, lambda = lambda, k = k, theta = theta, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model_2sp, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", 
                      "Sb", "Ib", "Lb")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib) / 
    (SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
       SILI_out$Ib + SILI_out$La + SILI_out$Lb)
  
  SILI_out$N<- SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
    SILI_out$Ib + SILI_out$La + SILI_out$Lb
  
  SILI_out$Na<- SILI_out$Sa + SILI_out$Ia + SILI_out$La 
  
  SILI_out$Nb<- SILI_out$Sb + SILI_out$Ib + SILI_out$Lb
  
  SILI_out$I<- SILI_out$Ia + SILI_out$Ib
  
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- SILI_out
  
  ##check that final equillibrium timesteps gives a difference of <0.001 individuals within the total pop
  SILI_out_tail[[i]]<- SILI_out$N[SILI_out$time == 109500] - SILI_out$N[SILI_out$time == 109470]
  
  print(i)
}

#check all parameter combos made it to equilibrium
SILI_out_tail_df<- do.call(rbind, SILI_out_tail)

#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)



SILIeq_df$period<- NA
SILIeq_df[1:3651, 13]<- 'short'
SILIeq_df[3652:7302, 13]<- 'long'


SILIeq_df$period<- factor(SILIeq_df$period, levels = c('short', 'long'))

#create variable for period facet
SILIeq_df$periodFacet[SILIeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~omega*" = 1/60"'
SILIeq_df$periodFacet[SILIeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'

SILIeq_df$periodFacet<- factor(SILIeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~omega*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'))


#SILI plot
SILI_2sp_InfTime<- ggplot(data = SILIeq_df, aes(x = time, y = I)) +
  facet_grid(~periodFacet, labeller = label_parsed) +
  geom_path(aes(color = 'I'), lwd = 1, linetype = 'solid') +
  geom_path(aes(x = time, y = Ia, color = 'Ia'), linetype = 'solid') +
  geom_path(aes(x = time, y = Ib, color = 'Ib'), linetype = 'dashed') +
  scale_color_manual(name = "Infected Species Population", 
                     values = c("I" = "black", "Ia" = '#C7E020FF', 'Ib' = '#287C8CFF')) +
  th +
  labs(y = 'Total Infected Population (I)', x = 'Time (days)', 
       color = "Infected Species Population", linetype = "Infected Species Population", 
       title = 'SILI 2 species')




#Full SILI 2 sp model

#make empty list to put output dataframes into
SILIeq_datalist<- list()
SILI_out_tail<- list()

#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_2sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp[i, 1]
  lambda<- params_combo_2sp[i, 5]
  psi<- params_combo_2sp[i, 2]
  theta<- beta * exp(lambda*(1 - psi))
  gamma<- params_combo_2sp[i, 3]
  omega<- params_combo_2sp[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, omega = omega, 
           mu = mu, lambda = lambda, k = k, theta = theta, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model_2sp, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", 
                      "Sb", "Ib", "Lb")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib) / 
    (SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
       SILI_out$Ib + SILI_out$La + SILI_out$Lb)
  
  SILI_out$N<- SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
    SILI_out$Ib + SILI_out$La + SILI_out$Lb
  
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- SILI_out[SILI_out$time == 109500, ]
  
  ##check that final equillibrium timesteps gives a difference of <0.001 individuals within the total pop
  SILI_out_tail[[i]]<- SILI_out$N[SILI_out$time == 109500] - SILI_out$N[SILI_out$time == 109470]
  
  print(i)
}

#check all parameter combos made it to equilibrium
SILI_out_tail_df<- do.call(rbind, SILI_out_tail)

#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset to create plots later
SILIeq_params_df_2sp<- data.frame(SILIeq_df, params_combo_2sp)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SILIeq_params_df_2sp$ips<- 1/SILIeq_params_df_2sp$gamma

SILIeq_params_df_2sp$ips_cat<- round(SILIeq_params_df_2sp$ips)

SILIeq_params_df_2sp$ips_cat[SILIeq_params_df_2sp$ips_cat == 5]<- '5 days infectious'
SILIeq_params_df_2sp$ips_cat[SILIeq_params_df_2sp$ips_cat == 548]<- '548 days infectious'

#rename epsilon in param df to omega
colnames(SILIeq_params_df_2sp)[13]<- 'omega'
SILIeq_params_df_2sp$lat_dur<- 1/SILIeq_params_df_2sp$omega 

SILIeq_params_df_2sp$lat_dur_cat<- round(SILIeq_params_df_2sp$lat_dur)

SILIeq_params_df_2sp$lat_dur_cat[SILIeq_params_df_2sp$lat_dur_cat == 60]<- '60 days latent'
SILIeq_params_df_2sp$lat_dur_cat[SILIeq_params_df_2sp$lat_dur_cat == 1643]<- '1643 days latent'


SILIeq_params_df_2sp$psi_perc<- SILIeq_params_df_2sp$psi*100


SILIeq_params_df_2sp$beta<- as.character(SILIeq_params_df_2sp$beta)
SILIeq_params_df_2sp$beta<- factor(SILIeq_params_df_2sp$beta, 
                                   levels = c('0.0001', '0.0005', '0.001'))

SILIeq_params_df_2sp$lambda<- factor(SILIeq_params_df_2sp$lambda, levels = c('-1', '-5', '-10'))

#SILIeq_params_df_2sp$period<- factor(SILIeq_params_df_2sp$period, levels = c('short', 'long'))

#create variable for lambda facet
SILIeq_params_df_2sp$lambdaFacet[SILIeq_params_df_2sp$lambda == -1]<- 'lambda*" = -1"'
SILIeq_params_df_2sp$lambdaFacet[SILIeq_params_df_2sp$lambda == -5]<- 'lambda*" = -5"'
SILIeq_params_df_2sp$lambdaFacet[SILIeq_params_df_2sp$lambda == -10]<- 'lambda*" = -10"'

SILIeq_params_df_2sp$lambdaFacet<- factor(SILIeq_params_df_2sp$lambdaFacet, 
                                          levels = c('lambda*" = -10"', 
                                                     'lambda*" = -5"', 'lambda*" = -1"'))
#create variable for period facet
SILIeq_params_df_2sp$periodFacet[SILIeq_params_df_2sp$period == 'short']<- 
  '"short:"~gamma*" = 1/5,"~omega*" = 1/60"'
SILIeq_params_df_2sp$periodFacet[SILIeq_params_df_2sp$period == 'long']<- 
  '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'

SILIeq_params_df_2sp$periodFacet<- factor(SILIeq_params_df_2sp$periodFacet, 
                                          levels = c('"short:"~gamma*" = 1/5,"~omega*" = 1/60"', 
                                                     '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'))


#SILI plot
SILI_2sp<- ggplot(data = SILIeq_params_df_2sp, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(periodFacet ~ lambdaFacet, scales = 'free', labeller = label_parsed) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(begin = 0.2, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = 'SILI 2 species') 


SILI_2sp




#revert beta to a continuous value
SILIeq_params_df_2sp$beta<- as.character(SILIeq_params_df_2sp$beta)
SILIeq_params_df_2sp$beta<- as.numeric(SILIeq_params_df_2sp$beta)

SILIeq_params_df_2sp$lambda<- as.character(SILIeq_params_df_2sp$lambda)
SILIeq_params_df_2sp$lambda<- as.numeric(SILIeq_params_df_2sp$lambda)

#redfine theta
SILIeq_params_df_2sp$theta<- SILIeq_params_df_2sp$beta * 
  exp(SILIeq_params_df_2sp$lambda*(1 - SILIeq_params_df_2sp$psi))

#quantify next-gen R0 for SILI 2 spp
SILIeq_params_df_2sp$R0<- (k*(SILIeq_params_df_2sp$beta + SILIeq_params_df_2sp$theta) *
                              (mu + SILIeq_params_df_2sp$omega)) / 
                          (2*(mu*(mu + SILIeq_params_df_2sp$gamma + SILIeq_params_df_2sp$omega)))






#rename expanded params space for SILI models
names(heatmap_df_2sp)<- c('beta', 'psi', 'gamma', 'omega', 'lambda', 'periods')

#run expanded param space through SILI 2sp model
SILIeq_datalist<- list() #for equilibria data
SILI_out_tail<- list() #for checking data reached equilibria

#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(heatmap_df_2sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- heatmap_df_2sp[i, 1]
  lambda<- heatmap_df_2sp[i, 5]
  psi<- heatmap_df_2sp[i, 2]
  theta<- beta * exp(lambda*(1 - psi))
  omega<- heatmap_df_2sp[i, 3]
  epsilon<- heatmap_df_2sp[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, omega = omega, 
           mu = mu, lambda = lambda, theta = theta, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model_2sp, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", 
                      "Sb", "Ib", "Lb")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib) / 
    (SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
       SILI_out$Ib + SILI_out$La + SILI_out$Lb)
  
  SILI_out$N<- SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
    SILI_out$Ib + SILI_out$La + SILI_out$Lb
  
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- SILI_out[SILI_out$time == 109500, ]
  
  ##check that final equillibrium timesteps gives a difference of <0.001 individuals within the total pop
  SILI_out_tail[[i]]<- SILI_out$N[SILI_out$time == 109500] - SILI_out$N[SILI_out$time == 109470]
  
  print(i)
}

#check all parameter combos made it to equilibrium
SILI_out_tail_df<- do.call(rbind, SILI_out_tail)


#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset 
SILI_heatmap_2sp<- data.frame(SILIeq_df, heatmap_df_2sp)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SILI_heatmap_2sp$ips<- 1/SILI_heatmap_2sp$gamma

SILI_heatmap_2sp$ips_cat<- round(SILI_heatmap_2sp$ips)

SILI_heatmap_2sp$ips_cat[SILI_heatmap_2sp$ips_cat == 5]<- '5 days infectious'
SILI_heatmap_2sp$ips_cat[SILI_heatmap_2sp$ips_cat == 548]<- '548 days infectious'


SILI_heatmap_2sp$lat_dur<- 1/SILI_heatmap_2sp$omega

SILI_heatmap_2sp$lat_dur_cat<- round(SILI_heatmap_2sp$lat_dur)

SILI_heatmap_2sp$lat_dur_cat[SILI_heatmap_2sp$lat_dur_cat == 60]<- '60 days latent'
SILI_heatmap_2sp$lat_dur_cat[SILI_heatmap_2sp$lat_dur_cat == 1643]<- '1643 days latent'


SILI_heatmap_2sp$psi_perc<- SILI_heatmap_2sp$psi*100

#SILI_heatmap_2sp$lambda<- factor(SILI_heatmap_2sp$lambda, levels = c('-1', '-5', '-10'))

#SILI_heatmap_2sp$period<- factor(SILI_heatmap_2sp$period, levels = c('short', 'long'))

#SILI_heatmap_2sp$lambda<- factor(SILI_heatmap_2sp$lambda, levels = c('-1', '-5', '-10'))

#create variable for lambda facet
SILI_heatmap_2sp$lambdaFacet[SILI_heatmap_2sp$lambda == -1]<- 'lambda*" = -1"'
SILI_heatmap_2sp$lambdaFacet[SILI_heatmap_2sp$lambda == -5]<- 'lambda*" = -5"'
SILI_heatmap_2sp$lambdaFacet[SILI_heatmap_2sp$lambda == -10]<- 'lambda*" = -10"'

SILI_heatmap_2sp$lambdaFacet<- factor(SILI_heatmap_2sp$lambdaFacet, 
                                      levels = c('lambda*" = -10"', 
                                                 'lambda*" = -5"', 'lambda*" = -1"'))
#create variable for period facet
SILI_heatmap_2sp$periodFacet[SILI_heatmap_2sp$periods == 'short']<- 
  '"short:"~gamma*" = 1/5,"~omega*" = 1/60"'
SILI_heatmap_2sp$periodFacet[SILI_heatmap_2sp$periods == 'long']<- 
  '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'

SILI_heatmap_2sp$periodFacet<- factor(SILI_heatmap_2sp$periodFacet, 
                                      levels = c('"short:"~gamma*" = 1/5,"~omega*" = 1/60"', 
                                                 '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'))


#redfine theta
SILI_heatmap_2sp$theta<- SILI_heatmap_2sp$beta * 
  exp(SILI_heatmap_2sp$lambda*(1 - SILI_heatmap_2sp$psi))

#quantify next-gen R0 for SILI 2 spp
SILI_heatmap_2sp$R0<-  (k*(SILIeq_params_df_2sp$beta + SILIeq_params_df_2sp$theta) *
                               (mu + SILIeq_params_df_2sp$omega)) / 
  (2*(mu*(mu + SILIeq_params_df_2sp$gamma + SILIeq_params_df_2sp$omega)))






SILI_heatmap_2sp<- as.data.frame(SILI_heatmap_2sp)


#need to make a df so ggplot doesn't read entries as duplicates
SILI_2sp_contour_df <- SILI_heatmap_2sp %>%
  group_by(periodFacet, lambdaFacet, periods, lambda, psi_perc, beta) %>%
  summarise(
    R0 = mean(R0),              # safe: one value per grid cell per facet
    prev_roost = mean(prev_roost),
    .groups = "drop")

#rework lambda and period facet levels
SILI_2sp_contour_df$lambdaFacet[SILI_2sp_contour_df$lambda == -1]<- 'lambda*" = -1"'
SILI_2sp_contour_df$lambdaFacet[SILI_2sp_contour_df$lambda == -5]<- 'lambda*" = -5"'
SILI_2sp_contour_df$lambdaFacet[SILI_2sp_contour_df$lambda == -10]<- 'lambda*" = -10"'

SILI_2sp_contour_df$lambdaFacet<- factor(SILI_2sp_contour_df$lambdaFacet, 
                                         levels = c('lambda*" = -10"', 
                                                    'lambda*" = -5"', 'lambda*" = -1"'))
#create variable for period facet
SILI_2sp_contour_df$periodFacet[SILI_2sp_contour_df$periods == 'short']<- 
  '"short:"~gamma*" = 1/5,"~omega*" = 1/60"'
SILI_2sp_contour_df$periodFacet[SILI_2sp_contour_df$periods == 'long']<- 
  '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'

SILI_2sp_contour_df$periodFacet<- factor(SILI_2sp_contour_df$periodFacet, 
                                         levels = c('"short:"~gamma*" = 1/5,"~omega*" = 1/60"', 
                                                    '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'))


SILI_heatmap_short<- ggplot(data = SILI_2sp_contour_df[SILI_2sp_contour_df$periods == 'short', ],
                            aes(x = psi_perc, y = beta)) +
  facet_grid(periodFacet ~ lambdaFacet, labeller = label_parsed) +
  geom_tile(aes(fill = prev_roost)) +
  geom_contour(aes(z = R0),color = "gray80", breaks = 1) +
  scale_fill_viridis_c(begin = 0.3, end = 0.9) +
  th +
  labs(x = "% Relatedness", y = expression(beta), fill = "Roost Prevalence", title = "SILI 2 species")


SILI_heatmap_long<- ggplot(data = SILI_2sp_contour_df[SILI_2sp_contour_df$periods == 'long', ],
                           aes(x = psi_perc, y = beta)) +
  facet_grid(periodFacet ~ lambdaFacet, labeller = label_parsed) +
  geom_tile(aes(fill = prev_roost)) +
  geom_contour(aes(z = R0),color = "gray80", breaks = c(0, 1)) +
  scale_fill_viridis_c(begin = 0.3, end = 0.9) +
  th +
  labs(x = "% Relatedness", y = expression(beta), fill = "Roost Prevalence")

#Figure 3B
SILI_2sp_heat<- SILI_heatmap_short / SILI_heatmap_long + plot_layout(axes = 'collect')






#Figure 3C and 3D
SIRS_2sp_prev<- (SIRS_2sp + SILI_2sp) + plot_layout(guides = 'collect', axes = 'collect') 



#Figure S2
(SIRS_2sp_DiseaseFree / SILI_2sp_DiseaseFree) + plot_layout(guides = 'collect')

#Figure S4
(SIRS_2sp_InfTime / SILI_2sp_InfTime) + plot_layout(guides = 'collect')
















#########################3 spp Models##################################################

#clear space for 3spp model data 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

#make general theme for plots later
th=theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(plot.title = element_text(hjust = 0.5))

#remove scientific notation
options(scipen=999)

#read in parameter space for 3 spp
params_combo_3sp<- read.csv('ParamCombos_v19_3sp.csv', header = T, sep = ',')
#make sure column names are correct3
names(params_combo_3sp)<- c('psi_ab', 'psi_ac', 'psi_bc', 'beta', 'gamma', 'epsilon', 'lambda', 'period')

sub_params_3sp<- as.data.frame(params_combo_3sp[rownames(params_combo_3sp[
  params_combo_3sp$psi_ab == 0.5 & params_combo_3sp$psi_ac == 0.5 & params_combo_3sp$psi_bc == 0.5 & 
    params_combo_3sp$beta == 0.0005,]), ])

#read in 3sp model data
load('OSCER_3sp_CST_models_20250820.Rdata')





#Make SIRS 3sp model
SIRS_model_3sp<- function(time, current_states, params){
  
  ## state the state variables
  Sa<- current_states[1] 
  Ia<- current_states[2]
  Ra<- current_states[3]
  
  Sb<- current_states[4]
  Ib<- current_states[5]
  Rb<- current_states[6]
  
  Sc<- current_states[7]
  Ic<- current_states[8]
  Rc<- current_states[9]
  
  ## write equations with parameters as a list
  with(as.list(params),{
    
    ## ODEs
    
    #rates of change for vampire bats
    dS_a<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + Ra + Rb + Rc))*(Sa + Ia + Ra)) - ((beta*Sa*Ia) + (theta_ab*Sa*Ib) + (theta_ac*Sa*Ic)) - (mu*Sa) + (epsilon*Ra)
    
    dI_a<- ((beta*Sa*Ia) + (theta_ab*Sa*Ib) + (theta_ac*Sa*Ic)) - (Ia*(mu + gamma))  
    
    dR_a<- (gamma*Ia) - (Ra*(mu + epsilon))
    
    #rates of change for other bat species 1
    dS_b<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + Ra + Rb + Rc))*(Sb + Ib + Rb)) - ((beta*Sb*Ib) + (theta_ab*Sb*Ia) + (theta_bc*Sb*Ic)) - (mu*Sb) + (epsilon*Rb)
    
    dI_b<- ((beta*Sb*Ib) + (theta_ab*Sb*Ia) + (theta_bc*Sb*Ic)) - (Ib*(mu + gamma))  
    
    dR_b<- (gamma*Ib) - (Rb*(mu + epsilon))
    
    #rates of change for other bat species 2
    dS_c<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + Ra + Rb + Rc))*(Sc + Ic + Rc)) - ((beta*Sc*Ic) + (theta_ac*Sc*Ia) + (theta_bc*Sc*Ib)) - (mu*Sc) + (epsilon*Rc) 
    
    dI_c<- ((beta*Sc*Ic) + (theta_ac*Sc*Ia) + (theta_bc*Sc*Ib)) - (Ic*(mu + gamma)) 
    
    dR_c<- (gamma*Ic) - (Rc*(mu + epsilon))
    
    ## make a vector of state variables, as list
    dy_SIRS<- c(dS_a, dI_a, dR_a, dS_b, dI_b, dR_b, dS_c, dI_c, dR_c)
    list(dy_SIRS)})
}



##demographic parameters for equal lifespans
mu<- (1/(12*365))
k<- 3000 #carrying capacity = 
b0<- 2/365 #2 bats per year
b1<- (b0-mu)/k


##define simulation run time in days (running for 20 years)
tmax<- 300*365
by<- 30 
time<- seq(0, tmax, by = by)



#disease-free models first
Sa0<- 1000 
Ia0<- 0
Ra0<- 0

Sb0<- 1000
Ib0<- 0
Rb0<- 0

Sc0<- 1000
Ic0<- 0
Rc0<- 0
Y0<- c(Sa0, Ia0, Ra0, Sb0, Ib0, Rb0, Sc0, Ic0, Rc0)


#make empty lists to put output dataframes into
SIRS_datalist<- list() #for equilibria data
SIRS_out_tail<- list() #for checking data reached equilibria

#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params_3sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_3sp[i, 4]
  lambda<- params_combo_3sp[i, 7]
  psi_ab<- params_combo_3sp[i, 1]
  psi_ac<- params_combo_3sp[i, 2]
  psi_bc<- params_combo_3sp[i, 3]
  theta_ab<- beta * exp(lambda*(1 - psi_ab))
  theta_ac<- beta * exp(lambda*(1 - psi_ac))
  theta_bc<- beta * exp(lambda*(1 - psi_bc))
  gamma<- params_combo_3sp[i, 5]
  epsilon<- params_combo_3sp[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, 
           gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, theta_ab = theta_ab,
           theta_ac = theta_ac, theta_bc = theta_bc, b0 = b0, b1 = b1)
  
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model_3sp, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", 
                      "Sb", "Ib", "Rb",
                      "Sc", "Ic", "Rc")
  
  #get total N
  SIRS_out$N<- SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Sc + 
    SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic + 
    SIRS_out$Ra + SIRS_out$Rb + SIRS_out$Rc
  
  SIRS_out$Na<- SIRS_out$Sa + SIRS_out$Ia + SIRS_out$Ra 
  
  SIRS_out$Nb<- SIRS_out$Sb + SIRS_out$Ib + SIRS_out$Rb 
  
  SIRS_out$Nc<- SIRS_out$Sc + SIRS_out$Ic + SIRS_out$Rc
  
  
  #get all data
  SIRS_datalist[[i]]<- SIRS_out
  
  print(i)
  
}

#check all parameter combos made it to equilibrium
SIRS_out_tail_df<- do.call(rbind, SIRS_out_tail)


#make time series data one df
SIRSeq_df<- do.call(rbind, SIRS_datalist)

SIRSeq_df$period<- NA
SIRSeq_df[1:10953, 15]<- 'short'
SIRSeq_df[10954:21906, 15]<- 'long'


SIRSeq_df$period<- factor(SIRSeq_df$period, levels = c('short', 'long'))


#create variable for period facet
SIRSeq_df$periodFacet[SIRSeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SIRSeq_df$periodFacet[SIRSeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SIRSeq_df$periodFacet<- factor(SIRSeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


#plot
SIRS_3sp_DiseaseFree<- ggplot(data = SIRSeq_df, aes(x = time, y = N)) +
  facet_grid(~periodFacet, labeller = label_parsed) +
  geom_path(aes(color = 'N'), lwd = 1, linetype = 'solid') +
  geom_path(aes(x = time, y = Na, color = 'Na'), linetype = 'solid') +
  geom_path(aes(x = time, y = Nb, color = 'Nb'), linetype = 'dashed') +
  geom_path(aes(x = time, y = Nc, color = 'Nc'), linetype = 'dotted') +
  scale_color_manual(name = "Species Population", values = c("N" = "black", "Na" = '#C7E020FF', 'Nb' = '#287C8CFF', 'Nc' = '#481F70FF')) +
  th +
  labs(y = 'Total Roost Population (N)', x = 'Time (days)', 
       color = 'Species Population', linetype = 'Species Population',
       title = 'SIRS 3 species')




#Example time series with infected individuals
## starting population values 
Sa0<- 999 
Ia0<- 1
Ra0<- 0

Sb0<- 999
Ib0<- 1
Rb0<- 0

Sc0<- 999
Ic0<- 1
Rc0<- 0
Y0<- c(Sa0, Ia0, Ra0, Sb0, Ib0, Rb0, Sc0, Ic0, Rc0)



#make empty lists to put output dataframes into
SIRS_datalist<- list() #for equilibria data
SIRS_out_tail<- list() #for checking data reached equilibria


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params_3sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_3sp[i, 4]
  lambda<- params_combo_3sp[i, 7]
  psi_ab<- params_combo_3sp[i, 1]
  psi_ac<- params_combo_3sp[i, 2]
  psi_bc<- params_combo_3sp[i, 3]
  theta_ab<- beta * exp(lambda*(1 - psi_ab))
  theta_ac<- beta * exp(lambda*(1 - psi_ac))
  theta_bc<- beta * exp(lambda*(1 - psi_bc))
  gamma<- params_combo_3sp[i, 5]
  epsilon<- params_combo_3sp[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, 
           gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, theta_ab = theta_ab,
           theta_ac = theta_ac, theta_bc = theta_bc, b0 = b0, b1 = b1)
  
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model_3sp, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", 
                      "Sb", "Ib", "Rb",
                      "Sc", "Ic", "Rc")
  
  #get total N
  SIRS_out$N<- SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Sc + 
    SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic + 
    SIRS_out$Ra + SIRS_out$Rb + SIRS_out$Rc
  
  SIRS_out$Na<- SIRS_out$Sa + SIRS_out$Ia + SIRS_out$Ra 
  
  SIRS_out$Nb<- SIRS_out$Sb + SIRS_out$Ib + SIRS_out$Rb 
  
  SIRS_out$Nc<- SIRS_out$Sc + SIRS_out$Ic + SIRS_out$Rc
  
  SIRS_out$I<- SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic 
  
  #get all data
  SIRS_datalist[[i]]<- SIRS_out
  
  print(i)
  
}

#check all parameter combos made it to equilibrium
SIRS_out_tail_df<- do.call(rbind, SIRS_out_tail)


#make time series data one df
SIRSeq_df<- do.call(rbind, SIRS_datalist)

SIRSeq_df$period<- NA
SIRSeq_df[1:10953, 16]<- 'short'
SIRSeq_df[10954:21906, 16]<- 'long'


SIRSeq_df$period<- factor(SIRSeq_df$period, levels = c('short', 'long'))

#create variable for period facet
SIRSeq_df$periodFacet[SIRSeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SIRSeq_df$periodFacet[SIRSeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SIRSeq_df$periodFacet<- factor(SIRSeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


#plot
SIRS_3sp_InfTime<- ggplot(data = SIRSeq_df, aes(x = time)) + 
  facet_grid(~periodFacet, labeller = label_parsed) +
  geom_line(aes(y = I,  color = 'I', group = 1), lwd = 1, linetype = 'solid') +
  geom_line(aes(y = Ia, color = 'Ia'), linetype = 'solid') +
  geom_line(aes(y = Ib, color = 'Ib'), linetype = 'dashed') +
  geom_line(aes(y = Ic, color = 'Ic'), linetype = 'dotted') +
  scale_color_manual(name = "Infected Species Population", 
                     values = c("I"  = "black", "Ia" = "#C7E020FF", "Ib" = "#287C8CFF", "Ic" = "#481F70FF")) +
  th +
  labs(y = 'Total Infected Population (I)', x = 'Time (days)', color = 'Infected Species Population', 
       linetype = 'Infected Species Population', title = 'SIRS 3 species')


#make beta numeric
SIRSeq_params_df_3sp$beta<- as.character(SIRSeq_params_df_3sp$beta)
SIRSeq_params_df_3sp$beta<- as.numeric(SIRSeq_params_df_3sp$beta)

#cross check  
unique(SIRSeq_params_df_3sp$beta)

#order periods
SIRSeq_params_df_3sp$period<- factor(SIRSeq_params_df_3sp$period, levels = c('short', 'long'))

#redifine thetas
SIRSeq_params_df_3sp$theta_ab<- SIRSeq_params_df_3sp$beta * 
  exp(SIRSeq_params_df_3sp$lambda*(1 - SIRSeq_params_df_3sp$psi_ab))

SIRSeq_params_df_3sp$theta_ac<- SIRSeq_params_df_3sp$beta * 
  exp(SIRSeq_params_df_3sp$lambda*(1 - SIRSeq_params_df_3sp$psi_ac))

SIRSeq_params_df_3sp$theta_bc<- SIRSeq_params_df_3sp$beta * 
  exp(SIRSeq_params_df_3sp$lambda*(1 - SIRSeq_params_df_3sp$psi_bc))


#quantify next gen R0
P<- ((SIRSeq_params_df_3sp$theta_ab * SIRSeq_params_df_3sp$theta_ab) + 
       (SIRSeq_params_df_3sp$theta_ac * SIRSeq_params_df_3sp$theta_ac) + 
       (SIRSeq_params_df_3sp$theta_bc * SIRSeq_params_df_3sp$theta_bc))/
  (SIRSeq_params_df_3sp$beta * SIRSeq_params_df_3sp$beta)


C<- (2*(SIRSeq_params_df_3sp$theta_ab*SIRSeq_params_df_3sp$theta_ac*SIRSeq_params_df_3sp$theta_bc))/
  ((SIRSeq_params_df_3sp$beta)^3)


Delta<- ((C/2)^2) - ((P/3)^3)
Delta<- as.numeric(Delta)
#Delta creates some negative numbers which will be problematic for calculating the next step using 
#a cubed root (calculating Lambda). Need to use some trig to fix this

#create cubed root function to calculate Lambda
cube_root<- function(x) sign(x) * abs(x)^(1/3)
  
## correct for calculating Lambda
Lambda<- rep(NA,length(Delta))

for (i in 1:length(Delta)) {
  if (Delta[i] >= 0) {
    
    Lambda[i]<- 1 + cube_root((C[i]/2) + sqrt(Delta[i])) + cube_root((C[i]/2) - sqrt(Delta[i]))
  } else {
    
    ## use trigonometric form (three real roots; choose the largest)
    ## ensure floating point inside acos is in [-1,1]
    arg<- (C[i]/2) / sqrt((P[i]/3)^3)
    
    ## numeric safety
    arg<- pmin(1, pmax(-1, arg))
    y<- 2 * sqrt(P[i]/3) * cos((1/3) * acos(arg))
    Lambda[i]<- 1 + y
  } 
}
#Now should have the same length of Lambdas as prev values
#calculate R0
SIRSeq_params_df_3sp$R0<- ((k*SIRSeq_params_df_3sp$beta) / 
                             (3*(mu + SIRSeq_params_df_3sp$gamma))) * Lambda




#sanity check using brute force
R0_num_datalist<- list()

for (i in 1:nrow(SIRSeq_params_df_3sp)) {
  
  S<- k/3 #disease free equal pop at eq
  
  M<- matrix(c(SIRSeq_params_df_3sp[i, 17] * S, SIRSeq_params_df_3sp[i, 29] * S, SIRSeq_params_df_3sp[i, 30] * S,
               SIRSeq_params_df_3sp[i, 29] * S, SIRSeq_params_df_3sp[i, 17] * S, SIRSeq_params_df_3sp[i, 31] * S,
               SIRSeq_params_df_3sp[i, 30] * S, SIRSeq_params_df_3sp[i, 31] * S, SIRSeq_params_df_3sp[i, 17] * S),
               nrow = 3, byrow = T)
  
  G<- M/(mu + SIRSeq_params_df_3sp[i, 18])
  
  Eigvals<- eigen(G)$values
  
  R0_num<- max(Re(Eigvals))
  
  R0_num_datalist[[i]]<- R0_num
  
}

R0_num<- do.call(rbind, R0_num_datalist)

SIRSeq_params_df_3sp$R0_num<- R0_num         
             

#all R0 should fall on 1:1 line with R0_num
ggplot(data = SIRSeq_params_df_3sp, aes(x = R0_num, y = R0, color = factor(beta))) +
  #facet_wrap(period ~ lambda, scales = 'free', labeller = label_parsed) +
  geom_jitter()+
  geom_abline(slope = 1, intercept = 0)+
  th 




##create variable for period facet for plotting
SIRSeq_params_df_3sp$periodFacet[SIRSeq_params_df_3sp$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SIRSeq_params_df_3sp$periodFacet[SIRSeq_params_df_3sp$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SIRSeq_params_df_3sp$periodFacet<- factor(SIRSeq_params_df_3sp$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


#create variable for lambda facet for plotting
SIRSeq_params_df_3sp$lambdaFacet[SIRSeq_params_df_3sp$lambda == -1]<- 'lambda*" = -1"'
SIRSeq_params_df_3sp$lambdaFacet[SIRSeq_params_df_3sp$lambda == -5]<- 'lambda*" = -5"'
SIRSeq_params_df_3sp$lambdaFacet[SIRSeq_params_df_3sp$lambda == -10]<- 'lambda*" = -10"'

SIRSeq_params_df_3sp$lambdaFacet<- factor(SIRSeq_params_df_3sp$lambdaFacet, levels = c('lambda*" = -10"', 
                                                         'lambda*" = -5"', 'lambda*" = -1"'))


#create variable for beta facet for plotting
SIRSeq_params_df_3sp$betaFacet[SIRSeq_params_df_3sp$beta == 0.0001]<- 'beta*" = 0.0001"'
SIRSeq_params_df_3sp$betaFacet[SIRSeq_params_df_3sp$beta == 0.0005]<- 'beta*" = 0.0005"'
SIRSeq_params_df_3sp$betaFacet[SIRSeq_params_df_3sp$beta == 0.001]<- 'beta*" = 0.001"'

SIRSeq_params_df_3sp$betaFacet<- factor(SIRSeq_params_df_3sp$betaFacet, levels = c('beta*" = 0.0001"', 
                                                                                       'beta*" = 0.0005"', 
                                                                                   'beta*" = 0.001"'))


#plot SIRS ternary plots and include R0

SIRSeq_params_df_3sp$R0_bin<- ifelse(SIRSeq_params_df_3sp$R0 < 1, 'R[0] < 1', 'R[0] > 1')

#Figure 4A
  ggtern(data = SIRSeq_params_df_3sp[SIRSeq_params_df_3sp$period == 'short', ], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, alpha = R0_bin, color = prev_roost)) +
    facet_grid(periodFacet*betaFacet ~ lambdaFacet, scales = 'free', labeller = label_parsed)+
    geom_point(size = 3) +
    #geom_hex_tern() +
    scale_color_viridis_c(begin = 0.2, end =  0.9) +
    scale_alpha_manual(values = c(0.05, 1))+
    th +
   theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalence', alpha = '', title = 'SIRS 3 species') +
    Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
    theme_arrowlong()
  
#Figure 4C
  ggtern(data = SIRSeq_params_df_3sp[SIRSeq_params_df_3sp$period == 'long', ], 
         aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
    facet_grid(periodFacet*betaFacet ~ lambdaFacet, scales = 'free', labeller = label_parsed)+
    geom_point(size = 3) +
    scale_color_viridis_c(begin = 0.2, end =  0.9) +
    th +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalence') +
    Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
    theme_arrowlong()


  
  
  
  
  




#SILI 3spp
#change epsilon to omega for SILI
names(params_combo_3sp)<- c('psi_ab', 'psi_ac', 'psi_bc', 'beta', 'gamma', 'omega', 'lambda', 'period')

names(sub_params_3sp)<- c('psi_ab', 'psi_ac', 'psi_bc', 'beta', 'gamma', 'omega', 'lambda', 'period')


#Create function for 3sp SILI model
SILI_model_3sp<- function(time, current_states, params){
  
  ## state the state variables
  Sa<- current_states[1] 
  Ia<- current_states[2]
  La<- current_states[3]
  
  Sb<- current_states[4]
  Ib<- current_states[5]
  Lb<- current_states[6]
  
  Sc<- current_states[7]
  Ic<- current_states[8]
  Lc<- current_states[9]
  
  ## write equations with parameters as a list
  with(as.list(params),{
    
    ## ODEs
    
    #rates of change for vampire bats
    dS_a<- ((b0 - (b1*(Sa + Sb + Sc +Ia + Ib + Ic + La + Lb + Lc)))*(Sa + Ia + La)) - ((beta*Sa*Ia) + (theta_ab*Sa*Ib) + (theta_ac*Sa*Ic)) - (mu*Sa)  
    
    dI_a<- ((beta*Sa*Ia) + (theta_ab*Sa*Ib) + (theta_ac*Sa*Ic)) - (Ia*(mu + gamma)) + (omega*La) 
    
    dL_a<- (gamma*Ia) - (La*(mu + omega))
    
    #rates of change for other bat species 1
    dS_b<- ((b0 - (b1*(Sa + Sb + Sc +Ia + Ib + Ic + La + Lb + Lc)))*(Sb + Ib + Lb)) - ((beta*Sb*Ib) + (theta_ab*Sb*Ia) + (theta_bc*Sb*Ic)) - (mu*Sb)
    
    dI_b<- ((beta*Sb*Ib) + (theta_ab*Sb*Ia) + (theta_bc*Sb*Ic)) - (Ib*(mu + gamma)) + (omega*Lb) 
    
    dL_b<- (gamma*Ib) - (Lb*(mu + omega))
    
    #rates of change for other bat species 2
    dS_c<- ((b0 - (b1*(Sa + Sb + Sc +Ia + Ib + Ic + La + Lb + Lc)))*(Sc + Ic + Lc)) - ((beta*Sc*Ic) + (theta_ac*Sc*Ia) + (theta_bc*Sc*Ib)) - (mu*Sc)
    
    dI_c<- ((beta*Sc*Ic) + (theta_ac*Sc*Ia) + (theta_bc*Sc*Ib)) - (Ic*(mu + gamma)) + (omega*Lc) 
    
    dL_c<- (gamma*Ic) - (Lc*(mu + omega))
    
    ## make a vector of state variables, as list
    dy_SILI<- c(dS_a, dI_a, dL_a, dS_b, dI_b, dL_b, dS_c, dI_c, dL_c)
    list(dy_SILI)})
}


#disease free starting pops
##starting population values
Sa0<- 1000 
Ia0<- 0
La0<- 0

Sb0<- 1000
Ib0<- 0
Lb0<- 0

Sc0<- 1000
Ic0<- 0
Lc0<- 0
Y0<- c(Sa0, Ia0, La0, Sb0, Ib0, Lb0, Sc0, Ic0, Lc0)


#same demographic parameters and time settings

#make empty lists to put output dataframes into
SILI_datalist<- list() #for equilibria data
SILI_out_tail<- list() #for checking data reached equilibria


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params_3sp)) {
  
  ## parameters identified by the row number in sub_params df
  beta<- sub_params_3sp[i, 4]
  lambda<- sub_params_3sp[i, 7]
  psi_ab<- sub_params_3sp[i, 1]
  psi_ac<- sub_params_3sp[i, 2]
  psi_bc<- sub_params_3sp[i, 3]
  theta_ab<- beta * exp(lambda*(1 - psi_ab))
  theta_ac<- beta * exp(lambda*(1 - psi_ac))
  theta_bc<- beta * exp(lambda*(1 - psi_bc))
  gamma<- sub_params_3sp[i, 5]
  omega<- sub_params_3sp[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, 
           gamma = gamma, omega = omega, 
           mu = mu, lambda = lambda, k = k, theta_ab = theta_ab,
           theta_ac = theta_ac, theta_bc = theta_bc, b0 = b0, b1 = b1)
  
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model_3sp, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", 
                      "Sb", "Ib", "Lb",
                      "Sc", "Ic", "Lc")
  
  #get total N
  SILI_out$N<- SILI_out$Sa + SILI_out$Sb + SILI_out$Sc + 
    SILI_out$Ia + SILI_out$Ib + SILI_out$Ic + 
    SILI_out$La + SILI_out$Lb + SILI_out$Lc
  
  SILI_out$Na<- SILI_out$Sa + SILI_out$Ia + SILI_out$La 
  
  SILI_out$Nb<- SILI_out$Sb + SILI_out$Ib + SILI_out$Lb 
  
  SILI_out$Nc<- SILI_out$Sc + SILI_out$Ic + SILI_out$Lc
  
  SILI_out$I<- SILI_out$Ia + SILI_out$Ib + SILI_out$Ic 
  
  #get all data
  SILI_datalist[[i]]<- SILI_out
  
  print(i)
  
}

#check all parameter combos made it to equilibrium
SILI_out_tail_df<- do.call(rbind, SILI_out_tail)


#make time series data one df
SILIeq_df<- do.call(rbind, SILI_datalist)

SILIeq_df$period<- NA
SILIeq_df[1:10953, 16]<- 'short'
SILIeq_df[10954:21906, 16]<- 'long'


SILIeq_df$period<- factor(SILIeq_df$period, levels = c('short', 'long'))

#create variable for period facet
SILIeq_df$periodFacet[SILIeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~omega*" = 1/60"'
SILIeq_df$periodFacet[SILIeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'

SILIeq_df$periodFacet<- factor(SILIeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~omega*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'))



SILI_3sp_DiseaseFree<- ggplot(data = SILIeq_df, aes(x = time, y = N)) +
  facet_grid(~periodFacet, labeller = label_parsed) +
  geom_path(aes(color = 'N'), lwd = 1, linetype = 'solid') +
  geom_path(aes(x = time, y = Na, color = 'Na'), linetype = 'solid') +
  geom_path(aes(x = time, y = Nb, color = 'Nb'), linetype = 'dashed') +
  geom_path(aes(x = time, y = Nc, color = 'Nc'), linetype = 'dotted') +
  scale_color_manual(name = "Species Population", values = c("N" = "black", "Na" = '#C7E020FF', 'Nb' = '#287C8CFF', 'Nc' = '#481F70FF')) +
  th +
  labs(y = 'Total Roost Population (N)', x = 'Time (days)', 
       color = 'Species Population', linetype = 'Species Population',
       title = 'SILI 3 species')







#Example time series with infected individuals
## starting population values 
Sa0<- 999 
Ia0<- 1
La0<- 0

Sb0<- 999
Ib0<- 1
Lb0<- 0

Sc0<- 999
Ic0<- 1
Lc0<- 0
Y0<- c(Sa0, Ia0, Ra0, Sb0, Ib0, Rb0, Sc0, Ic0, Rc0)


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params_3sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params_3sp[i, 4]
  lambda<- sub_params_3sp[i, 7]
  psi_ab<- sub_params_3sp[i, 1]
  psi_ac<- sub_params_3sp[i, 2]
  psi_bc<- sub_params_3sp[i, 3]
  theta_ab<- beta * exp(lambda*(1 - psi_ab))
  theta_ac<- beta * exp(lambda*(1 - psi_ac))
  theta_bc<- beta * exp(lambda*(1 - psi_bc))
  gamma<- sub_params_3sp[i, 5]
  omega<- sub_params_3sp[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, 
           gamma = gamma, omega = omega, 
           mu = mu, lambda = lambda, k = k, theta_ab = theta_ab,
           theta_ac = theta_ac, theta_bc = theta_bc, b0 = b0, b1 = b1)
  
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model_3sp, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", 
                      "Sb", "Ib", "Lb",
                      "Sc", "Ic", "Lc")
  
  #get total N
  SILI_out$N<- SILI_out$Sa + SILI_out$Sb + SILI_out$Sc + 
    SILI_out$Ia + SILI_out$Ib + SILI_out$Ic + 
    SILI_out$La + SILI_out$Lb + SILI_out$Lc
  
  SILI_out$Na<- SILI_out$Sa + SILI_out$Ia + SILI_out$La 
  
  SILI_out$Nb<- SILI_out$Sb + SILI_out$Ib + SILI_out$Lb 
  
  SILI_out$Nc<- SILI_out$Sc + SILI_out$Ic + SILI_out$Lc
  
  SILI_out$I<- SILI_out$Ia + SILI_out$Ib + SILI_out$Ic 
  
  #get all data
  SILI_datalist[[i]]<- SILI_out
  
  print(i)
  
}

#check all parameter combos made it to equilibrium
SILI_out_tail_df<- do.call(rbind, SILI_out_tail)


#make time series data one df
SILIeq_df<- do.call(rbind, SILI_datalist)

SILIeq_df$period<- NA
SILIeq_df[1:10953, 16]<- 'short'
SILIeq_df[10954:21906, 16]<- 'long'


SILIeq_df$period<- factor(SILIeq_df$period, levels = c('short', 'long'))

SILIeq_df$periodFacet[SILIeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~omega*" = 1/60"'
SILIeq_df$periodFacet[SILIeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'

SILIeq_df$periodFacet<- factor(SILIeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~omega*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'))

#SILI infectious pop example plot
SILI_3sp_InfTime<- ggplot(data = SILIeq_df, aes(x = time)) + 
  facet_grid(~periodFacet, labeller = label_parsed) +
  geom_line(aes(y = I,  color = 'I', group = 1), lwd = 1, linetype = 'solid') +
  geom_line(aes(y = Ia, color = 'Ia'), linetype = 'solid') +
  geom_line(aes(y = Ib, color = 'Ib'), linetype = 'dashed') +
  geom_line(aes(y = Ic, color = 'Ic'), linetype = 'dotted') +
  scale_color_manual(name = "Infected Species Population", 
                     values = c("I"  = "black", "Ia" = "#C7E020FF", "Ib" = "#287C8CFF", "Ic" = "#481F70FF")) +
  th +
  labs(y = 'Total Infected Population (I)', x = 'Time (days)', color = 'Infected Species Population', 
    linetype = 'Infected Species Population', title = 'SILI 3 species')
  
  
  


#make beta numeric
SILIeq_params_df_3sp$beta<- as.character(SILIeq_params_df_3sp$beta)
SILIeq_params_df_3sp$beta<- as.numeric(SILIeq_params_df_3sp$beta)

#cross check  
unique(SILIeq_params_df_3sp$beta)

#order periods
SILIeq_params_df_3sp$period<- factor(SILIeq_params_df_3sp$period, levels = c('short', 'long'))


#redifine thetas
SILIeq_params_df_3sp$theta_ab<- SILIeq_params_df_3sp$beta * 
  exp(SILIeq_params_df_3sp$lambda*(1 - SILIeq_params_df_3sp$psi_ab))

SILIeq_params_df_3sp$theta_ac<- SILIeq_params_df_3sp$beta * 
  exp(SILIeq_params_df_3sp$lambda*(1 - SILIeq_params_df_3sp$psi_ac))

SILIeq_params_df_3sp$theta_bc<- SILIeq_params_df_3sp$beta * 
  exp(SILIeq_params_df_3sp$lambda*(1 - SILIeq_params_df_3sp$psi_bc))



#quantify next gen R0
P<- ((SILIeq_params_df_3sp$theta_ab * SILIeq_params_df_3sp$theta_ab) + 
       (SILIeq_params_df_3sp$theta_ac * SILIeq_params_df_3sp$theta_ac) + 
       (SILIeq_params_df_3sp$theta_bc * SILIeq_params_df_3sp$theta_bc))/
  (SILIeq_params_df_3sp$beta * SILIeq_params_df_3sp$beta)

C<- (2*(SILIeq_params_df_3sp$theta_ab*SILIeq_params_df_3sp$theta_ac*SILIeq_params_df_3sp$theta_bc))/
  ((SILIeq_params_df_3sp$beta)^3)

Delta<- ((C/2)^2) - ((P/3)^3)
Delta<- as.numeric(Delta)
#Delta creates some negative numbers which will be problematic for calculating the next step using 
#a cubed root (calculating Lambda). Need to use some trig to fix this

#create cubed root function to calculate Lambda
cube_root<- function(x) sign(x) * abs(x)^(1/3)

## correct for calculating Lambda
Lambda<- rep(NA,length(Delta))

for (i in 1:length(Delta)) {
  if (Delta[i] >= 0) {
    
    Lambda[i]<- 1 + cube_root((C[i]/2) + sqrt(Delta[i])) + cube_root((C[i]/2) - sqrt(Delta[i]))
  } else {
    
    ## use trigonometric form (three real roots; choose the largest)
    ## ensure floating point inside acos is in [-1,1]
    arg<- (C[i]/2) / sqrt((P[i]/3)^3)
    
    ## numeric safety
    arg<- pmin(1, pmax(-1, arg))
    y<- 2 * sqrt(P[i]/3) * cos((1/3) * acos(arg))
    Lambda[i]<- 1 + y
  } 
}

#Now should have the same length of Lambdas as prev values
#calculate R0
SILIeq_params_df_3sp$R0<- ((k*SILIeq_params_df_3sp$beta) / 
                             (3*(mu + SILIeq_params_df_3sp$gamma))) * Lambda




#sanity check using brute force
R0_num_datalist<- list()

for (i in 1:nrow(SILIeq_params_df_3sp)) {
  
  S<- k/3 #disease free equal pop at eq
  
  M<- matrix(c(SILIeq_params_df_3sp[i, 17] * S, SILIeq_params_df_3sp[i, 29] * S, SILIeq_params_df_3sp[i, 30] * S,
               SILIeq_params_df_3sp[i, 29] * S, SILIeq_params_df_3sp[i, 17] * S, SILIeq_params_df_3sp[i, 31] * S,
               SILIeq_params_df_3sp[i, 30] * S, SILIeq_params_df_3sp[i, 31] * S, SILIeq_params_df_3sp[i, 17] * S),
             nrow = 3, byrow = T)
  
  G<- M/(mu + SILIeq_params_df_3sp[i, 18])
  
  Eigvals<- eigen(G)$values
  
  R0_num<- max(Re(Eigvals))
  
  R0_num_datalist[[i]]<- R0_num
  
}

R0_num<- do.call(rbind, R0_num_datalist)

SILIeq_params_df_3sp$R0_num<- R0_num         


#all R0 should fall on 1:1 line with R0_num
ggplot(data = SILIeq_params_df_3sp, aes(x = R0_num, y = R0, color = factor(beta))) +
  #facet_wrap(period ~ lambda, scales = 'free', labeller = label_parsed) +
  geom_jitter()+
  geom_abline(slope = 1, intercept = 0)+
  th 






#Figure S3
SIRS_3sp_DiseaseFree/ SILI_3sp_DiseaseFree + plot_layout(guides = 'collect')

#Figure S5
SIRS_3sp_InfTime / SILI_3sp_InfTime + plot_layout(guides = 'collect')



##create variable for period facet for plotting
SILIeq_params_df_3sp$periodFacet[SILIeq_params_df_3sp$period == 'short']<- '"short:"~gamma*" = 1/5,"~omega*" = 1/60"'
SILIeq_params_df_3sp$periodFacet[SILIeq_params_df_3sp$period == 'long']<- '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'

SILIeq_params_df_3sp$periodFacet<- factor(SILIeq_params_df_3sp$periodFacet, 
                                          levels = c('"short:"~gamma*" = 1/5,"~omega*" = 1/60"', 
                                                     '"long:"~gamma*" = 1/548,"~omega*" = 1/1643"'))


#create variable for lambda facet for plotting
SILIeq_params_df_3sp$lambdaFacet[SILIeq_params_df_3sp$lambda == -1]<- 'lambda*" = -1"'
SILIeq_params_df_3sp$lambdaFacet[SILIeq_params_df_3sp$lambda == -5]<- 'lambda*" = -5"'
SILIeq_params_df_3sp$lambdaFacet[SILIeq_params_df_3sp$lambda == -10]<- 'lambda*" = -10"'

SILIeq_params_df_3sp$lambdaFacet<- factor(SILIeq_params_df_3sp$lambdaFacet, levels = c('lambda*" = -10"', 
                                                                                       'lambda*" = -5"', 'lambda*" = -1"'))


#create variable for beta facet for plotting
SILIeq_params_df_3sp$betaFacet[SILIeq_params_df_3sp$beta == 0.0001]<- 'beta*" = 0.0001"'
SILIeq_params_df_3sp$betaFacet[SILIeq_params_df_3sp$beta == 0.0005]<- 'beta*" = 0.0005"'
SILIeq_params_df_3sp$betaFacet[SILIeq_params_df_3sp$beta == 0.001]<- 'beta*" = 0.001"'

SILIeq_params_df_3sp$betaFacet<- factor(SILIeq_params_df_3sp$betaFacet, levels = c('beta*" = 0.0001"', 
                                                                                   'beta*" = 0.0005"', 
                                                                                   'beta*" = 0.001"'))

#Prev figures
SILIeq_params_df_3sp$R0_bin<- ifelse(SILIeq_params_df_3sp$R0 < 1, 'R[0] < 1', 'R[0] > 1')

#Figure 4B
ggtern(data = SILIeq_params_df_3sp[SILIeq_params_df_3sp$period == 'short', ], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, alpha = R0_bin, color = prev_roost)) +
  facet_grid(periodFacet*betaFacet ~ lambdaFacet, scales = 'free', labeller = label_parsed)+
  geom_point(size = 3) +
  #geom_hex_tern() +
  scale_color_viridis_c(begin = 0.2, end =  0.9) +
  scale_alpha_manual(values = c(0.05, 1))+
  th +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalence', alpha = '', title = 'SILI 3 species') +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure 4D
ggtern(data = SILIeq_params_df_3sp[SILIeq_params_df_3sp$period == 'long', ], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(periodFacet*betaFacet ~ lambdaFacet, scales = 'free', labeller = label_parsed)+
  geom_point(size = 3) +
  scale_color_viridis_c(begin = 0.2, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalence') +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()





