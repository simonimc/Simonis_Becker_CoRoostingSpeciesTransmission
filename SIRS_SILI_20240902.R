#---
#title: SIRS_SILI_20240902
#author: Molly Simonis
#date: 2024-09-02
#---

#turn on packages
library(deSolve)
library(ggplot2)
library(patchwork)
library(ggtern)


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
params_combo_2sp<- read.csv('ParamCombos_v18_2sp.csv', header = T, sep = ',')
#make sure column names are correct3
names(params_combo_2sp)<- c('beta', 'psi', 'gamma', 'epsilon', 'lambda', 'period')

sub_params_2sp<- as.data.frame(rbind(params_combo_2sp[c(53, 179), ]))


#read in a tables of parameter combinations 3sp
params_combo_3sp<- read.csv('ParamCombos_v18_3sp.csv', header = T, sep = ',')
#make sure column names are correct3
names(params_combo_3sp)<- c('psi_ab', 'psi_ac', 'psi_bc', 'beta', 'gamma', 'epsilon', 'lambda', 'period')

sub_params_3sp<- as.data.frame(rbind(params_combo_3sp[c(23153, 78719), ]))

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
    dS_a<- ((b0 - b1*(Sa + Sb + Ia + Ib + Ra + Rb))*(Sa + Ia + Ra)) - ((beta*Sa*Ia) + (beta*psi*Sa*Ib)) - (mu*Sa) + (epsilon*Ra)
    
    dI_a<- ((beta*Sa*Ia) + (beta*psi*Sa*Ib)) - (Ia*(mu + gamma)) 
    
    dR_a<- (gamma*Ia) - (Ra*(mu + epsilon))
    
    #rates of change for other species 1
    dS_b<-((b0 - b1*(Sa + Sb + Ia + Ib + Ra + Rb))*(Sb + Ib + Rb)) - ((beta*Sb*Ib) + (beta*psi*Sb*Ia)) - (mu*Sb) + (epsilon*Rb)
    
    dI_b<-((beta*Sb*Ib) + (beta*psi*Sb*Ia)) - (Ib*(mu + gamma)) 
    
    dR_b<- (gamma*Ib) - (Rb*(mu + epsilon))
    
    
    ## make a vector of state variables, as list
    dy_SIRS<- c(dS_a, dI_a, dR_a, dS_b, dI_b, dR_b)
    list(dy_SIRS)})
}



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
    dS_a<- ((b0 - b1*(Sa + Sb + Ia + Ib + La + Lb))*(Sa + Ia + La)) - ((beta*Sa*Ia) + (beta*psi*Sa*Ib)) - (mu*Sa) 
    
    dI_a<- ((beta*Sa*Ia) + (beta*psi*Sa*Ib)) - (Ia*(mu + gamma)) + (epsilon*La)
    
    dL_a<- (gamma*Ia) - (La*(mu + epsilon))
    
    #rates of change for other species 1
    dS_b<-((b0 - b1*(Sa + Sb + Ia + Ib + La + Lb))*(Sb + Ib + Lb)) - ((beta*Sb*Ib) + (beta*psi*Sb*Ia)) - (mu*Sb) 
    
    dI_b<-((beta*Sb*Ib) + (beta*psi*Sb*Ia)) - (Ib*(mu + gamma)) + (epsilon*Lb)
    
    dL_b<- (gamma*Ib) - (Lb*(mu + epsilon))
    
    
    ## make a vector of state variables, as list
    dy_SILI<- c(dS_a, dI_a, dL_a, dS_b, dI_b, dL_b)
    list(dy_SILI)})
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


k<- 3000 #carrying capacity 
b0<- 2/365 #2 bats per year birth rate
b1<- b0/k

##define simulation run time in days (running for 85 years)
tmax<- 300*365
by<- 30 
time<- seq(0, tmax, by = by)


##lifespan
mu<- (1/(12*365)) #1 death every 12 years typically 


#make empty lists to put output dataframes into
SIRSeq_datalist<- list() #for equilibria data
SIRS_out_tail<- list() #for checking data reached equilibria

#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params_2sp)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params_2sp[i, 1]
  lambda<- sub_params_2sp[i, 5]
  psi<- exp(lambda*(1 - sub_params_2sp[i, 2]))
  gamma<- sub_params_2sp[i, 3]
  epsilon<- sub_params_2sp[i, 4]
  
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, b0 = b0, b1 = b1)
  
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
  beta<- sub_params_2sp[i, 1]
  lambda<- sub_params_2sp[i, 5]
  psi<- exp(lambda*(1 - sub_params_2sp[i, 2]))
  gamma<- sub_params_2sp[i, 3]
  epsilon<- sub_params_2sp[i, 4]
  
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, b0 = b0, b1 = b1)
  
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
  psi<- exp(lambda*(1 - params_combo_2sp[i, 2]))
  gamma<- params_combo_2sp[i, 3]
  epsilon<- params_combo_2sp[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, b0 = b0, b1 = b1)
  
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
                                               levels = c('0.0005', '0.001'))

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
  scale_color_viridis_d(begin = 0.2, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = 'SIRS 2 species') 








##simulation run time in days will remain the same

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
  psi<- exp(lambda*(1 - params_combo_2sp[i, 2]))
  gamma<- params_combo_2sp[i, 3]
  epsilon<- params_combo_2sp[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, b0 = b0, b1 = b1)
  
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
SILIeq_df$periodFacet[SILIeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SILIeq_df$periodFacet[SILIeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SILIeq_df$periodFacet<- factor(SILIeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


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
  psi<- exp(lambda*(1 - params_combo_2sp[i, 2]))
  gamma<- params_combo_2sp[i, 3]
  epsilon<- params_combo_2sp[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, b0 = b0, b1 = b1)
  
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
SILIeq_df$periodFacet[SILIeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SILIeq_df$periodFacet[SILIeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SILIeq_df$periodFacet<- factor(SILIeq_df$periodFacet, 
                                          levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                                     '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


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
  psi<- exp(lambda*(1 - params_combo_2sp[i, 2]))
  gamma<- params_combo_2sp[i, 3]
  epsilon<- params_combo_2sp[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, b0 = b0, b1 = b1)
  
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


SILIeq_params_df_2sp$lat_dur<- 1/SILIeq_params_df_2sp$epsilon

SILIeq_params_df_2sp$lat_dur_cat<- round(SILIeq_params_df_2sp$lat_dur)

SILIeq_params_df_2sp$lat_dur_cat[SILIeq_params_df_2sp$lat_dur_cat == 60]<- '60 days latent'
SILIeq_params_df_2sp$lat_dur_cat[SILIeq_params_df_2sp$lat_dur_cat == 1643]<- '1643 days latent'


SILIeq_params_df_2sp$psi_perc<- SILIeq_params_df_2sp$psi*100


SILIeq_params_df_2sp$beta<- as.character(SILIeq_params_df_2sp$beta)
SILIeq_params_df_2sp$beta<- factor(SILIeq_params_df_2sp$beta, 
                                   levels = c('0.0005', '0.001'))

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
  '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SILIeq_params_df_2sp$periodFacet[SILIeq_params_df_2sp$period == 'long']<- 
  '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SILIeq_params_df_2sp$periodFacet<- factor(SILIeq_params_df_2sp$periodFacet, 
                                          levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                                     '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


#SILI plot
SILI_2sp<- ggplot(data = SILIeq_params_df_2sp, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(periodFacet ~ lambdaFacet, scales = 'free', labeller = label_parsed) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(begin = 0.2, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = 'SILI 2 species') 






#combine full 2sp model plots
SIRS_2sp / SILI_2sp + plot_layout(guides = 'collect')

#combine 2sp disease-free plots
SIRS_2sp_DiseaseFree /SILI_2sp_DiseaseFree + plot_layout(guides = 'collect')

#combine 2sp infectious pop examples plots
SIRS_2sp_InfTime / SILI_2sp_InfTime + plot_layout(guides = 'collect')


################3 species Models#######################



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
    dS_a<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + Ra + Rb + Rc))*(Sa + Ia + Ra)) - ((beta*Sa*Ia) + (beta*psi_ab*Sa*Ib) + (beta*psi_ac*Sa*Ic)) - (mu*Sa) + (epsilon*Ra)
    
    dI_a<- ((beta*Sa*Ia) + (beta*psi_ab*Sa*Ib) + (beta*psi_ac*Sa*Ic)) - (Ia*(mu + gamma))  
    
    dR_a<- (gamma*Ia) - (Ra*(mu + epsilon))
    
    #rates of change for other bat species 1
    dS_b<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + Ra + Rb + Rc))*(Sb + Ib + Rb)) - ((beta*Sb*Ib) + (beta*psi_ab*Sb*Ia) + (beta*psi_bc*Sb*Ic)) - (mu*Sb) + (epsilon*Rb)
    
    dI_b<- ((beta*Sb*Ib) + (beta*psi_ab*Sb*Ia) + (beta*psi_bc*Sb*Ic)) - (Ib*(mu + gamma))  
    
    dR_b<- (gamma*Ib) - (Rb*(mu + epsilon))
    
    #rates of change for other bat species 2
    dS_c<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + Ra + Rb + Rc))*(Sc + Ic + Rc)) - ((beta*Sc*Ic) + (beta*psi_ac*Sc*Ia) + (beta*psi_bc*Sc*Ib)) - (mu*Sc) + (epsilon*Rc) 
    
    dI_c<- ((beta*Sc*Ic) + (beta*psi_ac*Sc*Ia) + (beta*psi_bc*Sc*Ib)) - (Ic*(mu + gamma)) 
    
    dR_c<- (gamma*Ic) - (Rc*(mu + epsilon))
    
    ## make a vector of state variables, as list
    dy_SIRS<- c(dS_a, dI_a, dR_a, dS_b, dI_b, dR_b, dS_c, dI_c, dR_c)
    list(dy_SIRS)})
}



##demographic parameters for equal lifespans
mu<- (1/(12*365))
k<- 3000 #carrying capacity = 
b0<- 2/365 #2 bats per year
b1<- b0/k


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
  beta<- sub_params_3sp[i, 4]
  lambda<- sub_params_3sp[i, 7]
  psi_ab<- exp(lambda*(1 - sub_params_3sp[i, 1]))
  psi_ac<- exp(lambda*(1 - sub_params_3sp[i, 2]))
  psi_bc<- exp(lambda*(1 - sub_params_3sp[i, 3]))
  gamma<- sub_params_3sp[i, 5]
  epsilon<- sub_params_3sp[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, 
           gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, b0 = b0, b1 = b1)
  
  
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
SIRSeq_df[1:3651, 15]<- 'short'
SIRSeq_df[3652:7302, 15]<- 'long'


SIRSeq_df$period<- factor(SIRSeq_df$period, levels = c('short', 'long'))


#create variable for period facet
SIRSeq_df$periodFacet[SIRSeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SIRSeq_df$periodFacet[SIRSeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SIRSeq_df$periodFacet<- factor(SIRSeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))

#SIRS plot
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
  beta<- sub_params_3sp[i, 4]
  lambda<- sub_params_3sp[i, 7]
  psi_ab<- exp(lambda*(1 - sub_params_3sp[i, 1]))
  psi_ac<- exp(lambda*(1 - sub_params_3sp[i, 2]))
  psi_bc<- exp(lambda*(1 - sub_params_3sp[i, 3]))
  gamma<- sub_params_3sp[i, 5]
  epsilon<- sub_params_3sp[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, 
           gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, b0 = b0, b1 = b1)
  
  
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
SIRSeq_df[1:3651, 16]<- 'short'
SIRSeq_df[3652:7302, 16]<- 'long'


SIRSeq_df$period<- factor(SIRSeq_df$period, levels = c('short', 'long'))

#create variable for period facet
SIRSeq_df$periodFacet[SIRSeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SIRSeq_df$periodFacet[SIRSeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SIRSeq_df$periodFacet<- factor(SIRSeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))


#SIRS plot
SIRS_3sp_InfTime<- ggplot(data = SIRSeq_df, aes(x = time, y = I)) +
  facet_grid(~periodFacet, labeller = label_parsed) +
  geom_path(aes(color = 'I'), lwd = 1, linetype = 'solid') +
  geom_path(aes(x = time, y = Ia, color = 'Ia'), linetype = 'solid') +
  geom_path(aes(x = time, y = Ib, color = 'Ib'), linetype = 'dashed') +
  geom_path(aes(x = time, y = Ic, color = 'Ic'), linetype = 'dotted') +
  scale_color_manual(name = "Infected Species Population", 
                     values = c("I" = "black", "Ia" = '#C7E020FF', 'Ib' = '#287C8CFF', 'Ic' = '#481F70FF')) +
  th +
  labs(y = 'Total Infected Population (I)', x = 'Time (days)', 
       color = "Infected Species Population", linetype = "Infected Species Population",
       title = 'SIRS 3 species')





#SILI 3 sp disease free
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
    dS_a<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + La + Lb + Lc))*(Sa + Ia + La)) - ((beta*Sa*Ia) + (beta*psi_ab*Sa*Ib) + (beta*psi_ac*Sa*Ic)) - (mu*Sa)  
    
    dI_a<- ((beta*Sa*Ia) + (beta*psi_ab*Sa*Ib) + (beta*psi_ac*Sa*Ic)) - (Ia*(mu + gamma)) + (epsilon*La) 
    
    dL_a<- (gamma*Ia) - (La*(mu + epsilon))
    
    #rates of change for other bat species 1
    dS_b<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + La + Lb + Lc))*(Sb + Ib + Lb)) - ((beta*Sb*Ib) + (beta*psi_ab*Sb*Ia) + (beta*psi_bc*Sb*Ic)) - (mu*Sb)
    
    dI_b<- ((beta*Sb*Ib) + (beta*psi_ab*Sb*Ia) + (beta*psi_bc*Sb*Ic)) - (Ib*(mu + gamma)) + (epsilon*Lb) 
    
    dL_b<- (gamma*Ib) - (Lb*(mu + epsilon))
    
    #rates of change for other bat species 2
    dS_c<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + La + Lb + Lc))*(Sc + Ic + Lc)) - ((beta*Sc*Ic) + (beta*psi_ac*Sc*Ia) + (beta*psi_bc*Sc*Ib)) - (mu*Sc)
    
    dI_c<- ((beta*Sc*Ic) + (beta*psi_ac*Sc*Ia) + (beta*psi_bc*Sc*Ib)) - (Ic*(mu + gamma)) + (epsilon*Lc) 
    
    dL_c<- (gamma*Ic) - (Lc*(mu + epsilon))
    
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
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params_3sp[i, 4]
  lambda<- sub_params_3sp[i, 7]
  psi_ab<- exp(lambda*(1 - sub_params_3sp[i, 1]))
  psi_ac<- exp(lambda*(1 - sub_params_3sp[i, 2]))
  psi_bc<- exp(lambda*(1 - sub_params_3sp[i, 3]))
  gamma<- sub_params_3sp[i, 5]
  epsilon<- sub_params_3sp[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, 
           gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, b0 = b0, b1 = b1)
  
  
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
SILIeq_df[1:3651, 16]<- 'short'
SILIeq_df[3652:7302, 16]<- 'long'


SILIeq_df$period<- factor(SILIeq_df$period, levels = c('short', 'long'))

#create variable for period facet
SILIeq_df$periodFacet[SILIeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SILIeq_df$periodFacet[SILIeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SILIeq_df$periodFacet<- factor(SILIeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))



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
  psi_ab<- exp(lambda*(1 - sub_params_3sp[i, 1]))
  psi_ac<- exp(lambda*(1 - sub_params_3sp[i, 2]))
  psi_bc<- exp(lambda*(1 - sub_params_3sp[i, 3]))
  gamma<- sub_params_3sp[i, 5]
  epsilon<- sub_params_3sp[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, 
           gamma = gamma, epsilon = epsilon, 
           mu = mu, lambda = lambda, k = k, b0 = b0, b1 = b1)
  
  
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
SILIeq_df[1:3651, 16]<- 'short'
SILIeq_df[3652:7302, 16]<- 'long'


SILIeq_df$period<- factor(SILIeq_df$period, levels = c('short', 'long'))

SILIeq_df$periodFacet[SILIeq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
SILIeq_df$periodFacet[SILIeq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'

SILIeq_df$periodFacet<- factor(SILIeq_df$periodFacet, 
                               levels = c('"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"', 
                                          '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'))

#SILI infectious pop example plot
SILI_3sp_InfTime<- ggplot(data = SILIeq_df, aes(x = time, y = I)) +
  facet_grid(~periodFacet, labeller = label_parsed) +
  geom_path(aes(color = 'I'), lwd = 1, linetype = 'solid') +
  geom_path(aes(x = time, y = Ia, color = 'Ia'), linetype = 'solid') +
  geom_path(aes(x = time, y = Ib, color = 'Ib'), linetype = 'dashed') +
  geom_path(aes(x = time, y = Ic, color = 'Ic'), linetype = 'dotted') +
  scale_color_manual(name = "Infected Species Population", values = c("I" = "black", "Ia" = '#C7E020FF', 'Ib' = '#287C8CFF', 'Ic' = '#481F70FF')) +
  th +
  labs(y = 'Total Infected Population (I)', x = 'Time (days)', 
       color = 'Infected Species Population', linetype = 'Infected Species Population',
       title = 'SILI 3 species')





#combine 3 sp plots 
SIRS_3sp_DiseaseFree / SILI_3sp_DiseaseFree + plot_layout(guides = 'collect')

SIRS_3sp_InfTime / SILI_3sp_InfTime + plot_layout(guides = 'collect')













#clear space for 3sp full param space model outputs
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

#read in 3sp model data
load('OSCER_3sp_CST_models_20240916.Rdata')


#fix lowest transmission values in data to not be NA
SIRSeq_params_df_3sp$beta<- replace_na(SIRSeq_params_df_3sp$beta, '0.0005')
SIRSeq_params_df_3sp$beta<- as.character(SIRSeq_params_df_3sp$beta)
SIRSeq_params_df_3sp$beta<- as.numeric(SIRSeq_params_df_3sp$beta)

SILIeq_params_df_3sp$beta<- replace_na(SILIeq_params_df_3sp$beta, '0.0005')
SILIeq_params_df_3sp$beta<- as.character(SILIeq_params_df_3sp$beta)
SILIeq_params_df_3sp$beta<- as.numeric(SILIeq_params_df_3sp$beta)

#cross check  
unique(SIRSeq_params_df_3sp$beta)
unique(SILIeq_params_df_3sp$beta)


#order periods
SIRSeq_params_df_3sp$period<- factor(SIRSeq_params_df_3sp$period, levels = c('short', 'long'))
SILIeq_params_df_3sp$period<- factor(SILIeq_params_df_3sp$period, levels = c('short', 'long'))


#combine SIRS and SILI eq data
#add a model type column
SIRSeq_params_df_3sp$model<- 'SIRS'
SILIeq_params_df_3sp$model<- 'SILI'

#subset for data needed
eq_SIRS<- SIRSeq_params_df_3sp[ , c(11, 13, 17, 20, 26:29)]
eq_SILI<- SILIeq_params_df_3sp[ , c(11, 13, 17, 20, 26:29)]

Eq_df<- as.data.frame(rbind(eq_SIRS, eq_SILI))
#order models as they come in the manuscript
Eq_df$model<- factor(Eq_df$model, levels = c('SIRS', 'SILI'))

#create variable for lambda facet
Eq_df$lambdaFacet[Eq_df$lambda == -1]<- 'lambda*" = -1"'
Eq_df$lambdaFacet[Eq_df$lambda == -5]<- 'lambda*" = -5"'
Eq_df$lambdaFacet[Eq_df$lambda == -10]<- 'lambda*" = -10"'

Eq_df$lambdaFacet<- factor(Eq_df$lambdaFacet, levels = c('lambda*" = -10"', 
                                                         'lambda*" = -5"', 'lambda*" = -1"'))
#create variable for period facet
Eq_df$periodFacet[Eq_df$period == 'short']<- '"short:"~gamma*" = 1/5,"~epsilon*" = 1/60"'
Eq_df$periodFacet[Eq_df$period == 'long']<- '"long:"~gamma*" = 1/548,"~epsilon*" = 1/1643"'


#faceted figs
ggtern(data = Eq_df[Eq_df$beta == 0.0005 & 
                                     Eq_df$period == 'short',], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(periodFacet*model ~ lambdaFacet, labeller = label_parsed)+
  geom_point(size = 3) +
  scale_color_viridis_c(begin = 0.2, end =  0.9, limits = c(0.05, 0.08)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalence', title = expression(~beta~'='~0.0005)) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()
  


ggtern(data = Eq_df[Eq_df$beta == 0.0005 & 
                      Eq_df$period == 'long',], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(periodFacet*model ~ lambdaFacet, labeller = label_parsed)+
  geom_point(size = 3) +
  scale_color_viridis_c(begin = 0.2, end =  0.9, limits = c(0.313, 0.315)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalence', title = expression(~beta~'='~0.0005)) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


ggtern(data = Eq_df[Eq_df$beta == 0.001 & 
                      Eq_df$period == 'short',], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(periodFacet*model ~ lambdaFacet, labeller = label_parsed)+
  geom_point(size = 3) +
  scale_color_viridis_c(begin = 0.2, end =  0.9, limits = c(0.05, 0.08)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalence', title = expression(~beta~'='~0.001)) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


ggtern(data = Eq_df[Eq_df$beta == 0.001 & 
                      Eq_df$period == 'long',], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(periodFacet*model ~ lambdaFacet, labeller = label_parsed)+
  geom_point(size = 3) +
  scale_color_viridis_c(begin = 0.2, end =  0.9, limits = c(0.313, 0.315)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalence', title = expression(~beta~'='~0.001)) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()












