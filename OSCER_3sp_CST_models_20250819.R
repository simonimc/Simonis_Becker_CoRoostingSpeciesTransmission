#---
#title: OSCER_3sp_CST_models_20250819
#author: Molly Simonis
#date: 2025-08-19
#---

#turn on packages
library(deSolve)


#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM


#read in a tables of parameter combinations 3sp
params_combo_3sp<- read.csv('ParamCombos_v19_3sp.csv', header = T, sep = ',')
#make sure column names are correct3
names(params_combo_3sp)<- c('psi_ab', 'psi_ac', 'psi_bc', 'beta', 'gamma', 'epsilon', 'lambda', 'period')



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
SIRSeq_datalist<- list() #for equilibria data
SIRS_out_tail<- list() #for checking data reached equilibria


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_3sp)) {
  
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
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic) / (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Sc + SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic + SIRS_out$Ra + SIRS_out$Rb + SIRS_out$Rc)
  
  SIRS_out$N<- SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Sc + 
    SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic + 
    SIRS_out$Ra + SIRS_out$Rb + SIRS_out$Rc
  
  SIRS_out$period<- params_combo_3sp[i, 8]
  
  #get all equilibrium only data
  SIRSeq_datalist[[i]]<- tail(SIRS_out, 1)
  
  ##check that final equillibrium timesteps gives a difference of <0.001 individuals within the total pop
  SIRS_out_tail[[i]]<- tail(SIRS_out$N)[2] - tail(SIRS_out$N)[1]
  
  print(i)
}



#check all parameter combos made it to equilibrium
SIRS_out_tail_df<- do.call(rbind, SIRS_out_tail)

#make equilibrium data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset to create plots later
SIRSeq_params_df_3sp<- data.frame(SIRSeq_df, params_combo_3sp)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SIRSeq_params_df_3sp$ips<- 1/SIRSeq_params_df_3sp$gamma

SIRSeq_params_df_3sp$ips_cat<- round(SIRSeq_params_df_3sp$ips)

SIRSeq_params_df_3sp$ips_cat[SIRSeq_params_df_3sp$ips_cat == 5]<- '5 days infectious'
SIRSeq_params_df_3sp$ips_cat[SIRSeq_params_df_3sp$ips_cat == 548]<- '548 days infectious'


SIRSeq_params_df_3sp$imm_dur<- 1/SIRSeq_params_df_3sp$epsilon

SIRSeq_params_df_3sp$imm_dur_cat<- round(SIRSeq_params_df_3sp$imm_dur)

SIRSeq_params_df_3sp$imm_dur_cat[SIRSeq_params_df_3sp$imm_dur_cat == 60]<- '60 days immune'
SIRSeq_params_df_3sp$imm_dur_cat[SIRSeq_params_df_3sp$imm_dur_cat == 1643]<- '1643 days immune'


SIRSeq_params_df_3sp$psi_ab_perc<- SIRSeq_params_df_3sp$psi_ab*100
SIRSeq_params_df_3sp$psi_ac_perc<- SIRSeq_params_df_3sp$psi_ac*100
SIRSeq_params_df_3sp$psi_bc_perc<- SIRSeq_params_df_3sp$psi_bc*100


SIRSeq_params_df_3sp$beta<- as.character(SIRSeq_params_df_3sp$beta)
SIRSeq_params_df_3sp$beta<- factor(SIRSeq_params_df_3sp$beta, 
                                   levels = c('0.0001', '0.0005', '0.001'))






#rename epsilon to omega for SILI parameter space
names(params_combo_3sp)<- c('psi_ab', 'psi_ac', 'psi_bc', 'beta', 'gamma', 'omega', 'lambda', 'period')

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


#all demographic parameters and time series the same

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



#make empty lists to put output dataframes into
SILIeq_datalist<- list() #for equilibria data
SILI_out_tail<- list() #for checking data reached equilibria


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_3sp)) {
  
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
  omega<- params_combo_3sp[i, 6]
  
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
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib + SILI_out$Ic) / 
    (SILI_out$Sa + SILI_out$Sb + SILI_out$Sc + 
       SILI_out$Ia + SILI_out$Ib + SILI_out$Ic + 
       SILI_out$La + SILI_out$Lb + SILI_out$Lc)
  
  SILI_out$N<- SILI_out$Sa + SILI_out$Sb + SILI_out$Sc + 
    SILI_out$Ia + SILI_out$Ib + SILI_out$Ic + 
    SILI_out$La + SILI_out$Lb + SILI_out$Lc
  
  SILI_out$period<- params_combo_3sp[i, 8]
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- tail(SILI_out, 1)
  
  ##check that final equillibrium timesteps gives a difference of <0.001 individuals within the total pop
  SILI_out_tail[[i]]<- tail(SILI_out$N)[2] - tail(SILI_out$N)[1]
  
  print(i)
}



#check all parameter combos made it to equilibrium
SILI_out_tail_df<- do.call(rbind, SILI_out_tail)

#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset to create plots later
SILIeq_params_df_3sp<- data.frame(SILIeq_df, params_combo_3sp)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SILIeq_params_df_3sp$ips<- 1/SILIeq_params_df_3sp$gamma

SILIeq_params_df_3sp$ips_cat<- round(SILIeq_params_df_3sp$ips)

SILIeq_params_df_3sp$ips_cat[SILIeq_params_df_3sp$ips_cat == 5]<- '5 days infectious'
SILIeq_params_df_3sp$ips_cat[SILIeq_params_df_3sp$ips_cat == 548]<- '548 days infectious'


SILIeq_params_df_3sp$lat_dur<- 1/SILIeq_params_df_3sp$omega

SILIeq_params_df_3sp$lat_dur_cat<- round(SILIeq_params_df_3sp$lat_dur)

SILIeq_params_df_3sp$lat_dur_cat[SILIeq_params_df_3sp$lat_dur_cat == 60]<- '60 days latent'
SILIeq_params_df_3sp$lat_dur_cat[SILIeq_params_df_3sp$lat_dur_cat == 1643]<- '1643 days latent'


SILIeq_params_df_3sp$psi_ab_perc<- SILIeq_params_df_3sp$psi_ab*100
SILIeq_params_df_3sp$psi_ac_perc<- SILIeq_params_df_3sp$psi_ac*100
SILIeq_params_df_3sp$psi_bc_perc<- SILIeq_params_df_3sp$psi_bc*100


SILIeq_params_df_3sp$beta<- as.character(SILIeq_params_df_3sp$beta)
SILIeq_params_df_3sp$beta<- factor(SILIeq_params_df_3sp$beta, 
                                   levels = c('0.0001', '0.0005', '0.001'))


save.image('OSCER_3sp_CST_models_20250820.Rdata')

