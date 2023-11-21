#---
#title: OSCER_SIRS_20231105_3sp_N1eqN2eqN3
#author: Molly Simonis
#date: 2023-11-05
#---

#turn on packages
library(deSolve)

#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
#graphics.off()##clear graphics
#gc()##free up some RAM


#determine start time of code
start<- Sys.time()



#read in a tables of parameter combinations
params_combo_3sp_short<- read.csv('ParamCombos_v17_3sp_short.csv', header = T, sep = ',')
params_combo_3sp_long<- read.csv('ParamCombos_v17_3sp_long.csv', header = T, sep = ',')
#make sure column names are correct3
names(params_combo_3sp_short)<- c('psi_ab', 'psi_ac', 'psi_bc', 'beta', 'gamma', 'epsilon')
names(params_combo_3sp_long)<- c('psi_ab', 'psi_ac', 'psi_bc', 'beta', 'gamma', 'epsilon')


#####################3 species models########################

#########################3 species SIRS######################
#Create function for the model
SIRS_model<- function(time, current_states, params){
  
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
    dS_a<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + Ra + Rb + Rc))*(Sa + Ia + Ra)) - ((beta*Sa*Ia) + (beta*psi_ab*Sa*Ib) + (beta*psi_ac*Sa*Ic)) - (mu_a*Sa) + (epsilon*Ra)
    
    dI_a<- ((beta*Sa*Ia) + (beta*psi_ab*Sa*Ib) + (beta*psi_ac*Sa*Ic)) - (Ia*(mu_a + gamma))  
    
    dR_a<- (gamma*Ia) - (Ra*(mu_a + epsilon))
    
    #rates of change for other bat species 1
    dS_b<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + Ra + Rb + Rc))*(Sb + Ib + Rb)) - ((beta*Sb*Ib) + (beta*psi_ab*Sb*Ia) + (beta*psi_bc*Sb*Ic)) - (mu_b*Sb) + (epsilon*Rb)
    
    dI_b<- ((beta*Sb*Ib) + (beta*psi_ab*Sb*Ia) + (beta*psi_bc*Sb*Ic)) - (Ib*(mu_b + gamma))  
    
    dR_b<- (gamma*Ib) - (Rb*(mu_b + epsilon))
    
    #rates of change for other bat species 2
    dS_c<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + Ra + Rb + Rc))*(Sc + Ic + Rc)) - ((beta*Sc*Ic) + (beta*psi_ac*Sc*Ia) + (beta*psi_bc*Sc*Ib)) - (mu_c*Sc) + (epsilon*Rc) 
    
    dI_c<- ((beta*Sc*Ic) + (beta*psi_ac*Sc*Ia) + (beta*psi_bc*Sc*Ib)) - (Ic*(mu_c + gamma)) 
    
    dR_c<- (gamma*Ic) - (Rc*(mu_c + epsilon))
    
    ## make a vector of state variables, as list
    dy_SIRS<- c(dS_a, dI_a, dR_a, dS_b, dI_b, dR_b, dS_c, dI_c, dR_c)
    list(dy_SIRS)})
}

##demographic parameters
mu_a<- (1/(6*365)) #1 death every 6 years typically (average lifespan of Desmodus Rotundus)
mu_b<- (1/(9*365)) #1 death every 9 years typically (average lifespan of Phyllostomus Discolor)
mu_c<- (1/(15*365)) #1 death every 15 years typically (average lifespan of many other bats)
k<- 450 #carrying capacity = 
b0<- 2/365 #2 bats per year
b1<- b0/k

##define simulation run time in days (running for 20 years)
tmax<- 365*20
by<- 10 
time<- seq(0, tmax, by = by)


#starting population values
#Na = Nb = Nc
Sa0<- 149 
Ia0<- 1
Ra0<- 0

Sb0<- 149
Ib0<- 1
Rb0<- 0

Sc0<- 149
Ic0<- 1
Rc0<- 0
Y0<- c(Sa0, Ia0, Ra0, Sb0, Ib0, Rb0, Sc0, Ic0, Rc0)

#for short infectious periods and waning immunity durations
#make empty list to put output dataframes into
SIRSeq_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_3sp_short)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_3sp_short[i, 4]
  psi_ab<- exp(-6*(1 - params_combo_3sp_short[i, 1]))
  psi_ac<- exp(-6*(1 - params_combo_3sp_short[i, 2]))
  psi_bc<- exp(-6*(1 - params_combo_3sp_short[i, 3]))
  gamma<- params_combo_3sp_short[i, 5]
  epsilon<- params_combo_3sp_short[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, gamma = gamma, epsilon = epsilon, mu_a = mu_a, mu_b = mu_b, mu_c = mu_c, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", "Sb", "Ib", "Rb", "Sc", "Ic", "Rc")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic) / (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Sc + SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic + SIRS_out$Ra + SIRS_out$Rb + SIRS_out$Rc)
  
  
  #get all equilibrium only data
  SIRSeq_datalist[[i]]<- SIRS_out[SIRS_out$time == 7300, ]
  
}


#make equilibrium data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset to create plots later
SIRSeq_params_df_3spSHORT_N1eqN2eqN3<- data.frame(SIRSeq_df, params_combo_3sp_short)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips<- 1/SIRSeq_params_df_3spSHORT_N1eqN2eqN3$gamma

SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat<- round(SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips)

SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat == 7]<- '7 days infectious'
SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat == 5]<- '5 days infectious'
SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat == 3]<- '3 days infectious'


SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur<- 1/SIRSeq_params_df_3spSHORT_N1eqN2eqN3$epsilon

SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur_cat<- round(SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur)

SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur_cat[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur_cat == 90]<- '90 days immune'
SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur_cat[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur_cat == 60]<- '60 days immune'
SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur_cat[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur_cat == 30]<- '30 days immune'


SIRSeq_params_df_3spSHORT_N1eqN2eqN3$psi_ab_perc<- SIRSeq_params_df_3spSHORT_N1eqN2eqN3$psi_ab*100
SIRSeq_params_df_3spSHORT_N1eqN2eqN3$psi_ac_perc<- SIRSeq_params_df_3spSHORT_N1eqN2eqN3$psi_ac*100
SIRSeq_params_df_3spSHORT_N1eqN2eqN3$psi_bc_perc<- SIRSeq_params_df_3spSHORT_N1eqN2eqN3$psi_bc*100


SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta<- as.character(SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta)
SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta<- factor(SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta, 
                                                   levels = c('0.0005', '0.0025', '0.005', '0.0075'))







##for long infectious periods and waning immunity durations
#make empty list to put output dataframes into
SIRSeq_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_3sp_long)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_3sp_long[i, 4]
  psi_ab<- exp(-6*(1 - params_combo_3sp_long[i, 1]))
  psi_ac<- exp(-6*(1 - params_combo_3sp_long[i, 2]))
  psi_bc<- exp(-6*(1 - params_combo_3sp_long[i, 3]))
  gamma<- params_combo_3sp_long[i, 5]
  epsilon<- params_combo_3sp_long[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, gamma = gamma, epsilon = epsilon, mu_a = mu_a, mu_b = mu_b, mu_c = mu_c, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", "Sb", "Ib", "Rb", "Sc", "Ic", "Rc")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic) / (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Sc + SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic + SIRS_out$Ra + SIRS_out$Rb + SIRS_out$Rc)
  
  
  #get all equilibrium only data
  SIRSeq_datalist[[i]]<- SIRS_out[SIRS_out$time == 7300, ]
  
}


#make equilibrium data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

#Add phylogenetic relatedness parameter values to equilibrium values for a full 
#equilibrium dataset to create plots later
SIRSeq_params_df_3spLONG_N1eqN2eqN3<- data.frame(SIRSeq_df, params_combo_3sp_long)

#For SIRS 3sp long ips and imm_dur Na = Nb
SIRSeq_params_df_3spLONG_N1eqN2eqN3$ips<- 1/SIRSeq_params_df_3spLONG_N1eqN2eqN3$gamma

SIRSeq_params_df_3spLONG_N1eqN2eqN3$ips_cat<- round(SIRSeq_params_df_3spLONG_N1eqN2eqN3$ips)

SIRSeq_params_df_3spLONG_N1eqN2eqN3$ips_cat[SIRSeq_params_df_3spLONG_N1eqN2eqN3$ips_cat == 730]<- '730 days infectious'
SIRSeq_params_df_3spLONG_N1eqN2eqN3$ips_cat[SIRSeq_params_df_3spLONG_N1eqN2eqN3$ips_cat == 548]<- '548 days infectious'
SIRSeq_params_df_3spLONG_N1eqN2eqN3$ips_cat[SIRSeq_params_df_3spLONG_N1eqN2eqN3$ips_cat == 365]<- '365 days infectious'

SIRSeq_params_df_3spLONG_N1eqN2eqN3$imm_dur<- 1/SIRSeq_params_df_3spLONG_N1eqN2eqN3$epsilon

SIRSeq_params_df_3spLONG_N1eqN2eqN3$imm_dur_cat<- round(SIRSeq_params_df_3spLONG_N1eqN2eqN3$imm_dur)

SIRSeq_params_df_3spLONG_N1eqN2eqN3$imm_dur_cat[SIRSeq_params_df_3spLONG_N1eqN2eqN3$imm_dur_cat == 2190]<- '2190 days immune'
SIRSeq_params_df_3spLONG_N1eqN2eqN3$imm_dur_cat[SIRSeq_params_df_3spLONG_N1eqN2eqN3$imm_dur_cat == 1643]<- '1643 days immune'
SIRSeq_params_df_3spLONG_N1eqN2eqN3$imm_dur_cat[SIRSeq_params_df_3spLONG_N1eqN2eqN3$imm_dur_cat == 1095]<- '1095 days immune'


SIRSeq_params_df_3spLONG_N1eqN2eqN3$psi_ab_perc<- SIRSeq_params_df_3spLONG_N1eqN2eqN3$psi_ab*100
SIRSeq_params_df_3spLONG_N1eqN2eqN3$psi_ac_perc<- SIRSeq_params_df_3spLONG_N1eqN2eqN3$psi_ac*100
SIRSeq_params_df_3spLONG_N1eqN2eqN3$psi_bc_perc<- SIRSeq_params_df_3spLONG_N1eqN2eqN3$psi_bc*100


SIRSeq_params_df_3spLONG_N1eqN2eqN3$beta<- as.character(SIRSeq_params_df_3spLONG_N1eqN2eqN3$beta)
SIRSeq_params_df_3spLONG_N1eqN2eqN3$beta<- factor(SIRSeq_params_df_3spLONG_N1eqN2eqN3$beta, 
                                                  levels = c('0.0005', '0.0025', '0.005', '0.0075'))




#Save environment for plotting later
save.image(paste(date(), '_SIRS3spN1eqN2eqN3eqN3.Rdata'))

#determine amount of time code took to run
Sys.time() - start