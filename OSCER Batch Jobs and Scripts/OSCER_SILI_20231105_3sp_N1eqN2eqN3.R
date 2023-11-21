#---
#title: OSCER_SILI_20231105_3sp_N1eqN2eqN3
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

#########################3 species SILI######################
#Create function for the model
SILI_model<- function(time, current_states, params){
  
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
    dS_a<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + La + Lb + Lc))*(Sa + Ia + La)) - ((beta*Sa*Ia) + (beta*psi_ab*Sa*Ib) + (beta*psi_ac*Sa*Ic)) - (mu_a*Sa)  
    
    dI_a<- ((beta*Sa*Ia) + (beta*psi_ab*Sa*Ib) + (beta*psi_ac*Sa*Ic)) - (Ia*(mu_a + gamma)) + (epsilon*La) 
    
    dL_a<- (gamma*Ia) - (La*(mu_a + epsilon))
    
    #rates of change for other bat species 1
    dS_b<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + La + Lb + Lc))*(Sb + Ib + Lb)) - ((beta*Sb*Ib) + (beta*psi_ab*Sb*Ia) + (beta*psi_bc*Sb*Ic)) - (mu_b*Sb)
    
    dI_b<- ((beta*Sb*Ib) + (beta*psi_ab*Sb*Ia) + (beta*psi_bc*Sb*Ic)) - (Ib*(mu_b + gamma)) + (epsilon*Lb) 
  
    dL_b<- (gamma*Ib) - (Lb*(mu_b + epsilon))
    
    #rates of change for other bat species 2
    dS_c<- ((b0 - b1*(Sa + Sb + Sc +Ia + Ib + Ic + La + Lb + Lc))*(Sc + Ic + Lc)) - ((beta*Sc*Ic) + (beta*psi_ac*Sc*Ia) + (beta*psi_bc*Sc*Ib)) - (mu_c*Sc)
    
    dI_c<- ((beta*Sc*Ic) + (beta*psi_ac*Sc*Ia) + (beta*psi_bc*Sc*Ib)) - (Ic*(mu_c + gamma)) + (epsilon*Lc) 
    
    dL_c<- (gamma*Ic) - (Lc*(mu_c + epsilon))
    
    ## make a vector of state variables, as list
    dy_SILI<- c(dS_a, dI_a, dL_a, dS_b, dI_b, dL_b, dS_c, dI_c, dL_c)
    list(dy_SILI)})
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
La0<- 0

Sb0<- 149
Ib0<- 1
Lb0<- 0

Sc0<- 149
Ic0<- 1
Lc0<- 0
Y0<- c(Sa0, Ia0, La0, Sb0, Ib0, Lb0, Sc0, Ic0, Lc0)

#for short infectious periods and latency durations
#make empty list to put output dataframes into
SILIeq_datalist<- list()


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
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", "Sb", "Ib", "Lb", "Sc", "Ic", "Lc")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib + SILI_out$Ic) / (SILI_out$Sa + SILI_out$Sb + SILI_out$Sc + SILI_out$Ia + SILI_out$Ib + SILI_out$Ic + SILI_out$La + SILI_out$Lb + SILI_out$Lc)
  
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- SILI_out[SILI_out$time == 7300, ]
  
}


#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset to create plots later
SILIeq_params_df_3spSHORT_N1eqN2eqN3<- data.frame(SILIeq_df, params_combo_3sp_short)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SILIeq_params_df_3spSHORT_N1eqN2eqN3$ips<- 1/SILIeq_params_df_3spSHORT_N1eqN2eqN3$gamma

SILIeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat<- round(SILIeq_params_df_3spSHORT_N1eqN2eqN3$ips)

SILIeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat[SILIeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat == 7]<- '7 days infectious'
SILIeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat[SILIeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat == 5]<- '5 days infectious'
SILIeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat[SILIeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat == 3]<- '3 days infectious'


SILIeq_params_df_3spSHORT_N1eqN2eqN3$lat_dur<- 1/SILIeq_params_df_3spSHORT_N1eqN2eqN3$epsilon

SILIeq_params_df_3spSHORT_N1eqN2eqN3$lat_dur_cat<- round(SILIeq_params_df_3spSHORT_N1eqN2eqN3$lat_dur)

SILIeq_params_df_3spSHORT_N1eqN2eqN3$lat_dur_cat[SILIeq_params_df_3spSHORT_N1eqN2eqN3$lat_dur_cat == 90]<- '90 days latent'
SILIeq_params_df_3spSHORT_N1eqN2eqN3$lat_dur_cat[SILIeq_params_df_3spSHORT_N1eqN2eqN3$lat_dur_cat == 60]<- '60 days latent'
SILIeq_params_df_3spSHORT_N1eqN2eqN3$lat_dur_cat[SILIeq_params_df_3spSHORT_N1eqN2eqN3$lat_dur_cat == 30]<- '30 days latent'


SILIeq_params_df_3spSHORT_N1eqN2eqN3$psi_ab_perc<- SILIeq_params_df_3spSHORT_N1eqN2eqN3$psi_ab*100
SILIeq_params_df_3spSHORT_N1eqN2eqN3$psi_ac_perc<- SILIeq_params_df_3spSHORT_N1eqN2eqN3$psi_ac*100
SILIeq_params_df_3spSHORT_N1eqN2eqN3$psi_bc_perc<- SILIeq_params_df_3spSHORT_N1eqN2eqN3$psi_bc*100


SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta<- as.character(SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta)
SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta<- factor(SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta, 
                                               levels = c('0.0005', '0.0025', '0.005', '0.0075'))







##for long infectious periods and latency durations
#make empty list to put output dataframes into
SILIeq_datalist<- list()


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
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", "Sb", "Ib", "Lb", "Sc", "Ic", "Lc")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib + SILI_out$Ic) / (SILI_out$Sa + SILI_out$Sb + SILI_out$Sc + SILI_out$Ia + SILI_out$Ib + SILI_out$Ic + SILI_out$La + SILI_out$Lb + SILI_out$Lc)
  
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- SILI_out[SILI_out$time == 7300, ]
  
}


#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)

#Add phylogenetic relatedness parameter values to equilibrium values for a full 
#equilibrium dataset to create plots later
SILIeq_params_df_3spLONG_N1eqN2eqN3<- data.frame(SILIeq_df, params_combo_3sp_long)

#For SILI 3sp long ips and lat_dur Na = Nb
SILIeq_params_df_3spLONG_N1eqN2eqN3$ips<- 1/SILIeq_params_df_3spLONG_N1eqN2eqN3$gamma

SILIeq_params_df_3spLONG_N1eqN2eqN3$ips_cat<- round(SILIeq_params_df_3spLONG_N1eqN2eqN3$ips)

SILIeq_params_df_3spLONG_N1eqN2eqN3$ips_cat[SILIeq_params_df_3spLONG_N1eqN2eqN3$ips_cat == 730]<- '730 days infectious'
SILIeq_params_df_3spLONG_N1eqN2eqN3$ips_cat[SILIeq_params_df_3spLONG_N1eqN2eqN3$ips_cat == 548]<- '548 days infectious'
SILIeq_params_df_3spLONG_N1eqN2eqN3$ips_cat[SILIeq_params_df_3spLONG_N1eqN2eqN3$ips_cat == 365]<- '365 days infectious'

SILIeq_params_df_3spLONG_N1eqN2eqN3$lat_dur<- 1/SILIeq_params_df_3spLONG_N1eqN2eqN3$epsilon

SILIeq_params_df_3spLONG_N1eqN2eqN3$lat_dur_cat<- round(SILIeq_params_df_3spLONG_N1eqN2eqN3$lat_dur)

SILIeq_params_df_3spLONG_N1eqN2eqN3$lat_dur_cat[SILIeq_params_df_3spLONG_N1eqN2eqN3$lat_dur_cat == 2190]<- '2190 days latent'
SILIeq_params_df_3spLONG_N1eqN2eqN3$lat_dur_cat[SILIeq_params_df_3spLONG_N1eqN2eqN3$lat_dur_cat == 1643]<- '1643 days latent'
SILIeq_params_df_3spLONG_N1eqN2eqN3$lat_dur_cat[SILIeq_params_df_3spLONG_N1eqN2eqN3$lat_dur_cat == 1095]<- '1095 days latent'


SILIeq_params_df_3spLONG_N1eqN2eqN3$psi_ab_perc<- SILIeq_params_df_3spLONG_N1eqN2eqN3$psi_ab*100
SILIeq_params_df_3spLONG_N1eqN2eqN3$psi_ac_perc<- SILIeq_params_df_3spLONG_N1eqN2eqN3$psi_ac*100
SILIeq_params_df_3spLONG_N1eqN2eqN3$psi_bc_perc<- SILIeq_params_df_3spLONG_N1eqN2eqN3$psi_bc*100


SILIeq_params_df_3spLONG_N1eqN2eqN3$beta<- as.character(SILIeq_params_df_3spLONG_N1eqN2eqN3$beta)
SILIeq_params_df_3spLONG_N1eqN2eqN3$beta<- factor(SILIeq_params_df_3spLONG_N1eqN2eqN3$beta, 
                                              levels = c('0.0005', '0.0025', '0.005', '0.0075'))




#Save environment for plotting later
save.image(paste(date(), '_SILI3spN1eqN2eqN3eqN3.Rdata'))

#determine amount of time code took to run
Sys.time() - start