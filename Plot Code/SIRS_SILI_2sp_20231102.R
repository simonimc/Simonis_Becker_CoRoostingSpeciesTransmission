#---
#title: SIRS_SILI_20231102
#author: Molly Simonis
#date: 2023-11-02
#---

#turn on packages
library(deSolve)
library(ggplot2)
library(patchwork)


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

#read in a tables of parameter combinations
params_combo_2sp_short<- read.csv('ParamCombos_v17_2sp_short.csv', header = T, sep = ',')
params_combo_2sp_long<- read.csv('ParamCombos_v17_2sp_long.csv', header = T, sep = ',')
#make sure column names are correct3
names(params_combo_2sp_short)<- c('beta', 'psi', 'gamma', 'epsilon')
names(params_combo_2sp_long)<- c('beta', 'psi', 'gamma', 'epsilon')


####################2 Species Models######################

################2 species Models SIRS#######################
#Create function for two species SIRS model
SIRS_model<- function(time, current_states, params){
  
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
    dS_a<- ((b0 - b1*(Sa + Sb + Ia + Ib + Ra + Rb))*(Sa + Ia + Ra)) - ((beta*Sa*Ia) + (beta*psi*Sa*Ib)) - (mu_a*Sa) + (epsilon*Ra)
    
    dI_a<- ((beta*Sa*Ia) + (beta*psi*Sa*Ib)) - (Ia*(mu_a + gamma)) 
    
    dR_a<- (gamma*Ia) - (Ra*(mu_a + epsilon))
    
    #rates of change for other species 1
    dS_b<-((b0 - b1*(Sa + Sb + Ia + Ib + Ra + Rb))*(Sb + Ib + Rb)) - ((beta*Sb*Ib) + (beta*psi*Sb*Ia)) - (mu_b*Sb) + (epsilon*Rb)
    
    dI_b<-((beta*Sb*Ib) + (beta*psi*Sb*Ia)) - (Ib*(mu_b + gamma)) 
    
    dR_b<- (gamma*Ib) - (Rb*(mu_b + epsilon))
    
    
    ## make a vector of state variables, as list
    dy_SIRS<- c(dS_a, dI_a, dR_a, dS_b, dI_b, dR_b)
    list(dy_SIRS)})
}

##demographic parameters that will remain the same
mu_a<- (1/(6*365)) #1 death every 6 years typically (average lifespan of Desmodus Rotundus)
mu_b<- (1/(9*365)) #1 death every 9 years typically (average lifespan of Phyllostomus Discolor)
k<- 450 #carrying capacity 
b0<- 2/365 #2 bats per year
b1<- b0/k

##define simulation run time in days (running for 20 years)
tmax<- 365*20
by<- 10
time<- seq(0, tmax, by = by)


## starting population values (Na = Nb)
Sa0<- 149 
Ia0<- 1
Ra0<- 0

Sb0<- 149
Ib0<- 1
Rb0<- 0

Y0<- c(Sa0, Ia0, Ra0, Sb0, Ib0, Rb0)



#for short infectious periods and waning immunity durations
#make empty list to put output dataframes into
SIRSeq_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_2sp_short)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp_short[i, 1]
  psi<- exp(-6*(1 - params_combo_2sp_short[i, 2]))
  gamma<- params_combo_2sp_short[i, 3]
  epsilon<- params_combo_2sp_short[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu_a = mu_a, mu_b = mu_b, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", 
                      "Sb", "Ib", "Rb")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib) / 
    (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
       SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb)
  
  
  #get all equilibrium only data
  SIRSeq_datalist[[i]]<- SIRS_out[SIRS_out$time == 7300, ]
  
}


#make equilibrium data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset 
SIRSeq_params_df_2spSHORT_N1eqN2<- data.frame(SIRSeq_df, params_combo_2sp_short)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SIRSeq_params_df_2spSHORT_N1eqN2$ips<- 1/SIRSeq_params_df_2spSHORT_N1eqN2$gamma

SIRSeq_params_df_2spSHORT_N1eqN2$ips_cat<- round(SIRSeq_params_df_2spSHORT_N1eqN2$ips)

SIRSeq_params_df_2spSHORT_N1eqN2$ips_cat[SIRSeq_params_df_2spSHORT_N1eqN2$ips_cat == 7]<- '7 days infectious'
SIRSeq_params_df_2spSHORT_N1eqN2$ips_cat[SIRSeq_params_df_2spSHORT_N1eqN2$ips_cat == 5]<- '5 days infectious'
SIRSeq_params_df_2spSHORT_N1eqN2$ips_cat[SIRSeq_params_df_2spSHORT_N1eqN2$ips_cat == 3]<- '3 days infectious'


SIRSeq_params_df_2spSHORT_N1eqN2$imm_dur<- 1/SIRSeq_params_df_2spSHORT_N1eqN2$epsilon

SIRSeq_params_df_2spSHORT_N1eqN2$imm_dur_cat<- round(SIRSeq_params_df_2spSHORT_N1eqN2$imm_dur)

SIRSeq_params_df_2spSHORT_N1eqN2$imm_dur_cat[SIRSeq_params_df_2spSHORT_N1eqN2$imm_dur_cat == 90]<- '90 days immune'
SIRSeq_params_df_2spSHORT_N1eqN2$imm_dur_cat[SIRSeq_params_df_2spSHORT_N1eqN2$imm_dur_cat == 60]<- '60 days immune'
SIRSeq_params_df_2spSHORT_N1eqN2$imm_dur_cat[SIRSeq_params_df_2spSHORT_N1eqN2$imm_dur_cat == 30]<- '30 days immune'


SIRSeq_params_df_2spSHORT_N1eqN2$psi_perc<- SIRSeq_params_df_2spSHORT_N1eqN2$psi*100


SIRSeq_params_df_2spSHORT_N1eqN2$beta<- as.character(SIRSeq_params_df_2spSHORT_N1eqN2$beta)
SIRSeq_params_df_2spSHORT_N1eqN2$beta<- factor(SIRSeq_params_df_2spSHORT_N1eqN2$beta, 
                                                 levels = c('0.0005', '0.0025', '0.005', '0.0075'))

#Supplemental Figure S6
ggplot(data = SIRSeq_params_df_2spSHORT_N1eqN2, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.9) +  
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = expression(SIRS~N[A]~'='~N[B])) +
  ylim(-0.01, 0.5)




#for long infectious periods and waning immunity durations
#make empty list to put output dataframes into
SIRSeq_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_2sp_long)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp_long[i, 1]
  psi<- exp(-6*(1 - params_combo_2sp_long[i, 2]))
  gamma<- params_combo_2sp_long[i, 3]
  epsilon<- params_combo_2sp_long[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu_a = mu_a, mu_b = mu_b, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", 
                      "Sb", "Ib", "Rb")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ia) / 
    (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
       SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb)
  
  
  #get all equilibrium only data
  SIRSeq_datalist[[i]]<- SIRS_out[SIRS_out$time == 7300, ]
  
}


#make equilibrium data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

#Add phylogenetic relatedness parameter values to equilibrium values for a full 
#equilibrium dataset to create plots later
SIRSeq_params_df_2spLONG_N1eqN2<- data.frame(SIRSeq_df, params_combo_2sp_long)

#For SIRS 2sp long ips and imm_dur N1 = N2
SIRSeq_params_df_2spLONG_N1eqN2$ips<- 1/SIRSeq_params_df_2spLONG_N1eqN2$gamma

SIRSeq_params_df_2spLONG_N1eqN2$ips_cat<- round(SIRSeq_params_df_2spLONG_N1eqN2$ips)

SIRSeq_params_df_2spLONG_N1eqN2$ips_cat[SIRSeq_params_df_2spLONG_N1eqN2$ips_cat == 730]<- '730 days infectious'
SIRSeq_params_df_2spLONG_N1eqN2$ips_cat[SIRSeq_params_df_2spLONG_N1eqN2$ips_cat == 548]<- '548 days infectious'
SIRSeq_params_df_2spLONG_N1eqN2$ips_cat[SIRSeq_params_df_2spLONG_N1eqN2$ips_cat == 365]<- '365 days infectious'

SIRSeq_params_df_2spLONG_N1eqN2$imm_dur<- 1/SIRSeq_params_df_2spLONG_N1eqN2$epsilon

SIRSeq_params_df_2spLONG_N1eqN2$imm_dur_cat<- round(SIRSeq_params_df_2spLONG_N1eqN2$imm_dur)

SIRSeq_params_df_2spLONG_N1eqN2$imm_dur_cat[SIRSeq_params_df_2spLONG_N1eqN2$imm_dur_cat == 2190]<- '2190 days immune'
SIRSeq_params_df_2spLONG_N1eqN2$imm_dur_cat[SIRSeq_params_df_2spLONG_N1eqN2$imm_dur_cat == 1643]<- '1643 days immune'
SIRSeq_params_df_2spLONG_N1eqN2$imm_dur_cat[SIRSeq_params_df_2spLONG_N1eqN2$imm_dur_cat == 1095]<- '1095 days immune'


SIRSeq_params_df_2spLONG_N1eqN2$psi_perc<- SIRSeq_params_df_2spLONG_N1eqN2$psi*100


SIRSeq_params_df_2spLONG_N1eqN2$beta<- as.character(SIRSeq_params_df_2spLONG_N1eqN2$beta)
SIRSeq_params_df_2spLONG_N1eqN2$beta<- factor(SIRSeq_params_df_2spLONG_N1eqN2$beta, 
                                                levels = c('0.0005', '0.0025', '0.005', '0.0075'))

#Supplemental Figure S7
ggplot(data = SIRSeq_params_df_2spLONG_N1eqN2, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = expression(SIRS~N[A]~'='~N[B])) +
  ylim(-0.01, 0.5)









## starting population values (Na > Nb)
Sa0<- 149 
Ia0<- 1
Ra0<- 0

Sb0<- 74
Ib0<- 1
Rb0<- 0

Y0<- c(Sa0, Ia0, Ra0, Sb0, Ib0, Rb0)


#for short infectious periods and waning immunity durations
#make empty list to put output dataframes into
SIRSeq_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_2sp_short)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp_short[i, 1]
  psi<- exp(-6*(1 - params_combo_2sp_short[i, 2]))
  gamma<- params_combo_2sp_short[i, 3]
  epsilon<- params_combo_2sp_short[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu_a = mu_a, mu_b = mu_b, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", 
                      "Sb", "Ib", "Rb")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib) / 
    (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
       SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb)
  
  
  #get all equilibrium only data
  SIRSeq_datalist[[i]]<- SIRS_out[SIRS_out$time == 7300, ]
  
}


#make equilibrium data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset to create plots later
SIRSeq_params_df_2spSHORT_N1grN2<- data.frame(SIRSeq_df, params_combo_2sp_short)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SIRSeq_params_df_2spSHORT_N1grN2$ips<- 1/SIRSeq_params_df_2spSHORT_N1grN2$gamma

SIRSeq_params_df_2spSHORT_N1grN2$ips_cat<- round(SIRSeq_params_df_2spSHORT_N1grN2$ips)

SIRSeq_params_df_2spSHORT_N1grN2$ips_cat[SIRSeq_params_df_2spSHORT_N1grN2$ips_cat == 7]<- '7 days infectious'
SIRSeq_params_df_2spSHORT_N1grN2$ips_cat[SIRSeq_params_df_2spSHORT_N1grN2$ips_cat == 5]<- '5 days infectious'
SIRSeq_params_df_2spSHORT_N1grN2$ips_cat[SIRSeq_params_df_2spSHORT_N1grN2$ips_cat == 3]<- '3 days infectious'


SIRSeq_params_df_2spSHORT_N1grN2$imm_dur<- 1/SIRSeq_params_df_2spSHORT_N1grN2$epsilon

SIRSeq_params_df_2spSHORT_N1grN2$imm_dur_cat<- round(SIRSeq_params_df_2spSHORT_N1grN2$imm_dur)

SIRSeq_params_df_2spSHORT_N1grN2$imm_dur_cat[SIRSeq_params_df_2spSHORT_N1grN2$imm_dur_cat == 90]<- '90 days immune'
SIRSeq_params_df_2spSHORT_N1grN2$imm_dur_cat[SIRSeq_params_df_2spSHORT_N1grN2$imm_dur_cat == 60]<- '60 days immune'
SIRSeq_params_df_2spSHORT_N1grN2$imm_dur_cat[SIRSeq_params_df_2spSHORT_N1grN2$imm_dur_cat == 30]<- '30 days immune'


SIRSeq_params_df_2spSHORT_N1grN2$psi_perc<- SIRSeq_params_df_2spSHORT_N1grN2$psi*100


SIRSeq_params_df_2spSHORT_N1grN2$beta<- as.character(SIRSeq_params_df_2spSHORT_N1grN2$beta)
SIRSeq_params_df_2spSHORT_N1grN2$beta<- factor(SIRSeq_params_df_2spSHORT_N1grN2$beta, 
                                               levels = c('0.0005', '0.0025', '0.005', '0.0075'))

#Supplemental Figure S8
ggplot(data = SIRSeq_params_df_2spSHORT_N1grN2, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.9) +  
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = expression(SIRS~N[A]~'>'~N[B])) +
  ylim(-0.01, 0.5)



#for long infectious periods and waning immunity durations
#make empty list to put output dataframes into
SIRSeq_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_2sp_long)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp_long[i, 1]
  psi<- exp(-6*(1 - params_combo_2sp_long[i, 2]))
  gamma<- params_combo_2sp_long[i, 3]
  epsilon<- params_combo_2sp_long[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu_a = mu_a, mu_b = mu_b, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", 
                      "Sb", "Ib", "Rb")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ia) / 
    (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Ia + 
       SIRS_out$Ib + SIRS_out$Ra + SIRS_out$Rb)
  
  
  #get all equilibrium only data
  SIRSeq_datalist[[i]]<- SIRS_out[SIRS_out$time == 7300, ]
  
}

#make equilibrium data one df
SIRSeq_df<- do.call(rbind, SIRSeq_datalist)

#Add phylogenetic relatedness parameter values to equilibrium values for a full 
#equilibrium dataset to create plots later
SIRSeq_params_df_2spLONG_N1grN2<- data.frame(SIRSeq_df, params_combo_2sp_long)

#For SIRS 2sp long ips and imm_dur N1 = N2
SIRSeq_params_df_2spLONG_N1grN2$ips<- 1/SIRSeq_params_df_2spLONG_N1grN2$gamma

SIRSeq_params_df_2spLONG_N1grN2$ips_cat<- round(SIRSeq_params_df_2spLONG_N1grN2$ips)

SIRSeq_params_df_2spLONG_N1grN2$ips_cat[SIRSeq_params_df_2spLONG_N1grN2$ips_cat == 730]<- '730 days infectious'
SIRSeq_params_df_2spLONG_N1grN2$ips_cat[SIRSeq_params_df_2spLONG_N1grN2$ips_cat == 548]<- '548 days infectious'
SIRSeq_params_df_2spLONG_N1grN2$ips_cat[SIRSeq_params_df_2spLONG_N1grN2$ips_cat == 365]<- '365 days infectious'

SIRSeq_params_df_2spLONG_N1grN2$imm_dur<- 1/SIRSeq_params_df_2spLONG_N1grN2$epsilon

SIRSeq_params_df_2spLONG_N1grN2$imm_dur_cat<- round(SIRSeq_params_df_2spLONG_N1grN2$imm_dur)

SIRSeq_params_df_2spLONG_N1grN2$imm_dur_cat[SIRSeq_params_df_2spLONG_N1grN2$imm_dur_cat == 2190]<- '2190 days immune'
SIRSeq_params_df_2spLONG_N1grN2$imm_dur_cat[SIRSeq_params_df_2spLONG_N1grN2$imm_dur_cat == 1643]<- '1643 days immune'
SIRSeq_params_df_2spLONG_N1grN2$imm_dur_cat[SIRSeq_params_df_2spLONG_N1grN2$imm_dur_cat == 1095]<- '1095 days immune'


SIRSeq_params_df_2spLONG_N1grN2$psi_perc<- SIRSeq_params_df_2spLONG_N1grN2$psi*100


SIRSeq_params_df_2spLONG_N1grN2$beta<- as.character(SIRSeq_params_df_2spLONG_N1grN2$beta)
SIRSeq_params_df_2spLONG_N1grN2$beta<- factor(SIRSeq_params_df_2spLONG_N1grN2$beta, 
                                              levels = c('0.0005', '0.0025', '0.005', '0.0075'))

#Supplemental Figure S9
ggplot(data = SIRSeq_params_df_2spLONG_N1grN2, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = expression(SIRS~N[A]~'>'~N[B])) +
  ylim(-0.01, 0.5)










################2 species Models SILI#######################
#Create function for two species SILI model
SILI_model<- function(time, current_states, params){
  
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
    dS_a<- ((b0 - b1*(Sa + Sb + Ia + Ib + La + Lb))*(Sa + Ia + La)) - ((beta*Sa*Ia) + (beta*psi*Sa*Ib)) - (mu_a*Sa) 
    
    dI_a<- ((beta*Sa*Ia) + (beta*psi*Sa*Ib)) - (Ia*(mu_a + gamma)) + (epsilon*La)
    
    dL_a<- (gamma*Ia) - (La*(mu_a + epsilon))
    
    #rates of change for other species 1
    dS_b<-((b0 - b1*(Sa + Sb + Ia + Ib + La + Lb))*(Sb + Ib + Lb)) - ((beta*Sb*Ib) + (beta*psi*Sb*Ia)) - (mu_b*Sb) 
    
    dI_b<-((beta*Sb*Ib) + (beta*psi*Sb*Ia)) - (Ib*(mu_b + gamma)) + (epsilon*Lb)
    
    dL_b<- (gamma*Ib) - (Lb*(mu_b + epsilon))
    
    
    ## make a vector of state variables, as list
    dy_SILI<- c(dS_a, dI_a, dL_a, dS_b, dI_b, dL_b)
    list(dy_SILI)})
}

##demographic parameters that will remain the same

##simulation run time in days will remain the same


## starting population values (Na = Nb)
Sa0<- 149 
Ia0<- 1
La0<- 0

Sb0<- 149
Ib0<- 1
Lb0<- 0

Y0<- c(Sa0, Ia0, La0, Sb0, Ib0, Lb0)


#for short infectious periods and latency durations
#make empty list to put output dataframes into
SILIeq_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_2sp_short)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp_short[i, 1]
  psi<- exp(-6*(1 - params_combo_2sp_short[i, 2]))
  gamma<- params_combo_2sp_short[i, 3]
  epsilon<- params_combo_2sp_short[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu_a = mu_a, mu_b = mu_b, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", 
                      "Sb", "Ib", "Lb")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib) / 
    (SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
       SILI_out$Ib + SILI_out$La + SILI_out$Lb)
  
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- SILI_out[SILI_out$time == 7300, ]
  
}


#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset to create plots later
SILIeq_params_df_2spSHORT_N1eqN2<- data.frame(SILIeq_df, params_combo_2sp_short)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SILIeq_params_df_2spSHORT_N1eqN2$ips<- 1/SILIeq_params_df_2spSHORT_N1eqN2$gamma

SILIeq_params_df_2spSHORT_N1eqN2$ips_cat<- round(SILIeq_params_df_2spSHORT_N1eqN2$ips)

SILIeq_params_df_2spSHORT_N1eqN2$ips_cat[SILIeq_params_df_2spSHORT_N1eqN2$ips_cat == 7]<- '7 days infectious'
SILIeq_params_df_2spSHORT_N1eqN2$ips_cat[SILIeq_params_df_2spSHORT_N1eqN2$ips_cat == 5]<- '5 days infectious'
SILIeq_params_df_2spSHORT_N1eqN2$ips_cat[SILIeq_params_df_2spSHORT_N1eqN2$ips_cat == 3]<- '3 days infectious'


SILIeq_params_df_2spSHORT_N1eqN2$lat_dur<- 1/SILIeq_params_df_2spSHORT_N1eqN2$epsilon

SILIeq_params_df_2spSHORT_N1eqN2$lat_dur_cat<- round(SILIeq_params_df_2spSHORT_N1eqN2$lat_dur)

SILIeq_params_df_2spSHORT_N1eqN2$lat_dur_cat[SILIeq_params_df_2spSHORT_N1eqN2$lat_dur_cat == 90]<- '90 days latent'
SILIeq_params_df_2spSHORT_N1eqN2$lat_dur_cat[SILIeq_params_df_2spSHORT_N1eqN2$lat_dur_cat == 60]<- '60 days latent'
SILIeq_params_df_2spSHORT_N1eqN2$lat_dur_cat[SILIeq_params_df_2spSHORT_N1eqN2$lat_dur_cat == 30]<- '30 days latent'


SILIeq_params_df_2spSHORT_N1eqN2$psi_perc<- SILIeq_params_df_2spSHORT_N1eqN2$psi*100


SILIeq_params_df_2spSHORT_N1eqN2$beta<- as.character(SILIeq_params_df_2spSHORT_N1eqN2$beta)
SILIeq_params_df_2spSHORT_N1eqN2$beta<- factor(SILIeq_params_df_2spSHORT_N1eqN2$beta, 
                                               levels = c('0.0005', '0.0025', '0.005', '0.0075'))

#Supplemental Figure S10
ggplot(data = SILIeq_params_df_2spSHORT_N1eqN2, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.9) +  
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = expression(SILI~N[A]~'='~N[B])) +
  ylim(-0.01, 0.5)




#for long infectious periods and latency durations
#make empty list to put output dataframes into
SILIeq_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_2sp_long)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp_long[i, 1]
  psi<- exp(-6*(1 - params_combo_2sp_long[i, 2]))
  gamma<- params_combo_2sp_long[i, 3]
  epsilon<- params_combo_2sp_long[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu_a = mu_a, mu_b = mu_b, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", 
                      "Sb", "Ib", "Lb")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ia) / 
    (SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
       SILI_out$Ib + SILI_out$La + SILI_out$Lb)
  
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- SILI_out[SILI_out$time == 7300, ]
  
}


#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)

#Add phylogenetic relatedness parameter values to equilibrium values for a full 
#equilibrium dataset to create plots later
SILIeq_params_df_2spLONG_N1eqN2<- data.frame(SILIeq_df, params_combo_2sp_long)

#For SILI 2sp long ips and lat_dur N1 = N2
SILIeq_params_df_2spLONG_N1eqN2$ips<- 1/SILIeq_params_df_2spLONG_N1eqN2$gamma

SILIeq_params_df_2spLONG_N1eqN2$ips_cat<- round(SILIeq_params_df_2spLONG_N1eqN2$ips)

SILIeq_params_df_2spLONG_N1eqN2$ips_cat[SILIeq_params_df_2spLONG_N1eqN2$ips_cat == 730]<- '730 days infectious'
SILIeq_params_df_2spLONG_N1eqN2$ips_cat[SILIeq_params_df_2spLONG_N1eqN2$ips_cat == 548]<- '548 days infectious'
SILIeq_params_df_2spLONG_N1eqN2$ips_cat[SILIeq_params_df_2spLONG_N1eqN2$ips_cat == 365]<- '365 days infectious'

SILIeq_params_df_2spLONG_N1eqN2$lat_dur<- 1/SILIeq_params_df_2spLONG_N1eqN2$epsilon

SILIeq_params_df_2spLONG_N1eqN2$lat_dur_cat<- round(SILIeq_params_df_2spLONG_N1eqN2$lat_dur)

SILIeq_params_df_2spLONG_N1eqN2$lat_dur_cat[SILIeq_params_df_2spLONG_N1eqN2$lat_dur_cat == 2190]<- '2190 days latent'
SILIeq_params_df_2spLONG_N1eqN2$lat_dur_cat[SILIeq_params_df_2spLONG_N1eqN2$lat_dur_cat == 1643]<- '1643 days latent'
SILIeq_params_df_2spLONG_N1eqN2$lat_dur_cat[SILIeq_params_df_2spLONG_N1eqN2$lat_dur_cat == 1095]<- '1095 days latent'


SILIeq_params_df_2spLONG_N1eqN2$psi_perc<- SILIeq_params_df_2spLONG_N1eqN2$psi*100


SILIeq_params_df_2spLONG_N1eqN2$beta<- as.character(SILIeq_params_df_2spLONG_N1eqN2$beta)
SILIeq_params_df_2spLONG_N1eqN2$beta<- factor(SILIeq_params_df_2spLONG_N1eqN2$beta, 
                                              levels = c('0.0005', '0.0025', '0.005', '0.0075'))

#Supplemental Figure S11
ggplot(data = SILIeq_params_df_2spLONG_N1eqN2, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = expression(SILI~N[A]~'='~N[B])) +
  ylim(-0.01, 0.5)











## starting population values (Na > Nb)
Sa0<- 149 
Ia0<- 1
La0<- 0

Sb0<- 74
Ib0<- 1
Lb0<- 0

Y0<- c(Sa0, Ia0, La0, Sb0, Ib0, Lb0)


#for short infectious periods and latency durations
#make empty list to put output dataframes into
SILIeq_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_2sp_short)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp_short[i, 1]
  psi<- exp(-6*(1 - params_combo_2sp_short[i, 2]))
  gamma<- params_combo_2sp_short[i, 3]
  epsilon<- params_combo_2sp_short[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu_a = mu_a, mu_b = mu_b, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", 
                      "Sb", "Ib", "Lb")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib) / 
    (SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
       SILI_out$Ib + SILI_out$La + SILI_out$Lb)
  
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- SILI_out[SILI_out$time == 7300, ]
  
}

#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)

#Add parameter values to equilibrium values for a full equilibrium dataset to create plots later
SILIeq_params_df_2spSHORT_N1grN2<- data.frame(SILIeq_df, params_combo_2sp_short)

#make ips and imm_durs/latency, and %relatedness, and beta categorical for equilibrium datasets
SILIeq_params_df_2spSHORT_N1grN2$ips<- 1/SILIeq_params_df_2spSHORT_N1grN2$gamma

SILIeq_params_df_2spSHORT_N1grN2$ips_cat<- round(SILIeq_params_df_2spSHORT_N1grN2$ips)

SILIeq_params_df_2spSHORT_N1grN2$ips_cat[SILIeq_params_df_2spSHORT_N1grN2$ips_cat == 7]<- '7 days infectious'
SILIeq_params_df_2spSHORT_N1grN2$ips_cat[SILIeq_params_df_2spSHORT_N1grN2$ips_cat == 5]<- '5 days infectious'
SILIeq_params_df_2spSHORT_N1grN2$ips_cat[SILIeq_params_df_2spSHORT_N1grN2$ips_cat == 3]<- '3 days infectious'


SILIeq_params_df_2spSHORT_N1grN2$lat_dur<- 1/SILIeq_params_df_2spSHORT_N1grN2$epsilon

SILIeq_params_df_2spSHORT_N1grN2$lat_dur_cat<- round(SILIeq_params_df_2spSHORT_N1grN2$lat_dur)

SILIeq_params_df_2spSHORT_N1grN2$lat_dur_cat[SILIeq_params_df_2spSHORT_N1grN2$lat_dur_cat == 90]<- '90 days latent'
SILIeq_params_df_2spSHORT_N1grN2$lat_dur_cat[SILIeq_params_df_2spSHORT_N1grN2$lat_dur_cat == 60]<- '60 days latent'
SILIeq_params_df_2spSHORT_N1grN2$lat_dur_cat[SILIeq_params_df_2spSHORT_N1grN2$lat_dur_cat == 30]<- '30 days latent'


SILIeq_params_df_2spSHORT_N1grN2$psi_perc<- SILIeq_params_df_2spSHORT_N1grN2$psi*100


SILIeq_params_df_2spSHORT_N1grN2$beta<- as.character(SILIeq_params_df_2spSHORT_N1grN2$beta)
SILIeq_params_df_2spSHORT_N1grN2$beta<- factor(SILIeq_params_df_2spSHORT_N1grN2$beta, 
                                               levels = c('0.0005', '0.0025', '0.005', '0.0075'))

#Supplemental Figure S12
ggplot(data = SILIeq_params_df_2spSHORT_N1grN2, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.9) +  
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = expression(SILI~N[A]~'>'~N[B])) +
  ylim(-0.01, 0.5)




#for short infectious periods and latency durations
#make empty list to put output dataframes into
SILIeq_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(params_combo_2sp_long)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- params_combo_2sp_long[i, 1]
  psi<- exp(-6*(1 - params_combo_2sp_long[i, 2]))
  gamma<- params_combo_2sp_long[i, 3]
  epsilon<- params_combo_2sp_long[i, 4]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi = psi, gamma = gamma, epsilon = epsilon, 
           mu_a = mu_a, mu_b = mu_b, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", 
                      "Sb", "Ib", "Lb")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ia) / 
    (SILI_out$Sa + SILI_out$Sb + SILI_out$Ia + 
       SILI_out$Ib + SILI_out$La + SILI_out$Lb)
  
  
  #get all equilibrium only data
  SILIeq_datalist[[i]]<- SILI_out[SILI_out$time == 7300, ]
  
}


#make equilibrium data one df
SILIeq_df<- do.call(rbind, SILIeq_datalist)

#equilibrium dataset to create plots later
SILIeq_params_df_2spLONG_N1grN2<- data.frame(SILIeq_df, params_combo_2sp_long)

#For SILI 2sp long ips and lat_dur N1 = N2
SILIeq_params_df_2spLONG_N1grN2$ips<- 1/SILIeq_params_df_2spLONG_N1grN2$gamma

SILIeq_params_df_2spLONG_N1grN2$ips_cat<- round(SILIeq_params_df_2spLONG_N1grN2$ips)

SILIeq_params_df_2spLONG_N1grN2$ips_cat[SILIeq_params_df_2spLONG_N1grN2$ips_cat == 730]<- '730 days infectious'
SILIeq_params_df_2spLONG_N1grN2$ips_cat[SILIeq_params_df_2spLONG_N1grN2$ips_cat == 548]<- '548 days infectious'
SILIeq_params_df_2spLONG_N1grN2$ips_cat[SILIeq_params_df_2spLONG_N1grN2$ips_cat == 365]<- '365 days infectious'

SILIeq_params_df_2spLONG_N1grN2$lat_dur<- 1/SILIeq_params_df_2spLONG_N1grN2$epsilon

SILIeq_params_df_2spLONG_N1grN2$lat_dur_cat<- round(SILIeq_params_df_2spLONG_N1grN2$lat_dur)

SILIeq_params_df_2spLONG_N1grN2$lat_dur_cat[SILIeq_params_df_2spLONG_N1grN2$lat_dur_cat == 2190]<- '2190 days latent'
SILIeq_params_df_2spLONG_N1grN2$lat_dur_cat[SILIeq_params_df_2spLONG_N1grN2$lat_dur_cat == 1643]<- '1643 days latent'
SILIeq_params_df_2spLONG_N1grN2$lat_dur_cat[SILIeq_params_df_2spLONG_N1grN2$lat_dur_cat == 1095]<- '1095 days latent'


SILIeq_params_df_2spLONG_N1grN2$psi_perc<- SILIeq_params_df_2spLONG_N1grN2$psi*100


SILIeq_params_df_2spLONG_N1grN2$beta<- as.character(SILIeq_params_df_2spLONG_N1grN2$beta)
SILIeq_params_df_2spLONG_N1grN2$beta<- factor(SILIeq_params_df_2spLONG_N1grN2$beta, 
                                              levels = c('0.0005', '0.0025', '0.005', '0.0075'))

#Supplemental Figure S13
ggplot(data = SILIeq_params_df_2spLONG_N1grN2, aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = expression(SILI~N[A]~'>'~N[B])) +
  ylim(-0.01, 0.5)







#Subplots for main text Figure 3
ggplot(data = SIRSeq_params_df_2spSHORT_N1eqN2[SIRSeq_params_df_2spSHORT_N1eqN2$ips_cat == '7 days infectious' & SIRSeq_params_df_2spSHORT_N1grN2$imm_dur_cat == '30 days immune', ], aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = expression(SIRS~N[A]~'='~N[B])) +
  ylim(-0.01, 0.5)


ggplot(data = SIRSeq_params_df_2spSHORT_N1grN2[SIRSeq_params_df_2spSHORT_N1grN2$ips_cat == '7 days infectious' & SIRSeq_params_df_2spSHORT_N1grN2$imm_dur_cat == '30 days immune', ], aes(x = psi_perc, y = prev_roost, color = beta, group = beta)) +
  geom_path(aes(color = beta), linewidth = 0.75) +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.9) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = '% Relatedness', y = 'Roost Prevalence', colour = expression(beta), title = expression(SIRS~N[A]~'>'~N[B])) +
  ylim(-0.01, 0.5)



