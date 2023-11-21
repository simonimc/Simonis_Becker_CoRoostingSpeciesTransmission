#---
#title: CSTmisc_plots_20231121
#author: Molly Simonis
#date: 2023-11-21
#---

## load packages
library(ape)
library(caper)
library(ggplot2)
library(reshape2)
library(deSolve)


#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

#disable scientific notation
options(scipen = 999)


## set theme
th=theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        strip.text=element_text(size=10),
        legend.text=element_text(size=10))+
  #theme(legend.position = "top")+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

## load in Upham's phylogeny and taxonomy from 
tree<- read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')
taxa<- read.csv('taxonomy_mamPhy_5911species.csv', header=T)

## clean up labels
tree$tip.label<- sapply(strsplit(tree$tip.label, '_'), function(x) paste(x[1], x[2], sep=' '))
taxa$tiplabel<- sapply(strsplit(taxa$tiplabel, '_'), function(x) paste(x[1], x[2], sep=' '))

## subset to bats
taxa<- taxa[taxa$ord == "CHIROPTERA", ]
tree<- keep.tip(tree, taxa$tiplabel)

## get matrix
cmat<- vcv.phylo(tree, cor = T)

## keep only half and nix the diagonal
cmat[lower.tri(cmat,diag = T)] = NA

## melt y
cset<- melt(cmat)
cset<- cset[!is.na(cset$value), ]

## view Desmodus specifically
dset<- cset[cset$Var1=="Desmodus rotundus", ]


#Main Text Figure 2
#plot distribution of relatedness between Neotropical bat species with lines for specific %s for vampire bat relatedness examples 
ggplot(cset,aes(value * 100))+
  geom_histogram(fill = 'slateblue4', color = 'slateblue4')+
  geom_vline(xintercept = 8, color = 'yellowgreen', lwd = 1)+ #this is Desmodus rotundus with Myotis velifer 
  geom_vline(xintercept = 55, color = 'springgreen4', lwd = 1)+ #this is Desmodus rotundus and Phylostomous discolor
  geom_vline(xintercept = 64, color = 'turquoise4', lwd = 1)+ #this is Desmodus rotundus and Diphylla ecaudata
  geom_vline(xintercept = 79, color = 'deepskyblue3', lwd = 1) + #this is Desmodus rotundus with Diaemus youngi
  th+
  scale_x_continuous(limits = c(-3, 103), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,360000), expand = c(0, 0)) +
  labs(x="Phylogenetic relatedness (%)",
       y="frequency")

#labeled specific species relationships and %s in powerpoint




#for exponential decay curve proof of concept
#make a variable for psi
psi<- rep(c(0.0001, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.9999), 
          times = 10)

#make a variable for lambda
lambda<- c(rep(-1, 21), rep(-2, 21), rep(-3, 21), rep(-4, 21), rep(-5, 21), rep(-6, 21), rep(-7, 21), rep(-8, 21), rep(-9, 21), rep(-10, 21))

#combine them into a dataset
exp_dec_df<- as.data.frame(cbind(psi, lambda))

#create actual decay value
exp_dec_df$y<- exp((1-psi)*lambda)

#make lambda categorical for plotting
exp_dec_df$lambda<- factor(exp_dec_df$lambda, levels = c('-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1'))

#Supplemental Figure S1
ggplot(data = exp_dec_df, aes(x = psi, y = y, color = lambda)) +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.9) +
  th +
  labs(x = 'Correlation Coefficient Between Host Species (psi)', y = 'Fraction of Intraspecific Transmission', 
       color = 'Shape of Exponential Decay (lambda)')



#########################2 species time series######################
#read in a tables of parameter combinations
params_combo_2sp_short<- read.csv('ParamCombos_v17_2sp_short.csv', header = T, sep = ',')
params_combo_2sp_long<- read.csv('ParamCombos_v17_2sp_long.csv', header = T, sep = ',')
#make sure column names are correct3
names(params_combo_2sp_short)<- c('beta', 'psi', 'gamma', 'epsilon')
names(params_combo_2sp_long)<- c('beta', 'psi', 'gamma', 'epsilon')


#pick parameter spaces with 99.99% relatedness across species, the longest infectious and shortest
#immune/latency durations for both short and long parameter spaces
sub_params<- as.data.frame(rbind(params_combo_2sp_short[189, ], params_combo_2sp_long[189, ]))

#add a column for 'short' vs 'long' parameter spaces
sub_params$param_cat<- c('short', 'long')
sub_params$param_cat<- factor(sub_params$param_cat, levels = c('short', 'long'))


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
SIRS_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params[i, 1]
  psi<- exp(-6*(1 - sub_params[i, 2]))
  gamma<- sub_params[i, 3]
  epsilon<- sub_params[i, 4]
  
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
  
  SIRS_out$param_cat<- sub_params[i, 5]
  
  #get all data
  SIRS_datalist[[i]]<- SIRS_out
  
}

SIRS_N1eqN2<- do.call(rbind, SIRS_datalist)


#melt I classes to long format
SIRS_N1eqN2_I<- melt(SIRS_N1eqN2, id.vars = c('time', 'param_cat'), measure.vars = c('Ia', 'Ib'), 
                     variable.name = 'I_spp', value.name = 'I_pop')

#Subplots for S2
SIRS_N1eqN2_Iplot<-ggplot(data = SIRS_N1eqN2_I, aes(x = time, y = I_pop, color = I_spp)) +
  facet_wrap(~param_cat, ncol = 1, nrow = 2, scales = 'free_y') +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.5) +
  th +
  labs(x = 'Time (days)', y = 'Infected Population in Roost', 
       title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'='~N[B]), color = 'Infected Species') +
  xlim(0, 2500)


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
SIRS_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params[i, 1]
  psi<- exp(-6*(1 - sub_params[i, 2]))
  gamma<- sub_params[i, 3]
  epsilon<- sub_params[i, 4]
  
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
  
  SIRS_out$param_cat<- sub_params[i, 5]
  
  #get all data
  SIRS_datalist[[i]]<- SIRS_out
  
}

SIRS_N1grN2<- do.call(rbind, SIRS_datalist)



#melt I classes to long format
SIRS_N1grN2_I<- melt(SIRS_N1grN2, id.vars = c('time', 'param_cat'), measure.vars = c('Ia', 'Ib'), 
                     variable.name = 'I_spp', value.name = 'I_pop')

#Subplot for S2
SIRS_N1grN2_Iplot<-ggplot(data = SIRS_N1grN2_I, aes(x = time, y = I_pop, color = I_spp)) +
  facet_wrap(~param_cat, ncol = 1, nrow = 2, scales = 'free_y') +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.5) +
  th +
  labs(x = 'Time (days)', y = 'Infected Population in Roost', 
       title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]), color = 'Infected Species') +
  xlim(0, 2500)

#Plot S2 
SIRS_N1eqN2_Iplot + SIRS_N1grN2_Iplot + plot_layout(guides = 'collect')




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
SILI_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params[i, 1]
  psi<- exp(-6*(1 - sub_params[i, 2]))
  gamma<- sub_params[i, 3]
  epsilon<- sub_params[i, 4]
  
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
  
  SILI_out$param_cat<- sub_params[i, 5]
  
  #get all data
  SILI_datalist[[i]]<- SILI_out
  
}

SILI_N1eqN2<- do.call(rbind, SILI_datalist)


#melt I classes to long format
SILI_N1eqN2_I<- melt(SILI_N1eqN2, id.vars = c('time', 'param_cat'), measure.vars = c('Ia', 'Ib'), 
                     variable.name = 'I_spp', value.name = 'I_pop')

#Subplot for S3
SILI_N1eqN2_Iplot<-ggplot(data = SILI_N1eqN2_I, aes(x = time, y = I_pop, color = I_spp)) +
  facet_wrap(~param_cat, ncol = 1, nrow = 2, scales = 'free_y') +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.5) +
  th +
  labs(x = 'Time (days)', y = 'Infected Population in Roost', 
       title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'='~N[B]), color = 'Infected Species') +
  xlim(0, 2500)



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
SILI_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params[i, 1]
  psi<- exp(-6*(1 - sub_params[i, 2]))
  gamma<- sub_params[i, 3]
  epsilon<- sub_params[i, 4]
  
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
  
  SILI_out$param_cat<- sub_params[i, 5]
  
  #get all data
  SILI_datalist[[i]]<- SILI_out
  
}

SILI_N1grN2<- do.call(rbind, SILI_datalist)


#melt I classes to long format
SILI_N1grN2_I<- melt(SILI_N1grN2, id.vars = c('time', 'param_cat'), measure.vars = c('Ia', 'Ib'), 
                     variable.name = 'I_spp', value.name = 'I_pop')

#Subplot for S3
SILI_N1grN2_Iplot<-ggplot(data = SILI_N1grN2_I, aes(x = time, y = I_pop, color = I_spp)) +
  facet_wrap(~param_cat, ncol = 1, nrow = 2, scales = 'free_y') +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.5) +
  th +
  labs(x = 'Time (days)', y = 'Infected Population in Roost', 
       title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]), color = 'Infected Species') +
  xlim(0, 2500)


#Plot S3
SILI_N1eqN2_Iplot + SILI_N1grN2_Iplot + plot_layout(guides = 'collect')



#########################3 species time series######################

#read in a tables of parameter combinations
params_combo_3sp_short<- read.csv('ParamCombos_v17_3sp_short.csv', header = T, sep = ',')
params_combo_3sp_long<- read.csv('ParamCombos_v17_3sp_long.csv', header = T, sep = ',')
#make sure column names are correct3
names(params_combo_3sp_short)<- c('psi_ab', 'psi_ac', 'psi_bc', 'beta', 'gamma', 'epsilon')
names(params_combo_3sp_long)<- c('psi_ab', 'psi_ac', 'psi_bc', 'beta', 'gamma', 'epsilon')

#pick parameter spaces with 99.99% relatedness across species, the longest infectious and shortest
#immune/latency durations for both short and long parameter spaces
sub_params<- as.data.frame(rbind(params_combo_3sp_short[27783, ], params_combo_3sp_long[27783, ]))

#add a column for 'short' vs 'long' parameter spaces
sub_params$param_cat<- c('short', 'long')
sub_params$param_cat<- factor(sub_params$param_cat, levels = c('short', 'long'))

#Create function for the SIRS model
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
SIRS_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params[i, 4]
  psi_ab<- exp(-6*(1 - sub_params[i, 1]))
  psi_ac<- exp(-6*(1 - sub_params[i, 2]))
  psi_bc<- exp(-6*(1 - sub_params[i, 3]))
  gamma<- sub_params[i, 5]
  epsilon<- sub_params[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, gamma = gamma, epsilon = epsilon, mu_a = mu_a, mu_b = mu_b, mu_c = mu_c, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", "Sb", "Ib", "Rb", "Sc", "Ic", "Rc")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic) / (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Sc + SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic + SIRS_out$Ra + SIRS_out$Rb + SIRS_out$Rc)
  
  SIRS_out$param_cat<- sub_params[i, 7]
  
  #get all equilibrium only data
  SIRS_datalist[[i]]<- SIRS_out
  
}
SIRS_N1eqN2eqN3<- do.call(rbind, SIRS_datalist)


#melt I classes to long format
SIRS_N1eqN2eqN3_I<- melt(SIRS_N1eqN2eqN3, id.vars = c('time', 'param_cat'), measure.vars = c('Ia', 'Ib', 'Ic'), 
                         variable.name = 'I_spp', value.name = 'I_pop')

#Subplot for S4
SIRS_N1eqN2eqN3_Iplot<-ggplot(data = SIRS_N1eqN2eqN3_I, aes(x = time, y = I_pop, color = I_spp)) +
  facet_wrap(~param_cat, ncol = 1, nrow = 2, scales = 'free_y') +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.8) +
  th +
  labs(x = 'Time (days)', y = 'Infected Population in Roost', 
       title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'='~N[B]~'='~N[C]), color = 'Infected Species') +
  xlim(0, 2500)




#starting population values
#Na > Nb > Nc
Sa0<- 149 
Ia0<- 1
Ra0<- 0

Sb0<- 74
Ib0<- 1
Rb0<- 0

Sc0<- 37
Ic0<- 1
Rc0<- 0
Y0<- c(Sa0, Ia0, Ra0, Sb0, Ib0, Rb0, Sc0, Ic0, Rc0)

#for short infectious periods and waning immunity durations
#make empty list to put output dataframes into
SIRS_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params[i, 4]
  psi_ab<- exp(-6*(1 - sub_params[i, 1]))
  psi_ac<- exp(-6*(1 - sub_params[i, 2]))
  psi_bc<- exp(-6*(1 - sub_params[i, 3]))
  gamma<- sub_params[i, 5]
  epsilon<- sub_params[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, gamma = gamma, epsilon = epsilon, mu_a = mu_a, mu_b = mu_b, mu_c = mu_c, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", "Sb", "Ib", "Rb", "Sc", "Ic", "Rc")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic) / (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Sc + SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic + SIRS_out$Ra + SIRS_out$Rb + SIRS_out$Rc)
  
  SIRS_out$param_cat<- sub_params[i, 7]
  
  #get all equilibrium only data
  SIRS_datalist[[i]]<- SIRS_out
  
}
SIRS_N1grN2grN3<- do.call(rbind, SIRS_datalist)

#melt I classes to long format
SIRS_N1grN2grN3_I<- melt(SIRS_N1grN2grN3, id.vars = c('time', 'param_cat'), measure.vars = c('Ia', 'Ib', 'Ic'), 
                         variable.name = 'I_spp', value.name = 'I_pop')

#Subplot for S4
SIRS_N1grN2grN3_Iplot<-ggplot(data = SIRS_N1grN2grN3_I, aes(x = time, y = I_pop, color = I_spp)) +
  facet_wrap(~param_cat, ncol = 1, nrow = 2, scales = 'free_y') +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.8) +
  th +
  labs(x = 'Time (days)', y = 'Infected Population in Roost', 
       title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'>'~N[C]), color = 'Infected Species') +
  xlim(0, 2500)



#starting population values
#Na > Nb < Nc
Sa0<- 149 
Ia0<- 1
Ra0<- 0

Sb0<- 37
Ib0<- 1
Rb0<- 0

Sc0<- 74
Ic0<- 1
Rc0<- 0
Y0<- c(Sa0, Ia0, Ra0, Sb0, Ib0, Rb0, Sc0, Ic0, Rc0)

#for short infectious periods and waning immunity durations
#make empty list to put output dataframes into
SIRS_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params[i, 4]
  psi_ab<- exp(-6*(1 - sub_params[i, 1]))
  psi_ac<- exp(-6*(1 - sub_params[i, 2]))
  psi_bc<- exp(-6*(1 - sub_params[i, 3]))
  gamma<- sub_params[i, 5]
  epsilon<- sub_params[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, gamma = gamma, epsilon = epsilon, mu_a = mu_a, mu_b = mu_b, mu_c = mu_c, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SIRS_out<- as.data.frame(lsoda(Y0, time, SIRS_model, pars))  
  
  ## name colmns in df
  names(SIRS_out)<- c("time", "Sa", "Ia", "Ra", "Sb", "Ib", "Rb", "Sc", "Ic", "Rc")
  
  #define prevalence in df
  SIRS_out$prev_roost<- (SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic) / (SIRS_out$Sa + SIRS_out$Sb + SIRS_out$Sc + SIRS_out$Ia + SIRS_out$Ib + SIRS_out$Ic + SIRS_out$Ra + SIRS_out$Rb + SIRS_out$Rc)
  
  SIRS_out$param_cat<- sub_params[i, 7]
  
  #get all equilibrium only data
  SIRS_datalist[[i]]<- SIRS_out
  
}
SIRS_N1grN2lsN3<- do.call(rbind, SIRS_datalist)


#melt I classes to long format
SIRS_N1grN2lsN3_I<- melt(SIRS_N1grN2lsN3, id.vars = c('time', 'param_cat'), measure.vars = c('Ia', 'Ib', 'Ic'), 
                         variable.name = 'I_spp', value.name = 'I_pop')

#Subplot for S4
SIRS_N1grN2lsN3_Iplot<-ggplot(data = SIRS_N1grN2lsN3_I, aes(x = time, y = I_pop, color = I_spp)) +
  facet_wrap(~param_cat, ncol = 1, nrow = 2, scales = 'free_y') +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.8) +
  th +
  labs(x = 'Time (days)', y = 'Infected Population in Roost', 
       title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'<'~N[C]), color = 'Infected Species') +
  xlim(0, 2500)


#Plot S4
SIRS_N1eqN2eqN3_Iplot + SIRS_N1grN2grN3_Iplot + SIRS_N1grN2lsN3_Iplot + plot_layout(guides = 'collect')






#Create function for the SILI model
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
SILI_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params[i, 4]
  psi_ab<- exp(-6*(1 - sub_params[i, 1]))
  psi_ac<- exp(-6*(1 - sub_params[i, 2]))
  psi_bc<- exp(-6*(1 - sub_params[i, 3]))
  gamma<- sub_params[i, 5]
  epsilon<- sub_params[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, gamma = gamma, epsilon = epsilon, mu_a = mu_a, mu_b = mu_b, mu_c = mu_c, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", "Sb", "Ib", "Lb", "Sc", "Ic", "Lc")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib + SILI_out$Ic) / (SILI_out$Sa + SILI_out$Sb + SILI_out$Sc + SILI_out$Ia + SILI_out$Ib + SILI_out$Ic + SILI_out$La + SILI_out$Lb + SILI_out$Lc)
  
  SILI_out$param_cat<- sub_params[i, 7]
  
  #get all equilibrium only data
  SILI_datalist[[i]]<- SILI_out
  
}
SILI_N1eqN2eqN3<- do.call(rbind, SILI_datalist)


#melt I classes to long format
SILI_N1eqN2eqN3_I<- melt(SILI_N1eqN2eqN3, id.vars = c('time', 'param_cat'), measure.vars = c('Ia', 'Ib', 'Ic'), 
                         variable.name = 'I_spp', value.name = 'I_pop')

#Subplot for S5 
SILI_N1eqN2eqN3_Iplot<-ggplot(data = SILI_N1eqN2eqN3_I, aes(x = time, y = I_pop, color = I_spp)) +
  facet_wrap(~param_cat, ncol = 1, nrow = 2, scales = 'free_y') +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.8) +
  th +
  labs(x = 'Time (days)', y = 'Infected Population in Roost', 
       title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'='~N[B]~'='~N[C]), color = 'Infected Species') +
  xlim(0, 2500)




#starting population values
#Na > Nb > Nc
Sa0<- 149 
Ia0<- 1
La0<- 0

Sb0<- 74
Ib0<- 1
Lb0<- 0

Sc0<- 37
Ic0<- 1
Lc0<- 0
Y0<- c(Sa0, Ia0, La0, Sb0, Ib0, Lb0, Sc0, Ic0, Lc0)

#for short infectious periods and latency durations
#make empty list to put output dataframes into
SILI_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params[i, 4]
  psi_ab<- exp(-6*(1 - sub_params[i, 1]))
  psi_ac<- exp(-6*(1 - sub_params[i, 2]))
  psi_bc<- exp(-6*(1 - sub_params[i, 3]))
  gamma<- sub_params[i, 5]
  epsilon<- sub_params[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, gamma = gamma, epsilon = epsilon, mu_a = mu_a, mu_b = mu_b, mu_c = mu_c, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", "Sb", "Ib", "Lb", "Sc", "Ic", "Lc")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib + SILI_out$Ic) / (SILI_out$Sa + SILI_out$Sb + SILI_out$Sc + SILI_out$Ia + SILI_out$Ib + SILI_out$Ic + SILI_out$La + SILI_out$Lb + SILI_out$Lc)
  
  SILI_out$param_cat<- sub_params[i, 7]
  
  #get all equilibrium only data
  SILI_datalist[[i]]<- SILI_out
  
}
SILI_N1grN2grN3<- do.call(rbind, SILI_datalist)

#melt I classes to long format
SILI_N1grN2grN3_I<- melt(SILI_N1grN2grN3, id.vars = c('time', 'param_cat'), measure.vars = c('Ia', 'Ib', 'Ic'), 
                         variable.name = 'I_spp', value.name = 'I_pop')

#Subplot for S5 
SILI_N1grN2grN3_Iplot<-ggplot(data = SILI_N1grN2grN3_I, aes(x = time, y = I_pop, color = I_spp)) +
  facet_wrap(~param_cat, ncol = 1, nrow = 2, scales = 'free_y') +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.8) +
  th +
  labs(x = 'Time (days)', y = 'Infected Population in Roost', 
       title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'>'~N[C]), color = 'Infected Species') +
  xlim(0, 2500)



#starting population values
#Na > Nb < Nc
Sa0<- 149 
Ia0<- 1
La0<- 0

Sb0<- 37
Ib0<- 1
Lb0<- 0

Sc0<- 74
Ic0<- 1
Lc0<- 0
Y0<- c(Sa0, Ia0, La0, Sb0, Ib0, Lb0, Sc0, Ic0, Lc0)

#for short infectious periods and latency durations
#make empty list to put output dataframes into
SILI_datalist<- list()


#Simulate across parameter space. The loop will take a few minutes to run
for (i in 1:nrow(sub_params)) {
  
  ## parameters identified by the row number in params_combo df
  beta<- sub_params[i, 4]
  psi_ab<- exp(-6*(1 - sub_params[i, 1]))
  psi_ac<- exp(-6*(1 - sub_params[i, 2]))
  psi_bc<- exp(-6*(1 - sub_params[i, 3]))
  gamma<- sub_params[i, 5]
  epsilon<- sub_params[i, 6]
  
  ## parameter vector identified by the row number in params_combo df
  pars<- c(beta = beta, psi_ab = psi_ab, psi_ac = psi_ac, psi_bc = psi_bc, gamma = gamma, epsilon = epsilon, mu_a = mu_a, mu_b = mu_b, mu_c = mu_c, k = k, b0 = b0, b1 = b1)
  
  ## simultion df identified by the row number in params_combo df
  SILI_out<- as.data.frame(lsoda(Y0, time, SILI_model, pars))  
  
  ## name colmns in df
  names(SILI_out)<- c("time", "Sa", "Ia", "La", "Sb", "Ib", "Lb", "Sc", "Ic", "Lc")
  
  #define prevalence in df
  SILI_out$prev_roost<- (SILI_out$Ia + SILI_out$Ib + SILI_out$Ic) / (SILI_out$Sa + SILI_out$Sb + SILI_out$Sc + SILI_out$Ia + SILI_out$Ib + SILI_out$Ic + SILI_out$La + SILI_out$Lb + SILI_out$Lc)
  
  SILI_out$param_cat<- sub_params[i, 7]
  
  #get all equilibrium only data
  SILI_datalist[[i]]<- SILI_out
  
}
SILI_N1grN2lsN3<- do.call(rbind, SILI_datalist)

#melt I classes to long format
SILI_N1grN2lsN3_I<- melt(SILI_N1grN2lsN3, id.vars = c('time', 'param_cat'), measure.vars = c('Ia', 'Ib', 'Ic'), 
                         variable.name = 'I_spp', value.name = 'I_pop')

#Subplot for S5 
SILI_N1grN2lsN3_Iplot<-ggplot(data = SILI_N1grN2lsN3_I, aes(x = time, y = I_pop, color = I_spp)) +
  facet_wrap(~param_cat, ncol = 1, nrow = 2, scales = 'free_y') +
  geom_line() +
  scale_color_viridis_d(begin = 0, end =  0.8) +
  th +
  labs(x = 'Time (days)', y = 'Infected Population in Roost', 
       title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'<'~N[C]), color = 'Infected Species') +
  xlim(0, 2500)



#Plot S5
SILI_N1eqN2eqN3_Iplot + SILI_N1grN2grN3_Iplot + SILI_N1grN2lsN3_Iplot + plot_layout(guides = 'collect')
