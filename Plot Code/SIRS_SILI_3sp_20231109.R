#---
#title: SIRS_SILI_3sp_20231109
#author: Molly Simonis
#date: 2023-11-09
#---

#turn on packages
library(ggplot2)
library(ggtern)
library(patchwork)


#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

#disable scientific notation
options(scipen = 999)

#make a theme for plotting
th<- theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(plot.title = element_text(hjust = 0.5))

#load in environments one at a time (to keep computer from crashing-- large data)
load('SIRS3spN1eqN2eqN3_20231109.Rdata')

#fix very small prevalence values for plotting
SIRSeq_params_df_3spSHORT_N1eqN2eqN3$prev_roost[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$prev_roost < 0]<- 0
SIRSeq_params_df_3spLONG_N1eqN2eqN3$prev_roost[SIRSeq_params_df_3spLONG_N1eqN2eqN3$prev_roost < 0]<- 0

####FOR ALL 3sp PLOTS: x vertex = species A, y vertex = species C, z vertex = species B
#beta = 0.0005
#Figure S14
ggtern(data = SIRSeq_params_df_3spSHORT_N1eqN2eqN3[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0005,], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S15
ggtern(data = SIRSeq_params_df_3spLONG_N1eqN2eqN3[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#beta = 0.0025
#Figure S16
ggtern(data = SIRSeq_params_df_3spSHORT_N1eqN2eqN3[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0025*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S17
ggtern(data = SIRSeq_params_df_3spLONG_N1eqN2eqN3[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0025*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#beta = 0.005
#Figure S18
ggtern(data = SIRSeq_params_df_3spSHORT_N1eqN2eqN3[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.005*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S19
ggtern(data = SIRSeq_params_df_3spLONG_N1eqN2eqN3[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.005*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#beta = 0.0075
#Figure S20
ggtern(data = SIRSeq_params_df_3spSHORT_N1eqN2eqN3[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0075*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S21
ggtern(data = SIRSeq_params_df_3spLONG_N1eqN2eqN3[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0075*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()



#Subset plot for Main Text Figure 3
ggtern(data = SIRSeq_params_df_3spSHORT_N1eqN2eqN3[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0005 & SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat == '7 days infectious' & SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur_cat == '30 days immune',], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Subset plot for Main Text Figure 3
ggtern(data = SIRSeq_params_df_3spSHORT_N1eqN2eqN3[SIRSeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0025 & SIRSeq_params_df_3spSHORT_N1eqN2eqN3$ips_cat == '7 days infectious' & SIRSeq_params_df_3spSHORT_N1eqN2eqN3$imm_dur_cat == '30 days immune',], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0025*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()




#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

#disable scientific notation
options(scipen = 999)

#make a theme for plotting
th<- theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(plot.title = element_text(hjust = 0.5))

#load in environments one at a time (to keep computer from crashing-- large data)
load('SIRS3spN1grN2grN3_20231109.Rdata')

#fix very small prevalence values for plotting
SIRSeq_params_df_3spSHORT_N1grN2grN3$prev_roost[SIRSeq_params_df_3spSHORT_N1grN2grN3$prev_roost < 0]<- 0
SIRSeq_params_df_3spLONG_N1grN2grN3$prev_roost[SIRSeq_params_df_3spLONG_N1grN2grN3$prev_roost < 0]<- 0

#beta = 0.0005
#Figure S22
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2grN3[SIRSeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0005,], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#Figure S23
ggtern(data = SIRSeq_params_df_3spLONG_N1grN2grN3[SIRSeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#beta = 0.0025
#Figure S24
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2grN3[SIRSeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0025*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S25
ggtern(data = SIRSeq_params_df_3spLONG_N1grN2grN3[SIRSeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0025*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#beta = 0.005
#Figure S26
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2grN3[SIRSeq_params_df_3spSHORT_N1grN2grN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.005*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S27
ggtern(data = SIRSeq_params_df_3spLONG_N1grN2grN3[SIRSeq_params_df_3spSHORT_N1grN2grN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.005*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#beta = 0.0075
#Figure S28
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2grN3[SIRSeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0075*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S29
ggtern(data = SIRSeq_params_df_3spLONG_N1grN2grN3[SIRSeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0075*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()



#Subset plot for Main Text Figure 3
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2grN3[SIRSeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0005 & SIRSeq_params_df_3spSHORT_N1grN2grN3$ips_cat == '7 days infectious' & SIRSeq_params_df_3spSHORT_N1grN2grN3$imm_dur_cat == '30 days immune',], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Subset plot for Main Text Figure 3
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2grN3[SIRSeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0025 & SIRSeq_params_df_3spSHORT_N1grN2grN3$ips_cat == '7 days infectious' & SIRSeq_params_df_3spSHORT_N1grN2grN3$imm_dur_cat == '30 days immune',], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0025*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()











#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

#disable scientific notation
options(scipen = 999)

#make a theme for plotting
th<- theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(plot.title = element_text(hjust = 0.5))

#load in environments one at a time (to keep computer from crashing-- large data)
load('SIRS3spN1grN2lsN3_20231109.Rdata')

#fix very small prevalence values for plotting
SIRSeq_params_df_3spSHORT_N1grN2lsN3$prev_roost[SIRSeq_params_df_3spSHORT_N1grN2lsN3$prev_roost < 0]<- 0
SIRSeq_params_df_3spLONG_N1grN2lsN3$prev_roost[SIRSeq_params_df_3spLONG_N1grN2lsN3$prev_roost < 0]<- 0

#beta = 0.0005
#Figure S230
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2lsN3[SIRSeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0005,], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#Figure S31
ggtern(data = SIRSeq_params_df_3spLONG_N1grN2lsN3[SIRSeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#beta = 0.0025
#Figure S32
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2lsN3[SIRSeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0025*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S33
ggtern(data = SIRSeq_params_df_3spLONG_N1grN2lsN3[SIRSeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0025*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#beta = 0.005
#Figure S34
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2lsN3[SIRSeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.005*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S35
ggtern(data = SIRSeq_params_df_3spLONG_N1grN2lsN3[SIRSeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.005*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#beta = 0.0075
#Figure S36
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2lsN3[SIRSeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0075*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S37
ggtern(data = SIRSeq_params_df_3spLONG_N1grN2lsN3[SIRSeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(imm_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0075*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()




#Subset plot for Main Text Figure 3
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2lsN3[SIRSeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0005 & SIRSeq_params_df_3spSHORT_N1grN2lsN3$ips_cat == '7 days infectious' & SIRSeq_params_df_3spSHORT_N1grN2lsN3$imm_dur_cat == '30 days immune',], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Subset plot for Main Text Figure 3
ggtern(data = SIRSeq_params_df_3spSHORT_N1grN2lsN3[SIRSeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0025 & SIRSeq_params_df_3spSHORT_N1grN2lsN3$ips_cat == '7 days infectious' & SIRSeq_params_df_3spSHORT_N1grN2lsN3$imm_dur_cat == '30 days immune',], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SIRS*';'~beta~'='~0.0025*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()














#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

#disable scientific notation
options(scipen = 999)

#make a theme for plotting
th<- theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(plot.title = element_text(hjust = 0.5))

#load in environments one at a time (to keep computer from crashing-- large data)
load('SILI3spN1eqN2eqN3_20231109.Rdata')


#beta = 0.0005
#Figure S38
ggtern(data = SILIeq_params_df_3spSHORT_N1eqN2eqN3[SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0005,], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#Figure S39
ggtern(data = SILIeq_params_df_3spLONG_N1eqN2eqN3[SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#beta = 0.0025
#Figure S40
ggtern(data = SILIeq_params_df_3spSHORT_N1eqN2eqN3[SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0025*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S41
ggtern(data = SILIeq_params_df_3spLONG_N1eqN2eqN3[SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0025*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#beta = 0.005
#Figure S42
ggtern(data = SILIeq_params_df_3spSHORT_N1eqN2eqN3[SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.005*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S43
ggtern(data = SILIeq_params_df_3spLONG_N1eqN2eqN3[SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.005*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#beta = 0.0075
#Figure S44
ggtern(data = SILIeq_params_df_3spSHORT_N1eqN2eqN3[SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0075*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S45
ggtern(data = SILIeq_params_df_3spLONG_N1eqN2eqN3[SILIeq_params_df_3spSHORT_N1eqN2eqN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0075*';'~N[A]~'='~N[B]~'='~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()













#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

#disable scientific notation
options(scipen = 999)

#make a theme for plotting
th<- theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(plot.title = element_text(hjust = 0.5))

#load in environments one at a time (to keep computer from crashing-- large data)
load('SILI3spN1grN2grN3_20231109.Rdata')


#beta = 0.0005
#Figure S46
ggtern(data = SILIeq_params_df_3spSHORT_N1grN2grN3[SILIeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0005,], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#Figure S47
ggtern(data = SILIeq_params_df_3spLONG_N1grN2grN3[SILIeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#beta = 0.0025
#Figure S48
ggtern(data = SILIeq_params_df_3spSHORT_N1grN2grN3[SILIeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0025*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S49
ggtern(data = SILIeq_params_df_3spLONG_N1grN2grN3[SILIeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0025*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#beta = 0.005
#Figure S50
ggtern(data = SILIeq_params_df_3spSHORT_N1grN2grN3[SILIeq_params_df_3spSHORT_N1grN2grN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.005*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S51
ggtern(data = SILIeq_params_df_3spLONG_N1grN2grN3[SILIeq_params_df_3spSHORT_N1grN2grN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.005*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#beta = 0.0075
#Figure S52
ggtern(data = SILIeq_params_df_3spSHORT_N1grN2grN3[SILIeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0075*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S53
ggtern(data = SILIeq_params_df_3spLONG_N1grN2grN3[SILIeq_params_df_3spSHORT_N1grN2grN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0075*';'~N[A]~'>'~N[B]~'>'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()











#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

#disable scientific notation
options(scipen = 999)

#make a theme for plotting
th<- theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(plot.title = element_text(hjust = 0.5))

#load in environments one at a time (to keep computer from crashing-- large data)
load('SILI3spN1grN2lsN3_20231109.Rdata')

#fix beta values
SILIeq_params_df_3spSHORT_N1grN2lsN3$beta<- params_combo_3sp_short$beta
SILIeq_params_df_3spLONG_N1grN2lsN3$beta<- params_combo_3sp_short$beta

#beta = 0.0005
#Figure S54
ggtern(data = SILIeq_params_df_3spSHORT_N1grN2lsN3[SILIeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0005,], 
       aes(x = psi_ab, y = psi_ac, z = psi_bc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#Figure S55
ggtern(data = SILIeq_params_df_3spLONG_N1grN2lsN3[SILIeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0005*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#beta = 0.0025
#Figure S56
ggtern(data = SILIeq_params_df_3spSHORT_N1grN2lsN3[SILIeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0025*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S57
ggtern(data = SILIeq_params_df_3spLONG_N1grN2lsN3[SILIeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0025,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0025*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#beta = 0.005
#Figure S58
ggtern(data = SILIeq_params_df_3spSHORT_N1grN2lsN3[SILIeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.005*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S59
ggtern(data = SILIeq_params_df_3spLONG_N1grN2lsN3[SILIeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.005,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.005*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()


#beta = 0.0075
#Figure S60
ggtern(data = SILIeq_params_df_3spSHORT_N1grN2lsN3[SILIeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0075*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()

#Figure S61
ggtern(data = SILIeq_params_df_3spLONG_N1grN2lsN3[SILIeq_params_df_3spSHORT_N1grN2lsN3$beta == 0.0075,], 
       aes(x = psi_ab_perc, y = psi_ac_perc, z = psi_bc_perc, color = prev_roost)) +  
  facet_grid(lat_dur_cat ~ ips_cat) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1, begin = 0, end =  0.9, limits = c(0, 0.5)) +
  th +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = 'A', y = 'C', z = 'B', colour = 'Roost Prevalance', title = expression(SILI*';'~beta~'='~0.0075*';'~N[A]~'>'~N[B]~'<'~N[C])) +
  Tarrowlab('% relatedness') + Larrowlab('% relatedness') + Rarrowlab('% relatedness') +
  theme_arrowlong()












