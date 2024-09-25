#---
#title: CSTmisc_plots_20240924
#author: Molly Simonis
#date: 2024-09-24
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


