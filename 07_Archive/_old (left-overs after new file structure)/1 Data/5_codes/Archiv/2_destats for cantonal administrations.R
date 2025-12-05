setwd("C:/SchmidLu/Dropbox/Projekt Nationalräte/1 Data")

library(foreign)
library(ggplot2)
library(grid)
library(scales)


rm(list=ls())

source("./5_codes/10_fine_grid.R") # ggplot layers
source("./5_codes/10_summarySE.R") # summarize function


data<-read.dta("./1_data/NRW_KANDIDATEN.dta")
#data<-read.dta("./1_data/data_federal_temp.dta")




ggplot(data=data,aes(y=tau1,x=windower)) +  
  geom_point(position=position_dodge(width=0.3),width=.1,height = 0.6,size=2.6,colour = "black",fill = "grey" )  +
  geom_errorbar(aes(ymin=tau_low, ymax=tau_up),width=.5,height = 0.6) +
  theme_bw_finegrid(base_size = 28) +                             
  xlab(" \n Time Difference to Leader of First Run  \n (in Hundredths of a Second)") + ylab("Effect of Low Left Digit\n") +
  scale_x_continuous(breaks=c(label1),limit=c(0,26.9),labels=label) +
  scale_y_continuous(limit=c(-0.04,0.02))



data<-read.dta("./1_data/nr_analysis.dta")



data$candidate_diff_marginal_mod <- data$candidate_diff_marginal+runif(length(data$candidate_diff_marginal),min = 0, max = 10)

data$elected_mod <- matrix(NA,length(data$candidate_diff_marginal),1)

data$elected_mod[data$candidate_diff_marginal_mod>=0] <- 1
data$elected_mod[data$candidate_diff_marginal_mod<0] <- 0

ggplot(data=data,aes(y=elected_mod,x=candidate_diff_marginal_mod)) +  
  geom_point( )  +

    theme_bw_finegrid(base_size = 28) +                             
  xlab(" \n Time Difference to Leader of First Run  \n (in Hundredths of a Second)") + ylab("Effect of Low Left Digit\n") +
  scale_x_continuous(breaks=c(label1),limit=c(0,26.9),labels=label) +
  scale_y_continuous(limit=c(-0.04,0.02))


