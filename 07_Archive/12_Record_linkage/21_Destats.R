library(scales)
library(haven)
library(tidyverse)
rm(list=ls())


setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte")
nr_data <- read_dta("./02_Processed_data/nationalraete_1931_2015.dta")

#---------------------------------------------
# Figure XY:
#---------------------------------------------

# Note: For this, we need the csv file results_wide_XY that includes the results from the 
#       RL algorithm


model_name <-  "Bisnode_only_1994_2017_Gen7_optimal_corrected" 

model_results <- data.table::fread(paste0("./02_Processed_data/12_Record_linkage/03_Data_analysis/results_wide_",model_name,".csv"),
                                   sep=";")

nr_data_all <- nr_data %>%
  left_join(model_results,by=c("year","ID")) %>%
  filter(year>=1994) %>%
  mutate(votemargin_rel_brackets=round(votemargin_rel ,1)) %>%
  rowwise() %>%
  mutate(mandate_lead_mean_1_4 = mean(c(mandate_lead_1,mandate_lead_2,mandate_lead_3, mandate_lead_4)),
         mandate_lag_mean_1_4 = mean(c(mandate_lag_1,mandate_lag_2,mandate_lag_3, mandate_lag_4)),
         mandate_relative=mandate_lead_mean_1_4/mandate_lag_mean_1_4) 
  
mandate_lead_mean_1_4_avg <- nr_data_all %>%
  group_by(votemargin_rel_brackets) %>%
  summarize(mandate_lead_mean_1_4=mean(mandate_lead_mean_1_4,na.rm=T),
            mandate_lag_mean_1_4=mean(mandate_lag_mean_1_4,na.rm=T),
            mandate_relative=mean(mandate_relative,na.rm=T))

ggplot(mandate_lead_mean_1_4_avg,aes(x=votemargin_rel_brackets,y=mandate_lead_mean_1_4)) +
  geom_point() +
  theme_bw()

ggplot(mandate_lead_mean_1_4_avg,aes(x=votemargin_rel_brackets,y=mandate_lag_mean_1_4)) +
  geom_point()+
  theme_bw()

ggplot(mandate_lead_mean_1_4_avg,aes(x=votemargin_rel_brackets,y=mandate_relative)) +
  geom_point()+
  theme_bw()




#---------------------------------------------
# Illustration table for presentation
#---------------------------------------------

View(nr_data %>% arrange(canton,name,ID,year))

nr_data_short <- nr_data %>% 
  arrange(canton,ID,year) %>%
  select(ID,year,name,firstname,listname_bfs,elected, votes,votemargin) %>%
  filter(ID %in% c("ZH-1995-9180","ZH-1979-0051","ZH-1999-0083")) 

print(xtable::xtable(nr_data_short,digits=c(0)), include.rownames=FALSE)
