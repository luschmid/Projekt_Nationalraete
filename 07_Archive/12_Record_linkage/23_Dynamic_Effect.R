library(scales)
library(haven)
library(tidyverse)
library(stringr)
rm(list=ls())


setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte")

Reshape_TexResults <- function(data_input){
  
  return(data_input %>% mutate(V2= as.numeric((str_replace_all(data_input$V2, "[\\(\\)\\\\]", ""))),
                               varname=rep(c("coef","se"),length(time_window)),
                               time=ifelse(!is.na(V1),V1,lag(V1))) %>%
           rename(value=V2) %>%
           select(varname,time,value) %>%
           pivot_wider(id_cols = time,names_from = varname))
}

#---------------------------------------------
# Figure XY:
#---------------------------------------------


# (i) Extensive margin

file_name <- "dynamic_effect_ext"
time_window <- c(-4:4)

model_results_ext <- data.table::fread(paste0("./05_Texts_and_presns/03_Political_Rents/tables/",file_name,".tex"),
                                       sep="&")


model_results_ext_wide <- Reshape_TexResults(data_input=model_results_ext)


ggplot(model_results_ext_wide,aes(x=time,y=coef)) +
  geom_point()+
  geom_errorbar(aes(ymin=coef-1.96*se, ymax=coef+1.96*se), width=.2,
                position=position_dodge(0.05)) +
  theme_bw(base_size=24) +
  scale_x_continuous(breaks=time_window) +
  xlab("Time relative to election year") + ylab("RD estimate")



# (ii) Intensive margin

file_name <- "dynamic_effect_int"
time_window <- c(-4:4)

model_results_int <- data.table::fread(paste0("./05_Texts_and_presns/03_Political_Rents/tables/",file_name,".tex"),
                                   sep="&")


model_results_int_wide <- Reshape_TexResults(data_input=model_results_int)

ggplot(model_results_int_wide,aes(x=time,y=coef)) +
  geom_point()+
  geom_errorbar(aes(ymin=coef-1.96*se, ymax=coef+1.96*se), width=.2,
                position=position_dodge(0.05)) +
  theme_bw(base_size=24) +
  scale_x_continuous(breaks=time_window) +
  xlab("Time relative to election year") + ylab("RD estimate")
  

                         