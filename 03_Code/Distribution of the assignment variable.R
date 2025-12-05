
# Figure 2: Distribution of the assignment variable for Switzerland

# 1. Load packages and set directory

library(scales)
library(tidyverse)
library(haven)
library(janitor)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))

# 2. Read in data on Switzerland 

rm(list = ls(all = TRUE))
data <- read_dta("./02_Processed_data/Politicians_Directorships_1931-2017.dta") %>%
  mutate(outcome_TL1=rowMeans(cbind(i_all_c1_L3,i_all_c1_L2, i_all_c1_L1, i_all_c1), na.rm = TRUE)) %>%
  mutate(outcome_T0=rowMeans(cbind(i_all_c1_F1,i_all_c1_F2, i_all_c1_F3, i_all_c1_F4), na.rm = TRUE)) %>%
  mutate(outcome_TF1=rowMeans(cbind(i_all_c1_F5,i_all_c1_F6, i_all_c1_F7, i_all_c1_F8), na.rm = TRUE)) 
data %>% tabyl(canton,year) 

data %>% filter((is.na(outcome_TL1)==F |
                  is.na(outcome_T0)==F | 
                  is.na(outcome_TF1)==F) &
                  is.na(votemargin_rel) ==F  ) 
  
#data$<- data$votemargin/data$eligible_cant

# 3. Histogram function

PlotHistogram <- function(datainput,window,binwidth_input=0.01,runningvariable,label=""){
# This function plots the distribution of the running variable votemargin_rel and saves it into our folder "05_Texts_and_presns/01_Running_variable"
  break_points <- unique(c(-rev(seq(0,window,binwidth_input)),seq(0,window,binwidth_input)))
    ggplot(datainput, aes_string(x = runningvariable)) +
    geom_histogram(aes(y = ..count.. / sum(..count..)), fill = "grey", colour = "black", 
                    breaks = break_points) +
    theme_bw(base_size = 42) +
    scale_x_continuous(limits = c(-window, window)) +
    geom_vline(xintercept = 0) +
    ylab("") + xlab("")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  ggsave(file = paste("./05_Texts_and_presns/03_Political_Rents/figures/votemargin_rel_polrents", 
                      window, "_",label,".pdf", sep = ""), width = 20, height = 12)
}

# 4. Histogram


PlotHistogram(datainput=data %>% filter(is.na(outcome_TL1)==F |
                                          is.na(outcome_T0)==F | 
                                          is.na(outcome_TF1)==F ),
              window=0.11,binwidth_input=0.002,
              runningvariable="votemargin_rel",label="")
