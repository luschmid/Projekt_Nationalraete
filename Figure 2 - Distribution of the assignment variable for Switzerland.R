
# Figure 2: Distribution of the assignment variable for Switzerland

# 1. Load packages and set directory

library(scales)
library(tidyverse)
library(readstata13)
library(janitor)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))

# 2. Read in data on Switzerland 

rm(list = ls(all = TRUE))
data <- read.dta13("./02_Processed_data/nationalraete_1931_2015.dta") 
data %>% tabyl(canton,year) 
  
data$votemargin_rel<- data$votemargin/data$eligible_cant

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

  ggsave(file = paste("./05_Texts_and_presns/01_Running_variable/figures/votemargin_rel_bw_", window, "_",label,".pdf", sep = ""), width = 20, height = 12)
}

# 4. Histogram

PlotHistogram(datainput=data %>% filter(is.na(elected_F1)==F),window=0.11,binwidth_input=0.002,
              runningvariable="votemargin_rel",label="ch")
