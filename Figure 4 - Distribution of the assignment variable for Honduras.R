# Figure 4: Distribution of the assignment variable for Honduras

# 1. Load packages and set directory

library(scales)
library(tidyverse)
library(readstata13)
library(janitor)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))


# 2. Read in data on Honduras 

rm(list = ls(all = TRUE))
data <- read.dta13("./02_Processed_data/15_Elections_Honduras/elections_hn_final.dta")

# 3. Histogram of running variable

ggplot(data, aes(x = votemargin_rel)) +
  geom_histogram(aes(y = ..count.. / sum(..count..)),
                 alpha=0.5,position = "stack",boundary=0,color="black",
                 binwidth=0.002) +
  theme_bw(base_size = 42) +
  scale_x_continuous(limits = c(-0.07, 0.07), 
                     breaks = seq(-0.06, 0.06, 0.02)) +
  geom_vline(xintercept = 0) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("") + xlab("")

ggsave(file = "./05_Texts_and_presns/01_Running_variable/figures/votemargin_rel_bw_hn_0_07.pdf", width = 20, height = 12)
