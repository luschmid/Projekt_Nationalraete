
# Figure A.1: Comparing our overall margin with the candidate margin

# 1. Load packages, set directory, and ggplot layers

library(tidyverse)
library(readstata13)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))
source("./03_Code/fine_grid.R") # ggplot layers

# 2. Read in data on Switzerland 

rm(list = ls(all = TRUE))
data <- read.dta13("./02_Processed_data/nationalraete_1931_2015.dta") 

data_elected <- data %>% filter(elected==1) %>% group_by(year, canton,list) %>% summarize(votes_h_min=min(votes))
data_notelected <- data %>% filter(elected==0) %>% group_by(year, canton,list) %>% summarize(votes_h_max=max(votes))

data_all <- data %>% select(year, canton,list,votes,votemargin,elected,eligible_cant) %>% 
                     left_join(data_elected,by=c("year", "canton","list")) %>% 
                     left_join(data_notelected,by=c("year", "canton","list"))%>%
                      mutate(votemargin_cand=ifelse(elected==1,(votes-votes_h_max),(votes-votes_h_min)))

# 3. Quality checks and destats

data_all %>% group_by(elected) %>% summarize(votemargin_cand_min=min(votemargin_cand,na.rm=T),
                                                                    votemargin_cand_max=max(votemargin_cand,na.rm=T))

data_all %>% group_by(elected) %>% summarize(votemargin_min=min(votemargin,na.rm=T),
                                             votemargin_max=max(votemargin,na.rm=T))

data_all <- data_all %>% mutate(votemargin_cand_rel=votemargin_cand/eligible_cant,
                                votemargin_rel=votemargin/eligible_cant,
                                votemargin_binding=ifelse(votemargin_rel==votemargin_cand_rel,1,0),
                                novotemargin=ifelse(is.na(votemargin_cand_rel),1,0))

data_all %>% tabyl(votemargin_binding,novotemargin)
data_all$votemargin_cand_rel[is.na(data_all$votemargin_cand_rel)==T] <- -0.5

# 4. Graph

ggplot(data_all,aes(x=votemargin_rel,y=votemargin_cand_rel))+
  geom_point() +
  theme_bw_finegrid(base_size=24)+
  xlab("") + ylab("") +ggtitle("") +
  theme(plot.caption = element_text(hjust = 0, face= "italic"), #Default is hjust=1
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot", 
        plot.title = element_text(size = 24)) +
  scale_y_continuous(breaks=seq(-0.5,0.5,0.25),labels=c("NA",seq(-0.25,0.5,0.25)))

ggsave(file = "./05_Texts_and_presns/01_Running_variable/figures/votemargin_comparison.png", dpi = 900, width = 20, height = 12)
