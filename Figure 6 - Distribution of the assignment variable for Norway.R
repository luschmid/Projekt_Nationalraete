# Figure 4: Distribution of the assignment variable for Honduras

# 1. Load packages and set directory

library(tidyverse)
library(readstata13)
library(janitor)
library(scales)
library(reprex)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalr?te"))

# 2. Read in data on Norway 

rm(list = ls(all = TRUE))

data <- read.dta13("./02_Processed_data/13_Running_variable/FivaSmithJune2019_processed.dta") %>% 
  select(year,districtid, district, party, partyname, pid, votes,voteshare,
         elected_F1,rank,votemargin_vote,votemargin_fiva) %>% 
  filter(year>=1953 & year<=1981) %>% 
  arrange(districtid,year, partyname,pid) %>%
  mutate(fiva_smith_sample=ifelse(is.na(votemargin_fiva)==F,"Marginal candidate",
                                  "Non-marginal candidate"))
data %>% tabyl(fiva_smith_sample,year)

# 3. Histogram options

theme_set(theme_bw(base_size = 42))
theme_update(legend.key = element_rect(size = 6, fill = "white", colour = "white"), 
             legend.key.size = unit(1.5, "cm")) 

# 4. Histogram

ggplot(data, aes( votemargin_vote,fill=factor(fiva_smith_sample))) +
  geom_histogram(aes(y = ..count.. / sum(..count..)),
                 alpha=0.5,position = "stack",boundary=0,color="black",
                 binwidth=0.005) +
  scale_fill_manual(values=c("black","grey")) +
  scale_x_continuous(limits=c(-0.36,0.36),
                     breaks=seq(-0.36,36,0.12)) +
  geom_vline(xintercept = 0) +
  ylab("") + xlab("")+
  theme(axis.text.y = element_text(angle = 0)) +
  theme(legend.position=c(0.78,0.87),
        legend.key.size = unit(3, 'lines')) +
  theme(legend.title = element_blank(),
        legend.box.just = "right",
        legend.key = element_rect(size = 6),
        legend.key.size = unit(2, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(file = "./05_Texts_and_presns/01_Running_variable/figures/hist_norway_36.pdf", width = 20, height = 12)

