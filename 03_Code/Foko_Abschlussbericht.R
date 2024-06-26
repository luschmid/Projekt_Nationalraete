#-------------------------------------------------
# (A) Libraries, working directory and load data
#-------------------------------------------------

rm(list = ls())

# (i) Libraries, working directory, and load functions

library(readxl)
library(tidyverse)
library(data.table)
library(haven)
library(xtable)
library(janitor)
library(vtable)

path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte"
setwd(path)

#-------------------------------------------------
# (A) Libraries, working directory and load data
#-------------------------------------------------


df_directors <- read_dta("02_Processed_data/Politicians_Directorships_1931-2017.dta")

df_directors_small <- df_directors %>% 
  filter(!is.na(elected)) 

df_directors_small %>% 
  count(elected,i_all_s1) %>%   # dplyr function
  adorn_totals()        # janitor function

out_op <- "latex" # latex or return

# (i) All candidates

sumtable(df_directors_small, vars = c('i_all_s1'),
         out=out_op,
         label=c("At least one mandate"),
                 summ=c('notNA(x)',
                        'mean(x)',
                        'sd(x)',
                        'min(x)',
                        'max(x)'))

sumtable(df_directors_small %>% filter(i_all_s1==1), 
         vars = c('n_all_sum_s1'),
         out=out_op,
         label=c("Number of mandates"),
         summ=c('notNA(x)',
                'mean(x)',
                'sd(x)',
                'min(x)',
                'max(x)'))

# (ii) Elected candidates

sumtable(df_directors_small %>% filter(elected==1), vars = c('i_all_s1'),
         out=out_op,
         label=c("At least one mandate"),
         summ=c('notNA(x)',
                'mean(x)',
                'sd(x)',
                'min(x)',
                'max(x)'))

sumtable(df_directors_small %>% filter(i_all_s1==1 & elected==1), 
         vars = c('n_all_sum_s1'),
         out=out_op,
         label=c("Number of mandates"),
         summ=c('notNA(x)',
                'mean(x)',
                'sd(x)',
                'min(x)',
                'max(x)'))

# (iii) Not elected candidates

sumtable(df_directors_small %>% filter(elected==0), vars = c('i_all_s1'),
         out=out_op,
         label=c("At least one mandate"),
         summ=c('notNA(x)',
                'mean(x)',
                'sd(x)',
                'min(x)',
                'max(x)'))

sumtable(df_directors_small %>% filter(i_all_s1==1 & elected==0), 
         vars = c('n_all_sum_s1'),
         out=out_op,
         label=c("Number of mandates"),
         summ=c('notNA(x)',
                'mean(x)',
                'sd(x)',
                'min(x)',
                'max(x)'))




df_directors_small %>% 
  select(ID,i_all_s1_F1, i_all_s1_F2,i_all_s1_F3, i_all_s1_F4)
