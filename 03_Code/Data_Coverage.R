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

#----------------------------------
# (B) Generate data coverage graph
#----------------------------------

# (i) Lookup function

Lookup <- function(year_pol,horizon){
# This function returns 1 if a given "horizon" is covered by the Sugarcube
# data for a specific election year ("year_pol"), and 0 otherwise. 
  
year<- c(1934,1943,1960,1962,1963,1964,1965,1966,1969,1972,1975,1979,1980,
              1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,
              1994,1995,1996,1997,1998,1999,2000,2001,2002,2003)
return(ifelse((year_pol+horizon) %in% year,1,0))
}

# (ii) Define all election years 

year_pol_all <- seq(1931,2003,4)

# (iii) Loop through all horizons for vector of all election years

out <- matrix(nrow=length(year_pol_all),ncol=13)

for(i in 1:13){
out[,i] <- unlist(lapply(year_pol_all,Lookup,horizon=(i-5)))
}

colnames(out) <- paste0(c(-4:8))
rownames(out) <- paste0(year_pol_all)

# (iv) Reshape matrix

out_df <-
  out %>%
  as_tibble() 
  
out_df$year_pol <- year_pol_all  

out_df_long <- out_df %>%
  pivot_longer(cols=c(1:13),names_to = "horizon")

# (v) Heat map

ggplot(out_df_long, aes(y=year_pol, x=horizon)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "darkgreen") + 
  theme_minimal() +
  scale_y_continuous(breaks=year_pol_all) +
  ylab("Year") + xlab("Horizon") + 
  theme(legend.position="none")
  
ggsave("./04_Results/04_Political_Rents/Data_Coverage.pdf",width = 8,height = 8)

out_horizon <- out_df_long %>%
  filter(value==1) %>%
  group_by(horizon) %>%
  summarize(year_mean=mean(year_pol)) %>%
  mutate(ordering=ifelse)

ggplot(out_horizon, aes(y=year_mean, x=horizon)) +
  geom_point()  +
  ylab("Year average") + xlab("Horizon") + 
  theme_minimal()

ggsave("./04_Results/04_Political_Rents/Data_Coverage_Mean.pdf",width = 8,height = 8)


#-------------------------------------
# (C) Generate overlap in regressions
#-------------------------------------

data_i_lrg_c1 <- read_dta(paste0(path, "./02_Processed_data/20_Analysis/sample_composition_i_lrg_c1.dta")) %>%
  select(starts_with("s_i"))

df_colnames <- colnames(data_i_lrg_c1)

data_i_lrg_c1$s_i_lrg_c1_L3_ob[data_i_lrg_c1$s_i_lrg_c1_L4_ob==1 & data_i_lrg_c1$s_i_lrg_c1_L3_ob==1 ]
data_i_lrg_c1$s_i_lrg_c1_L3_ob[data_i_lrg_c1$s_i_lrg_c1_L4_ob==1 & data_i_lrg_c1$s_i_lrg_c1_L3_ob==0 ]


