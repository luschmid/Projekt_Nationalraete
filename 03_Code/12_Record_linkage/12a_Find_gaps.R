library(tidyverse)
library(haven)
library(stringr)
library(readxl)

rm(list = ls())


path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte"
setwd(path)

# (A) Read in Sugarcube, rl data, and gaps ----

sugarcube_data <- data.table::fread("./02_Processed_data/10_Directors_1934_2003/Sugarcube_RLPostProcessing_Persons-Firmnames.csv",
                                    sep = ",", encoding = "UTF-8") %>%
  mutate(time_period = case_when(
    year_sug == 1933 ~ 1 ,
    year_sug == 1942 ~ 2 ,
    year_sug == 1959 ~ 3 ,
    year_sug == 1961 ~ 4 ,
    year_sug == 1962 ~ 5 ,
    year_sug == 1963 ~ 6 ,
    year_sug == 1965 ~ 7 ,
    year_sug == 1966 ~ 8 ,
    year_sug == 1969 ~ 9 ,
    year_sug == 1972 ~ 10,
    year_sug == 1975 ~ 11,
    year_sug == 1979 ~ 12,
    year_sug == 1980 ~ 13,
    year_sug == 1981 ~ 14,
    year_sug == 1982 ~ 15,
    year_sug == 1983 ~ 16,
    year_sug == 1984 ~ 17,
    year_sug == 1985 ~ 18,
    year_sug == 1986 ~ 19,
    year_sug == 1987 ~ 20,
    year_sug == 1988 ~ 21,
    year_sug == 1989 ~ 22,
    year_sug == 1990 ~ 23,
    year_sug == 1991 ~ 24,
    year_sug == 1992 ~ 25,
    year_sug == 1993 ~ 26,
    year_sug == 1994 ~ 27,
    year_sug == 1995 ~ 28,
    year_sug == 1996 ~ 29,
    year_sug == 1997 ~ 30,
    year_sug == 1998 ~ 31,
    year_sug == 1999 ~ 32,
    year_sug == 2000 ~ 33,
    year_sug == 2001 ~ 34,
    year_sug == 2002 ~ 35,
    year_sug == 2003 ~ 36),
    id_1_all=paste0(year_sug,"-",id_sug)
  )

rl_results <- data.table::fread("./02_Processed_data/12_Record_linkage/02_Sugarcube/RL-Results-15122022-100518bc-12f0-4d5d-8b42-1f4484952a2a.csv",
                                    sep = ",", encoding = "UTF-8")%>%
  mutate(time_period = case_when(
    year == 1933 ~ 1 ,
    year == 1942 ~ 2 ,
    year == 1959 ~ 3 ,
    year == 1961 ~ 4 ,
    year == 1962 ~ 5 ,
    year == 1963 ~ 6 ,
    year == 1965 ~ 7 ,
    year == 1966 ~ 8 ,
    year == 1969 ~ 9 ,
    year == 1972 ~ 10,
    year == 1975 ~ 11,
    year == 1979 ~ 12,
    year == 1980 ~ 13,
    year == 1981 ~ 14,
    year == 1982 ~ 15,
    year == 1983 ~ 16,
    year == 1984 ~ 17,
    year == 1985 ~ 18,
    year == 1986 ~ 19,
    year == 1987 ~ 20,
    year == 1988 ~ 21,
    year == 1989 ~ 22,
    year == 1990 ~ 23,
    year == 1991 ~ 24,
    year == 1992 ~ 25,
    year == 1993 ~ 26,
    year == 1994 ~ 27,
    year == 1995 ~ 28,
    year == 1996 ~ 29,
    year == 1997 ~ 30,
    year == 1998 ~ 31,
    year == 1999 ~ 32,
    year == 2000 ~ 33,
    year == 2001 ~ 34,
    year == 2002 ~ 35,
    year == 2003 ~ 36),
    id_1_all=paste0(year,"-",id)
  ) 


df_gaps <- read_excel("./02_Processed_data/12_Record_linkage/02_Sugarcube/13_RL_Output_Round_2_In/Erfassung_4.xlsx",skip=1)


df_after_post_processing <- read_dta("./02_Processed_data/12_Record_linkage/02_Sugarcube/11_RL_Output_Round1_In/RL_Output_Round1_fp_maj_IDs.dta") %>%
  mutate(time_period = case_when(
    year_1 == 1933 ~ 1 ,
    year_1 == 1942 ~ 2 ,
    year_1 == 1959 ~ 3 ,
    year_1 == 1961 ~ 4 ,
    year_1 == 1962 ~ 5 ,
    year_1 == 1963 ~ 6 ,
    year_1 == 1965 ~ 7 ,
    year_1 == 1966 ~ 8 ,
    year_1 == 1969 ~ 9 ,
    year_1 == 1972 ~ 10,
    year_1 == 1975 ~ 11,
    year_1 == 1979 ~ 12,
    year_1 == 1980 ~ 13,
    year_1 == 1981 ~ 14,
    year_1 == 1982 ~ 15,
    year_1 == 1983 ~ 16,
    year_1 == 1984 ~ 17,
    year_1 == 1985 ~ 18,
    year_1 == 1986 ~ 19,
    year_1 == 1987 ~ 20,
    year_1 == 1988 ~ 21,
    year_1 == 1989 ~ 22,
    year_1 == 1990 ~ 23,
    year_1 == 1991 ~ 24,
    year_1 == 1992 ~ 25,
    year_1 == 1993 ~ 26,
    year_1 == 1994 ~ 27,
    year_1 == 1995 ~ 28,
    year_1 == 1996 ~ 29,
    year_1 == 1997 ~ 30,
    year_1 == 1998 ~ 31,
    year_1 == 1999 ~ 32,
    year_1 == 2000 ~ 33,
    year_1 == 2001 ~ 34,
    year_1 == 2002 ~ 35,
    year_1 == 2003 ~ 36),
    id_1_all=paste0(year_1,"-",id_1)
  ) %>%
  mutate(year_pic = case_when(
    year_1 == 1933 ~ 1 ,
    year_1 == 1942 ~ 2 ,
    year_1 == 1959 ~ 3 ,
    year_1 == 1961 ~ 4 ,
    year_1 == 1962 ~ 5 ,
    year_1 == 1963 ~ 6 ,
    year_1 == 1965 ~ 7 ,
    year_1 == 1966 ~ 8 ,
    year_1 == 1969 ~ 9 ,
    year_1 == 1972 ~ 10,
    year_1 == 1975 ~ 11,
    year_1 == 1979 ~ 12,
    year_1 == 1980 ~ 13,
    year_1 == 1981 ~ 14,
    year_1 == 1982 ~ 15,
    year_1 == 1983 ~ 16,
    year_1 == 1984 ~ 17,
    year_1 == 1985 ~ 18,
    year_1 == 1986 ~ 19,
    year_1 == 1987 ~ 20,
    year_1 == 1988 ~ 21,
    year_1 == 1989 ~ 22,
    year_1 == 1990 ~ 23,
    year_1 == 1991 ~ 24,
    year_1 == 1992 ~ 25,
    year_1 == 1993 ~ 26,
    year_1 == 1994 ~ 27,
    year_1 == 1995 ~ 28,
    year_1 == 1996 ~ 29,
    year_1 == 1997 ~ 30,
    year_1 == 1998 ~ 31,
    year_1 == 1999 ~ 32,
    year_1 == 2000 ~ 33,
    year_1 == 2001 ~ 34,
    year_1 == 2002 ~ 35,
    year_1 == 2003 ~ 36),
    id_1_all=paste0(year_1,"-",id_1)
  )

# (B) Search for gaps ----


FindGaps <- function(rownumber,df_gaps,df_sugarcube,df_after_post_processing,rl_results){
# This function tries to find the gaps defined in the specific rownumber of
# the dataframe df_gaps. The gaps are search in df_sugarcube using exact name
# matches during the gap period and fuzzy firm mathches during the gap period. 


i= rownumber
gaps_time <- (df_gaps$time_period1[i]+1):(df_gaps$time_period2[i]-1)
  
out1 <- sugarcube_data %>% filter(firstname==df_gaps$firstname_sug1[i] &  
                            lastname==df_gaps$lastname_sug1[i] & 
                            time_period %in% gaps_time)

out2 <- sugarcube_data %>% filter(firstname==df_gaps$firstname_sug2[i] &  
                                    lastname==df_gaps$lastname_sug2[i] & 
                                    time_period %in% gaps_time)


out3 <-  sugarcube_data %>%
  filter(agrepl(df_gaps$firmnames1[i], firmnames,
                max.distance=0.2)  &
           time_period %in% gaps_time)

out4 <-  sugarcube_data %>%
  filter(agrepl(df_gaps$firmnames2[i], firmnames,
                max.distance=0.2)  &
           time_period %in% gaps_time)


out5 <- sugarcube_data %>% filter(firstname==df_gaps$firstname_sug1[i] &  
                                    lastname==df_gaps$lastname_sug1[i] & 
                                    time_period==df_gaps$time_period1[i] & 
                                    firmnames==df_gaps$firmnames1[i]) %>%
  mutate(gap="start")

out6 <- sugarcube_data %>% filter(firstname==df_gaps$firstname_sug2[i] &  
                                    lastname==df_gaps$lastname_sug2[i] & 
                                    time_period==df_gaps$time_period2[i] & 
                                    firmnames==df_gaps$firmnames2[i]) %>%
  mutate(gap="end")


id_1_postproc <- df_after_post_processing$id_1_all
id_1_rlresults <- rl_results$id_1_all

out_gap <- bind_rows(out1,out2,out3,out4) %>% 
  distinct(id_sug,year_sug,.keep_all = T) 

return(out5%>% 
  bind_rows(out_gap,out6) %>%
    mutate(odd=i%%2,
           id_in_rlresults=ifelse(id_1_all %in% id_1_rlresults,1,0),
           id_in_postproc=ifelse(id_1_all %in% id_1_postproc,1,0))
)

}


start_time <- Sys.time()
out_data <- FindGaps(rownumber=1,df_gaps=df_gaps,df_sugarcube=df_sugarcube)
end_time <- Sys.time()
end_time - start_time

# (C) Create outfile for analysis by Simon and Lukas regarind the source of the gaps

df_gaps %>% filter(is.na(time_period2))

out <- tibble()
# for (i in 1:length(df_gaps$id_0)){
for (i in 1:143){
  out <- bind_rows(out,df_gaps[i,] %>% select(id_0,ends_with("1")) %>% rename_with(~str_remove(., '1')))
  for (j in ((df_gaps[i,]$time_period1+1):(df_gaps[i,]$time_period2-1))){
    out <- bind_rows(out,tibble(time_period=j))  
  }
  out <- bind_rows(out,df_gaps[i,] %>% select(id_0,ends_with("2")) %>% rename_with(~str_remove(., '2')))
  
}


# 
# rl_results_cut <- rl_results %>% filter(firstname==df_gaps$firstname_sug2[i] &  
#                         name==df_gaps$lastname_sug2[i] & 
#                         time_period %in% gaps_time & 
#                         `source file`==1)
# print(rl_results_cut)
#print(rl_results_cut$year)
#print(out$year)
# rl_results_cut_notfound <- rl_results_cut %>% 
#   filter(!time_period %in% out_time_periods)

# 
# out_fuzzy <-  sugarcube_data %>%
#   filter(agrepl("Kathrein Electronic AG", firmnames,
#                 max.distance=0.1)  &
#            time_period %in% gaps_time)
# 
# max_distance_input_fn <- 0.1
# max_distance_input_ln <- 0.1
# out_fuzzy <-  sugarcube_data %>%
#   filter(agrepl(df_gaps$firstname_sug1[i], firstname,
#                 max.distance=max_distance_input_fn) &
#            agrepl(df_gaps$lastname_sug1[i], lastname,
#                   max.distance = max_distance_input_ln) &
#            time_period %in% gaps_time)
# 
# View(out_fuzzy)

# rl_results %>% filter(id==38398 & time_period %in% gaps_time )
# df_after_post_processing %>% filter(id_1==38398 & time_period %in% gaps_time)
