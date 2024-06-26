
# (A) Load packages, settings, and set working directory ----

library(tidyverse)
library(readstata13)
library(janitor)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))

#rm(list = ls(all = TRUE))
source("./03_Code/13_Running_variable/00_Vote_margin_analytical_highestaverage.R") 


# (B) Read in data and calculate vote margin ----


# (i) Read in data and look at districts and years

data <- read.dta13("./02_Processed_data/nationalraete_1931_2015.dta") %>%
  filter(year>=1971)%>%
  as_tibble() %>%
  mutate(pvotes=ifelse(is.na(pvotes)==F,pvotes,votes))

data %>% filter(year==1971) %>% tabyl(listname_bfs)




# (ii) Calculate running variable for current situation

data_swi_prepared <- PrepareData(data %>% select(-votemargin),
                                 votes_h_name="votes", 
                                 votes_j_name= "pvotes",
                                 alliance_name="alliance",
                                 suballiance_name="suballiance",
                                 party_name="list",
                                 districtname="canton",
                                 election_cyclename = "year",
                                 system="open",
                                 cand_id_name="ID_pers",
                                 alliances=TRUE) %>% 
  filter(is.na(votes_h) == F) # Note: filter out tacit elections


RV_actual <- CalculateMargins(data_input=data_swi_prepared,
                               system="open",
                               method="dHondt",
                               additional_vars=c("ID_pers","year"),
                               return_option=T) %>%
  rename(votemargin_actual=votemargin)


# (iv) Calculate running variable if there were no alliances

data_swi_prepared <- PrepareData(data %>% select(-votemargin),
                                 votes_h_name="votes", 
                                 votes_j_name= "pvotes",
                                 alliance_name="alliance",
                                 suballiance_name="suballiance",
                                 party_name="list",
                                 districtname="canton",
                                 election_cyclename = "year",
                                 system="open",
                                 cand_id_name="ID_pers",
                                 alliances=FALSE) %>% 
  filter(is.na(votes_h) == F) # Note: filter out tacit elections

RV_no_alliance1 <- CalculateMargins(data_input=data_swi_prepared,
                               system="open",
                               method="dHondt",
                               additional_vars=c("ID_pers","year"),
                               return_option=T)%>%
  rename(votemargin_no_alliance1=votemargin)


data_swi_prepared <- PrepareData(data %>% select(-votemargin),
                                 votes_h_name="votes", 
                                 votes_j_name= "pvotes",
                                 alliance_name="alliance",
                                 suballiance_name="suballiance",
                                 party_name="list",
                                 districtname="canton",
                                 election_cyclename = "year",
                                 system="open",
                                 cand_id_name="ID_pers",
                                 alliances=TRUE) %>% 
  filter(is.na(votes_h) == F) # Note: filter out tacit elections

RV_no_alliance2 <- CalculateMargins(data_input=data_swi_prepared ,
                                    system="open",
                                    method="dHondt",
                                    additional_vars=c("ID_pers","year"),
                                    return_option=T)%>%
  rename(votemargin_no_alliance2=votemargin)



# (v) Calculate running variable if there are only alliances in the same party (intra)
# 
# # Note: remove all alliance and suballiances for übrige Parteien 
# 
# df <- data %>%
#   mutate(alliance=ifelse(listname_bfs %in% c("Übrige/Autres"),"",alliance))%>%
#   mutate(suballiance=ifelse(listname_bfs %in% c("Übrige/Autres"),"",suballiance))%>%
#   mutate(party_alliance=ifelse(!alliance %in% c("","."),paste0(alliance,"-",listname_bfs),""))%>%
#   mutate(party_suballiance=ifelse(!suballiance %in% c("","."),paste0(suballiance,"-",listname_bfs),"")) 
# 
# df_sa <- df %>% 
#   group_by(party_alliance,canton,year) %>%
#   summarize(sd_sa=sd(as.numeric(as.factor(party_suballiance)),na.rm=T)) %>%
#   filter(party_alliance!="")
# 
# df_al <- df %>% 
#   group_by(list,canton,year) %>%
#   summarize(sd_al=sd(as.numeric(as.factor(party_alliance)),na.rm=T)) %>%
#   filter(list!="")
# 
# 
# df_all <- df %>% 
#   left_join(df_sa,by=c("year","canton","party_alliance")) %>% 
#   left_join(df_al,by=c("year","canton","list"))
# 
# #View(df %>% select(contains("alliance"),listname_bfs,list,canton,year))
# 
# data_swi_prepared <- PrepareData(df_all %>% select(-votemargin) ,
#                                  votes_h_name="votes", 
#                                  votes_j_name= "pvotes",
#                                  alliance_name="party_alliance",
#                                  suballiance_name="party_suballiance",
#                                  party_name="list",
#                                  districtname="canton",
#                                  election_cyclename = "year",
#                                  system="open",
#                                  cand_id_name="ID_pers",
#                                  alliances=TRUE) %>% 
#   filter(is.na(votes_h) == F)%>%
#   mutate(alliance_dummy=ifelse(sd_al>0 & !is.na(sd_al),alliance_dummy,0 ))%>%
#   mutate(suballiance_dummy=ifelse(sd_sa>0 & !is.na(sd_sa),suballiance_dummy,0 ))
# 
#   
# RV_only_intra <- CalculateMargins(data_input=data_swi_prepared %>% filter(year==2015 & canton=="ZH"),
#                                system="open",
#                                method="dHondt",
#                                additional_vars=c("ID_pers","year"),
#                                return_option=T)%>%
#   rename(votemargin_only_intra=votemargin)


data_swi_out <- data_swi_prepared %>% filter(year==2015) %>% 
  left_join(RV_actual %>%
              select(ID_pers,year,votemargin_actual), 
            by=c("ID_pers","year")) %>% 
  left_join(RV_no_alliance1 %>%
              select(ID_pers,year,votemargin_no_alliance1), 
            by=c("ID_pers","year")) %>% 
  left_join(RV_no_alliance2 %>%
              select(ID_pers,year,votemargin_no_alliance2), 
            by=c("ID_pers","year")) %>% 
  left_join(RV_only_intra %>%
              select(ID_pers,year,votemargin_only_intra), 
            by=c("ID_pers","year")) %>%
  mutate(elected_actual=ifelse(votemargin_actual>0,1,0),
         elected_no_alliance1=ifelse(votemargin_no_alliance1>0,1,0),
         elected_no_alliance2=ifelse(votemargin_no_alliance2>0,1,0),
         elected_only_intra=ifelse(votemargin_only_intra>0,1,0)) 


# Two tests

data_swi_out %>% filter(elected!=elected_actual)

data_swi_out %>% filter(elected_no_alliance1!=elected_no_alliance2)


# (C) Analyze data ----

data_swi_out_col <- data_swi_out %>%
  group_by(year,listname_bfs) %>%
  summarize(elected_actual=sum(elected_actual,na.rm=T),
            elected_no_alliance1=sum(elected_no_alliance1,na.rm=T),
            elected_only_intra=sum(elected_only_intra,na.rm=T))

sum(data_swi_out_col$elected_actual)
sum(data_swi_out_col$elected_only_intra)

data_swi_out  %>%
  group_by(vm_pos) %>%
  summarize(el_mean=mean(elected))

data_swi_out %>% tabyl(vm_pos,elected)

