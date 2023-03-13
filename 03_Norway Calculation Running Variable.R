#-------------------------------------------------------
# (A) Load packages, settings, and set working directory
#-------------------------------------------------------

library(tidyverse)
library(readstata13)
library(janitor)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))

#rm(list = ls(all = TRUE))
source("./03_Code/fine_grid.R") # ggplot layers
#source("./03_Code/13_Running_variable/04_Vote_margin_analytical_highestaverage.R") 

#--------------------------------------------
# (B) Read in data and calculate vote margin
#--------------------------------------------

# (i) Read in data and look at districts and years

data <- read.dta13("./01_Raw_data/13_Running_variable/FivaSmithJune2019.dta")%>% 
    select(year,districtid, district, electorate, party, partyname, rank, pid, 
           votes,elected) %>% 
    filter(year>=1953 & year<=1981) %>% 
    arrange(districtid,year, partyname,pid)

data %>% tabyl(district,year) 

# (ii) Prepare data and calculate vote margin

data_norway_prep <- PrepareData(data=data, 
                                votes_j_name="votes",
                                party_name="partyname", 
                                districtname="district", 
                                election_cyclename = "year",
                                alliances=F,
                                system="closed") 

norway_margins <- CalculateMargins(data_norway_prep,
                                   system="closed",
                                   method="SainteLague",
                                   rank="rank",
                                   margin_type="easiest",
                                   additional_vars=c("partyname",
                                                     "pid",
                                                     "district",
                                                     "districtid"))

norway_margins <- norway_margins %>% 
                  mutate(party_num=party) %>%
                  select(district,districtid,year,party_num,partyname,
                         votemargin,pid)

# (iii) Merge results on vote margin back to initial dataset

data_final_norway_elected <- data %>% 
  group_by(year,districtid )%>%
  summarize(elected_sum=sum(elected))

data_final_norway <- data %>% 
                     left_join(norway_margins %>% select(-districtid,),
                               by=c("district","year","partyname","pid")) %>%
                     left_join(data_final_norway_elected,
                               by=c("districtid","year"))

data.table::fwrite(data_final_norway, 
file = "./02_Processed_data/13_Running_variable/RV_Norway_1953_1981_Easiest.csv", 
sep = ";",row.names = F)

# -------------------
# (C) Quality Checks
#--------------------

# (i) Check whether all elected have positive margin and all non-elected have 
#     negative margin

summary(data_final_norway %>% filter(rank<=elected_sum & elected==1) %>% 
          select(votemargin))
summary(data_final_norway %>% filter(rank<=elected_sum & elected==0) %>% 
          select(votemargin))

data_final_norway %>% filter(rank<=elected_sum & elected==0 & votemargin>0)
data_final_norway %>% filter(rank<=elected_sum & elected==1 & votemargin<0)
# result: all fine

# (ii) Check whether candidates with a rank higher than the number of seats have 
#    a missing votemargin

summary(data_final_norway %>% filter(rank>elected_sum) %>% select(votemargin))
# result: all fine, for these cases we find only na's

# (iii) Check whether we can replicate the election results

seatsnorway_predicted <- GetSeatsHighestAverageLoop(data_input = 
                                                    data_norway_prep,
                           method="SainteLague") %>%
                           as_tibble() %>%
                           rename(seats_predicted=seats)
            
seatsnorway_predicted$seats_predicted <- as.numeric(
                                         seatsnorway_predicted$seats_predicted)

seatsnorway_official <-   data_norway_prep %>% 
                          group_by(year,districtid,party) %>% 
                          summarize(seats_official=sum(elected)) %>%
                          mutate(year=as.character(as.factor(year)),
                                 districtid=as.character(as.factor(districtid)),
                                 party=as.character(as.factor(party)))

seatsnorway_compare <- seatsnorway_predicted %>% left_join(seatsnorway_official,
                                    by=c("year","districtid","party"))

seatsnorway_compare %>% mutate(seats_diff=seats_official-seats_predicted) %>%
                        filter(seats_diff!=0)
# result: all fine, no difference

# (iv) Compare simulated vote margin to analytic vote margin
start_time <- Sys.time()
Simout <- GetRVSimulation_Candidate_Loop(data_input=data_norway_prep,
                                     convcrit=0.001,
                                     method="SainteLague",
                                     system="closed",
                                     type="Votes_Added")
end_time <- Sys.time()
end_time - start_time
Simout <- Simout %>% 
          mutate(rank=seat) %>%
          rename(party_num=party,
                 votemargin_sim=VS)
Sim_Compare <- Simout %>% left_join(data_final_norway,by=c("party_num","rank")) %>%
           mutate(votemargin_diff=votemargin_sim-votemargin)%>%
           select(year,districtid, district, partyname,pid,votemargin,
                  votemargin_sim,votemargin_diff)

summary(Sim_Compare$votemargin_diff)
data.table::fwrite(Sim_Compare, 
                   file = "./02_Processed_data/13_Running_variable/RV_Norway_1953_1981.csv", 
                   sep = ";",row.names = F)

# --------------------------------------------------------
# (D) Simulated margin for vote redistribution among parties
#---------------------------------------------------------

# (i) Run simulation

start_time <- Sys.time()
Simout1953 <- GetRVSimulation_Candidate_Loop(data_input=data_norway_prep %>%
                                           filter(year==1953),
                                         convcrit=0.001,
                                         method="SainteLague",
                                         system="closed",
                                         type="Party_redistribution")
end_time <- Sys.time()
end_time - start_time
Simout <- Simout %>% 
  mutate(rank=seat) %>%
  rename(party_num=party,
         votemargin_sim=VS)

# (ii) Quality checks

data_compare <- data_norway_prep %>% left_join(Simout,by=c("party","rank")) %>%
  arrange(districtid,party,rank) %>% 
  select(district,districtid,year,partyname,party,partyminusj,rank,
         votes_j,Votes_required,votemargin_sim,votes_range_max) %>%
  filter(Votes_required<=votes_range_max)

# data_norway_prep %>% left_join(Simout,by=c("party","rank")) %>%
#   arrange(districtid,party,rank) %>% 
#   select(district,districtid,year,partyname,party,partyminusj,rank,
#          votes_j,Votes_required,votemargin_sim,votes_range_max) %>%
#   filter(Votes_required>votes_range_max & districtid==1)

data_compare %>%  distinct(party,Votes_required,.keep_all = T) %>%
  filter(districtid==1 & votes_j==58120)

data_compare %>% distinct(districtid,year,partyname,party)%>%
  filter(districtid==1 & year==1953)

# (iii) Save original file

Simout <- Simout %>% mutate(party=as.numeric(as.character(party)))
data.table::fwrite(Simout, file = 
                     "./02_Processed_data/13_Running_variable/RV_Norway_1953_1981_Votes_redistr.csv", 
                  sep = ";",row.names = F)

# (iv) Aggregate to single margin and save

data_margin_final <- data_compare %>% 
                 group_by(districtid,year,partyname,party,rank) %>%
                 dplyr::summarize(votemargin_sim_min=min(votemargin_sim),
                                  votemargin_sim_max=max(votemargin_sim)) %>%
                 ungroup()
data_margin_final %>% filter(votemargin_sim_min<0 & votemargin_sim_max>0)  
# result: no case with positive and negative margin at the same time

data_margin_final <- data_margin_final %>% 
                     mutate(votemargin_redis=ifelse(votemargin_sim_min>=0,
                             votemargin_sim_min,
                             votemargin_sim_max)) %>%
                     select(districtid,year,partyname,rank,
                            votemargin_redis)


data.table::fwrite(data_margin_final, file = 
  "./02_Processed_data/13_Running_variable/RV_Norway_1953_1981_redis.csv", 
   sep = ";",row.names = F)


# --------------------------------------------------------
# (E) Simulated margin for easiest way as in Fiva and Smith (2018)
#---------------------------------------------------------

# (i) Run simulation

start_time <- Sys.time()
Simout1953 <- GetRVSimulation_Candidate_Loop(data_input=data_norway_prep %>%
                                               filter(year==1953),
                                             convcrit=0.001,
                                             method="SainteLague",
                                             system="closed",
                                             type="Votes_Subtracted_Or_Added")
end_time <- Sys.time()
end_time - start_time
Simout <- Simout %>% 
  mutate(rank=seat) %>%
  rename(party_num=party,
         votemargin_sim=VS)

# (ii) Quality checks

data_compare <- data_norway_prep %>% left_join(Simout,by=c("party","rank")) %>%
  arrange(districtid,party,rank) %>% 
  select(district,districtid,year,partyname,party,partyminusj,rank,
         votes_j,Votes_required,votemargin_sim,votes_range_max) %>%
  filter(Votes_required<=votes_range_max)

# data_norway_prep %>% left_join(Simout,by=c("party","rank")) %>%
#   arrange(districtid,party,rank) %>% 
#   select(district,districtid,year,partyname,party,partyminusj,rank,
#          votes_j,Votes_required,votemargin_sim,votes_range_max) %>%
#   filter(Votes_required>votes_range_max & districtid==1)

data_compare %>%  distinct(party,Votes_required,.keep_all = T) %>%
  filter(districtid==1 & votes_j==58120)

data_compare %>% distinct(districtid,year,partyname,party)%>%
  filter(districtid==1 & year==1953)

# (iii) Save original file

Simout <- Simout %>% mutate(party=as.numeric(as.character(party)))
data.table::fwrite(Simout, file = 
                     "./02_Processed_data/13_Running_variable/RV_Norway_1953_1981_Votes_subtracted_or_added.csv", 
                   sep = ";",row.names = F)

# (iv) Aggregate to single margin and save

data_margin_final <- data_compare %>% 
  group_by(districtid,year,partyname,party,rank) %>%
  dplyr::summarize(votemargin_sim_min=min(votemargin_sim),
                   votemargin_sim_max=max(votemargin_sim)) %>%
  ungroup()
data_margin_final %>% filter(votemargin_sim_min<0 & votemargin_sim_max>0)  
# result: no case with positive and negative margin at the same time

data_margin_final <- data_margin_final %>% 
  mutate(votemargin_redis=ifelse(votemargin_sim_min>=0,
                                 votemargin_sim_min,
                                 votemargin_sim_max)) %>%
  select(districtid,year,partyname,rank,
         votemargin_redis)


data.table::fwrite(data_margin_final, file = 
                     "./02_Processed_data/13_Running_variable/RV_Norway_1953_1981_redis.csv", 
                   sep = ";",row.names = F)

#--------------------------------------------------
# (F) Evaluate correlation of different votemargins 
# (our vote margin, easiest, redistribution)
#--------------------------------------------------

# (i) Method proposed by Lüchinger, Schelker, Schmid (2021)

df_std <- data.table::fread(
  file = "./02_Processed_data/13_Running_variable/RV_Norway_1953_1981.csv", 
  sep = ";")

# (ii) Simulated margin for easiest way as in Fiva and Smith (2018)

df_easiest <- data.table::fread(
  file = "./02_Processed_data/13_Running_variable/RV_Norway_1953_1981_Easiest.csv", 
  sep = ";")


# (iii) Simulated margin for vote distribution among parties

df_redis <- data.table::fread(
  file ="./02_Processed_data/13_Running_variable/RV_Norway_1953_1981_redis.csv",
  sep = ";")

# (iv) Mege all data and compare margins

data <- read.dta13("./01_Raw_data/13_Running_variable/FivaSmithJune2019.dta")%>% 
  select(year,districtid, district, electorate, party, partyname, rank, pid, 
         votes,elected) %>% 
  filter(year>=1953 & year<=1981) %>% 
  arrange(districtid,year, partyname,pid)

data %>% tabyl(district,year) 

df_all <- data %>%
  left_join(df_std, by=c("year","districtid","pid"))%>% 
  rename(votemargin_std=votemargin) %>%
  left_join(df_easiest %>% select(year,districtid,pid,votemargin), 
            by=c("year","districtid","pid")) %>%
  rename(votemargin_easiest=votemargin) %>%
  left_join(df_redis %>% select(year,districtid,pid,votemargin), 
            by=c("year","districtid","rank")) %>%
  rename(votemargin_redis=votemargin)

cor(df_all %>% select(starts_with("votemargin")))

# (i) Method proposed by Lüchinger, Schelker, Schmid (2021)

df_std <- data.table::fread(
  file = "./02_Processed_data/13_Running_variable/RV_Norway_1953_1981.csv", 
  sep = ";")

# (ii) Simulated margin for easiest way as in Fiva and Smith (2018)

df_easiest <- data.table::fread(
  file = "./02_Processed_data/13_Running_variable/RV_Norway_1953_1981_Easiest.csv", 
  sep = ";")


# (iii) Simulated margin for vote distribution among parties

df_redis <- data.table::fread(
  file ="./02_Processed_data/13_Running_variable/RV_Norway_1953_1981_redis.csv",
  sep = ";")

# (iv) Mege all data and compare margins

data <- read.dta13("./01_Raw_data/13_Running_variable/FivaSmithJune2019.dta")%>% 
  select(year,districtid, district, electorate, party, partyname, rank, pid, 
         votes,elected) %>% 
  filter(year>=1953 & year<=1981) %>% 
  arrange(districtid,year, partyname,pid)

data %>% tabyl(district,year) 

df_all <- data %>%
  left_join(df_std, by=c("year","districtid","pid"))%>% 
  rename(votemargin_std=votemargin) %>%
  left_join(df_easiest %>% select(year,districtid,pid,votemargin), 
            by=c("year","districtid","pid")) %>%
  rename(votemargin_easiest=votemargin) %>%
  left_join(df_redis %>% select(year,districtid,pid,votemargin), 
            by=c("year","districtid","rank")) %>%
  rename(votemargin_redis=votemargin)

cor(df_all %>% select(starts_with("votemargin")))



# check <- data_final_norway %>% filter(rank<=elected_sum & elected==0 & 
#                                        votemargin>0|
#                                         rank<=elected_sum & elected==1 & 
#                                       votemargin<0)%>% 
#   arrange(year,districtid,party)
# 
# check_2 <-data_final_norway %>% filter(year==1973&districtid==18)
# check_2 %>%  group_by(partyname) %>% 
#   summarize(elected_sum=sum(elected),
#             votes=first(votes)) 
# 
# View(check_2)
# data_final_norway%>%group_by(year,district,partyname,party) %>% 
#   summarize(elected_sum=sum(elected),
#             votes=first(votes)) %>%
#   filter(district=="nordland"& year==1973)
# 
# data_final_norway%>%filter(district=="nordland" & year==1973)%>% 
#   summarize(elected_sum=sum(elected),
#             votes=first(votes)) 
# data_final_norway%>%
#   filter(district=="nordland" & year==1973 & 
#          partyname=="det norske arbeiderparti")
# 
# data_final_norway%>%filter(partyname=="høyre og kristelig folkeparti")
# 
# Check simulation
# data_input <- PrepareData(data=data, votes_j_name="votes",
#                           party_name="partyname", districtname="district", 
#                           election_cyclename = "year",alliances=F,
#                           system="closed") %>%
#                           filter(district=="nordland" & year==1973)

# Simout <- GetRVSimulation_Candidate(data_input=data_norway_prep, 
#                                     n=sum(test$elected), 
#                                     convcrit=0.001,
#                                     method="SainteLague",
#                                     system="closed") 
# end_time <- Sys.time()
# end_time - start_time
# Simout <- Simout %>% 
#   mutate(rank=seat) %>%
#   rename(party_num=party,
#          votemargin_sim=VS)
# Sim_Compare <- Simout %>% left_join(data_final_norway,by=c("party_num","rank")) %>%
#   mutate(votemargin_diff=votemargin_sim-votemargin)%>%
#   select(year,districtid, district, partyname,pid,votemargin,
#          votemargin_sim,votemargin_diff)
