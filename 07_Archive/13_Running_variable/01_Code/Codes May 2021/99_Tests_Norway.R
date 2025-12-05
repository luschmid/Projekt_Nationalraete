library(tidyverse)
library(readstata13)
library(janitor)
library(scales)
library(ggExtra)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))

data <- read.dta13("./02_Processed_data/13_Running_variable/FivaSmithJune2019_processed.dta") %>% 
  select(year,districtid, district, party, partyname, pid, votes,voteshare,
         elected_F1,rank,votemargin_vote,margin,elected) %>% 
  filter(year>=1953 & year<=1981) %>% 
  arrange(districtid,year, partyname,pid) %>%
  mutate(fiva_smith_sample=ifelse(is.na(margin)==F,"Included in Fiva & Smith (2018)",
                                  "Excluded in Fiva & Smith (2018)"))
data_district_party <- data %>%
  group_by(year,districtid,party) %>%
  summarize(pvotes=first(votes),
            elected_sum=sum(elected))
votes <- data_district_party[c(1:6),c(3,4)]

method <- "SainteLague"; no_seats <- 8; 

data_input_test <- CalculateRatios(data_district_party[c(1:6),c(3,4)],8,method="SainteLague")
PartyMargin(party="v", data_input=data_input_test, no_seats=8,method="SainteLague")

data <- read.dta13("./01_Raw_data/13_Running_variable/FivaSmithJune2019.dta") %>% 
  select(year,districtid, district, party, partyname, pid, votes,voteshare,
         elected,rank) %>% 
  filter(year>=1953 & year<=1981) %>% 
  arrange(districtid,year, partyname,pid)


data_1_1953 <- data %>% filter(districtid == 1 & year == 1953) %>%
               mutate(voteshare=round(voteshare,3))
test <- PrepareData(data=data_1_1953, votes_j_name="votes",party_name="party", 
                    districtname="district", election_cyclename = "year",alliances=F,
                    system="closed") 

AllianceInformation <- GetAllianceSuballiance(1, test) # get info whether party is in sub-/alliance and in which
test2 <- CalculatePartyMargin(party_input=1, data_input=test, 8,method="SainteLague")

CalculatePartyMargins(data_input=test,no_seats=8,method="SainteLague")

AggregateMargins(data_input=test, no_seats=8,system="closed",method="SainteLague",rank="rank",
                 additional_vars=c("partyname","pid")) 
  

data_norway_prep <- PrepareData(data=data, votes_j_name="votes",party_name="party", districtname="district", election_cyclename = "year",alliances=F,
                    system="closed") 


norway_margins <- CalculateMargins(data=test,system="closed",method="SainteLague",rank="rank",
                                   additional_vars=c("partyname","pid"))