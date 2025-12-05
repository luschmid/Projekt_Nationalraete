#-----------------------------------------------------------
# (A) Libraries, paths, and load functions (from rl project)
#----------------------------------------------------------

library(tidyverse)
library(data.table)
library(dplyr, warn.conflict = FALSE, quietly = TRUE)
library(lubridate)
library(testit)
library(haven)
library(stringr)

rm(list = ls())
rm(list=setdiff(ls(all.names=TRUE), lsf.str(all.names=TRUE)))

setwd("E:/Projekt Nationalräte")
source("./03_Code/12_Record_linkage/00_Functions_Record_linkage.R")
source("./03_Code/12_Record_linkage/00_Functions.R")
path <- "./02_Processed_data/12_Record_linkage/03_Data_out_analysis/"

#--------------------------------
# B) Bisnode Gen 7: optimal cutoff
#--------------------------------

# (i) Set parameters

modelname="Generation_7_corrected";
year_range=c(1994:2017);
source="bisnode";
linkscore_cutoff="optimal";

# (ii) Load correspondence table

CorrespondenceTable <- GetCorrespondenceTable(modelname=modelname,
                                              source=source,
                                              linkscore_cutoff=linkscore_cutoff,
                                              corrected=T) %>%
  arrange(id_0) %>% 
  distinct(id_0, id_1)


# (iii) Load all elected candidates and bisnode data

nr_data_elected_1999_2015 <- read_dta("./02_Processed_data/02_Elections_1971_2015/nationalraete_1931_2015_Geo.dta") %>%
  filter((year==1999 | year==2015) & elected==1)

id_1_correspondencetable <- CorrespondenceTable %>% 
  filter(id_0 %in% nr_data_elected_1999_2015$ID) %>%
  select(id_1) %>%
  distinct() %>% 
  pull()

data_input <- data.table::fread("./02_Processed_data/11_Directors_1994_2018/Bisnode_Person-Firmen_Geo_Dataconstruct.csv",
                                    sep = ";", encoding = "UTF-8") %>%  
      filter(personenid %in% as.numeric(id_1_correspondencetable))


# (iv) Construct out data


year_1999_election <- 2002
year_2015_election <- 2016

data_input_cut <- BisnodeDataCut(data=data_input,
                                    id_bis=id_1_correspondencetable,
                                    startdate="1994-01-01",
                                    enddate="2017-12-31",
                                    gremium_categories=c("Verwaltungsrat","Bankrat"),
                                    cut_month=6,
                                    filter_conditions=filter_conditions,
                                    president=FALSE) %>%
  arrange(personenid,duns) %>%
  mutate(id_1=personenid)%>%
  left_join(CorrespondenceTable,by=c("id_1")) %>%
  left_join(nr_data_elected_1999_2015 %>% 
              dplyr::rename(id_0=ID,election_year=year) %>%
              select(name,firstname,canton,election_year,id_0,job,list),
            by=c("id_0")) %>%
  filter((eintrittjahr<=year_1999_election & austrittjahr>=year_1999_election & election_year==1999)| 
          (eintrittjahr<=year_2015_election & austrittjahr>=year_2015_election & election_year==2015)|
           is.na(eintrittjahr)) %>%
  left_join(data_input %>%
              mutate(id_1=personenid)%>%
              select(id_1,vorname,nachname,handelsregisternummer,duns,uid,firma) %>%
              distinct(),
            by=c("id_1","duns")) %>%
  mutate(year=ifelse(election_year==1999,year_1999_election,year_2015_election))%>%
  rename(firstname_pol=firstname,lastname_pol=name,lastname_bis=nachname,firstname_bis=vorname)%>%
  select(id_0,id_1,year,firstname_pol,lastname_pol,firstname_bis,lastname_bis,firma,duns,handelsregisternummer,uid) %>%
  distinct() %>%
  arrange(id_0,year,id_1,duns) %>%
  left_join(nr_data_elected_1999_2015 %>% dplyr::rename(id_0=ID), ,by=c("id_0"))


openxlsx::write.xlsx(data_input_cut, file = "./02_Processed_data/12_Record_linkage/03_Data_analysis/bisnode_registers_Generation7_optimal.xlsx")


# Quality check: check whether those with no mandate are not in record linkage correspondence table

nr_no_mandate <- data_input_cut %>% filter(is.na(firma)) %>% select(id_0) %>% pull()
  
no_mand <- CorrespondenceTable %>% filter(id_0 %in% nr_no_mandate)

data_input %>% filter(personenid %in% no_mand$id_1)

# result: these are all mandates with starting date in 2017 and later

# Quality check: 

id_0_correspondencetable <- CorrespondenceTable %>% 
  filter(id_0 %in% nr_data_elected_1999_2015$ID) %>%
  select(id_0) %>%
  distinct() %>% 
  pull()

data_input_cut %>%
  filter(!id_0 %in% id_0_correspondencetable)
