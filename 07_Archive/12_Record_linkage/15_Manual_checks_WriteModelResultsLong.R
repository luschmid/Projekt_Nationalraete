#-----------------------------------------------------------
# (A) Libraries, paths, and load functions (from rl project)
#----------------------------------------------------------

library(tidyverse)
library(data.table)
library(dplyr, warn.conflict = FALSE, quietly = TRUE)
library(lubridate)
library(testit)
library(haven)

rm(list = ls())

setwd("E:/Projekt Nationalräte")
source("./03_Code/12_Record_linkage/00_Functions_Record_linkage.R")
source("./03_Code/12_Record_linkage/00_Functions.R")

#--------------------------------------------------------
# (B) Bisnode Generation 7: optimal cutoff, November 2022
#--------------------------------------------------------

modelname_bisnode <- "Generation_7"
source="bisnode"
linkscore_cutoff <- "optimal"
sample_size <- 100

# (i) Load NR

nr_data <- read_dta("./02_Processed_data/nationalraete_1931_2015.dta")

id_all <- nr_data  %>%
  distinct(ID) %>%
  arrange(ID) %>%
  pull()

# (ii) Load Bisnode


# (a) Correspondence table


CorrespondenceTable <- GetCorrespondenceTable(modelname=modelname_bisnode,
                                              source=source,
                                              linkscore_cutoff=linkscore_cutoff) %>%
  arrange(id_0) %>% 
  distinct(id_0,id_1)


# (b) Random sample of correspondence table

id_0_correspondencetable <- CorrespondenceTable %>% select(id_0) %>%
  distinct()%>% 
  pull()


# (c) Correspondence table only for random sample

set.seed(123)

id_0_randomsample <- unique(sample(id_0_correspondencetable,sample_size))

CorrespondenceTable_Cut <- CorrespondenceTable %>%
  filter(id_0 %in% id_0_randomsample)

id_1_randomsample <- unique(CorrespondenceTable_Cut$id_1)

# (d) Load full bisnode data
# Note: As discussed, we focus on rechtsform=="Aktiengesellschaft" & gremium=="Verwaltungsrat"



data_input <- data.table::fread("./02_Processed_data/11_Directors_1994_2018/Bisnode_Person-Firmen_Geo_Dataconstruct.csv",
                                sep = ";", encoding = "UTF-8") %>%  
  dplyr::filter(rechtsform=="Aktiengesellschaft" & gremium=="Verwaltungsrat" ) 
  

data_input %>% filter(is.na(eintrittdatum) & is.na(austrittdatum))
# what to do about them?

# other question: duplicate entries

# (e) Make same transformation as in code (function: QueryPoliticianMandates)

director_data_cut <- BisnodeDataCut(data=data_input,id_bis=CorrespondenceTable_Cut$id_1,
                                    startdate="1994-01-01",
                                    enddate="2017-12-31",
                                    gremium_categories=c("Verwaltungsrat","Bankrat"),
                                    cut_month=6,
                                    filter_conditions="",
                                    president=FALSE) %>%
  mutate(personenid=as.character(personenid)) %>%
  left_join(CorrespondenceTable_Cut,by=c("personenid"="id_1"))

director_data_cut %>% 
  get_dupes(personenid,duns,eintrittjahr,austrittjahr)

# (f) Generate outfile

openxlsx::write.xlsx(director_data_cut, file = "./02_Processed_data/12_Record_linkage/01_Bisnode/Check_WriteModelResultsLong_bisnode_generation7.xlsx")


# (g) Save long format template for random sample

all_checks <- director_data_cut %>% expand(id_0, (1994:2017)) %>%
  dplyr::rename(year=`(1994:2017)`) %>%
  arrange(id_0,year)


openxlsx::write.xlsx(all_checks, file = "./02_Processed_data/12_Record_linkage/01_Bisnode/Check_WriteModelResultsLong_bisnode_generation7_template.xlsx")

  all_checks %>% 
  group_by(id_0) %>%
  summarize(nobs=n()) %>%
  filter(nobs!=24)

# Result: all fine. Only 55 obs b/c some of the 45 dropped have mandates that
# are out of our focus (Stiftungsrat etc)


  
  

#-----------------------------------------
# (B) Bisnode Generation 7, November 2022
#-----------------------------------------

modelname_bisnode <- "Generation_7"
source="bisnode"
linkscore_cutoff <- "optimal"
sample_size <- 100

# (i) Load NR

nr_data <- read_dta("./02_Processed_data/nationalraete_1931_2015.dta")

id_all <- nr_data  %>%
  distinct(ID) %>%
  arrange(ID) %>%
  pull()

# (ii) Load Bisnode


# (a) Correspondence table


CorrespondenceTable <- GetCorrespondenceTable(modelname=modelname_bisnode,
                                              source=source,
                                              linkscore_cutoff=linkscore_cutoff) %>%
  arrange(id_0) %>% 
  distinct(id_0,id_1)


# (b) Random sample of correspondence table

id_0_correspondencetable <- CorrespondenceTable %>% select(id_0) %>%
  distinct()%>% 
  pull()


# (c) Correspondence table only for random sample

set.seed(123)

id_0_randomsample <- unique(sample(id_0_correspondencetable,sample_size))

CorrespondenceTable_Cut <- CorrespondenceTable %>%
  filter(id_0 %in% id_0_randomsample)

id_1_randomsample <- unique(CorrespondenceTable_Cut$id_1)

# (d) Load full bisnode data
# Note: As discussed, we focus on rechtsform=="Aktiengesellschaft" & gremium=="Verwaltungsrat"



data_input <- data.table::fread("./02_Processed_data/11_Directors_1994_2018/Bisnode_Person-Firmen_Geo_Dataconstruct.csv",
                                sep = ";", encoding = "UTF-8") %>%  
  dplyr::filter(rechtsform=="Aktiengesellschaft" & gremium=="Verwaltungsrat" ) 


data_input %>% filter(is.na(eintrittdatum) & is.na(austrittdatum))
# what to do about them?

# other question: duplicate entries

# (e) Make same transformation as in code (function: QueryPoliticianMandates)

director_data_cut <- BisnodeDataCut(data=data_input,id_bis=CorrespondenceTable_Cut$id_1,
                                    startdate="1994-01-01",
                                    enddate="2017-12-31",
                                    gremium_categories=c("Verwaltungsrat","Bankrat"),
                                    cut_month=6,
                                    filter_conditions="",
                                    president=FALSE) %>%
  mutate(personenid=as.character(personenid)) %>%
  left_join(CorrespondenceTable_Cut,by=c("personenid"="id_1"))

director_data_cut %>% 
  get_dupes(personenid,duns,eintrittjahr,austrittjahr)

# (f) Generate outfile

openxlsx::write.xlsx(director_data_cut, file = "./02_Processed_data/12_Record_linkage/01_Bisnode/Check_WriteModelResultsLong_bisnode_generation7.xlsx")


# (g) Save long format template for random sample

all_checks <- director_data_cut %>% expand(id_0, (1994:2017)) %>%
  dplyr::rename(year=`(1994:2017)`) %>%
  arrange(id_0,year)


openxlsx::write.xlsx(all_checks, file = "./02_Processed_data/12_Record_linkage/01_Bisnode/Check_WriteModelResultsLong_bisnode_generation7_template.xlsx")

all_checks %>% 
  group_by(id_0) %>%
  summarize(nobs=n()) %>%
  filter(nobs!=24)

# Result: all fine. Only 55 obs b/c some of the 45 dropped have mandates that
# are out of our focus (Stiftungsrat etc)

# (ii) Read in manual checks by Oester and Marilley and consolidates by 
#      Mark and Simon, Email Mark 19 November 2022

library(readxl)
df_cor <- read_excel("./02_Processed_data/12_Record_linkage/04_Post_RL_QualChecks/bisnode_fp_firstcheck_Generation7_optimal_corrected.xlsx") %>%
  filter(is.na(CodierungFin)) %>%
  group_by(id_0,id_1) %>%
  mutate(age_firstmandate_max=max(age_firstmandate,na.rm=T),
         age_firstmandate_min=min(age_firstmandate,na.rm=T),
         age_lastmandate_max=max(age_lastmandate,na.rm=T))%>%
  select(id_0,id_1,starts_with("age")) %>%
  filter(age_firstmandate_max<=85 & age_lastmandate_max<=100 &
           age_firstmandate_min>=18  ) # see Email 19 November 2022

data.table::fwrite(
  x = df_cor,
  file = "./02_Processed_data/12_Record_linkage/01_Bisnode/RL-Results-generation_7_corrected.csv",
  sep = ";")
