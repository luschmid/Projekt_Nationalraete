#-----------------------------------------
# (A) Libraries, paths, and load functions
#-----------------------------------------

library(tidyverse)
library(data.table)
library(testit)
library(RSQLite)
library(odbc)
library(dbplyr)
library(dplyr, warn.conflict = FALSE, quietly = TRUE)
library(lubridate)
library(haven)

rm(list=setdiff(ls(all.names=TRUE), lsf.str(all.names=TRUE)))
setwd("E:/Projekt Nationalräte")
source("./03_Code/12_Record_linkage/00_Functions_Record_linkage.R")
db_path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte/02_Processed_data/database_political_rents.sqlite"

#----------------------
# B) Set parameters
#----------------------

modelname="Generation_7";
year_range=c(1994:2017);
source="bisnode";
batch_size=1000;
batch_number=1;
linkscore_cutoff="optimal";
path="./11_Directors_1994_2018/"

#--------------
# D) Load data
#--------------

CorrespondenceTable <- GetCorrespondenceTable(modelname=modelname,
                                              source=source,
                                              linkscore_cutoff=linkscore_cutoff) %>%
  arrange(id_0) %>% 
  distinct(id_0,id_1)

id_0_correspondencetable <- CorrespondenceTable %>% select(id_0) %>%
  distinct()%>% 
  pull()

id_1_correspondencetable <- CorrespondenceTable %>% select(id_1) %>%
  distinct() %>% 
  pull()

nr_data <- read_dta("./02_Processed_data/nationalraete_1931_2015.dta")

id_all <- nr_data  %>%
  distinct(ID) %>%
  arrange(ID) %>%
  pull()

data_input <- data.table::fread("./02_Processed_data/11_Directors_1994_2018/Bisnode_Person-Firmen_Geo_Dataconstruct.csv",
                                sep = ";", encoding = "UTF-8") %>%  
  mutate(id_0=as.character(personenid)) %>% 
  filter(rechtsform=="Aktiengesellschaft" & gremium=="Verwaltungsrat" ) %>%
  filter(personenid %in% id_1_correspondencetable)




###########################

# b) Manual checks Bisnode

CorrespondenceTable <- GetCorrespondenceTable(modelname="Generation_7",
                                              source=source,
                                              linkscore_cutoff=linkscore_cutoff) %>%
  arrange(id_0) %>% 
  distinct(id_0,id_1)

id_0_correspondencetable <- CorrespondenceTable %>% select(id_0) %>%
  distinct()%>% 
  pull()


N=100
set.seed(123)
out <- data.frame()
for (i in 1:N){
  out <- bind_rows(out,data.frame(id_0=sample(id_0_correspondencetable,1),year=sample(1994:2017,1))) 
}
out <- out %>% arrange(id_0,year)


CorrespondenceTable_sample <- CorrespondenceTable %>% 
  filter(id_0 %in% out$id_0)

id_1_correspondencetable <- CorrespondenceTable_sample
select(id_1) %>%
  distinct() %>% 
  pull()

my_db_connect <- dbConnect(RSQLite::SQLite(), db_path) 
data_to_check <- dbGetQuery(my_db_connect, "SELECT * FROM bisnode") %>% 
  filter(personenid %in% id_1_correspondencetable) %>%
  select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum)
dbDisconnect(my_db_connect)


path <- "C:/Schmidlu/Dropbox/Projekt NationalrÃ¤te/02_Processed_data/11_Directors_1994_2018/"
data.table::fwrite(x = out,file = paste0(path,"Manual_Check_Random_Sample",".csv"), sep = ";")
data.table::fwrite(x = CorrespondenceTable_sample,file = paste0(path,"Manual_Check_Correspondence_Table",".csv"),sep = ";")
data.table::fwrite(x = data_to_check %>% arrange(personenid),file = paste0(path,"Manual_Check_Bisnode_Data",".csv"),sep = ";")


# c) Check specific cases

pol_id <- "AG-2011-0185"
pol_id=c("ZH-1979-0051") # Christoph Blocher

results <- data.table::fread("C:/Schmidlu/Dropbox/Projekt NationalrÃ¤te/02_Processed_data/11_Directors_1994_2018/Bisnode_Results_Generation_4_optimal.csv",
                             sep = ";", encoding = "UTF-8")%>%
  filter(ID==pol_id)

CorrespondenceTable <- GetCorrespondenceTable(modelname="Generation_4",
                                              source="bisnode",
                                              linkscore_cutoff="optimal") %>%
  arrange(id_0) %>%
  filter(id_0==pol_id)

my_db_connect <- dbConnect(RSQLite::SQLite(), db_path) 
bisnode_data <- do.call(rbind.data.frame,lapply(CorrespondenceTable$id_1, function(y) {
  dbGetQuery(my_db_connect, paste0("SELECT * FROM bisnode WHERE personenid =='",y,"'"))}))
nr_data <- dbGetQuery(my_db_connect, paste0("SELECT * FROM nr_data")) %>%
  filter(ID==pol_id)
dbDisconnect(my_db_connect)

QueryPoliticianMandates(pol_id=pol_id,
                        correspondencetable=CorrespondenceTable,
                        data_input=bisnode_data,
                        source="bisnode",
                        year_range=c(1994:2017))

bisnode_data %>% select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum) %>%
  arrange(eintrittdatum) %>%
  distinct()

bisnode_data %>% select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum) %>%
  filter(duns==487025892) %>%
  arrange(eintrittdatum) %>%
  distinct()


# d) Check distribution of eintrittdatum, austrittdatum

my_db_connect <- dbConnect(RSQLite::SQLite(), db_path) 
bisnode_all <- dbGetQuery(my_db_connect, paste0("SELECT * FROM bisnode")) 
dbDisconnect(my_db_connect)

bisnode_all %>%
  janitor::tabyl(eintrittdatum) %>%
  arrange(eintrittdatum)

bisnode_all %>%
  janitor::tabyl(austrittdatum) %>%
  arrange(-percent)


#----------------------
# D) Sugarcube 
#----------------------

CorrespondenceTable <- GetCorrespondenceTable(modelname="Generation_1",
                                              source="sugarcube",
                                              linkscore_cutoff="optimal") %>%
  arrange(id_0)

test <- QueryPoliticianMandates(pol_id="AG-1931-0010",
                                correspondencetable=CorrespondenceTable,
                                year_range=c(1950:2015),
                                source="sugarcube")
View(test)  

CorrespondenceTable %>% filter(id_0=="AG-1931-0010") %>%
  arrange(year_1)


QueryPoliticianMandatesLoop(pol_ids=c("AG-1931-0010","AG-1931-0013"),
                            modelname="Generation_1",
                            year_range=c(1994:2017),
                            source="sugarcube",
                            linkscore_cutoff="optimal")


