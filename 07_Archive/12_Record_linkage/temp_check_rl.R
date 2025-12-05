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
library(readxl)

rm(list = ls())
rm(list=setdiff(ls(all.names=TRUE), lsf.str(all.names=TRUE)))

setwd("E:/Projekt Nationalräte")
source("./03_Code/12_Record_linkage/00_Functions_Record_linkage.R")
source("./03_Code/12_Record_linkage/00_Functions.R")
path <- "./02_Processed_data/12_Record_linkage/03_Data_out_analysis/"


#-------------------
# (B) Read in data
#-------------------

ids_post1994 <- read_dta("./02_Processed_data/nationalraete_1931_2015.dta") %>%
  filter(year>=1995) %>%
  distinct(ID) %>%
  pull()

length(ids_post1994)

df <- read_excel("./02_Processed_data/12_Record_linkage/04_Post_RL_QualChecks/bisnode_fp_firstcheck_Generation7_optimal_corrected.xlsx") %>%
  filter(is.na(CodierungFin)) %>%
  group_by(id_0,id_1) %>%
  mutate(age_firstmandate_max=max(age_firstmandate,na.rm=T),
         age_firstmandate_min=min(age_firstmandate,na.rm=T),
         age_lastmandate_max=max(age_lastmandate,na.rm=T))%>%
  select(id_0,id_1,starts_with("age")) %>%
  filter(age_firstmandate_max<=85 & age_lastmandate_max<=100 &
           age_firstmandate_min>=18 )%>% # see Email 19 November 2022 
  filter(id_0%in%ids_post1994)

df_correspondence_table <- df %>%
  distinct(id_0,id_1)


data_input <- data.table::fread("./02_Processed_data/11_Directors_1994_2018/Bisnode_Person-Firmen_Geo_Dataconstruct.csv",
                                sep = ";", encoding = "UTF-8")


data_wide <- BisnodeDataCut(data=data_input,id_bis=df_correspondence_table$id_1) %>%
  filter(eintrittjahr<austrittjahr) %>%
  mutate(spell_length=austrittjahr-eintrittjahr)

ggplot(data_wide,aes(spell_length)) +
  geom_histogram(breaks=c(0:25))

ggplot(data_wide,aes(austrittjahr)) +
  geom_histogram(breaks=c(1994:2017))

ggplot(data_wide,aes(eintrittjahr)) +
  geom_histogram(breaks=c(1994:2017))

ggplot(data_wide %>% filter(austrittjahr!=1995 & austrittjahr!=2016),aes(spell_length)) +
  geom_histogram(breaks=c(0:25))


# check for one specific person: AG-1975-9121, Ulrich Fischer

id_check <- df %>% filter(id_0=="AG-1975-9121") %>%
  distinct(id_1) %>%
  pull()

data_input %>% filter(personenid %in% id_check) %>%
  select(vorname,nachname,duns,personenid,eintrittdatum,austrittdatum,gremium,funktion,rechtsform)

ReShapeResults(pol_id="AG-1975-9121",
               nr_data=nr_data,
               input_data=df_all %>% filter(id_0=="AG-1975-9121"),
               max_period=(max(year_range_all)-min(year_range_all)+2))

