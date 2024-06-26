# Functions

GenerateCantonSpecificVariables <- function(data,unitvariable="ctn",timevariable="year"){
  
  data$ctn <- data[[unitvariable]]
  data$year <- data[[timevariable]]
  
  data <- fastDummies::dummy_cols(sugarcube_coll_by_year,
                                  select_columns = "ctn")
  
  ct_levels <- levels(as.factor(sugarcube_coll_by_year$ctn))
  
  for (i in 1:length(ct_levels)){
    data[[paste0("ctn_",ct_levels[i])]] <-data[[paste0("ctn_",ct_levels[i])]]*data$year
    data[[paste0("ctn_",ct_levels[i],"_sq")]] <-data[[paste0("ctn_",ct_levels[i])]]*data$year^2
  }
  
  ct_elections_coll_by_year <- fastDummies::dummy_cols(ct_elections_coll_by_year,
                                                       select_columns = "canton")
  
  ct_levels <- levels(as.factor(ct_elections_coll_by_year$canton))
  
  for (i in 1:length(ct_levels)){
    ct_elections_coll_by_year[[paste0("canton_",ct_levels[i])]] <- ct_elections_coll_by_year[[paste0("canton_",ct_levels[i])]]*ct_elections_coll_by_year$year
    ct_elections_coll_by_year[[paste0("canton_",ct_levels[i],"_sq")]] <- ct_elections_coll_by_year[[paste0("canton_",ct_levels[i])]]*ct_elections_coll_by_year$year^2
  }
  return(data)  
}


#-----------------------------------------
# (A) Libraries, paths, and load functions
#-----------------------------------------


library(tidyverse)
library(data.table)
library(testit)
library(dbplyr)
library(lubridate)
library(haven)
library(fixest)
library(modelsummary)
library(ggrepel)
library(haven)
library(pxR)
library(janitor)


rm(list=setdiff(ls(all.names=TRUE), lsf.str(all.names=TRUE)))
setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte")

out_path="C:/Schmidlu/Dropbox/Apps/Overleaf/Female_Suffrage_Directors"
#-----------------
# (B) Read in data
#-----------------


# (i) Sugarcube 

# (a) wide format (one observations is a person per year)

sugarcube <- data.table::fread("./02_Processed_data/10_Directors_1934_2003/Sugarcube_Person_CleanName-Gender-Geo.csv",
                               sep = ";", encoding = "UTF-8")

sugarcube_coll_by_year <- sugarcube %>%
  dplyr::group_by(year,ctn) %>%
  dplyr::summarize(mean_female_vr=1-mean(male,na.rm=T),
                   nobs=n()) %>% 
  mutate(female_suffrage_ct=0) %>%
  mutate(female_suffrage_ct=case_when(
    ctn %in% c("ag","fr","sh","zg") & year>=1971 ~ 1 ,
    ctn== "ar" & year>=1990 ~ 1 ,
    ctn== "ai" & year>=1991 ~ 1 ,
    ctn== "bl" & year>=1969 ~ 1 ,
    ctn== "bs" & year>=1967 ~ 1 ,
    ctn== "be" & year>=1972 ~ 1 ,
    ctn== "ge" & year>=1961 ~ 1 ,
    ctn== "gl" & year>=1971 ~ 1 ,
    ctn== "gr" & year>=1972 ~ 1 ,
    ctn== "ju" & year>=1977 ~ 1 ,
    ctn== "lu" & year>=1972 ~ 1 ,
    ctn== "ne" & year>=1960 ~ 1 ,
    ctn== "nw" & year>=1972 ~ 1 ,
    ctn== "ow" & year>=1972 ~ 1 ,
    ctn== "sz" & year>=1972 ~ 1 ,
    ctn== "so" & year>=1972 ~ 1 ,
    ctn== "sg" & year>=1972 ~ 1 ,
    ctn== "ti" & year>=1970 ~ 1 ,
    ctn== "tg" & year>=1972 ~ 1 ,
    ctn== "ur" & year>=1972 ~ 1 ,
    ctn== "vd" & year>=1959 ~ 1 ,
    ctn== "vs" & year>=1971 ~ 1 ,
    ctn== "zh" & year>=1971 ~ 1, 
    TRUE~0))%>%
  mutate(early_intro=case_when(
    ctn %in% c("ne","vd","bs","bl","ge")  ~ 1 ,
    TRUE~0),
    dummy_1971=ifelse(year>=1971,1,0)) %>%
  mutate(
    suffc_date=case_when(
      ctn %in% c("ag", "fr", "sh", "zg")~mdy("2-7-1971"),
      ctn=="ar"~mdy("4-30-1989"),
      ctn=="ai"~mdy("11-27-1990"),
      ctn=="bl"~mdy("6-23-1968"),
      ctn=="bs"~mdy("6-26-1966"),
      ctn=="be"~mdy("12-12-1971"),
      ctn=="ge"~mdy("3-6-1960"),
      ctn=="gl"~mdy("5-2-1971"),
      ctn=="gr"~mdy("3-5-1972"),
      ctn=="ju"~mdy("3-20-1977"),
      ctn=="lu"~mdy("10-25-1971"),
      ctn=="ne"~mdy("9-27-1959"),
      ctn=="nw"~mdy("4-30-1972"),
      ctn=="ow"~mdy("9-24-1972"),
      ctn=="sz"~mdy("3-5-1972"),
      ctn=="so"~mdy("6-6-1971"),
      ctn=="sg"~mdy("1-23-1972"),
      ctn=="ti"~mdy("10-19-1969"),
      ctn=="tg"~mdy("12-12-1971"),
      ctn=="ur"~mdy("1-30-1972"),
      ctn=="vd"~mdy("2-1-1959"),
      ctn=="vs"~mdy("4-12-1970"),
      ctn=="zh"~mdy("11-15-1970")),
    suff_year=ifelse(month(suffc_date)<=6,year(suffc_date),year(suffc_date)+1)) %>%
  mutate(year_date=paste0("6-1-",year)) %>%
  mutate(time_since_treatment=year(mdy(year_date))-year(suffc_date),
         time_to_treatment=-time_since_treatment)%>% 
  filter(ctn!="")

sugarcube_coll_by_year %>% select(year, ctn, time_since_treatment,time_to_treatment)

# (b) long format (one observations is a person per year)

sugarcube_long <- data.table::fread("./02_Processed_data/10_Directors_1934_2003/sugarcube_long.csv",
                                    sep = ";", encoding = "UTF-8")

sugarcube_long_by_year <- sugarcube_long %>%
  dplyr::group_by(year,ctn) %>%
  dplyr::summarize(mean_female_vr=1-mean(male,na.rm=T),
                   nobs=n()) %>% 
  mutate(female_suffrage_ct=0) %>%
  mutate(female_suffrage_ct=case_when(
    ctn %in% c("ag","fr","sh","zg") & year>=1971 ~ 1 ,
    ctn== "ar" & year>=1990 ~ 1 ,
    ctn== "ai" & year>=1991 ~ 1 ,
    ctn== "bl" & year>=1969 ~ 1 ,
    ctn== "bs" & year>=1967 ~ 1 ,
    ctn== "be" & year>=1972 ~ 1 ,
    ctn== "ge" & year>=1961 ~ 1 ,
    ctn== "gl" & year>=1971 ~ 1 ,
    ctn== "gr" & year>=1972 ~ 1 ,
    ctn== "ju" & year>=1977 ~ 1 ,
    ctn== "lu" & year>=1972 ~ 1 ,
    ctn== "ne" & year>=1960 ~ 1 ,
    ctn== "nw" & year>=1972 ~ 1 ,
    ctn== "ow" & year>=1972 ~ 1 ,
    ctn== "sz" & year>=1972 ~ 1 ,
    ctn== "so" & year>=1972 ~ 1 ,
    ctn== "sg" & year>=1972 ~ 1 ,
    ctn== "ti" & year>=1970 ~ 1 ,
    ctn== "tg" & year>=1972 ~ 1 ,
    ctn== "ur" & year>=1972 ~ 1 ,
    ctn== "vd" & year>=1959 ~ 1 ,
    ctn== "vs" & year>=1971 ~ 1 ,
    ctn== "zh" & year>=1971 ~ 1, 
    TRUE~0))%>%
  mutate(early_intro=case_when(
    ctn %in% c("ne","vd","bs","bl","ge")  ~ 1 ,
    TRUE~0),
    dummy_1971=ifelse(year>=1971,1,0)) %>%
  mutate(
    suffc_date=case_when(
      ctn %in% c("ag", "fr", "sh", "zg")~mdy("2-7-1971"),
      ctn=="ar"~mdy("4-30-1989"),
      ctn=="ai"~mdy("11-27-1990"),
      ctn=="bl"~mdy("6-23-1968"),
      ctn=="bs"~mdy("6-26-1966"),
      ctn=="be"~mdy("12-12-1971"),
      ctn=="ge"~mdy("3-6-1960"),
      ctn=="gl"~mdy("5-2-1971"),
      ctn=="gr"~mdy("3-5-1972"),
      ctn=="ju"~mdy("3-20-1977"),
      ctn=="lu"~mdy("10-25-1971"),
      ctn=="ne"~mdy("9-27-1959"),
      ctn=="nw"~mdy("4-30-1972"),
      ctn=="ow"~mdy("9-24-1972"),
      ctn=="sz"~mdy("3-5-1972"),
      ctn=="so"~mdy("6-6-1971"),
      ctn=="sg"~mdy("1-23-1972"),
      ctn=="ti"~mdy("10-19-1969"),
      ctn=="tg"~mdy("12-12-1971"),
      ctn=="ur"~mdy("1-30-1972"),
      ctn=="vd"~mdy("2-1-1959"),
      ctn=="vs"~mdy("4-12-1970"),
      ctn=="zh"~mdy("11-15-1970")),
    suff_year=ifelse(month(suffc_date)<=6,year(suffc_date),year(suffc_date)+1)) %>%
  mutate(year_date=paste0("6-1-",year)) %>%
  mutate(time_since_treatment=year(mdy(year_date))-year(suffc_date),
         time_to_treatment=-time_since_treatment)%>% 
  filter(ctn!="")

sugarcube_coll_by_year %>% select(year, ctn, time_since_treatment,time_to_treatment)

# (ii) Cantonal elections data

ct_elections <- read_dta("C:/Schmidlu/Dropbox/DataCompleteAll_1950_2017.dta")

ct_elections_coll_by_year <- ct_elections %>%
  mutate(canton=tolower(canton))%>%
  dplyr::group_by(year,canton) %>%
  dplyr::summarize(mean_female_nr=mean(sex,na.rm=T),
                   nobs=n(),
                   LastLegislation=first(LastLegislation)) %>%
  mutate(female_suffrage_ct=case_when(
    canton %in% c("ag","fr","sh","zg") & year>=1971 ~ 1 ,
    canton== "ar" & year>=1990 ~ 1 ,
    canton== "ai" & year>=1991 ~ 1 ,
    canton== "bl" & year>=1969 ~ 1 ,
    canton== "bs" & year>=1967 ~ 1 ,
    canton== "be" & year>=1972 ~ 1 ,
    canton== "ge" & year>=1961 ~ 1 ,
    canton== "gl" & year>=1971 ~ 1 ,
    canton== "gr" & year>=1972 ~ 1 ,
    canton== "ju" & year>=1977 ~ 1 ,
    canton== "lu" & year>=1972 ~ 1 ,
    canton== "ne" & year>=1960 ~ 1 ,
    canton== "nw" & year>=1972 ~ 1 ,
    canton== "ow" & year>=1972 ~ 1 ,
    canton== "sz" & year>=1972 ~ 1 ,
    canton== "so" & year>=1972 ~ 1 ,
    canton== "sg" & year>=1972 ~ 1 ,
    canton== "ti" & year>=1970 ~ 1 ,
    canton== "tg" & year>=1972 ~ 1 ,
    canton== "ur" & year>=1972 ~ 1 ,
    canton== "vd" & year>=1959 ~ 1 ,
    canton== "vs" & year>=1971 ~ 1 ,
    canton== "zh" & year>=1971 ~ 1, 
    TRUE~0))%>%
  mutate(
    suffc_date=case_when(
      canton %in% c("ag", "fr", "sh", "zg")~mdy("2-7-1971"),
      canton=="ar"~mdy("4-30-1989"),
      canton=="ai"~mdy("11-27-1990"),
      canton=="bl"~mdy("6-23-1968"),
      canton=="bs"~mdy("6-26-1966"),
      canton=="be"~mdy("12-12-1971"),
      canton=="ge"~mdy("3-6-1960"),
      canton=="gl"~mdy("5-2-1971"),
      canton=="gr"~mdy("3-5-1972"),
      canton=="ju"~mdy("3-20-1977"),
      canton=="lu"~mdy("10-25-1971"),
      canton=="ne"~mdy("9-27-1959"),
      canton=="nw"~mdy("4-30-1972"),
      canton=="ow"~mdy("9-24-1972"),
      canton=="sz"~mdy("3-5-1972"),
      canton=="so"~mdy("6-6-1971"),
      canton=="sg"~mdy("1-23-1972"),
      canton=="ti"~mdy("10-19-1969"),
      canton=="tg"~mdy("12-12-1971"),
      canton=="ur"~mdy("1-30-1972"),
      canton=="vd"~mdy("2-1-1959"),
      canton=="vs"~mdy("4-12-1970"),
      canton=="zh"~mdy("11-15-1970"))) %>%
  mutate(year_date=paste0("6-1-",year)) %>%
  mutate(time_since_treatment=year(mdy(year_date))-year(suffc_date),
         time_to_treatment=-time_since_treatment)

# (iii) Covariates from regulation project by Simon and Mark (email Mark 10.08.2022)

covariates <- read_dta("./01_Raw_data/16_Female_Suffrage_Directors/data_regulation.dta") %>%
  mutate(ctn=tolower(canton_abr))

# (iv) Read in vote data

female_votes_ct <- data.table::fread("./01_Raw_data/16_Female_Suffrage_Directors/female_suffrage_vote_1971_ct.csv",
                                     sep = ";", encoding = "UTF-8") %>%
  mutate(ctn=tolower(ctn)) %>%
  dplyr::rename(yes_percent_female_1971=yes_percent) %>%
  select(ctn,yes_percent_female_1971)


# (v) Generate canton dummy and canton-specific time trends var

sugarcube_coll_by_year <- GenerateCantonSpecificVariables(sugarcube_coll_by_year)
sugarcube_long_by_year <- GenerateCantonSpecificVariables(sugarcube_long_by_year)


# (vi) Merge data

sugarcube_all <- sugarcube_coll_by_year %>%
  left_join(female_votes_ct,by=c("ctn")) %>%
  left_join(covariates,by=c("ctn","year")) %>%
  mutate(yes_percent_female_1971_q=ntile(yes_percent_female_1971, 5),
         ctn_year=paste0(ctn,"_",year))

sugarcube_all %>% filter(is.na(yes_percent_female_1971)) %>% tabyl(ctn)

sugarcube_all %>% select(year,ctn,starts_with("yes_percent_female_1971")) %>%
  filter(year>2000)

# Result: Only JU has missing values for vote on female suffrage in 1971


sugarcube_all_long <- sugarcube_long_by_year %>%
  left_join(female_votes_ct,by=c("ctn")) %>%
  left_join(covariates,by=c("ctn","year")) %>%
  mutate(yes_percent_female_1971_q=ntile(yes_percent_female_1971, 5),
         ctn_year=paste0(ctn,"_",year))

save.image(file = "./02_Processed_data/16_Female_Suffrage_Directors/Female_Suffrage_Directors.RData")


