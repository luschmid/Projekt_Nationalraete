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

setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte")
source("C:/Schmidlu/Dropbox/Record Linkage/01_Code/03_Record_Linkage_Results/00_Functions.R")
#source("./03_Code/12_Record_linkage/00_Functions.R")
path <- "./02_Processed_data/12_Record_linkage/03_Data_out_analysis/"

#--------------------------------
# B) Bisnode Gen 7: optimal cutoff
#--------------------------------

# (i) Set parameters

modelname="Generation_7";
year_range=c(1994:2017);
source="bisnode";
linkscore_cutoff="optimal";

# (ii) Load correspondence table, nr data, and bisnode data


CorrespondenceTable <- GetCorrespondenceTable(modelname=modelname,
                                              source=source,
                                              linkscore_cutoff=linkscore_cutoff) %>%
  arrange(id_0) %>% 
  distinct(id_0, id_1)


id_1_correspondencetable <- CorrespondenceTable %>% select(id_1) %>%
  distinct() %>% 
  pull()

CorrespondenceTable %>% filter(id_1 %in% c(760,76920))


nr_data <- read_dta("./02_Processed_data/02_Elections_1971_2015/nationalraete_1931_2015_Geo.dta")
  
data_input <- fread("./02_Processed_data/11_Directors_1994_2018/Bisnode_Person_Geo.csv",
                                    sep = ";", encoding = "UTF-8") %>%  
      filter(personenid %in% as.numeric(id_1_correspondencetable))

data_input %>% filter(personenid %in% c(760,76920)) %>% select(vorname,nachname)
    
# (iii) Get minimal and maximal year per bisnode id to calculate implicit age

data_input_cut <- data_input %>% 
  select(personenid,anrede,vorname,nachname,geburtstag,eintrittdatum,austrittdatum,E_CNTR_w,N_CNTR_w) %>%
  #filter(is.na(eintrittdatum))
  mutate(eintrittdatum=as.character(eintrittdatum)) %>%
  mutate(austrittdatum=as.character(austrittdatum)) %>%
  mutate(eintrittdatum=ifelse(is.na(eintrittdatum), "1994-01-01",eintrittdatum)) %>%
  mutate(austrittdatum=ifelse(is.na(austrittdatum), "2017-12-31",austrittdatum)) %>%
  mutate(eintrittdatum=ymd(eintrittdatum)) %>%
  mutate(austrittdatum=ymd(austrittdatum)) %>%
  mutate(eintrittyear=as.numeric(year(eintrittdatum))) %>%
  mutate(austrittyear=as.numeric(year(austrittdatum))) %>%
  as_tibble() %>%
  dplyr::group_by(personenid) %>%
  dplyr::summarize(eintrittyear_min=min(eintrittyear,na.rm=T), 
                   austrittyear_max=max(austrittyear,na.rm=T))

data_input_cut %>% filter(personenid %in% c(760,76920)) 

    
# (iv) Prepare for merge
    
nr_data_to_merge <- nr_data %>%
  select(firstname,name,ID,birthyear,E_CNTR_w,N_CNTR_w,sex)  %>%
  mutate(sex_pol=ifelse(sex==1,"Male","Female")) %>%
  distinct() %>%
  dplyr::rename(lastname_pol=name,
                firstname_pol=firstname,
                birthyear_pol=birthyear,
                E_CNTR_w_pol=E_CNTR_w,
                N_CNTR_w_pol=N_CNTR_w,
                id_0=ID)

bisnode_to_merge <- data_input %>%   
  left_join(data_input_cut,by=c("personenid")) %>%    
  mutate(id_1=as.character(personenid),
         nchar_geb=nchar(geburtstag),
         geburtstag=ifelse(nchar_geb==4,paste0(geburtstag,"01-01"),geburtstag)) %>% 
  mutate(geburtstag=ymd(geburtstag)) %>%
  mutate(sex_bis=ifelse(anrede=="männlich","Male","Female"),
         birthyear_bis=year(geburtstag))%>%
  dplyr::rename(lastname_bis=nachname,
                firstname_bis=vorname,
                E_CNTR_w_bis=E_CNTR_w,
                N_CNTR_w_bis=N_CNTR_w,
                eintrittyear_min_bis= eintrittyear_min,
                austrittyear_max_bis=austrittyear_max) %>%
  select(id_1,ends_with("_bis"),gremium) 

data_input %>% filter(personenid=="813872")

bisnode_to_merge %>% filter(id_1 %in% c(760,76920)) 



# (v) Do merge and keep only those obs with rechtsform=="Aktiengesellschaft" & gremium=="Verwaltungsrat"

CorrespondenceTable %>% 
  distinct(id_0,id_1) %>% 
  left_join(nr_data_to_merge,by=c("id_0")) %>%
  filter(is.na(id_1)==FALSE) %>% 
  #left_join(bisnode_to_merge,by=c("id_1")) %>%
  filter(id_0=="AG-1935-0027") 

nr_data_to_merge %>% distinct(firstname_pol,lastname_pol,id_0,birthyear_pol,
                              E_CNTR_w_pol,N_CNTR_w_pol,sex_pol)


bisnode_to_merge %>% distinct(id_1)
data_input %>% filter(personenid==110)

# extract ids with multiple observations
nr_mult <- nr_data_to_merge %>% group_by(id_0) %>% mutate(nobs=n()) %>% filter(nobs>1 ) 
bisnode_mult <- bisnode_to_merge %>% group_by(id_1) %>% mutate(nobs=n()) %>% filter(nobs>1 ) 

data_input %>% filter(personenid==110) %>%
  select(duns,personenid,vorname,nachname,gremium,funktion)

CorrespondenceTable %>% 
  left_join(nr_mult %>% distinct(id_0,nobs),by=c("id_0")) %>%
  left_join(bisnode_mult %>% distinct(id_1,nobs),by=c("id_1")) %>%
  filter(!is.na(nobs.x) & !is.na(nobs.y))

nr_data_to_merge %>% filter(id_0=="AG-1935-0027") 
CorrespondenceTable %>% filter(id_0=="AG-1935-0027") 

ids_tolook <- CorrespondenceTable %>% filter(id_0=="AG-1935-0027") %>% select(id_1) %>% pull()

bisnode_to_merge %>% filter(id_1%in% ids_tolook)


data_out %>% filter(id_0=="AG-1935-0027") 



data_out <- nr_data_to_merge %>%
  left_join(CorrespondenceTable, 
            by=c("id_0")) %>%
  filter(is.na(id_1)==F) %>%
  left_join(bisnode_to_merge,by=c("id_1")) %>%
  mutate(distance=sqrt((E_CNTR_w_bis-E_CNTR_w_pol)^2+
                       (N_CNTR_w_bis-N_CNTR_w_pol)^2)/10000,
       age_diff=birthyear_pol-birthyear_bis,
       age_firstmandate=eintrittyear_min_bis-birthyear_pol,
       age_lastmandate=austrittyear_max_bis-birthyear_pol) %>%
  #filter(rechtsform=="Aktiengesellschaft" & gremium=="Verwaltungsrat" ) %>%
  select(id_0,id_1,firstname_pol,lastname_pol,firstname_bis, lastname_bis,
       sex_pol,sex_bis,
       distance,age_diff,age_firstmandate,age_lastmandate)

data_out %>% filter(id_1 %in% c(760,76920)) %>% tail()


data_fp_firstcheck <- data_out %>% 
  mutate(firstname_pol_str=tolower(StringReplace(firstname_pol)),
         firstname_bis_str=tolower(StringReplace(firstname_bis)),
         lastname_pol_str=tolower(StringReplace(lastname_pol)),
         lastname_bis_str=tolower(StringReplace(lastname_bis))) %>%
  filter(firstname_pol_str!=firstname_bis_str | 
                      lastname_pol_str!=lastname_bis_str) 

data_fp_firstcheck2 <- data_out %>% 
  mutate(firstname_pol_str=tolower(StringReplace(firstname_pol)),
         firstname_bis_str=tolower(StringReplace(firstname_bis)),
         lastname_pol_str=tolower(StringReplace(lastname_pol)),
         lastname_bis_str=tolower(StringReplace(lastname_bis))) %>%
  filter(firstname_pol_str==firstname_bis_str & 
           lastname_pol_str==lastname_bis_str)

data_fp_firstcheck2 %>% filter(id_1 %in% c(760,76920)) %>% tail()
data_fp_firstcheck2 %>% filter(id_1 %in% c(760,76920))

# (vii) Output for manual checks by Oester and Marilley (Oct/Nov 2022)

data.table::fwrite(
x = data_fp_firstcheck,
file = paste0("./02_Processed_data/12_Record_linkage/03_Data_analysis/",source,"_fp_firstcheck_",modelname,"_",
              gsub("\\.","\\_",linkscore_cutoff),".csv"),
sep = ";"
)

# (viii) Output for Bisnode bock 

# Note: For the output checks in (vii), we have only used matched observations
#       with differences either in the firstname or in the lastname. The 
#       following dataset includes all observations that have no differences
#       in the firstlame and in the lastname. The file 05_Setup_DataAnalysis.do
#       in the folder "03_Code\11_Directors_1994_2018" constructs this and
#       prepares the rl output for the analysis. 


data.table::fwrite(
  x = data_fp_firstcheck2,
  file = paste0("./02_Processed_data/12_Record_linkage/01_Bisnode/RL-Results-G7-perfectmatches.csv"),
  sep = ","
)


# 
# 
# # (ix) Read in manual checks by Oester and Marilley and consolidated by 
# #       Mark and Simon, Email Mark 19 November 2022
# 
# df_cor <- read_excel("./02_Processed_data/12_Record_linkage/04_Post_RL_QualChecks/bisnode_fp_firstcheck_Generation7_optimal_corrected.xlsx") %>%
#   filter(is.na(CodierungFin)) %>%
#   group_by(id_0,id_1) %>%
#   mutate(age_firstmandate_max=max(age_firstmandate,na.rm=T),
#          age_firstmandate_min=min(age_firstmandate,na.rm=T),
#          age_lastmandate_max=max(age_lastmandate,na.rm=T))%>%
#   select(id_0,id_1,starts_with("age")) %>%
#   filter(age_firstmandate_max<=85 & age_lastmandate_max<=100 &
#            age_firstmandate_min>=18  ) # see Email 19 November 2022
# 
# data.table::fwrite(
#   x = df_cor,
#   file = "./02_Processed_data/12_Record_linkage/01_Bisnode/RL-Results-generation_7_corrected.csv",
#   sep = ";")
# 
# # 
# 
# # (viii) check procedure 22.11.22
# 
# df_cor <- read_excel("./02_Processed_data/12_Record_linkage/04_Post_RL_QualChecks/bisnode_fp_firstcheck_Generation7_optimal_corrected.xlsx")
#                      
# # merge with outfile above
# data_fp_firstcheck %>%
#   anti_join(df_cor %>% mutate(id_1=as.character(id_1)) , by=c("id_0","id_1")) %>%
#   select(-firstname_pol_str,-firstname_bis_str,-lastname_pol_str,lastname_bis_str) %>%
#   print(n=50) 
# 
# # Result: Two connections different: 
# # id_0         id_1    firstname_pol lastname_pol firstname_bis lastname_bis sex_pol sex_bis distance age_diff age_firstmandate age_lastmand~1 lastn~2
# # TG-1991-0048 144716  Nicolo        Paganini     Nicol?        Paganini     Male    Male        1.11        0               32             51 pagani~
# # ZH-1987-0268 1134070 Stefan        Harangoz?    Stefan        Harangozo    Male    Male        1.26       NA               48             50 harang~
# #
# # Both connections seem to be correct connections
# 
# #df_cor %>% filter(id_0=="SG-2015-0131" & id_1== 144716)


#--------------------------------
# B) Bisnode Gen 7: cutoff 0.0
#--------------------------------

# (i) Set parameters

modelname="Generation_7";
year_range=c(1994:2017);
source="bisnode";
linkscore_cutoff=0;

# (ii) Load correspondence table, nr data, and bisnode data

CorrespondenceTable <- GetCorrespondenceTable(modelname=modelname,
                                              source=source,
                                              linkscore_cutoff=linkscore_cutoff) %>%
  arrange(id_0) %>% 
  distinct(id_0, id_1)

id_1_correspondencetable <- CorrespondenceTable %>% select(id_1) %>%
  distinct() %>% 
  pull()

nr_data <- read_dta("./02_Processed_data/02_Elections_1971_2015/nationalraete_1931_2015_Geo.dta")

data_input <- fread("./02_Processed_data/11_Directors_1994_2018/Bisnode_Person_Geo.csv",
                    sep = ";", encoding = "UTF-8") %>%  
  filter(personenid %in% as.numeric(id_1_correspondencetable))

# (iii) Get minimal and maximal year per bisnode id to calculate implicit age

data_input_cut <- data_input %>% 
  select(personenid,anrede,vorname,nachname,geburtstag,eintrittdatum,austrittdatum,E_CNTR_w,N_CNTR_w) %>%
  #filter(is.na(eintrittdatum))
  mutate(eintrittdatum=as.character(eintrittdatum)) %>%
  mutate(austrittdatum=as.character(austrittdatum)) %>%
  mutate(eintrittdatum=ifelse(is.na(eintrittdatum), "1994-01-01",eintrittdatum)) %>%
  mutate(austrittdatum=ifelse(is.na(austrittdatum), "2017-12-31",austrittdatum)) %>%
  mutate(eintrittdatum=ymd(eintrittdatum)) %>%
  mutate(austrittdatum=ymd(austrittdatum)) %>%
  mutate(eintrittyear=as.numeric(year(eintrittdatum))) %>%
  mutate(austrittyear=as.numeric(year(austrittdatum))) %>%
  as_tibble() %>%
  dplyr::group_by(personenid) %>%
  dplyr::summarize(eintrittyear_min=min(eintrittyear,na.rm=T), 
                   austrittyear_max=max(austrittyear,na.rm=T))

# (iv) Prepare for merge and then merge

nr_data_to_merge <- nr_data %>%
  select(firstname,name,ID,birthyear,E_CNTR_w,N_CNTR_w,sex)  %>%
  mutate(sex_pol=ifelse(sex==1,"Male","Female")) %>%
  distinct() %>%
  dplyr::rename(lastname_pol=name,
                firstname_pol=firstname,
                birthyear_pol=birthyear,
                E_CNTR_w_pol=E_CNTR_w,
                N_CNTR_w_pol=N_CNTR_w,
                id_0=ID)


bisnode_to_merge <- data_input %>%   
  left_join(data_input_cut,by=c("personenid")) %>%    
  mutate(id_1=as.character(personenid),
         nchar_geb=nchar(geburtstag),
         geburtstag=ifelse(nchar_geb==4,paste0(geburtstag,"01-01"),geburtstag)) %>% 
  mutate(geburtstag=ymd(geburtstag)) %>%
  mutate(sex_bis=ifelse(anrede=="männlich","Male","Female"),
         birthyear_bis=year(geburtstag))%>%
  dplyr::rename(lastname_bis=nachname,
                firstname_bis=vorname,
                E_CNTR_w_bis=E_CNTR_w,
                N_CNTR_w_bis=N_CNTR_w,
                eintrittyear_min_bis= eintrittyear_min,
                austrittyear_max_bis=austrittyear_max) %>%
  select(id_1,ends_with("_bis"),gremium) 


# (v) Do merge and keep only those obs with rechtsform=="Aktiengesellschaft" & gremium=="Verwaltungsrat"

data_out <- nr_data_to_merge %>%
  left_join(CorrespondenceTable, 
            by=c("id_0")) %>%
  filter(is.na(id_1)==F) %>%
  left_join(bisnode_to_merge,by=c("id_1")) %>%
  mutate(distance=sqrt((E_CNTR_w_bis-E_CNTR_w_pol)^2+
                         (N_CNTR_w_bis-N_CNTR_w_pol)^2)/10000,
         age_diff=birthyear_pol-birthyear_bis,
         age_firstmandate=eintrittyear_min_bis-birthyear_pol,
         age_lastmandate=austrittyear_max_bis-birthyear_pol) %>%
  filter(rechtsform=="Aktiengesellschaft" & gremium=="Verwaltungsrat" ) %>%
  select(id_0,id_1,firstname_pol,lastname_pol,firstname_bis, lastname_bis,
         sex_pol,sex_bis,
         distance,age_diff,age_firstmandate,age_lastmandate)

data_fp_firstcheck <- data_out %>% 
  mutate(firstname_pol_str=tolower(StringReplace(firstname_pol)),
         firstname_bis_str=tolower(StringReplace(firstname_bis)),
         lastname_pol_str=tolower(StringReplace(lastname_pol)),
         lastname_bis_str=tolower(StringReplace(lastname_bis))) %>%
  filter(firstname_pol_str!=firstname_bis_str | 
           lastname_pol_str!=lastname_bis_str) 


# (vi) remove already checked connections from optimal cutoff (see above)

fp_optimal <- data.table::fread(
  file = paste0("./02_Processed_data/12_Record_linkage/03_Data_analysis/",source,"_fp_firstcheck_",modelname,"_","optimal",".csv"),
  sep = ";") %>%
  select(id_0,id_1) %>%
  mutate(id_0=as.character(id_0),
         id_1=as.character(id_1)) %>%
  as_tibble() %>%
  distinct(id_0,id_1) %>%
  mutate(optimal=1)


# (vii) Focus on those below 85 for last start of a mandate and below 100 for the
#      last mandate end and those above 18 for the first start of a mandate

data_out <- data_fp_firstcheck %>%
  left_join(fp_optimal,by=c("id_0","id_1")) %>%
  filter(is.na(optimal)) %>%group_by(id_0,id_1) %>%
  mutate(age_firstmandate_max=max(age_firstmandate,na.rm=T),
         age_firstmandate_min=min(age_firstmandate,na.rm=T),
         age_lastmandate_max=max(age_lastmandate,na.rm=T))%>%
  filter(age_firstmandate_max<=85 & age_lastmandate_max<=100 &
           age_firstmandate_min>=18  )

data.table::fwrite(
  x = data_out,
  file = paste0("./02_Processed_data/12_Record_linkage/03_Data_analysis/",source,"_fp_firstcheck_",modelname,"_",
                gsub("\\.","\\_",linkscore_cutoff),".csv"),
  sep = ";"
)
# 
# #-----------------------------------------
# # C) Create outfile for data construction
# #-----------------------------------------
# 
# # (i) Read-in first batch and second batch
# 
# library(xlsx)
# 
# df_cor_firstbatch <- read_excel("./02_Processed_data/12_Record_linkage/01_Bisnode/04_Post_RL_QualChecks/bisnode_fp_firstcheck_Generation7_optimal_corrected.xlsx") %>%
#   filter(is.na(CodierungFin)) %>%
#   group_by(id_0,id_1) %>%
#   mutate(age_firstmandate_max=max(age_firstmandate,na.rm=T),
#          age_firstmandate_min=min(age_firstmandate,na.rm=T),
#          age_lastmandate_max=max(age_lastmandate,na.rm=T))%>%
#   select(id_0,id_1,starts_with("age")) %>%
#   filter(age_firstmandate_max<=85 & age_lastmandate_max<=100 &
#            age_firstmandate_min>=18  ) # see Email 19 November 2022
# 
# df_cor_secondbatch <- read_excel("./02_Processed_data/12_Record_linkage/01_Bisnode/04_Post_RL_QualChecks/bisnode_fp_firstcheck_Generation7_0_corrected.xlsx") %>%
#   filter(is.na(CodierungFin)) %>%
#   group_by(id_0,id_1) %>%
#   mutate(age_firstmandate_max=max(age_firstmandate,na.rm=T),
#          age_firstmandate_min=min(age_firstmandate,na.rm=T),
#          age_lastmandate_max=max(age_lastmandate,na.rm=T))%>%
#   select(id_0,id_1,starts_with("age")) %>%
#   filter(age_firstmandate_max<=85 & age_lastmandate_max<=100 &
#            age_firstmandate_min>=18) # see Email 19 November 2022
# 
# df_cor_firstbatch %>% distinct(id_0,id_1)
# df_cor_secondbatch %>% distinct(id_0,id_1)
# 
# # (ii) Combine first and second batch and create output file
# 
# df_out <- bind_rows(df_cor_firstbatch,df_cor_secondbatch)%>% 
#   distinct(id_0,id_1)
# 
# data.table::fwrite(
#   x = df_out,
#   file = paste0("./02_Processed_data/12_Record_linkage/01_Bisnode/04_Post_RL_QualChecks/bisnode_fp_firstcheck_Generation_7_0_final.csv"),
#   sep = ";"
# )

