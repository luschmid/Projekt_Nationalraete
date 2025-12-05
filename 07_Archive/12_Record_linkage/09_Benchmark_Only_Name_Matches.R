
#!/usr/bin/Rscript

# 1. Set path and load libraries and functions

path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte"
#path <- "E:/Projekt Nationalräte"

setwd(path)

library(tidyverse)
require(readxl)
library(data.table)

source("./03_Code/12_Record_linkage/00_Functions.R")


# 2. Read in nr data, bisnode data and test and validation set data

# (a) NR data incl. transformation of Bürgerorte


nr_data <- readstata13::read.dta13("./02_Processed_data/02_Elections_1971_2015/nationalraete_1931_2015_Geo.dta") %>%
  rename(id_0=ID)

nr_data_single_buergerort <- nr_data %>% 
  filter(is.na(E_CNTR_b2 )==T)  %>%
  rename(E_CNTR_b=E_CNTR_b1,N_CNTR_b=N_CNTR_b1 ) %>%
  select(id_0, canton, sex, name, firstname , birthyear , gdename, job, list, 
         E_CNTR_w,N_CNTR_w,E_CNTR_b,N_CNTR_b) %>%
  distinct(id_0,name, firstname,E_CNTR_w,N_CNTR_w,.keep_all = T) %>%
  arrange(id_0,name)

nr_data_multiple_buergerort <- nr_data %>% 
  filter(is.na(E_CNTR_b2 )==F)

nr_data_multiple_buergerort_long <- nr_data_multiple_buergerort %>%
  pivot_longer(starts_with(c("E_CNTR_b","N_CNTR_b")), names_to = "var", values_to = "test") %>%
  filter(is.na(test)==F) %>%
  distinct(id_0,name, firstname,E_CNTR_w,N_CNTR_w,var,.keep_all = T) %>%
  arrange(id_0,name)

nr_e_cntr <- nr_data_multiple_buergerort_long %>% select(id_0,var,test) %>% 
  filter(var %in% c("E_CNTR_b1","E_CNTR_b2","E_CNTR_b3","E_CNTR_b4","E_CNTR_b5","E_CNTR_b6")) %>%
  rename(E_CNTR_b=test) %>% 
  mutate(mergevar=substr(var,8,9))

nr_n_cntr <- nr_data_multiple_buergerort_long %>% select(id_0,var,test) %>% 
  filter(var %in% c("N_CNTR_b1","N_CNTR_b2","N_CNTR_b3","N_CNTR_b4","N_CNTR_b5","N_CNTR_b6")) %>%
  rename(N_CNTR_b=test) %>% 
  mutate(mergevar=substr(var,8,9))

nr_cntr_all <- nr_e_cntr %>% full_join(nr_n_cntr, by=c("id_0","mergevar")) %>%
  select(id_0,E_CNTR_b,N_CNTR_b)

nr_cntr_all %>% filter(is.na(E_CNTR_b) | is.na(N_CNTR_b))

df_multiple_buergerorte <- CreateDFAllBuergerorte(data=nr_data_multiple_buergerort,
                                                  data_buergerorte=nr_cntr_all) %>%
  select(id_0, canton, sex, name, firstname , birthyear , gdename, job, list, 
         E_CNTR_w,N_CNTR_w,E_CNTR_b,N_CNTR_b) 


# (b) Bisnode data

bisnode_data <-data.table::fread(file="./02_Processed_data/11_Directors_1994_2018/bisnode_export_record_linkage_oct2021.csv", 
                                 sep=",",encoding = "UTF-8")

# (c) Test and validation set

test_set <- ReadinGroundTruthFiles(file_name = "test_set_split_size_0_35_random_state_55.csv",
                                   separator=",",
                                   path_rel="./02_Processed_data/12_Record_linkage/01_Bisnode/")
validation_set <- ReadinGroundTruthFiles(file_name = "validation_set_split_size_0_35_random_state_55.csv",
                                         separator=",",
                                         path_rel="./02_Processed_data/12_Record_linkage/01_Bisnode/")

validation_set_unique <- validation_set %>% distinct(id_0)

ground_truth_sep <-bind_rows(test_set %>% mutate(tset=1),
                             validation_set %>% mutate(vset=1))%>%
  mutate(vset=ifelse(is.na(vset),0,1),
         tset=ifelse(is.na(tset),0,1))

nr_data_val_set <- nr_data %>% 
  right_join(validation_set,by=c("id_0")) 

# (d) Define dataset to search over
# Note: These are all observations in the NR data with all possible place of origins
#       and all possible residence municipalities who are in the validation set

df_tosearch <- bind_rows(nr_data_single_buergerort,
                         df_multiple_buergerorte) %>%
  distinct(id_0,name, firstname,E_CNTR_w,N_CNTR_w,E_CNTR_b,N_CNTR_b,.keep_all = T) %>%
  right_join(validation_set_unique,by=c("id_0"))


# 3. Implement Search
# Note: This was done on Unifr server on Oct 29, 2021. 

# run procedure (from command line with department name as input)

args = commandArgs(trailingOnly=TRUE)
#KotakorpiEtAl2007Loop(canton_input=args[1],data=data_ch,control_seat_allocation=0,check_convergence=1) 

#print(paste0("Arguments: ",args[1]," and", args[2]))

RL_StdSearch(data_target=bisnode_data,
             data_input=df_tosearch,
             max_distances=list(as.numeric(args[1]),as.numeric(args[2])),
             filename=paste0("RL_Std_2021_10_29","_",args[1],"_",args[2]),
             progress=T,
             append_choice=F,
             path="./02_Processed_data/12_Record_linkage/01_Bisnode/")


# 4. Read in results and calculate ROC Curve

path_files <- "//Srw-hpc5/e/Projekt Nationalräte/02_Processed_data/12_Record_linkage/01_Bisnode/"
filenames <- list.files(path=path_files,pattern="RL_Std_2021*")
filenames_fullpath=file.path(path_files,filenames)

data_std <- do.call("rbind",lapply(filenames_fullpath,FUN=function(files){ fread(files,encoding="UTF-8")}))

data_std <- data_std %>% rename(id_0=ID_Pol,
                                id_1=ID_VR,
                                firstname_0=Vorname_Pol,
                                name_0=Nachname_Pol,
                                firstname_1=Vorname_VR,
                                name_1=Nachname_VR,
                                w_distance=W_Distanz)

data_std_0_0 <- data_std[data_std$max_distance_input_name==0 &
                           data_std$max_distance_input_firstname==0 & 
                           name_0==name_1 & firstname_0==firstname_1 ,]

CalculateROCCurve(recordlinkagedata=data_std_0_0 %>%
                    select(id_0,id_1,w_distance,max_distance_input_name,max_distance_input_firstname),
                  groundtruthdata = ground_truth_sep,
                  idvars=c("id_0","id_1"),
                  tvvars=c("tset","vset"),
                  step_size = 0.02,
                  path_rel = "04_Results/02_Record_linkage",
                  file_name = "ROC_Curve_std_results_0_0" ,
                  source="std")

data_std_0.05_0.05 <- data_std[data_std$max_distance_input_name==0.05 &
                           data_std$max_distance_input_firstname==0.05 ,]

CalculateROCCurve(recordlinkagedata=data_std_0.05_0.05 %>%
                    select(id_0,id_1,w_distance,max_distance_input_name,max_distance_input_firstname),
                  groundtruthdata = ground_truth_sep,
                  idvars=c("id_0","id_1"),
                  tvvars=c("tset","vset"),
                  step_size = 0.02,
                  path_rel = "04_Results/02_Record_linkage",
                  file_name = "ROC_Curve_std_results_0.05_0.05" ,
                  source="std")
                  
                  
data_std_0.1_0.1 <- data_std[data_std$max_distance_input_name==0.1 &
                                 data_std$max_distance_input_firstname==0.1 ,]

CalculateROCCurve(recordlinkagedata=data_std_0.1_0.1  %>%
                    select(id_0,id_1,w_distance,max_distance_input_name,max_distance_input_firstname),
                  groundtruthdata = ground_truth_sep,
                  idvars=c("id_0","id_1"),
                  tvvars=c("tset","vset"),
                  step_size = 0.02,
                  path_rel = "04_Results/02_Record_linkage",
                  file_name = "ROC_Curve_std_results_0.1_0.1" ,
                  source="std")                  


data_std_0.05_0.1 <- data_std[data_std$max_distance_input_name==0.05 &
                               data_std$max_distance_input_firstname==0.1 ,]

CalculateROCCurve(recordlinkagedata=data_std_0.05_0.1  %>%
                    select(id_0,id_1,w_distance,max_distance_input_name,max_distance_input_firstname),
                  groundtruthdata = ground_truth_sep,
                  idvars=c("id_0","id_1"),
                  tvvars=c("tset","vset"),
                  step_size = 0.02,
                  path_rel = "04_Results/02_Record_linkage",
                  file_name = "ROC_Curve_std_results_0.05_0.1" ,
                  source="std")  
