
infile_csv_name <- "RL_Results_2021105.csv"
outfile_pdf_name = "ROC_Curve_results3_newmetric"

#-------------------------------------------------
# (A) Libraries, working directory and load data
#-------------------------------------------------

# (i) Libraries, working directory, and load functions

library(readxl)
library(tidyverse)
library(data.table)
library(haven)

path <- "C:/Schmidlu/Dropbox/Record Linkage"
setwd(path)

source("./01_Code/03_Record_Linkage_Results/00_Functions.R")


# (ii) Read in record linkage result files and reshape

rl_results1 <- ReadinResults(file_name = infile_csv_name,
                             path_rel="./02_Data/07_Output_RL_Bisnode/")


rl_results_wide1_all <- CreateWideDataset(data=rl_results1,
                                          generate_differences=T) 

rl_results_wide1_all %>% filter(!is.na(id_1)) 
rl_results_wide1_all %>% filter(!is.na(id_1)) %>%
  distinct(id_1)

rl_results_wide1_all %>% 
  filter(!is.na(id_1)) %>% 
  group_by(id_1) %>%
  mutate(nobs=n()) %>%
  ungroup() %>%
  filter(nobs>1) %>%
  arrange(id_1) %>%
  select(id_1,id_0,Cluster_ID,Link_Score, name_1,   name_0,  firstname_1, firstname_0)


# 2,122,759
# 2,103,645

# (iii) Read in ground truth files

validation_set <- ReadinGroundTruthFiles(file_name = "validation_set_randomized_buckets_0.35.csv",
                                   separator=",",
                                   path_rel="./02_Processed_data/12_Record_linkage/01_Bisnode/")


test_set <- ReadinGroundTruthFiles(file_name = "test_set_randomized_buckets_0.65.csv",
                                         separator=",",
                                         path_rel="./02_Processed_data/12_Record_linkage/01_Bisnode/")


ground_truth <- ReadinGroundTruthFiles(file_name = "bisnode_falsenegatives_final_july2021.csv",
                                   separator=";",
                                   path_rel="./02_Processed_data/12_Record_linkage/01_Bisnode/")
# commented out 7 February 2022
#
# test_set <- ReadinGroundTruthFiles(file_name = "test_set_split_size_0_35_random_state_55.csv",
#                                        separator=",",
#                                    path_rel="./02_Processed_data/12_Record_linkage/01_Bisnode/")
# validation_set <- ReadinGroundTruthFiles(file_name = "validation_set_split_size_0_35_random_state_55.csv",
#                                    separator=",",
#                                    path_rel="./02_Processed_data/12_Record_linkage/01_Bisnode/")

ground_truth_sep <-bind_rows(test_set %>% mutate(tset=1),
                             validation_set %>% mutate(vset=1))%>%
  mutate(vset=ifelse(is.na(vset),0,1),
         tset=ifelse(is.na(tset),0,1))

ground_truth_sep %>% dplyr::group_by(id_0) %>%
  dplyr::summarize(tset_max=max(tset,ra.rm=T),
            vset_max=max(vset,ra.rm=T)) 

# ground_truth %>% filter(id_0=="BEJU-2015-0182")
# ground_truth_sep %>% filter(id_0=="BEJU-2015-0182")
# ground_truth %>% filter(id_0=="BEJU-2011-0280")


# (iv) Check whether all obs are either in test or val set 

  
LookupDoubleEntries(data1=test_set,data2=validation_set)

validation_set %>% full_join(test_set,by=c("id_0"))  %>% print(n=50)
  
# Result: All fine, No. double entries: 0 

# (v) Check whether all obs in test/val set are in ground truth

(check <- ground_truth %>% left_join(ground_truth_sep,by=c("id_0","id_1")) %>%
  filter(is.na(tset) | is.na(vset) ))

# Result: 9 obs are only in ground truth

ground_truth %>% mutate(gt=1) %>% 
  full_join(ground_truth_sep %>% mutate(gts=1),by=c("id_0","id_1")) %>%
  filter(is.na(gt) | is.na(gts) )

# (vi) Read in Bisnode and NR date 

bisnode_data <- data.table::fread("./02_Processed_data/11_Directors_1994_2018/bisnode_export_record_linkage_oct2021.csv",
                             sep = ",", encoding = "UTF-8")

bisnode_data_nobs <- bisnode_data %>%
  group_by(name_bis, firstname_bis, e_cntr_w_bis, n_cntr_w_bis) %>%
  mutate(nobs=n()) %>%
  ungroup() %>%
  group_by(id_bis) %>%
  summarize(nobs_max_bis=max(nobs,na.rm=T)) 

janitor::tabyl(bisnode_data_nobs$nobs_max)

nr_data <- data.table::fread("C:/Schmidlu/Dropbox/Record Linkage/Data/Politikerdaten/nationalraete_export_record_linkage_july2021.csv", 
                             sep = ",", encoding = "UTF-8")

nr_data_nobs <- nr_data %>%
  group_by(name_polit, firstname_polit, e_cntr_w_polit, n_cntr_w_polit) %>%
  mutate(nobs=n()) %>%
  ungroup() %>%
  group_by(id_polit) %>%
  summarize(nobs_max_polit=max(nobs,na.rm=T)) 

janitor::tabyl(nr_data_nobs$nobs_max)

ground_truth_sep <- ground_truth_sep %>%
  left_join(bisnode_data_nobs %>% 
              dplyr::rename(id_1=id_bis) %>%
              mutate(id_1=as.character(id_1)),by=c("id_1"))%>%
  left_join(nr_data_nobs %>% dplyr::rename(id_0=id_polit),by=c("id_0"))

ggplot(ground_truth_sep,aes(nobs_max_bis)) +
  geom_histogram() +
  theme_bw()

janitor::tabyl(ground_truth_sep$nobs_max_bis)

ggplot(ground_truth_sep,aes(nobs_max_polit)) +
  geom_histogram() +
  theme_bw()

janitor::tabyl(ground_truth_sep$nobs_max_polit)



#------------------------
# (B) Calculate ROC Curve
#------------------------


# Note: Definitions based on https://en.wikipedia.org/wiki/Receiver_operating_characteristic

CalculateROCCurve(recordlinkagedata=rl_results_wide1_all %>%
                    select(Cluster_ID,Link_Score,id_0,id_1),
                  groundtruthdata = ground_truth_sep ,
                  idvars=c("id_0","id_1"),
                  tvvars=c("tset","vset"),
                  step_size = 0.02,
                  path=path,
                  path_rel = "04_Results/02_Record_linkage",
                  file_name = outfile_pdf_name) #,return_option=T)

outfile_pdf_name= "ROC_Curve_results3_only"

ground_truth_sep_uniques <- ground_truth_sep %>% 
  filter(nobs_max_bis==1 | is.na(id_1))

rl_results_wide1_all_uniques <- rl_results_wide1_all %>%
  select(Cluster_ID,Link_Score,id_0,id_1) %>%
  left_join(bisnode_data_nobs %>% 
              dplyr::rename(id_1=id_bis) %>%
              mutate(id_1=as.character(id_1)),by=c("id_1"))%>%
  left_join(nr_data_nobs %>% dplyr::rename(id_0=id_polit),by=c("id_0"))%>% 
  filter(nobs_max_bis==1| is.na(id_1))

CalculateROCCurve(recordlinkagedata=rl_results_wide1_all_uniques %>%
                    select(Cluster_ID,Link_Score,id_0,id_1),
                  groundtruthdata = ground_truth_sep_uniques ,
                  idvars=c("id_0","id_1"),
                  tvvars=c("tset","vset"),
                  step_size = 0.02,
                  path=path,
                  path_rel = "04_Results/02_Record_linkage",
                  file_name = outfile_pdf_name) #,return_option=T)




#--------------------
# (C) Explore results
#--------------------

# (i) Histogram of Score

View(rl_results)

ggplot(data=rl_results1,aes(x=Link_Score)) +
  geom_histogram()+
  theme_bw()

ggplot(data=rl_results2,aes(x=Link_Score)) +
  geom_histogram()+
  theme_bw()

# (ii) Create wide dataframe

rl_results_wide1 <- CreateWideDataset(data=KeepOnlyMatches(rl_results1),
                                      generate_differences=T)

rl_results_wide1 %>%   group_by(id_0) %>%
  mutate(nobs=n()) %>%
  filter(nobs>1) %>% 
  select(name_0 , name_1,  firstname_0, firstname_1,id_0,id_1,Link_Score, Cluster_ID) %>% 
  arrange(id_0)



# (iii) Histogram of link share among those with missing birthyear in Bisnode

# (a) All observations

ggplot(data=rl_results_wide1 %>% 
         filter(is.na(birthyear_1)),
       aes(x=Link_Score)) +
  geom_histogram()+
  theme_bw()


# (b) Observations with same name, firstname, and geo info

rl_results_wide1_sameinfo <- rl_results_wide1 %>%
  filter(name_0==name_1 & firstname_0== firstname_1 & 
           e_cntr_w_0==e_cntr_w_1 & n_cntr_w_0==n_cntr_w_1 &
           e_cntr_b_0==e_cntr_b_1 & n_cntr_b_0==n_cntr_b_1  )


ggplot(data=rl_results_wide1_sameinfo %>% 
         filter(is.na(birthyear_1)),
       aes(x=Link_Score)) +
  geom_histogram()+
  theme_bw()


# (iv) Relationship between geo distance and link score

# (a) All observations

ggplot(data=rl_results_wide1,
       aes(x=w_dist,y=Link_Score,color=birthyear_missing)) +
  geom_point()+
  theme_bw()



# (b) For those with the same firstname and lastname

rl_results_wide1_samename <- rl_results_wide1 %>%
  filter(name_0==name_1 & firstname_0== firstname_1)

ggplot(data=rl_results_wide1_samename,
       aes(x=w_dist,y=Link_Score,color=birthyear_missing)) +
  geom_point()+
  theme_bw()


#----------------------------------------------------------------------
# (D) Get precision and recall for specific cutoff value of Link_Score
#----------------------------------------------------------------------

# (i) Define specific cutoff value of Link_Score

ls_cutoff=0.001

# (ii) Calculate precision and recall

dataprep_precision_recall <- PreparationPrecisionRecall(recordlinkagedata=rl_results_wide2_all %>%
                                                          select(Cluster_ID,Link_Score,id_0,id_1),
                                                        #groundtruthdata = test_set,
                                                        groundtruthdata = ground_truth_sep,
                                                        idvars=c("id_0","id_1"),
                                                        tvvars=c("tset","vset"),
                                                        linkscore_cutoff = ls_cutoff)

# Quality check: look whether there are observations with more than one classi-
# fication or with not classification

dataprep_precision_recall %>% mutate(total = true_positives+true_negatives+false_positives+false_negatives) %>%
  filter(total>1)

dataprep_precision_recall %>% mutate(total = true_positives+true_negatives+false_positives+false_negatives) %>%
  filter(total==0)

dataprep_precision_recall %>% 
  filter((true_positives==0 & true_negatives==0 & false_positives==0 & false_negatives==0) |
           (is.na(true_positives) | is.na(true_negatives) | is.na(false_positives) | is.na(false_negatives)))

# Result: All fine. 

# check double entries

double_entries <- dataprep_precision_recall %>% 
  group_by(id_0,id_1)%>% 
  summarize(nobs=n()) %>%
  ungroup()%>%
  filter(nobs>1) 

double_entries %>% filter(!is.na(id_1))

# Result: All double entries are non-matches 

dataprep_precision_recall %>% 
  left_join(double_entries,by=c("id_0","id_1")) %>% 
  filter(nobs>1) %>%
  group_by(category)%>% 
  summarize(nobs_cat=n()) %>%
  ungroup()

# Result: All double entries are true negatives

# Conclusion: This test allows us to look at unique observations of id_0 and 
# id_1 to calculate the number ob observations in each of the four categories

GetPrecisionRecall(dataprep_precision_recall %>% filter(tset==1 | gt_data==0))
GetPrecisionRecall(dataprep_precision_recall %>% filter(vset==1 | gt_data==0))

# (iii) Look at four type of cases

# (a) True positives

df_true_pos <- LookUpConnection(data=rl_results_wide1_all %>% select(-Cluster_ID,-Link_Score),
                                entry=dataprep_precision_recall %>% 
                                  filter(true_positives==1),type="pair") %>%
  select(Cluster_ID,id_0,id_1,Link_Score,name_1,name_0,firstname_1,firstname_0,sex_1,sex_0,w_dist,
         b_dist,gt_data,tset,vset,rl_data,true_positives,true_negatives,false_positives,false_negatives)

# (b) True negatives

df_true_negs <-LookUpConnection(data=rl_results_wide1_all %>% select(-Cluster_ID,-Link_Score),
                                entry=dataprep_precision_recall %>% 
                                  distinct(id_0,id_1,.keep_all = T)%>%
                                  filter(true_negatives==1),type="pair") %>%
  select(Cluster_ID,id_0,id_1,Link_Score,name_1,name_0,firstname_1,firstname_0,sex_1,sex_0,w_dist,
         b_dist,gt_data,tset,vset,rl_data,true_positives,true_negatives,false_positives,false_negatives)

# (c) False positives

df_false_pos <-LookUpConnection(data=rl_results_wide1_all %>% select(-Cluster_ID,-Link_Score),
                                entry=dataprep_precision_recall %>% 
                                  filter(false_positives==1),type="pair") %>%
  select(Cluster_ID,id_0,id_1,Link_Score,name_1,name_0,firstname_1,firstname_0,sex_1,sex_0,w_dist,
         b_dist,gt_data,tset,vset,rl_data,true_positives,true_negatives,false_positives,false_negatives)

# (d) False negatives

df_false_negs <- LookUpConnection(data=rl_results_wide1_all ,
                                  entry=dataprep_precision_recall%>% 
                                    distinct(id_0,id_1,.keep_all = T)%>%
                                    filter(false_negatives==1)%>%
                                    select(id_1,id_0,rl_data,tset,vset,gt_data,
                                           true_positives,true_negatives,
                                           false_positives,false_negatives),
                                  idvars=c("id_0","id_1"),
                                  type="not_matched") %>%
  mutate(w_dist=NA,b_dist=NA) %>%
  select(Cluster_ID,id_0,id_1,Link_Score,name_1,name_0,firstname_1,firstname_0,sex_1,sex_0,w_dist,
         b_dist,gt_data,tset,vset,rl_data,true_positives,true_negatives,false_positives,false_negatives)

# out for Simon, Mark, and Sandro

df_out <- bind_rows(df_true_pos,df_true_negs,df_false_pos,df_false_negs) 

fwrite(x=df_out,
       file="04_Results/02_Record_linkage/Results_Wide_RL_20210921.csv", 
       sep=";")

# Quality check for false negatives

false_negs_check <- false_negs %>% select(id_1) %>% mutate(check=1)

ground_truth_sep %>% left_join(false_negs_check,by=c("id_1")) %>%
  filter(check==1 ) %>%
  summarize(sum_check=sum(check))

false_negs_check %>% left_join(ground_truth_sep %>% mutate(check2=1),by=c("id_1")) %>%
  filter(check==1 ) %>%
  summarize(sum_check=sum(check2))

# Result: all fine. 739 false negatives in ground truth data, all with id_0

# individual look ups

LookUpConnection(data=rl_results_wide1_all_fullvars %>%  select(-Cluster_ID,-Link_Score),
                 entry=precision_recall[[1]] %>% filter(id_0=="SG-2007-0049"),
                 type="pair")

LookUpConnection(data=rl_results_wide1_all_fullvars %>%  select(-Cluster_ID,-Link_Score),
                 entry=precision_recall[[1]] %>% filter(id_0=="SG-2007-0049"),
                 type="observations")

LookUpConnection(data=rl_results_wide1_all_fullvars,
                 entry=precision_recall[[1]] %>%filter(id_1=="1013799"),
                 type="observations")

## OLDER STUFF

# (W) Read in other data files

#bisnode_data <- data.table::fread("C:/Schmidlu/Dropbox/Record Linkage/Data/VRMandate/bisnode_falsenegatives_final_july2021.csv", 
#                                  sep = ";", encoding = "UTF-8")

#bisnode_data <- readstata13::read.dta13("./02_Processed_data/11_Directors_1994_2018/Bisnode_Person-Firmen_Geo.dta") %>%
#  select(ID,name, firstname, birthyear, sex, E_CNTR_w, N_CNTR_w,  E_CNTR_b, N_CNTR_b, language_w,eintrittdatum ,austrittdatum)

#nr_data <- data.table::fread("C:/Schmidlu/Dropbox/Record Linkage/Data/Politikerdaten/nationalraete_export_record_linkage_july2021.csv", 
#                             sep = ",", encoding = "UTF-8")
  
# (X) Check overlap test and validation set

# rl_results_wide1_all_fullvars %>% filter(id_0==overlap[1])
# ground_truth %>% filter(ID==overlap[1])
# 
# test_val_set %>% janitor::tabyl(vset,tset)
# 
# overlap <- test_val_set %>% filter(vset==1 & tset==1) %>% select(id_0) %>% pull()
# 
# test_set %>% filter(id_0==overlap[1])
# validation_set %>% filter(id_0==overlap[1])
# 
# 
# test_val_set <- validation_set %>% distinct(id_0) %>% mutate(vset=1) %>% 
#   full_join(test_set %>% distinct(id_0) %>% mutate(tset=1),by=c("id_0")) %>%
#   mutate(vset=ifelse(is.na(vset),0,1),
#          tset=ifelse(is.na(tset),0,1))
# 
# 
# CheckOverlap <- function(data1,data2,overlap){
#   for (i in 1:length(overlap)){
#     check1 <- data1 %>% filter(id_0==overlap[i]) 
#     check2 <- data2 %>% filter(id_0==overlap[i]) 
#     if (dim(check1)[1]==0){
#       print(paste0("Problem with id_0=",overlap[i])," in first dataset.")
#     }else if (dim(check2)[1]==0){
#       print(paste0("Problem with id_0=",overlap[i])," in second dataset.")
#     } else{
#       #print(paste0("No problem with id_0=",overlap[i],"."))
#     }
#   }
# }
# 
# CheckOverlap(data1=test_set,data2 = validation_set,overlap=overlap)
# 
# # Result: All as expected, all observations in overlap are in test and validation set. 

