rm(list = ls())

# (A) Choose recordlinkage data, groundtruth data and linkscore_cutoff  ----

# Note: rl_data data is either rl_results_gen3_wide or rl_results_postprocessing)
#        groundtruth data is either test_set or validation_set


rl_data <- "rl_results_gen3_wide"

gt_data <- "test_set" 

linkscore_cutoff <- 0.1


# (B) Libraries, working directory and load data ----

# (i) Libraries, working directory, and load functions

library(readxl)
library(tidyverse)
library(data.table)
library(haven)
library(xtable)

path <- "C:/Schmidlu/Dropbox/Record Linkage"
setwd(path)

nr_path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte"

source("./01_Code/03_Record_Linkage_Results/00_Functions.R")

infile_stata_name <- "RL_Output_Round1_fp_maj_IDs.dta"
source="sugarcube"
step_size=0.02
postprocessing=1


# (C) Read in datasets ----


# (i) Read in rl result files ----

rl_results_postprocessing <-  read_dta(paste0(nr_path,
                                              "/02_Processed_data/12_Record_linkage/02_Sugarcube/11_RL_Output_Round1_In/",infile_stata_name))%>%
  dplyr::rename(Link_Score=link_score,
                e_id_0=e_cntr_w_0, 
                n_id_0=n_cntr_w_0,
                e_id_1=e_cntr_w_1, 
                n_id_1=n_cntr_w_1) %>%
  mutate(Cluster_ID=row_number(),
         year_1=as.character(year_1))%>%
  mutate(time_period = case_when(
    year_1 == 1933 ~ 1 ,
    year_1 == 1942 ~ 2 ,
    year_1 == 1959 ~ 3 ,
    year_1 == 1961 ~ 4 ,
    year_1 == 1962 ~ 5 ,
    year_1 == 1963 ~ 6 ,
    year_1 == 1965 ~ 7 ,
    year_1 == 1966 ~ 8 ,
    year_1 == 1969 ~ 9 ,
    year_1 == 1972 ~ 10,
    year_1 == 1975 ~ 11,
    year_1 == 1979 ~ 12,
    year_1 == 1980 ~ 13,
    year_1 == 1981 ~ 14,
    year_1 == 1982 ~ 15,
    year_1 == 1983 ~ 16,
    year_1 == 1984 ~ 17,
    year_1 == 1985 ~ 18,
    year_1 == 1986 ~ 19,
    year_1 == 1987 ~ 20,
    year_1 == 1988 ~ 21,
    year_1 == 1989 ~ 22,
    year_1 == 1990 ~ 23,
    year_1 == 1991 ~ 24,
    year_1 == 1992 ~ 25,
    year_1 == 1993 ~ 26,
    year_1 == 1994 ~ 27,
    year_1 == 1995 ~ 28,
    year_1 == 1996 ~ 29,
    year_1 == 1997 ~ 30,
    year_1 == 1998 ~ 31,
    year_1 == 1999 ~ 32,
    year_1 == 2000 ~ 33,
    year_1 == 2001 ~ 34,
    year_1 == 2002 ~ 35,
    year_1 == 2003 ~ 36),
    id_1_all=paste0(year_1,"-",id_1)
  )

rl_results_gen3_long <-  ReadinResults("RL-Results-15122022-100518bc-12f0-4d5d-8b42-1f4484952a2a.csv", 
                                  path_rel="02_Data/08_Output_RL_Sugarcube/", source = "sugarcube")

rl_results_gen3_wide <- CreateWideDataset(data=rl_results_gen3_long ,
                                          source="sugarcube") %>%
  mutate(year_1=year) %>%
  mutate(time_period = case_when(
    year_1 == 1933 ~ 1 ,
    year_1 == 1942 ~ 2 ,
    year_1 == 1959 ~ 3 ,
    year_1 == 1961 ~ 4 ,
    year_1 == 1962 ~ 5 ,
    year_1 == 1963 ~ 6 ,
    year_1 == 1965 ~ 7 ,
    year_1 == 1966 ~ 8 ,
    year_1 == 1969 ~ 9 ,
    year_1 == 1972 ~ 10,
    year_1 == 1975 ~ 11,
    year_1 == 1979 ~ 12,
    year_1 == 1980 ~ 13,
    year_1 == 1981 ~ 14,
    year_1 == 1982 ~ 15,
    year_1 == 1983 ~ 16,
    year_1 == 1984 ~ 17,
    year_1 == 1985 ~ 18,
    year_1 == 1986 ~ 19,
    year_1 == 1987 ~ 20,
    year_1 == 1988 ~ 21,
    year_1 == 1989 ~ 22,
    year_1 == 1990 ~ 23,
    year_1 == 1991 ~ 24,
    year_1 == 1992 ~ 25,
    year_1 == 1993 ~ 26,
    year_1 == 1994 ~ 27,
    year_1 == 1995 ~ 28,
    year_1 == 1996 ~ 29,
    year_1 == 1997 ~ 30,
    year_1 == 1998 ~ 31,
    year_1 == 1999 ~ 32,
    year_1 == 2000 ~ 33,
    year_1 == 2001 ~ 34,
    year_1 == 2002 ~ 35,
    year_1 == 2003 ~ 36),
    id_1_all=paste0(year_1,"-",id_1)
  )

if (rl_data == "rl_results_gen3_wide"){
  recordlinkagedata <-   rl_results_gen3_wide
}

if (rl_data == "rl_results_postprocessing"){
  recordlinkagedata <-  rl_results_postprocessing
}



# (ii) Read in politician data and ground truth files ----

# Note: Keep only those in ground truth that are in NR data.

nr_data <- data.table::fread("./02_Data/04_Politicians_Sugarcube/nationalraete_export_record_linkage_sugarcube_july2021.csv",
                             sep = ",", encoding = "UTF-8"
)%>%
  mutate(
    id_0 = as.character(id_polit),
    e_id_0 = as.character(e_id_polit),
    n_id_0 = as.character(n_id_polit)
  ) 

id_0_in_nr_data <- nr_data %>%
  select(id_0) %>%
  distinct() %>%
  pull()

validation_set <-  ReadinGroundTruthFiles(file_name = "sugarcube_test_set_randomized_buckets_0.65.csv",
                                          separator=",",
                                          path_rel="./02_Data/02_Ground_Truth_Sugarcube/",
                                          source = "sugarcube")%>%
  filter(id_0 %in% id_0_in_nr_data)

#validation_set %>% filter(is.numeric(e_id_1)%%1>0)
#validation_set %>% filter(is.numeric(n_id_1)%%1>0)

test_set <-  ReadinGroundTruthFiles(file_name = "sugarcube_validation_set_randomized_buckets_0.35.csv",
                                    separator=",",
                                    path_rel="./02_Data/02_Ground_Truth_Sugarcube/",
                                    source = "sugarcube")%>%
  filter(id_0 %in% id_0_in_nr_data)

#test_set %>% filter(is.numeric(e_id_1)%%1>0)
#test_set %>% filter(is.numeric(n_id_1)%%1>0)


if (gt_data == "test_set"){
  groundtruthdata <-   test_set
}

if (gt_data == "validation_set"){
  groundtruthdata <-  validation_set
}


# (iii) Get id pairs of perfect matches ----

sugarcube_data <- data.table::fread("./02_Data/06_Directors_Sugarcube/sugarcube_export_record_linkage_july2022.csv",
                                    sep = ",", encoding = "UTF-8") %>%
  mutate(time_period = case_when(
    year_sug == 1933 ~ 1 ,
    year_sug == 1942 ~ 2 ,
    year_sug == 1959 ~ 3 ,
    year_sug == 1961 ~ 4 ,
    year_sug == 1962 ~ 5 ,
    year_sug == 1963 ~ 6 ,
    year_sug == 1965 ~ 7 ,
    year_sug == 1966 ~ 8 ,
    year_sug == 1969 ~ 9 ,
    year_sug == 1972 ~ 10,
    year_sug == 1975 ~ 11,
    year_sug == 1979 ~ 12,
    year_sug == 1980 ~ 13,
    year_sug == 1981 ~ 14,
    year_sug == 1982 ~ 15,
    year_sug == 1983 ~ 16,
    year_sug == 1984 ~ 17,
    year_sug == 1985 ~ 18,
    year_sug == 1986 ~ 19,
    year_sug == 1987 ~ 20,
    year_sug == 1988 ~ 21,
    year_sug == 1989 ~ 22,
    year_sug == 1990 ~ 23,
    year_sug == 1991 ~ 24,
    year_sug == 1992 ~ 25,
    year_sug == 1993 ~ 26,
    year_sug == 1994 ~ 27,
    year_sug == 1995 ~ 28,
    year_sug == 1996 ~ 29,
    year_sug == 1997 ~ 30,
    year_sug == 1998 ~ 31,
    year_sug == 1999 ~ 32,
    year_sug == 2000 ~ 33,
    year_sug == 2001 ~ 34,
    year_sug == 2002 ~ 35,
    year_sug == 2003 ~ 36),
    id_1_all=paste0(year_sug,"-",id_sug)
  )

sugarcube_data %>% filter(id_1_all %in% c("1991-112413", "1989-158857", "1995-2436", 
                          "1995-2436", "1986-3948", "1969-45255", "1969-45255", 
                          "1992-56200", "1995-77810", "1988-87936", "1989-93820", 
                          "1969-9682"))


perfectmatches <- nr_data %>%
  filter(id_polit %in% validation_set$id_0) %>%
  inner_join(sugarcube_data,
             by = c(
               "name_polit" = "name_sug",
               "firstname_polit" = "firstname_sug",
               "n_cntr_w_polit" = "gdenr_n_cntr_sug",
               "e_cntr_w_polit" = "gdenr_e_cntr_sug"
             )) %>%
  mutate(id_0 = id_polit, 
         id_1 = as.character(id_sug), 
         indi = 1) %>%
  mutate(n_id_1=n_id_0,
         e_id_1=e_id_0,
         year=as.character(year_sug),
         id_1 = as.character(floor(as.numeric(id_1)))) %>%
  select(id_0, id_1, year,e_id_0,n_id_0,e_id_1,n_id_1,indi)



# (D) Run classification by Lukas ----

# Note: This is after Lukas and Sandro have compared their classification
#       in March 2022

idvars = c("id_0", "id_1","year_1","e_id_0","n_id_0", "e_id_1","n_id_1")


# (i) Keep only matches above the linkscore cutoff


recordlinkagedata <- recordlinkagedata %>%
  select(Cluster_ID, Link_Score, id_0, id_1,year_1, e_id_0,n_id_0, e_id_1,n_id_1)

id_0_groundtruthdata <- groundtruthdata %>%
  distinct(id_0) %>%
  pull()


df_matches <- recordlinkagedata %>%
  filter(!is.na(id_1)) %>% # filter out only existing matches
  filter(is.element(id_0, id_0_groundtruthdata)) %>% # Note: focus only on politician ids in ground truth data
  filter(Link_Score >= linkscore_cutoff)%>%
  mutate(id_1 = as.character(id_1),
         e_id_0=as.character(e_id_0),
         n_id_0=as.character(n_id_0),
         e_id_1=as.character(e_id_1),
         n_id_1=as.character(n_id_1))


# (ii) Combine Record Linkage with Ground Truth Data and recode cases

df_all <- df_matches %>%
  mutate(rl_data = 1) %>%
  full_join(groundtruthdata %>% mutate(gt_data = 1, id_1 = as.character(id_1)),
            by = c(idvars)
  ) %>%
  mutate(
    gt_data = ifelse(is.na(gt_data), 0, 1),
    rl_data = ifelse(is.na(rl_data), 0, 1)
  ) %>%
  mutate(
    true_positives = ifelse(rl_data == 1 & gt_data == 1, 1, 0), # ids are in both datasets
    true_negatives = ifelse(rl_data == 0 & gt_data == 1 & is.na(id_1), 1, 0), # ids are only in gt data with missing id_1
    false_positives = ifelse(rl_data == 1 & gt_data == 0, 1, 0), # ids are only in rl data
    false_negatives = ifelse(rl_data == 0 & gt_data == 1 & !is.na(id_0) & !is.na(id_1), 1, 0), # ids are only in gt data with non-missing id_0 and id_1
    category = case_when(
      true_positives == 1 ~ "true positive",
      true_negatives == 1 ~ "true negative",
      false_positives == 1 ~ "false positive",
      false_negatives == 1 ~ "false negative"
    ),
    check = true_positives + true_negatives + false_positives + false_negatives
  ) %>%
  select(Cluster_ID, Link_Score, id_0, id_1,year_1, e_id_0,n_id_0, e_id_1,n_id_1,category,
         true_positives,true_negatives,false_positives,false_negatives)


# (iii) Merge nr and sugarcube information to rl results

df_final <- df_all %>% 
  left_join(sugarcube_data %>%
              select(id_sug,year_sug,e_id_sug,n_id_sug,name_sug,firstname_sug) %>%
              mutate(id_sug=as.character(id_sug),
                     year_sug=as.character(year_sug),
                     e_id_sug=as.character(e_id_sug),
                     n_id_sug=as.character(n_id_sug)), 
            by = c("id_1" = "id_sug",
                   "year_1" = "year_sug",
                   "e_id_1" = "e_id_sug",
                   "n_id_1" = "n_id_sug"
            )) %>% 
  left_join(nr_data %>%
              select(id_polit,e_id_polit,n_id_polit,name_polit, firstname_polit) %>%
              distinct(id_polit,e_id_polit,n_id_polit,.keep_all =T) %>%
              mutate(id_polit=as.character(id_polit),
                     e_id_polit=as.character(e_id_polit),
                     n_id_polit=as.character(n_id_polit)), 
            by = c("id_0" = "id_polit",
                   "e_id_0" = "e_id_polit",
                   "n_id_0" = "n_id_polit"
            )) %>%
  mutate(w_dist = sqrt((as.numeric(e_id_0) - as.numeric(e_id_1))^2 + (as.numeric(n_id_0) - as.numeric(n_id_1))^2))




# (E) Explore categories  ----


# (i) true positives ----

df_tp <- df_final %>% filter(category=="true positive")
df_tp_same <- df_tp %>% filter(w_dist==0 &name_polit == name_sug &  firstname_polit == firstname_sug )
df_tp_notsame <- df_tp %>% filter((w_dist!=0 | (name_polit != name_sug) | (firstname_polit != firstname_sug)))

out <- dim(df_tp_same)[1]/dim(df_tp)[1]



# (ii) false positives ----

df_fp <- df_final %>% filter(category=="false positive")
df_fp_same <- df_fp %>% filter(w_dist==0 &name_polit == name_sug &  firstname_polit == firstname_sug)
df_fp_notsame <- df_fp %>% filter((w_dist!=0 | (name_polit != name_sug) | (firstname_polit != firstname_sug)) )

#View(df_fp)

out <- c(out,dim(df_fp_same)[1]/dim(df_fp)[1])

df_fp_same <- df_fp %>% filter(name_polit == name_sug &  firstname_polit == firstname_sug)
dim(df_fp_same)[1]/dim(df_fp)[1]

mean(df_fp_same$w_dist,na.rm = T)

# RL results: 

# Validation set: 84% of false positives have the same first and the same lastname but a 
# different residence location. The average distance in the residence location is 
# 52.1 kilometers. 

# Test set: 85.5% of false positives have the same first and the same lastname but a 
# different residence location. The average distance in the residence location is 
# 51.67 kilometers. 


# Postprocessing results: 

# Validation set: 63% of false positives have the same first and the same lastname but a 
# different residence location. The average distance in the residence location is 
# 10.7 kilometers. 

# Test set: 71% of false positives have the same first and the same lastname but a 
# different residence location. The average distance in the residence location is 
# 19.17 kilometers. 


# (iii) true negatives ----

df_tn <- df_final %>% filter(category=="true negative")
df_tn_same <- df_tn %>% filter(w_dist==0 &name_polit == name_sug &  firstname_polit == firstname_sug)
df_tn_notsame <- df_tn %>% filter((w_dist!=0 | (name_polit != name_sug) | (firstname_polit != firstname_sug)) )

out <- c(out,dim(df_tn_same)[1]/dim(df_tn)[1])


# (iv) false negatives ----

df_fn <- df_final %>% filter(category=="false negative")
df_fn_same <- df_fn %>% filter(w_dist==0 &name_polit == name_sug &  firstname_polit == firstname_sug)
df_fn_notsame <- df_fn %>% filter((w_dist!=0 | (name_polit != name_sug) | (firstname_polit != firstname_sug)) )

out <- c(out,dim(df_fn_same)[1]/dim(df_fn)[1])

df_fn %>% head()

# Validation set: 8.8 % of all false negatives are identical in terms of name, firstname, and residence municipality. 
# Test set: 12.2% of all false negatives are identical in terms of name, firstname, and residence municipality. 


# (F) Output 

# (i) Relative share of categories ----


# df_out <- df_final %>% group_by(category) %>%
#   summarize(nobs=n()) %>%
#   mutate(share=nobs/dim(df_all)[1],
#          cat=case_when(category=="true positive" ~ 1,
#                        category=="false positive" ~ 2,
#                        category=="true negative" ~ 3,
#                        category=="false negative" ~ 4)) %>%
#   arrange(cat) %>%
#   select(share)
# 
# print(xtable::xtable(df_out), include.rownames=FALSE)

dataprep_vs <- df_final %>%
  group_by(id_0) %>%
  dplyr::summarize(
    true_positives = mean(true_positives),
    true_negatives = mean(true_negatives),
    false_positives = mean(false_positives),
    false_negatives = mean(false_negatives)
  ) %>%
  mutate(sum_all = true_positives + true_negatives + false_positives + false_negatives)


df <- dataprep_vs %>%
  summarize(
    true_positives = sum(true_positives, na.rm = T),
    false_positives = sum(false_positives, na.rm = T),
    true_negatives = sum(true_negatives, na.rm = T),
    false_negatives = sum(false_negatives, na.rm = T)
  )

df_sum <- df %>%
  summarize(n_all = true_positives + true_negatives + false_positives + false_negatives)


df/df_sum$n_all


precision_intensive <- df$true_positives / (df$true_positives + df$false_positives)
recall_intensive <- df$true_positives / (df$true_positives + df$false_negatives) # true positive rate (TPR)

# (ii) Share of observations identical in terms of name, firstname, and residence municipality  ----

print(xtable::xtable(data.frame(out)), include.rownames=FALSE)


# (iii) Sample output ----

# data.table::fwrite(
#   x = bind_rows(df_tp_same %>% mutate(same=1),df_tp_notsame %>% mutate(same=0),
#                 df_fp_same %>% mutate(same=1),df_fp_notsame %>% mutate(same=0),
#                 df_tn_same %>% mutate(same=1),df_tn_notsame %>% mutate(same=0),
#                 df_fn_same %>% mutate(same=1),df_fn_notsame %>% mutate(same=0)),
#   file = paste0(nr_path,"/12_Record_linkage/02_Sugarcube/20_Check_RL_Output/sample.csv"),
#   sep = ";",
#   append=F
# )


# (iv) Feedback Sandro (see email 28 August 2023) ----

# 2 a) False negatives: Check deduplication

ids_to_check <- df_fn_same %>% select(id_0) %>% pull() %>% unique()

id_to_check <- ids_to_check[12] # 8,9,15 

rl_results_gen3_wide %>%
  filter(id_0 %in% id_to_check) %>%
  arrange(year_1) %>%
  select(Link_Score,contains("name"),year_1,id_1,time_period,e_cntr_w_0,n_cntr_w_0,e_cntr_w_1,
         n_cntr_w_1) %>%
  print(n=30)
  
df_fn_same %>%
  filter(id_0 %in% id_to_check) %>%
  arrange(year_1) %>%
  select(id_0,contains("name"),year_1,id_1,e_id_0, n_id_0, e_id_1, n_id_1)

# 2 b) False negatives: Minuscule typos


id_to_check <-"BEJU-1955-0131" # "BEJU-1947-0032"  "BEJU-1975-9019", "BEJU-1951-0193"

rl_results_gen3_wide %>%
  filter(id_0 %in% id_to_check) %>%
  arrange(year_1) %>%
  select(Link_Score,contains("name"),year_1,id_1,time_period,e_cntr_w_0,n_cntr_w_0,e_cntr_w_1,
         n_cntr_w_1) %>%
  print(n=30)

sugarcube_data %>%
  filter(id_sug ==65430	& year_sug==1984) %>% 	
  arrange(year_sug) %>%
  select(firstname_sug,name_sug,year_sug,id_sug, e_id_sug, n_id_sug)

# 3) 

df_fn_same <- df_fn %>% filter(w_dist==0 &name_polit == name_sug &  firstname_polit == firstname_sug)%>%
  select(Link_Score,contains("name"),year_1,id_1,e_id_0,n_id_0,e_id_1,n_id_1) 


df_tp %>% filter(name_polit=="Delaloye")%>%
  select(Link_Score,contains("name"),year_1,id_1,e_id_0,n_id_0,e_id_1,n_id_1) %>%
  print(n=30)


rl_results_gen3_wide %>% filter(year_1==1996  & id_1== 36251)

rl_results_gen3_wide %>% 
