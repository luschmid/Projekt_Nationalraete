library(tidyverse) 

path <- "C:/Schmidlu/Dropbox/Record Linkage"

# (A) Test and validation set Bisnode ----

bis_test <- data.table::fread(paste0(path,"./02_Data/01_Ground_Truth_Bisnode/test_set_randomized_buckets_0.65.csv"),
                  sep = ",", encoding = "UTF-8") 
              
bis_test %>% distinct(ID)


bis_val <- data.table::fread(paste0(path,"./02_Data/01_Ground_Truth_Bisnode/validation_set_randomized_buckets_0.35.csv"),
                              sep = ",", encoding = "UTF-8") 

bis_val %>% distinct(ID)


# (B) Test and validation set Sugarcube ----

sug_test <- data.table::fread(paste0(path,"./02_Data/02_Ground_Truth_Sugarcube/sugarcube_test_set_randomized_buckets_0.65.csv"),
                              sep = ",", encoding = "UTF-8") 

sug_test %>% distinct(id_polit)


sug_val <- data.table::fread(paste0(path,"./02_Data/02_Ground_Truth_Sugarcube/sugarcube_validation_set_randomized_buckets_0.35.csv"),
                             sep = ",", encoding = "UTF-8") 

sug_val %>% distinct(id_polit)




