library(tidyverse)
library(janitor)

setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte")



# (i) Get Unique firstnames from Bisnode

sugarcube_data <- data.table::fread("C:/Schmidlu/Dropbox/Record Linkage/02_Data/06_Directors_Sugarcube/rl_sugarcube_input_data.csv",
                                  sep = ",", encoding = "UTF-8"
) 

sugarcube_data <-sugarcube_data %>% 
  group_by(firstname) %>%
  mutate(firstname_nobs=n()) %>%
  ungroup()

sugarcube_data %>% arrange(firstname_nobs)

sugarcube_distinct_names <- sugarcube_data %>% 
  distinct(firstname,firstname_nobs) %>%
  arrange(firstname)

data.table::fwrite(
  x = sugarcube_distinct_names,
  file ="./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/firstnames_unique.csv")

# Note: this file is used as an input in Python script "03_Classify_Gender_Sugarcube.py"


# (ii) Read in Python results for manual checks

infile_csv_names = list.files(path="./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/01_Name_Databases",pattern=paste0("firstnames_unique_notRecoded_"))

# Step 1: Get all files in Long format

out <- data.frame()
for (i in 1:length(infile_csv_names)){
new_file <- data.table::fread(paste0("./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/01_Name_Databases",infile_csv_names[i]),
                              sep = ";", encoding = "UTF-8") %>%
  as_tibble()
new_file$source <- substr(infile_csv_names[i],1,nchar(infile_csv_names[i])-28)
out <- bind_rows(out,new_file)%>%
  as_tibble()
}

# Step 2: Merge files in wide format

sources <- out %>% distinct(source) %>% pull()

for (i in 1:length(sources)){
out_temp <- out %>% filter(source==sources[i]) %>% dplyr::select(-source)
names(out_temp)[3:5] <- paste0(names(out_temp)[3:5],"_",i)
if (i==1){
out_all <- out_temp
} else{
out_temp <- out_temp %>% select(-CleanName)
out_all <- out_all %>% left_join(out_temp,by=c("OriginalName")) 
}
}

# Step 3: Recode gender

sugarcube_unique <- sugarcube_data %>% select(firstname,firstname_nobs) %>% 
  distinct(firstname,firstname_nobs) %>%
  arrange(firstname)

out_all <- out_all %>% 
  mutate(female=ifelse(female_1>0 | female_3>0,1,0)) %>% 
  mutate(male=ifelse(male_1>0 | male_3>0,1,0)) %>%
  mutate(conflict=ifelse(male==1 & female==1,1,0)) %>%
  mutate(undecided_1=ifelse((male_1>0 & female_1==0) | (female_1>0 & male_1==0),0,undecided_1)) %>%
  mutate(undecided_2=ifelse((male_2>0 & female_2==0) | (female_2>0 & male_2==0),0,undecided_2)) %>%
  mutate(undecided_3=ifelse((male_3>0 & female_3==0) | (female_3>0 & male_3==0),0,undecided_3)) %>%
  left_join(.,y=sugarcube_unique,
            by=c("OriginalName"="firstname")) %>%
  mutate(OriginalNameColl=StringReplace(OriginalName))

out_recoded <- out_all %>% 
  anti_join(out_tocheck %>% select(OriginalName), by=c("OriginalName"))

data.table::fwrite(out_tocheck %>% select(-male_2,-female_2,-undecided_2),
                   file ="./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/Gender_ToRecode.csv",
                   sep = ";") 

data.table::fwrite(out_recoded,
                   file ="./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/Gender_Recoded.csv",
                   sep = ";") 


# Step 4: Read-in manually recoded gender classification

# Note: This task was done by Daphne Caverzasio and Ladina Oester in June, July 2022 and supervised
# by Mark. The files are here: Dropbox\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\07_GenderCoding


Gender_Recoded <- data.table::fread("./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/Gender_Recoded.csv",
                              sep = ";", encoding = "UTF-8") %>%
  select(CleanName,OriginalName,male)

Gender_ToRecode <- data.table::fread("./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/Gender_ToRecode.csv",
                                    sep = ";", encoding = "UTF-8") %>%
  select(CleanName,OriginalName,male) %>%
  mutate(male=NA)

Gender_ManualRecode <-  read_dta("./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/Gender_ToRecoded_final.dta")%>%
  rename(CleanName=cleanname,
         male_manual=gender)

firstnames_unique <- data.table::fread(
  file ="./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/firstnames_unique.csv") %>%
  rename(OriginalName=firstname)


GenderAll  <- Gender_Recoded %>% 
  mutate(recoded=1) %>%
  bind_rows(Gender_ToRecode) %>%
  left_join(Gender_ManualRecode,by="CleanName")  %>%
  mutate(male=ifelse(is.na(recoded) & is.na(male),male_manual,male )) %>%
  select(OriginalName,male) %>%
  left_join(firstnames_unique,by="OriginalName") 
  

GenderAll %>% filter(OriginalName=="Adelheid gen. Heidi")


GenderAll %>% tabyl(male) 
GenderAll %>% group_by(male) %>% summarize(firstname_nobs_max=max(firstname_nobs,na.rm=T),
                                           firstname_nobs_mean=mean(firstname_nobs,na.rm=T)) 

View(GenderAll %>% filter(firstname_nobs>4 & is.na(male)) %>%
  arrange(firstname_nobs))


GenderNotRecoded <- GenderAll %>% filter(is.na(male))


# data.table::fwrite(GenderNotRecoded,
#                    file ="./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/GenderNotRecodedFinal.csv",
#                    sep = ";") 

data.table::fwrite(GenderAll %>% select(OriginalName,male),
                   file ="./02_Processed_data/10_Directors_1934_2003/07_GenderCoding/GenderRecodedFinal.csv",
                   sep = ";") 



# check for those that do not merge

test <- Gender_Recoded %>% 
  mutate(recoded=1) %>%
  bind_rows(Gender_ToRecode) %>%
  anti_join(Gender_ManualRecode,by="CleanName") %>%
  filter(is.na(recoded))

View(test)

# Result: These are only abbreviations (A.T., L.O) with no further information on the gender


# 
# out_all %>% distinct(OriginalNameColl)
# 
# out_all %>% arrange(-firstname_nobs) %>% print(n=100)
# 
# out_all %>% tabyl(undecided_1)
# out_all %>% tabyl(undecided_2)
# out_all %>% tabyl(undecided_3)
# 
# out_all %>% tabyl(undecided_1,undecided_2)
# out_all %>% tabyl(undecided_1,undecided_3)
# 
# out_tocheck <- out_all %>% filter(conflict==1 | (undecided_1>0 | undecided_3>0) | 
#                                     (male==0 & female==0)) 
# 
# out_tocheck %>% distinct(OriginalNameColl)
# 
# ggplot(data=out_tocheck %>% filter(firstname_nobs>1),aes(x=OriginalName,y=firstname_nobs)) +
#   geom_bar(stat="identity") +
#   coord_flip()
# 
# sum(out_tocheck$firstname_nobs,na.rm=T)
# sum(out_recoded$firstname_nobs,na.rm=T)



# 
# 
# vornamencom <- read_delim("./02_Data/06_Directors_Sugarcube/vornamencom_firstnames_unique_out.csv",
#                           delim = ";") 
# 
# 
# vornamen_firstnames_unique_notRecoded_out
# 
# vornamencom <- read_delim("./02_Data/06_Directors_Sugarcube/vornamencom_firstnames_unique_out.csv",
#   delim = ";") 
# names(vornamencom)[3:5] <- paste0(names(vornamencom)[3:5],"_vn")
# 
# 
# babyvornamen <- read_delim("./02_Data/06_Directors_Sugarcube/babyvornamen_firstnames_unique_out.csv",
#                         delim = ";") 
# names(babyvornamen)[3:5] <- paste0(names(babyvornamen)[3:5],"_ba")
# 
# gender_recoded <- vornamencom %>% left_join(babyvornamen,by=c("OriginalName")) %>%
#   select(-CleanName.y) %>%
#   rename(CleanName=CleanName.x)
# 
# 
# gender_recoded %>% tail(n=30) %>% print(n=30)
# 
# gender_recoded_unclear <- gender_recoded %>%
#   filter(undecided==1 | (male==0 & female==0))
# 
# gender_recoded_split <- gender_recoded %>%
#   filter(male>0 & female>0)
# 
# 
# bisnode_data %>% filter(str_detect(firstname,"Gampel AG"))
# bisnode_distinct_names %>% filter(str_detect(firstname,"Gampel AG"))
