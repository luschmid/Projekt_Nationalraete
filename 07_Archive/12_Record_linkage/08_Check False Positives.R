


# 1. Set path and load libraries and functions

path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte"
setwd(path)

library(tidyverse)
require(readxl)

source("./03_Code/12_Record_linkage/00_Functions.R")

# 2. Read in outfiles from round 1 and append them

path_rawdata1=paste0(path,"/02_Processed_data/12_Record_linkage/01_Bisnode/02_Round1_In_from_Coders/01_Sandro")
url_xlsx1 <- dir(path_rawdata1,pattern =".xlsx")
path_rawdata2=paste0(path,"/02_Processed_data/12_Record_linkage/01_Bisnode/02_Round1_In_from_Coders/02_Peter_Wallmueller")
url_xlsx2 <- dir(path_rawdata2,pattern =".xlsx")

df1 <- lapply(url_xlsx1, read_xlsx_files,path_rawdata=path_rawdata1) %>%
  bind_rows() %>%
  select(ID,Sicher,Unsicher,sourcefile)
df2 <- lapply(url_xlsx2, read_xlsx_files,path_rawdata=path_rawdata2 ) %>%
  bind_rows()%>%
  select(ID,Sicher,Unsicher,sourcefile)

df <- bind_rows(df1,df2)

# Correct special characters <U+FEFF>

df <- df %>% mutate(Sicher=ifelse(ID=="AG-2011-0042","3229636/3561882",Sicher),
                    Sicher=ifelse(ID=="BEJU-2003-0185","1831760/3266973/2487676",Sicher),
                    Unsicher=ifelse(ID=="BEJU-2003-0185","895126/380366/172062",Unsicher),
                    Sicher=ifelse(ID=="GE-1995-0089","1025330/3376737",Sicher),
                    Sicher=ifelse(ID=="SO-2011-0087","3395418/1346114",Sicher),
                    Sicher=ifelse(ID=="VD-2003-0059","2207836/2222851/2725321",Sicher),
                    Sicher=ifelse(ID=="ZH-1991-0020","572454/1142784",Sicher),
                    Sicher=ifelse(ID=="ZH-1999-0486","15463/197894/1400064/3300910",Sicher),
                    Sicher=ifelse(ID=="ZH-2015-0340","837579/3253006",Sicher),
                    Sicher=ifelse(ID=="ZH-2015-0340","837579/3253006",Sicher),
                    Sicher=ifelse(ID=="ZH-2015-0371","1172971/1468890/2245556/2529228",Sicher),
                    Sicher=ifelse(ID=="ZH-2015-0630","3037158/3543588",Sicher),
                    Unsicher=ifelse(ID=="ZH-2015-0630","3305822",Unsicher),
                    Unsicher=ifelse(ID=="SZ-2011-0034","1329429/464949",Unsicher),
                    Unsicher=ifelse(ID=="TG-1991-7080","211298/625250",Unsicher),
                    Unsicher=ifelse(ID=="TG-1995-0068","239778/1224184/1224186/1269231",Unsicher),
                    Sicher=ifelse(ID=="ZH-1999-0053","802346",Sicher),
                    Unsicher=ifelse(ID=="ZH-1999-0053","33182/671904/41826/748614/522718/972157/1347368",Unsicher),
                    Unsicher=ifelse(ID=="ZH-2007-0178","1838491/2277664",Unsicher),
                    Unsicher=ifelse(ID=="GE-1991-0093","1549418/3376771/2471700/1804103/1465568/1296616/1405574/1979334",Unsicher),
                    Unsicher=ifelse(ID=="ZH-1991-0565","252698/32046/839610/2346925/3185999",Unsicher),
                    Unsicher=ifelse(ID=="ZH-1995-0632","1234374/688299/3520451",Unsicher)
)

id_check = "TG-2015-0067"

df %>% filter(ID==id_check) %>% select(Unsicher)
df %>% filter(ID==id_check) %>% select(Sicher)

# 3. Recode Sicher/Unsicher delimiter

df$Sicher <- gsub('\\\\', '/',df$Sicher)
df$Sicher <- gsub('\\r\\n', '/',df$Sicher)

df$Unsicher <- gsub('\\\\', '/',df$Unsicher)
df$Unsicher <- gsub('\\r\\n', '/',df$Unsicher)


# 4. Split by Sicher/Unsicher by delimiter and reshape to long dataset
# Note: we keep only entries (connections), not non-entries

df_all_sicher <- bind_cols(df,as.data.frame(str_split_fixed(df$Sicher, pattern="/",n=100)))
df_all_unsicher <- bind_cols(df,as.data.frame(str_split_fixed(df$Unsicher, pattern="/",n=100)))

df_all_sicher_long <- df_all_sicher %>%   
  pivot_longer(!ID & !Sicher & !Unsicher & !sourcefile , names_to = "var", values_to = "id_1") %>%
  rename(id_0=ID) %>%
  filter(id_1!="") %>%
  select(id_0,id_1)

df_all_sicher_long %>% filter(id_0==df_all_sicher$ID[3])

df_all_unsicher_long <- df_all_unsicher %>%   
  pivot_longer(!ID & !Sicher & !Unsicher & !sourcefile, names_to = "var", values_to = "id_1") %>%
  rename(id_0=ID) %>%
  filter(id_1!="") %>%
  select(id_0,id_1)

df_final <- bind_rows(df_all_sicher_long,df_all_unsicher_long)
df_final$id_1_num <- as.numeric(df_final$id_1)
View(df_final %>% filter(is.na(id_1_num)))
df_final %>% filter(is.na(id_1_num)) %>% print(n=100) 

df_final <- df_final %>% filter(!is.na(id_1_num))%>%
  select(id_0,id_1_num) %>%
rename(id_1=id_1_num) %>%
  mutate(checked_firstround=1)


# 5. Merge with RL output

results_rl <- data.table::fread(file="04_Results/02_Record_linkage/Results_Wide_RL_20210921.csv", 
       sep=";",encoding = "UTF-8")

false_pos <- results_rl %>% filter(false_positives==1) %>% distinct(id_0)

df_check <- false_pos %>% left_join(df_final,by=c("id_0","id_1"))

sum(df_check$checked_firstround,na.rm=T)


# 6. Search for connections for false positives (as discussed with Simon, Oct 2021)

# (a) Prepare NR data


nr_data <- readstata13::read.dta13("./02_Processed_data/02_Elections_1971_2015/nationalraete_1931_2015_Geo.dta") %>%
  rename(id_0=ID)

nr_data_wide <- readstata13::read.dta13("./02_Processed_data/11_Directors_1994_2018/bisnode_out_wide.dta") %>%
  rename(id_0=ID)

nr_data_false_pos<- nr_data %>% left_join(false_pos %>% mutate(fp=1),by=c("id_0")) %>%
  filter(fp==1) %>%
  select(id_0, canton, sex, name, firstname , birthyear , gdename, job, list, 
         E_CNTR_w,N_CNTR_w, starts_with("E_CNTR_b"), starts_with("N_CNTR_b")) %>% 
  as_tibble()

nr_data_false_pos_single_buergerort <- nr_data_false_pos %>% 
  filter(is.na(E_CNTR_b2 )==T)  %>%
  rename(E_CNTR_b=E_CNTR_b1,N_CNTR_b=N_CNTR_b1 ) %>%
  select(id_0, canton, sex, name, firstname , birthyear , gdename, job, list, 
         E_CNTR_w,N_CNTR_w,E_CNTR_b,N_CNTR_b) %>%
  distinct(id_0,name, firstname,E_CNTR_w,N_CNTR_w,.keep_all = T) %>%
  arrange(id_0,name)
  

nr_data_false_pos_multiple_buergerort <- nr_data_false_pos %>% 
  filter(is.na(E_CNTR_b2 )==F)

nr_data_false_pos_multiple_buergerort_long <- nr_data_false_pos_multiple_buergerort %>%
  pivot_longer(!id_0 & !canton & !sex & !name & !firstname & !birthyear & !job & 
                 !gdename & !list & !E_CNTR_w & !N_CNTR_w, names_to = "var", values_to = "test") %>%
  filter(is.na(test)==F) %>%
  distinct(id_0,name, firstname,E_CNTR_w,N_CNTR_w,var,.keep_all = T) %>%
  arrange(id_0,name)

nr_e_cntr <- nr_data_false_pos_multiple_buergerort_long %>% select(id_0,var,test) %>% 
  filter(var %in% c("E_CNTR_b1","E_CNTR_b2","E_CNTR_b3","E_CNTR_b4","E_CNTR_b5","E_CNTR_b6")) %>%
  rename(E_CNTR_b=test)

nr_n_cntr <- nr_data_false_pos_multiple_buergerort_long %>% select(id_0,var,test) %>% 
  filter(var %in% c("N_CNTR_b1","N_CNTR_b2","N_CNTR_b3","N_CNTR_b4","N_CNTR_b5","N_CNTR_b6")) %>%
  rename(N_CNTR_b=test) %>% 
  mutate(mergevar=substr(var,8,9))


nr_cntr_all <- nr_e_cntr %>% full_join(nr_n_cntr, by=c("id_0","mergevar")) %>%
  select(id_0,E_CNTR_b,N_CNTR_b) %>% 
  mutate(mergevar=substr(var,8,9))

nr_cntr_all %>% filter(is.na(E_CNTR_b) | is.na(N_CNTR_b))

df_multiple_buergerorte <- CreateDFAllBuergerorte(data=nr_data_false_pos_multiple_buergerort,
                       data_buergerorte=nr_cntr_all) %>%
  select(-E_CNTR_b1,-E_CNTR_b2,-E_CNTR_b3,-E_CNTR_b4,-E_CNTR_b5,-E_CNTR_b6,
         -N_CNTR_b1,-N_CNTR_b2,-N_CNTR_b3,-N_CNTR_b4,-N_CNTR_b5,-N_CNTR_b6)


df_tosearch <- bind_rows(nr_data_false_pos_single_buergerort,
                         df_multiple_buergerorte) %>%
  distinct(id_0,name, firstname,E_CNTR_w,N_CNTR_w,E_CNTR_b,N_CNTR_b,.keep_all = T)

# (b) Search for connections in Bisnode data

bisnode_data <-data.table::fread(file="./02_Processed_data/11_Directors_1994_2018/bisnode_export_record_linkage_oct2021.csv", 
                                 sep=",",encoding = "UTF-8")

FalsePos_ToCheck <- GetResultBisnode(data_target=bisnode_data, data_input=df_tosearch,max_distance_input=0.2)


bisnode_data_firms <-data.table::fread(file="./02_Processed_data/11_Directors_1994_2018/bisnode_firms.csv", 
                                 sep=",",encoding = "UTF-8") %>% distinct(personenid,firma)
  ungroup()


# (c) Add all firms from Bisnode data

FalsePos_ToCheck_ID_VR <- FalsePos_ToCheck %>% distinct(ID_VR) %>% pull()

bisnode_data_firms_wide <- WriteAllVarsToWideFormat(data=bisnode_data_firms %>% 
                                                      filter(personenid %in% FalsePos_ToCheck_ID_VR),
                                                    idvar="personenid",
                                                    targetvars=c("firma"))

bisnode_data_firms_wide <- bisnode_data_firms_wide %>% rename(ID_VR=personenid)


# (d) Add all name information from NR data


FalsePos_ToCheck_ID_Pol <- FalsePos_ToCheck %>% distinct(ID_Pol) %>% pull()

nr_data_name <- readstata13::read.dta13(file="./02_Processed_data/nationalraete_1931_2015.dta") %>%
  select(ID,name) %>%
  distinct(ID,name)

nr_data_firstname <- readstata13::read.dta13(file="./02_Processed_data/nationalraete_1931_2015.dta") %>%
  select(ID,firstname) %>%
  distinct(ID,firstname)


nr_wide_name  <- WriteAllVarsToWideFormat(data=nr_data_name,
                                              idvar="ID",
                                              targetvars=c("name"))

nr_wide_firstname <- WriteAllVarsToWideFormat(data=nr_data_firstname,
                                           idvar="ID",
                                           targetvars=c("firstname"))

bisnode_data %>% filter(id_bis==2450156
)

# (e) Evaluate the connections

FalsePos_ToCheck_Unique <- FalsePos_ToCheck %>%
  distinct(ID_Pol,ID_VR, Vorname_Pol, Nachname_Pol, Vorname_VR, Nachname_VR,
           W_Distanz, B_Distanz, .keep_all = T)%>%
  arrange(ID_Pol,Vorname_Pol,Nachname_Pol)%>%
  left_join(bisnode_data_firms_wide %>% rename(ID_VR=personenid) ,by=c("ID_VR"))%>%
  left_join(nr_wide_firstname %>% rename(ID_Pol=ID),by=c("ID_Pol"))%>%
  left_join(nr_wide_name %>% rename(ID_Pol=ID),by=c("ID_Pol")) %>%
  select(-Vorname_Pol,-Nachname_Pol) %>%
  rename(Vorname_Pol=firstname,Nachname_Pol=name ) %>%
  select( ID_Pol,ID_VR,Vorname_VR,Nachname_VR,Vorname_Pol, Nachname_Pol,
          everything())

data.table::fwrite(x=FalsePos_ToCheck_Unique,
       file="02_Processed_data/12_Record_linkage/01_Bisnode/FalsePosCheck.csv", 
       ,sep=";")

# Note: This file was sent to Ladina Oester in October 2021.

# (f) Read-in first round results
# Note: This file was received from Ladina Oester in October 2021.


df <- read_excel("02_Processed_data/12_Record_linkage/01_Bisnode/FalsePosCheck_BackIn.xlsx")

library(lubridate)

age_df <- bisnode_data %>% select(id_bis,eintrittdatum_bis,austrittdatum_bis) %>%
  mutate(eintrittsdatum_bis_date=ymd(eintrittdatum_bis),
         austrittdatum_bis_date=ymd(austrittdatum_bis)) %>%
  group_by(id_bis) %>%
  summarize(eintrittdatum_bis_min=min(eintrittsdatum_bis_date,na.rm=T),
            austrittdatum_bis_max=max(austrittdatum_bis_date,na.rm=T))

bisnode_data %>% filter(id_bis==2)
bisnode_data %>% filter(id_bis==3)
bisnode_data %>% filter(id_bis==100)
age_df %>% filter(id_bis==100)
bisnode_data %>% filter(id_bis==101)
age_df %>% filter(id_bis==101)
bisnode_data %>% filter(id_bis==1001)
age_df %>% filter(id_bis==1001)

# Note: We create the minimum of the entry date of mandates per id (Eintrittsdatum)
#       and the maximum of the exit data of mandates per id (Austrittsdatum) to 
#       calculate the minimal implied entry age and the maximum implied exit
#       date of a possible connection.

df_second_round <- df %>% filter(is.na(Verbindung)) %>%
  left_join(nr_data %>% select(id_0,birthyear),by=c("ID_Pol"="id_0")) %>%
  left_join(bisnode_data %>% select(id_bis,birthyear_bis),by=c("ID_VR"="id_bis"))%>%
  left_join(age_df ,by=c("ID_VR"="id_bis"))%>%
  mutate(min_age_implied=year(eintrittdatum_bis_min)-birthyear,
         max_age_implied=year(austrittdatum_bis_max)-birthyear)
  rename(Geburtsjahr_Pol=birthyear,Geburtsjahr_VR=birthyear_bis)
  
hist(df_second_round$min_age_implied)
hist(df_second_round$max_age_implied)

  
  
# results_rl %>% filter(id_0=="AG-1983-0074")
# FalsePos_ToCheck_Unique %>% filter(ID_Pol=="AG-1983-0074" & ID_VR==1236988)
# FalsePos_ToCheck_Unique %>% filter(ID_Pol=="AG-2003-0050" & ID_VR==127947)
# 
# fp_all <- results_rl %>% filter(false_positives==T) %>% select(id_0,id_1)
# 
# df_check <- FalsePos_ToCheck_Unique %>% filter(ID_Pol %in% fp_all$id_0 & ID_VR %in% fp_all$id_1)
# 
# not_matched <- df_check %>% rename(id_0=ID_Pol,id_1=ID_VR) %>% mutate(onlycheck=1)%>%
#   full_join(fp_all %>% mutate(onlyfp=1), by=c("id_0","id_1")) %>%
#   filter(is.na(onlyfp) ) %>% select(id_0,id_1)
# 
# FalsePos_ToCheck_Unique %>% filter(ID_Pol %in% not_matched$id_0 & ID_VR %in% not_matched$id_1) %>%
#   distinct(ID_Pol,ID_VR)

# df_final %>% filter(id_0=="AG-2011-0042") 
# df %>% filter(ID=="AG-2011-0042") 
# df_final %>% filter(id_0=="AG-1995-9060") 
# old code: regex that did not work

#df <- replace_pattern_loop(data=df,variable="Sicher",
#                           pattern="[<U\\+FEFF>]",pattern_replace=""
#df$Sicher <- gsub('<U\\+FEFF>', '',df$Sicher, useBytes = T) 
#df$Sicher <- str_replace_all(df$Sicher[i], "[<U\\+FEFF>]", "")

#test <- gsub('[\\<U\\+FEFF\\>]', '',"3229636/<U+FEFF>3561882")
#df %>% filter(ID=="AG-2011-0042") %>% select(Sicher)
#gsub('[\\<U\\+FEFF\\>]', '',df$Sicher[1012])
#df[1012,]#df$Unsicher <- gsub('[<U\\+FEFF>]', '',df$Unsicher)
