
# (A) Libraries, paths, and working directory ----

library(tidyverse)
library(data.table)
library(dplyr, warn.conflict = FALSE, quietly = TRUE)
library(lubridate)
library(testit)
library(haven)
library(stringr)
library(readxl)
library(furrr)

rm(list = ls())
rm(list=setdiff(ls(all.names=TRUE), lsf.str(all.names=TRUE)))

setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte")


# (B) Read in data ----

data_sug <- fread("./02_Processed_data/10_Directors_1934_2003/Sugarcube_RLPostProcessing_Persons-Firmnames.csv",
                    sep = ",", encoding = "UTF-8") %>% 
  select(year_sug,firstname,lastname,firmnames,gdename,
         GdeNr_E_CNTR,GdeNr_N_CNTR,id_sug,male)%>%
  rename(lastname_1=lastname,firstname_1=firstname,gdename_1=gdename,
         firmnames_1=firmnames,year_1=year_sug,id_1=id_sug, male_1=male,
         e_cntr_w_1=GdeNr_E_CNTR,n_cntr_w_1=GdeNr_N_CNTR)

data_fn <- haven::read_dta("./02_Processed_data/12_Record_linkage/no_matches.dta") %>%
  filter(tocheck==1) %>% 
  rename(lastname_1=name_1) %>% 
  as_tibble() %>%
  select(id_0,firstname_0,name_0,time_period,male_0,e_cntr_w_0,n_cntr_w_0,gdename_0,
         firstname_1,lastname_1,id_1,year_1,firmnames_1) %>%
  rename(lastname_0=name_0)

data_rl_all <- haven::read_dta("./02_Processed_data/12_Record_linkage/rl_output_fn_check.dta") %>%
  rename(lastname_1=name_1) %>% 
  as_tibble() %>%
  select(id_0,firstname_0,name_0,time_period,male_0,e_cntr_w_0,n_cntr_w_0,gdename_0,
         firstname_1,lastname_1,id_1,year_1,firmnames_1) 

data3 <- openxlsx::read.xlsx("./02_Processed_data/12_Record_linkage/02_Sugarcube/10_RL_Output_Round1_Out/RL_Output_Round1_All.xlsx")


# (C) Functions


CutPolLinks <- function(data,pol_id){
# This function cuts the rl output data for a given politician id.    
return(data %>% filter(id_0==pol_id))  
}

GetUniqueSugarcube <- function(data){
# This function returns the unique sugarcube entries
  return(data %>% 
           distinct(firstname_0,lastname_0,firstname_1,lastname_1,e_cntr_w_0, n_cntr_w_0,gdename_0) %>% 
           filter(firstname_1!="" & lastname_1!=""))  
}


SugarSearch <- function(data1,data2,max_distance_input_fn=0,max_distance_input_ln=0){
# This function searches for entries in data1 using the firstname and lastname
# from data2. max_distance_input_fn and max_distance_input_ln are the allowed
# distances in terms of the firstname and lastname. 
  df_result <- data1 %>% 
    filter(agrepl(data2$firstname_1, firstname_1, 
                  max.distance=max_distance_input_fn) & 
             agrepl(data2$lastname_1, lastname_1, 
                    max.distance = max_distance_input_ln))
  return(df_result)  
}


SugarSearch_FalseNegatives <- function(pol_id,data1,data2,data3,
                                       max_distance_input_fn=0,
                                       max_distance_input_ln=0){
# This function looks for missed mathces in the Sugarcube data based on the 
# record linkage output for a specific politician. It proceeds as follows: 
# 1. It cuts all record linkage entries in a dataset (data2) for a 
# specific politician with an id (pol_id).
# 2. It gets the unique sugarcube entries in the record linkage data (across years)
# 3. It loops over all these entries and tries to find all entries in the
#    full Sugarcube dataset (data1). 
# 4. It labels previously found entries in the rl linkage output data (data3) as such and calculates
#    the distance in the residence municipality.
  
# data1=data_sug; data2=data_fn; pol_id="AG-1931-0002"; max_distance_input_fn=0;
#  max_distance_input_ln=0 

# 1. get unique politician list

data_pol <- CutPolLinks(data=data2,pol_id=pol_id)

# 2. get unique sugarcube entries that match to specific politician

data_pol_unique <- GetUniqueSugarcube(data=data_pol)


# 3. Search for further sugarcube entries

search_all <- data.frame()

for (i in 1:dim(data_pol_unique)[1]){

search_all <- bind_rows(search_all,SugarSearch(data1=data_sug,data2=data_pol_unique[i,],
                                               max_distance_input_fn=max_distance_input_fn,
                                               max_distance_input_ln=max_distance_input_ln))
search_all$e_cntr_w_0 <- data_pol_unique$e_cntr_w_0[i]
search_all$n_cntr_w_0 <- data_pol_unique$n_cntr_w_0[i]
search_all$gdename_0 <- data_pol_unique$gdename_0[i]
search_all$firstname_0 <- data_pol_unique$firstname_0[i]
search_all$lastname_0 <- data_pol_unique$lastname_0[i]
search_all$male_0 <- data_pol_unique$male_0[i]
}

# 4. Filter out only those that were not previously found and calculate distance

out <- search_all  %>%
  select(lastname_0,firstname_0,lastname_1,firstname_1,gdename_0, gdename_1,
         firmnames_1,year_1,id_1,e_cntr_w_0,n_cntr_w_0,e_cntr_w_1,n_cntr_w_1) %>%
  distinct(id_1,year_1,.keep_all = TRUE) %>%
  left_join(data3 %>% 
              select(id_0,id_1,year_1) %>%
              mutate(previously_found=1) %>%
              mutate(id_1=as.numeric(id_1)),
            by=c("id_1","year_1")) %>%
  mutate(distance=round(sqrt((n_cntr_w_0-n_cntr_w_1)^2+(e_cntr_w_0-e_cntr_w_1)^2)/1000,1)) %>%
  select(previously_found,year_1, id_1,starts_with("first"),starts_with("last"),starts_with("gde"),distance,everything()) %>%
  select(-starts_with("n_cntr"),-starts_with("e_cntr")) %>%
  arrange(lastname_0,firstname_0,year_1)

View(out %>% arrange(id_0,distance,year_1))

return(out)
}

# (D) Implement search

# (i) Simple search for one politician

test <- SugarSearch_FalseNegatives(pol_id="AG-1931-0002",
                                   data1=data_sug, 
                                   data2=data_fn, 
                                   max_distance_input_fn=0,
                                   max_distance_input_ln=0) 

View(test)

# (ii) Multiprocessing

ids_to_loop <- unique(data_fn$id_0)
no_cores <- availableCores() - 1
plan(multisession,workers=no_cores)

start_time <- Sys.time()
out_data <- future_map(ids_to_loop[1:10],
                       SugarSearch_FalseNegatives,
                       data1=data_sug, 
                       data2=data_fn, 
                       max_distance_input_fn=0,
                       max_distance_input_ln=0)
out_data_df <- as.data.frame(do.call(rbind, out_data))
end_time <- Sys.time()
end_time - start_time

View(out_data_df)

# computing time estimate on server: 6*length(ids_to_loop)/10*7/31*(1/60)


# (E) Implement search for specific examples (this is basis for email to Sandro, 3 August 2023)  ----


result <- data_sug %>% dplyr::filter(year_1 %in% c(1933,1942,1959) & firstname=="Roman" & lastname=="Abt")
View(result)

result <- data_sug %>% dplyr::filter(year_1 %in% c(1986) & firstname=="Josef" & lastname=="Ackermann")
View(result)


result <- data_sug %>% dplyr::filter(year_1 %in% c(1985) & firstname_1=="Josef" & lastname_1=="Ackermann")
View(result)

result <- data_sug %>% dplyr::filter(year_1 %in% c(1984) & firstname_1=="Josef" & lastname_1=="Ackermann")
View(result)

result <- data_sug %>% dplyr::filter(year_1 %in% c(1983) & firstname_1=="Josef" & lastname_1=="Ackermann")
View(result)

result <- data_sug %>% dplyr::filter(year_1 %in% c(1982) & firstname_1=="Josef" & lastname_1=="Ackermann")
View(result)

result <- data_sug %>% dplyr::filter(year_1 %in% c(1981) & firstname_1=="Josef" & lastname_1=="Ackermann")
View(result)

result <- data_sug %>% dplyr::filter(year_1 %in% c(1980) & firstname_1=="Josef" & lastname_1=="Ackermann")
View(result)

result <- data_sug %>% dplyr::filter(year_1 %in% c(1972:1986) & firstname_1=="Josef" & lastname_1=="Ackermann")
View(result)

result <- data_sug %>% dplyr::filter(year_sug %in% c(1981) & firstname_1=="Walter" & lastname_1=="Frei")
View(result)



# 
# #old code with specific commands
# # 2. get unique sugarcube entries that match to specific politician
# 
# data_pol_unique <- data_pol %>% 
#   distinct(firstname_1,name_1) %>% 
#   filter(firstname_1!="" & name_1!="")
# 
# # 3. get years with matches
# 
# data_pol_spec <- data_pol %>% 
#   right_join(data_pol_unique[2,],by=c("firstname_1","name_1")) %>%
#   select(firstname_1,name_1,year_1) 
# 
# years <- data_pol_spec %>% select(year_1) %>% pull()
# 
# # 4. Search for further sugarcube entries
# 
# 
# max_distance_input_fn <- 0
# max_distance_input_ln <- 0
# 
# 
# bigger_result <- data_sug %>% 
#   filter(agrepl(data_pol_unique$firstname_1[1], firstname, max.distance=max_distance_input_fn) & 
#            agrepl(data_pol_unique$name_1[1], lastname, max.distance = max_distance_input_ln))
# 
# View(bigger_result)
# 
# # (Z) Old Functions ----
# 
# CorrectName <- function(string_input){ # remove spacing after comma 
#   str_sub(string_input,2,-1)  
# }
# 
# StringReplace <- function(string_input){ # this is based on Mark Schelker's Email from 31 March 2020
#   return(str_replace_all(string_input,c("Ä"="ae","Ö"="oe","Ü"="ue","ä"="ae","ö"="oe","ü"="ue","È"="e","É"="e","Ê"="e","è"="e","é"="e","ê"="e","ë"="e","À"="a","Â"="a","à"="a","â"="a","Û"="u","Û"="u","û"="u","ù"="u","Ô"="o","Ô"="o","ô"="o","ò"="o","ó"="o" ,"Î"="i","Î"="i","î"="i","ì"="i","ï"="i","ç"="c","ß"="ss")))
# }
# 
# GetAllNames <- function(string_input){
#   splitted <- strsplit(string_input,split=",")
#   out <- list()
#   for (i in 1:length(splitted[[1]])){
#     if (i==1){out <- StringReplace(tolower(splitted[[1]])[1])}
#     else{out <- rbind(out,StringReplace(CorrectName(tolower(splitted[[1]])[i])))}
#   }
#   return(out)
# }
# 
# 
# 
# GetResult <- function(data_target, data_input,max_distance_input=0.2){
#   Names <- GetAllNames(data_input$Name)  
#   Firstnames <- GetAllNames(data_input$Vorname)  
#   GeoCodes <- GetAllNames(data_input$GdeNr_CNTR)
#   Municipality <- GetAllNames(data_input$Wohngemeinde)
#   out <- data.frame()
#   for (i in 1:length(Names)){
#     for (j in 1:length(Firstnames)){
#       for (k in 1:length(GeoCodes)){
#         GeoCodesSplit <- as.numeric(unlist(strsplit(as.character(GeoCodes[k]),split = "/")))
#         bigger_result <- data_target %>% filter(agrepl(Names[i], lastname, max.distance=max_distance_input) & agrepl(Firstnames[j], firstname, max.distance = 0.2))
#         bigger_result <- cbind(bigger_result,do.call(rbind, strsplit(as.character(bigger_result$GdeNr_CNTR), split = "/", fixed = FALSE)))
#         if (dim(bigger_result)[1]>0){
#           bigger_result <- bigger_result %>% mutate(Nachname_Pol=StringReplace(tolower(data_input$Name)),
#                                                     Vorname_Pol=StringReplace(tolower(data_input$Vorname)),
#                                                     Jahr_VR=as.numeric(as.character(year)),
#                                                     Geo_Code_Pol=GeoCodes[k],
#                                                     ID_Pol=data_input$ID,
#                                                     Alter=year-as.numeric(data_input$Geburtsjahr),
#                                                     Distanz=sqrt((as.numeric(as.character(`1`))-GeoCodesSplit[1])^2+(as.numeric(as.character(`2`))-GeoCodesSplit[2])^2)/1000,
#                                                     Wohngemeinde_Pol=Municipality[k],
#                                                     Beruf_Pol=data_input$Beruf) %>%
#             rename(Vorname_VR=firstname,
#                    Nachname_VR=lastname,
#                    ID_VR=PID,
#                    Wohngemeinde_VR=gdename,
#                    Firmen_VR=companies,
#                    Geo_Code_VR=GdeNr_CNTR)%>%
#             select(ID_Pol,ID_VR,
#                    Vorname_Pol,Nachname_Pol,
#                    Jahr_VR,Vorname_VR,Nachname_VR,
#                    Alter,Distanz,
#                    Wohngemeinde_Pol,Wohngemeinde_VR,
#                    Beruf_Pol,Firmen_VR,
#                    Geo_Code_Pol,Geo_Code_VR)%>%
#             arrange(ID_Pol,Vorname_VR,Nachname_VR,Wohngemeinde_Pol,Wohngemeinde_VR,Jahr_VR)
#           
#           
#           out <- rbind(out,bigger_result)
#         }
#         else{}
#       } 
#     } 
#   } 
#   if (dim(bigger_result)[1]>0){return(out %>% distinct(ID_Pol,ID_VR,Jahr_VR,Wohngemeinde_Pol,.keep_all = TRUE) )}
#   else{return(out)}
# }
# 
# GetResults <- function(data_target,data_input,max_distance_input=0.2){
#   out <- data.frame()
#   for (i in 1:dim(data_input)[1]){
#     print(paste("Iteration no.",i))
#     result <- GetResult(data_target,data_input[i,],max_distance_input) 
#     if (dim(result)[2]>0){out <- rbind(out,result)}
#   }
#   out <- out %>% group_by(ID_Pol) %>% 
#     mutate(CaseID = group_indices())
#   return(out)
# }
# 
# 
# GetResultsAllFiles <- function(data_target,path_input_files,
#                                path_output_files,max_distance_input){
#   file_names <- dir(path_input_files,pattern ="Sugarcube*")
#   print(paste("Files to generate: ",file_names,sep=""))
#   for (i in 1:length(file_names)){
#     print(paste("File name: ",as.character(file_names[i])))
#     my_data <- read_excel(paste(path_input_files,"/",file_names[i],sep="")) %>% filter(is.na(Geschlecht)==F)
#     Resultout <- GetResults(data_target=sugarcube,data_input=my_data,max_distance_input=max_distance_input)
#     write.table(Resultout,
#                 file=paste(path_output_files,"/Sugar_Checks",gsub("[^0-9]", "",file_names[i]),".csv",sep=""),
#                 sep=";")
#   }  
# }
# 
