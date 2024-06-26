
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

# (B) Functions ----

SugarSearch <- function(data,
                        searchterm_fn,
                        searchterm_ln,
                        max_distance_input_fn=0,
                        max_distance_input_ln=0,
                        max_distance_input_ln2=0,
                        year){
  # This function searches for entries in data using the firstname and lastname
  # specified in searchterm_fn and searchterm_ln. max_distance_input_fn and 
  # max_distance_input_ln are the allowed distances in terms of the firstname and 
  # lastname. All entries in the specified year will not be shown. 
  
  searchterm_fn1 <- str_split_i(searchterm_fn, "-", 1)
  searchterm_fn2 <- str_split_i(searchterm_fn, " ", 1)

  searchterm_ln1 <- str_split_i(searchterm_ln, "-", 1)
  searchterm_ln2 <- str_split_i(searchterm_ln, "-", 2)
  searchterm_ln3 <- str_split_i(searchterm_ln, " ", 1)
  searchterm_ln4 <- str_split_i(searchterm_ln, " ", 2)  
  
  return(bind_rows(data %>% 
           filter(agrepl(searchterm_fn1, firstname_stub, 
                         max.distance=max_distance_input_fn) & 
                    agrepl(searchterm_ln1, lastname_stub, 
                           max.distance = max_distance_input_ln) & 
                    year_sug!=year) %>%
    arrange(lastname,firstname), 
    data %>% 
      filter(agrepl(searchterm_fn1, firstname_stub, 
                    max.distance=max_distance_input_fn) & 
               agrepl(searchterm_ln2, lastname_stub, 
                      max.distance = max_distance_input_ln) & 
               year_sug!=year) %>%
      arrange(lastname,firstname), 
    data %>% 
      filter(agrepl(searchterm_fn2, firstname_stub, 
                    max.distance=max_distance_input_fn) & 
               agrepl(searchterm_ln1, lastname_stub, 
                      max.distance = max_distance_input_ln) & 
               year_sug!=year) %>%
      arrange(lastname,firstname), 
    data %>% 
      filter(agrepl(searchterm_fn2, firstname_stub, 
                    max.distance=max_distance_input_fn) & 
               agrepl(searchterm_ln2, lastname_stub, 
                      max.distance = max_distance_input_ln) & 
               year_sug!=year) %>%
      arrange(lastname,firstname), 
    data %>% 
      filter(agrepl(searchterm_fn1, firstname_stub, 
                    max.distance=max_distance_input_fn) & 
               agrepl(searchterm_ln1, lastname_split2_stub, 
                      max.distance = max_distance_input_ln2) & 
               year_sug!=year) %>%
      arrange(lastname,firstname), 
    data %>% 
      filter(agrepl(searchterm_fn2, firstname_stub, 
                    max.distance=max_distance_input_fn) & 
               agrepl(searchterm_ln1, lastname_split2_stub, 
                      max.distance = max_distance_input_ln2) & 
               year_sug!=year)) %>%
      arrange(lastname,firmnames) %>%
      distinct(year_sug,id_sug,.keep_all = TRUE) %>%
      arrange(lastname,firstname,year_sug,gdename)
  )
    
}



SugarSearchDF <- function(id,
                          data1,
                          data2,
                          max_distance_input_fn=0,
                          max_distance_input_ln=0,
                          max_distance_input_ln2=0,
                          year){
data1$target <- 1
  
return(bind_rows(data1[id,], SugarSearch(data=data2,
                   searchterm_ln=data1$lastname_stub[id],
                   searchterm_fn=data1$firstname_stub[id],
                   max_distance_input_fn=max_distance_input_fn,
                   max_distance_input_ln=max_distance_input_ln,
                   max_distance_input_ln2=max_distance_input_ln2,
                   year=data1$year[id]))%>%
         arrange(lastname,firstname,year_sug,gdename))
  
}

# (C) Implement search ----

# (i) Read in data ----

data_sug <- fread("./02_Processed_data/10_Directors_1934_2003/12_Panel_Round1/Sugarcube_Panel_start.csv",
                   sep = ",", encoding = "UTF-8") 

# (ii) Get sample of sugarcube data ----

set.seed(1234)
random_sample_size <- 1000
data_sug_sample <- data_sug[sample(nrow(data_sug),random_sample_size),] %>%
  mutate(id=row_number())


# (iii) Multiprocessing

ids_to_loop <- c(1:length(data_sug_sample$lastname)) 
no_cores <- availableCores() - 1
plan(multisession,workers=no_cores)

start_time <- Sys.time()
out_data <- future_map(ids_to_loop[1:8],
                       SugarSearchDF,
                       data1=data_sug_sample, 
                       data2=data_sug, 
                       max_distance_input_fn=0.1,
                       max_distance_input_ln=0.1,
                       max_distance_input_ln2=0.1)
out_data_df <- as.data.frame(do.call(rbind, out_data)) %>%
  select(target, everything()) 
end_time <- Sys.time()
end_time - start_time

View(out_data_df)


data_test <- fwrite(out_data_df,"./02_Processed_data/10_Directors_1934_2003/12_Panel_Round1/test_search.csv",
                   sep = ";") 

# (D) Test function

data_test <- fread("./02_Processed_data/10_Directors_1934_2003/12_Panel_Round1/test_searchfunction.csv",
                   sep = ";", encoding = "UTF-8") %>%
  mutate(firstname_stub=substr(firstname,1,15))%>%
  mutate(lastname_stub=substr(lastname,1,15))%>%
  mutate(lastname_split2_stub=substr(lastname_split2,1,15)) 


start_time <- Sys.time()
out <- SugarSearch(data=data_sug,
            searchterm_ln=substr(data_sug_sample$lastname_stub[1],1,15),
            searchterm_fn=substr(data_sug_sample$firstname_stub[1],1,15),
            max_distance_input_fn=0.05,
            max_distance_input_ln=0.05,
            max_distance_input_ln2=0.05,
            year=data_sug_sample$year_sug[1])
end_time <- Sys.time()
end_time - start_time




SugarSearch(data=data_test,
            searchterm_ln=substr("Abaecherli-Hu",1,15),
            searchterm_fn=substr("Alois",1,15),
            max_distance_input_fn=0.05,
            max_distance_input_ln=0.05,
            max_distance_input_ln2=0.05,
            year=1972)


SugarSearch(data=data_test,
            searchterm_ln="Alai de Grolimund",
            searchterm_fn="Christa",
            max_distance_input_fn=0.1,
            max_distance_input_ln=0.1,
            year=1972)


SugarSearch(data=data_sug,
            searchterm_ln="Schnyder",
            searchterm_fn="Isabelle Clara",
            max_distance_input_fn=0.2,
            max_distance_input_ln=0.2,
            year=1972)


out <-SugarSearch(data=data_sug,
            searchterm_ln="Schnyder",
            searchterm_fn="Isabelle Clara",
            max_distance_input_fn=0.4,
            max_distance_input_ln=0.4,
            year=1972)

View(out)

out <- SugarSearch(data=data_sug,
            searchterm_ln="Schny",
            searchterm_fn="Isabe",
            max_distance_input_fn=0.05,
            max_distance_input_ln=0.05,
            year=1972)
  
View(out)

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
