library(readxl)
library(tidyverse)
library(data.table)

## (A) Working directory and load data

path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte"
path_source_files <- "/02_Processed_data/12_Record_linkage/02_Sugarcube/01_Nationalrat_Random_Sample"
path_out_files <- "/02_Processed_data/12_Record_linkage/02_Sugarcube/02_Nationalrat_VR_Initial_Match"


setwd(path) 

my_data <- read_excel(paste(path,"/02_Processed_data/12_Record_linkage/01_Check_False_Negatives_out/Sugarcube_Excel_version/Check_False_Negatives_Sugarcube1.xlsx",sep="")) %>% filter(is.na(Geschlecht)==F)
sugarcube <- fread("./02_Processed_data/10_Directors_1934_2003/Sugarcube_Person-Geo.csv", sep = ",", encoding = "UTF-8") %>%
             select(PID,year,firstname,lastname,gdename,GdeNr_CNTR)
sugarcube_companies <- fread("./02_Processed_data/10_Directors_1934_2003/Sugarcube_Companies.csv", sep = ";", encoding = "UTF-8") 
sugarcube <- sugarcube %>% left_join(sugarcube_companies,by=c("PID","year"))

## (B) Functions

CorrectName <- function(string_input){ # remove spacing after comma 
str_sub(string_input,2,-1)  
}

StringReplace <- function(string_input){ # this is based on Mark Schelker's Email from 31 March 2020
return(str_replace_all(string_input,c("Ä"="ae","Ö"="oe","Ü"="ue","ä"="ae","ö"="oe","ü"="ue","È"="e","É"="e","Ê"="e","è"="e","é"="e","ê"="e","ë"="e","À"="a","Â"="a","à"="a","â"="a","Û"="u","Û"="u","û"="u","ù"="u","Ô"="o","Ô"="o","ô"="o","ò"="o","ó"="o" ,"Î"="i","Î"="i","î"="i","ì"="i","ï"="i","ç"="c","ß"="ss")))
}

GetAllNames <- function(string_input){
splitted <- strsplit(string_input,split=",")
out <- list()
for (i in 1:length(splitted[[1]])){
  if (i==1){out <- StringReplace(tolower(splitted[[1]])[1])}
  else{out <- rbind(out,StringReplace(CorrectName(tolower(splitted[[1]])[i])))}
}
return(out)
}



GetResult <- function(data_target, data_input,max_distance_input=0.2){
  Names <- GetAllNames(data_input$Name)  
  Firstnames <- GetAllNames(data_input$Vorname)  
  GeoCodes <- GetAllNames(data_input$GdeNr_CNTR)
  Municipality <- GetAllNames(data_input$Wohngemeinde)
  out <- data.frame()
  for (i in 1:length(Names)){
    for (j in 1:length(Firstnames)){
      for (k in 1:length(GeoCodes)){
        GeoCodesSplit <- as.numeric(unlist(strsplit(as.character(GeoCodes[k]),split = "/")))
        bigger_result <- data_target %>% filter(agrepl(Names[i], lastname, max.distance=max_distance_input) & agrepl(Firstnames[j], firstname, max.distance = 0.2))
        bigger_result <- cbind(bigger_result,do.call(rbind, strsplit(as.character(bigger_result$GdeNr_CNTR), split = "/", fixed = FALSE)))
        if (dim(bigger_result)[1]>0){
        bigger_result <- bigger_result %>% mutate(Nachname_Pol=StringReplace(tolower(data_input$Name)),
                                                  Vorname_Pol=StringReplace(tolower(data_input$Vorname)),
                                                  Jahr_VR=as.numeric(as.character(year)),
                                                  Geo_Code_Pol=GeoCodes[k],
                                                  ID_Pol=data_input$ID,
                                                  Alter=year-as.numeric(data_input$Geburtsjahr),
                                                  Distanz=sqrt((as.numeric(as.character(`1`))-GeoCodesSplit[1])^2+(as.numeric(as.character(`2`))-GeoCodesSplit[2])^2)/1000,
                                                  Wohngemeinde_Pol=Municipality[k],
                                                  Beruf_Pol=data_input$Beruf) %>%
                                          rename(Vorname_VR=firstname,
                                                 Nachname_VR=lastname,
                                                 ID_VR=PID,
                                                 Wohngemeinde_VR=gdename,
                                                 Firmen_VR=companies,
                                                 Geo_Code_VR=GdeNr_CNTR)%>%
                                          select(ID_Pol,ID_VR,
                                                 Vorname_Pol,Nachname_Pol,
                                                 Jahr_VR,Vorname_VR,Nachname_VR,
                                                 Alter,Distanz,
                                                 Wohngemeinde_Pol,Wohngemeinde_VR,
                                                 Beruf_Pol,Firmen_VR,
                                                 Geo_Code_Pol,Geo_Code_VR)%>%
                                         arrange(ID_Pol,Vorname_VR,Nachname_VR,Wohngemeinde_Pol,Wohngemeinde_VR,Jahr_VR)
        
        
      out <- rbind(out,bigger_result)
        }
      else{}
      } 
    } 
  } 
  if (dim(bigger_result)[1]>0){return(out %>% distinct(ID_Pol,ID_VR,Jahr_VR,Wohngemeinde_Pol,.keep_all = TRUE) )}
  else{return(out)}
}

GetResults <- function(data_target,data_input,max_distance_input=0.2){
out <- data.frame()
  for (i in 1:dim(data_input)[1]){
  print(paste("Iteration no.",i))
  result <- GetResult(data_target,data_input[i,],max_distance_input) 
  if (dim(result)[2]>0){out <- rbind(out,result)}
  }
out <- out %>% group_by(ID_Pol) %>% 
  mutate(CaseID = group_indices())
return(out)
}


GetResultsAllFiles <- function(data_target,path_input_files,
                               path_output_files,max_distance_input){
  file_names <- dir(path_input_files,pattern ="Sugarcube*")
  print(paste("Files to generate: ",file_names,sep=""))
  for (i in 1:length(file_names)){
    print(paste("File name: ",as.character(file_names[i])))
    my_data <- read_excel(paste(path_input_files,"/",file_names[i],sep="")) %>% filter(is.na(Geschlecht)==F)
    Resultout <- GetResults(data_target=sugarcube,data_input=my_data,max_distance_input=max_distance_input)
    write.table(Resultout,
                file=paste(path_output_files,"/Sugar_Checks",gsub("[^0-9]", "",file_names[i]),".csv",sep=""),
                sep=";")
  }  
}




## (C) Check functions

# (a) Transformation and correction of names

# out <- GetAllNames(my_data$Name[67])
# GetAllNames(my_data$GdeNr_CNTR[7])
# CorrectName(my_data$GdeNr_CNTR[7])
# 
# # (b) GetResult (one single candidate)
# tester_0.5 <- GetResult(data_target=sugarcube,data_input=my_data[1,],max_distance_input=0.5)
# tester_0.1 <- GetResult(data_target=sugarcube,data_input=my_data[1,],max_distance_input=0.1)
# tester_0.05 <- GetResult(data_target=sugarcube,data_input=my_data[1,],max_distance_input=0.05)
# test <- GetResult(data_target=sugarcube,data_input=my_data[24,],max_distance_input=0.1)
# # 881 obs
# test  %>% print(n=800)
# test  %>% distinct(Wohngemeinde_Pol)


# (c) GetResults (multiple candidates): 
#     Note: We explore different distances of name similarity (max_distance_input). 

#test_0.5 <- GetResults(data_target=sugarcube,data_input=my_data,max_distance_input=0.5)
#test_0.2 <- GetResults(data_target=sugarcube,data_input=my_data,max_distance_input=0.2)
#test_0.1 <- GetResults(data_target=sugarcube,data_input=my_data,max_distance_input=0.1)
#test_0.05 <- GetResults(data_target=sugarcube,data_input=my_data,max_distance_input=0.05)
#test_0.01 <- GetResults(data_target=sugarcube,data_input=my_data,max_distance_input=0.01)

# test_0.2 <- test_0.2 %>% mutate(test_0_2=1)
# test_0.1 <- test_0.1 %>% mutate(test_0_1=1) %>% select(ID_VR,Jahr_VR,test_0_1)
#test_0.05 <-  test_0.05 %>% mutate(test_0_05=1)%>% select(ID_VR,Jahr_VR,test_0_05)
#test_0.01 <- test_0.01 %>% mutate(test_0_01=1)%>% select(ID_VR,Jahr_VR,test_0_01)

# test_all <- test_0.2 %>% 
#             full_join(test_0.1,by=c("ID_VR","Jahr_VR"))%>% 
#             select(test_0_2,test_0_1,everything())
# 
# test_all %>% filter(test_0_2==1 & is.na(test_0_1)==T) %>% 
#              select(-test_0_2,-test_0_1 ,-ID_Pol ,-ID_VR,-Unsicher_Match,-CaseID)%>%
#              print(n=100)
# 
# write.table(test_all,paste(path,path_source_files,"/Sugar_Checks_1_test_all.csv",sep=""),sep=";")

# (d) Outfile for Laura Decet to test

#lauradecet_0.1 <- GetResults(data_target=sugarcube,data_input=my_data,max_distance_input=0.1)
#write.table(lauradecet_0.1,paste(path,path_source_files,"/Sugar_Checks_11.csv",sep=""),sep=";")

# (e) Generate all data files

GetResultsAllFiles(data_target=sugarcube,
                   path_input_files=paste(path,path_source_files,sep=""),
                   path_output_files=paste(path,path_out_files,sep=""),
                   max_distance_input=0.1)




######## old search method

s_nachname = "luescher"
s_vorname = "gottlieb"
small_result = sugarcube %>% filter(grepl(s_nachname,lastname) & grepl(s_vorname, firstname))
bigger_result = sugarcube %>% filter(agrepl(s_nachname, lastname, max.distance=0.2) & agrepl(s_vorname, firstname, max.distance = 0.2))
View(small_result)
View(bigger_result %>% dplyr::filter(!PID %in% small_result$PID))