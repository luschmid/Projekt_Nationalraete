
library(readstata13)
library(data.table)
library(fastLink)
library(tidyverse)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))
try(setwd("C:/Dropbox/Dropbox/Projekt Nationalräte")) 

rm(list = ls())

# (i) Read-in data

data_hn <- read.dta13("./02_Processed_data/15_Elections_Honduras/elections_hn.dta") %>% 
         as_tibble() %>% 
         group_by(Year,Department) %>% 
         mutate(rownumber=row_number()) %>%
         ungroup()
data_hn <- data_hn %>% mutate(ID=paste(substr(data_hn$Department,1,3),Year,rownumber,sep="-"))

data_hn$ID

# (ii) Functions (Unsupervised Classification per department)

DoFastLink <- function(data_input1,data_input2,varinput,cut.p.chosen){
  matches.out <- fastLink(
    dfA = data_input1,data_input2, 
    varnames = paste(varinput),
    stringdist.match = paste(varinput),
    partial.match = paste(varinput),
    cut.p=cut.p.chosen,
    n.cores=4
  )
  return(matches.out)
}

DoFastLinkLoop <- function(data,departmentname,years,varinput,cut.p.chosen){
  departments <- data %>% dplyr::select(!!departmentname) %>% unique() %>% pull()
  correspondenceall <- tibble()
  for (t in 1:(length(years)-1)){
    for (u in (t+1):(length(years))){
    for (i in 1:(length(departments))){
      print(paste("Year: ",years[t]))
      print(paste("Department: ",departments[i]))
      data1 <- data %>% filter(Department==as.character(departments[i])& Year==years[t])
      data2 <- data %>% filter(Department==as.character(departments[i])& Year==years[u])
      MatchOut <- DoFastLink(data_input1=data1,data_input2=data2,
                             varinput=varinput,
                             cut.p.chosen=cut.p.chosen) 
      
      correspondence<- tibble(ID1=data1$ID[MatchOut$matches$inds.a], 
                              ID2=data2$ID[MatchOut$matches$inds.b],
                              prob=MatchOut$posterior)
      
      if (i==1 & t==1 & u==2){
      write.table(correspondence,
                  file="./02_Processed_data/15_Elections_Honduras/MatchesHonduras.csv",
                  append=FALSE,sep=",",col.names=TRUE,row.names=TRUE)
      }
      else{
      write.table(correspondence,
                    file="./02_Processed_data/15_Elections_Honduras/MatchesHonduras.csv",
                  append=TRUE,sep=",",col.names=TRUE,row.names=TRUE)       
      }
    }
}
}
}

DisplayMatches <- function(data,ID1,ID2,varinput){
print(data %>% filter(ID==ID1) %>% select(!!varinput) %>% as.data.frame())
print(data %>% filter(ID==ID2) %>% select(!!varinput)%>% as.data.frame())
}

FindMatches <- function(IDinput,data,varinput){
  return(data %>% filter(ID==IDinput) %>% select(!!varinput) )
}

MatchesAll <- function(datainput,ID1,ID2,prob,varinput){
  data1 <- lapply(ID1,FindMatches,
                data=datainput,varinput=varinput)
  data1 <- data.frame(matrix(unlist(data1), nrow=length(data1), byrow=T))
  colnames(data1) <- paste(varinput,"1",sep="")
  data2 <- lapply(ID2,FindMatches,
                data=datainput,varinput=varinput)
  data2 <- data.frame(matrix(unlist(data2), nrow=length(data2), byrow=T))
  colnames(data2) <- paste(varinput,"2",sep="") 
  return(cbind(data1,data2,prob))
}

ReadInCsvFiles <- function(arg){
return(readr::read_csv2(paste("./02_Processed_data/15_Elections_Honduras/",arg,sep="")))
}

AppendFiles <- function(){
files <- list.files(path="./02_Processed_data/15_Elections_Honduras",pattern = "Matches*")
out <- lapply(files,ReadInCsvFiles)
return(do.call("rbind", out))
}


FindConnections <- function(data,ID1,ID2,ID){
  outfile <- data %>% dplyr::filter(ID1 %in% as.character(ID))  
  return(outfile$ID2)
}

FindConnectionsAll <- function(data,ID1,ID2,ID){
  # This function finds all possible connections for a certain ID in a corresponce table of IDs
  # Arguments: data: correspondence table (tibble or dataframe) 
  #            ID1: Name of first ID variable in data
  #            ID2: Name of second ID variable in data
  #            ID: ID(s) to be searched (string or vector)
  out <- ID # start with ID requested
  check=0
  while (check==0 ){
    outnew <- sort(levels(as.factor(c(ID,FindConnections(data=data,ID1="ID1", ID2="ID2",ID=c(out))))))
    if (identical(outnew,out)==F){
      out <- sort(levels(as.factor(c(out,outnew))))
    }
    else{check = 1}
  }
  return(out[out!=ID])
}

OrderDataFrame <- function(data,ID1,ID2){
  # This function orders the correspondence table by all possible connections 
  # Arguments: data: correspondence table (tibble or dataframe) 
  #            ID1: Name of first ID variable in data
  #            ID2: Name of second ID variable in data  
  data_ord <- tibble()
  IDs_all<- levels(as.factor(eval(data$ID1)))
  for (i in 1:length(IDs_all)){
    Connections <- FindConnections(data=data,
                                   ID1="ID1",
                                   ID2="ID2",
                                   ID=IDs_all[i])
    if (length(Connections)>=1){
      data_new <- data %>% filter(ID1 %in% Connections| ID1 %in% IDs_all[i]) %>% mutate(group=i)
      data_ord <- rbind(data_ord,data_new)
      data <- data %>% filter(!ID1 %in% Connections & !ID1%in%IDs_all[i])
    }
  } 
  return(data_ord)
}

LookUpID <- function(datainput,ID1input,ID2input){
# This function generates a matching ID (IDMatching) for one observation
# on a dataset (datainput) using IDs (ID1input,ID2input) from a correspondence table
  dataoutput <- datainput %>% 
    filter(ID==ID1input |ID==ID2input) %>% 
    mutate(IDMatching=ifelse(substr(ID1input,5,8)>substr(ID2input,5,8),ID2input,ID1input)) 
  return(dataoutput)  
}

LookUpIDs <- function(dataall,datamatching){
# This function generates a matching ID (IDMatching) for all observations 
  # (including those with not match)
  # Arguments: dataall: original dataset with an original ID
  #            datamatching : correspondence table
  outall <-tibble()
  for (i in 1:dim(datamatching)[1]){
    out <- LookUpID(datainput=dataall,ID1input=datamatching$ID1[i],
                    ID2input=datamatching$ID2[i])
    if (dim(out)[1]>=1){
    out$group <- datamatching$group[i]
    outall <- rbind(outall,out)
    dataall <- dataall %>% filter(ID!=datamatching$ID1[i]& ID!=datamatching$ID2[i])
    }
    else{}
  }
  dataall$IDMatching <- dataall$ID
  dataall$group <- max(outall$group)+as.numeric(as.factor(dataall$ID))
  return(rbind(outall,dataall))
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 



## (iii) Run program
# 
# DoFastLinkLoop(data=data_hn,
#                departmentname="Department",
#                years=c(2009,2013,2017),
#                varinput=c("Candidate","Year"),
#                cut.p.chosen=0.02)
# 
# MatchOutputHonduras %>% arrange(-prob)
# 
# MatchesHonduras <- MatchesAll(data=data_hn,
#            ID1=MatchOutputHonduras$ID1,
#            ID2=MatchOutputHonduras$ID2,
#            prob=MatchOutputHonduras$prob,
#            varinput=c("ID","Candidate","Party"))
# 
# write_delim(MatchesHonduras,
#             delim = ";", 
#             path="/MatchesHonduras.csv")


HondurasAllOut <-AppendFiles() %>% 
                 filter(prob>0.5) %>% 
                 select(ID1,ID2) %>%
                 arrange(ID1,ID2)
                
HondurasAllOutOrdered <- OrderDataFrame(HondurasAllOut,"ID1","ID2")
HondurasFinal <- LookUpIDs(dataall=data_hn,datamatching = HondurasAllOutOrdered)

HondurasFinal %>% select(Candidate,Year,Department,ID,IDMatching) %>%
                  filter(IDMatching =="CHO-2009-12")

write.table(HondurasFinal,
            file="./02_Processed_data/15_Elections_Honduras/MatchingHondurasAll.csv",
            append=FALSE,sep=",",col.names=TRUE,row.names=FALSE)

HondurasFinal %>% filter(Candidate=="MARCIO RENE ESPINAL CARDONA")












# test basic function ######################

data1 <- data %>% filter(Department=="FRANCISCO MORAZAN" & Year==2009)
data2 <- data %>% filter(Department=="ATLANTIDA" & Year==2013)
data_all <- data1 %>% bind_rows(data2) %>% arrange(Candidate) %>% print(n=112)
MatchOut <- DoFastLink(data_input1=data1,data_input2=data2,
                       vars=c("Candidate","Year"),
                       cut.p.chosen=0.15)
matchout1 <- data1 %>% filter(row_number()%in% MatchOut$matches$inds.a)
matchout2 <- data2 %>% filter(row_number()%in% MatchOut$matches$inds.b)
matchout1 %>% bind_rows(matchout2) %>% arrange(Candidate) 

matchout1 <- data1 %>% select(ID) %>%
             filter(row_number()%in% MatchOut$matches$inds.a)
matchout2 <- data2 %>%  select(ID) %>%
             filter(row_number()%in% MatchOut$matches$inds.b)

departmentname <- "Department"
cut.p.chosen <- 0.02
years <- c(2009,2013)
i=1;t=1;u=2;
data <- data_hn
varinput <- c("Candidate","Year")

DisplayMatches(data=data,
               ID1=MatchOutput$ID1[2],
               ID2=MatchOutput$ID2[2],
               vars=c("Candidate","Party")
)

FindMatches(data=data,
            IDinput=MatchOutputHonduras$ID1[2],
            vars=c("ID","Candidate","Party"))

LookUpID(datainput=data_hn,ID1input=HondurasAllOut$ID1[2],
         ID2input=HondurasAllOut$ID2[2])



test <- HondurasAllOut %>% filter(ID1=="ATL-2009-33" | ID1=="ATL-2013-40" | ID2=="ATL-2013-40" )
test <- rbind(test,tibble(ID1="ATL-2017-1",ID2="ATL-2020-40"))
test <- rbind(test,tibble(ID1="ATL-2009-33",ID2="ATL-2024-5"))

FindConnectionsAll(data=test,
                   ID1="ID1",
                   ID2="ID2",
                   ID=c("ATL-2009-33"))

FindConnections(data=HondurasAllOut,
                ID1="ID1",
                ID2="ID2",
                ID=c("ATL-2009-33"))
