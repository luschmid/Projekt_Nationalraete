#-----------------------------------------
# (A) Libraries, paths, and load functions
#-----------------------------------------


library(tidyverse)
library(data.table)
library(testit)
library(RSQLite)
library(odbc)
library(dbplyr)
library(dplyr, warn.conflict = FALSE, quietly = TRUE)
library(lubridate)read

rm(list=setdiff(ls(all.names=TRUE), lsf.str(all.names=TRUE)))
setwd("C:/Schmidlu/Dropbox/Record Linkage")
source("./01_Code/03_Record_Linkage_Results/00_Functions.R")
db_path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte/02_Processed_data/database_political_rents.sqlite"

#----------------------
# B) Specific functions
#----------------------

GetCorrespondenceTable <- function(modelname,source="bisnode",linkscore_cutoff="optimal"){
# This function returns the correspondence table of politician ids (id_0) and 
# director ids (id_1) for a specific model. The modelname can be specified as
# "Generation_1" and is then retransformed into Sandro's original name
  

# 1. Choose optimal link score

if (source=="bisnode"){
all_models <- data.table::fread("./03_Results/Bisnode-ModelResults.csv",
                                  sep = ";", encoding = "UTF-8"
)}
if (source=="sugarcube"){
    all_models <- data.table::fread("./03_Results/Sugarcube-ModelResults.csv",
                                    sep = ";", encoding = "UTF-8"
)}
  
assert("Model specified in modelname not found.", 
         dim(all_models %>% dplyr::filter(model==modelname))[1]>0)
  

if (linkscore_cutoff=="optimal"){
  linkscore_cutoff <- all_models %>% 
  filter(model==modelname) %>%
  select(linkscore_cutoff) %>%
  pull()
} 

# 2. Read in RL output data and return those connections with Link_Score higher than optimal linkscore cutoff


modelname_original <- ReplaceModelName(modelname=modelname,source=source,way="back")

path_relative <- ifelse(source=="bisnode","./02_Data/07_Output_RL_Bisnode/", "./02_Data/08_Output_RL_Sugarcube/")

rl_results1 <- ReadinResults(file_name = paste0("RL-Results-",modelname_original, ".csv"),
                             path_rel=path_relative,
                             source=source)

if (source=="bisnode"){
rl_results_wide1_all <- CreateWideDataset(
  data = rl_results1,
  source = source
)
return(rl_results_wide1_all %>% filter(Link_Score>=linkscore_cutoff) %>%
         select(id_0,id_1))
}

if (source=="sugarcube"){
rl_results_wide1_all <- CreateWideDataset(data=rl_results1,
                                          source="sugarcube",
                                          all_of = c(
                                            "name", "firstname", "sex", "year",
                                            "e_cntr_w", "n_cntr_w", "birthyear", "e_cntr_b",
                                            "n_cntr_b", "id"
                                          )) %>%
  select(id_0, id_1, year_1,Link_Score)

return(rl_results_wide1_all %>% filter(Link_Score>=linkscore_cutoff) %>%
         select(id_0,id_1,year_1))

}
}



BisnodeTransformToLong <- function(daf,id_underconsid){
# This function takes daf (director data) as an input and calculates for a given
# id (id_underconsid) the number range of years with at least one mandate. 
# The reference day is 1. January. For example, if someone has one mandate with 
# a starting date on 20 January 2000 and and end date on 30 March 2004, the function
# would return the years 2001,2002,2003, and 2004. 
  
daf_cut <- daf %>% filter(personenid==id_underconsid) %>%
  select(personenid,eintrittdatum,austrittdatum) %>%
  pivot_longer(cols=c(2:3)) %>%
  mutate(value_year=as.numeric(substr(value,1,4))) %>%
  mutate(value_day=as.numeric(substr(value,6,7))) %>%
  mutate(value_month=as.numeric(substr(value,9,10))) 

eintrittdatum_min <- daf_cut %>%
  filter(name=="eintrittdatum") %>% 
  summarize(value_year_min=min(value_year,na.rm=T)) %>%
  pull()

austrittdatum_max <- daf_cut %>%
  filter(name=="austrittdatum") %>% 
  summarize(value_year_max=max(value_year,na.rm=T)) %>%
  pull()

daf_cut_long <- data.frame(year=rep((eintrittdatum_min+1):austrittdatum_max),
                           id=id_underconsid)

# check if eintrittsdatum in first year is 1.1. (if yes, append this year)

min_eintrittdatum <- daf_cut %>% filter(value_year==eintrittdatum_min) %>%
  summarize(min_month=min(value_month,na.rm=T),
            min_day=min(value_day,na.rm=T))

if (min_eintrittdatum$min_day==1 & min_eintrittdatum$min_month==1){
daf_cut_long <- bind_rows(data.frame(year=eintrittdatum_min,
                                     id=id_underconsid),daf_cut_long)  
}
return(daf_cut_long)  
}


AddYearsWithNoMandates <- function(data,
                                   category_label="Verwaltungsrat",
                                   year_range){
  # This function adds years with no mandates to a given out_final dataframe for
  # a person
  # example: 
  # 1  1994 Verwaltungsrat        0
  # 2  1995 Verwaltungsrat       20
  # 3  1996 Verwaltungsrat       18
  # 4  1997 Verwaltungsrat       18
  # 5  1998 Verwaltungsrat       17
  # 6  1999 Verwaltungsrat        0
  
  if (dim(data)[1]>=1){
    year_missing <- year_range[!year_range %in% data$year]
  } else {
    year_missing <- year_range
  }
  
  if (length(year_missing)>0){
    
    out_final2 <- data.frame(year=year_missing,category=category_label,
                             mandates=0)
    
    out_final_all <- bind_rows(data,out_final2) %>%
      arrange(year)
  } else{
    out_final_all <- data
  }
  return(out_final_all)
}

MandatesYearly <- function(data,category_label,year_range=c(1994:2017)){
  
# This function aggregates bisnode input data (data) into the number of mandates
# per year
# Input example: 
# personenid eintrittdatum austrittdatum      duns        gremium eintrittjahr austrittjahr
# 1         760    1994-01-01    1994-07-11 482337987 Verwaltungsrat         1994         1994
# 2         760    1994-01-01    1995-08-30 487802795 Verwaltungsrat         1994         1995
# 3         760    1994-01-01    1997-03-07 480006188 Verwaltungsrat         1994         1996
# 4         760    1994-01-01    1997-03-12 480001759 Verwaltungsrat         1994         1996
# 
# Output example:
# 
# 1  1994 Verwaltungsrat       20
# 2  1995 Verwaltungsrat       20
# 3  1996 Verwaltungsrat       18
# 4  1997 Verwaltungsrat       18
# 5  1998 Verwaltungsrat       17
# 6  1999 Verwaltungsrat       10
# 7  2000 Verwaltungsrat       12
# 8  2001 Verwaltungsrat       11
# 9  2002 Verwaltungsrat       11
# 10  2003 Verwaltungsrat      12
  
out <- data.frame()

# 1. Generate outfile for those with no mandates


if (dim(data)[1]==0){
out_final <- data.frame(year=year_range,
                        category=category_label,
                        mandates=0) 

} else {
  
# 2. Reshape to Long dataframe and append for those with positive mandates
  
for (j in 1:dim(data)[1]){
  out <- bind_rows(out, BisnodeTransformToLong(daf=data[j,],
                                               id_underconsid=data$personenid[j]) %>%
                     mutate(duns=data$duns[j],
                            funktion=data$funktion[j],
                            gremium=data$gremium[j])) 
}

out_collapsed_vr <- out %>% distinct(year,id,duns,gremium)

  
  # 2. Aggregate at yearly level
  
  out_final <- out_collapsed_vr %>% 
    group_by(year,gremium) %>% 
    summarize(mandates=n()) %>%
    filter(year %in% year_range) %>% 
    mutate(category=category_label) %>%
    select(year,category,mandates)
}
  return(out_final)
  
}

QueryPoliticianMandates <- function(pol_id,
                                    data_input,
                                    correspondencetable,
                                    source="bisnode",
                                    year_range=c(1994:2017)){
  
# This function uses the correspondence table and the director data (data_input) 
# as an input and returns the number of connections for a certain year range 
# for the categories Verwaltungsrat or Verwaltungsrat_Praesident. 
# The source is either "bisnode" or "sugarcube. 
# test id: pol_id="ZH-1979-0051" (Blocher)
  
assert("Correspondence table has no entries.",dim(correspondencetable)[1]>0)
#assert("Director data has no entries.",dim(director_data)[1]>0)
  

# 1. Restrict correspondencetable to pol_id of interest
  
  
correspondencetable_id <- correspondencetable %>% 
  filter(id_0==pol_id) %>%
  distinct(id_0,id_1)


if (source=="sucarcube"){
  out_final1 <-correspondencetable_id %>% 
    mutate(year=as.numeric(as.character(year_1)))%>% 
    group_by(year) %>% 
    summarize(mandates=n()) %>%
    filter(year %in% year_range)
}

  
# 2. (Bisnode)  Create out file for those pol_ids with no connection

if (source=="bisnode"){
if (dim(correspondencetable_id)[1]==0){
  out_final <- data.frame(year=year_range,
                          Verwaltungsrat=0,
                          Verwaltungsrat_Praesident=0) %>% 
    pivot_longer(cols=2:3,names_to="category",
                 values_to="mandates") %>%
    mutate(id_0=pol_id)
}
else{


# 3. (Bisnode) Data manipulation for eintrittdatum and austrittdatum
  
# a) All mandates in Verwaltungsrat

director_data_cut <- data_input %>%
  filter(personenid %in% correspondencetable_id$id_1) %>%
  select(personenid,eintrittdatum,vorname, nachname,austrittdatum,duns,gremium,funktion) %>%
  filter(gremium%in% c("Verwaltungsrat","Bankrat")) %>%
  mutate(eintrittdatum=as.character(eintrittdatum)) %>%  
  mutate(austrittdatum=as.character(austrittdatum)) %>%  
  mutate(eintrittdatum=ifelse(is.na(eintrittdatum), "1994-01-01",eintrittdatum)) %>%
  mutate(austrittdatum=ifelse(is.na(austrittdatum), "2017-12-31",austrittdatum)) %>%
  mutate(eintrittdatum=ymd(eintrittdatum)) %>%  
  mutate(austrittdatum=ymd(austrittdatum)) %>% 
  filter(austrittdatum>"1994-01-01")%>% 
  mutate(eintrittdatum=ifelse(eintrittdatum<"1994-01-01",as.character("1994-01-01"),as.character(eintrittdatum))) %>% # first entry date: 1994-01-01
  mutate(eintrittdatum=ymd(eintrittdatum)) %>%  
  distinct(duns,personenid,gremium,eintrittdatum,austrittdatum) %>% # remove duplicates in variable "funktion" and multiple entries
  arrange(eintrittdatum,austrittdatum) %>%
  mutate(eintrittjahr=ifelse(month(eintrittdatum)>=6,year(eintrittdatum)+1,year(eintrittdatum))) %>% # Stichtag: 1. June of the respective year
  mutate(austrittjahr=ifelse(month(austrittdatum)>=6,year(austrittdatum),year(austrittdatum)-1)) %>% # dito
  filter(eintrittjahr<=austrittjahr) # remove those mandates that were started and ended within a year, where year starts on Stichtag defined above (see Blocher example with duns 486686525)


# b) All presidencies

director_data_cut_president <- data_input %>%
  filter(personenid %in% correspondencetable_id$id_1) %>%
  select(personenid,eintrittdatum,vorname, nachname,austrittdatum,duns,gremium,funktion) %>%
  filter(gremium%in% c("Verwaltungsrat","Bankrat")) %>%
  filter(funktion %in% c("Präsident/in","Präsident des Bankrates"))%>%
  mutate(eintrittdatum=as.character(eintrittdatum)) %>%  
  mutate(austrittdatum=as.character(austrittdatum)) %>%  
  mutate(eintrittdatum=ifelse(is.na(eintrittdatum), "1994-01-01",eintrittdatum)) %>%
  mutate(austrittdatum=ifelse(is.na(austrittdatum), "2017-12-31",austrittdatum)) %>%
  mutate(eintrittdatum=ymd(eintrittdatum)) %>%  
  mutate(austrittdatum=ymd(austrittdatum)) %>% 
  filter(austrittdatum>"1994-01-01")%>% 
  mutate(eintrittdatum=ifelse(eintrittdatum<"1994-01-01",as.character("1994-01-01"),as.character(eintrittdatum))) %>% # first entry date: 1994-01-01
  mutate(eintrittdatum=ymd(eintrittdatum)) %>%  
  distinct(duns,personenid,gremium,eintrittdatum,austrittdatum) %>% 
  arrange(eintrittdatum,austrittdatum)%>%
  mutate(eintrittjahr=ifelse(month(eintrittdatum)>=6,year(eintrittdatum)+1,year(eintrittdatum))) %>% 
  mutate(austrittjahr=ifelse(month(austrittdatum)>=6,year(austrittdatum),year(austrittdatum)-1)) %>%
  filter(eintrittjahr<=austrittjahr) 

  
# 4. (Bisnode) Count mandates per year (for those with positive number of mandates)

out_final <- MandatesYearly(director_data_cut,category_label="Verwaltungsrat")
out_final_president <- MandatesYearly(director_data_cut_president,category_label="Verwaltungsrat_Praesident")

# 5. Add years with no mandates

out_final <- AddYearsWithNoMandates(data=out_final,category_label="Verwaltungsrat",
                                    year_range=year_range)
out_final_president <- AddYearsWithNoMandates(data=out_final_president,
                                              category_label="Verwaltungsrat_Praesident",
                                              year_range=year_range)
out_final <- out_final %>% bind_rows(out_final_president) %>% mutate(id_0=pol_id)
}
}
return(out_final)
}



ReShapeBisnodeResults <- function(pol_id,nr_data,input_data,max_period=86){
# This function reshapes the input_data for a specific pol_id to an output
# data format that can be easily transformed in Stata and merged to our nr_data
# for the RDD estimations. The output data is at the year-pol_id observations
# and reports lags and leads for all categories of mandates (e.g., Verwaltungsrat
# and Verwaltungsrat_Praesident)
#
# example input data  
# 1 ZH-1979-0051  1994 Verwaltungsrat       20
# 2 ZH-1979-0051  1995 Verwaltungsrat       20
# 3 ZH-1979-0051  1996 Verwaltungsrat       18
# 4 ZH-1979-0051  1997 Verwaltungsrat       18
# 5 ZH-1979-0051  1998 Verwaltungsrat       17
# 6 ZH-1979-0051  1999 Verwaltungsrat       10
#
# example output data
# category                    year          ID          mandate_lag_86 mandate_lag_85 
# 1            Verwaltungsrat 1979 ZH-1979-0051             NA             NA           
# 2 Verwaltungsrat_Praesident 1979 ZH-1979-0051             NA             NA            
# 3            Verwaltungsrat 1983 ZH-1979-0051             NA             NA           
# 4 Verwaltungsrat_Praesident 1983 ZH-1979-0051             NA             NA             

  # 1. Get all years from NR data for a specific pol_id
  
  nr_years <- nr_data %>%
    filter(ID==pol_id) %>%
    select(ID,year)
  
  #min <- min(input_data$year)-min(nr_years$year)
  #max <- max(input_data$year)-min(nr_years$year)
  
  # 2. Reshape results and loop over all years and all mandate types
  
  results_reduced <- input_data %>% 
    select(year,everything()) 
  
  #levels_category <- colnames(results_reduced)[2:dim(results_reduced)[2]]
  levels_category <- levels(as.factor(results_reduced$category))
  
  out_df_final <- data.frame()
  
  for (yr in nr_years$year){
    out_df <- data.frame()
    #print(paste0("Year: ",yr))
    
    for (cat in levels_category){
      
      #print(paste0("Gremium/Funktion: ",gf))
      
      #results_reduced_temp <- results_reduced[,(j+1)]
      results_reduced_temp <- results_reduced %>%
        filter(category==cat) %>%
        select(mandates,year)

      out_df <- rbind(out_df,t(results_reduced_temp$mandates))
    }
    
    colnames(out_df) <- results_reduced_temp$year-yr
    out_df$category <- levels_category
    
    out_df$year <- yr
    out_df$ID <- pol_id
    out_df_final <- bind_rows(out_df_final,out_df)

  }
  
  # 3. Add variables for those years that are missing in previous query
  
  row.names(out_df_final) <- c(1:length(row.names(out_df_final)))
  
  out_df_final <- out_df_final%>%  select(category,year,ID,everything())
  
  years_all <- seq(-max_period,max_period,1) 
  years_covered <- sort(as.numeric(as.character(colnames(out_df_final[,4:dim(out_df_final)[2]]))))
  years_to_add <- setdiff(years_all,years_covered)
  
  cols_to_replace <- colnames(out_df_final[,4:dim(out_df_final)[2]])
  cols_to_replace_numeric <- as.numeric(cols_to_replace)
  
  cols_to_replace[cols_to_replace_numeric<0] <- paste0("_lag_",abs(cols_to_replace_numeric[cols_to_replace_numeric<0]))
  cols_to_replace[cols_to_replace_numeric>0] <- paste0("_lead_",abs(cols_to_replace_numeric[cols_to_replace_numeric>0]))
  cols_to_replace[cols_to_replace_numeric==0] <- paste0("_",abs(cols_to_replace_numeric[cols_to_replace_numeric==0]))
  
  colnames(out_df_final) <- c("category","year","ID",paste0("mandate",cols_to_replace))

  for (yr in years_to_add){
  if (yr<0){
  out_df_final[[paste0("mandate_lag_",abs(yr))]] <- NA
  }
  if (yr==0){
      out_df_final[[paste0("mandate_",yr)]] <- NA
  } 
  if (yr>0){
    out_df_final[[paste0("mandate_lead_",yr)]] <- NA
  }  
  }
  
  return(  
    out_df_final %>% select(category,year,ID,paste0("mandate_lag_",seq(max_period,1,-1)),
                            "mandate_0",
                            paste0("mandate_lead_",seq(1,max_period,1)))
    )
  
}

WriteModelResults <- function(modelname="Generation_7",
                              year_range=c(1994:2017),
                              source="bisnode",
                              batch_size=1000,
                              batch_number=1,
                              linkscore_cutoff="optimal",
                              path="C:/Schmidlu/Dropbox/Projekt Nationalräte/02_Processed_data/11_Directors_1994_2018"){
  
  
  # 1. Load correspondence table and save vector of id_0 and id_1 there in 
  
  CorrespondenceTable <- GetCorrespondenceTable(modelname=modelname,
                                                source=source,
                                                linkscore_cutoff=linkscore_cutoff) %>%
    arrange(id_0) %>% 
    distinct(id_0,id_1)
  
  id_0_correspondencetable <- CorrespondenceTable %>% select(id_0) %>%
    distinct()%>% 
    pull()
  
  id_1_correspondencetable <- CorrespondenceTable %>% select(id_1) %>%
    distinct() %>% 
    pull()
  
  # 2. Load NR data, get all unique IDs from NR data, and load director data and 
  
  my_db_connect <- dbConnect(RSQLite::SQLite(), db_path) 
  nr_data <- dbGetQuery(my_db_connect, paste0("SELECT * FROM nr_data")) %>%
    select(ID,year)
  
  id_all <- nr_data  %>%
    distinct(ID) %>%
    arrange(ID) %>%
    pull()
  
  if (source=="bisnode"){ 
    data_input <- dbGetQuery(my_db_connect, "SELECT * FROM bisnode") %>% 
      filter(rechtsform: "Aktiengesellschaft" & gremium=="Verwaltungsrat" ) %>%
      filter(personenid %in% id_1_correspondencetable)
    }
  
  if (source=="sugarcube"){ 
    data_input <- dbGetQuery(my_db_connect, "SELECT * FROM sugarcube")%>%
      filter(personenid %in% id_1_correspondencetable)
    }

  dbDisconnect(my_db_connect)
  
  
  # 3. Write results for every ID in csv file
  
  id_all_batch <- id_all[((batch_number-1)*batch_size+1):(batch_number*batch_size)]
  id_all_batch <- id_all_batch[!is.na(id_all_batch)]

  for (i in 1:length(id_all_batch)){
    
    print(paste0("New ID: ", id_all_batch[i]))
    
    out <- QueryPoliticianMandates(pol_id=id_all_batch[i],
                                   data_input=data_input,
                                   correspondencetable=CorrespondenceTable,
                                   year_range=year_range,
                                   source=source)
                                
    out_right_format <- ReShapeBisnodeResults(pol_id=id_all_batch[i],
                                              nr_data=nr_data,
                                              input_data=out,
                                              max_period=86) 
    
    data.table::fwrite(
      x = out_right_format,
      file = paste0(path,"/Bisnode_Results_",modelname,"_",linkscore_cutoff,"_batch_",batch_number,".csv"),
      sep = ";",
      append=T
    )
  }
}



#----------------------
# C) Bisnode 
#----------------------

# a) Lookup

WriteModelResults(modelname="Generation_7",
                  year_range=c(1994:2017),
                  source="bisnode",
                  batch_number=1,
                  batch_size=100000,
                  linkscore_cutoff="optimal",
                  path="C:/Schmidlu/Dropbox/Projekt Nationalräte/02_Processed_data/11_Directors_1994_2018")


# b) Manual checks

CorrespondenceTable <- GetCorrespondenceTable(modelname="Generation_4",
                                              source=source,
                                              linkscore_cutoff=linkscore_cutoff) %>%
  arrange(id_0) %>% 
  distinct(id_0,id_1)

id_0_correspondencetable <- CorrespondenceTable %>% select(id_0) %>%
  distinct()%>% 
  pull()


N=100
set.seed(123)
out <- data.frame()
for (i in 1:N){
out <- bind_rows(out,data.frame(id_0=sample(id_0_correspondencetable,1),year=sample(1994:2017,1))) 
}
out <- out %>% arrange(id_0,year)


CorrespondenceTable_sample <- CorrespondenceTable %>% 
  filter(id_0 %in% out$id_0)

id_1_correspondencetable <- CorrespondenceTable_sample
  select(id_1) %>%
  distinct() %>% 
  pull()

my_db_connect <- dbConnect(RSQLite::SQLite(), db_path) 
data_to_check <- dbGetQuery(my_db_connect, "SELECT * FROM bisnode") %>% 
  filter(personenid %in% id_1_correspondencetable) %>%
  select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum)
dbDisconnect(my_db_connect)


path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte/02_Processed_data/11_Directors_1994_2018/"
data.table::fwrite(x = out,file = paste0(path,"Manual_Check_Random_Sample",".csv"), sep = ";")
data.table::fwrite(x = CorrespondenceTable_sample,file = paste0(path,"Manual_Check_Correspondence_Table",".csv"),sep = ";")
data.table::fwrite(x = data_to_check %>% arrange(personenid),file = paste0(path,"Manual_Check_Bisnode_Data",".csv"),sep = ";")


# c) Check specific cases

pol_id <- "AG-2011-0185"
pol_id=c("ZH-1979-0051") # Christoph Blocher

results <- data.table::fread("C:/Schmidlu/Dropbox/Projekt Nationalräte/02_Processed_data/11_Directors_1994_2018/Bisnode_Results_Generation_4_optimal.csv",
                                   sep = ";", encoding = "UTF-8")%>%
  filter(ID==pol_id)

CorrespondenceTable <- GetCorrespondenceTable(modelname="Generation_4",
                                              source="bisnode",
                                              linkscore_cutoff="optimal") %>%
  arrange(id_0) %>%
  filter(id_0==pol_id)

my_db_connect <- dbConnect(RSQLite::SQLite(), db_path) 
bisnode_data <- do.call(rbind.data.frame,lapply(CorrespondenceTable$id_1, function(y) {
  dbGetQuery(my_db_connect, paste0("SELECT * FROM bisnode WHERE personenid =='",y,"'"))}))
nr_data <- dbGetQuery(my_db_connect, paste0("SELECT * FROM nr_data")) %>%
  filter(ID==pol_id)
dbDisconnect(my_db_connect)

QueryPoliticianMandates(pol_id=pol_id,
                        correspondencetable=CorrespondenceTable,
                        data_input=bisnode_data,
                        source="bisnode",
                        year_range=c(1994:2017))

bisnode_data %>% select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum) %>%
  arrange(eintrittdatum) %>%
  distinct()

bisnode_data %>% select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum) %>%
  filter(duns==487025892) %>%
  arrange(eintrittdatum) %>%
  distinct()


# d) Check distribution of eintrittdatum, austrittdatum

my_db_connect <- dbConnect(RSQLite::SQLite(), db_path) 
bisnode_all <- dbGetQuery(my_db_connect, paste0("SELECT * FROM bisnode")) 
dbDisconnect(my_db_connect)

bisnode_all %>%
  janitor::tabyl(eintrittdatum) %>%
  arrange(eintrittdatum)

bisnode_all %>%
  janitor::tabyl(austrittdatum) %>%
  arrange(-percent)


#----------------------
# D) Sugarcube 
#----------------------

CorrespondenceTable <- GetCorrespondenceTable(modelname="Generation_1",
                                              source="sugarcube",
                                              linkscore_cutoff="optimal") %>%
  arrange(id_0)

test <- QueryPoliticianMandates(pol_id="AG-1931-0010",
                                correspondencetable=CorrespondenceTable,
                                year_range=c(1950:2015),
                                source="sugarcube")
View(test)  

CorrespondenceTable %>% filter(id_0=="AG-1931-0010") %>%
  arrange(year_1)


QueryPoliticianMandatesLoop(pol_ids=c("AG-1931-0010","AG-1931-0013"),
                            modelname="Generation_1",
                            year_range=c(1994:2017),
                            source="sugarcube",
                            linkscore_cutoff="optimal")



# Check functions 
# 
# out <- QueryPoliticianMandates(pol_ids=pol_id,
#                                    modelname="Generation_4",
#                                    year_range=c(1960:2017),
#                                    source="bisnode",
#                                    linkscore_cutoff="optimal")
# 
# out_reshaped <- out %>%
#   pivot_longer(cols=c(2:5),names_to="gremium_funktion",values_to="mandate")
# 
# Blocher <- ReShapeBisnodeResults(pol_id=pol_id,out=out_reshaped)
# 
# 
#  CorrespondenceTable <- GetCorrespondenceTable(modelname="Generation_4",
#                                                source="bisnode",
#                                                linkscore_cutoff="optimal") %>%
#  arrange(id_0)

# BisnodeTransformToLong(daf=bisnode_data,id_underconsid="2666115")

# test <- QueryPoliticianMandates(pol_id="AG-1931-0010",
#                                 correspondencetable=CorrespondenceTable,
#                                 source="bisnode",
#                                 year_range=c(1950:2015))
# test %>% tail(n=20)
# 
# ids_test <- CorrespondenceTable %>% 
#   filter(id_0=="AG-1931-0010") %>%  
#   select(id_1) %>% 
#   pull()
# 
# director_data %>% filter(personenid %in% ids_test) %>%
#   select(vorname,nachname,gremium,funktion,unterschrift,geburtstag)
# out <- do.call(rbind.data.frame,lapply(c("AG-1931-0010","AG-1931-0013"), function(y) {
#   dbGetQuery(my_db_connect, paste0("SELECT * FROM nr_data WHERE ID =='",y,"'"))}))
# 
# 
# data <- lapply(c("AG-1931-0013"), function(y) {
#   dbGetQuery(my_db_connect, paste0("SELECT * FROM nr_data WHERE ID =='",y,"'"))})
# 

