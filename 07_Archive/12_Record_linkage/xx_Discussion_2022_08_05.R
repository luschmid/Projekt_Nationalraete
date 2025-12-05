# Idee Mark on dataset

CorrespondenceTable <- GetCorrespondenceTable(modelname="Generation_4",
                                              source="bisnode",
                                              linkscore_cutoff="optimal") %>%
  arrange(id_0) %>%
  distinct(id_0,id_1)


id_0_correspondencetable <- CorrespondenceTable %>% select(id_0) %>%
  distinct()

id_1_correspondencetable <- CorrespondenceTable %>% select(id_1) %>%
  distinct() %>% 
  pull()

my_db_connect <- dbConnect(RSQLite::SQLite(), db_path) 
bisnode_data_all <- dbGetQuery(my_db_connect, "SELECT * FROM bisnode" )
dbDisconnect(my_db_connect)
bisnode_data_all_cut <- bisnode_data_all %>% filter(personenid)

# decision: save bisnode data before going into loop per person and keep only
# those observations in correspondence table
  
### check funktion and gremium

bisnode_data_all %>%
  janitor::tabyl(gremium) 

bisnode_data_all %>% filter(gremium=="Verwaltungsrat") %>%
  janitor::tabyl(funktion)  

bisnode_data_all %>% 
  filter(gremium=="Verwaltungsrat" ) %>%
  filter(personenid %in% id_1_correspondencetable) %>%
  janitor::tabyl(funktion)  


# decision: Mark looks whether gremium=="Verwaltungsrat" is enough for 
# restricting our sample to include only AGs (he will check the rechtsform variable)


### check president and member

presidents <- bisnode_data_all %>% 
  filter(gremium=="Verwaltungsrat" & funktion=="Präsident/in") %>%
  filter(personenid %in% id_1_correspondencetable) %>%
  select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum) %>%
  distinct()


presidents[c(80:140),]

bisnode_data_all %>% 
  filter(personenid==110 & duns==480003251)%>%
  select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum)


bisnode_data_all %>% 
  filter(personenid==499 & duns==488867599)%>%
  select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum)

bisnode_data_all %>% 
  filter(personenid==1512 & duns==480538842 )%>%
  select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum)


bisnode_data_all %>% 
  filter(personenid==1417 & duns==480065978 )%>%
  select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum)


bisnode_data_all %>% 
  filter(personenid==1157 & duns==481479590 )%>%
  select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum)


bisnode_data_all %>% 
  filter(personenid==1781 & duns==481692036 )%>%
  select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum)

director_data_cut <- bisnode_data_all %>% 
  filter(personenid==2130 & duns==482620648 )%>%
  select(duns,personenid,vorname,nachname,gremium,funktion,eintrittdatum,austrittdatum) 


# decision: not consistent coding (sometimes overlap, sometimes not), we need to
# address this directly in the code# 

