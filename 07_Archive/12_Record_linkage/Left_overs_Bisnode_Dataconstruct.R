%>%
  mutate(funktion=str_replace(funktion,"\\/in",""))

out_final1_wide <- out_final1 %>% pivot_wider(names_from=c(gremium,funktion),
                                              id_cols=year,
                                              values_from=mandates,
                                              values_fill=0)


my_db_connect <- dbConnect(RSQLite::SQLite(), db_path) 
dbListTables(my_db_connect)
data <- lapply(correspondencetable_id$id_1, function(y) {
  dbGetQuery(my_db_connect, paste0("SELECT * FROM bisnode WHERE personenid ==",y))})
dbDisconnect(my_db_connect)


as.data.frame(do.call(rbind, data_input))

out_reshaped <- out %>%
  pivot_longer(cols=c(3:dim(out)[2]),names_to="category",values_to="mandates") 



# QueryPoliticianMandatesLoop <- function(pol_ids,
#                                     modelname,
#                                     data_input,
#                                     year_range,
#                                     source,
#                                     linkscore_cutoff){
# 
# # This function returns the number of connections per year for several 
# # politician ids (pol_ids) for a certain year range (year_range). The modelname 
# # is the RL model (e.g., "Generation_1"). The source is either "bisnode" or 
# # "sugarcube. 
# 
#   CorrespondenceTable <- GetCorrespondenceTable(modelname=modelname,
#                                                 source=source,
#                                                 linkscore_cutoff=linkscore_cutoff) %>%
#     arrange(id_0) 
# 
# 
# out <- data.frame()
# for (i in 1:length(pol_ids)){
# out<- bind_rows(out,QueryPoliticianMandates(pol_id=pol_ids[i],
#                                             data_input=data_input,
#                                   correspondencetable=CorrespondenceTable,
#                                   year_range=year_range,
#                                   source=source))
# }
# return(out)
# }


# 
#WriteModelResultsout
#QueryPoliticianMandatesLoop(pol_ids=id_all[i],
#                                    data_input=data_input,
#                                    modelname=modelname,
#                                    year_range=year_range,
#                                    source=source,
#                                    linkscore_cutoff=linkscore_cutoff)%>%
#   select(id_0,year,everything())