
data_minyear <-    data %>%
  group_by(ID)%>%
  summarize(year_min=min(year))

data_canton_votes <-    data %>%
  group_by(canton,year)%>%
  summarize(votes_sum=sum(votes),
            no_seats=sum(elected))

data <- data %>% inner_join(data_minyear,by=c("ID")) %>% 
  inner_join(data_canton_votes,by=c("canton","year")) %>% 
  mutate(votes_rel=votes/(votes_sum/no_seats), 
         votemargin_rel=votemargin/(votes_sum/no_seats),
         age=year-birthyear)  %>%
  mutate(sex=as.numeric(as.factor(data$sex))-1)



data_forward <- data  %>%
  group_by(ID) %>% mutate(votes_f1=lead(votes),
                          votes_rel_f1=lead(votes_rel),
                          elected_f1=lead(elected))  %>% 
  select(ID,year,ends_with("f1")) %>% 
  mutate(participation_f1=ifelse(elected_f1%in%c(0,1),1,0), 
         elected_f1=ifelse(is.na(elected_f1),0,elected_f1))


data_lag <- data %>% inner_join(data_canton_votes,by=c("canton","year")) %>%
  group_by(ID) %>% 
  mutate(votes_l1=lag(votes),
         votes_rel_l1=lag(votes_rel),
         elected_l1=lag(elected))  %>% 
  select(ID,year,ends_with("l1"))


data_working<- data %>% inner_join(data_forward,by=c("ID","year")) %>%
  inner_join(data_lag,by=c("ID","year")) %>%
  select(name, firstname, ID,year,year_min,canton, leftright,sex,age, votemargin, votemargin_rel,participation_f1,starts_with("votes"),starts_with("elected"),cand_before1931)%>%
  filter(elected==0 & votemargin_rel<0|elected==1 & votemargin_rel>0) # temporary: remove four cases FR-1931-9003 and FR-1935-0013 in 1939 as well as ZH-1943-0030 and in 1967


data_working_ct<-data_working%>%select(canton)

dummy <- dummyVars(~ ., data =data_working_ct ) # make dummy for canton
data_working <- cbind(data_working,predict(dummy, data_working_ct))

rm(data_minyear,data_canton_votes,data_forward,data_lag,data,nr_1931_2015,rv_analytical_all) 

write_csv(data_working,"Running variable/Data/data_working.csv",na=".")