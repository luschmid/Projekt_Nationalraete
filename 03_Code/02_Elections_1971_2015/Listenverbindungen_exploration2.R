
# (A) Load packages, settings, and set working directory ----

library(tidyverse)
library(readstata13)
library(janitor)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))

#rm(list = ls(all = TRUE))
source("./03_Code/13_Running_variable/00_Vote_margin_analytical_highestaverage.R") 


# (B) Read in data and calculate vote margin ----


# (i) Read in data and look at districts and years ----

data <- read.dta13("./02_Processed_data/nationalraete_1931_2015.dta") %>%
  filter(year>=1971)%>%
  as_tibble() %>%
  mutate(pvotes=ifelse(is.na(pvotes)==F,pvotes,votes))

# (ii) Calculate seat allocation for current situation ----

data_swi_prepared <- PrepareData(data %>% select(-votemargin),
                                 votes_h_name="votes", 
                                 votes_j_name= "pvotes",
                                 alliance_name="alliance",
                                 suballiance_name="suballiance",
                                 party_name="list",
                                 districtname="canton",
                                 election_cyclename = "year",
                                 system="open",
                                 cand_id_name="ID_pers",
                                 alliances=TRUE) %>% 
  filter(is.na(votes_h) == F) # Note: filter out tacit elections


seat_allocation_currenct <- GetSeatsHighestAverageLoop(data_input = data_swi_prepared ,
                                       method="dHondt") %>%
  as_tibble()

# quality check: compare our calculated seat allocation with seat allocation in the data

data_check <- data_swi_prepared %>% 
  group_by(district,year,party,listname_bfs,list) %>% 
  summarize(elected_actual=sum(elected,na.rm=TRUE)) %>% 
  left_join(seat_allocation_currenct %>% 
              rename(district=districtid) %>% 
              mutate(year=as.numeric(as.character(year)))%>% 
              mutate(party=as.numeric(as.character(party))),
            by=c("district","year","party"))

data_check %>% filter(elected_actual!=seats)

# result: no difference in seat allocation 

# (iii) Calculate seat allocation if there were no alliances ----

data_swi_prepared <- PrepareData(data %>% select(-votemargin),
                                 votes_h_name="votes", 
                                 votes_j_name= "pvotes",
                                 alliance_name="alliance",
                                 suballiance_name="suballiance",
                                 party_name="list",
                                 districtname="canton",
                                 election_cyclename = "year",
                                 system="open",
                                 cand_id_name="ID_pers",
                                 alliances=FALSE) %>% 
  filter(is.na(votes_h) == F) # Note: filter out tacit elections


seat_allocation_no_alliance1<- GetSeatsHighestAverageLoop(data_input = data_swi_prepared ,
                                                       method="dHondt") %>%
  as_tibble()

seat_allocation_no_alliance1 <- seat_allocation_no_alliance1 %>%
  rename(seats_no_alliance1=seats)

data_swi_prepared <- PrepareData(data %>% select(-votemargin) %>% mutate(alliance="",suballiance=""),
                                 votes_h_name="votes", 
                                 votes_j_name= "pvotes",
                                 alliance_name="alliance",
                                 suballiance_name="suballiance",
                                 party_name="list",
                                 districtname="canton",
                                 election_cyclename = "year",
                                 system="open",
                                 cand_id_name="ID_pers",
                                 alliances=TRUE) %>% 
  filter(is.na(votes_h) == F) # Note: filter out tacit elections


seat_allocation_no_alliance2<- GetSeatsHighestAverageLoop(data_input = data_swi_prepared ,
                                                          method="dHondt") %>%
  as_tibble()

seat_allocation_no_alliance2 <- seat_allocation_no_alliance2 %>%
  rename(seats_no_alliance2=seats)

# quality check: compare our seat allocation of the two versions with no alliance

seat_allocation_no_alliance1 %>% 
  left_join(seat_allocation_no_alliance2, 
            by=c("districtid","year","party")) %>%
  filter(seats_no_alliance1!=seats_no_alliance2)


# result: no difference in seat allocation

# (iv) Construct one dataset

data_all <- data_swi_prepared %>% 
  group_by(district,year,party,listname_bfs,list) %>% 
  summarize(elected_actual=sum(elected,na.rm=TRUE)) %>% 
  left_join(seat_allocation_currenct %>% 
              rename(district=districtid) %>% 
              mutate(year=as.numeric(as.character(year)))%>% 
              mutate(party=as.numeric(as.character(party))),
            by=c("district","year","party"))%>% 
  left_join(seat_allocation_no_alliance1 %>% 
              rename(district=districtid) %>% 
              mutate(year=as.numeric(as.character(year)))%>% 
              mutate(party=as.numeric(as.character(party))),
            by=c("district","year","party"))%>% 
  left_join(seat_allocation_no_alliance2 %>% 
              rename(district=districtid) %>% 
              mutate(year=as.numeric(as.character(year)))%>% 
              mutate(party=as.numeric(as.character(party))),
            by=c("district","year","party")) 



# (C) Plot results

data_years <- data_all %>%
  group_by(year) %>%
  summarize(seats_total=sum(seats,na.rm = TRUE)) %>%
  ungroup()  

data_perparty_years <- data_all %>%
  left_join(data_years,by=c("year")) %>%
  mutate(listname_bfs=ifelse(listname_bfs=="LPS/PLS","FDP/PLR (PRD)",listname_bfs))%>% 
  group_by(listname_bfs,year) %>%
  summarize(seats=sum(seats,na.rm = TRUE),
            seats_no_alliance1=sum(seats_no_alliance1,na.rm = TRUE),
            seats_total=mean(seats_total)) %>%
  mutate(seats_share=seats/seats_total,
         seats_no_alliance1_share=seats_no_alliance1/seats_total)
  ungroup() 


data_party <- data_perparty_years %>%
  group_by(listname_bfs) %>%
  summarize(seats=mean(seats),
            seats_no_alliance1=mean(seats_no_alliance1),
            seats_share=mean(seats_share),
            seats_no_alliance1_share=mean(seats_no_alliance1_share)) %>% 
  filter(seats>5) %>%
  arrange(seats) %>%
  mutate(id=row_number()) %>%
  mutate(colorvar=case_when(listname_bfs=="SP/PS" ~ rgb(255, 0, 0, maxColorValue = 255),
                            listname_bfs=="FDP/PLR (PRD)" ~ rgb(0, 68, 212, maxColorValue = 255),
                            listname_bfs=="SVP/UDC" ~ rgb(0, 120, 50, maxColorValue = 255),
                            listname_bfs=="CVP/PDC" ~ rgb(247, 113, 3, maxColorValue = 255),
                            listname_bfs=="GPS/PES" ~ rgb(3, 242, 26, maxColorValue = 255),
                            listname_bfs=="BDP/PBD" ~ rgb(255, 220, 0, maxColorValue = 255),
                            listname_bfs=="GLP/PVL" ~ rgb(125, 55, 115, maxColorValue = 255))) %>%
  mutate(seats_loc=seats_share-sign(seats-seats_no_alliance1)*0.002)



ggplot(data=data_perparty_years , 
       aes(x=as.factor(listname_bfs),
           y=seats_share,
           color=as.factor(listname_bfs)
       )) +
  geom_point(size=4,shape = 1, stroke=3) + 
  geom_point(size=4,data=data_party, 
             aes(x=as.factor(listname_bfs),
                 y=seats_no_alliance1_share,
                 color=as.factor(listname_bfs)), shape  = 19, stroke=3)  +
  scale_x_discrete(breaks=as.factor(data_party$listname_bfs),
                   labels=as.factor(data_party$listname_bfs)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5L)) +
  scale_color_manual(breaks=as.factor(data_party$listname_bfs),
                   values=data_party$colorvar) +
  coord_flip() + 
  geom_segment(aes(x = as.factor(listname_bfs), y = seats_loc, xend = as.factor(listname_bfs), 
                   yend =  seats_no_alliance1_share),size=1) +
  theme_bw(base_size=24) +
  theme(legend.position = "none") +
  xlab("") + ylab("Anteil Sitze") +
  theme(axis.text.y = element_text(angle = 0)) +
  theme(legend.title = element_blank(),
        legend.box.just = "right",
        legend.key = element_rect(size = 6),
        legend.key.size = unit(2, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "white"), 
        axis.line.y = element_line(colour = "black")) + 
  facet_grid(~year)


# per party over time

ggplot(data=data_party , 
       aes(x=as.factor(id),
           y=seats_share,
           color=as.factor(id)
       )) +
  geom_point(size=4,shape = 1, stroke=3) + 
  geom_point(size=4,data=data_party, 
             aes(x=as.factor(id),
                 y=seats_no_alliance1_share,
                 color=as.factor(id)), shape  = 19, stroke=3)  +
  scale_x_discrete(breaks=as.factor(data_party$id),
                   labels=as.factor(data_party$listname_bfs)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5L)) +
  scale_color_manual(breaks=as.factor(data_party$id),
                     values=data_party$colorvar) +
  coord_flip() + 
  geom_segment(aes(x = id, y = seats_loc, xend = id, 
                   yend =  seats_no_alliance1_share),size=1) +
  theme_bw(base_size=24) +
  theme(legend.position = "none") +
  xlab("") + ylab("Anteil Sitze") +
  theme(axis.text.y = element_text(angle = 0)) +
  theme(legend.title = element_blank(),
        legend.box.just = "right",
        legend.key = element_rect(size = 6),
        legend.key.size = unit(2, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "white"), 
        axis.line.y = element_line(colour = "black")) 



# over time


listname_bfs_unique <- data_party %>% distinct(listname_bfs) %>% pull()

data_party_long <- data_perparty_years %>%
  pivot_longer(
    cols = c("seats_share","seats_no_alliance1_share"),
    names_to = "scenario",
    names_prefix = "wk",
    values_to = "seat_share",
    values_drop_na = TRUE
  ) %>%
  filter(listname_bfs %in% listname_bfs_unique)%>%
  mutate(colorvar=case_when(listname_bfs=="SP/PS" ~ rgb(255, 0, 0, maxColorValue = 255),
                            listname_bfs=="FDP/PLR (PRD)" ~ rgb(0, 68, 212, maxColorValue = 255),
                            listname_bfs=="SVP/UDC" ~ rgb(0, 120, 50, maxColorValue = 255),
                            listname_bfs=="CVP/PDC" ~ rgb(247, 113, 3, maxColorValue = 255),
                            listname_bfs=="GPS/PES" ~ rgb(3, 242, 26, maxColorValue = 255),
                            listname_bfs=="BDP/PBD" ~ rgb(255, 220, 0, maxColorValue = 255),
                            listname_bfs=="GLP/PVL" ~ rgb(125, 55, 115, maxColorValue = 255))) 



ggplot(data=data_party_long %>% filter(scenario=="seats_share") , 
       aes(x=year,
           y=seat_share,
           color=factor(listname_bfs))) +
  geom_point(size=3,shape=1)+
  geom_line(size=2)+
  geom_point(size=3,data=data_party_long %>% filter(scenario=="seats_no_alliance1_share") , 
             aes(x=year,
                 y=seat_share,
                 color=factor(listname_bfs)),shape=19)+
  geom_line(size=2,data=data_party_long %>% filter(scenario=="seats_no_alliance1_share") , 
            aes(x=year,
                y=seat_share,
                color=factor(listname_bfs)),linetype="dashed")  +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5L)) +
  scale_color_manual(breaks=as.factor(data_party_long$listname_bfs),
                     values=data_party_long$colorvar) +
  theme_bw(base_size=24) +
  
  theme(legend.title = element_blank(),
        legend.box.just = "right",
        legend.key = element_rect(size = 6),
        legend.key.size = unit(2, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black")) 


