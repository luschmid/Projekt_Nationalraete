#-------------------------------------------------------
# (A) Load packages, settings, and set working directory
#-------------------------------------------------------

library(tidyverse)
library(readstata13)
library(janitor)
library(scales)
library(ggExtra)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))

#rm(list = ls(all = TRUE))

# -------------------------------------------
# (B) Read in data and calculate vote margin
#--------------------------------------------

# (i) Read in data and look at districts and years

data <- read.dta13("./02_Processed_data/nationalraete_1931_2015.dta") %>% 
  select(year,canton, list, ID,eligible_cant,elected,votes,pvotes,
         alliance, suballiance) %>%
  group_by(year,canton, list) %>%
  mutate(rank_h = 
           dense_rank(-as.numeric(votes+elected+row_number()/1000))) 
  

data_ch_prep <- PrepareData(data=data, 
                                votes_j_name="pvotes",
                                party_name="list", 
                                districtname="canton", 
                                election_cyclename = "year",
                                alliances=TRUE,
                                system="closed",
                                alliance_name="alliance", 
                                suballiance_name="suballiance") 


switzerland_margins_closedlist <- CalculateMargins(data_ch_prep,
                                   system="closed",
                                   method="dHondt",
                                   rank="rank_h",
                                   margin_type="add",
                                   additional_vars=c("ID",
                                                     "year",
                                                     "canton"))

data.table::fwrite(switzerland_margins_closedlist, 
          file = "./02_Processed_data/13_Running_variable/RV_CH_Closedlist.csv", 
          sep = ";",row.names = F)


# (ii) Merge datasets and recode marginal candidates

switzerland_margins_closedlist <- data.table::fread(
             file = "./02_Processed_data/13_Running_variable/RV_CH_Closedlist.csv", 
             sep = ";")


data_elected_perparty <- data %>%
                         group_by(year,canton, list) %>%
                         summarize(elected_sum=sum(elected,na.rm=T))

data_all <- data %>% left_join(switzerland_margins_closedlist %>%
                                select(ID,year,votemargin),
                               by=c("ID","year")) %>%
                     left_join(data_elected_perparty,
                               by=c("year","canton","list"))%>%
                     arrange(year,canton, list,rank_h) %>%
                     ungroup() %>%
                     mutate(sample=
                            ifelse(rank_h==elected_sum | rank_h==elected_sum+1 ,
                                   "Marginal candidate",
                                  "Not marginal candidate"),
                            votemargin_rel=votemargin/eligible_cant)


# (iii) Plausibility checks

data_all %>% group_by(elected) %>%
  dplyr::summarize(margin_min=min(votemargin,na.rm=T), 
                   margin_max=max(votemargin,na.rm=T))

# Result: looks fine, all elected candidates have positive margin, all not 
#         elected candidates have negative margin

data_all %>% filter(sample=="Marginal candidate") %>% print(n=100) %>%
             select(year,canton,list,votes,pvotes,rank_h,votemargin,
                    elected_sum,sample )

# Result: looks fine, votemargin changes from positive to negative for all 
#         parties with elected candidates and is negative for all parties with
#         no elected candidate
  
data_all %>% filter(is.na(votemargin))

# Result: 63 observations, all in small cantons

data_nobs <- data_all %>% group_by(canton,year) %>% summarize(nobs=n()) %>%
             ungroup()

data_check <- data_all %>% left_join(data_nobs,
                      by=c("year","canton")) %>%
              filter(is.na(votemargin) & nobs>1) %>% 
              select(year,canton,nobs)
# Note: observations with more than one observation and missing votemargin

data_all %>% filter(year==1943 & canton=="NW")                         

# (iv) Plot

ggplot(data_all, aes(votemargin_rel,fill=factor(sample))) +
  geom_histogram(aes(y = ..count.. / sum(..count..)),
                 alpha=0.5,position = "stack",boundary=0,color="black") +
  theme_bw(base_size = 28) +
  scale_fill_manual(values=c("black","grey50")) +
  scale_x_continuous(labels=percent) +
  geom_vline(xintercept = 0) +
  ylab("") + xlab("")+
  scale_x_continuous(limits=c(-0.5,0.5))+
  theme(axis.text.y = element_text(angle = 0)) +
  theme(legend.position="bottom")+
  theme(legend.title = element_blank())

ggsave(file = "./05_Texts_and_presns/01_Running_variable/figures/hist_switzerland_50.pdf", width = 20, height = 12)


# check functions
# data_ag_1931 <- data_ch_prep %>% filter(district=="AG" & year==1931)
# data_input <- data_ch_prep %>% filter(district=="AG" & year==1931)
# data_ag_1931 %>% distinct(votes_j,.keep_all = T)
# data_BE_1935 <- data_ch_prep %>% filter(district=="BE" & year==1935) 
# data_BE_1935 %>% distinct(votes_j,.keep_all = T)
# 
# data <- data_BE_1935 %>% distinct(pvotes,.keep_all = T)
# 
# CalculatePartyMargin(party_input=341,
#                      data_input=data_BE_1935,
#                      no_seats=31,
#                      method="dHondt",
#                      margin_type="add")

