
library(readstata13)
library(data.table)
library(tidyverse)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))

#rm(list = ls())

# (i) Read-in data

data <- read_dta("./01_Raw_data/Electoral data/Honduras/elections_hn_fuzzy_merge_final.dta")

MatchingHondurasAll <-readstata13::read.dta13("./02_Processed_data/15_Elections_Honduras/elections_hn_final.dta")%>% 
                      mutate(Partyname=Party)%>% 
                      dplyr::select(id_Stata,Candidate,Department,Year,Partyname, Votes,Department,Partyname,Elected)%>% 
                      rename(votes_h=Votes,
                             year=Year,
                             district=Department,
                             elected=Elected) 



honduras_seats <-readr::read_delim("./01_Raw_data/15_Elections_Honduras/honduras_seats.csv",
                                  delim=";")
                                        
MatchingHondurasAllByParty <- MatchingHondurasAll %>% 
                              group_by(district,year,Partyname) %>%
                              summarize(votes=sum(votes_h)) %>% 
                              mutate(party=row_number())
MatchingHondurasAll <- MatchingHondurasAll %>% left_join(MatchingHondurasAllByParty,by=c("district","year","Partyname"))


# Example: Atlantida 2009

MatchingHondurasPartyATL <- MatchingHondurasAllByParty %>% filter(district=="ATLANTIDA" & year==2009) 
MatchingHondurasATL <- MatchingHondurasAll %>% filter(district=="ATLANTIDA" & year==2009) 
MatchingHonduras_ATL_2009_PN <- MatchingHondurasAll %>% filter(district=="ATLANTIDA" & year==2009 & 
                                                                 Party=="PARTIDO NACIONAL DE HONDURAS") 

MatchingHondurasATLNACIONAL <-MatchingHondurasAll %>%filter(district=="ATLANTIDA" & year==2009 & Party=="PARTIDO NACIONAL DE HONDURAS") 

GetVSHare(votes=MatchingHondurasPartyATL$votes, n=8) %>% 
        left_join(MatchingHondurasPartyATL,by=c("party")) %>%
        select(party,seat,VS,party)

CalculateCandidateMargin(party_input="PARTIDO NACIONAL DE HONDURAS", 
                         data_input=MatchingHondurasATLNACIONAL)

MatchingHondurasAll %>% filter(district=="ATLANTIDA" & year==2009)%>% arrange(Party,-votes_h) %>% print(n=100)



# Scale up to all years, districts, parties, and candidates



CalculateMargins <- function(data_input, data_input_seats) {
  data_year_district <- data_input %>%
    group_by(year, district) %>%
    summarize(nobs = n()) %>%
    select(year, district)%>%
    ungroup()%>%
    mutate(id=row_number())
  Out <- data.frame()
  for (i in c(data_year_district$id)){
    #print(paste("Loop number: ",i,", Dimension of Outfile:",dim(Out)[1],",",dim(Out)[2],sep=""))
      dataworking <- data_input[data_input$year == data_year_district$year[i] & data_input$district == data_year_district$district[i], ]
      no_seats <- data_input_seats$no_seats[data_input_seats$year==data_year_district$year[i] & data_input_seats$district==data_year_district$district[i]]
      Out <- rbind(Out, as.data.frame(AggregateMarginsLargestRemainder(data_input=dataworking, no_seats=no_seats)))
    }
  return(Out)
}

# test functions
CalculateCandidateMargins(MatchingHonduras_ATL_2009_PN)
MatchingHonduras_ATL_2009_PN$votes_h

# function that loops over all years, districts, and parties (takes around 30-45minutes on Lukas' PC)

start_time <- Sys.time()
Margins <- CalculateMargins(data_input=MatchingHondurasAll , data_input_seats=honduras_seats)
end_time <- Sys.time()
end_time - start_time
write.table(Margins %>% select(rank_h ,votemargin,candmargin, Candidate,votes_h,year,district,Partyname ), 
            file = "./02_Processed_data/15_Elections_Honduras/Margins_Honduras.csv", sep = ";")



# check some cases if running variable is correct (see Excel file "Bsp Largest Remainder - Honduras.xlsx" in foler "07 Archive")


Margins <-read.table(file = "./02_Processed_data/15_Elections_Honduras/Margins_Honduras.csv", sep = ";",header=T)
Margins_old <-read.table(file = "./02_Processed_data/15_Elections_Honduras/Margins_Honduras_old.csv", sep = ";",header=T) %>%
                        select(Candidate,votes_h,year,district,votemargin)%>%
                       rename(votemargin_old=votemargin)

dim(Margins)
dim(Margins_old)

MarginsCompare <-Margins %>%
                  left_join(Margins_old,by=c("Candidate","votes_h","year","district")) %>%
                  mutate(votemargin_diff=votemargin_old-votemargin)
summary(MarginsCompare$votemargin_diff)
MarginsCompare %>% filter(is.na(votemargin_old)) 

# (a) Atlantida 2009


Margins %>% select(Candidate,year,Party,votemargin,candmargin,district,votes_h) %>%
  filter(year==2009 & district=="ATLANTIDA" ) %>% 
  arrange(Party,-votes_h)

PartyVotesATL2009 <- MatchingHondurasAllByParty %>%  dplyr::filter(year=="2009" & district=="ATLANTIDA")

GetVSHare(PartyVotesATL2009$votes, 8)
# result: party margin appears to be correct. 

#  Candidate year             Party                         votemargin candmargin  district votes_h
#  RODOLFO  IRIAS NAVAS 2009  PARTIDO NACIONAL DE HONDURAS      13395       6241  ATLANTIDA   47929
#47929-34534=13395 # this is the difference of the first ranked RODOLFO  IRIAS NAVAS to the fifth-ranked candidate
# result: vote margin is correct



# (b) Cases for which party margin (and not candidate margin) is binding

Margins%>%select(id_Stata,Candidate,year,Party,votemargin,candmargin,district) %>%
  filter(votemargin!=candmargin ) %>% head(n=20)

# ATLANTIDA-2009-0009, DIOGENES ORLANDO ALVAREZ RODAS 2009 PARTIDO INNOVACION Y UNIDAD SOCIAL DEMOCRATA

Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==2009 & district=="ATLANTIDA" & Party=="PARTIDO INNOVACION Y UNIDAD SOCIAL DEMOCRATA" ) %>% 
  arrange(Party,-votes_h)

MatchingHondurasAllByParty %>%  dplyr::filter(year=="2009" & district=="ATLANTIDA")
GetVSHare(PartyVotesATL2009$votes, 8)

#######################
# (c) checks with Simon
#######################


#rm(list = setdiff(ls(), lsf.str())) # remove all objects except functions, load objects at top after this command
#set.seed(1234)
#sort(sample(c(1:3025),20))
# 29  342  557  693  697  796  845  855  873 1539 1629 1827 1867 1869 1918 1993 2074 2500 2580 2759

#  DARIO ALEJANDRO MUNGUIA QUEZADA (Margins[29,])

Case <- Margins[29,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==2013 & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS1 <- GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)



#  WENDY VANESA FLORES SAUCEDA (Margins[342,])

Case <- Margins[342,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==2013 & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS2 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)



#  OSCAR HUMBERTO MEJIA HERNANDEZ (Margins[557,])

Case <- Margins[557,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS3 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)



#  IVETH OBDULIA MATUTE BETANCOURTH (Margins[693,])

Case <- Margins[693,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS4 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)



#  PEDRO GUEVARA (Margins[697,])


Case <- Margins[697,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS5 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)



#  LUIS ALONSO MEMBRENO CASTILLO (Margins[796,])


Case <- Margins[796,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS6 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)


#  BERNARDINO ESPINOZA PALMA (Margins[845,])

Case <- Margins[845,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS7 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)


#  GLORIA ARGENTINA BONILLA BONILLA (Margins[855,])

Case <- Margins[855,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS8 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)

#  OSCAR EMILIO CRUZ PINEDA (Margins[873,])

Case <- Margins[873,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS9 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)


#  SINDI CAROLINA VALDEZ MATUTE (Margins[1539,])

Case <- Margins[1539,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS10 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)


#  MARIO ALBERTO GARCIA MARTINEZ (Margins[1827,])

Case <- Margins[1827,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS11 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)



#  GREYDI MERARI LAGOS GUTIERREZ (Margins[1867,])

Case <- Margins[1867,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS12 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)


#  DAVID ARMANDO REYES OSORTO  (Margins[1869,])

Case <- Margins[1869,,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS13 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)

#  IRIS YAMILETH MENDOZA SANCHEZ (Margins[1918,])

Case <- Margins[1918,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS14 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)

#  MODESTO BANEGAS EUCEDA (Margins[1993,])

Case <- Margins[1993,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS15 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)


#  GILMA ESPERANZA ARDON ALBERTO (Margins[2074,])

Case <- Margins[2074,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS16 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)

#  MAURICIO ROBERTO DOMINGUEZ HERNANDEZ (Margins[2500,])

Case <- Margins[2500,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS17 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)


#  DINI ELIZABETH VASQUEZ PALMA (Margins[2580,])

Case <- Margins[2580,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS18 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)


#  FREDY GERARDO FLORES VALLECILLO (Margins[2759,])

Case <- Margins[2759,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS19 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)



#  RODRIGO AVILA LEVAS (EINERBEZIRK)


Case <- Margins[Margins$Candidate=="RODRIGO AVILA LEVAS" & Margins$year==2017,]

CandVotesExample <- Margins %>% select(Candidate,year,Party,district,votes_h,votemargin,candmargin) %>%
  filter(year==Case$year & district==Case$district & Party==Case$Party) %>% 
  arrange(Party,-votes_h)

PartyVotesExample <- MatchingHondurasAllByParty %>%  dplyr::filter(year==Case$year& district==Case$district)
(no_seats <- honduras_seats$no_seats[honduras_seats$district==Case$district  & honduras_seats$year==Case$year ])
VS20 <-GetVSHare(PartyVotesExample$votes, no_seats)
df <- as.data.frame(rbind(PartyVotesExample$votes))
colnames(df) <- PartyVotesExample$Partyname

write.table(df, 
            file = "./02_Processed_data/15_Elections_Honduras/temp.csv", sep = ";",row.names=F)

