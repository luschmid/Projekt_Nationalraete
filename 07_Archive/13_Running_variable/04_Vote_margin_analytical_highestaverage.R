
library(tidyverse, quietly = TRUE)


try(setwd("D:/SchmidLu/Dropbox/Projekt Nationalr?te"))
try(setwd("C:/Users/Lukas/Dropbox/Projekt Nationalr?te"))

source("./03_Code/fine_grid.R") # ggplot layers


###################
# (A) Functions
###################

#rm(list = ls(all = TRUE))
options("encoding" = "UTF-8")

PrepareData <- function(data, votes_h_name, votes_j_name, alliance_name=party_name, suballiance_name=party_name, party_name, districtname, election_cyclename = "year") {
  if (length(data[[alliance_name]]) == 0) stop("Error: Alliance variable not found. ")
  if (length(data[[suballiance_name]]) == 0) stop("Error: Suballiance variable not found. ")
  if (length(data[[party_name]]) == 0) stop("Error: Party variable not found. ")
  if (is.numeric(data[[votes_h_name]]) == F) stop("Error: Variable for individual votes (votes_h) contains non-numeric values.")
  if (is.numeric(data[[votes_j_name]]) == F) stop("Error: Variable for party votes (votes_j) contains non-numeric values. ")
  data$year <- data[[election_cyclename]]
  data$district <- data[[districtname]]
  data$alliance_orig <- as.character(data[[alliance_name]]) # write party name into free alliance/suballiance column for parties that are not in alliance/suballiance
  data$suballiance_orig <- as.character(data[[suballiance_name]])
  data$alliance <- data[[party_name]]
  data$suballiance <- data[[party_name]]
  data$party <- paste(as.character(data$district), as.character(data$year), data[[party_name]])
  data$alliance[data$alliance_orig != "" & data$alliance_orig != "."] <- data$alliance_orig[data$alliance_orig != "" & data$alliance_orig != "."]
  data$suballiance[data$suballiance_orig != "" & data$suballiance_orig != "."] <- data$suballiance_orig[data$suballiance_orig != "" & data$suballiance_orig != "."]
  data$alliance <- as.numeric(as.factor(data$alliance)) # make alliance, suballiance and party variable numeric
  data$suballiance <- as.numeric(as.factor(data$suballiance))
  data$party <- as.numeric(as.factor(data$party))
  data$alliance_dummy <- 0
  data$suballiance_dummy <- 0
  data$alliance_dummy[data$alliance_orig != "" & data$alliance_orig != "."] <- 1
  data$suballiance_dummy[data$suballiance_orig != "" & data$suballiance_orig != "."] <- 1
  data$votes_h <- data[[votes_h_name]]
  data$votes_j <- data[[votes_j_name]]
  data <- data %>%
    group_by(party) %>%
    mutate(rank_h_inv = dense_rank(as.numeric(votes_h))) %>%
    mutate(rank_h_inv_max = max(rank_h_inv)) %>%
    mutate(rank_h = rank_h_inv_max - rank_h_inv + 1)
  dataout <- data %>%
    arrange(alliance, suballiance, party) %>%
    select(-ends_with("_orig"), -rank_h_inv_max, -rank_h_inv)
  return(dataout)
}


CollapseVotes <- function(data_input, aggregation_level, aggregation_unit = NA) { #
  aggregation_level <- enquo(aggregation_level)
  aggregation_unit <- enquo(aggregation_unit)
  check_aggregation_level <- names(data_input %>% ungroup() %>% select(!!aggregation_level))
  if (check_aggregation_level == "suballiance") {
    data_input <- data_input %>% filter(alliance == !!aggregation_unit)
  }
  if (check_aggregation_level == "party") {
    data_input <- data_input %>% filter(suballiance == !!aggregation_unit)
  }
  data_out <- data_input %>%
    group_by(party) %>%
    filter(row_number() == 1) %>%
    group_by(!!aggregation_level) %>%
    summarize(votes_j = sum(votes_j)) %>%
    as.matrix()
  return(data_out)
}


CalculateDHondtNumbers <- function(votes, no_seats) {
  div <- matrix(rep(c(1:no_seats), length(votes[, 2])), nrow = length(votes[, 2]), byrow = T)
  DHondt_wide <- as.data.frame(cbind(votes[, 1], matrix(rep(as.numeric(votes[, 2]), no_seats), nrow = length(votes[, 2])) / div))
  colnames(DHondt_wide) <- c("party", paste(1:no_seats, sep = ""))
  DHondt_long <- gather(DHondt_wide, paste(1:no_seats, sep = ""), key = "seat", value = "dHondtno")
  return(DHondt_long)
}


PartyMargin <- function(party, data_input, no_seats) {
  party <- enquo(party)
  data_others <- data_input %>%
    filter(party != !!party) %>%
    mutate(dHondtno_rank = dim(data_input)[1] - dense_rank(as.numeric(dHondtno)) - no_seats + 1) %>%
    filter(dHondtno_rank <= no_seats) %>%
    mutate(seat = no_seats - dHondtno_rank + 1) %>%
    select(dHondtno, seat)
  data_out <- data_input %>%
    filter(party == !!party) %>%
    mutate(seat = as.numeric(seat)) %>%
    left_join(data_others, by = c("seat")) %>%
    mutate(margin = (as.numeric(dHondtno.x) - as.numeric(dHondtno.y)) * seat) %>%
    select(party, seat, margin)
  return(data_out)
}

GetAllianceSuballiance <- function(party_input, data_input) {
  party_input <- enquo(party_input)
  alliance_dummy <- data_input %>%
    filter(party == !!party_input) %>%
    ungroup() %>%
    select(alliance_dummy) %>%
    filter(row_number() == 1) %>%
    pull()
  suballiance_dummy <- data_input %>%
    filter(party == !!party_input) %>%
    ungroup() %>%
    select(suballiance_dummy) %>%
    filter(row_number() == 1) %>%
    pull()
  alliance <- data_input %>%
    filter(party == !!party_input) %>%
    ungroup() %>%
    select(alliance) %>%
    filter(row_number() == 1) %>%
    pull()
  suballiance <- data_input %>%
    filter(party == !!party_input) %>%
    ungroup() %>%
    select(suballiance) %>%
    filter(row_number() == 1) %>%
    pull()
  return(data.frame(alliance_dummy = alliance_dummy, suballiance_dummy = suballiance_dummy, alliance = alliance, suballiance = suballiance))
}


CalculatePartyMargin <- function(party_input, data_input, no_seats) {
  print(paste("Party:", party_input))
  AllianceInformation <- GetAllianceSuballiance(party_input, data_input) # get info whether party is in sub-/alliance and in which
  # 1. Calculation of margin across alliances
  VotesAlliance <- CollapseVotes(data_input, aggregation_level = alliance)
  AllianceDHondtNumbers <- CalculateDHondtNumbers(VotesAlliance, no_seats)
  Alliancemargin <- PartyMargin(AllianceInformation$alliance, AllianceDHondtNumbers, no_seats)
  colnames(Alliancemargin) <- c("alliance", "seat_a", "margin_a")
  # 2. Calculation of margin across suballiances (only for parties in alliance)
  if (AllianceInformation$alliance_dummy == 1) {
    Suballiancemargin_all <- data.frame()
    for (i in c(1:no_seats)) {
      VotesSuballiance <- CollapseVotes(data_input, aggregation_level = suballiance, aggregation_unit = AllianceInformation$alliance)
      SuballianceDHondtNumbers <- CalculateDHondtNumbers(VotesSuballiance, i)
      Suballiancemargin <- PartyMargin(AllianceInformation$suballiance, SuballianceDHondtNumbers, i)
      Suballiancemargin$seat_a <- i
      colnames(Suballiancemargin) <- c("suballiance", "seat_s", "margin_s", "seat_a")
      Suballiancemargin_all <- rbind(Suballiancemargin_all, Suballiancemargin)
    }
    # 3. Calculation of margin across parties (only for parties in suballiance)
    if (AllianceInformation$alliance_dummy == 1 & AllianceInformation$suballiance_dummy == 1) {
      Partymargin_all <- data.frame()
      for (i in c(1:no_seats)) {
        VotesParty <- CollapseVotes(data_input, aggregation_level = party, aggregation_unit = AllianceInformation$suballiance)
        PartyDHondtNumbers <- CalculateDHondtNumbers(VotesParty, i)
        Partymargin <- PartyMargin(party_input, PartyDHondtNumbers, i)
        Partymargin$seat_s <- i
        colnames(Partymargin) <- c("party", "seat_p", "margin_p", "seat_s")
        Partymargin_all <- rbind(Partymargin_all, Partymargin)
      }
    }
  }
  # 4. Aggregation
  # Case 1: No alliance
  if (AllianceInformation$alliance_dummy == 0) {
    colnames(Alliancemargin) <- c("party", "seat_p", "partymargin")
    Alliancemargin$party <- party_input
    AllpartymarginsOut <- Alliancemargin
    AllpartymarginsOut <- AllpartymarginsOut %>% rename(seat = seat_p)
  }
  # Case 2: Only alliance but no suballiance
  if (AllianceInformation$alliance_dummy == 1 & AllianceInformation$suballiance_dummy == 0) {
    Allpartymargins <- Suballiancemargin_all %>% left_join(Alliancemargin, by = c("seat_a"))
    Allpartymargins <- Allpartymargins %>% mutate(margin_min = pmin(margin_a, margin_s))
    Allpartymargins$party <- party_input
    AllpartymarginsOut <- Allpartymargins %>% filter(seat_s <= seat_a) # keep only seats that are possible to make for suballiance
    AllpartymarginsOut <- AllpartymarginsOut %>%
      group_by(seat_s) %>%
      summarize(partymargin = max(margin_min), party = mean(party))
    AllpartymarginsOut <- AllpartymarginsOut %>% rename(seat = seat_s)
  }
  # Case 3: Alliance and suballiance
  if (AllianceInformation$alliance_dummy == 1 & AllianceInformation$suballiance_dummy == 1) {
    Allpartymargins <- Partymargin_all %>%
      left_join(Suballiancemargin_all, by = c("seat_s")) %>%
      left_join(Alliancemargin, by = c("seat_a"))
    Allpartymargins <- Allpartymargins %>% mutate(margin_min = pmin(margin_a, margin_s, margin_p))
    Allpartymargins$party <- party_input
    AllpartymarginsOut <- Allpartymargins %>% filter(seat_p <= seat_s & seat_s <= seat_a) # keep only seats that are possible to make for suballiance and party
    AllpartymarginsOut <- AllpartymarginsOut %>%
      group_by(seat_s) %>%
      summarize(partymargin = max(margin_min), party = mean(party))
    AllpartymarginsOut <- AllpartymarginsOut %>% rename(seat = seat_s)
  }
  return(AllpartymarginsOut)
}

CalculatePartyMargins <- function(data_input, no_seats) {
  Out <- data.frame()
  all_parties <- as.numeric(levels(as.factor(data_input$party)))
  for (j in 1:length(all_parties)) {
    Out <- rbind(Out, CalculatePartyMargin(all_parties[j], data_input, no_seats))
  }
  return(Out)
}

CalculateCandidateMargin <- function(party_input, data_input) {
  data <- data_input[data_input$party == party_input, ]
  out <- data.frame()
  out_all <- data.frame()
  if (dim(data)[1] > 1) { # calculate candidate margin only for parties with more than one candidate
    candidates <- max(data$rank_h)
    # print(paste("Party:",party_input))
    for (h in 1:candidates) {
      out_h <- data.frame()
      for (hcom in 1:(candidates + 1)) {
        if (h != hcom & hcom < (candidates + 1)) {
          candmargin <- data$votes_h[data$rank_h == h] - data$votes_h[data$rank_h == hcom]
          out <- data.frame(candmargin = candmargin, rank_h = h, rank_h_compare = hcom)
          out_h <- rbind(out_h, out)
        }
        if (h != hcom & hcom == (candidates + 1)) {
          out <- data.frame(candmargin = data$votes_h[data$rank_h == h], rank_h = h, rank_h_compare = hcom)
          out_h <- rbind(out_h, out)
        }
      }
      out_all <- rbind(out_h, out_all)
    }
  } else {
    out_all <- data.frame(candmargin = NA, rank_h = 1, rank_h_compare = 1)
  } # candidate margin is missing for parties with only one candidate
  out_all$party <- party_input
  return(out_all)
}

CalculateCandidateMargins <- function(data_input) {
  party_inputs <- as.numeric(levels(as.factor(data_input$party)))
  ncand <- dim(data_input)[1]
  out_all <- data.frame()
  for (j in party_inputs) {
    out_all <- rbind(CalculateCandidateMargin(data_input = data_input, party_input = j), out_all)
  }
  return(out_all)
}

AggregateMargins <- function(data_input, no_seats) {
  PartyMarginOut <- CalculatePartyMargins(data_input, no_seats) %>% mutate(mergevar = seat)
  CandidateMarginsOut <- CalculateCandidateMargins(data_input) %>% mutate(mergevar = ifelse(rank_h < rank_h_compare, rank_h_compare - 1, rank_h_compare))
  MarginsMerged <- CandidateMarginsOut %>%
    left_join(PartyMarginOut, by = c("party", "mergevar")) %>%
    mutate(margin_min = ifelse(is.na(candmargin) == F, pmin(candmargin, partymargin))) # minimum of candidate and partymargin for all canidates for whom both are available, partymargin for those with no candidate margin (single-candidate lists)
  MarginsFinal <- MarginsMerged %>%
    group_by(party, rank_h) %>%
    summarize(votemargin = max(margin_min),candmargin=min(candmargin)) # maximum of all potential seat configurations
  #PartyMargintoMerge <- MarginsMerged %>% filter(rank_h == seat)
  return(MarginsFinal %>% left_join(data_input, by = c("party", "rank_h")) %>% 
           #left_join(PartyMargintoMerge, by = c("party", "rank_h")) %>%
             select(party, rank_h, canton, year, firstname, name, votemargin,candmargin, pvotes, votes_h, elected, everything()))
}

CalculateMargins <- function(data_input) {
  data_singlecandidate <- data_input %>%
    group_by(year, district) %>%
    summarize(nobs = n()) %>%
    filter(nobs == 1)
  data_input <- data_input %>%
    left_join(data_singlecandidate, by = c("year", "district")) %>%
    filter(is.na(nobs) == T) # do calculation only for districts with more than one candidate
  years <- as.numeric(levels(as.factor(data_input$year)))
  Out <- data.frame()
  for (yr in years) {
    print(paste("New year:", yr))
    dataperyear <- data_input[data_input$year == yr, ]
    districts <- levels(as.factor(dataperyear$district))
    for (di in districts) {
      print(paste("New district:", di))
      dataworking <- dataperyear[dataperyear$district == di, ]
      no_seats <- sum(dataworking$elected)
      Out <- rbind(Out, as.data.frame(AggregateMargins(dataworking, no_seats)))
    }
  }
  return(Out)
}


###############################
# (B) Test of Functions
###############################

# as.factor(data_sz_out$alliance)
# as.factor(data_sz_out$suballiance)
# as.factor(data_sz_out$party)

# votes <- CollapseVotes(data_sz_out,pvotes,aggregation_level=alliance)
# dHondtout <- CalculateDHondtNumbers(votes,4)
# PartyMargin(4,dHondtout,4)


# Test of Allianceout
#
# GetAllianceSuballiance(13,data)
# data %>% filter(party==1)
# data %>% filter(suballiance==2)

# Test of CalculatePartyMargin
# CalculatePartyMargin(11,data_sz_out,4)
# data %>% filter(party==11)
#
# for (i in 1:11){
# print(paste("Party:",i))
# print(CalculatePartyMargin(i,data_sz_out,4))
# }
#
# (PartyMarginOut <- CalculatePartyMargins(data_sz_out,4))
#
# (CandidateMarginOut <- CalculateCandidateMargin(13,data_sz_out))
# (CalculateCandidateMarginsOut <- CalculateCandidateMargins(data_sz_out))
# min(CalculateCandidateMarginsOut$candmargin[CalculateCandidateMarginsOut$rank_h<CalculateCandidateMarginsOut$rank_h_compare])
# max(CalculateCandidateMarginsOut$candmargin[CalculateCandidateMarginsOut$rank_h>CalculateCandidateMarginsOut$rank_h_compare])
#
# PartyMarginOut <- PartyMarginOut %>% rename(mergevar=seat)
# CalculateCandidateMarginsOut <- CalculateCandidateMarginsOut %>% mutate(mergevar=ifelse(rank_h<rank_h_compare,rank_h_compare-1,rank_h_compare))
#
# MarginsMerged <- CalculateCandidateMarginsOut %>%
#                               left_join(PartyMarginOut,by=c("party","mergevar")) %>%
#                               left_join(data_sz_out,by=c("party","rank_h")) %>%
#                               select(firstname,name,candmargin,partymargin,votes,party,rank_h, rank_h_compare)
#
# (MarginsMerged <- MarginsMerged  %>%  mutate(margin_min = pmin(candmargin,partymargin)))
#
#
# (MarginsFinal <- MarginsMerged %>% group_by(party,rank_h) %>% summarize(votemargin=max(margin_min))  )%>% print(n=100)
#


###############################
# (C) Test of RV Functions
###############################

df <- data.frame(
  party = c(rep("A", 3), rep("B", 3), rep("C", 3)),
  votes = c(22, 8, 5,20, 15, 10, 12, 5, 3),
  firstname = c( "Brenda", "Charles", "Diana", "Z", "X", "Y","V", "W", "U"),
  name = rep("", 9),
  elected = c( 1, 0, 0, 1, 1, 0,0, 0, 0)
) %>%
  group_by(party) %>%
  mutate(pvotes = sum(votes), alliance = "", suballiance = "", canton = 1, year = 2013)

df_out <- PrepareData(df, "votes", "pvotes", "alliance", "suballiance", "party", "canton", "year")
AggregateMargins(df_out, 3)

# for aggregation table in presentation (from AggregateMargins function)

Margins_out <- MarginsMerged  %>% arrange(rank_h,seat) %>% 
                  select(rank_h,seat,partymargin,candmargin)%>%
                  mutate(rank_h=factor(rank_h,labels=c("Brenda","Charles","Diana")))
library(xtable)
print(xtable(Margins_out,digits=0),include.rownames=FALSE)

PartyMarginOut <- PartyMarginOut %>% mutate( party= c(rep("A", 3), rep("B", 3), rep("C", 3)))%>%
                  select(-mergevar)

print(xtable(PartyMarginOut,digits=c(0,0,0,1)),include.rownames=FALSE)





# (i) Run program for selected canton-years (SZ 2015 with all possible alliance/suballiance combinations, AG 1931 only with alliances )

data_all <- read.dta13("./02_Processed_data/nationalraete_1931_2015.dta") %>%
  as_tibble() %>%
  select(ID, canton, year, firstname, name, votes, pvotes, list, alliance, suballiance, elected) %>%
  mutate(party_orig = list, alliance_orig = alliance, suballiance_orig = suballiance)

data_sz_2015 <- data_all %>% filter(canton == "SZ" & year == 2015)
data_sz_2015 %>% print(n = 100)

data_ag_1931 <- data_all %>% filter(canton == "AG" & year == 1931)
data_ar_1931 <- data_all %>% filter(canton == "AR" & year == 1931)


# data <- data_sz_2015 # temporary

data_sz_2015_out <- PrepareData(data_sz_2015, "votes", "pvotes", "alliance_orig", "suballiance_orig", "list", "canton", "year")
data_ag_1931_out <- PrepareData(data_ag_1931, "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year")
data_ar_1931_out <- PrepareData(data_ar_1931, "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year")


AggregateMargins(data_sz_2015_out, 4)
AggregateMargins(data_ag_1931_out, 12)
AggregateMargins(data_ar_1931_out, 2)



# (ii) Run program for full dataset

data_prepared_1931_2015 <- PrepareData(data_all, "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year") %>% filter(is.na(votes_h) == F) # filter out "stille Wahlen"
data_prepared_1931_2015 %>%
  group_by(year, canton) %>%
  summarize(nobs = n()) %>%
  filter(nobs == 1)

start_time <- Sys.time()
RV_out <- CalculateMargins(data_prepared_1931_2015)
end_time <- Sys.time()
end_time - start_time


write.table(RV_out, file = "./02_Processed_data/13_Running_variable/RV_Analytical_R.csv", sep = ";")


# (iii) Check aggregate values

RV_out %>% filter(is.na(votemargin))

RV_out %>%
  filter(is.na(votemargin) == F) %>%
  group_by(year, elected) %>%
  summarize(votemargin_min = min(votemargin), votemargin_max = max(votemargin)) %>%
  print(n = 100)




# (iii) Check R implementation with previous implementation from Stata

RV_out <- read.table(file = "./02_Processed_data/13_Running_variable/RV_Analytical_R.csv", sep = ";")
head(RV_out)

data_rv_stata <- read.dta13("./02_Processed_data/13_Running_variable/RV_Analytical_all.dta") %>%
  as_tibble() %>%
  arrange(year, party_num, ID) %>%
  select(year, ID, votemargin) %>%
  rename(votemargin_stata = votemargin)

RV_compare_R_Stata <- RV_out %>%
  left_join(data_rv_stata, by = c("ID", "year")) %>%
  mutate(votemargin_diff_rstata = votemargin - votemargin_stata) %>%
  select(year, canton, ID, name, firstname, elected, votemargin, votemargin_stata, votemargin_diff_rstata, partymargin, alliance, suballiance, alliance_dummy, suballiance_dummy) %>%
  as_tibble()

RV_compare_R_Stata_diff <- RV_compare_R_Stata %>% filter(votemargin_diff_rstata > 0.00001)

options("encoding" = "UTF-8")
write.table(RV_compare_R_Stata_diff, file = "./02_Processed_data/13_Running_variable/RV_compare_R_Stata_diff.csv", sep = ";")
