
#-------------------------------------------------------
# Functions for vote margin for highest average methods
#-------------------------------------------------------

#------------------------------------------------------
# (A) Load packages, settings, and set working directory
#-------------------------------------------------------

library(tidyverse, quietly = TRUE)
options(dplyr.summarise.inform = FALSE)
options(warn = 1)
options("encoding" = "UTF-8")
library(furrr, quietly = TRUE)


#-------------------------------------------------------------------------------
# (B) Prepare data and collapse votes
# 
# Description: PrepareData prepares the data for the later functions 
#              GetAllianceSuballiance gets information on alliances/suballiances
#              CollapseVotes collapses party votes from the  candidate dataset
#              CollapseVotesCand collapses candidate votes from the  candidate dataset
#              CalculateRatios calculates Ratios (dHondt numbers) in long format
#-------------------------------------------------------------------------------

PrepareData <- function(data,
                        seat_name,
                        votes_j_name,
                        votes_h_name=votes_j_name, 
                        alliance_name=party_name, 
                        suballiance_name=party_name, 
                        party_name, 
                        districtname, 
                        election_cyclename = "year",
                        system="open",
                        cand_id_name,
                        alliances=FALSE) {
  # -----------------------------------------------------------------
  # This function translates the variable names from a dataset to the 
  # corresponding variables names used in the functions
  # system:     open or closed list pr system
  # alliances:  are alliances allowed (TRUE) or not (FALSE)
  # -----------------------------------------------------------------
  if (length(data[[party_name]]) == 0) stop("Error: Party variable not found. ")
  if (is.numeric(data[[votes_h_name]]) == F) 
    stop("Error: Variable for individual votes (votes_h) 
         contains non-numeric values.")
  if (is.numeric(data[[votes_j_name]]) == F) 
    stop("Error: Variable for party votes (votes_j) contains 
         non-numeric values. ")
  data$year <- data[[election_cyclename]]
  data$district <- data[[districtname]]
  data$votes_j <- data[[votes_j_name]]
  data$party_orig <- paste(as.character(data$district), 
                      as.character(data$year), 
                      data[[party_name]])
  if (alliances==TRUE){ # dummy whether alliances are possible
    if (length(data[[alliance_name]]) == 0) 
      stop("Error: Alliance variable not found. ")
    if (length(data[[suballiance_name]]) == 0) 
      stop("Error: Suballiance variable not found. ")
 data$alliance_orig <- as.character(data[[alliance_name]]) 
    # note: 
 data$suballiance_orig <- as.character(data[[suballiance_name]]) 
 data <- data %>% mutate(alliance=ifelse(alliance_orig%in%c("","."),
                                          party_orig,alliance_orig),
                         suballiance=ifelse(suballiance_orig%in%c("","."),
                                            party_orig,suballiance_orig),
                         alliance_dummy=ifelse(alliance_orig%in%c("","."),
                                               0,1),
                         suballiance_dummy=ifelse(suballiance_orig%in%c("","."),
                                               0,1))
  data$party <- as.numeric(as.factor(as.character(data$party_orig)))
  data$alliance <- as.numeric(as.factor(as.character(data$alliance)))
  data$suballiance <- as.numeric(as.factor(as.character(data$suballiance)))
  
  # note: write party name into free alliance/suballiance column for parties 
  #       that are not in alliance/suballiance; make alliance, suballiance and 
  #       party variable numeric and generate  dummy for whether a party is in 
  #       an alliance and on for suballiance
  }
  if (alliances==FALSE){
  data$party <- as.numeric(as.factor(as.character(data$party_orig)))
  data$alliance <- as.numeric(as.factor(as.character(data$party_orig)))
  data$suballiance <- as.numeric(as.factor(as.character(data$party_orig)))
  data$alliance_dummy <- 0
  data$suballiance_dummy <- 0
  }
    
  
  if (system=="open"){ # for open list pr systems
  data$votes_h <- data[[votes_h_name]]
  data$id_cand <- as.numeric(as.factor(paste0(data$year,data$district,data$party,data$alliance,data[[cand_id_name]])))
  }
  
  
  if (system=="closed"){ # for open list pr systems
    data$seat <- data[[seat_name]]
  }

  dataout <- data %>%
    select(-ends_with("_orig"))
  return(dataout)
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
  return(data.frame(alliance_dummy = alliance_dummy, 
                    suballiance_dummy = suballiance_dummy, 
                    alliance = alliance, suballiance = suballiance))
}

CollapseVotes <- function(data_input, aggregation_level, aggregation_unit = NA){ 
  aggregation_level <- enquo(aggregation_level)
  aggregation_unit <- enquo(aggregation_unit)
  check_aggregation_level <- names(data_input %>% ungroup() %>% 
                                    select(!!aggregation_level))
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


CollapseVotesCand <-  function(data_input, aggregation_level, aggregation_unit = NA){ 
  aggregation_level <- enquo(aggregation_level)
  aggregation_unit <- enquo(aggregation_unit)
  check_aggregation_level <- names(data_input %>% ungroup() %>% 
                                     select(!!aggregation_level))
  if (check_aggregation_level == "suballiance") {
    data_input <- data_input %>% filter(alliance == !!aggregation_unit)
  }
  if (check_aggregation_level == "party") {
    data_input <- data_input %>% filter(suballiance == !!aggregation_unit)
  }
  data_out <- data_input %>%
    group_by(!!aggregation_level) %>%
    summarize(votes_h = sum(votes_h)) %>%
    as.matrix()
  return(data_out)
}

# 
# CalculateDHondtNumbers <- function(votes, n) {
#   votes <- as.matrix(votes)
#   div <- matrix(rep(c(1:n), length(votes[, 2])), 
#                 nrow = length(votes[, 2]), byrow = T)
#   DHondt_wide <- as.data.frame(cbind(votes[, 1], 
#                   matrix(rep(as.numeric(votes[, 2]), n),
#                   nrow = length(votes[, 2])) / div))
#   colnames(DHondt_wide) <- c("party", paste(1:n, sep = ""))
#   DHondt_long <- gather(DHondt_wide, paste(1:n, sep = ""), 
#                   key = "seat", value = "No")
#   return(DHondt_long)
# }
# 
# CalculateSainteLagueNumbers <- function(votes, n) {
#   votes <- as.matrix(votes)
#   SainteLague_wide <- as.data.frame(cbind(votes[, 1], 
#                       matrix(rep(as.numeric(votes[, 2]), n), 
#                       nrow = length(votes[, 2])) / div))
#   colnames(SainteLague_wide) <- c("party", paste(1:n, sep = ""))
#   SainteLague_long <- gather(SainteLague_wide, paste(1:n, sep = ""), 
#                             key = "seat", value = "No")
#   return(SainteLague_long)
# }

CalculateRatios <- function(votes, n,method="dHondt") {
  votes <- as.matrix(votes)
  if (method=="dHondt"){
  div <- matrix(rep(c(1:(n)), length(votes[, 2])), 
                nrow = length(votes[, 2]), byrow = T)
  }
  if (method=="SainteLague"){
  div <- matrix(rep(c(1.4,seq(3,n*2,2)), 
                    length(votes[,2])), 
                nrow = length(votes[,2]), 
                byrow = T)
  }
  Votes_wide <- as.data.frame(cbind(votes[, 1], 
                                    matrix(rep(as.numeric(votes[, 2]), 
                                               n), 
                                           nrow = length(votes[, 2])) / div))
  colnames(Votes_wide) <- c("party", paste(1:n, sep = ""))
  Votes_long <- gather(Votes_wide, paste(1:n, sep = ""), 
                       key = "seat", 
                       value = "No")
  return(Votes_long)
}

AppendDataFrames <- function(data_input,parties){
  df_out <- data.frame()
  for (j in parties){
    df_out <- bind_rows(df_out,data_input %>% mutate(party=j)) 
  }
  return(df_out)
}

#-------------------------------------------------------------------------------
# (C) Calculate party margin
# 
# Description:  PartyMargin calculates the party margin for a single alliance/
#                party/suballiance
#               CalculatePartyMargin calculates the party margin for alliance/
#                party/suballiance and aggregates this margin
#               CalculatePartyMargins calculates party margin for different 
#                parties
#-------------------------------------------------------------------------------

PartyMargin <- function(party, data_input, n,method="dHondt") {

  party <- enquo(party)

  if (method=="dHondt"){
    df_mult <- data.frame(mult=c(1:n),
                          seat=c(1:n))
  }
  if (method=="SainteLague"){
    df_mult <- data.frame(mult=c(1.4,seq(3,n*2,2)),
                          seat=c(1:n))
  }
  
  data_input <- data_input %>% mutate(seat=as.numeric(seat),
                                      No=as.numeric(No))
  

  data_out <- data_input %>%
    filter(party == !!party) %>%
    mutate(seat_merge= as.numeric(seat)) %>%  
    left_join(df_mult, by = c("seat"))

  data_others <- data_input %>%
    filter(party != !!party) %>%
    mutate(No_rank = dense_rank(-as.numeric(No))) %>%
    filter(No_rank <= n) %>%
    mutate(seat_merge = n - No_rank + 1) %>%
    select(No, seat_merge,seat) %>%  
    left_join(df_mult, by = c("seat"))

    data_out <- data_out  %>%
      left_join(data_others, by = c("seat_merge")) %>%
      mutate(mult=mult.x) %>%
      mutate(margin = (as.numeric(No.x) - as.numeric(No.y)) * mult) %>%
      rename(seat=seat.x)%>%
      select(party, seat, margin)

  return(data_out)
}

CalculatePartyMargin <- function(party_input, 
                                 data_input, 
                                 n,
                                 method="dHondt",
                                 treshold="none",
                                 rank="rank") {
  #----------------------------------------------------------------------------
  # This function calculates the party margin aggregated (if necessary) over 
  # party, suballiance, alliance
  # Input: party_input: choice of party
  #        data_input:  votes data with cols party, alliance, suballiance, votes
  #        n:           number of seats in a district
  #        method:      allocation method (dHondt,SainteLague)
  #-----------------------------------------------------------------------------
  #print(paste("Party:", party_input))
  AllianceInformation <- GetAllianceSuballiance(party_input, data_input) 
  # note: get info whether party is in sub-/alliance and in which
  #print(paste("Alliance information:", AllianceInformation))
  # 1. Calculation of margin across alliances
  VotesAlliance <- CollapseVotes(data_input, aggregation_level = alliance)
  AllianceNumbers <- CalculateRatios(votes=VotesAlliance, 
                                     n=n,
                                     method=method)
  Alliancemargin <- PartyMargin(party=AllianceInformation$alliance, 
                                data_input=AllianceNumbers, 
                                n=n,
                                method=method)
  colnames(Alliancemargin) <- c("alliance", "seat_a", "margin_a")
  
  # 2. Calculation of margin across suballiances (only for parties in alliance)
  if (AllianceInformation$alliance_dummy == 1) {
    Suballiancemargin_all <- data.frame()
    for (i in c(1:n)) {
      VotesSuballiance <- CollapseVotes(data_input, 
                                aggregation_level = suballiance, 
                                aggregation_unit = AllianceInformation$alliance)
      SuballianceNumbers <- CalculateRatios(votes=VotesSuballiance, 
                                            i,
                                            method=method)
      Suballiancemargin <- PartyMargin(party=AllianceInformation$suballiance,
                                       data_input= SuballianceNumbers, 
                                       n=i,
                                       method=method)
      Suballiancemargin$seat_a <- i
      colnames(Suballiancemargin) <- c("suballiance", "seat_s", 
                                       "margin_s", "seat_a")
      Suballiancemargin_all <- rbind(Suballiancemargin_all, Suballiancemargin)
    }
    # 3. Calculation of margin across parties (only for parties in suballiance)
    if (AllianceInformation$alliance_dummy == 1 & 
        AllianceInformation$suballiance_dummy == 1) {
      Partymargin_all <- data.frame()
      for (i in c(1:n)) {
        VotesParty <- CollapseVotes(data_input, 
                             aggregation_level = party, 
                             aggregation_unit = AllianceInformation$suballiance)
        PartyNumbers <- CalculateRatios(votes=VotesParty, i,method=method)
        Partymargin <- PartyMargin(party=party_input, 
                                   data_input=PartyNumbers, 
                                   n=i,
                                   method=method)
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
  if (AllianceInformation$alliance_dummy == 1 & 
      AllianceInformation$suballiance_dummy == 0) {
    Allpartymargins <- Suballiancemargin_all %>% 
                       left_join(Alliancemargin, by = c("seat_a"))
    Allpartymargins <- Allpartymargins %>% 
                       mutate(margin_min = pmin(margin_a, margin_s))
    Allpartymargins$party <- party_input
    AllpartymarginsOut <- Allpartymargins %>% filter(seat_s <= seat_a) 
    # note: keep only seats that are possible to make for suballiance
    AllpartymarginsOut <- AllpartymarginsOut %>%
      group_by(seat_s) %>%
      summarize(partymargin = max(margin_min), party = mean(party))
    AllpartymarginsOut <- AllpartymarginsOut %>% rename(seat = seat_s)
  }
  # Case 3: Alliance and suballiance
  if (AllianceInformation$alliance_dummy == 1 & 
      AllianceInformation$suballiance_dummy == 1) {
    Allpartymargins <- Partymargin_all %>%
      left_join(Suballiancemargin_all, by = c("seat_s")) %>%
      left_join(Alliancemargin, by = c("seat_a"))
    Allpartymargins <- Allpartymargins %>% 
                       mutate(margin_min = pmin(margin_a, margin_s, margin_p))
    Allpartymargins$party <- party_input
    AllpartymarginsOut <- Allpartymargins %>% 
                          filter(seat_p <= seat_s & seat_s <= seat_a) 
    # note: keep only seats that are possible to make for suballiance and party
    AllpartymarginsOut <- AllpartymarginsOut %>%
      group_by(seat_p) %>%
      summarize(partymargin = max(margin_min), party = mean(party))
    AllpartymarginsOut <- AllpartymarginsOut %>% rename(seat = seat_p)
  }
  return(AllpartymarginsOut)
}

CalculatePartyMargins <- function(data_input, n,method="dHondt") {
  Out <- data.frame()
  all_parties <- as.numeric(levels(as.factor(data_input$party)))
  for (j in 1:length(all_parties)) {
    #no_candidates <- data_input[as.factor(data_input$party)==all_parties[j]]
    Out <- rbind(Out, CalculatePartyMargin(party=all_parties[j], data_input,
                                           n=n,method=method))
  }
  return(Out)
}

#-------------------------------------------------------------------------------
# (D) Calculate candidate margin
# 
# Description:  CalculateCandidateMargin calculates the candidate margin for all
#                candidates in a party
#               CalculateCandidateMargins calculates candidate margin for 
#                different parties
#-------------------------------------------------------------------------------

CalculateCandidateMargin <- function(party_input, data_input) {
  data <- data_input[data_input$party == party_input, ]
  data <- data %>%
    mutate(rank_h_inv = dense_rank(as.numeric(votes_h+elected+
                                              row_number()/1000))) %>% 
    # note: + elected is for coin flip cases where candidates have the same 
    #         number of votes
    #       + row_number()/1000 is for cases in which both candidates have the 
    #         same elected status (e.g., Liste 12 FDP 1967 in ZH)
    mutate(rank_h_inv_max = max(rank_h_inv)) %>%
    mutate(rank_h = rank_h_inv_max - rank_h_inv + 1)
  out <- data.frame()
  out_all <- data.frame()
  if (dim(data)[1] > 1) { 
    # calculate candidate margin only for parties with more than one candidate
    candidates <- max(data$rank_h)
    # print(paste("Party:",party_input))
    for (h in 1:candidates) {
      #print(h)
      out_h <- data.frame()
      for (hcom in 1:(candidates + 1)) {
        if (h != hcom & hcom < (candidates + 1)) {
          candmargin <- min(data$votes_h[data$rank_h == h]) - 
                        min(data$votes_h[data$rank_h == hcom]) 
          # note: min is for cases with more than one canidate 
          out <- data.frame(candmargin = candmargin, rank_h = h, 
                            rank_h_compare = hcom)
          out_h <- rbind(out_h, out)
        }
        if (h != hcom & hcom == (candidates + 1)) {
          out <- data.frame(candmargin = data$votes_h[data$rank_h == h], 
                            rank_h = h, 
                            rank_h_compare = hcom)
          out_h <- rbind(out_h, out)
        }
      }
      out_all <- rbind(out_h, out_all)
    }
  } else {
    out_all <- data.frame(candmargin = data$votes_h, 
                          rank_h = 1, 
                          rank_h_compare = 1)
  } # candidate margin is missing for parties with only one candidate
  out_all$party <- party_input
  return(out_all)
}

CalculateCandidateMargins <- function(data_input) {
  party_inputs <- levels(as.factor(data_input$party))
  ncand <- dim(data_input)[1]
  out_all <- data.frame()
  for (j in 1:length(party_inputs)) {
    out_all <- rbind(CalculateCandidateMargin(data_input = data_input, 
                                              party_input = party_inputs[j]), 
                                              out_all)
  }
  return(out_all)
}



#-------------------------------------------------------------------------------
# (D) Aggregate candidate and party margin
# 
# Description:  AggregateMargins aggregates the party and the candidate margin  
#                for all candidates in a party
#-------------------------------------------------------------------------------

AggregateMargins <- function(data_input, 
                             n,
                             system="open",
                             method="dHondt",
                             rank="rank",
                             additional_vars="",
                             threshold="none"){
  #--------------------------------------------------------------------------
  # This function aggregates party and candidate margin for open-list systems
  #--------------------------------------------------------------------------

  
  if (threshold=="spain_2004_2019"){
    data_input <- data_input %>% 
      mutate(p_votesRelative=p_votes/(p_votesPro_other + p_votes + blank_votes)) %>%
      filter(p_votesRelative>=0.03)
  }
  
  PartyMarginOut <- CalculatePartyMargins(data_input=data_input, 
                                          n=n,method=method)%>% 
                    mutate(mergevar = seat,party=as.character(party))
    
  
  if (system=="open"){
  CandidateMarginsOut <- CalculateCandidateMargins(data_input) %>% 
                         mutate(mergevar = ifelse(rank_h < rank_h_compare, 
                                           rank_h_compare - 1, rank_h_compare),
                                           party=as.character(party))
  MarginsMerged <- CandidateMarginsOut %>%
    left_join(PartyMarginOut, by = c("party", "mergevar")) %>%
    mutate(margin_min = ifelse(is.na(candmargin) == F, 
                               pmin(candmargin, partymargin),
                               partymargin)) 
  # note: minimum of candidate and partymargin for all canidates for whom both 
  #       are available, partymargin for those with no candidate margin 
  #       (single-candidate lists)
  MarginsFinal <- MarginsMerged %>%
    group_by(party, rank_h) %>%
    summarize(votemargin = max(margin_min,na.rm=T),candmargin=min(candmargin,na.rm=T)) 
  # note: maximum of all potential seat configurations
  MarginsFinal <- MarginsFinal %>%
           mutate(votemargin=ifelse(is.na(votemargin)==T,candmargin,votemargin)) 
  # note: replace votemargin as candmargin for all single candidates in which on  
  #       other list exists (e.g., D?hler	Edmund, AI 1931)
  data_input <- data_input %>%
    mutate(party=as.character(party))%>%
    group_by(party,year) %>%
    mutate(rank_h_inv = 
            dense_rank(as.numeric(votes_h+elected+row_number()/1000))) %>%
    mutate(rank_h_inv_max = max(rank_h_inv)) %>%
    mutate(rank_h = rank_h_inv_max - rank_h_inv + 1)
  }
  
  if (system=="closed"){
   MarginsFinal <-PartyMarginOut %>% rename(rank_h=mergevar,
                                            votemargin=partymargin)
   data_input <- data_input %>%
      mutate(party=as.character(party))
   data_input$rank_h <- data_input[[rank]]
  }

  return(MarginsFinal %>% left_join(data_input, by = c("party", "rank_h")) %>% 
             select(year, party,votemargin,all_of(additional_vars))) 
}

#-------------------------------------------------------------------------------
# (E) Calculate vote margin for several years, districts, and parties
# 
# Description:  CalculateMargins calculates the vote margin for several  
#                years, districts, and parties
#-------------------------------------------------------------------------------



CalculateMargins <- function(data_input,
                             system="open",
                             method="dHondt",
                             rank="rank", 
                             convcrit=0.001,
                             calculate_all_seats=TRUE,
                             total_seats_name="none",
                             additional_vars=c("alliance"),
                             threshold="none",
                             multiproc=TRUE,
                             rv_type="analytical",
                             return_option=T,
                             outfile_name,
                             print_seat=FALSE,
                             print_party=FALSE) {
  
  
  # 1) Check for single candidates
  
  data_singlecandidate <- data_input %>%
    group_by(year, district) %>%
    summarize(nobs = n()) %>%
    filter(nobs == 1)
 
  # 2) Loop over years and districts
  
  years <- as.numeric(levels(as.factor(data_input$year)))
  Out <- data.frame()
  for (yr in years) {
    print(paste("New year:", yr))
    dataperyear <- data_input[data_input$year == yr, ]
    districts <- levels(as.factor(dataperyear$district))
    for (di in districts) {
      print(paste("New district:", di))
      dataworking <- dataperyear[dataperyear$district == di, ]
      
      
      # 3) Define number of seats as sum of elected candidates if not provided
      #    by variable total_seats_name (string)
      
      if (total_seats_name=="none"){
      n <- sum(dataworking$elected)
      } else{
      n <- dataworking[[total_seats_name]][1]  
      }
      
      # 4) Calculate analytical vote margins as in Luechinger, Schelker, Schmid (2024, PA)
      
      if (rv_type=="analytical" & method %in% c("dHondt","SainteLague")){
      Out <- rbind(Out, as.data.frame(AggregateMargins(dataworking, 
                                              n=n,
                                              system=system,
                                              method=method,
                                              rank=rank,
                                              threshold=threshold,
                                              additional_vars=additional_vars)))
      }
      
      if (rv_type=="analytical" & method %in% c("Hare")){
        Out <- rbind(Out, as.data.frame(AggregateMarginsLargestRemainder(data_input=dataworking, 
                                                                         no_seats=n,
                                                                         system=system,
                                                                         method=method,
                                                                         rank=rank,
                                                                         threshold=threshold,
                                                                         additional_vars=additional_vars)))
      }
      
      
      # 5) Calculate simulated vote margins 
      
      if (rv_type=="simulation"){
        if (multiproc==TRUE){
        out_temp <- as.data.frame(GetRVSimulation_Candidate_MP(data_input=dataworking, 
                                                              system=system,
                                                              n=n,
                                                              calculate_all_seats=calculate_all_seats,
                                                              method=method,
                                                              convcrit=convcrit,
                                                              threshold=threshold,
                                                              additional_vars=additional_vars,
                                                              print_seat=print_seat,
                                                              print_party=print_party))
          
          
        } else{
        out_temp <- as.data.frame(GetRVSimulation_Candidate(data_input=dataworking, 
                                                                 system=system,
                                                                 n=n,
                                                                 calculate_all_seats=calculate_all_seats,
                                                                 method=method,
                                                                 convcrit=convcrit,
                                                                 threshold=threshold,
                                                                 print_seat=print_seat,
                                                                 print_party=print_party))
        }

        if (return_option==T){Out <- rbind(Out,out_temp)
        } else{
        data.table::fwrite(out_temp, file = outfile_name, 
                           sep = ";",row.names = F,append=T)  
        }
      }
    }
  }
  if (return_option==T){return(Out)}
  }

#-------------------------------------------------------------------------------
# (F) Calculate seat allocation
# 
# Description:  GetSeatsHighestAverage calculates the number of seats for a  
#                single party
#               GetSeatsHighestAverageLoop calculates the number of seats for   
#                several years, districts, and parties 
#-------------------------------------------------------------------------------


CutDataByQuorum <- function(data_input, threshold){
if (threshold=="spain_2004_2019"){
    
    blankv <- data_input$blank_votes[1]
    p_votes_total <- data_input %>%
      filter(list_pos==1) %>%
      select(votes_j) %>%
      summarize(p_votes_total=sum(votes_j)) %>%
      pull()
    
    data_input_out <- data_input %>%
      mutate(votes_j_rel=votes_j/(p_votes_total+blankv))%>% 
      filter(votes_j_rel>=0.03) %>%
      select(-votes_j_rel)
    

}
if (threshold=="israel_2009_2013"){
  
  p_votes_total <- data_input %>%
    filter(seat==1) %>%
    select(votes_j) %>%
    summarize(p_votes_total=sum(votes_j)) %>%
    pull()
  
  data_input_out <- data_input %>%
    mutate(votes_j_rel=votes_j/p_votes_total)%>% 
    filter(votes_j_rel>=0.02) %>%
    mutate(total_valid_votes=sum(votes_j))%>%
    select(-votes_j_rel)
  
  
}
  
if (threshold=="israel_2014_2022"){
  
  
  p_votes_total <- data_input %>%
    filter(seat==1) %>%
    select(votes_j) %>%
    summarize(p_votes_total=sum(votes_j)) %>%
    pull()
    
    data_input_out <- data_input %>%
      mutate(votes_j_rel=votes_j/p_votes_total)%>% 
      filter(votes_j_rel>=0.0325) %>%
      mutate(total_valid_votes=sum(votes_j))%>%
      select(-votes_j_rel)
    
}
    
    return(data_input_out)
  
}


GetSeatsHighestAverage <- function(party_input, 
                                   data_input, 
                                   n,
                                   method="dHondt",
                                   threshold="none") {
  #----------------------------------------------------------------------------
  # This function calculates the number of seats for a party using highest 
  # average methods
  # Input: party_input: choice of party
  #        data_input:  votes data with columns party, alliance, suballiance, 
  #                     votes
  #        n:           number of seats in a district
  #        method:      allocation method (dHondt,SainteLague)
  #        threshold:   threshold for parties to be considered
  #-----------------------------------------------------------------------------
  #n <- sum(data_input$elected)
  
  quorum=TRUE # set quorum to true as default
  
  if (threshold!="none"){
  data_input <- CutDataByQuorum(data_input,threshold)
  quorum <- is.element(party_input,levels(as.factor(data_input$party)))
  
  if (quorum==FALSE){ out <- 0  # calculation only for parties who made quorum
  }
  }
  
  if (quorum!=FALSE | threshold=="none"){
  
  AllianceInformation <- GetAllianceSuballiance(party_input, data_input) 
  
  # 1. Calculation of number of seats across alliances
  
  VotesAlliance <- CollapseVotes(data_input, aggregation_level = alliance)
  AllianceSeats <- CalculateRatios(votes=VotesAlliance, 
                                   n=n,
                                   method=method) %>% 
    mutate(No_rank = row_number(-as.numeric(No))) %>%
    filter(No_rank <= n) %>%
    group_by(party) %>%
    summarize(seats_alliance=n())%>%
    rename(alliance=party)
  
  out <- max(AllianceSeats$seats_alliance[AllianceSeats$alliance==
                                            AllianceInformation$alliance],0)

  # 2. Calculation of suballiance seats
  if (AllianceInformation$alliance_dummy == 1 & out>0) {
    VotesSuballiance <- CollapseVotes(data_input, 
                                aggregation_level = suballiance, 
                                aggregation_unit = AllianceInformation$alliance)
    SuballianceSeats <- CalculateRatios(votes=VotesSuballiance, 
                         n= AllianceSeats$seats_alliance[AllianceSeats$alliance==AllianceInformation$alliance],
                         method=method) %>% 
      mutate(No_rank = row_number(-as.numeric(No))) %>%
      filter(No_rank <= AllianceSeats$seats_alliance[AllianceSeats$alliance==AllianceInformation$alliance]) %>%
      group_by(party) %>%
      summarize(seats_suballiance=n())%>%
      rename(suballiance=party)
    
    out <- max(SuballianceSeats$seats_suballiance[SuballianceSeats$suballiance==
                                                 AllianceInformation$suballiance],0)

    # 3. Calculation of margin across parties (only for parties in suballiance)
    if (AllianceInformation$alliance_dummy == 1 & 
        AllianceInformation$suballiance_dummy == 1 & out>0) {
      VotesParty <- CollapseVotes(data_input, aggregation_level = party, 
                             aggregation_unit = AllianceInformation$suballiance)
      PartySeats <- CalculateRatios(votes=VotesParty, 
                                    n=SuballianceSeats$seats_suballiance[SuballianceSeats$suballiance==AllianceInformation$suballiance],method=method) %>% 
        mutate(No_rank = row_number(-as.numeric(No))) %>%
        filter(No_rank <= SuballianceSeats$seats_suballiance[SuballianceSeats$suballiance==AllianceInformation$suballiance]) %>%
        group_by(party) %>%
        summarize(seats_party=n())
      out <- max(PartySeats$seats_party[PartySeats$party==party_input],0)  
    }
  }
  }
  return(out)
}

GetSeatsHighestAverageLoop<- function(data_input,method="dHondt"){
data_input$seats_check <- NA
df_out <- data.frame()
years_all <- levels(as.factor(as.character(data_input$year)))   
print(paste0("Years: ", paste0(years_all,collapse=", ")))
for (yr in years_all) {
df1 <- data_input[data_input$year==yr,]
districts_all <- levels(as.factor(as.character(df1$district)))
#print(districts_all)
for (district in districts_all) {
df2 <- df1[df1$district==district,]
parties_all <- levels(as.factor(as.character(df2$party)))
#print(parties_all)
for (party in parties_all) {
df_out <- rbind(df_out,data.frame(year=yr,
                       districtid=district,
                       party=party,
                       seats=GetSeatsHighestAverage(party_input=party, 
                                                    data_input=df2, 
                                                    n=sum(df2$elected),
                                                    method=method)))
}  
}
}
return(df_out)
}




GetSeatsIsrael <- function(party_input,
                           data_input,
                           n,
                           method,
                           threshold="none") {
  #----------------------------------------------------------------------------
  # This function calculates the number of seats for a party using the seat 
  # distribution method for Israel 2009-2023. 
  # Input: party_input: choice of party
  #        data_input:  votes data with columns party, alliance, suballiance, 
  #                     votes
  #        n:           number of seats in a district
  #        method:      allocation method
  #        threshold:   threshold for parties to be considered
  #-----------------------------------------------------------------------------
  #n <- sum(data_input$elected)
  
  quorum <- TRUE # set quorum to true as default

  
  
  # 1. Check whether party is above threshold
  
  if (threshold!="none"){
    data_input <- CutDataByQuorum(data_input,threshold)
    quorum <- is.element(party_input,levels(as.factor(data_input$party)))
    
    if (quorum==FALSE){ out <- 0  # calculation only for parties who made quorum
    }
  }
  
  if (quorum==TRUE | threshold=="none"){
    
    AllianceInformation <- GetAllianceSuballiance(party_input, data_input) 
    
    # 2. Calculation of number of seats according to first round of Hare system
    
    VotesParty <- CollapseVotes(data_input, aggregation_level = list_name) %>%
      as_tibble() %>%
      mutate(votes_j=as.numeric(votes_j))
    
    df_party_names <- data_input %>% distinct(list_name,party)
    
    Seats_Round1 <- GetIntegerHare(VotesParty$votes_j,sum(data_input$seats_j[data_input$seat==1]))
    
    df_seats_round1 <- tibble(seats_round1=as.numeric(Seats_Round1),
                              list_name=VotesParty$list_name)
    
    RemSeats <- as.numeric(sum(data_input$seats_j[data_input$seat==1])-sum(Seats_Round1))
    
    
    # 3. Distribute seats to alliances using dHondt method
    
    data_input <- data_input %>% 
      left_join(df_seats_round1,by=c("list_name")) %>%
      mutate(seats_round1=ifelse(is.na(seats_round1),0,seats_round1))
    
    data_input_seats <- data_input %>% 
      group_by(list_name) %>%
      summarize(seats_round1=first(seats_round1),
                alliance=first(alliance)) %>% 
      group_by(alliance) %>%
      summarize(seats_round1=sum(seats_round1))
      
    
    VotesAlliance <- CollapseVotes(data_input, aggregation_level = alliance) %>% # continue here 
      as_tibble() %>%
      mutate(votes_l=as.numeric(votes_j)) %>%
      left_join(data_input_seats,by=c("alliance")) %>%
      mutate(HA_seat_l=0) %>%
      select(-votes_j)
    
    for (s in 1:RemSeats){
    VotesAlliance <- VotesAlliance %>%
        mutate(indic_l=votes_l/(seats_round1+1+HA_seat_l))%>%
        mutate(rank = dense_rank(as.numeric(-indic_l))) %>%
        mutate(HA_seat_l=ifelse(rank==1,HA_seat_l+1,HA_seat_l))

    }
    
    
    
    # 4. Distribute alliances seats from step 2 to parties using DHondt method
    
    data_input_seats <- data_input %>% 
      group_by(list_name) %>%
      summarize(seats_round1=first(seats_round1),
                alliance=first(alliance),
                votes_j=first(votes_j)) %>%
      left_join(VotesAlliance %>% select(HA_seat_l,alliance),
                by=c("alliance"))%>%
      mutate(HA_seat_j=0) %>%
      ungroup() %>%
      arrange(alliance,list_name)
    
    HA_seat_l_max <- max(data_input_seats$HA_seat_l)
    
    for (s in 1:HA_seat_l_max){
      data_input_seats <- data_input_seats %>%
        mutate(indic_j=votes_j/(seats_round1+1+HA_seat_j))%>%
        group_by(alliance) %>%
        mutate(rank = dense_rank(as.numeric(-indic_j))) %>%
        ungroup() %>%
        mutate(HA_seat_j=ifelse(rank==1 & HA_seat_l>=s,HA_seat_j+1,HA_seat_j)) 
      
    }
    
    out <- data_input_seats %>% 
      mutate(seats=seats_round1+HA_seat_j) %>% 
      select(list_name,seats) %>%
      left_join(df_party_names,by=c("list_name"))
    out <- out %>% filter(party==party_input) %>% select(seats) %>% pull()
  }   
  
  return(out)
}



GetElectionStatusBrazil<- function(data_input,
                                   n,
                                   method,
                                   threshold="none",
                                   precision=0.000000000000001) {
  #----------------------------------------------------------------------------
  # This function calculates the election status for candidates in a district
  # in a specific year for Brazil 2006 until 2018
  # Input: data_input:  votes data with columns party, alliance, suballiance, 
  #                     votes
  #        n:           number of seats in a district
  #        method:      allocation method
  #        threshold:   threshold for parties to be considered
  #-----------------------------------------------------------------------------
  #n <- sum(data_input$elected)
  
  # 1) Sum up vote_j and votes_h over all alliances and
  
  votes_j_alliance  <- CollapseVotes(data_input = data_input,aggregation_level=alliance) %>% 
    as_tibble()
  votes_h_alliance  <- CollapseVotesCand(data_input = data_input,aggregation_level=alliance)%>% 
    as_tibble()
  
  QE <- round((sum(votes_j_alliance$votes_j) + sum(votes_h_alliance$votes_h))/n,1)
  
  votes_all_alliance <- votes_j_alliance %>%
    left_join(votes_h_alliance, by=c("alliance")) %>%
    mutate(votes_l=votes_j+votes_h,
           QP=floor(votes_l/QE)) %>%
    mutate(coalition_threshold=ifelse(votes_l >=QE,1,0))%>%
    dplyr::rename(votes_h_alliance=votes_h)%>%
    select(alliance,votes_l,QP,votes_h_alliance,coalition_threshold)%>%
    ungroup()%>%
    mutate(total_seats=QP)
  
  
  # 2) Distribution Step 1: Assign seats to alliances based on largest remainder
  
  year <- as.numeric(as.character(data_input$year))[1]
  
  if (year<2018){
    
    votes_all_alliance <- votes_all_alliance %>%  
      filter(coalition_threshold==1) 
    
    data_input <- data_input %>%
      left_join(votes_all_alliance, by=c("alliance")) %>%
      mutate(coalition_threshold=ifelse(votes_l >=QE,1,0))%>%
      filter(coalition_threshold==1) %>%
      group_by(district,year,alliance)  %>%
      mutate(votes_h_rank= dense_rank(as.numeric(-votes_h))) %>%
      ungroup() %>%
      mutate(elected_sim=ifelse(votes_h_rank<=QP & coalition_threshold==1,1,0),
             eligible_cand=1-elected_sim,
             eligible_coal_threshold=(coalition_threshold*eligible_cand)) %>%
      group_by(district,year,alliance)  %>%
      mutate(eligible_coal=max(eligible_cand),
             coal_threshold=max(coalition_threshold),
             threshold_eligible_coal = max(eligible_coal_threshold)) %>%
      ungroup()%>%
      mutate(total_seats=QP)
  }
  
  if (year==2018){
    data_input <- data_input %>%
      left_join(votes_all_alliance, by=c("alliance")) %>%
      mutate(cand_threshold=ifelse(votes_h >= 0.1*QE,1,0))%>%
      filter(cand_threshold==1) %>%
      group_by(district,year,alliance)  %>%
      mutate(votes_h_rank= dense_rank(as.numeric(-votes_h))) %>%
      ungroup() %>%
      mutate(elected_sim=ifelse(votes_h_rank<=QP & cand_threshold==1,1,0),
             eligible_cand=1-elected_sim,
             eligible_cand_threshold=(cand_threshold*eligible_cand)) %>%
      group_by(district,year,alliance)  %>%
      mutate(eligible_coal=max(eligible_cand),
             coal_threshold=max(cand_threshold),
             threshold_eligible_coal = max(eligible_cand_threshold)) %>%
      ungroup()%>%
      mutate(total_seats=QP)
  }
  
  
  # 3) Distribution Step 2: Assign remaining seats based on highest average method
  
  # 3a) Initialization of eligible coalitions (eligible_state) and remaining seats
  #     (remaining_seats)
  
  eligible_state <- max(data_input$threshold_eligible_coal)
  remaining_seats <- n - sum(data_input$elected_sim)
  
  # 3b) Start of while loop
  
  while (remaining_seats > 0 & eligible_state==1){
    
    #  3c) Distribution according to largest remainder method
    #      Note: the term HA stands for highest average but it is not correctly 
    #            chosen as the distribution follows the largest remainder method. 
    #            However, for consistency with the Stata implementation, we chose
    #            to keep the same names as in the Stata code. 
    
    votes_all_alliance <- votes_all_alliance %>%
      mutate(HA=round(votes_l/(total_seats + 1),precision)) %>%
      mutate(HA_rank= dense_rank(as.numeric(-HA)),
             votes_l_rank= dense_rank(as.numeric(-votes_l)),
             votes_h_rank= dense_rank(as.numeric(-votes_h_alliance)))
    
    votes_all_alliance_cut <- votes_all_alliance %>% 
      filter(HA_rank==1)
    
    # 3d) Alliance tie breakers if two alliances have the same HA
    
    if (dim(votes_all_alliance_cut)[1]>1){
      votes_l_rank_min <- min(votes_all_alliance_cut$votes_l_rank)
      votes_all_alliance_cut <- votes_all_alliance_cut %>% 
        filter(votes_l_rank==votes_l_rank_min)
    }
    
    if (dim(votes_all_alliance_cut)[1]>1){
      votes_h_rank_min <- min(votes_all_alliance_cut$votes_h_rank)
      votes_all_alliance_cut <- votes_all_alliance_cut %>% 
        filter(votes_h_rank==votes_h_rank_min)
    }
    
    # 3e) Update election status of candidates (elected_new)
    
    data_input <- data_input %>%
      left_join(votes_all_alliance_cut %>%
                  select(alliance) %>%
                  mutate(additional_seat=1), by=c("alliance")) %>%
      mutate(elected_new=ifelse(additional_seat==1 & votes_h_rank==total_seats+1,1,0))
    
    # 3f) Candidate tie breakers: older candidate wins
    
    if (dim(data_input %>% filter(elected_new==1))[1]>1){
      birthdate_min <- min(data_input %>% filter(elected_new==1) %>% select(birthdate_enc) %>% pull())
      data_input <- data_input %>%
        mutate(elected_new=ifelse(birthdate_enc==birthdate_min,1,0))
    }
    
    # 3g) Update after seat distribution
    
    # (i) Election status of candidates
    
    if (year<2018){
      
      data_input <- data_input  %>%
        mutate(elected_sim=ifelse(is.na(elected_new),0,elected_sim+elected_new),
               eligible_cand=1-elected_sim,
               eligible_cand_threshold=(coal_threshold*eligible_cand)) %>%
        group_by(district,year,alliance)  %>%
        mutate(eligible_coal=max(eligible_cand),
               coal_threshold=max(coal_threshold,na.rm=T),
               threshold_eligible_coal = max(eligible_cand_threshold,na.rm=T),
               total_seats=sum(elected_sim)) %>%
        ungroup() %>%
        select(-elected_new,-additional_seat)
      
    }
    
    if (year==2018){
      
      data_input <- data_input  %>%
        mutate(elected_sim=ifelse(is.na(elected_new),0,elected_sim+elected_new),
               eligible_cand=1-elected_sim,
               eligible_cand_threshold=(cand_threshold*eligible_cand)) %>%
        group_by(district,year,alliance)  %>%
        mutate(eligible_coal=max(eligible_cand),
               coal_threshold=max(cand_threshold,na.rm=T),
               threshold_eligible_coal = max(eligible_cand_threshold,na.rm=T),
               total_seats=sum(elected_sim)) %>%
        ungroup() %>%
        select(-elected_new,-additional_seat)
      
    }
    
    # (ii) Seats of alliances for new round of while loop
    
    alliance_new_seat <- votes_all_alliance_cut %>% 
      select(alliance) %>% 
      pull()
    
    votes_all_alliance <- votes_all_alliance %>%
      mutate(total_seats=ifelse(alliance==alliance_new_seat,total_seats+1,total_seats))
    
    eligible_state <- max(data_input$threshold_eligible_coal)
    remaining_seats <- n - sum(data_input$elected_sim)
  }
  
  # 4) Return results
  
  return(data_input %>% 
           mutate(elected=elected_sim) %>%
           select(id_cand,elected))
}




CheckSeatAllocation <- function(data_input,
                                n,
                                method,
                                threshold,
                                system,
                                print_year=FALSE,
                                print_district=FALSE){
  # This function checks whether the actual seat distribution (seats_j) is the same
  # as the seat distribution obtained by our algorithm for a specific district
  # in a specific year. 
  #
  #if (method=="SainteLague" | method=="dHondt" ){
  
  # 1. Calculate seat distribution of algorithm 
  
  if (method=="Israel"){
    out_algo  <- lapply(c(1:max(data_input$party)),
                        GetSeatsIsrael,
                        data_input=data_input,
                        n=n,
                        method=method,
                        threshold=threshold)
  }
  
  if (method=="Brazil"){
    
    df_out <- tibble()
    
    year_all <- levels(as.factor(data_input$year))
    district_all <- levels(as.factor(data_input$district))
   for (yr in year_all) {
      if (print_year==TRUE){print(paste0("Year: ",yr))}
      data_yr<- data_input %>% filter(year==yr)
      for (di in district_all) {
        if (print_district==TRUE){print(paste0("District: ",di))}
        data_yr_di <- data_yr %>% filter(district==di)
        df_sim <- GetElectionStatusBrazil(data_input=data_yr_di,
                                                         n=sum(data_yr_di$elected),
                                                         method=method,
                                                         threshold=threshold) %>%
          dplyr::rename(elected_sim=elected) 
        
        df_out <- bind_rows(df_out,data_yr_di %>% 
                              left_join(df_sim,by=c("id_cand")) %>%
                              select(id_cand,year,district,elected_sim,elected))
        
      }
    }
        
  }
  
  # 2. Merge seats distribution of algorithm with real seat distribution and 
  #    return result
  
  # 2a) Closed list systems
  
  if (system=="closed" ){
  
  distribution_algo <- data.frame(party=c(1:max(data_input$party)),
                                  seats_sim=as.numeric(out_algo))
  
  df_out <- data_input %>%
    left_join(distribution_algo, by=c("party")) %>%
    mutate(elected_sim=ifelse(seat<=seats_sim,1,0)) %>%
    select(party,seat,elected,elected_sim)
  }
    
return(df_out %>%
         mutate(elected_diff=elected_sim-elected))
 
}
    
  

#-------------------------------------------------------------------------------
# (G) Simulate vote margin
# 
# Description: GetMinVotesQuorum calculates the minimum number of party votes required
#                   for party j to achieve a seat in systems with a quorum
#              GetRVSimulation_Candidate calculates the vote margin for all
#                parties and candidates in a district using a simulation that 
#                adds votes  until the respective candidate is elected
#              GetRVSimulation_Candidate_Loop: calls GetRVSimulation_Candidate
#-------------------------------------------------------------------------------

GetMinVotesQuorum <- function(threshold,data_input,party_j){
  # This function calculates the minimum number of party votes required
  # for party j to achieve a seat in systems with a quorum
  
  if (threshold=="spain_2004_2019"){
    
    
    blankvotes <- data_input$blank_votes[1]
    p_votes_others <- data_input %>%
      filter(party!=j & list_pos==1) %>%
      select(votes_j) %>%
      summarize(p_votes_others=sum(votes_j)) %>%
      pull()
    return(1/0.97*(blankvotes+p_votes_others)*0.03)
    
  }
}

GetRVSimulation_Candidate <- function(data_input, 
                                      n, 
                                      convcrit=0.001,
                                      method="dHondt",
                                      system="open",
                                      threshold="none",
                                      print_party=FALSE,
                                      print_seat=FALSE,
                                      calculate_all_seats) {
  #---------------------------------------------------------------------------
  # data_input: data after PrepareData with votes_j as party votes and votes_h 
  #             as candidate votes
  # n:          number of seats
  # convcrit:   convergence criterion
  # method:     method to distribute seats
  # system:     close-list or open-list
  # 
  #---------------------------------------------------------------------------
  out <- data.frame()
  if (system=="closed"){
  parties_all <- levels(as.factor(data_input$party))
  for (j in parties_all) {
    if (print_party==TRUE){print(paste0("Party: ", j))}
    if (calculate_all_seats==TRUE){no_seat <- n
    } else{
    no_elected <- dim(data_input %>% filter(party==j & elected==1))[1] 
    no_seat <- max(no_elected,1)
    }
    for (i in 1:no_seat) {
    if (print_seat==TRUE){ print(paste0("Seat: ", i))}
      
      votessim <- GetVotesRequieredClosedList(data_input=data_input,
                                    i=i,
                                    j=j,
                                    n=n,
                                    convcrit=convcrit,
                                    method=method,
                                    threshold=threshold)

      out_ji <- data.frame(party = j, seat = i, 
                      votemargin = round(data_input$votes_j[data_input$party==j][1]
                      - votessim, 4))
      out <- rbind(out, out_ji)
    }
  }
  }
  
  
  if (system=="open"){
    parties_all <- levels(as.factor(data_input$party))
    
    for (j in parties_all) {
      ids_all <- levels(as.factor(data_input$id_cand[data_input$party==j]))
      if (print_party==TRUE){print(paste0("Party:", j))}
      
      for (i in ids_all) {
      #print(paste0("ID: ",i))
        
        votessim <- GetVotesRequieredOpenList(data_input=data_input,
                                      i=i,
                                      j=j,
                                      n=n,
                                      method=method,
                                      threshold=threshold,
                                      convcrit=convcrit)

        
        out_ji <- data.frame(party = j, id_cand = data_input$id_cand[data_input$id_cand==i], year=data_input$year[1],
                             district=data_input$district[1],
                             votemargin = round(data_input$votes_h[data_input$id_cand==i & data_input$party==j]
                                                - votessim, 4))
        out <- rbind(out, out_ji)
      }
    }
  }
  
  return(out)
}


GetRVSimulation_Candidate_MP <- function(data_input, 
                                      n, 
                                      convcrit=0.001,
                                      method="dHondt",
                                      system="open",
                                      threshold="none",
                                      additional_vars,
                                      calculate_all_seats=TRUE,
                                      print_party=FALSE,
                                      print_seat=FALSE) {
  #---------------------------------------------------------------------------
  # data_input: data after PrepareData with votes_j as party votes and votes_h 
  #             as candidate votes
  # n:          number of seats
  # convcrit:   convergence criterion
  # method:     method to distribute seats
  # system:     close-list or open-list
  # 
  #---------------------------------------------------------------------------
  out <- data.frame()
  
  if (system=="closed"){
    
    # 1. Expand dataset to all seats for a specific party
    
    data_input_exp <- expand.grid(party=c(1:max(data_input$party)),
                              seat=c(1:n)) %>%
      mutate(id=row_number()) %>%
      left_join(data_input %>% distinct(party,.keep_all=TRUE) %>% select(-seat) 
                , by=c("party")) %>%
      mutate(id=as.factor(id))%>%
      arrange(id) %>%
      as_tibble() 
    
    ids <- levels(as.factor(data_input_exp$id))
    
    # 2. Run MP calculation of GetVotesRequieredClosedList_MP
    
    no_cores <- availableCores() - 1
    plan(multisession,workers=no_cores)
    out_list <- future_map(ids,
                           GetVotesRequieredClosedList_MP,
                           data_input=data_input_exp,
                           n=n,
                           convcrit=convcrit,
                           method=method,
                           threshold=threshold)
   
    out_df <- as.data.frame(do.call(rbind, out_list)) %>%
      mutate(id=ids) %>%
      dplyr::rename(votes_j_req=V1) %>%
      right_join(data_input_exp, by=c("id")) %>%
      mutate(votemargin=round(votes_j-votes_j,4)) %>%
      select(year,party,seat,votemargin, all_of(additional_vars))

}

  
  if (system=="open"){
    
    # 1. Run MP calculation of GetVotesRequieredOpenList_MP

      ids <- levels(as.factor(data_input$id_cand))
      
      no_cores <- availableCores() - 1
      plan(multisession,workers=no_cores)
      out_list <- future_map(ids,
                             GetVotesRequieredOpenList_MP,
                             data_input=data_input,
                             n=n,
                             convcrit=convcrit,
                             method=method,
                             threshold=threshold)
      
      out_df <- as.data.frame(do.call(rbind, out_list)) %>%
        mutate(id=ids) %>%
        dplyr::rename(votes_j_req=V1) %>%
        right_join(data_input_exp, by=c("id")) %>%
        mutate(votemargin=round(votes_h-votessim,4))
  }

  return(out_df)
}

GetVotesRequieredClosedList <- function(data_input,
                              j,
                              n,
                              i,
                              method,
                              threshold,
                              convcrit=0.001){
  
    votessim <- 0.5  # note: start with 1 vote in first iteration (first line after while loop)
    seats <- 0
    
    while (seats < i) {
      votessim <- votessim * 2
      data_input$votes_j[data_input$party==j] <- 
        votessim
      if (method=="SainteLague" | method=="dHondt" ){
        seats <- GetSeatsHighestAverage(party_input=j,
                                        data_input=data_input, 
                                        n=n,
                                        method=method,
                                        threshold=threshold)
      }
      if (method=="Israel"){
      seats <- GetSeatsIsrael(party_input=j,
                                        data_input=data_input, 
                                        n=n,
                                        method=method,
                                        threshold=threshold)  
        
      #print(seats)
      #print(votessim)
      }
      }
    
    xlow <- votessim / 2
    xhigh <- votessim
    

    while ((xhigh - xlow) > convcrit ) { #& i_new<10
      votessim <- (xlow + xhigh) / 2
      
      data_input$votes_j[data_input$party==j] <- 
        votessim
      if (method=="SainteLague" | method=="dHondt"){
        seats <- GetSeatsHighestAverage(party_input=j,
                                        data_input=data_input, 
                                        n=n,
                                        method=method,
                                        threshold=threshold)
      }
      if (method=="Israel"){
        seats <- GetSeatsIsrael(party_input=j,
                                data_input=data_input, 
                                n=n,
                                method=method,
                                threshold=threshold)  
        
      }
      if (seats < i) {
        xlow <- votessim
      }
      if (seats >= i) {
        xhigh <- votessim
      }
    }
    return(votessim)
}



GetVotesRequieredClosedList_MP <- function(id,
                                           data_input,
                                           n,
                                           method,
                                           threshold,
                                           convcrit=0.001){
  
  votessim <- 0.5  # note: start with 1 vote in first iteration (first line after while loop)
  seats <- 0
  
  j <- data_input$party[data_input$id==id]
  i <- data_input$seat[data_input$id==id]
  
  while (seats < i) {
    votessim <- votessim * 2
    data_input$votes_j[data_input$party==j] <- 
      votessim
    if (method=="SainteLague" | method=="dHondt" ){
      seats <- GetSeatsHighestAverage(party_input=j,
                                      data_input=data_input, 
                                      n=n,
                                      method=method,
                                      threshold=threshold)
    }
    if (method=="Israel"){
      seats <- GetSeatsIsrael(party_input=j,
                              data_input=data_input, 
                              n=n,
                              method=method,
                              threshold=threshold)  
      
      #print(seats)
      #print(votessim)
    }
  }
  
  xlow <- votessim / 2
  xhigh <- votessim
  
  
  while ((xhigh - xlow) > convcrit ) { #& i_new<10
    votessim <- (xlow + xhigh) / 2
    
    data_input$votes_j[data_input$party==j] <- 
      votessim
    if (method=="SainteLague" | method=="dHondt"){
      seats <- GetSeatsHighestAverage(party_input=j,
                                      data_input=data_input, 
                                      n=n,
                                      method=method,
                                      threshold=threshold)
    }
    if (method=="Israel"){
      seats <- GetSeatsIsrael(party_input=j,
                              data_input=data_input, 
                              n=n,
                              method=method,
                              threshold=threshold)  
      
    }
    if (seats < i) {
      xlow <- votessim
    }
    if (seats >= i) {
      xhigh <- votessim
    }
  }
  return(votessim)
}


GetVotesRequieredOpenList <- function(data_input,
                                        i,
                                        j,
                                        n,
                                        method="dHondt",
                                        threshold,
                                        convcrit=0.001){
  
  
  # 1. Define initial values of candidate votes and party votes
  
  votes_h_initial <- data_input$votes_h[data_input$id_cand==i]
  votes_j_initial <- data_input$votes_j[data_input$party==j][1]
  
  
  # 2. Iteration 1: Find interval that changes the election status of candidate i
  
  # a) Define initial values of votes for candidate i (votessim), number of seats
  #    for party j (seats),  the rank of a candidate on their party list (rank_h), 
  #    and the election status (elected)
  
  votessim <- 0.5  # note: start with 1 vote in first iteration (first line after while loop)
  seats <- 0
  rank_h <- length(data_input$votes_j[data_input$party==j])
  elected <- 0
  
  while (elected==0) {
    
  # b) Double candidate votes (votes_h) and party votes (for all countries other
  #    than Brazil) 
  # Note: In Brazil, voters can either vote for a party (changing votes_j) or
  #       a candidate (votes_h). In our simulation, we only change votes_h and
  #       leave votes_j unchanged. 

    votessim <- votessim * 2
    data_input$votes_h[data_input$id_cand==i] <- votessim
    if (method!="Brazil"){
    data_input$votes_j[data_input$party==j] <- votes_j_initial-votes_h_initial+votessim
    }
    
    # c) Calculate a party's number of seats (for all countries other than Brazil)
    #    or election status of a candidate i (for Brazil)
    
    if (method=="SainteLague" | method=="dHondt"){
      seats <- GetSeatsHighestAverage(party_input=j,
                                      data_input=data_input, 
                                      n=n,
                                      method=method,
                                      threshold=threshold)
      
    }
    
    if (method=="Hare"){ #continue here for Honduras
      seats <- GetSeatsHare(party_input=j,
                            data_input=data_input,
                            n=n,
                            method=method,
                            threshold=threshold)
      
    }
    
    if (method=="Brazil"){ 
      
      elected <- GetElectionStatusBrazil(data_input=data_input,
                              n=sum(data_input$elected),
                              method=method,
                              threshold=threshold)%>% 
        arrange(-elected) %>%
        filter(id_cand==i) %>%
        select(elected) %>%
        pull()
      
      elected <- ifelse(length(elected)==1,elected,0)
      # Note: replace elected status for those who do not fulfill threshold
      
    }
    
    # d) Calculate election status of candidate i (for all countries other than Brazil)
    
    if (method!="Brazil"){ 
    data_input_party <- data_input[data_input$party==j,]
    row_number <- row_number(data_input$votes_h[data_input$party==j])
    rank_h <- rank(-data_input$votes_h[data_input$party==j]+row_number/1000)[data_input_party$party==j & data_input_party$id_cand==i]
    elected <- ifelse(rank_h > seats,1,0)
    }
  }
  
  # 3. Iteration 2: Find exact votessim that changes the election status of candidate i
  
  # a) Define new interval and update votes
  
  xlow <- votessim / 2
  xhigh <- votessim
  
  while ((xhigh - xlow) > convcrit) {
    votessim <- (xlow + xhigh) / 2
    data_input$votes_h[data_input$id_cand==i] <- votessim
    #print(paste0("votessim:",votessim))
    
    if (method!="Brazil"){
      data_input$votes_j[data_input$party==j] <- votes_j_initial-votes_h_initial+votessim
    }
    
    # b) Calculate a party's number of seats (for all countries other than Brazil)
    #    or election status of a candidate i (for Brazil)
    

    if (method=="SainteLague" | method=="dHondt"){
      seats <- GetSeatsHighestAverage(party_input=j,
                                      data_input=data_input, 
                                      n=n,
                                      method=method,
                                      threshold=threshold)
    }
    
    if (method=="Brazil"){ 
      
      elected <- GetElectionStatusBrazil(data_input=data_input,
                                         n=sum(data_input$elected),
                                         method=method,
                                         threshold=threshold)%>% 
        filter(id_cand==i) %>%
        select(elected) %>%
        pull()
      
      #print(paste0("elected:",elected))
      
      
      elected <- ifelse(length(elected)==1,elected,0)
      
      
    }
    
    # c) Calculate election status of candidate i (for all countries other than Brazil)
    
    if (method!="Brazil"){ 
      data_input_party <- data_input[data_input$party==j,]
      row_number <- row_number(data_input$votes_h[data_input$party==j])
      rank_h <- rank(-data_input$votes_h[data_input$party==j]+row_number/1000)[data_input_party$party==j & data_input_party$id_cand==i]
      elected <- ifelse(rank_h > seats,1,0)
    }
  
    # d) Update interval for search

    if (elected==0) {
      xlow <- votessim
    }
    if (elected==1) {
      xhigh <- votessim
    }
  }
 
  return(votessim)
}



GetVotesRequieredOpenList_MP <- function(id, 
                                      data_input,
                                      n,
                                      method="dHondt",
                                      threshold,
                                      convcrit=0.001){
  
  j <- data_input$party[data_input$id==id]
  i <- data_input$seat[data_input$id==id]
  
  votessim <- 0.5  # note: start with 1 vote in first iteration (first line after while loop)
  seats <- 0
  rank_h <- length(data_input$votes_j[data_input$party==j])
  
  votes_h_initial <- data_input$votes_h[data_input$id_cand==i]
  votes_j_initial <- data_input$votes_j[data_input$party==j][1]
  
  while (rank_h > seats) {
    #print(votessim)
    votessim <- votessim * 2
    data_input$votes_j[data_input$party==j] <- votes_j_initial-votes_h_initial+votessim
    data_input$votes_h[data_input$id_cand==i] <- votessim
    
    if (method=="SainteLague" | method=="dHondt"){
      seats <- GetSeatsHighestAverage(party_input=j,
                                      data_input=data_input, 
                                      n=n,
                                      method=method,
                                      threshold=threshold)
      
    }
    
    if (method=="Hare"){ #continue here
      seats <- GetSeatsHare(party_input=j,
                            data_input=data_input,
                            n=n,
                            method=method,
                            threshold=threshold)
      
    }
    data_input_party <- data_input[data_input$party==j,]
    row_number <- row_number(data_input$votes_h[data_input$party==j])
    rank_h <- rank(-data_input$votes_h[data_input$party==j]+row_number/1000)[data_input_party$party==j & data_input_party$id_cand==i]
  }
  
  xlow <- votessim / 2
  xhigh <- votessim
  
  while ((xhigh - xlow) > convcrit) {
    votessim <- (xlow + xhigh) / 2
    data_input$votes_j[data_input$party==j] <- votes_j_initial-votes_h_initial+votessim
    data_input$votes_h[data_input$id_cand==i] <- votessim
    
    if (method=="SainteLague" | method=="dHondt"){
      seats <- GetSeatsHighestAverage(party_input=j,
                                      data_input=data_input, 
                                      n=n,
                                      method=method,
                                      threshold=threshold)
    }
    data_input_party <- data_input[data_input$party==j,]
    row_number <- row_number(data_input$votes_h[data_input$party==j])
    rank_h <- rank(-data_input$votes_h[data_input$party==j]+row_number/1000)[data_input_party$party==j & data_input_party$id_cand==i]
    
    if (rank_h > seats) {
      xlow <- votessim
    }
    if (rank_h <= seats) {
      xhigh <- votessim
    }
  }
  
  return(votessim)
}

GetRVSimulation_Candidate_Loop <- function(data_input, 
                                      convcrit=0.001,
                                      method="dHondt",
                                      system="open",
                                      threshold="none",
                                      rank="rank",
                                      additional_vars,
                                      print_seat=FALSE,
                                      print_party=FALSE) {  
  #---------------------------------------------------------------------------
  # data_input: data after PrepareData with votes_j as party votes and votes_h 
  #             as candidate votes
  # n:          number of seats
  # convcrit:   convergence criterion
  # method:     method to distribute seats
  # system:     close-list or open-list
  #---------------------------------------------------------------------------  
  df_out <- data.frame()
  data_input$seat <- data_input[[rank]]
  
      df_out <- rbind(df_out, GetRVSimulation_Candidate(data_input=data_input, 
                                                          n=sum(data_input$elected), 
                                                          convcrit=convcrit,
                                                          method=method,
                                                          system=system,
                                                          threshold=threshold,
                                                          print_seat=print_seat,
                                                          print_party=print_party))

  
  if (system=="open"){
    df_out <- df_out %>%
      group_by(party) %>%
      mutate(seat = dense_rank(as.numeric(-votemargin+
                                                row_number()/1000))) 
    
    data_input <- data_input%>%
      group_by(party) %>%
      mutate(seat = dense_rank(as.numeric(-votes_h+
                                            row_number()/1000))) 
    
  }


  return(df_out %>% left_join(data_input %>% mutate(party=as.character(party)) , by = c("party", "seat")) %>% 
         select( year, party,votemargin, all_of(additional_vars))
         )
}



#-----------------------
# (Z) Test of Functions
#-----------------------  

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
# (PartyMarginOut <- CalculatePartyMargins(data_sz_2015_out ,4))
#
# (CandidateMarginOut <- CalculateCandidateMargin(13,data_sz_2015_out ))
# (CalculateCandidateMarginsOut <- CalculateCandidateMargins(data_sz_out))
# min(CalculateCandidateMarginsOut$candmargin[CalculateCandidateMarginsOut$rank_h<
# CalculateCandidateMarginsOut$rank_h_compare])
# max(CalculateCandidateMarginsOut$candmargin[CalculateCandidateMarginsOut$rank_h>
#CalculateCandidateMarginsOut$rank_h_compare])
#
# PartyMarginOut <- PartyMarginOut %>% rename(mergevar=seat)
# CalculateCandidateMarginsOut <- CalculateCandidateMarginsOut %>% 
# mutate(mergevar=ifelse(rank_h<rank_h_compare,rank_h_compare-1,rank_h_compare))
#
# MarginsMerged <- CalculateCandidateMarginsOut %>%
#                               left_join(PartyMarginOut,by=c("party","mergevar")) %>%
#                               left_join(data_sz_out,by=c("party","rank_h")) %>%
#                               select(firstname,name,candmargin,partymargin,votes,party,
#                                       rank_h, rank_h_compare)
#
# (MarginsMerged <- MarginsMerged  %>%  mutate(margin_min = pmin(candmargin,partymargin)))
#
#
# (MarginsFinal <- MarginsMerged %>% group_by(party,rank_h) %>% summarize(votemargin=max(margin_min))  )%>% print(n=100)
#


###############################
# (C) Test of RV Functions
###############################
# 
# df <- data.frame(
#   party = c(rep("A", 3), rep("B", 3), rep("C", 3)),
#   votes = c(22, 8, 5,20, 15, 10, 12, 5, 3),
#   firstname = c( "Brenda", "Charles", "Diana", "Z", "X", "Y","V", "W", "U"),
#   name = rep("", 9),
#   elected = c( 1, 0, 0, 1, 1, 0,0, 0, 0)
# ) %>%
#   group_by(party) %>%
#   mutate(pvotes = sum(votes), alliance = "", suballiance = "", canton = 1, year = 2013)
# 
# df_out <- PrepareData(df, "votes", "pvotes", "alliance", "suballiance", "party", "canton", "year")
# AggregateMargins(df_out, 3)
# 
# # for aggregation table in presentation (from AggregateMargins function)
# 
# Margins_out <- MarginsMerged  %>% arrange(rank_h,seat) %>% 
#                   select(rank_h,seat,partymargin,candmargin)%>%
#                   mutate(rank_h=factor(rank_h,labels=c("Brenda","Charles",
#                                                        "Diana")))
# library(xtable)
# print(xtable(Margins_out,digits=0),include.rownames=FALSE)
# 
# PartyMarginOut <- PartyMarginOut %>% mutate( party= c(rep("A", 3), rep("B", 3), rep("C", 3)))%>%
#                   select(-mergevar)
# 
# print(xtable(PartyMarginOut,digits=c(0,0,0,1)),include.rownames=FALSE)
# 
# 
# 

# 
# # (i) Run program for selected canton-years (SZ 2015 with all possible alliance/suballiance combinations, AG 1931 only with alliances )
# 
# data_all <- read.dta13("./02_Processed_data/nationalraete_1931_2015.dta") %>%
#   as_tibble() %>%
#   select(ID, canton, year, firstname, name, votes, pvotes, list, alliance, suballiance, elected) %>%
#   mutate(party=as.character(list),
#          party_orig = list, 
#          alliance_orig = alliance, 
#          suballiance_orig = suballiance,
#          pvotes=ifelse(is.na(pvotes)==F,pvotes,votes))  # replace party votes for single-party candidates (e.g. Beat Graf, AI 1987)
#          
# 
# data_sz_2015 <- data_all %>% filter(canton == "SZ" & year == 2015)
# data_sz_2015 %>% print(n = 100)
# 
# data_ag_1931 <- data_all %>% filter(canton == "AG" & year == 1931)
# data_ar_1931 <- data_all %>% filter(canton == "AR" & year == 1931)
# data_be_1931 <- data_all %>% filter(canton == "BE" & year == 1931)
# data_ti_2011 <- data_all %>% filter(canton == "TI" & year == 2011)
# data_vs_1967<- data_all %>% filter(canton == "VS" & year == 1967)
# 
# # data <- data_sz_2015 # temporary
# 
# data_sz_2015_out <- PrepareData(data_sz_2015, "votes", "pvotes", "alliance_orig", "suballiance_orig", "list", "canton", "year")
# data_ag_1931_out <- PrepareData(data_ag_1931, "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year")
# data_ar_1931_out <- PrepareData(data_ar_1931, "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year")
# data_be_1931_out <- PrepareData(data_be_1931, "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year")
# data_ti_2011_out <- PrepareData(data_ti_2011, "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year")
# data_vs_1967_out <- PrepareData(data_vs_1967, "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year")
# data_zh_1967_out <- PrepareData( data_all %>% filter(canton == "ZH" & year == 1967), "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year")
# data_ai_1931_out <- PrepareData(data_all %>% filter(canton == "AI" & year == 1931), "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year")
# data_ai_1987_out <- PrepareData(data_all %>% filter(canton == "AI" & year == 1987), "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year")
# 
# 
# CalculateCandidateMargin(party_input = 4,data_input = data_zh_1967_out) %>% filter(rank_h<3)
# CalculatePartyMargin(party_input = 4,data_input = data_zh_1967_out,n=sum(data_zh_1967_out$elected))
# 
# CalculateCandidateMargin(party_input =1,data_input = data_ai_1931_out)
# CalculatePartyMargin(party_input =1,data_input = data_ai_1931_out,n=sum(data_ai_1931_out$elected))
# 
# CalculateCandidateMargin(party_input = 2,data_input = data_ai_1987_out) 
# CalculatePartyMargin(party_input = 2,data_input = data_ai_1987_out,n=1)
# 
# 
# AggregateMargins(data_ag_1931_out, 12)
# AggregateMargins(data_ar_1931_out, 2)
# AggregateMargins(data_ai_1931_out, 2)
# AggregateMargins(data_ai_1987_out, 1)
# 
# 
# print(AggregateMargins(data_be_1931_out, sum(data_be_1931_out$elected)),n=145)
# aggmargin_ti_2011 <- AggregateMargins(data_ti_2011_out, 8) 
# aggmargin_ti_2011%>% filter(party==6)
# aggmargin_vs_1967 <- AggregateMargins(data_vs_1967_out, sum(data_vs_1967$elected)) 
# aggmargin_vs_1967%>% filter(name=="Dellberg")
# aggmargin_zh_1967 <-AggregateMargins(data_zh_1967_out, 35)
# aggmargin_zh_1967%>% arrange(-votes_h) %>% filter(party==4)
# %>% arrange(-votes_h) %>% filter(party==4)
# 
# 
# # (ii) Run program for full dataset
# 
# data_prepared_1931_2015 <- PrepareData(data_all, "votes", "pvotes", "alliance", "suballiance", "list", "canton", "year") %>% filter(is.na(votes_h) == F) # filter out "stille Wahlen"
# data_prepared_1931_2015 %>%
#   group_by(year, canton) %>%
#   summarize(nobs = n()) %>%
#   filter(nobs == 1)
# 
# start_time <- Sys.time()
# RV_out <- CalculateMargins(data_prepared_1931_2015)
# end_time <- Sys.time()
# end_time - start_time
# 
# RV_out  %>% filter(canton=="ZH" & year==1967 & list=="Liste 12. Freisinnige Liste Land") %>% arrange(-votes_h)
# RV_out  %>% filter(canton=="AI" & year==1931) %>% arrange(-votes_h)
# 
# data.table::fwrite(RV_out, file = "./02_Processed_data/13_Running_variable/RV_Analytical_R.csv", sep = ";",row.names = F)
# 
# 
# # (iii) Check aggregate values
# 
# RV_out %>% filter(is.na(votemargin))
# 
# RV_out %>%
#   filter(is.na(votemargin) == F) %>%
#   group_by(year, elected) %>%
#   summarize(votemargin_min = min(votemargin), votemargin_max = max(votemargin)) %>%
#   print(n = 100)
# 
# 
# RV_out %>% filter(votemargin>0 & elected==0) %>% select(ID,name,firstname, year,canton, votemargin ,elected, starts_with("vote"))
# RV_out %>% filter(votemargin<0 & elected==1) %>% select(ID,name,firstname, year,canton, votemargin ,elected, starts_with("vote"))
# 
# 
# 
# 
# # (iii) Check R implementation with previous implementation from Stata
# 
# RV_out <- read.table(file = "./02_Processed_data/13_Running_variable/rv_analytical_r.csv", sep = ";")
# #head(RV_out)
# 
# data_rv_stata <- read.dta13("./02_Processed_data/13_Running_variable/RV_Analytical_all_old.dta") %>%
#   as_tibble() %>%
#   arrange(year, party_num, ID) %>%
#   select(year, ID, votemargin) %>%
#   rename(votemargin_stata = votemargin)
# 
# RV_compare_R_Stata <- RV_out %>%
#   left_join(data_rv_stata, by = c("ID", "year")) %>%
#   mutate(votemargin_diff_rstata = votemargin - votemargin_stata) %>%
#   select(year, canton, ID, name, firstname, elected, votemargin, votemargin_stata, votemargin_diff_rstata,  alliance, suballiance, alliance_dummy, suballiance_dummy) %>%
#   as_tibble()
# 
# RV_compare_R_Stata_diff <- RV_compare_R_Stata %>% filter(votemargin_diff_rstata > 1) 
# RV_compare_R_Stata_diff %>%tabyl(canton,year)
# 
# RV_compare_R_Stata %>% filter(canton=="ZH" & year==1967 & name=="B?hler") %>% print(n=100) 
# 
# options("encoding" = "UTF-8")
# write.table(RV_compare_R_Stata_diff, file = "./02_Processed_data/13_Running_variable/RV_compare_R_Stata_diff.csv", sep = ";")
