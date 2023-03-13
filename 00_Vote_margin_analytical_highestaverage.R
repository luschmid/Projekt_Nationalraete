
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

#-------------------------------------------------------------------------------
# (B) Prepare data and collapse votes
# 
# Description: PrepareData prepares the data for the later functions 
#              GetAllianceSuballiance gets information on alliances/suballiances
#              CollapseVotes gets party votes from the  candidate dataset
#              CalculateRatios calculates Ratios (dHondt numbers) in long format
#-------------------------------------------------------------------------------

PrepareData <- function(data,votes_j_name,
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

PartyMargin <- function(party, data_input, n,method="dHondt",
                        margin_type="Votes_Added") {
# margin_type: "Votes_Added" for only adding votes (as in L?chinger, Schelker, 
#               Schmid 2021), "easiest" for adding or subtracting votes
#               as in Fiva and Smith (2018)
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
  
  if (margin_type=="easiest"){
    
    data_out <- data_out  %>%
      left_join(data_others, by = c("seat_merge")) %>%
      mutate(mult=pmin(mult.x , mult.y)) %>%
      mutate(margin = (as.numeric(No.x) - as.numeric(No.y)) * mult) %>%
      rename(seat=seat.x)%>%
      select(party, seat, margin)
  }
  
  if (margin_type=="Votes_Added"){
    
    data_out <- data_out  %>%
      left_join(data_others, by = c("seat_merge")) %>%
      mutate(mult=mult.x) %>%
      mutate(margin = (as.numeric(No.x) - as.numeric(No.y)) * mult) %>%
      rename(seat=seat.x)%>%
      select(party, seat, margin)
  }
  return(data_out)
}

CalculatePartyMargin <- function(party_input, 
                                 data_input, 
                                 n,
                                 method="dHondt",
                                 margin_type="Votes_Added",
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
                                method=method,
                                margin_type=margin_type)
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
                                       method=method,
                                       margin_type=margin_type)
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
                                   method=method,
                                   margin_type=margin_type)
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

CalculatePartyMargins <- function(data_input, n,method="dHondt",
                                  margin_type="Votes_Added") {
  Out <- data.frame()
  all_parties <- as.numeric(levels(as.factor(data_input$party)))
  for (j in 1:length(all_parties)) {
    #no_candidates <- data_input[as.factor(data_input$party)==all_parties[j]]
    Out <- rbind(Out, CalculatePartyMargin(party=all_parties[j], data_input,
                                           margin_type=margin_type,
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
                             margin_type="Votes_Added",
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
                                          n=n,method=method,
                                          margin_type=margin_type)%>% 
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
    summarize(votemargin = max(margin_min),candmargin=min(candmargin)) 
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
             select(year, party,votemargin)) #, all_of(additional_vars)
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
                             margin_type="Votes_Added",
                             additional_vars="",
                             threshold="none",
                             rv_type="analytical",
                             return_option=T,
                             outfile_name) {
  data_singlecandidate <- data_input %>%
    group_by(year, district) %>%
    summarize(nobs = n()) %>%
    filter(nobs == 1)
  #data_input <- data_input %>%
   # left_join(data_singlecandidate, by = c("year", "district")) %>%
    #filter(is.na(nobs) == T) 
  # note: do calculation only for districts with more than one candidate
  years <- as.numeric(levels(as.factor(data_input$year)))
  Out <- data.frame()
  for (yr in years) {
    print(paste("New year:", yr))
    dataperyear <- data_input[data_input$year == yr, ]
    districts <- levels(as.factor(dataperyear$district))
    for (di in districts) {
      print(paste("New district:", di))
      dataworking <- dataperyear[dataperyear$district == di, ]
      n <- sum(dataworking$elected)
      if (rv_type=="analytical"){
      Out <- rbind(Out, as.data.frame(AggregateMargins(dataworking, 
                                              n,system=system,
                                              method=method,
                                              rank=rank,
                                              margin_type=margin_type,
                                              threshold=threshold,
                                              additional_vars=additional_vars)))
      }
      
      if (rv_type=="simulation"){
        
        out_temp <- as.data.frame(GetRVSimulation_Candidate_Loop(data_input=dataworking, 
                                                                 system=system,
                                                                 method=method,
                                                                 rank = rank,
                                                                 convcrit=convcrit,
                                                                 margin_type=margin_type,
                                                                 threshold=threshold,
                                                                 additional_vars=additional_vars))

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
  
  quorum=T
  
  if (threshold!="none"){
  data_input <- CutDataByQuorum(data_input,threshold)
  quorum <- is.element(party_input,levels(as.factor(data_input$party)))
  
  if (quorum==F){ out <- 0  # calculation only for parties who made quorum
  }
  }
  
  if (quorum!=F | threshold=="none"){
  
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
      filter(No_rank <= n) %>%
      group_by(party) %>%
      summarize(seats_suballiance=n())%>%
      rename(suballiance=party)
    
    out <- max(SuabllianceSeats$seats_alliance[SuballianceSeats$suballiance==
                                          SuballianceInformation$suballiance],0)

    # 3. Calculation of margin across parties (only for parties in suballiance)
    if (AllianceInformation$alliance_dummy == 1 & 
        AllianceInformation$suballiance_dummy == 1 & out>0) {
      VotesParty <- CollapseVotes(data_input, aggregation_level = party, 
                             aggregation_unit = AllianceInformation$suballiance)
      PartySeats <- CalculateRatios(votes=VotesParty, 
                                    n=n,method=method) %>% 
        mutate(No_rank = row_number(-as.numeric(No))) %>%
        filter(No_rank <= n) %>%
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
districts_all <- levels(as.factor(as.character(df1$districtid)))
#print(districts_all)
for (district in districts_all) {
df2 <- df1[df1$districtid==district,]
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

#-------------------------------------------------------------------------------
# (G) Simulate vote margin
# 
# Description: GetMinVotesQuorum calculates the minimum number of party votes required
#                   for party j to achieve a seat in systems with a quorum
#             GetRVSimulation_Candidate calculates the vote margin for all
#                candidates in a district using a simulation that adds votes 
#                until the respective candidate is elected
#              GetRVSimulation_Candidate_Redistribution  
#                calculates the vote margin for all candidates in a district 
#                using a simulation that redistributes party votes from a party
#                until the respective candidate is elected
#               GetRVSimulation_Candidate_Loop calculates simulated vote margin
#                for multiple years and districts
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
                                      threshold="none") {
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
    
    for (i in 1:n) {
  
      
      votessim <- GetVotesRequieredClosedList(data_input=data_input,
                                    seat=i,
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
      print(paste0("Party: ",j))
      
      for (i in ids_all) {
      print(paste0("ID: ",i))
        
        #print(paste("Party: ",j))
        
        votessim <- GetVotesRequieredOpenList(data_input=data_input,
                                      i=i,
                                      j=j,
                                      n=n,
                                      method=method,
                                      threshold=threshold,
                                      convcrit=convcrit)

        
        out_ji <- data.frame(party = j, seat = i, 
                             votemargin = round(data_input$votes_h[data_input$id_cand==i & data_input$party==j]
                                                - votessim, 4))
        out <- rbind(out, out_ji)
      }
    }
  }
  
  return(out)
}

GetVotesRequieredClosedList <- function(data_input,
                              j,
                              n,
                              threshold,
                              convcrit=0.001){
  
    votessim <- 0.5  # note: start with 1 vote in first iteration (first line after while loop)
    seats <- 0
    
    while (seats < i) {
      votessim <- votessim * 2
      data_input$votes_j[data_input$party==j] <- 
        votessim
      if (method=="SainteLague" | method=="dHondt"){
        seats <- GetSeatsHighestAverage(party_input=j,
                                        data_input=data_input, 
                                        n=n,
                                        method=method,
                                        threshold=threshold)
      }
    }
    xlow <- votessim / 2
    xhigh <- votessim
    

    while ((xhigh - xlow) > convcrit & i_new<10) {
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
                                      margin_type="Votes_Added") {  
  #---------------------------------------------------------------------------
  # data_input: data after PrepareData with votes_j as party votes and votes_h 
  #             as candidate votes
  # n:          number of seats
  # convcrit:   convergence criterion
  # method:     method to distribute seats
  # system:     close-list or open-list
  # type:       "Party_redistribution" if redistribution of party votes 
  #             "Votes_Added" if votes are added (our standard case)
  #             "Votes_Subtracted_Or_Added" for easiest ways as in 
  #             Fiva & Smith (2018)
  #---------------------------------------------------------------------------  
  df_out <- data.frame()
  data_input$seat <- data_input[[rank]]
  
      if (margin_type== "Votes_Added"){
      df_out <- rbind(df_out, GetRVSimulation_Candidate(data_input=data_input, 
                                                          n=sum(data_input$elected), 
                                                          convcrit=convcrit,
                                                          method=method,
                                                          system=system,
                                                          threshold=threshold))
      }
  
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

#---------------
# Rep package PA
#---------------


library(readstata13)
path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte"

data_all <- read.dta13(paste0(path,"./05_Texts_and_presns/01_Running_variable/submissions/17 pol analysis final/replication package/data_switzerland.dta")) %>%
  as_tibble() %>%
  mutate(pvotes=ifelse(is.na(pvotes)==F,pvotes,votes))  # replace party votes for single-party candidates (e.g. Beat Graf, AI 1987)

data_prepared_1931_2015 <- PrepareData(data_all,  
                                       votes_h_name="votes", 
                                       votes_j_name= "pvotes",
                                       alliance_name="alliance",
                                       suballiance_name="suballiance", 
                                       party_name="list", 
                                       districtname="canton", 
                                       election_cyclename = "year",
                                       system="open",
                                       cand_id_name="id",
                                       alliances=TRUE) %>% 
  filter(is.na(votes_h) == F) # Note: filter out tacit elections

start_time <- Sys.time()
RV_out <- CalculateMargins(data_prepared_1931_2015)
end_time <- Sys.time()
end_time - start_time


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
