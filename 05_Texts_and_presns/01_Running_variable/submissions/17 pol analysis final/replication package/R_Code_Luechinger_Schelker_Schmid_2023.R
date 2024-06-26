#--------------------
# A) Install packages
#--------------------

install.packages("devtools")
library(devtools)
install_version("data.table", version = "1.14.6")
install_version("furrr", version = "0.3.1")
install_version("haven", version = "2.5.1")
install_version("RColorBrewer", version = "1.1.3")
install_version("Rmpfr", version = "0.8-9")
install_version("tidyverse", version = "1.3.2")

#--------------------------------------
# B) General settings and load packages
#--------------------------------------

rm(list = ls(all = TRUE))
path <- "..."

library(data.table)
library(furrr)
library(haven)
library(RColorBrewer)
library(Rmpfr)
library(tidyverse)

options(dplyr.summarise.inform = FALSE,
        dplyr.summarise.inform = FALSE)
options(warn = 1)
options("encoding" = "UTF-8")


#---------------------
# C) Functions 
#---------------------
# Description: The functions are organized in the following five sections: 
#             (i) Prepare data and collapse votes
#             (ii) Calculate party margin for highest average systems
#             (iii) Calculate party margin for largest remainder systems
#             (iv) Calculate candidate margin
#             (v) Aggregate candidate and party margin
#             (vi) Simulation approaches

# (i) Prepare data and collapse votes
# 
# Description: PrepareData prepares the data for the later functions 
#              GetAllianceSuballiance gets information on alliances/suballiances
#              CollapseVotes gets party votes from the  candidate dataset
#              CalculateRatios calculates Ratios (dHondt numbers) in long format

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

# (ii) Calculate party margin for highest average systems

# Description:  PartyMargin calculates the party margin for a single alliance/
#                party/suballiance
#               CalculatePartyMargin calculates the party margin for alliance/
#                party/suballiance and aggregates this margin
#               CalculatePartyMargins calculates party margin for different 
#                parties


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
    Allpartymargins <- Alliancemargin %>% 
      left_join(Suballiancemargin_all, by = c("seat_a")) %>%
      full_join(Partymargin_all, by = c("seat_s"),
                relationship = "many-to-many") 
    # note: This "many-to-many" merge is needed b/c of nested structure of
    #       the difference margins at the alliance, suballiance, and party level.
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


# (iii) Calculate party margin for largest remainder systems


CheckNA <- function(x) {
  if (any(is.na(x))) {
    print("The vector of votes includes missing values.")
    return(1)
  }
  else {
    return(0)
  }
}


GetRatioHare <- function(votes, n) {
  if (CheckNA(votes) == 0) {
    quota <- mpfr(sum(votes) / n,128)
    return(votes / quota)
  }
}


GetIntegerHare <- function(votes, n) {
  if (CheckNA(votes) == 0) {
    ratio <- GetRatioHare(votes,n)
    return(floor(ratio))
  }
}

GetRemainderHare <- function(votes, n) {
  return(GetRatioHare(votes,n) - GetIntegerHare(votes, n))
}

GetHighestRemainder <- function(remainder, freeseats) {
  return(ifelse(rank(-remainder, ties.method = "random") <= freeseats, 1, 0))
}

GetSeatsHare <- function(votes, n) {
  Seats_Round1 <- GetIntegerHare(votes, n)
  Remainders_Round1 <- GetRemainderHare(votes, n)
  freeseats <- n - sum(Seats_Round1)
  Seats_Round2 <- GetHighestRemainder(Remainders_Round1, freeseats)
  return(Seats_Round1 + Seats_Round2)
}


GetVotesReqParty <- function(deltai,votes_minusj,votes_others,n){
  return(as.numeric(mpfr((votes_minusj* n + deltai * votes_others) / (n - deltai),128)))
}

GetIntegerHareModVotes <- function(votes_req_j,votes,j,n){
  votes[j] <- votes_req_j
  return(as.numeric(GetIntegerHare(c(votes),n)[j])+1)
}


GetVotesReq <- function(j, minusj,votes, n) {
  votes_others <- sum(votes[row(as.matrix(votes)) != j])
  l1 <- unlist(lapply(c(-n:(n-1)),GetVotesReqParty,votes[minusj],votes_others,n))
  l2 <- unlist(lapply(l1,GetIntegerHareModVotes,votes,j,n))
  out_df <- matrix(c(l1,l2,rep(j,(2*n)),rep(minusj,(2*n)),c(-(n-1):(n))),ncol=5,nrow=(2*n))
  colnames(out_df) <- c("votes_req_all","i","j","minusj","deltai")
  return(out_df)
}



GetVotesReqAll <- function(votes, n) {
  for (j in 1:length(votes)) {
    for (minusj in c(1:length(votes))[-j]) {
      VRi <- GetVotesReq(j, minusj,votes, n)
      if (j == 1 & minusj == c(1:length(votes))[-j][1]) {
        VRi_all <- VRi
      } else {
        VRi_all <- rbind(VRi_all, VRi)
      }
    }
  }
  
  VR <- cbind(VRi_all,NA,NA,NA,NA)
  for (u in 1:dim(VR)[1]) {
    minusj <- as.numeric(VR[u,4])
    j <- as.numeric(VR[u,3])
    votes_mod <- mpfr(votes,128)
    votes_mod[j] <- VR[u,1]
    minusj_revelant <- ifelse(j>minusj,minusj,minusj-1)
    VR[u,6] <- as.numeric(n - sum(GetIntegerHare(votes_mod, n))) # free seats after round 1
    VR[u,7] <- as.numeric(rank(-GetRemainderHare(votes_mod, n)[-j])[minusj_revelant]) # rank -j
    VR[u,8] <- as.numeric(GetRatioHare(votes_mod, n)[j]-floor(GetRatioHare(votes_mod, n)[j])) # remainderj
    VR[u,9] <- as.numeric(GetRatioHare(votes_mod, n)[minusj]-floor(GetRatioHare(votes_mod, n)[minusj])) # remainderminusj
  }
  colnames(VR) <- c("votes_req_all","i","j","minusj","deltai","freeseats","rank_minusj","remainderj","remainderminusj")
  return(VR)
}




GetVSHare <- function(data_input, n) {
  if (is.data.frame(data_input)==TRUE){
    votes <- data_input$votes_j
    VR <- GetVotesReqAll(votes, n) 
    VR <- VR[which(VR[,6]==VR[,7]  & abs(VR[,8]-VR[,9])<0.01 ),]  
    VR_df <- data.frame(votes_req = as.numeric(VR[,1]),party= as.numeric(VR[,3]),seat=as.numeric(VR[,2]))
    votes_df <- data.frame(votes = votes, party=c(1:length(votes)), party_original= data_input$party)
    VR_out <- VR_df %>%
      left_join(votes_df, by = c("party")) %>%
      mutate(VS = votes - as.numeric(votes_req)) %>%
      arrange(party_original, seat) %>%
      filter(seat>0) %>%
      select(-party) %>%
      dplyr::rename(party=party_original)
  } else{
    votes <- data_input
    VR <- GetVotesReqAll(votes, n) 
    VR <- VR[which(VR[,6]==VR[,7]  & abs(VR[,8]-VR[,9])<0.01 ),]  
    VR_df <- data.frame(votes_req = as.numeric(VR[,1]),party= as.numeric(VR[,3]),seat=as.numeric(VR[,2]))
    votes_df <- data.frame(votes = votes, party=c(1:length(votes)))
    VR_out <- VR_df %>%
      left_join(votes_df, by = c("party")) %>%
      mutate(VS = votes - as.numeric(votes_req)) %>%
      arrange(party_original, seat) %>%
      filter(seat>0) 
  }
  
  return(VR_out)
}


# (iv) Calculate candidate margin

# Description:  CalculateCandidateMargin calculates the candidate margin for all
#                candidates in a party
#               CalculateCandidateMargins calculates candidate margin for 
#                different parties

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
          # note: min is for cases with more than one candidate 
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



# (v) Aggregate candidate and party margin
#
# Description:  AggregateMargins aggregates the party and the candidate margin  
#                for all candidates in a party

AggregateMargins <- function(data_input, 
                             n,
                             system="open",
                             method="dHondt",
                             rank="rank",
                             candmargin_return=0,
                             additional_vars="",
                             threshold="none"){
  #--------------------------------------------------------------------------
  # This function aggregates party and candidate margin for open-list systems
  #--------------------------------------------------------------------------

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
      summarize(votemargin = max(margin_min),candmargin=min(candmargin)) 
    # note: maximum of all potential seat configurations
    MarginsFinal <- MarginsFinal %>%
      mutate(votemargin=ifelse(is.na(votemargin)==T,candmargin,votemargin)) 
    # note: replace votemargin as candmargin for all single candidates in which on  
    #       other list exists (e.g., Daehler	Edmund, AI 1931)
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
  
  if(candmargin_return==1){
  return(MarginsFinal %>% left_join(data_input, by = c("party", "rank_h")) %>% 
           select(year, party,votemargin,candmargin,all_of(additional_vars)))
  }else{
    return(MarginsFinal %>% left_join(data_input, by = c("party", "rank_h")) %>% 
             select(year, party,votemargin,all_of(additional_vars)))
  }
}



AggregateMarginsLargestRemainder <- function(data_input, 
                                             no_seats,
                                             system="open",
                                             method="Hare",
                                             candmargin_return=0,
                                             rank="rank",
                                             additional_vars="",
                                             threshold="none"){ 
  
  
  # (i) Aggregate party votes
  if (any(colnames(data_input) == "votes_j")==F){ # a) calculate party votes as sum of individual votes if party votes are not present
    data_inputparty <- data_input %>% 
      group_by(district,year,party) %>%
      summarize(votes_j=sum(votes_h)) 
    data_input <- data_input %>% left_join(data_inputparty,by=c("district","year","party")) %>%
      group_by(district,year,party)  %>%
      mutate(rank_h_inv = dense_rank(as.numeric(votes_h))) %>%
      mutate(rank_h_inv_max = max(rank_h_inv)) %>%
      mutate(rank_h = rank_h_inv_max - rank_h_inv + 1)
  }  else { # b) calculate party votes as mean of party votes if party votes are present
    data_inputparty <- data_input %>% 
      group_by(district,year,party) %>%
      summarize(votes_j=mean(pvotes)) 
    data_input <- data_input %>%  group_by(district,year,party)  %>%
      mutate(rank_h_inv = dense_rank(as.numeric(votes_h))) %>%
      mutate(rank_h_inv_max = max(rank_h_inv)) %>%
      mutate(rank_h = rank_h_inv_max - rank_h_inv + 1)  %>% 
      mutate(party=as.character(party))
  } 
  
  if (method=="Hare" & system=="open") { 
    # (ii) Calculate party and candidate margin
    PartyMarginOut <- GetVSHare(data_inputparty, no_seats) %>% mutate(mergevar = seat) %>% 
      rename(partymargin=VS) # temp: either VS or partymargin
    CandidateMarginsOut <- CalculateCandidateMargins(data_input) %>% 
      mutate(mergevar = ifelse(rank_h < rank_h_compare, rank_h_compare - 1, rank_h_compare))
    # (iii) Aggregate party and candidate margin
    MarginsMerged <- CandidateMarginsOut %>%
      left_join(PartyMarginOut %>% mutate(party=as.character(party)) , by = c("party", "mergevar")) %>%
      mutate(margin_min = ifelse(is.na(candmargin) == F, pmin(candmargin, partymargin),partymargin)) # minimum of candidate and partymargin for all canidates for whom both are available, partymargin for those with no candidate margin (single-candidate lists)
    MarginsFinal <- MarginsMerged %>%
      group_by(party, rank_h) %>%
      summarize(votemargin = max(margin_min)) # maximum of all potential seat configurations
    
  }
  
  if(candmargin_return==1){
    return(MarginsFinal %>% left_join(data_input, by = c("party", "rank_h")) %>% 
             select(year, party,votemargin,candmargin,all_of(additional_vars)))
  }else{
    return(MarginsFinal %>% left_join(data_input, by = c("party", "rank_h")) %>% 
             select(year, party,votemargin,all_of(additional_vars)))
  }

}

# (v) Calculate vote margin for several years, districts, and parties
# 
# Description:  CalculateMargins calculates the vote margin for several  
#                years, districts, and parties

CalculateMargins <- function(data_input,
                             system="open",
                             method="dHondt",
                             rank="rank", 
                             candmargin_return=0,
                             convcrit=0.001,
                             additional_vars="",
                             threshold="none",
                             return_option=T,
                             outfile_name) {
  data_singlecandidate <- data_input %>%
    group_by(year, district) %>%
    summarize(nobs = n()) %>%
    filter(nobs == 1)
  
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
      if (method %in% c("dHondt","SainteLague")){
        Out <- rbind(Out, as.data.frame(AggregateMargins(dataworking, 
                                                         n,system=system,
                                                         method=method,
                                                         candmargin_return=candmargin_return,
                                                         rank=rank,
                                                         threshold=threshold,
                                                         additional_vars=additional_vars)))
      }
      
      if (method %in% c("Hare")){
        Out <- rbind(Out, as.data.frame(AggregateMarginsLargestRemainder(data_input=dataworking, 
                                                                         no_seats=n,
                                                                         candmargin_return=candmargin_return,
                                                                         system=system,
                                                                         method=method,
                                                                         rank=rank,
                                                                         threshold=threshold,additional_vars=additional_vars)))
      }

    }
  }
  if (return_option==T){return(Out)}
}

# (vi) Simulation approaches

RandomizeElectionResults <- function(data,m){
  # This function draws random votes
  # Input: data: dataframe with votes_h_share that are used as probabilities for the random draw
  #         m: number of random votes
  # Output: vector of length m with random votes
  rmultinom(1, size = m, prob = data$votes_h_share)
}

GetPartyVotes <- function(data){
  # This function aggregates candidate votes to party votes
  # Input: data: dataframe with variables year,Department,party_num, votes_j, votes_h
  data_votes_j <- data %>% group_by(party_num) %>%
    select(-votes_j) %>%
    dplyr::summarize(votes_j=sum(votes_h,na.rm=T),seats=first(seats))  %>% ungroup()
  return(data_votes_j)
}

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


GetNumberSeats <- function(data){
  # This function calculates the number of seats of an alliance/suballianc/party based on D'Hondt numbers
  # Input: data with party name and party votes
  #DHondtNumbers <- CalculateDHondtNumbers(data,no_seats=max(data$seats))
  DHondtNumbers <- CalculateRatios(votes=data, n=max(data$seats),method="dHondt") 
  data_out <- DHondtNumbers %>%  
    left_join(data,by=c("party"))%>% 
    mutate(dHondtno_rank=rank(-as.numeric(No),ties.method = "random" )) %>%
    ungroup %>%
    filter(dHondtno_rank<=seats)  
  
  return(data_out %>%  group_by(party) %>%
           summarize(seats=n()) %>% ungroup())
}

GetElected_DHondt <- function(m=10000,data,control_seat_allocation=0){
  # This function calculates one random election result at the candidate level (if control_seat_allocation=0)
  # or the elected status for all candidates (if control_seat_allocation=1)
  # Input: m: number of votes to sample
  #        data: input data with variables id_Stata,party_num,year,elected, votes_h, votes_j, party_num
  #        control_seat_allocation: indicator that is 1 for check of seat distribution and 0 for random elections
  # Output: dataframe that includes id_Stata,party_num,year,elected
  if (control_seat_allocation==0){
    data$votes_h <- RandomizeElectionResults(data,m=m)
  }
  
  # (a) Sum up individual votes to party votes (incl. "additional candidates")
  partyvotes <- GetPartyVotes(data)%>% rename(party=party_num)
  
  # (b) Distribute seats to parties
  
  PartySeats <- GetNumberSeats(partyvotes) %>% rename(party_num=party)
  
  # (c) Get individual election result
  data_final <- data %>% select(-seats) %>% filter(ID!="Additional candidate") %>% 
    left_join(PartySeats,by=c("party_num"))%>%
    mutate(seats=ifelse(is.na(seats),0,seats)) %>%
    group_by(party_num)%>%
    mutate(votes_h_rank=rank(-votes_h,ties.method = "random" )) %>%
    ungroup() %>%
    mutate(elected=ifelse(votes_h_rank<=seats,1,0)) %>%
    select(ID,party_num,elected) %>%
    arrange(party_num,ID)
  return(data_final)
}

GetElectedSum<- function(m,data,M=20000,control_seat_allocation=0){
  # This function calculates M random elections result at the candidate level (if control_seat_allocation=0)
  # or the elected status for all candidates (if control_seat_allocation=1)
  # Input: m: number of votes to sample
  #        M: number of random elections (only if control_seat_allocation=1)
  #        data: input data with variables id_Stata,party_num,year,elected, votes_h, votes_j, party_num
  #        control_seat_allocation: indicator that is 1 for check of seat distribution and 0 for random elections
  # Output: dataframe that includes id_Stata,party_num,year,elected  
  if (control_seat_allocation==1){
    GetElected_DHondt(data=data,control_seat_allocation=1) %>% mutate(elected_sum=elected,M=1)
  }
  else{
    l <- lapply(rep(m,M),GetElected_DHondt,data,control_seat_allocation=0) 
    df <- do.call(rbind.data.frame, l) %>%
      group_by(ID,party_num) %>%
      dplyr::summarize(elected_sum=sum(elected))%>%
      mutate(M=M)
    return(df)
  }
}


KotakorpiEtAl2017<- function(m=10000,data,M=20000,control_seat_allocation=0,check_convergence=1){
  # This function calculates the procedure by Kotakorpi et al. 2017 for one district in a year
  # Input: data: original data
  #        m: number of votes to sample
  #        M: number of random elections (only if control_seat_allocation=1)
  #        control_seat_allocation: indicator that is 1 for check of seat distribution and 0 for random elections
  #        check_convergence: indicator that is 1 if procedure 2 of Kotakorpi et al. 2017 is chosen
  #                           and convergence is checked (Kotakorpi et al. 2017, p. 419)
  # Output: df_out_all that lists the number of times a candidate is elected (elected_sum), the totel 
  #         number of elections (M), and te final running variable p=elected_sum/M
  if (control_seat_allocation==1){
    df_out_all <- GetElectedSum(data=data,control_seat_allocation=1)
  }
  else{
    df_out_all <- GetElectedSum(m,data,M,control_seat_allocation=0)
  }
  return(df_out_all %>% mutate(m=m))
}



KotakorpiEtAl2017Loop<- function(m_seq,data,M=20000,control_seat_allocation=0,check_convergence=1){
  # This function loops over KotakorpiEtAl using parallel
  # computing on multiple cores. 
  no_cores <- availableCores() - 1
  plan(multisession,workers=no_cores)
    out_data <- future_map(m_seq,KotakorpiEtAl2017,
                           data=data,M=M,
                           control_seat_allocation=control_seat_allocation,
                           check_convergence=check_convergence)
return(as.data.frame(do.call(rbind, out_data))
)
}




#---------------------------------------------------
# D) Calculation of running variable for Switzerland
#---------------------------------------------------

# (i) Read in data and prepare it for calculation of running variable

data_swi <- read_dta(paste0(path,"/data_switzerland.dta")) %>%
  as_tibble() %>%
  mutate(pvotes=ifelse(is.na(pvotes)==F,pvotes,votes))  # replace party votes for single-party candidates (e.g. Beat Graf, AI 1987)

data_swi_prepared <- PrepareData(data_swi,
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


# (ii) Calculate running variable

RV_out_swi <- CalculateMargins(data_input=data_swi_prepared,
                           system="open",
                           method="dHondt",
                           additional_vars=c("ID_pers","year"),
                           return_option=T)

RV_out_swi %>% arrange(ID_pers)

# (iii) Merge votemargin with initial dataset and write new Stata dataset

data_swi_out <- data_swi %>% 
  left_join(RV_out_swi %>%
              select(ID_pers,year,votemargin), 
            by=c("ID_pers","year")) %>%
  mutate(votemargin_rel=votemargin/eligible_cant)

write_dta(data_swi_out, paste0(path,"/data_switzerland.dta"))



#-------------------------------------------------
# E) Calculation of running variable for Honduras
#------------------------------------------------

# (i) Read in data and prepare it for calculation of running variable

data_hon <-read_dta(paste0(path,"/data_honduras.dta"))%>% 
  as_tibble() 

data_hon_prepared <- PrepareData(data_hon,
                                 votes_h_name="votes", 
                                 votes_j_name= "pvotes",
                                 party_name="party",
                                 districtname="department",
                                 election_cyclename = "year",
                                 system="open",
                                 cand_id_name="id_Stata",
                                 alliances=FALSE) %>% 
  filter(is.na(votes_h) == F) # Note: filter out tacit elections


# (ii) Calculate running variable


RV_out_hon <- CalculateMargins(data_input=data_hon_prepared,
                           system="open",
                           method="Hare",
                           additional_vars=c("ID_pers","year"),
                           return_option=T) %>%
  as_tibble()


# (iii) Merge votemargin with initial dataset and write new Stata dataset

data_hon_out <- data_hon %>% 
  left_join(RV_out_hon %>%
              select(ID_pers,year,votemargin), 
            by=c("ID_pers","year")) %>% 
  mutate(votemargin_rel=votemargin/total_electoral)

write_dta(data_hon_out, paste0(path,"/data_honduras.dta"))


#-----------------------------------------------
# F) Calculation of running variable for Norway
#----------------------------------------------

# (i) Read in data and look at districts and years


data_nor <-read_dta(paste0(path,"/data_norway.dta"))%>% 
  as_tibble() 

# (ii) Prepare data and calculate vote margin

data_norway_prep <- PrepareData(data=data_nor, 
                                votes_j_name="pvotes",
                                party_name="party", 
                                districtname="districtid", 
                                election_cyclename = "year",
                                alliances=F,
                                system="closed") 

RV_out_nor <- CalculateMargins(data_norway_prep %>% filter(year<=1981),
                                   system="closed",
                                   method="SainteLague",
                                   rank="rank",
                                   additional_vars=c("ID_pers",
                                                     "year"))

# (iii) Merge votemargin with initial dataset and write new Stata dataset

data_nor_out <- data_nor %>% 
  left_join(RV_out_nor %>%
              select(ID_pers,year,votemargin), by=c("ID_pers","year")) %>%
  mutate(votemargin_rel=votemargin/electorate)

write_dta(data_nor_out, paste0(path,"/data_norway.dta"))



#--------------
# G) Figure 4
#--------------

data_elected <- data_swi %>% filter(elected==1) %>% group_by(year, canton,list) %>% summarize(votes_h_min=min(votes))
data_notelected <- data_swi %>% filter(elected==0) %>% group_by(year, canton,list) %>% summarize(votes_h_max=max(votes))

data_all <- data_swi_out %>% 
  select(year, canton,list,votes,votemargin,elected,eligible_cant) %>% 
  left_join(data_elected,by=c("year", "canton","list")) %>% 
  left_join(data_notelected,by=c("year", "canton","list"))%>%
  mutate(votemargin_cand=ifelse(elected==1,(votes-votes_h_max),(votes-votes_h_min)))%>%
  mutate(votemargin_cand_rel=votemargin_cand/eligible_cant,
         votemargin_rel=votemargin/eligible_cant,
         votemargin_binding=ifelse(votemargin_rel==votemargin_cand_rel,1,0),
         novotemargin=ifelse(is.na(votemargin_cand_rel),1,0))

data_all$votemargin_cand_rel[is.na(data_all$votemargin_cand_rel)==T] <- -0.5


ggplot(data_all,aes(x=votemargin_rel,y=votemargin_cand_rel))+
  geom_point() +
  theme_bw(base_size=24)+
  xlab("") + ylab("") +ggtitle("") +
  theme(plot.caption = element_text(hjust = 0, face= "italic"), #Default is hjust=1
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot", 
        plot.title = element_text(size = 24)) +
  scale_y_continuous(breaks=seq(-0.5,0.5,0.25),labels=c("NA",seq(-0.25,0.5,0.25)))

ggsave(file = paste0(path,"/Figures/Figure4.png"), dpi = 900, width = 20, height = 12)

#--------------
# H) Figure 5
#--------------

# (i) Read in data on running example from paper

df <- data.frame(votes_h=c(22,8,5,25,15,5,8,7,5),
                 votes_j=c(rep(35,3),rep(45,3),rep(20,3)),
                 elected=c(1,0,0,1,1,0,0,0,0),
                 party_num=c(rep(1,3),rep(2,3),rep(3,3)),
                 seats=rep(3,9),
                 ID=as.character(1:9))

df$votes_h_share <- df$votes_h/sum(df$votes_h)
data=df

# (ii)  Run simulation

set.seed(05071975) 
m_all <- c(seq(1,10,1),seq(20,100,10),seq(200,1000,100),seq(2000,10000,1000))
df_sim <- KotakorpiEtAl2017Loop(m_seq=m_all,
                                data=df,
                                M=20000,
                                control_seat_allocation=0,
                                check_convergence=0)




# (iii) Read in data from simulation

df_all <- df_sim %>%
  mutate(ID=as.character(ID),
         elected_sum_share=elected_sum/20000) %>%
  left_join(df,by=c("ID","party_num")) %>%
  select(votes_h,votes_j,elected,party_num,ID,votes_h_share,elected_sum_share,m) %>%
  as_tibble() %>%
  mutate(candidate_description=ifelse(ID%in%c(1,4,7),"Candidate 1",ifelse(ID%in%c(2,5,8),"Candidate 2","Candidate 3")))

# (iv) Figure

cols <- brewer.pal(3, "Set1")

ggplot(df_all,aes(x=m,y=elected_sum_share,group=ID,color=as.factor(party_num),
                  linetype=as.factor(candidate_description))) +
  geom_line(size=1.5) +
  theme_bw(base_size=32) +
  xlab("") + ylab("") + 
  scale_color_manual(labels=paste0("Party " ,c(1:3)),values=cols)+
  scale_linetype_manual(values=c("solid","dashed","dotted")) +
  theme(legend.position="bottom",legend.title = element_blank()) +
  scale_x_continuous(trans='log10',limits=c(0.22,10000)) +
  annotate("text", x=0.78,y=0.35, label="0.35", size=8) +
  annotate("text", x=0.78,y=0.45, label="0.45", size=8) +
  annotate("text", x=0.78,y=0.2, label="0.20", size=8) +
  geom_segment(aes(x=0.98, y=0.2, xend=1.05, yend=0.2),  size=2,color="black") +
  geom_segment(aes(x=0.98, y=0.35, xend=1.05, yend=0.35),  size=2,color="black") +
  geom_segment(aes(x=0.98, y=0.45, xend=1.05, yend=0.45),  size=2,color="black") +
  guides(shape = guide_legend(override.aes = list(size = 3)),
         linetype = guide_legend(override.aes = list(shape = NA)),
         size = "none") +
  theme(legend.position = "none")

ggsave(paste0(path,"/Figures/Figure5.pdf"),width=20,height=12)



#--------------
# I) Figure B.1
#--------------

break_points <- unique(c(-rev(seq(0,0.11,0.002)),seq(0,0.11,0.002)))

ggplot(data_swi_out, aes(x = votemargin_rel)) +
  geom_histogram(aes(y = after_stat(density)), fill = "grey", colour = "black", 
                 breaks = break_points) +
  theme_bw(base_size = 42) +
  scale_x_continuous(limits = c(-0.11, 0.11)) +
  geom_vline(xintercept = 0) +
  ylab("") + xlab("")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggsave(file = paste0(path,"/Figures/FigureB1.pdf"), width = 20, height = 12)


#--------------
# J) Figure B.2
#--------------

ggplot(data_hon_out, aes(x = votemargin_rel)) +
  geom_histogram(aes(y = after_stat(density)),
                 alpha=0.5,position = "stack",boundary=0,color="black",
                 binwidth=0.002) +
  theme_bw(base_size = 42) +
  scale_x_continuous(limits = c(-0.07, 0.07), 
                     breaks = seq(-0.06, 0.06, 0.02)) +
  geom_vline(xintercept = 0) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("") + xlab("")

ggsave(file = paste0(path,"/Figures/FigureB2.pdf"), width = 20, height = 12)


#--------------
# K) Figure B.3
#--------------

theme_set(theme_bw(base_size = 42))
theme_update(legend.key = element_rect(size = 6, fill = "white", colour = "white"), 
             legend.key.size = unit(1.5, "cm")) 

ggplot(data_nor_out, aes(x=votemargin_rel,fill=factor(marginal_candidate))) +
  geom_histogram(aes(y = after_stat(density)),
                 alpha=0.5,position = "stack",boundary=0,color="black",
                 binwidth=0.005) +
  scale_fill_manual(values=c("black","grey")) +
  scale_x_continuous(limits=c(-0.36,0.36),
                     breaks=seq(-0.36,36,0.12)) +
  geom_vline(xintercept = 0) +
  ylab("") + xlab("")+
  theme(axis.text.y = element_text(angle = 0)) +
  theme(legend.position=c(0.78,0.87),
        legend.key.size = unit(3, 'lines')) +
  theme(legend.title = element_blank(),
        legend.box.just = "right",
        legend.key = element_rect(size = 6),
        legend.key.size = unit(2, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(file = paste0(path,"/Figures/FigureB3.pdf"), width = 20, height = 12)