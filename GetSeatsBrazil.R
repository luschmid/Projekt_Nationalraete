

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
