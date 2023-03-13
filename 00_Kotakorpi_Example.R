
# 2. Functions

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

GetNumberSeats <- function(data){
# This function calculates the number of seats of an alliance/suballianc/party based on D'Hondt numbers
# Input: data with party name and party votes
  DHondtNumbers <- CalculateDHondtNumbers(data,no_seats=max(data$seats))
  data_out <- DHondtNumbers %>%  
    left_join(data,by=c("party"))%>% 
    mutate(dHondtno_rank=rank(-as.numeric(dHondtno),ties.method = "random" )) %>%
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

CheckConvergence <- function(data1,data2){
  # This function checks whether the election probability has converged, meaning
  # that within party j, the ranking of p and the ranking of votes_h are equal 
  # (see Kotakorpi et al. 2017, footnote 3)
  # Input: data1: original data
  #        data2: data after M random votes of size m have been taken
  #               variable elected_sum is the number of times candidate h is elected
  # Output: 1 if convergence was achieved, 0 otherwise
  data_check <- data1 %>% left_join(data2,by=c("ID","party_num")) %>%
    select(ID,starts_with("votes_h"),starts_with("elected_sum"),party_num,elected)
  party_num_all <- levels(as.factor(data_check$party_num))
  
  data_check_out <- data.frame()
  for (k in 1:length(party_num_all)){
    data_check_party <- data_check %>% filter(party_num==party_num_all[k])
    data_check_party$votes_h_rank <- data.table::frank(data_check_party,-votes_h,-elected)  
    data_check_party$elected_sum_rank <- data.table::frank(data_check_party,votes_h_rank)  
    data_check_out <- rbind(data_check_out,data_check_party)
  }
  data_check_out$rank_diff <- data_check_out$votes_h_rank-data_check_out$elected_sum_rank
  check_diff <- summary(data_check_out$rank_diff)
  return(ifelse(abs(check_diff[1])>0 | abs(check_diff[6]),0,1 ))  
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
  return(df_out_all )
}

KotakorpiEtAl2017Loop<- function(m_min=1,m_max=20000,data,M=20000,control_seat_allocation=0,check_convergence=1){
  # This function loops over KotakorpiEtAl for different m the procedure by Kotakorpi et al. 2017 for one district in a year
for (i in m_min:m_max){
out_data <- KotakorpiEtAl2017(m=i,data=df,M=M,control_seat_allocation=control_seat_allocation,check_convergence=check_convergence)
out_data$m <- i
if (i==m_min){
  fwrite(out_data, file=paste0("./02_Processed_data/13_Running_variable/03_kotakorpietal2017_example/kotakorpietal2017_example.csv"),  append=F)
}
else{
  fwrite(out_data, file=paste0("./02_Processed_data/13_Running_variable/03_kotakorpietal2017_example/kotakorpietal2017_example.csv"),  append=T)
}
}
}

CheckProcedure <- function(data){
  cans <- as.character(levels(as.factor(data$canton)))
  lapply(cans,KotakorpiEtAl2017Loop,data,control_seat_allocation=1,check_convergence=0)
  files <-  list.files(path='./02_Processed_data/13_Running_variable',pattern="kotakorpietal2017_ch",full.names = F, recursive = TRUE)
  df <-  do.call(rbind, lapply(paste0("./02_Processed_data/13_Running_variable/",files), function(x) read.csv(x, stringsAsFactors = FALSE))) %>%
    select(ID,elected_sum)
  lapply(paste0("./02_Processed_data/13_Running_variable/",files), function(x) file.remove(x))
  
  return(data %>% left_join(df,by=c("ID")) %>% select(ID,elected,elected_sum))
}
