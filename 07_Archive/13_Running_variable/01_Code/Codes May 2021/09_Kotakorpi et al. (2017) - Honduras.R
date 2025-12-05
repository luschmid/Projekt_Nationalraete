
suppressMessages(library(tidyverse))
suppressMessages(library(readstata13))

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))
try(setwd("E:/Projekt Nationalräte"))
options(dplyr.summarise.inform = FALSE)

source("./03_Code/13_Running_variable/05_Vote_margin_analytical_largestremainder.R")

#------------------
# Set parameters
#------------------

set.seed(05071975)

#------------------
# Read in datas
#------------------

data_hn <- readstata13::read.dta13("./02_Processed_data/15_Elections_Honduras/elections_hn_final.dta") %>%
  mutate(Sex_num=ifelse(Sex=="H",1,0),
         Department = sub(" ","",Department),
         Department = sub(" ","",Department),
         Department = sub(" ","",Department)) %>%
  rename(year=Year,elected=Elected,votes_h=Votes) %>%
  select(id_Stata,Candidate,year, Department,Party,votes_h,elected) 

data_hn_votes_j <- data_hn %>% group_by(year,Department,Party) %>% 
  dplyr::summarize(votes_j=sum(votes_h,na.rm=T)) %>%
  mutate(party_num=row_number())
data_hn_seats <- data_hn %>% group_by(year,Department) %>%
  dplyr::summarize(seats=sum(elected,na.rm=T),
                   votes_dep=sum(votes_h,na.rm=T))

data_hn <- data_hn %>% left_join(data_hn_votes_j,by=c("Department","Party","year")) %>%
  left_join(data_hn_seats,by=c("Department","year")) %>%
  mutate(votes_h_share=votes_h/votes_dep) %>%
  arrange(Department,year,party_num)

data_hn %>% filter(Department=="ATLANTIDA" & year==2009) %>% select(Candidate,votes_h,votes_h_share,Party)

# Check: do votes_h_share sum up to 1?
check <- data_hn %>% group_by(Department,year) %>% summarize(test=sum(votes_h_share))
summary(check$test)
# Result: yes!

# Functions --------------------------

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
  data_votes_j <- data %>% group_by(year,Department,party_num) %>%
    select(-votes_j) %>%
    dplyr::summarize(votes_j=sum(votes_h,na.rm=T))
  data <- data %>% select(-votes_j) %>% left_join(data_votes_j,by=c("Department","party_num","year")) 
  return(list(data_votes_j,data))
}

GetElected <- function(m=10000,data,control_seat_allocation=0){
  # This function calculates one random election result at the candidate level (if control_seat_allocation=0)
  # or the elected status for all candidates (if control_seat_allocation=1)
  # Input: m: number of votes to sample
  #        data: input data with variables id_Stata,party_num,year,elected, votes_h, votes_j, party_num
  #        control_seat_allocation: indicator that is 1 for check of seat distribution and 0 for random elections
  # Output: dataframe that includes id_Stata,party_num,year,elected
  if (control_seat_allocation==0){
    data$votes_h <- RandomizeElectionResults(data,m=m)
  }
  partyvotes <- GetPartyVotes(data)
  seats <- data.frame(party_num=partyvotes[[1]]$party_num,
                      votes_j=partyvotes[[1]]$votes_j,
                      seats_obtained=as.numeric(GetSeatsHare(votes=partyvotes[[1]]$votes_j,n=data[1,]$seats)))
  data_final <- data %>% left_join(seats,by=c("party_num"))%>%
    group_by(year,Department,party_num)%>%
    mutate(votes_h_rank=rank(-votes_h,ties.method = "random" )) %>%
    ungroup() %>%
    mutate(elected=ifelse(votes_h_rank<=seats_obtained,1,0)) %>%
    select(id_Stata,party_num,year,elected) %>%
    arrange(party_num,id_Stata)
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
    GetElected(data=data,control_seat_allocation=1) %>% mutate(elected_sum=elected,M=1)
  }
  else{
    l <- lapply(rep(m,M),GetElected,data,control_seat_allocation=0) 
    df <- do.call(rbind.data.frame, l) %>%
      group_by(id_Stata,year,party_num) %>%
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
  data_check <- data1 %>% left_join(data2,by=c("id_Stata","year","party_num")) %>%
    select(id_Stata,starts_with("votes_h"),starts_with("elected_sum"),party_num,elected,Candidate)
  party_num_all <- levels(as.factor(data_check$party_num))
  
  data_check_out <- data.frame()
  for (k in 1:length(party_num_all)){
    data_check_party <- data_check %>% filter(party_num==party_num_all[k])
    data_check_party$votes_h_rank <- data.table::frank(data_check_party,-votes_h,-elected,Candidate)  
    data_check_party$elected_sum_rank <- data.table::frank(data_check_party,votes_h_rank)  
    data_check_out <- rbind(data_check_out,data_check_party)
  }
  data_check_out$rank_diff <- data_check_out$votes_h_rank-data_check_out$elected_sum_rank
  check_diff <- summary(data_check_out$rank_diff)
  return(ifelse(abs(check_diff[1])>0 | abs(check_diff[6]),0,1 ))  
}


KotakorpiEtAl2007<- function(m=10000,data,M=20000,control_seat_allocation=0,check_convergence=1){
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
    if (check_convergence==1){
      m=data$seats[1]*20
      M=2000
      df_out_all <- GetElectedSum(m,data,M,control_seat_allocation=0)
      convergence_indi=CheckConvergence(data1=data,data2=df_out_all)
      i=1
      while (convergence_indi!=1){
        print(paste0("Round number: ",i))
        M=10000
        df_out <- GetElectedSum(m,data,M) %>% dplyr::rename(elected_sum_new=elected_sum,
                                                            M_new=M)
        df_out_all <-df_out_all %>% left_join(df_out,by=c("id_Stata","year", "party_num")) %>%
          mutate(M=M+M_new,
                 elected_sum=elected_sum+elected_sum_new)%>%
          select(-ends_with("_new")) %>%
          mutate(p=elected_sum/M)
        convergence_indi=CheckConvergence(data1=data,data2=df_out_all)
        i=i+1
      }
    }
  }
  return(df_out_all )
}

KotakorpiEtAl2007Loop <- function(department,data,control_seat_allocation=0,check_convergence=1) {
  # This function applies KotakorpiEtAl2007 for all years the procedure by one district and more than one year
  # Input: data: original data
  #        control_seat_allocation: indicator that is 1 for check of seat distribution and 0 for random elections
  #        check_convergence: indicator that is 1 if procedure 2 of Kotakorpi et al. 2017 is chosen
  #                           and convergence is checked (Kotakorpi et al. 2017, p. 419)
  # Output: csv file with output of KotakorpiEtAl2007 for all years 
  data_dep <- data %>% filter(Department==department)
  years <- as.character(levels(as.factor(data_dep$year)))
  for (yr in years){
    data_dep_year <- data_dep %>% filter(year==yr)
    print(paste0("Year: ",yr))
    temp <- KotakorpiEtAl2007(data=data_dep_year,control_seat_allocation=control_seat_allocation,check_convergence=check_convergence)
    if (yr==years[1]){
      fwrite(temp, file=paste0("./02_Processed_data/13_Running_variable/kotakorpietal2017",department,".csv"),  append=F)
    }
    else{
      fwrite(temp, file=paste0("./02_Processed_data/13_Running_variable/kotakorpietal2017",department,".csv"),  append=T)
    }
  }
}

CheckProcedure <- function(data){
  deps <- as.character(levels(as.factor(data$Department)))
  lapply(deps,KotakorpiEtAl2007Loop,data,control_seat_allocation=1,check_convergence=0)
  files <-  list.files(path='./02_Processed_data/13_Running_variable',pattern="kotakorpietal2017")
  df <-  do.call(rbind, lapply(paste0("./02_Processed_data/13_Running_variable/",files), function(x) read.csv(x, stringsAsFactors = FALSE))) %>%
    select(id_Stata,year,elected_sum)
  lapply(paste0("./02_Processed_data/13_Running_variable/",files), function(x) file.remove(x))
  
  return(data %>% left_join(df,by=c("id_Stata","year")) %>% select(id_Stata,year,Candidate,elected,elected_sum))
}

# run procedure (from command line with department name as input)

args = commandArgs(trailingOnly=TRUE)
KotakorpiEtAl2007Loop(department=args[1],data=data_hn,control_seat_allocation=0,check_convergence=1) 


# Bind all output files together  --------------------------


files = list.files(path='./02_Processed_data/13_Running_variable',pattern="kota")

# First apply read.csv, then rbind
df = do.call(rbind, lapply(paste0("./02_Processed_data/13_Running_variable/",files), function(x) read.csv(x, stringsAsFactors = FALSE))) %>%
  mutate(p=elected_sum/M)
summary(df$M)

data.table::fwrite(df, file=paste0("./02_Processed_data/13_Running_variable/kotakorpietal2017_honduras.csv"),  append=F)



# test functions
# 
# data_test <- data_hn %>% filter(year==2009 & Department=="ATLANTIDA")
# data_test$votes_h <- RandomizeElectionResults(m=100,data=data_test)
# partyvotes_test <- GetPartyVotes(data_test)
# 
# seats_test <- data.frame(party_num=partyvotes_test[[1]]$party_num,
#                          seats_obtained=as.numeric(GetSeatsHare(votes=partyvotes_test[[1]]$votes_j,n=data_test[1,]$seats)))
# 
# data_final_test <- data_test %>% left_join(seats_test,by=c("party_num"))%>%
#                                  group_by(year,Department,party_num)%>%
#                                  mutate(votes_h_rank=rank(-votes_h,ties.method = "random" )) %>%
#                                  ungroup() %>%
#                                  mutate(elected=ifelse(votes_h_rank<=seats_obtained,1,0)) %>%
#                                  select(id_Stata,year,Department,party_num,votes_h,votes_h_rank,elected,seats_obtained) %>%
#                                  arrange(year,Department,party_num,votes_h,)
# 
# test <- RandomizeSeats(m=1000,data=data_test)
# sum(test$elected)
# 
# out_kota <- KotakorpiEtAl2007(data=data_test,control_seat_allocation=1)

#deps <- as.character(levels(as.factor(data_hn$Department)))
# 
# check whether we can replicate seat allocation
# df_out <- CheckProcedure(data_hn) %>% mutate(diff=elected-elected_sum)
# summary(df_out$diff)
# result: no difference in real election status and election status after our procedure

