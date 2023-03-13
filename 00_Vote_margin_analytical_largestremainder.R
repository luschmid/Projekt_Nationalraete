#rm(list = ls())
#try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))
#setwd("//Srw-hpc5/e//Projekt Nationalräte")
#setwd("C:/Users/lukas/Dropbox/Projekt Nationalräte")

library(dplyr)
library(tidyr)
library(Rmpfr)
library(data.table)


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




GetVSHare <- function(votes, n) {
  VR <- GetVotesReqAll(votes, n) 
  VR <- VR[which(VR[,6]==VR[,7]  & abs(VR[,8]-VR[,9])<0.01 ),]  
  VR_df <- data.frame(votes_req = as.numeric(VR[,1]),party= as.numeric(VR[,3]),seat=as.numeric(VR[,2]))
  votes_df <- data.frame(votes = votes, party= c(1:length(votes)))
  VR_out <- VR_df %>%
    left_join(votes_df, by = c("party")) %>%
    mutate(VS = votes - as.numeric(votes_req)) %>%
    arrange(party, seat) %>%
    filter(seat>0)
  return(VR_out)
}




GetRVSimulation <- function(votes, n, convcrit) {
  out <- data.frame()
  for (j in 1:length(votes)) {
    # print(paste("Party: ",j))
    for (i in 1:n) {
      # print(paste("Seat: ",i))
      votessim <- votes
      votessim[j] <- 0.5 # start with 1 vote in first iteration (first line after while loop)
      # print(votessim)
      seats <- 0
      while (seats < i) {
        votessim[j] <- votessim[j] * 2
        seats <- GetSeatsHare(votessim, n)[j]
        # print(votessim)
      }
      xlow <- votessim[j] / 2
      xhigh <- votessim[j]
      # print(paste("Seat change at: ",votessim[j]))

      while ((xhigh - xlow) > convcrit) {
        votessim[j] <- (xlow + xhigh) / 2
        seats <- GetSeatsHare(votessim, n)[j]
        if (seats < i) {
          xlow <- votessim[j]
        }
        if (seats >= i) {
          xhigh <- votessim[j]
        }
        # print(paste("xhigh: ",xhigh))
        # print(paste("xlow: ",xlow))
      }
      out_ji <- data.frame(party = j, seat = i, VS = round(votes[j] - votessim[j], 4))
      out <- rbind(out, out_ji)
    }
  }
  return(out)
}


CompareAnalayticalSimulation <- function(nparties, nseats, niterations, seedstart = 1, print_every_iteration = 100) {
  it <- seedstart
  out <- data.frame()
  check <- 0
  set.seed(seedstart)
  while (it < (niterations + seedstart)) {
    if (it %% print_every_iteration == 0 & it>seedstart) {
      print(paste("Iteration:", it-seedstart+1))
    }
    votes <- round(runif(nparties, 0, 100000))
    OutSimulation <- GetRVSimulation(votes, nseats, 0.000001)
    OutAnalytical <- GetVSHare(votes, nseats)
    OutCompare <- OutAnalytical %>% full_join(OutSimulation, by = c("party", "seat"))
    check <- max(OutCompare$VS.x - OutCompare$VS.y)
    if (check>0.1 | length(OutCompare$VS.x[is.na(OutCompare$VS.x)]==T)>0 | length(OutCompare$VS.x[is.na(OutCompare$VS.y)]==T)>0){
        difference <- 1
        }
    else{difference <- 0}
    if (it==seedstart) {results <- c(nparties, nseats,(it-seedstart+1),difference,votes)}
    else{results <- rbind(results,c(nparties, nseats,(it-seedstart+1),difference,votes))}
    it <- it + 1
  }
  df_results <- as.data.frame(results)
  colnames(df_results) <- c("nparties", "nseats","iteration","difference",paste("votes", 1:nparties, sep=""))
  rownames(df_results) <- c(1:niterations)
  return(df_results)
}


CompareAnalayticalSimulationAll <- function(nparties_max, seats, niterations, print_every_iteration = 100) {
  for (np in 2:nparties_max) {
    print(paste("Number of parties: ", np))
    #for (ns in 1:nseats_max) {
    for (ns in seats) {
      print(paste("Number of seats: ", ns))
      seedstart_loop <- (np-1)*length(seats)+which(seats==ns)
      Out<- CompareAnalayticalSimulation(nparties = np, nseats = ns, niterations = 100, seedstart = seedstart_loop, print_every_iteration = print_every_iteration)
      Out <- data.frame(cbind(Out,matrix(NA,niterations,nparties_max-np)))
      colnames(Out) <-  c("nparties", "nseats","iteration","difference",paste("votes", 1:nparties_max, sep=""))
      fwrite(Out, file="./02_Processed_data/13_Running_variable/Simulation_RV_LargestRemainder.csv",  append=T)
      }
  }
}

# 
# # test functions
# votes <- c(47000, 16000, 15800, 12000, 6100, 3100)
# i <- 2
# j <- 1
# minusj <- 2
# n <- 10
# GetVotesReq( j, minusj,votes, n)
# 
# V1 <- GetVotesReqAll(votes, n)
# GetVSHare(votes, n)
# 
# 
# # test functions for example in paper
# votes <- c(45,35,15)
# i <- 1
# j <- 2
# minusj <- 2
# n <- 10
# GetVotesReq(votes, n, j, minusj)
# V1 <- GetVotesReqAll(votes, n)
# GetVSHare(votes, n)
# 
# # presentation
# 
# library(xtable)
# 
# votes <- c(35,45,20); n=3
# GetSeatsHare(votes, n)
# GetIntegerHare(votes, n)
# GetRemainderHare(votes, n)
# 
# votes_sol <- GetVotesReqAll(votes, n)
# votes_sol_colnames=colnames(votes_sol)
# votes_sol_df <- as.data.frame(matrix(as.numeric(votes_sol),ncol=dim(votes_sol)[2],nrow=dim(votes_sol)[1]))  
# colnames(votes_sol_df) <- votes_sol_colnames
# votes_sol_df <- votes_sol_df %>%
#                 mutate(j=factor(j,labels=c("A","B","C")),
#                        minusj=factor(minusj,labels=c("A","B","C")), 
#                        i_floor=i-1,
#                        deltai=deltai-1) %>% 
#                 select(deltai,j,minusj,votes_req_all,rank_minusj,freeseats,i,i_floor) %>%
#                 arrange(j,minusj,deltai)
# 
# print(xtable(votes_sol_df[c(1:12),c(1:4,7)],digits=c(0,0,0,0,1,0)),include.rownames=FALSE)
# print(xtable(votes_sol_df,digits=c(0,0,0,0,1,0,0,0,0)),include.rownames=FALSE)
# 
# print(xtable(GetVSHare(votes, n) %>% 
#                mutate(party=factor(party,labels=c("A","B","C"))) %>% 
#                select(party,seat,votes,votes_req,VS),digits=c(0,0,0,0,1,1)),include.rownames=FALSE)
# 
# 
# 
# GetRatioHare(c(45,45,20), n)
# votes <- c(35,45,20); n=3
# GetRatioHare(votes, n)
# print(xtable(V[c(1:10),],digits=1),include.rownames=FALSE)
# 
# 
# # new example with two parties
# 
# GetSeatsHare(votes, n)
# GetIntegerHare(votes, n)
# GetRemainderHare(votes, n)
# 
# 
# votes <- c(60,40); n=2
# GetSeatsHare(votes, n)
# GetIntegerHare(votes, n)
# (VSHare <- GetVSHare(votes, n))
# 
# print(xtable(VSHare,digits=1),include.rownames=FALSE)
# 
# 
# 
# # test analytical solution vs. simulation
# 
# nseats <- 10 # wikipedia example
# votes <- c(47000, 16000, 15800, 12000, 6100, 3100)
# OutSimulation <- GetRVSimulation(votes, nseats, 0.000001)
# OutSimulation <- OutSimulation %>% rename(VS_Sim = VS)
# OutAnalytical <- GetVSHare(votes, nseats) %>% rename(VS_Ana = VS)
# OutCompare <- OutAnalytical %>% full_join(OutSimulation, by = c("party", "seat"))
# check <- max(OutCompare$VS_Ana - OutCompare$VS_Sim)
# 
# max(abs(OutCompare$VS_Ana - OutCompare$VS_Sim))
# 
# CompareAnalayticalSimulationAll(nparties_max = 20, seats=c(1:30,seq(40,100,10)), niterations = 1000, print_every_iteration = 1000)
# 
# 
# CompareAnalayticalSimulationAllDifference <- read.csv(file="./02_Processed_data/13_Running_variable/Simulation_RV_LargestRemainder.csv",sep=",") %>%
#                                              distinct(votes1,votes2,votes3,votes4,votes5,votes6,votes7,votes8,votes9,votes10,
#                                                       votes11,votes12,votes13,votes14,votes15,votes16,votes17,votes18,votes19,votes20,.keep_all = T) %>%
#                                              filter(difference==1)
# 
# 


# votes <- c(60728,85355)
# nseats <- 7
# OutSimulation <- GetRVSimulation(votes, nseats, 0.000001)
# OutSimulation <- OutSimulation %>% rename(VS_Sim = VS)
# OutAnalytical <- GetVSHare(votes, nseats) %>% rename(VS_Ana = VS)
# OutCompare <- OutAnalytical %>% full_join(OutSimulation, by = c("party", "seat"))
# (check <- max(OutCompare$VS_Ana - OutCompare$VS_Sim))


# nseats <- 22
# votes <- c(16808,26037,5786)
# OutSimulation <- GetRVSimulation(votes, nseats, 0.000001)
# OutSimulation <- OutSimulation %>% rename(VS_Sim = VS)
# OutAnalytical <- GetVSHare(votes, nseats) %>% rename(VS_Ana = VS)
# OutCompare <- OutAnalytical %>% full_join(OutSimulation, by = c("party", "seat")) %>% print(n=100) %>%
# (check <- max(OutCompare$VS_Ana - OutCompare$VS_Sim))
