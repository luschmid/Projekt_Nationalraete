rm(list = ls())
setwd("D:/SchmidLu/Dropbox/Projekt Nationalräte")
library(dplyr)
library(tidyr)
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
    quota <- sum(votes) / n
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
  return( GetRatioHare(votes,n) - GetIntegerHare(votes, n))
}

GetHighestRemainder <- function(remainder, freeseats) {
  return(ifelse(rank(-remainder, ties.method = "average") <= freeseats, 1, 0))
}

GetSeatsHare <- function(votes, n) {
  Seats_Round1 <- GetIntegerHare(votes, n)
  Remainders_Round1 <- GetRemainderHare(votes, n)
  freeseats <- n - sum(Seats_Round1)
  Seats_Round2 <- GetHighestRemainder(Remainders_Round1, freeseats)
  return(Seats_Round1 + Seats_Round2)
}



GetVotesReq <- function(votes, n, j, minusj) {
  votes_others <- sum(votes[row(as.matrix(votes)) != j])
  for (deltai in (-n:(n-1))) {
    votes_req_j <- (votes[minusj] * n + deltai * votes_others) / (n - deltai)
    votes_mod <- votes
    votes_mod[j] <- votes_req_j
    i <- GetIntegerHare(votes_mod,n)[j]
    # print(paste("votes req:",votes_req_j,sep=" "))
    if (deltai == (-n)) {
      votes_req_all <- votes_req_j
      i_all <- i 
    }
    else {
      votes_req_all <- append(votes_req_all, votes_req_j)
      i_all <- append(i_all,i)
    }
  }
  out_df <- data.frame(deltai = c(-n:(n-1)), j, minusj, votes_req_all,i=i_all)
  return(out_df)
}


GetVotesReqAll <- function(votes, n) {
  for (j in 1:length(votes)) {
    for (minusj in c(1:length(votes))[-j]) {
      VRi <- GetVotesReq(votes, n, j, minusj)
      if (j == 1 & minusj == c(1:length(votes))[-j][1]) {
        VRi_all <- VRi
      } else {
        VRi_all <- rbind(VRi_all, VRi)
      }
    }
  }
  return(VRi_all)
}




GetRankFreeSeats <- function(votes, n, VR) {
  I <- c(1:length(votes))
  J <- c(0:(n-1))
  VR_out <- VR %>%  expand(I, I, J)
  colnames(VR_out) <- c("j", "minusj", "i")
  VR_out <- VR_out %>%
    filter(j != minusj) %>%
    left_join(VR, by = c("j", "minusj", "i")) 
  VR_out[["rank"]] <- NA
  VR_out[["freeseats"]] <- NA
  VR_out[["remainderj"]] <- NA
  VR_out[["remainderminusj"]] <- NA
  for (u in 1:dim(VR_out)[1]) {
    minusj <- VR_out[["minusj"]][u]
    j <- VR_out[["j"]][u]
    votes_mod <- votes
    votes_mod[j] <- VR_out[["votes_req_all"]][u]
    minusj_revelant <- ifelse(j>minusj,minusj,minusj-1)
    VR_out[["freeseats"]][u] <- n - sum(GetIntegerHare(votes_mod, n))
    VR_out[["rank"]][u] <- rank(-GetRemainderHare(votes_mod, n)[-j])[minusj_revelant]
    VR_out[["remainderj"]][u] <- GetRatioHare(votes_mod, n)[j]-floor(GetRatioHare(votes_mod, n)[j])
    VR_out[["remainderminusj"]][u] <- GetRatioHare(votes_mod, n)[minusj]-floor(GetRatioHare(votes_mod, n)[minusj])
  }
  return(VR_out)
}




GetVSHare <- function(votes, n) {
  VRi_all <- GetVotesReqAll(votes, n) 
  votes_df <- data.frame(votes = votes, j = c(1:length(votes)))
  VR <- GetRankFreeSeats(votes, n, VRi_all) %>%
    filter(rank == freeseats ) %>%
    mutate(remainder_diff=remainderj-remainderminusj) %>%
    filter(abs(remainder_diff) <=0.1 ) %>% # this is a workaround b/c of the rounding procedure of R
    left_join(votes_df, by = c("j")) %>%
    mutate(VS = votes - votes_req_all) %>%
    rename(party = j, seat = i) %>%
    mutate(seat=seat+1) %>%
    select(party, seat, VS) %>%
    arrange(party, seat)
  return(VR)
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


CompareAnalayticalSimulationAll <- function(nparties_max, nseats_max, niterations, seedstart = 1, print_every_iteration = 100) {
  for (np in 2:nparties_max) {
    print(paste("Number of parties: ", np))
    #for (ns in 1:nseats_max) {
    for (ns in c(1:30,seq(40,100,10))) {
        print(paste("Number of seats: ", ns))
      Out<- CompareAnalayticalSimulation(nparties = np, nseats = ns, niterations = 100, seedstart = 1200, print_every_iteration = print_every_iteration)
      Out <- data.frame(cbind(Out,matrix(NA,niterations,nparties_max-np)))
      colnames(Out) <-  c("nparties", "nseats","iteration","difference",paste("votes", 1:nparties_max, sep=""))
      fwrite(Out, file="./02_Processed_data/13_Running_variable/Simulation_RV_LargestRemainder.csv",  append=T)
      }
  }
}


# test functions
votes <- c(47000, 16000, 15800, 12000, 6100, 3100)
i <- 2
j <- 1
minusj <- 2
n <- 10
GetVotesReq(votes, n, j, minusj)
V1 <- GetVotesReqAll(votes, n)
GetRankFreeSeats(votes, n, V1)
GetVSHare(votes, n)




votes <- c(35,45,20); n=3
GetSeatsHare(votes, n)
GetIntegerHare(votes, n)
GetRemainderHare(votes, n)
votes_req <- GetVotesReqAll(votes, n) %>% mutate(j=factor(j,labels=c("A","B","C")),
                                                 minusj=factor(j,labels=c("B","C","A")))
print(xtable(votes_req[c(1:10),],digits=1),include.rownames=FALSE)

votes_sol <- GetRankFreeSeats(votes,n,GetVotesReqAll(votes, n) ) %>% 
mutate(j=factor(j,labels=c("A","B","C")),
       minusj=factor(j,labels=c("B","C","A"))) %>% 
       select(deltai,i,itilde,j,minusj,votes_req_all,rank,freeseats,i_new)

print(xtable(votes_sol[c(1:10),],digits=c(0,0,0,0,0,0,1,0,0,0)),include.rownames=FALSE)


GetVSHare(votes, n)
votes <- c(35,45,20); n=3
GetRatioHare(votes, n)
print(xtable(V[c(1:10),],digits=1),include.rownames=FALSE)


GetSeatsHare(votes, n)
GetIntegerHare(votes, n)
GetRemainderHare(votes, n)


# test analytical solution vs. simulation

nseats <- 10
votes <- c(47000, 16000, 15800, 12000, 6100, 3100)
OutSimulation <- GetRVSimulation(votes, nseats, 0.000001)
OutSimulation <- OutSimulation %>% rename(VS_Sim = VS)
OutAnalytical <- GetVSHare(votes, nseats) %>% rename(VS_Ana = VS)
OutCompare <- OutAnalytical %>% full_join(OutSimulation, by = c("party", "seat"))
check <- max(OutCompare$VS_Ana - OutCompare$VS_Sim)

max(abs(OutCompare$VS_Ana - OutCompare$VS_Sim))







votes <- c(60728,85355)
nseats <- 7
OutSimulation <- GetRVSimulation(votes, nseats, 0.000001)
OutSimulation <- OutSimulation %>% rename(VS_Sim = VS)
OutAnalytical <- GetVSHare(votes, nseats) %>% rename(VS_Ana = VS)
OutCompare <- OutAnalytical %>% full_join(OutSimulation, by = c("party", "seat"))
(check <- max(OutCompare$VS_Ana - OutCompare$VS_Sim))

CompareAnalayticalSimulationAll(nparties_max = 20, nseats_max = 20, niterations = 1000, seedstart = 1200, print_every_iteration = 1000)


CompareAnalayticalSimulationAllDifference <- read.csv(file="./02_Processed_data/13_Running_variable/Simulation_RV_LargestRemainder.csv",sep=",") %>%
                                             distinct(votes1,votes2,votes3,votes4,votes5,votes6,votes7,votes8,votes9,votes10,
                                                      votes11,votes12,votes13,votes14,votes15,votes16,votes17,votes18,votes19,votes20,.keep_all = T) %>%
                                             filter(difference==1)


nseats <- 22
votes <- c(16808,26037,5786)
OutSimulation <- GetRVSimulation(votes, nseats, 0.000001)
OutSimulation <- OutSimulation %>% rename(VS_Sim = VS)
OutAnalytical <- GetVSHare(votes, nseats) %>% rename(VS_Ana = VS)
OutCompare <- OutAnalytical %>% full_join(OutSimulation, by = c("party", "seat")) %>% print(n=100) %>%
(check <- max(OutCompare$VS_Ana - OutCompare$VS_Sim))
