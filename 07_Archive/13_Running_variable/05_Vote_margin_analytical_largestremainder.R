rm(list = ls())
library(dplyr)

CheckNA <- function(x) {
  if (any(is.na(x))) {
    print("The vector of votes includes missing values.")
    return(1)
  }
  else {
    return(0)
  }
}

GetIntegerHare <- function(votes, n) {
  if (CheckNA(votes) == 0) {
    quota <- sum(votes) / n
    return(floor(votes / quota))
  }
}

GetRemainderHare <- function(votes, n) {
  IntegerHare <- GetIntegerHare(votes, n)
  quota <- sum(votes) / n
  return((votes / quota) - IntegerHare)
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



GetVotesReq <- function(votes, i, j, n, minusj) {
  votes_others <- sum(votes[row(as.matrix(votes)) != j])
  for (itilde in 0:(n - i)) {
    votes_req_j <- (votes[minusj] * n + (i - itilde) * votes_others) / (n - i + itilde)
    # print(paste("votes req:",votes_req_j,sep=" "))
    if (itilde == 0) {
      votes_req_all <- votes_req_j
    }
    else {
      votes_req_all <- append(votes_req_all, votes_req_j)
    }
  }
  out_df <- data.frame(i, itilde = c(0:(n - i)), j, minusj, votes_req_all)
  return(out_df)
}



GetVotesReqAll <- function(votes, n) {
  for (j in 1:length(votes)) {
    for (i in 1:n) {
      for (minusj in c(1:length(votes))[-j]) {
        VSi <- GetVotesReq(votes, i, j, n, minusj)
        if (j == 1 & i == 1 & minusj == c(1:length(votes))[-j][1]) {
          VSi_all <- VSi
        } else {
          VSi_all <- rbind(VSi_all, VSi)
        }
      }
    }
  }
  return(VSi_all)
}


GetRankFreeSeats <- function(votes, VS, n) {
  VS[["rank"]] <- NA
  VS[["freeseats"]] <- NA
  VS[["i_new"]] <- NA
  for (u in 1:dim(VS)[1]) {
    minusj <- VS[["minusj"]][u]
    j <- VS[["j"]][u]
    votes_mod <- votes
    votes_mod[j] <- VS[["votes_req_all"]][u] - 0.0000001
    VS[["freeseats"]][u] <- n - sum(GetIntegerHare(votes_mod, n))
    VS[["rank"]][u] <- rank(-GetRemainderHare(votes_mod, n))[minusj]
    VS[["i_new"]][u] <- GetIntegerHare(votes_mod, n)[j] + 1
  }
  return(VS)
}

GetVSHare <- function(votes, n) {
  VSi_all <- GetVotesReqAll(votes, n) %>% filter(votes_req_all >= 0)
  votes_df <- data.frame(votes = votes, j = c(1:length(votes)))
  VS <- GetRankFreeSeats(votes, VSi_all, n) %>%
    filter(rank <= freeseats) %>%
    group_by(j, i_new) %>%
    summarize(votes_req_all_min = min(votes_req_all)) %>%
    left_join(votes_df, by = c("j")) %>%
    mutate(VS = votes - votes_req_all_min) %>%
    rename(party = j, seat = i_new) %>%
    select(-votes_req_all_min)
  return(VS)
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


CompareAnalayticalSimulation <- function(nparties, nseats, niterations, seedstart = 1) {
  it <- seedstart
  out <- data.frame()
  check <- 0
  while (abs(check) <= 0.1 & it < (niterations + seedstart)) {
    set.seed(it)
    if (it %% 100 == 0) {
      print(paste("Iteration:", it))
    }
    votes <- round(runif(nparties, 0, 100000))
    OutSimulation <- GetRVSimulation(votes, nseats, 0.000001)
    OutAnalytical <- GetVSHare(votes, nseats)
    check <- max(OutSimulation$VS - OutAnalytical$VS)
    it <- it + 1
  }
  if (it < (niterations + seedstart)) {
    return(as.vector(votes))
  }
}


nseats <- 10
votes <- c(47000, 16000, 15800, 12000, 6100, 3100)
OutSimulation <- GetRVSimulation(votes, nseats, 0.000001)
OutSimulation <- OutSimulation %>% rename(VS_Sim = VS)
OutAnalytical <- GetVSHare(votes, nseats)
OutAnalytical <- OutAnalytical %>% rename(VS_Ana = VS)
(OutCompare <- OutAnalytical %>% left_join(OutSimulation, by = c("party", "seat")))

max(abs(OutCompare$VS_Ana - OutCompare$VS_Sim))

(Out2_partyvotes <- CompareAnalayticalSimulation(nparties = 4, nseats = 10, niterations = 1000, seedstart = 1200))


# test overall function
votes <- c(47000, 16000, 15800, 12000, 6100, 3100)
i=2;j=1;jminus=2;n=10
GetVotesReq(votes, i, j, n, jminus)

# GetVSHare(votes,n)
#

#
# # test separate functions
VSi_all <- GetVotesReqAll(votes,n)
# head(VSi_all)
#
# library(tidyverse)
# VSi_all %>% filter(j==1 & minusj==2) %>% distinct(votes_req_all,.keep_all=T)
# VSi_all_positive <- VSi_all %>% filter(votes_req_all>=0) # %>% distinct(votes_req_all,.keep_all=T)
# #head(VSi_all_positive)
#
# VSi_all_out <- GetRankFreeSeats(votes,VSi_all_positive,n)
# head(VSi_all_out)
#
# VSi_all_out_possible <- VSi_all_out %>% filter(rank<=freeseats)
# head(VSi_all_out_possible) %>% as_tibble %>% filter(j==1 & i==1)  %>% print(n=100)
#
#
# VSi_all_min <- VSi_all_out_possible %>% group_by(j,i_new) %>%
#       summarize(votes_req_all_min=min(votes_req_all))
#
#
# VSi_all_min %>% as_tibble()    %>% print(n=100)
# VSi_all_out %>% filter(j==5 & i_new==2)
