
# (i) Sample data (https://en.wikipedia.org/wiki/Largest_remainder_method)

partyvotes <- c(47000, 16000, 15800, 12000, 6100, 3100)


# (ii) Calculate integer and remainder part for seat distribution

CheckNA <- function(x){
if(any(is.na(x))){print("The vector of votes includes missing values.")
return(1)
}
else {return(0)}
}  

GetIntegerHare<- function(votes,n){
if (CheckNA(votes)==0){
quota<-sum(votes)/n
return(floor(votes/quota))  
}
}

GetIntegerHare(partyvotes,10)

GetRemainderHare<- function(votes,n){
IntegerHare <- GetIntegerHare(votes,n)
quota<-sum(votes)/n
return((votes/quota)-IntegerHare)
}

GetHighestRemainder<- function(remainder,freeseats){
return(ifelse(rank(-remainder,ties.method="average")<=freeseats,1,0))
}

GetSeatsHare<- function(votes,n){
Seats_Round1 <- GetIntegerHare(votes,n) 
Remainders_Round1 <- GetRemainderHare(votes,n)
freeseats <- n-sum(Seats_Round1)
Seats_Round2 <-GetHighestRemainder(Remainders_Round1,freeseats)
return(Seats_Round1+Seats_Round2)
}

GetNumberFreeSeats <- function(votes,j,i,n){
seatsround1 <- floor(votes*(n-i+1)/sum(votes[row(as.matrix(votes))!=j]))
return(n-(i-1)-sum(seatsround1[row(as.matrix(seatsround1))!=j]))
}


GetVotesPlusVS1 <- function(votes,j,i,n){
out <- ((i-1)/(n-i+1))*sum(votes[row(as.matrix(votes))!=j])
return(out)
}

GetVotesPlusVS1Corrected <- function(votes,j,i,n){
votes[j] <- GetVotesPlusVS1(votes,j,i,n)
freeseats<- GetNumberFreeSeats(votes,j,i,n)
if (freeseats>0){  
out <- ((i-1)/(n-i+1))*sum(votes[row(as.matrix(votes))!=j])
}
if (freeseats==0){  
out <- ((i-1)/(n-i+1))*sum(votes[row(as.matrix(votes))!=j])+0.001
}
return(out)
}

GetVS2minusj <- function(votes,j,i,n,minusj,ishifted){
votes_others <- sum(votes[row(as.matrix(votes))!=j])
votes_minusj <- votes[row(as.matrix(votes))==minusj]
freeseats<- GetNumberFreeSeats(votes,j,i,n)
#print(paste("No free seats:",freeseats),sep=" ")
if (freeseats>0){  
out <- (votes_minusj*n-floor(votes_minusj*(n-(i-1))/votes_others)*(n/(n-(i-1)))*votes_others)/(n-(i-1)+floor(votes_minusj*(n-(i-1))/votes_others))
}
if (freeseats==0){  
out <- (votes_minusj*n-floor(votes_minusj*(n-(i-1))/(votes_others+(n-(i-1))/n*0.001))*(n/(n-(i-1)))*(votes_others+0.001)-(n-(i-1))*0.001)/(n-(i-1)+floor(votes_minusj*(n-(i-1))/(votes_others+(n-(i-1))/n*0.001)))
}
#if (length(votes)==2){out <- ceiling(votes[j]-(2*i-1)/(2*n-2*i+1)*votes[row(as.matrix(votes))!=j])}
#if (length(ishifted)>0){out <- (votes_minusj*n-floor(votes_minusj*(n-(ishifted-1))/votes_others)*(n/(n-(i-1)))*votes_others)/(n-(i-1)+floor(votes_minusj*(n-(ishifted-1))/votes_others))} # this is for the case when total seats change
return(out)
}


GetVS2s <- function(votes,j,i,n,ishifted=i){
minusj_all <- row(as.matrix(votes))[c(row(as.matrix(votes))!=j)]
VS2_all <- mapply(GetVS2minusj,minusj=c(minusj_all),MoreArgs=list(votes=votes,j=j,i=i,n=n,ishifted=ishifted))
return(VS2_all)
}


GetRatioRound2<- function(votes,j,i,n,VS2){
votes[j] <- GetVotesPlusVS1Corrected(votes,j,i,n)+VS2
return(votes*n/sum(votes))
}

GetIntegerRound2<- function(votes,j,i,n,VS2){
return(floor(GetRatioRound2(votes,j,i,n,VS2)))
}

GetRemainderRound2<- function(votes,j,i,n,VS2){
return(GetRatioRound2(votes,j,i,n,VS2)-floor(GetRatioRound2(votes,j,i,n,VS2)))
}

GetFreeSeatsRound2<- function(votes,j,i,n,VS2){
return(n-sum(GetIntegerRound2(votes,j,i,n,VS2)))
}

GetRemainderRankRound2<- function(votes,j,i,n,VS2){
return(rank(-GetRemainderRound2(votes,j,i,n,VS2)[row(as.matrix(votes))!=j]))
}


GetRemainderRound2s<- function(votes,j,i,n){
VS2s <- GetVS2s(votes,j,i,n)
RemaindersAll <- mapply(GetRemainderRound2,VS2=c(VS2s),MoreArgs=list(votes=votes,j=j,i=i,n=n))
return(RemaindersAll)
}



GetRemainderRankRound2s<- function(votes,j,i,n){
VS2s <- GetVS2s(votes,j,i,n)
RemainderRanksAll <- mapply(GetRemainderRankRound2,VS2=c(VS2s),MoreArgs=list(votes=votes,j=j,i=i,n=n))
return(RemainderRanksAll)
}

GetFreeSeatsRound2s<- function(votes,j,i,n){
VS2s <- GetVS2s(votes,j,i,n)
FreeSeats <- mapply(GetFreeSeatsRound2,VS2=c(VS2s),MoreArgs=list(votes=votes,j=j,i=i,n=n))
return(FreeSeats)
}

GetRelevantRemainderRank<- function(votes,j,i,n){
Remainderanks <- GetRemainderRankRound2s(votes,j,i,n)
return(diag(Remainderanks))
}

GetRatioRound2(votes,j,i,n,7692.31)
GetIntegerRound2(votes,j,i,n,7692.31)
GetRemainderRound2(votes,j,i,n,7692.31)
GetFreeSeatsRound2(votes,j,i,n,7692.31)
GetFreeSeatsRound2s(votes,j,i,n)
GetRemainderRound2s(votes,j,i,n)
GetRemainderRankRound2s(votes,j,i,n)
GetVS2s(votes,j,i,n)

GetRelevantVS2<- function(votes,j,i,n){
Rank_minusj <- diag(GetRemainderRankRound2s(votes,j,i,n))
FreeSeatsRound2<- GetFreeSeatsRound2s(votes,j,i,n)
VS2_All <- GetVS2s(votes,j,i,n)
VS2_Final <- VS2_All[Rank_minusj==FreeSeatsRound2]
return(VS2_Final)
}

GetVS <- function(votes,j,i,n){
out <- votes[row(as.matrix(votes))==j]-(GetVotesPlusVS1Corrected(votes,j,i,n)+GetRelevantVS2(votes,j,i,n))
return(out)
}


GetRelevantVS2_OLD<- function(votes,j,i,n){
if(length(votes)==2){ # case with only two parties
  out <- GetVS2s(votes,j,i,n)
}

if(length(votes)>2){ # case with more than two parties
FreeSeats <- GetFreeSeatsRound2s(votes,j,i,n)
votes_afterround1 <- votes
votes_afterround1[j] <- GetVotesPlusVS1(votes,j,i,n)
SeatsRound1 <- GetIntegerHare(votes,n)
SeatsRound2<- GetIntegerHare(votes_afterround1,n)
RemainderRank <- GetRelevantRemainderRank(votes,j,i,n)
VS2s <- GetVS2s(votes,j,i,n)
relevant_j <- row(as.matrix(FreeSeats))[FreeSeats==(length(votes)-RemainderRank)] # note: we do not substract 1 on the right side of the expression in brackets b/c votes has length J and Free Seats has length (J-1)
out <- VS2s[relevant_j]
if(length(out)==0 & length(which(FreeSeats>(length(votes-1)-RemainderRank)))>0){ # case for which party j is better ranked than necessary to achieve seat
out <-min(VS2s[which(FreeSeats>(length(votes-1)-RemainderRank))])
}
votes_aftervs2 <- votes
votes_aftervs2[j] <- ifelse(length(out)>0,out+votes[j],0)
SeatsAfterVS2 <- GetIntegerHare(votes_aftervs2,n)
SeatChange <- (max(abs(SeatsRound1-SeatsRound2),abs(SeatsRound1-SeatsAfterVS2)))
#}
if (SeatChange!=0){ # case in which seat allocation changes after round 1 or after calculation of vs2
jchange <- c(row(as.matrix(votes))[(SeatsRound2-SeatsRound1)!=0],row(as.matrix(votes))[(SeatsRound2-SeatsAfterVS2)!=0]) # party for which seat change happens
ishifted <- c(1:n)
VS2s <-unique(matrix(mapply(GetVS2s,ishifted=ishifted, MoreArgs=list(votes=votes,j=j,i=i,n=n)), ncol = 1))
VS2s <- VS2s[VS2s>0]
#print(paste("VS2s=",VS2s))
VS2_potential <- NA
for (jmod in 1:length(VS2s)){
FreeSeats <-GetFreeSeatsRound2(votes,j,i,n,VS2s[jmod])
RemainderOthers <- GetRemainderRound2(votes,j,i,n,VS2s[jmod])[row(as.matrix(votes))!=j]
Remainderj <- GetRemainderRound2(votes,j,i,n,VS2s[jmod])[row(as.matrix(votes))==j]
jequalremainder <- which(round(RemainderOthers%%10,7)==round(Remainderj%%10,7))
#print(paste("VS2s[jmod]:",VS2s[jmod]))
RemainderRank <- GetRemainderRankRound2(votes,j,i,n,VS2s[jmod])
if (FreeSeats>=(length(votes)-RemainderRank[jequalremainder])){
out <- rbind(out,VS2s[jmod])
}
else{}
}
out <- min(out,na.rm=T)
}
}
return(out)
}

GetVS_OLD <- function(votes,j,i,n){
if (length(GetRelevantVS2(votes,j,i,n))>0 ){
out <- votes[row(as.matrix(votes))==j]-(GetVotesPlusVS1(votes,j,i,n)+GetRelevantVS2(votes,j,i,n))
}
if (length(votes)==2){
out <-  GetRelevantVS2(votes,j,i,n)  
}
if (length(GetRelevantVS2(votes,j,i,n))==0){
out <- votes[row(as.matrix(votes))==j]-GetVotesPlusVS1(votes,j,i,n)
}
#votes_required[which(votes_required<0)] <- 0
return(out)
#return(VS)
}


GetVSs <- function(votes,n){
out <- data.frame()
for(j in 1:length(votes)){
  #print(paste("Party:",j))
VS=mapply(GetVS,i=c(1:n),MoreArgs=list(votes=votes,j=j,n=n))
out_j <- data.frame(party =j,seat = c(1:n))
out_j$VS <- VS
out <- rbind(out,out_j)
}
return(out)
}

GetVSs(votes,10)


CheckGetVS <- function(votes,n){
  # This function computes the vote surplus/shortfall (VS) for a vector of votes and n seats. 
  # In a next step, the function then calculates the seats for the modified vote vector votes_mod (votes-VS) and based on that all remainders and the relevant remainder (Remainder_relevant_j). 
  # First check: The function calculates check1 which indicates whether Remainder_relevant_j has length 1 or whether its rank is zero (if free seats>0 and remainder<1). 
  # Second check: The function computes whether the party with the relevant remainder actually gets a seat (if free seats>0 and remainder<1). 
  # Third check: The function checks whether we have a situation with no free seat and party j has an remainder of 1. 
  # Return: The function only returns problematic cases. As long as an empty dataframe is returned, everything works fine. 
  VS_all <- GetVSs(votes,n)
  #print(VS_all)
  for (j in 1:length(votes)){
    for (i in 1:n){
      #print(paste("Party:",j))
      #print(paste("Seat:",i))
      votes_mod <- votes
      votes_mod[j] <- votes[j]-VS_all$VS[VS_all$seat==i & VS_all$party==j] 
      GetSeatsHare(votes_mod,n)
      SeatsRound1 <- GetIntegerHare(votes_mod,n)
      FreeSeats <- n-sum(SeatsRound1)
      RemainderAll <- GetRemainderHare(votes_mod,n)
      Remainderj <- RemainderAll[row(as.matrix(votes))==j]
      if (FreeSeats>0){
        #print("new calculation")
        if (Remainderj<1){
          Remainder_others <- RemainderAll[row(as.matrix(votes))!=j]
          Remainder_relevant_j <- Remainder_others[round(Remainder_others%%10,10)==round(Remainderj%%10,10)]
          RankRemainder_relevant_j <-length(votes)-rank(Remainder_others)[round(Remainder_others%%10,10)==round(Remainder_relevant_j%%10,10)]
          VS_all$check_lengthremainder[VS_all$seat==i & VS_all$party==j] <- length(Remainder_relevant_j) 
          VS_all$RankRemainder_relevant_j[VS_all$seat==i & VS_all$party==j] <- ifelse(length(RankRemainder_relevant_j)>0,RankRemainder_relevant_j,0)
          VS_all$FreeSeats[VS_all$seat==i & VS_all$party==j] <- FreeSeats
        }
        if (Remainderj==1){} # remainder of 1 is fine as long as there are not 0 free seats (see if condition above)
        else{}
      }
      if (FreeSeats>0 & Remainderj==1){
        VS_all$check_remainder_freeseats[VS_all$seat==i & VS_all$party==j] <- 1
      }
    }
  }
  #return(VS_all)
  return(VS_all[VS_all$check_lengthremainder!=1 | VS_all$RankRemainder_relevant_j==0 | RankRemainder_relevant_j | VS_all$RankRemainder_relevant_j>VS_all$FreeSeats | VS_all$check_remainder_freeseats!=1,])
}



CheckGetVSs <- function(nparties,nseats,niterations){
  it=1
  set.seed(it)
  out <- data.frame()
  while (dim(out)[1]==0 & it <niterations){
    votes <- round(runif(nparties,0,100000))
    #print(votes)
    out <- CheckGetVS(votes,nseats)
    it=it+1
    #print(dim(out)[1])
  }
  if(it<niterations){return(out)}
}



CheckGetVSsAll<- function(max_nparties,max_nseats,max_niterations){
  for (j in 2:max_nparties){
    print(paste("Number of parties:",j))
    for (i in 2:max_nseats){
      print(paste("Number of seats:",i))
      for (it in 1:max_niterations){
        if (it%%10==0){
          print(paste("Number of iteration:",it))
        }
        CheckGetVSs(j,i,it)
      }
    }
  }
}



GetRVSimulation <- function(votes,n,convcrit){
  out <- data.frame()
  for (j in 1:length(votes)){
    #print(paste("Party: ",j))  
    for (i in 1:n){
      #print(paste("Seat: ",i))  
      votessim <- votes
      votessim[j] <- 0.5 # start with 1 vote in first iteration (first line after while loop)
      #print(votessim)
      seats <- 0
      while (seats<i){
        votessim[j] <- votessim[j]*2
        seats <- GetSeatsHare(votessim,n)[j]
        #print(votessim)
      }
      xlow=votessim[j]/2 
      xhigh=votessim[j]
      #print(paste("Seat change at: ",votessim[j]))  
      
      while ((xhigh-xlow)>convcrit){
        votessim[j]<- (xlow+xhigh)/2
        seats <- GetSeatsHare(votessim,n)[j]
        if (seats<i){
          xlow=votessim[j]  
        }
        if (seats>=i){
          xhigh=votessim[j]
        }
        #print(paste("xhigh: ",xhigh))  
        #print(paste("xlow: ",xlow))  
        
      }
      out_ji <- data.frame(party=j,seat=i,VS=round(votes[j]-votessim[j],4))
      out <- rbind(out,out_ji)
    }
  }
  return(out)
}


library(tidyverse)
votes <- c(0,50000,40000)
nseats=10
i=5;j=5;
votes <- c(47000, 16000, 15800, 12000, 6100, 3100)
OutSimulation <- GetRVSimulation(votes,nseats,0.000001)
OutSimulation <- OutSimulation %>% rename(VS_Sim=VS)
OutAnalytical <- GetVSs(votes,nseats)
OutAnalytical <- OutAnalytical %>% rename(VS_Ana=VS)

(OutCompare <- OutAnalytical %>%  left_join(OutSimulation,by=c("party","seat")))
GetRemainderRankRound2s(votes,5,5,10)
GetRemainderRound2s(votes,5,5,10)

diag(GetRemainderRankRound2s(votes,5,5,10))
GetFreeSeatsRound2s(votes,5,5,10)
# GetRelevantRemainderRank(partyvotes,1,5,10)
# GetRelevantVS2(partyvotes,1,5,10)
# GetVotesPlusVS1(partyvotes,1,5,10)






CompareAnalayticalSimulation <- function(nparties,nseats,niterations,seedstart=1){
  it=seedstart
  out <- data.frame()
  check=0
  while (abs(check)<=0.1 & it<(niterations+seedstart)){
    set.seed(it)
    if(it%%100==0){print(paste("Iteration:",it))}
    votes <- round(runif(nparties,0,100000))
    OutSimulation <- GetRVSimulation(votes,nseats,0.000001)
    OutAnalytical <- GetVSs(votes,nseats)
    check<-max(OutSimulation$VS-OutAnalytical$VS)
    it=it+1
  }
  if(it<(niterations+seedstart)){return(as.vector(votes))}
}


# (iii) Test of seat allocation
# 
GetIntegerHare(partyvotes,10)
(wiki_remainders <- GetRemainderHare(partyvotes,10))
GetHighestRemainder(wiki_remainders,3)
GetSeatsHare(partyvotes,10)
partyvotes_mod <- partyvotes
GetSeatsHare(partyvotes,10)

# (iv) Test of calcuation of vote shortfall


# V1 <- GetVS(partyvotes,1,5,10)
# partyvotes_mod=c(V1,partyvotes[2:6])
# GetNumberFreeSeats(partyvotes,1,5,10)
# VS2minusj <- GetVS2minusj(partyvotes_mod,1,5,10,5)
# GetVS2s(partyvotes,1,5,10)
# GetRatioRound2(partyvotes,1,5,10,VS2minusj)
# GetIntegerRound2(partyvotes,1,5,10,VS2minusj)
# GetRemainderRound2(partyvotes,1,5,10,VS2minusj)
# GetFreeSeatsRound2(partyvotes,1,5,10,VS2minusj)
# GetRemainderRankRound2(partyvotes,1,5,10,VS2minusj) 
# GetRemainderRound2s(partyvotes,1,5,10) # check this matrix with Simon: it is easier for party 1 to catch up with party 5 in the second round 
# GetRemainderRankRound2s(partyvotes,1,5,10)
# diag(GetRemainderRankRound2s(partyvotes,1,5,10))
# GetFreeSeatsRound2s(partyvotes,1,5,10)
# GetRelevantRemainderRank(partyvotes,1,5,10)
# GetRelevantVS2(partyvotes,1,5,10)
# GetVotesPlusVS1(partyvotes,1,5,10)



# (v) check GetVS for different seats

# (a) Wikipedia example 
# 
# (partyvotes <- c(47000, 16000, 15800, 12000, 6100, 3100)) # Wikipedia example
# GetSeatsHare(partyvotes,10)
# 
# # (a1) focus on party 1
# 
# 
# GetVS(partyvotes,1,5,10)
# partyvotes_mod <- c(partyvotes[1]-GetVS(partyvotes,1,5,10)+1,partyvotes[2:6])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# GetVS(partyvotes,1,4,10)
# partyvotes_mod <- c(partyvotes[1]-GetVS(partyvotes,1,4,10)+1,partyvotes[2:6])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# 
# GetVS(partyvotes,1,3,10)
# partyvotes_mod <- c(partyvotes[1]-GetVS(partyvotes,1,3,10)+1,partyvotes[2:6])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# 
# GetVS(partyvotes,1,2,10)
# partyvotes_mod <- c(partyvotes[1]-GetVS(partyvotes,1,2,10)+1,partyvotes[2:6])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# GetVS(partyvotes,1,1,10)
# partyvotes_mod <- c(partyvotes[1]-GetVS(partyvotes,1,1,10)+1,partyvotes[2:6])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# # (a2) focus on party 5
# 
# GetVS(partyvotes,5,5,10)
# partyvotes_mod <- c(partyvotes[1:4],partyvotes[5]-GetVS(partyvotes,5,5,10),partyvotes[6])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# GetRemainderRound2s(partyvotes,5,5,10)
# GetRemainderRankRound2s(partyvotes,5,5,10)
# diag(GetRemainderRankRound2s(partyvotes,5,5,10))
# GetFreeSeatsRound2s(partyvotes,5,5,10)
# GetRelevantVS2(partyvotes,5,5,10)
# GetVS2s(partyvotes,5,5,10)
# 
# GetVS(partyvotes,5,8,10)
# partyvotes_mod <- c(partyvotes[1:4],partyvotes[5]-GetVS(partyvotes,5,8,10),partyvotes[6])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# GetVS(partyvotes,5,10,10)
# partyvotes_mod <- c(partyvotes[1:4],partyvotes[5]-GetVS(partyvotes,5,10,10),partyvotes[6])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# # (b) Example with four parties and different votes
# 
# (partyvotes <- c(40000, 16000, 45800, 12000 ))
# GetSeatsHare(partyvotes,10)
# 
# 
# # (b1) focus on party 1
# 
# GetVS(partyvotes,1,5,10)
# partyvotes_mod <- c(partyvotes[1]-GetVS(partyvotes,1,5,10)+1,partyvotes[2:4])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# GetVS(partyvotes,1,4,10)
# partyvotes_mod <- c(partyvotes[1]-GetVS(partyvotes,1,4,10)+1,partyvotes[2:4])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# 
# GetVS(partyvotes,1,3,10)
# partyvotes_mod <- c(partyvotes[1]-GetVS(partyvotes,1,3,10)+1,partyvotes[2:4])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# 
# GetVS(partyvotes,1,2,10)
# partyvotes_mod <- c(partyvotes[1]-GetVS(partyvotes,1,2,10)+1,partyvotes[2:4])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# GetVS(partyvotes,1,1,10)
# partyvotes_mod <- c(partyvotes[1]-GetVS(partyvotes,1,1,10)+1,partyvotes[2:4])
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 


# (b2) focus on party 4
# 
# GetVS(partyvotes,4,10,10)
# partyvotes_mod <- c(partyvotes[1:3],partyvotes[4]-GetVS(partyvotes,4,10,10))
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# GetVS(partyvotes,4,5,10)
# partyvotes_mod <- c(partyvotes[1:3],partyvotes[4]-GetVS(partyvotes,4,5,10))
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# GetVS(partyvotes,4,2,10)
# partyvotes_mod <- c(partyvotes[1:3],partyvotes[4]-GetVS(partyvotes,4,2,10))
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# 
# GetVS(partyvotes,4,1,10)
# partyvotes_mod <- c(partyvotes[1:3],partyvotes[4]-GetVS(partyvotes,4,2,10))
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)


# (vi) check GetVSs  

#GetVSs(partyvotes,10)



CheckGetVS(c(26552,37212),10)
GetVSs(c(26551,37212),10)
GetVS(c(26551,37212),1,1,10)
CheckGetVSs(5,10,50) 
CheckGetVSsAll(10,10,100)  





  

# (vii) check RV using simulation approach  



# (Out <- GetRVSimulation(c(9,20),2,0.000001))
# GetVSs(c(9,20),2)
# 
# GetSeatsHare(c(9,3),2)
# GetSeatsHare(c(60.5,20),2)
# 
# GetIntegerHare(c(60,20),2)
# GetRemainderHare(c(60,20),2)
# GetIntegerHare(c(9,27),2)
# GetRemainderHare(c(9,27),2)
# 
# 
# (Out <- GetRVSimulation(c(9,20,11),2,0.000001))
# GetVSs(c(9,20,11),2)




(Out2_partyvotes <- CompareAnalayticalSimulation(nparties=4,nseats=10,niterations=1000,seedstart=1200))


for (np in 4:20){
print(paste("Number of parties: ",np))
for (ns in 1:20){
print(paste("Number of seats: ",ns))
(Out_Check <- CompareAnalayticalSimulation(nparties=np,nseats=ns,niterations=1000,seedstart=1))
}
}

#partyvotes <- c(26551,37212,57285)
#partyvotes <- c(54998,55267,23889)
#partyvotes <- c(65769,24986,30005)
#partyvotes <- c(12739, 76672,5442)
partyvotes <- c(879,22891,52524)

(OutSim <- GetRVSimulation(Out2_partyvotes,10,0.000001))
(OutAna <- GetVSs(Out2_partyvotes,10))

as.numeric(Out1$VS-Out2$VS)


library(tidyverse)
OutSim %>% rename(VSSim=VS) %>% left_join(OutAna,by=c("party","seat"))  %>% rename(VSAna=VS)%>% 
           select(party,seat,VSAna,VSSim) %>%
           mutate(VSdiff=round(VSAna-VSSim,2))


votes <- c(26551,37212,57285); j=2; i=1;n=3
partyvotes <- votes
V1 <- GetVS(partyvotes,j,i,n)
GetRatioRound2(partyvotes,j,i,n,836.6) # Note: here party 3 gets two seats in the second round and party 2's remainder is the same as party 3's remainder with two (not one!) seats
GetIntegerRound2(partyvotes,j,i,n,836.6)
GetRemainderRound2(partyvotes,j,i,n,836.6)

for (j in 1:3){
for (i in 1:10){
GetRelevantVS2(c(73103,61960,56252,22619),j,i,10)
}
}


# # Old stuff
# 
# # debuging GetVS(5,5,10) which returns nun-numeric value
# # reason: no solution for which rank(remainder_minusjtilde)==number of free seats, however in one solution, there are more seats available than the relevant party minusjtilde achieves and thus this 
# # challenge: what happens if for all VS2_minusjtilde there are less seats than the relevant party minusjtilde achieves?
# 
# V1 <- (partyvotes,5,5,10)
# partyvotes_mod <- partyvotes
# partyvotes_mod[5] <- V1
# GetSeatsHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)
# 
# 
# GetNumberFreeSeats(partyvotes,5,5,10)
# GetVS2s(partyvotes,5,5,10)
# GetRemainderRound2s(partyvotes,5,5,10)
# GetRemainderRankRound2s(partyvotes,5,5,10)
# GetFreeSeatsRound2s(partyvotes,5,5,10)
# 5-diag(GetRemainderRankRound2s(partyvotes,5,5,10))
# VS2s_all <- GetVS2s(partyvotes,5,5,10)
# 
# for (j in c(1:length(VS2s_all))){
# print(GetFreeSeatsRound2(partyvotes,5,5,10,VS2s_all[j]))
# print(GetRatioRound2(partyvotes,5,5,10,VS2s_all[j]))
# #print(GetRemainderRankRound2s(partyvotes,5,5,10))  
# }
# 
# GetRelevantVS2(partyvotes,5,5,10)
# 
# VS <- GetVS(partyvotes,5,5,10)
# partyvotes_mod <- partyvotes
# partyvotes_mod[5] <- partyvotes[5]-VS
# GetSeatsHare(partyvotes_mod,10)
# GetIntegerHare(partyvotes_mod,10)
# GetRemainderHare(partyvotes_mod,10)

