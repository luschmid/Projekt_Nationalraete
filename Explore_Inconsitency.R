source("C:/Schmidlu/Dropbox/Projekt Nationalräte/03_Code/13_Running_variable/00_Vote_margin_analytical_highestaverage.R")


# Functions

GetSeatsHighestAverageAll <- function(data,n,method="dHondt",threshold="none"){
seat_out <- rep(NA,dim(data)[1])
for (j in c(1:dim(data)[1])){
seat_out[j] <- GetSeatsHighestAverage(party_input=j, data_input=data,n=n,
                                      method=method,threshold=threshold)
}
return(seat_out)
}

DataPrepHelp <- function(votes){
  
  df <- data.frame(votes,party=c(1:length(votes)),
                   district=rep(1,length(votes)),
                   year=rep(1,length(votes)))
  
  return(PrepareData(data=df,votes_j_name="votes",
                     party_name="party", 
                     districtname="district", 
                     election_cyclename = "year",
                     alliances=F,
                     system="closed"))
}

SeatDistribution_Indirect <- function(votes,n=nseats,correction_unit=0.00001){
# 1. Start 
votes_j <- votes[1]
votes_minusj <- votes[2:length(votes)]
ci_all <-  rep(NA,n)
party_seat <-  rep(NA,n)
# 2. Search c for each i
for (i in 1:n){
df_prep <- DataPrepHelp(c(votes_j,votes_minusj))
seats <- GetSeatsHighestAverageAll(df_prep,n=n)
seats_others <- seats[2:length(votes)]
partys_underconsid <- which(seats_others>=1)
if (length(partys_underconsid)>0){
ci <- rep(NA,length(partys_underconsid))
for (minusj in partys_underconsid){
ci[minusj] <- (votes_minusj[minusj]*i-votes_j*seats[minusj+1])/(i+seats[minusj+1]) # Note: seats[minusj+1] is itilde
}
ci_binding <- min(ci)
minusj_ci_binding <- which(ci==ci_binding)
# 3. Update votes for next iteration
votes_j <- votes_j + ci_binding + correction_unit
votes_minusj[which(ci==ci_binding)] <- votes_minusj[which(ci==ci_binding)] -ci_binding -correction_unit
if (i>1) {ci_all[i] <- ci_all[(i-1)]+ci_binding
} else{
ci_all[i] <- ci_binding  
party_seat[i] <- which(ci==ci_binding) + 1
}
}else{
ci_all[i] <- NA 
party_seat[i] <- NA
}

}
return(list(ci_all,party_seat))
}

SeatDistribution_Direct <- function(votes,n=nseats,correction_unit=0.00001){
  # 1. Start 
  votes_j <- votes[1]
  votes_minusj <- votes[2:length(votes)]
  ci_all <-  rep(NA,n)
  party_seat <-  rep(NA,n)
  df_prep <- DataPrepHelp(c(votes_j,votes_minusj))
  seats <- GetSeatsHighestAverageAll(df_prep,n=n)
  seats_others <- seats[2:length(votes)]
  
  
  # 2. Search c for each i
  for (i in 1:n){
    partys_underconsid <- which(seats_others>=i)
    
    if (length(partys_underconsid)>0){
    ci <- rep(NA,length(partys_underconsid))
    
    for (minusj in partys_underconsid){
      #print(partys_underconsid)
      ci[minusj] <- (votes_minusj[minusj]*i-votes_j*(seats_others[minusj]-i+1))/(i+(seats_others[minusj]-i+1)) # Note: seats[minusj+1] is itilde
    }
    ci_binding <- min(ci,na.rm=T)
    minusj_ci_binding <- which(ci==ci_binding)
    ci_all[i] <- ci_binding
    party_seat[i] <- which(ci==ci_binding)
    } else{
    ci_all[i] <- NA
    }
    
  }
  return(list(ci_all,party_seat))
}

RandomizeElections <- function(nparties=3,nseats=3,runif_min=10,runif_max=20){

set.seed(123)
min_diff <- -5
while (min_diff<=0){
# 1. Randomize elections with first party as the one the as minimal no. votes
votes <- cumsum(round(runif(nparties,min=runif_min,max=runif_max)))

# 2. Calculate reshuffling using direct and indirect rule
out_indirect <- SeatDistribution_Indirect(votes,n=nseats)[[1]]
out_direct <- SeatDistribution_Direct(votes,n=nseats)[[1]]

# 3. Calculate difference between indirect and direct and take minimum
diff <- out_indirect-out_direct
min_diff <- min(diff,na.rm=T)
}
return(list(votes,out_direct,out_indirect))

}

# B) Look at one example

# votes <- c(20,45,35)
# 
# SeatDistribution_Indirect(votes,n=3)
# SeatDistribution_Direct(votes,n=3)

