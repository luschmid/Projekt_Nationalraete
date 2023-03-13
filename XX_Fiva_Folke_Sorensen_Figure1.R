

df_example <- data.frame(votes=c(35,45,20),party=c(1,2,3),
                         district=rep(1,3),year=rep(2021,3))

df_example_prep <- PrepareData(data=df_example, 
                                votes_j_name="votes",
                                party_name="party", 
                                districtname="district", 
                                election_cyclename = "year",
                                alliances=F,
                                system="closed") 



SimulateVotesAdd <- function(party_input,
                             data_input,
                             votes_min=0,
                             votes_max=1000,
                             votes_increase=1,
                             no_seats,
                             method="dHondt"){

df_out <- NULL

for (i in seq(votes_min,votes_max,votes_increase)){
#print(i)
data_input$votes_j[data_input$party==party_input] <- i
df_out<- rbind(df_out,data.frame(votes=i,
                             seats=GetSeatsHighestAverage(party_input=party_input,
                                                         data_input=data_input, 
                                                         no_seats=no_seats,
                                                         method=method)))
}
return(df_out)
}

result <- SimulateVotesAdd(party_input=1,
                          data_input=df_example_prep,
                          votes_min=0,
                          votes_max=1000,
                          votes_increase=1,
                          no_seats=3,
                          method="dHondt")

result <- result %>% mutate(votes_total=65+votes,
                            vote_share=votes/votes_total)

result %>% head(n=30)

ggplot(data=result,aes(x=vote_share,y=seats)) +
  geom_point() +
  theme_bw()

result$seat_diff<- c(0,result$seats[2:(length(result$seats))]-
                       result$seats[1:(length(result$seats)-1)])             

result %>% filter(seat_diff>0)
