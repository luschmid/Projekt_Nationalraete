rm(list=ls())

data_example <- data.frame(district=rep(1,9),
                           year=rep(2024,9),
                           id=c(1:9),
                           party=c(rep("P1",3),rep("P2",3),rep("P3",3)),
                           elected=c(1,1,0,1,0,0,0,0,0),
                           pvotes=c(rep(45,3),rep(35,3),rep(20,3)),
                           cvotes=c(c(25,11,9),c(22,8,5),c(8,7,5))
)
print(data_example) 

data_example_prep <- PrepareData(data=data_example, 
                                 district_name="district", 
                                 election_cycle_name = "year",
                                 system="open",
                                 cand_id_name="id",
                                 party_name="party", 
                                 votes_j_name="pvotes",
                                 votes_h_name="cvotes",
                                 alliances=F) 

print(data_example_prep) 


CalculateMargins(data_input=data_example_prep ,
                 method="dHondt",
                 convcrit=0.001,
                 system="open",
                 return_option = TRUE)



