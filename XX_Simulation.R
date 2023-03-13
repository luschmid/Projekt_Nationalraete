### RANDOM ELECTIONS ################## 

nseats_district=3;
eligibles_district=1000;
turnout_district=0.5;
voter_preferences_district=c(0.45,0.35,0.2);
individual_votes_distribution=c(0.5,0.3,0.2);



RandomElectionResult <- function(nseats_district,
                                 eligibles_district,
                                 turnout_district,
                                 voter_preferences_district,
                                 individual_votes_distribution){
  # 1. Calculate voters and votes per district 
  voters_district <- eligibles_district*turnout_district
  votes_district <- voters_district*nseats_district  
  out <- data.frame()
  # 2. Calculate voter preferences (with random component) and party votes 
  #    (with random component) 
  voter_preferences_district <- voter_preferences_district+rnorm(length(voter_preferences_district),0,0.01)
  votes_parties=as.vector(table(sample(c(1:length(voter_preferences_district)), size=votes_district, replace=TRUE, 
                                       prob=voter_preferences_district)))
  #print(votes_parties)
  # 3. Calclate individual votes
  for (i in 1:length(voter_preferences_district)){
    #print(i)
    votes_h=as.vector(table(sample(c(1:nseats_district), size=votes_parties[i], replace=TRUE, prob=individual_votes_distribution))) 
    if (length(votes_h)<nseats_district){votes_h <- c(votes_h,rep(0,nseats_district-length(votes_h)))}# fill up candidates with zero votes with a 0
    out <- rbind(out,data.frame(party=rep(i,nseats_district),
                                candidate=c(1:nseats_district),
                                votes_h=votes_h)
    )
  }
  return(out)
}

RandomElectionResult(
  nseats_district=3,
  eligibles_district=1000,
  turnout_district=0.5,
  voter_preferences_district=c(0.45,0.35,0.2),
  individual_votes_distribution=c(0.5,0.3,0.2)
)



RandomElectionResult(nseats_district,
                     eligibles_district,
                     turnout_district,
                     voter_preferences_district,
                     individual_votes_distribution)


nseats_districts=5;eligibles_district=100;turnout_district=0.2;voter_preferences_district=c(0.1,0.3,0.6);
individual_votes_distribution=rep(1/5,5);ndistricts=3;

SimulateElection<- function(ndistricts,
                            nparties, 
                            nseats, 
                            eligibles, 
                            turnout,voter_preferences,
                            seedstart = 1){
  #################################################################################################################
  # This function simulates the results of one election for the investigation of the imbalance problem regarding small/big districts
  # (e.g. Francisco Morazan in Honduras). 
  # Arguments:    ndistricts: number of districts (scalar)
  #               nparties: number of parties within district (list)
  #               nseats:  number of seats in district (list)
  #               eligibles: eligible voters in district (list)
  #               turnout: turnout in district(list)
  #               voter_preferences: preferences for voters for parties (list [ndistricts] of lists [number of parties])
  #               seedstart: start of seed
  #################################################################################################################
  ElectionResultAll <- data.frame()
  set.seed(seedstart)
  for (j in c(1:ndistricts)){
    #print(j)
    ElectionResult <- RandomElectionResult(nseats_district=nseats[j],
                                           eligibles_district=eligibles[j],
                                           turnout_district=turnout[j],
                                           voter_preferences_district=voter_preferences[[j]],
                                           individual_votes_distribution=rep(1/nseats[j],nseats[j]))
    ElectionResult$district <- j
    ElectionResult$year <- 2020
    ElectionResultAll <- rbind(ElectionResultAll,ElectionResult)
  }
  
  df_input_seats <- data.frame(no_seats=nseats,year=2020, district=c(1:ndistricts),eligibles=eligibles)
  
  Margins <-CalculateMargins(data_input=ElectionResultAll , data_input_seats=df_input_seats) %>%
    left_join(df_input_seats,by=c("year","district")) %>%
    mutate(votemargin_rel=votemargin/eligibles)
  Margins <-   MakeDummies(Margins,"district",prefix="district_") 
  
  
  out <- GenerateRDDTable(datainput=Margins,
                          outcomevariables=paste("district_",1:length(levels(as.factor(Margins$district))),sep=""),
                          runningvariable="votemargin_rel",
                          districtvariable = "district",
                          fuzzyvariable=0,
                          pvector=c(2, 3, 2, 3, 2, 3),
                          hvector= c(1, 1, 0.05, 0.05, 0.01, 0.01),
                          latex=0)
  return(out)
}


SimulateElections<- function(niterations=100,ndistricts,nparties, nseats, eligibles, turnout,voter_preferences,
                             seedstart = 1){
  # This function simulates the results of multiple elections for the investigation of the imbalance problem regarding small/big districts
  # (e.g. Francisco Morazan in Honduras). 
  out <- list()
  i=1
  while(i < niterations){
    print(paste("Iteration number:",i,sep=""))
    out <- list(out,SimulateElection(nparties=nparties, nseats=nseats, ndistricts=ndistricts,
                                     eligibles=eligibles, 
                                     turnout=turnout,voter_preferences=voter_preferences,
                                     seedstart = i)  )
    i=i+1
  }
  return(out)
}

nparties_case <- c(3,7,4)

SimulateElections(niterations=100,nparties=nparties_case, nseats=nparties_case, ndistricts=length(nparties_case),
                  eligibles=c(10,20,5)*10000, 
                  turnout=rep(0.5,ndistricts),voter_preferences=list(c(rep(1/nparties_case[1],nparties_case[1])),
                                                                     c(rep(1/nparties_case[2],nparties_case[2])),
                                                                     c(rep(1/nparties_case[3],nparties_case[3]))), # uniform party preference distribution,
                  seedstart = 1)


# ndistricts=3;
# nparties=c(3,7,2); nseats=round(nparties);eligibles=c(10,20,5)*10000; turnout <- rep(0.5,ndistricts)
# voter_preferences <-c(0.1,0.5,0.4)
# 
### END OF RANDOM ELECTIONS ################## 