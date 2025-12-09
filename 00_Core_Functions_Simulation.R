
CalculateMargins <- function(data_input,
                             system="open",
                             method="dHondt",
                             convcrit=0.001,
                             return_option=FALSE,
                             outfile_name=NULL,
                             threshold="none",
                             calculate_all_seats=TRUE,
                             total_seats_name="none",
                             append=FALSE,
                             additional_vars=c("alliance"),
                             print_party=FALSE) {
  
  
  # 0. Check arguments
  
  if (return_option==FALSE & is.null(outfile_name)) {
    stop("Outfile name is missing. Please specify outfile_name")
  }
  
  # 1. Check for single candidates
  
  data_singlecandidate <- data_input %>%
    group_by(year, district) %>%
    summarize(nobs = n()) %>%
    filter(nobs == 1)
  
  # 2. Loop over years and districts
  
  years <- as.numeric(levels(as.factor(data_input$year)))
  Out <- data.frame()
  for (yr in years) {
    print(paste("New year:", yr))
    dataperyear <- data_input[data_input$year == yr, ]
    districts <- levels(as.factor(dataperyear$district))
    for (di in districts) {
      print(paste("New district:", di))
      dataworking <- dataperyear[dataperyear$district == di, ]
      
      
      # 3. Define number of seats as sum of elected candidates if not provided
      #    by variable total_seats_name (string)
      
      if (total_seats_name=="none"){
        n <- sum(dataworking$elected)
      } else{
        n <- dataworking[[total_seats_name]][1]  
      }
      
      
      # 4. Calculate simulated vote margins 
      
      out_temp <- as.data.frame(GetRVSimulation_Candidate(data_input=dataworking, 
                                                          system=system,
                                                          n=n,
                                                          calculate_all_seats=calculate_all_seats,
                                                          method=method,
                                                          convcrit=convcrit,
                                                          threshold=threshold,
                                                          print_party=print_party))
      
      
      if (return_option==TRUE){
        Out <- bind_rows(Out,out_temp)
      } else{
        if (append==FALSE){
          data.table::fwrite(out_temp, file = outfile_name, 
                             sep = ";",row.names = F,append=F)  
        } else{
          data.table::fwrite(out_temp, file = outfile_name, 
                             sep = ";",row.names = F,append=T)  
        }
      }
    }
  }
  if (return_option==TRUE){
    return(Out)
  }
  
}

CalculateRatios <- function(votes, n,method="dHondt") {
  votes <- as.matrix(votes)
  
  # 1. Calculate ratios for different highest average methods
  
  if (method=="dHondt"){
    div <- matrix(rep(c(1:(n)), length(votes[, 2])), 
                  nrow = length(votes[, 2]), byrow = T)
  }
  if (method=="SainteLague"){
    div <- matrix(rep(c(1.4,seq(3,n*2,2)), 
                      length(votes[,2])), 
                  nrow = length(votes[,2]), 
                  byrow = T)
  }
  
  # 2. Write in long dataframe and change column names to party names
  
  Votes_wide <- as.data.frame(cbind(votes[, 1], 
                                    matrix(rep(as.numeric(votes[, 2]), 
                                               n), 
                                           nrow = length(votes[, 2])) / div))
  colnames(Votes_wide) <- c("party", paste(1:n, sep = ""))
 
   # 3. Reshape to wide and return
  
  Votes_long <- gather(Votes_wide, paste(1:n, sep = ""), 
                       key = "seat", 
                       value = "No")
  return(Votes_long)
}



CollapseVotes <- function(data_input, aggregation_level, aggregation_unit = NA){ 

  # 1. Define aggregation level and aggregation units
  
  aggregation_level <- enquo(aggregation_level)
  aggregation_unit <- enquo(aggregation_unit)
  check_aggregation_level <- names(data_input %>% ungroup() %>% 
                                     select(!!aggregation_level))
  
  # 2. Select only relevant alliance or suballiance
  
  if (check_aggregation_level == "suballiance") {
    data_input <- data_input %>% filter(alliance == !!aggregation_unit)
  }
  if (check_aggregation_level == "party") {
    data_input <- data_input %>% filter(suballiance == !!aggregation_unit)
  }

  # 3. Collapse by aggregation_level
  
    data_out <- data_input %>%
    group_by(party) %>%
    filter(row_number() == 1) %>%
    group_by(!!aggregation_level) %>%
    summarize(votes_j = sum(votes_j)) %>%
    as.matrix()
  
  return(data_out)
}

CutDataByQuorum <- function(data_input, threshold){
  
  # 1. Spain
  
  if (threshold=="spain_2004_2023"){
    
    blankv <- data_input$blank_votes[1]
    p_votes_total <- data_input %>%
      filter(seat==1) %>%
      select(votes_j) %>%
      summarize(p_votes_total=sum(votes_j)) %>%
      pull()
    
    data_input_out <- data_input %>%
      mutate(votes_j_rel=votes_j/(p_votes_total+blankv))%>% 
      filter(votes_j_rel>=0.03) %>%
      select(-votes_j_rel)
    
    
  }
  
  # 2. Israel
  
  if (threshold=="israel_2009_2013"){
    
    p_votes_total <- data_input %>%
      filter(seat==1) %>%
      select(votes_j) %>%
      summarize(p_votes_total=sum(votes_j)) %>%
      pull()
    
    data_input_out <- data_input %>%
      mutate(votes_j_rel=votes_j/p_votes_total)%>% 
      filter(votes_j_rel>=0.02) %>%
      mutate(total_valid_votes=sum(votes_j))%>%
      select(-votes_j_rel)
    
    
  }
  
  if (threshold=="israel_2014_2022"){
    
    
    p_votes_total <- data_input %>%
      filter(seat==1) %>%
      select(votes_j) %>%
      summarize(p_votes_total=sum(votes_j)) %>%
      pull()
    
    data_input_out <- data_input %>%
      mutate(votes_j_rel=votes_j/p_votes_total)%>% 
      filter(votes_j_rel>=0.0325) %>%
      mutate(total_valid_votes=sum(votes_j))%>%
      select(-votes_j_rel)
    
  }
  
  return(data_input_out)
  
}

GetAllianceSuballiance <- function(party_input, data_input) {
  party_input <- enquo(party_input)
  
  # 1. Generate alliance and suballiance dummies
  
  alliance_dummy <- data_input %>%
    filter(party == !!party_input) %>%
    ungroup() %>%
    select(alliance_dummy) %>%
    filter(row_number() == 1) %>%
    pull()
  suballiance_dummy <- data_input %>%
    filter(party == !!party_input) %>%
    ungroup() %>%
    select(suballiance_dummy) %>%
    filter(row_number() == 1) %>%
    pull()
  
  # 2. Get relevant alliance and suballiance
  
  alliance <- data_input %>%
    filter(party == !!party_input) %>%
    ungroup() %>%
    select(alliance) %>%
    filter(row_number() == 1) %>%
    pull()
  suballiance <- data_input %>%
    filter(party == !!party_input) %>%
    ungroup() %>%
    select(suballiance) %>%
    filter(row_number() == 1) %>%
    pull()
  return(data.frame(alliance_dummy = alliance_dummy, 
                    suballiance_dummy = suballiance_dummy, 
                    alliance = alliance, suballiance = suballiance))
}

GetRVSimulation_Candidate <- function(data_input, 
                                      n, 
                                      convcrit=0.001,
                                      method="dHondt",
                                      system="open",
                                      threshold="none",
                                      calculate_all_seats=TRUE,
                                      additional_vars="",
                                      print_party=FALSE) {
  
  # data_input: data after PrepareData with votes_j as party votes and votes_h 
  #             as candidate votes
  # n:          number of seats
  # convcrit:   convergence criterion
  # method:     method to distribute seats
  # system:     close-list or open-list
  
  # 0. Check arguments
  
  if (is.null(n)) {
    stop("Number of seats is missing. Please specify argument n in the function")
  }
  
  out <- data.frame()
  
  # 1. Code for closed-list systems
  
  if (system=="closed"){
    parties_all <- as.numeric(levels(as.factor(data_input$party)))
    for (j in parties_all) {
      if (calculate_all_seats==TRUE){no_seat <- n
      } else{
        no_elected <- dim(data_input %>% filter(party==j & elected==1))[1] 
        no_seat <- max(no_elected,1)
      }
      for (i in 1:no_seat) {
        
        votessim <- GetVotesRequiredClosedList(data_input=data_input,
                                                i=i,
                                                j=j,
                                                n=n,
                                                convcrit=convcrit,
                                                method=method,
                                                threshold=threshold)
        
        out_ji <- data.frame(party = j, seat = i, 
                             votemargin = round(data_input$votes_j[data_input$party==j][1]
                                                - votessim, 4))
        out <- rbind(out, out_ji)
      }
    }
    out <- data_input %>% left_join(out, by=c("party","seat"))
    
  }
  
  # 2. Code for open-list systems
  
  if (system=="open"){
    parties_all <- as.numeric(levels(as.factor(data_input$party)))
    
    for (j in parties_all) {
      ids_all <- levels(as.factor(data_input$id_cand[data_input$party==j]))
      if (print_party==TRUE){print(paste0("Party:", j))}
      for (i in ids_all) {
        
        votessim <- GetVotesRequiredOpenList(data_input=data_input,
                                              i=i,
                                              j=j,
                                              n=n,
                                              method=method,
                                              threshold=threshold,
                                              convcrit=convcrit)
        
        
        out_ji <- data.frame(party = j, id_cand = data_input$id_cand[data_input$id_cand==i], year=data_input$year[1],
                             district=data_input$district[1],
                             votemargin = round(data_input$votes_h[data_input$id_cand==i & data_input$party==j]
                                                - votessim, 4))
        out <- rbind(out, out_ji)
      }
    }
    out <- data_input %>% 
      select(-district,-year) %>%
      left_join(out, by=c("party","id_cand"))
  }
  
  return(out)
}


GetSeatsHighestAverage <- function(party_input, 
                                   data_input, 
                                   n,
                                   method="dHondt",
                                   threshold="none") {
  # This function calculates the number of seats for a party using highest 
  # average methods
  # Input: party_input: choice of party
  #        data_input:  votes data with columns party, alliance, suballiance, 
  #                     votes
  #        n:           number of seats in a district
  #        method:      allocation method (dHondt,SainteLague)
  #        threshold:   threshold for parties to be considered
  
  # 0. Check whether each party makes threshold and include only parties that
  #    fulfill threshold
  
  
  quorum=TRUE # set quorum (binary indicator whether threshold is passed) to true as default: 
  
  if (threshold!="none"){
    data_input <- CutDataByQuorum(data_input,threshold)
    quorum <- is.element(party_input,levels(as.factor(data_input$party)))
    
    if (quorum==FALSE){ out <- 0  # calculation only for parties who made quorum
    }
  }
  
  if (quorum!=FALSE | threshold=="none"){
    
    AllianceInformation <- GetAllianceSuballiance(party_input, data_input) 
    
    # 1. Calculation of number of seats across alliances
    
    VotesAlliance <- CollapseVotes(data_input, aggregation_level = alliance)
    AllianceSeats <- CalculateRatios(votes=VotesAlliance, 
                                     n=n,
                                     method=method) %>% 
      mutate(No_rank = row_number(-as.numeric(No))) %>%
      filter(No_rank <= n) %>%
      group_by(party) %>%
      summarize(seats_alliance=n())%>%
      rename(alliance=party)
    
    out <- max(AllianceSeats$seats_alliance[AllianceSeats$alliance==
                                              AllianceInformation$alliance],0)
    
    # 2. Calculation of suballiance seats
    
    if (AllianceInformation$alliance_dummy == 1 & out>0) {
      VotesSuballiance <- CollapseVotes(data_input, 
                                        aggregation_level = suballiance, 
                                        aggregation_unit = AllianceInformation$alliance)
      SuballianceSeats <- CalculateRatios(votes=VotesSuballiance, 
                                          n= AllianceSeats$seats_alliance[AllianceSeats$alliance==AllianceInformation$alliance],
                                          method=method) %>% 
        mutate(No_rank = row_number(-as.numeric(No))) %>%
        filter(No_rank <= AllianceSeats$seats_alliance[AllianceSeats$alliance==AllianceInformation$alliance]) %>%
        group_by(party) %>%
        summarize(seats_suballiance=n())%>%
        rename(suballiance=party)
      
      out <- max(SuballianceSeats$seats_suballiance[SuballianceSeats$suballiance==
                                                      AllianceInformation$suballiance],0)
      
      # 3. Calculation of margin across parties (only for parties in suballiance)
      if (AllianceInformation$alliance_dummy == 1 & 
          AllianceInformation$suballiance_dummy == 1 & out>0) {
        VotesParty <- CollapseVotes(data_input, aggregation_level = party, 
                                    aggregation_unit = AllianceInformation$suballiance)
        PartySeats <- CalculateRatios(votes=VotesParty, 
                                      n=SuballianceSeats$seats_suballiance[SuballianceSeats$suballiance==AllianceInformation$suballiance],method=method) %>% 
          mutate(No_rank = row_number(-as.numeric(No))) %>%
          filter(No_rank <= SuballianceSeats$seats_suballiance[SuballianceSeats$suballiance==AllianceInformation$suballiance]) %>%
          group_by(party) %>%
          summarize(seats_party=n())
        out <- max(PartySeats$seats_party[PartySeats$party==party_input],0)  
      }
    }
  }
  return(out)
}



GetVotesRequiredClosedList <- function(data_input,
                                          j,
                                          n,
                                          i,
                                          method,
                                          threshold,
                                          convcrit=0.001){
    
    # 1. Iteration 1: Find interval that changes the election status of seat i
  
    votessim <- 0.5  # note: start with 1 vote in first iteration (first line after while loop)
    seats <- 0
    
    while (seats < i) {
      votessim <- votessim * 2
      data_input$votes_j[data_input$party==j] <- votessim
      
      if (method=="SainteLague" | method=="dHondt" ){
        seats <- GetSeatsHighestAverage(party_input=j,
                                        data_input=data_input, 
                                        n=n,
                                        method=method,
                                        threshold=threshold)
      }
      
    }
    
    
    # 2. Iteration 2: Find exact votessim that changes the election status of seat i
  
    xlow <- votessim / 2
    xhigh <- votessim
    
    
    while ((xhigh - xlow) > convcrit ) { #& i_new<10
      votessim <- (xlow + xhigh) / 2
      
      data_input$votes_j[data_input$party==j] <- 
        votessim
      if (method=="SainteLague" | method=="dHondt"){
        seats <- GetSeatsHighestAverage(party_input=j,
                                        data_input=data_input, 
                                        n=n,
                                        method=method,
                                        threshold=threshold)
      }
      
      if (seats < i) {
        xlow <- votessim
      }
      if (seats >= i) {
        xhigh <- votessim
      }
    }
    return(votessim)
  }




GetVotesRequiredOpenList <- function(data_input,
                                      i,
                                      j,
                                      n,
                                      method="dHondt",
                                      threshold,
                                      convcrit=0.001){
  
  
  # 1. Define initial values of candidate votes and party votes
  
  votes_h_initial <- data_input$votes_h[data_input$id_cand==i]
  votes_j_initial <- data_input$votes_j[data_input$party==j][1]
  
  
  # 2. Iteration 1: Find interval that changes the election status of candidate i
  
  # a) Define initial values of votes for candidate i (votessim), number of seats
  #    for party j (seats),  the rank of a candidate on their party list (rank_h), 
  #    and the election status (elected)
  
  votessim <- 0.5  # note: start with 1 vote in first iteration (first line after while loop)
  seats <- 0
  rank_h <- length(data_input$votes_j[data_input$party==j])
  elected <- 0
  out <- NULL
  
  while (elected==0) {
    
    # b) Double candidate votes (votes_h) and party votes (for all countries other
    #    than Brazil) 
    # Note: In Brazil, voters can either vote for a party (changing votes_j) or
    #       a candidate (votes_h). In our simulation, we only change votes_h and
    #       leave votes_j unchanged. 
    
    votessim <- votessim * 2
    data_input$votes_h[data_input$id_cand==i] <- votessim
    if (method!="Brazil"){
      data_input$votes_j[data_input$party==j] <- votes_j_initial-votes_h_initial+votessim
    }
    
    # c) Calculate a party's number of seats (for all countries other than Brazil)
    #    or election status of a candidate i (for Brazil)
    
    if (method=="SainteLague" | method=="dHondt"){
      seats <- GetSeatsHighestAverage(party_input=j,
                                      data_input=data_input, 
                                      n=n,
                                      method=method,
                                      threshold=threshold)
      
    }
    
    if (method=="Hare"){ #continue here for Honduras
      seats <- GetSeatsHare(party_input=j,
                            data_input=data_input,
                            n=n,
                            method=method,
                            threshold=threshold)
      
    }
    
    if (method=="Brazil"){ 
      
      elected <- GetElectionStatusBrazil(data_input=data_input,
                                         n=sum(data_input$elected),
                                         method=method,
                                         threshold=threshold)%>% 
        arrange(-elected) %>%
        filter(id_cand==i) %>%
        select(elected) %>%
        pull()
      
      elected <- ifelse(length(elected)==1,elected,0)
      # Note: replace elected status for those who do not fulfill threshold
      
    }
    
    # d) Calculate election status of candidate i (for all countries other than Brazil)
    
    if (method %in% c("SainteLague","dHondt","Hare")){ 
      data_input_party <- data_input[data_input$party==j,]
      row_number <- row_number(data_input$votes_h[data_input$party==j])
      rank_h <- rank(-data_input$votes_h[data_input$party==j]+row_number/1000)[data_input_party$party==j & data_input_party$id_cand==i]
      elected <- ifelse(rank_h > seats,0,1)
    }
    out <- c(out,votessim)
  }
  
  # 3. Iteration 2: Find exact votessim that changes the election status of candidate i
  
  # a) Define new interval and update votes
  
  xlow <- votessim / 2
  xhigh <- votessim
  
  while ((xhigh - xlow) > convcrit) {
    votessim <- (xlow + xhigh) / 2
    data_input$votes_h[data_input$id_cand==i] <- votessim
    #print(paste0("votessim:",votessim))
    
    if (method %in% c("SainteLague","dHondt","Hare")){
      data_input$votes_j[data_input$party==j] <- votes_j_initial-votes_h_initial+votessim
    }
    
    # b) Calculate a party's number of seats (for all countries other than Brazil)
    #    or election status of a candidate i (for Brazil)
    
    
    if (method %in% c("SainteLague","dHondt","Hare")){ 
      seats <- GetSeatsHighestAverage(party_input=j,
                                      data_input=data_input, 
                                      n=n,
                                      method=method,
                                      threshold=threshold)
    }
    
    if (method=="Brazil"){ 
      
      elected <- GetElectionStatusBrazil(data_input=data_input,
                                         n=sum(data_input$elected),
                                         method=method,
                                         threshold=threshold)%>% 
        filter(id_cand==i) %>%
        select(elected) %>%
        pull()
      
      #print(paste0("elected:",elected))
      
      
      elected <- ifelse(length(elected)==1,elected,0)
      
      
    }
    
    # c) Calculate election status of candidate i (for all countries other than Brazil)
    
    if (method %in% c("SainteLague","dHondt","Hare")){ 
      data_input_party <- data_input[data_input$party==j,]
      row_number <- row_number(data_input$votes_h[data_input$party==j])
      rank_h <- rank(-data_input$votes_h[data_input$party==j]+row_number/1000)[data_input_party$party==j & data_input_party$id_cand==i]
      elected <- ifelse(rank_h > seats,0,1)
    }
    
    # d) Update interval for search
    
    if (elected==0) {
      xlow <- votessim
    }
    if (elected==1) {
      xhigh <- votessim
    }
    out <- c(out,votessim)
    
  }
  
  return(votessim)
}


PrepareData <- function(data,
                        election_cycle_name = "year",
                        district_name, 
                        system,
                        party_name = NULL,
                        votes_j_name = NULL,
                        alliances=FALSE,
                        additional_vars="",
                        votes_h_name=votes_j_name, 
                        cand_id_name = NULL, 
                        rank_name = NULL,
                        alliance_name=party_name, 
                        suballiance_name=party_name) {
  
  # This function translates the variable names from a dataset to the 
  # corresponding variables names used in the functions
  # system:     open or closed list pr system
  # alliances:  are alliances allowed (TRUE) or not (FALSE)
  
  # 1. Check whether necessary variable names are specified
  
  # a) Variable names that always need to be specified
  
  if (is.null(system)) {
    stop("Information whether it is an open or closed-list system is missing. Please specify system")
  }
  
  if (is.null(district_name) || is.null(data[[district_name]])) {
    stop("District name variable is missing or not found in the dataset. Please specify district_name")
  }
  
  if (is.null(party_name) || is.null(data[[party_name]])) {
    stop("Party name variable is missing or not found in the dataset. Please specify party_name.")
  }
  
  if (is.null(votes_j_name) || is.null(data[[votes_j_name]])) {
    stop("Votes_j variable is missing or not found in the dataset. Please specify votes_j_name.")
  }
  
  # b) Variable names that need to be specified only in open- or closed list systems
  
  if (system == "open") {
    if (is.null(cand_id_name) || is.null(data[[cand_id_name]])) {
      stop("Candidate ID is required for open list systems. Please specify cand_id_name.")
      if (is.null(votes_h_name) || is.null(data[[votes_h_name]])) {
        stop("Candidate votes variable is required for closed list systems. Please specify votes_h_name.")
      }
    }
  }
  
  if (system == "closed") {
    if (is.null(rank_name) || is.null(data[[rank_name]])) {
      stop("Rank name variable is required for closed list systems. Please specify votes_h_name.")
    }
  }
  
  # 2. Rename general variables (year, district, party votes)
  
  
  data$year <- data[[election_cycle_name]]
  data$district <- data[[district_name]]
  data$votes_j <- data[[votes_j_name]]
  
  # 3. Rename party, alliance, suballiances variables
  
  data$party_orig <- paste(as.character(data$district), 
                           as.character(data$year), 
                           data[[party_name]])
  
  data$alliance_orig <- as.character(data[[alliance_name]]) 
  data$suballiance_orig <- as.character(data[[suballiance_name]]) 
  data <- data %>% mutate(alliance=ifelse(alliance_orig%in%c("","."),
                                          party_orig,alliance_orig),
                          suballiance=ifelse(suballiance_orig%in%c("","."),
                                             party_orig,suballiance_orig),
                          alliance_dummy=ifelse(alliance_orig%in%c("","."),
                                                0,1),
                          suballiance_dummy=ifelse(suballiance_orig%in%c("","."),
                                                   0,1))
  data$party <- as.numeric(as.factor(as.character(data$party_orig)))
  data$alliance <- as.numeric(as.factor(as.character(data$alliance)))
  data$suballiance <- as.numeric(as.factor(as.character(data$suballiance)))
  
  
  if (alliances==FALSE){
    data$party <- as.numeric(as.factor(as.character(data$party_orig)))
    data$alliance <- as.numeric(as.factor(as.character(data$party_orig)))
    data$suballiance <- as.numeric(as.factor(as.character(data$party_orig)))
    data$alliance_dummy <- 0
    data$suballiance_dummy <- 0
  }
  
  # 4. Rename variables for open-list systems
  
  if (system=="open"){ 
    data$votes_h <- data[[votes_h_name]]
    data$id_cand <- as.numeric(as.factor(paste0(data$year,data$district,data$party,data$alliance,data[[cand_id_name]])))
  }
  
  # 5. Select core variables for calculation of closeness measure and additional
  #    variables for open-list systems and return data
  
  if (system=="open"){
    
    if (additional_vars[1]!=""){
      dataout <- data %>%
        select(district, year,party, id_cand, elected, votes_j, votes_h, alliance, 
               suballiance, alliance_dummy, suballiance_dummy, ends_with(additional_vars))
    }
    if (additional_vars[1]==""){
      dataout <- data %>%
        select(district, year, party, id_cand, elected, votes_j, votes_h, alliance,  
               suballiance, alliance_dummy, suballiance_dummy)
    }
    return(dataout %>% arrange(year,district,party,id_cand))  
  }
  
  # 5.  Rename ranking variable for closed-list systems, select core variables 
  #     for calculation of closeness measure, select additional variables  
  #     and return data
  
  if (system=="closed"){
    data$seat <- data[[rank_name]]
    
    if (additional_vars[1]!=""){
      dataout <- data %>%
        select(district, year,party, seat, elected, votes_j, alliance, suballiance, 
               alliance_dummy, suballiance_dummy, ends_with(additional_vars))
    }
    if (additional_vars[1]==""){
      dataout <- data %>%
        select(district, year, party, seat, elected, votes_j, alliance, suballiance, 
               alliance_dummy, suballiance_dummy)
    }
    return(dataout %>% arrange(year,district,party,seat))  
  }  
}

