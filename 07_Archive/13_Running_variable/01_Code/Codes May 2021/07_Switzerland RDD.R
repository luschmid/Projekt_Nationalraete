
library(scales)
library(tidyverse)
library(xtable)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(readstata13)
library(rdrobust)
library(ggrepel)
library(dummies)
library(janitor)
library(ggExtra)

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))
try(setwd("E:/Projekt Nationalräte"))

rm(list = ls(all = TRUE))
source("./03_Code/fine_grid.R") # ggplot layers


###################
# (A) Read-in data
###################


data <- read.dta13("./02_Processed_data/nationalraete_1931_2015.dta") 
cantons <- as.character(paste("canton_",levels(as.factor(data$canton)),sep=""))
#data <- data %>% filter((elected == 1 & votemargin_rel > 0) | (elected == 0 & votemargin_rel < 0))
data %>% tabyl(canton,year) 
  
data$votemargin_rel1<- data$votemargin/data$eligible_cant
data$votemargin_rel2<- data$votemargin/data$voters_cant 
data$votemargin_rel3<- data$votemargin*data$no_seats/data$eligible_cant
data$votemargin_rel4<- data$votemargin*data$no_seats/data$voters_cant

data$votemargin_rel1_L1<- data$votemargin_L1/data$eligible_cant_L1
data$votemargin_rel2_L1<- data$votemargin_L1/data$voters_cant_L1
data$votemargin_rel3_L1<- data$votemargin_L1*data$no_seats_L1/data$eligible_cant_L1
data$votemargin_rel4_L1<- data$votemargin_L1*data$no_seats_L1/data$voters_cant_L1

cor(data$votemargin_rel3,data$votemargin_rel)

data %>% select(voters_cant,eligible_cant,votemargin_rel,votemargin_rel1)

data_hn <- read.dta13("./02_Processed_data/15_Elections_Honduras/elections_hn_final.dta") %>%
                  mutate(Sex_num=ifelse(Sex=="H",1,0),
                         Department = sub(" ","",Department),
                         Department = sub(" ","",Department),
                         Department = sub(" ","",Department)) %>%
                  rename(year=Year,elected=Elected)
deps <- as.character(paste("Department_",levels(as.factor(data_hn$Department)),sep=""))

data_hn$votemargin_rel1<- data_hn$Votemargin/data_hn$total_electoral
data_hn$votemargin_rel2<- data_hn$Votemargin/data_hn$votaron
data_hn$votemargin_rel3<- data_hn$Votemargin*data_hn$no_seats/data_hn$total_electoral
data_hn$votemargin_rel4<- data_hn$Votemargin*data_hn$no_seats/data_hn$votaron

data_hn$votemargin_rel1_L1<- data_hn$Votemargin_L1/data_hn$total_electoral_L1
data_hn$votemargin_rel2_L1<- data_hn$Votemargin_L1/votaron_L1 
data_hn$votemargin_rel3_L1<- data_hn$Votemargin_L1*data_hn$no_seats_L1/data_hn$total_electoral_L1 
data_hn$votemargin_rel4_L1<- data_hn$Votemargin_L1*data_hn$no_seats_L1/data_hn$votaron_L1



data_hn %>% tabyl(Sex,Department)

data$votemargin_rel1_rounded <- round(data$votemargin_rel1,3)
data %>% tabyl(elected,votemargin_rel1_rounded)

data %>% filter(elected==0 & votemargin_rel1>=0)%>% select(name,firstname,year,votemargin_rel1,votemargin_rel3,party,votes_h,pvotes)
data %>% filter(elected==1 & votemargin_rel1<0)%>% select(name,firstname,year,votemargin_rel1,votemargin_rel3,party,votes_h,pvotes)

data %>% filter(name=="Musy")%>% select(year,name,firstname,elected,votemargin_rel1,votemargin_rel3)

###################
# (B) Functions
###################

GetMissingValues <- function(datainput,varnames){
  for (i in 1:length(varnames)){
    datainput$outcome <- datainput[[varnames[i]]]
    print(paste0("No. non-missing:",varnames[i],": ",datainput %>%summarise(count = sum(is.na(outcome)==F))%>%pull()))
  }  
}

MakeDummyVar <- function(v, prefix = "") {
  #########################################################
  # This function creates a dummy for a certain variable v
  #########################################################
  s <- sort(unique(v))
  d <- outer(v, s, function(v, s) 1L * (v == s))
  colnames(d) <- paste0(prefix, s)
  d
}

MakeDummies <- function(data,v, prefix = "") {
  #############################################################################
  # This function creates a dummy for a certain variable v in a dataframe data
  #############################################################################
  return(cbind(data,MakeDummyVar(data[[v]],prefix=prefix)))
}

AddStars <- function(datainput) {
  ######################################################################################
  # This function adds stars to coefficients 
  # Input: datainput: vector of estimate (first element) and std. error (second element)
  # Outuput: Vector with added stars
  ########################################################################################
  if (datainput[1]==0){coef <-datainput[1]}
  else{coef <- paste0(sprintf("%.3f",datainput[1]), "$^{",symnum(2*pnorm(-abs(datainput[1]/datainput[2])), corr = FALSE,
                                                                     cutpoints = c(0, .01,.05, .1, 1),
                                                                     symbols = c("***","**","*","")),"}$")}
  return(c(coef,paste0("(",sprintf("%.3f",datainput[2]),")")))
}

AddSquares <- function(datainput) {
  ######################################################################################
  # This function adds stars to coefficients 
  # Input: datainput: vector of estimate (first element) and std. error (second element)
  # Outuput: Vector with added stars
  ########################################################################################
return(paste0("[",sprintf("%.3f",datainput),"]"))
}


SummarizeData <- function(datainput,vars, districtvariable,  runningvariable,
                          bin_digits = 2, year_start = 1931, year_end = 2015, only_startyear = 0, only_participation = 0) {
  ###############################################################################################################################################
  # This function aggregates the raw data from candidates in several districts into aggregates in a certain bin of the running variable votemargin_rel.
  # Output: List with first element: raw data plus dummies for district
  #                    second element: means of raw data
  # Arguments:   bin_digits are the number of digits
  #              only_startyear=1 means that only candidates are considered who run the first time
  #              only_participation=1 means that means are only calculated for candidates who run in the next election
  ###############################################################################################################################################
  if (only_startyear == 1) {
    datainput <- datainput %>% filter(year == year_min & cand_before1931 == 0)
  }
  if (only_participation == 1) {
    datainput <- datainput %>% filter(participation_F1 == 1)
  }
  
  datainput <- MakeDummies(data=datainput,v=districtvariable, prefix = paste(districtvariable,"_",sep=""))
  datainput$runningvariable <- as.numeric(as.character(datainput[[runningvariable]]))
  
  data_estimation <- datainput %>%
    mutate(votemargin_rel_rounded = floor(runningvariable * 10^bin_digits)/10^bin_digits) %>%
    filter(year >= year_start & year <= year_end) %>%
    select(c(runningvariable,votemargin_rel_rounded,vars, starts_with( paste(districtvariable,"_",sep=""))))%>%
    group_by(votemargin_rel_rounded) %>%
    summarise_all(mean, na.rm = TRUE)%>%
    mutate(elected = ifelse(votemargin_rel_rounded >= 0, 1, 0))
  return(list(datainput,data_estimation))
}

PlotHistogram <- function(datainput,window,binwidth_input=0.01,runningvariable,label=""){
  ##################################################################################################################################################
  # This function plots the distribution of the running variable votemargin_rel and saves it into our folder "05_Texts_and_presns/01_Running_variable"
  ##################################################################################################################################################
  break_points <- unique(c(-rev(seq(0,window,binwidth_input)),seq(0,window,binwidth_input)))
    ggplot(datainput, aes_string(x = runningvariable)) +
    geom_histogram(aes(y = ..count.. / sum(..count..)), fill = "grey", colour = "black", 
                    breaks = break_points) +
    theme_bw(base_size = 42) +
    scale_x_continuous(limits = c(-window, window), label = percent) +
    geom_vline(xintercept = 0) +
    ylab("") + xlab("")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    

  
  ggsave(file = paste("./05_Texts_and_presns/01_Running_variable/figures/votemargin_rel_bw_", window, "_",label,".pdf", sep = ""), width = 20, height = 12)
}



PlotDiscontinuity <- function(data_all, y, ylab=y, kernel="triangular",binselect="qsmv",runningvariable,window = 1, polynomials = 4,label="") {
  #################################################################################################################
  # This function generates an RDD plot for one variable and saves it  into our folder "05_Texts_and_presns/01_Running_variable"
  #################################################################################################################
  data_all$runningvariable <- as.numeric(as.character(data_all[[runningvariable]]))
  data_all$y <- as.numeric(as.character(data_all[[y]]))
  data_all<- data_all[is.na(data_all$elected_F1)==F,]
  rdplot(data_all$y,data_all$runningvariable,
         kernel = kernel, 
         h=window,
         support=c(-window,window),
         x.lim=c(-window,window),
         p=polynomials,
         binselect=binselect,
         col.dots="grey",
         col.lines="black",
         title = "",
         x.label="\n  Relative vote margin",
         y.label=ylab) 
    
  ggsave(file = paste("./05_Texts_and_presns/01_Running_variable/figures/", y, "_",label,".pdf", sep = ""), width = 30, height = 20)
}


PlotDiscontinuities <- function(data,vars,varlabels,districtvariable,runningvariable,
                                kernel="triangular",binselect="qsmv", 
                                year_start , year_end , only_startyear = 0, polynomials = 4,
                                only_participation = 0,label="",window=0.01){
  #################################################################################################################
  # This function generates RDD plots for multiple variables by calling the function PlotDiscontinuity. 
  # Important note: The dots of the plot are calculated based on binned data (generated in SummarizeData), 
  #                  the blue line is based on all observations and not on the binned dots. 
  # Arguments:    data: raw data
  #               vars: variables to plot
  #               varlabels: label of variable in plot (y-axis), should be the same order as argument vars
  #               districtvariable: name of the district variable
  #               bin_digits: how many digits should be used for binned dots  
  #               year_start: start of sample period
  #               year_end = end of sample period
  #               only_startyear=1 means that only candidates are considered who run the first time
  #               only_participation=1 means that means are only calculated for candidates who run in the next election
  #               label: use if for additional label in name of pdf plot (example hn for Honduras)
  #################################################################################################################  
    for (i in 1:length(vars)) {
    PlotDiscontinuity(data_all=data, y = vars[i], runningvariable=runningvariable,
                      ylab = varlabels[i],  window = window, polynomials = polynomials ,
                      label=label,binselect = binselect,kernel=kernel)
  }
}



RDEstimation <- function(datainput,outcomevariable,runningvariable,districtvariable,fuzzyvariable=0,pvector,hvector,latex=0){
  #################################################################################################################
  # This function calculates the RD estimate using the RD robust package. 
  # Arguments:    datainput: raw data
  #               outcomevariable: name of outcome variable
  #               runningvariable: name of running variable
  #               fuzzyvariable: name of fuzzy treatment variable (if missing, sharp RD estimate is calculated)
  #               districtvariable: name of district variable (string)
  #               pvector: vector of polynomials
  #               hvector: vector of bandwidths (should be the same lenght as pvector)
  #               latex: a latex output is generated if 1, a table otherwise
  #################################################################################################################
  
  datainput <- MakeDummies(data=datainput,v=districtvariable, prefix = paste(districtvariable,"_",sep=""))
  
  result_file <- as.data.frame(matrix(NA, 3, length(pvector)))
  
 for (i in 1:length(pvector)) {
    if (fuzzyvariable==0){ # Case 1: Sharp RDD
      if (hvector[i] == 1) {rd_model <- rdrobust(datainput[[outcomevariable]], datainput[[runningvariable]], p = pvector[i],masspoints="off",cluster=datainput$ID_pers)}
      else {rd_model <- rdrobust(datainput[[outcomevariable]], datainput[[runningvariable]], p = pvector[i], h = hvector[i],masspoints="off",cluster=datainput$ID_pers)}
      result_file[, i] <- c(AddStars(c(rd_model$coef[3],rd_model$se[3])), sprintf("%.3f", rd_model$bws[1, 1]))
    }
    if (fuzzyvariable!=0){ # Case 2: Fuzzy RDD
      if (hvector[i] == 1) {rd_model <- rdrobust(datainput[[outcomevariable]], datainput[[runningvariable]], fuzzy = datainput[[fuzzyvariable]], p = pvector[i],masspoints="off",cluster=datainput$ID_pers)}
      else {rd_model <- rdrobust(datainput[[outcomevariable]], datainput[[runningvariable]], p = pvector[i], fuzzy = datainput[[fuzzyvariable]], h = hvector[i],masspoints="off",cluster=datainput$ID_pers)}
      result_file[, i] <- c(AddStars(c(rd_model$coef[1],rd_model$se[1])), sprintf("%.3f", rd_model$bws[1, 1]))
    }
  } 
  
  rownames(result_file) <- c("Estimate", "", "Bw")
  colnames(result_file) <- paste0("(",c(1:length(pvector)),")")
  if (latex==1){return(xtable(result_file))}
  else{return(result_file)}
}


GenerateRDDTable <- function(datainput,outcomevariables,outcomevariables_label=outcomevariables,runningvariable,
                             districtvariable,fuzzyvariable=0,pvector,hvector,latex=0){
  #################################################################################################################
  # This function generates the table for the paper using the function RDEstimation. 
  # Arguments:    datainput: raw data
  #               outcomevariables: name of outcome variables (vector)
  #               outcomevariables_labels: label of outcome variables (vector)
  #               runningvariable: string name of running variable 
  #               fuzzyvariable: string name of fuzzy treatment variable (if missing, sharp RD estimate is calculated)
  #               districtvariable: name of district variable (string)
  #               pvector: vector of polynomials
  #               hvector: vector of bandwidths (should be the same lenght as pvector)
  #               latex: a latex output is generated if 1, a table otherwise
  #################################################################################################################

output <- data.frame()  
for (i in 1:length(outcomevariables)){
  print(i)
  if (i<length(outcomevariables)){
  result <- RDEstimation(datainput=datainput,
                         outcome=outcomevariables[i],
                         runningvariable=runningvariable,
                         districtvariable=districtvariable,
                         pvector=pvector,
                         hvector= hvector)[c(1:2),]
  rownames(result) <- c(outcomevariables_label[i], strrep(" ",i))
  output <- rbind(output,result)
  }
  else{
  result <- RDEstimation(datainput=datainput,
                           outcome=outcomevariables[i],
                           runningvariable=runningvariable,
                           districtvariable=districtvariable,
                           pvector=pvector,
                           hvector= hvector)
  rownames(result) <- c(outcomevariables_label[i], "","Bandwidth")
  output <- rbind(output,result,as.vector(pvector))
  rownames(output) <- c(rownames(output)[1:(dim(output)[1]-1)],"Polynomials")
}
}
if (latex==1) return(xtable(output))
else{return(output)}
}




# test functions
#data_ave <- SummarizeData(data, vars=vars,
#                          districtvariable = districtvariable, bin_digits = 3, 
#                          year_start = 1931, year_end = 2015, only_startyear = 0, only_participation = 0)
#PlotDiscontinuity(data_all=data_ave[[1]],data_means=data_ave[[2]], y = "canton_ZH", ylab = "Canton of Zurich", group = "elected", window = 0.01, polynomials = 2)




PlotEstimateDifferentBandwidths <- function(outcomevariable,datainput,runningvariable,districtvariable,fuzzyvariable=0,pvector,hvector,label=""){
  #################################################################################################################
  # This function calculates the RD for different bandwidths using the RD robust package. 
  # Arguments:    datainput: raw data
  #               outcomevariable: name of outcome variable
  #               runningvariable: name of running variable
  #               fuzzyvariable: name of fuzzy treatment variable (if missing, sharp RD estimate is calculated)
  #               districtvariable: name of district variable (string)
  #               pvector: vector of polynomials
  #               hvector: vector of bandwidths
  #               label: use if for additional label in name of pdf plot (example hn for Honduras)
  #################################################################################################################
  
for (i in 1:length(pvector)){
results <- RDEstimation(datainput=data,outcome=outcomevariable,
                        runningvariable=runningvariable,
                        districtvariable = districtvariable,
                        pvector=rep(pvector[i],length(hvector)),
                        hvector= hvector,
                        latex=0)
df <- data.frame(beta= as.numeric(gsub("\\*|\\$|\\{|\\}|\\^", "", results[1,])),
                 beta_se=as.numeric(gsub("\\(|\\)", "", results[2,])),
                 bandwidth=as.numeric(results[3,]))

ggplot(data=df,aes(x=bandwidth,y=beta))+
  geom_point(size=6)+
  geom_errorbar(aes(ymax = beta+1.96*beta_se, ymin = beta-1.96*beta_se))+
  theme_bw(base_size = 82)+
  scale_x_continuous(label = percent) +
  xlab("Bandwidth") + ylab("\n  Estimate")

ggsave(file = paste("./05_Texts_and_presns/01_Running_variable/figures/bws_",outcomevariable, "_pol",pvector[i],"_",label,".pdf", sep = ""), width = 30, height = 20)


}
}


PlotsAll <- function(outcomevariables,datainput,runningvariables,varlabel=outcomevariables,districtvariable,year_start , year_end,fuzzyvariable=0,hvector,label=""){
  #################################################################################################################
  # This function creates plots for (A) RD estimate for different bandwidths; (B) Classical RD plot for different windows and possibly different polynomials
  # Arguments:    datainput: raw data
  #               outcomevariable: name of outcome variable (note: first variable is main outcome variable that is used for bw choice in part (BS))
  #               runningvariables: name of different running variables
  #               fuzzyvariable: name of fuzzy treatment variable (if missing, sharp RD estimate is calculated)
  #               districtvariable: name of district variable (string)
  #               label: use if for additional label in name of pdf plot (example hn for Honduras)
  #################################################################################################################
  
  
  
  
  # (A) Plot RD estimate for different bandwidths
  
  for (i in 1:length(runningvariables)){
    rv_min <-summary(datainput[[runningvariables[i]]])[1]
    rv_max <-summary(datainput[[runningvariables[i]]])[6]
    rv_0.25q <- quantile(datainput[[runningvariables[i]]][datainput$elected==1],probs=0.025,na.rm=T)
    
    print(paste("RD bw choice, running variable: ",runningvariables[i]))
    lapply(outcomevariables,PlotEstimateDifferentBandwidths,
           datainput=datainput,
           runningvariable=runningvariables[i],
           districtvariable = districtvariable,
           pvector=2,
           label=paste(label,"_",runningvariables[i],sep=""),
           hvector= seq(rv_0.25q,max(abs(rv_min),rv_max),length.out=50)) # 50 bins between min and max
  }
  
  
  
  # (B) Classical RD plot for different windows and possibly different polynomials
  
  for (i in 1:length(runningvariables)){
    #print(paste("Classical RD plot, RV: ",runningvariables[i]))
    rdbw <- rdbwselect(datainput[[outcomevariables[1]]], 
                       datainput[[runningvariables[i]]], 
                       p = 1,masspoints = F)
    rv_min <-summary(datainput[[runningvariables[i]]])[1]
    rv_max <-summary(datainput[[runningvariables[i]]])[6]
    
    windows <- c(max(abs(rv_min),rv_max),rdbw$bws[1]/4)
    windows_label <- c("full","bw_divided_by_4")
    pvector <- c(2,1)
    for (k in 1:length(windows)){
      print(paste("Classical RD plot, windows: ",windows[k]))
      PlotDiscontinuities(data=datainput,vars=outcomevariables,
                          runningvariable=runningvariables[i],
                          varlabels=varlabel, districtvariable = districtvariable,
                          bin_digits = 3, 
                          window=windows[k],
                          year_start = year_start, year_end = year_end, 
                          only_startyear = 0, only_participation = 0,
                          polynomials = pvector[k],
                          label=paste(label,"_",runningvariables[i],"_win",windows_label[k],"_pol",pvector[k],sep=""))
    }
  }
  
}


GenerateLongTable <- function(outcomevariables,datainput,runningvariables,varlabel=outcomevariables,districtvariable,year_start,year_end,fuzzyvariable=0,label="",latex=0){
  #################################################################################################################
  # This function creates plots for (A) RD estimate for different bandwidths; (B) Classical RD plot for different windows and possibly different polynomials
  # Arguments:    datainput: raw data
  #               outcomevariables: name of outcome variables (note: main outcome variable should be the first in this vector b/c bw is chosen based on mainoutcome variable)
  #               runningvariables: name of different running variables
  #               fuzzyvariable: name of fuzzy treatment variable (if missing, sharp RD estimate is calculated)
  #               districtvariable: name of district variable (string)
  #               windows: vector of windows for classical RD plot
  #               label: use if for additional label in name of outtable (example hn for Honduras)
  #################################################################################################################
  datainput <- datainput[is.na(datainput$elected_F1)==F,]
  datainput <- MakeDummies(data=datainput,v=districtvariable, prefix = paste(districtvariable,"_",sep=""))
  
  for (i in 1:length(runningvariables)){
    output_all <- data.frame()  
    output <- data.frame()  
    print(paste0("Running variable:", runningvariables[i]))
    datainput$runningvariable <- datainput[[runningvariables[i]]]
    for (j in 1:length(outcomevariables)){
      print(paste0("Outcome variable:", outcomevariables[j]))
      datainput$outcome <- datainput[[outcomevariables[j]]]
      if (j==1){ # bandwidth is chosen based on optimal bandwidth of outcome variable which shall be the first in the vector outcomevariables
     rdbw1 <- rdbwselect(datainput$outcome, datainput$runningvariable,
                          kernel="triangular", bwselect="mserd",p = 1,cluster=datainput$ID_pers)
      rdbw2 <- rdbwselect(datainput$outcome, datainput$runningvariable, 
                          kernel="triangular", bwselect="mserd",p = 2,cluster=datainput$ID_pers)
      if (rdbw1$bws[1]!=rdbw1$bws[2]){print("Bandwidths (left and right) are not the same (pol=1)")}
      if (rdbw2$bws[1]!=rdbw2$bws[2]){print("Bandwidths (left and right) are not the same (pol=2)")}
      window <- c(rdbw1$bws[2],rdbw1$bws[2]/2,rdbw2$bws[2],rdbw2$bws[2]/2)
      }
      
      if (j<length(outcomevariables)){
        result <- GenerateRDDTable(datainput=datainput,
                                         outcomevariables=outcomevariables[j],
                                         outcomevariables_label=varlabel[j],
                                         runningvariable="runningvariable",
                                         districtvariable = districtvariable,
                                         fuzzyvariable=fuzzyvariable,
                                         pvector=c(1,1,2,2),
                                         hvector= window,
                                         latex=0)[c(1:2),]
        output <- cbind(c(sprintf("%.3f", mean(datainput$outcome,na.rm=T)),AddSquares(sd(datainput$outcome,na.rm=T))),result)
        rownames(output) <- c(paste0(gsub("\\_","",varlabel[j]),strrep(" ",(i-1)*length(outcomevariables)+j)), strrep(" ",(i-1)*length(outcomevariables)+j))
        colnames(output) <- c(paste0("(",1:dim(output)[2],")")) 
        output_all <- rbind(output_all,output)
      }
      else{
        result <- GenerateRDDTable(datainput=datainput,
                                    outcomevariables=outcomevariables[j],
                                   outcomevariables_label=varlabel[j],
                                   runningvariable="runningvariable",
                                    districtvariable = districtvariable,
                                    fuzzyvariable=fuzzyvariable,
                                    pvector=c(1,1,2,2),
                                    hvector= c(rep(window,2)),
                                    latex=0)
        if (i<length(runningvariables[i])) result <- result[c(1:3),]
        rownames(result) <- c(paste0(sub("_","",varlabel[j]),strrep(" ",(i-1)*length(outcomevariables)+j)), strrep(" ",(i-1)*length(outcomevariables)+j),paste0("Bandwidth",strrep(" ",(i-1)*length(outcomevariables)+j)),paste0("Polynomials",strrep(" ",(i))))
        df_label <- data.frame(var1="",var2="",var3="",var4="",var5="")
        rownames(df_label) <- ""
        #rownames(df_label) <- paste0("\\textbf{Running variable: $",gsub("_","\\_",runningvariables[i], fixed=TRUE),"$}")
        colnames(df_label) <- paste0("(",c(1:dim(df_label)[2]),")")
        output <- cbind(c(sprintf("%.3f", mean(datainput$outcome,na.rm=T)),AddSquares(sd(datainput$outcome,na.rm=T)),"",""),result)
        #rownames(output) <- c(paste0(gsub("\\_","",varlabel[j]),strrep(" ",(i-1)*length(outcomevariables)+j)), strrep(" ",(i-1)*length(outcomevariables)+j))
        colnames(output) <- c(paste0("(",1:dim(output)[2],")")) 
        output <- rbind(df_label,output)
        output_all <- rbind(output_all,output)
      }

      
    }
 if (latex==1) {print(xtable(output_all), 
                        file = paste0("./05_Texts_and_presns/01_Running_variable/tables/balance_tables_",runningvariables[i],"_",label,".tex"),
                        only.contents = T,
                        booktabs=T,
                        sanitize.text.function = function(x) {x})}
   else{return(output_all)
}
}
}

### RANDOM ELECTIONS ################## 

RandomElectionResult <- function(nseats_district,eligibles_district,turnout_district,voter_preferences_district,
                                 individual_votes_distribution){
  voters_district <- eligibles_district*turnout_district
  votes_district <- voters_district*nseats_district  
  out <- data.frame()
  voter_preferences_district <- voter_preferences_district+rnorm(length(voter_preferences_district),0,0.01)
  votes_parties=as.vector(table(sample(c(1:length(voter_preferences_district)), size=votes_district, replace=TRUE, 
                                       prob=voter_preferences_district)))
  #print(votes_parties)
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



RandomElectionResult(nseats_district,eligibles_district,turnout_district,voter_preferences_district,
                                 individual_votes_distribution)

# nseats_district=5;eligibles_district=100;turnout_district=0.2;voter_preferences_district=c(0.1,0.3,0.6);
# individual_votes_distribution=rep(1/5,5)

SimulateElection<- function(ndistricts,nparties, nseats, eligibles, turnout,voter_preferences,
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
  ElectionResult <- RandomElectionResult(nseats_district=nseats[j],eligibles_district=eligibles[j],turnout_district=turnout[j],
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



##############################
# (C) RD estimation
##############################


# (i) Histogram of running variable


PlotHistogram(datainput=data %>% filter(is.na(elected_F1)==F),window=0.11,binwidth_input=0.002,
              runningvariable="votemargin_rel1",label="ch")
PlotHistogram(datainput=data_hn %>% filter(is.na(Elected_F1)==F),window=0.1,binwidth_input=0.002,
              runningvariable="votemargin_rel1",label="hn")

PlotHistogram(datainput=data,binwidth_input=0.005,window=2*0.112,runningvariable="votemargin_rel1",label="ch")


# (ii) RDD graphs for balance tests and outcome variables

varname <- c("elected_F1", "participation_F1", "sex", "year","no_seats","age","votemargin_L1",cantons)
varlabel <- c("Elected in t+1", "Participation in t+1", "Male", "Year","Seats of canton","Age","votemargin_L1",cantons)

PlotDiscontinuities(data=data,vars="age",
                    runningvariable="votemargin_rel1",
                    districtvariable = "canton",
                    varlabels = "age",
                    bin_digits = 3, 
                    year_start = 1931, year_end = 2015, only_startyear = 0, only_participation = 0,
                    polynomials = 1,
                    label="0.018",
                    window=0.039/2 )

GetMissingValues(datainput=data,varnames=c("elected_F1", "participation_F1", "sex","canton", "year","no_seats","age","votemargin_rel1","votemargin_rel2","votemargin_rel3","votemargin_rel4","votemargin_rel1_L1","votemargin_rel2_L1","votemargin_rel3_L1","votemargin_rel4_L1"))

varname <- c("elected_F1", "participation_F1", "sex", "year","no_seats","age","votemargin_rel1_L1","votemargin_rel2_L1","votemargin_rel3_L1","votemargin_rel4_L1",cantons)
varlabel <- c("Elected in t+1", "Participation in t+1", "Male", "Year","Seats of canton","Age","votemargin_rel1_L1","votemargin_rel2_L1","votemargin_rel3_L1","votemargin_rel4_L1",cantons)
  
PlotsAll(datainput=data,outcome=varname,
                 runningvariables=c("votemargin_rel1","votemargin_rel2","votemargin_rel3","votemargin_rel4"),
                 districtvariable = "canton",
                 pvector=c(2),
                 windows= c(1,0.01),
                 hvector=c(1),
                 year_start=1931,year_end=2015,
                 label="ch") 


GenerateLongTable(outcomevariables=varname,
                    datainput=data,
                    runningvariables=c("votemargin_rel1"),
                    varlabel=varlabel,districtvariable="canton",
                    year_start=1931,year_end=2015,
                    fuzzyvariable=0,label="ch",
                    latex=1)

PlotDiscontinuities(data=data,vars="age",varlabels=varlabel, 
                    runningvariable="votemargin_rel1",
                    districtvariable = "canton",
                    bin_digits = 3, 
                    year_start = 1931, year_end = 2015, only_startyear = 0, only_participation = 0,
                    polynomials = 2,
                    label="ch_100",
                    window=0.015)

# Honduras

varname <- c("Elected_F1", "participation_F1", "Sex_num", "year","no_seats","votemargin_rel1_L1","votemargin_rel2_L1","votemargin_rel3_L1","votemargin_rel4_L1",deps)
varlabel <- c("Elected in t+1", "Participation in t+1", "Male", "Year","no_seats","votemargin_rel1_L1","votemargin_rel2_L1","votemargin_rel3_L1","votemargin_rel4_L1",deps)
  
PlotsAll(datainput=data_hn,outcome="elected_F1",
         runningvariables=c("votemargin_rel1","votemargin_rel2","votemargin_rel3","votemargin_rel4"),
         districtvariable = "Department",
         pvector=c(2),
         windows= c(1,0.01),
         hvector=c(1),
         year_start=2009,year_end=2017,
         label="ch") 


GenerateLongTable(outcomevariables=varname,
                  datainput=data_hn,
                  runningvariables=c("votemargin_rel1","votemargin_rel2","votemargin_rel3","votemargin_rel4"),
                  varlabel=varlabel,districtvariable="Department",
                  year_start=2009,year_end=2019,
                  fuzzyvariable=0,
                  latex=1,
                  label="hn")


data %>% filter(votemargin>0 & elected==0) %>% select(ID,name,firstname, year,canton, votemargin ,elected, starts_with("vote"))
data %>% filter(votemargin<0 & elected==1) %>% select(ID,name,firstname, year,canton, votemargin ,elected, starts_with("vote"))





## Additional investigations

# investigate relationship between no_seats and turnout

data_ct <- data %>% group_by(canton,year) %>% select(canton,year,seats_cant,eligible_cant,voters_cant) %>%
        summarize_all(mean) %>%
        mutate(turnout=voters_cant/eligible_cant*100)

cor(data_ct$seats_cant,data_ct$turnout)


ggplot(data_ct,aes(x=seats_cant,y=turnout)) +
  geom_jitter()

summary(lm(turnout~seats_cant,data=data_ct))

# result: very small negative relationship between no_seats and turnout (correlation: -0.01)


# investigate the imbalance in Francisco Morazan

data_hn_pairs <- data_hn %>% mutate(votemargin_rel_abs=abs(votemargin_rel))%>%
          group_by(Department,Year,votemargin_rel_abs) %>% 
          mutate(id_pair=group_indices(),
                 nobs_pair=max(row_number())) %>%
          ungroup() %>%
          select(Candidate,Department,Year,votemargin_rel,Elected,id_pair,nobs_pair,no_seats) %>%
          arrange(Department,Year,abs(votemargin_rel))

data_pos <- data_hn_pairs %>% mutate(votemargin_rel_pos=ifelse(votemargin_rel>=0,1,0))  %>% 
        filter(nobs_pair==1 & votemargin_rel_pos==1) %>%
        group_by(no_seats) %>%
        summarize(nobs_pos=n()) 

data_neg <- data_hn_pairs %>% mutate(votemargin_rel_pos=ifelse(votemargin_rel>=0,1,0))  %>% 
  filter(nobs_pair==1 & votemargin_rel_pos==0) %>%
  group_by(no_seats) %>%
  summarize(nobs_neg=n())


data_nobs <- data_hn %>% group_by(no_seats) %>%
  summarize(nobs_all=n())

data_pos %>% full_join(data_neg,by=c("no_seats")) %>%  full_join(data_nobs,by=c("no_seats")) %>% 
  arrange(no_seats) %>% mutate(ratio=nobs_pos/nobs_neg,ratio_binding=(nobs_pos+nobs_neg)/nobs_all)  %>%
  ggplot(aes(no_seats,ratio_binding)) +
  geom_point()

# (iii) Sharp RD estimation

rdbw <- rdbwselect(data$elected_F1, data$votemargin_rel, p = 1)
rdbw$bws

rdbw <- rdbwselect(data_hn$Elected_F1, data_hn$votemargin_rel, p = 2)
rdbw$bws

# (a) Switzerland

RDEstimation(datainput=data,outcome="elected_F1",
             runningvariable="votemargin_rel1",
             districtvariable = "canton",
             pvector=c(1,1,1,2,2,2),
             hvector= c(0.078,0.039,0.078/4,0.112,0.056,0.112/4)) 

RDEstimation(datainput=datainput,
             outcome=outcomevariables[j],
                 runningvariable="votemargin_rel1",
                 districtvariable = districtvariable,
                 pvector=c(1,1,1,2,2,2),
                 hvector= c(0.078,0.039,0.078/4,0.112,0.056,0.112/4),
                 latex=0)

RDEstimation(datainput=data%>%filter(participation_F1==1),outcome="elected_F1",
             runningvariable="votemargin_rel",
             districtvariable = "canton",
             pvector=c(2, 3, 2, 3, 2, 3),
             hvector= c(1, 1, 0.05, 0.05, 0.01, 0.01)) 


GenerateRDDTable(datainput=data,
                 outcomevariables=c("elected_F1","sex","year","canton_ZH"),
                 outcomevariables_label=c("Elected in next election","Year","Male","Zurich"),
                 runningvariable="votemargin_rel",
                 districtvariable = "canton",
                 fuzzyvariable=0,
                 pvector=c(2, 3, 2, 3, 2, 3),
                 hvector= c(1, 1, 0.05, 0.05, 0.01, 0.01),
                 latex=1)









# (b) Honduras

RDEstimation(datainput=data_hn,outcome="Elected_F1",
             runningvariable="votemargin_rel",
             districtvariable = "Department",
             pvector=c(2, 3, 2, 3, 2, 3),
             hvector= c(1, 1, 0.05, 0.05, 0.01, 0.01)) 

RDEstimation(datainput=data_hn,outcome="Department_FRANCISCOMORAZAN",
             runningvariable="votemargin_rel",
             districtvariable = "Department",
             pvector=c(2, 3, 2, 3, 2, 3),
             hvector= c(1, 1, 0.05, 0.05, 0.01, 0.01)) 


GenerateRDDTable(datainput=data_hn,
                 outcomevariables=c("Elected_F1","Sex_num","Year","no_seats"),
                 outcomevariables_label=c("Elected in next election","Male","Year","no_seats"),
                 runningvariable="votemargin_rel",
                 districtvariable = "Department",
                 fuzzyvariable=0,
                 pvector=c(2, 3, 2, 3, 2, 3),
                 hvector= c(1, 1, 0.05, 0.05, 0.01, 0.01),
                 latex=0)



  


# (iv) Fuzzy RDD

# (a) Switzerland

RDEstimation(datainput=data,outcome="elected_F1",
             runningvariable="votemargin_rel",
             fuzzy="incumbent_F1",
             pvector=c(2, 3, 2, 3, 2, 3),
             hvector= c(1, 1, 0.05, 0.05, 0.01, 0.01)) 




# rd3 <- rdrobust(data$elected_F1, data$votemargin_rel, p = 2, bwselect = "CV")
# result_file
# 
# data %>%
#   group_by(elected) %>%
#   filter(abs(votemargin_rel) <= 0.11) %>%
#   summarize(n = n()) # get nobs
# data %>%
#   group_by(elected) %>%
#   filter(abs(votemargin_rel) <= 0.15) %>%
#   summarize(n = n()) # get nobs
# data %>%
#   group_by(elected) %>%
#   filter(abs(votemargin_rel) <= 0.05) %>%
#   summarize(n = n()) # get nobs
# data %>%
#   group_by(elected) %>%
#   filter(abs(votemargin_rel) <= 0.01) %>%
#   summarize(n = n()) # get nobs




# (v) Additinal information on data


data %>%
  arrange(votes) %>%
  tail()
data %>%
  filter(abs(votemargin_rel) <= window) %>%
  arrange(abs(votemargin_rel)) %>%
  print(n = Inf)


data_closest <- data %>%
  filter(abs(votemargin_rel) <= 0.01) %>%
  select(name, firstname, year, canton, votes, votemargin_rel, votemargin) %>%
  arrange(abs(votemargin))

print(xtable(data_closest[c(1:20), ], digits = c(0, 0, 0, 0, 0, 0, 4, 0)), include.rownames = FALSE)


data_farthest <- data %>%
  filter(abs(votemargin_rel) <= 0.01) %>%
  select(name, firstname, year, canton, votes, votemargin_rel, votemargin) %>%
  arrange(-abs(votemargin))

print(xtable(data_farthest[c(1:20), ], digits = c(0, 0, 0, 0, 0, 0, 4, 0)), include.rownames = FALSE)

data_farthest <- data %>%
  select(name, firstname, year, canton, votes, votemargin_rel, votemargin) %>%
  arrange(-abs(votemargin))

print(xtable(data_farthest[c(1:20), ], digits = c(0, 0, 0, 0, 0, 0, 4, 0)), include.rownames = FALSE)
