
library(scales)
library(tidyverse)
library(xlsx)
library(openintro)
library(xtable)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(readstata13)
library(rdrobust)
library(ggrepel)
library(caret)
library(rdrobust)

try(setwd("D:/SchmidLu/Dropbox/Projekt Nationalräte")) # ggplot layers

rm(list=ls(all=TRUE))
source("D:/SchmidLu/Dropbox/Various/Software/fine_grid.R") # ggplot layers


#################
# (A) Merge data
#################

nr_1931_2015 <- read.dta13("1 Data/1_data/nationalraete_1931_2015.dta")
rv_analytical_all <- read.dta13("Running variable/Data/RV_analytical_all.dta") %>%
                     select(ID,year,votemargin)

data <- nr_1931_2015 %>% inner_join(rv_analytical_all,by=c("ID","year")) %>% 
            select(firstname, name, ID,year,canton, elected, votes, birthyear, leftright,  votemargin, sex,cand_before1931)%>% 
            as_tibble() %>%
            arrange(ID,year) %>%
            as_tibble()

data_minyear <-    data %>%
                       group_by(ID)%>%
                       summarize(year_min=min(year))

data_canton_votes <-    data %>%
  group_by(canton,year)%>%
  summarize(votes_sum=sum(votes),
            no_seats=sum(elected))

data <- data %>% inner_join(data_minyear,by=c("ID")) %>% 
        inner_join(data_canton_votes,by=c("canton","year")) %>% 
        mutate(votes_rel=votes/(votes_sum/no_seats), 
        votemargin_rel=votemargin/(votes_sum/no_seats),
        age=year-birthyear)  %>%
        mutate(sex=as.numeric(as.factor(data$sex))-1)
 


data_forward <- data  %>%
             group_by(ID) %>% mutate(votes_f1=lead(votes),
                                     votes_rel_f1=lead(votes_rel),
                                     elected_f1=lead(elected))  %>% 
             select(ID,year,ends_with("f1")) %>% 
             mutate(participation_f1=ifelse(elected_f1%in%c(0,1),1,0), 
                    elected_f1=ifelse(is.na(elected_f1),0,elected_f1))

            
data_lag <- data %>% inner_join(data_canton_votes,by=c("canton","year")) %>%
                      group_by(ID) %>% 
                      mutate(votes_l1=lag(votes),
                          votes_rel_l1=lag(votes_rel),
                          elected_l1=lag(elected))  %>% 
                      select(ID,year,ends_with("l1"))


data_working<- data %>% inner_join(data_forward,by=c("ID","year")) %>%
                        inner_join(data_lag,by=c("ID","year")) %>%
                        select(name, firstname, ID,year,year_min,canton, leftright,sex,age, votemargin, votemargin_rel,participation_f1,starts_with("votes"),starts_with("elected"),cand_before1931)%>%
                        filter(elected==0 & votemargin_rel<0|elected==1 & votemargin_rel>0) # temporary: remove four cases FR-1931-9003 and FR-1935-0013 in 1939 as well as ZH-1943-0030 and in 1967

                        
data_working_ct<-data_working%>%select(canton)

dummy <- dummyVars(~ ., data =data_working_ct ) # make dummy for canton
data_working <- cbind(data_working,predict(dummy, data_working_ct))

rm(data_minyear,data_canton_votes,data_forward,data_lag,data,nr_1931_2015,rv_analytical_all) 

write_csv(data_working,"Running variable/Data/data_working.csv",na=".")

data_working %>% distinct(ID) %>% summarize(n=n()) # get nobs


##############################
# (B) RD estimation
##############################

# (i) Histogram of running variable


  ggplot(data_working, aes(x = votemargin_rel)) +
    geom_histogram(aes(y=..count../sum(..count..)),fill="grey",colour="black",breaks= seq(-0.05,0.05,0.005))+
    theme_bw_finegrid(base_size=42) +
    scale_x_continuous(limits = c(-0.05,0.05),label=percent) +
    geom_vline(xintercept=0) +
    ylab("Density") + xlab("\n  Relative vote difference")

ggsave(file=paste("./1 Data/2_figures/votemargin_rel_bw_",0.05,".pdf",sep=""),width=20,height=12 )



ggplot(data_working, aes(x = votemargin_rel)) +
  geom_histogram(aes(y=..count../sum(..count..)),fill="grey",colour="black",breaks= seq(-0.01,0.01,0.001),binwidth =window/10)+
  theme_bw_finegrid(base_size=42) +
  scale_x_continuous(limits = c(-window,window),label=percent) +
  geom_vline(xintercept=0) +
  ylab("Density") + xlab("\n  Relative vote difference")

ggsave(file=paste("./1 Data/2_figures/votemargin_rel_bw_",window,".pdf",sep=""),width=20,height=12 )



# (ii) RDD graphs for balance tests and outcome variables

CutData <- function(data,y,bin_digits=2,year_start=1931,year_end=2015,only_startyear=0,only_participation=0){

    if (only_startyear==1){data <- data %>%  filter(year==year_min & cand_before1931==0)}
    if (only_participation==1){data <- data %>%  filter(participation_f1==1)}
  
    data_estimation <- data %>% 
    mutate(votemargin_rel_rounded=floor(votemargin_rel*10^bin_digits)/10^bin_digits+(5/(10^(bin_digits+1))))%>% 
    group_by(votemargin_rel_rounded)%>%
    filter(year>=year_start & year<=year_end)%>%
    summarize(age_mean=mean(age),sex_mean=mean(sex),
              left_right_mean=mean(leftright,na.rm=T),
              participation_f1_mean=mean(participation_f1,na.rm=T),
              elected_f1_mean=mean(elected_f1,na.rm=T),
              year_mean=mean(year,na.rm=T),AG_mean=mean(cantonAG,na.rm=T),
              cantonAI_mean=mean(cantonAI,na.rm=T),
              cantonAR_mean=mean(cantonAR,na.rm=T),
              cantonBE_mean=mean(cantonBE,na.rm=T),
              cantonJU_mean=mean(cantonJU,na.rm=T),
              cantonBL_mean=mean(cantonBL,na.rm=T),
              cantonBS_mean=mean(cantonBS,na.rm=T),
              cantonFR_mean=mean(cantonFR,na.rm=T),
              cantonGE_mean=mean(cantonGE,na.rm=T),
              cantonGL_mean=mean(cantonGL,na.rm=T),
              cantonGR_mean=mean(cantonGR,na.rm=T),
              cantonLU_mean=mean(cantonLU,na.rm=T),
              cantonNE_mean=mean(cantonNE,na.rm=T),
              cantonNW_mean=mean(cantonNW,na.rm=T),
              cantonOW_mean=mean(cantonOW,na.rm=T),
              cantonSG_mean=mean(cantonSG,na.rm=T),
              cantonSH_mean=mean(cantonSH,na.rm=T),
              cantonSO_mean=mean(cantonSO,na.rm=T),
              cantonSZ_mean=mean(cantonSZ,na.rm=T),
              cantonTG_mean=mean(cantonTG,na.rm=T),
              cantonTI_mean=mean(cantonTI,na.rm=T),
              cantonUR_mean=mean(cantonUR,na.rm=T),
              cantonVD_mean=mean(cantonVD,na.rm=T),
              cantonVS_mean=mean(cantonVS,na.rm=T),
              cantonZG_mean=mean(cantonZG,na.rm=T),
              cantonZH_mean=mean(cantonZH,na.rm=T),
              nobs=n()) %>%
    mutate(elected=ifelse(votemargin_rel_rounded>=0,1,0))
    return(data_estimation)
}



PlotDiscontinuity = function(data,y,y_mean,ylab,group="elected",window=0.05,polynomials=2,only_participation=0) {
  
  if(only_participation==1){data_working_plot <- data_working  %>%  filter(participation_f1==1)}
  else{data_working_plot <- data_working}
  
  ggplot(data , aes_string(y=y_mean,x="votemargin_rel_rounded",group=group))+
    geom_point(size=4)+
    #geom_text_repel(aes(label=nobs))+
    geom_smooth(data=data_working_plot ,aes_string(y=y ,x="votemargin_rel",group=group),method="lm",formula=y ~ poly(x, 2, raw=TRUE))+ 
    theme_bw_finegrid(base_size=52) +
    scale_x_continuous(limits = c(-window,window),label=percent) +
    geom_vline(xintercept=0) +
    ylab(ylab) + xlab("\n  Relative vote difference") 
    ggsave(file=paste("./1 Data/2_figures/",y,".pdf",sep=""),width=30,height=20)
}


data_estimation <- CutData(data_working,y,bin_digits=3,year_start=1931,year_end=2015,only_startyear=0,only_participation=0)
varname=as.character(list("elected_f1","participation_f1","age","sex","cantonZH","cantonTI","year"))
varname2 <- map(varname, ~ paste(.x,"_mean",sep="")) 
varlabel <- as.character(list("Elected in t+1","Participation in t+1","Age","Male","Canton of Zurich","Canton of Ticino","Year"))


for (i in 1:length(varname)){
PlotDiscontinuity(data_estimation,y_mean=varname2[[i]],y=varname[i],ylab=varlabel[i],group="elected",window=0.01, polynomials = 2)
}

data_estimation <- CutData(data_working,y,bin_digits=3,year_start=1931,year_end=2015,only_startyear=0,only_participation=1)
PlotDiscontinuity(data_estimation,y_mean="elected_f1_mean",y="elected_f1",ylab="Participation in t+1",group="elected",window=0.01, polynomials = 2,only_participation=1)





# (iii) RD estimation 

rdbwselect(data_working$elected_f1,data_working$votemargin_rel,p=2)


rd1 <- rdrobust(data_working$elected_f1,data_working$votemargin_rel,p=2)
rd2 <- rdrobust(data_working$elected_f1,data_working$votemargin_rel,p=3)
rd3 <- rdrobust(data_working$elected_f1,data_working$votemargin_rel,p=2,h=0.05)
rd4 <- rdrobust(data_working$elected_f1,data_working$votemargin_rel,p=3,h=0.05)
rd5 <- rdrobust(data_working$elected_f1,data_working$votemargin_rel,p=2,h=0.01)
rd6 <- rdrobust(data_working$elected_f1,data_working$votemargin_rel,p=3,h=0.01)


summary(rd3)
rd6$bws

for (i in 1:6){
summary(paste("rd",i,sep=""))
}

rd2 <- rdrobust(data_working$elected_f1,data_working$votemargin_rel,p=2,bwselect="IK")
rd3 <- rdrobust(data_working$elected_f1,data_working$votemargin_rel,p=2,bwselect="CV")


data_working %>% group_by(elected) %>% filter(abs(votemargin_rel)<=0.11) %>% summarize(n=n()) # get nobs
data_working %>% group_by(elected) %>% filter(abs(votemargin_rel)<=0.15) %>% summarize(n=n()) # get nobs
data_working %>% group_by(elected) %>% filter(abs(votemargin_rel)<=0.05) %>% summarize(n=n()) # get nobs
data_working %>% group_by(elected) %>% filter(abs(votemargin_rel)<=0.01) %>% summarize(n=n()) # get nobs


# (iv) Additinal information on data

data_working %>% arrange(votes)%>% tail()
data_working %>% filter( abs(votemargin_rel)<=window) %>% arrange(abs(votemargin_rel)) %>%print(n=Inf)


data_closest <- data_working %>% filter(abs(votemargin_rel)<=0.01 ) %>% 
                 select(name,firstname,year, canton, votes, votemargin_rel,votemargin) %>% 
                 arrange(abs(votemargin))
                 
print(xtable(data_closest[c(1:20),],digits=c(0,0,0,0,0,0,4,0)),include.rownames=FALSE)


data_farthest <- data_working %>% filter(abs(votemargin_rel)<=0.01 ) %>% 
  select(name,firstname,year, canton, votes, votemargin_rel,votemargin) %>% 
  arrange(-abs(votemargin))

print(xtable(data_farthest[c(1:20),],digits=c(0,0,0,0,0,0,4,0)),include.rownames=FALSE)

data_farthest <- data_working %>% 
  select(name,firstname,year, canton, votes, votemargin_rel,votemargin,party) %>% 
  arrange(-abs(votemargin))

print(xtable(data_farthest[c(1:20),],digits=c(0,0,0,0,0,0,4,0)),include.rownames=FALSE)


### OLD CODE ###################
# 
# 
# for (i in 1:length(covariates)){
#   data_estimation <- data_working  %>% 
#     filter(year==year_min ) %>% 
#     select(paste(covariates[i],sep=""),votemargin_rel)  %>%
#     arrange(votemargin_rel) 
#   colnames(data_estimation) <- c("y","x") 
#   
#   rd <- rdrobust(data_estimation$y,data_estimation$x,p=2)
#   summary(rd)
#   attributes(rd)
#   balance_results[1,i] <- rd$coef[1]
#   balance_results[2,i] <-rd$se[1]
# }
# 
# print(xtable(balance_results,digits=2),include.rownames=FALSE)
# 
# rdplot(data_estimation$y,data_estimation$x,p=2,x.label="Relative vote distance",nbins=20,cex=2,cex.axis=2,cex.lab=2,cex.main=2,ci=T)
# 
# 
#   data_working  %>% 
#   filter(year==year_min ) %>% 
#   group_by(elected)%>% 
#   #filter(elected==0 & votemargin_rel>-window| elected==1 & votemargin_rel<window)%>% 
#   summarize(age_mean=mean(age),sex_mean=mean(sex),left_right_mean=mean(leftright,na.rm=T),nobs=n()) 
#   
# 
# 
