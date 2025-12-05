library(arm)
library(foreign)
library(plotrix)
library(ggplot2)
library(graphics)
library(xtable)
library(RColorBrewer)
library(grid)
library(scales)

rm(list=ls())




  setwd("C:/Dropbox/Projekt Nationalräte/1 Data")
  source("./5_codes/10_fine_grid.R") # ggplot layers
  source("./5_codes/10_multi_plot.R") # multiplot (requires grid library)


  #########################################################
  # (A) Number of candidates per canton by elected status
  #########################################################


  data<-read.dta("./1_data/NRW_KANDIDATEN.dta")

  data$gewaehlt[data$gewaehlt=="G"]<-"Elected   "
  data$gewaehlt[data$gewaehlt=="N"]<-"Not elected"

  data$knr_numeric <- NA
  data$knr_numeric[!is.na(data$kantonsname)]<-   transform[as.numeric(data$kantonsname[!is.na(data$kantonsname)])]   
  

  kantonsnamen <- levels(factor(data$kantonsname))

  data$knr_numeric <- array(match(data$kantonsname, kantonsnamen), dim=length(data$kantonsname), dimnames=length(data$kantonsname))

  data$kantonsnummer_mod <- 27-data$knr_numeric


 hist_cut <- ggplot(data, aes(x=factor(kantonsnummer_mod), fill=factor(gewaehlt)))

(p1 <-  hist_cut + geom_bar() +
  theme_bw_finegrid(base_size = 56) +
  xlab("") + ylab("\n Number of candidates") +
  scale_y_continuous(labels = comma) +
  scale_x_discrete(breaks=c(1:26),labels = kantonsnamen[26:1]) +
  coord_flip() +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("grey20","grey70")) +
  theme(legend.direction = "horizontal", legend.position = "bottom")  )
  
  
  ggsave(file="./2_figures/histogram.ps",width=20,height=20 )


  #########################################################
  # (B) Number of candidates per year
  #########################################################


  total_candidates=matrix(NA,15,1)
  year_range=matrix(NA,15,1)
  index=matrix(NA,15,1)

  data$year <- as.numeric(substr(data$nationalratswahl, 4, nchar(data$nationalratswahl)))

  for (i in 1:length(total_candidates)){
    total_candidates[i]=length(data$gewaehlt[is.element(data$year,c(1971+(i-1)*4):(1971+(i)*4)-1)])
    year_range[i]=(1971+(i-1)*4)
    index[i]=i
  }




  df=data.frame(total_candidates=total_candidates,year_range=as.numeric(year_range),index=index)

  df <- df[df$year_range<2015,]

#df=df[df$index<15] # drop last period because we have less data on this time

  
  ts_cut <- ggplot(df, aes(x=index,y=total_candidates))
(p2 <-ts_cut + geom_point(size=5.7) +geom_line(size=1.3) +
  theme_bw_finegrid(base_size = 58) +
  xlab("") + ylab("Number of candidates \n") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = index,label=year_range) +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("grey20","grey70")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) ) 

ggsave(file="./2_figures/time_series.jpg",width=20,height=15 )


#########################################################
# (C) Candidates across space and time together
#########################################################

jpeg(filename = "./2_figures/candidates_all.jpg", width = 1200, height = 880)
multiplot(p1, p2, cols=2)
dev.off()


ggsave(file="./2_figures/candidates_all.jpg",width=10,height=20 )

#########################################################
# (D) Destats
#########################################################

df = data.frame(yscale=c(1:9), 
                value=c(0.2945,(1-0.2945),0.08,(1-0.08),0.1179,0.1131,0.1288,0.1147,(1- sum(0.1179,0.1131,0.1288,0.1147))), 
                group=(c(3,3,2,2,1,1,1,1,1)),
                cat=c(rep(c(1,2),2),c(1:5)),
                label=c("Female","Male","Yes", "No","SP","CVP","FDP","SVP","Other parties"))

df$percent <- paste("(",round(df$value*100,0),"%)",sep="")


df$xpos <-matrix(NA,length(df$yscale),1) 
df$xpos[df$cat==1] <- df$value[df$cat==1]/2
df$xpos[df$cat==2] <- df$value[df$cat==1]+df$value[df$cat==2]/2

df$xpos[df$cat==3 & df$group==1] <- sum(df$value[is.element(df$cat,c(1:2))& df$group==1])+df$value[df$cat==3& df$group==1]/2
df$xpos[df$cat==4 & df$group==1] <- sum(df$value[is.element(df$cat,c(1:3))& df$group==1])+df$value[df$cat==4& df$group==1]/2
df$xpos[df$cat==5 & df$group==1] <- sum(df$value[is.element(df$cat,c(1:4))& df$group==1])+df$value[df$cat==5& df$group==1]/2

col.set <- c(brewer.pal(8,"Pastel1")[c(1:5)])


ggplot(df, aes(x = group, y = value,group=factor(group),fill=factor(cat))) +
  geom_bar(stat='identity',color="black") + 
  geom_text(aes(x= group+0.1, y = xpos,label=label),size=10.1) +
  geom_text(aes(x= group-0.1, y = xpos,label=percent),size=7.1) +
  coord_flip() +
  theme_bw_finegrid(base_size = 36)  + 
  xlab("") + ylab("\n Share") +
  scale_y_continuous(limits = c(0, 1),labels=percent) +
  scale_x_continuous(breaks=c(1:3),labels=c("Party","Elected","Gender")) +
  scale_fill_manual(values=col.set ,guide=F) +
  theme(legend.direction = "horizontal", legend.position = "bottom") +
  theme(panel.grid.minor.y =  element_line(colour = "white", size = 0.5,linetype="dotted")) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="black", fill="grey90" ) )+ 
  theme(panel.margin = unit(0.5, "lines"))

ggsave(file="./2_figures/destat_graphs.jpg",width=20,height=15 )


