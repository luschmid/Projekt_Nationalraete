# Comparison highest average vs. largest remainder


#setwd("D:/SchmidLu/Dropbox/Projekt Nationalräte/Running variable")
setwd("C:/Users/StempfeF/Dropbox/Projekt Nationalräte/Running variable")
setwd("/Users/Florence/Dropbox/Projekt Nationalräte/Running variable")

install.packages("tidyverse")
install.packages("readstata13")
install.packages("ggplot2")
install.packages("magrittr")
install.packages("dplyr")
install.packages("tidyr")

library(tidyverse)
library(readstata13)
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyr)

####################################
# (A) Generate file with party votes
####################################

votes <- c(100,80,30,20)
#votes_alt <- c(112500,62500,5000,20000)
#votes <- votes_alt

df <- data.frame(votes=votes,
                 party=LETTERS[1:4])



####################################
# (B) Highest average method
####################################


for (j in 1:8){
  df[,j+2] <- votes/j
}

colnames(df)[3:(j+2)] = 1:j
df_long <-  df %>%
  gather(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,key="seat",value="bruchzahl")


df_long <- df_long %>%
           mutate(rank=rank(bruchzahl)) %>%
           mutate(elected_highest_average=ifelse(33-rank<=8,1,0)) %>%
           arrange(party,seat) 

df_long$seat <- as.numeric(df_long$seat )
df_long$elected_highest_average <- as.numeric(df_long$elected_highest_average )

threshold <- df_long %>%
             filter(rank==25)%>%
             pull(bruchzahl)

ggplot(df_long,aes(x=seat,y=bruchzahl,color=factor(party)))+
  geom_point(aes(shape=factor(elected_highest_average)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  xlab("Seat number") + ylab("Quotient") +
  theme(legend.position="bottom") +
  geom_hline(yintercept=threshold,linetype="dashed")+
  ggtitle("D'Hondt")
  
ggsave(file="./Figures/DHondt.pdf",width=20,height=15 )



# Sainte Lague
df_long$bruchzahlSL=df_long$votes/(2*df_long$seat-1)

df_long <- df_long %>%
  mutate(rank=rank(bruchzahlSL)) %>%
  mutate(elected_SL=ifelse(33-rank<=8,1,0)) %>%
  arrange(party,seat) 

thresholdSL <- df_long %>%
  filter(rank==25)%>%
  pull(bruchzahlSL)


ggplot(df_long,aes(x=seat,y=bruchzahlSL,color=factor(party)))+
  geom_point(aes(shape=factor(elected_SL)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  xlab("Seat number") + ylab("Quotient") +
  theme(legend.position="bottom") +
  geom_hline(yintercept=thresholdSL,linetype="dashed") +
  ggtitle("Sainte-Laguë")

ggsave(file="./Figures/SaintLague.pdf",width=20,height=15 )





####################################
# (C) Largest remainder method
####################################

df_long <- df_long %>%
           left_join(df)%>%
           select(party, seat, bruchzahl,votes,elected_highest_average)

Quotient <- sum(df$votes)/8

df_long <- df_long %>%
           mutate(largest_remainder=votes-(seat-1)*Quotient) %>%
           mutate(rank=rank(largest_remainder)) %>%
           mutate(elected_largest_remainder=ifelse(33-rank<=8,1,0)) %>%
           arrange(party,seat) %>%
           select(-rank) 
           

ggplot(df_long,aes(x=seat,y=largest_remainder,color=factor(party)))+
  geom_point(aes(shape=factor(elected_largest_remainder)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  xlab("Seat number") + ylab("Remaining votes") +
  scale_y_continuous(limits=c(-200,100))+
  geom_hline(yintercept=Quotient,linetype="dashed") +
  theme(legend.position="bottom")+
  ggtitle("Hare")

ggsave(file="./Figures/Hare.pdf",width=20,height=15 )


####################################
# (D) Plot votes to lose/gain a seat
####################################

df_Hare = read.dta13("./ExamplesMS/PartymarginHare.dta")%>%
          select(party, number,ActualVotes, Partymargin, SumVotes)

df_Hare$Partymargin = df_Hare$Partymargin*-1
df_Hare$elected = ifelse(df_Hare$Partymargin > 0, c(1), c(0)) 

ggplot(df_Hare,aes(x=number,y=Partymargin,color=factor(party)))+
  geom_point(aes(shape=factor(elected)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  xlab("Seat number") + ylab("Partymargin") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme(legend.position="bottom")+
  ggtitle("Hare")

ggsave(file="./Figures/partymargin_hare.pdf",width=20,height=15 )


df_DHondt = read.dta13("./ExamplesMS/PartymarginDHondt.dta") %>%
            select(party, number,partymargin_ji, SumVotes)

df_DHondt$elected = ifelse(df_DHondt$partymargin_ji > 0, c(1), c(0)) 

ggplot(df_DHondt,aes(x=number,y=partymargin_ji,color=factor(party)))+
  geom_point(aes(shape=factor(elected)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  xlab("Seat number") + ylab("Partymargin") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme(legend.position="bottom")+
  ggtitle("D'Hondt")

ggsave(file="./Figures/partymargin_dhondt.pdf",width=20,height=15 )


df_SaintLague = read.dta13("./ExamplesMS/PartymarginSaintLague.dta") %>%
  select(party, number,partymargin_ji, SumVotes)

df_SaintLague$elected = ifelse(df_SaintLague$partymargin_ji > 0, c(1), c(0)) 

ggplot(df_SaintLague,aes(x=number,y=partymargin_ji,color=factor(party)))+
  geom_point(aes(shape=factor(elected)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  xlab("Seat number") + ylab("Partymargin") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme(legend.position="bottom") +
  ggtitle("Sainte-Laguë")

ggsave(file="./Figures/partymargin_saintlague.pdf",width=20,height=15 )



# relative to SumVotes

df_Hare = read.dta13("./ExamplesMS/PartymarginHare.dta")%>%
  select(party, number,ActualVotes, Partymargin, SumVotes)

df_Hare$Partymargin = df_Hare$Partymargin*-1
df_Hare$elected = ifelse(df_Hare$Partymargin > 0, c(1), c(0)) 
df_Hare$hypotheticalVoteShare = (-df_Hare$Partymargin+df_Hare$ActualVotes)/df_Hare$SumVotes
df_Hare$ActualVoteShare = df_Hare$ActualVotes/230
df_Hare$PartymarginRelative = df_Hare$ActualVoteShare-df_Hare$hypotheticalVoteShare
df_Hare$PartymarginRelative2 = df_Hare$Partymargin/(230-df_Hare$Partymargin)

 
ggplot(df_Hare,aes(x=number,y=PartymarginRelative,color=factor(party)))+
  geom_point(aes(shape=factor(elected)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  scale_y_continuous(limits=c(-1,0.5))+
  xlab("Seat number") + ylab("Difference vote share") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme(legend.position="bottom") +
  ggtitle("Hare")

ggsave(file="./Figures/partymargin_hare_relative.pdf",width=20,height=15)


# Alternativ
ggplot(df_Hare,aes(x=number,y=PartymarginRelative2,color=factor(party)))+
  geom_point(aes(shape=factor(elected)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  scale_y_continuous(limits=c(-1,1))+
  xlab("Seat number") + ylab("Difference vote share") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme(legend.position="bottom") +
  ggtitle("Hare")

ggsave(file="./Figures/partymargin_hare_relative2.pdf",width=20,height=15)



df_DHondt = read.dta13("./ExamplesMS/PartymarginDHondt.dta") %>%
  select(party, number,partymargin_ji, SumVotes)


df_DHondt$ActualVotes[df_DHondt$party == 1] <- 100
df_DHondt$ActualVotes[df_DHondt$party == 2] <- 80
df_DHondt$ActualVotes[df_DHondt$party == 3] <- 30
df_DHondt$ActualVotes[df_DHondt$party == 4] <- 20


df_DHondt$elected = ifelse(df_DHondt$partymargin_ji > 0, c(1), c(0)) 
df_DHondt$hypotheticalVoteShare = (-df_DHondt$partymargin_ji+df_DHondt$ActualVotes)/df_DHondt$SumVotes
df_DHondt$ActualVoteShare = df_DHondt$ActualVotes/230
df_DHondt$PartymarginRelative = df_DHondt$ActualVoteShare-df_DHondt$hypotheticalVoteShare
df_DHondt$PartymarginRelative2 = df_DHondt$partymargin_ji/(230-df_DHondt$partymargin_ji)

df_DHondt$elected = ifelse(df_DHondt$partymargin_ji > 0, c(1), c(0)) 

ggplot(df_DHondt,aes(x=number,y=PartymarginRelative,color=factor(party)))+
  geom_point(aes(shape=factor(elected)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  scale_y_continuous(limits=c(-1,0.5))+
  xlab("Seat number") + ylab("Difference vote share") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme(legend.position="bottom") +
  ggtitle("D'Hondt")

ggsave(file="./Figures/partymargin_dhondt_relative.pdf",width=20,height=15)


# Alternativ
ggplot(df_DHondt,aes(x=number,y=PartymarginRelative2,color=factor(party)))+
  geom_point(aes(shape=factor(elected)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  scale_y_continuous(limits=c(-1,1))+
  xlab("Seat number") + ylab("Difference vote share") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme(legend.position="bottom") +
  ggtitle("D'Hondt")

ggsave(file="./Figures/partymargin_dhondt_relative2.pdf",width=20,height=15)


df_SaintLague = read.dta13("./ExamplesMS/PartymarginSaintLague.dta")%>%
  select(party, number, partymargin_ji, SumVotes)

df_SaintLague$ActualVotes[df_SaintLague$party == 1] <- 100
df_SaintLague$ActualVotes[df_SaintLague$party == 2] <- 80
df_SaintLague$ActualVotes[df_SaintLague$party == 3] <- 30
df_SaintLague$ActualVotes[df_SaintLague$party == 4] <- 20


df_SaintLague$elected = ifelse(df_SaintLague$partymargin_ji > 0, c(1), c(0)) 
df_SaintLague$hypotheticalVoteShare = (-df_SaintLague$partymargin_ji+df_SaintLague$ActualVotes)/df_SaintLague$SumVotes
df_SaintLague$ActualVoteShare = df_SaintLague$ActualVotes/230
df_SaintLague$PartymarginRelative = df_SaintLague$ActualVoteShare-df_SaintLague$hypotheticalVoteShare
df_SaintLague$PartymarginRelative2 = df_SaintLague$partymargin_ji/(230-df_SaintLague$partymargin_ji)

df_SaintLague$elected = ifelse(df_SaintLague$partymargin_ji > 0, c(1), c(0)) 


ggplot(df_SaintLague,aes(x=number,y=PartymarginRelative,color=factor(party)))+
  geom_point(aes(shape=factor(elected)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  scale_y_continuous(limits=c(-1,0.5))+
  xlab("Seat number") + ylab("Difference vote share") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme(legend.position="bottom") +
  ggtitle("Sainte-Laguë")

ggsave(file="./Figures/partymargin_saintlague_relative.pdf",width=20,height=15)


# Alternativ
ggplot(df_SaintLague,aes(x=number,y=PartymarginRelative2,color=factor(party)))+
  geom_point(aes(shape=factor(elected)),size=4)+
  geom_line()+
  theme_bw(base_size=40) +
  scale_shape_manual("Elected",values=c(1,16))+
  scale_color_brewer("Party",palette="Set1") +
  scale_x_continuous(breaks=c(1:8)) +
  scale_y_continuous(limits=c(-1,1))+
  xlab("Seat number") + ylab("Difference vote share") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme(legend.position="bottom") +
  ggtitle("Sainte-Laguë")

ggsave(file="./Figures/partymargin_saintlague_relative2.pdf",width=20,height=15)


# Close elected

# Df_total with all three methods
df_total = merge(df_Hare,df_DHondt,by=c("party","number")) %>%
           select(party, number, PartymarginRelative.x, PartymarginRelative.y, PartymarginRelative2.x, PartymarginRelative2.y)

names(df_total)[names(df_total)=="PartymarginRelative.x"] <- "Hare"
names(df_total)[names(df_total)=="PartymarginRelative.y"] <- "DHondt"
names(df_total)[names(df_total)=="PartymarginRelative2.x"] <- "Method2Hare"
names(df_total)[names(df_total)=="PartymarginRelative2.y"] <- "Method2DHondt"

df_total = merge(df_total,df_SaintLague,by=c("party","number")) %>%
           select(party, number, Hare, DHondt, Method2Hare, Method2DHondt, PartymarginRelative, PartymarginRelative2)

names(df_total)[names(df_total)=="PartymarginRelative"] <- "SaintLague"
names(df_total)[names(df_total)=="PartymarginRelative2"] <- "Method2SaintLague"


df_total_long1 = gather(df_total, Variable, Method1, c(Hare, DHondt, SaintLague)) %>%
                 select(party, number, Variable, Method1)

df_total_long1$close = abs(df_total_long1$Method1)
df_total_long1 = df_total_long1[order(df_total_long1$close),] 



df_total_long2 = gather(df_total, Variable, Method2, c(Method2Hare, Method2DHondt, Method2SaintLague)) %>%
                 select(party, number, Variable, Method2)
df_total_long2$close = abs(df_total_long2$Method2)
df_total_long2 = df_total_long2[order(df_total_long2$close),] 



df_total_long = merge(df_total_long1,df_total_long2,by=c("party","number"))


####################################
# (E) Distribution of partymargins
####################################





# Standartisation (mean 0, standard deviation 1)
df_total$z_PartymarginHare=scale(df_total$PartymarginHare)
df_total$z_PartymarginDHondt=scale(df_total$PartymarginDHondt)

# Means
meanHare=mean(df_total$PartymarginHare)
meanDHondt=mean(df_total$PartymarginDHondt)

# Long format
names(df_total)[names(df_total)=="PartymarginHare"] <- "Hare"
names(df_total)[names(df_total)=="PartymarginDHondt"] <- "DHondt"
df_total_long <- gather(df_total, Method, Partymargin, c(Hare, DHondt))
df_total_long$elected = ifelse((df_total_long$electedHare==1 & df_total_long$Method=="Hare")| (df_total_long$electedDHondt==1 & df_total_long$Method=="DHondt"), 1, 0)
df_total_long1
df_total_long$SumVotes = ifelse(df_total_long$Method=="Hare", df_total_long$SumVotesHare, df_total_long$SumVotesDHondt)

df_total_long=df_total_long %>%
              select(party, number, ActualVotes, Method, Partymargin, elected, z_Partymargin,SumVotes)


# Histogram z_Partymargin
ggplot(df_total_long, aes(x = z_Partymargin)) +
  geom_histogram(aes(fill=Method), binwidth = 0.1, position="dodge") +
  theme_bw(base_size=40) +
  theme(legend.position="bottom")

ggsave(file="./Figures/z_Partymargin.pdf",width=20,height=15 )


# Histogram Partymargin
ggplot(df_total_long, aes(x = Partymargin)) +
  geom_histogram(aes(fill=Method), binwidth = 100, position="dodge") +
  scale_y_continuous(breaks=c(seq(0,10,2))) +
  theme_bw(base_size=40) +
  theme(legend.position="bottom")+
  geom_vline(xintercept=meanHare, color="turquoise3", size=1.5)+
  geom_vline(xintercept=meanDHondt, color="coral", size=1.5)+


ggsave(file="./Figures/Partymargin.pdf",width=20,height=15 )

