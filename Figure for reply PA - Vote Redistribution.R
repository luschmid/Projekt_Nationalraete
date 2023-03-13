rm(list = ls())

setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte")
source("./03_Code/13_Running_variable/00_Vote_margin_analytical_highestaverage.R")
library(tidyverse) # for the pipe
library(RColorBrewer)
library(latex2exp)


RedistributeVotes <- function(data,
                              votes_changed,
                              party_of_interest,
                              parties_to_steal,
                              type="add"){
  # This function redistributes votes to our party of interest from a set of other
  # parties (parties_to_steal) if type is "add". If type is "subtract", 
  # the votes are redistributed from our party of interest to the target parties. 
  
  if (type=="add"){
    for (j in 1:length(parties_to_steal)){
      data$votes_j[data$party==parties_to_steal[j]] <- data$votes_j[data$party==parties_to_steal[j]]-as.numeric(votes_changed[j])
    }  
    data$votes_j[data$party==party_of_interest] <- data$votes_j[data$party==party_of_interest] + sum(votes_changed) 
  }
  
  if (type=="subtract"){
    for (j in 1:length(parties_to_steal)){
      data$votes_j[data$party==parties_to_steal[j]] <- data$votes_j[data$party==parties_to_steal[j]]+as.numeric(votes_changed[j])
    }  
    data$votes_j[data$party==party_of_interest] <- data$votes_j[data$party==party_of_interest] - sum(votes_changed) 
  }
  return(data)
}




df_example <- data.frame(party=c("A","B","C"),
                         votes_j=c(45,35,20),
                         district="example",
                         year=2022)

n_seats <- 3
method <- "dHondt"


dp_example <- PrepareData(data=df_example,
                               votes_j_name="votes_j",
                               party_name="party", 
                               districtname="district", 
                               election_cyclename = "year",
                               alliances=F,
                               system="closed") 

out_all <- out <- CalculateRatios(votes=dp_example %>% select(party,votes_j), 
                                  n=n_seats,
                                  method=method) %>%
  arrange(-No)%>% 
  filter((party==1 & seat==3) | (party==2 & seat==1))%>% 
  mutate(votes_redistr=0)

data_mod <- dp_example

i=1
while (i <=20){
print(paste0("iteration:",i))
data_mod <- RedistributeVotes(data=data_mod,
                              votes_changed=1,
                              party_of_interest=1,
                              parties_to_steal=2,
                              type="add") 
out <- CalculateRatios(votes=data_mod %>% select(party,votes_j), 
               n=n_seats,
               method=method) %>%
  arrange(-No)
print(out)
out_all <- out_all %>% bind_rows(out %>% 
                                   filter((party==1 & seat==3) | (party==2 & seat==1))%>% 
                                 mutate(votes_redistr=i))
i <- i +1 
}

cols <- c(brewer.pal(3,"Set1")[2],brewer.pal(3,"Set1")[1],brewer.pal(3,"Set1")[3])
# Note: same color as in Figure 5 of the paper: P1 is blue, P2 is red, P3 is green


ggplot(data=out_all,aes(x=votes_redistr,
                        y=No,
                        group=factor(party),
                        color=factor(party))) +
         geom_line(linewidth=2) +
  geom_hline(yintercept=21,color=cols[3],linewidth=2) +
  theme_classic(base_size=32) +
  scale_y_continuous(breaks=seq(10,40,5),limits=c(10,40)) +
  scale_color_manual(values=cols) +
  xlab(TeX("Number of votes redistributed from $P_2$ to $P_1$")) + ylab("D'Hondt ratio") +
  scale_x_continuous(breaks=c(0:30)) +
  annotate("text",x=4,y=35,
           label=TeX("D'Hondt ratio seat 1 of $P_2$"), 
           parse=T,
           size=10,
           color=cols[2])+
  annotate("text",x=4,y=18,
           label=TeX("D'Hondt ratio seat 3 of $P_1$"), 
           parse=T,
           size=10,
           color=cols[1])+
  annotate("text",x=4,y=22,
           label=TeX("D'Hondt ratio seat 1 of $P_3$"), 
           parse=T,
           size=10,
           color=cols[3])+ 
  theme(legend.position = "none") +
  geom_segment(aes(x = 15, y = 10, xend = 15, yend = 20),color="black",
               linetype="dashed",
               size=2) +
  geom_segment(aes(x = 16, y = 25, xend = 15, yend = 20),
               arrow = arrow(length = unit(0.9, "cm")),
               color="black",
               size=1.5)+
  annotate("text",x=13,y=28,
           label=TeX("$P_1$ has the same D'Hondt ratio "), 
           parse=T,
           size=10,
           color="black",
           hjust = 0)+
  annotate("text",x=13,y=26.5,
           label=TeX("as $P_2$ but $P_3$ obtains seat"), 
           parse=T,
           size=10,
           color="black",
           hjust = 0) + 
  coord_cartesian(expand = FALSE)

ggsave("./04_Results/01_Running_variable/Redistribute_Votes_Example.pdf",width=20,height=12)


data_mod <- dp_example
method <- "dHondt"
out_all <-   out <- CalculateRatios(votes=data_mod %>% select(party,votes_j), 
                                    n=n_seats,
                                    method=method)%>% 
  filter((party==3 & seat==1) | (party==2 & seat==1))%>% 
  mutate(votes_redistr=0)
#print(out)

i=1
while (i <=8){
  print(paste0("iteration:",i))
  data_mod <- RedistributeVotes(data=data_mod,
                                votes_changed=1,
                                party_of_interest=2,
                                parties_to_steal=3,
                                type="subtract") 
  out <- CalculateRatios(votes=data_mod %>% select(party,votes_j), 
                         n=n_seats,
                         method=method) %>%
    arrange(-No)
  #print(out)
  out_all <- out_all %>% bind_rows(out %>% 
                                     filter((party==3 & seat==1) | (party==2 & seat==1))%>% 
                                     mutate(votes_redistr=i))
  i <- i +1
}

cols <- c(brewer.pal(3,"Set1")[2],brewer.pal(3,"Set1")[1],brewer.pal(3,"Set1")[3])
# Note: same color as in Figure 5 of the paper: P1 is blue, P2 is red, P3 is green


ggplot(data=out_all,aes(x=votes_redistr,
                        y=No,
                        group=factor(party),
                        color=factor(party))) +
  geom_line(linewidth=2) +
  geom_hline(yintercept=22.5,color=cols[1],linewidth=2) +
  theme_classic(base_size=32) +
  scale_color_manual(values=cols[2:3]) +
  scale_x_continuous(breaks=c(0:9),limits=c(0,9)) +
  scale_y_continuous(breaks=seq(0,35,5)) +
  xlab(TeX("Number of votes redistributed from $P_2$ to $P_3$")) + ylab("D'Hondt ratio") +
  annotate("text",x=3.8,y=29,
           label=TeX("D'Hondt ratio seat 1 of $P_2$"), 
           parse=T,
           size=10,
           color=cols[2])+
  annotate("text",x=3.8,y=26,
           label=TeX("D'Hondt ratio seat 1 of $P_3$"), 
           parse=T,
           size=10,
           color=cols[3])+
  annotate("text",x=1.3,y=24,
           label=TeX("D'Hondt ratio seat 2 of $P_1$"), 
           parse=T,
           size=10,
           color=cols[1])+ 
  theme(legend.position = "none") +
  geom_segment(aes(x = 7.5, y = 0, xend = 7.5, yend = 27.5),color="black",
               linetype="dashed",
               size=2) +
  geom_segment(aes(x = 6, y = 31, xend = 7.5, yend = 27.5),
               arrow = arrow(length = unit(0.9, "cm")),
               color="black",
               size=1.5)+
  annotate("text",x=5,y=34,
           label="Solution according to equation (5)", 
           size=10,
           color="black",
           hjust = 0)+
  annotate("text",x=5,y=32,
           label="in Grofman and Selb (2009)", 
           size=10,
           color="black",
           hjust = 0) +
  geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 22.5),color="black",
               linetype="dashed",
               size=2) +
  geom_segment(aes(x = 3, y = 18, xend = 2.5, yend = 22.5),
               arrow = arrow(length = unit(0.9, "cm")),
               color="black",
               size=1.5)+
  annotate("text",x=4,y=16,
           label=TeX("$P_1$ loses second seat to $P_3$"), 
           parse=T,
           size=10)+ 
  coord_cartesian(expand = FALSE)
ggsave("./04_Results/01_Running_variable/Redistribute_Votes_Example_Selb.pdf",width=20,height=12)

