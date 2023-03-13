# Figure A.2: Arbitrariness of simulated election probabilities

# 1. Load packages and set directory

suppressMessages(library(tidyverse))
suppressMessages(library(readstata13))
suppressMessages(library(data.table))

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))

# 2. Example data from paper

df <- data.frame(votes_h=c(22,8,5,25,15,5,8,7,5),
                 votes_j=c(rep(35,3),rep(45,3),rep(20,3)),
                 elected=c(1,0,0,1,1,0,0,0,0),
                 party_num=c(rep(1,3),rep(2,3),rep(3,3)),
                 seats=rep(3,9),
                 ID=as.character(1:9))

df$votes_h_share <- df$votes_h/sum(df$votes_h)
data=df

# 3. Run simulation
# Note: done on unifr server (Feb 2021)
# set.seed(05071975) 
# Functions are here: 00_Kotakorpi_Example.R
# KotakorpiEtAl2017Loop(m_min=1,m_max=200000,data=df,M=20000,control_seat_allocation=0,check_convergence=0)

# 4. Read in data from simulation

df_sim <- fread(file=paste0("./02_Processed_data/13_Running_variable/03_kotakorpietal2017_example/kotakorpietal2017_example.csv")) %>%
          mutate(ID=as.character(ID),
                 elected_sum_share=elected_sum/20000)

df_all <- df %>% left_join(df_sim,by=c("ID","party_num")) %>%
  select(votes_h,votes_j,elected,party_num,ID,votes_h_share,elected_sum_share,m) %>%
  as_tibble() %>%
  mutate(candidate_description=ifelse(ID%in%c(1,4,7),"Candidate 1",ifelse(ID%in%c(2,5,8),"Candidate 2","Candidate 3")))

# 5. Figure

cols <- RColorBrewer::brewer.pal(3, "Set1")

ggplot(df_all,aes(x=m,y=elected_sum_share,group=ID,color=as.factor(party_num),
                  linetype=as.factor(candidate_description))) +
  geom_line(size=1.5) +
  theme_bw(base_size=32) +
  xlab("") + ylab("") + 
  scale_color_manual(labels=paste0("Party " ,c(1:3)),values=cols)+
  scale_linetype_manual(values=c("solid","dashed","dotted")) +
  theme(legend.position="bottom",legend.title = element_blank()) +
  scale_x_continuous(trans='log10',limits=c(0.22,10000)) +
  annotate("text", x=0.78,y=0.35, label="0.35", size=8) +
  annotate("text", x=0.78,y=0.45, label="0.45", size=8) +
  annotate("text", x=0.78,y=0.2, label="0.20", size=8) +
  geom_segment(aes(x=0.98, y=0.2, xend=1.05, yend=0.2),  size=2,color="black") +
  geom_segment(aes(x=0.98, y=0.35, xend=1.05, yend=0.35),  size=2,color="black") +
  geom_segment(aes(x=0.98, y=0.45, xend=1.05, yend=0.45),  size=2,color="black") +
  guides(shape = guide_legend(override.aes = list(size = 3)),
         linetype = guide_legend(override.aes = list(shape = NA)),
         size = FALSE) +
  theme(legend.position = "none")

ggsave("./04_Results/01_Running_variable/KotakorpiEtAl_Example.pdf",width=20,height=12)

