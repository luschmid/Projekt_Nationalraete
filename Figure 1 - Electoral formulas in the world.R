
# Figure 1: Electoral formulas in the world

# 1. Load packages and set directory

library(maps)
library(countrycode)
library(tidyverse)
library(mapdata)
library(RColorBrewer)

setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte")

# 2. Read in and transform IDEA data

rm(list=ls())
df <- read_delim("./01_Raw_data/13_Running_variable/IDEAS data 20191125 answers-export.csv",delim=";") 

df <- df%>%mutate(Value=ifelse(Value %in% c("Mixed","PR","Plurality/Majority and PR"),"PR and mixed",Value),
                  Value=ifelse(Value %in% c("Not applicable","NA","Other","In transition"),"Not applicable, other, and in transition",Value))


# 3. Load map data and merge IDEA to map

world <- map_data("world")
world$iso3c <- countrycode(world$region,origin='country.name',destination = 'iso3c')

df$iso3c <- countrycode(df$Country,origin='country.name',destination = 'iso3c')

world_electoralsystems = right_join(world,df, by = "iso3c")
head(world_electoralsystems)

# 4. Black and white plot for paper

cols=c(brewer.pal(3,"Set1")[c(2,3)],"white")

theme_map <- function(...) {
  theme_void() +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.margin = margin(6, 6, 6, 6)
    )
}


cols <- c("gray30","gray80","white")

ggplot(world_electoralsystems) +
  geom_polygon(aes(x = long, y = lat, fill = Value, group = group), 
               color = "black") +
  coord_equal() +
  scale_fill_manual(values=cols,limits=c("PR and mixed","Plurality/Majority","Not applicable, other, and in transition")) +
  theme_map() +
  theme(legend.position = "bottom") +
  xlab("") + ylab("") + theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=rel(3)),
        legend.spacing.x = unit(0.5,"cm")) 

ggsave(file="./05_Texts_and_presns/01_Running_variable/figures/Ideas_Map_bw.pdf",width=25,height=11 )


# 5. Colored plot for presentation 

cols <- c(rgb(128/255,179/255,255/255),rgb(173/255,222/255,173/255),"white")

ggplot(world_electoralsystems) +
  geom_polygon(aes(x = long, y = lat, fill = Value, group = group), 
               color = "black") +
  coord_equal() +
  scale_fill_manual(values=cols,limits=c("PR and mixed","Plurality/Majority","Not applicable, other, and in transition")) +
  theme_map() +
  theme(legend.position = "bottom") +
  xlab("") + ylab("") + theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=rel(3)),
        legend.spacing.x = unit(0.5,"cm")) 

ggsave(file="./05_Texts_and_presns/01_Running_variable/figures/Ideas_Map_col.pdf",width=25,height=11 )

