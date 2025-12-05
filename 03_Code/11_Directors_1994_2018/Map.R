
# Figure 1: Electoral formulas in the world

# 1. Load packages and set directory

library(maps)
library(countrycode)
library(tidyverse)
library(mapdata)
library(RColorBrewer)

setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte")
path_overleaf <- "C:/Schmidlu/Dropbox/Apps/Overleaf/Political_Rents"


# 2. Read in and transform IDEA data

df <- tibble::tribble(
  ~Country,        ~Value,
  "Argentina",     1,
  "Australia",     0,
  "Austria",       1,
  "Belgium",       0,
  "Brazil",        1,
  "Canada",        0,
  "Chile",         1,
  "Colombia",      1,
  "Czech Rep.",    0,
  "Denmark",       0,
  "Finland",       0,
  "France",        1,
  "Germany",       0,
  "Greece",        1,
  "Hong Kong",     0,
  "Hungary",       0,
  "India",         0,
  "Indonesia",     0,
  "Ireland",       1,
  "Israel",        1,
  "Italy",         0,
  "Japan",         0,
  "Luxembourg",    0,
  "Malaysia",      0,
  "Mexico",        0,
  "Netherlands",   0,
  "New Zealand",   0,
  "Norway",        0,
  "Peru",          1,
  "Philippines",   1,
  "Poland",        0,
  "Portugal",      1,
  "Russia",        1,
  "Singapore",     0,
  "South Africa",  0,
  "South Korea",   1,
  "Spain",         0,
  "Sri Lanka",     1,
  "Sweden",        0,
  "Switzerland",   0,
  "Taiwan",        0,
  "Thailand",      0,
  "Turkey",        1,
  "UK",            0,
  "US",            1,
  "Venezuela",     0,
  "Zimbabwe",      0)

df %>% summarize(mean_value=mean(Value))

# 3. Load map data and merge IDEA to map

world <- map_data("world") %>% 
  filter(region != "Antarctica")  # Remove Antarctica
world$iso3c <- countrycode(world$region,origin='country.name',destination = 'iso3c')

df$iso3c <- countrycode(df$Country,origin='country.name',destination = 'iso3c')


world_regulations = world %>% left_join(df, by = "iso3c") %>%
  mutate(col= case_when(Value==1 ~ "gray30",
                        Value==0 ~ "gray80",
                        is.na(Value) ~"white"))
head(world_regulations)

world_regulations %>% filter(!is.na(Value)) %>% distinct(Country, Value)

# 4. Black and white plot for paper


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


cols <- c("gray30", "gray80", "white" )
ggplot(world_regulations) +
  geom_polygon(aes(x = long, y = lat, fill = col, group = group), 
               color = "black") +
  coord_equal() +
  scale_fill_manual(values=cols,
                    limits=cols,
                    labels=c("Probibited or major restrictions  ",
                             "Allowed, possibly with minor restrictions on government-controlled firms  ",
                             "No data")) +
  theme_map() +
  theme(legend.position = "bottom") +
  xlab("") + ylab("") + theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=rel(3)),
        legend.spacing.x = unit(0.5,"cm")) 

ggsave(file = paste0(path_overleaf,"/figures/map_bw.pdf"), 
       width = 25, height = 10)


# 
# # 5. Colored plot for presentation 
# 
# cols <- c(rgb(128/255,179/255,255/255),rgb(173/255,222/255,173/255),"white")
# 
# ggplot(world_electoralsystems) +
#   geom_polygon(aes(x = long, y = lat, fill = Value, group = group), 
#                color = "black") +
#   coord_equal() +
#   scale_fill_manual(values=cols,limits=c("PR and mixed","Plurality/Majority","Not applicable, other, and in transition")) +
#   theme_map() +
#   theme(legend.position = "bottom") +
#   xlab("") + ylab("") + theme(legend.title = element_blank()) +
#   theme(legend.text=element_text(size=rel(3)),
#         legend.spacing.x = unit(0.5,"cm")) 
# 
# ggsave(file="./05_Texts_and_presns/01_Running_variable/figures/Ideas_Map_col.pdf",width=25,height=11 )
# 
