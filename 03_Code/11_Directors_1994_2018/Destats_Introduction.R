library(tidyverse)
library(scales)
library(ggpubr)
library(ggtext)  # Needed for markdown in axis labels


path <- "C:/Schmidlu/Dropbox/Projekt Nationalräte"
path_overleaf <- "C:/Schmidlu/Dropbox/Apps/Overleaf/Political_Rents"

setwd(path)

rescale <- 2.48/0.516

# (A) Presentation ----

df <- tibble(mandates=c(0.05,2.48,NA,0.032,0.516,NA),
             group=c(1:6)) %>%
  mutate(mandates_rescaled = ifelse(group %in% c(4,5),mandates*rescale,mandates))


ggplot(data = df, aes(x = group, y = mandates_rescaled)) + 
  geom_bar(stat = "identity", width = c(0.2, 0.2)) + 
  coord_flip() + 
  theme_classic(base_size = 48) + 
  labs(x = NULL, y = NULL) + 
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_markdown()  # Enable bold text
  ) +
  geom_text(
    aes(label = mandates),
    hjust = -0.2, vjust = 0.4, size = 12) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c(
      "General population",
      "National councillors",
      "**Number of directorships:**",
      "General population",
      "National councillors",
      "**Share with directorship:**"
    )
  ) + 
  scale_y_continuous(limits=c(0,3.5))

ggsave(file = paste0(path_overleaf,"/figures/intro.pdf"), 
       width = 25, height = 10)
