
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

try(setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte"))
try(setwd("C:/Dropbox/Projekt Nationalräte")) 

#rm(list = ls())
source("./03_Code/fine_grid.R") # ggplot layers



###################
# (A) Read-in data
###################


data <- read.dta13("./02_Processed_data/15_Elections_Honduras/elections_hn_final.dta")

##############################
# (B) RD estimation
##############################

# (i) Histogram of running variable

ggplot(data, aes(x = votemargin_rel)) +
  geom_histogram(aes(y = ..count.. / sum(..count..)),
                 alpha=0.5,position = "stack",boundary=0,color="black",
                 binwidth=0.002) +
  theme_bw(base_size = 42) +
  scale_x_continuous(limits = c(-0.07, 0.07), 
                     breaks = seq(-0.06, 0.06, 0.02), 
                     label = percent) +
  geom_vline(xintercept = 0) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("") + xlab("")

ggsave(file = "./05_Texts_and_presns/01_Running_variable/figures/votemargin_rel_bw_hn_0_07.pdf", width = 20, height = 12)

window <- 0.01

ggplot(data, aes(x = votemargin_rel)) +
  geom_histogram(aes(y = ..count.. / sum(..count..)), fill = "grey", colour = "black", breaks = seq(-0.01, 0.01, 0.001), binwidth = window / 10) +
  theme_bw_finegrid(base_size = 42) +
  scale_x_continuous(limits = c(-window, window), label = percent) +
  geom_vline(xintercept = 0) +
  ylab("Density") + xlab("\n  Relative vote difference")

ggsave(file = paste("./05_Texts_and_presns/01_Running_variable/figures/votemargin_rel_bw_hn_", window, ".pdf", sep = ""), width = 20, height = 12)




ggplot(data, aes(y = Elected_F1, x = votemargin_rel, group = Elected)) +
  geom_point(size = 4) +
  # geom_text_repel(aes(label=nobs))+
  geom_smooth(data = data, aes(y = Elected_F1, x = votemargin_rel, group = Elected), method = "lm", formula = y ~ poly(x, 2, raw = TRUE)) +
  theme_bw_finegrid(base_size = 52) +
  scale_x_continuous(limits = c(-window, window), label = percent) +
  geom_vline(xintercept = 0) +
  ylab(ylab) + xlab("\n  Relative vote difference")
ggsave(file = paste("./05_Texts_and_presns/01_Running_variable/figures/", y, ".pdf", sep = ""), width = 30, height = 20)



rdbwselect(data$Elected_F1, data$votemargin_rel, p = 2)


pvector <- c(2, 3, 2, 3, 2, 3)
hvector <- c(1, 1, 0.05, 0.05, 0.01, 0.01)
result_file <- as.data.frame(matrix(NA, 3, length(pvector)))

for (i in 1:length(pvector)) {
  if (hvector[i] == 1) {
    rd_model <- rdrobust(data$Elected_F1, data$votemargin_rel, p = pvector[i])
  }
  else {
    rd_model <- rdrobust(data$Elected_F1, data$votemargin_rel, p = pvector[i], h = hvector[i])
  }
  summary(rd_model)
  result_file[, i] <- c(sprintf("%.3f", rd_model$coef[1]), sprintf("(%.3f)", rd_model$se[1]), sprintf("%.2f", rd_model$bws[1, 1]))
}

rownames(result_file) <- c("Estimate", "", "Bw")
print(xtable(result_file))

rd2 <- rdrobust(data$Elected_F1, data$votemargin_rel, p = 2, bwselect = "IK")
rd3 <- rdrobust(data$Elected_F1, data$votemargin_rel, p = 2, bwselect = "CV")