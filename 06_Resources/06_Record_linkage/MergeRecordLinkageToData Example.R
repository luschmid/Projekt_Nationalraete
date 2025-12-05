

library(ggplot2)
library(dplyr)
library(tidyverse) ## dplyr provides the join functions
library(AER) 
library(xtable) 
library(gcookbook) # For the data set
library(RColorBrewer) # For the data set
library(tidyverse)

rm(list = ls())


# (0) Mock connection file and mock original name file

data_out_all <- tibble(UniqueID1=c(1,1,2,2,4,4,4),
                  UniqueID2=c(2,3,14,15,17,18,19),
                  Group=c(c(1,1,2,2,3,3,3)))

df <- tibble(UniqueID=c(1,2,3,14,15,17,18,19),
             Name=c("Schaltegger","Schaltegger","Schaltegger","Lüchinger","Lüchinger","Brandes","Branding","Branding"),
             Vorname=c("Chris","Christoph A.","Christoph","Simon","Simon D.","Leif","Leif","Live"))


# (i) Reshape data to have groups in wide format

data_out_all1 <- data_out_all %>% select(UniqueID1,Group)
data_out_all2 <- data_out_all %>% select(UniqueID2,Group)
data_out_all1 <- rename(data_out_all1, UniqueID=UniqueID1)      
data_out_all2 <- rename(data_out_all2, UniqueID=UniqueID2)      
out <- rbind(data_out_all1,data_out_all2)
out$FixCol <- out$UniqueID
out_new <- as.data.frame(dcast(setDT(out), UniqueID~ Group, value.var=c("FixCol"), fun=mean))


# (ii) Generate IDs

# ID_Final: Take first group match and neglect further connections
# IDFinalAlt: Recode all possible connections as the same ID
# Note: In between those two extremes there are a lot more possible connections. We should
#       deal with this in a later stage. 

no_groups <- dim(out_new)[2]-1
out_new$IDFinal <- NA
out_new$IDFinalAlt <- NA

for (i in 1:no_groups){
  out_new$IDFinal[!is.na(out_new[[1+i]]) & is.na(out_new$IDFinal)]  <- i
  out_new$IDFinalAlt[!is.na(out_new[[1+i]]) ]  <- min(out_new$IDFinal[!is.na(out_new[[1+i]])])
}
out_new


# (iii) Merge ID-File to Original File


(df_final <- left_join(df,out_new,key="UniqueID"))

