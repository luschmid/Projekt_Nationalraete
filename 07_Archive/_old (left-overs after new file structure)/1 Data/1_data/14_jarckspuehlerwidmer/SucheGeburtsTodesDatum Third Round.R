
library(readxl)
library(xlsx)

try(setwd("D:/SchmidLu/Dropbox/Projekt Nationalräte/1 Data/1_data/14_jarckspuehlerwidmer"))
try(setwd("C:/Users/Lukas/Dropbox/Projekt Nationalräte/1 Data/1_data/14_jarckspuehlerwidmer"))


# (i) Read-in data


df <-  read_excel("Zivis_out_wide.xlsx")

head(df)
dim(df)

set.seed(12345)

# (ii) Define common block

df$id <- c(1:dim(df)[1])
common_sample <- sample(1:dim(df)[1],276)
df_common_block <- df[is.element(df$id,common_sample),]
df_common_block <- df_common_block[order(df_common_block$ID),]
df_common_block$sampleordering <- c(1:276)
df_notcommon_block <- df[is.element(df$id,common_sample)==FALSE,]
df_notcommon_block$sampleordering <- c(1:dim(df_notcommon_block)[1])

write.csv2(df_common_block, "./Outfiles3/AllCommonblock.csv",na = "")  
write.csv2(df_notcommon_block, "./Outfiles3/AllNotCommonblock.csv",na = "")  


# (iii) Get packages for Silo 1


for (i in 1:12){
  notcommon_sample <- sample(df_notcommon_block$id,977)
  df_out1 <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==TRUE,]
  df_out2 <- df_common_block[is.element(df_common_block$sampleordering,c(((i-1)*23+1):(i*23))),]  
  df_out <- rbind(df_out1,df_out2) 
  df_out <- df_out[order(df_out$ID),]
  file <-  paste("./Outfiles3/GeburtsTodesdaten3_A_",i,".csv",sep="")
  write.csv2(df_out, file,na = "")  
  df_notcommon_block <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==FALSE,]
}



# (iv) Get packages for Silo 2


for (i in 1:11){
  notcommon_sample <- sample(df_notcommon_block$id,977)
  df_out1 <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==TRUE,]
  df_out2 <- df_common_block[is.element(df_common_block$sampleordering,c(((i-1)*23+1):(i*23))),]  
  df_out <- rbind(df_out1,df_out2) 
  df_out <- df_out[order(df_out$ID),]
  file <-  paste("./Outfiles3/GeburtsTodesdaten3_B_",i,".csv",sep="")
  write.csv2(df_out, file,na = "")  
  df_notcommon_block <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==FALSE,]
}


df_out1 <- df_notcommon_block
df_out2 <- df_common_block[is.element(df_common_block$sampleordering,c(((i-1)*23+1):(i*23))),]  
df_out <- rbind(df_out1,df_out2) 
df_out <- df_out[order(df_out$ID),]
file <-  paste("./Outfiles3/GeburtsTodesdaten3_B_",12,".csv",sep="")
write.csv2(df_out, file,na = "")  
df_notcommon_block <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==FALSE,]

