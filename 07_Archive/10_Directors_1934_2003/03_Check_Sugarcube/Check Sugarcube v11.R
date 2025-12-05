

setwd("E:/Corrected (Final Scans)")

no.folders <- length(list.dirs())
name.folders <- list.dirs()
name.folders <- name.folders[name.folders!="."]


# first draw: seeds <- c(1234,12348,12345658,123947)
# second draw: seeds <- c(12345,125348,123445658,1523947)
seeds <- c(123434,12344228,1234123125658,12391211247)



N=150

df <- data.frame(jahrgang=matrix(NA,4*N),dateiname=matrix(NA,4*N),	eintrag=matrix(NA,4*N),batch=matrix(NA,4*N))


for (batch in 1:4){
set.seed(seeds[batch])  
for (i in 1:N){
name.folder <- sample(name.folders,1)
setwd(paste("E:/Corrected (Final Scans)/",name.folder,sep=""))
df$jahrgang[(batch-1)*N+i] <- name.folder
df$dateiname[(batch-1)*N+i] <- sample(list.files(pattern=".tif"),1)
df$eintrag[(batch-1)*N+i] <- floor(runif(1,1,200))
df$batch[(batch-1)*N+i] <- batch
}
}


df

df <- df[order(df$batch,df$jahrgang,df$dateiname,df$eintrag),]


setwd("D:/SchmidLu/Dropbox/Projekt Nationalräte/2 Data Collection/Verwaltungsratsmandate")
write.table(df, file = "Check Sugarcube v11_Zufallsstichprobe3.csv",sep = ";", row.names = F)

