
library(readxl)
library(xlsx)

try(setwd("D:/SchmidLu/Dropbox/Projekt Nationalräte/1 Data/1_data/14_jarckspuehlerwidmer"))
try(setwd("C:/Users/Lukas/Dropbox/Projekt Nationalräte/1 Data/1_data/14_jarckspuehlerwidmer"))


# (i) Read-in data


df <-  read_excel("JarckSpuehlerWidmer_out_wide.xlsx")

head(df)
dim(df)

set.seed(12345)

# (ii) Define common block

df$id <- c(1:dim(df)[1])
common_sample <- sample(1:dim(df)[1],320)
df_common_block <- df[is.element(df$id,common_sample),]
df_common_block <- df_common_block[order(df_common_block$ID),]
df_common_block$sampleordering <- c(1:320)
df_notcommon_block <- df[is.element(df$id,common_sample)==FALSE,]
df_notcommon_block$sampleordering <- c(1:dim(df_notcommon_block)[1])

write.csv2(df_common_block, "./Outfiles/AllCommonblock.csv",na = "")  
write.csv2(df_notcommon_block, "./Outfiles/AllNotCommonblock.csv",na = "")  


# (iii) Get packages for Jarck and Widmer


for (i in 1:4){
  notcommon_sample <- sample(df_notcommon_block$id,920)
  df_out1 <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==TRUE,]
  df_out2 <- df_common_block[is.element(df_common_block$sampleordering,c(((i-1)*80+1):(i*80))),]  
  df_out <- rbind(df_out1,df_out2) 
  df_out <- df_out[order(df_out$ID),]
  file <-  paste("./Outfiles/GeburtsTodesdaten_Jana_Jarck_",i,".csv",sep="")
  write.csv2(df_out, file,na = "")  
  df_notcommon_block <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==FALSE,]
}



for (i in 1:4){
  notcommon_sample <- sample(df_notcommon_block$id,920)
  df_out1 <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==TRUE,]
  df_out2 <- df_common_block[is.element(df_common_block$sampleordering,c(((i-1)*80+1):(i*80))),]  
  df_out <- rbind(df_out1,df_out2) 
  df_out <- df_out[order(df_out$ID),]
  file <-  paste("./Outfiles/GeburtsTodesdaten_Livia_Widmer_",i,".csv",sep="")
  write.csv2(df_out, file,na = "")  
  df_notcommon_block <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==FALSE,]
}


for (i in 1:4){
  notcommon_sample <- sample(df_notcommon_block$id,920)
  df_out1 <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==TRUE,]
  df_out2 <- df_common_block[is.element(df_common_block$sampleordering,c(((i-1)*80+1):(i*80))),]  
  df_out <- rbind(df_out1,df_out2) 
  df_out <- df_out[order(df_out$ID),]
  file <-  paste("./Outfiles/GeburtsTodesdaten_Yves_Spuehler_",i,".csv",sep="")
  write.csv2(df_out, file,na = "")  
  df_notcommon_block <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==FALSE,]
}

# (iv) Get packages for Zivi

for (i in 1:15){
  notcommon_sample <- sample(df_notcommon_block$id,980)
  df_out1 <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==TRUE,]
  df_out2 <- df_common_block[is.element(df_common_block$sampleordering,c(((i-1)*20+1):(i*20))),]  
  df_out <- rbind(df_out1,df_out2) 
  df_out <- df_out[order(df_out$ID),]
  file <-  paste("./Outfiles/GeburtsTodesdaten_Zivi_",i,".csv",sep="")
  write.csv2(df_out, file,na = "")  
  df_notcommon_block <- df_notcommon_block[is.element(df_notcommon_block$id,notcommon_sample)==FALSE,]
}

i=16
df_out1 <- df_notcommon_block
df_out2 <- df_common_block[is.element(df_common_block$sampleordering,c(((i-1)*20+1):(i*20))),]  
df_out <- rbind(df_out1,df_out2) 
file <-  paste("./Outfiles/GeburtsTodesdaten_Zivi_",i,".csv",sep="")
write.csv2(df_out, file,na = "") 




# Jarck/Widmer: ca. 4 Tausender-Päckli
# jedes Päckli: 76 gemeinsame Kandidaten, jeweils dieselben; dazu 924 unterschiedliche Kandidaten
# 
#
# Zivi: ca. 19 Tausenderpäckli
# jedes Päckli: 16 gemeinsame Kandidaten, jeweils dieselben; 984 unterschiedliche
# 
# Vorgehen: 304 ordnen, beim Zivi in der ersten Tranche Personen 1-16, zweiten Tranche Personen
# 17-32 etc. einfügen, bei Jarck/Widmer 1-76 in ersten Tranche, 77-152 in zweiten Tranche etc. 

# Am Schluss nach ID ordnen

# Alle Pakete (das 304er und die 9xxer) werden zufällig und ohne zurücklegen gezogen

26646+304*2

# Letztes Paket wird grösser!! (1254 Personen statt 1000)

# Daten basieren auf JarckSpuehlerWidmer_out_wide, zu ergänzen jeweils um 
# Spalte für Quelle Geburtsdatum und Quelle Todesdatum, löschen: Wohnorte/Bemerkungen (letzte Spalten)

# Vorgehen: Suche Geburts- und Todesdatum
# 1. HLS
# 2. Todesanzeigenportal
# 3. Wikipedia
# sobald ein Datum gefunden: Eintrag Eingabe und Quellennummer, weiter zur nächsten Person ohne Sichtung weiterer Quellen
