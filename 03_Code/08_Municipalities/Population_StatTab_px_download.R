
# load libraries
library(pxR) # this is needed to read px-files, the format used on Stat-Tab
library(tidyverse)
library(stringr)
library(haven)
try(setwd("D:/SchmidLu/Dropbox/Projekt Nationalräte"))

# first part of the url (this is pasted into the part where the file is read, next step)
bfs.url <- "https://www.pxweb.bfs.admin.ch/DownloadFile.aspx?"

# read the px-file:
bevoelkerung.px <- read.px(paste0(bfs.url, "file=px-x-0102020000_201"),
                           encoding="cp1252")  # I needed this encoding option 

detach(package:pxR, unload = TRUE)
detach(package:plyr, unload = TRUE)
bevoelkerung.df <- as.data.frame(bevoelkerung.px)  #  save as data frame; note: for Windows, you may need to modify this  

bevoelkerung.df %>% distinct(Kanton.......Bezirk........Gemeinde.........) %>% print(max=50)

rm(bevoelkerung.px)
str(bevoelkerung.df)

bevoelkerung.df.final <- bevoelkerung.df %>% 
  filter(Demographische.Komponente=="Bestand am 1. Januar", 
         Geschlecht=="Geschlecht - Total",
         Staatsangehörigkeit..Kategorie.=="Staatsangehörigkeit - Total",
         !Kanton.......Bezirk........Gemeinde.........%in% c("Schweiz"),  
         !grepl("-",Kanton.......Bezirk........Gemeinde.........), 
         !grepl(">>",Kanton.......Bezirk........Gemeinde.........))%>%
         #Jahr %in% seq(1983,2015,4))%>%
  mutate(municipalityno=gsub(pattern = "[^0-9]",replacement = "",x = Kanton.......Bezirk........Gemeinde.........), 
         municipality=gsub(pattern = "[0-9.]",replacement = "",Kanton.......Bezirk........Gemeinde.........)) %>%
  rename(year=Jahr,
         population=value) %>%
  select(municipality,municipalityno,year,population)%>%
  as_tibble()

bevoelkerung.df.final$municipalityno <- as.numeric(bevoelkerung.df.final$municipalityno)
bevoelkerung.df.final$year <- as.numeric(as.character(bevoelkerung.df.final$year))

bevoelkerung.df.final %>% distinct(municipality) %>% print(n=50)


write_dta(bevoelkerung.df.final, "./02_Processed_data/08_Municipalities/Population_1981_2017.dta")
