#install.packages("data.table")

library(dplyr)
library(readr)
library(data.table)

setwd("C:/Schmidlu/Dropbox/Record Linkage/Data/VRMandate")
bisnode <- fread("Bisnode_Person-Firmen_Geo.csv", sep = ";", encoding = "UTF-8")

bisnode <- bisnode %>% select(personenid, anrede, titel, vorname, nachname, W_PLZ4, wohnort,
                              bürgerort,geburtstag,firma,gremium, plzdomiziladresse, 
                              ortdomiziladresse, eintrittdatum, austrittdatum)




##SO-2007-0011

unsicher <- bisnode %>% filter(personenid==2778923)
View(unsicher)
# https://business-monitor.ch/de/p/daniel-bleuer-5105016


# BEJU-1995-0573

unsicher <- bisnode %>% filter(personenid==1515414)
View(unsicher)
# falsches Geburtsjahr, https://www.bielertagblatt.ch/nachrichten/biel/martin-zesiger-ist-der-neue-leiter-der-abteilung-bildung-und-kultur

## BL-2007-0028
unsicher <- bisnode %>% filter(personenid==1111362)
View(unsicher)


## GE-1991-0058
unsicher <- bisnode %>% filter(personenid==598809|personenid==425696)
View(unsicher)

## GR-1999-0020
unsicher <- bisnode %>% filter(personenid %in% c(1387722, 1575206, 2113035, 2245166, 2831076, 283019))
View(unsicher)
# andere Person: Prof. Willi Geiger

# NE-1995-0001
unsicher <- bisnode %>% filter(personenid %in% c(111726, 394896, 796882, 401527, 997662))
View(unsicher)
# https://business-monitor.ch/de/p/jean-pierre-authier-1109382

## ZH-1995-0489	
unsicher <- bisnode %>% filter(personenid %in% c(3268564))
View(unsicher)

## BEJU-1995-0185
unsicher <- bisnode %>% filter(personenid %in% c(9231, 42933, 3350704, 1249729))
View(unsicher)
# Merkmale des "1249729"-Eintrags stimmen mit anderne Einträgen überein

## BL-1987-0021	
unsicher <- bisnode %>% filter(personenid %in% c(6758, 3139524, 676403))
View(unsicher)
# Merkmale stimmen überein

## GE-2011-9152
unsicher <- bisnode %>% filter(personenid %in% c(2501893))
View(unsicher)
# Gemäss Youtube nicht gleiche Person (unterschiedliches Alter)

##  LU-1975-0050	
unsicher <- bisnode %>% filter(personenid %in%c(187887, 1561834, 309882, 3098231))
View(unsicher)
# Alexander Willi: gleiche Person

## LU-1999-0050
unsicher <- bisnode %>% filter(personenid %in%c(545164, 3512262, 316773))
View(unsicher)
# Otto Laubacher: gleiche Person

## ZH-2003-0423	
unsicher <- bisnode %>% filter(personenid %in%c(87350))
View(unsicher)
# Gleiche Adresse auf Monetas wie auf einer Seite zum Beruf 
# https://www.monetas.ch/de/671/SHAB-HR-Meldungen-14-07-1995-Kanton-Zuerich.htm?year=1995&shabnumber=135&canton=ZH+&date=14-07-1995
# https://www.yumpu.com/de/document/view/3376680/behorden-verwaltung-adressen-2007-stadt-wil
# Akeret René, Leiter Wiler Integrations- und Präventionsprojekte, ... Jung René, Musiklehrperson, Brühlbergstrasse 10, 8400 Winterthur (Vertreter ...

# GE-1995-0059
unsicher <- bisnode %>% filter(personenid %in%c(1528119))
View(unsicher)
# Politiker ist in Genf, VR-Mandatsträger in der Waadt
# Keine Anhaltspunkte, dass dies die gleiche Person ist
	
# VD-1995-0164
unsicher <- bisnode %>% filter(personenid %in% c(361966,391590,392420,877239,1365387,1371468,1527451,1964616,1968804,2838933,2968615,3338184,846674))
View(unsicher)
# 846674: Gleiches Firmendomizil wie andere Fälle
# 3338184: Gleiche Gegend des Firmendomizils wie anderere Fälle


# ZH-1995-0453	
unsicher <- bisnode %>% filter(personenid %in% c(2997534))
View(unsicher)
# Gleicher Job: http://andreas-meierhofer.ch/firma/ansprechpartner/index.html

# ZH-1995-0595	
unsicher <- bisnode %>% filter(personenid %in% c(1158222))
View(unsicher)
# Gleicher Wohnort

# ZH-2011-0225	
unsicher <- bisnode %>% filter(personenid %in% c(3279963))
View(unsicher)
# Philostudent http://gbs-schweiz.org/blog/author/lukas-gloor/

# BEJU-1991-0247
unsicher <- bisnode %>% filter(personenid %in% c(5526,511998,515667,164907,941225))
View(unsicher)
# Gleicher Wohnort, gleicher Heimatort, Mandate im Baubereich

# BEJU-2003-0266
unsicher <- bisnode %>% filter(personenid %in% c(143341,2315298,3009381))
View(unsicher)
# Jus -> Rechtsanwältin, gleicher Wohn- und Heimatort

# BEJU-2007-0251
unsicher <- bisnode %>% filter(personenid %in% c(705251, 1221557))
View(unsicher)
# Anderes Business, andere Gemeinde

# BL-2003-0049	
unsicher <- bisnode %>% filter(personenid %in% c(231455))
View(unsicher)
# Gleiche Wohngemeinde wie Firma (neu)

#GR-1999-0011-5	
unsicher <- bisnode %>% filter(personenid %in% c(2512100,79946))
View(unsicher)

## LU-2003-0074
unsicher <- bisnode %>% filter(personenid %in% c(1263451,1640691))
View(unsicher)

## SG-1999-0090
unsicher <- bisnode %>% filter(personenid %in% c(2162376,2543929))
View(unsicher)
# andere Heimat- und Wohnort

## ZH-2003-0391-5
unsicher <- bisnode %>% filter(personenid %in% c(1477286,2493258,2628224,3043589,3460250,1402276))
View(unsicher)
# Falscher Geburtstag in Bisnodedaten
# https://www.schadenanwaelte.ch/schadenanwaelte/david-husmann/

## BEJU-1991-0450	
unsicher <- bisnode %>% filter(personenid %in% c(994901,2486263,7969))
View(unsicher) 
# Beide Mandate im Immobilienbereich

## VD-1991-0041
unsicher <- bisnode %>% filter(personenid %in% c(633917))
View(unsicher) 
# Beruf passt zu Mandaten


## VS-2011-9030
unsicher <- bisnode %>% filter(personenid %in% c(1679445,1906854,2468150,1966187))
View(unsicher) 
# Gleicher Wohn- und Bürgerort

## ZG-2007-0027
unsicher <- bisnode %>% filter(personenid %in% c(310671,310674,558650,1504074,1958270,2047165,700935))
View(unsicher) 
# Gleiches Business
# https://www.yumpu.com/de/document/read/1107772/unser-neuer-internetauftritt-aktuelle-projekte-jego-ag

## ZH-1999-0768
unsicher <- bisnode %>% filter(personenid %in% c(1041924))
View(unsicher) 
# Anderer Heimatort, kein Beruf angegeben in Politikerdaten, Alter

## AG-2003-0079
unsicher <- bisnode %>% filter(personenid %in% c(554396,3091985))
View(unsicher) 
# Andere Person als Aargauer Regierungsrätin
# https://www.hogradirekt.ch/ueber-uns/2017_kleiner/


## TG-1995-0068
unsicher <- bisnode %>% filter(personenid %in% c(1224184,1224186,1269231))
View(unsicher) 
# Anderer Wohnort

# TI-2011-9046
unsicher <- bisnode %>% filter(personenid %in% c(2808878))
View(unsicher) 
# https://www.linkedin.com/in/alessandro-lucchini/?originalSubdomain=ch
	

# ZH-1999-0053
unsicher <- bisnode %>% filter(personenid %in% c(802346,748614))
View(unsicher) 
# anderer Wohnort, anderer Beruf

# ZH-2007-0178
unsicher <- bisnode %>% filter(personenid %in% c(1838491,2277664))
View(unsicher) 


# ZH-2007-0662-7	
unsicher <- bisnode %>% filter(personenid %in% c(541055,3539119))
View(unsicher) 

# AG-2011-0144	
unsicher <- bisnode %>% filter(personenid %in% c(576841))
View(unsicher)
bisnode %>% filter(firma %in% c("Robinex AG"))
# Ähnlicher Wohnort wie Firma, Heimatort	

# BEJU-1995-0073	
unsicher <- bisnode %>% filter(personenid %in% c(14957, 	3159242))
View(unsicher)
# Ähnliche Region

# BEJU-2011-0160
unsicher <- bisnode %>% filter(personenid %in% c(1262364,675039,1250458,1274291,1274366))
View(unsicher)
# Garten ist eine andere Person


# FR-2007-0057	
unsicher <- bisnode %>% filter(personenid %in% c(1409473))
View(unsicher)
bisnode %>% filter(firma %in% c("Garage Nicoli R. Sàrl"))
# Anderes Business, eher kein Familienunternehmen

# ZH-1991-0416-8	
unsicher <- bisnode %>% filter(personenid %in% c(1898618, 73565))
View(unsicher)
# Anderer Wohnort, anderes Business


# AG-2015-0005	
unsicher <- bisnode %>% filter(personenid %in% c(2182945, 2744073))
View(unsicher)
# Anderer Wohnort, anderes Business
# https://www.aargauerzeitung.ch/aargau/fricktal/bei-der-jakem-ag-in-muenchwilen-verlieren-80-mitarbeiter-ihren-job-128850665

# VS-1991-0020
unsicher <- bisnode %>% filter(personenid %in% c(342284, 2511366, 683413))
View(unsicher)
# Gleiche Region

# VS-2003-0047
unsicher <- bisnode %>% filter(personenid %in% c(3131028))
View(unsicher)
# Anderes Business als Beruf

# ZH-1995-0452	
unsicher <- bisnode %>% filter(personenid %in% c(813289))
View(unsicher)
# Schweizerzeit und SVP
	
# BEJU-1999-0306 
unsicher <- bisnode %>% filter(personenid %in% c(478828, 1486029, 3148398))
View(unsicher)
# https://home.kpmg/ch/de/home/contacts/r/martin-rohrbach.html
# https://www.linkedin.com/in/martin-rohrbach-bb0b2474/?originalSubdomain=ch


# BEJU-2015-0356
unsicher <- bisnode %>% filter(personenid %in% c(755674, 2163416, 2394953, 3328752, 3353513, 1250932))
View(unsicher)
# Dr-Titel über Mandate hinweg, Heimatorte sind angrenzend

# SG-1991-0034
unsicher <- bisnode %>% filter(personenid %in% c(918820))
View(unsicher)
# Nur Namen übereinstimmend	


# ZH-1991-0308
unsicher <- bisnode %>% filter(personenid %in% c(3209695))
View(unsicher)
# Anderer Wohnort und Business

# AG-1995-0068
unsicher <- bisnode %>% filter(personenid %in% c(712554,1092632))
View(unsicher)
# UBS Mandat ist andere Person

# BEJU-1991-0082
unsicher <- bisnode %>% filter(personenid %in% c(1087176))
View(unsicher)
# PdA und linker Verlag, andere Merkmale stimmen auch überein

# BEJU-1991-0177
unsicher <- bisnode %>% filter(personenid %in% c(14319,1172878))
View(unsicher)
# Andere Region, anderes Business, anderes Geburtsdatum

# BEJU-1995-0324
unsicher <- bisnode %>% filter(personenid %in% c(14323,19332 ,168045,964359,1219890 ,1460900,1525762 ,1557123 ,1631752,1962277,1969334,2624486,2681143,2709383,2994149,398881,895783,398881895783))
View(unsicher)
# alles im Schuhbusiness


# BEJU-1999-0305
unsicher <- bisnode %>% filter(personenid %in% c(1236988,42742))
View(unsicher)
# Polizist ist eher nicht im Agrarbereich tätig

# GE-2011-0030
unsicher <- bisnode %>% filter(personenid %in% c(1923448))
View(unsicher)
# Merkmale stimmen überein

# LU-1999-0087
unsicher <- bisnode %>% filter(personenid %in% c(2552928))
View(unsicher)
# Merkmale stimmen nicht überein


# LU-2003-0021
unsicher <- bisnode %>% filter(personenid %in% c(318715,3096958,232900))
View(unsicher)
# Andere Mandate, andere Region

# SG-1991-0010
unsicher <- bisnode %>% filter(personenid %in% c(1000870, 1025655))
View(unsicher)
# https://www.trauerportal-ostschweiz.ch/traueranzeige/urs-bernhardsgruetter
# https://krj.ch/ansprechpartner/urs-bernhardsgruetter/


# SG-2003-0116
unsicher <- bisnode %>% filter(personenid %in% c(1546352,485396,3292133))
View(unsicher)
# https://www.linkedin.com/in/markus-ritter-0ab01866/?originalSubdomain=ch

# SG-2011-0174
unsicher <- bisnode %>% filter(personenid %in% c(291878,742999,1658922,1305972))
View(unsicher)
# Andere Region

# SZ-1995-0012
unsicher <- bisnode %>% filter(personenid %in% c(1255659))
View(unsicher)
# Andererer Beruf, SDler hat wohl eher kein Buchhandlungsmandat

# TG-2003-0052
unsicher <- bisnode %>% filter(personenid %in% c(3582764))
View(unsicher)
# Andererer Beruf, andere Region	

# VD-1995-0126
unsicher <- bisnode %>% filter(personenid %in% c(1490466, 575429))
View(unsicher)
# Andererer Beruf, andere Region	



# Arbeitsteilung Mitte/Ende Batch 11, 3. März 2020

setwd("C:/Schmidlu/Dropbox/Projekt Nationalräte/02_Processed_data/12_Record_linkage")
check_batch_11_15 <- fread("help_check_batch_11_15.csv", sep = ";", encoding = "UTF-8") %>% rename(personenid=Sicher)

unsicher <- bisnode %>% 
            filter(personenid %in% check_batch_11_15$personenid) %>% 
            left_join(check_batch_11_15,by=c("personenid"))%>%
            select(ID,personenid,everything())%>%
            arrange(batch,ID)
View(unsicher)

write_excel_csv(unsicher, "bisnode_tocheck_batch_11_15.csv")














	





