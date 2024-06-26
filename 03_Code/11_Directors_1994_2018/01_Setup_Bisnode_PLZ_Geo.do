clear 	
cap log close
set more off	  

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"

clear


**** Read-in Bisnode data ****
******************************

import delimited "$path\01_Raw_data\11_Directors_1994_2018\Uni Fribourg Personendaten 20-07-2018.csv", delimiters(";") varnames(1) clear
save "$path\02_Processed_data\11_Directors_1994_2018\Personendaten.dta", replace

import delimited "$path\01_Raw_data\11_Directors_1994_2018\Uni Fribourg Firmendaten 20-07-2018.csv", delimiters(";") varnames(1) clear
save "$path\02_Processed_data\11_Directors_1994_2018\Firmendaten.dta", replace



**** Add geo-coordinates to Bisnode data ****
*********************************************


*** Merge PLZ-Geo data with Bisnode data ***

** Prepare Bisnode datset **
use "$path\02_Processed_data\11_Directors_1994_2018\Personendaten.dta", clear
gen PLZWmiss = 1 if plzwohnort==""
destring plzwohnort, g(PLZ4) force
sort PLZ4
// correct PLZ in Bisnode for those that do not originally merge with Post-Data
replace PLZ4 = 1304 if PLZ4 == 1118
replace PLZ4 = 2123 if PLZ4 == 1221
replace PLZ4 = 1237 if PLZ4 == 1238
replace PLZ4 = 1282 if PLZ4 == 1282
replace PLZ4 = 1537 if PLZ4 == 1487
replace PLZ4 = 1538 if PLZ4 == 1488
replace PLZ4 = 1607 if PLZ4 == 1501
replace PLZ4 = 1607 if PLZ4 == 1501
replace PLZ4 = 1610 if PLZ4 == 1502
replace PLZ4 = 1610 if PLZ4 == 1502
replace PLZ4 = 1589 if PLZ4 == 1594
replace PLZ4 = 1656 if PLZ4 == 1655
replace PLZ4 = 1695 if PLZ4 == 1743
replace PLZ4 = 1867 if PLZ4 == 1855
replace PLZ4 = 1863 if PLZ4 == 1861
replace PLZ4 = 1882 if PLZ4 == 1881
replace PLZ4 = 1893 if PLZ4 == 1883
replace PLZ4 = 1934 if PLZ4 == 1935
replace PLZ4 = 1938 if PLZ4 == 1939
replace PLZ4 = 1955 if PLZ4 == 1956
replace PLZ4 = 1690 if PLZ4 == 1960
replace PLZ4 = 1690 if PLZ4 == 1960
replace PLZ4 = 1981 if PLZ4 == 1989
replace PLZ4 = 1063 if PLZ4 == 2032
replace PLZ4 = 2024 if PLZ4 == 2040
replace PLZ4 = 1566 if PLZ4 == 2040
replace PLZ4 = 1566 if PLZ4 == 2041
replace PLZ4 = 2046 if PLZ4 == 2045
replace PLZ4 = 2046 if PLZ4 == 2047
replace PLZ4 = 2087 if PLZ4 == 2078
replace PLZ4 = 1740 if PLZ4 == 2211
replace PLZ4 = 1786 if PLZ4 == 2281
replace PLZ4 = 2950 if PLZ4 == 2650
replace PLZ4 = 2046 if PLZ4 == 2837
replace PLZ4 = 3110 if PLZ4 == 3310
replace PLZ4 = 3615 if PLZ4 == 3519
replace PLZ4 = 3762 if PLZ4 == 3759
replace PLZ4 = 1865 if PLZ4 == 5411
replace PLZ4 = 6598 if PLZ4 == 5698
replace PLZ4 = 4153 if PLZ4 == 5743
replace PLZ4 = 1462 if PLZ4 == 5939
replace PLZ4 = 3946 if PLZ4 == 6119
replace PLZ4 = 6276 if PLZ4 == 6272
replace PLZ4 = 6280 if PLZ4 == 6282
replace PLZ4 = 5605 if PLZ4 == 6505
replace PLZ4 = 6805 if PLZ4 == 6508
replace PLZ4 = 8152 if PLZ4 == 8150
replace PLZ4 = 9326 if PLZ4 == 8326
replace PLZ4 = 8664 if PLZ4 == 8557
replace PLZ4 = 8645 if PLZ4 == 8647
replace PLZ4 = 8645 if PLZ4 == 8648
replace PLZ4 = 8645 if PLZ4 == 8649
replace PLZ4 = 6340 if PLZ4 == 8944
replace PLZ4 = 8135 if PLZ4 == 8944
replace PLZ4 = 8330 if PLZ4 == 9330
replace PLZ4 = 9434 if PLZ4 == 9433
replace PLZ4 = 8370 if PLZ4 == 9572

replace PLZ4 = 8564 if wohnort == "Lipperswil" & PLZ4 == 8664

/*// missing PLZ for Wohnort in Schweiz (merge on "Wohnort" using uniquePLZ observations from Swiss Post)
preserve // empty PLZ4 for Wohnort
keep if PLZ4==. 
gen ch = 0
replace ch = 1 if wohnland=="Schweiz"
drop PLZ4
gen count = 1
collapse (sum) count, by(wohnort ch)
duplicates tag wohnort, gen(dup)
rename wohnort gdename
sort gdename
save "${path}02_Processed_data\11_Directors_1994_2018\Bisnode_missingWohnPLZ_temp1.dta", replace
use "$path\02_Processed_data\08_Municipalities/SwissPost-PLZ_GdeNr2018-Geo.dta", clear // generate a list of Gemeindenamen with PLZ to merge on missings
keep if yPLZout>=2018
duplicates drop gdename, force
sort gdename
save "$path\02_Processed_data\11_Directors_1994_2018/Gde-PLZ-List_temp.dta", replace
use "$path\02_Processed_data\11_Directors_1994_2018/Bisnode_missingWohnPLZ_temp1.dta", clear
merge m:1 gdename using "$path\02_Processed_data\11_Directors_1994_2018/Gde-PLZ-List_temp.dta"
keep if _merge==3
keep gdename PLZ4
duplicates drop gdename PLZ4, force
save "$path\02_Processed_data\11_Directors_1994_2018/Bisnode_missingWohnPLZ-merge_temp1.dta", replace
restore
*/

replace PLZ4= 5646 if wohnort == "Abtwil" & PLZWmiss == 1
replace PLZ4= 6043 if wohnort == "Adligenswil" & PLZWmiss == 1
replace PLZ4= 8452 if wohnort == "Adlikon" & PLZWmiss == 1
replace PLZ4= 3703 if wohnort == "Aeschi bei Spiez" & PLZWmiss == 1
replace PLZ4= 3672 if wohnort == "Aeschlen" & PLZWmiss == 1
replace PLZ4= 6927 if wohnort == "Agra" & PLZWmiss == 1
replace PLZ4= 1860 if wohnort == "Aigle" & PLZWmiss == 1
replace PLZ4= 7472 if wohnort == "Albula/Alvra" & PLZWmiss == 1
replace PLZ4= 3112 if wohnort == "Allmendingen" & PLZWmiss == 1
replace PLZ4= 6010 if wohnort == "Alpnach" & PLZWmiss == 1
replace PLZ4= 1715 if wohnort == "Alterswil" & PLZWmiss == 1
replace PLZ4= 9450 if wohnort == "Altstätten" & PLZWmiss == 1
replace PLZ4= 8873 if wohnort == "Amden" & PLZWmiss == 1
replace PLZ4= 5600 if wohnort == "Ammerswil" & PLZWmiss == 1
replace PLZ4= 1247 if wohnort == "Anières" & PLZWmiss == 1
replace PLZ4= 3961 if wohnort == "Anniviers" & PLZWmiss == 1
replace PLZ4= 6532 if wohnort == "Arbedo-Castione" & PLZWmiss == 1
replace PLZ4= 6939 if wohnort == "Arosio" & PLZWmiss == 1
replace PLZ4= 1273 if wohnort == "Arzier-Le Muids" & PLZWmiss == 1
replace PLZ4= 6612 if wohnort == "Ascona" & PLZWmiss == 1
replace PLZ4= 1170 if wohnort == "Aubonne" & PLZWmiss == 1
replace PLZ4= 5105 if wohnort == "Auenstein" & PLZWmiss == 1
replace PLZ4= 4302 if wohnort == "Augst" & PLZWmiss == 1
replace PLZ4= 1742 if wohnort == "Autigny" & PLZWmiss == 1
replace PLZ4= 1285 if wohnort == "Avusy" & PLZWmiss == 1
replace PLZ4= 6346 if wohnort == "Baar" & PLZWmiss == 1
replace PLZ4= 8164 if wohnort == "Bachs" & PLZWmiss == 1
replace PLZ4= 5400 if wohnort == "Baden" & PLZWmiss == 1
replace PLZ4= 1947 if wohnort == "Bagnes" & PLZWmiss == 1
replace PLZ4= 4525 if wohnort == "Balm bei Günsberg" & PLZWmiss == 1
replace PLZ4= 1257 if wohnort == "Bardonnex" & PLZWmiss == 1
replace PLZ4= 1666 if wohnort == "Bas-Intyamon" & PLZWmiss == 1
replace PLZ4= 4089 if wohnort == "Basel" & PLZWmiss == 1
replace PLZ4= 1268 if wohnort == "Begnins" & PLZWmiss == 1
replace PLZ4= 1293 if wohnort == "Bellevue" & PLZWmiss == 1
replace PLZ4= 6582 if wohnort == "Bellinzona" & PLZWmiss == 1
replace PLZ4= 1773 if wohnort == "Belmont-Broye" & PLZWmiss == 1
replace PLZ4= 7482 if wohnort == "Bergün Filisur" & PLZWmiss == 1
replace PLZ4= 8222 if wohnort == "Beringen" & PLZWmiss == 1
replace PLZ4= 3376 if wohnort == "Berken" & PLZWmiss == 1
replace PLZ4= 3341 if wohnort == "Bern" & PLZWmiss == 1
replace PLZ4= 1233 if wohnort == "Bernex" & PLZWmiss == 1
replace PLZ4= 2544 if wohnort == "Bettlach" & PLZWmiss == 1
replace PLZ4= 1880 if wohnort == "Bex" & PLZWmiss == 1
replace PLZ4= 4562 if wohnort == "Biberist" & PLZWmiss == 1
replace PLZ4= 4105 if wohnort == "Biel-Benken" & PLZWmiss == 1
replace PLZ4= 4102 if wohnort == "Binningen" & PLZWmiss == 1
replace PLZ4= 9220 if wohnort == "Bischofszell" & PLZWmiss == 1
replace PLZ4= 6816 if wohnort == "Bissone" & PLZWmiss == 1
replace PLZ4= 3919 if wohnort == "Blatten" & PLZWmiss == 1
replace PLZ4= 1807 if wohnort == "Blonay" & PLZWmiss == 1
replace PLZ4= 1279 if wohnort == "Bogis-Bossey" & PLZWmiss == 1
replace PLZ4= 4556 if wohnort == "Bolken" & PLZWmiss == 1
replace PLZ4= 7606 if wohnort == "Bondo" & PLZWmiss == 1
replace PLZ4= 1172 if wohnort == "Bougy-Villars" & PLZWmiss == 1
replace PLZ4= 1091 if wohnort == "Bourg-en-Lavaux" & PLZWmiss == 1
replace PLZ4= 1034 if wohnort == "Boussens" & PLZWmiss == 1
replace PLZ4= 6932 if wohnort == "Breganzona" & PLZWmiss == 1
replace PLZ4= 6837 if wohnort == "Breggia" & PLZWmiss == 1
replace PLZ4= 7159 if wohnort == "Breil/Brigels" & PLZWmiss == 1
replace PLZ4= 5620 if wohnort == "Bremgarten (AG)" & PLZWmiss == 1
replace PLZ4= 3047 if wohnort == "Bremgarten bei Bern" & PLZWmiss == 1
replace PLZ4= 3900 if wohnort == "Brig-Glis" & PLZWmiss == 1
replace PLZ4= 5200 if wohnort == "Brugg" & PLZWmiss == 1
replace PLZ4= 6827 if wohnort == "Brusino Arsizio" & PLZWmiss == 1
replace PLZ4= 2555 if wohnort == "Brügg" & PLZWmiss == 1
replace PLZ4= 4587 if wohnort == "Buchegg" & PLZWmiss == 1
replace PLZ4= 3615 if wohnort == "Buchholterberg" & PLZWmiss == 1
replace PLZ4= 1630 if wohnort == "Bulle" & PLZWmiss == 1
replace PLZ4= 6374 if wohnort == "Buochs" & PLZWmiss == 1
replace PLZ4= 3400 if wohnort == "Burgdorf" & PLZWmiss == 1
replace PLZ4= 3323 if wohnort == "Bäriswil" & PLZWmiss == 1
replace PLZ4= 5225 if wohnort == "Bözberg" & PLZWmiss == 1
replace PLZ4= 3274 if wohnort == "Bühl" & PLZWmiss == 1
replace PLZ4= 9615 if wohnort == "Bütschwil-Ganterschwil" & PLZWmiss == 1
replace PLZ4= 5619 if wohnort == "Büttikon" & PLZWmiss == 1
replace PLZ4= 6814 if wohnort == "Cadempino" & PLZWmiss == 1
replace PLZ4= 6543 if wohnort == "Calanca" & PLZWmiss == 1
replace PLZ4= 6760 if wohnort == "Campello" & PLZWmiss == 1
replace PLZ4= 6914 if wohnort == "Carona" & PLZWmiss == 1
replace PLZ4= 7422 if wohnort == "Cazis" & PLZWmiss == 1
replace PLZ4= 6330 if wohnort == "Cham" & PLZWmiss == 1
replace PLZ4= 1637 if wohnort == "Charmey" & PLZWmiss == 1
replace PLZ4= 1474 if wohnort == "Cheyres-Châbles" & PLZWmiss == 1
replace PLZ4= 7075 if wohnort == "Churwalden" & PLZWmiss == 1
replace PLZ4= 1474 if wohnort == "Châbles" & PLZWmiss == 1
replace PLZ4= 1224 if wohnort == "Chêne-Bougeries" & PLZWmiss == 1
replace PLZ4= 1893 if wohnort == "Collombey-Muraz" & PLZWmiss == 1
replace PLZ4= 1903 if wohnort == "Collonges" & PLZWmiss == 1
replace PLZ4= 1291 if wohnort == "Commugny" & PLZWmiss == 1
replace PLZ4= 1296 if wohnort == "Coppet" & PLZWmiss == 1
replace PLZ4= 2036 if wohnort == "Corcelles-Cormondrèche" & PLZWmiss == 1
replace PLZ4= 1782 if wohnort == "Cormagens" & PLZWmiss == 1
replace PLZ4= 2087 if wohnort == "Cornaux" & PLZWmiss == 1
replace PLZ4= 1304 if wohnort == "Cossonay" & PLZWmiss == 1
replace PLZ4= 1796 if wohnort == "Courgevaux" & PLZWmiss == 1
replace PLZ4= 1791 if wohnort == "Courtaman" & PLZWmiss == 1
replace PLZ4= 1096 if wohnort == "Cully" & PLZWmiss == 1
replace PLZ4= 7270 if wohnort == "Davos" & PLZWmiss == 1
replace PLZ4= 8253 if wohnort == "Diessenhofen" & PLZWmiss == 1
replace PLZ4= 1304 if wohnort == "Dizy" & PLZWmiss == 1
replace PLZ4= 7418 if wohnort == "Domleschg" & PLZWmiss == 1
replace PLZ4= 1041 if wohnort == "Dommartin" & PLZWmiss == 1
replace PLZ4= 4143 if wohnort == "Dornach" & PLZWmiss == 1
replace PLZ4= 4558 if wohnort == "Drei Höfe" & PLZWmiss == 1
replace PLZ4= 4658 if wohnort == "Däniken" & PLZWmiss == 1
replace PLZ4= 8114 if wohnort == "Dänikon" & PLZWmiss == 1
replace PLZ4= 8132 if wohnort == "Egg" & PLZWmiss == 1
replace PLZ4= 9315 if wohnort == "Egnach" & PLZWmiss == 1
replace PLZ4= 6020 if wohnort == "Emmen" & PLZWmiss == 1
replace PLZ4= 5304 if wohnort == "Endingen" & PLZWmiss == 1
replace PLZ4= 6373 if wohnort == "Ennetbürgen" & PLZWmiss == 1
replace PLZ4= 1066 if wohnort == "Epalinges" & PLZWmiss == 1
replace PLZ4= 8586 if wohnort == "Erlen" & PLZWmiss == 1
replace PLZ4= 8726 if wohnort == "Ernetschwil" & PLZWmiss == 1
replace PLZ4= 6472 if wohnort == "Erstfeld" & PLZWmiss == 1
replace PLZ4= 8360 if wohnort == "Eschlikon" & PLZWmiss == 1
replace PLZ4= 6182 if wohnort == "Escholzmatt-Marbach" & PLZWmiss == 1
replace PLZ4= 1489 if wohnort == "Estavayer" & PLZWmiss == 1
replace PLZ4= 4107 if wohnort == "Ettingen" & PLZWmiss == 1
replace PLZ4= 3617 if wohnort == "Fahrni" & PLZWmiss == 1
replace PLZ4= 8320 if wohnort == "Fehraltorf" & PLZWmiss == 1
replace PLZ4= 8552 if wohnort == "Felben" & PLZWmiss == 1
replace PLZ4= 6938 if wohnort == "Fescoggia" & PLZWmiss == 1
replace PLZ4= 1044 if wohnort == "Fey" & PLZWmiss == 1
replace PLZ4= 6145 if wohnort == "Fischbach" & PLZWmiss == 1
replace PLZ4= 8376 if wohnort == "Fischingen" & PLZWmiss == 1
replace PLZ4= 2114 if wohnort == "Fleurier" & PLZWmiss == 1
replace PLZ4= 7017 if wohnort == "Flims" & PLZWmiss == 1
replace PLZ4= 6173 if wohnort == "Flühli" & PLZWmiss == 1
replace PLZ4= 3636 if wohnort == "Forst" & PLZWmiss == 1
replace PLZ4= 8808 if wohnort == "Freienbach" & PLZWmiss == 1
replace PLZ4= 8428 if wohnort == "Freienstein-Teufen" & PLZWmiss == 1
replace PLZ4= 1700 if wohnort == "Fribourg" & PLZWmiss == 1
replace PLZ4= 1532 if wohnort == "Fétigny" & PLZWmiss == 1
replace PLZ4= 6579 if wohnort == "Gambarogno" & PLZWmiss == 1
replace PLZ4= 1272 if wohnort == "Genolier" & PLZWmiss == 1
replace PLZ4= 1200 if wohnort == "Genève" & PLZWmiss == 1
replace PLZ4= 4563 if wohnort == "Gerlafingen" & PLZWmiss == 1
replace PLZ4= 1728 if wohnort == "Gibloux" & PLZWmiss == 1
replace PLZ4= 1762 if wohnort == "Givisiez" & PLZWmiss == 1
replace PLZ4= 3207 if wohnort == "Golaten" & PLZWmiss == 1
replace PLZ4= 3989 if wohnort == "Goms" & PLZWmiss == 1
replace PLZ4= 8274 if wohnort == "Gottlieben" & PLZWmiss == 1
replace PLZ4= 1376 if wohnort == "Goumoëns" & PLZWmiss == 1
replace PLZ4= 3376 if wohnort == "Graben" & PLZWmiss == 1
replace PLZ4= 3989 if wohnort == "Grafschaft" & PLZWmiss == 1
replace PLZ4= 1091 if wohnort == "Grandvaux" & PLZWmiss == 1
replace PLZ4= 1686 if wohnort == "Grangettes" & PLZWmiss == 1
replace PLZ4= 6929 if wohnort == "Gravesano" & PLZWmiss == 1
replace PLZ4= 3784 if wohnort == "Gsteig" & PLZWmiss == 1
replace PLZ4= 7545 if wohnort == "Guarda" & PLZWmiss == 1
replace PLZ4= 8523 if wohnort == "Hagenbuch" & PLZWmiss == 1
replace PLZ4= 3419 if wohnort == "Hasle bei Burgdorf" & PLZWmiss == 1
replace PLZ4= 8773 if wohnort == "Haslen" & PLZWmiss == 1
replace PLZ4= 5212 if wohnort == "Hausen bei Brugg" & PLZWmiss == 1
replace PLZ4= 2855 if wohnort == "Haute-Sorne" & PLZWmiss == 1
replace PLZ4= 2068 if wohnort == "Hauterive" & PLZWmiss == 1
replace PLZ4= 8908 if wohnort == "Hedingen" & PLZWmiss == 1
replace PLZ4= 9514 if wohnort == "Heiligkreuz" & PLZWmiss == 1
replace PLZ4= 8444 if wohnort == "Henggart" & PLZWmiss == 1
replace PLZ4= 6133 if wohnort == "Hergiswil bei Willisau" & PLZWmiss == 1
replace PLZ4= 9102 if wohnort == "Herisau" & PLZWmiss == 1
replace PLZ4= 6024 if wohnort == "Hildisrieden" & PLZWmiss == 1
replace PLZ4= 6281 if wohnort == "Hochdorf" & PLZWmiss == 1
replace PLZ4= 8242 if wohnort == "Hofen" & PLZWmiss == 1
replace PLZ4= 4112 if wohnort == "Hofstetten-Flüh" & PLZWmiss == 1
replace PLZ4= 3622 if wohnort == "Homberg" & PLZWmiss == 1
replace PLZ4= 8507 if wohnort == "Homburg" & PLZWmiss == 1
replace PLZ4= 6038 if wohnort == "Honau" & PLZWmiss == 1
replace PLZ4= 9326 if wohnort == "Horn" & PLZWmiss == 1
replace PLZ4= 8181 if wohnort == "Höri" & PLZWmiss == 1
replace PLZ4= 8553 if wohnort == "Hüttlingen" & PLZWmiss == 1
replace PLZ4= 8307 if wohnort == "Illnau-Effretikon" & PLZWmiss == 1
replace PLZ4= 6440 if wohnort == "Ingenbohl" & PLZWmiss == 1
replace PLZ4= 3232 if wohnort == "Ins" & PLZWmiss == 1
replace PLZ4= 6993 if wohnort == "Iseo" & PLZWmiss == 1
replace PLZ4= 8640 if wohnort == "Jona" & PLZWmiss == 1
replace PLZ4= 1083 if wohnort == "Jorat-Mézières" & PLZWmiss == 1
replace PLZ4= 1254 if wohnort == "Jussy" & PLZWmiss == 1
replace PLZ4= 5466 if wohnort == "Kaiserstuhl" & PLZWmiss == 1
replace PLZ4= 5083 if wohnort == "Kaisten" & PLZWmiss == 1
replace PLZ4= 3273 if wohnort == "Kappelen" & PLZWmiss == 1
replace PLZ4= 3126 if wohnort == "Kaufdorf" & PLZWmiss == 1
replace PLZ4= 4468 if wohnort == "Kienberg" & PLZWmiss == 1
replace PLZ4= 8956 if wohnort == "Killwangen" & PLZWmiss == 1
replace PLZ4= 3038 if wohnort == "Kirchlindach" & PLZWmiss == 1
replace PLZ4= 3213 if wohnort == "Kleinbösingen" & PLZWmiss == 1
replace PLZ4= 5313 if wohnort == "Klingnau" & PLZWmiss == 1
replace PLZ4= 5322 if wohnort == "Koblenz" & PLZWmiss == 1
replace PLZ4= 6010 if wohnort == "Kriens" & PLZWmiss == 1
replace PLZ4= 8314 if wohnort == "Kyburg" & PLZWmiss == 1
replace PLZ4= 6403 if wohnort == "Küssnacht am Rigi" & PLZWmiss == 1
replace PLZ4= 1148 if wohnort == "L'Isle" & PLZWmiss == 1
replace PLZ4= 2300 if wohnort == "La Chaux-de-Fonds" & PLZWmiss == 1
replace PLZ4= 2405 if wohnort == "La Chaux-du-Milieu" & PLZWmiss == 1
replace PLZ4= 1278 if wohnort == "La Rippe" & PLZWmiss == 1
replace PLZ4= 1634 if wohnort == "La Roche" & PLZWmiss == 1
replace PLZ4= 2314 if wohnort == "La Sagne" & PLZWmiss == 1
replace PLZ4= 1814 if wohnort == "La Tour-de-Peilz" & PLZWmiss == 1
replace PLZ4= 1635 if wohnort == "La Tour-de-Trême" & PLZWmiss == 1
replace PLZ4= 2075 if wohnort == "La Tène" & PLZWmiss == 1
replace PLZ4= 7031 if wohnort == "Laax" & PLZWmiss == 1
replace PLZ4= 8853 if wohnort == "Lachen" & PLZWmiss == 1
replace PLZ4= 6814 if wohnort == "Lamone" & PLZWmiss == 1
replace PLZ4= 1200 if wohnort == "Lancy" & PLZWmiss == 1
replace PLZ4= 7302 if wohnort == "Landquart" & PLZWmiss == 1
replace PLZ4= 4242 if wohnort == "Laufen" & PLZWmiss == 1
replace PLZ4= 8212 if wohnort == "Laufen-Uhwiesen" & PLZWmiss == 1
replace PLZ4= 5085 if wohnort == "Laufenburg" & PLZWmiss == 1
replace PLZ4= 3177 if wohnort == "Laupen" & PLZWmiss == 1
replace PLZ4= 1000 if wohnort == "Lausanne" & PLZWmiss == 1
replace PLZ4= 1689 if wohnort == "Le Châtelard" & PLZWmiss == 1
replace PLZ4= 1052 if wohnort == "Le Mont-sur-Lausanne" & PLZWmiss == 1
replace PLZ4= 1261 if wohnort == "Le Vaud" & PLZWmiss == 1
replace PLZ4= 8574 if wohnort == "Lengwil" & PLZWmiss == 1
replace PLZ4= 3775 if wohnort == "Lenk" & PLZWmiss == 1
replace PLZ4= 1607 if wohnort == "Les Thioleyres" & PLZWmiss == 1
replace PLZ4= 3952 if wohnort == "Leuk" & PLZWmiss == 1
replace PLZ4= 8310 if wohnort == "Lindau" & PLZWmiss == 1
replace PLZ4= 3673 if wohnort == "Linden" & PLZWmiss == 1
replace PLZ4= 4514 if wohnort == "Lommiswil" & PLZWmiss == 1
replace PLZ4= 6962 if wohnort == "Lugano" & PLZWmiss == 1
replace PLZ4= 1184 if wohnort == "Luins" & PLZWmiss == 1
replace PLZ4= 7147 if wohnort == "Lumnezia" & PLZWmiss == 1
replace PLZ4= 4542 if wohnort == "Luterbach" & PLZWmiss == 1
replace PLZ4= 9405 if wohnort == "Lutzenberg" & PLZWmiss == 1
replace PLZ4= 6000 if wohnort == "Luzern" & PLZWmiss == 1
replace PLZ4= 8224 if wohnort == "Löhningen" & PLZWmiss == 1
replace PLZ4= 4574 if wohnort == "Lüsslingen" & PLZWmiss == 1
replace PLZ4= 4574 if wohnort == "Lüsslingen-Nennigkofen" & PLZWmiss == 1
replace PLZ4= 3435 if wohnort == "Lützelflüh" & PLZWmiss == 1
replace PLZ4= 4464 if wohnort == "Maisprach" & PLZWmiss == 1
replace PLZ4= 7208 if wohnort == "Malans" & PLZWmiss == 1
replace PLZ4= 6928 if wohnort == "Manno" & PLZWmiss == 1
replace PLZ4= 1920 if wohnort == "Martigny" & PLZWmiss == 1
replace PLZ4= 3800 if wohnort == "Matten bei Interlaken" & PLZWmiss == 1
replace PLZ4= 4713 if wohnort == "Matzendorf" & PLZWmiss == 1
replace PLZ4= 8122 if wohnort == "Maur" & PLZWmiss == 1
replace PLZ4= 6045 if wohnort == "Meggen" & PLZWmiss == 1
replace PLZ4= 5274 if wohnort == "Mettauertal" & PLZWmiss == 1
replace PLZ4= 4115 if wohnort == "Metzerlen-Mariastein" & PLZWmiss == 1
replace PLZ4= 2013 if wohnort == "Milvignes" & PLZWmiss == 1
replace PLZ4= 8756 if wohnort == "Mitlödi" & PLZWmiss == 1
replace PLZ4= 9123 if wohnort == "Mogelsberg" & PLZWmiss == 1
replace PLZ4= 1148 if wohnort == "Moiry" & PLZWmiss == 1
replace PLZ4= 1973 if wohnort == "Mont-Noble" & PLZWmiss == 1
replace PLZ4= 1786 if wohnort == "Mont-Vully" & PLZWmiss == 1
replace PLZ4= 6926 if wohnort == "Montagnola" & PLZWmiss == 1
replace PLZ4= 3963 if wohnort == "Montana" & PLZWmiss == 1
replace PLZ4= 1515 if wohnort == "Montanaire" & PLZWmiss == 1
replace PLZ4= 6802 if wohnort == "Monteceneri" & PLZWmiss == 1
replace PLZ4= 6996 if wohnort == "Monteggio" & PLZWmiss == 1
replace PLZ4= 2362 if wohnort == "Montfaucon" & PLZWmiss == 1
replace PLZ4= 1870 if wohnort == "Monthey" & PLZWmiss == 1
replace PLZ4= 1824 if wohnort == "Montreux" & PLZWmiss == 1
replace PLZ4= 6295 if wohnort == "Mosen" & PLZWmiss == 1
replace PLZ4= 1510 if wohnort == "Moudon" & PLZWmiss == 1
replace PLZ4= 2740 if wohnort == "Moutier" & PLZWmiss == 1
replace PLZ4= 3074 if wohnort == "Muri bei Bern" & PLZWmiss == 1
replace PLZ4= 5103 if wohnort == "Möriken-Wildegg" & PLZWmiss == 1
replace PLZ4= 9402 if wohnort == "Mörschwil" & PLZWmiss == 1
replace PLZ4= 8555 if wohnort == "Müllheim" & PLZWmiss == 1
replace PLZ4= 4719 if wohnort == "Mümliswil-Ramiswil" & PLZWmiss == 1
replace PLZ4= 3110 if wohnort == "Münsingen" & PLZWmiss == 1
replace PLZ4= 2008 if wohnort == "Neuchâtel" & PLZWmiss == 1
replace PLZ4= 2560 if wohnort == "Nidau" & PLZWmiss == 1
replace PLZ4= 8172 if wohnort == "Niederglatt" & PLZWmiss == 1
replace PLZ4= 3853 if wohnort == "Niederried bei Interlaken" & PLZWmiss == 1
replace PLZ4= 4421 if wohnort == "Nuglar-St. Pantaleon" & PLZWmiss == 1
replace PLZ4= 1260 if wohnort == "Nyon" & PLZWmiss == 1
replace PLZ4= 4625 if wohnort == "Oberbuchsiten" & PLZWmiss == 1
replace PLZ4= 9245 if wohnort == "Oberbüren" & PLZWmiss == 1
replace PLZ4= 5420 if wohnort == "Oberehrendingen" & PLZWmiss == 1
replace PLZ4= 8154 if wohnort == "Oberglatt" & PLZWmiss == 1
replace PLZ4= 5062 if wohnort == "Oberhof" & PLZWmiss == 1
replace PLZ4= 6208 if wohnort == "Oberkirch" & PLZWmiss == 1
replace PLZ4= 8942 if wohnort == "Oberrieden" & PLZWmiss == 1
replace PLZ4= 7134 if wohnort == "Obersaxen" & PLZWmiss == 1
replace PLZ4= 7138 if wohnort == "Obersaxen Mundaun" & PLZWmiss == 1
replace PLZ4= 5415 if wohnort == "Obersiggenthal" & PLZWmiss == 1
replace PLZ4= 9242 if wohnort == "Oberuzwil" & PLZWmiss == 1
replace PLZ4= 3114 if wohnort == "Oberwichtrach" & PLZWmiss == 1
replace PLZ4= 3298 if wohnort == "Oberwil bei Büren" & PLZWmiss == 1
replace PLZ4= 1867 if wohnort == "Ollon" & PLZWmiss == 1
replace PLZ4= 4305 if wohnort == "Olsberg" & PLZWmiss == 1
replace PLZ4= 6664 if wohnort == "Onsernone" & PLZWmiss == 1
replace PLZ4= 1430 if wohnort == "Orges" & PLZWmiss == 1
replace PLZ4= 1865 if wohnort == "Ormont-Dessus" & PLZWmiss == 1
replace PLZ4= 1607 if wohnort == "Oron" & PLZWmiss == 1
replace PLZ4= 8913 if wohnort == "Ottenbach" & PLZWmiss == 1
replace PLZ4= 6902 if wohnort == "Paradiso" & PLZWmiss == 1
replace PLZ4= 1303 if wohnort == "Penthaz" & PLZWmiss == 1
replace PLZ4= 1258 if wohnort == "Perly-Certoux" & PLZWmiss == 1
replace PLZ4= 2748 if wohnort == "Petit-Val" & PLZWmiss == 1
replace PLZ4= 5735 if wohnort == "Pfeffikon" & PLZWmiss == 1
replace PLZ4= 4148 if wohnort == "Pfeffingen" & PLZWmiss == 1
replace PLZ4= 8330 if wohnort == "Pfäffikon" & PLZWmiss == 1
replace PLZ4= 2515 if wohnort == "Plateau de Diesse" & PLZWmiss == 1
replace PLZ4= 3638 if wohnort == "Pohlern" & PLZWmiss == 1
replace PLZ4= 6988 if wohnort == "Ponte Tresa" & PLZWmiss == 1
replace PLZ4= 1897 if wohnort == "Port-Valais" & PLZWmiss == 1
replace PLZ4= 1197 if wohnort == "Prangins" & PLZWmiss == 1
replace PLZ4= 1724 if wohnort == "Praroman" & PLZWmiss == 1
replace PLZ4= 1211 if wohnort == "Pregny-Chambésy" & PLZWmiss == 1
replace PLZ4= 1008 if wohnort == "Prilly" & PLZWmiss == 1
replace PLZ4= 1428 if wohnort == "Provence" & PLZWmiss == 1
replace PLZ4= 1028 if wohnort == "Préverenges" & PLZWmiss == 1
replace PLZ4= 1009 if wohnort == "Pully" & PLZWmiss == 1
replace PLZ4= 8640 if wohnort == "Rapperswil (SG)" & PLZWmiss == 1
replace PLZ4= 8640 if wohnort == "Rapperswil-Jona" & PLZWmiss == 1
replace PLZ4= 8462 if wohnort == "Rheinau" & PLZWmiss == 1
replace PLZ4= 4310 if wohnort == "Rheinfelden" & PLZWmiss == 1
replace PLZ4= 9532 if wohnort == "Rickenbach bei Wil" & PLZWmiss == 1
replace PLZ4= 8739 if wohnort == "Rieden" & PLZWmiss == 1
replace PLZ4= 6703 if wohnort == "Riviera" & PLZWmiss == 1
replace PLZ4= 2019 if wohnort == "Rochefort" & PLZWmiss == 1
replace PLZ4= 4938 if wohnort == "Rohrbach" & PLZWmiss == 1
replace PLZ4= 1032 if wohnort == "Romanel-sur-Lausanne" & PLZWmiss == 1
replace PLZ4= 6547 if wohnort == "Rossa" & PLZWmiss == 1
replace PLZ4= 6000 if wohnort == "Rothenburg" & PLZWmiss == 1
replace PLZ4= 8964 if wohnort == "Rudolfstetten-Friedlisberg" & PLZWmiss == 1
replace PLZ4= 3422 if wohnort == "Rüdtligen-Alchenflüh" & PLZWmiss == 1
replace PLZ4= 5235 if wohnort == "Rüfenach" & PLZWmiss == 1
replace PLZ4= 5464 if wohnort == "Rümikon" & PLZWmiss == 1
replace PLZ4= 3154 if wohnort == "Rüschegg" & PLZWmiss == 1
replace PLZ4= 3295 if wohnort == "Rüti bei Büren" & PLZWmiss == 1
replace PLZ4= 7247 if wohnort == "Saas" & PLZWmiss == 1
replace PLZ4= 7104 if wohnort == "Safiental" & PLZWmiss == 1
replace PLZ4= 2072 if wohnort == "Saint-Blaise" & PLZWmiss == 1
replace PLZ4= 1265 if wohnort == "Saint-Cergue" & PLZWmiss == 1
replace PLZ4= 1898 if wohnort == "Saint-Gingolph" & PLZWmiss == 1
replace PLZ4= 3961 if wohnort == "Saint-Jean" & PLZWmiss == 1
replace PLZ4= 1890 if wohnort == "Saint-Maurice" & PLZWmiss == 1
replace PLZ4= 1162 if wohnort == "Saint-Prex" & PLZWmiss == 1
replace PLZ4= 1450 if wohnort == "Sainte-Croix" & PLZWmiss == 1
replace PLZ4= 7563 if wohnort == "Samnaun" & PLZWmiss == 1
replace PLZ4= 6575 if wohnort == "San Nazzaro" & PLZWmiss == 1
replace PLZ4= 6534 if wohnort == "San Vittore" & PLZWmiss == 1
replace PLZ4= 6577 if wohnort == "Sant'Abbondio" & PLZWmiss == 1
replace PLZ4= 7536 if wohnort == "Santa Maria im Münstertal" & PLZWmiss == 1
replace PLZ4= 6417 if wohnort == "Sattel" & PLZWmiss == 1
replace PLZ4= 2537 if wohnort == "Sauge" & PLZWmiss == 1
replace PLZ4= 1073 if wohnort == "Savigny" & PLZWmiss == 1
replace PLZ4= 8208 if wohnort == "Schaffhausen" & PLZWmiss == 1
replace PLZ4= 5107 if wohnort == "Schinznach" & PLZWmiss == 1
replace PLZ4= 5107 if wohnort == "Schinznach Dorf" & PLZWmiss == 1
replace PLZ4= 8418 if wohnort == "Schlatt" & PLZWmiss == 1
replace PLZ4= 6231 if wohnort == "Schlierbach" & PLZWmiss == 1
replace PLZ4= 3855 if wohnort == "Schwanden bei Brienz" & PLZWmiss == 1
replace PLZ4= 6215 if wohnort == "Schwarzenbach" & PLZWmiss == 1
replace PLZ4= 6103 if wohnort == "Schwarzenberg" & PLZWmiss == 1
replace PLZ4= 3150 if wohnort == "Schwarzenburg" & PLZWmiss == 1
replace PLZ4= 8762 if wohnort == "Schwändi" & PLZWmiss == 1
replace PLZ4= 8577 if wohnort == "Schönholzerswilen" & PLZWmiss == 1
replace PLZ4= 7545 if wohnort == "Scuol/Schuls" & PLZWmiss == 1
replace PLZ4= 8607 if wohnort == "Seegräben" & PLZWmiss == 1
replace PLZ4= 6377 if wohnort == "Seelisberg" & PLZWmiss == 1
replace PLZ4= 4206 if wohnort == "Seewen" & PLZWmiss == 1
replace PLZ4= 3662 if wohnort == "Seftigen" & PLZWmiss == 1
replace PLZ4= 1623 if wohnort == "Semsales" & PLZWmiss == 1
replace PLZ4= 9468 if wohnort == "Sennwald" & PLZWmiss == 1
replace PLZ4= 6713 if wohnort == "Serravalle" & PLZWmiss == 1
replace PLZ4= 7514 if wohnort == "Sils im Engadin/Segl" & PLZWmiss == 1
replace PLZ4= 2577 if wohnort == "Siselen" & PLZWmiss == 1
replace PLZ4= 4334 if wohnort == "Sisseln" & PLZWmiss == 1
replace PLZ4= 4503 if wohnort == "Solothurn" & PLZWmiss == 1
replace PLZ4= 9020 if wohnort == "St. Gallen" & PLZWmiss == 1
replace PLZ4= 9430 if wohnort == "St. Margrethen" & PLZWmiss == 1
replace PLZ4= 7116 if wohnort == "St. Martin" & PLZWmiss == 1
replace PLZ4= 3927 if wohnort == "St. Niklaus" & PLZWmiss == 1
replace PLZ4= 7028 if wohnort == "St. Peter" & PLZWmiss == 1
replace PLZ4= 1736 if wohnort == "St. Silvester" & PLZWmiss == 1
replace PLZ4= 8175 if wohnort == "Stadel" & PLZWmiss == 1
replace PLZ4= 6370 if wohnort == "Stans" & PLZWmiss == 1
replace PLZ4= 5603 if wohnort == "Staufen" & PLZWmiss == 1
replace PLZ4= 3940 if wohnort == "Steg" & PLZWmiss == 1
replace PLZ4= 9323 if wohnort == "Steinach" & PLZWmiss == 1
replace PLZ4= 6422 if wohnort == "Steinen" & PLZWmiss == 1
replace PLZ4= 6312 if wohnort == "Steinhausen" & PLZWmiss == 1
replace PLZ4= 3632 if wohnort == "Stocken-Höfen" & PLZWmiss == 1
replace PLZ4= 2557 if wohnort == "Studen" & PLZWmiss == 1
replace PLZ4= 1036 if wohnort == "Sullens" & PLZWmiss == 1
replace PLZ4= 7452 if wohnort == "Surses" & PLZWmiss == 1
replace PLZ4= 2572 if wohnort == "Sutz-Lattrigen" & PLZWmiss == 1
replace PLZ4= 1554 if wohnort == "Sédeilles" & PLZWmiss == 1
replace PLZ4= 2720 if wohnort == "Tavannes" & PLZWmiss == 1
replace PLZ4= 6654 if wohnort == "Terre di Pedemonte" & PLZWmiss == 1
replace PLZ4= 4106 if wohnort == "Therwil" & PLZWmiss == 1
replace PLZ4= 2075 if wohnort == "Thielle-Wavre" & PLZWmiss == 1
replace PLZ4= 3600 if wohnort == "Thun" & PLZWmiss == 1
replace PLZ4= 6808 if wohnort == "Torricella-Taverne" & PLZWmiss == 1
replace PLZ4= 8219 if wohnort == "Trasadingen" & PLZWmiss == 1
replace PLZ4= 9043 if wohnort == "Trogen" & PLZWmiss == 1
replace PLZ4= 1872 if wohnort == "Troistorrents" & PLZWmiss == 1
replace PLZ4= 3946 if wohnort == "Turtmann-Unterems" & PLZWmiss == 1
replace PLZ4= 2513 if wohnort == "Twann-Tüscherz" & PLZWmiss == 1
replace PLZ4= 1423 if wohnort == "Tévenon" & PLZWmiss == 1
replace PLZ4= 8707 if wohnort == "Uetikon" & PLZWmiss == 1
replace PLZ4= 8142 if wohnort == "Uitikon" & PLZWmiss == 1
replace PLZ4= 3944 if wohnort == "Unterbäch" & PLZWmiss == 1
replace PLZ4= 5726 if wohnort == "Unterkulm" & PLZWmiss == 1
replace PLZ4= 8476 if wohnort == "Unterstammheim" & PLZWmiss == 1
replace PLZ4= 3322 if wohnort == "Urtenen" & PLZWmiss == 1
replace PLZ4= 8592 if wohnort == "Uttwil" & PLZWmiss == 1
replace PLZ4= 9240 if wohnort == "Uzwil" & PLZWmiss == 1
replace PLZ4= 6833 if wohnort == "Vacallo" & PLZWmiss == 1
replace PLZ4= 2829 if wohnort == "Val Terbi" & PLZWmiss == 1
replace PLZ4= 1873 if wohnort == "Val-d'Illiez" & PLZWmiss == 1
replace PLZ4= 1654 if wohnort == "Val-de-Charmey" & PLZWmiss == 1
replace PLZ4= 2057 if wohnort == "Val-de-Ruz" & PLZWmiss == 1
replace PLZ4= 2114 if wohnort == "Val-de-Travers" & PLZWmiss == 1
replace PLZ4= 2735 if wohnort == "Valbirse" & PLZWmiss == 1
replace PLZ4= 1536 if wohnort == "Valbroye" & PLZWmiss == 1
replace PLZ4= 6951 if wohnort == "Valcolla" & PLZWmiss == 1
replace PLZ4= 1337 if wohnort == "Vallorbe" & PLZWmiss == 1
replace PLZ4= 7558 if wohnort == "Valsot" & PLZWmiss == 1
replace PLZ4= 1219 if wohnort == "Vernier" & PLZWmiss == 1
replace PLZ4= 1290 if wohnort == "Versoix" & PLZWmiss == 1
replace PLZ4= 1255 if wohnort == "Veyrier" & PLZWmiss == 1
replace PLZ4= 1752 if wohnort == "Villars-sur-Glâne" & PLZWmiss == 1
replace PLZ4= 1091 if wohnort == "Villette (Lavaux)" & PLZWmiss == 1
replace PLZ4= 2057 if wohnort == "Villiers" & PLZWmiss == 1
replace PLZ4= 8604 if wohnort == "Volketswil" & PLZWmiss == 1
replace PLZ4= 7149 if wohnort == "Vrin" & PLZWmiss == 1
replace PLZ4= 1418 if wohnort == "Vuarrens" & PLZWmiss == 1
replace PLZ4= 1509 if wohnort == "Vucherens" & PLZWmiss == 1
replace PLZ4= 1134 if wohnort == "Vufflens-le-Château" & PLZWmiss == 1
replace PLZ4= 1696 if wohnort == "Vuisternens-en-Ogoz" & PLZWmiss == 1
replace PLZ4= 1585 if wohnort == "Vully-les-Lacs" & PLZWmiss == 1
replace PLZ4= 3152 if wohnort == "Wahlern" & PLZWmiss == 1
replace PLZ4= 4437 if wohnort == "Waldenburg" & PLZWmiss == 1
replace PLZ4= 9205 if wohnort == "Waldkirch" & PLZWmiss == 1
replace PLZ4= 4323 if wohnort == "Wallbach" & PLZWmiss == 1
replace PLZ4= 3380 if wohnort == "Walliswil bei Niederbipp" & PLZWmiss == 1
replace PLZ4= 3377 if wohnort == "Walliswil bei Wangen" & PLZWmiss == 1
replace PLZ4= 4612 if wohnort == "Wangen bei Olten" & PLZWmiss == 1
replace PLZ4= 9479 if wohnort == "Wartau" & PLZWmiss == 1
replace PLZ4= 6485 if wohnort == "Wassen" & PLZWmiss == 1
replace PLZ4= 3251 if wohnort == "Wengi" & PLZWmiss == 1
replace PLZ4= 8907 if wohnort == "Wettswil am Albis" & PLZWmiss == 1
replace PLZ4= 8489 if wohnort == "Wildberg" & PLZWmiss == 1
replace PLZ4= 3428 if wohnort == "Wiler bei Utzenstorf" & PLZWmiss == 1
replace PLZ4= 8185 if wohnort == "Winkel" & PLZWmiss == 1
replace PLZ4= 8400 if wohnort == "Winterthur" & PLZWmiss == 1
replace PLZ4= 5064 if wohnort == "Wittnau" & PLZWmiss == 1
replace PLZ4= 3043 if wohnort == "Wohlen bei Bern" & PLZWmiss == 1
replace PLZ4= 3077 if wohnort == "Worb" & PLZWmiss == 1
replace PLZ4= 9175 if wohnort == "Wünnewil-Flamatt" & PLZWmiss == 1
replace PLZ4= 1169 if wohnort == "Yens" & PLZWmiss == 1
replace PLZ4= 1432 if wohnort == "Yverdon-les-Bains" & PLZWmiss == 1
replace PLZ4= 1853 if wohnort == "Yvorne" & PLZWmiss == 1
replace PLZ4= 6301 if wohnort == "Zug" & PLZWmiss == 1
replace PLZ4= 5330 if wohnort == "Zurzach" & PLZWmiss == 1
replace PLZ4= 8004 if wohnort == "Zürich" & PLZWmiss == 1

/*// (still) Missing Wohnort PLZ & Wohnland=="Schweiz"
gen c_persID = 1
bysort personenid: egen obspersID = total(c_persID)
gen c_missPLZ = 0
replace c_missPLZ = 1 if PLZ4 == .
bysort personenid: egen obsmissPLZ = total(c_missPLZ)
*br * if obsmissPLZ < obspersID & obsmissPLZ != 0  // are there observations within Person ID that have non-missing PLZ -- > 0 Obs.
drop c_persID c_missPLZ obspersID obsmissPLZ

preserve
keep if wohnland == "Schweiz" & PLZWmiss==1 & PLZ4==.  // Obs. with Swiss "Wohnland" which are still missing.
gen count = 1
collapse (sum) count, by(wohnort)
sort wohnort
save "$path\02_Processed_data\11_Directors_1994_2018/Bisnode_missingWohnPLZ_temp2.dta", replace
restore
*/

// Provide PLZ options for those observations, where PLZ4 is missing, Wohnland=="Schweiz", and No. of observations per Wohnort >= 4
gen PLZ4_1 =.
replace PLZ4_1 = 8883 if wohnort == "Murg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6516 if wohnort == "Gerra" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6950 if wohnort == "Gola di Lago" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2543 if wohnort == "Lengnau" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9428 if wohnort == "Lachen AR" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5628 if wohnort == "Birri" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6418 if wohnort == "Rothenturn" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8508 if wohnort == "Unterhörstetten" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8910 if wohnort == "Affoltern" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3203 if wohnort == "Allenlüften" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3203 if wohnort == "Mauss" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4616 if wohnort == "Kappel" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8640 if wohnort == "Rapperswi-Jona" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4665 if wohnort == "Küngoldingen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5212 if wohnort == "Hausen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1997 if wohnort == "Sornard" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5524 if wohnort == "Niederwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8733 if wohnort == "Eschenbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3255 if wohnort == "Rapperswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8904 if wohnort == "Aesch bei Birmensdorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3047 if wohnort == "Stuckishaus" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8614 if wohnort == "Gossau" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3985 if wohnort == "Münster" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2022 if wohnort == "La Grande-Béroche" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3636 if wohnort == "Forst bei Langenbühl" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8703 if wohnort == "Erlenbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1218 if wohnort == "Grand-Saconnex" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1787 if wohnort == "Mur" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8854 if wohnort == "Siebnen-Galgenen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8370 if wohnort == "Wiezikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7536 if wohnort == "Santa Maria" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4462 if wohnort == "Rickenbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4153 if wohnort == "Reinach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1946 if wohnort == "Bourg-St-Bernard" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3203 if wohnort == "Spengelried" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5015 if wohnort == "Obererlinsbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3114 if wohnort == "Niederwichtrach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5085 if wohnort == "Sulz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8802 if wohnort == "Kilchberg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8305 if wohnort == "Dietilikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8717 if wohnort == "Benken" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3860 if wohnort == "Willigen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3065 if wohnort == "Habstetten" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6433 if wohnort == "Stoos" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3963 if wohnort == "Crans" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7082 if wohnort == "Lain" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6981 if wohnort == "Biogno" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6106 if wohnort == "Schachen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3365 if wohnort == "Niedergrasswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6535 if wohnort == "Roveredo" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5415 if wohnort == "Nussbaumen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3035 if wohnort == "Aspi b. Seedorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1724 if wohnort == "Essert" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5620 if wohnort == "Bremgarten" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3018 if wohnort == "Bümpliz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7494 if wohnort == "Davos Wiesen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5736 if wohnort == "Burg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1941 if wohnort == "Chemin-Dessus" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9127 if wohnort == "Wald bei St. Peterzell" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3074 if wohnort == "Melchenbühl" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9100 if wohnort == "Hersiau" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8525 if wohnort == "Kefikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6317 if wohnort == "Oberwil bei Zug" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8806 if wohnort == "Bäch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9533 if wohnort == "Dietschwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9053 if wohnort == "Teufen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3065 if wohnort == "Geristein" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3095 if wohnort == "Spiegel bei Bern" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2555 if wohnort == "Brügg bei Biel" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4114 if wohnort == "Hofstetten" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8630 if wohnort == "Rüti" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6982 if wohnort == "Serocca" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9322 if wohnort == "Neukirch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3508 if wohnort == "Arni bei Biglen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9535 if wohnort == "Wilen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1614 if wohnort == "Granges" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6130 if wohnort == "Willisau-Stadt" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1227 if wohnort == "Carouge" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1897 if wohnort == "Le Bouveret" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9463 if wohnort == "Kobelwald" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8514 if wohnort == "Bänikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4105 if wohnort == "Benken BL" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4147 if wohnort == "Aesch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8050 if wohnort == "Oerlikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3294 if wohnort == "Büren" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6962 if wohnort == "Albonago" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3960 if wohnort == "Muraz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3775 if wohnort == "Lenk im Simmenthal" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9411 if wohnort == "Mohren" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3508 if wohnort == "Arni" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9533 if wohnort == "Kirchberg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3098 if wohnort == "Schliern bei Köniz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6196 if wohnort == "Welsikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6950 if wohnort == "Lelgio" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9656 if wohnort == "Wildhaus- Alt St. Johann" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3038 if wohnort == "Oberlindach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8903 if wohnort == "Birmensdorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2500 if wohnort == "Biel" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5630 if wohnort == "Muri" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3983 if wohnort == "Breiten VS" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8268 if wohnort == "Mannenbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8332 if wohnort == "Rumlikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8362 if wohnort == "Ifwil bei Balterswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9466 if wohnort == "Haag" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9470 if wohnort == "BuchsSG" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8352 if wohnort == "Schottikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3203 if wohnort == "Buttenried" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1092 if wohnort == "Belmont" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3045 if wohnort == "Weissenstein b. Meikirch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3267 if wohnort == "Aspi bei Seedorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5112 if wohnort == "Thalheim" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5113 if wohnort == "Holderbank" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4332 if wohnort == "Stein" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9542 if wohnort == "Oberhofen TG" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3116 if wohnort == "Kirchdorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3033 if wohnort == "Innerberg bei Säriswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6515 if wohnort == "Progero" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1964 if wohnort == "Conthey-Bourg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8514 if wohnort == "Bissegg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1740 if wohnort == "Neyruz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3203 if wohnort == "Juchlishaus" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8306 if wohnort == "Wangen-Brüttisellen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8620 if wohnort == "Wetzikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7278 if wohnort == "Monstein" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5420 if wohnort == "Unterehrendingen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3924 if wohnort == "Mattsand" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3918 if wohnort == "Wiler" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8046 if wohnort == "Affoltern b. Zürich" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3854 if wohnort == "Oberried" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8000 if wohnort == "Zürichhorn" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3754 if wohnort == "Entschwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4438 if wohnort == "Bärenwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1787 if wohnort == "Haut-Vully" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3075 if wohnort == "Vielbringen b. Worb" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8452 if wohnort == "Niederwil ZH" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5413 if wohnort == "Birmenstorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6052 if wohnort == "Hergiswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2022 if wohnort == "La Grande Béroche" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6460 if wohnort == "Altdorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8700 if wohnort == "Küsnacht" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3653 if wohnort == "Oberhofen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6436 if wohnort == "Ried SZ" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8824 if wohnort == "Schönenberg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8309 if wohnort == "Birchwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6513 if wohnort == "Carasso" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5428 if wohnort == "Schönbühl" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3900 if wohnort == "Briga" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1814 if wohnort == "La Tour-de- Peilz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8580 if wohnort == "Auenhofen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1253 if wohnort == "Vandoeuvre" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3292 if wohnort == "Busswil BE" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9430 if wohnort == "St. Margarethen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9422 if wohnort == "Staad" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4914 if wohnort == "Roggwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8486 if wohnort == "Rikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6535 if wohnort == "Roverdo GR" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5610 if wohnort == "Wohlen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4515 if wohnort == "Oberdorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1246 if wohnort == "Corsier" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8107 if wohnort == "Buchs" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6516 if wohnort == "Gerra Piano" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8635 if wohnort == "Oberdürnten" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8580 if wohnort == "Hatswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3852 if wohnort == "Ringgenberg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9533 if wohnort == "Wolfikon SG" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4853 if wohnort == "Balzenwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5332 if wohnort == "Rekingen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1700 if wohnort == "Villars-les-Joncs" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8104 if wohnort == "Weiningen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6574 if wohnort == "Vira Gambarogno" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5417 if wohnort == "Ennetturgi" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6781 if wohnort == "Ossasco" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5733 if wohnort == "Leimbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1024 if wohnort == "Ecublens" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5412 if wohnort == "Vogelsang AG" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1741 if wohnort == "Cottens" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3145 if wohnort == "Oberscherli" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8252 if wohnort == "Schlatt bei Diessenhofen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3047 if wohnort == "Münchwilen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8340 if wohnort == "Hadlikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5723 if wohnort == "Teufenthal" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4104 if wohnort == "Oberwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9322 if wohnort == "Hegi" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5608 if wohnort == "Stetten" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6900 if wohnort == "Cassarate" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9437 if wohnort == "Marbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8524 if wohnort == "Horben b. Frauenfeld" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8362 if wohnort == "Ifwil b. Balterswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8737 if wohnort == "Ernetschwil Gommiswald" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8315 if wohnort == "Kemptthal-Grafstal" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1950 if wohnort == "Molignon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8252 if wohnort == "Mett-Oberschlatt" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8357 if wohnort == "Maischhausen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3766 if wohnort == "Reidenbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8558 if wohnort == "Helsighausen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1934 if wohnort == "Fontenelle" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8498 if wohnort == "Gibswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9562 if wohnort == "Weingarten b. Märwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6600 if wohnort == "Locarno Monti" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8640 if wohnort == "Jona-Rapperswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1806 if wohnort == "St.Légier-La Chiésaz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1731 if wohnort == "Ependes" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8890 if wohnort == "Flumserberg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3065 if wohnort == "Flugbrunnen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2572 if wohnort == "Lattrigen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1680 if wohnort == "Romont" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6802 if wohnort == "Monte Ceneri" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8525 if wohnort == "Wilen-Neunforn" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3855 if wohnort == "Schwanden" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3086 if wohnort == "Wald" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3250 if wohnort == "Busswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3907 if wohnort == "Simplon Kulm" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5453 if wohnort == "Busslingen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7522 if wohnort == "La Punt Chamues-ch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8500 if wohnort == "Erzenholz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8964 if wohnort == "Friedlisberg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7433 if wohnort == "Patzen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2336 if wohnort == "Le Peu-Claude" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7017 if wohnort == "Flims-Dorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8570 if wohnort == "Weerswilen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8954 if wohnort == "Geroldwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8132 if wohnort == "Egg ZH" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8152 if wohnort == "Optikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8890 if wohnort == "Grossberg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2000 if wohnort == "Neuenburg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7260 if wohnort == "Frauenkirch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3186 if wohnort == "Guin" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8472 if wohnort == "Oberohringen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1869 if wohnort == "Daviaz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3172 if wohnort == "Niederwangen bei Bern" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6950 if wohnort == "Tessserete" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3538 if wohnort == "Emmental" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1290 if wohnort == "Chavannes-de-Bois" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2720 if wohnort == "Le Reussilles" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3806 if wohnort == "Bönigen bei Interlaken" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4463 if wohnort == "Buss" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3257 if wohnort == "Kosthofen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5316 if wohnort == "Gippingen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1025 if wohnort == "St-Sulpice" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3945 if wohnort == "Getwing" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8967 if wohnort == "Widen b. Bremgarten" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6005 if wohnort == "St. Niklausen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1782 if wohnort == "Sonnaz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 5106 if wohnort == "Veltheim" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8000 if wohnort == "Zürch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8400 if wohnort == "Sennhof" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1898 if wohnort == "St-Gingolph" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4532 if wohnort == "Freldbrunnen-St. Niklaus" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1700 if wohnort == "Freiburg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9470 if wohnort == "Räfis" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1908 if wohnort == "Tzoumaz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3235 if wohnort == "Cerlier" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8855 if wohnort == "Wangen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8574 if wohnort == "Dettighofen bei Oberhofen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8574 if wohnort == "Dettighofen b. Oberhofen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9204 if wohnort == "Andwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1661 if wohnort == "Le Pâquier" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6988 if wohnort == "Ponte" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6310 if wohnort == "Willisau-Land" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 4500 if wohnort == "Soletta" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1965 if wohnort == "Diolly" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8552 if wohnort == "Wellhausen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7250 if wohnort == "Mezzaselva" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1934 if wohnort == "Montagnier" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9500 if wohnort == "Bronschofen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3854 if wohnort == "Oberried BE" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8606 if wohnort == "Nänikon-Greifensee" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3963 if wohnort == "Bluche" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2053 if wohnort == "Les Cerniers" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8738 if wohnort == "Uetliburg bei Gommiswald" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9434 if wohnort == "Aus" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2500 if wohnort == "Biel/Bienne" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8307 if wohnort == "Luckhausen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2500 if wohnort == "Bienne" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8505 if wohnort == "Dettighofen bei Lanzenneunforn" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1482 if wohnort == "Cugy" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8545 if wohnort == "Sulz ZH" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9422 if wohnort == "Buechen bei Staad" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1273 if wohnort == "Arzier-Le-Muids" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6600 if wohnort == "Locarno-Solduno" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1934 if wohnort == "Châble VS" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8532 if wohnort == "Kartause Ittingen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1908 if wohnort == "La Tsoumaz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8309 if wohnort == "Breite bei Nürensdorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3173 if wohnort == "Oberwangen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3172 if wohnort == "Niederwangen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8966 if wohnort == "Oberwil AG" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1470 if wohnort == "Lully" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1934 if wohnort == "Médières" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3780 if wohnort == "Gruben bei Gstaad" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9500 if wohnort == "Wil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6512 if wohnort == "Giubisco" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8967 if wohnort == "Widen AG" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 9325 if wohnort == "Freidorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1844 if wohnort == "Villeneuve" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6979 if wohnort == "Bré sopra Lugano" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7743 if wohnort == "Zalende" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7063 if wohnort == "Prada" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 7014 if wohnort == "Trin-Mulin" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 6545 if wohnort == "Landarenca" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3257 if wohnort == "Ammerzwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1475 if wohnort == "Vernay FR" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1700 if wohnort == "Frbourg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3088 if wohnort == "Vorderfultigen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3556 if wohnort == "Fankhaus" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 2520 if wohnort == "Schafis" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8487 if wohnort == "Zell" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 1785 if wohnort == "Cressier" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8904 if wohnort == "Aesch b. B." & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 3255 if wohnort == "Seewil bei Dieterswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_1 = 8903 if wohnort == "Landikon" & PLZ4 == . & wohnland == "Schweiz"
gen PLZ4_2 =.
replace PLZ4_2 = 5273 if wohnort == "Oberhofen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 4413 if wohnort == "Büren" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 3422 if wohnort == "Kirchberg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 1673 if wohnort == "Ecublens" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8782 if wohnort == "Rüti" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8832 if wohnort == "Bäch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 9063 if wohnort == "Stein" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8000 if wohnort == "Leimbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 1702 if wohnort == "Oberwil bei Zug" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 4333 if wohnort == "Münchwilen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8302 if wohnort == "Bänikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8542 if wohnort == "Kefikon" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 5426 if wohnort == "Lengnau" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 6950 if wohnort == "Roveredo" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 3047 if wohnort == "Bremgarten" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 1432 if wohnort == "Belmont" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 9200 if wohnort == "Gossau" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 3033 if wohnort == "Wohlen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 6284 if wohnort == "Sulz" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 3762 if wohnort == "Erlenbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 4436 if wohnort == "Oberdorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 6287 if wohnort == "Aesch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8478 if wohnort == "Thalheim" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8966 if wohnort == "Oberwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8537 if wohnort == "Nussbaumen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 1116 if wohnort == "Cottens" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 4523 if wohnort == "Niederwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 6934 if wohnort == "Biogno" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 5734 if wohnort == "Reinach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 9325 if wohnort == "Roggwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 9470 if wohnort == "Buchs" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8463 if wohnort == "Benken" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 3074 if wohnort == "Muri" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 6274 if wohnort == "Eschenbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 6221 if wohnort == "Rickenbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 9100 if wohnort == "Wilen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8640 if wohnort == "Rapperswil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 2316 if wohnort == "Affoltern" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 6482 if wohnort == "Wiler" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 4117 if wohnort == "Burg" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8354 if wohnort == "Hofstetten" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8905 if wohnort == "Arni" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8234 if wohnort == "Stetten" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8926 if wohnort == "Kappel" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 1804 if wohnort == "Corsier" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 4718 if wohnort == "Holderbank" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8915 if wohnort == "Hausen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8240 if wohnort == "Altdorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 2540 if wohnort == "Granges" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8585 if wohnort == "Andwil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 1527 if wohnort == "Villeneuve" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8196 if wohnort == "Wil" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 6144 if wohnort == "Zell" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 8636 if wohnort == "Wald" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 1132 if wohnort == "Lully" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 6064 if wohnort == "St. Niklausen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 2123 if wohnort == "St-Sulpice" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 1434 if wohnort == "Ependes" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 2088 if wohnort == "Cressier" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 2538 if wohnort == "Romont" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_2 = 1053 if wohnort == "Cugy" & PLZ4 == . & wohnland == "Schweiz"
gen PLZ4_3 =.
replace PLZ4_3 = 5272 if wohnort == "Büren" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 3960 if wohnort == "Granges" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 5033 if wohnort == "Buchs" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 4252 if wohnort == "Wiler" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 6370 if wohnort == "Oberdorf" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 6222 if wohnort == "Bäch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 8180 if wohnort == "Nussbaumen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 4613 if wohnort == "Rickenbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 8904 if wohnort == "Aesch" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 6900 if wohnort == "Biogno" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 9655 if wohnort == "Stein" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 9428 if wohnort == "Wilen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 8046 if wohnort == "Affoltern" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 2564 if wohnort == "Belmont" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 9044 if wohnort == "Wald" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_3 = 1233 if wohnort == "Lully" & PLZ4 == . & wohnland == "Schweiz"
gen PLZ4_4 =.
replace PLZ4_4 = 6060 if wohnort == "Wilen" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_4 = 5024 if wohnort == "Büren" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_4 = 8427 if wohnort == "Wiler" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_4 = 1563 if wohnort == "Belmont" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_4 = 8260 if wohnort == "Stein" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_4 = 9532 if wohnort == "Rickenbach" & PLZ4 == . & wohnland == "Schweiz"
gen PLZ4_5 =.
replace PLZ4_5 = 8558 if wohnort == "Büren" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_5 = 8545 if wohnort == "Rickenbach" & PLZ4 == . & wohnland == "Schweiz"
replace PLZ4_5 = 9213 if wohnort == "Wilen" & PLZ4 == . & wohnland == "Schweiz"
gen PLZ4_6 =.
replace PLZ4_6 = 3312 if wohnort == "Büren" & PLZ4 == . & wohnland == "Schweiz"
gen PLZ4_7 =.
replace PLZ4_7 = 6370 if wohnort == "Büren" & PLZ4 == . & wohnland == "Schweiz"
gen PLZ4_8 =.
replace PLZ4_8 = 6386 if wohnort == "Büren" & PLZ4 == . & wohnland == "Schweiz"

rename PLZ4_* hPLZ4_*
gen help1=_n
reshape long hPLZ4_, i(help1) j(obs)
rename hPLZ4_ hPLZ4
drop if hPLZ4 == . & obs > 1
replace PLZ4 = hPLZ4 if PLZ4 == . & hPLZ4 != . 
rename obs wohnort_PLZopt
drop hPLZ4 help1

// We also check observations for which Wohnland != "Schweiz" -- > 
// we found cases where the Wohnort is actually in Switzerland.
// However, the PLZ was not missing and was from the Wohnort in Switzerland 
// We checked the first 2000 Obs (number of obs per wohnort < 30).
// 3003 observations out of a total of 293389 wohnorte (indicated non-Swiss) where corrected
replace wohnland = "Schweiz" if wohnort == "Steinen" & wohnland == "Deutschland" & PLZ4 == 6422
replace wohnland = "Schweiz" if wohnort == "Egg" & wohnland == "" & PLZ4 == 8132
replace wohnland = "Schweiz" if wohnort == "Villars-Epeney" & wohnland == "" & PLZ4 == 1404
replace wohnland = "Schweiz" if wohnort == "Müllheim" & wohnland == "Deutschland" & PLZ4 == 8555
replace wohnland = "Schweiz" if wohnort == "Pfäffikon" & wohnland == "" & PLZ4 == 8330
replace wohnland = "Schweiz" if wohnort == "Montreux" & wohnland == "Frankreich" & PLZ4 == 1824
replace wohnland = "Schweiz" if wohnort == "Binningen" & wohnland == "Deutschland" & PLZ4 == 4102
replace wohnland = "Schweiz" if wohnort == "Rheinfelden" & wohnland == "" & PLZ4 == 4310
replace wohnland = "Schweiz" if wohnort == "Genève" & wohnland == "Frankreich" & PLZ4 == 1200
replace wohnland = "Schweiz" if wohnort == "Lindau" & wohnland == "Deutschland" & PLZ4 == 8315
replace wohnland = "Schweiz" if wohnort == "Saint-Gingolph" & wohnland == "Frankreich" & PLZ4 == 1898
replace wohnland = "Schweiz" if wohnort == "Murg" & wohnland == "Deutschland" & PLZ4 == 8877
replace wohnland = "Schweiz" if wohnort == "Laufenburg" & wohnland == "Deutschland" & PLZ4 == 5080
replace wohnland = "Schweiz" if wohnort == "Zürich" & wohnland == "Deutschland" & PLZ4 == 8000
replace wohnland = "Schweiz" if wohnort == "Bellevue" & wohnland == "Vereinigte Staaten" & PLZ4 == 1293
replace wohnland = "Schweiz" if wohnort == "Lugano" & wohnland == "Italien" & PLZ4 == 6900
replace wohnland = "Schweiz" if wohnort == "Vandoeuvres" & wohnland == "" & PLZ4 == 1253
replace wohnland = "Schweiz" if wohnort == "Pringy" & wohnland == "Frankreich" & PLZ4 == 1663
replace wohnland = "Schweiz" if wohnort == "Niederglatt" & wohnland == "" & PLZ4 == 8172
replace wohnland = "Schweiz" if wohnort == "Burgdorf" & wohnland == "Deutschland" & PLZ4 == 3400
replace wohnland = "Schweiz" if wohnort == "San Nazzaro" & wohnland == "Italien" & PLZ4 == 6575
replace wohnland = "Schweiz" if wohnort == "Fribourg" & wohnland == "" & PLZ4 == 1700
replace wohnland = "Schweiz" if wohnort == "Collonges" & wohnland == "Frankreich" & PLZ4 == 1903
replace wohnland = "Schweiz" if wohnort == "Staufen" & wohnland == "Deutschland" & PLZ4 == 5603
replace wohnland = "Schweiz" if wohnort == "Schlierbach" & wohnland == "Frankreich" & PLZ4 == 6231
replace wohnland = "Schweiz" if wohnort == "Stadel" & wohnland == "" & PLZ4 == 8175
replace wohnland = "Schweiz" if wohnort == "Koblenz" & wohnland == "Deutschland" & PLZ4 == 5322
replace wohnland = "Schweiz" if wohnort == "Meggen" & wohnland == "" & PLZ4 == 6045
replace wohnland = "Schweiz" if wohnort == "Blonay" & wohnland == "Frankreich" & PLZ4 == 1807
replace wohnland = "Schweiz" if wohnort == "Fischingen" & wohnland == "Deutschland" & PLZ4 == 8376
replace wohnland = "Schweiz" if wohnort == "Sant'Abbondio" & wohnland == "Italien" & PLZ4 == 6577
replace wohnland = "Schweiz" if wohnort == "Martigny" & wohnland == "Frankreich" & PLZ4 == 1920
replace wohnland = "Schweiz" if wohnort == "Esslingen" & wohnland == "Deutschland" & PLZ4 == 8133
replace wohnland = "Schweiz" if wohnort == "Meggen" & wohnland == "Deutschland" & PLZ4 == 6045
replace wohnland = "Schweiz" if wohnort == "Müllheim" & wohnland == "" & PLZ4 == 8555
replace wohnland = "Schweiz" if wohnort == "Mörschwil" & wohnland == "Deutschland" & PLZ4 == 9402
replace wohnland = "Schweiz" if wohnort == "Porza" & wohnland == "Österreich" & PLZ4 == 6948
replace wohnland = "Schweiz" if wohnort == "Grandson" & wohnland == "" & PLZ4 == 1422
replace wohnland = "Schweiz" if wohnort == "Rheinau" & wohnland == "Deutschland" & PLZ4 == 8462
replace wohnland = "Schweiz" if wohnort == "Oberrieden" & wohnland == "Deutschland" & PLZ4 == 8942
replace wohnland = "Schweiz" if wohnort == "Epagny" & wohnland == "Frankreich" & PLZ4 == 1663
replace wohnland = "Schweiz" if wohnort == "Monthey" & wohnland == "Vereinigte Staaten" & PLZ4 == 1870

// Non-Swiss observations: weed-out foreign Wohnorte in non-missing PLZ (only those that do not originally merge with Post-data)
gen foreignWO = .   
replace foreignWO = 1 if PLZ4 == 2070 & wohnort == "La Marsa-Tunis"
replace foreignWO = 1 if PLZ4 == 2440 & wohnort == "Geel"
replace foreignWO = 1 if PLZ4 == 2627 & wohnort == "Weissenberg"
replace foreignWO = 1 if PLZ4 == 3500 & wohnort == "Krems an der Donau"
replace foreignWO = 1 if PLZ4 == 3580 & wohnort == "Alfaz del Pi"
replace foreignWO = 1 if PLZ4 == 4150 & wohnort == "Porto"
replace foreignWO = 1 if PLZ4 == 5960 & wohnort == "Itzig"
replace foreignWO = 1 if PLZ4 == 6100 & wohnort == "Seefeld"
replace foreignWO = 1 if PLZ4 == 6800 & wohnort == "Feldkirch"
replace foreignWO = 1 if PLZ4 == 6811 & wohnort == "Göfis"
replace foreignWO = 1 if PLZ4 == 6812 & wohnort == "Meiningen"
replace foreignWO = 1 if PLZ4 == 6840 & wohnort == "Götzis"
replace foreignWO = 1 if PLZ4 == 6841 & wohnort == "Mäder"
replace foreignWO = 1 if PLZ4 == 6845 & wohnort == "Hohenems"
replace foreignWO = 1 if PLZ4 == 6861 & wohnort == "Alberschwende"
replace foreignWO = 1 if PLZ4 == 6911 & wohnort == "Campione d'Italia"
replace foreignWO = 1 if PLZ4 == 6973 & wohnort == "Höchst"
replace foreignWO = 1 if PLZ4 == 7600 & wohnort == "Stellenbosch"
replace foreignWO = 1 if PLZ4 == 8238 & wohnort == "Büsingen"
replace foreignWO = 1 if PLZ4 == 8326 & wohnort == "Horn"
replace foreignWO = 1 if PLZ4 == 9130 & wohnort == "Winnetka"
replace foreignWO = 1 if PLZ4 == 9131 & wohnort == "Grafenstein"
replace foreignWO = 1 if PLZ4 == 9485 & wohnort == "Nendeln"
replace foreignWO = 1 if PLZ4 == 9486 & wohnort == "Schaanwald"
replace foreignWO = 1 if PLZ4 == 9487 & wohnort == "Bendern"
replace foreignWO = 1 if PLZ4 == 9487 & wohnort == "Gamprin-Bendern"
replace foreignWO = 1 if PLZ4 == 9488 & wohnort == "Schellenberg"
replace foreignWO = 1 if PLZ4 == 9490 & wohnort == "Schellenberg"
replace foreignWO = 1 if PLZ4 == 9490 & wohnort == "Vaduz"
replace foreignWO = 1 if PLZ4 == 9491 & wohnort == "Ruggell"
replace foreignWO = 1 if PLZ4 == 9492 & wohnort == "Balzers"
replace foreignWO = 1 if PLZ4 == 9492 & wohnort == "Eschen"
replace foreignWO = 1 if PLZ4 == 9493 & wohnort == "Mauren"
replace foreignWO = 1 if PLZ4 == 9493 & wohnort == "Mauren FL"
replace foreignWO = 1 if PLZ4 == 9494 & wohnort == "Schaan"
replace foreignWO = 1 if PLZ4 == 9494 & wohnort == "Vaduz"
replace foreignWO = 1 if PLZ4 == 9495 & wohnort == "Triesen"
replace foreignWO = 1 if PLZ4 == 9496 & wohnort == "Balzers"
replace foreignWO = 1 if PLZ4 == 9497 & wohnort == "Triesenberg"
replace foreignWO = 1 if PLZ4 == 9498 & wohnort == "Planken"
replace foreignWO = 1 if PLZ4 == 10000 & wohnort == "Pristina"
replace foreignWO = 1 if PLZ4 == 10019 & wohnort == "New York"
replace foreignWO = 1 if PLZ4 == 10431 & wohnort == "Athen"
replace foreignWO = 1 if PLZ4 == 10555 & wohnort == "Berlin"
replace foreignWO = 1 if PLZ4 == 10629 & wohnort == "Berlin"
replace foreignWO = 1 if PLZ4 == 11523 & wohnort == "Palmbrush Trail"
replace foreignWO = 1 if PLZ4 == 12685 & wohnort == "Berlin"
replace foreignWO = 1 if PLZ4 == 13001 & wohnort == "Marseille"
replace foreignWO = 1 if PLZ4 == 13040 & wohnort == "Borgo D'ale"
replace foreignWO = 1 if PLZ4 == 14547 & wohnort == "Dionysos"
replace foreignWO = 1 if PLZ4 == 20122 & wohnort == "Milvignes"
replace foreignWO = 1 if PLZ4 == 22060 & wohnort == "Campione dItalia"
replace foreignWO = 1 if PLZ4 == 25870 & wohnort == "Bonnay"
replace foreignWO = 1 if PLZ4 == 30659 & wohnort == "Hannover"
replace foreignWO = 1 if PLZ4 == 35210 & wohnort == "Izmir"
replace foreignWO = 1 if PLZ4 == 37085 & wohnort == "Göttingen Geismar"
replace foreignWO = 1 if PLZ4 == 37201 & wohnort == "Nashville"
replace foreignWO = 1 if PLZ4 == 39400 & wohnort == "Morez"
replace foreignWO = 1 if PLZ4 == 40237 & wohnort == "Düsseldorf"
replace foreignWO = 1 if PLZ4 == 40629 & wohnort == "Düsseldorf"
replace foreignWO = 1 if PLZ4 == 43100 & wohnort == "Raanana"
replace foreignWO = 1 if PLZ4 == 44240 & wohnort == "La Chapelle-sur-Erdre"
replace foreignWO = 1 if PLZ4 == 50100 & wohnort == "Florence"
replace foreignWO = 1 if PLZ4 == 61350 & wohnort == "Bad homburg"
replace foreignWO = 1 if PLZ4 == 64560 & wohnort == "Riedstadt"
replace foreignWO = 1 if PLZ4 == 65817 & wohnort == "Eppstein"
replace foreignWO = 1 if PLZ4 == 67435 & wohnort == "Neustadt an der Weinstrasse"
replace foreignWO = 1 if PLZ4 == 68220 & wohnort == "Hegenheim"
replace foreignWO = 1 if PLZ4 == 68390 & wohnort == "Sausheim"
replace foreignWO = 1 if PLZ4 == 70184 & wohnort == "Stuttgart"
replace foreignWO = 1 if PLZ4 == 71723 & wohnort == "Grossbottwar"
replace foreignWO = 1 if PLZ4 == 74100 & wohnort == "Ambilly"
replace foreignWO = 1 if PLZ4 == 74130 & wohnort == "Brison"
replace foreignWO = 1 if PLZ4 == 74140 & wohnort == "Chens-sur-Léman"
replace foreignWO = 1 if PLZ4 == 74140 & wohnort == "Messery"
replace foreignWO = 1 if PLZ4 == 74140 & wohnort == "Nernier"
replace foreignWO = 1 if PLZ4 == 74140 & wohnort == "Sciez"

// Foreign cities if wohnland == "Schweiz" and PLZ missing, and No. of observations per Wohnort >= 4
replace foreignWO = 1 if wohnort == "Dubai" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "London" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Singapur" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Athen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Hamburg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Bangkok" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Monaco" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "München" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Londres" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Abu Dhabi" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Dubaï" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Grand Cayman" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Paris" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Beyrouth" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Reutlingen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Mauren" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Stockholm" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Hong Kong" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "St. John" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Casablanca" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Milano" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "New York" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Berlin" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Habère-Lullin" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "San José" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Dubai Marina" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Les Adrets" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Gibraltar" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Montevideo" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Weil am Rhein" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Frankfurt am Main" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Gjakove" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Singapour" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Dornbirn" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Graz" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Lörrach" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Monte-Carlo" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Toronto" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Annecy" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Düsseldorf" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Helsinki" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Pristina" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Archamps" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Bruxelles" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Como" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Monte Carlo" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Moskau" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Schoten" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Shanghai" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Tel Aviv" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Tokyo" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Vancouver" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Voëns" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Washington D" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Barcelona" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Buenos Aires" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Detmold" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Divonne-les-Bains" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kümmersbruck" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Teneriffa" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Büsingen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Dresden" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Lunden" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Bad Säckingen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Challex" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Douvaine" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Guernsey" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Hohenems" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Istanbul" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Izlake" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "San Bartolomeo V" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Santeramo" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Aloha" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Beijing" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Bregenz" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Dalian" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Hirschberg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Hod" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Konstanz" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Madrid" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Moscou" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Nairobi" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Nürnberg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Sergy" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "St. Louis" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Stuttgart" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Séoul" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Thoiry" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Viry" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Amsterdam" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Castelnau-le-Lez" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Cottbus" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Datchet" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Djursholm" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Gex" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Hall" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Herzliya Pituach" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kapstadt" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kiev" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kingwood" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Lugrin" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Prag" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "S. Bartolomeo V" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Sofia" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Täby" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Waldegg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Wiesbaden" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Andratx" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Bad Homburg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Bisingen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Divonnes-les-Bains" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Doussard" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Grandola" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kopenhagen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "La Muraz" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Laveno Mombello" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Maldonado" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Malmö" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Nassau" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Ohringen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Ornex" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Phuket" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Tel-Aviv" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Weingarten" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Annemasse" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Augsburg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Collonges-sous-Salève" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Den Haag" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Gaillard" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Gjilan" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Gravina" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Heilbronn" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Holzgerlingen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Huningue" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Istog" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Jakarta" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Jávea" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kelkheim" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Krimpen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "La Haye" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Le Sapey" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Luxembourg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Luxemburg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Mexico D" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Milan" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Montréal" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Oberursel" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Ouagadougou" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Pers-Jussy" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Prague" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Prévessin-Moëns" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "San Francisco" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Sechelt" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Shtime" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Southlake" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "St-Gildas de Rhuys" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Varsovie" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Vienne" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Wittelsheim" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Armação de Pêra" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Baie d'Urfe" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Barasso" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Besançon" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Blizikuce/Vrba-Tudorovici" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Bloomington" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Cabiate" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Calice Ligure" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Chiesa" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Chonburi" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Clermont" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Cserszegtomaj" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Dietingen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Dortmund" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Feigères" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Feldkirch" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Fowler" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Halle" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Hong-Kong" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Hongkong" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kalmar" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kiew" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kleines Wiesenthal" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kuala Lumpur" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Les Houches" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Limassol" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Lund" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Lyon" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Maccagno con Pino" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Mannheim" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Meppen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Midland" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Neumünster" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Oltingue" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Ottawa" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Poole" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Poway" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Raleigh" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Ronago" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Rougement-le-Château" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Schifferstadt" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Sciez" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "St-Cergues" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Todtmoos" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Tortola" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Wellington" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Acapulco" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Aelesund" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Alland" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Antibes" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Argra" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Aschaffenburg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Athènes" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Attenschwiler" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Au am Rhein" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Badenweiler" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Barcellona" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Berlin DE" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Bingen am Rhein" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Boca Raton" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Buchen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Budapest" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Bürg-Neuhaus" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Cartagena" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Cartagena des Indias" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Castell'Alfero" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Cernex" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Chens-sur-Léman" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Cipressa" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Cluses" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Colleyville" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Collonge-sous-Salève" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Dangio" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Dienersdorf" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Douglas" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Elmira N" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Emmendingen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Excenevex" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Falmenta" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Ferney-Voltaire" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Fürstenzell" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Gailingen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "George Town" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Gjakovë" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Golfe-Juan" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Hegenheim" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Herceg Novi" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Hinterbrühl" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Jeddah" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kharkiv" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Kozakli-Akent" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Köln" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Lacanau" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Lagos" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Laguna Niguel" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Le Caire" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Leymen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Liberty Hill" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Marburg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Messery" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Mettendorf" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Mils" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Neuilly" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "New Jersey" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "North Norfolk" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Novo" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Olgiate Comasco" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Oslo" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Panama City" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Pierre de Bresse" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Poisson" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Prishtina" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Providenciales" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Raanana" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Rayong" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Reno" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Reuthe" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Rio de Janeiro" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Roma" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Rosengarten" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Rosà" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Santiago du Chili" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Smith" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Soresina" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Spring Lake" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Sprockhövel" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "St-Genis-Pouilly" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "St-Julien-en-Genevois" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "St-Louis" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "St. Maarten" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Strasbourg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Sulzbach" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Sydney" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Tannenheim" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Taufkirchen" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Tehran-Iran" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Thessaloniki" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Tianjin" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Trinity" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Uppsala" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Venlo" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Vescia di Foligno" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Ville-la-Grand" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Villers-le-Lac" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Wartberg" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Wehr" & wohnland == "Schweiz"
replace foreignWO = 1 if wohnort == "Winzer" & wohnland == "Schweiz"

// Define foreign Wohnort if Wohnland != "Schweiz", PLZ4 == missing
replace foreignWO = 1 if wohnland!="Schweiz" & PLZ4 == .
// Correct foreign Wohnort if Wohnland != "Schweiz", PLZ4 == missing, and No. of observations per Wohnort >= 30  
// --> Zero corrections to be made.

rename PLZ4 W_PLZ4

// missing PLZ for Bürgerort in Schweiz (merge on "Wohnort" using Post-PLZ observations)
gen B_PLZ4 = plzbürgerort

/*
preserve // empty PLZ4 for Bürgerort
gen PLZBmiss = 1 if B_PLZ4==.
keep if B_PLZ4==.
gen nation_ch = 1 if nationalität == "Schweiz"
drop plzbürgerort B_PLZ4
gen count = 1
collapse (sum) count, by(bürgerort)
duplicates drop bürgerort, force
rename bürgerort gdename // just for merge with List of Gdename
sort gdename
save "$path\02_Processed_data\11_Directors_1994_2018/Bisnode_missingBuergerPLZ_temp1.dta", replace
merge m:1 gdename using "$path\02_Processed_data\11_Directors_1994_2018/Gde-PLZ-List_temp.dta"
keep if _merge==3
keep gdename PLZ4
duplicates drop gdename PLZ4, force
save "$path\02_Processed_data\11_Directors_1994_2018/Bisnode_missingBuergerPLZ-merge_temp1.dta", replace
restore
*/

gen PLZBmiss = 1 if B_PLZ4==.
replace B_PLZ4 = 9547 if bürgerort == "Aadorf" & PLZBmiss == 1
replace B_PLZ4 = 5000 if bürgerort == "Aarau" & PLZBmiss == 1
replace B_PLZ4 = 3270 if bürgerort == "Aarberg" & PLZBmiss == 1
replace B_PLZ4 = 4663 if bürgerort == "Aarburg" & PLZBmiss == 1
replace B_PLZ4 = 4912 if bürgerort == "Aarwangen" & PLZBmiss == 1
replace B_PLZ4 = 5646 if bürgerort == "Abtwil" & PLZBmiss == 1
replace B_PLZ4 = 6722 if bürgerort == "Acquarossa" & PLZBmiss == 1
replace B_PLZ4 = 3715 if bürgerort == "Adelboden" & PLZBmiss == 1
replace B_PLZ4 = 6043 if bürgerort == "Adligenswil" & PLZBmiss == 1
replace B_PLZ4 = 8452 if bürgerort == "Adlikon" & PLZBmiss == 1
replace B_PLZ4 = 8134 if bürgerort == "Adliswil" & PLZBmiss == 1
replace B_PLZ4 = 4714 if bürgerort == "Aedermannsdorf" & PLZBmiss == 1
replace B_PLZ4 = 3426 if bürgerort == "Aefligen" & PLZBmiss == 1
replace B_PLZ4 = 2558 if bürgerort == "Aegerten" & PLZBmiss == 1
replace B_PLZ4 = 3703 if bürgerort == "Aeschi bei Spiez" & PLZBmiss == 1
replace B_PLZ4 = 3672 if bürgerort == "Aeschlen" & PLZBmiss == 1
replace B_PLZ4 = 4588 if bürgerort == "Aetingen" & PLZBmiss == 1
replace B_PLZ4 = 8914 if bürgerort == "Aeugst" & PLZBmiss == 1
replace B_PLZ4 = 9562 if bürgerort == "Affeltrangen" & PLZBmiss == 1
replace B_PLZ4 = 8910 if bürgerort == "Affoltern am Albis" & PLZBmiss == 1
replace B_PLZ4 = 3462 if bürgerort == "Affoltern im Emmental" & PLZBmiss == 1
replace B_PLZ4 = 3951 if bürgerort == "Agarn" & PLZBmiss == 1
replace B_PLZ4 = 6982 if bürgerort == "Agno" & PLZBmiss == 1
replace B_PLZ4 = 6780 if bürgerort == "Airolo" & PLZBmiss == 1
replace B_PLZ4 = 3955 if bürgerort == "Albinen" & PLZBmiss == 1
replace B_PLZ4 = 7472 if bürgerort == "Albula/Alvra" & PLZBmiss == 1
replace B_PLZ4 = 3473 if bürgerort == "Alchenstorf" & PLZBmiss == 1
replace B_PLZ4 = 2942 if bürgerort == "Alle" & PLZBmiss == 1
replace B_PLZ4 = 4123 if bürgerort == "Allschwil" & PLZBmiss == 1
replace B_PLZ4 = 6147 if bürgerort == "Altbüron" & PLZBmiss == 1
replace B_PLZ4 = 8852 if bürgerort == "Altendorf" & PLZBmiss == 1
replace B_PLZ4 = 1715 if bürgerort == "Alterswil" & PLZBmiss == 1
replace B_PLZ4 = 6246 if bürgerort == "Altishofen" & PLZBmiss == 1
replace B_PLZ4 = 8595 if bürgerort == "Altnau" & PLZBmiss == 1
replace B_PLZ4 = 6939 if bürgerort == "Alto Malcantone" & PLZBmiss == 1
replace B_PLZ4 = 9450 if bürgerort == "Altstätten" & PLZBmiss == 1
replace B_PLZ4 = 6286 if bürgerort == "Altwis" & PLZBmiss == 1
replace B_PLZ4 = 7492 if bürgerort == "Alvaneu" & PLZBmiss == 1
replace B_PLZ4 = 8873 if bürgerort == "Amden" & PLZBmiss == 1
replace B_PLZ4 = 5600 if bürgerort == "Ammerswil" & PLZBmiss == 1
replace B_PLZ4 = 8580 if bürgerort == "Amriswil" & PLZBmiss == 1
replace B_PLZ4 = 3633 if bürgerort == "Amsoldingen" & PLZBmiss == 1
replace B_PLZ4 = 7442 if bürgerort == "Andeer" & PLZBmiss == 1
replace B_PLZ4 = 7159 if bürgerort == "Andiast" & PLZBmiss == 1
replace B_PLZ4 = 9545 if bürgerort == "Anetswil" & PLZBmiss == 1
replace B_PLZ4 = 1247 if bürgerort == "Anières" & PLZBmiss == 1
replace B_PLZ4 = 3961 if bürgerort == "Anniviers" & PLZBmiss == 1
replace B_PLZ4 = 9050 if bürgerort == "Appenzell" & PLZBmiss == 1
replace B_PLZ4 = 1143 if bürgerort == "Apples" & PLZBmiss == 1
replace B_PLZ4 = 1974 if bürgerort == "Arbaz" & PLZBmiss == 1
replace B_PLZ4 = 4424 if bürgerort == "Arboldswil" & PLZBmiss == 1
replace B_PLZ4 = 9320 if bürgerort == "Arbon" & PLZBmiss == 1
replace B_PLZ4 = 7546 if bürgerort == "Ardez" & PLZBmiss == 1
replace B_PLZ4 = 1957 if bürgerort == "Ardon" & PLZBmiss == 1
replace B_PLZ4 = 4422 if bürgerort == "Arisdorf" & PLZBmiss == 1
replace B_PLZ4 = 5628 if bürgerort == "Aristau" & PLZBmiss == 1
replace B_PLZ4 = 4144 if bürgerort == "Arlesheim" & PLZBmiss == 1
replace B_PLZ4 = 3508 if bürgerort == "Arni (BE)" & PLZBmiss == 1
replace B_PLZ4 = 6823 if bürgerort == "Arogno" & PLZBmiss == 1
replace B_PLZ4 = 7027 if bürgerort == "Arosa" & PLZBmiss == 1
replace B_PLZ4 = 6410 if bürgerort == "Arth" & PLZBmiss == 1
replace B_PLZ4 = 6543 if bürgerort == "Arvigo" & PLZBmiss == 1
replace B_PLZ4 = 1273 if bürgerort == "Arzier-Le Muids" & PLZBmiss == 1
replace B_PLZ4 = 6864 if bürgerort == "Arzo" & PLZBmiss == 1
replace B_PLZ4 = 1042 if bürgerort == "Assens" & PLZBmiss == 1
replace B_PLZ4 = 6999 if bürgerort == "Astano" & PLZBmiss == 1
replace B_PLZ4 = 2954 if bürgerort == "Asuel" & PLZBmiss == 1
replace B_PLZ4 = 1616 if bürgerort == "Attalens" & PLZBmiss == 1
replace B_PLZ4 = 5056 if bürgerort == "Attelwil" & PLZBmiss == 1
replace B_PLZ4 = 4536 if bürgerort == "Attiswil" & PLZBmiss == 1
replace B_PLZ4 = 4302 if bürgerort == "Augst" & PLZBmiss == 1
replace B_PLZ4 = 1484 if bürgerort == "Aumont" & PLZBmiss == 1
replace B_PLZ4 = 6661 if bürgerort == "Auressio" & PLZBmiss == 1
replace B_PLZ4 = 3938 if bürgerort == "Ausserberg" & PLZBmiss == 1
replace B_PLZ4 = 4944 if bürgerort == "Auswil" & PLZBmiss == 1
replace B_PLZ4 = 1475 if bürgerort == "Autavaux" & PLZBmiss == 1
replace B_PLZ4 = 1742 if bürgerort == "Autigny" & PLZBmiss == 1
replace B_PLZ4 = 2012 if bürgerort == "Auvernier" & PLZBmiss == 1
replace B_PLZ4 = 5644 if bürgerort == "Auw" & PLZBmiss == 1
replace B_PLZ4 = 6670 if bürgerort == "Avegno Gordevio" & PLZBmiss == 1
replace B_PLZ4 = 1580 if bürgerort == "Avenches" & PLZBmiss == 1
replace B_PLZ4 = 7447 if bürgerort == "Avers" & PLZBmiss == 1
replace B_PLZ4 = 1643 if bürgerort == "Avry-devant-Pont" & PLZBmiss == 1
replace B_PLZ4 = 1754 if bürgerort == "Avry-sur-Matran" & PLZBmiss == 1
replace B_PLZ4 = 1237 if bürgerort == "Avully" & PLZBmiss == 1
replace B_PLZ4 = 1285 if bürgerort == "Avusy" & PLZBmiss == 1
replace B_PLZ4 = 1966 if bürgerort == "Ayent" & PLZBmiss == 1
replace B_PLZ4 = 3961 if bürgerort == "Ayer" & PLZBmiss == 1
replace B_PLZ4 = 6346 if bürgerort == "Baar" & PLZBmiss == 1
replace B_PLZ4 = 8184 if bürgerort == "Bachenbülach" & PLZBmiss == 1
replace B_PLZ4 = 7310 if bürgerort == "Bad Ragaz" & PLZBmiss == 1
replace B_PLZ4 = 5400 if bürgerort == "Baden" & PLZBmiss == 1
replace B_PLZ4 = 1947 if bürgerort == "Bagnes" & PLZBmiss == 1
replace B_PLZ4 = 6828 if bürgerort == "Balerna" & PLZBmiss == 1
replace B_PLZ4 = 9436 if bürgerort == "Balgach" & PLZBmiss == 1
replace B_PLZ4 = 1338 if bürgerort == "Ballaigues" & PLZBmiss == 1
replace B_PLZ4 = 6275 if bürgerort == "Ballwil" & PLZBmiss == 1
replace B_PLZ4 = 3254 if bürgerort == "Balm bei Messen" & PLZBmiss == 1
replace B_PLZ4 = 4710 if bürgerort == "Balsthal" & PLZBmiss == 1
replace B_PLZ4 = 3256 if bürgerort == "Bangerten" & PLZBmiss == 1
replace B_PLZ4 = 1257 if bürgerort == "Bardonnex" & PLZBmiss == 1
replace B_PLZ4 = 1666 if bürgerort == "Bas-Intyamon" & PLZBmiss == 1
replace B_PLZ4 = 1788 if bürgerort == "Bas-Vully" & PLZBmiss == 1
replace B_PLZ4 = 8255 if bürgerort == "Basadingen-Schlattingen" & PLZBmiss == 1
replace B_PLZ4 = 4089 if bürgerort == "Basel" & PLZBmiss == 1
replace B_PLZ4 = 2925 if bürgerort == "Basse-Allaine" & PLZBmiss == 1
replace B_PLZ4 = 2854 if bürgerort == "Bassecourt" & PLZBmiss == 1
replace B_PLZ4 = 8303 if bürgerort == "Bassersdorf" & PLZBmiss == 1
replace B_PLZ4 = 1446 if bürgerort == "Baulmes" & PLZBmiss == 1
replace B_PLZ4 = 8493 if bürgerort == "Bauma" & PLZBmiss == 1
replace B_PLZ4 = 3803 if bürgerort == "Beatenberg" & PLZBmiss == 1
replace B_PLZ4 = 6375 if bürgerort == "Beckenried" & PLZBmiss == 1
replace B_PLZ4 = 8228 if bürgerort == "Beggingen" & PLZBmiss == 1
replace B_PLZ4 = 5637 if bürgerort == "Beinwil (Freiamt)" & PLZBmiss == 1
replace B_PLZ4 = 5712 if bürgerort == "Beinwil am See" & PLZBmiss == 1
replace B_PLZ4 = 1782 if bürgerort == "Belfaux" & PLZBmiss == 1
replace B_PLZ4 = 4512 if bürgerort == "Bellach" & PLZBmiss == 1
replace B_PLZ4 = 5454 if bürgerort == "Bellikon" & PLZBmiss == 1
replace B_PLZ4 = 6582 if bürgerort == "Bellinzona" & PLZBmiss == 1
replace B_PLZ4 = 3997 if bürgerort == "Bellwald" & PLZBmiss == 1
replace B_PLZ4 = 1773 if bürgerort == "Belmont-Broye" & PLZBmiss == 1
replace B_PLZ4 = 1432 if bürgerort == "Belmont-sur-Yverdon" & PLZBmiss == 1
replace B_PLZ4 = 3123 if bürgerort == "Belp" & PLZBmiss == 1
replace B_PLZ4 = 4431 if bürgerort == "Bennwil" & PLZBmiss == 1
replace B_PLZ4 = 7482 if bürgerort == "Bergün Filisur" & PLZBmiss == 1
replace B_PLZ4 = 7484 if bürgerort == "Bergün/Bravuogn" & PLZBmiss == 1
replace B_PLZ4 = 8267 if bürgerort == "Berlingen" & PLZBmiss == 1
replace B_PLZ4 = 3341 if bürgerort == "Bern" & PLZBmiss == 1
replace B_PLZ4 = 9442 if bürgerort == "Berneck" & PLZBmiss == 1
replace B_PLZ4 = 1233 if bürgerort == "Bernex" & PLZBmiss == 1
replace B_PLZ4 = 6215 if bürgerort == "Beromünster" & PLZBmiss == 1
replace B_PLZ4 = 8543 if bürgerort == "Bertschikon" & PLZBmiss == 1
replace B_PLZ4 = 1609 if bürgerort == "Besencens" & PLZBmiss == 1
replace B_PLZ4 = 3991 if bürgerort == "Betten" & PLZBmiss == 1
replace B_PLZ4 = 3366 if bürgerort == "Bettenhausen" & PLZBmiss == 1
replace B_PLZ4 = 9553 if bürgerort == "Bettwiesen" & PLZBmiss == 1
replace B_PLZ4 = 5618 if bürgerort == "Bettwil" & PLZBmiss == 1
replace B_PLZ4 = 1880 if bürgerort == "Bex" & PLZBmiss == 1
replace B_PLZ4 = 6710 if bürgerort == "Biasca" & PLZBmiss == 1
replace B_PLZ4 = 4562 if bürgerort == "Biberist" & PLZBmiss == 1
replace B_PLZ4 = 4105 if bürgerort == "Biel-Benken" & PLZBmiss == 1
replace B_PLZ4 = 3507 if bürgerort == "Biglen" & PLZBmiss == 1
replace B_PLZ4 = 1681 if bürgerort == "Billens" & PLZBmiss == 1
replace B_PLZ4 = 3996 if bürgerort == "Binn" & PLZBmiss == 1
replace B_PLZ4 = 4102 if bürgerort == "Binningen" & PLZBmiss == 1
replace B_PLZ4 = 6992 if bürgerort == "Bioggio" & PLZBmiss == 1
replace B_PLZ4 = 5244 if bürgerort == "Birrhard" & PLZBmiss == 1
replace B_PLZ4 = 9220 if bürgerort == "Bischofszell" & PLZBmiss == 1
replace B_PLZ4 = 3982 if bürgerort == "Bitsch" & PLZBmiss == 1
replace B_PLZ4 = 1145 if bürgerort == "Bière" & PLZBmiss == 1
replace B_PLZ4 = 3919 if bürgerort == "Blatten" & PLZBmiss == 1
replace B_PLZ4 = 3368 if bürgerort == "Bleienbach" & PLZBmiss == 1
replace B_PLZ4 = 3674 if bürgerort == "Bleiken bei Oberdiessbach" & PLZBmiss == 1
replace B_PLZ4 = 6717 if bürgerort == "Blenio" & PLZBmiss == 1
replace B_PLZ4 = 1807 if bürgerort == "Blonay" & PLZBmiss == 1
replace B_PLZ4 = 3638 if bürgerort == "Blumenstein" & PLZBmiss == 1
replace B_PLZ4 = 6743 if bürgerort == "Bodio" & PLZBmiss == 1
replace B_PLZ4 = 3065 if bürgerort == "Bolligen" & PLZBmiss == 1
replace B_PLZ4 = 3766 if bürgerort == "Boltigen" & PLZBmiss == 1
replace B_PLZ4 = 7402 if bürgerort == "Bonaduz" & PLZBmiss == 1
replace B_PLZ4 = 7606 if bürgerort == "Bondo" & PLZBmiss == 1
replace B_PLZ4 = 2944 if bürgerort == "Bonfol" & PLZBmiss == 1
replace B_PLZ4 = 1615 if bürgerort == "Bossonnens" & PLZBmiss == 1
replace B_PLZ4 = 5623 if bürgerort == "Boswil" & PLZBmiss == 1
replace B_PLZ4 = 1041 if bürgerort == "Bottens" & PLZBmiss == 1
replace B_PLZ4 = 4814 if bürgerort == "Bottenwil" & PLZBmiss == 1
replace B_PLZ4 = 8598 if bürgerort == "Bottighofen" & PLZBmiss == 1
replace B_PLZ4 = 2017 if bürgerort == "Boudry" & PLZBmiss == 1
replace B_PLZ4 = 1172 if bürgerort == "Bougy-Villars" & PLZBmiss == 1
replace B_PLZ4 = 1699 if bürgerort == "Bouloz" & PLZBmiss == 1
replace B_PLZ4 = 1091 if bürgerort == "Bourg-en-Lavaux" & PLZBmiss == 1
replace B_PLZ4 = 1035 if bürgerort == "Bournens" & PLZBmiss == 1
replace B_PLZ4 = 1932 if bürgerort == "Bovernier" & PLZBmiss == 1
replace B_PLZ4 = 3533 if bürgerort == "Bowil" & PLZBmiss == 1
replace B_PLZ4 = 9502 if bürgerort == "Braunau" & PLZBmiss == 1
replace B_PLZ4 = 7516 if bürgerort == "Bregaglia" & PLZBmiss == 1
replace B_PLZ4 = 6837 if bürgerort == "Breggia" & PLZBmiss == 1
replace B_PLZ4 = 7159 if bürgerort == "Breil/Brigels" & PLZBmiss == 1
replace B_PLZ4 = 5620 if bürgerort == "Bremgarten (AG)" & PLZBmiss == 1
replace B_PLZ4 = 3047 if bürgerort == "Bremgarten bei Bern" & PLZBmiss == 1
replace B_PLZ4 = 3856 if bürgerort == "Brienzwiler" & PLZBmiss == 1
replace B_PLZ4 = 3900 if bürgerort == "Brig-Glis" & PLZBmiss == 1
replace B_PLZ4 = 4805 if bürgerort == "Brittnau" & PLZBmiss == 1
replace B_PLZ4 = 1636 if bürgerort == "Broc" & PLZBmiss == 1
replace B_PLZ4 = 9512 if bürgerort == "Bronschhofen" & PLZBmiss == 1
replace B_PLZ4 = 2149 if bürgerort == "Brot-Dessous" & PLZBmiss == 1
replace B_PLZ4 = 5200 if bürgerort == "Brugg" & PLZBmiss == 1
replace B_PLZ4 = 3307 if bürgerort == "Brunnenthal" & PLZBmiss == 1
replace B_PLZ4 = 6827 if bürgerort == "Brusino Arsizio" & PLZBmiss == 1
replace B_PLZ4 = 7744 if bürgerort == "Brusio" & PLZBmiss == 1
replace B_PLZ4 = 2555 if bürgerort == "Brügg" & PLZBmiss == 1
replace B_PLZ4 = 1719 if bürgerort == "Brünisried" & PLZBmiss == 1
replace B_PLZ4 = 3237 if bürgerort == "Brüttelen" & PLZBmiss == 1
replace B_PLZ4 = 8414 if bürgerort == "Buch am Irchel" & PLZBmiss == 1
replace B_PLZ4 = 8524 if bürgerort == "Buch bei Frauenfeld" & PLZBmiss == 1
replace B_PLZ4 = 4587 if bürgerort == "Buchegg" & PLZBmiss == 1
replace B_PLZ4 = 3615 if bürgerort == "Buchholterberg" & PLZBmiss == 1
replace B_PLZ4 = 6033 if bürgerort == "Buchrain" & PLZBmiss == 1
replace B_PLZ4 = 1630 if bürgerort == "Bulle" & PLZBmiss == 1
replace B_PLZ4 = 1452 if bürgerort == "Bullet" & PLZBmiss == 1
replace B_PLZ4 = 6374 if bürgerort == "Buochs" & PLZBmiss == 1
replace B_PLZ4 = 3400 if bürgerort == "Burgdorf" & PLZBmiss == 1
replace B_PLZ4 = 1183 if bürgerort == "Bursins" & PLZBmiss == 1
replace B_PLZ4 = 1030 if bürgerort == "Bussigny-près-Lausanne" & PLZBmiss == 1
replace B_PLZ4 = 1608 if bürgerort == "Bussigny-sur-Oron" & PLZBmiss == 1
replace B_PLZ4 = 9565 if bürgerort == "Bussnang" & PLZBmiss == 1
replace B_PLZ4 = 4917 if bürgerort == "Busswil bei Melchnau" & PLZBmiss == 1
replace B_PLZ4 = 2116 if bürgerort == "Buttes" & PLZBmiss == 1
replace B_PLZ4 = 6018 if bürgerort == "Buttisholz" & PLZBmiss == 1
replace B_PLZ4 = 5632 if bürgerort == "Buttwil" & PLZBmiss == 1
replace B_PLZ4 = 4463 if bürgerort == "Buus" & PLZBmiss == 1
replace B_PLZ4 = 8344 if bürgerort == "Bäretswil" & PLZBmiss == 1
replace B_PLZ4 = 3323 if bürgerort == "Bäriswil" & PLZBmiss == 1
replace B_PLZ4 = 3315 if bürgerort == "Bätterkinden" & PLZBmiss == 1
replace B_PLZ4 = 2014 if bürgerort == "Bôle" & PLZBmiss == 1
replace B_PLZ4 = 3806 if bürgerort == "Bönigen" & PLZBmiss == 1
replace B_PLZ4 = 5315 if bürgerort == "Böttstein" & PLZBmiss == 1
replace B_PLZ4 = 5225 if bürgerort == "Bözberg" & PLZBmiss == 1
replace B_PLZ4 = 3215 if bürgerort == "Büchslen" & PLZBmiss == 1
replace B_PLZ4 = 3263 if bürgerort == "Büetigen" & PLZBmiss == 1
replace B_PLZ4 = 5624 if bürgerort == "Bünzen" & PLZBmiss == 1
replace B_PLZ4 = 3935 if bürgerort == "Bürchen" & PLZBmiss == 1
replace B_PLZ4 = 3294 if bürgerort == "Büren an der Aare" & PLZBmiss == 1
replace B_PLZ4 = 6233 if bürgerort == "Büron" & PLZBmiss == 1
replace B_PLZ4 = 4227 if bürgerort == "Büsserach" & PLZBmiss == 1
replace B_PLZ4 = 9615 if bürgerort == "Bütschwil" & PLZBmiss == 1
replace B_PLZ4 = 9615 if bürgerort == "Bütschwil-Ganterschwil" & PLZBmiss == 1
replace B_PLZ4 = 8236 if bürgerort == "Büttenhardt" & PLZBmiss == 1
replace B_PLZ4 = 5619 if bürgerort == "Büttikon" & PLZBmiss == 1
replace B_PLZ4 = 6965 if bürgerort == "Cadro" & PLZBmiss == 1
replace B_PLZ4 = 6543 if bürgerort == "Calanca" & PLZBmiss == 1
replace B_PLZ4 = 6528 if bürgerort == "Camorino" & PLZBmiss == 1
replace B_PLZ4 = 6760 if bürgerort == "Campello" & PLZBmiss == 1
replace B_PLZ4 = 6720 if bürgerort == "Campo (Blenio)" & PLZBmiss == 1
replace B_PLZ4 = 6684 if bürgerort == "Campo (Vallemaggia)" & PLZBmiss == 1
replace B_PLZ4 = 7113 if bürgerort == "Camuns" & PLZBmiss == 1
replace B_PLZ4 = 6837 if bürgerort == "Caneggio" & PLZBmiss == 1
replace B_PLZ4 = 6825 if bürgerort == "Capolago" & PLZBmiss == 1
replace B_PLZ4 = 6987 if bürgerort == "Caslano" & PLZBmiss == 1
replace B_PLZ4 = 6875 if bürgerort == "Castel San Pietro" & PLZBmiss == 1
replace B_PLZ4 = 7433 if bürgerort == "Casti-Wergenstein" & PLZBmiss == 1
replace B_PLZ4 = 7027 if bürgerort == "Castiel" & PLZBmiss == 1
replace B_PLZ4 = 7126 if bürgerort == "Castrisch" & PLZBmiss == 1
replace B_PLZ4 = 7422 if bürgerort == "Cazis" & PLZBmiss == 1
replace B_PLZ4 = 6655 if bürgerort == "Centovalli" & PLZBmiss == 1
replace B_PLZ4 = 6683 if bürgerort == "Cerentino" & PLZBmiss == 1
replace B_PLZ4 = 6959 if bürgerort == "Certara" & PLZBmiss == 1
replace B_PLZ4 = 1589 if bürgerort == "Chabrey" & PLZBmiss == 1
replace B_PLZ4 = 3966 if bürgerort == "Chalais" & PLZBmiss == 1
replace B_PLZ4 = 6330 if bürgerort == "Cham" & PLZBmiss == 1
replace B_PLZ4 = 1955 if bürgerort == "Chamoson" & PLZBmiss == 1
replace B_PLZ4 = 2735 if bürgerort == "Champoz" & PLZBmiss == 1
replace B_PLZ4 = 1537 if bürgerort == "Champtauroz" & PLZBmiss == 1
replace B_PLZ4 = 1443 if bürgerort == "Champvent" & PLZBmiss == 1
replace B_PLZ4 = 1874 if bürgerort == "Champéry" & PLZBmiss == 1
replace B_PLZ4 = 1284 if bürgerort == "Chancy" & PLZBmiss == 1
replace B_PLZ4 = 3961 if bürgerort == "Chandolin" & PLZBmiss == 1
replace B_PLZ4 = 1801 if bürgerort == "Chardonne" & PLZBmiss == 1
replace B_PLZ4 = 1637 if bürgerort == "Charmey" & PLZBmiss == 1
replace B_PLZ4 = 2947 if bürgerort == "Charmoille" & PLZBmiss == 1
replace B_PLZ4 = 1464 if bürgerort == "Chavannes-le-Chêne" & PLZBmiss == 1
replace B_PLZ4 = 1435 if bürgerort == "Chavornay" & PLZBmiss == 1
replace B_PLZ4 = 3971 if bürgerort == "Chermignon" & PLZBmiss == 1
replace B_PLZ4 = 2906 if bürgerort == "Chevenez" & PLZBmiss == 1
replace B_PLZ4 = 1474 if bürgerort == "Cheyres-Châbles" & PLZBmiss == 1
replace B_PLZ4 = 6830 if bürgerort == "Chiasso" & PLZBmiss == 1
replace B_PLZ4 = 6764 if bürgerort == "Chiggiogna" & PLZBmiss == 1
replace B_PLZ4 = 7001 if bürgerort == "Chur" & PLZBmiss == 1
replace B_PLZ4 = 7075 if bürgerort == "Churwalden" & PLZBmiss == 1
replace B_PLZ4 = 1660 if bürgerort == "Château-d'Oex" & PLZBmiss == 1
replace B_PLZ4 = 1610 if bürgerort == "Châtillens" & PLZBmiss == 1
replace B_PLZ4 = 1744 if bürgerort == "Chénens" & PLZBmiss == 1
replace B_PLZ4 = 1720 if bürgerort == "Chésopelloz" & PLZBmiss == 1
replace B_PLZ4 = 1224 if bürgerort == "Chêne-Bougeries" & PLZBmiss == 1
replace B_PLZ4 = 1225 if bürgerort == "Chêne-Bourg" & PLZBmiss == 1
replace B_PLZ4 = 6702 if bürgerort == "Claro" & PLZBmiss == 1
replace B_PLZ4 = 1595 if bürgerort == "Clavaleyres" & PLZBmiss == 1
replace B_PLZ4 = 2882 if bürgerort == "Clos du Doubs" & PLZBmiss == 1
replace B_PLZ4 = 6877 if bürgerort == "Coldrerio" & PLZBmiss == 1
replace B_PLZ4 = 1239 if bürgerort == "Collex-Bossy" & PLZBmiss == 1
replace B_PLZ4 = 6925 if bürgerort == "Collina d'Oro" & PLZBmiss == 1
replace B_PLZ4 = 1893 if bürgerort == "Collombey-Muraz" & PLZBmiss == 1
replace B_PLZ4 = 1222 if bürgerort == "Collonge-Bellerive" & PLZBmiss == 1
replace B_PLZ4 = 1903 if bürgerort == "Collonges" & PLZBmiss == 1
replace B_PLZ4 = 1223 if bürgerort == "Cologny" & PLZBmiss == 1
replace B_PLZ4 = 6949 if bürgerort == "Comano" & PLZBmiss == 1
replace B_PLZ4 = 1535 if bürgerort == "Combremont-le-Grand" & PLZBmiss == 1
replace B_PLZ4 = 1291 if bürgerort == "Commugny" & PLZBmiss == 1
replace B_PLZ4 = 1426 if bürgerort == "Concise" & PLZBmiss == 1
replace B_PLZ4 = 1232 if bürgerort == "Confignon" & PLZBmiss == 1
replace B_PLZ4 = 1976 if bürgerort == "Conthey" & PLZBmiss == 1
replace B_PLZ4 = 1296 if bürgerort == "Coppet" & PLZBmiss == 1
replace B_PLZ4 = 2036 if bürgerort == "Corcelles-Cormondrèche" & PLZBmiss == 1
replace B_PLZ4 = 1082 if bürgerort == "Corcelles-le-Jorat" & PLZBmiss == 1
replace B_PLZ4 = 1562 if bürgerort == "Corcelles-près-Payerne" & PLZBmiss == 1
replace B_PLZ4 = 1374 if bürgerort == "Corcelles-sur-Chavornay" & PLZBmiss == 1
replace B_PLZ4 = 1720 if bürgerort == "Corminboeuf" & PLZBmiss == 1
replace B_PLZ4 = 2087 if bürgerort == "Cornaux" & PLZBmiss == 1
replace B_PLZ4 = 2952 if bürgerort == "Cornol" & PLZBmiss == 1
replace B_PLZ4 = 1727 if bürgerort == "Corpataux" & PLZBmiss == 1
replace B_PLZ4 = 1802 if bürgerort == "Corseaux" & PLZBmiss == 1
replace B_PLZ4 = 1747 if bürgerort == "Corserey" & PLZBmiss == 1
replace B_PLZ4 = 1808 if bürgerort == "Corsier-sur-Vevey" & PLZBmiss == 1
replace B_PLZ4 = 2016 if bürgerort == "Cortaillod" & PLZBmiss == 1
replace B_PLZ4 = 1304 if bürgerort == "Cossonay" & PLZBmiss == 1
replace B_PLZ4 = 2853 if bürgerort == "Courfaivre" & PLZBmiss == 1
replace B_PLZ4 = 2950 if bürgerort == "Courgenay" & PLZBmiss == 1
replace B_PLZ4 = 1796 if bürgerort == "Courgevaux" & PLZBmiss == 1
replace B_PLZ4 = 2823 if bürgerort == "Courroux" & PLZBmiss == 1
replace B_PLZ4 = 2738 if bürgerort == "Court" & PLZBmiss == 1
replace B_PLZ4 = 2905 if bürgerort == "Courtedoux" & PLZBmiss == 1
replace B_PLZ4 = 2108 if bürgerort == "Couvet" & PLZBmiss == 1
replace B_PLZ4 = 1299 if bürgerort == "Crans-près-Céligny" & PLZBmiss == 1
replace B_PLZ4 = 1526 if bürgerort == "Cremin" & PLZBmiss == 1
replace B_PLZ4 = 1023 if bürgerort == "Crissier" & PLZBmiss == 1
replace B_PLZ4 = 6980 if bürgerort == "Croglio" & PLZBmiss == 1
replace B_PLZ4 = 1406 if bürgerort == "Cronay" & PLZBmiss == 1
replace B_PLZ4 = 2746 if bürgerort == "Crémines" & PLZBmiss == 1
replace B_PLZ4 = 1653 if bürgerort == "Crésuz" & PLZBmiss == 1
replace B_PLZ4 = 1588 if bürgerort == "Cudrefin" & PLZBmiss == 1
replace B_PLZ4 = 1096 if bürgerort == "Cully" & PLZBmiss == 1
replace B_PLZ4 = 7142 if bürgerort == "Cumbels" & PLZBmiss == 1
replace B_PLZ4 = 6963 if bürgerort == "Cureggia" & PLZBmiss == 1
replace B_PLZ4 = 6981 if bürgerort == "Curio" & PLZBmiss == 1
replace B_PLZ4 = 1521 if bürgerort == "Curtilles" & PLZBmiss == 1
replace B_PLZ4 = 1298 if bürgerort == "Céligny" & PLZBmiss == 1
replace B_PLZ4 = 6252 if bürgerort == "Dagmersellen" & PLZBmiss == 1
replace B_PLZ4 = 6383 if bürgerort == "Dallenwil" & PLZBmiss == 1
replace B_PLZ4 = 2933 if bürgerort == "Damphreux" & PLZBmiss == 1
replace B_PLZ4 = 7270 if bürgerort == "Davos" & PLZBmiss == 1
replace B_PLZ4 = 9113 if bürgerort == "Degersheim" & PLZBmiss == 1
replace B_PLZ4 = 3053 if bürgerort == "Deisswil bei Münchenbuchsee" & PLZBmiss == 1
replace B_PLZ4 = 1568 if bürgerort == "Delley-Portalban" & PLZBmiss == 1
replace B_PLZ4 = 2800 if bürgerort == "Delémont" & PLZBmiss == 1
replace B_PLZ4 = 1410 if bürgerort == "Denezy" & PLZBmiss == 1
replace B_PLZ4 = 1026 if bürgerort == "Denges" & PLZBmiss == 1
replace B_PLZ4 = 5026 if bürgerort == "Densbüren" & PLZBmiss == 1
replace B_PLZ4 = 4551 if bürgerort == "Derendingen" & PLZBmiss == 1
replace B_PLZ4 = 2802 if bürgerort == "Develier" & PLZBmiss == 1
replace B_PLZ4 = 8157 if bürgerort == "Dielsdorf" & PLZBmiss == 1
replace B_PLZ4 = 3756 if bürgerort == "Diemtigen" & PLZBmiss == 1
replace B_PLZ4 = 9444 if bürgerort == "Diepoldsau" & PLZBmiss == 1
replace B_PLZ4 = 6036 if bürgerort == "Dierikon" & PLZBmiss == 1
replace B_PLZ4 = 3264 if bürgerort == "Diessbach bei Büren" & PLZBmiss == 1
replace B_PLZ4 = 8253 if bürgerort == "Diessenhofen" & PLZBmiss == 1
replace B_PLZ4 = 8953 if bürgerort == "Dietikon" & PLZBmiss == 1
replace B_PLZ4 = 8305 if bürgerort == "Dietlikon" & PLZBmiss == 1
replace B_PLZ4 = 8474 if bürgerort == "Dinhard" & PLZBmiss == 1
replace B_PLZ4 = 7182 if bürgerort == "Disentis/Mustèr" & PLZBmiss == 1
replace B_PLZ4 = 4243 if bürgerort == "Dittingen" & PLZBmiss == 1
replace B_PLZ4 = 1304 if bürgerort == "Dizy" & PLZBmiss == 1
replace B_PLZ4 = 7013 if bürgerort == "Domat/Ems" & PLZBmiss == 1
replace B_PLZ4 = 2056 if bürgerort == "Dombresson" & PLZBmiss == 1
replace B_PLZ4 = 1564 if bürgerort == "Domdidier" & PLZBmiss == 1
replace B_PLZ4 = 7418 if bürgerort == "Domleschg" & PLZBmiss == 1
replace B_PLZ4 = 7433 if bürgerort == "Donath" & PLZBmiss == 1
replace B_PLZ4 = 1408 if bürgerort == "Donneloye" & PLZBmiss == 1
replace B_PLZ4 = 8458 if bürgerort == "Dorf" & PLZBmiss == 1
replace B_PLZ4 = 4143 if bürgerort == "Dornach" & PLZBmiss == 1
replace B_PLZ4 = 5605 if bürgerort == "Dottikon" & PLZBmiss == 1
replace B_PLZ4 = 3293 if bürgerort == "Dotzigen" & PLZBmiss == 1
replace B_PLZ4 = 8582 if bürgerort == "Dozwil" & PLZBmiss == 1
replace B_PLZ4 = 4558 if bürgerort == "Drei Höfe" & PLZBmiss == 1
replace B_PLZ4 = 4202 if bürgerort == "Duggingen" & PLZBmiss == 1
replace B_PLZ4 = 1266 if bürgerort == "Duillier" & PLZBmiss == 1
replace B_PLZ4 = 4657 if bürgerort == "Dulliken" & PLZBmiss == 1
replace B_PLZ4 = 1195 if bürgerort == "Dully" & PLZBmiss == 1
replace B_PLZ4 = 8108 if bürgerort == "Dällikon" & PLZBmiss == 1
replace B_PLZ4 = 4658 if bürgerort == "Däniken" & PLZBmiss == 1
replace B_PLZ4 = 8114 if bürgerort == "Dänikon" & PLZBmiss == 1
replace B_PLZ4 = 5312 if bürgerort == "Döttingen" & PLZBmiss == 1
replace B_PLZ4 = 8600 if bürgerort == "Dübendorf" & PLZBmiss == 1
replace B_PLZ4 = 3186 if bürgerort == "Düdingen" & PLZBmiss == 1
replace B_PLZ4 = 8632 if bürgerort == "Dürnten" & PLZBmiss == 1
replace B_PLZ4 = 3465 if bürgerort == "Dürrenroth" & PLZBmiss == 1
replace B_PLZ4 = 6245 if bürgerort == "Ebersecken" & PLZBmiss == 1
replace B_PLZ4 = 6031 if bürgerort == "Ebikon" & PLZBmiss == 1
replace B_PLZ4 = 9642 if bürgerort == "Ebnat-Kappel" & PLZBmiss == 1
replace B_PLZ4 = 1040 if bürgerort == "Echallens" & PLZBmiss == 1
replace B_PLZ4 = 1646 if bürgerort == "Echarlens" & PLZBmiss == 1
replace B_PLZ4 = 1024 if bürgerort == "Echichens" & PLZBmiss == 1
replace B_PLZ4 = 1311 if bürgerort == "Eclépens" & PLZBmiss == 1
replace B_PLZ4 = 5078 if bürgerort == "Effingen" & PLZBmiss == 1
replace B_PLZ4 = 8132 if bürgerort == "Egg" & PLZBmiss == 1
replace B_PLZ4 = 3939 if bürgerort == "Eggerberg" & PLZBmiss == 1
replace B_PLZ4 = 9036 if bürgerort == "Eggersriet" & PLZBmiss == 1
replace B_PLZ4 = 3536 if bürgerort == "Eggiwil" & PLZBmiss == 1
replace B_PLZ4 = 8193 if bürgerort == "Eglisau" & PLZBmiss == 1
replace B_PLZ4 = 5704 if bürgerort == "Egliswil" & PLZBmiss == 1
replace B_PLZ4 = 9315 if bürgerort == "Egnach" & PLZBmiss == 1
replace B_PLZ4 = 6243 if bürgerort == "Egolzwil" & PLZBmiss == 1
replace B_PLZ4 = 5420 if bürgerort == "Ehrendingen" & PLZBmiss == 1
replace B_PLZ4 = 6205 if bürgerort == "Eich" & PLZBmiss == 1
replace B_PLZ4 = 9452 if bürgerort == "Eichberg" & PLZBmiss == 1
replace B_PLZ4 = 5074 if bürgerort == "Eiken" & PLZBmiss == 1
replace B_PLZ4 = 8840 if bürgerort == "Einsiedeln" & PLZBmiss == 1
replace B_PLZ4 = 3943 if bürgerort == "Eischoll" & PLZBmiss == 1
replace B_PLZ4 = 8353 if bürgerort == "Elgg" & PLZBmiss == 1
replace B_PLZ4 = 8548 if bürgerort == "Ellikon an der Thur" & PLZBmiss == 1
replace B_PLZ4 = 8767 if bürgerort == "Elm" & PLZBmiss == 1
replace B_PLZ4 = 8352 if bürgerort == "Elsau" & PLZBmiss == 1
replace B_PLZ4 = 3922 if bürgerort == "Embd" & PLZBmiss == 1
replace B_PLZ4 = 8424 if bürgerort == "Embrach" & PLZBmiss == 1
replace B_PLZ4 = 6020 if bürgerort == "Emmen" & PLZBmiss == 1
replace B_PLZ4 = 6376 if bürgerort == "Emmetten" & PLZBmiss == 1
replace B_PLZ4 = 5304 if bürgerort == "Endingen" & PLZBmiss == 1
replace B_PLZ4 = 6390 if bürgerort == "Engelberg" & PLZBmiss == 1
replace B_PLZ4 = 8765 if bürgerort == "Engi" & PLZBmiss == 1
replace B_PLZ4 = 8755 if bürgerort == "Ennenda" & PLZBmiss == 1
replace B_PLZ4 = 5408 if bürgerort == "Ennetbaden" & PLZBmiss == 1
replace B_PLZ4 = 6373 if bürgerort == "Ennetbürgen" & PLZBmiss == 1
replace B_PLZ4 = 1667 if bürgerort == "Enney" & PLZBmiss == 1
replace B_PLZ4 = 6162 if bürgerort == "Entlebuch" & PLZBmiss == 1
replace B_PLZ4 = 1066 if bürgerort == "Epalinges" & PLZBmiss == 1
replace B_PLZ4 = 5012 if bürgerort == "Eppenberg-Wöschnau" & PLZBmiss == 1
replace B_PLZ4 = 3272 if bürgerort == "Epsach" & PLZBmiss == 1
replace B_PLZ4 = 3947 if bürgerort == "Ergisch" & PLZBmiss == 1
replace B_PLZ4 = 4952 if bürgerort == "Eriswil" & PLZBmiss == 1
replace B_PLZ4 = 3619 if bürgerort == "Eriz" & PLZBmiss == 1
replace B_PLZ4 = 3235 if bürgerort == "Erlach" & PLZBmiss == 1
replace B_PLZ4 = 3758 if bürgerort == "Erlenbach im Simmental" & PLZBmiss == 1
replace B_PLZ4 = 5017 if bürgerort == "Erlinsbach" & PLZBmiss == 1
replace B_PLZ4 = 8273 if bürgerort == "Ermatingen" & PLZBmiss == 1
replace B_PLZ4 = 8726 if bürgerort == "Ernetschwil" & PLZBmiss == 1
replace B_PLZ4 = 3957 if bürgerort == "Erschmatt" & PLZBmiss == 1
replace B_PLZ4 = 4228 if bürgerort == "Erschwil" & PLZBmiss == 1
replace B_PLZ4 = 3424 if bürgerort == "Ersigen" & PLZBmiss == 1
replace B_PLZ4 = 6472 if bürgerort == "Erstfeld" & PLZBmiss == 1
replace B_PLZ4 = 8264 if bürgerort == "Eschenz" & PLZBmiss == 1
replace B_PLZ4 = 8360 if bürgerort == "Eschlikon" & PLZBmiss == 1
replace B_PLZ4 = 6192 if bürgerort == "Escholzmatt" & PLZBmiss == 1
replace B_PLZ4 = 6182 if bürgerort == "Escholzmatt-Marbach" & PLZBmiss == 1
replace B_PLZ4 = 1670 if bürgerort == "Esmonts" & PLZBmiss == 1
replace B_PLZ4 = 1186 if bürgerort == "Essertines-sur-Rolle" & PLZBmiss == 1
replace B_PLZ4 = 1417 if bürgerort == "Essertines-sur-Yverdon" & PLZBmiss == 1
replace B_PLZ4 = 1489 if bürgerort == "Estavayer" & PLZBmiss == 1
replace B_PLZ4 = 1695 if bürgerort == "Estavayer-le-Gibloux" & PLZBmiss == 1
replace B_PLZ4 = 1473 if bürgerort == "Estavayer-le-Lac" & PLZBmiss == 1
replace B_PLZ4 = 1687 if bürgerort == "Estévenens" & PLZBmiss == 1
replace B_PLZ4 = 1163 if bürgerort == "Etoy" & PLZBmiss == 1
replace B_PLZ4 = 6217 if bürgerort == "Ettiswil" & PLZBmiss == 1
replace B_PLZ4 = 3306 if bürgerort == "Etzelkofen" & PLZBmiss == 1
replace B_PLZ4 = 1985 if bürgerort == "Evolène" & PLZBmiss == 1
replace B_PLZ4 = 1262 if bürgerort == "Eysins" & PLZBmiss == 1
replace B_PLZ4 = 3617 if bürgerort == "Fahrni" & PLZBmiss == 1
replace B_PLZ4 = 5615 if bürgerort == "Fahrwangen" & PLZBmiss == 1
replace B_PLZ4 = 2916 if bürgerort == "Fahy" & PLZBmiss == 1
replace B_PLZ4 = 6763 if bürgerort == "Faido" & PLZBmiss == 1
replace B_PLZ4 = 1726 if bürgerort == "Farvagny-le-Grand" & PLZBmiss == 1
replace B_PLZ4 = 8320 if bürgerort == "Fehraltorf" & PLZBmiss == 1
replace B_PLZ4 = 8552 if bürgerort == "Felben" & PLZBmiss == 1
replace B_PLZ4 = 7153 if bürgerort == "Fellers" & PLZBmiss == 1
replace B_PLZ4 = 2063 if bürgerort == "Fenin-Vilars-Saules" & PLZBmiss == 1
replace B_PLZ4 = 3206 if bürgerort == "Ferenbalm" & PLZBmiss == 1
replace B_PLZ4 = 7444 if bürgerort == "Ferrera" & PLZBmiss == 1
replace B_PLZ4 = 8246 if bürgerort == "Feuerthalen" & PLZBmiss == 1
replace B_PLZ4 = 8834 if bürgerort == "Feusisberg" & PLZBmiss == 1
replace B_PLZ4 = 7235 if bürgerort == "Fideris" & PLZBmiss == 1
replace B_PLZ4 = 3984 if bürgerort == "Fieschertal" & PLZBmiss == 1
replace B_PLZ4 = 7477 if bürgerort == "Filisur" & PLZBmiss == 1
replace B_PLZ4 = 8757 if bürgerort == "Filzbach" & PLZBmiss == 1
replace B_PLZ4 = 1925 if bürgerort == "Finhaut" & PLZBmiss == 1
replace B_PLZ4 = 6145 if bürgerort == "Fischbach" & PLZBmiss == 1
replace B_PLZ4 = 8498 if bürgerort == "Fischenthal" & PLZBmiss == 1
replace B_PLZ4 = 8376 if bürgerort == "Fischingen" & PLZBmiss == 1
replace B_PLZ4 = 5442 if bürgerort == "Fislisbach" & PLZBmiss == 1
replace B_PLZ4 = 8416 if bürgerort == "Flaach" & PLZBmiss == 1
replace B_PLZ4 = 9230 if bürgerort == "Flawil" & PLZBmiss == 1
replace B_PLZ4 = 2114 if bürgerort == "Fleurier" & PLZBmiss == 1
replace B_PLZ4 = 7017 if bürgerort == "Flims" & PLZBmiss == 1
replace B_PLZ4 = 8894 if bürgerort == "Flums" & PLZBmiss == 1
replace B_PLZ4 = 8247 if bürgerort == "Flurlingen" & PLZBmiss == 1
replace B_PLZ4 = 7306 if bürgerort == "Fläsch" & PLZBmiss == 1
replace B_PLZ4 = 6454 if bürgerort == "Flüelen" & PLZBmiss == 1
replace B_PLZ4 = 6173 if bürgerort == "Flühli" & PLZBmiss == 1
replace B_PLZ4 = 1423 if bürgerort == "Fontanezier" & PLZBmiss == 1
replace B_PLZ4 = 3636 if bürgerort == "Forst" & PLZBmiss == 1
replace B_PLZ4 = 3312 if bürgerort == "Fraubrunnen" & PLZBmiss == 1
replace B_PLZ4 = 8500 if bürgerort == "Frauenfeld" & PLZBmiss == 1
replace B_PLZ4 = 8808 if bürgerort == "Freienbach" & PLZBmiss == 1
replace B_PLZ4 = 8428 if bürgerort == "Freienstein-Teufen" & PLZBmiss == 1
replace B_PLZ4 = 5423 if bürgerort == "Freienwil" & PLZBmiss == 1
replace B_PLZ4 = 3510 if bürgerort == "Freimettigen" & PLZBmiss == 1
replace B_PLZ4 = 2027 if bürgerort == "Fresens" & PLZBmiss == 1
replace B_PLZ4 = 1700 if bürgerort == "Fribourg" & PLZBmiss == 1
replace B_PLZ4 = 5070 if bürgerort == "Frick" & PLZBmiss == 1
replace B_PLZ4 = 1055 if bürgerort == "Froideville" & PLZBmiss == 1
replace B_PLZ4 = 3714 if bürgerort == "Frutigen" & PLZBmiss == 1
replace B_PLZ4 = 7533 if bürgerort == "Fuldera" & PLZBmiss == 1
replace B_PLZ4 = 1926 if bürgerort == "Fully" & PLZBmiss == 1
replace B_PLZ4 = 6696 if bürgerort == "Fusio" & PLZBmiss == 1
replace B_PLZ4 = 8118 if bürgerort == "Fällanden" & PLZBmiss == 1
replace B_PLZ4 = 1173 if bürgerort == "Féchy" & PLZBmiss == 1
replace B_PLZ4 = 1532 if bürgerort == "Fétigny" & PLZBmiss == 1
replace B_PLZ4 = 4414 if bürgerort == "Füllinsdorf" & PLZBmiss == 1
replace B_PLZ4 = 8547 if bürgerort == "Gachnang" & PLZBmiss == 1
replace B_PLZ4 = 3863 if bürgerort == "Gadmen" & PLZBmiss == 1
replace B_PLZ4 = 9056 if bürgerort == "Gais" & PLZBmiss == 1
replace B_PLZ4 = 9030 if bürgerort == "Gaiserwald" & PLZBmiss == 1
replace B_PLZ4 = 8854 if bürgerort == "Galgenen" & PLZBmiss == 1
replace B_PLZ4 = 6579 if bürgerort == "Gambarogno" & PLZBmiss == 1
replace B_PLZ4 = 3945 if bürgerort == "Gampel" & PLZBmiss == 1
replace B_PLZ4 = 3945 if bürgerort == "Gampel-Bratsch" & PLZBmiss == 1
replace B_PLZ4 = 9473 if bürgerort == "Gams" & PLZBmiss == 1
replace B_PLZ4 = 6978 if bürgerort == "Gandria" & PLZBmiss == 1
replace B_PLZ4 = 5272 if bürgerort == "Gansingen" & PLZBmiss == 1
replace B_PLZ4 = 9608 if bürgerort == "Ganterschwil" & PLZBmiss == 1
replace B_PLZ4 = 5412 if bürgerort == "Gebenstorf" & PLZBmiss == 1
replace B_PLZ4 = 6284 if bürgerort == "Gelfingen" & PLZBmiss == 1
replace B_PLZ4 = 4460 if bürgerort == "Gelterkinden" & PLZBmiss == 1
replace B_PLZ4 = 1294 if bürgerort == "Genthod" & PLZBmiss == 1
replace B_PLZ4 = 6925 if bürgerort == "Gentilino" & PLZBmiss == 1
replace B_PLZ4 = 1200 if bürgerort == "Genève" & PLZBmiss == 1
replace B_PLZ4 = 4563 if bürgerort == "Gerlafingen" & PLZBmiss == 1
replace B_PLZ4 = 8954 if bürgerort == "Geroldswil" & PLZBmiss == 1
replace B_PLZ4 = 6410 if bürgerort == "Gersau" & PLZBmiss == 1
replace B_PLZ4 = 3115 if bürgerort == "Gerzensee" & PLZBmiss == 1
replace B_PLZ4 = 6142 if bürgerort == "Gettnau" & PLZBmiss == 1
replace B_PLZ4 = 6232 if bürgerort == "Geuensee" & PLZBmiss == 1
replace B_PLZ4 = 1728 if bürgerort == "Gibloux" & PLZBmiss == 1
replace B_PLZ4 = 1735 if bürgerort == "Giffers" & PLZBmiss == 1
replace B_PLZ4 = 1182 if bürgerort == "Gilly" & PLZBmiss == 1
replace B_PLZ4 = 1188 if bürgerort == "Gimel" & PLZBmiss == 1
replace B_PLZ4 = 1276 if bürgerort == "Gingins" & PLZBmiss == 1
replace B_PLZ4 = 5073 if bürgerort == "Gipf-Oberfrick" & PLZBmiss == 1
replace B_PLZ4 = 6038 if bürgerort == "Gisikon" & PLZBmiss == 1
replace B_PLZ4 = 6074 if bürgerort == "Giswil" & PLZBmiss == 1
replace B_PLZ4 = 1762 if bürgerort == "Givisiez" & PLZBmiss == 1
replace B_PLZ4 = 1196 if bürgerort == "Gland" & PLZBmiss == 1
replace B_PLZ4 = 8750 if bürgerort == "Glarus" & PLZBmiss == 1
replace B_PLZ4 = 8752 if bürgerort == "Glarus Nord" & PLZBmiss == 1
replace B_PLZ4 = 8756 if bürgerort == "Glarus Süd" & PLZBmiss == 1
replace B_PLZ4 = 8192 if bürgerort == "Glattfelden" & PLZBmiss == 1
replace B_PLZ4 = 1544 if bürgerort == "Gletterens" & PLZBmiss == 1
replace B_PLZ4 = 9403 if bürgerort == "Goldach" & PLZBmiss == 1
replace B_PLZ4 = 8638 if bürgerort == "Goldingen" & PLZBmiss == 1
replace B_PLZ4 = 1124 if bürgerort == "Gollion" & PLZBmiss == 1
replace B_PLZ4 = 8738 if bürgerort == "Gommiswald" & PLZBmiss == 1
replace B_PLZ4 = 3989 if bürgerort == "Goms" & PLZBmiss == 1
replace B_PLZ4 = 4955 if bürgerort == "Gondiswil" & PLZBmiss == 1
replace B_PLZ4 = 5728 if bürgerort == "Gontenschwil" & PLZBmiss == 1
replace B_PLZ4 = 6672 if bürgerort == "Gordevio" & PLZBmiss == 1
replace B_PLZ4 = 6518 if bürgerort == "Gorduno" & PLZBmiss == 1
replace B_PLZ4 = 2025 if bürgerort == "Gorgier" & PLZBmiss == 1
replace B_PLZ4 = 8274 if bürgerort == "Gottlieben" & PLZBmiss == 1
replace B_PLZ4 = 9213 if bürgerort == "Gottshaus" & PLZBmiss == 1
replace B_PLZ4 = 1376 if bürgerort == "Goumoens-la-Ville" & PLZBmiss == 1
replace B_PLZ4 = 3376 if bürgerort == "Graben" & PLZBmiss == 1
replace B_PLZ4 = 9470 if bürgerort == "Grabs" & PLZBmiss == 1
replace B_PLZ4 = 3989 if bürgerort == "Grafschaft" & PLZBmiss == 1
replace B_PLZ4 = 1543 if bürgerort == "Grandcour" & PLZBmiss == 1
replace B_PLZ4 = 1422 if bürgerort == "Grandson" & PLZBmiss == 1
replace B_PLZ4 = 1091 if bürgerort == "Grandvaux" & PLZBmiss == 1
replace B_PLZ4 = 1763 if bürgerort == "Granges-Paccot" & PLZBmiss == 1
replace B_PLZ4 = 1686 if bürgerort == "Grangettes" & PLZBmiss == 1
replace B_PLZ4 = 8606 if bürgerort == "Greifensee" & PLZBmiss == 1
replace B_PLZ4 = 4203 if bürgerort == "Grellingen" & PLZBmiss == 1
replace B_PLZ4 = 2540 if bürgerort == "Grenchen" & PLZBmiss == 1
replace B_PLZ4 = 3993 if bürgerort == "Grengiols" & PLZBmiss == 1
replace B_PLZ4 = 6404 if bürgerort == "Greppen" & PLZBmiss == 1
replace B_PLZ4 = 1432 if bürgerort == "Gressy" & PLZBmiss == 1
replace B_PLZ4 = 5014 if bürgerort == "Gretzenbach" & PLZBmiss == 1
replace B_PLZ4 = 1971 if bürgerort == "Grimisuat" & PLZBmiss == 1
replace B_PLZ4 = 3823 if bürgerort == "Grindelwald" & PLZBmiss == 1
replace B_PLZ4 = 3257 if bürgerort == "Grossaffoltern" & PLZBmiss == 1
replace B_PLZ4 = 6146 if bürgerort == "Grossdietwil" & PLZBmiss == 1
replace B_PLZ4 = 3506 if bürgerort == "Grosshöchstetten" & PLZBmiss == 1
replace B_PLZ4 = 6022 if bürgerort == "Grosswangen" & PLZBmiss == 1
replace B_PLZ4 = 1882 if bürgerort == "Gryon" & PLZBmiss == 1
replace B_PLZ4 = 5722 if bürgerort == "Gränichen" & PLZBmiss == 1
replace B_PLZ4 = 3979 if bürgerort == "Grône" & PLZBmiss == 1
replace B_PLZ4 = 8627 if bürgerort == "Grüningen" & PLZBmiss == 1
replace B_PLZ4 = 7213 if bürgerort == "Grüsch" & PLZBmiss == 1
replace B_PLZ4 = 3784 if bürgerort == "Gsteig" & PLZBmiss == 1
replace B_PLZ4 = 1738 if bürgerort == "Guggisberg" & PLZBmiss == 1
replace B_PLZ4 = 6222 if bürgerort == "Gunzwil" & PLZBmiss == 1
replace B_PLZ4 = 3208 if bürgerort == "Gurbrü" & PLZBmiss == 1
replace B_PLZ4 = 6482 if bürgerort == "Gurtnellen" & PLZBmiss == 1
replace B_PLZ4 = 3663 if bürgerort == "Gurzelen" & PLZBmiss == 1
replace B_PLZ4 = 3956 if bürgerort == "Guttet-Feschel" & PLZBmiss == 1
replace B_PLZ4 = 1251 if bürgerort == "Gy" & PLZBmiss == 1
replace B_PLZ4 = 8214 if bürgerort == "Gächlingen" & PLZBmiss == 1
replace B_PLZ4 = 8594 if bürgerort == "Güttingen" & PLZBmiss == 1
replace B_PLZ4 = 3804 if bürgerort == "Habkern" & PLZBmiss == 1
replace B_PLZ4 = 8523 if bürgerort == "Hagenbuch" & PLZBmiss == 1
replace B_PLZ4 = 7023 if bürgerort == "Haldenstein" & PLZBmiss == 1
replace B_PLZ4 = 8215 if bürgerort == "Hallau" & PLZBmiss == 1
replace B_PLZ4 = 5705 if bürgerort == "Hallwil" & PLZBmiss == 1
replace B_PLZ4 = 3419 if bürgerort == "Hasle bei Burgdorf" & PLZBmiss == 1
replace B_PLZ4 = 8773 if bürgerort == "Haslen" & PLZBmiss == 1
replace B_PLZ4 = 6083 if bürgerort == "Hasliberg" & PLZBmiss == 1
replace B_PLZ4 = 8925 if bürgerort == "Hausen am Albis" & PLZBmiss == 1
replace B_PLZ4 = 1669 if bürgerort == "Haut-Intyamon" & PLZBmiss == 1
replace B_PLZ4 = 2855 if bürgerort == "Haute-Sorne" & PLZBmiss == 1
replace B_PLZ4 = 1648 if bürgerort == "Hauteville" & PLZBmiss == 1
replace B_PLZ4 = 8580 if bürgerort == "Hefenhofen" & PLZBmiss == 1
replace B_PLZ4 = 9410 if bürgerort == "Heiden" & PLZBmiss == 1
replace B_PLZ4 = 3625 if bürgerort == "Heiligenschwendi" & PLZBmiss == 1
replace B_PLZ4 = 9514 if bürgerort == "Heiligkreuz" & PLZBmiss == 1
replace B_PLZ4 = 3412 if bürgerort == "Heimiswil" & PLZBmiss == 1
replace B_PLZ4 = 1714 if bürgerort == "Heitenried" & PLZBmiss == 1
replace B_PLZ4 = 9633 if bürgerort == "Hemberg" & PLZBmiss == 1
replace B_PLZ4 = 8231 if bürgerort == "Hemmental" & PLZBmiss == 1
replace B_PLZ4 = 4715 if bürgerort == "Herbetswil" & PLZBmiss == 1
replace B_PLZ4 = 8535 if bürgerort == "Herdern" & PLZBmiss == 1
replace B_PLZ4 = 6133 if bürgerort == "Hergiswil bei Willisau" & PLZBmiss == 1
replace B_PLZ4 = 9102 if bürgerort == "Herisau" & PLZBmiss == 1
replace B_PLZ4 = 5027 if bürgerort == "Herznach" & PLZBmiss == 1
replace B_PLZ4 = 4577 if bürgerort == "Hessigkofen" & PLZBmiss == 1
replace B_PLZ4 = 8442 if bürgerort == "Hettlingen" & PLZBmiss == 1
replace B_PLZ4 = 6024 if bürgerort == "Hildisrieden" & PLZBmiss == 1
replace B_PLZ4 = 8340 if bürgerort == "Hinwil" & PLZBmiss == 1
replace B_PLZ4 = 8335 if bürgerort == "Hittnau" & PLZBmiss == 1
replace B_PLZ4 = 6281 if bürgerort == "Hochdorf" & PLZBmiss == 1
replace B_PLZ4 = 8182 if bürgerort == "Hochfelden" & PLZBmiss == 1
replace B_PLZ4 = 8242 if bürgerort == "Hofen" & PLZBmiss == 1
replace B_PLZ4 = 3858 if bürgerort == "Hofstetten bei Brienz" & PLZBmiss == 1
replace B_PLZ4 = 4112 if bürgerort == "Hofstetten-Flüh" & PLZBmiss == 1
replace B_PLZ4 = 6276 if bürgerort == "Hohenrain" & PLZBmiss == 1
replace B_PLZ4 = 2000 if bürgerort == "Hohentannen" & PLZBmiss == 1
replace B_PLZ4 = 3949 if bürgerort == "Hohtenn" & PLZBmiss == 1
replace B_PLZ4 = 3622 if bürgerort == "Homberg" & PLZBmiss == 1
replace B_PLZ4 = 8634 if bürgerort == "Hombrechtikon" & PLZBmiss == 1
replace B_PLZ4 = 8507 if bürgerort == "Homburg" & PLZBmiss == 1
replace B_PLZ4 = 8135 if bürgerort == "Horgen" & PLZBmiss == 1
replace B_PLZ4 = 5075 if bürgerort == "Hornussen" & PLZBmiss == 1
replace B_PLZ4 = 3623 if bürgerort == "Horrenbach-Buchen" & PLZBmiss == 1
replace B_PLZ4 = 4557 if bürgerort == "Horriwil" & PLZBmiss == 1
replace B_PLZ4 = 6048 if bürgerort == "Horw" & PLZBmiss == 1
replace B_PLZ4 = 8457 if bürgerort == "Humlikon" & PLZBmiss == 1
replace B_PLZ4 = 9064 if bürgerort == "Hundwil" & PLZBmiss == 1
replace B_PLZ4 = 4950 if bürgerort == "Huttwil" & PLZBmiss == 1
replace B_PLZ4 = 4615 if bürgerort == "Hägendorf" & PLZBmiss == 1
replace B_PLZ4 = 9308 if bürgerort == "Häggenschwil" & PLZBmiss == 1
replace B_PLZ4 = 5607 if bürgerort == "Hägglingen" & PLZBmiss == 1
replace B_PLZ4 = 4640 if bürgerort == "Härkingen" & PLZBmiss == 1
replace B_PLZ4 = 3510 if bürgerort == "Häutligen" & PLZBmiss == 1
replace B_PLZ4 = 1987 if bürgerort == "Hérémence" & PLZBmiss == 1
replace B_PLZ4 = 3429 if bürgerort == "Höchstetten" & PLZBmiss == 1
replace B_PLZ4 = 3631 if bürgerort == "Höfen" & PLZBmiss == 1
replace B_PLZ4 = 4434 if bürgerort == "Hölstein" & PLZBmiss == 1
replace B_PLZ4 = 8194 if bürgerort == "Hüntwangen" & PLZBmiss == 1
replace B_PLZ4 = 8825 if bürgerort == "Hütten" & PLZBmiss == 1
replace B_PLZ4 = 8115 if bürgerort == "Hüttikon" & PLZBmiss == 1
replace B_PLZ4 = 1977 if bürgerort == "Icogne" & PLZBmiss == 1
replace B_PLZ4 = 7302 if bürgerort == "Igis" & PLZBmiss == 1
replace B_PLZ4 = 7130 if bürgerort == "Ilanz/Glion" & PLZBmiss == 1
replace B_PLZ4 = 6434 if bürgerort == "Illgau" & PLZBmiss == 1
replace B_PLZ4 = 8307 if bürgerort == "Illnau-Effretikon" & PLZBmiss == 1
replace B_PLZ4 = 6440 if bürgerort == "Ingenbohl" & PLZBmiss == 1
replace B_PLZ4 = 3232 if bürgerort == "Ins" & PLZBmiss == 1
replace B_PLZ4 = 6034 if bürgerort == "Inwil" & PLZBmiss == 1
replace B_PLZ4 = 6461 if bürgerort == "Isenthal" & PLZBmiss == 1
replace B_PLZ4 = 6993 if bürgerort == "Iseo" & PLZBmiss == 1
replace B_PLZ4 = 6661 if bürgerort == "Isorno" & PLZBmiss == 1
replace B_PLZ4 = 1914 if bürgerort == "Isérables" & PLZBmiss == 1
replace B_PLZ4 = 4452 if bürgerort == "Itingen" & PLZBmiss == 1
replace B_PLZ4 = 3063 if bürgerort == "Ittigen" & PLZBmiss == 1
replace B_PLZ4 = 1656 if bürgerort == "Jaun" & PLZBmiss == 1
replace B_PLZ4 = 7233 if bürgerort == "Jenaz" & PLZBmiss == 1
replace B_PLZ4 = 8640 if bürgerort == "Jona" & PLZBmiss == 1
replace B_PLZ4 = 8916 if bürgerort == "Jonen" & PLZBmiss == 1
replace B_PLZ4 = 9536 if bürgerort == "Jonschwil" & PLZBmiss == 1
replace B_PLZ4 = 1061 if bürgerort == "Jorat-Menthue" & PLZBmiss == 1
replace B_PLZ4 = 1083 if bürgerort == "Jorat-Mézières" & PLZBmiss == 1
replace B_PLZ4 = 1008 if bürgerort == "Jouxtens-Mézery" & PLZBmiss == 1
replace B_PLZ4 = 1254 if bürgerort == "Jussy" & PLZBmiss == 1
replace B_PLZ4 = 4303 if bürgerort == "Kaiseraugst" & PLZBmiss == 1
replace B_PLZ4 = 5466 if bürgerort == "Kaiserstuhl" & PLZBmiss == 1
replace B_PLZ4 = 5625 if bürgerort == "Kallern" & PLZBmiss == 1
replace B_PLZ4 = 8722 if bürgerort == "Kaltbrunn" & PLZBmiss == 1
replace B_PLZ4 = 4535 if bürgerort == "Kammersrohr" & PLZBmiss == 1
replace B_PLZ4 = 3717 if bürgerort == "Kandergrund" & PLZBmiss == 1
replace B_PLZ4 = 8926 if bürgerort == "Kappel am Albis" & PLZBmiss == 1
replace B_PLZ4 = 3122 if bürgerort == "Kehrsatz" & PLZBmiss == 1
replace B_PLZ4 = 8573 if bürgerort == "Kemmental" & PLZBmiss == 1
replace B_PLZ4 = 3309 if bürgerort == "Kernenried" & PLZBmiss == 1
replace B_PLZ4 = 6064 if bürgerort == "Kerns" & PLZBmiss == 1
replace B_PLZ4 = 3210 if bürgerort == "Kerzers" & PLZBmiss == 1
replace B_PLZ4 = 8593 if bürgerort == "Kesswil" & PLZBmiss == 1
replace B_PLZ4 = 4703 if bürgerort == "Kestenholz" & PLZBmiss == 1
replace B_PLZ4 = 4468 if bürgerort == "Kienberg" & PLZBmiss == 1
replace B_PLZ4 = 3917 if bürgerort == "Kippel" & PLZBmiss == 1
replace B_PLZ4 = 5054 if bürgerort == "Kirchleerau" & PLZBmiss == 1
replace B_PLZ4 = 3038 if bürgerort == "Kirchlindach" & PLZBmiss == 1
replace B_PLZ4 = 4245 if bürgerort == "Kleinlützel" & PLZBmiss == 1
replace B_PLZ4 = 5313 if bürgerort == "Klingnau" & PLZBmiss == 1
replace B_PLZ4 = 7247 if bürgerort == "Klosters" & PLZBmiss == 1
replace B_PLZ4 = 8058 if bürgerort == "Kloten" & PLZBmiss == 1
replace B_PLZ4 = 8934 if bürgerort == "Knonau" & PLZBmiss == 1
replace B_PLZ4 = 6213 if bürgerort == "Knutwil" & PLZBmiss == 1
replace B_PLZ4 = 5322 if bürgerort == "Koblenz" & PLZBmiss == 1
replace B_PLZ4 = 3510 if bürgerort == "Konolfingen" & PLZBmiss == 1
replace B_PLZ4 = 6217 if bürgerort == "Kottwil" & PLZBmiss == 1
replace B_PLZ4 = 9214 if bürgerort == "Kradolf-Schönenberg" & PLZBmiss == 1
replace B_PLZ4 = 3704 if bürgerort == "Krattigen" & PLZBmiss == 1
replace B_PLZ4 = 3325 if bürgerort == "Krauchthal" & PLZBmiss == 1
replace B_PLZ4 = 8280 if bürgerort == "Kreuzlingen" & PLZBmiss == 1
replace B_PLZ4 = 3179 if bürgerort == "Kriechenwil" & PLZBmiss == 1
replace B_PLZ4 = 4566 if bürgerort == "Kriegstetten" & PLZBmiss == 1
replace B_PLZ4 = 6010 if bürgerort == "Kriens" & PLZBmiss == 1
replace B_PLZ4 = 9622 if bürgerort == "Krinau" & PLZBmiss == 1
replace B_PLZ4 = 9643 if bürgerort == "Krummenau" & PLZBmiss == 1
replace B_PLZ4 = 4447 if bürgerort == "Känerkinden" & PLZBmiss == 1
replace B_PLZ4 = 5742 if bürgerort == "Kölliken" & PLZBmiss == 1
replace B_PLZ4 = 3098 if bürgerort == "Köniz" & PLZBmiss == 1
replace B_PLZ4 = 7240 if bürgerort == "Küblis" & PLZBmiss == 1
replace B_PLZ4 = 5444 if bürgerort == "Künten" & PLZBmiss == 1
replace B_PLZ4 = 6403 if bürgerort == "Küssnacht am Rigi" & PLZBmiss == 1
replace B_PLZ4 = 5024 if bürgerort == "Küttigen" & PLZBmiss == 1
replace B_PLZ4 = 1342 if bürgerort == "L'Abbaye" & PLZBmiss == 1
replace B_PLZ4 = 1355 if bürgerort == "L'Abergement" & PLZBmiss == 1
replace B_PLZ4 = 1148 if bürgerort == "L'Isle" & PLZBmiss == 1
replace B_PLZ4 = 2947 if bürgerort == "La Baroche" & PLZBmiss == 1
replace B_PLZ4 = 1756 if bürgerort == "La Brillaz" & PLZBmiss == 1
replace B_PLZ4 = 2406 if bürgerort == "La Brévine" & PLZBmiss == 1
replace B_PLZ4 = 2300 if bürgerort == "La Chaux-de-Fonds" & PLZBmiss == 1
replace B_PLZ4 = 2405 if bürgerort == "La Chaux-du-Milieu" & PLZBmiss == 1
replace B_PLZ4 = 2117 if bürgerort == "La Côte-aux-Fées" & PLZBmiss == 1
replace B_PLZ4 = 2520 if bürgerort == "La Neuveville" & PLZBmiss == 1
replace B_PLZ4 = 1634 if bürgerort == "La Roche" & PLZBmiss == 1
replace B_PLZ4 = 2314 if bürgerort == "La Sagne" & PLZBmiss == 1
replace B_PLZ4 = 1814 if bürgerort == "La Tour-de-Peilz" & PLZBmiss == 1
replace B_PLZ4 = 1635 if bürgerort == "La Tour-de-Trême" & PLZBmiss == 1
replace B_PLZ4 = 7031 if bürgerort == "Laax" & PLZBmiss == 1
replace B_PLZ4 = 8853 if bürgerort == "Lachen" & PLZBmiss == 1
replace B_PLZ4 = 1287 if bürgerort == "Laconnex" & PLZBmiss == 1
replace B_PLZ4 = 3931 if bürgerort == "Lalden" & PLZBmiss == 1
replace B_PLZ4 = 6814 if bürgerort == "Lamone" & PLZBmiss == 1
replace B_PLZ4 = 4432 if bürgerort == "Lampenberg" & PLZBmiss == 1
replace B_PLZ4 = 1200 if bürgerort == "Lancy" & PLZBmiss == 1
replace B_PLZ4 = 2525 if bürgerort == "Landeron-Combes" & PLZBmiss == 1
replace B_PLZ4 = 3434 if bürgerort == "Landiswil" & PLZBmiss == 1
replace B_PLZ4 = 7302 if bürgerort == "Landquart" & PLZBmiss == 1
replace B_PLZ4 = 4438 if bürgerort == "Langenbruck" & PLZBmiss == 1
replace B_PLZ4 = 4513 if bürgerort == "Langendorf" & PLZBmiss == 1
replace B_PLZ4 = 4900 if bürgerort == "Langenthal" & PLZBmiss == 1
replace B_PLZ4 = 8135 if bürgerort == "Langnau am Albis" & PLZBmiss == 1
replace B_PLZ4 = 6260 if bürgerort == "Langnau bei Reiden" & PLZBmiss == 1
replace B_PLZ4 = 3550 if bürgerort == "Langnau im Emmental" & PLZBmiss == 1
replace B_PLZ4 = 8585 if bürgerort == "Langrickenbach" & PLZBmiss == 1
replace B_PLZ4 = 7083 if bürgerort == "Lantsch/Lenz" & PLZBmiss == 1
replace B_PLZ4 = 3782 if bürgerort == "Lauenen" & PLZBmiss == 1
replace B_PLZ4 = 6424 if bürgerort == "Lauerz" & PLZBmiss == 1
replace B_PLZ4 = 4242 if bürgerort == "Laufen" & PLZBmiss == 1
replace B_PLZ4 = 8212 if bürgerort == "Laufen-Uhwiesen" & PLZBmiss == 1
replace B_PLZ4 = 5085 if bürgerort == "Laufenburg" & PLZBmiss == 1
replace B_PLZ4 = 3177 if bürgerort == "Laupen" & PLZBmiss == 1
replace B_PLZ4 = 4712 if bürgerort == "Laupersdorf" & PLZBmiss == 1
replace B_PLZ4 = 3543 if bürgerort == "Lauperswil" & PLZBmiss == 1
replace B_PLZ4 = 1000 if bürgerort == "Lausanne" & PLZBmiss == 1
replace B_PLZ4 = 3822 if bürgerort == "Lauterbrunnen" & PLZBmiss == 1
replace B_PLZ4 = 4426 if bürgerort == "Lauwil" & PLZBmiss == 1
replace B_PLZ4 = 6633 if bürgerort == "Lavertezzo" & PLZBmiss == 1
replace B_PLZ4 = 1892 if bürgerort == "Lavey-Morcles" & PLZBmiss == 1
replace B_PLZ4 = 1175 if bürgerort == "Lavigny" & PLZBmiss == 1
replace B_PLZ4 = 6692 if bürgerort == "Lavizzara" & PLZBmiss == 1
replace B_PLZ4 = 2414 if bürgerort == "Le Cerneux-Péquignot" & PLZBmiss == 1
replace B_PLZ4 = 1341 if bürgerort == "Le Chenit" & PLZBmiss == 1
replace B_PLZ4 = 1689 if bürgerort == "Le Châtelard" & PLZBmiss == 1
replace B_PLZ4 = 1611 if bürgerort == "Le Crêt" & PLZBmiss == 1
replace B_PLZ4 = 1699 if bürgerort == "Le Flon" & PLZBmiss == 1
replace B_PLZ4 = 1211 if bürgerort == "Le Grand-Saconnex" & PLZBmiss == 1
replace B_PLZ4 = 1345 if bürgerort == "Le Lieu" & PLZBmiss == 1
replace B_PLZ4 = 2400 if bürgerort == "Le Locle" & PLZBmiss == 1
replace B_PLZ4 = 2340 if bürgerort == "Le Noirmont" & PLZBmiss == 1
replace B_PLZ4 = 5325 if bürgerort == "Leibstadt" & PLZBmiss == 1
replace B_PLZ4 = 4935 if bürgerort == "Leimiswil" & PLZBmiss == 1
replace B_PLZ4 = 8574 if bürgerort == "Lengwil" & PLZBmiss == 1
replace B_PLZ4 = 3775 if bürgerort == "Lenk" & PLZBmiss == 1
replace B_PLZ4 = 1978 if bürgerort == "Lens" & PLZBmiss == 1
replace B_PLZ4 = 5600 if bürgerort == "Lenzburg" & PLZBmiss == 1
replace B_PLZ4 = 2127 if bürgerort == "Les Bayards" & PLZBmiss == 1
replace B_PLZ4 = 2326 if bürgerort == "Les Bois" & PLZBmiss == 1
replace B_PLZ4 = 2416 if bürgerort == "Les Brenets" & PLZBmiss == 1
replace B_PLZ4 = 2206 if bürgerort == "Les Geneveys-sur-Coffrane" & PLZBmiss == 1
replace B_PLZ4 = 1483 if bürgerort == "Les Montets" & PLZBmiss == 1
replace B_PLZ4 = 2325 if bürgerort == "Les Planchettes" & PLZBmiss == 1
replace B_PLZ4 = 2316 if bürgerort == "Les Ponts-de-Martel" & PLZBmiss == 1
replace B_PLZ4 = 2126 if bürgerort == "Les Verrières" & PLZBmiss == 1
replace B_PLZ4 = 5316 if bürgerort == "Leuggern" & PLZBmiss == 1
replace B_PLZ4 = 3952 if bürgerort == "Leuk" & PLZBmiss == 1
replace B_PLZ4 = 3954 if bürgerort == "Leukerbad" & PLZBmiss == 1
replace B_PLZ4 = 1911 if bürgerort == "Leytron" & PLZBmiss == 1
replace B_PLZ4 = 1945 if bürgerort == "Liddes" & PLZBmiss == 1
replace B_PLZ4 = 3213 if bürgerort == "Liebistorf" & PLZBmiss == 1
replace B_PLZ4 = 4254 if bürgerort == "Liesberg" & PLZBmiss == 1
replace B_PLZ4 = 4410 if bürgerort == "Liestal" & PLZBmiss == 1
replace B_PLZ4 = 2514 if bürgerort == "Ligerz" & PLZBmiss == 1
replace B_PLZ4 = 1357 if bürgerort == "Lignerolle" & PLZBmiss == 1
replace B_PLZ4 = 2523 if bürgerort == "Lignières" & PLZBmiss == 1
replace B_PLZ4 = 6853 if bürgerort == "Ligornetto" & PLZBmiss == 1
replace B_PLZ4 = 8310 if bürgerort == "Lindau" & PLZBmiss == 1
replace B_PLZ4 = 3673 if bürgerort == "Linden" & PLZBmiss == 1
replace B_PLZ4 = 6682 if bürgerort == "Linescio" & PLZBmiss == 1
replace B_PLZ4 = 6014 if bürgerort == "Littau" & PLZBmiss == 1
replace B_PLZ4 = 6600 if bürgerort == "Locarno" & PLZBmiss == 1
replace B_PLZ4 = 6527 if bürgerort == "Lodrino" & PLZBmiss == 1
replace B_PLZ4 = 6616 if bürgerort == "Losone" & PLZBmiss == 1
replace B_PLZ4 = 6558 if bürgerort == "Lostallo" & PLZBmiss == 1
replace B_PLZ4 = 4654 if bürgerort == "Lostorf" & PLZBmiss == 1
replace B_PLZ4 = 1682 if bürgerort == "Lovatens" & PLZBmiss == 1
replace B_PLZ4 = 2732 if bürgerort == "Loveresse" & PLZBmiss == 1
replace B_PLZ4 = 1683 if bürgerort == "Lucens" & PLZBmiss == 1
replace B_PLZ4 = 8426 if bürgerort == "Lufingen" & PLZBmiss == 1
replace B_PLZ4 = 6962 if bürgerort == "Lugano" & PLZBmiss == 1
replace B_PLZ4 = 7147 if bürgerort == "Lumnezia" & PLZBmiss == 1
replace B_PLZ4 = 6078 if bürgerort == "Lungern" & PLZBmiss == 1
replace B_PLZ4 = 5242 if bürgerort == "Lupfig" & PLZBmiss == 1
replace B_PLZ4 = 4419 if bürgerort == "Lupsingen" & PLZBmiss == 1
replace B_PLZ4 = 3215 if bürgerort == "Lurtigen" & PLZBmiss == 1
replace B_PLZ4 = 4542 if bürgerort == "Luterbach" & PLZBmiss == 1
replace B_PLZ4 = 6154 if bürgerort == "Luthern" & PLZBmiss == 1
replace B_PLZ4 = 1095 if bürgerort == "Lutry" & PLZBmiss == 1
replace B_PLZ4 = 7141 if bürgerort == "Luven" & PLZBmiss == 1
replace B_PLZ4 = 7245 if bürgerort == "Luzein" & PLZBmiss == 1
replace B_PLZ4 = 6000 if bürgerort == "Luzern" & PLZBmiss == 1
replace B_PLZ4 = 3250 if bürgerort == "Lyss" & PLZBmiss == 1
replace B_PLZ4 = 3636 if bürgerort == "Längenbühl" & PLZBmiss == 1
replace B_PLZ4 = 1773 if bürgerort == "Léchelles" & PLZBmiss == 1
replace B_PLZ4 = 8224 if bürgerort == "Löhningen" & PLZBmiss == 1
replace B_PLZ4 = 2576 if bürgerort == "Lüscherz" & PLZBmiss == 1
replace B_PLZ4 = 4574 if bürgerort == "Lüsslingen-Nennigkofen" & PLZBmiss == 1
replace B_PLZ4 = 9604 if bürgerort == "Lütisburg" & PLZBmiss == 1
replace B_PLZ4 = 3435 if bürgerort == "Lützelflüh" & PLZBmiss == 1
replace B_PLZ4 = 4935 if bürgerort == "Madiswil" & PLZBmiss == 1
replace B_PLZ4 = 4312 if bürgerort == "Magden" & PLZBmiss == 1
replace B_PLZ4 = 6983 if bürgerort == "Magliaso" & PLZBmiss == 1
replace B_PLZ4 = 7304 if bürgerort == "Maienfeld" & PLZBmiss == 1
replace B_PLZ4 = 4464 if bürgerort == "Maisprach" & PLZBmiss == 1
replace B_PLZ4 = 7208 if bürgerort == "Malans" & PLZBmiss == 1
replace B_PLZ4 = 7074 if bürgerort == "Malix" & PLZBmiss == 1
replace B_PLZ4 = 2735 if bürgerort == "Malleray" & PLZBmiss == 1
replace B_PLZ4 = 6102 if bürgerort == "Malters" & PLZBmiss == 1
replace B_PLZ4 = 6713 if bürgerort == "Malvaglia" & PLZBmiss == 1
replace B_PLZ4 = 8265 if bürgerort == "Mammern" & PLZBmiss == 1
replace B_PLZ4 = 5318 if bürgerort == "Mandach" & PLZBmiss == 1
replace B_PLZ4 = 1775 if bürgerort == "Mannens-Grandsivaz" & PLZBmiss == 1
replace B_PLZ4 = 1613 if bürgerort == "Maracon" & PLZBmiss == 1
replace B_PLZ4 = 1261 if bürgerort == "Marchissy" & PLZBmiss == 1
replace B_PLZ4 = 6817 if bürgerort == "Maroggia" & PLZBmiss == 1
replace B_PLZ4 = 1633 if bürgerort == "Marsens" & PLZBmiss == 1
replace B_PLZ4 = 8464 if bürgerort == "Marthalen" & PLZBmiss == 1
replace B_PLZ4 = 1920 if bürgerort == "Martigny" & PLZBmiss == 1
replace B_PLZ4 = 1921 if bürgerort == "Martigny-Combe" & PLZBmiss == 1
replace B_PLZ4 = 8933 if bürgerort == "Maschwanden" & PLZBmiss == 1
replace B_PLZ4 = 7425 if bürgerort == "Masein" & PLZBmiss == 1
replace B_PLZ4 = 6900 if bürgerort == "Massagno" & PLZBmiss == 1
replace B_PLZ4 = 1871 if bürgerort == "Massongex" & PLZBmiss == 1
replace B_PLZ4 = 1692 if bürgerort == "Massonnens" & PLZBmiss == 1
replace B_PLZ4 = 7433 if bürgerort == "Mathon" & PLZBmiss == 1
replace B_PLZ4 = 1753 if bürgerort == "Matran" & PLZBmiss == 1
replace B_PLZ4 = 3800 if bürgerort == "Matten bei Interlaken" & PLZBmiss == 1
replace B_PLZ4 = 3322 if bürgerort == "Mattstetten" & PLZBmiss == 1
replace B_PLZ4 = 4713 if bürgerort == "Matzendorf" & PLZBmiss == 1
replace B_PLZ4 = 9548 if bürgerort == "Matzingen" & PLZBmiss == 1
replace B_PLZ4 = 1453 if bürgerort == "Mauborget" & PLZBmiss == 1
replace B_PLZ4 = 8122 if bürgerort == "Maur" & PLZBmiss == 1
replace B_PLZ4 = 1148 if bürgerort == "Mauraz" & PLZBmiss == 1
replace B_PLZ4 = 6809 if bürgerort == "Medeglia" & PLZBmiss == 1
replace B_PLZ4 = 7185 if bürgerort == "Medel (Lucmagn)" & PLZBmiss == 1
replace B_PLZ4 = 6045 if bürgerort == "Meggen" & PLZBmiss == 1
replace B_PLZ4 = 8706 if bürgerort == "Meilen" & PLZBmiss == 1
replace B_PLZ4 = 1252 if bürgerort == "Meinier" & PLZBmiss == 1
replace B_PLZ4 = 3860 if bürgerort == "Meiringen" & PLZBmiss == 1
replace B_PLZ4 = 4917 if bürgerort == "Melchnau" & PLZBmiss == 1
replace B_PLZ4 = 6815 if bürgerort == "Melide" & PLZBmiss == 1
replace B_PLZ4 = 5507 if bürgerort == "Mellingen" & PLZBmiss == 1
replace B_PLZ4 = 8879 if bürgerort == "Mels" & PLZBmiss == 1
replace B_PLZ4 = 6825 if bürgerort == "Mendrisio" & PLZBmiss == 1
replace B_PLZ4 = 5737 if bürgerort == "Menziken" & PLZBmiss == 1
replace B_PLZ4 = 6313 if bürgerort == "Menzingen" & PLZBmiss == 1
replace B_PLZ4 = 6122 if bürgerort == "Menznau" & PLZBmiss == 1
replace B_PLZ4 = 2827 if bürgerort == "Mervelier" & PLZBmiss == 1
replace B_PLZ4 = 6563 if bürgerort == "Mesocco" & PLZBmiss == 1
replace B_PLZ4 = 3254 if bürgerort == "Messen" & PLZBmiss == 1
replace B_PLZ4 = 5274 if bürgerort == "Mettauertal" & PLZBmiss == 1
replace B_PLZ4 = 8932 if bürgerort == "Mettmenstetten" & PLZBmiss == 1
replace B_PLZ4 = 1217 if bürgerort == "Meyrin" & PLZBmiss == 1
replace B_PLZ4 = 6805 if bürgerort == "Mezzovico-Vira" & PLZBmiss == 1
replace B_PLZ4 = 1295 if bürgerort == "Mies" & PLZBmiss == 1
replace B_PLZ4 = 6986 if bürgerort == "Miglieglia" & PLZBmiss == 1
replace B_PLZ4 = 2013 if bürgerort == "Milvignes" & PLZBmiss == 1
replace B_PLZ4 = 3532 if bürgerort == "Mirchel" & PLZBmiss == 1
replace B_PLZ4 = 1721 if bürgerort == "Misery-Courtion" & PLZBmiss == 1
replace B_PLZ4 = 1565 if bürgerort == "Missy" & PLZBmiss == 1
replace B_PLZ4 = 3972 if bürgerort == "Miège" & PLZBmiss == 1
replace B_PLZ4 = 9123 if bürgerort == "Mogelsberg" & PLZBmiss == 1
replace B_PLZ4 = 1148 if bürgerort == "Moiry" & PLZBmiss == 1
replace B_PLZ4 = 8753 if bürgerort == "Mollis" & PLZBmiss == 1
replace B_PLZ4 = 1415 if bürgerort == "Molondin" & PLZBmiss == 1
replace B_PLZ4 = 1125 if bürgerort == "Monnaz" & PLZBmiss == 1
replace B_PLZ4 = 1973 if bürgerort == "Mont-Noble" & PLZBmiss == 1
replace B_PLZ4 = 1786 if bürgerort == "Mont-Vully" & PLZBmiss == 1
replace B_PLZ4 = 1185 if bürgerort == "Mont-sur-Rolle" & PLZBmiss == 1
replace B_PLZ4 = 1440 if bürgerort == "Montagny-près-Yverdon" & PLZBmiss == 1
replace B_PLZ4 = 2027 if bürgerort == "Montalchez" & PLZBmiss == 1
replace B_PLZ4 = 3963 if bürgerort == "Montana" & PLZBmiss == 1
replace B_PLZ4 = 1515 if bürgerort == "Montanaire" & PLZBmiss == 1
replace B_PLZ4 = 1669 if bürgerort == "Montbovon" & PLZBmiss == 1
replace B_PLZ4 = 1475 if bürgerort == "Montbrelloz" & PLZBmiss == 1
replace B_PLZ4 = 1354 if bürgerort == "Montcherand" & PLZBmiss == 1
replace B_PLZ4 = 6513 if bürgerort == "Monte Carasso" & PLZBmiss == 1
replace B_PLZ4 = 6802 if bürgerort == "Monteceneri" & PLZBmiss == 1
replace B_PLZ4 = 6996 if bürgerort == "Monteggio" & PLZBmiss == 1
replace B_PLZ4 = 1674 if bürgerort == "Montet (Glâne)" & PLZBmiss == 1
replace B_PLZ4 = 1174 if bürgerort == "Montherod" & PLZBmiss == 1
replace B_PLZ4 = 1870 if bürgerort == "Monthey" & PLZBmiss == 1
replace B_PLZ4 = 2924 if bürgerort == "Montignez" & PLZBmiss == 1
replace B_PLZ4 = 1043 if bürgerort == "Montilliez" & PLZBmiss == 1
replace B_PLZ4 = 1824 if bürgerort == "Montreux" & PLZBmiss == 1
replace B_PLZ4 = 2828 if bürgerort == "Montsevelier" & PLZBmiss == 1
replace B_PLZ4 = 6834 if bürgerort == "Morbio Inferiore" & PLZBmiss == 1
replace B_PLZ4 = 1541 if bürgerort == "Morens (FR)" & PLZBmiss == 1
replace B_PLZ4 = 1110 if bürgerort == "Morges" & PLZBmiss == 1
replace B_PLZ4 = 6433 if bürgerort == "Morschach" & PLZBmiss == 1
replace B_PLZ4 = 9612 if bürgerort == "Mosnang" & PLZBmiss == 1
replace B_PLZ4 = 1510 if bürgerort == "Moudon" & PLZBmiss == 1
replace B_PLZ4 = 2740 if bürgerort == "Moutier" & PLZBmiss == 1
replace B_PLZ4 = 5037 if bürgerort == "Muhen" & PLZBmiss == 1
replace B_PLZ4 = 4322 if bürgerort == "Mumpf" & PLZBmiss == 1
replace B_PLZ4 = 3903 if bürgerort == "Mund" & PLZBmiss == 1
replace B_PLZ4 = 7138 if bürgerort == "Mundaun" & PLZBmiss == 1
replace B_PLZ4 = 9313 if bürgerort == "Muolen" & PLZBmiss == 1
replace B_PLZ4 = 4856 if bürgerort == "Murgenthal" & PLZBmiss == 1
replace B_PLZ4 = 3074 if bürgerort == "Muri bei Bern" & PLZBmiss == 1
replace B_PLZ4 = 2338 if bürgerort == "Muriaux" & PLZBmiss == 1
replace B_PLZ4 = 1793 if bürgerort == "Murten" & PLZBmiss == 1
replace B_PLZ4 = 1428 if bürgerort == "Mutrux" & PLZBmiss == 1
replace B_PLZ4 = 7431 if bürgerort == "Mutten" & PLZBmiss == 1
replace B_PLZ4 = 4132 if bürgerort == "Muttenz" & PLZBmiss == 1
replace B_PLZ4 = 6933 if bürgerort == "Muzzano" & PLZBmiss == 1
replace B_PLZ4 = 8708 if bürgerort == "Männedorf" & PLZBmiss == 1
replace B_PLZ4 = 8561 if bürgerort == "Märstetten" & PLZBmiss == 1
replace B_PLZ4 = 1533 if bürgerort == "Ménières" & PLZBmiss == 1
replace B_PLZ4 = 4313 if bürgerort == "Möhlin" & PLZBmiss == 1
replace B_PLZ4 = 8617 if bürgerort == "Mönchaltorf" & PLZBmiss == 1
replace B_PLZ4 = 5103 if bürgerort == "Möriken-Wildegg" & PLZBmiss == 1
replace B_PLZ4 = 9402 if bürgerort == "Mörschwil" & PLZBmiss == 1
replace B_PLZ4 = 5642 if bürgerort == "Mühlau" & PLZBmiss == 1
replace B_PLZ4 = 3204 if bürgerort == "Mühleberg" & PLZBmiss == 1
replace B_PLZ4 = 3127 if bürgerort == "Mühlethurnen" & PLZBmiss == 1
replace B_PLZ4 = 8555 if bürgerort == "Müllheim" & PLZBmiss == 1
replace B_PLZ4 = 4719 if bürgerort == "Mümliswil-Ramiswil" & PLZBmiss == 1
replace B_PLZ4 = 3053 if bürgerort == "Münchenbuchsee" & PLZBmiss == 1
replace B_PLZ4 = 4142 if bürgerort == "Münchenstein" & PLZBmiss == 1
replace B_PLZ4 = 3110 if bürgerort == "Münsingen" & PLZBmiss == 1
replace B_PLZ4 = 3985 if bürgerort == "Münster-Geschinen" & PLZBmiss == 1
replace B_PLZ4 = 3225 if bürgerort == "Müntschemier" & PLZBmiss == 1
replace B_PLZ4 = 7537 if bürgerort == "Müstair" & PLZBmiss == 1
replace B_PLZ4 = 6289 if bürgerort == "Müswangen" & PLZBmiss == 1
replace B_PLZ4 = 3914 if bürgerort == "Naters" & PLZBmiss == 1
replace B_PLZ4 = 1973 if bürgerort == "Nax" & PLZBmiss == 1
replace B_PLZ4 = 6244 if bürgerort == "Nebikon" & PLZBmiss == 1
replace B_PLZ4 = 9125 if bürgerort == "Neckertal" & PLZBmiss == 1
replace B_PLZ4 = 8412 if bürgerort == "Neftenbach" & PLZBmiss == 1
replace B_PLZ4 = 1669 if bürgerort == "Neirivue" & PLZBmiss == 1
replace B_PLZ4 = 1996 if bürgerort == "Nendaz" & PLZBmiss == 1
replace B_PLZ4 = 9650 if bürgerort == "Nesslau" & PLZBmiss == 1
replace B_PLZ4 = 8759 if bürgerort == "Netstal" & PLZBmiss == 1
replace B_PLZ4 = 2008 if bürgerort == "Neuchâtel" & PLZBmiss == 1
replace B_PLZ4 = 4623 if bürgerort == "Neuendorf" & PLZBmiss == 1
replace B_PLZ4 = 3176 if bürgerort == "Neuenegg" & PLZBmiss == 1
replace B_PLZ4 = 5432 if bürgerort == "Neuenhof" & PLZBmiss == 1
replace B_PLZ4 = 6206 if bürgerort == "Neuenkirch" & PLZBmiss == 1
replace B_PLZ4 = 8212 if bürgerort == "Neuhausen am Rheinfall" & PLZBmiss == 1
replace B_PLZ4 = 8525 if bürgerort == "Neunforn" & PLZBmiss == 1
replace B_PLZ4 = 4704 if bürgerort == "Niederbipp" & PLZBmiss == 1
replace B_PLZ4 = 4626 if bürgerort == "Niederbuchsiten" & PLZBmiss == 1
replace B_PLZ4 = 9246 if bürgerort == "Niederbüren" & PLZBmiss == 1
replace B_PLZ4 = 5015 if bürgerort == "Niedererlinsbach" & PLZBmiss == 1
replace B_PLZ4 = 3942 if bürgerort == "Niedergesteln" & PLZBmiss == 1
replace B_PLZ4 = 8172 if bürgerort == "Niederglatt" & PLZBmiss == 1
replace B_PLZ4 = 5013 if bürgerort == "Niedergösgen" & PLZBmiss == 1
replace B_PLZ4 = 8155 if bürgerort == "Niederhasli" & PLZBmiss == 1
replace B_PLZ4 = 9526 if bürgerort == "Niederhelfenschwil" & PLZBmiss == 1
replace B_PLZ4 = 3087 if bürgerort == "Niedermuhlern" & PLZBmiss == 1
replace B_PLZ4 = 3283 if bürgerort == "Niederried bei Kallnach" & PLZBmiss == 1
replace B_PLZ4 = 5443 if bürgerort == "Niederrohrdorf" & PLZBmiss == 1
replace B_PLZ4 = 3424 if bürgerort == "Niederösch" & PLZBmiss == 1
replace B_PLZ4 = 3116 if bürgerort == "Noflen" & PLZBmiss == 1
replace B_PLZ4 = 2149 if bürgerort == "Noiraigue" & PLZBmiss == 1
replace B_PLZ4 = 1757 if bürgerort == "Noréaz" & PLZBmiss == 1
replace B_PLZ4 = 6207 if bürgerort == "Nottwil" & PLZBmiss == 1
replace B_PLZ4 = 6986 if bürgerort == "Novaggio" & PLZBmiss == 1
replace B_PLZ4 = 6883 if bürgerort == "Novazzano" & PLZBmiss == 1
replace B_PLZ4 = 1845 if bürgerort == "Noville" & PLZBmiss == 1
replace B_PLZ4 = 4421 if bürgerort == "Nuglar-St. Pantaleon" & PLZBmiss == 1
replace B_PLZ4 = 4208 if bürgerort == "Nunningen" & PLZBmiss == 1
replace B_PLZ4 = 4453 if bürgerort == "Nusshof" & PLZBmiss == 1
replace B_PLZ4 = 1260 if bürgerort == "Nyon" & PLZBmiss == 1
replace B_PLZ4 = 8752 if bürgerort == "Näfels" & PLZBmiss == 1
replace B_PLZ4 = 3096 if bürgerort == "Oberbalm" & PLZBmiss == 1
replace B_PLZ4 = 4538 if bürgerort == "Oberbipp" & PLZBmiss == 1
replace B_PLZ4 = 3414 if bürgerort == "Oberburg" & PLZBmiss == 1
replace B_PLZ4 = 9245 if bürgerort == "Oberbüren" & PLZBmiss == 1
replace B_PLZ4 = 3672 if bürgerort == "Oberdiessbach" & PLZBmiss == 1
replace B_PLZ4 = 9413 if bürgerort == "Oberegg" & PLZBmiss == 1
replace B_PLZ4 = 5420 if bürgerort == "Oberehrendingen" & PLZBmiss == 1
replace B_PLZ4 = 8425 if bürgerort == "Oberembrach" & PLZBmiss == 1
replace B_PLZ4 = 5036 if bürgerort == "Oberentfelden" & PLZBmiss == 1
replace B_PLZ4 = 5108 if bürgerort == "Oberflachs" & PLZBmiss == 1
replace B_PLZ4 = 3988 if bürgerort == "Obergesteln" & PLZBmiss == 1
replace B_PLZ4 = 8154 if bürgerort == "Oberglatt" & PLZBmiss == 1
replace B_PLZ4 = 3999 if bürgerort == "Obergoms" & PLZBmiss == 1
replace B_PLZ4 = 9126 if bürgerort == "Oberhelfenschwil" & PLZBmiss == 1
replace B_PLZ4 = 3653 if bürgerort == "Oberhofen am Thunersee" & PLZBmiss == 1
replace B_PLZ4 = 8843 if bürgerort == "Oberiberg" & PLZBmiss == 1
replace B_PLZ4 = 6208 if bürgerort == "Oberkirch" & PLZBmiss == 1
replace B_PLZ4 = 5727 if bürgerort == "Oberkulm" & PLZBmiss == 1
replace B_PLZ4 = 3616 if bürgerort == "Oberlangenegg" & PLZBmiss == 1
replace B_PLZ4 = 4324 if bürgerort == "Obermumpf" & PLZBmiss == 1
replace B_PLZ4 = 8942 if bürgerort == "Oberrieden" & PLZBmiss == 1
replace B_PLZ4 = 5452 if bürgerort == "Oberrohrdorf" & PLZBmiss == 1
replace B_PLZ4 = 7134 if bürgerort == "Obersaxen" & PLZBmiss == 1
replace B_PLZ4 = 7138 if bürgerort == "Obersaxen Mundaun" & PLZBmiss == 1
replace B_PLZ4 = 1716 if bürgerort == "Oberschrot" & PLZBmiss == 1
replace B_PLZ4 = 5415 if bürgerort == "Obersiggenthal" & PLZBmiss == 1
replace B_PLZ4 = 8477 if bürgerort == "Oberstammheim" & PLZBmiss == 1
replace B_PLZ4 = 3531 if bürgerort == "Oberthal" & PLZBmiss == 1
replace B_PLZ4 = 8868 if bürgerort == "Oberurnen" & PLZBmiss == 1
replace B_PLZ4 = 9242 if bürgerort == "Oberuzwil" & PLZBmiss == 1
replace B_PLZ4 = 3114 if bürgerort == "Oberwichtrach" & PLZBmiss == 1
replace B_PLZ4 = 3298 if bürgerort == "Oberwil bei Büren" & PLZBmiss == 1
replace B_PLZ4 = 3765 if bürgerort == "Oberwil im Simmental" & PLZBmiss == 1
replace B_PLZ4 = 6315 if bürgerort == "Oberägeri" & PLZBmiss == 1
replace B_PLZ4 = 3363 if bürgerort == "Oberönz" & PLZBmiss == 1
replace B_PLZ4 = 3424 if bürgerort == "Oberösch" & PLZBmiss == 1
replace B_PLZ4 = 8912 if bürgerort == "Obfelden" & PLZBmiss == 1
replace B_PLZ4 = 8758 if bürgerort == "Obstalden" & PLZBmiss == 1
replace B_PLZ4 = 4566 if bürgerort == "Oekingen" & PLZBmiss == 1
replace B_PLZ4 = 4702 if bürgerort == "Oensingen" & PLZBmiss == 1
replace B_PLZ4 = 4943 if bürgerort == "Oeschenbach" & PLZBmiss == 1
replace B_PLZ4 = 8618 if bürgerort == "Oetwil am See" & PLZBmiss == 1
replace B_PLZ4 = 8955 if bürgerort == "Oetwil an der Limmat" & PLZBmiss == 1
replace B_PLZ4 = 4665 if bürgerort == "Oftringen" & PLZBmiss == 1
replace B_PLZ4 = 1867 if bürgerort == "Ollon" & PLZBmiss == 1
replace B_PLZ4 = 4600 if bürgerort == "Olten" & PLZBmiss == 1
replace B_PLZ4 = 1213 if bürgerort == "Onex" & PLZBmiss == 1
replace B_PLZ4 = 6664 if bürgerort == "Onsernone" & PLZBmiss == 1
replace B_PLZ4 = 8152 if bürgerort == "Opfikon" & PLZBmiss == 1
replace B_PLZ4 = 3629 if bürgerort == "Oppligen" & PLZBmiss == 1
replace B_PLZ4 = 1350 if bürgerort == "Orbe" & PLZBmiss == 1
replace B_PLZ4 = 1430 if bürgerort == "Orges" & PLZBmiss == 1
replace B_PLZ4 = 6945 if bürgerort == "Origlio" & PLZBmiss == 1
replace B_PLZ4 = 4466 if bürgerort == "Ormalingen" & PLZBmiss == 1
replace B_PLZ4 = 1866 if bürgerort == "Ormont-Dessous" & PLZBmiss == 1
replace B_PLZ4 = 1865 if bürgerort == "Ormont-Dessus" & PLZBmiss == 1
replace B_PLZ4 = 1607 if bürgerort == "Oron" & PLZBmiss == 1
replace B_PLZ4 = 1610 if bürgerort == "Oron-la-Ville" & PLZBmiss == 1
replace B_PLZ4 = 1944 if bürgerort == "Orsières" & PLZBmiss == 1
replace B_PLZ4 = 6763 if bürgerort == "Osco" & PLZBmiss == 1
replace B_PLZ4 = 8475 if bürgerort == "Ossingen" & PLZBmiss == 1
replace B_PLZ4 = 3070 if bürgerort == "Ostermundigen" & PLZBmiss == 1
replace B_PLZ4 = 8913 if bürgerort == "Ottenbach" & PLZBmiss == 1
replace B_PLZ4 = 1377 if bürgerort == "Oulens-sous-Echallens" & PLZBmiss == 1
replace B_PLZ4 = 6915 if bürgerort == "Pambio-Noranco" & PLZBmiss == 1
replace B_PLZ4 = 1142 if bürgerort == "Pampigny" & PLZBmiss == 1
replace B_PLZ4 = 6902 if bürgerort == "Paradiso" & PLZBmiss == 1
replace B_PLZ4 = 7433 if bürgerort == "Patzen-Fardün" & PLZBmiss == 1
replace B_PLZ4 = 1530 if bürgerort == "Payerne" & PLZBmiss == 1
replace B_PLZ4 = 1059 if bürgerort == "Peney-le-Jorat" & PLZBmiss == 1
replace B_PLZ4 = 1375 if bürgerort == "Penthéréaz" & PLZBmiss == 1
replace B_PLZ4 = 2034 if bürgerort == "Peseux" & PLZBmiss == 1
replace B_PLZ4 = 2748 if bürgerort == "Petit-Val" & PLZBmiss == 1
replace B_PLZ4 = 4915 if bürgerort == "Pfaffnau" & PLZBmiss == 1
replace B_PLZ4 = 5735 if bürgerort == "Pfeffikon" & PLZBmiss == 1
replace B_PLZ4 = 4148 if bürgerort == "Pfeffingen" & PLZBmiss == 1
replace B_PLZ4 = 8505 if bürgerort == "Pfyn" & PLZBmiss == 1
replace B_PLZ4 = 7315 if bürgerort == "Pfäfers" & PLZBmiss == 1
replace B_PLZ4 = 8330 if bürgerort == "Pfäffikon" & PLZBmiss == 1
replace B_PLZ4 = 1716 if bürgerort == "Plaffeien" & PLZBmiss == 1
replace B_PLZ4 = 1228 if bürgerort == "Plan-les-Ouates" & PLZBmiss == 1
replace B_PLZ4 = 1737 if bürgerort == "Plasselb" & PLZBmiss == 1
replace B_PLZ4 = 2515 if bürgerort == "Plateau de Diesse" & PLZBmiss == 1
replace B_PLZ4 = 6742 if bürgerort == "Pollegio" & PLZBmiss == 1
replace B_PLZ4 = 1318 if bürgerort == "Pompaples" & PLZBmiss == 1
replace B_PLZ4 = 1649 if bürgerort == "Pont-la-Ville" & PLZBmiss == 1
replace B_PLZ4 = 2900 if bürgerort == "Porrentruy" & PLZBmiss == 1
replace B_PLZ4 = 1897 if bürgerort == "Port-Valais" & PLZBmiss == 1
replace B_PLZ4 = 7745 if bürgerort == "Poschiavo" & PLZBmiss == 1
replace B_PLZ4 = 1724 if bürgerort == "Praroman" & PLZBmiss == 1
replace B_PLZ4 = 6772 if bürgerort == "Prato (Leventina)" & PLZBmiss == 1
replace B_PLZ4 = 4133 if bürgerort == "Pratteln" & PLZBmiss == 1
replace B_PLZ4 = 1211 if bürgerort == "Pregny-Chambésy" & PLZBmiss == 1
replace B_PLZ4 = 6523 if bürgerort == "Preonzo" & PLZBmiss == 1
replace B_PLZ4 = 1243 if bürgerort == "Presinge" & PLZBmiss == 1
replace B_PLZ4 = 1746 if bürgerort == "Prez-vers-Noréaz" & PLZBmiss == 1
replace B_PLZ4 = 1677 if bürgerort == "Prez-vers-Siviriez" & PLZBmiss == 1
replace B_PLZ4 = 1008 if bürgerort == "Prilly" & PLZBmiss == 1
replace B_PLZ4 = 1624 if bürgerort == "Progens" & PLZBmiss == 1
replace B_PLZ4 = 1428 if bürgerort == "Provence" & PLZBmiss == 1
replace B_PLZ4 = 7424 if bürgerort == "Präz" & PLZBmiss == 1
replace B_PLZ4 = 1682 if bürgerort == "Prévonloup" & PLZBmiss == 1
replace B_PLZ4 = 1070 if bürgerort == "Puidoux" & PLZBmiss == 1
replace B_PLZ4 = 1009 if bürgerort == "Pully" & PLZBmiss == 1
replace B_PLZ4 = 2603 if bürgerort == "Péry-La Heutte" & PLZBmiss == 1
replace B_PLZ4 = 8885 if bürgerort == "Quarten" & PLZBmiss == 1
replace B_PLZ4 = 6775 if bürgerort == "Quinto" & PLZBmiss == 1
replace B_PLZ4 = 3036 if bürgerort == "Radelfingen" & PLZBmiss == 1
replace B_PLZ4 = 8197 if bürgerort == "Rafz" & PLZBmiss == 1
replace B_PLZ4 = 6026 if bürgerort == "Rain" & PLZBmiss == 1
replace B_PLZ4 = 7557 if bürgerort == "Ramosch" & PLZBmiss == 1
replace B_PLZ4 = 8262 if bürgerort == "Ramsen" & PLZBmiss == 1
replace B_PLZ4 = 6862 if bürgerort == "Rancate" & PLZBmiss == 1
replace B_PLZ4 = 3975 if bürgerort == "Randogne" & PLZBmiss == 1
replace B_PLZ4 = 8558 if bürgerort == "Raperswilen" & PLZBmiss == 1
replace B_PLZ4 = 8640 if bürgerort == "Rapperswil-Jona" & PLZBmiss == 1
replace B_PLZ4 = 3942 if bürgerort == "Raron" & PLZBmiss == 1
replace B_PLZ4 = 9445 if bürgerort == "Rebstein" & PLZBmiss == 1
replace B_PLZ4 = 4565 if bürgerort == "Recherswil" & PLZBmiss == 1
replace B_PLZ4 = 1718 if bürgerort == "Rechthalten" & PLZBmiss == 1
replace B_PLZ4 = 2732 if bürgerort == "Reconvilier" & PLZBmiss == 1
replace B_PLZ4 = 8105 if bürgerort == "Regensdorf" & PLZBmiss == 1
replace B_PLZ4 = 9038 if bürgerort == "Rehetobel" & PLZBmiss == 1
replace B_PLZ4 = 3713 if bürgerort == "Reichenbach im Kandertal" & PLZBmiss == 1
replace B_PLZ4 = 8864 if bürgerort == "Reichenburg" & PLZBmiss == 1
replace B_PLZ4 = 6263 if bürgerort == "Reiden" & PLZBmiss == 1
replace B_PLZ4 = 4418 if bürgerort == "Reigoldswil" & PLZBmiss == 1
replace B_PLZ4 = 4919 if bürgerort == "Reisiswil" & PLZBmiss == 1
replace B_PLZ4 = 5057 if bürgerort == "Reitnau" & PLZBmiss == 1
replace B_PLZ4 = 5453 if bürgerort == "Remetschwil" & PLZBmiss == 1
replace B_PLZ4 = 5236 if bürgerort == "Remigen" & PLZBmiss == 1
replace B_PLZ4 = 3647 if bürgerort == "Reutigen" & PLZBmiss == 1
replace B_PLZ4 = 8462 if bürgerort == "Rheinau" & PLZBmiss == 1
replace B_PLZ4 = 9424 if bürgerort == "Rheineck" & PLZBmiss == 1
replace B_PLZ4 = 4310 if bürgerort == "Rheinfelden" & PLZBmiss == 1
replace B_PLZ4 = 7403 if bürgerort == "Rhäzüns" & PLZBmiss == 1
replace B_PLZ4 = 6263 if bürgerort == "Richenthal" & PLZBmiss == 1
replace B_PLZ4 = 8805 if bürgerort == "Richterswil" & PLZBmiss == 1
replace B_PLZ4 = 9532 if bürgerort == "Rickenbach bei Wil" & PLZBmiss == 1
replace B_PLZ4 = 3216 if bürgerort == "Ried bei Kerzers" & PLZBmiss == 1
replace B_PLZ4 = 8739 if bürgerort == "Rieden" & PLZBmiss == 1
replace B_PLZ4 = 4533 if bürgerort == "Riedholz" & PLZBmiss == 1
replace B_PLZ4 = 4125 if bürgerort == "Riehen" & PLZBmiss == 1
replace B_PLZ4 = 6452 if bürgerort == "Riemenstalden" & PLZBmiss == 1
replace B_PLZ4 = 5323 if bürgerort == "Rietheim" & PLZBmiss == 1
replace B_PLZ4 = 1097 if bürgerort == "Riex" & PLZBmiss == 1
replace B_PLZ4 = 3132 if bürgerort == "Riggisberg" & PLZBmiss == 1
replace B_PLZ4 = 6826 if bürgerort == "Riva San Vitale" & PLZBmiss == 1
replace B_PLZ4 = 1071 if bürgerort == "Rivaz" & PLZBmiss == 1
replace B_PLZ4 = 6703 if bürgerort == "Riviera" & PLZBmiss == 1
replace B_PLZ4 = 6265 if bürgerort == "Roggliswil" & PLZBmiss == 1
replace B_PLZ4 = 4938 if bürgerort == "Rohrbach" & PLZBmiss == 1
replace B_PLZ4 = 4938 if bürgerort == "Rohrbachgraben" & PLZBmiss == 1
replace B_PLZ4 = 1180 if bürgerort == "Rolle" & PLZBmiss == 1
replace B_PLZ4 = 8590 if bürgerort == "Romanshorn" & PLZBmiss == 1
replace B_PLZ4 = 6167 if bürgerort == "Romoos" & PLZBmiss == 1
replace B_PLZ4 = 6622 if bürgerort == "Ronco sopra Ascona" & PLZBmiss == 1
replace B_PLZ4 = 8427 if bürgerort == "Rorbas" & PLZBmiss == 1
replace B_PLZ4 = 9400 if bürgerort == "Rorschach" & PLZBmiss == 1
replace B_PLZ4 = 9404 if bürgerort == "Rorschacherberg" & PLZBmiss == 1
replace B_PLZ4 = 6547 if bürgerort == "Rossa" & PLZBmiss == 1
replace B_PLZ4 = 2842 if bürgerort == "Rossemaison" & PLZBmiss == 1
replace B_PLZ4 = 1658 if bürgerort == "Rossinière" & PLZBmiss == 1
replace B_PLZ4 = 6000 if bürgerort == "Rothenburg" & PLZBmiss == 1
replace B_PLZ4 = 6418 if bürgerort == "Rothenthurm" & PLZBmiss == 1
replace B_PLZ4 = 4852 if bürgerort == "Rothrist" & PLZBmiss == 1
replace B_PLZ4 = 8919 if bürgerort == "Rottenschwil" & PLZBmiss == 1
replace B_PLZ4 = 1659 if bürgerort == "Rougemont" & PLZBmiss == 1
replace B_PLZ4 = 1463 if bürgerort == "Rovray" & PLZBmiss == 1
replace B_PLZ4 = 8964 if bürgerort == "Rudolfstetten-Friedlisberg" & PLZBmiss == 1
replace B_PLZ4 = 1542 if bürgerort == "Rueyres-les-Prés" & PLZBmiss == 1
replace B_PLZ4 = 5102 if bürgerort == "Rupperswil" & PLZBmiss == 1
replace B_PLZ4 = 3251 if bürgerort == "Ruppoldsried" & PLZBmiss == 1
replace B_PLZ4 = 8322 if bürgerort == "Russikon" & PLZBmiss == 1
replace B_PLZ4 = 1281 if bürgerort == "Russin" & PLZBmiss == 1
replace B_PLZ4 = 6017 if bürgerort == "Ruswil" & PLZBmiss == 1
replace B_PLZ4 = 6028 if bürgerort == "Römerswil" & PLZBmiss == 1
replace B_PLZ4 = 3373 if bürgerort == "Röthenbach bei Herzogenbuchsee" & PLZBmiss == 1
replace B_PLZ4 = 3538 if bürgerort == "Röthenbach im Emmental" & PLZBmiss == 1
replace B_PLZ4 = 3439 if bürgerort == "Rüderswil" & PLZBmiss == 1
replace B_PLZ4 = 3422 if bürgerort == "Rüdtligen-Alchenflüh" & PLZBmiss == 1
replace B_PLZ4 = 3089 if bürgerort == "Rüeggisberg" & PLZBmiss == 1
replace B_PLZ4 = 3415 if bürgerort == "Rüegsau" & PLZBmiss == 1
replace B_PLZ4 = 5235 if bürgerort == "Rüfenach" & PLZBmiss == 1
replace B_PLZ4 = 5464 if bürgerort == "Rümikon" & PLZBmiss == 1
replace B_PLZ4 = 8153 if bürgerort == "Rümlang" & PLZBmiss == 1
replace B_PLZ4 = 3154 if bürgerort == "Rüschegg" & PLZBmiss == 1
replace B_PLZ4 = 8803 if bürgerort == "Rüschlikon" & PLZBmiss == 1
replace B_PLZ4 = 3295 if bürgerort == "Rüti bei Büren" & PLZBmiss == 1
replace B_PLZ4 = 3421 if bürgerort == "Rüti bei Lyssach" & PLZBmiss == 1
replace B_PLZ4 = 3099 if bürgerort == "Rüti bei Riggisberg" & PLZBmiss == 1
replace B_PLZ4 = 4933 if bürgerort == "Rütschelen" & PLZBmiss == 1
replace B_PLZ4 = 4522 if bürgerort == "Rüttenen" & PLZBmiss == 1
replace B_PLZ4 = 3783 if bürgerort == "Saanen" & PLZBmiss == 1
replace B_PLZ4 = 7247 if bürgerort == "Saas" & PLZBmiss == 1
replace B_PLZ4 = 3908 if bürgerort == "Saas Balen" & PLZBmiss == 1
replace B_PLZ4 = 3906 if bürgerort == "Saas Fee" & PLZBmiss == 1
replace B_PLZ4 = 3905 if bürgerort == "Saas-Almagell" & PLZBmiss == 1
replace B_PLZ4 = 6073 if bürgerort == "Sachseln" & PLZBmiss == 1
replace B_PLZ4 = 5745 if bürgerort == "Safenwil" & PLZBmiss == 1
replace B_PLZ4 = 7109 if bürgerort == "Safien" & PLZBmiss == 1
replace B_PLZ4 = 7104 if bürgerort == "Safiental" & PLZBmiss == 1
replace B_PLZ4 = 7152 if bürgerort == "Sagogn" & PLZBmiss == 1
replace B_PLZ4 = 2732 if bürgerort == "Saicourt" & PLZBmiss == 1
replace B_PLZ4 = 2354 if bürgerort == "Saignelégier" & PLZBmiss == 1
replace B_PLZ4 = 2364 if bürgerort == "Saint-Brais" & PLZBmiss == 1
replace B_PLZ4 = 1188 if bürgerort == "Saint-George" & PLZBmiss == 1
replace B_PLZ4 = 1898 if bürgerort == "Saint-Gingolph" & PLZBmiss == 1
replace B_PLZ4 = 2610 if bürgerort == "Saint-Imier" & PLZBmiss == 1
replace B_PLZ4 = 3961 if bürgerort == "Saint-Jean" & PLZBmiss == 1
replace B_PLZ4 = 3961 if bürgerort == "Saint-Luc" & PLZBmiss == 1
replace B_PLZ4 = 1806 if bürgerort == "Saint-Légier-La Chiésaz" & PLZBmiss == 1
replace B_PLZ4 = 1958 if bürgerort == "Saint-Léonard" & PLZBmiss == 1
replace B_PLZ4 = 1890 if bürgerort == "Saint-Maurice" & PLZBmiss == 1
replace B_PLZ4 = 1187 if bürgerort == "Saint-Oyens" & PLZBmiss == 1
replace B_PLZ4 = 1113 if bürgerort == "Saint-Saphorin-sur-Morges" & PLZBmiss == 1
replace B_PLZ4 = 2882 if bürgerort == "Saint-Ursanne" & PLZBmiss == 1
replace B_PLZ4 = 1450 if bürgerort == "Sainte-Croix" & PLZBmiss == 1
replace B_PLZ4 = 6954 if bürgerort == "Sala Capriasca" & PLZBmiss == 1
replace B_PLZ4 = 8269 if bürgerort == "Salenstein" & PLZBmiss == 1
replace B_PLZ4 = 8599 if bürgerort == "Salmsach" & PLZBmiss == 1
replace B_PLZ4 = 1922 if bürgerort == "Salvan" & PLZBmiss == 1
replace B_PLZ4 = 7503 if bürgerort == "Samedan" & PLZBmiss == 1
replace B_PLZ4 = 7563 if bürgerort == "Samnaun" & PLZBmiss == 1
replace B_PLZ4 = 6575 if bürgerort == "San Nazzaro" & PLZBmiss == 1
replace B_PLZ4 = 6534 if bürgerort == "San Vittore" & PLZBmiss == 1
replace B_PLZ4 = 6577 if bürgerort == "Sant'Abbondio" & PLZBmiss == 1
replace B_PLZ4 = 6592 if bürgerort == "Sant'Antonino" & PLZBmiss == 1
replace B_PLZ4 = 6583 if bürgerort == "Sant'Antonio" & PLZBmiss == 1
replace B_PLZ4 = 6541 if bürgerort == "Santa Maria in Calanca" & PLZBmiss == 1
replace B_PLZ4 = 7320 if bürgerort == "Sargans" & PLZBmiss == 1
replace B_PLZ4 = 5614 if bürgerort == "Sarmenstorf" & PLZBmiss == 1
replace B_PLZ4 = 6060 if bürgerort == "Sarnen" & PLZBmiss == 1
replace B_PLZ4 = 1683 if bürgerort == "Sarzens" & PLZBmiss == 1
replace B_PLZ4 = 1242 if bürgerort == "Satigny" & PLZBmiss == 1
replace B_PLZ4 = 6417 if bürgerort == "Sattel" & PLZBmiss == 1
replace B_PLZ4 = 2537 if bürgerort == "Sauge" & PLZBmiss == 1
replace B_PLZ4 = 2873 if bürgerort == "Saulcy" & PLZBmiss == 1
replace B_PLZ4 = 1073 if bürgerort == "Savigny" & PLZBmiss == 1
replace B_PLZ4 = 1965 if bürgerort == "Savièse" & PLZBmiss == 1
replace B_PLZ4 = 7460 if bürgerort == "Savognin" & PLZBmiss == 1
replace B_PLZ4 = 1907 if bürgerort == "Saxon" & PLZBmiss == 1
replace B_PLZ4 = 8208 if bürgerort == "Schaffhausen" & PLZBmiss == 1
replace B_PLZ4 = 6197 if bürgerort == "Schangnau" & PLZBmiss == 1
replace B_PLZ4 = 7412 if bürgerort == "Scharans" & PLZBmiss == 1
replace B_PLZ4 = 6467 if bürgerort == "Schattdorf" & PLZBmiss == 1
replace B_PLZ4 = 2827 if bürgerort == "Schelten" & PLZBmiss == 1
replace B_PLZ4 = 6214 if bürgerort == "Schenkon" & PLZBmiss == 1
replace B_PLZ4 = 3305 if bürgerort == "Scheunen" & PLZBmiss == 1
replace B_PLZ4 = 2556 if bürgerort == "Scheuren" & PLZBmiss == 1
replace B_PLZ4 = 7228 if bürgerort == "Schiers" & PLZBmiss == 1
replace B_PLZ4 = 5107 if bürgerort == "Schinznach" & PLZBmiss == 1
replace B_PLZ4 = 5107 if bürgerort == "Schinznach Dorf" & PLZBmiss == 1
replace B_PLZ4 = 7168 if bürgerort == "Schlans" & PLZBmiss == 1
replace B_PLZ4 = 8418 if bürgerort == "Schlatt" & PLZBmiss == 1
replace B_PLZ4 = 9050 if bürgerort == "Schlatt-Haslen" & PLZBmiss == 1
replace B_PLZ4 = 8226 if bürgerort == "Schleitheim" & PLZBmiss == 1
replace B_PLZ4 = 7151 if bürgerort == "Schleuis" & PLZBmiss == 1
replace B_PLZ4 = 6231 if bürgerort == "Schlierbach" & PLZBmiss == 1
replace B_PLZ4 = 8010 if bürgerort == "Schlieren" & PLZBmiss == 1
replace B_PLZ4 = 5044 if bürgerort == "Schlossrued" & PLZBmiss == 1
replace B_PLZ4 = 3082 if bürgerort == "Schlosswil" & PLZBmiss == 1
replace B_PLZ4 = 8716 if bürgerort == "Schmerikon" & PLZBmiss == 1
replace B_PLZ4 = 5046 if bürgerort == "Schmiedrued" & PLZBmiss == 1
replace B_PLZ4 = 5425 if bürgerort == "Schneisingen" & PLZBmiss == 1
replace B_PLZ4 = 3253 if bürgerort == "Schnottwil" & PLZBmiss == 1
replace B_PLZ4 = 6288 if bürgerort == "Schongau" & PLZBmiss == 1
replace B_PLZ4 = 3855 if bürgerort == "Schwanden bei Brienz" & PLZBmiss == 1
replace B_PLZ4 = 6215 if bürgerort == "Schwarzenbach" & PLZBmiss == 1
replace B_PLZ4 = 6103 if bürgerort == "Schwarzenberg" & PLZBmiss == 1
replace B_PLZ4 = 3150 if bürgerort == "Schwarzenburg" & PLZBmiss == 1
replace B_PLZ4 = 9103 if bürgerort == "Schwellbrunn" & PLZBmiss == 1
replace B_PLZ4 = 8603 if bürgerort == "Schwerzenbach" & PLZBmiss == 1
replace B_PLZ4 = 6423 if bürgerort == "Schwyz" & PLZBmiss == 1
replace B_PLZ4 = 8762 if bürgerort == "Schwändi" & PLZBmiss == 1
replace B_PLZ4 = 8718 if bürgerort == "Schänis" & PLZBmiss == 1
replace B_PLZ4 = 5040 if bürgerort == "Schöftland" & PLZBmiss == 1
replace B_PLZ4 = 8577 if bürgerort == "Schönholzerswilen" & PLZBmiss == 1
replace B_PLZ4 = 6143 if bürgerort == "Schötz" & PLZBmiss == 1
replace B_PLZ4 = 8863 if bürgerort == "Schübelbach" & PLZBmiss == 1
replace B_PLZ4 = 3054 if bürgerort == "Schüpfen" & PLZBmiss == 1
replace B_PLZ4 = 6170 if bürgerort == "Schüpfheim" & PLZBmiss == 1
replace B_PLZ4 = 7545 if bürgerort == "Scuol/Schuls" & PLZBmiss == 1
replace B_PLZ4 = 3475 if bürgerort == "Seeberg" & PLZBmiss == 1
replace B_PLZ4 = 8607 if bürgerort == "Seegräben" & PLZBmiss == 1
replace B_PLZ4 = 2747 if bürgerort == "Seehof" & PLZBmiss == 1
replace B_PLZ4 = 6377 if bürgerort == "Seelisberg" & PLZBmiss == 1
replace B_PLZ4 = 5707 if bürgerort == "Seengen" & PLZBmiss == 1
replace B_PLZ4 = 4206 if bürgerort == "Seewen" & PLZBmiss == 1
replace B_PLZ4 = 3662 if bürgerort == "Seftigen" & PLZBmiss == 1
replace B_PLZ4 = 1525 if bürgerort == "Seigneux" & PLZBmiss == 1
replace B_PLZ4 = 2545 if bürgerort == "Selzach" & PLZBmiss == 1
replace B_PLZ4 = 1933 if bürgerort == "Sembrancher" & PLZBmiss == 1
replace B_PLZ4 = 6204 if bürgerort == "Sempach" & PLZBmiss == 1
replace B_PLZ4 = 1623 if bürgerort == "Semsales" & PLZBmiss == 1
replace B_PLZ4 = 1304 if bürgerort == "Senarclens" & PLZBmiss == 1
replace B_PLZ4 = 9468 if bürgerort == "Sennwald" & PLZBmiss == 1
replace B_PLZ4 = 1724 if bürgerort == "Senèdes" & PLZBmiss == 1
replace B_PLZ4 = 5703 if bürgerort == "Seon" & PLZBmiss == 1
replace B_PLZ4 = 1355 if bürgerort == "Sergey" & PLZBmiss == 1
replace B_PLZ4 = 6713 if bürgerort == "Serravalle" & PLZBmiss == 1
replace B_PLZ4 = 1080 if bürgerort == "Servion" & PLZBmiss == 1
replace B_PLZ4 = 6997 if bürgerort == "Sessa" & PLZBmiss == 1
replace B_PLZ4 = 8472 if bürgerort == "Seuzach" & PLZBmiss == 1
replace B_PLZ4 = 9475 if bürgerort == "Sevelen" & PLZBmiss == 1
replace B_PLZ4 = 7157 if bürgerort == "Siat" & PLZBmiss == 1
replace B_PLZ4 = 8225 if bürgerort == "Siblingen" & PLZBmiss == 1
replace B_PLZ4 = 3976 if bürgerort == "Sierre" & PLZBmiss == 1
replace B_PLZ4 = 3535 if bürgerort == "Signau" & PLZBmiss == 1
replace B_PLZ4 = 3658 if bürgerort == "Sigriswil" & PLZBmiss == 1
replace B_PLZ4 = 7514 if bürgerort == "Sils im Engadin/Segl" & PLZBmiss == 1
replace B_PLZ4 = 3907 if bürgerort == "Simplon" & PLZBmiss == 1
replace B_PLZ4 = 5643 if bürgerort == "Sins" & PLZBmiss == 1
replace B_PLZ4 = 1950 if bürgerort == "Sion" & PLZBmiss == 1
replace B_PLZ4 = 2577 if bürgerort == "Siselen" & PLZBmiss == 1
replace B_PLZ4 = 4334 if bürgerort == "Sisseln" & PLZBmiss == 1
replace B_PLZ4 = 1677 if bürgerort == "Siviriez" & PLZBmiss == 1
replace B_PLZ4 = 4503 if bürgerort == "Solothurn" & PLZBmiss == 1
replace B_PLZ4 = 1688 if bürgerort == "Sommentier" & PLZBmiss == 1
replace B_PLZ4 = 7176 if bürgerort == "Somvix" & PLZBmiss == 1
replace B_PLZ4 = 2605 if bürgerort == "Sonceboz-Sombeval" & PLZBmiss == 1
replace B_PLZ4 = 6637 if bürgerort == "Sonogno" & PLZBmiss == 1
replace B_PLZ4 = 2615 if bürgerort == "Sonvilier" & PLZBmiss == 1
replace B_PLZ4 = 1286 if bürgerort == "Soral" & PLZBmiss == 1
replace B_PLZ4 = 1642 if bürgerort == "Sorens" & PLZBmiss == 1
replace B_PLZ4 = 2716 if bürgerort == "Sornetan" & PLZBmiss == 1
replace B_PLZ4 = 1062 if bürgerort == "Sottens" & PLZBmiss == 1
replace B_PLZ4 = 2887 if bürgerort == "Soubey" & PLZBmiss == 1
replace B_PLZ4 = 2748 if bürgerort == "Souboz" & PLZBmiss == 1
replace B_PLZ4 = 9042 if bürgerort == "Speicher" & PLZBmiss == 1
replace B_PLZ4 = 3700 if bürgerort == "Spiez" & PLZBmiss == 1
replace B_PLZ4 = 8751 if bürgerort == "Spiringen" & PLZBmiss == 1
replace B_PLZ4 = 7435 if bürgerort == "Splügen" & PLZBmiss == 1
replace B_PLZ4 = 8957 if bürgerort == "Spreitenbach" & PLZBmiss == 1
replace B_PLZ4 = 1713 if bürgerort == "St. Antoni" & PLZBmiss == 1
replace B_PLZ4 = 7245 if bürgerort == "St. Antönien Ascharina" & PLZBmiss == 1
replace B_PLZ4 = 9020 if bürgerort == "St. Gallen" & PLZBmiss == 1
replace B_PLZ4 = 8727 if bürgerort == "St. Gallenkappel" & PLZBmiss == 1
replace B_PLZ4 = 9430 if bürgerort == "St. Margrethen" & PLZBmiss == 1
replace B_PLZ4 = 7116 if bürgerort == "St. Martin" & PLZBmiss == 1
replace B_PLZ4 = 7500 if bürgerort == "St. Moritz" & PLZBmiss == 1
replace B_PLZ4 = 3927 if bürgerort == "St. Niklaus" & PLZBmiss == 1
replace B_PLZ4 = 7028 if bürgerort == "St. Peter-Pagig" & PLZBmiss == 1
replace B_PLZ4 = 9127 if bürgerort == "St. Peterzell" & PLZBmiss == 1
replace B_PLZ4 = 1736 if bürgerort == "St. Silvester" & PLZBmiss == 1
replace B_PLZ4 = 3773 if bürgerort == "St. Stephan" & PLZBmiss == 1
replace B_PLZ4 = 1717 if bürgerort == "St. Ursen" & PLZBmiss == 1
replace B_PLZ4 = 6854 if bürgerort == "Stabio" & PLZBmiss == 1
replace B_PLZ4 = 8175 if bürgerort == "Stadel" & PLZBmiss == 1
replace B_PLZ4 = 5053 if bürgerort == "Staffelbach" & PLZBmiss == 1
replace B_PLZ4 = 3933 if bürgerort == "Staldenried" & PLZBmiss == 1
replace B_PLZ4 = 6370 if bürgerort == "Stans" & PLZBmiss == 1
replace B_PLZ4 = 4656 if bürgerort == "Starrkirch-Wil" & PLZBmiss == 1
replace B_PLZ4 = 3613 if bürgerort == "Steffisburg" & PLZBmiss == 1
replace B_PLZ4 = 3940 if bürgerort == "Steg" & PLZBmiss == 1
replace B_PLZ4 = 8260 if bürgerort == "Stein am Rhein" & PLZBmiss == 1
replace B_PLZ4 = 9323 if bürgerort == "Steinach" & PLZBmiss == 1
replace B_PLZ4 = 6422 if bürgerort == "Steinen" & PLZBmiss == 1
replace B_PLZ4 = 6416 if bürgerort == "Steinerberg" & PLZBmiss == 1
replace B_PLZ4 = 6312 if bürgerort == "Steinhausen" & PLZBmiss == 1
replace B_PLZ4 = 4556 if bürgerort == "Steinhof" & PLZBmiss == 1
replace B_PLZ4 = 8162 if bürgerort == "Steinmaur" & PLZBmiss == 1
replace B_PLZ4 = 8499 if bürgerort == "Sternenberg" & PLZBmiss == 1
replace B_PLZ4 = 5608 if bürgerort == "Stetten (AG)" & PLZBmiss == 1
replace B_PLZ4 = 7459 if bürgerort == "Stierva" & PLZBmiss == 1
replace B_PLZ4 = 3632 if bürgerort == "Stocken-Höfen" & PLZBmiss == 1
replace B_PLZ4 = 4802 if bürgerort == "Strengelbach" & PLZBmiss == 1
replace B_PLZ4 = 2557 if bürgerort == "Studen" & PLZBmiss == 1
replace B_PLZ4 = 8712 if bürgerort == "Stäfa" & PLZBmiss == 1
replace B_PLZ4 = 1433 if bürgerort == "Suchy" & PLZBmiss == 1
replace B_PLZ4 = 5034 if bürgerort == "Suhr" & PLZBmiss == 1
replace B_PLZ4 = 1036 if bürgerort == "Sullens" & PLZBmiss == 1
replace B_PLZ4 = 3457 if bürgerort == "Sumiswald" & PLZBmiss == 1
replace B_PLZ4 = 7472 if bürgerort == "Surava" & PLZBmiss == 1
replace B_PLZ4 = 6210 if bürgerort == "Sursee" & PLZBmiss == 1
replace B_PLZ4 = 7452 if bürgerort == "Surses" & PLZBmiss == 1
replace B_PLZ4 = 1437 if bürgerort == "Suscévaz" & PLZBmiss == 1
replace B_PLZ4 = 1510 if bürgerort == "Syens" & PLZBmiss == 1
replace B_PLZ4 = 1712 if bürgerort == "Tafers" & PLZBmiss == 1
replace B_PLZ4 = 7015 if bürgerort == "Tamins" & PLZBmiss == 1
replace B_PLZ4 = 7553 if bürgerort == "Tarasp" & PLZBmiss == 1
replace B_PLZ4 = 7422 if bürgerort == "Tartar" & PLZBmiss == 1
replace B_PLZ4 = 1180 if bürgerort == "Tartegnin" & PLZBmiss == 1
replace B_PLZ4 = 2720 if bürgerort == "Tavannes" & PLZBmiss == 1
replace B_PLZ4 = 5306 if bürgerort == "Tegerfelden" & PLZBmiss == 1
replace B_PLZ4 = 6598 if bürgerort == "Tenero-Contra" & PLZBmiss == 1
replace B_PLZ4 = 1734 if bürgerort == "Tentlingen" & PLZBmiss == 1
replace B_PLZ4 = 6654 if bürgerort == "Terre di Pedemonte" & PLZBmiss == 1
replace B_PLZ4 = 9423 if bürgerort == "Thal" & PLZBmiss == 1
replace B_PLZ4 = 8800 if bürgerort == "Thalwil" & PLZBmiss == 1
replace B_PLZ4 = 8236 if bürgerort == "Thayngen" & PLZBmiss == 1
replace B_PLZ4 = 4106 if bürgerort == "Therwil" & PLZBmiss == 1
replace B_PLZ4 = 2075 if bürgerort == "Thielle-Wavre" & PLZBmiss == 1
replace B_PLZ4 = 1410 if bürgerort == "Thierrens" & PLZBmiss == 1
replace B_PLZ4 = 3600 if bürgerort == "Thun" & PLZBmiss == 1
replace B_PLZ4 = 4922 if bürgerort == "Thunstetten" & PLZBmiss == 1
replace B_PLZ4 = 7430 if bürgerort == "Thusis" & PLZBmiss == 1
replace B_PLZ4 = 1226 if bürgerort == "Thônex" & PLZBmiss == 1
replace B_PLZ4 = 3367 if bürgerort == "Thörigen" & PLZBmiss == 1
replace B_PLZ4 = 7453 if bürgerort == "Tinizong" & PLZBmiss == 1
replace B_PLZ4 = 9555 if bürgerort == "Tobel-Tägerschen" & PLZBmiss == 1
replace B_PLZ4 = 3125 if bürgerort == "Toffen" & PLZBmiss == 1
replace B_PLZ4 = 1131 if bürgerort == "Tolochenaz" & PLZBmiss == 1
replace B_PLZ4 = 7418 if bürgerort == "Tomils" & PLZBmiss == 1
replace B_PLZ4 = 1749 if bürgerort == "Torny" & PLZBmiss == 1
replace B_PLZ4 = 1748 if bürgerort == "Torny-le-Grand" & PLZBmiss == 1
replace B_PLZ4 = 6808 if bürgerort == "Torricella-Taverne" & PLZBmiss == 1
replace B_PLZ4 = 3456 if bürgerort == "Trachselwald" & PLZBmiss == 1
replace B_PLZ4 = 2722 if bürgerort == "Tramelan" & PLZBmiss == 1
replace B_PLZ4 = 2105 if bürgerort == "Travers" & PLZBmiss == 1
replace B_PLZ4 = 1733 if bürgerort == "Treyvaux" & PLZBmiss == 1
replace B_PLZ4 = 6234 if bürgerort == "Triengen" & PLZBmiss == 1
replace B_PLZ4 = 4632 if bürgerort == "Trimbach" & PLZBmiss == 1
replace B_PLZ4 = 9043 if bürgerort == "Trogen" & PLZBmiss == 1
replace B_PLZ4 = 1256 if bürgerort == "Troinex" & PLZBmiss == 1
replace B_PLZ4 = 1872 if bürgerort == "Troistorrents" & PLZBmiss == 1
replace B_PLZ4 = 3556 if bürgerort == "Trub" & PLZBmiss == 1
replace B_PLZ4 = 3555 if bürgerort == "Trubschachen" & PLZBmiss == 1
replace B_PLZ4 = 7166 if bürgerort == "Trun" & PLZBmiss == 1
replace B_PLZ4 = 8465 if bürgerort == "Trüllikon" & PLZBmiss == 1
replace B_PLZ4 = 8856 if bürgerort == "Tuggen" & PLZBmiss == 1
replace B_PLZ4 = 7189 if bürgerort == "Tujetsch" & PLZBmiss == 1
replace B_PLZ4 = 8488 if bürgerort == "Turbenthal" & PLZBmiss == 1
replace B_PLZ4 = 5300 if bürgerort == "Turgi" & PLZBmiss == 1
replace B_PLZ4 = 3946 if bürgerort == "Turtmann-Unterems" & PLZBmiss == 1
replace B_PLZ4 = 2513 if bürgerort == "Twann-Tüscherz" & PLZBmiss == 1
replace B_PLZ4 = 5522 if bürgerort == "Tägerig" & PLZBmiss == 1
replace B_PLZ4 = 8274 if bürgerort == "Tägerwilen" & PLZBmiss == 1
replace B_PLZ4 = 2575 if bürgerort == "Täuffelen" & PLZBmiss == 1
replace B_PLZ4 = 1423 if bürgerort == "Tévenon" & PLZBmiss == 1
replace B_PLZ4 = 3923 if bürgerort == "Törbel" & PLZBmiss == 1
replace B_PLZ4 = 3182 if bürgerort == "Ueberstorf" & PLZBmiss == 1
replace B_PLZ4 = 5028 if bürgerort == "Ueken" & PLZBmiss == 1
replace B_PLZ4 = 3661 if bürgerort == "Uetendorf" & PLZBmiss == 1
replace B_PLZ4 = 5619 if bürgerort == "Uezwil" & PLZBmiss == 1
replace B_PLZ4 = 6153 if bürgerort == "Ufhusen" & PLZBmiss == 1
replace B_PLZ4 = 8142 if bürgerort == "Uitikon" & PLZBmiss == 1
replace B_PLZ4 = 3988 if bürgerort == "Ulrichen" & PLZBmiss == 1
replace B_PLZ4 = 3944 if bürgerort == "Unterbäch" & PLZBmiss == 1
replace B_PLZ4 = 9033 if bürgerort == "Untereggen" & PLZBmiss == 1
replace B_PLZ4 = 8845 if bürgerort == "Unteriberg" & PLZBmiss == 1
replace B_PLZ4 = 5726 if bürgerort == "Unterkulm" & PLZBmiss == 1
replace B_PLZ4 = 3614 if bürgerort == "Unterlangenegg" & PLZBmiss == 1
replace B_PLZ4 = 6465 if bürgerort == "Unterschächen" & PLZBmiss == 1
replace B_PLZ4 = 3800 if bürgerort == "Unterseen" & PLZBmiss == 1
replace B_PLZ4 = 8476 if bürgerort == "Unterstammheim" & PLZBmiss == 1
replace B_PLZ4 = 4916 if bürgerort == "Untersteckholz" & PLZBmiss == 1
replace B_PLZ4 = 7204 if bürgerort == "Untervaz" & PLZBmiss == 1
replace B_PLZ4 = 6314 if bürgerort == "Unterägeri" & PLZBmiss == 1
replace B_PLZ4 = 7114 if bürgerort == "Uors-Peiden" & PLZBmiss == 1
replace B_PLZ4 = 8901 if bürgerort == "Urdorf" & PLZBmiss == 1
replace B_PLZ4 = 9107 if bürgerort == "Urnäsch" & PLZBmiss == 1
replace B_PLZ4 = 4937 if bürgerort == "Ursenbach" & PLZBmiss == 1
replace B_PLZ4 = 1670 if bürgerort == "Ursy" & PLZBmiss == 1
replace B_PLZ4 = 3322 if bürgerort == "Urtenen" & PLZBmiss == 1
replace B_PLZ4 = 8616 if bürgerort == "Uster" & PLZBmiss == 1
replace B_PLZ4 = 3628 if bürgerort == "Uttigen" & PLZBmiss == 1
replace B_PLZ4 = 8592 if bürgerort == "Uttwil" & PLZBmiss == 1
replace B_PLZ4 = 3427 if bürgerort == "Utzenstorf" & PLZBmiss == 1
replace B_PLZ4 = 8730 if bürgerort == "Uznach" & PLZBmiss == 1
replace B_PLZ4 = 9240 if bürgerort == "Uzwil" & PLZBmiss == 1
replace B_PLZ4 = 6833 if bürgerort == "Vacallo" & PLZBmiss == 1
replace B_PLZ4 = 7535 if bürgerort == "Val Müstair" & PLZBmiss == 1
replace B_PLZ4 = 2829 if bürgerort == "Val Terbi" & PLZBmiss == 1
replace B_PLZ4 = 1873 if bürgerort == "Val-d'Illiez" & PLZBmiss == 1
replace B_PLZ4 = 1654 if bürgerort == "Val-de-Charmey" & PLZBmiss == 1
replace B_PLZ4 = 2057 if bürgerort == "Val-de-Ruz" & PLZBmiss == 1
replace B_PLZ4 = 2114 if bürgerort == "Val-de-Travers" & PLZBmiss == 1
replace B_PLZ4 = 2735 if bürgerort == "Valbirse" & PLZBmiss == 1
replace B_PLZ4 = 1536 if bürgerort == "Valbroye" & PLZBmiss == 1
replace B_PLZ4 = 6951 if bürgerort == "Valcolla" & PLZBmiss == 1
replace B_PLZ4 = 1441 if bürgerort == "Valeyres-sous-Montagny" & PLZBmiss == 1
replace B_PLZ4 = 1358 if bürgerort == "Valeyres-sous-Rances" & PLZBmiss == 1
replace B_PLZ4 = 1565 if bürgerort == "Vallon" & PLZBmiss == 1
replace B_PLZ4 = 1337 if bürgerort == "Vallorbe" & PLZBmiss == 1
replace B_PLZ4 = 7116 if bürgerort == "Vals" & PLZBmiss == 1
replace B_PLZ4 = 7558 if bürgerort == "Valsot" & PLZBmiss == 1
replace B_PLZ4 = 3953 if bürgerort == "Varen" & PLZBmiss == 1
replace B_PLZ4 = 2537 if bürgerort == "Vauffelin" & PLZBmiss == 1
replace B_PLZ4 = 1325 if bürgerort == "Vaulion" & PLZBmiss == 1
replace B_PLZ4 = 1627 if bürgerort == "Vaulruz" & PLZBmiss == 1
replace B_PLZ4 = 2028 if bürgerort == "Vaumarcus" & PLZBmiss == 1
replace B_PLZ4 = 7082 if bürgerort == "Vaz/Obervaz" & PLZBmiss == 1
replace B_PLZ4 = 3068 if bürgerort == "Vechigen" & PLZBmiss == 1
replace B_PLZ4 = 2943 if bürgerort == "Vendlincourt" & PLZBmiss == 1
replace B_PLZ4 = 3973 if bürgerort == "Venthône" & PLZBmiss == 1
replace B_PLZ4 = 6664 if bürgerort == "Vergeletto" & PLZBmiss == 1
replace B_PLZ4 = 6992 if bürgerort == "Vernate" & PLZBmiss == 1
replace B_PLZ4 = 1904 if bürgerort == "Vernayaz" & PLZBmiss == 1
replace B_PLZ4 = 1219 if bürgerort == "Vernier" & PLZBmiss == 1
replace B_PLZ4 = 1290 if bürgerort == "Versoix" & PLZBmiss == 1
replace B_PLZ4 = 1800 if bürgerort == "Vevey" & PLZBmiss == 1
replace B_PLZ4 = 1988 if bürgerort == "Vex" & PLZBmiss == 1
replace B_PLZ4 = 1255 if bürgerort == "Veyrier" & PLZBmiss == 1
replace B_PLZ4 = 1993 if bürgerort == "Veysonnaz" & PLZBmiss == 1
replace B_PLZ4 = 1820 if bürgerort == "Veytaux" & PLZBmiss == 1
replace B_PLZ4 = 6943 if bürgerort == "Vezia" & PLZBmiss == 1
replace B_PLZ4 = 6921 if bürgerort == "Vico Morcote" & PLZBmiss == 1
replace B_PLZ4 = 1679 if bürgerort == "Villaraboud" & PLZBmiss == 1
replace B_PLZ4 = 1583 if bürgerort == "Villarepos" & PLZBmiss == 1
replace B_PLZ4 = 1685 if bürgerort == "Villariaz" & PLZBmiss == 1
replace B_PLZ4 = 1682 if bürgerort == "Villars-Bramard" & PLZBmiss == 1
replace B_PLZ4 = 1423 if bürgerort == "Villars-Burquin" & PLZBmiss == 1
replace B_PLZ4 = 1040 if bürgerort == "Villars-le-Terroir" & PLZBmiss == 1
replace B_PLZ4 = 1752 if bürgerort == "Villars-sur-Glâne" & PLZBmiss == 1
replace B_PLZ4 = 1651 if bürgerort == "Villarvolard" & PLZBmiss == 1
replace B_PLZ4 = 1690 if bürgerort == "Villaz-Saint-Pierre" & PLZBmiss == 1
replace B_PLZ4 = 1091 if bürgerort == "Villette (Lavaux)" & PLZBmiss == 1
replace B_PLZ4 = 5234 if bürgerort == "Villigen" & PLZBmiss == 1
replace B_PLZ4 = 5613 if bürgerort == "Villmergen" & PLZBmiss == 1
replace B_PLZ4 = 1694 if bürgerort == "Villorsonnens" & PLZBmiss == 1
replace B_PLZ4 = 7324 if bürgerort == "Vilters" & PLZBmiss == 1
replace B_PLZ4 = 3234 if bürgerort == "Vinelz" & PLZBmiss == 1
replace B_PLZ4 = 1895 if bürgerort == "Vionnaz" & PLZBmiss == 1
replace B_PLZ4 = 3930 if bürgerort == "Visp" & PLZBmiss == 1
replace B_PLZ4 = 3932 if bürgerort == "Visperterminen" & PLZBmiss == 1
replace B_PLZ4 = 3961 if bürgerort == "Vissoie" & PLZBmiss == 1
replace B_PLZ4 = 6632 if bürgerort == "Vogorno" & PLZBmiss == 1
replace B_PLZ4 = 8459 if bürgerort == "Volken" & PLZBmiss == 1
replace B_PLZ4 = 8604 if bürgerort == "Volketswil" & PLZBmiss == 1
replace B_PLZ4 = 1927 if bürgerort == "Vollèges" & PLZBmiss == 1
replace B_PLZ4 = 4803 if bürgerort == "Vordemwald" & PLZBmiss == 1
replace B_PLZ4 = 8857 if bürgerort == "Vorderthal" & PLZBmiss == 1
replace B_PLZ4 = 1896 if bürgerort == "Vouvry" & PLZBmiss == 1
replace B_PLZ4 = 7149 if bürgerort == "Vrin" & PLZBmiss == 1
replace B_PLZ4 = 1628 if bürgerort == "Vuadens" & PLZBmiss == 1
replace B_PLZ4 = 1509 if bürgerort == "Vucherens" & PLZBmiss == 1
replace B_PLZ4 = 1134 if bürgerort == "Vufflens-le-Château" & PLZBmiss == 1
replace B_PLZ4 = 1431 if bürgerort == "Vugelles-La Mothe" & PLZBmiss == 1
replace B_PLZ4 = 1486 if bürgerort == "Vuissens" & PLZBmiss == 1
replace B_PLZ4 = 1687 if bürgerort == "Vuisternens-devant-Romont" & PLZBmiss == 1
replace B_PLZ4 = 1696 if bürgerort == "Vuisternens-en-Ogoz" & PLZBmiss == 1
replace B_PLZ4 = 1085 if bürgerort == "Vulliens" & PLZBmiss == 1
replace B_PLZ4 = 1585 if bürgerort == "Vully-les-Lacs" & PLZBmiss == 1
replace B_PLZ4 = 1963 if bürgerort == "Vétroz" & PLZBmiss == 1
replace B_PLZ4 = 3618 if bürgerort == "Wachseldorn" & PLZBmiss == 1
replace B_PLZ4 = 4246 if bürgerort == "Wahlen" & PLZBmiss == 1
replace B_PLZ4 = 3152 if bürgerort == "Wahlern" & PLZBmiss == 1
replace B_PLZ4 = 6318 if bürgerort == "Walchwil" & PLZBmiss == 1
replace B_PLZ4 = 4437 if bürgerort == "Waldenburg" & PLZBmiss == 1
replace B_PLZ4 = 9205 if bürgerort == "Waldkirch" & PLZBmiss == 1
replace B_PLZ4 = 9104 if bürgerort == "Waldstatt" & PLZBmiss == 1
replace B_PLZ4 = 8880 if bürgerort == "Walenstadt" & PLZBmiss == 1
replace B_PLZ4 = 3512 if bürgerort == "Walkringen" & PLZBmiss == 1
replace B_PLZ4 = 4323 if bürgerort == "Wallbach" & PLZBmiss == 1
replace B_PLZ4 = 1784 if bürgerort == "Wallenried" & PLZBmiss == 1
replace B_PLZ4 = 8304 if bürgerort == "Wallisellen" & PLZBmiss == 1
replace B_PLZ4 = 3380 if bürgerort == "Walliswil bei Niederbipp" & PLZBmiss == 1
replace B_PLZ4 = 3377 if bürgerort == "Walliswil bei Wangen" & PLZBmiss == 1
replace B_PLZ4 = 3272 if bürgerort == "Walperswil" & PLZBmiss == 1
replace B_PLZ4 = 8468 if bürgerort == "Waltalingen" & PLZBmiss == 1
replace B_PLZ4 = 7158 if bürgerort == "Waltensburg/Vuorz" & PLZBmiss == 1
replace B_PLZ4 = 5622 if bürgerort == "Waltenschwil" & PLZBmiss == 1
replace B_PLZ4 = 9428 if bürgerort == "Walzenhausen" & PLZBmiss == 1
replace B_PLZ4 = 3380 if bürgerort == "Wangen an der Aare" & PLZBmiss == 1
replace B_PLZ4 = 4612 if bürgerort == "Wangen bei Olten" & PLZBmiss == 1
replace B_PLZ4 = 3374 if bürgerort == "Wangenried" & PLZBmiss == 1
replace B_PLZ4 = 3372 if bürgerort == "Wanzwil" & PLZBmiss == 1
replace B_PLZ4 = 9479 if bürgerort == "Wartau" & PLZBmiss == 1
replace B_PLZ4 = 8532 if bürgerort == "Warth-Weiningen" & PLZBmiss == 1
replace B_PLZ4 = 6485 if bürgerort == "Wassen" & PLZBmiss == 1
replace B_PLZ4 = 3665 if bürgerort == "Wattenwil" & PLZBmiss == 1
replace B_PLZ4 = 9631 if bürgerort == "Wattwil" & PLZBmiss == 1
replace B_PLZ4 = 8872 if bürgerort == "Weesen" & PLZBmiss == 1
replace B_PLZ4 = 6353 if bürgerort == "Weggis" & PLZBmiss == 1
replace B_PLZ4 = 8570 if bürgerort == "Weinfelden" & PLZBmiss == 1
replace B_PLZ4 = 8484 if bürgerort == "Weisslingen" & PLZBmiss == 1
replace B_PLZ4 = 4716 if bürgerort == "Welschenrohr" & PLZBmiss == 1
replace B_PLZ4 = 3251 if bürgerort == "Wengi" & PLZBmiss == 1
replace B_PLZ4 = 4493 if bürgerort == "Wenslingen" & PLZBmiss == 1
replace B_PLZ4 = 6106 if bürgerort == "Werthenstein" & PLZBmiss == 1
replace B_PLZ4 = 5430 if bürgerort == "Wettingen" & PLZBmiss == 1
replace B_PLZ4 = 8907 if bürgerort == "Wettswil am Albis" & PLZBmiss == 1
replace B_PLZ4 = 3114 if bürgerort == "Wichtrach" & PLZBmiss == 1
replace B_PLZ4 = 8967 if bürgerort == "Widen" & PLZBmiss == 1
replace B_PLZ4 = 9443 if bürgerort == "Widnau" & PLZBmiss == 1
replace B_PLZ4 = 4537 if bürgerort == "Wiedlisbach" & PLZBmiss == 1
replace B_PLZ4 = 8546 if bürgerort == "Wiesendangen" & PLZBmiss == 1
replace B_PLZ4 = 3053 if bürgerort == "Wiggiswil" & PLZBmiss == 1
replace B_PLZ4 = 4806 if bürgerort == "Wikon" & PLZBmiss == 1
replace B_PLZ4 = 8492 if bürgerort == "Wila" & PLZBmiss == 1
replace B_PLZ4 = 8218 if bürgerort == "Wilchingen" & PLZBmiss == 1
replace B_PLZ4 = 8489 if bürgerort == "Wildberg" & PLZBmiss == 1
replace B_PLZ4 = 9658 if bürgerort == "Wildhaus" & PLZBmiss == 1
replace B_PLZ4 = 9657 if bürgerort == "Wildhaus-Alt St. Johann" & PLZBmiss == 1
replace B_PLZ4 = 3918 if bürgerort == "Wiler (Lötschen)" & PLZBmiss == 1
replace B_PLZ4 = 3428 if bürgerort == "Wiler bei Utzenstorf" & PLZBmiss == 1
replace B_PLZ4 = 5058 if bürgerort == "Wiliberg" & PLZBmiss == 1
replace B_PLZ4 = 6130 if bürgerort == "Willisau" & PLZBmiss == 1
replace B_PLZ4 = 6130 if bürgerort == "Willisau Stadt" & PLZBmiss == 1
replace B_PLZ4 = 3752 if bürgerort == "Wimmis" & PLZBmiss == 1
replace B_PLZ4 = 5210 if bürgerort == "Windisch" & PLZBmiss == 1
replace B_PLZ4 = 8185 if bürgerort == "Winkel" & PLZBmiss == 1
replace B_PLZ4 = 4451 if bürgerort == "Wintersingen" & PLZBmiss == 1
replace B_PLZ4 = 8400 if bürgerort == "Winterthur" & PLZBmiss == 1
replace B_PLZ4 = 5463 if bürgerort == "Wislikofen" & PLZBmiss == 1
replace B_PLZ4 = 9300 if bürgerort == "Wittenbach" & PLZBmiss == 1
replace B_PLZ4 = 4108 if bürgerort == "Witterswil" & PLZBmiss == 1
replace B_PLZ4 = 5064 if bürgerort == "Wittnau" & PLZBmiss == 1
replace B_PLZ4 = 5610 if bürgerort == "Wohlen (AG)" & PLZBmiss == 1
replace B_PLZ4 = 3043 if bürgerort == "Wohlen bei Bern" & PLZBmiss == 1
replace B_PLZ4 = 6387 if bürgerort == "Wolfenschiessen" & PLZBmiss == 1
replace B_PLZ4 = 9427 if bürgerort == "Wolfhalden" & PLZBmiss == 1
replace B_PLZ4 = 4704 if bürgerort == "Wolfisberg" & PLZBmiss == 1
replace B_PLZ4 = 4628 if bürgerort == "Wolfwil" & PLZBmiss == 1
replace B_PLZ4 = 6110 if bürgerort == "Wolhusen" & PLZBmiss == 1
replace B_PLZ4 = 3077 if bürgerort == "Worb" & PLZBmiss == 1
replace B_PLZ4 = 3252 if bürgerort == "Worben" & PLZBmiss == 1
replace B_PLZ4 = 4923 if bürgerort == "Wynau" & PLZBmiss == 1
replace B_PLZ4 = 3474 if bürgerort == "Wynigen" & PLZBmiss == 1
replace B_PLZ4 = 4954 if bürgerort == "Wyssachen" & PLZBmiss == 1
replace B_PLZ4 = 8820 if bürgerort == "Wädenswil" & PLZBmiss == 1
replace B_PLZ4 = 9175 if bürgerort == "Wünnewil-Flamatt" & PLZBmiss == 1
replace B_PLZ4 = 5303 if bürgerort == "Würenlingen" & PLZBmiss == 1
replace B_PLZ4 = 1432 if bürgerort == "Yverdon-les-Bains" & PLZBmiss == 1
replace B_PLZ4 = 1462 if bürgerort == "Yvonand" & PLZBmiss == 1
replace B_PLZ4 = 1853 if bürgerort == "Yvorne" & PLZBmiss == 1
replace B_PLZ4 = 3309 if bürgerort == "Zauggenried" & PLZBmiss == 1
replace B_PLZ4 = 4495 if bürgerort == "Zeglingen" & PLZBmiss == 1
replace B_PLZ4 = 5079 if bürgerort == "Zeihen" & PLZBmiss == 1
replace B_PLZ4 = 3920 if bürgerort == "Zermatt" & PLZBmiss == 1
replace B_PLZ4 = 7543 if bürgerort == "Zernez" & PLZBmiss == 1
replace B_PLZ4 = 4417 if bürgerort == "Ziefen" & PLZBmiss == 1
replace B_PLZ4 = 8588 if bürgerort == "Zihlschlacht-Sitterdorf" & PLZBmiss == 1
replace B_PLZ4 = 7432 if bürgerort == "Zillis-Reischen" & PLZBmiss == 1
replace B_PLZ4 = 7205 if bürgerort == "Zizers" & PLZBmiss == 1
replace B_PLZ4 = 4800 if bürgerort == "Zofingen" & PLZBmiss == 1
replace B_PLZ4 = 8702 if bürgerort == "Zollikon" & PLZBmiss == 1
replace B_PLZ4 = 4528 if bürgerort == "Zuchwil" & PLZBmiss == 1
replace B_PLZ4 = 6301 if bürgerort == "Zug" & PLZBmiss == 1
replace B_PLZ4 = 1719 if bürgerort == "Zumholz" & PLZBmiss == 1
replace B_PLZ4 = 8126 if bürgerort == "Zumikon" & PLZBmiss == 1
replace B_PLZ4 = 4455 if bürgerort == "Zunzgen" & PLZBmiss == 1
replace B_PLZ4 = 5330 if bürgerort == "Zurzach" & PLZBmiss == 1
replace B_PLZ4 = 3770 if bürgerort == "Zweisimmen" & PLZBmiss == 1
replace B_PLZ4 = 3907 if bürgerort == "Zwischbergen" & PLZBmiss == 1
replace B_PLZ4 = 3532 if bürgerort == "Zäziwil" & PLZBmiss == 1
replace B_PLZ4 = 8004 if bürgerort == "Zürich" & PLZBmiss == 1

/*
// (still) Missing Buergerort PLZ & Nationalität==Schweiz
preserve
gen nat_ch = 0
replace nat_ch =1 if nationalität == "Schweiz"
keep if nat_ch==1 & PLZBmiss==1 & B_PLZ4==.
gen count = 1
collapse (sum) count, by(bürgerort)
sort bürgerort
save "$path\02_Processed_data\11_Directors_1994_2018/Bisnode_missingBuergerPLZ_temp2.dta", replace
restore

gen c_persID = 1
bysort personenid: egen obspersID_B = total(c_persID)
gen c_missPLZ = 0
replace c_missPLZ = 1 if B_PLZ4 == .
bysort personenid: egen obsmissPLZ = total(c_missPLZ)
*br * if obsmissPLZ < obspersID_B & obsmissPLZ != 0  // are there observations within Person ID that have non-missing PLZ -- > 0 Obs.
drop c_persID c_missPLZ obspersID_B obsmissPLZ
*/

// Replace PLZ of missing Bürgerort-PLZ if nationalität=="Schweiz" for observations with >= 8 occurances
gen hPLZ4_1 =.
replace hPLZ4_1 = 1717 if bürgerort == "St-Ours" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4114 if bürgerort == "Hofstetten" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5452 if bürgerort == "Oberrohrdorf und Gossau SG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4000 if bürgerort == "Bâle" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7482 if bürgerort == "Bergün" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3906 if bürgerort == "Saas-Fee" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3415 if bürgerort == "Hasle BE" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8802 if bürgerort == "Kilchberg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5610 if bürgerort == "Wohlen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5417 if bürgerort == "Untersiggent" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1611 if bürgerort == "Crêt" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3185 if bürgerort == "Schmitten" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5413 if bürgerort == "Birmenstorf AG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3512 if bürgerort == "Walkringen/BE" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9000 if bürgerort == "St-Gall et Stein" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9463 if bürgerort == "Oberriet SG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8104 if bürgerort == "Weiningen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1228 if bürgerort == "Plan-les-Quates" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3625 if bürgerort == "Schwendi / Heiligenschwendi" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6418 if bürgerort == "Rothenthu" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3178 if bürgerort == "Bösingen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7175 if bürgerort == "Sumvigt GR" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6052 if bürgerort == "Heriswil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8802 if bürgerort == "Kilchberg ZH" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3271 if bürgerort == "Radelfingen bei Aarberg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6926 if bürgerort == "Collinda d'Oro" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8855 if bürgerort == "Wangen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9411 if bürgerort == "Reute" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8418 if bürgerort == "Unterschlatt" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9001 if bürgerort == "St. Gallen-Tablat" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2500 if bürgerort == "Biel" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5445 if bürgerort == "Heggenswil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8303 if bürgerort == "Basserdorf" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9606 if bürgerort == "Bütschwil-Ganterswil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1712 if bürgerort == "Tavel" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1482 if bürgerort == "Cugy" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9204 if bürgerort == "Andwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1684 if bürgerort == "Mézières" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1694 if bürgerort == "Orsonnen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9533 if bürgerort == "Kirchberg SG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9200 if bürgerort == "Gossau" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5723 if bürgerort == "Teufenthal" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8251 if bürgerort == "Kleinandelfing" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6026 if bürgerort == "Rain/LU" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5733 if bürgerort == "Leimbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7320 if bürgerort == "Sargen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4853 if bürgerort == "Murgenthaler" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1218 if bürgerort == "Grand-Saconnex" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1786 if bürgerort == "Haut-Vully" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8508 if bürgerort == "Gündelhart-Hörhausen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6541 if bürgerort == "Santa Maria" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8200 if bürgerort == "Sciaffusa" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6460 if bürgerort == "Altdorf" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1025 if bürgerort == "Saint-Sulpice VD" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1024 if bürgerort == "Ecublens" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4533 if bürgerort == "Niederwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8570 if bürgerort == "Weerswilen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6949 if bürgerort == "Camano/TI" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9444 if bürgerort == "Diepoldsau-Schmitter" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4952 if bürgerort == "Eriswil BE" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4633 if bürgerort == "Ifenthal" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1653 if bürgerort == "Châtel-sur-Montsalv" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8580 if bürgerort == "Räuchlisberg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3186 if bürgerort == "Guin" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9127 if bürgerort == "St.Peterzell" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8718 if bürgerort == "Schänis-Rufi" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1344 if bürgerort == "L' Abbaye" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1700 if bürgerort == "Friburgo" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7184 if bürgerort == "Medel" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1713 if bürgerort == "St-Antoni" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2500 if bürgerort == "Bözingen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1000 if bürgerort == "Losanna" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1864 if bürgerort == "Ormond-Dessus" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6130 if bürgerort == "Willisau-Stadt" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5054 if bürgerort == "Mooslerau" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8400 if bürgerort == "Winterth" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9514 if bürgerort == "Wuppenau" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1687 if bürgerort == "Vuisternens-dev-Romont" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2500 if bürgerort == "Bienne" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3974 if bürgerort == "Mollens" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1734 if bürgerort == "Tinterin" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1731 if bürgerort == "Ependes" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8344 if bürgerort == "Baeretswil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4553 if bürgerort == "Subinge" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8890 if bürgerort == "Flums-Kleinberg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1091 if bürgerort == "Villette" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9053 if bürgerort == "Teufen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1032 if bürgerort == "Romanel" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3550 if bürgerort == "Langnau" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5332 if bürgerort == "Rekingen AG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8707 if bürgerort == "Uetikon am See" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1473 if bürgerort == "Châtillon" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8630 if bürgerort == "Rüti" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1724 if bürgerort == "Le Mouret" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3297 if bürgerort == "Leuzingen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1718 if bürgerort == "Dirlaret" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4657 if bürgerort == "Dullikon" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8756 if bürgerort == "Glarona Sud" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6463 if bürgerort == "Bürglen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2616 if bürgerort == "Renan" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5734 if bürgerort == "Reinach AG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4914 if bürgerort == "Roggwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9200 if bürgerort == "Gossau SG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8585 if bürgerort == "Dünnershaus" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4573 if bürgerort == "Lohn" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6166 if bürgerort == "Hasle" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4942 if bürgerort == "Walterswil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1864 if bürgerort == "Les Ormonts" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3953 if bürgerort == "Loèche" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4900 if bürgerort == "Schoren" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5085 if bürgerort == "Sulz" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2022 if bürgerort == "La Grande Béroche" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1345 if bürgerort == "Lieu" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1717 if bürgerort == "St-Ursen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7463 if bürgerort == "Riom-Parsonz" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7180 if bürgerort == "Disentis" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8306 if bürgerort == "Wangen-Brüttiselle" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9035 if bürgerort == "Grub" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6130 if bürgerort == "Willisau-Land" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8306 if bürgerort == "Wangen-Brüttisellen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6182 if bürgerort == "Eschholzmatt" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7324 if bürgerort == "Vilters-Wangs" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9500 if bürgerort == "Wil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9053 if bürgerort == "Teufen AR" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8717 if bürgerort == "Benken" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5453 if bürgerort == "Remetschwil AG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1614 if bürgerort == "Granges" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3600 if bürgerort == "Strättligen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6634 if bürgerort == "Brione Verzasca" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1844 if bürgerort == "Villeneuve" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5106 if bürgerort == "Weltheim" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4656 if bürgerort == "Starrkirch" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3855 if bürgerort == "Brienz" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2360 if bürgerort == "Bémont" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1025 if bürgerort == "St-Sulpice" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9205 if bürgerort == "Waldkirch-Bernhardzell" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9533 if bürgerort == "Kirchberg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8890 if bürgerort == "Flums-Grossberg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6106 if bürgerort == "Werthestein" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5112 if bürgerort == "Thalheim AG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3236 if bürgerort == "Champion" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7250 if bürgerort == "Klosters-Serneus" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8585 if bürgerort == "Birwinken" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1563 if bürgerort == "Dompierre" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6083 if bürgerort == "Reuti" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4663 if bürgerort == "Arburg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6144 if bürgerort == "Zell" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1727 if bürgerort == "Corpataux-Magnedens" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1489 if bürgerort == "Montborget" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8600 if bürgerort == "Stettbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1072 if bürgerort == "Forel" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8524 if bürgerort == "Uesslingen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1934 if bürgerort == "Bagnens" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1717 if bürgerort == "St-Ours et Fribourg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1754 if bürgerort == "Avry" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6052 if bürgerort == "Hergiswil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1031 if bürgerort == "Mex" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3942 if bürgerort == "Niedergestein" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8733 if bürgerort == "Eschenbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5112 if bürgerort == "Thalheim" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3508 if bürgerort == "Arni" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2512 if bürgerort == "Tüscherz" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3178 if bürgerort == "Bösigen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2855 if bürgerort == "Sceut" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1934 if bürgerort == "Bagne" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6300 if bürgerort == "Zugo" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7525 if bürgerort == "Scanf" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8733 if bürgerort == "Eschenbach SG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7180 if bürgerort == "Mustér" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5413 if bürgerort == "Birmenstorf" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2336 if bürgerort == "Bois" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1717 if bürgerort == "St.Ours" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1724 if bürgerort == "Montécu" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4000 if bürgerort == "Basel-Stadt" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4500 if bürgerort == "Soleure" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6500 if bürgerort == "Bellinzone" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3512 if bürgerort == "Walringen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3086 if bürgerort == "Wald" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7175 if bürgerort == "Sumvigt" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2103 if bürgerort == "Le Val-de-Travers" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5736 if bürgerort == "Burg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2543 if bürgerort == "Lengnau" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9524 if bürgerort == "Zuzwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1341 if bürgerort == "Chenit" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1400 if bürgerort == "Yverdon Les Bains" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9642 if bürgerort == "Ebnat" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8887 if bürgerort == "Mels-Weisstannen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4656 if bürgerort == "Starrkirch SO" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7418 if bürgerort == "Feldis" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4153 if bürgerort == "Reinach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3998 if bürgerort == "Reckingen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1740 if bürgerort == "Neyruz" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5620 if bürgerort == "Bremgarten" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8824 if bürgerort == "Schönenberg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7175 if bürgerort == "Somvix GR" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8750 if bürgerort == "Glaris et Amden" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3984 if bürgerort == "Fischertal" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4500 if bürgerort == "Soletta" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6133 if bürgerort == "Lucerna e Hergiswil/LU" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1784 if bürgerort == "Wagenried" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1741 if bürgerort == "Cottens" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1774 if bürgerort == "Montagny" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3216 if bürgerort == "Ried" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3114 if bürgerort == "Niederwichtrach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2504 if bürgerort == "Bönzingen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8000 if bürgerort == "Zurich" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 2022 if bürgerort == "La Grande-Béroche" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8750 if bürgerort == "Glarona" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4616 if bürgerort == "Kappel" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8355 if bürgerort == "Aardorf TG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4000 if bürgerort == "Bâle-Ville" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6808 if bürgerort == "Torricella- Taverne" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8620 if bürgerort == "Wetzikon" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6221 if bürgerort == "Rickenbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9463 if bürgerort == "Oberriet" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3255 if bürgerort == "Rapperswil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9437 if bürgerort == "Marbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3852 if bürgerort == "Ringgenberg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1716 if bürgerort == "Planfayon" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1091 if bürgerort == "Cully-Lutry-Puidoux" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9524 if bürgerort == "Zuziwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1085 if bürgerort == "Vuillens" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4229 if bürgerort == "Beinwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3792 if bürgerort == "Gessenay" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8107 if bürgerort == "Buchs" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 9542 if bürgerort == "Münchwilen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1625 if bürgerort == "Sales" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1969 if bürgerort == "St-Martin" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6112 if bürgerort == "Doppelschwand" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8585 if bürgerort == "Dünnershaus TG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4147 if bürgerort == "Aesch" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5015 if bürgerort == "Obererlinsbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 4436 if bürgerort == "Oberdorf" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1680 if bürgerort == "Romont" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 7110 if bürgerort == "Suraua" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5630 if bürgerort == "Muri" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 3267 if bürgerort == "Seedorf" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8700 if bürgerort == "Küsnacht ZH" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8307 if bürgerort == "Les Kybourg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 8762 if bürgerort == "Schwanden" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 5106 if bürgerort == "Veltheim" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1227 if bürgerort == "Carouge" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6345 if bürgerort == "Heuheim" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1786 if bürgerort == "Vully-le-Bas" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 6802 if bürgerort == "Montecenri" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_1 = 1784 if bürgerort == "Esserts" & B_PLZ4 == . & nationalität == "Schweiz"
gen hPLZ4_2 =.
replace hPLZ4_2 = 3911 if bürgerort == "Ried" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 9655 if bürgerort == "St-Gall et Stein" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1700 if bürgerort == "St-Ours et Fribourg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 5426 if bürgerort == "Lengnau" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 9312 if bürgerort == "Heggenswil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 4462 if bürgerort == "Rickenbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8463 if bürgerort == "Benken" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 5746 if bürgerort == "Walterswil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 3038 if bürgerort == "Bremgarten" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1673 if bürgerort == "Ecublens" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8478 if bürgerort == "Thalheim" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8000 if bürgerort == "Leimbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 7493 if bürgerort == "Schmitten" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1053 if bürgerort == "Cugy" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 7536 if bürgerort == "Santa Maria" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 7084 if bürgerort == "Brienz" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 2843 if bürgerort == "Châtillon" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 3303 if bürgerort == "Zuzwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 6370 if bürgerort == "Oberdorf" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 6284 if bürgerort == "Sulz" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 3422 if bürgerort == "Kirchberg" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8782 if bürgerort == "Rüti" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1609 if bürgerort == "St-Martin" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8243 if bürgerort == "Altdorf" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 3033 if bürgerort == "Wohlen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1731 if bürgerort == "Sales" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 4333 if bürgerort == "Münchwilen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8196 if bürgerort == "Wil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 5033 if bürgerort == "Buchs" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 3303 if bürgerort == "Zuziwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1434 if bürgerort == "Ependes" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8905 if bürgerort == "Arni" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8873 if bürgerort == "Glaris et Amden" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 6000 if bürgerort == "Lucerna e Hergiswil/LU" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1247 if bürgerort == "Villette" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1682 if bürgerort == "Dompierre" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 5524 if bürgerort == "Niederwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 3380 if bürgerort == "Wangen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1818 if bürgerort == "Tavel" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1527 if bürgerort == "Villeneuve" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 9325 if bürgerort == "Roggwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8575 if bürgerort == "Bürglen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1122 if bürgerort == "Romanel" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 7433 if bürgerort == "Lohn" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1095 if bürgerort == "Cully-Lutry-Puidoux" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 9200 if bürgerort == "Oberrohrdorf und Gossau SG" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 6196 if bürgerort == "Marbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 3074 if bürgerort == "Muri" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1862 if bürgerort == "Les Ormonts" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1083 if bürgerort == "Mézières" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 6287 if bürgerort == "Aesch" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 2538 if bürgerort == "Romont" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 2540 if bürgerort == "Granges" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8625 if bürgerort == "Gossau" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8487 if bürgerort == "Zell" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 9034 if bürgerort == "Grub" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1890 if bürgerort == "Mex" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 5734 if bürgerort == "Reinach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8640 if bürgerort == "Rapperswil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 2123 if bürgerort == "St-Sulpice" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8636 if bürgerort == "Wald" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1116 if bürgerort == "Cottens" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 6274 if bürgerort == "Eschenbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 6462 if bürgerort == "Seedorf" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 1475 if bürgerort == "Forel" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 5637 if bürgerort == "Beinwil" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8354 if bürgerort == "Hofstetten" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_2 = 8585 if bürgerort == "Andwil" & B_PLZ4 == . & nationalität == "Schweiz"
gen hPLZ4_3 =.
replace hPLZ4_3 = 4613 if bürgerort == "Rickenbach" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 3960 if bürgerort == "Granges" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 6211 if bürgerort == "Buchs" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 3986 if bürgerort == "Ried" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 4612 if bürgerort == "Wangen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 3858 if bürgerort == "Hofstetten" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 8904 if bürgerort == "Aesch" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 1226 if bürgerort == "Villette" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 9044 if bürgerort == "Wald" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 4515 if bürgerort == "Oberdorf" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 1070 if bürgerort == "Cully-Lutry-Puidoux" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 8235 if bürgerort == "Lohn" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_3 = 1818 if bürgerort == "Sales" & B_PLZ4 == . & nationalität == "Schweiz"
gen hPLZ4_4 =.
replace hPLZ4_4 = 9470 if bürgerort == "Buchs" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_4 = 8306 if bürgerort == "Wangen" & B_PLZ4 == . & nationalität == "Schweiz"
replace hPLZ4_4 = 9532 if bürgerort == "Rickenbach" & B_PLZ4 == . & nationalität == "Schweiz"
gen hPLZ4_5 =.
replace hPLZ4_5 = 8545 if bürgerort == "Rickenbach" & B_PLZ4 == . & nationalität == "Schweiz"
gen help1=_n
reshape long hPLZ4_, i(help1) j(obs)
rename hPLZ4_ hPLZ4
drop if hPLZ4 == . & obs > 1
replace B_PLZ4 = hPLZ4 if B_PLZ4 == . & hPLZ4 != . 
rename obs Bürgerort_PLZopt1
drop hPLZ4 help1

// Replace missing Bürgerort PLZ if Nationalität != "Schweiz" but Bürgerort is in Switzerland (for observations with >= 10 occurances)
// Missing Bürgerort occurs with nationalität == ""
gen hPLZ4_1 = .
replace hPLZ4_1 = 3153 if bürgerort == "Rüeschegg" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 4942 if bürgerort == "Walterswil" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 6274 if bürgerort == "Eschenbach" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1892 if bürgerort == "Lavey" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 2360 if bürgerort == "Bémont" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 9001 if bürgerort == "Tablat" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1083 if bürgerort == "Mézières" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1731 if bürgerort == "Ependes" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 3953 if bürgerort == "Loèche" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 3267 if bürgerort == "Seedorf" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 2088 if bürgerort == "Cressier" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 9533 if bürgerort == "Kirchberg SG" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 5630 if bürgerort == "Muri" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 5610 if bürgerort == "Wohlen" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 8890 if bürgerort == "Flums-Grossberg" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 5085 if bürgerort == "Sulz" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1417 if bürgerort == "Essertines" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 6196 if bürgerort == "Marbach" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 3282 if bürgerort == "Bargen" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1734 if bürgerort == "Tinterin" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 3615 if bürgerort == "Bucholterberg" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1031 if bürgerort == "Mex" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1844 if bürgerort == "Villeneuve" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 5610 if bürgerort == "Wohlen AG" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1341 if bürgerort == "Chenit" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1864 if bürgerort == "Ormond-Dessus" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1862 if bürgerort == "Ormond-Dessous" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1611 if bürgerort == "Crêt" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 7418 if bürgerort == "Feldis" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 4332 if bürgerort == "Stein" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1321 if bürgerort == "Arnex" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1625 if bürgerort == "Sales" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1741 if bürgerort == "Cottens" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1482 if bürgerort == "Cugy" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1754 if bürgerort == "Avry" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1969 if bürgerort == "Saint-Martin" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 3538 if bürgerort == "Rothenbach" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 3960 if bürgerort == "Granges" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_1 = 1682 if bürgerort == "Dompierre" & B_PLZ4 == . & nationalität == ""
gen hPLZ4_2 = .
replace hPLZ4_2 = 1731 if bürgerort == "Sales" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 8733 if bürgerort == "Eschenbach" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1277 if bürgerort == "Arnex" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 3033 if bürgerort == "Wohlen" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1186 if bürgerort == "Essertines" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1053 if bürgerort == "Cugy" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1609 if bürgerort == "Saint-Martin" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 6284 if bürgerort == "Sulz" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1890 if bürgerort == "Mex" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 6462 if bürgerort == "Seedorf" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1563 if bürgerort == "Dompierre" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 8233 if bürgerort == "Bargen" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1434 if bürgerort == "Ependes" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1684 if bürgerort == "Mézières" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 5746 if bürgerort == "Walterswil" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1116 if bürgerort == "Cottens" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 9063 if bürgerort == "Stein" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1614 if bürgerort == "Granges" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 3074 if bürgerort == "Muri" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 9437 if bürgerort == "Marbach" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1785 if bürgerort == "Cressier" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_2 = 1527 if bürgerort == "Villeneuve" & B_PLZ4 == . & nationalität == ""
gen hPLZ4_3 = .
replace hPLZ4_3 = 2540 if bürgerort == "Granges" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_3 = 1818 if bürgerort == "Sales" & B_PLZ4 == . & nationalität == ""
replace hPLZ4_3 = 9655 if bürgerort == "Stein" & B_PLZ4 == . & nationalität == ""

gen nat_ch = 1 if nationalität == "Schweiz" // define nationality indicator Schweiz
replace nat_ch = 1 if hPLZ4_1 != . // correct nationality indicator (swiss), if nationalität != "" (but Bürgerort is in CH --> as above)

gen help1=_n
reshape long hPLZ4_, i(help1) j(obs)
rename hPLZ4_ hPLZ4
drop if hPLZ4 == . & obs > 1
replace B_PLZ4 = hPLZ4 if B_PLZ4 == . & hPLZ4 != . 
rename obs Bürgerort_PLZopt2
gen Bürgerort_PLZopt = Bürgerort_PLZopt1
replace Bürgerort_PLZopt = Bürgerort_PLZopt2 if Bürgerort_PLZopt1 == 1 // there is no overlap of _PLZopt1 and 2
drop hPLZ4 help1 Bürgerort_PLZopt1 Bürgerort_PLZopt2

order duns personenid anrede titel vorname nachname W_PLZ4 wohnort wohnland foreign B_PLZ4 bürgerort nationalität ///
geburtstag gremium funktion unterschrift eintrittdatum austrittdatum austrittindikator kapitalanteil PLZWmiss PLZBmiss ///
wohnort_PLZopt Bürgerort_PLZopt plzwohnort plzbürgerort
save "$path\02_Processed_data\11_Directors_1994_2018/Personendaten_corr.dta", replace

* Merge Geo-Info on Bisnode PLZ Wohnort & PLZ Bürgerort
use "$path\02_Processed_data\11_Directors_1994_2018/Personendaten_corr.dta", clear
rename W_PLZ4 PLZ4
merge m:1 PLZ4 using "$path\02_Processed_data\08_Municipalities\uniquePLZ_(Avg)GdeNr2018-Geo.dta" // merge on Wohnort PLZ   GdeNr-Geo_GdeNr-PLZ_1931-2018.dta
drop if _merge == 2
rename PLZ4 W_PLZ4
rename GdeNr_E_CNTR E_CNTR_w
rename GdeNr_N_CNTR N_CNTR_w
rename _merge _mergePLZ_w
rename gdenr_2018 gdenr_2018_w
label var W_PLZ4 "PLZ Wohnort"
label var E_CNTR_w "E-Coordinate Wohnort"
label var N_CNTR_w "N-Coordinate Wohnort"
label var _mergePLZ_w "_merge PLZ Wohnort"
label var gdenr_2018_w "GdeNr. 2018 Wohnort"

rename B_PLZ4 PLZ4
sort PLZ4
merge m:1 PLZ4 using "$path\02_Processed_data\08_Municipalities\uniquePLZ_(Avg)GdeNr2018-Geo.dta" // merge on Bürgerort PLZ
drop if _merge == 2
rename PLZ4 B_PLZ4
rename GdeNr_E_CNTR E_CNTR_b
rename GdeNr_N_CNTR N_CNTR_b
rename _merge _mergePLZ_b
rename gdenr_2018 gdenr_2018_b
label var B_PLZ4 "PLZ Buergerort"
label var E_CNTR_b "E-Coordinate Buergerort"
label var N_CNTR_b "N-Coordinate Buergerort"
label var _mergePLZ_b "_merge PLZ Buergerort"
label var gdenr_2018_b "GdeNr. 2018 Buergerort"

gen avggeo = 0
replace avggeo = 1 if countPLZ4>1
label var avggeo "Average Gemeinde-Centroid for all GdeNr pro PLZ"
drop countPLZ4

label var PLZWmiss "Initially missing PLZ Wohnort"
label var PLZBmiss "Initially missing PLZ Bürgerort"
label var wohnort_PLZopt "missing PLZ: Wohnort-PLZ Options"
label var Bürgerort_PLZopt "missing PLZ: Bürgerort-PLZ Options"

sort gdenr_2018_w
rename gdenr_2018_w gdenr_2018

preserve
use "$path\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018.dta", clear
duplicates drop gdenr_2018, force
sort gdenr_2018
save "$path\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018_temp.dta", replace
restore

merge m:1 gdenr_2018 using "$path\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018_temp.dta"
drop if _merge==2
rename _merge _mergelanguage_w
rename gdenr_2018 gdenr_2018_w
rename language language_w
drop gem_name gem_cd
label var language_w "majority language in Wohnort"

sort gdenr_2018_b
rename gdenr_2018_b gdenr_2018
merge m:1 gdenr_2018 using "$path\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018_temp.dta"
drop if _merge==2
rename _merge _mergelanguage_b
rename gdenr_2018 gdenr_2018_b
rename language language_b
drop gem_name gem_cd
label var language_b "majority language in Bürgerort"

drop _mergelanguage_w _mergelanguage_b

sort personenid
save "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Person_Geo.dta", replace
erase "$path\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018_temp.dta"



/*** Merge Person data with Firmen data ***

use "$path\02_Processed_data\11_Directors_1994_2018\Firmendaten.dta", clear
sort duns
duplicates report duns
save "$path\02_Processed_data\11_Directors_1994_2018\Firmendaten_tmp.dta", replace

use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Person_Geo.dta", clear
sort duns
merge m:1 duns using "$path\02_Processed_data\11_Directors_1994_2018\Firmendaten_tmp.dta"

/*
    Result                           # of obs.
    -----------------------------------------
    not matched                            12
        from master                        12  (_merge==1)
        from using                          0  (_merge==2)

    matched                         5,898,177  (_merge==3)
    -----------------------------------------

br personenid anrede titel vorname nachname nationalität funktion duns if _merge==1
personenid	anrede	titel	vorname	nachname	nationalität	funktion	duns
3544532	männlich		Bernard	Terrier		Leiter der Zweigniederlassung	480223555
3544533	männlich		Ralph	Bobbiá		Leiter der Zweigniederlassung	480223558
3544534	weiblich		Iris	Klausmann		Leiter der Zweigniederlassung	480223559
3544535	weiblich		Nadine	Meyer		Leiter der Zweigniederlassung	480223560
3546801	weiblich		Erika	Peter		Leiter der Zweigniederlassung	480224704
3546808	männlich		Bruno	Kuster		Betriebsleiter/in	480224714
1914931	männlich		Werner	Muster	Schweiz	Präsident/in	481295843
1914933	männlich		Hans	Muster	Schweiz	Mitglied	481295843
1915058	männlich	Dipl. Ing.	Vaclav	Sobriecz	Polen	Direktor/in	481295843
1914932	männlich		Peter	Muster	Schweiz	Mitglied	481295843
2008142	weiblich		Sonja	Kvas	Deutschland	Inhaber/in	481996440
3038797	männlich		Federico	Nannini	Italien	Prokurist/in	488252610
other 1013 cases are all missing obseravations.
*/

drop _merge
sort personenid duns
save "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Person-Firmen_Geo.dta",replace

cap erase "$path\02_Processed_data\11_Directors_1994_2018/Personendaten_corr.dta"
cap erase "$path\02_Processed_data\11_Directors_1994_2018/Firmendaten_tmp.dta"
*/
