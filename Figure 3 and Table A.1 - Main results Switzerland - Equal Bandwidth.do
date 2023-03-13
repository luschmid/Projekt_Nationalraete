* 1. Set directories 

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

set more off

* 2. Define variables (votemargin, covariates) and options (kernel and cluster)

global outcomes elected_F1 participation_F1 

global covs votemargin_rel_L1 age sex year no_seats ///
	canton_1 canton_2 canton_3 canton_4 canton_5 canton_6 canton_7 canton_8 ///
	canton_9 canton_10 canton_11 canton_12 canton_13 canton_14 canton_15 ///
	canton_16 canton_17 canton_18 canton_19 canton_20 canton_21 canton_22 ///
	canton_23 canton_24 canton_25 canton_26 party_sp party_cvp party_fdp ///
	party_svp

global vars_all $outcomes $covs
	
global kernel triangular // Options: triangular (default), epanechnikov, uniform
global secl cluster // Options: cluster nncluster

* 3. Read in data, generate running variable, election in next election,
* 	 party variables, and cantonal dummy variable

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear

g votemargin_rel_L1 = L.votemargin_rel

gen party_sp= listname_bfs=="SP/PS" if year>=1971 
gen party_cvp= listname_bfs=="CVP/PDC" if year>=1971 
gen party_fdp= listname_bfs=="FDP/PLR (PRD)" if year>=1971 
gen party_svp= listname_bfs=="SVP/UDC" if year>=1971 
replace party_svp=1 if listname_bfs=="BDP/PBD" & year>=1971 

tab canton, gen(canton_)

* 4. Variable labeling

label var votemargin_rel_L1 "Vote margin in t-1"
label var votemargin_rel "Vote margin (eligibles)"

label var canton_1 "Aargau"
label var canton_2 "Appenzell Innerrhoden"
label var canton_3 "Appenzell Ausserrhoden"
label var canton_4 "Bern"
label var canton_5 "Basel Landschaft"
label var canton_6 "Basel Stadt"
label var canton_7 "Fribourg"
label var canton_8 "Geneva"
label var canton_9 "Glarus"
label var canton_10 "Graubünden"
label var canton_11 "Jura"
label var canton_12 "Lucerne"
label var canton_13 "Neuchâtel"
label var canton_14 "Nidwalden"
label var canton_15 "Obwalden"
label var canton_16 "St. Gallen"
label var canton_17 "Schaffhausen"
label var canton_18 "Solothurn"
label var canton_19 "Schwyz"
label var canton_20 "Thurgau"
label var canton_21 "Ticino"
label var canton_22 "Uri"
label var canton_23 "Vaud"
label var canton_24 "Valais"
label var canton_25 "Zug"
label var canton_26 "Zurich"

label var party_sp "Social democrats (SP)"
label var party_cvp "Christian democratats (CVP)"
label var party_fdp "Liberals (FDP)"
label var party_svp "National convervatives (SVP)"

* 5. RDD estimation

cap program drop rdd_estimation
do "$path\03_Code\13_Running_variable\rdd_estimation.do"

rdd_estimation votemargin_rel "ch_all"
