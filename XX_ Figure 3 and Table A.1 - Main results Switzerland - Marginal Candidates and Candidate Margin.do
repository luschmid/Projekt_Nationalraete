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

* 5. Identify marginal candidates

* Note: We define a marginal candidate as either the candidate with the lowest
*		number of votes among elected candidate or candidate with the highest
*		number of votes among not elected candidate in those parties with at
*		least one seat achieved and at least one candidate who is not elected. 

* a) Marginal candidates in terms of candidate votes

bysort year canton list: egen rank_elected=rank(votes) if elected==1, track
bysort year canton list: egen mrank_not_elected=rank(-votes) if elected==0, track
bysort year canton list: egen no_seats_list=sum(elected) 
bysort year canton list: egen elected_min=min(elected) 

gen marginal_candidate=0
replace marginal_candidate=1 if (rank_elected==1 | mrank_not_elected==1) & (no_seats_list>0 & elected_min==0)

bysort year canton list: egen sum_marginal_candidates=sum(marginal_candidate)
tab sum_marginal_candidates

sort year canton list votemargin elected
br list year canton  first* votes elected *rank_* marginal_candidate if sum_marginal_candidates==3
* Result: All fine. Two situations where marginal candidates have the same number of votes.

* b) Marginal candidates in terms of vote margin
* Note: Multiple marginal candidates are possible if the party margin is binding.

bysort year canton list: egen rank_elected_vm=rank(votemargin) if elected==1, track
bysort year canton list: egen mrank_not_elected_vm=rank(-votemargin) if elected==0, track

gen marginal_candidate_vm=0
replace marginal_candidate_vm=1 if (rank_elected_vm==1 | mrank_not_elected_vm==1) & (no_seats_list>0 & elected_min==0)

* 6. Calculate candidate margin

bysort year canton list: egen temp1=min(votes) if elected==1
bysort year canton list: egen temp2=max(votes) if elected==0
bysort year canton list: egen votes_min=max(temp1) 
bysort year canton list: egen votes_max=max(temp2) 
drop temp1 temp2 

gen candidate_margin=.
replace candidate_margin=votes-votes_max if elected==1
replace candidate_margin=votes-votes_min if elected==0

bysort elected: sum candidate_margin
corr votemargin candidate_margin

gen candidate_margin_rel=candidate_margin/eligible_cant

* 7. RD Estimation for candidate margin

rdd_estimation candidate_margin_rel "ch_all"


* 8. RD Estimation for marginal candidates

* a) Marginal candidates in terms of votes

preserve
keep if marginal_candidate==1
rdd_estimation votemargin_rel "ch_marg"
restore

* b) Marginal candidates in terms of vote margin

preserve
keep if marginal_candidate_vm==1
rdd_estimation votemargin_rel "ch_marg_vm"
restore
