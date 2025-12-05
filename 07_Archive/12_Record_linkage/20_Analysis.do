* 1. Set directories 

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

set more off

* 2. Define model name, variables (votemargin, covariates), 
* 	 and options (kernel and cluster)

global model_name "Bisnode_only_1994_2017_Gen7_optimal_corrected" 
*global model_name "Bisnode_only_1994_2017_Gen7_optimal_high_capital"

global outcomes mandates_lead4 

global covs votemargin_rel_L1 age sex year no_seats ///
	canton_1 canton_2 canton_3 canton_4 canton_5 canton_6 canton_7 canton_8 ///
	canton_9 canton_10 canton_11 canton_12 canton_13 canton_14 canton_15 ///
	canton_16 canton_17 canton_18 canton_19 canton_20 canton_21 canton_22 ///
	canton_23 canton_24 canton_25 canton_26 party_sp party_cvp party_fdp ///
	party_svp

global vars_all $outcomes $covs
	
global kernel triangular // Options: triangular (default), epanechnikov, uniform
global secl cluster // Options: cluster nncluster

* 3. Read in record linkage data and save it 

insheet using "$path\02_Processed_data\12_Record_linkage\03_Data_analysis\results_wide_$model_name.csv", ///
	delimiter(";") clear
rename id ID
keep ID year mandate_lag_4 mandate_lag_3 mandate_lag_2 mandate_lag_1 mandate_0 mandate_lead_1 mandate_lead_2 mandate_lead_3 mandate_lead_4
save "$path\02_Processed_data\12_Record_linkage\03_Data_analysis\results_wide_$model_name.dta", replace


* 4. Read in nr data and merge record linkage data, generate running variable, 
* election in next election, party variables, and cantonal dummy variable

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear

sort ID year
bysort ID: gen time_period=_n

merge 1:1 ID year using "$path\02_Processed_data\12_Record_linkage\03_Data_analysis\results_wide_$model_name.dta", ///
 gen(merge_model)
 
bysort ID: egen year_min=min(year)

keep if year>=1995

egen mandate_pre=rowmean(mandate_lag_1 mandate_lag_2 mandate_lag_3 mandate_lag_4)
egen mandate_post=rowmean(mandate_lead_1 mandate_lead_2 mandate_lead_3 mandate_lead_4)
*br ID year mandate_lead* mandates_lead4

forv i=1(1)4 {
g dummy_mandate_lag_`i'=.
replace dummy_mandate_lag_`i'=1 if mandate_lag_`i'>0  & mandate_lag_`i'!=.
replace dummy_mandate_lag_`i'=0 if mandate_lag_`i'==0 
g dummy_mandate_lead_`i'=.
replace dummy_mandate_lead_`i'=1 if mandate_lead_`i'>0 & mandate_lag_`i'!=.
replace dummy_mandate_lead_`i'=0 if mandate_lead_`i'==0	
}

g dummy_mandate_0=.
replace dummy_mandate_0=1 if mandate_0>0
replace dummy_mandate_0=0 if mandate_0==0

order ID year canton name firstname birthyear age sex job votes elected ///
dummy_mandate_lag_4 dummy_mandate_lag_3 dummy_mandate_lag_2 dummy_mandate_lag_1 ///
dummy_mandate_0 dummy_mandate_lead_1 dummy_mandate_lead_2 dummy_mandate_lead_3 ///
dummy_mandate_lead_4 

g dummy_mandate_post=.
replace dummy_mandate_post=1 if mandate_post>0
replace dummy_mandate_post=0 if mandate_post==0

g dummy_mandate_pre=.
replace dummy_mandate_pre=1 if mandate_pre>0
replace dummy_mandate_pre=0 if mandate_pre==0

egen ID_num=group(ID)
sort ID_num ID_time 
xtset ID_num ID_time  
g votemargin_rel_L1 = L.votemargin_rel

gen party_sp= listname_bfs=="SP/PS" if year>=1971 
gen party_cvp= listname_bfs=="CVP/PDC" if year>=1971 
gen party_fdp= listname_bfs=="FDP/PLR (PRD)" if year>=1971 
gen party_svp= listname_bfs=="SVP/UDC" if year>=1971 
replace party_svp=1 if listname_bfs=="BDP/PBD" & year>=1971 

tab canton, gen(canton_)

* 5. Variable labeling

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


* 6. Dynamic effect

local counter=-4

foreach var of varlist dummy_mandate_lag_4-dummy_mandate_lead_4 {
rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year>=1995 
est store p1
if "`var'" == "dummy_mandate_lag_4" {
	esttab p1  using "$path\05_Texts_and_presns\03_Political_Rents\tables\\dynamic_effect_ext.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`counter'") varwidth(21) star(* 0.10 ** 0.05 *** 0.01) fragment ///
	sca("h_l Bandwidth" "N_h_l Left obs" "N_h_r Right obs")  ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes replace
	}
else {
	esttab p1 using "$path\05_Texts_and_presns\03_Political_Rents\tables\\dynamic_effect_ext.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`counter'") varwidth(21) fragment star(* 0.10 ** 0.05 *** 0.01) ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes append
	}
local ++counter	
}

local counter=-4

foreach var of varlist mandate_lag_4-mandate_lead_4 {
rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year>=1995 
est store p1
if "`var'" == "mandate_lag_4" {
	esttab p1 using "$path\05_Texts_and_presns\03_Political_Rents\tables\\dynamic_effect_int.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`counter'") varwidth(21) star(* 0.10 ** 0.05 *** 0.01) fragment ///
	sca("h_l Bandwidth" "N_h_l Left obs" "N_h_r Right obs")  ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes replace
	}
else {
	esttab p1 using "$path\05_Texts_and_presns\03_Political_Rents\tables\\dynamic_effect_int.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`counter'") varwidth(21) fragment star(* 0.10 ** 0.05 *** 0.01) ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes append
	}
local ++counter	
}

/*esttab p1  using "$path\05_Texts_and_presns\03_Political_Rents\tables\\dynamic_effect.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`label'") varwidth(21) star(* 0.10 ** 0.05 *** 0.01) fragment ///
	sca("h_l Bandwidth" "N_h_l Left obs" "N_h_r Right obs")  ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes replace
*/
	

*** 

* 7. Mean effect

* a) Lead variables 

rdrobust dummy_mandate_post votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
	
rdrobust dummy_mandate_post  votemargin_rel, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
	
rdrobust mandate_post  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
	
rdrobust mandate_post  votemargin_rel, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
	
* b) Lagged variables 
	
rdrobust dummy_mandate_pre  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
	
rdrobust dummy_mandate_pre  votemargin_rel, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
	
rdrobust mandate_pre  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
	
rdrobust mandate_pre  votemargin_rel, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
	
* 8. Mean effect for first candidacy

rdrobust dummy_mandate_post votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year==year_min

rdrobust dummy_mandate_post  votemargin_rel, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year==year_min
	
rdrobust mandate_post  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year==year_min
	
rdrobust mandate_post  votemargin_rel, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year==year_min
	
* b) Lagged variables 
	
rdrobust dummy_mandate_pre  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year==year_min
	
rdrobust dummy_mandate_pre  votemargin_rel, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year==year_min
	
rdrobust mandate_pre  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year==year_min
	
rdrobust mandate_pre  votemargin_rel, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year==year_min
	
* 9. Mean effect by party

* a) extensive margin

rdrobust dummy_mandate_post  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if  party_fdp==1	
	 
rdrobust dummy_mandate_post  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if  party_sp==1	
	
rdrobust dummy_mandate_post  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if party_cvp==1	
	
rdrobust dummy_mandate_post  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if party_svp==1		

* b) intensive margin
	
rdrobust mandate_post votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if  party_fdp==1	
	 
rdrobust mandate_post  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year>=1995 & party_sp==1	
	
rdrobust mandate_post  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year>=1995 & party_cvp==1	
	
rdrobust mandate_post  votemargin_rel, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1), if year>=1995 & party_svp==1		

* 10. Binscatter
* Install package: net install binsreg, from(https://raw.githubusercontent.com/nppackages/binsreg/master/stata) replace


binsreg dummy_mandates_lead4 votemargin_rel // w i.t, at(`evalcovar´)

binsreg mandates_lead4 votemargin_rel if inrange(votemargin_rel,-0.059,0.059)

rdplot dummy_mandates_lead4 votemargin_rel ///
	if -0.059 <= votemargin_rel & votemargin_rel <= 0.059, ///
	p(1) kernel($kernel) binselect(qs) 
	
rdplot dummy_mandates_lag4 votemargin_rel ///
	if -0.059 <= votemargin_rel & votemargin_rel <= 0.059, ///
	p(1) kernel($kernel) binselect(qs) 	
	
binsreg mandates_lead4 votemargin_rel , nbins(10) polyreg(2)


/*
* 6. RDD estimation

cap program drop rdd_estimation
do "$path\03_Code\13_Running_variable\rdd_estimation.do"

rdd_estimation votemargin_rel "ch_all"
