* 1. Set directories 

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

set more off

* 2. Define variables (votemargin, covariates) and options (kernel and cluster)

global margin votemargin_rel
local margin $margin // only for file saving

global covs elected_F1 participation_F1 votemargin_L1 age sex year no_seats ///
	canton_1 canton_2 canton_3 canton_4 canton_5 canton_6 canton_7 canton_8 ///
	canton_9 canton_10 canton_11 canton_12 canton_13 canton_14 canton_15 ///
	canton_16 canton_17 canton_18 canton_19 canton_20 canton_21 canton_22 ///
	canton_23 canton_24 canton_25 canton_26 party_sp party_cvp party_fdp ///
	party_svp

global kernel triangular // Options: triangular (default), epanechnikov, uniform
global secl cluster // Options: cluster nncluster

* 3. Read in data, generate running variable, election in next election,
* 	 party variables, and cantonal dummy variable

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear

capture drop votemargin_rel // fix this!!
gen votemargin_rel= votemargin/eligible_cant
drop votemargin_L1
xtset ID_pers ID_time
g votemargin_L1 = L.$margin

gen party_sp= listname_bfs=="SP/PS" if year>=1971 
gen party_cvp= listname_bfs=="CVP/PDC" if year>=1971 
gen party_fdp= listname_bfs=="FDP/PLR (PRD)" if year>=1971 
gen party_svp= listname_bfs=="SVP/UDC" if year>=1971 
replace party_svp=1 if listname_bfs=="BDP/PBD" & year>=1971 

tab canton, gen(canton_)

* 4. Variable labeling

label var votemargin_L1 "Vote margin in t-1"
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

* 5. RD Estimation

* (a) Summary statistics and bandwidth definition

sum $covs

rdrobust elected_F1 $margin, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all 
global bw_opt_p1 =e(h_l)
global bbw_opt_p1 =e(b_l)
di $bw_opt_p1
di $bbw_opt_p1

rdrobust elected_F1 $margin, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all
global bw_opt_p2 =e(h_l)
global bbw_opt_p2 =e(b_l)
di $bw_opt_p2
di $bbw_opt_p2

* (b) RD plot
/*
foreach z of varlist $covs {

local label: var label `z' 

cap rdplot `z' $margin if elected_F1 < ., ///
	kernel($kernel) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off) ///
	ti(`label', position(11) orientation(horizontal) si(huge))) ///
	plotregion(lcolor(black))
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig1.pdf", ///
	as(pdf) replace

cap rdplot `z' $margin ///
	if -$bw_opt_p1 <= $margin & $margin <= $bw_opt_p1 & elected_F1 < ., ///
	p(1) kernel($kernel) h($bw_opt_p1 $bw_opt_p1) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig2.pdf", ///
	as(pdf) replace
	
cap rdplot `z' $margin ///
	if -$bw_opt_p2 <= $margin & $margin <= $bw_opt_p2 & elected_F1 < ., ///
	p(2) kernel($kernel) h($bw_opt_p2 $bw_opt_p2) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off))
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig4.pdf", ///
	as(pdf) replace
}*/

* (c) Balance table

local cp1_5 = 0
local cp1_10 = 0
local cp1_div2_5 = 0
local cp1_div2_10 = 0
local cp2_5 = 0
local cp2_10 = 0
local cp2_div2_5 = 0
local cp2_div2_10 = 0

foreach z of varlist $covs {
qui {
local label: var label `z' 
rdrobust  `z' $margin if elected_F1 < ., p(1) h($bw_opt_p1) b($bbw_opt_p1) ///
	kernel($kernel) vce($secl ID_pers) all
est store p1
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1" local ++cp1_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_10
local bw_opt_p1_div2=$bw_opt_p1/2
local bbw_opt_p1_div2=$bbw_opt_p1/2
rdrobust  `z' $margin if elected_F1 < ., p(1) h(`bw_opt_p1_div2') b(`bbw_opt_p1_div2') ///
	kernel($kernel) vce($secl ID_pers) all 
est store p1_div2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_div2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_div2_10
rdrobust  `z' $margin if elected_F1 < ., p(2) h($bw_opt_p2) b($bbw_opt_p2) ///
	kernel($kernel) vce($secl ID_pers) all
est store p2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_10
local bw_opt_p2_div2=$bw_opt_p2/2
local bbw_opt_p2_div2=$bbw_opt_p2/2
rdrobust  `z' $margin if elected_F1 < ., p(2) h(`bw_opt_p2_div2') b(`bbw_opt_p2_div2') ///
	kernel($kernel) vce($secl ID_pers) all
est store p2_div2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_div2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_div2_10

if "`z'" == "elected_F1" {
	esttab p1 p1_div2 p2 p2_div2 using "$path\05_Texts_and_presns\01_Running_variable\tables\\`margin'_balance.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`label'") varwidth(21) star(* 0.10 ** 0.05 *** 0.01) fragment ///
	sca("h_l Bandwidth" "N_h_l Left obs" "N_h_r Right obs")  ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes replace

	}

else {
	esttab p1 p1_div2 p2 p2_div2 using "$path\05_Texts_and_presns\01_Running_variable\tables\\`margin'_balance.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`label'") varwidth(21) fragment star(* 0.10 ** 0.05 *** 0.01) ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes append
	}
}
}

mat rrates = J(2,4,.)
mat rrates[1,1] =`cp1_5'
mat rrates[2,1] =`cp1_10'
mat rrates[1,2] =`cp1_div2_5'
mat rrates[2,2] =`cp1_div2_10'
mat rrates[1,3] =`cp2_5'
mat rrates[2,3] =`cp2_10'
mat rrates[1,4] =`cp2_div2_5'
mat rrates[2,4] =`cp2_div2_10'
mat rownames rrates = "5\%" "10\%"
mat colnames rrates = "\hspace{0.1cm}" "\hspace{0.2cm}" "\hspace{0.3cm}" "\hspace{0.4cm}"
*outtable using "$path\05_Texts_and_presns\01_Running_variable\tables\\`margin'_balance", ///
*	mat(rrates) append nobox
esttab matrix(rrates) using "$path\05_Texts_and_presns\01_Running_variable\tables\\`margin'_balance.tex",  ///
	nomti nonum append fragment

* (d) Bandwidth figures

foreach z of varlist $covs {
qui {
mat pointestimate = J(1,15,.)
mat CI = J(2,15,.)
local col=1
forv i = 0.1(0.1)1.6 {
local bw_opt_p1_tmp=$bw_opt_p1*`i'
local bbw_opt_p1_tmp=$bbw_opt_p1*`i'
rdrobust  `z' $margin if elected_F1 < ., p(1) h(`bw_opt_p1_tmp') b(`bbw_opt_p1_tmp') ///
	kernel($kernel) vce($secl ID_pers) all
mat pointestimate[1,`col']=e(tau_bc)
mat CI[1,`col'] = e(ci_l_rb)
mat CI[2,`col'] = e(ci_r_rb)
local ++col
}
mat colnames pointestimate = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
mat colnames CI = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
coefplot matrix(pointestimate), ci(CI) vertical nolabel yline(0) ///
	mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig3.pdf", ///
	as(pdf) replace
clear matrix

mat pointestimate = J(1,15,.)
mat CI = J(2,15,.)
local col=1
forv i = 0.1(0.1)1.6 {
local bw_opt_p2_tmp=$bw_opt_p2*`i'
local bbw_opt_p2_tmp=$bbw_opt_p2*`i'
rdrobust  `z' $margin if elected_F1 < ., p(2) h(`bw_opt_p2_tmp') b(`bbw_opt_p2_tmp') ///
	kernel($kernel) vce($secl ID_pers) all
mat pointestimate[1,`col']=e(tau_bc)
mat CI[1,`col'] = e(ci_l_rb)
mat CI[2,`col'] = e(ci_r_rb)
local ++col
}
mat colnames pointestimate = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
mat colnames CI = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
coefplot matrix(pointestimate), ci(CI) vertical nolabel yline(0) ///
	mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig5.pdf", ///
	as(pdf) replace
}
}

