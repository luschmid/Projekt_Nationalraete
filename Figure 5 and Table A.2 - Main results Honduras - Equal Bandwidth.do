* 1. Set directories 

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

set more off

* 2. Define variables (votemargin, covariates) and options (kernel and cluster)

global margin votemargin_rel
local margin $margin // only for file saving

global covs elected_F1 participation_F1 votemargin_L1 sex year no_seats ///
	Department_1 Department_2 Department_3 Department_4 Department_5 Department_6 Department_7 Department_8 ///
	Department_9 Department_10 Department_11 Department_12 Department_13 Department_14 Department_15 ///
	Department_16 Department_17 Department_18 Party_pn Party_pl Party_libre

global vars_all $outcomes $covs
	
global kernel triangular // Options: triangular (default), epanechnikov, uniform
global secl cluster // Options: cluster nncluster

* 3. Read in data, generate running variable and variable renaming

use "$path\02_Processed_data\15_Elections_Honduras\elections_hn_final.dta", clear

drop votemargin_rel
gen votemargin_rel= Votemargin/total_electoral
capture drop votemargin_L1
xtset ID_pers ID_time
g votemargin_L1 = L.$margin
label var votemargin_L1 "Vote margin in t-1"

rename Elected elected
rename Elected_F1 elected_F1 
rename Elected_L1 elected_L1
rename Year year
gen sex=.
replace sex=1 if Sex=="H"
replace sex=0 if Sex=="M"
drop Sex

gen Party_pn= Party=="PARTIDO NACIONAL DE HONDURAS" 
gen Party_pl= Party=="PARTIDO LIBERAL DE HONDURAS" 
gen Party_libre= Party=="PARTIDO LIBERTAD Y REFUNDACION" 

sort Department
tab Department, gen(Department_)

* 4. Variable labeling

label var votemargin_rel "Vote margin (eligibles)"
label var elected_F1 "Elected in t+1"
label var participation_F1 "Participation in t+1"
label var no_seats "Number of seats"
label var year "Year"
label var sex "Sex"

label var Department_1 "ATLANTIDA"
label var Department_2 "CHOLUTECA"
label var Department_3 "COLON"
label var Department_4 "COMAYAGUA"
label var Department_5 "COPAN"
label var Department_6 "CORTES"
label var Department_7 "EL PARAISO"
label var Department_8 "FRANCISCO MORAZAN"
label var Department_9 "GRACIAS A DIOS"
label var Department_10 "INTIBUCA"
label var Department_11 "ISLAS DE LA BAHIA "
label var Department_12 "LA PAZ"
label var Department_13 "LEMPIRA"
label var Department_14 "OCOTEPEQUE"
label var Department_15 "OLANCHO"
label var Department_16 "SANTA BARBARA"
label var Department_17 "VALLE"
label var Department_18 "YORO"

label var Party_pn "National conservatives (PN)"
label var Party_pl "Social liberals (PL)"
label var Party_libre "Social democrats (Libre)"

* 5. Manipulation test

rddensity $margin, kernel($kernel) plot


* 5. RDD estimation

cap program drop rdd_estimation
do "$path\03_Code\13_Running_variable\rdd_estimation.do"

rdd_estimation votemargin_rel "hn_all"

/*

* 6. RD Estimation

* (a) Summary statistics and bandwidth definition

sum $covs

rdrobust elected_F1 $margin, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
global bw_opt_p1 =e(h_l)
global bbw_opt_p1 =e(b_l)
di $bw_opt_p1
di $bbw_opt_p1

rdrobust elected_F1 $margin, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
global bw_opt_p2 =e(h_l)
global bbw_opt_p2 =e(b_l)
di $bw_opt_p2
di $bbw_opt_p2

* (b) (RD plot)
/*
foreach z of varlist $covs {

local label: var label `z' 

cap rdplot `z' $margin if elected_F1 < ., ///
	kernel($kernel) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off) ///
	ti(`label', position(11) orientation(horizontal) si(huge)))
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig1_hn.pdf", ///
	as(pdf) replace

cap rdplot `z' $margin ///
	if -$bw_opt_p1 <= $margin & $margin <= $bw_opt_p1 & elected_F1 < ., ///
	p(1) kernel($kernel) h($bw_opt_p1 $bw_opt_p1) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off))
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig2_hn.pdf", ///
	as(pdf) replace
	
cap rdplot `z' $margin ///
	if -$bw_opt_p2 <= $margin & $margin <= $bw_opt_p2 & elected_F1 < ., ///
	p(2) kernel($kernel) h($bw_opt_p2 $bw_opt_p2) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off))
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig4_hn.pdf", ///
	as(pdf) replace
} */

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
sum `z' if elected_F1 < .

rdrobust  `z' $margin if elected_F1 < ., p(1) h($bw_opt_p1) rho(1) ///
	kernel($kernel) vce($secl ID_pers) all
est store p1
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1" local ++cp1_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_10
local bw_opt_p1_div2=$bw_opt_p1/2
local bbw_opt_p1_div2=$bbw_opt_p1/2
rdrobust  `z' $margin if elected_F1 < ., p(1) h(`bw_opt_p1_div2') rho(1) ///
	kernel($kernel) vce($secl ID_pers) all 
est store p1_div2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_div2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_div2_10
rdrobust  `z' $margin if elected_F1 < ., p(2) h($bw_opt_p2) rho(1) ///
	kernel($kernel) vce($secl ID_pers) all
est store p2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_10
local bw_opt_p2_div2=$bw_opt_p2/2
local bbw_opt_p2_div2=$bbw_opt_p2/2
rdrobust  `z' $margin if elected_F1 < ., p(2) h(`bw_opt_p2_div2') rho(1)  ///
	kernel($kernel) vce($secl ID_pers) all
est store p2_div2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_div2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_div2_10

if "`z'" == "elected_F1" {
	noisily esttab p1 p1_div2 p2 p2_div2 using "$path\05_Texts_and_presns\01_Running_variable\tables\\`margin'_balance_hn.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`label'") varwidth(21) ///
	sca("h_l Bandwidth" "N_h_l Left obs" "N_h_r Right obs") ///
	star(* 0.10 ** 0.05 *** 0.01) fragment ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes replace
}

else{
	noisily esttab p1 p1_div2 p2 p2_div2 using "$path\05_Texts_and_presns\01_Running_variable\tables\\`margin'_balance_hn.tex", ///
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
esttab matrix(rrates) using "$path\05_Texts_and_presns\01_Running_variable\tables\\`margin'_balance_hn.tex",  ///
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
cap rdrobust  `z' $margin if elected_F1 < ., p(1) h(`bw_opt_p1_tmp') b(`bbw_opt_p1_tmp') ///
	kernel($kernel) vce($secl ID_pers) all
if _rc==0{
mat pointestimate[1,`col']=e(tau_bc)
mat CI[1,`col'] = e(ci_l_rb)
mat CI[2,`col'] = e(ci_r_rb)
}
else{
}
local ++col
}
mat colnames pointestimate = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
mat colnames CI = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
coefplot matrix(pointestimate), ci(CI) vertical nolabel yline(0) ///
	mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig3_hn.pdf", ///
	as(pdf) replace
clear matrix

mat pointestimate = J(1,15,.)
mat CI = J(2,15,.)
local col=1
forv i = 0.1(0.1)1.6 {
local bw_opt_p2_tmp=$bw_opt_p2*`i'
local bbw_opt_p2_tmp=$bbw_opt_p2*`i'
cap rdrobust  `z' $margin if elected_F1 < ., p(2) h(`bw_opt_p2_tmp') b(`bbw_opt_p2_tmp') ///
	kernel($kernel) vce($secl ID_pers) all
if _rc==0{
mat pointestimate[1,`col']=e(tau_bc)
mat CI[1,`col'] = e(ci_l_rb)
mat CI[2,`col'] = e(ci_r_rb)
}
else{
}
local ++col
}
mat colnames pointestimate = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
mat colnames CI = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
coefplot matrix(pointestimate), ci(CI) vertical nolabel yline(0) ///
	mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig5_hn.pdf", ///
	as(pdf) replace
}
}
