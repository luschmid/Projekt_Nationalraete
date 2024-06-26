
*-------------------------------------------
* (A) Set path and install package ietoolkit 
*-------------------------------------------

global path "..."
cd "$path"

ssc install ietoolkit
ieboilstart, version(17.0) adopath("$path/ado", strict)

*----------------------
* (B) Load RDD function
*----------------------

cap program drop rdd_estimation

program rdd_estimation

* (a) Bandwidth definition

rdrobust elected_F1 `1', p(1) bwselect(mserd) kernel(triangular)  ///
	vce(cluster ID_pers) all rho(1) 
global bw_opt_p1 =e(h_l)
di $bw_opt_p1

rdrobust elected_F1 `1', p(2) bwselect(mserd) kernel(triangular)  ///
	vce(cluster ID_pers) all rho(1)
global bw_opt_p2 =e(h_l)
di $bw_opt_p2

* (b) RD plot

local label: var label elected_F1 

cap rdplot elected_F1 `1'  ///
	if -$bw_opt_p1 <= `1' & `1' <= $bw_opt_p1 & elected_F1 < ., ///
	p(1) kernel(triangular) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "./Figures/`2'.pdf", ///
	as(pdf) replace

* (c) Bandwidth figures

qui {
mat pointestimate = J(1,15,.)
mat CI = J(2,15,.)
local col=1
forv i = 0.1(0.1)1.6 {
local bw_opt_p1_tmp=$bw_opt_p1*`i'
rdrobust elected_F1 `1' if elected_F1 < ., p(1) h(`bw_opt_p1_tmp')  ///
	kernel(triangular) vce(cluster ID_pers) all rho(1)
mat pointestimate[1,`col']=e(tau_cl)
mat CI[1,`col'] = e(ci_l_cl)
mat CI[2,`col'] = e(ci_r_cl)
local ++col
}
mat colnames pointestimate = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
mat colnames CI = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
coefplot matrix(pointestimate), ci(CI) vertical nolabel yline(0) ///
	mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)
graph export "./Figures/`3'.pdf", ///
	as(pdf) replace
clear matrix
}

* (d) Balance table

local cp1_5 = 0
local cp1_10 = 0
local cp1_div2_5 = 0
local cp1_div2_10 = 0
local cp2_5 = 0
local cp2_10 = 0
local cp2_div2_5 = 0
local cp2_div2_10 = 0

global len_covs: word count $covs

foreach z of varlist $vars_all {
qui {
local label: var label `z' 
rdrobust `z' `1' if elected_F1 < ., p(1) h($bw_opt_p1) rho(1) ///
	kernel(triangular) vce(cluster ID_pers) all
est store p1
if `e(pv_rb)' < .05 & "`z'" != "elected_F1"   local ++cp1_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1"    local ++cp1_10
local bw_opt_p1_div2=$bw_opt_p1/2
rdrobust `z' `1' if elected_F1 < ., p(1) h(`bw_opt_p1_div2') rho(1) ///
	kernel(triangular) vce(cluster ID_pers) all 
est store p1_div2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1"   local ++cp1_div2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1"    local ++cp1_div2_10
rdrobust `z' `1' if elected_F1 < ., p(2) h($bw_opt_p2) rho(1) ///
	kernel(triangular) vce(cluster ID_pers) all
est store p2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1"   local ++cp2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1"    local ++cp2_10
local bw_opt_p2_div2=$bw_opt_p2/2
rdrobust `z' `1' if elected_F1 < ., p(2) h(`bw_opt_p2_div2') rho(1) ///
	kernel(triangular) vce(cluster ID_pers) all
est store p2_div2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1"   local ++cp2_div2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1"    local ++cp2_div2_10

if "`z'" == "elected_F1" {
	esttab p1 p1_div2 p2 p2_div2 using "./Tables/`4'.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`label'") varwidth(21) nostar fragment ///
	sca("h_l Bandwidth" "N_h_l Left obs" "N_h_r Right obs")  ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes replace

	}

else {
	esttab p1 p1_div2 p2 p2_div2 using "./Tables/`4'.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`label'") varwidth(21) fragment nostar ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes append
	}
}
}

mat rrates = J(2,4,.)
mat rrates[1,1] =`cp1_5'/$len_covs*100
mat rrates[2,1] =`cp1_10'/$len_covs*100
mat rrates[1,2] =`cp1_div2_5'/$len_covs*100
mat rrates[2,2] =`cp1_div2_10'/$len_covs*100
mat rrates[1,3] =`cp2_5'/$len_covs*100
mat rrates[2,3] =`cp2_10'/$len_covs*100
mat rrates[1,4] =`cp2_div2_5'/$len_covs*100
mat rrates[2,4] =`cp2_div2_10'/$len_covs*100
mat rownames rrates = "5\%" "10\%"
mat colnames rrates = "\hspace{0.1cm}" "\hspace{0.2cm}" "\hspace{0.3cm}" "\hspace{0.4cm}"
esttab matrix(rrates) using "./Tables/`4'.tex",  ///
	nomti nonum append fragment

end


*----------------------------
* (C) Figure 1 and Table B.1
*----------------------------

* (i) Define variables (votemargin, covariates)

global covs votemargin_rel_L1 age sex year no_seats ///
	canton_1 canton_2 canton_3 canton_4 canton_5 canton_6 canton_7 canton_8 ///
	canton_9 canton_10 canton_11 canton_12 canton_13 canton_14 canton_15 ///
	canton_16 canton_17 canton_18 canton_19 canton_20 canton_21 canton_22 ///
	canton_23 canton_24 canton_25 canton_26 party_sp party_cvp party_fdp ///
	party_svp

global vars_all elected_F1 $covs
	
* (ii) Read in data and perform RD estimations

use "data_switzerland.dta", clear
xtset ID_pers ID_time
sort ID_pers ID_time
g votemargin_rel_L1 = L.votemargin_rel
gen elected_F1=F1.elected
replace elected_F1=0 if elected_F1==. & year!=2015
label var votemargin_rel "Vote margin"
label var votemargin_rel_L1 "Vote margin in t-1"
label var elected_F1 "Elected in t+1"

rdd_estimation votemargin_rel "Figure1A" "Figure1B" "TableB1"

*----------------------------
* (D) Figure 2 and Table B.2
*----------------------------

* (i) Define variables (votemargin, covariates)

global covs votemargin_rel_L1 sex year no_seats ///
	department_1 department_2 department_3 department_4 department_5 department_6 department_7 department_8 ///
	department_9 department_10 department_11 department_12 department_13 department_14 department_15 ///
	department_16 department_17 department_18 party_pn party_pl party_libre

global vars_all elected_F1 $covs
	
* (ii) Read in data and perform RD estimations

use "data_honduras.dta", clear

xtset ID_pers ID_time
sort ID_pers ID_time
g votemargin_rel_L1 = L.votemargin_rel 
gen elected_F1=F1.elected
replace elected_F1=0 if elected_F1==. & year!=2017
label var votemargin_rel "Vote margin"
label var votemargin_rel_L1 "Vote margin in t-1"
label var elected_F1 "Elected in t+1"

rdd_estimation votemargin_rel "Figure2A" "Figure2B" "TableB2"

*--------------
* (E) Figure 3 
*---------------

* (i) Read in data and perform RD estimations

use "data_norway.dta", clear


xtset ID_pers ID_time
sort ID_pers ID_time
gen elected_F1=F1.elected
replace elected_F1=0 if elected_F1==.
drop if year>=1985
* Note: Electoral system changes in 1985.  
label var votemargin_rel "Vote margin"
label var elected_F1 "Elected in t+1"

rdrobust elected_F1 votemargin_rel, p(1) bwselect(mserd) kernel(triangular)  ///
	vce(cluster ID_pers) all rho(1), if marginal_candidate=="Marginal candidate"
global bw_sample_marginals =e(h_l)

rdrobust elected_F1 votemargin_rel, p(1) bwselect(mserd) kernel(triangular)  ///
	vce(cluster ID_pers) all rho(1)
global bw_full_sample =e(h_l)

* (ii) RD plots

rdplot elected_F1 votemargin_rel ///
	if -$bw_sample_marginals <= votemargin_rel & votemargin_rel<= $bw_sample_marginals ///
	& marginal_candidate=="Marginal candidate" & elected_F1 < ., ///
	p(1) kernel(triangular) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "./Figures/Figure3A.pdf", ///
	as(pdf) replace	
	
rdplot elected_F1 votemargin_rel ///
	if -$bw_full_sample <= votemargin_rel & votemargin_rel<= $bw_full_sample & elected_F1 < ., ///
	p(1) kernel(triangular) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)  ) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "./Figures/Figure3B.pdf", ///
	as(pdf) replace
