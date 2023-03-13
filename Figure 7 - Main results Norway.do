* 1. Set directories 

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

set more off

* 2. Read in data on running variable

import delimited using ///
	"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981.csv", ///
	encoding("utf-8") clear 
save "$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981.dta", ///
	replace
	
* 3. Read in data on from Fiva and Smith (2018), small correction,
*	 focus on sample period, and merge data on running variable

use "$path\01_Raw_data\13_Running_variable\FivaSmithJune2019.dta", clear

* a) Small correction of data

replace pid=22893 if candidatename_orig=="Lars Fagerland" & year==1973
* Note: The original pid refers to the pid of another person. 
*		Fiva and Smith (2018) made a mistake by editing the original candidate 
*		name. It changed from Lars Fagerland to Styrk Lothe. We have to make 
*		this change, otherwise the candidates run in the same district and year 
*		in different parties

* b) Focus on years 1953-1985

keep if inrange(year,1953,1985)
*drop if candidatename_orig=="Knudsen, A. G." & year==1912
* note: drop duplicate

* c) Merge vote margins for Norway and generate running variable

merge 1:1 year partyname districtid pid using ///
	"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981.dta", ///
	gen(merge1)
drop if merge1 == 2
* Drop one votemargin w/out candidate
drop merge1

* d) Generate vote margin

gen votemargin_vote=votemargin/approvedvotesoverall
gen votemargin_elig=votemargin/electorate
rename margin votemargin_fiva


* 5. Generate lead election variable 

* (a) Gen unique pid
* Note: We do not know how and when Jon drops duplicate pids. Instead, 
* 		we use our way and assign a unique pid per district

gen double pid_distid=districtid*100000+pid
egen time_id=group(year)
xtset pid_distid time_id

* (b) Gen unique pid

sort pid_distid time_id
gen elected_F1=F1.elected
replace elected_F1=0 if elected_F1==.

save "$path\02_Processed_data\13_Running_variable\FivaSmithJune2019_processed.dta", replace
* Note: For R file Figure 6 - Histogram Norway

* 6. Restrict sample period and define options (margin, kernel, cluster, polynomials)

keep if inrange(year,1953,1981)
global margin votemargin_elig
local margin $margin
global kernel triangular // Options: triangular (default), epanechnikov, uniform
global secl cluster // Options: cluster nncluster
global p 1
tab party, gen(party_)
tab district, gen(district_)

global covs age female first_year last_year ever_seat party_1 party_2 party_3 party_4 party_5 party_6 party_7 party_8 party_9 party_10  party_11 party_12 party_13 party_14 party_15 party_16 party_17 party_18 party_19 district_1 district_2 district_3 district_4 district_5 district_6 district_7 district_8 district_9 district_10  district_11 district_12 district_13 district_14 district_15 district_16 district_17 district_18 district_19 district_20 district_21 

* 7. RD estimation

rdrobust elected_F1 votemargin_fiva, p($p) bwselect(mserd) kernel($kernel)  ///
	vce($secl pid_distid) all

rdrobust elected_F1 $margin, p($p) bwselect(mserd) kernel($kernel)  ///
	vce($secl pid_distid) all, if votemargin_fiva!=.
global bw_sample_fiva =e(h_l)
global bbw_sample_fiva =e(b_l)

rdrobust elected_F1 $margin, p($p) bwselect(mserd) kernel($kernel)  ///
	vce($secl pid_distid) all
global bw_full_sample =e(h_l)
global bbw_full_sample =e(b_l)
global bw_opt_p1 =e(h_l)
global bbw_opt_p1 =e(b_l)

rdrobust elected_F1 $margin, p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl pid_distid) all
global bw_full_sample =e(h_l)
global bbw_full_sample =e(b_l)
global bw_opt_p2 =e(h_l)
global bbw_opt_p2 =e(b_l)
	


* 8. RD plots

rdplot elected_F1 $margin ///
	if -$bw_sample_fiva <= $margin & $margin<= $bw_sample_fiva ///
	& votemargin_fiva < . & elected_F1 < ., ///
	p($p) kernel($kernel) h($bw_sample_fiva $bw_sample_fiva) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\norway_fiva_sample.pdf", ///
	as(pdf) replace	
	
rdplot elected_F1 $margin ///
	if -$bw_full_sample <= $margin & $margin<= $bw_full_sample & elected_F1 < ., ///
	p($p) kernel($kernel) h($bw_full_sample $bw_full_sample) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\norway_full_sample.pdf", ///
	as(pdf) replace
	
* 9. Bandwidth figures


*qui {
mat pointestimate = J(1,15,.)
mat CI = J(2,15,.)
local col=1
forv i = 0.1(0.1)1.6 {
local bw_opt_p1_tmp=$bw_full_sample*`i'
local bbw_opt_p1_tmp=$bbw_full_sample*`i'
rdrobust  elected_F1 $margin if elected_F1 < ., p(1) h(`bw_opt_p1_tmp') b(`bbw_opt_p1_tmp') ///
	kernel($kernel) vce($secl pid_distid) all
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
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\norway_bw.pdf", ///
	as(pdf) replace
clear matrix
*}



* 10. Balance tables

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
	kernel($kernel) vce($secl pid) all
est store p1
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1" local ++cp1_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_10
local bw_opt_p1_div2=$bw_opt_p1/2
local bbw_opt_p1_div2=$bbw_opt_p1/2
rdrobust  `z' $margin if elected_F1 < ., p(1) h(`bw_opt_p1_div2') b(`bbw_opt_p1_div2') ///
	kernel($kernel) vce($secl pid) all 
est store p1_div2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_div2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_div2_10
rdrobust  `z' $margin if elected_F1 < ., p(2) h($bw_opt_p2) b($bbw_opt_p2) ///
	kernel($kernel) vce($secl pid) all
est store p2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_10
local bw_opt_p2_div2=$bw_opt_p2/2
local bbw_opt_p2_div2=$bbw_opt_p2/2
rdrobust  `z' $margin if elected_F1 < ., p(2) h(`bw_opt_p2_div2') b(`bbw_opt_p2_div2') ///
	kernel($kernel) vce($secl pid) all
est store p2_div2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_div2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_div2_10

if "`z'" == "elected_F1" {
	esttab p1 p1_div2 p2 p2_div2 using "$path\05_Texts_and_presns\01_Running_variable\tables\\no_`margin'_balance.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`label'") varwidth(21) star(* 0.10 ** 0.05 *** 0.01) fragment ///
	sca("h_l Bandwidth" "N_h_l Left obs" "N_h_r Right obs")  ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes replace

	}

else {
	esttab p1 p1_div2 p2 p2_div2 using "$path\05_Texts_and_presns\01_Running_variable\tables\\no_`margin'_balance.tex", ///
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
esttab matrix(rrates) using "$path\05_Texts_and_presns\01_Running_variable\tables\\no_`margin'_balance.tex",  ///
	nomti nonum append fragment
