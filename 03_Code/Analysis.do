********************************************************************************
*                                                                              *
*           Politicians as Directors - Evidence from Close Elections           *
*                                                                              *
********************************************************************************

clear
set more off
version 17

*** Set paths

global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
global path_ol "C:\Schmidlu\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"

*global path "C:\Current\Dropbox\Projekt Nationalräte"
*path_ol "C:\Current\Dropbox\Apps\Overleaf\Political_Rents"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

*global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path_ol "E:\12. Cloud\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"

***********************************************************************
* A) Estimation of treatment effect over time using rdrobust package
***********************************************************************

* Note: These results are in the files results_round_1 and results_round2 in the
*		following dropbox folder: Projekt Nationalräte\04_Results\04_Political_Rents

* (i) Main effects

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

foreach outcome of varlist ///
	i_all_c1 i_all_c2 i_all_s1 i_all_s2 ///
	i_lrg_c1 i_lrg_c2 i_lrg_s1 i_lrg_s2 ///
	i_sml_c1 i_sml_c2 i_sml_s1 i_sml_s2 ///
	i_prs_c1 i_prs_c2 i_prs_s1 i_prs_s2 ///
	n_all_sum_c1 n_all_sum_c2 n_all_sum_s1 n_all_sum_s2 ///
	n_lrg_sum_c1 n_lrg_sum_c2 n_lrg_sum_s1 n_lrg_sum_s2 ///
	n_sml_sum_c1 n_sml_sum_c2 n_sml_sum_s1 n_sml_sum_s2 ///
	n_prs_sum_c1 n_prs_sum_c2 n_prs_sum_s1 n_prs_sum_s2 ///
	n_all_avg_c1 n_all_avg_c2 n_all_avg_s1 n_all_avg_s2 ///
	n_lrg_avg_c1 n_lrg_avg_c2 n_lrg_avg_s1 n_lrg_avg_s2 ///
	n_sml_avg_c1 n_sml_avg_c2 n_sml_avg_s1 n_sml_avg_s2 ///
	n_prs_avg_c1 n_prs_avg_c2 n_prs_avg_s1 n_prs_avg_s2 {
* prc_lk_s1 rcl_lk_s1 f1_lk_s1 prc_mn_s1 rcl_mn_s1 f1_mn_s1 

*	qui {
	clear matrix

foreach k in tri uni epa {
	
	foreach m in ///
		rb_ob_p1_ay cl_ob_p1_ay rb_hb_p1_ay cl_hb_p1_ay ///
		rb_ob_p2_ay cl_ob_p2_ay rb_hb_p2_ay cl_hb_p2_ay ///
		rb_ob_p1_1y cl_ob_p1_1y rb_hb_p1_1y cl_hb_p1_1y ///
		rb_ob_p2_1y cl_ob_p2_1y rb_hb_p2_1y cl_hb_p2_1y {
	mat coef_`m' = J(1,20,.)
	mat CI_`m' = J(2,20,.)
	}
	local i = 1

	foreach var of varlist `outcome'_L4 `outcome'_L3 `outcome'_L2 `outcome'_L1 ///
		`outcome' `outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4 ///
		`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8 {

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_ay = e(h_l)
		local bw_half_p1_ay = `bw_opt_p1_ay'/2
		
		* Robust
		mat coef_rb_ob_p1_ay[1,`i']=e(tau_bc)
		mat CI_rb_ob_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p1_ay[2,`i'] = e(ci_r_rb)

		* Conventional
		mat coef_cl_ob_p1_ay[1,`i']=e(tau_cl)
		mat CI_cl_ob_p1_ay[1,`i'] = e(ci_l_cl)
		mat CI_cl_ob_p1_ay[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1_ay') rho(1) 
			
		* Robust
		mat coef_rb_hb_p1_ay[1,`i']=e(tau_bc)
		mat CI_rb_hb_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p1_ay[2,`i'] = e(ci_r_rb)

		* Conventional
		mat coef_cl_hb_p1_ay[1,`i']=e(tau_cl)
		mat CI_cl_hb_p1_ay[1,`i'] = e(ci_l_cl)
		mat CI_cl_hb_p1_ay[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_ay = e(h_l)
		local bw_half_p2_ay = `bw_opt_p2_ay'/2
		
		* Robust
		mat coef_rb_ob_p2_ay[1,`i']=e(tau_bc)
		mat CI_rb_ob_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p2_ay[2,`i'] = e(ci_r_rb)

		* Conventional
		mat coef_cl_ob_p2_ay[1,`i']=e(tau_cl)
		mat CI_cl_ob_p2_ay[1,`i'] = e(ci_l_cl)
		mat CI_cl_ob_p2_ay[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p2_ay') rho(1) 
			
		* Robust
		mat coef_rb_hb_p2_ay[1,`i']=e(tau_bc)
		mat CI_rb_hb_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p2_ay[2,`i'] = e(ci_r_rb)

		* Conventional
		mat coef_cl_hb_p2_ay[1,`i']=e(tau_cl)
		mat CI_cl_hb_p2_ay[1,`i'] = e(ci_l_cl)
		mat CI_cl_hb_p2_ay[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_1y = e(h_l)
		local bw_half_p1_1y = `bw_opt_p1_1y'/2
		
		* Robust
		mat coef_rb_ob_p1_1y[1,`i']=e(tau_bc)
		mat CI_rb_ob_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p1_1y[2,`i'] = e(ci_r_rb)

		* Conventional
		mat coef_cl_ob_p1_1y[1,`i']=e(tau_cl)
		mat CI_cl_ob_p1_1y[1,`i'] = e(ci_l_cl)
		mat CI_cl_ob_p1_1y[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p1_1y') rho(1) 
			
		* Robust
		mat coef_rb_hb_p1_1y[1,`i']=e(tau_bc)
		mat CI_rb_hb_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p1_1y[2,`i'] = e(ci_r_rb)

		* Conventional
		mat coef_cl_hb_p1_1y[1,`i']=e(tau_cl)
		mat CI_cl_hb_p1_1y[1,`i'] = e(ci_l_cl)
		mat CI_cl_hb_p1_1y[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_1y = e(h_l)
		local bw_half_p2_1y = `bw_opt_p2_1y'/2
		
		* Robust
		mat coef_rb_ob_p2_1y[1,`i']=e(tau_bc)
		mat CI_rb_ob_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p2_1y[2,`i'] = e(ci_r_rb)

		* Conventional
		mat coef_cl_ob_p2_1y[1,`i']=e(tau_cl)
		mat CI_cl_ob_p2_1y[1,`i'] = e(ci_l_cl)
		mat CI_cl_ob_p2_1y[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p2_1y') rho(1) 
			
		* Robust
		mat coef_rb_hb_p2_1y[1,`i']=e(tau_bc)
		mat CI_rb_hb_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p2_1y[2,`i'] = e(ci_r_rb)

		* Conventional
		mat coef_cl_hb_p2_1y[1,`i']=e(tau_cl)
		mat CI_cl_hb_p2_1y[1,`i'] = e(ci_l_cl)
		mat CI_cl_hb_p2_1y[2,`i'] = e(ci_r_cl)
				
		local ++i
		
		}
	local title: var label `outcome'

	foreach m in ///
		rb_ob_p1_ay cl_ob_p1_ay rb_hb_p1_ay cl_hb_p1_ay ///
		rb_ob_p2_ay cl_ob_p2_ay rb_hb_p2_ay cl_hb_p2_ay ///
		rb_ob_p1_1y cl_ob_p1_1y rb_hb_p1_1y cl_hb_p1_1y ///
		rb_ob_p2_1y cl_ob_p2_1y rb_hb_p2_1y cl_hb_p2_1y {

		mat colnames coef_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" ///
			"0 yr" "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
		mat colnames CI_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" "0 yr" ///
			"+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" "+6 yrs" ///
			"+7 yrs" "+8 yrs"
	}
	
	coefplot matrix(coef_rb_ob_p1_ay), ci(CI_rb_ob_p1_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Robust, optimal bandwidth, 1st order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_ob_p1_`k'_ay_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_cl_ob_p1_ay), ci(CI_cl_ob_p1_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Conventional, optimal bandwidth, 1st order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_cl_ob_p1_`k'_ay_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p1_ay), ci(CI_rb_hb_p1_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Robust, half optimal bandwidth, 1st order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_hb_p1_`k'_ay_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_cl_hb_p1_ay), ci(CI_cl_hb_p1_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Conventional, half optimal bandwidth, 1st order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_cl_hb_p1_`k'_ay_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_ob_p2_ay), ci(CI_rb_ob_p2_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Robust, optimal bandwidth, 2nd order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_ob_p2_`k'_ay_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_cl_ob_p2_ay), ci(CI_cl_ob_p2_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Conventional, optimal bandwidth, 2nd order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_cl_ob_p2_`k'_ay_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p2_ay), ci(CI_rb_hb_p2_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Robust, half optimal bandwidth, 2nd order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_hb_p2_`k'_ay_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_cl_hb_p2_ay), ci(CI_cl_hb_p2_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Conventional, half optimal bandwidth, 2nd order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_cl_hb_p2_`k'_ay_`outcome'.pdf", as(pdf) replace
	
	coefplot matrix(coef_rb_ob_p1_1y), ci(CI_rb_ob_p1_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Robust, optimal bandwidth, 1st order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_ob_p1_`k'_1y_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_cl_ob_p1_1y), ci(CI_cl_ob_p1_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Conventional, optimal bandwidth, 1st order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_cl_ob_p1_`k'_1y_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p1_1y), ci(CI_rb_hb_p1_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Robust, half optimal bandwidth, 1st order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_hb_p1_`k'_1y_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_cl_hb_p1_1y), ci(CI_cl_hb_p1_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Conventional, half optimal bandwidth, 1st order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_cl_hb_p1_`k'_1y_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_ob_p2_1y), ci(CI_rb_ob_p2_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Robust, optimal bandwidth, 2nd order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_ob_p2_`k'_1y_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_cl_ob_p2_1y), ci(CI_cl_ob_p2_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Conventional, optimal bandwidth, 2nd order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_cl_ob_p2_`k'_1y_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p2_1y), ci(CI_rb_hb_p2_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Robust, half optimal bandwidth, 2nd order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_hb_p2_`k'_1y_`outcome'.pdf", as(pdf) replace

	coefplot matrix(coef_cl_hb_p2_1y), ci(CI_cl_hb_p2_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "Conventional, half optimal bandwidth, 2nd order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_cl_hb_p2_`k'_1y_`outcome'.pdf", as(pdf) replace
	
	clear matrix
*	}
}
}

* (ii) Explore sample composition between optimal and half-optimal bandwidth

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

mat coef_rb_ob = J(1,20,.)
mat CI_rb_ob = J(2,20,.)
mat coef_rb_hb = J(1,20,.)
mat CI_rb_hb = J(2,20,.)

local outcome i_lrg_c1
local i 1
foreach var of varlist `outcome'_L4 `outcome'_L3 `outcome'_L2 `outcome'_L1 ///
		`outcome' `outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4 ///
		`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8 {	
	
rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(tri) vce(cluster ID_num) all rho(1)
* Bandwidth
local bw_opt = e(h_l)
local bw_half = `bw_opt'/2

* Mark sample

gen s_`var'_ob=0
replace s_`var'_ob=1 if inrange(votemargin_rel,- `bw_opt', `bw_opt')

* Out
mat coef_rb_ob[1,`i']=e(tau_bc)
mat CI_rb_ob[1,`i'] = e(ci_l_rb)
mat CI_rb_ob[2,`i'] = e(ci_r_rb)

rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(tri) vce(cluster ID_num) all rho(1) h(`bw_half')
			
* Mark sample

gen s_`var'_hb=0
replace s_`var'_hb=1 if inrange(votemargin_rel,- `bw_half', `bw_half')

* Out
mat coef_rb_hb[1,`i']=e(tau_bc)
mat CI_rb_hb[1,`i'] = e(ci_l_rb)
mat CI_rb_hb[2,`i'] = e(ci_r_rb)
local ++i
}

mat colnames coef_rb_ob = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" ///
			"0 yr" "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
mat colnames CI_rb_ob = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" "0 yr" ///
			"+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" "+6 yrs" ///
			"+7 yrs" "+8 yrs"
			
mat colnames coef_rb_hb = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" ///
			"0 yr" "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
mat colnames CI_rb_hb = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" "0 yr" ///
			"+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" "+6 yrs" ///
			"+7 yrs" "+8 yrs"


coefplot matrix(coef_rb_ob), ci(CI_rb_ob) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'" "Robust, optimal bandwidth, 1st order polynomial, tri, 1st election")
		
coefplot matrix(coef_rb_hb), ci(CI_rb_hb) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'" "Robust, optimal bandwidth, 1st order polynomial, tri, 1st election")

preserve
keep year list* birthyear sex canton votemargin_rel s_*
save "$path\02_Processed_data\20_Analysis\sample_composition_`outcome'.dta", replace
restore




* (ii) Analysis of differences

use "$path\02_Processed_data\20_Analysis\sample_composition_i_lrg_c1.dta", clear

* a) Recode variables

gen leftist=0
replace leftist=1 if listname_bfs=="SP/PS" | listname_bfs=="GPS/PES" | ///
	 listname_bfs=="PdA/PST"

gen centrist=0
replace centrist=1 if listname_bfs=="CVP/PDC" | listname_bfs=="FDP/PLR (PRD)" | ///
	 listname_bfs=="LPS/PLS"

gen rightist=0
replace rightist=1 if listname_bfs=="SVP/UDC" | listname_bfs=="Rep./Rép."  | /// 
	listname_bfs=="Lega"  | listname_bfs=="FPS/PSL" | listname_bfs=="MCR"
	
gen canton_lrg=0 
replace canton_lrg=1 if canton=="BE" | canton=="VD"  | canton=="ZH"  


* b) Loop through outcome variables and check means of covariates



local outcome i_lrg_c1

foreach covar in canton_lrg year birthyear sex leftist centrist rightist  { 
mat mean_ob = J(1,20,.)
mat ci_ob = J(2,20,.)
mat mean_hb = J(1,20,.)
mat ci_hb = J(2,20,.)
di "variable: `outcome'"
local i 1
foreach var in s_`outcome'_L4 s_`outcome'_L3 s_`outcome'_L2 s_`outcome'_L1 ///
		s_`outcome' s_`outcome'_F1 s_`outcome'_F2 s_`outcome'_F3 s_`outcome'_F4 ///
		s_`outcome'_F5 s_`outcome'_F6 s_`outcome'_F7 s_`outcome'_F8 {
di "variable: `var'"
sum  `covar'  if `var'_ob==1
mat mean_ob[1,`i']=r(mean)
mat ci_ob[1,`i'] =r(N)

mat ci_ob[1,`i'] =r(mean)-1.96*sqrt(r(sd)^2/r(N))
mat ci_ob[2,`i'] =r(mean)+1.96*sqrt(r(sd)^2/r(N))

sum `covar'  if `var'_hb==1 
mat mean_hb[1,`i']=r(mean)
mat ci_hb[1,`i'] =r(mean)-1.96*sqrt(r(sd)^2/r(N))
mat ci_hb[2,`i'] =r(mean)+1.96*sqrt(r(sd)^2/r(N))

*sum year birthyear sex leftist centrist rightist canton_lrg if `var'_ob==1
*sum year birthyear sex leftist centrist rightist canton_lrg if `var'_hb==1 
*tab listname_bfs if `var'_ob==1
*tab listname_bfs if `var'_hb==1 
*tab canton if `var'_ob==1
*tab canton if `var'_hb==1 	
local ++i
}

mat colnames mean_ob = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" ///
			"0 yr" "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
mat colnames ci_ob = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" "0 yr" ///
			"+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" "+6 yrs" ///
			"+7 yrs" "+8 yrs"
			
mat colnames mean_hb = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" ///
			"0 yr" "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
mat colnames ci_hb = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" "0 yr" ///
			"+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" "+6 yrs" ///
			"+7 yrs" "+8 yrs"
			
coefplot matrix(mean_ob), ci(ci_ob) vertical nolabel  ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("Mean of `covar', Outcome: `outcome' optimal bandwidth" )
		
graph export "$path_ol\figures\fig_mean_`outcome'_`covar'_ob.pdf", as(pdf) replace


coefplot matrix(mean_hb), ci(ci_hb) vertical nolabel  ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("Mean of `covar', Outcome: `outcome' half-optimal bandwidth" )

graph export "$path_ol\figures\fig_mean_`outcome'_`covar'_hb.pdf", as(pdf) replace
clear matrix		
}

***********************************************************************
* B) Estimation of treatment effect over time using rdrobust package
***********************************************************************

* Note: These results are in the files results_round_1 and results_round2 in the
*		following dropbox folder: Projekt Nationalräte\04_Results\04_Political_Rents

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .