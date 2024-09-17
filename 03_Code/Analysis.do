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

global path "C:\Current\Dropbox\Projekt Nationalräte"
global path_ol "C:\Current\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"

global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"
global path_ol "E:\12. Cloud\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"

***********************************************************************
* A) Estimation of treatment effect over time using rdrobust package
***********************************************************************

* Note: These results are in the files results_round_1 and results_round2 in the
*		following dropbox folder: Projekt Nationalräte\04_Results\04_Political_Rents

* (i) Main effects (round 1 and 2)

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

* a) Drop those with missing vote margin
keep if votemargin_rel < .

* b) Generate variable for first election and votemargin left/right

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

gen votemargin_rel_l=(1-elected)*votemargin_rel
gen votemargin_rel_r=elected*votemargin_rel

* c) Set up matrices

mat coef_rb_ob = J(1,20,.)
mat CI_rb_ob = J(2,20,.)
mat coef_rb_hb = J(1,20,.)
mat CI_rb_hb = J(2,20,.)

* d) Loop over all outcome variables

local outcome i_lrg_c1
local i 1
foreach var of varlist `outcome'_L4 `outcome'_L3 `outcome'_L2 `outcome'_L1 ///
		`outcome' `outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4 ///
		`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8 {
			
* d1) Optimal bandwidth estimation			
	
rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(tri) vce(cluster ID_num) all rho(1)
local bw_opt = e(h_l)
local bw_half = `bw_opt'/2

* d2) Mark sample

*gen s_`var'_ob=0
*replace s_`var'_ob=1 if inrange(votemargin_rel,- `bw_opt', `bw_opt') & `var'!=.

gen s_`var'_ob=0
replace s_`var'_ob=1-abs(votemargin_rel/`bw_opt') if inrange(votemargin_rel,-`bw_opt', `bw_opt') & `var'!=.



* d3) Results to matrices 

mat coef_rb_ob[1,`i']=e(tau_bc)
mat CI_rb_ob[1,`i'] = e(ci_l_rb)
mat CI_rb_ob[2,`i'] = e(ci_r_rb)


* d4) Halfoptimal bandwidth

rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(tri) vce(cluster ID_num) all rho(1) h(`bw_half')
			
* d4) Mark sample

*gen s_`var'_hb=0
*replace s_`var'_hb=1 if inrange(votemargin_rel,- `bw_half', `bw_half') & `var'!=.

gen s_`var'_hb=0
replace s_`var'_hb=1-abs(votemargin_rel/`bw_half') if inrange(votemargin_rel,-`bw_half', `bw_half') & `var'!=.

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
keep year list* birthyear sex canton votemargin_rel s_* ID_num `outcome'*
save "$path\02_Processed_data\20_Analysis\sample_composition_`outcome'.dta", replace
restore


* (ii) Analysis of covariate differences across samples (round 2 mean comparison)

use "$path\02_Processed_data\20_Analysis\sample_composition_i_lrg_c1.dta", clear

* a) Recode variables

gen leftist=0 if year>=1971
replace leftist=1 if listname_bfs=="SP/PS" | listname_bfs=="GPS/PES" | ///
	 listname_bfs=="PdA/PST"

gen centrist=0 if year>=1971
replace centrist=1 if listname_bfs=="CVP/PDC" | listname_bfs=="FDP/PLR (PRD)" | ///
	 listname_bfs=="LPS/PLS"

gen rightist=0 if year>=1971
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
foreach var in `outcome'_L4 `outcome'_L3 `outcome'_L2 `outcome'_L1 ///
		`outcome' `outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4 ///
		`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8 {
di "variable: `var'"

reg  `covar' [aweight=s_`var'_ob] , vce(cluster ID_num)

mat mean_ob[1,`i']=_b[_cons]
local ser =_se[_cons]
mat ci_ob[1,`i'] =mean_ob[1,`i']-1.96*`ser'
mat ci_ob[2,`i'] =mean_ob[1,`i']+1.96*`ser'

reg  `covar' [aweight=s_`var'_hb] , vce(cluster ID_num)

mat mean_hb[1,`i']=_b[_cons]
local ser =_se[_cons]
mat ci_hb[1,`i'] =mean_hb[1,`i']-1.96*`ser'
mat ci_hb[2,`i'] =mean_hb[1,`i']+1.96*`ser'

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
* B) Estimation of treatment effect over time using regression function
***********************************************************************

* Note: These results are in the files results_round_1 and results_round2 in the
*		following dropbox folder: Projekt Nationalräte\04_Results\04_Political_Rents


use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear

* a) Drop those with missing vote margin
keep if votemargin_rel < .

* b) Generate variable for first election and votemargin left/right

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

gen votemargin_rel_l=(1-elected)*votemargin_rel
gen votemargin_rel_r=elected*votemargin_rel

*drop if votemargin_rel==0
replace elected=1 if votemargin_rel==0

* c) Set up matrices

	mat coef_rb_ob = J(1,20,.)
	mat CI_rb_ob = J(2,20,.)
	mat coef_rb_hb = J(1,20,.)
	mat CI_rb_hb = J(2,20,.)
	mat coef_rg_ob = J(1,20,.)
	mat CI_rg_ob = J(2,20,.)
	mat coef_rg_hb = J(1,20,.)
	mat CI_rg_hb = J(2,20,.)

	* d) Loop over all outcome variables

	local outcome i_lrg_c1
	local k tri
	local i 1
	foreach var of varlist `outcome' {
				
	* d1) Optimal bandwidth estimation			
		
	rdrobust `var' votemargin_rel , p(1) bwselect(mserd) ///
				kernel(`k') all rho(1) vce(hc2) 
	local bw_opt = e(h_l)
	di `bw_opt'
	*local bw_half = `bw_opt'/2

	* d3) Results to matrices 

	mat coef_rb_ob[1,`i']=e(tau_bc)
	mat CI_rb_ob[1,`i'] = e(ci_l_rb)
	mat CI_rb_ob[2,`i'] = e(ci_r_rb)

	* d4) Implement regression version

	gen weight_ob=1-abs(votemargin_rel/`bw_opt') if inrange(votemargin_rel,-`bw_opt', `bw_opt')

	reg  `var' elected votemargin_rel_l votemargin_rel_r /// 
		if inrange(votemargin_rel,-`bw_opt', `bw_opt') [aweight=weight_ob] , ///
		 vce(hc2)

	mat coef_rg_ob[1,`i']=_b[elected]
	mat CI_rg_ob[1,`i'] = _b[elected] - invttail(e(df_r),.025)*_se[elected]
	mat CI_rg_ob[2,`i'] = _b[elected] + invttail(e(df_r),.025)*_se[elected]

	/*
	* d5) Halfoptimal bandwidth

	rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
				kernel(`k') vce(cluster ID_num) all rho(1) h(`bw_half')
				
	mat coef_rb_hb[1,`i']=e(tau_bc)
	mat CI_rb_hb[1,`i'] = e(ci_l_rb)
	mat CI_rb_hb[2,`i'] = e(ci_r_rb)

	* d6) Implement regression version

	*gen weigth_ob=1-abs(votemargin_rel/`bw_opt')

	reg  `var' elected votemargin_rel_l votemargin_rel_r, /// 
		cluster(ID_num), if inrange(votemargin_rel,- `bw_half', `bw_half')
		
	mat coef_rg_hb[1,`i']=_b[elected]
	mat CI_rg_hb[1,`i'] = _b[elected] - invttail(e(df_r),.025)*_se[elected]
	mat CI_rg_hb[2,`i'] = _b[elected] + invttail(e(df_r),.025)*_se[elected]
	*/
	local ++i
	}

******************************************
* C) Subsamples 
******************************************

* Note: These results are in the files results_round_4 in the
*		following dropbox folder: Projekt Nationalräte\04_Results\04_Political_Rents

* (i) Early and late

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {

*	qui {
	clear matrix

foreach k in tri {
foreach s in early late {
	
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
		
		local new_era_year 1975
		preserve
		if "`s'"=="early" {
		keep if year< `new_era_year'	
		} 
		else {
		keep if year>= `new_era_year'		
		}

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_ay = e(h_l)
		local bw_half_p1_ay = `bw_opt_p1_ay'/2
		
		* Robust
		mat coef_rb_ob_p1_ay[1,`i']=e(tau_bc)
		mat CI_rb_ob_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1_ay') rho(1) 
			
		* Robust
		mat coef_rb_hb_p1_ay[1,`i']=e(tau_bc)
		mat CI_rb_hb_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_ay = e(h_l)
		local bw_half_p2_ay = `bw_opt_p2_ay'/2
		
		* Robust
		mat coef_rb_ob_p2_ay[1,`i']=e(tau_bc)
		mat CI_rb_ob_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p2_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p2_ay') rho(1) 
			
		* Robust
		mat coef_rb_hb_p2_ay[1,`i']=e(tau_bc)
		mat CI_rb_hb_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p2_ay[2,`i'] = e(ci_r_rb)
		
		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_1y = e(h_l)
		local bw_half_p1_1y = `bw_opt_p1_1y'/2
		
		* Robust
		mat coef_rb_ob_p1_1y[1,`i']=e(tau_bc)
		mat CI_rb_ob_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p1_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p1_1y') rho(1) 
			
		* Robust
		mat coef_rb_hb_p1_1y[1,`i']=e(tau_bc)
		mat CI_rb_hb_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p1_1y[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_1y = e(h_l)
		local bw_half_p2_1y = `bw_opt_p2_1y'/2
		
		* Robust
		mat coef_rb_ob_p2_1y[1,`i']=e(tau_bc)
		mat CI_rb_ob_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p2_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p2_1y') rho(1) 
			
		* Robust
		mat coef_rb_hb_p2_1y[1,`i']=e(tau_bc)
		mat CI_rb_hb_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p2_1y[2,`i'] = e(ci_r_rb)
		
		restore
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
		title("`outcome'_`s': `title'" "Robust, optimal bandwidth, 1st order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_ob_p1_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p1_ay), ci(CI_rb_hb_p1_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, half optimal bandwidth, 1st order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_hb_p1_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_ob_p2_ay), ci(CI_rb_ob_p2_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, optimal bandwidth, 2nd order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_ob_p2_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p2_ay), ci(CI_rb_hb_p2_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, half optimal bandwidth, 2nd order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_hb_p2_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_ob_p1_1y), ci(CI_rb_ob_p1_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, optimal bandwidth, 1st order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_ob_p1_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p1_1y), ci(CI_rb_hb_p1_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, half optimal bandwidth, 1st order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_hb_p1_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_ob_p2_1y), ci(CI_rb_ob_p2_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, optimal bandwidth, 2nd order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_ob_p2_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p2_1y), ci(CI_rb_hb_p2_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, half optimal bandwidth, 2nd order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_hb_p2_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace
	
	clear matrix
*	}
}
}
}


* (ii) Parties

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

* a) Recode variables


gen leftist=0 if year>=1971
replace leftist=1 if listname=="CSP/PCS"
replace leftist=1 if listname=="DSP/PSD"
replace leftist=1 if listname=="F-GA/VA-F"
replace leftist=1 if listname=="GA/ASV"
replace leftist=1 if listname=="GPS/PES"
replace leftist=1 if listname=="POCH"
replace leftist=1 if listname=="PSA"
replace leftist=1 if listname=="PdA/PST"
replace leftist=1 if listname=="SAP/PSO"
replace leftist=1 if listname=="SP/PS"
replace leftist=1 if listname=="Sol."

gen centrist=0 if year>=1971
replace centrist=1 if listname=="BDP/PBD"
replace centrist=1 if listname=="CVP/PDC"
replace centrist=1 if listname=="EVP/PEV"
replace centrist=1 if listname=="FDP/PLR (PRD)"
replace centrist=1 if listname=="GLP/PVL"
replace centrist=1 if listname=="LPS/PLS"
replace centrist=1 if listname=="LdU/AdI"

gen rightist=0 if year>=1971
replace rightist=1 if listname=="EDU/UDF"
replace rightist=1 if listname=="FPS/PSL"
replace rightist=1 if listname=="Lega"
replace rightist=1 if listname=="MCR"
replace rightist=1 if listname=="Rep./Rép."
replace rightist=1 if listname=="SD/DS"
replace rightist=1 if listname=="SVP/UDC"

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {

*	qui {
	clear matrix

foreach k in tri {
foreach s of varlist leftist centrist rightist{
	
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
		
		local new_era_year 1975
		preserve
		keep if `s'==1

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_ay = e(h_l)
		local bw_half_p1_ay = `bw_opt_p1_ay'/2
		
		* Robust
		mat coef_rb_ob_p1_ay[1,`i']=e(tau_bc)
		mat CI_rb_ob_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1_ay') rho(1) 
			
		* Robust
		mat coef_rb_hb_p1_ay[1,`i']=e(tau_bc)
		mat CI_rb_hb_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_ay = e(h_l)
		local bw_half_p2_ay = `bw_opt_p2_ay'/2
		
		* Robust
		mat coef_rb_ob_p2_ay[1,`i']=e(tau_bc)
		mat CI_rb_ob_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p2_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p2_ay') rho(1) 
			
		* Robust
		mat coef_rb_hb_p2_ay[1,`i']=e(tau_bc)
		mat CI_rb_hb_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p2_ay[2,`i'] = e(ci_r_rb)
		
		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_1y = e(h_l)
		local bw_half_p1_1y = `bw_opt_p1_1y'/2
		
		* Robust
		mat coef_rb_ob_p1_1y[1,`i']=e(tau_bc)
		mat CI_rb_ob_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p1_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p1_1y') rho(1) 
			
		* Robust
		mat coef_rb_hb_p1_1y[1,`i']=e(tau_bc)
		mat CI_rb_hb_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p1_1y[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_1y = e(h_l)
		local bw_half_p2_1y = `bw_opt_p2_1y'/2
		
		* Robust
		mat coef_rb_ob_p2_1y[1,`i']=e(tau_bc)
		mat CI_rb_ob_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_ob_p2_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p2_1y') rho(1) 
			
		* Robust
		mat coef_rb_hb_p2_1y[1,`i']=e(tau_bc)
		mat CI_rb_hb_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_rb_hb_p2_1y[2,`i'] = e(ci_r_rb)
		
		restore
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
		title("`outcome'_`s': `title'" "Robust, optimal bandwidth, 1st order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_ob_p1_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p1_ay), ci(CI_rb_hb_p1_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, half optimal bandwidth, 1st order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_hb_p1_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_ob_p2_ay), ci(CI_rb_ob_p2_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, optimal bandwidth, 2nd order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_ob_p2_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p2_ay), ci(CI_rb_hb_p2_ay) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, half optimal bandwidth, 2nd order polynomial, `k', all elections")
	graph export "$path_ol\figures\fig_rb_hb_p2_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_ob_p1_1y), ci(CI_rb_ob_p1_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, optimal bandwidth, 1st order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_ob_p1_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p1_1y), ci(CI_rb_hb_p1_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, half optimal bandwidth, 1st order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_hb_p1_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_ob_p2_1y), ci(CI_rb_ob_p2_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, optimal bandwidth, 2nd order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_ob_p2_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace

	coefplot matrix(coef_rb_hb_p2_1y), ci(CI_rb_hb_p2_1y) vertical nolabel yline(0) ///
		mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "Robust, half optimal bandwidth, 2nd order polynomial, `k', 1st election")
	graph export "$path_ol\figures\fig_rb_hb_p2_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace
	
	clear matrix
*	}
}
}
}


*************************************************
* D) Bias-corrected robust CIs and undersmoothing 
*************************************************

* Note: These results are in the files results_round_5 in the following dropbox
*       folder: Projekt Nationalräte\04_Results\04_Political_Rents
*       The estimates take seriously the points of Calonico et al. (2014, Stata
*       Journal, p. 918) that a) bias correction or b) undersmoothing are rele-
*       vant to correct biases in CIs - not point estimates - and that a) and b)
*       are alternative approaches to address this issue.

* Round 5

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {

*	qui {
	clear matrix

foreach k in tri {
	
	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {
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
		mat coef_est_ob_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1_ay') rho(1) 
			
		* Conventional
		mat coef_est_hb_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_ay[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_ay = e(h_l)
		local bw_half_p2_ay = `bw_opt_p2_ay'/2
		
		* Robust
		mat coef_est_ob_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p2_ay') rho(1) 
			
		* Conventional
		mat coef_est_hb_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_ay[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_1y = e(h_l)
		local bw_half_p1_1y = `bw_opt_p1_1y'/2
		
		* Robust
		mat coef_est_ob_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p1_1y') rho(1) 
			
		* Conventional
		mat coef_est_hb_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_1y[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_1y = e(h_l)
		local bw_half_p2_1y = `bw_opt_p2_1y'/2
		
		* Robust
		mat coef_est_ob_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p2_1y') rho(1) 
			
		* Conventional
		mat coef_est_hb_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_1y[2,`i'] = e(ci_r_cl)
		
		local ++i
		
		}
	local title: var label `outcome'

	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {

		mat colnames coef_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" ///
			"0 yr" "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
		mat colnames CI_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" "0 yr" ///
			"+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" "+6 yrs" ///
			"+7 yrs" "+8 yrs"
	}
	
		coefplot ///
		(matrix(coef_est_ob_p1_ay), ci(CI_est_ob_p1_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_ay), ci(CI_est_hb_p1_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p1_`k'_ay_`outcome'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_ay), ci(CI_est_ob_p2_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_ay), ci(CI_est_hb_p2_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p2_`k'_ay_`outcome'.pdf", as(pdf) replace
		
				coefplot ///
		(matrix(coef_est_ob_p1_1y), ci(CI_est_ob_p1_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_1y), ci(CI_est_hb_p1_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_p1_`k'_1y_`outcome'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_1y), ci(CI_est_ob_p2_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_1y), ci(CI_est_hb_p2_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_p2_`k'_1y_`outcome'.pdf", as(pdf) replace

	clear matrix
	*}
}
}




*****************************************************************
* E) Bias-corrected robust CIs and undersmoothing for subsamples
*****************************************************************

* Note: round 6

* (i) Early and late years

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {

*	qui {
	clear matrix

foreach k in tri {
foreach s in early late {
	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {
	mat coef_`m' = J(1,20,.)
	mat CI_`m' = J(2,20,.)
	}
	local i = 1

	foreach var of varlist `outcome'_L4 `outcome'_L3 `outcome'_L2 `outcome'_L1 ///
		`outcome' `outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4 ///
		`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8 {
			
		local new_era_year 1975
		preserve
		if "`s'"=="early" {
		keep if year< `new_era_year'	
		} 
		else {
		keep if year>= `new_era_year'		
		}			

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_ay = e(h_l)
		local bw_half_p1_ay = `bw_opt_p1_ay'/2
		
		* Robust
		mat coef_est_ob_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1_ay') rho(1) 
			
		* Conventional
		mat coef_est_hb_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_ay[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_ay = e(h_l)
		local bw_half_p2_ay = `bw_opt_p2_ay'/2
		
		* Robust
		mat coef_est_ob_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p2_ay') rho(1) 
			
		* Conventional
		mat coef_est_hb_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_ay[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_1y = e(h_l)
		local bw_half_p1_1y = `bw_opt_p1_1y'/2
		
		* Robust
		mat coef_est_ob_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p1_1y') rho(1) 
			
		* Conventional
		mat coef_est_hb_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_1y[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_1y = e(h_l)
		local bw_half_p2_1y = `bw_opt_p2_1y'/2
		
		* Robust
		mat coef_est_ob_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p2_1y') rho(1) 
			
		* Conventional
		mat coef_est_hb_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_1y[2,`i'] = e(ci_r_cl)
		
		local ++i
		restore
		}
	local title: var label `outcome'

	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {

		mat colnames coef_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" ///
			"0 yr" "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
		mat colnames CI_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" "0 yr" ///
			"+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" "+6 yrs" ///
			"+7 yrs" "+8 yrs"
	}
	
		coefplot ///
		(matrix(coef_est_ob_p1_ay), ci(CI_est_ob_p1_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_ay), ci(CI_est_hb_p1_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s', : `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p1_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_ay), ci(CI_est_ob_p2_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_ay), ci(CI_est_hb_p2_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "2nd order polynomial, `k',  all elections")
		graph export "$path_ol\figures\fig_p2_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p1_1y), ci(CI_est_ob_p1_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_1y), ci(CI_est_hb_p1_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "1st order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_p1_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_1y), ci(CI_est_ob_p2_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_1y), ci(CI_est_hb_p2_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "2nd order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_p2_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace

	clear matrix
	
	*}
}
}
}



* (ii) Parties

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

* a) Recode variables


gen leftist=0 if year>=1971
replace leftist=1 if listname=="CSP/PCS"
replace leftist=1 if listname=="DSP/PSD"
replace leftist=1 if listname=="F-GA/VA-F"
replace leftist=1 if listname=="GA/ASV"
replace leftist=1 if listname=="GPS/PES"
replace leftist=1 if listname=="POCH"
replace leftist=1 if listname=="PSA"
replace leftist=1 if listname=="PdA/PST"
replace leftist=1 if listname=="SAP/PSO"
replace leftist=1 if listname=="SP/PS"
replace leftist=1 if listname=="Sol."

gen centrist=0 if year>=1971
replace centrist=1 if listname=="BDP/PBD"
replace centrist=1 if listname=="CVP/PDC"
replace centrist=1 if listname=="EVP/PEV"
replace centrist=1 if listname=="FDP/PLR (PRD)"
replace centrist=1 if listname=="GLP/PVL"
replace centrist=1 if listname=="LPS/PLS"
replace centrist=1 if listname=="LdU/AdI"

gen rightist=0 if year>=1971
replace rightist=1 if listname=="EDU/UDF"
replace rightist=1 if listname=="FPS/PSL"
replace rightist=1 if listname=="Lega"
replace rightist=1 if listname=="MCR"
replace rightist=1 if listname=="Rep./Rép."
replace rightist=1 if listname=="SD/DS"
replace rightist=1 if listname=="SVP/UDC"


bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {

*	qui {
	clear matrix

foreach k in tri {
foreach s of varlist leftist centrist rightist{
	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {
	mat coef_`m' = J(1,20,.)
	mat CI_`m' = J(2,20,.)
	}
	local i = 1

	foreach var of varlist `outcome'_L4 `outcome'_L3 `outcome'_L2 `outcome'_L1 ///
		`outcome' `outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4 ///
		`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8 {
			
		preserve
		keep if `s'==1


		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_ay = e(h_l)
		local bw_half_p1_ay = `bw_opt_p1_ay'/2
		
		* Robust
		mat coef_est_ob_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1_ay') rho(1) 
			
		* Conventional
		mat coef_est_hb_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_ay[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_ay = e(h_l)
		local bw_half_p2_ay = `bw_opt_p2_ay'/2
		
		* Robust
		mat coef_est_ob_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p2_ay') rho(1) 
			
		* Conventional
		mat coef_est_hb_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_ay[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_1y = e(h_l)
		local bw_half_p1_1y = `bw_opt_p1_1y'/2
		
		* Robust
		mat coef_est_ob_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p1_1y') rho(1) 
			
		* Conventional
		mat coef_est_hb_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_1y[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_1y = e(h_l)
		local bw_half_p2_1y = `bw_opt_p2_1y'/2
		
		* Robust
		mat coef_est_ob_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p2_1y') rho(1) 
			
		* Conventional
		mat coef_est_hb_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_1y[2,`i'] = e(ci_r_cl)
		
		local ++i
		restore
		}
	local title: var label `outcome'

	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {

		mat colnames coef_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" ///
			"0 yr" "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
		mat colnames CI_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" "0 yr" ///
			"+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" "+6 yrs" ///
			"+7 yrs" "+8 yrs"
	}
	
		coefplot ///
		(matrix(coef_est_ob_p1_ay), ci(CI_est_ob_p1_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_ay), ci(CI_est_hb_p1_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s', : `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p1_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_ay), ci(CI_est_ob_p2_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_ay), ci(CI_est_hb_p2_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "2nd order polynomial, `k',  all elections")
		graph export "$path_ol\figures\fig_p2_`k'_ay_`outcome'_`s'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p1_1y), ci(CI_est_ob_p1_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_1y), ci(CI_est_hb_p1_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "1st order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_p1_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_1y), ci(CI_est_ob_p2_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_1y), ci(CI_est_hb_p2_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome'_`s': `title'" "2nd order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_p2_`k'_1y_`outcome'_`s'.pdf", as(pdf) replace

	clear matrix
	*}
}
}
}

************************************
* (F) Analysis at term level
************************************
* Note: round 7

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 {
	egen `outcome'_TL2=rowmax(`outcome'_L7 `outcome'_L6 `outcome'_L5 `outcome'_L4)
	egen `outcome'_TL1=rowmax(`outcome'_L3 `outcome'_L2 `outcome'_L1 `outcome')
	egen `outcome'_T0=rowmax(`outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4)
	egen `outcome'_TF1=rowmax(`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8)
}
	
foreach outcome of varlist ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {
	egen `outcome'_TL2=rowmean(`outcome'_L7 `outcome'_L6 `outcome'_L5 `outcome'_L4)
	egen `outcome'_TL1=rowmean(`outcome'_L3 `outcome'_L2 `outcome'_L1 `outcome')
	egen `outcome'_T0=rowmean(`outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4)
	egen `outcome'_TF1=rowmean(`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8)
}


bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {

*	qui {
	clear matrix

foreach k in tri {
	
	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {
	mat coef_`m' = J(1,20,.)
	mat CI_`m' = J(2,20,.)
	}
	local i = 1

	foreach var of varlist `outcome'_TL2 `outcome'_TL1 `outcome'_T0 `outcome'_TF1 {

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_ay = e(h_l)
		local bw_half_p1_ay = `bw_opt_p1_ay'/2
		
		* Robust
		mat coef_est_ob_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1_ay') rho(1) 
			
		* Conventional
		mat coef_est_hb_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_ay[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_ay = e(h_l)
		local bw_half_p2_ay = `bw_opt_p2_ay'/2
		
		* Robust
		mat coef_est_ob_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p2_ay') rho(1) 
			
		* Conventional
		mat coef_est_hb_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_ay[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p1_1y = e(h_l)
		local bw_half_p1_1y = `bw_opt_p1_1y'/2
		
		* Robust
		mat coef_est_ob_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p1_1y') rho(1) 
			
		* Conventional
		mat coef_est_hb_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_1y[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1)

		* Bandwidth
		local bw_opt_p2_1y = e(h_l)
		local bw_half_p2_1y = `bw_opt_p2_1y'/2
		
		* Robust
		mat coef_est_ob_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p2_1y') rho(1) 
			
		* Conventional
		mat coef_est_hb_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_1y[2,`i'] = e(ci_r_cl)
		
		local ++i
		
		}
	local title: var label `outcome'

	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {

		mat colnames coef_`m' = "-2 term" "-1 term" "0 term" "+1 term"
		mat colnames CI_`m' = "-2 term" "-1 term" "0 term" "+1 term"
	}
	
		coefplot ///
		(matrix(coef_est_ob_p1_ay), ci(CI_est_ob_p1_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_ay), ci(CI_est_hb_p1_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p1_`k'_ay_`outcome'_term.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_ay), ci(CI_est_ob_p2_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_ay), ci(CI_est_hb_p2_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p2_`k'_ay_`outcome'_term.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p1_1y), ci(CI_est_ob_p1_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_1y), ci(CI_est_hb_p1_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_p1_`k'_1y_`outcome'_term.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_1y), ci(CI_est_ob_p2_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_1y), ci(CI_est_hb_p2_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_p2_`k'_1y_`outcome'_term.pdf", as(pdf) replace

	clear matrix
	*}
}
}



***************
* G) Fuzzy RDD 
***************

* Round 8

* Note: For the definition of incumbent status: We opt for a simple solution with a 
*		common reference date instead of  year-specific reference date. We guess 
*		that this only marginally  affects the first stage and has no effect on 
* 		the reduced form. 

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0


foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {

*	qui {

	clear matrix

foreach k in tri {
	
	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {
	mat coef_`m' = J(1,20,.)
	mat CI_`m' = J(2,20,.)
	}
	local i = 1
	
	
	foreach var in F1 F2 F3 F4 F5 F6 F7 F8 {

		rdrobust `outcome'_`var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) fuzzy(inoffice_`var' sharpbw)

		* Bandwidth
		local bw_opt_p1_ay = e(h_l)
		local bw_half_p1_ay = `bw_opt_p1_ay'/2
		
		* Robust
		mat coef_est_ob_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `outcome'_`var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1_ay') rho(1) fuzzy(inoffice_`var' sharpbw) 
			
		* Conventional
		mat coef_est_hb_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_ay[2,`i'] = e(ci_r_cl)

		rdrobust `outcome'_`var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) fuzzy(inoffice_`var' sharpbw)

		* Bandwidth
		local bw_opt_p2_ay = e(h_l)
		local bw_half_p2_ay = `bw_opt_p2_ay'/2
		
		* Robust
		mat coef_est_ob_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_ay[2,`i'] = e(ci_r_rb)

		rdrobust `outcome'_`var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p2_ay') rho(1) fuzzy(inoffice_`var' sharpbw)
			
		* Conventional
		mat coef_est_hb_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_ay[2,`i'] = e(ci_r_cl)
		
		rdrobust `outcome'_`var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1) fuzzy(inoffice_`var' sharpbw)

		* Bandwidth
		local bw_opt_p1_1y = e(h_l)
		local bw_half_p1_1y = `bw_opt_p1_1y'/2
		
		* Robust
		mat coef_est_ob_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_1y[2,`i'] = e(ci_r_rb)


		rdrobust `outcome'_`var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p1_1y') rho(1) fuzzy(inoffice_`var' sharpbw)
			
		* Conventional
		mat coef_est_hb_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_1y[2,`i'] = e(ci_r_cl)

		rdrobust `outcome'_`var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1) fuzzy(inoffice_`var' sharpbw)

		* Bandwidth
		local bw_opt_p2_1y = e(h_l)
		local bw_half_p2_1y = `bw_opt_p2_1y'/2
		
		* Robust
		mat coef_est_ob_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_1y[2,`i'] = e(ci_r_rb)


		rdrobust `outcome'_`var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p2_1y') rho(1) ///
			fuzzy(inoffice_`var' sharpbw)
			
		* Conventional
		mat coef_est_hb_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_1y[2,`i'] = e(ci_r_cl)
		
		local ++i
		
		}
	local title: var label `outcome'

	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {

		mat colnames coef_`m' = "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
		mat colnames CI_`m' =  "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///  
			"+6 yrs" "+7 yrs" "+8 yrs"
	}
	
		coefplot ///
		(matrix(coef_est_ob_p1_ay), ci(CI_est_ob_p1_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_ay), ci(CI_est_hb_p1_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_frd_p1_`k'_ay_`outcome'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_ay), ci(CI_est_ob_p2_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_ay), ci(CI_est_hb_p2_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_frd_p2_`k'_ay_`outcome'.pdf", as(pdf) replace
		
				coefplot ///
		(matrix(coef_est_ob_p1_1y), ci(CI_est_ob_p1_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_1y), ci(CI_est_hb_p1_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_frd_p1_`k'_1y_`outcome'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_1y), ci(CI_est_ob_p2_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_1y), ci(CI_est_hb_p2_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_frd_p2_`k'_1y_`outcome'.pdf", as(pdf) replace

	clear matrix
	*}
}
}


************************
* H) RDD with covariates
************************

* Round 9

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

tab canton, gen(canton_)
tab year, gen(year_)
tostring year, gen(year_str)
gen ct_yr=canton+year_str
tab ct_yr, gen(ct_yr_)

local covs "ct_yr_*"

*local covs "canton_2 canton_3 canton_4 canton_5 canton_6 canton_7 canton_8 canton_9 canton_10 canton_11 canton_12 canton_13 canton_14 canton_15 canton_16 canton_17 canton_18 canton_19 canton_20 canton_21 canton_22 canton_23 canton_24 canton_25 canton_26 year_2 year_3 year_4 year_5 year_6 year_7 year_8 year_9 year_10 year_11 year_12 year_13 year_14 year_15 year_16 year_17 year_18 year_19 year_20 year_21 year_22"

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {

*	qui {
	clear matrix

foreach k in tri {
	
	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {
	mat coef_`m' = J(1,20,.)
	mat CI_`m' = J(2,20,.)
	}
	local i = 1

	foreach var of varlist `outcome'_L4 `outcome'_L3 `outcome'_L2 `outcome'_L1 ///
		`outcome' `outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4 ///
		`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8 {

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) covs(`covs')

		* Bandwidth
		local bw_opt_p1_ay = e(h_l)
		local bw_half_p1_ay = `bw_opt_p1_ay'/2
		
		* Robust
		mat coef_est_ob_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1_ay') rho(1) covs(`covs') 
			
		* Conventional
		mat coef_est_hb_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_ay[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) covs(`covs')

		* Bandwidth
		local bw_opt_p2_ay = e(h_l)
		local bw_half_p2_ay = `bw_opt_p2_ay'/2
		
		* Robust
		mat coef_est_ob_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p2_ay') rho(1) covs(`covs') 
			
		* Conventional
		mat coef_est_hb_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_ay[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1) covs(`covs')

		* Bandwidth
		local bw_opt_p1_1y = e(h_l)
		local bw_half_p1_1y = `bw_opt_p1_1y'/2
		
		* Robust
		mat coef_est_ob_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p1_1y') rho(1) covs(`covs') 
			
		* Conventional
		mat coef_est_hb_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_1y[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1) covs(`covs')

		* Bandwidth
		local bw_opt_p2_1y = e(h_l)
		local bw_half_p2_1y = `bw_opt_p2_1y'/2
		
		* Robust
		mat coef_est_ob_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p2_1y') rho(1) covs(`covs') 
			
		* Conventional
		mat coef_est_hb_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_1y[2,`i'] = e(ci_r_cl)
		
		local ++i
		
		}
	local title: var label `outcome'

	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {

		mat colnames coef_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" ///
			"0 yr" "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
		mat colnames CI_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" "0 yr" ///
			"+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" "+6 yrs" ///
			"+7 yrs" "+8 yrs"
	}
	
		coefplot ///
		(matrix(coef_est_ob_p1_ay), ci(CI_est_ob_p1_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_ay), ci(CI_est_hb_p1_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_cov_p1_`k'_ay_`outcome'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_ay), ci(CI_est_ob_p2_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_ay), ci(CI_est_hb_p2_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_cov_p2_`k'_ay_`outcome'.pdf", as(pdf) replace
		
				coefplot ///
		(matrix(coef_est_ob_p1_1y), ci(CI_est_ob_p1_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_1y), ci(CI_est_hb_p1_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_cov_p1_`k'_1y_`outcome'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_1y), ci(CI_est_ob_p2_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_1y), ci(CI_est_hb_p2_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_cov_p2_`k'_1y_`outcome'.pdf", as(pdf) replace

	clear matrix
	*}
}
}



**************************
* I) Fixed small bandwidth
**************************

* Round 10

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

* define fixed bandwidth
local fb 0.005

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {

*	qui {
	clear matrix

foreach k in tri {
	
	foreach m in ///
		est_fb_p1_ay ///
		est_fb_p2_ay ///
		est_fb_p1_1y ///
		est_fb_p2_1y {
	mat coef_`m' = J(1,20,.)
	mat CI_`m' = J(2,20,.)
	}
	local i = 1

	foreach var of varlist `outcome'_L4 `outcome'_L3 `outcome'_L2 `outcome'_L1 ///
		`outcome' `outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4 ///
		`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8 {


		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`fb') rho(1) 
			
		mat coef_est_fb_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_fb_p1_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_fb_p1_ay[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`fb') rho(1) 
			
		mat coef_est_fb_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_fb_p2_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_fb_p2_ay[2,`i'] = e(ci_r_cl)
		

		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`fb') rho(1) 
			
		mat coef_est_fb_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_fb_p1_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_fb_p1_1y[2,`i'] = e(ci_r_cl)


		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`fb') rho(1) 
			
		mat coef_est_fb_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_fb_p2_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_fb_p2_1y[2,`i'] = e(ci_r_cl)
		
		local ++i
		
		}
	local title: var label `outcome'

	foreach m in ///
		est_fb_p1_ay ///
		est_fb_p2_ay ///
		est_fb_p1_1y ///
		est_fb_p2_1y {

		mat colnames coef_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" ///
			"0 yr" "+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" ///
			"+6 yrs" "+7 yrs" "+8 yrs"
			
		mat colnames CI_`m' = "-4 yrs" "-3 yrs" "-2 yrs" "-1 yr" "0 yr" ///
			"+1 yr" "+2 yrs" "+3 yrs" "+4 yrs" "+5 yrs" "+6 yrs" ///
			"+7 yrs" "+8 yrs"
	}
	
		coefplot ///
		(matrix(coef_est_fb_p1_ay), ci(CI_est_fb_p1_ay) mcolor(blue) ciopts(lc(blue)) ), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_fb_p1_`k'_ay_`outcome'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_fb_p2_ay), ci(CI_est_fb_p2_ay) mcolor(blue) ciopts(lc(blue))), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_fb_p2_`k'_ay_`outcome'.pdf", as(pdf) replace
		
				coefplot ///
		(matrix(coef_est_fb_p1_1y), ci(CI_est_fb_p1_1y) mcolor(blue) ciopts(lc(blue))), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_fb_p1_`k'_1y_`outcome'.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_fb_p2_1y), ci(CI_est_fb_p2_1y) mcolor(blue) ciopts(lc(blue))), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_fb_p2_`k'_1y_`outcome'.pdf", as(pdf) replace

	clear matrix
	*}
}
}


*********************************************
* (J) Analysis at term level with covariates
*********************************************

* Note: round 11

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

tab canton, gen(canton_)
tab year, gen(year_)
tostring year, gen(year_str)
gen ct_yr=canton+year_str
tab ct_yr, gen(ct_yr_)

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

local covs "ct_yr_*"

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 {
	egen `outcome'_TL1=rowmax(`outcome'_L3 `outcome'_L2 `outcome'_L1 `outcome')
	egen `outcome'_T0=rowmax(`outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4)
	egen `outcome'_TF1=rowmax(`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8)
}
	
foreach outcome of varlist ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {
	egen `outcome'_TL1=rowmean(`outcome'_L3 `outcome'_L2 `outcome'_L1 `outcome')
	egen `outcome'_T0=rowmean(`outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4)
	egen `outcome'_TF1=rowmean(`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8)
}

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {

*	qui {
	clear matrix

foreach k in tri {
	
	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {
	mat coef_`m' = J(1,20,.)
	mat CI_`m' = J(2,20,.)
	}
	local i = 1

	foreach var of varlist `outcome'_TL1 `outcome'_T0 `outcome'_TF1 {

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) covs(`covs')

		* Bandwidth
		local bw_opt_p1_ay = e(h_l)
		local bw_half_p1_ay = `bw_opt_p1_ay'/2
		
		* Robust
		mat coef_est_ob_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1_ay') rho(1) covs(`covs')
			
		* Conventional
		mat coef_est_hb_p1_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_ay[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) covs(`covs')

		* Bandwidth
		local bw_opt_p2_ay = e(h_l)
		local bw_half_p2_ay = `bw_opt_p2_ay'/2
		
		* Robust
		mat coef_est_ob_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_ay[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_ay[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(2) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p2_ay') rho(1) covs(`covs')
			
		* Conventional
		mat coef_est_hb_p2_ay[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_ay[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_ay[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1) covs(`covs')

		* Bandwidth
		local bw_opt_p1_1y = e(h_l)
		local bw_half_p1_1y = `bw_opt_p1_1y'/2
		
		* Robust
		mat coef_est_ob_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p1_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p1_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(1) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p1_1y') rho(1) covs(`covs')
			
		* Conventional
		mat coef_est_hb_p1_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p1_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p1_1y[2,`i'] = e(ci_r_cl)

		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all rho(1) covs(`covs')

		* Bandwidth
		local bw_opt_p2_1y = e(h_l)
		local bw_half_p2_1y = `bw_opt_p2_1y'/2
		
		* Robust
		mat coef_est_ob_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_ob_p2_1y[1,`i'] = e(ci_l_rb)
		mat CI_est_ob_p2_1y[2,`i'] = e(ci_r_rb)


		rdrobust `var' votemargin_rel if first_yr == 1, p(2) bwselect(mserd) ///
			kernel(`k') vce(cluster ID_num) all h(`bw_half_p2_1y') rho(1) covs(`covs')
			
		* Conventional
		mat coef_est_hb_p2_1y[1,`i']=e(tau_cl)
		mat CI_est_hb_p2_1y[1,`i'] = e(ci_l_cl)
		mat CI_est_hb_p2_1y[2,`i'] = e(ci_r_cl)
		
		local ++i
		
		}
	local title: var label `outcome'

	foreach m in ///
		est_ob_p1_ay est_hb_p1_ay ///
		est_ob_p2_ay est_hb_p2_ay ///
		est_ob_p1_1y est_hb_p1_1y ///
		est_ob_p2_1y est_hb_p2_1y {

		mat colnames coef_`m' =  "-1 term" "0 term" "+1 term"
		mat colnames CI_`m' =  "-1 term" "0 term" "+1 term"
	}
	
		coefplot ///
		(matrix(coef_est_ob_p1_ay), ci(CI_est_ob_p1_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_ay), ci(CI_est_hb_p1_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p1_`k'_ay_`outcome'_term_cov.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_ay), ci(CI_est_ob_p2_ay) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_ay), ci(CI_est_hb_p2_ay) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p2_`k'_ay_`outcome'_term_cov.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p1_1y), ci(CI_est_ob_p1_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p1_1y), ci(CI_est_hb_p1_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_p1_`k'_1y_`outcome'_term_cov.pdf", as(pdf) replace
		
		coefplot ///
		(matrix(coef_est_ob_p2_1y), ci(CI_est_ob_p2_1y) mcolor(black) ciopts(lc(black)) offset(-0.1) ) ///
		(matrix(coef_est_hb_p2_1y), ci(CI_est_hb_p2_1y) mcolor(blue) ciopts(lc(blue)) offset(0.1)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "2nd order polynomial, `k', first election")
		graph export "$path_ol\figures\fig_p2_`k'_1y_`outcome'_term_cov.pdf", as(pdf) replace

	clear matrix
	*}
}
}



*********************************************
* (K) Combined plot
*********************************************

* Note: round 12

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

tostring year, gen(year_str)
gen ct_yr=canton+year_str
tab ct_yr, gen(ct_yr_)

local covs "ct_yr_*"
local fb 0.01

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 {
	egen `outcome'_TL1=rowmax(`outcome'_L3 `outcome'_L2 `outcome'_L1 `outcome')
	egen `outcome'_T0=rowmax(`outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4)
	egen `outcome'_TF1=rowmax(`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8)
}
	
foreach outcome of varlist ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {
	egen `outcome'_TL1=rowmean(`outcome'_L3 `outcome'_L2 `outcome'_L1 `outcome')
	egen `outcome'_T0=rowmean(`outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4)
	egen `outcome'_TF1=rowmean(`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8)
}

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	i_prs_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 ///
	n_prs_sum_c1 {


*	qui {
	clear matrix

foreach k in tri {
	
	foreach m in est_ob_p1 ///
				 est_hb_p1 est_hb_p1_cov ///
				 est_fb_p1 est_fb_p1_cov {
	mat coef_`m' = J(1,20,.)
	mat CI95_`m' = J(2,20,.)
	mat CI90_`m' = J(2,20,.)
	}
	local i = 1

	foreach var of varlist `outcome'_TL1 `outcome'_T0 `outcome'_TF1 {

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) level(95) 

		* Bandwidth
		local bw_opt_p1 = e(h_l)
		local bw_half_p1 = `bw_opt_p1'/2
		
		* Robust
		mat coef_est_ob_p1[1,`i']=e(tau_rb)
		mat CI95_est_ob_p1[1,`i'] = e(ci_l_rb)
		mat CI95_est_ob_p1[2,`i'] = e(ci_r_rb)
		
		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) level(90)
		
		mat CI90_est_ob_p1[1,`i'] = e(ci_l_rb)
		mat CI90_est_ob_p1[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1') rho(1) level(95) 
			
		* Conventional
		mat coef_est_hb_p1[1,`i']=e(tau_cl)
		mat CI95_est_hb_p1[1,`i'] = e(ci_l_cl)
		mat CI95_est_hb_p1[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1') rho(1) level(90)
				
		mat CI90_est_hb_p1[1,`i'] = e(ci_l_cl)
		mat CI90_est_hb_p1[2,`i'] = e(ci_r_cl)

		
		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) covs(`covs')

		* Bandwidth
		local bw_opt_p1 = e(h_l)
		local bw_half_p1 = `bw_opt_p1'/2
		
		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1') rho(1) covs(`covs') level(95)
			
		* Conventional
		mat coef_est_hb_p1_cov[1,`i']=e(tau_cl)
		mat CI95_est_hb_p1_cov[1,`i'] = e(ci_l_cl)
		mat CI95_est_hb_p1_cov[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1') rho(1) covs(`covs') level(90)
			
		mat CI90_est_hb_p1_cov[1,`i'] = e(ci_l_cl)
		mat CI90_est_hb_p1_cov[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) kernel(`k') ///
			vce(cluster ID_num) all h(`fb') rho(1) level(95)
			
		* Conventional
		mat coef_est_fb_p1[1,`i']=e(tau_cl)
		mat CI95_est_fb_p1[1,`i'] = e(ci_l_cl)
		mat CI95_est_fb_p1[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) kernel(`k') ///
			vce(cluster ID_num) all h(`fb') rho(1) level(90)	
		
		mat CI90_est_fb_p1[1,`i'] = e(ci_l_cl)
		mat CI90_est_fb_p1[2,`i'] =  e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) kernel(`k') ///
			vce(cluster ID_num) all h(`fb') rho(1) covs(`covs') level(95)
			
		* Conventional
		mat coef_est_fb_p1_cov[1,`i']=e(tau_cl)
		mat CI95_est_fb_p1_cov[1,`i'] = e(ci_l_cl)
		mat CI95_est_fb_p1_cov[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) kernel(`k') ///
			vce(cluster ID_num) all h(`fb') rho(1) covs(`covs') level(90)
			
		mat CI90_est_fb_p1_cov[1,`i'] = e(ci_l_cl)
		mat CI90_est_fb_p1_cov[2,`i'] = e(ci_r_cl)	
		local ++i
			}
				
	local title: var label `outcome'

	foreach m in est_ob_p1 ///
				 est_hb_p1 est_hb_p1_cov ///
				 est_fb_p1 est_fb_p1_cov {
		mat colnames coef_`m' =  "-1 term" "0 term" "+1 term"
		mat colnames CI95_`m' =  "-1 term" "0 term" "+1 term" 
		mat colnames CI90_`m' =  "-1 term" "0 term" "+1 term" 
		}
	
		coefplot ///
		(matrix(coef_est_ob_p1), ci(CI95_est_ob_p1) mcolor(black) ciopts(lc(black) lw(thin)) offset(-0.2)) ///
		(matrix(coef_est_ob_p1), ci(CI90_est_ob_p1) mcolor(black) ciopts(lc(black) lw(thick)) offset(-0.2)) ///		
		(matrix(coef_est_hb_p1), ci(CI95_est_hb_p1) mcolor(black) ciopts(lc(black) lw(thin)) offset(-0.1)) ///
		(matrix(coef_est_hb_p1), ci(CI90_est_hb_p1) mcolor(black) ciopts(lc(black) lw(thick)) offset(-0.1)) ///		
		(matrix(coef_est_hb_p1_cov), ci(CI95_est_hb_p1_cov) mcolor(black) ciopts(lc(black) lw(thin)) offset(0.0)) ///
		(matrix(coef_est_hb_p1_cov), ci(CI90_est_hb_p1_cov) mcolor(black) ciopts(lc(black) lw(thick)) offset(0.0)) ///		
		(matrix(coef_est_fb_p1), ci(CI95_est_fb_p1) mcolor(black) ciopts(lc(black) lw(thin)) offset(0.1)) ///
		(matrix(coef_est_fb_p1), ci(CI90_est_fb_p1) mcolor(black) ciopts(lc(black) lw(thick)) offset(0.1)) ///		
		(matrix(coef_est_fb_p1_cov), ci(CI95_est_fb_p1_cov) mcolor(black) ciopts(lc(black) lw(thin)) offset(0.2)) ///
		(matrix(coef_est_fb_p1_cov), ci(CI90_est_fb_p1_cov) mcolor(black) ciopts(lc(black) lw(thick)) offset(0.2)), ///							
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p1_`k'_`outcome'_term_combined.pdf", as(pdf) replace		
		
		coefplot ///
		(matrix(coef_est_ob_p1), ci(CI95_est_ob_p1) mcolor(black) ciopts(lc(black) lw(thin)) offset(-0.2)) ///
		(matrix(coef_est_ob_p1), ci(CI90_est_ob_p1) mcolor(black) ciopts(lc(black) lw(thick)) offset(-0.2)) ///		
		(matrix(coef_est_hb_p1), ci(CI95_est_hb_p1) mcolor(blue) ciopts(lc(blue) lw(thin)) offset(-0.1)) ///
		(matrix(coef_est_hb_p1), ci(CI90_est_hb_p1) mcolor(blue) ciopts(lc(blue) lw(thick)) offset(-0.1)) ///		
		(matrix(coef_est_hb_p1_cov), ci(CI95_est_hb_p1_cov) mcolor(blue) ciopts(lc(blue) lw(thin)) offset(0.0)) ///
		(matrix(coef_est_hb_p1_cov), ci(CI90_est_hb_p1_cov) mcolor(blue) ciopts(lc(blue) lw(thick)) offset(0.0)) ///		
		(matrix(coef_est_fb_p1), ci(CI95_est_fb_p1) mcolor(ebblue) ciopts(lc(ebblue) lw(thin)) offset(0.1)) ///
		(matrix(coef_est_fb_p1), ci(CI90_est_fb_p1) mcolor(ebblue) ciopts(lc(ebblue) lw(thick)) offset(0.1)) ///		
		(matrix(coef_est_fb_p1_cov), ci(CI95_est_fb_p1_cov) mcolor(ebblue) ciopts(lc(ebblue) lw(thin)) offset(0.2)) ///
		(matrix(coef_est_fb_p1_cov), ci(CI90_est_fb_p1_cov) mcolor(ebblue) ciopts(lc(ebblue) lw(thick)) offset(0.2)), ///							
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p1_`k'_`outcome'_term_combined_color.pdf", as(pdf) replace	
	clear matrix
}
}

*********************************************
* (L) Combined plot including population
*********************************************

* Note: round 13
* Use "rowmean" instead of "rowmax" for extensive margin variable. 
* Interpretation: not probability of having a mandate in at least one year of 
* the legislative period, but rather share of years with at least one mandate. 
* Reason: We cannot create an analogous variable for the general population 
* w/out having panel data for the general population.

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

tostring year, gen(year_str)
gen ct_yr=canton+year_str
tab ct_yr, gen(ct_yr_)

local covs "ct_yr_*"
local fb 0.01

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1  {
	egen `outcome'_TL1=rowmean(`outcome'_L3 `outcome'_L2 `outcome'_L1 `outcome')
	egen `outcome'_T0=rowmean(`outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4)
	egen `outcome'_TF1=rowmean(`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8)
}
	
foreach outcome of varlist ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1  {
	egen `outcome'_TL1=rowmean(`outcome'_L3 `outcome'_L2 `outcome'_L1 `outcome')
	egen `outcome'_T0=rowmean(`outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4)
	egen `outcome'_TF1=rowmean(`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8)
}

foreach outcome of varlist ///
	i_all_c1 ///
	i_lrg_c1 ///
	i_sml_c1 ///
	n_all_sum_c1 ///
	n_lrg_sum_c1 ///
	n_sml_sum_c1 {

	clear matrix	
	mat coef_cand = J(1,20,.)
	mat CI95_cand = J(2,20,.)
	mat CI90_cand = J(2,20,.)
	
	reg `outcome'_T0 elected, vce(cluster ID_num)
	
	mat coef_cand[1,2] = _b[elected]
	mat CI95_cand[1,2] = _b[elected] - invttail(e(df_r),0.025)*_se[elected]
	mat CI95_cand[2,2] = _b[elected] + invttail(e(df_r),0.025)*_se[elected]
	mat CI90_cand[1,2] = _b[elected] - invttail(e(df_r),0.05)*_se[elected]
	mat CI90_cand[2,2] = _b[elected] + invttail(e(df_r),0.05)*_se[elected]
	
foreach k in tri {
	
	foreach m in est_ob_p1 ///
				 est_hb_p1 est_hb_p1_cov ///
				 est_fb_p1 est_fb_p1_cov {
	mat coef_`m' = J(1,20,.)
	mat CI95_`m' = J(2,20,.)
	mat CI90_`m' = J(2,20,.)
	}
	local i = 1

	foreach var of varlist `outcome'_TL1 `outcome'_T0 `outcome'_TF1 {

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) level(95) 

		* Bandwidth
		local bw_opt_p1 = e(h_l)
		local bw_half_p1 = `bw_opt_p1'/2
		
		* Robust
		mat coef_est_ob_p1[1,`i'] = e(tau_rb)
		mat CI95_est_ob_p1[1,`i'] = e(ci_l_rb)
		mat CI95_est_ob_p1[2,`i'] = e(ci_r_rb)
		
		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) level(90)
		
		mat CI90_est_ob_p1[1,`i'] = e(ci_l_rb)
		mat CI90_est_ob_p1[2,`i'] = e(ci_r_rb)

		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1') rho(1) level(95) 
			
		* Conventional
		mat coef_est_hb_p1[1,`i'] = e(tau_cl)
		mat CI95_est_hb_p1[1,`i'] = e(ci_l_cl)
		mat CI95_est_hb_p1[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1') rho(1) level(90)
				
		mat CI90_est_hb_p1[1,`i'] = e(ci_l_cl)
		mat CI90_est_hb_p1[2,`i'] = e(ci_r_cl)

		
		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all rho(1) covs(`covs')

		* Bandwidth
		local bw_opt_p1 = e(h_l)
		local bw_half_p1 = `bw_opt_p1'/2
		
		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1') rho(1) covs(`covs') level(95)
			
		* Conventional
		mat coef_est_hb_p1_cov[1,`i'] = e(tau_cl)
		mat CI95_est_hb_p1_cov[1,`i'] = e(ci_l_cl)
		mat CI95_est_hb_p1_cov[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(`k') ///
			vce(cluster ID_num) all h(`bw_half_p1') rho(1) covs(`covs') level(90)
			
		mat CI90_est_hb_p1_cov[1,`i'] = e(ci_l_cl)
		mat CI90_est_hb_p1_cov[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) kernel(`k') ///
			vce(cluster ID_num) all h(`fb') rho(1) level(95)
			
		* Conventional
		mat coef_est_fb_p1[1,`i'] = e(tau_cl)
		mat CI95_est_fb_p1[1,`i'] = e(ci_l_cl)
		mat CI95_est_fb_p1[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) kernel(`k') ///
			vce(cluster ID_num) all h(`fb') rho(1) level(90)	
		
		mat CI90_est_fb_p1[1,`i'] = e(ci_l_cl)
		mat CI90_est_fb_p1[2,`i'] =  e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) kernel(`k') ///
			vce(cluster ID_num) all h(`fb') rho(1) covs(`covs') level(95)
			
		* Conventional
		mat coef_est_fb_p1_cov[1,`i'] = e(tau_cl)
		mat CI95_est_fb_p1_cov[1,`i'] = e(ci_l_cl)
		mat CI95_est_fb_p1_cov[2,`i'] = e(ci_r_cl)
		
		rdrobust `var' votemargin_rel, p(1) kernel(`k') ///
			vce(cluster ID_num) all h(`fb') rho(1) covs(`covs') level(90)
			
		mat CI90_est_fb_p1_cov[1,`i'] = e(ci_l_cl)
		mat CI90_est_fb_p1_cov[2,`i'] = e(ci_r_cl)	
		local ++i
			}
				
	local title: var label `outcome'

	foreach m in est_ob_p1 ///
				 est_hb_p1 est_hb_p1_cov ///
				 est_fb_p1 est_fb_p1_cov {
		mat colnames coef_`m' =  "-1 term" "0 term" "+1 term"
		mat colnames CI95_`m' =  "-1 term" "0 term" "+1 term" 
		mat colnames CI90_`m' =  "-1 term" "0 term" "+1 term" 
		}
	
	preserve
	use "$path\02_Processed_data\21_Population\population_data.dta", clear

	replace ellegperiod = 0 if ellegperiod == .
		
		mat coef_pop = J(1,20,.)
		mat CI95_pop = J(2,20,.)
		mat CI90_pop = J(2,20,.)
		
		reg `outcome' ellegperiod, robust
		
		mat coef_pop[1,2] = _b[ellegperiod]
		mat CI95_pop[1,2] = _b[ellegperiod] - invttail(e(df_r),0.025)*_se[ellegperiod]
		mat CI95_pop[2,2] = _b[ellegperiod] + invttail(e(df_r),0.025)*_se[ellegperiod]
		mat CI90_pop[1,2] = _b[ellegperiod] - invttail(e(df_r),0.05)*_se[ellegperiod]
		mat CI90_pop[2,2] = _b[ellegperiod] + invttail(e(df_r),0.05)*_se[ellegperiod]
	restore
	
		mat colnames coef_pop =  "-1 term" "0 term" "+1 term"
		mat colnames CI95_pop =  "-1 term" "0 term" "+1 term" 
		mat colnames CI90_pop =  "-1 term" "0 term" "+1 term" 
		mat colnames coef_cand =  "-1 term" "0 term" "+1 term"
		mat colnames CI95_cand =  "-1 term" "0 term" "+1 term" 
		mat colnames CI90_cand =  "-1 term" "0 term" "+1 term" 	
	
		coefplot ///
		(matrix(coef_est_ob_p1), ci(CI95_est_ob_p1) mcolor(black) ciopts(lc(black) lw(medthin)) offset(-0.2)) ///
		(matrix(coef_est_ob_p1), ci(CI90_est_ob_p1) mcolor(black) ciopts(lc(black) lw(medthick)) offset(-0.2)) ///		
		(matrix(coef_est_hb_p1), ci(CI95_est_hb_p1) mcolor(black) ciopts(lc(black) lw(medthin)) offset(-0.1)) ///
		(matrix(coef_est_hb_p1), ci(CI90_est_hb_p1) mcolor(black) ciopts(lc(black) lw(medthick)) offset(-0.1)) ///		
		(matrix(coef_est_hb_p1_cov), ci(CI95_est_hb_p1_cov) mcolor(black) ciopts(lc(black) lw(medthin)) offset(0.0)) ///
		(matrix(coef_est_hb_p1_cov), ci(CI90_est_hb_p1_cov) mcolor(black) ciopts(lc(black) lw(medthick)) offset(0.0)) ///		
		(matrix(coef_est_fb_p1), ci(CI95_est_fb_p1) mcolor(black) ciopts(lc(black) lw(medthin)) offset(0.1)) ///
		(matrix(coef_est_fb_p1), ci(CI90_est_fb_p1) mcolor(black) ciopts(lc(black) lw(medthick)) offset(0.1)) ///		
		(matrix(coef_est_fb_p1_cov), ci(CI95_est_fb_p1_cov) mcolor(black) ciopts(lc(black) lw(medthin)) offset(0.2)) ///
		(matrix(coef_est_fb_p1_cov), ci(CI90_est_fb_p1_cov) mcolor(black) ciopts(lc(black) lw(medthick)) offset(0.2)) ///		
		(matrix(coef_cand), ci(CI95_cand) mcolor(black) ciopts(lc(black) lw(medthin)) offset(-0.3)) ///
		(matrix(coef_cand), ci(CI90_cand) mcolor(black) ciopts(lc(black) lw(medthick)) offset(-0.3)) ///	
		(matrix(coef_pop), ci(CI95_pop) mcolor(black) ciopts(lc(black) lw(medthin)) offset(-0.4)) ///
		(matrix(coef_pop), ci(CI90_pop) mcolor(black) ciopts(lc(black) lw(medthick)) offset(-0.4)), ///	
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p1_`k'_`outcome'_term_pop.pdf", as(pdf) replace		
		
		coefplot ///
		(matrix(coef_est_ob_p1), ci(CI95_est_ob_p1) mcolor(black) ciopts(lc(black) lw(medthin)) offset(-0.2)) ///
		(matrix(coef_est_ob_p1), ci(CI90_est_ob_p1) mcolor(black) ciopts(lc(black) lw(medthick)) offset(-0.2)) ///		
		(matrix(coef_est_hb_p1), ci(CI95_est_hb_p1) mcolor(blue) ciopts(lc(blue) lw(medthin)) offset(-0.1)) ///
		(matrix(coef_est_hb_p1), ci(CI90_est_hb_p1) mcolor(blue) ciopts(lc(blue) lw(medthick)) offset(-0.1)) ///		
		(matrix(coef_est_hb_p1_cov), ci(CI95_est_hb_p1_cov) mcolor(blue) ciopts(lc(blue) lw(medthin)) offset(0.0)) ///
		(matrix(coef_est_hb_p1_cov), ci(CI90_est_hb_p1_cov) mcolor(blue) ciopts(lc(blue) lw(medthick)) offset(0.0)) ///		
		(matrix(coef_est_fb_p1), ci(CI95_est_fb_p1) mcolor(ebblue) ciopts(lc(ebblue) lw(medthin)) offset(0.1)) ///
		(matrix(coef_est_fb_p1), ci(CI90_est_fb_p1) mcolor(ebblue) ciopts(lc(ebblue) lw(medthick)) offset(0.1)) ///		
		(matrix(coef_est_fb_p1_cov), ci(CI95_est_fb_p1_cov) mcolor(ebblue) ciopts(lc(ebblue) lw(medthin)) offset(0.2)) ///
		(matrix(coef_est_fb_p1_cov), ci(CI90_est_fb_p1_cov) mcolor(ebblue) ciopts(lc(ebblue) lw(medthick)) offset(0.2)) ///
		(matrix(coef_cand), ci(CI95_cand) mcolor(lavender) ciopts(lc(lavender) lw(medthin)) offset(-0.3)) ///
		(matrix(coef_cand), ci(CI90_cand) mcolor(lavender) ciopts(lc(lavender) lw(medthick)) offset(-0.3)) ///
		(matrix(coef_pop), ci(CI95_pop) mcolor(purple) ciopts(lc(purple) lw(medthin)) offset(-0.4)) ///
		(matrix(coef_pop), ci(CI90_pop) mcolor(purple) ciopts(lc(purple) lw(medthick)) offset(-0.4)), ///
		vertical nolabel yline(0) ylabel(, angle(horizontal) gsty(dot)) ///
		graphregion(fcolor(white) lcolor(white)) legend(off) ///
		title("`outcome': `title'" "1st order polynomial, `k', all elections")
		graph export "$path_ol\figures\fig_p1_`k'_`outcome'_term_pop_color.pdf", as(pdf) replace	
	clear matrix
}
}



***************
* Z) RDD plot
***************


* Note: For the definition of incumbent status: We opt for a simple solution with a 
*		common reference date instead of  year-specific reference date. We guess 
*		that this only marginally  affects the first stage and has no effect on 
* 		the reduced form. 

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

bysort ID: egen min_yr=min(year)
g first_yr = 0
replace first_yr = 1 if year == min_yr & cand_before1931 == 0

local outcome i_all_c1_L4

rdrobust `outcome' votemargin_rel if first_yr==1, p(1) bwselect(mserd) kernel(tri) ///
	vce(cluster ID_num) all rho(1)
	
local bw_opt_p1_1y = e(h_l)
local bw_half_p1_1y = `bw_opt_p1_1y'/2

cap rdplot `outcome' votemargin_rel  ///
	if -`bw_half_p1_1y' <= votemargin_rel & votemargin_rel <= `bw_half_p1_1y' & first_yr==1, ///
	p(1) kernel(triangular) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 