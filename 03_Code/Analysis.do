********************************************************************************
*                                                                              *
*           Politicians as Directors - Evidence from Close Elections           *
*                                                                              *
********************************************************************************

clear
set more off
version 17

*** Set paths

*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "C:\Current\Dropbox\Projekt Nationalräte"
*path_ol "C:\Current\Dropbox\Apps\Overleaf\Political_Rents"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
global path_ol "E:\12. Cloud\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"

/*** Estimation variant 1
* Note: Both variants produce the same results.

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear


local outcome i_lrg_c1 // prc_lk_s1 rcl_lk_s1 f1_lk_s1 prc_mn_s1 rcl_mn_s1 f1_mn_s1 

foreach var of varlist votemargin_rel_F4 votemargin_rel_F3 votemargin_rel_F2 ///
	votemargin_rel_F1 votemargin_rel votemargin_rel_L1 votemargin_rel_L2 ///
	votemargin_rel_L3 votemargin_rel_L4 {

	rdrobust `outcome' `var', p(1) bwselect(mserd) kernel(triangular) ///
		vce(cluster ID_num) all rho(1)
}
*/

*** Estimation variant 2
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