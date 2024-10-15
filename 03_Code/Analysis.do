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
*global path_ol "C:\Schmidlu\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"
global path_ol "C:\Schmidlu\Dropbox\Apps\Overleaf\Political_Rents"

global path "C:\Current\Dropbox\Projekt Nationalräte"
global path_ol "C:\Current\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"

global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"
*global path_ol "E:\12. Cloud\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"
global path_ol "E:\12. Cloud\Dropbox\Apps\Overleaf\Political_Rents"


*********************************************
* (A) Data preparation and macros
*********************************************

use "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", clear
keep if votemargin_rel < .

tostring year, gen(year_str)
gen ct_yr=canton+year_str
tab ct_yr, gen(ct_yr_)

foreach outcome of varlist ///
	n_all_sum_c1 i_all_c1 ///
	n_lrg_sum_c1 i_lrg_c1 ///
	n_sml_sum_c1 i_sml_c1 { 
	egen `outcome'_TL1=rowmean(`outcome'_L3 `outcome'_L2 `outcome'_L1 `outcome')
	egen `outcome'_T0=rowmean(`outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4)
	egen `outcome'_TF1=rowmean(`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8)
}

global covs "ct_yr_*"
global fb 0.01

*********************************************
* (B) Baseline figure
*********************************************

foreach outcome of varlist ///
	n_all_sum_c1 i_all_c1 ///
	n_lrg_sum_c1 i_lrg_c1 ///
	n_sml_sum_c1 i_sml_c1 {  

	clear matrix	
	mat coef_cand = J(1,20,.)
	mat CI95_cand = J(2,20,.)
	mat CI90_cand = J(2,20,.)

	reg `outcome'_T0 elected, vce(cluster ID_num)
	* Control mean: sum `outcome'_T0  if elected == 0
		
	mat coef_cand[1,2] = _b[elected]
	mat CI95_cand[1,2] = _b[elected] - invttail(e(df_r),0.025)*_se[elected]
	mat CI95_cand[2,2] = _b[elected] + invttail(e(df_r),0.025)*_se[elected]
	mat CI90_cand[1,2] = _b[elected] - invttail(e(df_r),0.05)*_se[elected]
	mat CI90_cand[2,2] = _b[elected] + invttail(e(df_r),0.05)*_se[elected]
			
	foreach m in est_ob_p1 ///
		est_hb_p1 est_hb_p1_cov ///
		est_fb_p1 est_fb_p1_cov {
		mat coef_`m' = J(1,20,.)
		mat CI95_`m' = J(2,20,.)
		mat CI90_`m' = J(2,20,.)
		}
		local i = 1

		foreach var of varlist `outcome'_TL1  `outcome'_T0 `outcome'_TF1 { 

			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all rho(1) level(95)
				
			* Bandwidth
			local bw_opt_p1 = e(h_l)
			local bw_half_p1 = `bw_opt_p1'/2

			* Robust
			mat coef_est_ob_p1[1,`i'] = e(tau_rb)
			mat CI95_est_ob_p1[1,`i'] = e(ci_l_rb)
			mat CI95_est_ob_p1[2,`i'] = e(ci_r_rb)			

			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all rho(1) level(90)

			mat CI90_est_ob_p1[1,`i'] = e(ci_l_rb)
			mat CI90_est_ob_p1[2,`i'] = e(ci_r_rb)

			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all h(`bw_half_p1') rho(1) level(95) 

			* Conventional
			mat coef_est_hb_p1[1,`i'] = e(tau_cl)
			mat CI95_est_hb_p1[1,`i'] = e(ci_l_cl)
			mat CI95_est_hb_p1[2,`i'] = e(ci_r_cl)

			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all h(`bw_half_p1') rho(1) level(90)
					
			mat CI90_est_hb_p1[1,`i'] = e(ci_l_cl)
			mat CI90_est_hb_p1[2,`i'] = e(ci_r_cl)

			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all rho(1) covs($covs)

			* Bandwidth
			local bw_opt_p1 = e(h_l)
			local bw_half_p1 = `bw_opt_p1'/2
			
			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all h(`bw_half_p1') rho(1) covs($covs) level(95)
				
			* Conventional
			mat coef_est_hb_p1_cov[1,`i'] = e(tau_cl)
			mat CI95_est_hb_p1_cov[1,`i'] = e(ci_l_cl)
			mat CI95_est_hb_p1_cov[2,`i'] = e(ci_r_cl)

			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all h(`bw_half_p1') rho(1) covs($covs) level(90)

			mat CI90_est_hb_p1_cov[1,`i'] = e(ci_l_cl)
			mat CI90_est_hb_p1_cov[2,`i'] = e(ci_r_cl)

			rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
				vce(cluster ID_num) all h($fb) rho(1) level(95)

			* Conventional
			mat coef_est_fb_p1[1,`i'] = e(tau_cl)
			mat CI95_est_fb_p1[1,`i'] = e(ci_l_cl)
			mat CI95_est_fb_p1[2,`i'] = e(ci_r_cl)

			rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
				vce(cluster ID_num) all h($fb) rho(1) level(90)	

			mat CI90_est_fb_p1[1,`i'] = e(ci_l_cl)
			mat CI90_est_fb_p1[2,`i'] =  e(ci_r_cl)
			
			rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
				vce(cluster ID_num) all h($fb) rho(1) covs($covs) level(95)

			* Conventional
			mat coef_est_fb_p1_cov[1,`i'] = e(tau_cl)
			mat CI95_est_fb_p1_cov[1,`i'] = e(ci_l_cl)
			mat CI95_est_fb_p1_cov[2,`i'] = e(ci_r_cl)

			rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
				vce(cluster ID_num) all h($fb) rho(1) covs($covs) level(90)

			mat CI90_est_fb_p1_cov[1,`i'] = e(ci_l_cl)
			mat CI90_est_fb_p1_cov[2,`i'] = e(ci_r_cl)	
			local ++i
		}
		
		
		* Export graphs		

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
			* Control mean: sum `outcome'  if ellegperiod == 0

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
			graphregion(fcolor(white) lcolor(white)) legend(off)  title("")
			graph export "$path_ol\figures\fig_p1_`k'_`outcome'_term_pop_color.pdf", as(pdf) replace

		clear matrix

}

*********************************************
* (C) Baseline and robustness tables
*********************************************

cap erase  "$path_ol/tables/app_baseline.tex"
gen votemargin = votemargin_rel
gen n_all_TL1 = n_all_sum_c1_TL1
gen i_all_TL1 = i_all_c1_TL1
gen n_lrg_TL1 = n_lrg_sum_c1_TL1
gen i_lrg_TL1 = i_lrg_c1_TL1
gen n_sml_TL1 = n_sml_sum_c1_TL1
gen i_sml_TL1 = i_sml_c1_TL1
gen n_all_T0 = n_all_sum_c1_T0
gen i_all_T0 = i_all_c1_T0
gen n_lrg_T0 = n_lrg_sum_c1_T0
gen i_lrg_T0 = i_lrg_c1_T0
gen n_sml_T0 = n_sml_sum_c1_T0
gen i_sml_T0 = i_sml_c1_T0
gen n_all_TF1 = n_all_sum_c1_TF1
gen i_all_TF1 = i_all_c1_TF1
gen n_lrg_TF1 = n_lrg_sum_c1_TF1
gen i_lrg_TF1 = i_lrg_c1_TF1
gen n_sml_TF1 = n_sml_sum_c1_TF1
gen i_sml_TF1 = i_sml_c1_TF1

clear programs

program tables 

replace votemargin = votemargin_rel + `1'

local trimmed = `2'

if `trimmed' == 1 {

	foreach outcome in ///
		n_all ///
		n_lrg ///
		n_sml {
		sum `outcome'_TL1, d
		replace `outcome'_TL1 = . if `outcome'_TL1 > `r(p99)'
		sum `outcome'_T0, d
		replace `outcome'_T0 = . if `outcome'_T0 > `r(p99)'
		sum `outcome'_TF1, d
		replace `outcome'_TF1 = . if `outcome'_TF1 > `r(p99)'
		}
}

foreach outcome in ///
	n_all i_all ///
	n_lrg i_lrg ///
	n_sml i_sml {
		
	reg `outcome'_T0 elected, vce(cluster ID_num)
	* Control mean: sum `outcome'_T0  if elected == 0

	foreach k in tri {
		
		foreach m in est_ob_p1 ///
					 est_hb_p1 est_hb_p1_cov ///
					 est_fb_p1 est_fb_p1_cov {
		}
		local i = 1

		foreach var of varlist `outcome'_TL1  `outcome'_T0 `outcome'_TF1 { 
			
			local e_first = (`i'-1)*5+1
			local e_second = (`i'-1)*5+2
			local e_third = (`i'-1)*5+3
			local e_fourth = (`i'-1)*5+4
			local e_fifth = (`i'-1)*5+5

			rdrobust `var' votemargin, p(`3') bwselect(mserd) kernel(`4') ///
				vce(cluster ID_num) all rho(1) level(95)
				
			* Bandwidth
			local bw_opt_p1 = e(h_l)
			local bw_half_p1 = `bw_opt_p1'/2
			
			* Robust
			local ci_l_rb=string(e(ci_l_rb), "%12.2f")
			local ci_r_rb=string(e(ci_r_rb), "%12.2f")
			estadd local estimate ="[`ci_l_rb',`ci_r_rb']"
			estadd local se ="-"
			estadd scalar control_mean=e(tau_bc_l)
			estadd local covariates ="No"
			estadd local bw_type ="Optimal"
			estadd scalar bw = e(h_l)
			estadd scalar po = e(p)
			estadd local ke = "Tri."
			estadd scalar nobs_l = e(N_h_l)
			estadd scalar nobs_r = e(N_h_r)
			est store e`e_first'

			rdrobust `var' votemargin, p(`3') bwselect(mserd) kernel(`4') ///
				vce(cluster ID_num) all h(`bw_half_p1') rho(1) level(95) 

			* Conventional
			local estimate = string(e(tau_cl), "%12.3f")
			if e(pv_cl)>=0.1 estadd local estimate ="`estimate'"
			if e(pv_cl)>=0.05 & e(pv_cl)<0.1 estadd local estimate ="`estimate'\sym{*}"
			if e(pv_cl)>=0.01 & e(pv_cl)<0.05 estadd local estimate ="`estimate'\sym{**}"
			if e(pv_cl)<0.01 estadd local estimate ="`estimate'\sym{***}"
			local se =string(e(se_tau_cl), "%12.3f") 
			estadd local se= "(`se')" 
			estadd scalar control_mean=e(tau_bc_l)
			estadd local covariates ="No"
			estadd local bw_type ="Half"
			estadd scalar bw = e(h_l)
			estadd scalar po = e(p)
			estadd local ke = "Tri."
			estadd scalar nobs_l = e(N_h_l)
			estadd scalar nobs_r = e(N_h_r)
			est store e`e_second'

			rdrobust `var' votemargin, p(`3') bwselect(mserd) kernel(`4') ///
				vce(cluster ID_num) all rho(1) covs($covs)

			* Bandwidth
			local bw_opt_p1 = e(h_l)
			local bw_half_p1 = `bw_opt_p1'/2

			rdrobust `var' votemargin, p(`3') bwselect(mserd) kernel(`4') ///
				vce(cluster ID_num) all h(`bw_half_p1') rho(1) covs($covs) level(95)

			* Conventional
			local estimate = string(e(tau_cl), "%12.3f")
			if e(pv_cl)>=0.1 estadd local estimate ="`estimate'"
			if e(pv_cl)>=0.05 & e(pv_cl)<0.1 estadd local estimate ="`estimate'\sym{*}"
			if e(pv_cl)>=0.01 & e(pv_cl)<0.05 estadd local estimate ="`estimate'\sym{**}"
			if e(pv_cl)<0.01 estadd local estimate ="`estimate'\sym{***}"
			local se =string(e(se_tau_cl), "%12.3f") 
			estadd local se= "(`se')" 
			estadd scalar control_mean=e(tau_bc_l)
			estadd local covariates ="Yes"
			estadd local bw_type ="Half"
			estadd scalar bw = e(h_l)
			estadd scalar po = e(p)
			estadd local ke = "Tri."
			estadd scalar nobs_l = e(N_h_l)
			estadd scalar nobs_r = e(N_h_r)
			est store e`e_third'
			
			rdrobust `var' votemargin, p(`3') kernel(`4') ///
				vce(cluster ID_num) all h($fb) rho(1) level(95)

			* Conventional
			local estimate = string(e(tau_cl), "%12.3f")
			if e(pv_cl)>=0.1 estadd local estimate ="`estimate'"
			if e(pv_cl)>=0.05 & e(pv_cl)<0.1 estadd local estimate ="`estimate'\sym{*}"
			if e(pv_cl)>=0.01 & e(pv_cl)<0.05 estadd local estimate ="`estimate'\sym{**}"
			if e(pv_cl)<0.01 estadd local estimate ="`estimate'\sym{***}"
			local se =string(e(se_tau_cl), "%12.3f") 
			estadd local se= "(`se')" 
			estadd scalar control_mean=e(tau_bc_l)
			estadd local covariates ="No"
			estadd local bw_type ="Fixed"
			estadd scalar bw = e(h_l)
			estadd scalar po = e(p)
			estadd local ke = "Tri."
			estadd scalar nobs_l = e(N_h_l)
			estadd scalar nobs_r = e(N_h_r)
			est store e`e_fourth'
			
			rdrobust `var' votemargin, p(`3') kernel(`4') ///
				vce(cluster ID_num) all h($fb) rho(1) covs($covs) level(95)

			* Conventional
			local estimate = string(e(tau_cl), "%12.3f")
			if e(pv_cl)>=0.1 estadd local estimate ="`estimate'"
			if e(pv_cl)>=0.05 & e(pv_cl)<0.1 estadd local estimate ="`estimate'\sym{*}"
			if e(pv_cl)>=0.01 & e(pv_cl)<0.05 estadd local estimate ="`estimate'\sym{**}"
			if e(pv_cl)<0.01 estadd local estimate ="`estimate'\sym{***}"
			local se =string(e(se_tau_cl), "%12.3f") 
			estadd local se= "(`se')" 
			estadd scalar control_mean=e(tau_bc_l)
			estadd local covariates ="Yes"
			estadd local bw_type ="Fixed"
			estadd scalar bw = e(h_l)
			estadd scalar po = e(p)
			estadd local ke = "Tri."
			estadd scalar nobs_l = e(N_h_l)
			estadd scalar nobs_r = e(N_h_r)
			est store e`e_fifth'

			local ++i
				}
		* Export tables
		
		if "`outcome'"=="i_all" | "`outcome'"=="i_lrg" | "`outcome'"=="i_sml" {
			
		esttab e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 e11 e12 e13 e14 e15 ///
		using "$path_ol/tables/app_plac_`1'_trim_`2'_poly_`3'_kern_`4'.tex", ///
		varwidth(21) fragment nogap noli nomti nonum nonotes append drop(Conventional Bias-corrected Robust) ///
		stats(estimate se control_mean bw nobs_l nobs_r covariates bw_type po ke, ///
		fmt(%12.3f %12.3f %12.3f %12.3f %12.0f %12.0f %12.0f %12.0f %12.0f %12.0f) ///
			labels("Estimate" "Standard error" "Control mean" "Bandwidth" ///
			"Observations left" "Observations right" "Covariates" ///
			"Bandwidth type" "Polynomial order" "Kernel type" ))
			
		} 
		else{
			
		esttab e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 e11 e12 e13 e14 e15 ///
		using "$path_ol/tables/app_plac_`1'_trim_`2'_poly_`3'_kern_`4'.tex", ///
		varwidth(21) fragment nogap noli nomti nonum nonotes append drop(Conventional Bias-corrected Robust) ///
		stats(estimate se control_mean bw nobs_l nobs_r, ///
		fmt(%12.3f %12.3f %12.3f %12.3f %12.0f %12.0f) ///
			labels("Estimate" "Standard error" "Control mean" "Bandwidth" ///
			"Observations left" "Observations right"))	
			
		}

	}
	}
end

tables	0		0		1		tri
tables	0		0		2		tri
tables	0		0		1		uni
tables	-0.06	0		1		tri
tables	0.06	0		1		tri
tables	0		1		1		tri

!rename "$path_ol/tables/app_plac_0_trim_0_poly_1_kern_tri.tex" "$path_ol/tables/app_baseline.tex"




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