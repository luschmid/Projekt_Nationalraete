********************************************************************************
*                                                                              *
*           Politicians as Directors - Evidence from Close Elections           *
*                                                                              *
********************************************************************************

clear
set more off
version 18

*** Set paths

global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path_ol "C:\Schmidlu\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"
global path_ol "C:\Schmidlu\Dropbox\Apps\Overleaf\Political_Rents"

*global path "C:\Current\Dropbox\Projekt Nationalräte"
*global path_ol "C:\Current\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"

global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"
*global path_ol "E:\12. Cloud\Dropbox\Projekt Nationalräte\04_Results\04_Political_Rents"
global path_ol "E:\12. Cloud\Dropbox\Apps\Overleaf\Political_Rents"

* ssc install sumstats


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
	n_sml_sum_c1 i_sml_c1 ///
	n_six_sum_c1 i_six_c1 { 
	egen `outcome'_TL1=rowmean(`outcome'_L3 `outcome'_L2 `outcome'_L1 `outcome')
	egen `outcome'_T0=rowmean(`outcome'_F1 `outcome'_F2 `outcome'_F3 `outcome'_F4)
	egen `outcome'_TF1=rowmean(`outcome'_F5 `outcome'_F6 `outcome'_F7 `outcome'_F8)
}

global covs "ct_yr_*"
global fb 0.01

gen age=year-birthyear
sum elected votemargin_rel i_all_c1_T0 n_all_sum_c1_T0 age sex

** Destats

preserve
keep if (i_all_c1_T0 != . | i_all_c1_TL1  != . |  i_all_c1_TF1 != .)  & votemargin_rel!= . 
unique ID 

sumstats (elected votemargin_rel i_all_c1_T0 i_lrg_c1_T0 i_sml_c1_T0 ///
	n_all_sum_c1_T0 n_lrg_sum_c1_T0 n_sml_sum_c1_T0) ///
    using "$path_ol/tables/destats.xlsx", replace stats(mean sd min max n)
restore 


*********************************************
* (B) Baseline figure
*********************************************

foreach outcome of varlist ///
	n_all_sum_c1 i_all_c1 ///
	n_lrg_sum_c1 i_lrg_c1 ///
	n_sml_sum_c1 i_sml_c1  {  // n_six_sum_c1 i_six_c1

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

			* Bias-corrected robust
			
			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all rho(1) level(95)
				
			local bw_opt_p1 = e(h_l)
			local bw_half_p1 = `bw_opt_p1'/2

			mat coef_est_ob_p1[1,`i'] = e(tau_rb)
			mat CI95_est_ob_p1[1,`i'] = e(ci_l_rb)
			mat CI95_est_ob_p1[2,`i'] = e(ci_r_rb)			

			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all rho(1) level(90)

			mat CI90_est_ob_p1[1,`i'] = e(ci_l_rb)
			mat CI90_est_ob_p1[2,`i'] = e(ci_r_rb)
			

			* Conventional with no covariates
			
			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all h(`bw_half_p1') rho(1) level(95) 

			mat coef_est_hb_p1[1,`i'] = e(tau_cl)
			mat CI95_est_hb_p1[1,`i'] = e(ci_l_cl)
			mat CI95_est_hb_p1[2,`i'] = e(ci_r_cl)

			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all h(`bw_half_p1') rho(1) level(90)
					
			mat CI90_est_hb_p1[1,`i'] = e(ci_l_cl)
			mat CI90_est_hb_p1[2,`i'] = e(ci_r_cl)

			* Conventional with covariates
			
			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all rho(1) covs($covs)

			local bw_opt_p1 = e(h_l)
			local bw_half_p1 = `bw_opt_p1'/2
						
			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all h(`bw_half_p1') rho(1) covs($covs) level(95)	
			
			mat coef_est_hb_p1_cov[1,`i'] = e(tau_cl)
			mat CI95_est_hb_p1_cov[1,`i'] = e(ci_l_cl)
			mat CI95_est_hb_p1_cov[2,`i'] = e(ci_r_cl)

			rdrobust `var' votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
				vce(cluster ID_num) all h(`bw_half_p1') rho(1) covs($covs) level(90)

			mat CI90_est_hb_p1_cov[1,`i'] = e(ci_l_cl)
			mat CI90_est_hb_p1_cov[2,`i'] = e(ci_r_cl)

			* Conventional with fixed bandwidth
			
			rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
				vce(cluster ID_num) all h($fb) rho(1) level(95)

			mat coef_est_fb_p1[1,`i'] = e(tau_cl)
			mat CI95_est_fb_p1[1,`i'] = e(ci_l_cl)
			mat CI95_est_fb_p1[2,`i'] = e(ci_r_cl)

			rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
				vce(cluster ID_num) all h($fb) rho(1) level(90)	

			mat CI90_est_fb_p1[1,`i'] = e(ci_l_cl)
			mat CI90_est_fb_p1[2,`i'] =  e(ci_r_cl)
			
			rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
				vce(cluster ID_num) all h($fb) rho(1) covs($covs) level(95)

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
		
		gen candidate = 0
		replace candidate = 1 if ellegperiod < .
		replace ellegperiod = 0 if ellegperiod == .
		gen weight = 4/count
			
			mat coef_pop = J(1,20,.)
			mat CI95_pop = J(2,20,.)
			mat CI90_pop = J(2,20,.)
			
			reg `outcome' ellegperiod [iweight=weight], robust  // elected against population, identical results with aweight
			* Control mean: sum `outcome' [iweight=weight] if ellegperiod == 0
			
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
*gen n_six_TF1 = n_six_sum_c1_TF1
*gen i_six_TF1 = i_six_c1_TF1

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
	i_all n_all ///
	i_lrg n_lrg ///
	i_sml n_sml {  // i_six n_six
		
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
			estadd local ke = "`4'"
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
			estadd local ke = "`4'"
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
			estadd local ke = "`4'"
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
			estadd local ke = "`4'"
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
			estadd local ke = "`4'"
			estadd scalar nobs_l = e(N_h_l)
			estadd scalar nobs_r = e(N_h_r)
			est store e`e_fifth'

			local ++i
				}
		* Export tables
		
		if "`outcome'"=="n_all" | "`outcome'"=="n_lrg" | "`outcome'"=="n_sml" | "`outcome'"=="n_six" {
			
		esttab e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 e11 e12 e13 e14 e15 ///
		using "$path_ol/tables/app_plac_`1'_trim_`2'_poly_`3'_kern_`4'`5'.tex", ///
		varwidth(21) fragment nogap noli nomti nonum nonotes append drop(Conventional Bias-corrected Robust) ///
		stats(estimate se control_mean bw nobs_l nobs_r covariates bw_type po ke, ///
		fmt(%12.3f %12.3f %12.3f %12.3f %12.0f %12.0f %12.0f %12.0f %12.0f %12.0f) ///
			labels("Estimate" "Standard error" "Control mean" "Bandwidth" ///
			"Observations left" "Observations right" "Covariates" ///
			"Bandwidth type" "Polynomial order" "Kernel type" ))
			
		} 
		else{
			
		esttab e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 e11 e12 e13 e14 e15 ///
		using "$path_ol/tables/app_plac_`1'_trim_`2'_poly_`3'_kern_`4'`5'.tex", ///
		varwidth(21) fragment nogap noli nomti nonum nonotes append drop(Conventional Bias-corrected Robust) ///
		stats(estimate se control_mean bw nobs_l nobs_r, ///
		fmt(%12.3f %12.3f %12.3f %12.3f %12.0f %12.0f) ///
			labels("Estimate" "Standard error" "Control mean" "Bandwidth" ///
			"Observations left" "Observations right"))	
			
		}

	}
	}
end

cap erase  "$path_ol/tables/app_plac_0_trim_0_poly_1_kern_uni.tex"
cap erase  "$path_ol/tables/app_plac_0_trim_0_poly_2_kern_tri.tex"

* Main results
tables	0		0		1		tri

* Robustness
tables	0		0		2		tri
tables	0		0		1		uni

* Heterogeneity
preserve
keep if party_left_wing==1
tables	0		0		1		tri 	_left
restore

preserve
keep if party_center==1
tables	0		0		1		tri 	_center
restore

preserve
keep if party_right_wing==1
tables	0		0		1		tri 	_right
restore

preserve
keep if sex==1
tables	0		0		1		tri 	_men
restore

preserve
keep if year < 1971
tables	0		0		1		tri 	_early
restore

preserve
keep if year >= 1971
tables	0		0		1		tri 	_late
restore


tables	-0.06	0		1		tri
tables	0.06	0		1		tri
tables	0		1		1		tri
tables	0		1		1		tri

!rename "$path_ol/tables/app_plac_0_trim_0_poly_1_kern_tri.tex" "$path_ol/tables/app_baseline.tex"


*******************
* (D) Balance tests
*******************

* (i) Define and create covariates
* Note: we use the four-year lag for votemargin_rel (votemargin_rel_L4) 
* 	    because previous election took place four years before. 

global covs_bal votemargin_rel_L4 ///
	party_left_wing party_center party_right_wing ///
	age sex ///
	year no_seats /// 
	canton_1 canton_2 canton_3 canton_4 canton_5 canton_6 canton_7 canton_8 ///
	canton_9 canton_10 canton_11 canton_12 canton_13 canton_14 canton_15 ///
	canton_16 canton_17 canton_18 canton_19 canton_20 canton_21 canton_22 ///
	canton_23 canton_24 canton_25 canton_26  
	
tab canton, gen(canton_)
label var votemargin_rel_L4  "Previous vote margin"
label var age   "Age"
label var canton_1   "Aargau"				
label var canton_2   "Appenzell Innerrhoden"
label var canton_3   "Appenzell Ausserrhoden"
label var canton_4   "Bern"
label var canton_5   "Basel Landschaft"
label var canton_6   "Basel Stadt"
label var canton_7   "Fribourg"
label var canton_8   "Geneva"
label var canton_9   "Glarus"
label var canton_10  "Graubünden"
label var canton_11  "Jura"
label var canton_12  "Lucerne"
label var canton_13  "Neuchâtel"
label var canton_14  "Nidwalden"
label var canton_15  "Obwalden"
label var canton_16  "St. Gallen"
label var canton_17  "Schaffhausen"
label var canton_18  "Solothurn"
label var canton_19  "Schwyz"
label var canton_20  "Thurgau"
label var canton_21  "Ticino"
label var canton_22  "Uri"
label var canton_23  "Vaud"
label var canton_24  "Valais"
label var canton_25  "Zug"
label var canton_26  "Zürich"

* (ii) Read out bandwidths for both outcome variables "at least one directorship"
* 	   and "number of directorships"

rdrobust i_all_c1_T0 votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
		vce(cluster ID_num) all rho(1) level(95)
		
local bw_opt_i = e(h_l)
local bw_half_i = `bw_opt_i'/2

rdrobust n_all_sum_c1_T0 votemargin_rel, p(1) bwselect(mserd) kernel(tri) ///
		vce(cluster ID_num) all rho(1) level(95)
		
local bw_opt_n = e(h_l)
local bw_half_n = `bw_opt_n'/2

* (iii) Set up rejection rates

local i_ob_5  = 0
local i_hb_5  = 0
local n_ob_5  = 0
local n_hb_5  = 0
local fb_5    = 0
local i_ob_10 = 0
local i_hb_10 = 0
local n_ob_10 = 0
local n_hb_10 = 0
local fb_10   = 0
local i_ob_sum  = 0
local i_hb_sum  = 0
local n_ob_sum  = 0
local n_hb_sum  = 0
local fb_sum 	  = 0

* (iv) Loop over all covariates

capture erase "$path_ol/tables/balance_table.tex"
capture erase "$path_ol/tables/balance_table_rr.tex"

foreach var of varlist $covs_bal {   
	
* (a) Robust for at least one directorship
rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
	vce(cluster ID_num) all rho(1) level(95) h(`bw_opt_i')

if e(pv_rb)< .{ // output for those variables with variation within bandwidth
	local ci_l_rb=string(e(ci_l_rb), "%12.2f")
	local ci_r_rb=string(e(ci_r_rb), "%12.2f")
	estadd local estimate ="[`ci_l_rb',`ci_r_rb']"
	if e(pv_rb)<0.1 local ++i_ob_10
	if e(pv_rb)<0.05 local ++i_ob_5
	local ++i_ob_sum
	estadd local se =""
} 
else { // output for those variables without variation within bandwidth
	estadd local estimate ="-"
	estadd local se =""
}

estadd scalar bw = e(h_l)
estadd scalar nobs_l = e(N_h_l)
estadd scalar nobs_r = e(N_h_r)
est store e1	

* (b) Conventional for at least one directorship
rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
	vce(cluster ID_num) all h(`bw_half_i') rho(1) level(95) 
	
if e(pv_cl)< .{ // output for those variables with variation within bandwidth
	local estimate = string(e(tau_cl), "%12.3f")
	if e(pv_cl)>=0.1 estadd local estimate ="`estimate'"
	if e(pv_cl)>=0.05 & e(pv_cl)<0.1 estadd local estimate ="`estimate'\sym{*}"
	if e(pv_cl)>=0.01 & e(pv_cl)<0.05 estadd local estimate ="`estimate'\sym{**}"
	if e(pv_cl)<0.01 estadd local estimate ="`estimate'\sym{***}"
	if e(pv_cl)<0.1 local ++i_hb_10
	if e(pv_cl)<0.05 local ++i_hb_5
	local ++i_hb_sum
	local se =string(e(se_tau_cl), "%12.3f") 
	estadd local se= "(`se')" 
} 
else { // output for those variables without variation within bandwidth
	estadd local estimate ="-"
	estadd local se =""
}	

estadd scalar bw = e(h_l)
estadd scalar nobs_l = e(N_h_l)
estadd scalar nobs_r = e(N_h_r)
est store e2

* (c) Robust for number of directorships
rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
	vce(cluster ID_num) all rho(1) level(95) h(`bw_opt_n')
	
if e(pv_rb)< .{ // output for those variables with variation within bandwidth
	local ci_l_rb=string(e(ci_l_rb), "%12.2f")
	local ci_r_rb=string(e(ci_r_rb), "%12.2f")
	estadd local estimate ="[`ci_l_rb',`ci_r_rb']"
	if e(pv_rb)<0.1 local ++n_ob_10
	if e(pv_rb)<0.05 local ++n_ob_5
	local ++n_ob_sum
	estadd local se =""
} 
else { // output for those variables without variation within bandwidth
	estadd local estimate ="-"
	estadd local se =""
}	
	
estadd scalar bw = e(h_l)
estadd scalar nobs_l = e(N_h_l)
estadd scalar nobs_r = e(N_h_r)
est store e3	

* (d) Conventional for number of directorships
rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
	vce(cluster ID_num) all h(`bw_half_n') rho(1) level(95) 
	
if e(pv_cl)< .{ // output for those variables with variation within bandwidth
	local estimate = string(e(tau_cl), "%12.3f")
	if e(pv_cl)>=0.1 estadd local estimate ="`estimate'"
	if e(pv_cl)>=0.05 & e(pv_cl)<0.1 estadd local estimate ="`estimate'\sym{*}"
	if e(pv_cl)>=0.01 & e(pv_cl)<0.05 estadd local estimate ="`estimate'\sym{**}"
	if e(pv_cl)<0.01 estadd local estimate ="`estimate'\sym{***}"
	if e(pv_cl)<0.1 local ++n_hb_10
	if e(pv_cl)<0.05 local ++n_hb_5
	local ++n_hb_sum
	local se =string(e(se_tau_cl), "%12.3f") 
	estadd local se= "(`se')" 
} 
else { // output for those variables without variation within bandwidth
	estadd local estimate ="-"
	estadd local se =""
}	

estadd scalar bw = e(h_l)
estadd scalar nobs_l = e(N_h_l)
estadd scalar nobs_r = e(N_h_r)
est store e4

* (e) Conventional for fixed bandwidth
rdrobust `var' votemargin_rel, p(1) kernel(tri) ///
	vce(cluster ID_num) all h(0.01) rho(1) level(95) 
	
if e(pv_cl)< .{ // output for those variables with variation within bandwidth
	local estimate = string(e(tau_cl), "%12.3f")
	if e(pv_cl)>=0.1 estadd local estimate ="`estimate'"
	if e(pv_cl)>=0.05 & e(pv_cl)<0.1 estadd local estimate ="`estimate'\sym{*}"
	if e(pv_cl)>=0.01 & e(pv_cl)<0.05 estadd local estimate ="`estimate'\sym{**}"
	if e(pv_cl)<0.01 estadd local estimate ="`estimate'\sym{***}"
	if e(pv_cl)<0.1 local ++fb_10
	if e(pv_cl)<0.05 local ++fb_5
	local ++fb_sum
	local se =string(e(se_tau_cl), "%12.3f") 
	estadd local se= "(`se')" 
} 
else { // output for those variables without variation within bandwidth
	estadd local estimate ="-"
	estadd local se =""
}

estadd scalar bw = e(h_l)
estadd scalar nobs_l = e(N_h_l)
estadd scalar nobs_r = e(N_h_r)
est store e5

* (v) Export tables

local var_label: variable label `var'

if "`var'"=="canton_26" | {
	
esttab e1 e2 e3 e4 e5  ///
using "$path_ol/tables/balance_table.tex", ///
varwidth(21) fragment nogap noli nomti nonum nonotes append drop(Conventional Bias-corrected Robust) ///
stats(estimate se bw nobs_l nobs_r, ///
fmt(%12.3f %12.3f %12.3f %12.0f %12.0f) ///
	labels("`var_label'" " " "Bandwidth" ///
	"Observations left" "Observations right"))
	
} 
else{
	
esttab e1 e2 e3 e4 e5  ///
using "$path_ol/tables/balance_table.tex", ///
varwidth(21) fragment nogap noli nomti nonum nonotes append drop(Conventional Bias-corrected Robust) ///
stats(estimate se , ///
fmt(%12.3f %12.3f) ///
	labels("`var_label'" " " ))	
	
}
}

* (vi) Export rejection rates

mat rrates = J(2,5,.)
mat rrates[1,1] =`i_ob_5'/`i_ob_sum'*100
mat rrates[2,1] =`i_ob_10'/`i_ob_sum'*100
mat rrates[1,2] =`i_hb_5'/`i_hb_sum'*100
mat rrates[2,2] =`i_hb_10'/`i_hb_sum'*100
mat rrates[1,3] =`n_ob_5'/`n_ob_sum'*100
mat rrates[2,3] =`n_ob_10'/`n_ob_sum'*100
mat rrates[1,4] =`n_hb_5'/`n_hb_sum'*100
mat rrates[2,4] =`n_hb_10'/`n_hb_sum'*100
mat rrates[1,5] =`fb_5'/`fb_sum'*100
mat rrates[2,5] =`fb_10'/`fb_sum'*100

mat rownames rrates = "5\%" "10\%"
mat colnames rrates = "\hspace{0.1cm}" "\hspace{0.2cm}" "\hspace{0.3cm}" "\hspace{0.4cm}" "\hspace{0.5cm}"
esttab matrix(rrates) using "$path_ol/tables/balance_table_rr.tex",  ///
	nomti nonum append fragment



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