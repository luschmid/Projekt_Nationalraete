******************************************************
* (A) Set path
******************************************************

clear
set more off
version 17

*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
global path "C:\Current\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

******************************************************
* (B) Read in and merge data for periods 1934-2003 
*     (source 1) and 1994-2018 (source 2)
******************************************************

use "$path\02_Processed_data\10_Directors_1934_2003\Politicians_Directorships_1934-2003.dta", clear

merge 1:1 ID year canton name firstname birthyear sex job elected votemargin_rel ///
	incumbent tenure list listname_bfs cand_before1931 municipality ///
	municipality_num language using ///
	"$path\02_Processed_data\11_Directors_1994_2018\Politicians_Directorships_1994-2017.dta"
	
* Quality checks
* _merge == 3: 266,290 obs = 10 * 26,629
sum year if _merge == 1 // Only years before 1994
sum year if _merge == 2 // Only years after 2003
sum year if _merge == 3 // Only years 1994-2003
table year ///
	if inlist(year, 1931, 1935, 1939, 1943, 1947, 1951, 1955, 1959, 1963, 1967) | ///
	inlist(year, 1971, 1975, 1979, 1983, 1987, 1991, 1995, 1999, 2003, 2007) | ///
	inlist(year, 2011, 2015), statistic(total elected) nototals
// Election information complete? Yes, discrepancies are tacit elections.
table year ///
	if inlist(year, 1931, 1935, 1939, 1943, 1947, 1951, 1955, 1959, 1963, 1967) | ///
	inlist(year, 1971, 1975, 1979, 1983, 1987, 1991, 1995, 1999, 2003, 2007) | ///
	inlist(year, 2011, 2015), statistic(count votemargin_rel) nototals
duplicates report ID year
drop _merge

* Comparison of data sources
foreach var in n_all_sum n_lrg_sum n_sml_sum n_prs_sum ///
	n_all_avg n_lrg_avg n_sml_avg n_prs_avg ///
	i_all i_lrg i_sml i_prs {
	qui pwcorr `var'_s1 `var'_s2
	di "`var'" _col(20) %9.2f `r(rho)'
}

sum n_all_sum_s1 n_all_sum_s2 if inrange(year,1994,2003)
sum n_lrg_sum_s1 n_lrg_sum_s2 if inrange(year,1994,2003)
sum n_sml_sum_s1 n_sml_sum_s2 if inrange(year,1994,2003)
sum n_prs_sum_s1 n_prs_sum_s2 if inrange(year,1994,2003)

sum n_all_avg_s1 n_all_avg_s2 if inrange(year,1994,2003)
sum n_lrg_avg_s1 n_lrg_avg_s2 if inrange(year,1994,2003)
sum n_sml_avg_s1 n_sml_avg_s2 if inrange(year,1994,2003)
sum n_prs_avg_s1 n_prs_avg_s2 if inrange(year,1994,2003)

sum i_all_s1 i_all_s2 if inrange(year,1994,2003)
sum i_lrg_s1 i_lrg_s2 if inrange(year,1994,2003)
sum i_sml_s1 i_sml_s2 if inrange(year,1994,2003)
sum i_prs_s1 i_prs_s2 if inrange(year,1994,2003)

foreach var in n_all_sum n_lrg_sum n_sml_sum n_prs_sum ///
	n_all_avg n_lrg_avg n_sml_avg n_prs_avg ///
	i_all i_lrg i_sml i_prs {
		gen diff_`var' = `var'_s1 - `var'_s2
}

sum diff_*
* The large differences in the between -92 and 90 (all firms) and -33 and 30
* (large firms) are a mistake in the Bisnode data: They mixed up Alfred Heer
* from ZH with Alfred Heer from GL.

/*
br ID year canton name firstname elected n_all_sum* if inrange(year,1994,2003)
br ID year canton name firstname elected n_all_sum* if (ID == "LU-2003-0035" | ID == "ZH-1979-0051") & inrange(year,1994,2003)
*/

foreach var in n_all_sum n_lrg_sum n_sml_sum n_prs_sum n_all_avg n_lrg_avg ///
	n_sml_avg n_prs_avg i_all i_lrg i_sml i_prs {
	gen `var'_c1 = `var'_s1
	replace `var'_c1 = `var'_s2 if `var'_s1 == . & `var'_s2 != .
}

foreach var in n_all_sum n_lrg_sum n_sml_sum n_prs_sum n_all_avg n_lrg_avg ///
	n_sml_avg n_prs_avg i_all i_lrg i_sml i_prs {
	gen `var'_c2 = `var'_s2
	replace `var'_c2 = `var'_s1 if `var'_s2 == . & `var'_s1 != .
}

label var i_all_c1 "At least one mandate, all companies"
label var i_lrg_c1 "At least one mandate, large companies"
label var i_sml_c1 "At least one mandate, small companies"
label var i_prs_c1 "At least one mandate as president"
label var n_all_sum_c1 "Number of mandates, all companies"
label var n_all_avg_c1 "Avg. number of mandates, all companies"
label var n_lrg_sum_c1 "Number of mandates, large companies"
label var n_lrg_avg_c1 "Avg. number of mandates, large companies"
label var n_sml_sum_c1 "Number of mandates, small companies"
label var n_sml_avg_c1 "Avg. number of mandates, small companies"
label var n_prs_sum_c1 "Number of mandates as president"
label var n_prs_avg_c1 "Avg. number of mandates as president"
label var i_all_c2 "At least one mandate, all companies"
label var i_lrg_c2 "At least one mandate, large companies"
label var i_sml_c2 "At least one mandate, small companies"
label var i_prs_c2 "At least one mandate as president"
label var n_all_sum_c2 "Number of mandates, all companies"
label var n_all_avg_c2 "Avg. number of mandates, all companies"
label var n_lrg_sum_c2 "Number of mandates, large companies"
label var n_lrg_avg_c2 "Avg. number of mandates, large companies"
label var n_sml_sum_c2 "Number of mandates, small companies"
label var n_sml_avg_c2 "Avg. number of mandates, small companies"
label var n_prs_sum_c2 "Number of mandates as president"
label var n_prs_avg_c2 "Avg. number of mandates as president"
*label var prc_lk "Precision at link-level"
*label var rcl_lk "Recall at link-level"
*label var f1_lk "F1 at link-level"
*label var prc_mn "Precision at mandate-level"
*label var rcl_mn "Recakll at mandate-level"
*label var f1_mn "F1 at mandate-level"


******************************************************
* (C) Create lags and leads
******************************************************

encode ID, gen(ID_num)
xtset ID_num year

foreach var of varlist elected votemargin_rel n_* i_* prc_* rcl_* ///
	f1_* {
	forv i = 1(1)8 {
	gen `var'_L`i' = L`i'.`var'
	gen `var'_F`i' = F`i'.`var'
}
}

order ID_num ID-job tenure-dir_year  ///
	elected_F8 elected_F7 elected_F6 elected_F5 elected_F4 elected_F3 elected_F2 ///
	elected_F1 elected elected_L1 elected_L2 elected_L3 elected_L4 elected_L5 ///
	elected_L6 elected_L7 elected_L8 ///
	votemargin_rel_F8 votemargin_rel_F7 votemargin_rel_F6 votemargin_rel_F5 ///
	votemargin_rel_F4 votemargin_rel_F3 votemargin_rel_F2 votemargin_rel_F1 ///
	votemargin_rel votemargin_rel_L1 votemargin_rel_L2 votemargin_rel_L3 ///
	votemargin_rel_L4 votemargin_rel_L5 votemargin_rel_L6 votemargin_rel_L7 ///
	votemargin_rel_L8 ///
	incumbent_F8 incumbent_F7 incumbent_F6 incumbent_F5 incumbent_F4 ///
	incumbent_F3 incumbent_F2 incumbent_F1 incumbent incumbent_L1 incumbent_L2 ///
	incumbent_L3 incumbent_L4 incumbent_L5 incumbent_L6 incumbent_L7 ///
	incumbent_L8 ///
	n_all_sum_c1_F8 n_all_sum_c1_F7 n_all_sum_c1_F6 n_all_sum_c1_F5 ///
	n_all_sum_c1_F4 n_all_sum_c1_F3 n_all_sum_c1_F2 n_all_sum_c1_F1 n_all_sum_c1 ///
	n_all_sum_c1_L1 n_all_sum_c1_L2 n_all_sum_c1_L3 n_all_sum_c1_L4 ///
	n_all_sum_c1_L5 n_all_sum_c1_L6 n_all_sum_c1_L7 n_all_sum_c1_L8 ///
	n_all_sum_c2_F8 n_all_sum_c2_F7 n_all_sum_c2_F6 n_all_sum_c2_F5 ///
	n_all_sum_c2_F4 n_all_sum_c2_F3 n_all_sum_c2_F2 n_all_sum_c2_F1 n_all_sum_c1 ///
	n_all_sum_c2_L1 n_all_sum_c2_L2 n_all_sum_c2_L3 n_all_sum_c2_L4 ///
	n_all_sum_c2_L5 n_all_sum_c2_L6 n_all_sum_c2_L7 n_all_sum_c2_L8 ///	
	n_lrg_sum_c1_F8 n_lrg_sum_c1_F7 n_lrg_sum_c1_F6 n_lrg_sum_c1_F5 ///
	n_lrg_sum_c1_F4 n_lrg_sum_c1_F3 n_lrg_sum_c1_F2 n_lrg_sum_c1_F1 n_lrg_sum_c1 ///
	n_lrg_sum_c1_L1 n_lrg_sum_c1_L2 n_lrg_sum_c1_L3 n_lrg_sum_c1_L4 ///
	n_lrg_sum_c1_L5 n_lrg_sum_c1_L6 n_lrg_sum_c1_L7 n_lrg_sum_c1_L8 ///
	n_lrg_sum_c2_F8 n_lrg_sum_c2_F7 n_lrg_sum_c2_F6 n_lrg_sum_c2_F5 ///
	n_lrg_sum_c2_F4 n_lrg_sum_c2_F3 n_lrg_sum_c2_F2 n_lrg_sum_c2_F1 n_lrg_sum_c2 ///
	n_lrg_sum_c2_L1 n_lrg_sum_c2_L2 n_lrg_sum_c2_L3 n_lrg_sum_c2_L4 ///
	n_lrg_sum_c2_L5 n_lrg_sum_c2_L6 n_lrg_sum_c2_L7 n_lrg_sum_c2_L8 ///
	n_sml_sum_c1_F8 n_sml_sum_c1_F7 n_sml_sum_c1_F6 n_sml_sum_c1_F5 ///
	n_sml_sum_c1_F4 n_sml_sum_c1_F3 n_sml_sum_c1_F2 n_sml_sum_c1_F1 n_sml_sum_c1 ///
	n_sml_sum_c1_L1 n_sml_sum_c1_L2 n_sml_sum_c1_L3 n_sml_sum_c1_L4 ///
	n_sml_sum_c1_L5 n_sml_sum_c1_L6 n_sml_sum_c1_L7 n_sml_sum_c1_L8 ///
	n_sml_sum_c2_F8 n_sml_sum_c2_F7 n_sml_sum_c2_F6 n_sml_sum_c2_F5 ///
	n_sml_sum_c2_F4 n_sml_sum_c2_F3 n_sml_sum_c2_F2 n_sml_sum_c2_F1 n_sml_sum_c2 ///
	n_sml_sum_c2_L1 n_sml_sum_c2_L2 n_sml_sum_c2_L3 n_sml_sum_c2_L4 ///
	n_sml_sum_c2_L5 n_sml_sum_c2_L6 n_sml_sum_c2_L7 n_sml_sum_c2_L8 ///
	n_prs_sum_c1_F8 n_prs_sum_c1_F7 n_prs_sum_c1_F6 n_prs_sum_c1_F5 ///
	n_prs_sum_c1_F4 n_prs_sum_c1_F3 n_prs_sum_c1_F2 n_prs_sum_c1_F1 n_prs_sum_c1 ///
	n_prs_sum_c1_L1 n_prs_sum_c1_L2 n_prs_sum_c1_L3 n_prs_sum_c1_L4 ///
	n_prs_sum_c1_L5 n_prs_sum_c1_L6 n_prs_sum_c1_L7 n_prs_sum_c1_L8 ///
	n_prs_sum_c2_F8 n_prs_sum_c2_F7 n_prs_sum_c2_F6 n_prs_sum_c2_F5 ///
	n_prs_sum_c2_F4 n_prs_sum_c2_F3 n_prs_sum_c2_F2 n_prs_sum_c2_F1 n_prs_sum_c2 ///
	n_prs_sum_c2_L1 n_prs_sum_c2_L2 n_prs_sum_c2_L3 n_prs_sum_c2_L4 ///
	n_prs_sum_c2_L5 n_prs_sum_c2_L6 n_prs_sum_c2_L7 n_prs_sum_c2_L8 ///	
	n_all_avg_c1_F8 n_all_avg_c1_F7 n_all_avg_c1_F6 n_all_avg_c1_F5 ///
	n_all_avg_c1_F4 n_all_avg_c1_F3 n_all_avg_c1_F2 n_all_avg_c1_F1 n_all_avg_c1 ///
	n_all_avg_c1_L1 n_all_avg_c1_L2 n_all_avg_c1_L3 n_all_avg_c1_L4 ///
	n_all_avg_c1_L5 n_all_avg_c1_L6 n_all_avg_c1_L7 n_all_avg_c1_L8 ///
	n_all_avg_c2_F8 n_all_avg_c2_F7 n_all_avg_c2_F6 n_all_avg_c2_F5 ///
	n_all_avg_c2_F4 n_all_avg_c2_F3 n_all_avg_c2_F2 n_all_avg_c2_F1 n_all_avg_c1 ///
	n_all_avg_c2_L1 n_all_avg_c2_L2 n_all_avg_c2_L3 n_all_avg_c2_L4 ///
	n_all_avg_c2_L5 n_all_avg_c2_L6 n_all_avg_c2_L7 n_all_avg_c2_L8 ///	
	n_lrg_avg_c1_F8 n_lrg_avg_c1_F7 n_lrg_avg_c1_F6 n_lrg_avg_c1_F5 ///
	n_lrg_avg_c1_F4 n_lrg_avg_c1_F3 n_lrg_avg_c1_F2 n_lrg_avg_c1_F1 n_lrg_avg_c1 ///
	n_lrg_avg_c1_L1 n_lrg_avg_c1_L2 n_lrg_avg_c1_L3 n_lrg_avg_c1_L4 ///
	n_lrg_avg_c1_L5 n_lrg_avg_c1_L6 n_lrg_avg_c1_L7 n_lrg_avg_c1_L8 ///
	n_lrg_avg_c2_F8 n_lrg_avg_c2_F7 n_lrg_avg_c2_F6 n_lrg_avg_c2_F5 ///
	n_lrg_avg_c2_F4 n_lrg_avg_c2_F3 n_lrg_avg_c2_F2 n_lrg_avg_c2_F1 n_lrg_avg_c2 ///
	n_lrg_avg_c2_L1 n_lrg_avg_c2_L2 n_lrg_avg_c2_L3 n_lrg_avg_c2_L4 ///
	n_lrg_avg_c2_L5 n_lrg_avg_c2_L6 n_lrg_avg_c2_L7 n_lrg_avg_c2_L8 ///
	n_sml_avg_c1_F8 n_sml_avg_c1_F7 n_sml_avg_c1_F6 n_sml_avg_c1_F5 ///
	n_sml_avg_c1_F4 n_sml_avg_c1_F3 n_sml_avg_c1_F2 n_sml_avg_c1_F1 n_sml_avg_c1 ///
	n_sml_avg_c1_L1 n_sml_avg_c1_L2 n_sml_avg_c1_L3 n_sml_avg_c1_L4 ///
	n_sml_avg_c1_L5 n_sml_avg_c1_L6 n_sml_avg_c1_L7 n_sml_avg_c1_L8 ///
	n_sml_avg_c2_F8 n_sml_avg_c2_F7 n_sml_avg_c2_F6 n_sml_avg_c2_F5 ///
	n_sml_avg_c2_F4 n_sml_avg_c2_F3 n_sml_avg_c2_F2 n_sml_avg_c2_F1 n_sml_avg_c2 ///
	n_sml_avg_c2_L1 n_sml_avg_c2_L2 n_sml_avg_c2_L3 n_sml_avg_c2_L4 ///
	n_sml_avg_c2_L5 n_sml_avg_c2_L6 n_sml_avg_c2_L7 n_sml_avg_c2_L8 ///
	n_prs_avg_c1_F8 n_prs_avg_c1_F7 n_prs_avg_c1_F6 n_prs_avg_c1_F5 ///
	n_prs_avg_c1_F4 n_prs_avg_c1_F3 n_prs_avg_c1_F2 n_prs_avg_c1_F1 n_prs_avg_c1 ///
	n_prs_avg_c1_L1 n_prs_avg_c1_L2 n_prs_avg_c1_L3 n_prs_avg_c1_L4 ///
	n_prs_avg_c1_L5 n_prs_avg_c1_L6 n_prs_avg_c1_L7 n_prs_avg_c1_L8 ///
	n_prs_avg_c2_F8 n_prs_avg_c2_F7 n_prs_avg_c2_F6 n_prs_avg_c2_F5 ///
	n_prs_avg_c2_F4 n_prs_avg_c2_F3 n_prs_avg_c2_F2 n_prs_avg_c2_F1 n_prs_avg_c2 ///
	n_prs_avg_c2_L1 n_prs_avg_c2_L2 n_prs_avg_c2_L3 n_prs_avg_c2_L4 ///
	n_prs_avg_c2_L5 n_prs_avg_c2_L6 n_prs_avg_c2_L7 n_prs_avg_c2_L8 ///
	i_all_c1_F8 i_all_c1_F7 i_all_c1_F6 i_all_c1_F5 i_all_c1_F4 i_all_c1_F3 ///
	i_all_c1_F2 i_all_c1_F1 i_all_c1 i_all_c1_L1 i_all_c1_L2 i_all_c1_L3 ///
	i_all_c1_L4 i_all_c1_L5 i_all_c1_L6 i_all_c1_L7 i_all_c1_L8 ///
	i_all_c2_F8 i_all_c2_F7 i_all_c2_F6 i_all_c2_F5 i_all_c2_F4 i_all_c2_F3 ///
	i_all_c2_F2 i_all_c2_F1 i_all_c2 i_all_c2_L1 i_all_c2_L2 i_all_c2_L3 ///
	i_all_c2_L4 i_all_c2_L5 i_all_c2_L6 i_all_c2_L7 i_all_c2_L8 ///
	i_lrg_c1_F8 i_lrg_c1_F7 i_lrg_c1_F6 i_lrg_c1_F5 i_lrg_c1_F4 i_lrg_c1_F3 ///
	i_lrg_c1_F2 i_lrg_c1_F1 i_lrg_c1 i_lrg_c1_L1 i_lrg_c1_L2 i_lrg_c1_L3 ///
	i_lrg_c1_L4 i_lrg_c1_L5 i_lrg_c1_L6 i_lrg_c1_L7 i_lrg_c1_L8 ///
	i_lrg_c2_F8 i_lrg_c2_F7 i_lrg_c2_F6 i_lrg_c2_F5 i_lrg_c2_F4 i_lrg_c2_F3 ///
	i_lrg_c2_F2 i_lrg_c2_F1 i_lrg_c2 i_lrg_c2_L1 i_lrg_c2_L2 i_lrg_c2_L3 ///
	i_lrg_c2_L4 i_lrg_c2_L5 i_lrg_c2_L6 i_lrg_c2_L7 i_lrg_c2_L8 ///
	i_sml_c1_F8 i_sml_c1_F7 i_sml_c1_F6 i_sml_c1_F5 i_sml_c1_F4 i_sml_c1_F3 ///
	i_sml_c1_F2 i_sml_c1_F1 i_sml_c1 i_sml_c1_L1 i_sml_c1_L2 i_sml_c1_L3 ///
	i_sml_c1_L4 i_sml_c1_L5 i_sml_c1_L6 i_sml_c1_L7 i_sml_c1_L8 ///
	i_sml_c2_F8 i_sml_c2_F7 i_sml_c2_F6 i_sml_c2_F5 i_sml_c2_F4 i_sml_c2_F3 ///
	i_sml_c2_F2 i_sml_c2_F1 i_sml_c2 i_sml_c2_L1 i_sml_c2_L2 i_sml_c2_L3 ///
	i_sml_c2_L4 i_sml_c2_L5 i_sml_c2_L6 i_sml_c2_L7 i_sml_c2_L8 ///	
	i_prs_c1_F8 i_prs_c1_F7 i_prs_c1_F6 i_prs_c1_F5 i_prs_c1_F4 i_prs_c1_F3 ///
	i_prs_c1_F2 i_prs_c1_F1 i_prs_c1 i_prs_c1_L1 i_prs_c1_L2 i_prs_c1_L3 ///
	i_prs_c1_L4 i_prs_c1_L5 i_prs_c1_L6 i_prs_c1_L7 i_prs_c1_L8 ///
	i_prs_c2_F8 i_prs_c2_F7 i_prs_c2_F6 i_prs_c2_F5 i_prs_c2_F4 i_prs_c2_F3 ///
	i_prs_c2_F2 i_prs_c2_F1 i_prs_c2 i_prs_c2_L1 i_prs_c2_L2 i_prs_c2_L3 ///
	i_prs_c2_L4 i_prs_c2_L5 i_prs_c2_L6 i_prs_c2_L7 i_prs_c2_L8 ///	
	n_all_sum_s1_F8 n_all_sum_s1_F7 n_all_sum_s1_F6 n_all_sum_s1_F5 ///
	n_all_sum_s1_F4 n_all_sum_s1_F3 n_all_sum_s1_F2 n_all_sum_s1_F1 n_all_sum_s1 ///
	n_all_sum_s1_L1 n_all_sum_s1_L2 n_all_sum_s1_L3 n_all_sum_s1_L4 ///
	n_all_sum_s1_L5 n_all_sum_s1_L6 n_all_sum_s1_L7 n_all_sum_s1_L8 ///
	n_all_sum_s2_F8 n_all_sum_s2_F7 n_all_sum_s2_F6 n_all_sum_s2_F5 ///
	n_all_sum_s2_F4 n_all_sum_s2_F3 n_all_sum_s2_F2 n_all_sum_s2_F1 n_all_sum_s1 ///
	n_all_sum_s2_L1 n_all_sum_s2_L2 n_all_sum_s2_L3 n_all_sum_s2_L4 ///
	n_all_sum_s2_L5 n_all_sum_s2_L6 n_all_sum_s2_L7 n_all_sum_s2_L8 ///	
	n_lrg_sum_s1_F8 n_lrg_sum_s1_F7 n_lrg_sum_s1_F6 n_lrg_sum_s1_F5 ///
	n_lrg_sum_s1_F4 n_lrg_sum_s1_F3 n_lrg_sum_s1_F2 n_lrg_sum_s1_F1 n_lrg_sum_s1 ///
	n_lrg_sum_s1_L1 n_lrg_sum_s1_L2 n_lrg_sum_s1_L3 n_lrg_sum_s1_L4 ///
	n_lrg_sum_s1_L5 n_lrg_sum_s1_L6 n_lrg_sum_s1_L7 n_lrg_sum_s1_L8 ///
	n_lrg_sum_s2_F8 n_lrg_sum_s2_F7 n_lrg_sum_s2_F6 n_lrg_sum_s2_F5 ///
	n_lrg_sum_s2_F4 n_lrg_sum_s2_F3 n_lrg_sum_s2_F2 n_lrg_sum_s2_F1 n_lrg_sum_s2 ///
	n_lrg_sum_s2_L1 n_lrg_sum_s2_L2 n_lrg_sum_s2_L3 n_lrg_sum_s2_L4 ///
	n_lrg_sum_s2_L5 n_lrg_sum_s2_L6 n_lrg_sum_s2_L7 n_lrg_sum_s2_L8 ///
	n_sml_sum_s1_F8 n_sml_sum_s1_F7 n_sml_sum_s1_F6 n_sml_sum_s1_F5 ///
	n_sml_sum_s1_F4 n_sml_sum_s1_F3 n_sml_sum_s1_F2 n_sml_sum_s1_F1 n_sml_sum_s1 ///
	n_sml_sum_s1_L1 n_sml_sum_s1_L2 n_sml_sum_s1_L3 n_sml_sum_s1_L4 ///
	n_sml_sum_s1_L5 n_sml_sum_s1_L6 n_sml_sum_s1_L7 n_sml_sum_s1_L8 ///
	n_sml_sum_s2_F8 n_sml_sum_s2_F7 n_sml_sum_s2_F6 n_sml_sum_s2_F5 ///
	n_sml_sum_s2_F4 n_sml_sum_s2_F3 n_sml_sum_s2_F2 n_sml_sum_s2_F1 n_sml_sum_s2 ///
	n_sml_sum_s2_L1 n_sml_sum_s2_L2 n_sml_sum_s2_L3 n_sml_sum_s2_L4 ///
	n_sml_sum_s2_L5 n_sml_sum_s2_L6 n_sml_sum_s2_L7 n_sml_sum_s2_L8 ///
	n_prs_sum_s1_F8 n_prs_sum_s1_F7 n_prs_sum_s1_F6 n_prs_sum_s1_F5 ///
	n_prs_sum_s1_F4 n_prs_sum_s1_F3 n_prs_sum_s1_F2 n_prs_sum_s1_F1 n_prs_sum_s1 ///
	n_prs_sum_s1_L1 n_prs_sum_s1_L2 n_prs_sum_s1_L3 n_prs_sum_s1_L4 ///
	n_prs_sum_s1_L5 n_prs_sum_s1_L6 n_prs_sum_s1_L7 n_prs_sum_s1_L8 ///
	n_prs_sum_s2_F8 n_prs_sum_s2_F7 n_prs_sum_s2_F6 n_prs_sum_s2_F5 ///
	n_prs_sum_s2_F4 n_prs_sum_s2_F3 n_prs_sum_s2_F2 n_prs_sum_s2_F1 n_prs_sum_s2 ///
	n_prs_sum_s2_L1 n_prs_sum_s2_L2 n_prs_sum_s2_L3 n_prs_sum_s2_L4 ///
	n_prs_sum_s2_L5 n_prs_sum_s2_L6 n_prs_sum_s2_L7 n_prs_sum_s2_L8 ///	
	n_all_avg_s1_F8 n_all_avg_s1_F7 n_all_avg_s1_F6 n_all_avg_s1_F5 ///
	n_all_avg_s1_F4 n_all_avg_s1_F3 n_all_avg_s1_F2 n_all_avg_s1_F1 n_all_avg_s1 ///
	n_all_avg_s1_L1 n_all_avg_s1_L2 n_all_avg_s1_L3 n_all_avg_s1_L4 ///
	n_all_avg_s1_L5 n_all_avg_s1_L6 n_all_avg_s1_L7 n_all_avg_s1_L8 ///
	n_all_avg_s2_F8 n_all_avg_s2_F7 n_all_avg_s2_F6 n_all_avg_s2_F5 ///
	n_all_avg_s2_F4 n_all_avg_s2_F3 n_all_avg_s2_F2 n_all_avg_s2_F1 n_all_avg_s1 ///
	n_all_avg_s2_L1 n_all_avg_s2_L2 n_all_avg_s2_L3 n_all_avg_s2_L4 ///
	n_all_avg_s2_L5 n_all_avg_s2_L6 n_all_avg_s2_L7 n_all_avg_s2_L8 ///	
	n_lrg_avg_s1_F8 n_lrg_avg_s1_F7 n_lrg_avg_s1_F6 n_lrg_avg_s1_F5 ///
	n_lrg_avg_s1_F4 n_lrg_avg_s1_F3 n_lrg_avg_s1_F2 n_lrg_avg_s1_F1 n_lrg_avg_s1 ///
	n_lrg_avg_s1_L1 n_lrg_avg_s1_L2 n_lrg_avg_s1_L3 n_lrg_avg_s1_L4 ///
	n_lrg_avg_s1_L5 n_lrg_avg_s1_L6 n_lrg_avg_s1_L7 n_lrg_avg_s1_L8 ///
	n_lrg_avg_s2_F8 n_lrg_avg_s2_F7 n_lrg_avg_s2_F6 n_lrg_avg_s2_F5 ///
	n_lrg_avg_s2_F4 n_lrg_avg_s2_F3 n_lrg_avg_s2_F2 n_lrg_avg_s2_F1 n_lrg_avg_s2 ///
	n_lrg_avg_s2_L1 n_lrg_avg_s2_L2 n_lrg_avg_s2_L3 n_lrg_avg_s2_L4 ///
	n_lrg_avg_s2_L5 n_lrg_avg_s2_L6 n_lrg_avg_s2_L7 n_lrg_avg_s2_L8 ///
	n_sml_avg_s1_F8 n_sml_avg_s1_F7 n_sml_avg_s1_F6 n_sml_avg_s1_F5 ///
	n_sml_avg_s1_F4 n_sml_avg_s1_F3 n_sml_avg_s1_F2 n_sml_avg_s1_F1 n_sml_avg_s1 ///
	n_sml_avg_s1_L1 n_sml_avg_s1_L2 n_sml_avg_s1_L3 n_sml_avg_s1_L4 ///
	n_sml_avg_s1_L5 n_sml_avg_s1_L6 n_sml_avg_s1_L7 n_sml_avg_s1_L8 ///
	n_sml_avg_s2_F8 n_sml_avg_s2_F7 n_sml_avg_s2_F6 n_sml_avg_s2_F5 ///
	n_sml_avg_s2_F4 n_sml_avg_s2_F3 n_sml_avg_s2_F2 n_sml_avg_s2_F1 n_sml_avg_s2 ///
	n_sml_avg_s2_L1 n_sml_avg_s2_L2 n_sml_avg_s2_L3 n_sml_avg_s2_L4 ///
	n_sml_avg_s2_L5 n_sml_avg_s2_L6 n_sml_avg_s2_L7 n_sml_avg_s2_L8 ///
	n_prs_avg_s1_F8 n_prs_avg_s1_F7 n_prs_avg_s1_F6 n_prs_avg_s1_F5 ///
	n_prs_avg_s1_F4 n_prs_avg_s1_F3 n_prs_avg_s1_F2 n_prs_avg_s1_F1 n_prs_avg_s1 ///
	n_prs_avg_s1_L1 n_prs_avg_s1_L2 n_prs_avg_s1_L3 n_prs_avg_s1_L4 ///
	n_prs_avg_s1_L5 n_prs_avg_s1_L6 n_prs_avg_s1_L7 n_prs_avg_s1_L8 ///
	n_prs_avg_s2_F8 n_prs_avg_s2_F7 n_prs_avg_s2_F6 n_prs_avg_s2_F5 ///
	n_prs_avg_s2_F4 n_prs_avg_s2_F3 n_prs_avg_s2_F2 n_prs_avg_s2_F1 n_prs_avg_s2 ///
	n_prs_avg_s2_L1 n_prs_avg_s2_L2 n_prs_avg_s2_L3 n_prs_avg_s2_L4 ///
	n_prs_avg_s2_L5 n_prs_avg_s2_L6 n_prs_avg_s2_L7 n_prs_avg_s2_L8 ///
	i_all_s1_F8 i_all_s1_F7 i_all_s1_F6 i_all_s1_F5 i_all_s1_F4 i_all_s1_F3 ///
	i_all_s1_F2 i_all_s1_F1 i_all_s1 i_all_s1_L1 i_all_s1_L2 i_all_s1_L3 ///
	i_all_s1_L4 i_all_s1_L5 i_all_s1_L6 i_all_s1_L7 i_all_s1_L8 ///
	i_all_s2_F8 i_all_s2_F7 i_all_s2_F6 i_all_s2_F5 i_all_s2_F4 i_all_s2_F3 ///
	i_all_s2_F2 i_all_s2_F1 i_all_s2 i_all_s2_L1 i_all_s2_L2 i_all_s2_L3 ///
	i_all_s2_L4 i_all_s2_L5 i_all_s2_L6 i_all_s2_L7 i_all_s2_L8 ///
	i_lrg_s1_F8 i_lrg_s1_F7 i_lrg_s1_F6 i_lrg_s1_F5 i_lrg_s1_F4 i_lrg_s1_F3 ///
	i_lrg_s1_F2 i_lrg_s1_F1 i_lrg_s1 i_lrg_s1_L1 i_lrg_s1_L2 i_lrg_s1_L3 ///
	i_lrg_s1_L4 i_lrg_s1_L5 i_lrg_s1_L6 i_lrg_s1_L7 i_lrg_s1_L8 ///
	i_lrg_s2_F8 i_lrg_s2_F7 i_lrg_s2_F6 i_lrg_s2_F5 i_lrg_s2_F4 i_lrg_s2_F3 ///
	i_lrg_s2_F2 i_lrg_s2_F1 i_lrg_s2 i_lrg_s2_L1 i_lrg_s2_L2 i_lrg_s2_L3 ///
	i_lrg_s2_L4 i_lrg_s2_L5 i_lrg_s2_L6 i_lrg_s2_L7 i_lrg_s2_L8 ///
	i_sml_s1_F8 i_sml_s1_F7 i_sml_s1_F6 i_sml_s1_F5 i_sml_s1_F4 i_sml_s1_F3 ///
	i_sml_s1_F2 i_sml_s1_F1 i_sml_s1 i_sml_s1_L1 i_sml_s1_L2 i_sml_s1_L3 ///
	i_sml_s1_L4 i_sml_s1_L5 i_sml_s1_L6 i_sml_s1_L7 i_sml_s1_L8 ///
	i_sml_s2_F8 i_sml_s2_F7 i_sml_s2_F6 i_sml_s2_F5 i_sml_s2_F4 i_sml_s2_F3 ///
	i_sml_s2_F2 i_sml_s2_F1 i_sml_s2 i_sml_s2_L1 i_sml_s2_L2 i_sml_s2_L3 ///
	i_sml_s2_L4 i_sml_s2_L5 i_sml_s2_L6 i_sml_s2_L7 i_sml_s2_L8 ///	
	i_prs_s1_F8 i_prs_s1_F7 i_prs_s1_F6 i_prs_s1_F5 i_prs_s1_F4 i_prs_s1_F3 ///
	i_prs_s1_F2 i_prs_s1_F1 i_prs_s1 i_prs_s1_L1 i_prs_s1_L2 i_prs_s1_L3 ///
	i_prs_s1_L4 i_prs_s1_L5 i_prs_s1_L6 i_prs_s1_L7 i_prs_s1_L8 ///
	i_prs_s2_F8 i_prs_s2_F7 i_prs_s2_F6 i_prs_s2_F5 i_prs_s2_F4 i_prs_s2_F3 ///
	i_prs_s2_F2 i_prs_s2_F1 i_prs_s2 i_prs_s2_L1 i_prs_s2_L2 i_prs_s2_L3 ///
	i_prs_s2_L4 i_prs_s2_L5 i_prs_s2_L6 i_prs_s2_L7 i_prs_s2_L8 ///
	prc_lk_s1_F8 prc_lk_s1_F7 prc_lk_s1_F6 prc_lk_s1_F5 prc_lk_s1_F4 ///
	prc_lk_s1_F3 prc_lk_s1_F2 prc_lk_s1_F1 prc_lk_s1 prc_lk_s1_L1 prc_lk_s1_L2 ///
	prc_lk_s1_L3 prc_lk_s1_L4 prc_lk_s1_L5 prc_lk_s1_L6 prc_lk_s1_L7 ///
	prc_lk_s1_L8 ///
	rcl_lk_s1_F8 rcl_lk_s1_F7 rcl_lk_s1_F6 rcl_lk_s1_F5 rcl_lk_s1_F4 ///
	rcl_lk_s1_F3 rcl_lk_s1_F2 rcl_lk_s1_F1 rcl_lk_s1 rcl_lk_s1_L1 rcl_lk_s1_L2 ///
	rcl_lk_s1_L3 rcl_lk_s1_L4 rcl_lk_s1_L5 rcl_lk_s1_L6 rcl_lk_s1_L7 ///
	rcl_lk_s1_L8 ///
	f1_lk_s1_F8 f1_lk_s1_F7 f1_lk_s1_F6 f1_lk_s1_F5 f1_lk_s1_F4 ///
	f1_lk_s1_F3 f1_lk_s1_F2 f1_lk_s1_F1 f1_lk_s1 f1_lk_s1_L1 f1_lk_s1_L2 ///
	f1_lk_s1_L3 f1_lk_s1_L4 f1_lk_s1_L5 f1_lk_s1_L6 f1_lk_s1_L7 ///
	f1_lk_s1_L8 ///
	prc_mn_s1_F8 prc_mn_s1_F7 prc_mn_s1_F6 prc_mn_s1_F5 prc_mn_s1_F4 ///
	prc_mn_s1_F3 prc_mn_s1_F2 prc_mn_s1_F1 prc_mn_s1 prc_mn_s1_L1 prc_mn_s1_L2 ///
	prc_mn_s1_L3 prc_mn_s1_L4 prc_mn_s1_L5 prc_mn_s1_L6 prc_mn_s1_L7 ///
	prc_mn_s1_L8 ///
	rcl_mn_s1_F8 rcl_mn_s1_F7 rcl_mn_s1_F6 rcl_mn_s1_F5 rcl_mn_s1_F4 ///
	rcl_mn_s1_F3 rcl_mn_s1_F2 rcl_mn_s1_F1 rcl_mn_s1 rcl_mn_s1_L1 rcl_mn_s1_L2 ///
	rcl_mn_s1_L3 rcl_mn_s1_L4 rcl_mn_s1_L5 rcl_mn_s1_L6 rcl_mn_s1_L7 ///
	rcl_mn_s1_L8 ///
	f1_mn_s1_F8 f1_mn_s1_F7 f1_mn_s1_F6 f1_mn_s1_F5 f1_mn_s1_F4 ///
	f1_mn_s1_F3 f1_mn_s1_F2 f1_mn_s1_F1 f1_mn_s1 f1_mn_s1_L1 f1_mn_s1_L2 ///
	f1_mn_s1_L3 f1_mn_s1_L4 f1_mn_s1_L5 f1_mn_s1_L6 f1_mn_s1_L7 ///
	f1_mn_s1_L8 ///
/*
	prc_lk_s2_F8 prc_lk_s2_F7 prc_lk_s2_F6 prc_lk_s2_F5 prc_lk_s2_F4 ///
	prc_lk_s2_F3 prc_lk_s2_F2 prc_lk_s2_F1 prc_lk_s2 prc_lk_s2_L1 prc_lk_s2_L2 ///
	prc_lk_s2_L3 prc_lk_s2_L4 prc_lk_s2_L5 prc_lk_s2_L6 prc_lk_s2_L7 ///
	prc_lk_s2_L8 ///
	rcl_lk_s2_F8 rcl_lk_s2_F7 rcl_lk_s2_F6 rcl_lk_s2_F5 rcl_lk_s2_F4 ///
	rcl_lk_s2_F3 rcl_lk_s2_F2 rcl_lk_s2_F1 rcl_lk_s2 rcl_lk_s2_L1 rcl_lk_s2_L2 ///
	rcl_lk_s2_L3 rcl_lk_s2_L4 rcl_lk_s2_L5 rcl_lk_s2_L6 rcl_lk_s2_L7 ///
	rcl_lk_s2_L8 ///
	f1_lk_s2_F8 f1_lk_s2_F7 f1_lk_s2_F6 f1_lk_s2_F5 f1_lk_s2_F4 ///
	f1_lk_s2_F3 f1_lk_s2_F2 f1_lk_s2_F1 f1_lk_s2 f1_lk_s2_L1 f1_lk_s2_L2 ///
	f1_lk_s2_L3 f1_lk_s2_L4 f1_lk_s2_L5 f1_lk_s2_L6 f1_lk_s2_L7 ///
	f1_lk_s2_L8
	prc_mn_s2_F8 prc_mn_s2_F7 prc_mn_s2_F6 prc_mn_s2_F5 prc_mn_s2_F4 ///
	prc_mn_s2_F3 prc_mn_s2_F2 prc_mn_s2_F1 prc_mn_s2 prc_mn_s2_L1 prc_mn_s2_L2 ///
	prc_mn_s2_L3 prc_mn_s2_L4 prc_mn_s2_L5 prc_mn_s2_L6 prc_mn_s2_L7 ///
	prc_mn_s2_L8 ///
	rcl_mn_s2_F8 rcl_mn_s2_F7 rcl_mn_s2_F6 rcl_mn_s2_F5 rcl_mn_s2_F4 ///
	rcl_mn_s2_F3 rcl_mn_s2_F2 rcl_mn_s2_F1 rcl_mn_s2 rcl_mn_s2_L1 rcl_mn_s2_L2 ///
	rcl_mn_s2_L3 rcl_mn_s2_L4 rcl_mn_s2_L5 rcl_mn_s2_L6 rcl_mn_s2_L7 ///
	rcl_mn_s2_L8 ///
	f1_mn_s2_F8 f1_mn_s2_F7 f1_mn_s2_F6 f1_mn_s2_F5 f1_mn_s2_F4 ///
	f1_mn_s2_F3 f1_mn_s2_F2 f1_mn_s2_F1 f1_mn_s2 f1_mn_s2_L1 f1_mn_s2_L2 ///
	f1_mn_s2_L3 f1_mn_s2_L4 f1_mn_s2_L5 f1_mn_s2_L6 f1_mn_s2_L7 ///
	f1_mn_s2_L8
*/
save "$path\02_Processed_data\Politicians_Directorships_1931-2017.dta", replace
