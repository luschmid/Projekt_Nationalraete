clear
cap log close
set more 1
version 17

global data "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\12_Record_linkage\01_Bisnode\04_Post_RL_QualChecks\"

/** Post-Record Linkage cleaning - Optimal (F1 score)
* Eliminate false positives among RL matches with optimal F1 score (conservative approach: eliminate clear cases only)

* Read-in Ladina Oester/Nathan Marilley

import excel using "$data\bisnode_fp_firstcheck_Generation7_optimal_oester_in.xlsx", firstrow clear
gen obs = _n
save "$data\bisnode_fp_firstcheck_Generation7_optimal_oester_in.dta", replace

import excel using "$data\bisnode_fp_firstcheck_Generation7_optimal_marilley_in.xlsx", firstrow clear
gen obs = _n
foreach var in id_0 id_1 Codierung firstname_pol lastname_pol firstname_bis lastname_bis sex_pol sex_bis distance age_diff age_firstmandate age_lastmandate {
	rename `var' m_`var'
}
save "$data\bisnode_fp_firstcheck_Generation7_optimal_marilley_in.dta", replace

merge 1:1 obs using "$data\bisnode_fp_firstcheck_Generation7_optimal_oester_in.dta"
drop _merge

rename Codierung CodierungOester
rename m_Codierung CodierungMarilley

list * if  m_id_0 != id_0 & m_id_1 != id_1
list * if  m_firstname_pol != firstname_pol & m_lastname_bis != lastname_bis

drop m_*
gen flag1 = 1 if CodierungOester != CodierungMarilley
gen flag2 = 1 if CodierungOester == 2 | CodierungMarilley == 2

preserve
keep if flag1 == 1
gen CodierungCheck1 = .
order obs id_0 id_1 CodierungCheck1 Codierung*
export excel using "$data\bisnode_fp_firstcheck_Generation7_optimal_OeMaDiff1.xlsx", firstrow(variables) replace
// manual checking/recoding
import excel using "$data\bisnode_fp_firstcheck_Generation7_optimal_OeMaDiff1Fin.xlsx", firstrow clear
keep obs CodierungCheck1
save "$data\bisnode_fp_firstcheck_Generation7_optimal_Coding1Fin.dta", replace
restore

preserve
keep if flag2 == 1 & flag1 != 1
gen CodierungCheck2 = .
order obs id_0 id_1 CodierungCheck2 Codierung*
export excel using "$data\bisnode_fp_firstcheck_Generation7_optimal_OeMaDiff2.xlsx", firstrow(variables) replace
// manual checking/recoding
import excel using "$data\bisnode_fp_firstcheck_Generation7_optimal_OeMaDiff2Fin.xlsx", firstrow clear
keep obs CodierungCheck2
save "$data\bisnode_fp_firstcheck_Generation7_optimal_Coding2Fin.dta", replace
restore

merge 1:1 obs using "$data\bisnode_fp_firstcheck_Generation7_optimal_Coding1Fin.dta"
drop _merge
merge 1:1 obs using "$data\bisnode_fp_firstcheck_Generation7_optimal_Coding2Fin.dta"
drop _merge

gen CodierungFin = CodierungCheck1  // checked all cases in which CodierungOester != CodierungMarilley
replace CodierungFin = 1 if CodierungOester == 1 & CodierungOester == CodierungMarilley & CodierungCheck1 == .
replace CodierungFin = 1 if CodierungCheck2 == 1  // checked all cases where Marilley and Oester where unsure
label var CodierungFin "1: False Positive"
drop CodierungCheck1 CodierungCheck2 CodierungMarilley CodierungOester
order obs id_0 id_1 CodierungFin
export excel using "$data\bisnode_fp_firstcheck_Generation7_optimal_corrected.xlsx", firstrow(variables) replace
export delimited using "$data\bisnode_fp_firstcheck_Generation7_optimal_corrected.csv",  delimiter(",")  replace

erase "$data\bisnode_fp_firstcheck_Generation7_optimal_Coding1Fin.dta"
erase "$data\bisnode_fp_firstcheck_Generation7_optimal_Coding2Fin.dta"

*/


** Post-Record Linkage cleaning - "0"
* Eliminate false positives among RL matches with optimal F1 score (conservative approach: eliminate clear cases only)

* Read-in Ladina Oester/Nathan Marilley

import excel using "$data\bisnode_fp_firstcheck_Generation7_0_oester_in.xlsx", firstrow clear
*gen obs = _n
save "$data\bisnode_fp_firstcheck_Generation7_0_oester_in.dta", replace

import excel using "$data\bisnode_fp_firstcheck_Generation7_0_marilley_in.xlsx", firstrow clear
*gen obs = _n
foreach var in id_0 id_1 Codierung firstname_pol lastname_pol firstname_bis lastname_bis sex_pol sex_bis distance age_diff age_firstmandate age_lastmandate {
	rename `var' m_`var'
}
save "$data\bisnode_fp_firstcheck_Generation7_0_marilley_in.dta", replace

merge 1:1 obs using "$data\bisnode_fp_firstcheck_Generation7_0_oester_in.dta"
drop _merge

rename Codierung CodierungOester
rename m_Codierung CodierungMarilley

list * if  m_id_0 != id_0 & m_id_1 != id_1
list * if  m_firstname_pol != firstname_pol & m_lastname_bis != lastname_bis

drop m_*
gen flag1 = 1 if CodierungOester != CodierungMarilley
gen flag2 = 1 if CodierungOester == 2 | CodierungMarilley == 2

preserve
keep if flag1 == 1
gen CodierungCheck1 = .
order obs id_0 id_1 CodierungCheck1 Codierung*
export excel using "$data\bisnode_fp_firstcheck_Generation7_0_OeMaDiff1.xlsx", firstrow(variables) replace
// manual checking/recoding
import excel using "$data\bisnode_fp_firstcheck_Generation7_0_OeMaDiff1Fin.xlsx", firstrow clear
keep obs CodierungCheck1
save "$data\bisnode_fp_firstcheck_Generation7_0_Coding1Fin.dta", replace
restore

preserve
keep if flag2 == 1 & flag1 != 1
gen CodierungCheck2 = .
order obs id_0 id_1 CodierungCheck2 Codierung*
export excel using "$data\bisnode_fp_firstcheck_Generation7_0_OeMaDiff2.xlsx", firstrow(variables) replace
// manual checking/recoding
import excel using "$data\bisnode_fp_firstcheck_Generation7_0_OeMaDiff2Fin.xlsx", firstrow clear
keep obs CodierungCheck2
drop if obs==.
save "$data\bisnode_fp_firstcheck_Generation7_0_Coding2Fin.dta", replace
restore

merge 1:1 obs using "$data\bisnode_fp_firstcheck_Generation7_0_Coding1Fin.dta"
drop _merge
merge 1:1 obs using "$data\bisnode_fp_firstcheck_Generation7_0_Coding2Fin.dta"
drop _merge

gen CodierungFin = CodierungCheck1  // checked all cases in which CodierungOester != CodierungMarilley
replace CodierungFin = 1 if CodierungOester == 1 & CodierungOester == CodierungMarilley & CodierungCheck1 == .
replace CodierungFin = 1 if CodierungCheck2 == 1  // checked all cases where Marilley and Oester where unsure
label var CodierungFin "1: False Positive"
drop CodierungCheck1 CodierungCheck2 CodierungMarilley CodierungOester
order obs id_0 id_1 CodierungFin
export excel using "$data\bisnode_fp_firstcheck_Generation7_0_corrected.xlsx", firstrow(variables) replace
export delimited using "$data\bisnode_fp_firstcheck_Generation7_0_corrected.csv",  delimiter(",")  replace

erase "$data\bisnode_fp_firstcheck_Generation7_0_Coding1Fin.dta"
erase "$data\bisnode_fp_firstcheck_Generation7_0_Coding2Fin.dta"

