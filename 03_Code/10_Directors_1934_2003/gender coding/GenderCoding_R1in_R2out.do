version 17
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"


*****************
* Gender Coding *
*****************

* Read-in: First round of gender-coding

forv i = 1/2 {  // Caverzasio == 1, Oester == 2
import excel using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R1_in`i'.xlsx", firstrow clear
gen code_id = _n
tab check_R1
tab male_R1 
tab undecided_R1 
tab dontknow_R1

foreach j in cleanname male_R1 undecided_R1 dontknow_R1 check_R1 Comment {
	rename `j' `j'_`i'
}
save "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R1_in`i'.dta", replace
}

use "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R1_in1.dta", clear
merge 1:1 code_id using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R1_in2.dta"
order code_id cleanname* male_R1* undecided_R1* dontknow_R1* check_R1* Comment*
drop cleanname_2 check_R1* _merge
rename cleanname_1 cleanname
drop H I J

gen Diff_R1 = .
replace Diff_R1 = 1 if male_R1_1 != male_R1_2
replace Diff_R1 = 1 if undecided_R1_1 != undecided_R1_2
replace Diff_R1 = 1 if dontknow_R1_1 != dontknow_R1_2

order nobs code_id cleanname Diff_R1 male_R1_1 male_R1_2 undecided_R1_1 undecided_R1_2 dontknow_R1_1 dontknow_R1_2 Comment_1 Comment_2
save "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R1.dta", replace

/* Produce files for second round (R2) of coding
gen male_R2 = .
gen undecided_R2 = .
gen dontknow_R2 = .
gen check_R2 = ""
gen Comment_R2 = ""

order code_id cleanname Diff_R1 male_R2 undecided_R2 dontknow_R2 check_R2 male_R1_1 male_R1_2 undecided_R1_1 undecided_R1_2 dontknow_R1_1 dontknow_R1_2 Comment_1 Comment_2 Comment_R2 

keep if Diff_R1 == 1
drop Diff_R1

export excel using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R2_out1.xlsx", firstrow(variables) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R2_out2.xlsx", firstrow(variables) replace


