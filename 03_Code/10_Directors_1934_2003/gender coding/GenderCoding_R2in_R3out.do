version 17
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"


*****************
* Gender Coding *
*****************

* Read-in: Second round of gender-coding

forv i = 1/2 {  // Caverzasio == 1, Oester == 2
import excel using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R2_in`i'.xlsx", firstrow cellrange(A1:O1630) clear
tab check_R2
tab male_R2 
tab undecided_R2 
tab dontknow_R2

foreach j in cleanname male_R2 undecided_R2 dontknow_R2 check_R2 Comment_R2 {
	rename `j' `j'_`i'
}
save "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R2_in`i'.dta", replace
}

use "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R2_in1.dta", clear
merge 1:1 code_id using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R2_in2.dta"
order code_id cleanname* male_R2* undecided_R2* dontknow_R2* check_R2* Comment_R2*
drop cleanname_2 check_R2* _merge
rename cleanname_1 cleanname

sort cleanname

// merge with BfS: Vornamen der Bevölkerung nach Geschlecht, Schweiz, 2015-2019  (https://www.bfs.admin.ch/asset/de/13927395)
preserve
forv i = 2015(1)2019 {
import excel using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\su-q-01.04.00.12.xlsx", sheet("`i'") clear
rename A cleanname
rename B female
rename C male 
replace cleanname = usubinstr(cleanname, "-", "",.)
gen rows = _n
drop if rows < 7
drop if cleanname == ""
replace female = "0" if female == "*"  // flag! most probably from BfS perspective missing != 0  -- > but this is irrelevant for this application
replace male = "0" if male == "*"
destring female, gen(femaleBfS)
destring male, gen(maleBfS)
gen year = `i'
drop male female
save "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\BfS_`i'.dta", replace
}

forv i = 2018(-1)2015 {
	append using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\BfS_`i'.dta"
}
collapse (mean) maleBfS femaleBfS, by(cleanname)
save "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\BfS_Vornamen2015-19.dta", replace
forv i = 2015(1)2019 {
	erase "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\BfS_`i'.dta"
}	
restore

merge 1:1 cleanname using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\BfS_Vornamen2015-19.dta"
drop if _merge == 2
drop _merge

egen genBfS_tot = rowtotal(maleBfS femaleBfS)
replace genBfS_tot = . if genBfS_tot == 0
gen genBfS_mbal = maleBfS / genBfS_tot
gen genderBfS = 1 if genBfS_mbal >= 0.9 & genBfS_mbal != .
replace genderBfS = 2 if genBfS_mbal <= 0.1 & genBfS_mbal != .
gen genderBfS_und = 1 if genBfS_mbal != . & (genBfS_mbal > 0.1 & genBfS_mbal < 0.9)  // gender proportions too asymmetric (90% cutoff)

gen Diff_R2 = 0
replace Diff_R2 = 1 if male_R2_1 != male_R2_2
replace Diff_R2 = 1 if undecided_R2_1 != undecided_R2_2
replace Diff_R2 = 1 if dontknow_R2_1 != dontknow_R2_2

gen male1 = 1 if male_R1_1 == 1
gen male2 = 1 if male_R1_2 == 1
gen male3 = 1 if male_R2_1 == 1
gen male4 = 1 if male_R2_2 == 1

gen female1 = 1 if male_R1_1 == 2
gen female2 = 1 if male_R1_2 == 2
gen female3 = 1 if male_R2_1 == 2
gen female4 = 1 if male_R2_2 == 2

gen male_R3 = male_R2_1 if male_R2_1 == male_R2_2
gen undecided_R3 = undecided_R2_1 if undecided_R2_1 == undecided_R2_2
gen dontknow_R3 = dontknow_R2_1 if dontknow_R2_1 == dontknow_R2_2

gen Flag_BfS = 1 if male_R3 != . & male_R3 != genderBfS	& genderBfS != .	// flag difference in coding and BfS
replace Flag_BfS = 1 if undecided_R3 == 1 & genderBfS != .
replace Flag_BfS = 1 if genderBfS_H == 1 & male_R3 != .

replace male_R3 = genderBfS if male_R3 == . & male_R2_1 != male_R2_2 & 	undecided_R3 != 1	// fill in BfS if no coding consensus
replace undecided_R3 = genderBfS_und if genderBfS_und == 1 & undecided_R3 == .   // take BfS indifferent indicator if proportions too are asymmetric

gen check_R3 = ""
gen Comment_R3 = ""

egen R_male = rowtotal(male1 male2 male3 male4)
replace R_male = R_male / 4
replace R_male = . if R_male == 0
egen R_female = rowtotal(female1 female2 female3 female4)
replace R_female = R_female / 4
replace R_female = . if R_female == 0
egen R_undecided = rowtotal(undecided_R2_1 undecided_R2_2 undecided_R1_1 undecided_R1_2)
replace R_undecided = R_undecided / 4
replace R_undecided = . if R_undecided == 0
egen R_dontknow = rowtotal(dontknow_R2_1 dontknow_R2_2 dontknow_R1_1 dontknow_R1_2)
replace R_dontknow = R_dontknow / 4
replace R_dontknow = . if R_dontknow == 0

order code_id cleanname Diff_R2 check_R3 Flag_BfS male_R3 undecided_R3 dontknow_R3 genderBfS genderBfS_und maleBfS femaleBfS R_male R_female R_undecided R_dontknow male_R2* undecided_R2* dontknow_R2* Comment_R2*
drop male1 male2 male3 male4 female1 female2 female3 female4


/* Produce files for third round (R3) of coding
export excel using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R3_out.xlsx", firstrow(variables) replace



