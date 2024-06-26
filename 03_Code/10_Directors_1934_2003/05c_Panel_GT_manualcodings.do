clear
cap log close
set more 1
version 17

*global path_Panel_GT "C:/Schmidlu/Dropbox/Projekt Nationalräte"
*global path_Panel_GT_rl "C:/Schmidlu/Dropbox/Record Linkage"
global path_Panel_GT "E:/12. Cloud\Dropbox/Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\12_Panel_ID_GT"


**** Establish Ground Truth for Panel ID creation in Sugarcube

/** 1) Read-in manual GT codings of first round

* a) read in all manually linked files by coder

* Ehrengruber

foreach n in 1 6 12 { 
	import excel using "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.xlsx", first clear
	gen coder = "Ehrengruber"
	gen batchnr = `n'
	tostring target, replace
	replace pid_new = "" if pid_new == " " | pid_new == "."
	keep if pid_new != ""
	save "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.dta", replace	
}

* Gaberthüel

foreach n in 2 {
	import excel using "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.xlsx", first clear
	gen coder = "Gaberthuel"
	gen batchnr = `n'
	tostring target, replace
	replace pid_new = "" if pid_new == " " | pid_new == "."
	keep if pid_new != ""
	save "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.dta", replace	
}

* Kotatis

foreach n in 3 7 15 18 19 { 
	import excel using "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.xlsx", first clear
	gen coder = "Kotatis"
	gen batchnr = `n'
	tostring target, replace
	replace pid_new = "" if pid_new == " " | pid_new == "."
	keep if pid_new != ""
	save "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.dta", replace	
}

* Haid

foreach n in 4 17 { 
	import excel using "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.xlsx", first clear
	gen coder = "Haid"
	gen batchnr = `n'
	tostring target, replace
	replace pid_new = "" if pid_new == " " | pid_new == "."
	keep if pid_new != ""
	save "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.dta", replace	
}

* Muenger

foreach n in 5 8 14 { 
	import excel using "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.xlsx", first clear
	gen coder = "Muenger"
	gen batchnr = `n'
	tostring target, replace
	replace pid_new = "" if pid_new == " " | pid_new == "."
	keep if pid_new != ""
	save "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.dta", replace	
}

* Heisenberg

foreach n in 9 16 { 
	import excel using "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.xlsx", first clear
	gen coder = "Heisenberg"
	gen batchnr = `n'
	tostring target, replace
	replace pid_new = "" if pid_new == " " | pid_new == "."
	keep if pid_new != ""
	save "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.dta", replace	
}

* Beck

foreach n in 11 13 { 
	import excel using "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.xlsx", first clear
	gen coder = "Beck"
	gen batchnr = `n'
	tostring target, replace
	replace pid_new = "" if pid_new == " " | pid_new == "."
	keep if pid_new != ""
	save "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'_in.dta", replace	
}

* b) Merge and append

forv m = 1(1)9 { 
	use "$path_Panel_GT\Panel_Round1\Panel_Round1_`m'_in.dta", clear
	keep pid year id_dupl gde_options block_id pid_new coder batchnr target lastname firstname gdename firmnamelocation
	order pid year id_dupl gde_options block_id pid_new coder batchnr target lastname firstname gdename firmnamelocation
	rename pid_new pid_new1
	rename coder coder1
	rename batchnr batchnr1
	rename target target1 
	rename lastname lastname1
	rename firstname firstname1
	rename gdename gdename1
	rename firmnamelocation firmnamelocation1
	save "$path_Panel_GT\Panel_Round1\Panel_Round1_`m'B1_in.dta", replace 
}

forv m = 11(1)19 {
	use "$path_Panel_GT\Panel_Round1\Panel_Round1_`m'_in.dta", clear
	keep pid year id_dupl gde_options block_id pid_new coder batchnr target lastname firstname gdename firmnamelocation
	order pid year id_dupl gde_options block_id pid_new coder batchnr target lastname firstname gdename firmnamelocation
	rename pid_new pid_new2
	rename coder coder2
	rename batchnr batchnr2
	rename target target2
	rename lastname lastname2
	rename firstname firstname2
	rename gdename gdename2
	rename firmnamelocation firmnamelocation2
	local n = `m'-10
	save "$path_Panel_GT\Panel_Round1\Panel_Round1_`n'B2_in.dta", replace 
}

forv i = 1(1)9 { 
	use "$path_Panel_GT\Panel_Round1\Panel_Round1_`i'B1_in.dta", clear
	merge 1:1 pid year id_dupl gde_options block_id using "$path_Panel_GT\Panel_Round1\Panel_Round1_`i'B2_in.dta"
	save "$path_Panel_GT\Panel_Round1\Panel_Round1_`i'B1_2_in.dta", replace
}

use "$path_Panel_GT\Panel_Round1\Panel_Round1_1B1_2_in.dta", clear
forv i = 2(1)9 {
	append using "$path_Panel_GT\Panel_Round1\Panel_Round1_`i'B1_2_in.dta", force
}

sort block_id year

* c) compare and prepare for 2nd Round checks 

destring target1, replace
destring target2, replace

gen pid_check = pid_new1 == pid_new2 if _merge == 3
tab pid_check  

foreach var in pid_new target lastname firstname gdename firmnamelocation {
	gen `var' = `var'1 if _merge == 3 & pid_check == 1
	replace `var' = `var'1 if _merge == 1
	replace `var' = `var'2 if _merge == 2
	}
replace pid_new = "" if _merge != 3 | pid_check == 0
gen tocheck = 0 if _merge==3
replace tocheck = 1 if _merge != 3 | pid_check == 0 

sort block_id year
order pid id_dupl gde_options block_id tocheck target pid_new lastname firstname year gdename firmnamelocation

* d) Export for 2nd Round Check

keep pid id_dupl gde_options block_id tocheck target pid_new lastname firstname year gdename firmnamelocation
export excel using "$path_Panel_GT\Panel_Round2\Panel_Round2.xlsx", first(var) replace 

* e) erase tmp files

forv m = 1(1)9 { 
	erase "$path_Panel_GT\Panel_Round1\Panel_Round1_`m'_in.dta"
	*erase "$path_Panel_GT\Panel_Round1\Panel_Round1_`m'B1_in.dta"
	*erase "$path_Panel_GT\Panel_Round1\Panel_Round1_`m'B2_in.dta"
	*erase "$path_Panel_GT\Panel_Round1\Panel_Round1_`m'B1_2_in.dta"
}
forv m = 11(1)19 { 
	erase "$path_Panel_GT\Panel_Round1\Panel_Round1_`m'_in.dta"
}

*/

** 2) Read-in manual GT codings of second round

* a) read in manually linked files by two coder

foreach n in 1 2 { 
	import excel using "$path_Panel_GT\Panel_Round2\Panel_Round2_`n'_in.xlsx", first clear
	gen coder = ""
	replace coder = "Beck" if `n' == 1
	replace coder = "Haid" if `n' == 2
	gen batchnr = `n'
	tostring target, replace
	replace pid_new = "" if pid_new == " " | pid_new == "."
	rename pid_new pid_`n'
	save "$path_Panel_GT\Panel_Round2\Panel_Round2_`n'_in.dta", replace	
}

* b) merge coders results

use "$path_Panel_GT\Panel_Round2\Panel_Round2_1_in.dta", clear
merge 1:1 pid year id_dupl gde_options block_id using "$path_Panel_GT\Panel_Round2\Panel_Round2_2_in.dta"
drop _merge

* c) prepare output for final round of coding
gen pid_diff = 1 if pid_1 != pid_2
sort block_id year
order pid id_dupl gde_options block_id target pid_diff pid_1 pid_2 year lastname firstname year gdename firmnamelocation
export excel using "$path_Panel_GT\Panel_Round3\Round_3.xlsx", first(var) replace 




















	