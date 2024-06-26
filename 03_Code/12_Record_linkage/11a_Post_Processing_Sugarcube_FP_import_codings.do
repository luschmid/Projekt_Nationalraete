clear
cap log close
set more 1
version 17

*global path "C:/Schmidlu/Dropbox/Projekt Nationalräte"
*global path_rl "C:/Schmidlu/Dropbox/Record Linkage"
global path "E:/12. Cloud\Dropbox/Projekt Nationalräte"
global path_rl "E:/12. Cloud\Dropbox/Record Linkage"


* Read in all datasets from manual FP checking (team of coders)

* a) read in all manually corrected files by coder

* Brunschwiler

foreach n in 1 9 12 15 33 36 49 52 62 73 {  
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force
	gen coder = "Brunschweiler"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}

* Ehrengruber

foreach n in 2 3 4 5 6 7 10 11 21 39 72 79 { 
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force
	gen coder = "Ehrengruber"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}

* Gaberthüel

foreach n in 31 35 40 42 43 44 45 46 47 69 82 { 
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force
	gen coder = "Gaberthuel"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}

* Kotatis

foreach n in 16 17 18 19 20 22 32 34 51 71 74 { 
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force
	gen coder = "Kotatis"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}

* Haid

foreach n in 8 13 48 64 { 
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force
	gen coder = "Haid"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}

* Muenger

foreach n in 37 68 75 80 { 
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force
	gen coder = "Muenger"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}

* Heisenberg

foreach n in 14	38 41 50 76 81 { 
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force
	gen coder = "Heisenberg"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}

* Beck

foreach n in 61 78 { 
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force
	gen coder = "Beck"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}

* Reichmuth
/*
foreach n in  { // 62
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	gen coder = "Reichmuth"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}
*/
* Stern

foreach n in 63 65 66 67 70 77 { 
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force
	gen coder = "Stern"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}
/*
* Rychener

foreach n in  { // 65
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force	
	gen coder = "Rychener"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}

* Bossard

foreach n in { // 65
	import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.xlsx", first clear
	tostring fp, replace force
	gen coder = "Bossard"
	gen batchnr = `n'
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'_in.dta", replace	
}
*/

* b) Check fp = 2 & comments

use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_1_in.dta", clear
foreach m in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ///
			31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 ///
			61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 { 
	append using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`m'_in.dta", force
} 
drop if id_0 == "" // empty lines
tab fp
replace fp = "0" if fp == ""
replace fp = "0" if fp == " "
replace fp = "0" if fp == "."
replace fp = "1" if fp == "11"
replace fp = "1" if fp == "old"
replace fp = "1" if fp == "young"
replace fp = "0" if fp == " "
destring fp, replace
tab fp

*br * if fp == 2
*br * if comment != ""
/*
preserve
keep if fp == "2"
export excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_Code2.xlsx", first(var) replace
restore

preserve
keep if comment != ""
export excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_Comments.xlsx", first(var) replace
restore
*/

* c) Coder evaluation

gen n = 1
bysort coder: egen cases = total(n)
bysort coder: egen totfp = total(fp)
gen pctfp = totfp/cases*100
tab coder, sum(pctfp)

* d) Merge and compare

forv m = 1(1)22 { 
	use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`m'_in.dta", clear
	keep link_score id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 firstname_0 name_0 firstname_1 name_1 fp coder batchnr
	order link_score id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 firstname_0 name_0 firstname_1 name_1 fp coder batchnr
	drop if id_0 == "" & id_1 == ""
	rename fp fp1
	rename coder coder1
	rename batchnr batchnr1
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`m'B1_in.dta", replace 
}

forv m = 31(1)52 {
	use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`m'_in.dta", clear
	keep id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 firstname_0 name_0 firstname_1 name_1 fp coder batchnr
	order id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 firstname_0 name_0 firstname_1 name_1 fp coder batchnr
	drop if id_0 == "" & id_1 == ""
	rename fp fp2
	rename coder coder2
	rename batchnr batchnr2
	local n = `m'-30
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'B2_in.dta", replace 
}

forv m = 61(1)82 { 
use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`m'_in.dta", clear
	keep id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 firstname_0 name_0 firstname_1 name_1 fp coder batchnr
	order id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 firstname_0 name_0 firstname_1 name_1 fp coder batchnr
	drop if id_0 == "" & id_1 == "" 
	rename fp fp3
	rename coder coder3
	rename batchnr batchnr3
	local n = `m'-60
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`n'B3_in.dta", replace 
}

forv i = 1(1)22  { 
	use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B1_in.dta", clear
	merge 1:1 id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B2_in.dta"
	rename _merge merge1
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B1_2_in.dta", replace
}

forv i = 1(1)22 {
	use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B1_2_in.dta", clear
	merge 1:1 id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B3_in.dta"
	rename _merge merge2
	save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B1_3_in.dta", replace
}

use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_1B1_3_in.dta", clear
forv i = 2(1)22  {
	append using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B1_3_in.dta", force
}

forv j = 1(1)3 {
	replace fp`j' = "0" if fp`j' == ""
	replace fp`j' = "0" if fp`j' == " "
	replace fp`j' = "0" if fp`j' == "."
	replace fp`j' = "1" if fp`j' == "11"
	replace fp`j' = "1" if fp`j' == "old"
	replace fp`j' = "1" if fp`j' == "young"
	replace fp`j' = "0" if fp`j' == " "
	tab fp`j'
	destring fp`j', replace
}

gen fptotal = fp1 + fp2 + fp3
gen fp_simp = 0
replace fp_simp = 1 if fptotal >= 1
gen fp_maj = 0
replace fp_maj = 1 if fptotal >= 2
gen fp_unani = 0
replace fp_unani = 1 if fptotal == 3

tab fp_simp
tab fp_maj
tab fp_unani

* e) Coders performance

foreach coder in Brunschweiler Ehrengruber Gaberthuel Kotatis Haid Muenger Heisenberg Beck Stern {
gen `coder' = 1 if coder1 == "`coder'" | coder2 == "`coder'" | coder3 == "`coder'"
egen tot_`coder' = total(`coder')
	forv x = 1(1)3 {
		gen fp_`coder'_`x' = fp`x' if coder`x' == "`coder'"	
	}
egen fp_`coder' = rowtotal(fp_`coder'_1 fp_`coder'_2 fp_`coder'_3), missing
drop fp_`coder'_1 fp_`coder'_2 fp_`coder'_3
gen maj_`coder' = 1 if fp_`coder' == fp_maj
replace maj_`coder' = 0 if fp_`coder' != fp_maj & fp_`coder' != .
egen majtot_`coder' = total(maj_`coder')
gen majpct_`coder' = majtot_`coder' / tot_`coder' * 100
}
sum majtot_* majpct*

* f) Keep only dataset after majority voting on fp

drop if fp_maj == 1  // drop false positives
keep link_score id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_fp_maj_IDs.dta", replace

* g) erase creates temp files

foreach m in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ///
			31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 ///
			61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 { 
	erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`m'_in.dta"
} 
forv i = 1(1)22 {
	erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B1_in.dta"
	erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B2_in.dta"
	erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B3_in.dta"
	erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B1_2_in.dta"
	erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'B1_3_in.dta"
}	
