clear
cap log close
set more 1
version 17

*global path "C:/Schmidlu/Dropbox/Projekt Nationalräte"
*global path_rl "C:/Schmidlu/Dropbox/Record Linkage"
global path "E:/12. Cloud\Dropbox/Projekt Nationalräte"
global path_rl "E:/12. Cloud\Dropbox/Record Linkage"
global sugarcube_model_name "15122022-100518bc-12f0-4d5d-8b42-1f4484952a2a"


* (i) Read in Sugarcube results by Sandro and transform them into wide format

insheet using ///
	"$path_rl\02_Data\08_Output_RL_Sugarcube\RL-Results-15122022-100518bc-12f0-4d5d-8b42-1f4484952a2a.csv", ///
	delim(",") clear
	
* a) Keep politician data

preserve
keep if sourcefile==0
keep cluster_id link_score id e_cntr_w n_cntr_w
rename id id_0 
rename e_cntr_w e_cntr_w_0 
rename n_cntr_w n_cntr_w_0 
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\rl_sugarcube_nr.dta", replace
restore

* b) Keep Sugarcube data

keep if sourcefile==1
keep cluster_id	e_id n_id id year firstname name

rename id id_1 
rename year year_1 
rename e_id e_cntr_w_1 
rename n_id n_cntr_w_1 
rename firstname firstname_1 
rename name name_1 

* c) Merge

merge m:1 cluster_id using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\rl_sugarcube_nr.dta"
drop _merge
order cluster_id link_score id_0 e_cntr_w_0 n_cntr_w_0  id_1 year e_cntr_w_1 n_cntr_w_1 
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\rl_sugarcube_cortable.dta", replace

	
*reshape wide name firstname sex e_cntr_w n_cntr_w birthyear e_cntr_b n_cntr_b id, i(cluster_id) j(sourcefile)
	
	
* (ii) Prepare NR files
	
use "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta",clear

sort ID year
bysort ID: gen timevar=_n
rename sex male

* a) Information variables (only used for manual checks)

foreach var of varlist firstname name job list male gdename{
preserve 
keep ID `var' year timevar
duplicates drop ID `var', force
keep ID `var' timevar
reshape wide `var', i(ID) j(timevar)
gen `var'_0=`var'1
forvalues i = 2(1)11 {
cap replace `var'_0=`var'_0+"; " + `var'`i' if `var'`i'!=""
}

keep ID `var'_0 
rename ID id_0

save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_`var'.dta", replace
restore 		
}

* b) Birthyear

preserve
bysort ID: egen birthyear_mean=mean(birthyear)
keep ID birthyear_mean
duplicates drop ID birthyear_mean, force
gen birthyear=floor(birthyear_mean)
keep ID birthyear
rename ID id_0
rename birthyear birthyear_0
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_birthyear.dta", replace
restore

/*
* c) Geocodes

preserve
keep ID E_CNTR_w N_CNTR_w
duplicates drop ID E_CNTR_w N_CNTR_w, force
rename E_CNTR_w e_cntr_w_1 
rename N_CNTR_w n_cntr_w_1 
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_geocodes.dta", replace
restore
*/

* (iii) Prepare Sugarcube files

use "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_RLPostProcessing_Persons-Firmnames.dta",clear   // original/early version from May 2023

* a) Keep only selected variables

keep id_sug year_sug e_id_sug n_id_sug firstname lastname male firmnames gdename
rename id_sug id_1 
rename year_sug year_1
rename n_id_sug n_cntr_w_1
rename e_id_sug e_cntr_w_1
rename firstname firstname_1 
rename lastname name_1
rename male male_1
rename gdename gdename_1 
rename firmnames firmnames_1
tostring id_1, replace
tostring year_1, replace

* b) Drop duplicates in terms of id_1 year_1 n_cntr_w_1 e_cntr_w_1 firstname_1 name_1

duplicates drop id_1 year_1 n_cntr_w_1 e_cntr_w_1 firstname_1 name_1, force

save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\sug_to_merge.dta", replace



* (iv) Merge all datasets for FP check

use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\rl_sugarcube_cortable.dta", clear

*duplicates tag id_0 e_cntr_w_0 n_cntr_w_0 id_1 year_1 e_cntr_w_1 n_cntr_w_1, gen(dups)
*sort id_0 e_cntr_w_0 n_cntr_w_0 id_1 year_1 n_cntr_w_1 e_cntr_w_1
duplicates drop id_0 e_cntr_w_0 n_cntr_w_0 id_1 year_1 e_cntr_w_1 n_cntr_w_1, force

* a) Politician data

local nr_vars "firstname name job list male gdename birthyear"
foreach var of local nr_vars {
merge m:1 id_0 using  "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_`var'.dta", gen(merge_`var')
}

drop if merge_firstname==2
drop merge*

* b) Sugarcube data

merge m:1 id_1 year_1 n_cntr_w_1 e_cntr_w_1 firstname_1 name_1 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\sug_to_merge.dta", gen(merge_sug)
keep if merge_sug==3

* c) NR data (candidacy years)
preserve 
use "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta", clear
keep ID year
bysort ID: egen maxcandyear = max(year)
tostring maxcandyear, replace
bysort ID: egen mincandyear = min(year)
tostring mincandyear, replace
gen candyears = mincandyear + "-" + maxcandyear
drop maxcandyear mincandyear year
duplicates drop ID candyears, force
rename ID id_0
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\NR_candyears.dta", replace
restore

* d) Additional variables (distance, implicit age, time in office)

gen fp=""
gen comment=""
gen w_dist=sqrt((e_cntr_w_1- e_cntr_w_0)^2+(n_cntr_w_1-n_cntr_w_0)^2)/1000
destring year_1, replace
gen age_implicit=year_1-birthyear
sort id_0
merge m:1 id_0 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\NR_candyears.dta"
keep if _merge == 3
drop _merge


* e) Indicate FP based on age

replace fp = "old" if age_implicit > 99  // cut-off according to Bisnode decision (see Email Lukas: 19.11.2022)
replace fp = "young" if age_implicit < 18  // cut-off according to Bisnode decision (see Email Lukas: 19.11.2022)

* f) Sorting

egen id_0_num=group(id_0)

gen firstname_1_first3=substr(firstname_1,1,3)
gen name_1_first3=substr(name_1,1,3)

sort id_0_num firstname_1_first3 name_1_first3 n_cntr_w_1 e_cntr_w_1 year_1

* Ground Truth Generation files: ID_Pol,Vorname_VR,Nachname_VR,Wohngemeinde_Pol,Wohngemeinde_VR,Jahr_VR

drop merge_sug

* g) Ordering of variables

order link_score cluster_id id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 fp comment id_0_num firstname_0 name_0 firstname_1 name_1 year_1 w_dist candyears age_implicit  firmnames_1 male_0 male_1 gdename_0 gdename_1 birthyear_0

export excel using 	"$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1_All.xlsx", first(var) replace

* h) Cutting dataset into 22 blocks (527 ids per block)

forvalues i=1(1)22{
preserve
keep if inrange(id_0_num,((`i'-1)*527)+1,`i'*527)
local j = `i' + 30  // all files recieve a twin numbered with i+30
local k = `i' + 60  // all files recieve a twin numbered with i+60
export excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1_`i'.xlsx", first(var) replace
export excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1_`j'.xlsx", first(var) replace
export excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1_`k'.xlsx", first(var) replace
restore	
}

** erase temporary files
/*foreach var in firstname name job list male gdename {
	erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_`var'.dta"
}*/
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\rl_sugarcube_nr.dta"
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\rl_sugarcube_cortable.dta"
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_birthyear.dta"
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\sug_to_merge.dta"
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\NR_candyears.dta"

