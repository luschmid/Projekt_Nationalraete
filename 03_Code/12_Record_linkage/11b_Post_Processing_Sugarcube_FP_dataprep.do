clear
cap log close
set more 1
version 17

global path "C:/Schmidlu/Dropbox/Projekt Nationalräte"
global path_rl "C:/Schmidlu/Dropbox/Record Linkage"
*global path "E:/12. Cloud\Dropbox/Projekt Nationalräte"
*global path_rl "E:/12. Cloud\Dropbox/Record Linkage"
global sugarcube_model_name "15122022-100518bc-12f0-4d5d-8b42-1f4484952a2a"


* (i) Read in Sugarcube results by Sandro and transform them into wide format

* a) Old file with deduplication error

insheet using ///
	"$path_rl\02_Data\08_Output_RL_Sugarcube\RL-Results-15122022-100518bc-12f0-4d5d-8b42-1f4484952a2a.csv", ///
	delim(",") clear
	
replace year="." if year=="nan"
	
*duplicates drop id year e_id n_id sourcefile , force

save "$path_rl\02_Data\08_Output_RL_Sugarcube\RL-Results-15122022-100518bc-12f0-4d5d-8b42-1f4484952a2a.dta",replace

* b) New file after correction deduplication error (31 August 2023)

insheet using ///
"$path_rl\02_Data\08_Output_RL_Sugarcube\extended_result_set_01092023\extended_result_set_15122022-100518bc-12f0-4d5d-8b42-1f4484952a2a-01092023.csv", ///
	delim(",") clear

rename clusterid cluster_id
rename linkscore link_score

duplicates report cluster_id id year e_id n_id sourcefile

tostring year, replace

merge 1:1 cluster_id id year e_id n_id sourcefile using "$path_rl\02_Data\08_Output_RL_Sugarcube\RL-Results-15122022-100518bc-12f0-4d5d-8b42-1f4484952a2a.dta"

gen newobservation=0
replace newobservation=1 if _merge==1
 

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
keep cluster_id	e_id n_id id year firstname name newobservation

rename id id_1 
rename year year_1 
rename e_id e_cntr_w_1 
rename n_id n_cntr_w_1 
rename firstname firstname_1 
rename name name_1 

* c) Merge

merge m:1 cluster_id using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\rl_sugarcube_nr.dta"
drop _merge
order newobservation cluster_id link_score id_0 e_cntr_w_0 n_cntr_w_0  id_1 year_1 e_cntr_w_1 n_cntr_w_1 
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

use "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_RLPostProcessing_Persons-Firmnames.dta",clear

* Compare newest version of data-file with wide-firmnames with previous ones (from 3.5.2023)
preserve  //  old file (Sugarcube_RLPost..._230503/ contained 12 missing geo-codes
use "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_RLPostProcessing_Persons-Firmnames_230503.dta", clear
gen geo_flag = 0
replace geo_flag = 1 if id_sug == 106906 & year_sug == 1992
replace geo_flag = 1 if id_sug == 112413 & year_sug == 1991
replace geo_flag = 1 if id_sug == 158857 & year_sug == 1989
replace geo_flag = 1 if id_sug == 2436 & year_sug == 1995
replace geo_flag = 1 if id_sug == 3948 & year_sug == 1986
replace geo_flag = 1 if id_sug == 45255 & year_sug == 1969
replace geo_flag = 1 if id_sug == 56200 & year_sug == 1992
replace geo_flag = 1 if id_sug == 77810 & year_sug == 1995
replace geo_flag = 1 if id_sug == 87936 & year_sug == 1988
replace geo_flag = 1 if id_sug == 93820 & year_sug == 1989
replace geo_flag = 1 if id_sug == 9682 & year_sug == 1969
save "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_RLPostProcessing_Persons-Firmnames_230503_tmp.dta", replace
restore

merge 1:1 id_sug year_sug e_id_sug n_id_sug ID_dupl firstname lastname male using "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_RLPostProcessing_Persons-Firmnames_230503_tmp.dta"
drop if _merge==2 


* a) Keep only selected variables

keep id_sug year_sug e_id_sug n_id_sug firstname lastname male firmnames gdename geo_flag
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

unique id_0
unique id_0 if newobservation==0


*duplicates tag id_0 e_cntr_w_0 n_cntr_w_0 id_1 year_1 e_cntr_w_1 n_cntr_w_1, gen(dups)
*sort id_0 e_cntr_w_0 n_cntr_w_0 id_1 year_1 n_cntr_w_1 e_cntr_w_1
*bysort id_0 e_cntr_w_0 n_cntr_w_0 id_1 year_1 e_cntr_w_1 n_cntr_w_1: egen newobs=sd(newobservation)
duplicates drop id_0 e_cntr_w_0 n_cntr_w_0 id_1 year_1 e_cntr_w_1 n_cntr_w_1 newobservation, force



* a) Politician data

local nr_vars "firstname name job list male gdename birthyear"
foreach var of local nr_vars {
merge m:1 id_0 using  "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_`var'.dta", gen(merge_`var')
}

drop if merge_firstname==2 // drop all nr candidates with no match
drop merge*

* b) Sugarcube data

merge m:1 id_1 year_1 n_cntr_w_1 e_cntr_w_1 name_1 firstname_1 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\sug_to_merge.dta", gen(merge_sug)
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
gen age_implicit=year_1-birthyear_0
sort id_0
merge m:1 id_0 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\NR_candyears.dta"
keep if _merge == 3
drop _merge


* e) Indicate FP based on age

replace fp = "old" if age_implicit > 99  // cut-off according to Bisnode decision (see Email Lukas: 19.11.2022)
replace fp = "young" if age_implicit < 18  // cut-off according to Bisnode decision (see Email Lukas: 19.11.2022)

* f) Sorting

preserve
use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1_All.dta", clear
keep id_0 id_0_num
duplicates drop
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\id_0_num.dta", replace
restore

merge m:1 id_0 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\id_0_num.dta"


*egen id_0_num=group(id_0)

gen firstname_1_first3=substr(firstname_1,1,3)
gen name_1_first3=substr(name_1,1,3)

sort id_0_num firstname_1_first3 name_1_first3 n_cntr_w_1 e_cntr_w_1 year_1

* Ground Truth Generation files: ID_Pol,Vorname_VR,Nachname_VR,Wohngemeinde_Pol,Wohngemeinde_VR,Jahr_VR

drop merge_sug

* g) Ordering of variables

order newobservation geo_flag link_score cluster_id id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 fp comment id_0_num firstname_0 name_0 firstname_1 name_1 year_1 w_dist candyears age_implicit firmnames_1 male_0 male_1 gdename_0 gdename_1 birthyear_0


* h) Save complete dataset of Round 1b (after deduplication error by Sandro, 31 August 2023)

*export excel using 	"$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_All.xlsx", first(var) replace

save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_All.dta", replace

/*
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_1_in.xlsx", clear first

preserve
keep id_0
duplicates drop id_0, force
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1_1_temp.dta", replace
restore

tostring comment, replace

append using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_All.dta"

merge m:1 id_0 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1_1_temp.dta"
*/

* h) Cutting dataset into 22 blocks (527 ids per block)

matrix to_check=J(100,1,0)

use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_All.dta", clear

forvalues i=1(1)22{
* Part (A): all files from 1 to 22
* step 1: read in old recoded files (May to August 2023)
preserve
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'_in.xlsx", first clear
drop if cluster_id==. & id_0==""
keep id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 fp comment
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'_in.dta", replace
restore
* step 2: read in new files (September 2023) and merge old recoded files
preserve
keep if inrange(id_0_num,((`i'-1)*527)+1,`i'*527)
drop comment
rename fp fp_new

merge 1:1 id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`i'_in.dta", gen(mergevar)
tab newobservation mergevar
list id_0 id_1 year_1 firstname_1 name_1 if newobservation==0 & mergevar==1
gen obs_tocheck = mergevar==1
*sum obs_tocheck if obs_tocheck==1
*matrix to_check[`i',1]=`r(n)'
sort id_0_num firstname_1_first3 name_1_first3 n_cntr_w_1 e_cntr_w_1 year_1
order link_score cluster_id id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 fp comment id_0_num firstname_0 name_0 firstname_1 name_1 year_1 w_dist candyears age_implicit firmnames_1 male_0 male_1 gdename_0 gdename_1 birthyear_0 job_0 list_0 obs_tocheck
drop newobservation	geo_flag _merge	firstname_1_first3 name_1_first3 fp_new	mergevar
*export excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`i'_out.xlsx", first(var) replace
restore
* Part (B): all files from 31 to 52
local j = `i' + 30  
* step 1: read in old recoded files (May to August 2023)
preserve
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`j'_in.xlsx", first clear
drop if cluster_id==. & id_0==""
keep id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 fp comment
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`j'_in.dta", replace
restore
* step 2: read in new files (September 2023) and merge old recoded files
preserve
keep if inrange(id_0_num,((`i'-1)*527)+1,`i'*527)
drop comment
rename fp fp_new
merge 1:1 id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`j'_in.dta", gen(mergevar)
tab newobservation mergevar
list id_0 id_1 year_1 firstname_1 name_1 if newobservation==0 & mergevar==1
gen obs_tocheck = mergevar==1
*sum obs_tocheck if obs_tocheck==1
*matrix to_check[`j',1]=`r(n)'
sort id_0_num firstname_1_first3 name_1_first3 n_cntr_w_1 e_cntr_w_1 year_1
order link_score cluster_id id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 fp comment id_0_num firstname_0 name_0 firstname_1 name_1 year_1 w_dist candyears age_implicit firmnames_1 male_0 male_1 gdename_0 gdename_1 birthyear_0 job_0 list_0 obs_tocheck
drop newobservation	geo_flag _merge	firstname_1_first3 name_1_first3 fp_new	mergevar
*export excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`j'_out.xlsx", first(var) replace
restore

* Part (C): all files from 61 to 82
local k = `i' + 60  
* step 1: read in old recoded files (May to August 2023)
preserve
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`k'_in.xlsx", first clear
drop if cluster_id==. & id_0==""
keep id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 fp comment
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`k'_in.dta", replace
restore
* step 2: read in new files (September 2023) and merge old recoded files
preserve
keep if inrange(id_0_num,((`i'-1)*527)+1,`i'*527)
drop comment
rename fp fp_new
merge 1:1 id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 year_1 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1_`k'_in.dta", gen(mergevar)
tab newobservation mergevar
list id_0 id_1 year_1 firstname_1 name_1 if newobservation==0 & mergevar==1
gen obs_tocheck = mergevar==1
*sum obs_tocheck if obs_tocheck==1
*matrix to_check[`k',1]=`r(n)'
sort id_0_num firstname_1_first3 name_1_first3 n_cntr_w_1 e_cntr_w_1 year_1
order link_score cluster_id id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 fp comment id_0_num firstname_0 name_0 firstname_1 name_1 year_1 w_dist candyears age_implicit firmnames_1 male_0 male_1 gdename_0 gdename_1 birthyear_0 job_0 list_0 obs_tocheck
drop newobservation	geo_flag _merge	firstname_1_first3 name_1_first3 fp_new	mergevar
*export excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`k'_out.xlsx", first(var) replace
restore
}


/*
* i) Quality checks: are Simon's and Lukas's old and young codings identical?


forvalues i=20(1)20{
local j = `i' + 30  // all files receive a twin numbered with i+30
local k = `i' + 60  // all files recevie a twin numbered with i+60
display "file number: " `i'
clear
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`i'_out.xlsx", first
order fp
set seed 1234
sort link_score-obs_tocheck
gen row=_n
keep if obs_tocheck==1
keep row fp firstname_0 name_0 firstname_1 name_1 n_cntr_w_1 e_cntr_w_1 year_1 firmnames_1 id_1
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\temp_`i'.dta", replace
clear
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`j'_out.xlsx", first
order fp
set seed 1234
sort link_score-obs_tocheck
gen row=_n
keep if obs_tocheck==1
keep row fp firstname_0 name_0 firstname_1 name_1 n_cntr_w_1 e_cntr_w_1 year_1 firmnames_1 id_1
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\temp_`j'.dta", replace
clear
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`k'_out.xlsx", first
order fp
set seed 1234
sort link_score-obs_tocheck
gen row=_n
keep if obs_tocheck==1
keep row fp firstname_0 name_0 firstname_1 name_1 n_cntr_w_1 e_cntr_w_1 year_1 firmnames_1 id_1
merge 1:1 row fp using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\temp_`i'.dta", gen(m1)
merge 1:1 row fp using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\temp_`j'.dta", gen(m2)
}

* Yes, after some iterations!

* j) Quality checks: are files of Round1b the same as of Round1?

forvalues i=1(1)22{
display "file number: " `i'
quietly{
clear
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`i'.xlsx", clear first
drop if newobservation==1
drop if cluster_id==.
replace link_score=round(link_score,5)
replace w_dist=round(w_dist,5)
sort cluster_id id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 year_1 e_cntr_w_1 firstname_0 name_0 firstname_1	name_1
gen id_totest=_n 
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`i'.dta", replace
clear
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1_`i'_out.xlsx", clear first
replace link_score=round(link_score,5)
replace w_dist=round(w_dist,5)
drop if cluster_id==.
sort cluster_id link_score id_0 e_cntr_w_0 n_cntr_w_0 id_1 year_1 n_cntr_w_1 e_cntr_w_1 firstname_0	name_0	firstname_1	name_1
gen id_totest=_n 
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1_`i'.dta", replace
}
use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1_`i'.dta", clear
cf3 _all using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`i'.dta", id(id_totest) noverbose
*local j = `i' + 30  // all files recieve a twin numbered with i+30
*local k = `i' + 60  // all files recieve a twin numbered with i+60
*export excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`j'.xlsx", first(var) replace
*export excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\10_RL_Output_Round1_Out\RL_Output_Round1b_`k'.xlsx", first(var) replace
}
*/


* result: no difference between the files everything fine. This test was done before
*		  Mark found out that we have problems with 13 geocodes (Zoom 4 September 2023).


** erase temporary files
/*foreach var in firstname name job list male gdename {
	erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_`var'.dta"
}
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\rl_sugarcube_nr.dta"
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\rl_sugarcube_cortable.dta"
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_birthyear.dta"
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\sug_to_merge.dta"
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\NR_candyears.dta"
