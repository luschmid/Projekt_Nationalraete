clear
cap log close
set more 1
version 17

global path "C:/Schmidlu/Dropbox/Projekt Nationalräte"
global path_rl "C:/Schmidlu/Dropbox/Record Linkage"
*global path "E:/12. Cloud\Dropbox/Projekt Nationalräte"
*global path_rl "E:/12. Cloud\Dropbox/Record Linkage"

**************************************************************
* A) Read in output by coders completed by November 28, 2023
**************************************************************

* Note: The coers were Anna-Lena Beck, Aristomenis Zeus Kotatis, 
*		Jasper Ehrengruber, Julien Gaberthüel Michael Stern, 
*		Paul Vincent Heisenberg. 

* (i) Read in Excel files

forvalues i=1(1)27{
* a) Read in all batches
local j = `i' + 30  // all files receive a twin numbered with i+30
local k = `i' + 60  // all files receive a twin numbered with i+60
di "j:" "`j'"
di "k:" "`k'"
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`i'.xlsx", clear firstrow
rename Coding Coding1
recast str1000 firmnames, force
gen obs_id=_n
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`i'.dta", replace
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`j'.xlsx", clear firstrow
rename Coding Coding2
recast str1000 firmnames, force
gen obs_id=_n
*drop if obs_id==100
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`j'.dta", replace
import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`k'.xlsx", clear firstrow
rename Coding Coding3
recast str1000 firmnames, force
gen obs_id=_n
*drop if obs_id==100
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`k'.dta", replace
* b) Merge all files
* b1) Merge command for those between 30 and 57
use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`i'.dta", replace
merge 1:1 ID id_sug year_sug e_id_sug n_id_sug odd gap firstname lastname time_period firmnames obs_id using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`j'.dta",gen(mergevar)
* b2) Exit if not matched
quietly sum mergevar
if r(min)<3 {
	di "Not the same observations"
	sum xy
}
drop mergevar
* b3) Merge command for those between 60 and 87
merge 1:1 ID id_sug year_sug e_id_sug n_id_sug odd gap firstname lastname time_period firmnames obs_id using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`k'.dta",gen(mergevar)
* b4) Exit if not matched
quietly sum mergevar 
if r(min)<3 {
	di "Not the same observations"
	sum xy
}
* c) Generate one final coding variable
replace Coding1=0 if Coding1==.
replace Coding2=0 if Coding2==.
replace Coding3=0 if Coding3==.
gen Coding=Coding1 + Coding2 + Coding3
drop mergevar odd 
gen Group= `i'
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`i'_all.dta", replace
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`i'.dta"
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`j'.dta"
erase "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`k'.dta"
}



* (ii) Append all files 

use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_1_all.dta", clear

forvalues i=2(1)27{
append using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_`i'_all.dta"
}

tab Coding if Coding>0
keep if Coding>=2


* (iii) Keep only those with at least one coding, remove duplicates, and merge to 
*	   sorted out false positives in previous round (see Email Mark November 24, 2023)


rename ID id_0
rename id_sug id_1
rename year_sug year_1
rename e_id_sug e_cntr_w_1
rename n_id_sug n_cntr_w_1
tostring id_1, replace
duplicates drop id_0 id_1 year_1 e_cntr_w_1 n_cntr_w_1, force

preserve
use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1b_IDs_eliminatedFalsePositives.dta", clear
duplicates drop id_0 id_1 year_1 e_cntr_w_1 n_cntr_w_1, force
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1b_IDs_eliminatedFalsePositives_nodups.dta", replace
restore

merge 1:1 id_0 id_1 year_1 e_cntr_w_1 n_cntr_w_1 using  "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1b_IDs_eliminatedFalsePositives_nodups.dta"

gen ID=id_0

* (iv) Keep gap observations that were NOT sorted out as false positives in 
* 	   previous FP round (summer 2023)

preserve
keep if _merge==1
keep id_0 id_1 year_1 e_cntr_w_1 n_cntr_w_1
gen gap_consistent=1
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_consistent.dta", replace
restore

* (v) Keep gap observations that were sorted out as false positives in previous 
*	   FP round (summer 2023)

preserve
keep if _merge==3
keep  e_cntr_w_0 n_cntr_w_0  id_0 id_1 year_1 e_cntr_w_1 n_cntr_w_1
gen gap_inconsistent=1
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_inconsistent.dta", replace
restore

***********************************************************
* B) Prepare data for manual checks of inconsistent codings
***********************************************************

* (i) Prepare politician data for merge
* Note: One complication is that for the FN checks (gaps control) done in October
*		and November 2023, we identified political candidates at the ID level
*		and not at the ID-id_0-id_1 level as in entire machine learning step  
* 		and as during our FP checks. Since candidates may move (different geo 
*		locations), change firstnames and lastnames, and give different birthyears, 
*		we first need to prepare the politicians data to match them to the data
*		with the inconsistent coding. 

use "$path/02_Processed_data/02_Elections_1971_2015/nationalraete_1931_2015_Geo.dta",clear

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

* c) Candidacy years
preserve 
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

* (ii) Merge all datasets

* a) Start with Mark's correspondence table after removing FP

use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1b_fp_maj_IDs.dta", clear

* b) Append gaps (consistent and inconsistent)

append using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_consistent.dta"
append using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_inconsistent.dta"

* c) Flag all political candidates with at least one inconsistent gap 
* 	 observations and keep only those observations

replace gap_inconsistent=0 if gap_inconsistent==.
bysort id_0: egen gap_inconsistent_max=max(gap_inconsistent)

keep if gap_inconsistent_max==1


* (d) Merge Sugarcube data

rename e_cntr_w_1 e_id_sug
rename n_cntr_w_1 n_id_sug
rename year_1 year_sug
rename id_1 id_sug
rename e_cntr_w_0 e_id_polit
rename n_cntr_w_0 n_id_polit

destring id_sug, replace

gen ID_dupl=1

merge m:1 id_sug year_sug e_id_sug n_id_sug ID_dupl using "$path/02_Processed_data/10_Directors_1934_2003/Sugarcube_RLPostProcessing_Persons-Firmnames.dta", gen(merge_all)
keep if merge_all==3

rename e_id_sug   e_cntr_w_1 	
rename n_id_sug   n_cntr_w_1 	
rename year_sug   year_1 		
rename id_sug     id_1 		
rename e_id_polit e_cntr_w_0 	
rename n_id_polit n_cntr_w_0 	
rename lastname name_1
rename firstname firstname_1
rename gdename gdename_1
rename male male_1
rename ctn ctn_1
rename language_w language_w_1
rename firmnames firmnames_1

* e) Merge politician data

local nr_vars "firstname name job list male gdename birthyear"
foreach var of local nr_vars {
merge m:1 id_0 using  "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_`var'.dta", gen(merge_`var')
keep if merge_`var'==3
}

merge m:1 id_0 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\NR_candyears.dta"
keep if _merge==3


* f) Additional variables (distance, implicit age, time in office)

gen fp=""
gen comment=""
gen w_dist=sqrt((e_cntr_w_1- e_cntr_w_0)^2+(n_cntr_w_1-n_cntr_w_0)^2)/1000
destring year_1, replace
gen age_implicit=year_1-birthyear
sort id_0
keep if _merge == 3
drop _merge


* g) Indicate FP based on age

replace fp = "old" if age_implicit > 99  // cut-off according to Bisnode decision (see Email Lukas: 19.11.2022)
replace fp = "young" if age_implicit < 18  // cut-off according to Bisnode decision (see Email Lukas: 19.11.2022)

* h) Sorting (according to first 4 letters of firmnames)

egen id_0_num=group(id_0)
gen firmnames_1_first3=substr(firmnames_1,1,4)
sort id_0_num firmnames_1_first3 year_1 
gen odd=mod(id_0_num,2)
drop merge_*

* g) Ordering of variables

order odd link_score id_0 e_cntr_w_0 n_cntr_w_0 id_1 n_cntr_w_1 e_cntr_w_1 gap_consistent gap_inconsistent gap_inconsistent_max fp comment id_0_num firstname_0 name_0 firstname_1 name_1 year_1 w_dist candyears age_implicit  firmnames_1 male_0 male_1 gdename_0 gdename_1 birthyear_0

export excel using 	"$path\02_Processed_data\12_Record_linkage\02_Sugarcube\16_RL_Output_Round_4_Out\RL_Output_Round5_All.xlsx", first(var) replace

sum age_implicit if gap_inconsistent==1
hist age_implicit if gap_inconsistent==1


*********************************************
* C) Oufile after FP and FN checks (Dec 2023)
*********************************************

* a) Start with Mark's correspondence table after removing FP

use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1b_fp_maj_IDs.dta", clear

* b) Append gaps (consistent and inconsistent)

append using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_consistent.dta"
append using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_inconsistent.dta"

drop gap_consistent gap_inconsistent

save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\RL_Output_after_FP_and_FN.dta", replace





/* old
* (ii) Read in coder files

import excel using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\14_RL_Output_Round_3_Out\0_out_list.xlsx", /// 
 clear firstrow  cellrange("G2:U29") 
rename FirstBatch Batch1
rename Coder Coder1
rename SecondBatch Batch2
rename M Coder2
rename ThirdBatch Batch3
rename S Coder3
keep Coder* Batch*
gen Group=_n 
reshape long Coder Batch, i(Group) j(new)
gen work = subinstr(Batch, ",", "", .)
moss work, match("([0-9]+)")  regex
drop Batch
destring _match1, gen(Batch)
keep Coder Batch Group

save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\15_RL_Output_Round_3_In\RL_gaps_Round3_coder_overview.dta", replace
