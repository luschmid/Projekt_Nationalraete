clear
cap log close
set more 1
version 17

global path "C:/Schmidlu/Dropbox/Projekt Nationalräte"
global path_rl "C:/Schmidlu/Dropbox/Record Linkage"
*global path "E:/12. Cloud\Dropbox/Projekt Nationalräte"
*global path_rl "E:/12. Cloud\Dropbox/Record Linkage"

*******************************
* A) Prepare election data
*******************************

* (i) NR elected information

use "$path\02_Processed_data\nationalraete_1931_2015.dta",clear

collapse (max) elected, by(ID)
rename ID id_0

save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_elected.dta", replace

* (ii) NR panel in wide format

use "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta",clear
	
local vars "name firstname birthyear job gdename"
* Note: Create a file with all different name/firstname spellings and birthyears
* which we then append to original dataset.

foreach var in `vars'{	
preserve
bysort ID `var': gen indi_`var'=_n		
keep if indi_`var'==1
* Note: Keep only one observation per ID-var group.

bysort ID: gen observation_number =_n
keep ID `var' observation_number
reshape wide `var', i(ID) j(observation_number) 
save "$path\02_Processed_data\01_Elections_1931_1975\\`var'", replace
restore
}

keep ID year 

merge m:1 ID using "$path\02_Processed_data\01_Elections_1931_1975\name.dta", ///
	gen(merge_firstname)
merge m:1 ID using "$path\02_Processed_data\01_Elections_1931_1975\firstname.dta", ///
	gen(merge_name)
merge m:1 ID using "$path\02_Processed_data\01_Elections_1931_1975\birthyear.dta", ///
	gen(merge_birthyear)
merge m:1 ID using "$path\02_Processed_data\01_Elections_1931_1975\job.dta", ///
	gen(merge_job)
merge m:1 ID ///
	using "$path\02_Processed_data\01_Elections_1931_1975\gdename.dta", ///
	gen(merge_municipality)
	
	
gen name=name1
replace name=name+ ", " + name2 if name2!=""
replace name=name+ ", " + name3 if name3!=""

gen firstname=firstname1
replace firstname=firstname+ ", " + firstname2 if firstname2!=""
replace firstname=firstname+ ", " + firstname3 if firstname3!=""

tostring birthyear1, replace
tostring birthyear2, replace

gen birthyear=birthyear1
replace birthyear=birthyear+ ", " + birthyear2 if birthyear2!="."

gen gdename=gdename1
replace gdename=gdename+ ", " + gdename2 if gdename2!=""
replace gdename=gdename+ ", " + gdename3 if gdename3!=""
replace gdename=gdename+ ", " + gdename4 if gdename4!=""

keep ID year name firstname birthyear gdename

duplicates drop ID, force
rename ID id_0

erase "$path\02_Processed_data\01_Elections_1931_1975\name.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\firstname.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\birthyear.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\job.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\gdename.dta"

save "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_wide.dta", replace

**********************
* B) Read in Sugarcube
**********************
/*
use PID year year_sug firstname lastname gdename CID firmname using "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", clear
duplicates drop PID year_sug firstname lastname gdename CID firmname, force
compress
rename PID id_1
rename year_sug year_1
label var year_1 "=year_sug: ONLY ID-Variable for RL purpose"
rename firstname firstname_1
rename lastname name_1
rename gdename gdename_1
tostring id_1, replace
save "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person_CleanName-Gender-Geo_Companies_small.dta", replace
*/

*******************************
* C) Read in ground truth files
*******************************

insheet using ///
	"$path\02_Processed_data\12_Record_linkage\01_Bisnode\test_set_randomized_buckets_0.65.csv", ///
	delim(",") clear
* a) Rename variables

rename id id_0 
rename personenid id_1

* b) Save all politician IDs in test set
preserve
duplicates drop id_0, force
keep id_0
save "$path\02_Processed_data\12_Record_linkage\01_Bisnode\bisnode_test_set_polids.dta", replace
restore
* c) Drop all those with no match in GT data
drop if id_1==.

* d) Check for duplicates
duplicates report id_0 id_1 

* e) Keep only relevant variables
keep id_0 id_1 
save "$path\02_Processed_data\12_Record_linkage\01_Bisnode\bisnode_test_set_randomized_buckets_0.65.dta", replace



*********************************
* D) Read in record linkage files
*********************************

* (i) Read in RL output after postprocessing

import delim using "$path/02_Processed_data/12_Record_linkage/01_Bisnode/RL-Results-G7-perfectmatches.csv", delimiters(",") encoding("UTF-8") clear
gen source = "perfect"
preserve
import delim using "$path/02_Processed_data/12_Record_linkage/01_Bisnode/04_Post_RL_QualChecks/bisnode_fp_firstcheck_Generation7_optimal_corrected.csv", delimiters(",") encoding("UTF-8") clear
gen source = "optimal"
save "$path/02_Processed_data/12_Record_linkage/01_Bisnode/bisnode_fp_firstcheck_Generation7_optimal_corrected.dta", replace
restore
preserve
import delim using "$path/02_Processed_data/12_Record_linkage/01_Bisnode/04_Post_RL_QualChecks/bisnode_fp_firstcheck_Generation7_0_corrected.csv", delimiters(",") encoding("UTF-8") clear
gen source = "0"
save "$path/02_Processed_data/12_Record_linkage/01_Bisnode/bisnode_fp_firstcheck_Generation7_0_corrected.dta", replace
restore

append using "$path/02_Processed_data/12_Record_linkage/01_Bisnode/bisnode_fp_firstcheck_Generation7_optimal_corrected.dta"
append using "$path/02_Processed_data/12_Record_linkage/01_Bisnode/bisnode_fp_firstcheck_Generation7_0_corrected.dta"

drop if codierungfin == 1 // False positives
drop if age_firstmandate < 18  // see Email Lukas: 19.11.2022
drop if age_lastmandate < 18
drop if age_firstmandate > 85
drop if age_lastmandate > 99	

keep id_0 id_1
duplicates report id_0 id_1   // duplicates because  manual post-processing was at spell-level, not  personID-level (i.e. various mandates per person)
duplicates drop id_0 id_1, force


save "$path\02_Processed_data\12_Record_linkage\01_Bisnode\RL_Output_after_FP_and_FN_no_duplicates_temp.dta", replace

erase "$path\02_Processed_data\12_Record_linkage\01_Bisnode\bisnode_fp_firstcheck_Generation7_0_corrected.dta"
erase "$path\02_Processed_data\12_Record_linkage\01_Bisnode\bisnode_fp_firstcheck_Generation7_optimal_corrected.dta"



*****************************************************************************
* E) Merge final RL output (after FP and FN checks) to ground truth files and 
*	 generate metrics															
*****************************************************************************

* (i) Test set at link level

* a) Read in test set

use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\sugarcube_test_set_randomized_buckets_0.65.dta", clear

tostring id_1, replace

* b) Merge RL output

merge 1:1 id_0 id_1 year_1 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\RL_Output_after_FP_and_FN_no_duplicates_temp.dta", gen(merge_testset)

*merge 1:1 id_0 id_1 year_1 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\rl_sugarcube_cortable.dta", gen(merge_testset)


* c) Merge politician IDs in test set and drop all those not in test set

merge m:1 id_0 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\sugarcube_test_set_polids.dta"

keep if _merge==3

* d) Merge politician ID with coding for elected (at least once)

merge m:1 id_0 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\nr_elected.dta", gen(merge_elected)

drop if merge_elected==2

* e) Merge data about when politicians were in office

merge m:1 id_0 year_1 using "$path\02_Processed_data\02_Elections_1971_2015\OfficeYears_GT1.dta", gen(merge_offyears)
drop if merge_offyears == 2

gen in_office=0 if elected==1
replace in_office=1 if merge_offyears==3 & elected==1

* f) Recode TP, FP, and FN

gen TP=0 
replace TP=1 if merge_testset==3
gen FP=0 
replace FP=1 if merge_testset==2
gen FN=0 
replace FN=1 if merge_testset==1

* g) Output for entire test set and by elected status

preserve
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1  F1 TP FP FN
restore 

preserve
keep if in_office==1
tab in_office,missing
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1  F1 TP FP FN
restore 

preserve
keep if in_office!=1
tab in_office,missing
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1  F1 TP FP FN
restore 

preserve
keep if elected==1
tab elected,missing
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall  F1 TP FP FN
restore

preserve
keep if elected==0
tab elected,missing
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1  F1 TP FP FN
restore

* h) Collapse TP, FP, and FN by politician ID and year for RDD test

gen year=year_1
replace year=year_1+1 if year_1<=1963
* Note: replace year_sug (for GT and RL) with true years (June 24). 

preserve
collapse (sum) TP FP FN , by(id_0 year)

gen prc_lk=TP/(FP+TP)
gen rcl_lk=TP/(FN+TP)
gen f1_lk = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
rename id_0 ID
keep ID year prc_lk rcl_lk f1_lk

save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\pr_link.dta", replace
restore

sum prc_lk if year==1939 | year==1943 | year==1947 | year==1951 | year==1955 | /// 
	year==1959 | year==1963 | year==1967 | year==1971 | year==1975 | ///
	year==1979 | year==1983 | year==1987 | year==1991 | ///
	year==1995 | year==1999 | year==2003


* (ii) Test set at mandate level

* a) Merge Sugarcube and politician data for visual inspection

merge 1:m id_1 year_1 using "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person_CleanName-Gender-Geo_Companies_small.dta", gen(merge_sug)
drop if merge_sug == 2

merge m:1 id_0 using "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_wide.dta", gen(merge_nr)
drop if merge_nr == 2

sort id_0 year_1 id_1 firmname

order id_0 year_1 id_1 firstname firstname_1  name name_1  gdename  gdename_1 firmname TP FP FN in_office
br id_0 year_1 id_1 firstname firstname_1  name name_1  gdename  gdename_1 firmname TP FP FN in_office

* b) Metrics at mandate level (comparable to GT 2)

preserve
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1  TP FP FN
restore 

preserve
keep if in_office==1
tab in_office,missing
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1 TP FP FN
restore 

preserve
keep if in_office!=1
tab in_office,missing
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1 TP FP FN
restore 

preserve
keep if elected==1
tab elected,missing
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1  TP FP FN
restore

preserve
keep if elected==0
tab elected,missing
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1 TP FP FN
restore

* c) Collapse TP, FP, and FN by politician ID and year for RDD test

collapse (sum) TP FP FN , by(id_0 year)

gen prc_mn=TP/(FP+TP)
gen rcl_mn=TP/(FN+TP)
gen f1_mn = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
rename id_0 ID
keep ID year prc_mn rcl_mn f1_mn

save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\pr_mandate.dta", replace

