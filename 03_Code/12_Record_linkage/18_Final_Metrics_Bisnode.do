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

save "$path\02_Processed_data\12_Record_linkage\01_Bisnode\nr_elected.dta", replace

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


*******************************
* B) Read in ground truth files
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

* f) Merge NR information

preserve
use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
keep ID birthyear
duplicates drop ID, force
rename ID id_0
save "$path\02_Processed_data\nr_id_birthyear.dta", replace
restore

merge m:1 id_0 using "$path\02_Processed_data\nr_id_birthyear.dta"
erase "$path\02_Processed_data\nr_id_birthyear.dta"
keep if _merge==3
drop _merge

* g) Merge Bisnode information

preserve
use personenid duns gremium rechtsform funktion eintrittdatum austrittdatum ///
	using "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Person-Firmen_Geo.dta", clear

gen entrydate = date(eintrittdatum, "YMD")
format entrydate %td
gen exitdate = date(austrittdatum, "YMD")
format exitdate %td
gen flag = 1 if exitdate<entrydate & entrydate != . & exitdate != .
gen exitdate_tmp = exitdate
replace exitdate = entrydate if flag == 1 // assumption: confused entry/exit date
replace entrydate = exitdate_tmp if flag == 1
drop exitdate_tmp flag
gen entryear=year(entrydate)
gen exityear=year(exitdate)

replace entryear = 1994 if eintrittdatum == ""
replace exityear = 2017 if austrittdatum == ""

rename personenid id_1
bysort id_1: egen entryear_min=min(entryear)
bysort id_1: egen exityear_max=max(exityear)
keep id_1 entryear_min exityear_max entryear exityear duns gremium rechtsform funktion

save "$path\02_Processed_data\bisnode_id_mandates.dta", replace
restore

merge 1:m id_1 using "$path\02_Processed_data\bisnode_id_mandates.dta",gen(merge_id1)
erase "$path\02_Processed_data\bisnode_id_mandates.dta"
keep if merge_id1==3

* h) Keep only Verwaltungsräte and apply age restrictions

gen age_firstmandate=entryear_min-birthyear
gen age_lastmandate=exityear_max-birthyear

drop if age_firstmandate < 18  // see Email Lukas: 19.11.2022
drop if age_lastmandate < 18
drop if age_firstmandate > 85
drop if age_lastmandate > 99	

keep if gremium == "Verwaltungsrat" & rechtsform == "Aktiengesellschaft"
drop if inlist(funktion, "Aktuar/in (nicht Mitglied)", "Sekretär/in (nicht Mitglied)", "Beisitzer/in", "Ausseramtliche/r Konkursverwalter/in", "Generalsekretär/in (nicht Mitglied)", "Liquidator/in", "Kassier/in (nicht Mitglied)", "Protokollführer/in (nicht Mitglied)", "Verwalter/in (nicht Mitglied)")

* i) Construct dataset at id_0-id_1-mandate-year level 

bysort id_1 duns: gen mandid = _n
expand 25, gen(expy)
bysort id_1 duns mandid: gen year=_n + 1993
drop expy

drop if year<entryear | year>exityear
drop if year>2017

duplicates drop id_0 id_1 duns year, force
keep id_0 id_1 duns year

save "$path\02_Processed_data\12_Record_linkage\01_Bisnode\bisnode_test_set_final.dta", replace


*********************************
* C) Read in record linkage files
*********************************

* (i) Read in RL output after postprocessing
* Note: Split in three parts: (i) perfect matches (no post-processing), 
* 		(ii) optimal cutoff, (iii) cutoff of 0. For (ii), we look at all *
* 		connections. For (iii), we focus only on connections within 
*       rechtsform=="Aktiengesellschaft" & gremium=="Verwaltungsrat". 
* 		Important: This dataset is generated in 05_Setup_DataAnalysis.do. 
*		Corrections for only Verwaltungsräte and age are already implemented. 

use "$path\02_Processed_data\11_Directors_1994_2018\RL_G7_NR-Bisnode_1994-2017.dta", clear
rename NRid id_0
rename VRid id_1 

* (ii) Keep only id_0 that are in test dataset

merge m:1 id_0 using "$path\02_Processed_data\12_Record_linkage\01_Bisnode\bisnode_test_set_polids.dta"
keep if _merge==3
drop _merge

* (iii) Keep only those observations with a mandate (drop true negatives)

unique id_0 if id_1!=.

keep if id_1!=.

keep id_0 id_1 duns year

save "$path\02_Processed_data\12_Record_linkage\01_Bisnode\rl_test_set.dta", replace


*****************************************************************************
* D) Merge final RL output (after FP and FN checks) to ground truth files and 
*	 generate metrics															
*****************************************************************************

* (i) Test set at link level

* a) Read in test set

use "$path\02_Processed_data\12_Record_linkage\01_Bisnode\bisnode_test_set_final.dta", clear

* b) Merge RL output

merge 1:1 id_0 id_1 duns year using "$path\02_Processed_data\12_Record_linkage\01_Bisnode\rl_test_set.dta", gen(merge_testset)

* c) Recode TP, FP, and FN

gen TP=0 
replace TP=1 if merge_testset==3
gen FP=0 
replace FP=1 if merge_testset==2
gen FN=0 
replace FN=1 if merge_testset==1

* d) Output for entire test set and by elected status

preserve
collapse (sum) TP FP FN
gen precision=TP/(FP+TP)
gen recall=TP/(FN+TP)
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum  TP FP FN precision recall F1
restore 


