clear
clear matrix
clear mata
cap log close
set more 1


**** Recordlinkage: Checking (wrong?) False Positives which seem to be suspicious according to Sandro (email 02.01.2023) ****

global data "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\12_Record_linkage\02_Sugarcube\"


*** import wrong "false positives" from Sandros checks after Record Linkage (Email 02.01.2023)
// what exactly is this?
// why is year and birthyear identical?
// why are there duplicates?
// why do 459 out of the remaining 674 match with obs in GT? Wasn't the problem that they do NOT match ("wrong" FPs)?

import delimited using  "${data}wrong_fps_100518bc-12f0-4d5d-8b42-1f4484952a2a.csv", encoding("utf-8") clear
rename id id_sug
rename year year_sug
rename e_id e_id_sug
rename n_id n_id_sug
rename name name_wfps
rename firstname firstname_wfps
rename birthyear birthyear_wfps
rename sex sex_wfps
rename gdename gdename_wfps
replace name_wfps = ustrlower(name_wfps)
replace firstname_wfps = ustrlower(firstname_wfps)
keep id_sug year_sug e_id_sug n_id_sug *_wfps
sort id_sug year_sug e_id_sug n_id_sug
sum id_sug
duplicates report id_sug year_sug e_id_sug n_id_sug
duplicates drop id_sug year_sug e_id_sug n_id_sug, force // why are there duplicates?
sum id_sug
save "${data}wrong_fps_100518bc-12f0-4d5d-8b42-1f4484952a2a.dta", replace

*** Compare to Ground truth delivered to Sandro

use "${data}Sugarcube_Ground_Truth_210714.dta", clear
sum finalmatch

// renaming of variables for "Sugarcube Ground Thruth ID-file" for Sandro
rename ID_Pol id_polit
rename Geo_Code_E_Pol e_id_polit
rename Geo_Code_N_Pol n_id_polit
rename ID_VR id_sug
rename Jahr_VR year_sug
rename Geo_Code_E_VR e_id_sug 
rename Geo_Code_N_VR n_id_sug

order id_polit e_id_polit n_id_polit id_sug year_sug e_id_sug n_id_sug

*** merge

merge m:1 id_sug year_sug e_id_sug n_id_sug using "${data}wrong_fps_100518bc-12f0-4d5d-8b42-1f4484952a2a.dta"

drop if _merge == 1 // focus on problem cases _merge == 2 (and compare with _merge == 3)
sort firstname_wfps name_wfps e_id_sug n_id_sug year
br Vorname_Pol Nachname_Pol Vorname_VR Nachname_VR Alter firstname_wfps name_wfps e_id_sug n_id_sug _merge year_sug
order Vorname_Pol Nachname_Pol Vorname_VR Nachname_VR Alter firstname_wfps name_wfps e_id_sug n_id_sug _merge year_sug
keep Vorname_Pol Nachname_Pol Vorname_VR Nachname_VR Alter firstname_wfps name_wfps e_id_sug n_id_sug _merge year_sug id_sug
export excel using "${data}wrong_fps_merge2-3_temp.xlsx", firstrow(variables) replace

*** Look up non-matches: compare to initial sample taken by Lukas, before Students worked in GT
/*
Comparisons made (manually):
 - against Lukas' initial sample
 - against "Sugar_Checks1_final" sample (== Ground Truth before dropping non-matches)
 - against "Sugarcube_Ground_Truth"
*/

* Restore Lukas' initial match NR-VR match (some fuzzy-merge procedure)
preserve
forv v = 1(1)19 {
		import delimited using "${data}02_Nationalrat_VR_Initial_Match\Sugar_Checks`v'.csv", encoding("utf-8") clear
		save "${data}02_Nationalrat_VR_Initial_Match\Sugar_Checks_NR_initialMatch`v'.dta", replace
	}
use "${data}02_Nationalrat_VR_Initial_Match\Sugar_Checks_NR_initialMatch1.dta", clear
forv v = 2(1)19 {
	append using "${data}02_Nationalrat_VR_Initial_Match\Sugar_Checks_NR_initialMatch`v'.dta", force
	}	
save "${data}02_Nationalrat_VR_Initial_Match\Sugar_Checks_NR_initialMatch.dta", replace
forv v = 1(1)19 {
	erase "${data}02_Nationalrat_VR_Initial_Match\Sugar_Checks_NR_initialMatch`v'.dta"
	}
restore

** Example matches/non-matches:
*  What are the reasons for the supposed (according to Sandro) "wrong" false positives (269 obs)?
*  (1) Not a False Positive: very likely other candidate with the same first-/lastname in same location: 168 obs
*  (2) The observation was not chosen in the initial sampling by Lukas: 40 obs
*  (3) wrong False Positive: seems to be error in GT coding: 7 obs
* see "${data}wrong_fps_merge2-3_checkMS.xlsx"


erase "${data}02_Nationalrat_VR_Initial_Match\Sugar_Checks_NR_initialMatch.dta"
erase "${data}wrong_fps_merge2-3_temp.xlsx"



