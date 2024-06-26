version 17
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\"


*** RL Post-Processing: Prepare company info for post-processing procedures  

* Prepare dataset for RL post-processing and manual FP cleaning (Nationalrat-VR RL)  -- > prepared for Sandro Cilurzo May 2023

use "$path\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", clear
sort PID year ID_dupl gde_option mandID CID 
egen UPID = group(PID year ID_dupl gde_options)
bysort UPID: gen n = _n

preserve
keep UPID mandID firmname
replace firmname = strproper(firmname)
replace firmname = usubinstr(firmname, "Ag", "AG",.)
replace firmname = usubinstr(firmname, "Sa", "SA",.)
replace firmname = usubinstr(firmname, "Si", "SI",.)
sort UPID firmname
replace firmname = usubinstr(firmname, ";", "",.)
bysort UPID (firmname): gen firmstr = firmname[1]  // alphabetical order of firm names
by UPID: replace firmstr = firmstr[_n-1] + "; " + firmname if _n > 1 
by UPID: replace firmstr = firmstr[_N]
bysort UPID: gen n = _n
keep if n == 1
drop n firmname
rename firmstr firmnames
save "$path\Sugarcube_Person_Company_Names_wide_alpha.dta", replace
restore

drop if n != 1
drop CID compFunct mandate_flag firmname firmlocation owners capital function signature n // these are variables for the long-format
sort UPID
merge 1:1 UPID using "$path\Sugarcube_Person_Company_Names_wide_alpha.dta"
keep id_sug year_sug firstname lastname gdename male GdeNr_E_CNTR GdeNr_N_CNTR ctn language_w e_id_sug n_id_sug firmnames ID_dupl gde_options geo_merge

duplicates report id_sug year_sug e_id_sug n_id_sug
order id_sug year_sug e_id_sug n_id_sug
sort id_sug year_sug e_id_sug n_id_sug ID_dupl
save "$path\Sugarcube_RLPostProcessing_Persons-Firmnames.dta", replace
export delimited using "$path\Sugarcube_RLPostProcessing_Persons-Firmnames.csv", delimiter(",") replace

*erase "$path\Sugarcube_Person_Company_Names_wide_alpha.dta"


/*** Original: had error in dropping command. It used keep if mandID == 1 (!= dropping the first entry!)

keep if mandID == 1
drop CID compFunct mandate_flag firmname firmlocation owners capital function signature // these are variables for the long-format
sort UPID
merge 1:1 UPID using "$path\Sugarcube_Person_Company_Names_wide_alpha.dta"
drop if _merge == 2
keep id_sug year_sug firstname lastname gdename male GdeNr_E_CNTR GdeNr_N_CNTR ctn language_w e_id_sug n_id_sug firmnames ID_dupl gde_options geo_merge

duplicates report id_sug year_sug e_id_sug n_id_sug
order id_sug year_sug e_id_sug n_id_sug
sort id_sug year_sug e_id_sug n_id_sug ID_dupl
save "$path\Sugarcube_RLPostProcessing_Persons-Firmnames_old.dta", replace
export delimited using "$path\Sugarcube_RLPostProcessing_Persons-Firmnames.csv", delimiter(",") replace
