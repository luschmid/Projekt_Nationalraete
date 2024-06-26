version 17
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\"


*** Panel ID: Prepare company info for procedures to establish Ground Truth

* Prepare Ground Truth for Person Panel ID (August 2023)

use "$path\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", clear
sort PID year ID_dupl gde_option mandID CID 
egen UPID = group(PID year ID_dupl gde_options)

preserve
keep UPID mandID firmname firmlocation
replace firmname = strproper(firmname)
replace firmname = usubinstr(firmname, "&" , " & ",.)
replace firmname = " " + firmname + " "
replace firmname = usubinstr(firmname, " Ag ", " AG ",.)
replace firmname = usubinstr(firmname, " Sa ", " SA ",.)
replace firmname = usubinstr(firmname, " Si ", " SI ",.)
replace firmname = usubinstr(firmname, "«" , "",.)
replace firmname = usubinstr(firmname, "»" , "",.)
replace firmname = usubinstr(firmname, "  " , " ",.)
replace firmname = usubinstr(firmname, "  " , " ",.)
replace firmname = ustrtrim(firmname)
replace firmlocation = strproper(firmlocation)
foreach ctn in Ag Ai Ar Be Bl Bs Fr Ge Gl Gr Ju Lu Ne Nw Ow Sg Sh So Sz Tg Ti Ur Vd Vs Zg Zh {
	replace firmlocation = usubinstr(firmlocation, "(`ctn')" , ustrupper("(`ctn')"),.)
}
gen firmnamelocation = firmname + " " + "(" + firmlocation + ")"
replace firmnamelocation = usubinstr(firmnamelocation, "  " , " ",.)
drop firmname firmlocation
sort UPID firmnamelocation
replace firmnamelocation = usubinstr(firmnamelocation, ";", "",.)
bysort UPID (firmnamelocation): gen firmstr = firmnamelocation[1]  // alphabetical order of firm names
by UPID: replace firmstr = firmnamelocation[_n-1] + "; " + firmnamelocation if _n > 1 
by UPID: replace firmstr = firmstr[_N]
bysort UPID: gen n = _n
keep if n == 1
drop n firmnamelocation mandID
rename firmstr firmnamelocation
save "$path\Sugarcube_Person_CompanyNames-Location_wide_alpha.dta", replace
restore

duplicates drop PID year ID_dupl gde_option, force // just keep one observation
drop CID compFunct mandate_flag firmname firmlocation owners capital function signature // these are variables for the long-format
sort UPID
merge 1:1 UPID using "$path\Sugarcube_Person_CompanyNames-Location_wide_alpha.dta"
keep PID year firstname lastname city gdename male GdeNr_E_CNTR GdeNr_N_CNTR ctn language_w firmnamelocation ID_dupl gde_options geo_merge

duplicates report PID year ID_dupl gde_options
order PID year ID_dupl gde_options
sort PID year ID_dupl gde_options
*save "$path\Sugarcube_PersPanelID_Firmnames_wide.dta", replace
export delimited using "$path\Sugarcube_PersPanelID_Firmnames_wide.csv", delimiter(",") replace

erase "$path\Sugarcube_Person_CompanyNames-Location_wide_alpha.dta"