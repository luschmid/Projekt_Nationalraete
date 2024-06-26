version 17
clear all
set more off
cap log close


global raw "E:\12. Cloud\Dropbox\Projekt Nationalräte\01_Raw_data\18_Interessenbindungen\IG_1995-2012_MarcoPortmann\"
global processed "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\18_Interessenbindungen\"

import excel "${raw}tab_affiliations.xlsx", first clear
destring Mandate, replace

sort mp_id year

preserve
import excel "${raw}tab_mitglieder.xlsx", first clear
sort mp_id
save "${processed}tab_mitglieder.dta", replace
restore

merge m:1 mp_id using "${processed}tab_mitglieder.dta"

*br * if _merge == 2 // (empty inst_id's, no value in keeping those)
drop if _merge == 2
drop _merge

sort inst_id

preserve
import excel "${raw}tab_inst.xlsx", first clear
sort inst_id
save "${processed}tab_inst.dta", replace
restore

merge m:1 inst_id using "${processed}tab_inst.dta"

br * if _merge == 1
drop _merge
tab year

sort mp_id year
save "${processed}Interessenbindungen_1996-2014.dta", replace

erase "${processed}tab_inst.dta"
erase "${processed}tab_mitglieder.dta"