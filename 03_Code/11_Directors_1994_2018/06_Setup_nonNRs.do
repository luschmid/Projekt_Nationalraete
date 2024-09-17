clear
cap log close
set more 1
version 17

global dataNR "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\"
global dataRL "E:\12. Cloud\Dropbox\Record Linkage\02_Data\"
global tempBis "C:\Users\schelkem\BisnodeSetup_tmp\"  

*global dataNR "C:\Current\Dropbox\Projekt Nationalräte\02_Processed_data\"
*global dataRL "C:\Current\Dropbox\Record Linkage\02_Data\"
*global tempBis "C:\BisnodeSetup_tmp\"

mkdir "C:\Users\schelkem\BisnodeSetup_tmp"

********************************************************************************
*** Prepare Bisnode dataset for summary statistics in event-time (elections) ***
********************************************************************************

*** Steps
* Step 1: Get information on which Bisnode-IDs are matched to NR candidates
* Step 2: Eliminate all NR-candidate matches from Bisnode database
* Step 3: Collapse at ID level


*** Step 1: Read-in Bisnode_RL dataset, identify Bisnode-ID receiving a match with NR data 
use "$dataNR\11_Directors_1994_2018\RL_G7_NR-Bisnode_1994-2017.dta", clear
keep VRid 
duplicates drop VRid, force
drop if VRid == .
save "$tempBis\RL_bisnode_ids.dta", replace


* Step 2: Drop all Bisnode IDs with link to NR candidates
use "$dataNR\11_Directors_1994_2018\Bisnode_Person-Firmen_Geo.dta", clear
rename personenid VRid
merge m:1 VRid using "$tempBis\RL_bisnode_ids.dta"
drop if _merge == 3

// drop all non-AG observations
keep if gremium == "Verwaltungsrat" & rechtsform == "Aktiengesellschaft"
drop if inlist(funktion, "Aktuar/in (nicht Mitglied)", "Sekretär/in (nicht Mitglied)", "Beisitzer/in", "Ausseramtliche/r Konkursverwalter/in", "Generalsekretär/in (nicht Mitglied)", "Liquidator/in", "Kassier/in (nicht Mitglied)", "Protokollführer/in (nicht Mitglied)", "Verwalter/in (nicht Mitglied)")

keep VRid duns funktion eintrittdatum austrittdatum kapitalnominalag

// restructure from spell-level to panel (as in "05_Setup_DataAnalysis.do")
bysort VRid duns funktion: gen mandid = _n
expand 25, gen(expy)
bysort VRid duns funktion mandid: gen year=_n + 1993
drop expy
sort VRid duns funktion mandid year
order year VRid duns funktion mandid

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
tab entryear
tab exityear
replace entryear = 1993 if eintrittdatum == ""
replace exityear = 2019 if austrittdatum == ""
tab entryear
tab exityear
drop if year<entryear | year>exityear

duplicates report VRid duns funktion year 
duplicates drop VRid duns funktion year, force  // drop "true" duplicates ("funktion" included)

* Overlapping spells as President and Member --- > keep highest function (e.g., President)
bysort VRid duns year: gen no_duplMand = _n  
bysort VRid duns year: egen avg_duplMand = mean(no_duplMand)
gen duplMand = 0
bysort VRid duns year: replace duplMand = 1 if avg_duplMand>1
drop avg_duplMand

gen MandH = 1 if funktion=="Präsident/in"
replace MandH = 2 if funktion=="Vizepräsident/in"
replace MandH = 3 if funktion=="Delegierte/r"
replace MandH = 4 if funktion=="Mitglied"
replace MandH = 5 if funktion=="Einziges Mitglied"
replace MandH = 6 if funktion=="Verwalter/in"
replace MandH = 7 if funktion=="Aktuar/in"
replace MandH = 8 if funktion=="Sekretär/in"
bysort VRid duns year: egen minMandH = min(MandH)
keep if MandH == minMandH  // keep the "highest" function
duplicates report VRid duns year 
drop duplMand no_duplMand MandH minMandH

sort VRid year
keep if year>=1994 & year<2018
rename VRid PID
rename duns CID
keep PID year CID kapitalnominalag 


* Step 3: Import capital threshold, collapse at ID level
** import capital thresholds
merge m:1 year using "$dataNR\11_Directors_1994_2018\bisnode_cap90.dta"
keep if _merge == 3
drop _merge
gen all = 1 
gen lrg = 1 if CID != . & kapitalnominalag > cap90  & kapitalnominalag < . 
gen sml = 1 if CID != . & kapitalnominalag <= cap90  & kapitalnominalag < .

** collapse
collapse (sum) n_all_sum=all n_lrg_sum=lrg n_sml_sum=sml, by(PID year)  

foreach var in all lrg sml { 
	gen i_`var' = .
	replace i_`var' = 0 if n_`var'_sum == 0
	replace i_`var' = 1 if n_`var'_sum > 0 & n_`var'_sum < .
}
foreach var of varlist n_all_sum-i_sml { 
	rename `var' `var'_c1  // version with sugarcube from 1934-2003, bisnode afterwards
}
keep if year > 2003
sort year PID 
save "$dataNR\11_Directors_1994_2018\Bisnode_nonNRs.dta", replace


rmdir "$tempBis"

