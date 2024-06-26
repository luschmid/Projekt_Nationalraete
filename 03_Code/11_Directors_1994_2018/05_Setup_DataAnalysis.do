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

mkdir "C:\Users\schelkem\BisnodeSetup_tmp", public  // reduce traffic on dropbox: create directory for all the tmp files.
* mkdir "C:\BisnodeSetup_tmp", public  // reduce traffic on dropbox: create directory for all the tmp files.


**** Setup Bisnode data for analyses ****
* Step 1: read-in Generation 7 after post-processing (map of IDs)
* Step 2: merge Bisnode data on Bisnode ID
* Step 3: establish panel structure
* Step 4: merge NR data on NR ID
* Step 5: collapse data to create dataset for analysis


*** Step 1: Read-in Generation 7 after post-processing (map of IDs) 
import delim using "$dataNR\12_Record_linkage\01_Bisnode\RL-Results-G7-perfectmatches.csv", delimiters(",") encoding("UTF-8") clear
gen source = "perfect"
preserve
import delim using "$dataNR\12_Record_linkage\01_Bisnode\04_Post_RL_QualChecks\bisnode_fp_firstcheck_Generation7_optimal_corrected.csv", delimiters(",") encoding("UTF-8") clear
gen source = "optimal"
save "$tempBis\bisnode_fp_firstcheck_Generation7_optimal_corrected.dta", replace
restore
preserve
import delim using "$dataNR\12_Record_linkage\01_Bisnode\04_Post_RL_QualChecks\bisnode_fp_firstcheck_Generation7_0_corrected.csv", delimiters(",") encoding("UTF-8") clear
gen source = "0"
save "$tempBis\bisnode_fp_firstcheck_Generation7_0_corrected.dta", replace
restore

append using "$tempBis\bisnode_fp_firstcheck_Generation7_optimal_corrected.dta"
append using "$tempBis\bisnode_fp_firstcheck_Generation7_0_corrected.dta"

drop if codierungfin == 1 // False positives
drop if age_firstmandate < 18  // see Email Lukas: 19.11.2022
drop if age_lastmandate < 18
drop if age_firstmandate > 85
drop if age_lastmandate > 99	

keep id_0 id_1
duplicates report id_0 id_1   // duplicates because  manual post-processing was at spell-level, not  personID-level (i.e. various mandates per person)
duplicates drop id_0 id_1, force

duplicates report id_1 // create dataset unique in bisnode ID
sort id_1
bysort id_1: gen n = _n
reshape wide id_0, i(id_1) j(n)
rename id_1 VRid
rename id_01 NRid_1
rename id_02 NRid_2
rename id_03 NRid_3
rename id_04 NRid_4
gen NRlinks = 1  // some Bisnode IDs have up to 4 NR links
replace NRlinks = 2 if NRid_2 != "" & NRid_3 == "" & NRid_4 == ""
replace NRlinks = 3 if NRid_2 != "" & NRid_3 != "" & NRid_4 == ""
replace NRlinks = 4 if NRid_2 != "" & NRid_3 != "" & NRid_4 != ""
save "$tempBis\RL_bisnodeIDs_NRids.dta", replace

*** Step 2: merge Bisnode data on Bisnode IDs from RL
use "$dataNR\11_Directors_1994_2018\Bisnode_Person-Firmen_Geo.dta", clear
rename personenid VRid
merge m:1 VRid using "$tempBis\RL_bisnodeIDs_NRids.dta"
keep if _merge == 3
drop _merge

keep VRid duns vorname nachname W_PLZ4 wohnort wohnland nationalität geburtstag gremium rechtsform funktion eintrittdatum austrittdatum E_CNTR_w N_CNTR_w uid handelsregisternummer firma handelsname plzdomiziladresse ortdomiziladresse kantondomiziladresse nogacode unterschrift kapitalnominalag kapitaleinbezahltag NRid_* NRlinks
save "$tempBis\RL_G7_Bisnode_Person-Firmen_Geo.dta", replace

*** Step 3: Panel structure
*use "$tempBis\RL_G7_Bisnode_Person-Firmen_Geo.dta", clear
duplicates report VRid duns handelsregisternummer rechtsform gremium funktion eintrittdatum austrittdatum
duplicates report VRid duns handelsregisternummer rechtsform gremium funktion eintrittdatum austrittdatum kapitalnominalag kapitaleinbezahltag plzdomiziladresse firma // (NRid_* NRlinks ) it looks as if these are genuine duplicates in Bisnode 
*duplicates tag VRid duns handelsregisternummer rechtsform gremium funktion eintrittdatum austrittdatum kapitalnominalag kapitaleinbezahltag plzdomiziladresse firma, gen(dupl) // visual inspection  
*br * if dupl>0   // visual inspection suggests genuine duplicates
*drop dupl
duplicates drop VRid duns handelsregisternummer rechtsform gremium funktion eintrittdatum austrittdatum kapitalnominalag kapitaleinbezahltag plzdomiziladresse firma, force


* keep only AG's and VR mandates (post-processing only considered those)
keep if gremium == "Verwaltungsrat" & rechtsform == "Aktiengesellschaft"
drop if inlist(funktion, "Aktuar/in (nicht Mitglied)", "Sekretär/in (nicht Mitglied)", "Beisitzer/in", "Ausseramtliche/r Konkursverwalter/in", "Generalsekretär/in (nicht Mitglied)", "Liquidator/in", "Kassier/in (nicht Mitglied)", "Protokollführer/in (nicht Mitglied)", "Verwalter/in (nicht Mitglied)")

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

duplicates report VRid duns rechtsform gremium funktion year 
duplicates drop VRid duns rechtsform gremium funktion year, force  // drop "true" duplicates ("funktion" included)

* Overlapping spells as President and Member --- > keep highest function (e.g., President)
bysort VRid duns rechtsform gremium year: gen no_duplMand = _n  
bysort VRid duns rechtsform gremium year: egen avg_duplMand = mean(no_duplMand)
gen duplMand = 0
bysort VRid duns rechtsform gremium year: replace duplMand = 1 if avg_duplMand>1
drop avg_duplMand
/*
Funktion			Freq.	Percent	Cum.
			
Aktuar/in			36		0.06	0.06
Delegierte/r		1,620	2.71	2.77
Einziges Mitglied	9,213	15.42	18.19
Mitglied			27,348	45.77	63.96
Präsident/in		16,239	27.18	91.14
Sekretär/in			1,847	3.09	94.23
Verwalter/in		25		0.04	94.28
Vizepräsident/in	3,420	5.72	100.00
*/
gen MandH = 1 if funktion=="Präsident/in"
replace MandH = 2 if funktion=="Vizepräsident/in"
replace MandH = 3 if funktion=="Delegierte/r"
replace MandH = 4 if funktion=="Mitglied"
replace MandH = 5 if funktion=="Einziges Mitglied"
replace MandH = 6 if funktion=="Verwalter/in"
replace MandH = 7 if funktion=="Aktuar/in"
replace MandH = 8 if funktion=="Sekretär/in"
bysort VRid duns rechtsform gremium year: egen minMandH = min(MandH)
keep if MandH == minMandH  // keep the "highest" function
duplicates report VRid duns rechtsform gremium year 
drop duplMand no_duplMand MandH minMandH

sort VRid year
rename vorname firstname_bis 
rename nachname name_bis
rename wohnort residence_bis
keep if year>1993 & year<2019

save "$tempBis\RL_Bisnode-NRids_1994-2018_tmp.dta", replace

*** Step 4: merge NR info on NRids
* prepare NR dataset
clear
set obs 105 // create year dataset 1915 - 2020
gen n = _n-1
gen year = 1915 + n
drop n
save "$tempBis\years.dta", replace

use "$dataNR\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta", clear  // create NRid dataset
keep ID
duplicates drop ID, force
save "$tempBis\NRids.dta", replace
 
*use "$dataNR\11_Directors_1994_2018\NRids.dta", clear
cross using "$tempBis\years.dta"  // panel ids: NR-years
rename ID NRid
sort NRid year
save "$tempBis\NRid-years.dta", replace

use "$dataNR\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta", clear  // prepare NR dataset
gen eleyear = year 
rename ID NRid
sort year
merge 1:1 NRid year using "$tempBis\NRid-years.dta"
drop _merge
gen candidate = 0 // years with candidacy
replace candidate = 1 if elected != .
label var candidate "Was a political candidate in latest election"
order NRid canton year name firstname birthyear sex job candidate elected
sort NRid year
save "$tempBis\RL_NRinfo_ypanel.dta", replace 

* prepare bisnode datasets (different NR links)
use "$tempBis\RL_Bisnode-NRids_1994-2018_tmp.dta", clear

preserve
keep if NRid_1 != ""
gen NRid = NRid_1
drop NRid_*
sort NRid year
save "$tempBis\RL_Bisnode-NRids_1994-2018_tmpNR1.dta", replace
restore

preserve
keep if NRid_2 != ""
gen NRid = NRid_2
drop NRid_*
sort NRid year
save "$tempBis\RL_Bisnode-NRids_1994-2018_tmpNR2.dta", replace
restore

preserve
keep if NRid_3 != ""
gen NRid = NRid_3
drop NRid_*
sort NRid year
save "$tempBis\RL_Bisnode-NRids_1994-2018_tmpNR3.dta", replace
restore

preserve
keep if NRid_4 != ""
gen NRid = NRid_4
drop NRid_*
sort NRid year
save "$tempBis\RL_Bisnode-NRids_1994-2018_tmpNR4.dta", replace
restore

* merge Bisnode & NR data
clear
forv i = 1(1)4 {
	use "$tempBis\RL_NRinfo_ypanel.dta", clear
	merge 1:m NRid year using "$tempBis\RL_Bisnode-NRids_1994-2018_tmpNR`i'.dta"
	keep if _merge == 3
	drop _merge
	save "$tempBis\RL_G7_NR-Bisnode_tmp`i'.dta", replace
}
use "$tempBis\RL_G7_NR-Bisnode_tmp1.dta", clear
forv i = 2(1)4 {
	append using "$tempBis\RL_G7_NR-Bisnode_tmp`i'.dta"
}	
duplicates report NRid VRid year duns rechtsform gremium funktion
sort NRid year VRid duns rechtsform gremium funktion
save "$tempBis\RL_G7_NR-Bisnode_tmp.dta", replace

* merge with universe of NRs
use "$tempBis\RL_NRinfo_ypanel.dta", clear
merge 1:m NRid year using "$tempBis\RL_G7_NR-Bisnode_tmp.dta"
drop _merge


* keep only observations for the period 1994-2017 (the 2018 contains data only until end of june and is thus incomplete)
keep if year>=1994 & year<2018 
/*
bysort NRid: egen NRcand = total(votes)
keep if NRcand > 0  // keep those observations that were involved in an election in this period (1980-2018), drop candidates other candidates (earlier)
*/

label var mandid "VR mandate ID"
label var NRlinks "No of NRid per VRid"
label var entrydate "VR entry date"
label var exitdate "VR exit date"
label var entryear "VR entry year"
label var exityear "VR exit year"
label var eleyear "Election year"

duplicates report NRid VRid year duns funktion // are there duplicate entries for the same person (NRid & VRid) in the same company (duns) and same function per year?  (no!)

sort NRid year VRid 
save "$dataNR\11_Directors_1994_2018\RL_G7_NR-Bisnode_1994-2017.dta", replace


*** Step 5: collapse data to create dataset for analysis
*use "$dataNR\11_Directors_1994_2018\RL_G7_NR-Bisnode_1994-2017.dta", clear

** Variable creation and naming
rename NRid ID
rename gdename municipality
rename gdenr_2018_w municipality_num
rename language_w language
rename VRid PID
rename duns CID

** locals
local first canton name firstname birthyear sex job elected ///
votemargin_rel incumbent tenure list listname_bfs cand_before1931 ///
municipality municipality_num language 

keep `first' ID year PID CID funktion kapitalnominalag // potentially: branch and real-estate company indicator

** collapse
gen all = 1 if CID != .

bysort year: egen cap75 = pctile(kapitalnominalag), p(75)
gen lrg = 1 if CID != . & kapitalnominalag > cap75  & kapitalnominalag < . 
gen sml = 1 if CID != . & kapitalnominalag <= cap75  & kapitalnominalag < .
gen prs = 1 if CID != . & funktion == "Präsident/in"

collapse (first) `first' (sum) all lrg sml prs, by(ID year PID)
collapse (first) `first' (sum) n_all_sum=all n_lrg_sum=lrg ///
	n_sml_sum=sml n_prs_sum=prs (mean) n_all_avg=all n_lrg_avg=lrg ///
	n_sml_avg=sml n_prs_avg=prs, by(ID year)  
	
foreach var in all lrg sml prs {
	gen i_`var' = .
	replace i_`var' = 0 if n_`var'_sum == 0
	replace i_`var' = 1 if n_`var'_sum > 0 & n_`var'_sum < .
}
label var ID "Candidate ID"
label var canton "Canton"
label var cand_before1931 "Candidate participated in 1925 and/or 1928 election"
label var tenure "Number of days in office at election date"
label var incumbent "Incumbent dummy at election date"
label var elected "Elected"
label var votemargin_rel "Running variable (relative vote margin)"
label var firstname "First name of candidate"
label var name "Surname of candidate"
label var list "Party list"
label var listname_bfs "Party names"
label var birthyear "Year of birth"
label var sex "Male=1, Female=0"
label var job "Job title"
label var municipality "Municipality name 2018 (residence)"
label var municipality_num "Municipality number 2018 (residence)"
label var language "Majority language in residence municipality"
label var i_all "At least one mandate, all companies"
label var i_lrg "At least one mandate, large companies"
label var i_sml "At least one mandate, small companies"
label var i_prs "At least one mandate as president"
label var n_all_sum "Number of mandates, all companies"
label var n_all_avg "Avg. number of mandates, all companies"
label var n_lrg_sum "Number of mandates, large companies"
label var n_sml_sum "Number of mandates, small companies"
label var n_lrg_avg "Avg. number of mandates, large companies"
label var n_sml_avg "Avg. number of mandates, small companies"
label var n_prs_sum "Number of mandates as president"
label var n_prs_avg "Avg. number of mandates as president"

foreach var of varlist i_* n_* {
	rename `var' `var'_s2
}

duplicates report ID year
sort ID year
save "$dataNR\11_Directors_1994_2018\Politicians_Directorships_1994-2017.dta", replace

erase "$tempBis\bisnode_fp_firstcheck_Generation7_0_corrected.dta"
erase "$tempBis\RL_G7_Bisnode_Person-Firmen_Geo.dta"
erase "$tempBis\years.dta"
erase "$tempBis\NRids.dta"
erase "$tempBis\NRid-years.dta"
forv i = 1(1)4 {
	erase "$tempBis\RL_Bisnode-NRids_1994-2018_tmpNR`i'.dta"
	erase "$tempBis\RL_G7_NR-Bisnode_tmp`i'.dta"
}
erase "$tempBis\RL_G7_NR-Bisnode_tmp.dta"
erase "$tempBis\RL_Bisnode-NRids_1994-2018_tmp.dta"
erase "$tempBis\RL_NRinfo_ypanel.dta"
erase "$tempBis\RL_bisnodeIDs_NRids.dta"
rmdir "$tempBis"

