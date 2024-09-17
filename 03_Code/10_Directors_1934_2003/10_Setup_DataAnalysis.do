clear
cap log close
set more 1
version 17

global dataNR "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\"
global dataRL "E:\12. Cloud\Dropbox\Record Linkage\02_Data\"
global tempSug "C:\Users\schelkem\SugarcubeSetup_tmp\"  

*global dataNR "C:\Current\Dropbox\Projekt Nationalräte\02_Processed_data\"
*global dataRL "C:\Current\Dropbox\Record Linkage\02_Data\"
*global tempSug "C:\SugarcubeSetup_tmp\"  

mkdir "C:\Users\schelkem\SugarcubeSetup_tmp", public  // reduce traffic on dropbox: create directory for all the tmp files.
*mkdir "C:\SugarcubeSetup_tmp", public  // reduce traffic on dropbox: create directory for all the tmp files.


**** Setup Sugarcube data for analyses ****
* Step 1: read-in Generation 3 after manual post-processing (FP & FN)
* Step 2: merge Sugarcube data on Sugarcube-IDs from RL
* Step 3: merge NR data on NR-IDs from RL
* Step 4: clean capital & functions
* Step 5: collapse data to create dataset for analysis


*** Step 1: Read-in Generation 3 after manual post-processing 
use "$dataNR\12_Record_linkage\02_Sugarcube\RL_Output_after_FP_and_FN.dta", clear
rename id_0 id_pol
rename e_cntr_w_0 e_id_polit
rename n_cntr_w_0 n_id_polit
rename id_1 id_sug
destring id_sug, replace
rename n_cntr_w_1 n_id_sug
rename e_cntr_w_1 e_id_sug
rename year_1 year_sug  // year_sug, n_id_sug, e_id_sug are pure ID variable! Not to be used as actual geo-location or year indication. Use the proper variables in the Sugarcube dataset

duplicates report id_pol n_id_polit e_id_polit id_sug n_id_sug e_id_sug year_sug  // no duplicates

duplicates tag id_pol id_sug n_id_sug e_id_sug year_sug, gen(duplNR)  // look for geo-duplicates in NRid
bysort id_pol id_sug n_id_sug e_id_sug year_sug: egen maxlink = max(link_score)  // keep only link with highest link score 
drop if duplNR>0 & maxlink!=link_score
duplicates report id_pol id_sug n_id_sug e_id_sug year_sug
drop n_id_polit e_id_polit maxlink duplNR

duplicates report id_sug n_id_sug e_id_sug year_sug // create dataset unique in Sugarcube ID
sort id_sug n_id_sug e_id_sug year_sug
bysort id_sug n_id_sug e_id_sug year_sug: gen n = _n
reshape wide id_pol link_score, i(id_sug n_id_sug e_id_sug year_sug) j(n)
gen NRlinks = 1  // some Sugarcube IDs have up to 3 NR links
replace NRlinks = 2 if id_pol2 != "" & id_pol3 == ""
replace NRlinks = 3 if id_pol2 != "" & id_pol3 != ""
sort year_sug id_sug n_id_sug e_id_sug 
save "$tempSug\RL_sugarcubeIDs_NRids.dta", replace

*** Step 2: merge Sugarcube data on Sugarcube IDs from RL
use "$dataNR\10_Directors_1934_2003\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", clear
sort year_sug id_sug n_id_sug e_id_sug
merge m:1 year_sug id_sug n_id_sug e_id_sug using "$tempSug\RL_sugarcubeIDs_NRids.dta"

keep if _merge == 3
drop _merge
rename titlejob titlejob_bus
rename firstname firstname_bus
rename lastname name_bus
rename male male_bus
rename gdename gdename_bus
rename gdenr_2018 gdenr_2018_bus
rename ctn ctn_bus
sort year id_pol1  // note (!): year != year_sug  -- > year_sug is purely ID variable in RL. year_sug is before year correction for early years!!!
save "$tempSug\RL_PP_Sugarcube_Person-Firmen_NRids.dta", replace

*** Step 3: merge NR data on NR-IDs from RL
* prepare NR dataset
clear
set obs 105 // create year dataset 1915 - 2020
gen n = _n-1
gen year = 1915 + n
drop n
save "$tempSug\years.dta", replace

use "$dataNR\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta", clear  // create NRid dataset
keep ID
duplicates drop ID, force
save "$tempSug\NRids.dta", replace
 
cross using "$tempSug\years.dta"  // panel ids: NR-years
rename ID NRid
sort NRid year
save "$tempSug\NRid-years.dta", replace

use "$dataNR\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta", clear  // prepare NR dataset
gen eleyear = year 
rename ID NRid
sort year
merge 1:1 NRid year using "$tempSug\NRid-years.dta"
drop _merge
gen id_pol = NRid
order NRid canton year name firstname birthyear sex job elected
sort year id_pol
save "$tempSug\RL_NRinfo_ypanel.dta", replace 

* prepare Sugarcube datasets (different NR links)
use "$tempSug\RL_PP_Sugarcube_Person-Firmen_NRids.dta", clear

preserve
keep if id_pol1 != ""
gen id_pol = id_pol1
drop id_pol1 id_pol2 id_pol3
gen linkscore = link_score1
drop link_score1 link_score2 link_score3
sort year id_pol
save "$tempSug\RL_Sugarcube-NRids_1934-2003_tmpNR1.dta", replace
restore

preserve
keep if id_pol2 != ""
gen id_pol = id_pol2
drop id_pol1 id_pol2 id_pol3 
gen linkscore = link_score2
drop link_score1 link_score2 link_score3
sort year id_pol
save "$tempSug\RL_Sugarcube-NRids_1934-2003_tmpNR2.dta", replace
restore

preserve
keep if id_pol3 != ""
gen id_pol = id_pol3
drop id_pol1 id_pol2 id_pol3
gen linkscore = link_score3
drop link_score1 link_score2 link_score3 
sort year id_pol
save "$tempSug\RL_Sugarcube-NRids_1934-2003_tmpNR3.dta", replace
restore

* merge Sugarcube & NR data
clear
forv i = 1(1)3 {
	use "$tempSug\RL_NRinfo_ypanel.dta", clear
	merge 1:m year id_pol using "$tempSug\RL_Sugarcube-NRids_1934-2003_tmpNR`i'.dta"
	keep if _merge == 3
	drop _merge
	save "$tempSug\RL_NR_Sugarcube_1934-2003_tmpNR`i'.dta", replace
}
use "$tempSug\RL_NR_Sugarcube_1934-2003_tmpNR1.dta", clear
forv i = 2(1)3 {
	append using "$tempSug\RL_NR_Sugarcube_1934-2003_tmpNR`i'.dta"
}	
duplicates report id_pol id_sug n_id_sug e_id_sug year CID
sort year id_pol id_sug n_id_sug e_id_sug CID 
save "$tempSug\RL_NR_Sugarcube_1934-2003_tmpNR.dta", replace

* merge with universe of NRs
use "$tempSug\RL_NRinfo_ypanel.dta", clear
merge 1:m year id_pol using "$tempSug\RL_NR_Sugarcube_1934-2003_tmpNR.dta"
drop _merge

* keep only observations for the period 1931-2003 
keep if year>=1931 & year<=2003 

* drop unnecessary or redundant variables
drop ID_dupl id_sug year_sug e_id_sug n_id_sug


* Deduplicate due to gde_options in Sugarcube-Persons (e.g., Rüti...)
duplicates report NRid year e_id_polit n_id_polit PID gde_options CID
duplicates tag NRid year e_id_polit n_id_polit PID CID, gen(geodupl)
tab geodupl
duplicates drop NRid year e_id_polit n_id_polit PID CID, force  // drop duplicate information, erase geo-information of remaining match (we do not know which match is correct)
replace gdename_bus = city if geodupl>0
replace gdenr_2018_bus = . if geodupl>0
replace GdeNr_E_CNTR = . if geodupl>0
replace GdeNr_N_CNTR = . if geodupl>0 
replace ctn_bus = "" if geodupl>0
drop geodupl gde_options


*** Step 4: clean capital & functions

** clean capital variable
/*
*ssc install egenmore
recast str capital // the egenmore program cannot handle StrL, therefore recast the string
egen capital_strings = sieve(capital), omit(.0123456789)
preserve
gen n=1
collapse (sum) n, by(capital_strings)
save "$tempSug\capital_strings.dta", replace
restore
*/
rename capital capital_orig
gen capital = capital_orig
replace capital = ustrlower(capital)
replace capital = ustrregexra(capital,"[a-z]","")
replace capital = usubinstr(capital, ")", "",.)
replace capital = usubinstr(capital, "í", "",.)
replace capital = usubinstr(capital, "'", "",.)
replace capital = usubinstr(capital, "-", "",.)
replace capital = usubinstr(capital, ",", "",.)
replace capital = usubinstr(capital, "!", "",.)
replace capital = usubinstr(capital, "`", "",.)
replace capital = usubinstr(capital, "�", "",.)
replace capital = usubinstr(capital, "/", "",.)
replace capital = usubinstr(capital, "?", "",.)
replace capital = usubinstr(capital, ":", "",.)
replace capital = usubinstr(capital, "è", "",.)
replace capital = usubinstr(capital, "ó", "",.)
replace capital = usubinstr(capital, "ö", "",.)
replace capital = usubinstr(capital, "ü", "",.)
replace capital = usubinstr(capital, "ä", "",.)
replace capital = usubinstr(capital, "   ", " ",.)
replace capital = usubinstr(capital, "  ", " ",.)
replace capital = usubinstr(capital, " ", "",.)
replace capital = ustrtrim(capital)
destring capital, replace force 
replace capital = 0.5 if capital_orig == "0.S"
replace capital = 0.5 if capital_orig == "0.5)."
replace capital = 0.05 if capital_orig == "0:05"
replace capital = 0.3 if capital_orig == "0:3"
replace capital = 0.3 if capital_orig == "0.3) 5"
replace capital = 0.1 if capital_orig == "0-1"
replace capital = 0.1 if capital_orig == "0.'"
replace capital = 0.05 if capital_orig == "0.05l"
replace capital = 32.66 if capital_orig == "32.6602).5"
replace capital = 0.1 if capital_orig == "/0,1"
replace capital = 145.0 if capital_orig == "145:0"
replace capital = 0.1 if capital_orig == "0.1) 5"
replace capital = 0.98 if capital_orig == "O,98"
replace capital = 0.06 if capital_orig == "0.C6"
replace capital = . if capital_orig == "Liesta140,1"
replace capital = . if capital_orig == ",4,0"
replace capital = . if capital_orig == "41,"
replace capital = 0.1 if capital_orig == "0.i"
replace capital = 0.05 if capital_orig == "0..05"
replace capital = 0.05 if capital_orig == "0.0.5"
replace capital = 0.05 if capital_orig == "0.05) Pr.-"
replace capital = 0.1 if capital_orig == "0.1) Vr."
replace capital = 0.1 if capital_orig == "0.1."
replace capital = 0.15 if capital_orig == "0.15."
replace capital = 0.1 if capital_orig == "0.25)' t..."
replace capital = 140 if capital_orig == "FMV). Sion 140.0"
replace capital = 0.05 if capital_orig == "OISA). Vevey 10.051"

sum capital, d

sort NRid year PID CID
save "$dataNR\10_Directors_1934_2003\RL_NR-Sugarcube_1934-2003.dta", replace


*** Step 5: collapse data to create dataset for analysis
use "$dataNR\10_Directors_1934_2003\RL_NR-Sugarcube_1934-2003.dta", clear

** Variable creation and naming
rename NRid ID
rename gdename municipality
rename gdenr_2018_w municipality_num
rename language_w language
gen dir_year = 0
replace dir_year = 1 if inlist(year, 1934, 1943, 1960, 1962, 1963, 1964, 1965, ///
	1966, 1969, 1972) | year == 1975 | inrange(year, 1979, 2003) 
	
** locals
local first canton name firstname birthyear sex job elected ///
votemargin_rel incumbent tenure list listname_bfs cand_before1931 ///
municipality municipality_num language dir_year ///
EDateJoining1 EDateLeaving1 EDateJoining2 EDateLeaving2 EDateJoining3 ///
EDateLeaving3 EDateJoining4 EDateLeaving4 EDateJoining5 EDateLeaving5 ///
EDateJoining6 EDateLeaving6

keep `first' ID year PID CID capital compFunct 

** collapse
gen all = 1 if CID != .
*bysort year: egen cap75 = pctile(capital), p(75)
*gen lrg = 1 if CID != . & capital > cap75  & capital < . 
*gen sml = 1 if CID != . & capital <= cap75  & capital < .
bysort year: egen cap90 = pctile(capital), p(90)

// read out and save capital thresholds
preserve
keep cap90 year
collapse (first) cap90, by(year)
save "$dataNR\10_Directors_1934_2003\sugar_cap90.dta", replace
restore

gen lrg = 1 if CID != . & capital > cap90  & capital < . 
gen sml = 1 if CID != . & capital <= cap90  & capital < .
gen prs = 1 if CID != . & compFunct == "president" | compFunct == "pre" | compFunct == "president/administrator"

collapse (first) `first' (sum) all lrg sml prs, by(ID year PID) 
collapse (first) `first' (sum) n_all_sum=all n_lrg_sum=lrg ///
	n_sml_sum=sml n_prs_sum=prs (mean) n_all_avg=all n_lrg_avg=lrg ///
	n_sml_avg=sml n_prs_avg=prs, by(ID year)  

foreach var of varlist n_* {  
	replace `var' = . if dir_year == 0
}
replace n_prs_sum = . if year < 1975  
replace n_prs_avg = . if year < 1975
// Information on function would be available in 1934. 
// We ignore this year to avoid composition effects for individual lags and leads.


foreach var in all lrg sml prs { 
	gen i_`var' = .
	replace i_`var' = 0 if dir_year == 1 & n_`var'_sum == 0
	replace i_`var' = 1 if dir_year == 1 & n_`var'_sum > 0 & n_`var'_sum < .
}


*** Step 6: Definitions of who is in office and when

gen inoffice=0
forvalues i=1(1)6 { 
gen EYearJoining`i' = year(EDateJoining`i')
gen EYearLeaving`i' = year(EDateLeaving`i')
bysort ID: egen EYearJoining_temp=min(EYearJoining`i')
bysort ID: egen EYearLeaving_temp=min(EYearLeaving`i')
replace inoffice=1 if year>=EYearJoining_temp & year<= EYearLeaving_temp
drop EYearJoining_temp EYearLeaving_temp
}
drop EDate* EYear*


*** Step 7: Import information on precision & recall at the link and mandate level

merge 1:1 ID year using "$dataNR\12_Record_linkage\02_Sugarcube\pr_link.dta"
drop _merge
merge 1:1 ID year using "$dataNR\12_Record_linkage\02_Sugarcube\pr_mandate.dta"
drop _merge


*** Labelling

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
label var dir_year "Year available in Directory of Directors"
label var municipality "Municipality name 2018 (residence)"
label var municipality_num "Municipality number 2018 (residence)"
label var language "Majority language in residence municipality"
label var inoffice "Person is currently in National Council"
label var i_all "At least one mandate, all companies"
label var i_lrg "At least one mandate, large companies"
label var i_sml "At least one mandate, small companies"
label var i_prs "At least one mandate as president"
label var n_all_sum "Number of mandates, all companies"
label var n_all_avg "Avg. number of mandates, all companies"
label var n_lrg_sum "Number of mandates, large companies"
label var n_lrg_avg "Avg. number of mandates, large companies"
label var n_sml_sum "Number of mandates, small companies"
label var n_sml_avg "Avg. number of mandates, small companies"
label var n_prs_sum "Number of mandates as president"
label var n_prs_avg "Avg. number of mandates as president"
label var prc_lk "Precision at link-level"
label var rcl_lk "Recall at link-level"
label var f1_lk "F1 at link-level"
label var prc_mn "Precision at mandate-level"
label var rcl_mn "Recall at mandate-level"
label var f1_mn "F1 at mandate-level"

foreach var of varlist i_* n_* prc_* rcl_* f1_* {
	rename `var' `var'_s1
}

duplicates report ID year
sort ID year
save "$dataNR\10_Directors_1934_2003\Politicians_Directorships_1934-2003.dta", replace

