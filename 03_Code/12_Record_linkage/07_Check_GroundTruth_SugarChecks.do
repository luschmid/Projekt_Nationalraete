clear
clear matrix
clear mata
cap log close
set more 1


**** Recordlinkage: Sugarchecks ****

global data "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\12_Record_linkage\02_Sugarcube\"

/* First round coding
forv v = 1(1)19 {
	forv u = 1(1)2 {
		import excel using "${data}04_Round1_In_from_Coders\Sugar_Checks`v'_out_`u'.xlsx", sheet("Sugar_Checks`v'") first clear
		rename Match Match`u'
		rename Comment Comment`u'
		cap rename A Ini_sort
		drop if Ini_sort == . & ID_Pol == "" // there are missings at the end of the data in Sugar_Checks 4 & 8
		save "${data}04_Round1_In_from_Coders\Sugar_Checks`v'_out_`u'.dta", replace
	}
	preserve
	use "${data}04_Round1_In_from_Coders\Sugar_Checks`v'_out_2.dta", clear // keep only coding & comments
	keep Ini_sort ID_Pol ID_VR Jahr_VR Comment2 Match2
	merge 1:1 Ini_sort ID_Pol ID_VR Jahr_VR using "${data}04_Round1_In_from_Coders\Sugar_Checks`v'_out_1.dta"
	drop _merge
	replace Match1 = 1 if Match1 != .  // some have coded matches with CaseID
	replace Match2 = 1 if Match2 != .  // some have coded matches with CaseID
	gen match_check = 1 if Match1 == Match2
	replace match_check = 0 if Match1 != Match2
	gen match_final = .
	gen comment = ""
	order comment match_final match_check Match1 Match2 CaseID Vorname_Pol Nachname_Pol Jahr_VR Vorname_VR Nachname_VR Alter Distanz Wohngemeinde_Pol Wohngemeinde_VR ///
	Beruf_Pol Firmen_VR Geo_Code_Pol Geo_Code_VR Ini_sort ID_Pol ID_VR Comment1 Comment2 
	save "${data}04_Round1_In_from_Coders\Sugar_Checks`v'_comp.dta", replace
	rename match_final match_final1
	rename comment comment_1
	export excel using "${data}05_Round2_Out_to_Coders\Sugar_Checks`v'_comp1.xlsx", sheet("Sugar_Checks`v'") first(var) replace
	rename match_final1 match_final2
	rename comment_1 comment_2
	export excel using "${data}05_Round2_Out_to_Coders\Sugar_Checks`v'_comp2.xlsx", sheet("Sugar_Checks`v'") first(var) replace
	restore
	erase "${data}04_Round1_In_from_Coders\Sugar_Checks`v'_comp.dta"
	forv u = 1(1)2 {
		erase "${data}04_Round1_In_from_Coders\Sugar_Checks`v'_out_`u'.dta"
	}	
}

* Second round coding
forv v = 1(1)19 {
	forv u = 1(1)2 {
		import excel using "${data}06_Round2_In_from_Coders\Sugar_Checks`v'_comp`u'_out.xlsx", sheet("Sugar_Checks`v'") first clear
		drop if Ini_sort == . & ID_Pol == "" 
		save "${data}06_Round2_In_from_Coders\Sugar_Checks`v'_comp`u'_2nd.dta", replace
	}
	use "${data}06_Round2_In_from_Coders\Sugar_Checks`v'_comp1_2nd.dta", clear
	merge 1:1 Ini_sort ID_Pol ID_VR Jahr_VR using "${data}06_Round2_In_from_Coders\Sugar_Checks`v'_comp2_2nd.dta"
	gen match_check2 = 1 if match_final1 == match_final2
	replace match_check2 = 0 if match_final1 != match_final2
	replace match_final1 = . if match_final1 == 0 & comment_1 != ""  // coders where not supposed to indicate "0" (leave empty)
	replace match_final2 = . if match_final2 == 0 & comment_2 != ""
	gen finalmatch = 77777
	replace finalmatch = 1 if Match1 == 1 & Match2 == 1  // both first-round coders agree
	replace finalmatch = 0 if Match1 == . & Match2 == .  // both first-round coders agree
	replace finalmatch = 1 if match_final1 == 1 & match_final2 == 1  // both second-round coders agree on match (typically in cases where disagreement in first round; in few cases: BOTH agree against first)
	replace finalmatch = 0 if match_final1 == . & match_final2 == . & match_check == 0  // both second-round coders agree in case the first coders disagreed
	replace finalmatch = 99999 if match_final1 == . & match_final2 == . & match_check == 1 & (comment_1 != "" | comment_2 != "") // BOTH second-round coders DISagree with BOTH first round coders: 99999 = check!
	replace finalmatch = 99999 if match_final1 == 1 & match_final2 == . & match_check == 1 & (comment_1 != "" | comment_2 != "") // one second-round coder DISagrees with BOTH first round coders: 99999 = check!
	replace finalmatch = 99999 if match_final1 == . & match_final2 == 1 & match_check == 1 & (comment_1 != "" | comment_2 != "") // one second-round coder DISagrees with BOTH first round coders: 99999 = check!
	replace finalmatch = 66666 if match_final1 == 1 & match_final2 == . & match_check == 1 & (comment_1 == "" | comment_2 != "") // one second-round coder DISagrees with BOTH first round coders, no comment: 66666 = check!
	replace finalmatch = 66666 if match_final1 == . & match_final2 == 1 & match_check == 1 & (comment_1 == "" | comment_2 != "") // one second-round coder DISagrees with BOTH first round coders, no comment: 66666 = check!
	// finalmatch = 77777 --> to be coded in last (3rd) round (disagreement of second-round coders and first-round coders) 
	gen comment_ms = ""
	drop _merge
	order comment_ms comment_1 comment_2 finalmatch match_check2 match_final1 match_final2 Match1 Match2 match_check
	export excel using "${data}07_Round3_Out_to_Coders\Sugar_Checks`v'_compfin.xlsx", sheet("Sugar_Checks`v'") first(var) replace
	
	forv u = 1(1)2 {
		erase "${data}06_Round2_In_from_Coders\Sugar_Checks`v'_comp`u'_2nd.dta"
	}
}

*/

* Read-in final codings
forv u = 1(1)19 {
		import excel using "${data}09_Final_Coding\Sugar_Checks`u'_final.xlsx", sheet("Sugar_Checks`u'") first clear
		gen fileID = `u'
		save "${data}09_Final_Coding\Sugar_Checks`u'_final.dta", replace
	}
use "${data}09_Final_Coding\Sugar_Checks1_final.dta", clear
forv u = 2(1)19 {
 append using "${data}09_Final_Coding\Sugar_Checks`u'_final.dta", force
}
drop if Nachname_Pol == "" // there are lots of missings from excel-import
replace Distanz = "" if Distanz=="NA"
destring Distanz, replace

// compare final matches with perfect matches
gen perfectmatch = 0 
replace perfectmatch = 1 if Vorname_Pol==Vorname_VR & Nachname_Pol==Nachname_VR & Distanz==0 & Alter<90 & Alter>17
tab perfectmatch finalmatch
replace finalmatch = 1 if perfectmatch == 1 & finalmatch == 0  // this is a convention decided by Lukas & Mark (01.07.2021): it affects 94 cases in which the perfect match is unlikely to be correct because the VR mandates suggest a person to be too old in later years (e.g., beyond 100 years) or the profession is very unlikely (or other). The convention was chosen because the record linkage algorithm has no way to deal with this ambiguity (perfect match on attributes but not coded as match).
drop perfectmatch

/*
** Final manual checks by Mark: corrections directly made in "Sugar_Checks`u'final.xlsx" files 
// more than one match per ID_VR - Jahr_VR pair
preserve
duplicates tag ID_Pol ID_VR Jahr_VR if finalmatch == 1, gen(duplIDs)
keep if duplIDs>0 & duplIDs<10
export excel using "${data}09_Final_Coding\Sugar_Checks_duplicate-matches.xlsx", first(var) replace 
restore

// check if ID_Pol - ID_VR - Jahr_VR combinations have all been found
bysort ID_Pol ID_VR Jahr_VR: egen matchsum = total(finalmatch)
gen help = 1
bysort ID_Pol ID_VR Jahr_VR: egen matchsumpot = total(help) // sum up all potential unique ID_Pol - ID_VR - Jahr_VR combinations (lots of combinations that remain irrelevant)
replace matchsumpot = . if matchsum<1 // we are not interested in all combinations, but only those who have been identified at least once as match.
drop help
tab matchsumpot matchsum // -- > if found once, checked and made consistent over various files
*/

** corrections January 2023 (mostly they become too old when linking over time through VR mandates)
replace finalmatch = 1 if Vorname_VR=="alfred" & Nachname_VR=="diezi" &  Geo_Code_VR == "2677500/1253800" & Jahr_VR == 1993
replace finalmatch = 0 if Vorname_VR=="heinrich" & Nachname_VR=="walder" &  Geo_Code_VR == "2683100/1247100"
replace finalmatch = 1 if Vorname_VR=="ernst" & Nachname_VR=="weber" &  Geo_Code_VR == "2683400/1241300"
replace finalmatch = 0 if Vorname_VR=="henri" & Nachname_VR=="pasche" &  Geo_Code_VR == "2538200/1152400"
replace finalmatch = 0 if Vorname_VR=="max" & Nachname_VR=="keller" &  Geo_Code_VR == "2697200/1261700"
replace finalmatch = 0 if Vorname_VR=="paul" & Nachname_VR=="brunner" &  Geo_Code_VR == "2683100/1247100"
replace finalmatch = 0 if Vorname_VR=="pierre" & Nachname_VR=="rochat" &  Geo_Code_VR == "2538200/1152400" & inrange(Jahr_VR, 1961, 1963)
replace finalmatch = 0 if Vorname_VR=="walter" & Nachname_VR=="bosshard" &  Geo_Code_VR == "2671100/1246800"
replace finalmatch = 0 if Vorname_VR=="walter" & Nachname_VR=="mueller" &  Geo_Code_VR == "2672500/1251000"
replace finalmatch = 0 if Vorname_VR=="werner" & Nachname_VR=="hadorn" &  Geo_Code_VR == "2618600/1170600"

save "${data}09_Final_Coding\Sugar_Checks_final.dta", replace

forv u = 1(1)19 {
 erase "${data}09_Final_Coding\Sugar_Checks`u'_final.dta"
}


* Prepare ground truth and export ground truth's ID's for Sandro
use "${data}09_Final_Coding\Sugar_Checks_final.dta", clear

// generate geo-ID variable (normalize for merge)
split Geo_Code_Pol, parse(/)
rename Geo_Code_Pol1 Geo_Code_E_Pol
rename Geo_Code_Pol2 Geo_Code_N_Pol
destring Geo_Code_E_Pol, replace
destring Geo_Code_N_Pol, replace
replace Geo_Code_E_Pol = int(Geo_Code_E_Pol) 
replace Geo_Code_N_Pol = int(Geo_Code_N_Pol)
drop Geo_Code_Pol
split Geo_Code_VR, parse(/)
rename Geo_Code_VR1 Geo_Code_E_VR
rename Geo_Code_VR2 Geo_Code_N_VR
destring Geo_Code_E_VR, replace
destring Geo_Code_N_VR, replace
replace Geo_Code_E_VR = int(Geo_Code_E_VR)
replace Geo_Code_N_VR = int(Geo_Code_N_VR)
drop Geo_Code_VR

// reduce Sugar_Check_final to only contain the relevant match-info for NR candiates  -- > ground truth
foreach var in ID_VR Jahr_VR Geo_Code_E_VR Geo_Code_N_VR {  // only keep information for actual matches
	replace `var' = . if finalmatch == 0
}
// eliminate all duplicate information per relevant NR candidate 
// (do not use duplicates drop, as there may be multiple identical (w.r.t. ID-year) but independent VR observations within the same geo-location (e.g., gottlieb luescher))
// -- > 1) if there are NO VR matches (help1 == 0): keep at least one NR candidate observation (help2 == 1)
// -- > 2) if there are VR matches (help1 > 0): keep all non-missing VR matches

bysort ID_Pol Geo_Code_E_Pol Geo_Code_N_Pol: egen help1 = total(finalmatch)  // number of matches per NR
bysort ID_Pol Geo_Code_E_Pol Geo_Code_N_Pol: gen help2 = _n

gen groundtruth1 = 1 if help1 == 0 & help2 == 1  // 1) no VR matches: keep first observation for all NR without a VR match
gen groundtruth2 = 1 if ID_VR != .  // 2) keep all observations with VR matches
gen groundtruth = 0
replace groundtruth = 1 if groundtruth1 == 1 | groundtruth2 == 1
drop if groundtruth == 0
drop help1 help2 groundtruth*
duplicates report ID_Pol Geo_Code_E_Pol Geo_Code_N_Pol ID_VR Jahr_VR Geo_Code_E_VR Geo_Code_N_VR
duplicates drop ID_Pol Geo_Code_E_Pol Geo_Code_N_Pol ID_VR Jahr_VR Geo_Code_E_VR Geo_Code_N_VR, force

// geo-corrections for match with newest sugarcube geo-referencing

replace Geo_Code_E_VR = 2716824 if Geo_Code_E_VR == 2716825  // 54 cases of slight miscoding of Lugano
replace Geo_Code_N_VR = 1095774 if Geo_Code_N_VR == 1095775

replace Geo_Code_E_VR = 2504100 if ID_VR == 11724 & Jahr_VR == 1942
replace Geo_Code_E_VR = 2540500 if ID_VR == 953 & Jahr_VR == 1933
replace Geo_Code_E_VR = 2564300 if ID_VR == 6602 & Jahr_VR == 1933
replace Geo_Code_E_VR = 2626400 if ID_VR == 12491 & Jahr_VR == 1961
replace Geo_Code_E_VR = 2626400 if ID_VR == 13158 & Jahr_VR == 1962
replace Geo_Code_E_VR = 2626400 if ID_VR == 10933 & Jahr_VR == 1959
replace Geo_Code_E_VR = 2567100 if ID_VR == 1426 & Jahr_VR == 1933
replace Geo_Code_E_VR = 2685800 if ID_VR == 8317 & Jahr_VR == 1933

replace Geo_Code_N_VR = 1132100 if ID_VR == 11724 & Jahr_VR == 1942
replace Geo_Code_N_VR = 1187000 if ID_VR == 953 & Jahr_VR == 1933
replace Geo_Code_N_VR = 1192000 if ID_VR == 6602 & Jahr_VR == 1933
replace Geo_Code_N_VR = 1198900 if ID_VR == 12491 & Jahr_VR == 1961
replace Geo_Code_N_VR = 1198900 if ID_VR == 13158 & Jahr_VR == 1962
replace Geo_Code_N_VR = 1198900 if ID_VR == 10933 & Jahr_VR == 1959
replace Geo_Code_N_VR = 1229100 if ID_VR == 1426 & Jahr_VR == 1933
replace Geo_Code_N_VR = 1243800 if ID_VR == 8317 & Jahr_VR == 1933

replace Geo_Code_E_VR = . if ID_VR == 135352 & Jahr_VR == 1987
replace Geo_Code_N_VR = . if ID_VR == 135352 & Jahr_VR == 1987

sort ID_Pol Geo_Code_E_Pol Geo_Code_N_Pol ID_VR Jahr_VR Geo_Code_E_VR Geo_Code_N_VR
save "${data}Sugarcube_Ground_Truth.dta", replace

// renaming of variables for "Sugarcube Ground Thruth ID-file" for Sandro
rename ID_Pol id_polit
rename Geo_Code_E_Pol e_id_polit
rename Geo_Code_N_Pol n_id_polit
rename ID_VR id_sug
rename Jahr_VR year_sug
rename Geo_Code_E_VR e_id_sug 
rename Geo_Code_N_VR n_id_sug

keep id_polit e_id_polit n_id_polit id_sug year_sug e_id_sug n_id_sug
order id_polit e_id_polit n_id_polit id_sug year_sug e_id_sug n_id_sug

save "${data}Sugarcube_Ground_Truth_ID_Sandro.dta", replace
export delimited using "${data}Sugarcube_Ground_Truth_ID_Sandro.csv", delimiter(",")  replace

