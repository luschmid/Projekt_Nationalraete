
********************************
* (A) SET OVERALL PARAMETERS
********************************

/*** DEFINE MAIN PATH ***/
global hauptpfad "C:\Schmidlu\Dropbox\Projekt Nationalräte"


* PATH OF FUZZY MERGE FILES
local fuzzyPath .\fuzzy_merge

//

/*** SETTINGS: SET LOCALS ***/
* the cantons name
* global canton BE

cd "$hauptpfad"
use "$hauptpfad\02_Processed_Data\15_Elections_Honduras\elections_hn_fuzzy_merge_inptut.dta", clear

drop Party_abbr

save "$hauptpfad\02_Processed_Data\15_Elections_Honduras\elections_hn_mod.dta", replace


levelsof Department, local(departnmentlevels) 
display `departnmentlevels'
foreach x of local departnmentlevels {
*foreach x in ATLANTIDA{
global department `x'
display "$department"




**************************************
* (B) DO FUZZY MERGE
*****************************************


*if "$department"=="ATLANTIDA"{
* cantonal years
global years 2009 2013 2017
* the matching variables
global vars Candidate 
* the weights for the matching variables (for matches)
global yesweights 15
* the weights for the matching variables (for non-matches)
global noweights 5
* define thresholds for merging at the end (colour cells to simplify task of recoders)
global threshold_opt 0.9
global threshold_min 0.6
global threshold_max 0.9
*}

do "$hauptpfad\03_Code\15_Elections_Honduras\fuzzymergePROCEDURE"

}


**************************************
* (C) APPEND ALL FUZZY MERGE FILES
*****************************************

local thresholds 0.6 0.9


foreach x of local departnmentlevels {
local t_help=6
foreach t of local thresholds {
if "`x'"=="ATLANTIDA"{
use "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`t'_`x'.dta",clear
keep id_Stata candidateid_ score_
rename id_Stata id_Stata_`t_help'
*rename score_ score_`t_help'
save "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`t'_all.dta",replace
}
else{
use "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`t'_all.dta",clear
append using "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`t'_`x'.dta"
keep id_Stata* candidateid_ score_
replace id_Stata_`t_help'=id_Stata if id_Stata_`t_help'==""
*replace score_`t_help'=score_ if score_`t_help'==.

drop id_Stata
save "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`t'_all.dta",replace
}
erase "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`t'_`x'.dta"
local t_help=`t_help'+3
}
erase "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/`x'_all.dta"
erase "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_out_`x'.dta"
}

******************************************************
* (D) COMBINE FUZZY MERGE FILES WITH INITIAL DATASET
******************************************************

use "$hauptpfad\02_Processed_Data\15_Elections_Honduras\elections_hn_fuzzy_merge_input.dta", clear

sort Department Year Party Candidate Votes
gen candidateid_=_n

local thresholds 0.6 0.9
local t_help=6
foreach t of local thresholds {
merge 1:1 candidateid_ using "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`t'_all.dta", gen(merge_`t_help')
local t_help=`t_help'+3
}

drop merge*
order id_Stata* Candidate Year Party
sort Department Candidate Year Party

* generate a single id_Stata b/c id_Stata_6 and id_Stata_9 are exactly equal

gen id_Stata=id_Stata_9
drop id_Stata_6 id_Stata_9
order id_Stata* Candidate Year Party


save "$hauptpfad\02_Processed_Data\15_Elections_Honduras\elections_hn_fuzzy_merge_ouptut.dta", replace

erase "$hauptpfad\02_Processed_Data\15_Elections_Honduras\elections_hn_mod.dta"
display "CALCULATIONS COMPLETE; QUITTING"

*exit, STATA clear
