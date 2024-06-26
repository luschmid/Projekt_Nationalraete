clear
timer clear
timer on 1
set more off


disp "WELCOME TO FUZZY MERGE :)"
disp "MAIN PATH IS: $hauptpfad"
disp "DEPARTMENT IS: $department"
disp "LIST OF ALL YEARS: $years"
disp "LIST OF MATCHING VARIABLES: $vars"
disp "LIST OF YES-WEIGHTS: $yesweights"
disp "LIST OF NO-WEIGHTS: $noweights"

/**********************
	(A) SETUP
**********************/
* The canton
local department $department
* list of all years
local years $years
* the matching variables
local vars $vars
* the weights for the matching variables (for matches)
local yesweights $yesweights
* the weights for the matching variables (for non-matches)
local noweights $noweights
* define thresholds for merging at the end (colour cells to simplify task of recoders)
local threshold_opt $threshold_opt
local threshold_min $threshold_min
local threshold_max $threshold_max
* list of all desired thresholds the a fuzzy merge to be considered a match
local thresholds 0.6 0.9
* the variables needed for the reshaping (basically, the variables for which we want to keep the information for all years)
local vars1 candidateid* score* Year* Department* Party* Candidate* ID_TSE* Votes* Elected* 
local vars2 candidateid_ score_ Year_ Department_ Party_ Candidate_ ID_TSE_ Votes_ Elected_ 
* Drop unnecessary globals
macro drop canton years vars yesweights noweights threshold_min threshold_max
** STILL TO SOLVE / IMPROVE HERE: BEZIRK-THING!
** ALSO, DO NOT MERGE WITH MANUAL FILE FROM ERFASSER IF NOT CANTON BE; SEE "/*** MERGE ID_RECODED INFORMATION TO THE DATASETS ***/"





/******************************
	(B) DATASET PREPARATION
******************************/

use "$hauptpfad\02_Processed_Data\15_Elections_Honduras\elections_hn_mod.dta"
sort Department Year Party Candidate Votes
gen candidateid=_n
keep if Department=="`department'"
save "$hauptpfad\02_Processed_Data\15_Elections_Honduras\1_Reclink/`department'_all", replace


foreach year in `years'{
	preserve
	keep if Year == `year'
	save "$hauptpfad\02_Processed_Data\15_Elections_Honduras\1_Reclink/`department'_`year'.dta", replace
	restore
}



/*******************************************************
	(C) MAIN PROGRAMM: DO FUZZY MERGES AND SAVE DATASETS
********************************************************/




foreach threshold in `thresholds'{


	* (i) Use master dataset including all observations (department_all) and subsequently do fuzzy merge for each year seperately

	foreach year in `years'{
		use "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/`department'_`year'.dta"
		//destring occupations_first, force replace
		replace Year = Year+1
		gen idmatch = _n
		rename * *_`year'
		save "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/`department'new_`year'.dta", replace
	}


	foreach year in `years'{
		clear
		use "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/`department'_all.dta"
		gen idmatch = _n

		local vars_u
		foreach var in `vars'{
			local vars_u `vars_u' `var'_`year'
		}

		replace idmatch =_n
		disp "`vars'"
		disp "`year'"
		disp "`vars_u'"

		reclink `vars' using "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/`department'new_`year'.dta", idm(idmatch) idu(idmatch_`year') gen(score_`year') wmatch(`yesweights') wnomatch(`noweights') _merge(_fuzzymerge_`year') uvarlist(`vars_u') minscore(`threshold')
		capture rename UCandidate Candidate_`year'
		capture rename UParty Party_`year'
		save "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/inter_`year'", replace
	}	
	* (ii) Reshape temporary dataset (inter_`year') into wide format


	foreach year in `years'{
		clear
		use "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/inter_`year'.dta", replace
		*gen districtname_recoded_`year'=districtname_recoded
		keep `vars1'
		bysort candidateid: g id_indi=_n	
		reshape wide `vars2' , i( candidateid) j(id_indi)
		save "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/inter_`year'_wide.dta", replace
	}

	local vars_all
	foreach year in `years'{
		foreach var in `vars2'{
			local vars_all `vars_all' `var'`year'
		}
	}
	
	* (iii) Merge all fuzzy merge and bring it to long format
		    // note:  1. reshape long format: multiple observations (saved in wide format)
			// 		  2. reshape long format: years


	use "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/`department'_all.dta", clear
	foreach year in `years'{
		merge 1:1 candidateid using "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/inter_`year'_wide.dta", nogen
	}

	reshape long `vars_all', i( candidateid) j(check) // 1. reshape

	egen TotalScore = rowmean(score_*)
	drop if TotalScore==.
	drop TotalScore


	reshape long `vars2' , i( candidateid check) j(yearNew) // 2. reshape
	drop if score_==.
	
	sort candidateid candidateid_ yearNew // 1. step: keep only blocks of observations
	duplicates drop candidateid_, force
										  
	bysort candidateid_  yearNew: egen score_max=max(score_) if year!=yearNew  // 2. step: keep only observations with the maximal score
	bysort candidateid_: gen indi_candid=_N
	keep if score_max==score_ | score_max==.
	

	

	
	* (iv) Generate id that makes sense for recoders (id_Stata: Department-Year-Number)
	
	bysort candidateid: egen yearNewMin = min(yearNew)
	
	sort yearNew candidateid_
	
	bysort yearNew: gen temp = _n
	
	replace temp = 0 if yearNew!=yearNewMin
	bysort candidateid: egen number = max(temp)
	drop temp*
	
	tostring number, replace
	gen temp1=substr(number,1,1) 
	gen temp2=substr(number,2,1)
	gen temp3=substr(number,3,1)
	gen temp4=substr(number,4,1)
	replace number="0"+number if temp2==""
	replace number="0"+number if temp3==""
	replace number="0"+number if temp4==""
	
	gen id_Stata=Department+"-"+string(yearNewMin)+"-"+number
	
	
	* (v) Sort and order data and save dataset (uzzy_`threshold'_`canton'.dta)

	
	
	sort Department Candidate Year
	order id_Stata Department Candidate Year
		
	save "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`threshold'_`department'.dta", replace
	display "fuzzy_`threshold' saved:  $S_TIME  $S_DATE"
	clear
}

* (vi) Generate indicator whether merging differs between threshold_min and threshold_max

use "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`threshold_opt'_`department'.dta", clear // start with optimal threshold level
gen candidateid_threshold_opt=candidateid
drop candidateid
merge 1:1 candidateid_ using "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`threshold_min'_`department'.dta", gen(merge_thresholdsmin) // merge minimal threshold level
gen candidateid_threshold_min=candidateid
drop candidateid
merge 1:1 candidateid_ using "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_`threshold_max'_`department'.dta", gen(merge_thresholdsmax) // merge maximal threshold level
gen candidateid_threshold_max=candidateid
drop candidateid

* (vi-a) check whether ids differ between threshold_opt and threshold_max

gen indicator=0
replace indicator=1 if candidateid_threshold_opt!= candidateid_threshold_max //indi_min1!=indi_min2 | indi_max1!=indi_max2

	

	
* (vii) Generate outfile for typists

*export excel "$hauptpfad/02_Processed_Data/15_Elections_Honduras/using erf_`department'.xlsx", firstrow(variables) replace


save "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/fuzzy_out_`department'.dta", replace



* (vii) Erase temporary datafiles

/*** ERASE TEMPORARY DATAFILES ***/
foreach year in `years'{
	capture erase "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/inter_`year'.dta"
	capture erase "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/inter_`year'_wide.dta"
	capture erase "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/`department'new_`year'.dta"
	capture erase "$hauptpfad/02_Processed_Data/15_Elections_Honduras/1_Reclink/`department'_`year'.dta"
}


timer off 1
timer list 1

*exit, STATA clear
