clear
timer clear
timer on 1
set more off


disp "WELCOME TO FUZZY MERGE :)"
disp "MAIN PATH IS: $hauptpfad"
disp "CANTON IS: $canton"
disp "LIST OF ALL YEARS: $years"
disp "LIST OF MATCHING VARIABLES: $vars"
disp "LIST OF YES-WEIGHTS: $yesweights"
disp "LIST OF NO-WEIGHTS: $noweights"

/**********************
	(A) SETUP
**********************/
* The canton
local canton $canton
* list of all years
local years $years
* the matching variables
local vars $vars
* the weights for the matching variables (for matches)
local yesweights $yesweights
* the weights for the matching variables (for non-matches)
local noweights $noweights
* list of all desired thresholds the a fuzzy merge to be considered a match
local thresholds 0.6 0.7 0.8 0.9
* define thresholds for merging at the end (colour cells to simplify task of recoders)
local threshold_opt $threshold_opt
local threshold_min $threshold_min
local threshold_max $threshold_max
* the variables needed for the reshaping (basically, the variables for which we want to keep the information for all years)
//local vars1 candidateid lastname_ firstname_ party_ gender_ birthyear_ LeftRight_ canton_
local vars1 candidateid_ lastname_ firstname_ listname_ list_nr_ score_ party_ town_ gender_ birthyear_ LeftRight_ canton_
local vars2 candidateid_*  lastname* firstname* listname* list_nr*  score* party* town*  gender* birthyear* LeftRight* canton*
// variables to keep (equal var1 and var2 but with candidateid)
local vars3  candidateid* score* lastname* firstname* listname* list_nr* party* town* gender* birthyear* LeftRight* canton*

* Drop unnecessary globals
//macro drop canton years vars yesweights noweights
** STILL TO SOLVE / IMPROVE HERE: BEZIRK-THING!
** ALSO, DO NOT MERGE WITH MANUAL FILE FROM ERFASSER IF NOT CANTON BE; SEE "/*** MERGE ID_RECODED INFORMATION TO THE DATASETS ***/"






/******************************
	(B) DATASET PREPARATION
******************************/

cd "$path_fuzzy_data"
use `canton'_all.dta



foreach year in `years'{
	preserve
	capture gen gender=geschlecht
	keep if year == `year'
	save `canton'_`year'.dta, replace
	restore
}



/*******************************************************
	(C) MAIN PROGRAMM: DO FUZZY MERGES AND SAVE DATASETS
********************************************************/

foreach threshold in `thresholds'{

//local threshold 0.6


	foreach year in `years'{
		use `canton'_`year'.dta
		replace year = year+1
		gen idmatch = _n
		rename * *_`year'
		save `canton'new_`year', replace
	}


	foreach year in `years'{
		clear
		* (i) Use master dataset including all observations (canton_all) and subsequently do fuzzy merge for each year seperately
		use `canton'_all.dta
		gen idmatch = _n

		local vars_u
		foreach var in `vars'{
			local vars_u `vars_u' `var'_`year'
		}

		replace idmatch =_n

		reclink `vars' using `canton'new_`year'.dta, idm(idmatch) idu(idmatch_`year') gen(score_`year') wmatch(`yesweights') wnomatch(`noweights') _merge(_fuzzymerge_`year') uvarlist(`vars_u') minscore(`threshold')

		capture rename Ulastname lastname_`year'
		capture rename Ufirstname firstname_`year'
		capture rename ULeftRight LeftRight_`year'
		capture rename Ubirthyear birthyear_`year'
		capture rename Ugender gender_`year'	
		capture rename Utown_nr town_nr_`year'	
		

		//rename canton canton_`year'

		//drop firstname_start* lastname_start*
		save inter_`year', replace
	}	
	
	* (ii) Reshape temporary dataset (inter_`year') into wide format
		
	foreach year in `years'{
		clear
		use inter_`year', replace
		//drop year_string*
		//gen bezirk_`year'=districtname_`year'
		keep `vars3'
		bysort candidateid: g id_indi=_n	
		reshape wide `vars2' , i( candidateid) j(id_indi)
		save inter_`year'_wide.dta, replace
	}
	

	local vars_all 
	foreach year in `years'{
		foreach var in `vars1'{
			local vars_all `vars_all' `var'`year'
		}
	}
	
	display "`vars_all'"
	
	* (iii) Merge all fuzzy merge and bring it to long format
		    // note:  1. reshape long format: multiple observations (saved in wide format)
			// 		  2. reshape long format: years


	use `canton'_all, clear
	foreach year in `years'{
		merge 1:1 candidateid using inter_`year'_wide, nogen
	}
	
	save temp.dta, replace
	
	reshape long `vars_all', i( candidateid) j(check)

	egen TotalScore = rowmean(score_*)
	drop if TotalScore==.
	drop TotalScore


	reshape long `vars1' , i( candidateid check) j(yearNew)
	drop if score_==.
	
	sort candidateid candidateid_ yearNew // 1. step: keep only blocks of observations
	duplicates drop candidateid_, force
										  
	bysort candidateid_  yearNew: egen score_max=max(score_) if year!=yearNew  // 2. step: keep only observations with the maximal score
	bysort candidateid_: gen indi_candid=_N
	keep if score_max==score_ | score_max==.
	
	order score_ candidateid candidateid_ yearNew firstname_ lastname_ birthyear   yearNew party_ 
	
	
	
	* (iv) Generate id that makes sense for recoders (id_Stata: Canton-Year-Number)

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
	
	gen id_Stata=canton+"-"+string(yearNewMin)+"-"+number
	
	gen firstname_start_ = strlower(substr(firstname_,1,3))
	gen lastname_start_ = strlower(substr(lastname_,1,3))



	save fuzzy_`threshold'_`canton'.dta, replace
	display "fuzzy_`threshold' saved:  $S_TIME  $S_DATE"
	clear
}



	* (vi) Generate indicator whether merging differs between threshold_min and threshold_max

	use fuzzy_`threshold_opt'_`canton'.dta, clear // start with optimal threshold level
	gen candidateid_threshold_opt=candidateid
	drop candidateid
	merge 1:1 candidateid_ using fuzzy_`threshold_min'_`canton'.dta, gen(merge_thresholdsmin) // merge minimal threshold level
	gen candidateid_threshold_min=candidateid
	drop candidateid
	merge 1:1 candidateid_ using fuzzy_`threshold_max'_`canton'.dta, gen(merge_thresholdsmax) // merge maximal threshold level
	gen candidateid_threshold_max=candidateid
	drop candidateid
	/*bysort candidateid_threshold_min: egen indi_min1=min(candidateid_threshold_opt)
	bysort candidateid_threshold_max: egen indi_min2=min(candidateid_threshold_opt)
	bysort candidateid_threshold_min: egen indi_max1=max(candidateid_threshold_opt)
	bysort candidateid_threshold_max: egen indi_max2=max(candidateid_threshold_opt)
	*/
	* (vi-a) check whether ids differ between threshold_min and threshold_max

	gen indicator=0
	replace indicator=1 if candidateid_threshold_opt!= candidateid_threshold_max // indi_min1!=indi_min2 | indi_max1!=indi_max2

	order indicator candidateid_threshold_opt candidateid_threshold_min candidateid_threshold_max firstname_ lastname_  yearNew  birthyear_  canton   town_    candidateid_ 
	sort canton firstname_start_ lastname_start_ town_ yearNew candidateid_threshold_opt candidateid_

	
	* (vi-b) check whether some persons are at very different positions in the dataset 

	
	preserve

	gen numbering=_n
	sort yearNew
	bysort candidateid_threshold_min: gen yearNew_ind=_n
	xtset candidateid_threshold_min yearNew_ind
	
	gen numbering_l=abs(D1.numbering)
	bysort candidateid_threshold_min: egen numbering_diff=max(numbering_l)
	
	// keep if  inrange(numbering_diff,8,100000) // save all those with a substantial distance between two equal observations 
	order indicator candidateid_threshold_min candidateid_threshold_min candidateid_threshold_max canton firstname_start_ lastname_start_ yearNew candidateid_
	sort canton firstname_start_ lastname_start_ town_ yearNew  candidateid_threshold_opt candidateid_ 
	
	save fuzzy_check_numbering_`canton'.dta, replace // USE THIS DATA TO MERGE LATER TO CHECKED LISTS
	restore
	
	
	* (vii) Generate outfile for typists
	
	sort firstname_start_ lastname_start_  candidateid_threshold_opt
	egen candidateid_threshold_min_sorted=group(candidateid_threshold_min)
	
	order indicator  candidateid_threshold_opt candidateid_threshold_min candidateid_threshold_max canton firstname_start_ lastname_start_ yearNew candidateid_
	sort canton  firstname_start_ lastname_start_  town_ yearNew candidateid_threshold_opt candidateid_

	// for erfasser in germany: 
	/*
	order indicator candidateid_threshold_min_sorted candidateid_threshold_opt candidateid_threshold_min candidateid_threshold_max canton firstname_start_ lastname_start_ yearNew candidateid_
	sort canton candidateid_threshold_min_sorted firstname_start_ lastname_start_  town_ yearNew candidateid_threshold_opt candidateid_
	*/
	preserve

	capture gen id_recoded = ""
	gen falsch = ""
	gen unsicher = ""
	gen Fehler = ""
	gen Kommentar = ""
	rename firstname_ Vorname
	rename lastname_ Nachname
	rename birthyear_ Geburtsjahr
	rename listname_ Listenbezeichnung
	rename gender Geschlecht
	//rename occupations_first_ Beruf
	//rename occupations_second_ Beruf2
	//rename title_ Titel
	rename town_ Gemeinde
	rename yearNew Wahljahr
	//rename bezirk_ Bezirk


	local keepErfasser candidateid_threshold_min_sorted indicator candidateid_threshold_opt candidateid_threshold_min candidateid_threshold_max candidateid_ Vorname Nachname id_Stata id_recoded unsicher Kommentar Geschlecht Geburtsjahr Listenbezeichnung Gemeinde Wahljahr
	order `keepErfasser'
	keep `keepErfasser'
	export excel using erf_`canton'.xlsx, firstrow(variables) replace

	restore

	save fuzzy_out_`canton'.dta, replace

//////////////////////////////////////////////////////////////////////////////////////////////////////
// INCLUDE HERE IF NEEDED: MERGE ID_RECODED INFORMATION TO THE DATASETS (DO-FILE: COMPARE MEN VS. MACHINE)***///
//////////////////////////////////////////////////////////////////////////////////////////////////////

	* (vii) Erase temporary datafiles


	/*** ERASE TEMPORARY DATAFILES ***/
	foreach year in `years'{
	capture erase inter_`year'.dta
	capture erase inter_`year'_wide.dta
	capture erase `canton'new_`year'.dta
	capture erase `canton'_`year'.dta
	}



	timer off 1
	timer list 1


*exit, STATA clear









