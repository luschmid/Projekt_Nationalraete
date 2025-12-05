******************************************************
* D'Hondt running variable with Listenverbindung *
******************************************************

// Potential issues: 1. Start values for votes (local Votes_h=1 // start value for votes)
// 					 2. Rounding of votes in loop 


capture cd "D:\SchmidLu\Dropbox\Projekt Nationalräte"
	if _rc==0{
		global path "D:\SchmidLu\Dropbox\Projekt Nationalräte"
}

capture cd "C:\Current\02. Academic\Dropbox\Projekt Nationalräte""
	if _rc==0{
		global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
}


******************************************
*  Set parameters							 *
******************************************

timer on 1
global convcrit 0.00005

do "$path\1 Data\5_codes\DHondt_Program.do" 

cd "$path\Running variable\Data"

//local can "SZ"
//local yr 2011

*********************************************************
* A) Read in seats data for loop over cantons and years *
*********************************************************

quietly{
preserve
clear
import excel "$path\1 Data\0_original_data\Verteilung der Nationalratssitze auf die Kantone.xlsx",  first sheet("out")
reshape long Seats_ct, i(canton) j(year)
keep if year<2019 
save "$path\1 Data\0_original_data\Seats_perCtYear.dta", replace
restore

use "$path\1 Data\1_data\nationalraete_1931_2015.dta", clear
merge m:1 canton year using "$path\1 Data\0_original_data\Seats_perCtYear.dta", gen(merge_seats)
keep if votes!=. // exclude all elections with "stillen Wahlen"
collapse (mean) Seats_ct , by(canton year)
save "Seats_perCtYear", replace
}

//local cans "SZ" //temporary
//local cans "AG AI AR BE BL BS FR GE GL GR JU LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH" //temporary
local cans "VD VS ZG ZH" //temporary

//qui levelsof canton, local(cans) // 1. Loop over cantons
foreach can of local cans{

qui use "Seats_perCtYear", clear
qui keep if canton=="`can'" 
di "NEW CANTON:" "`can'"

//local yrs 2015 //temporary
qui levelsof year, local(yrs) // 2. Loop over years (might be different for each canton depending on when stille Wahlen took place)
foreach yr of local yrs{

qui use "Seats_perCtYear", clear
qui sum Seats_ct if canton=="`can'" & year==`yr'
global n=r(max) // define total number of seats in district

di "NEW YEAR:" "`yr'"




******************************************
* B) Read-in Data 						 *
******************************************
quietly{
use "$path\1 Data\1_data\nationalraete_1931_2015.dta", clear

drop if name==""
replace list=name if list==""
keep if canton=="`can'" & year==`yr'
rename votes votes_h
rename pvotes votes_j
rename list party

replace alliance=party if alliance=="." | alliance==""
replace suballiance=party if suballianc=="." | suballiance==""

sort alliance suballiance party
egen party_num=group(party) // numeric version of party
egen alliance_num=group(alliance) // numeric version of alliance
egen suballiance_num=group(alliance suballiance) // numeric version of suballiance

gsort alliance party -votes_h
bysort party: gen id=_n // indicator for alliance

keep ID canton year name firstname alliance* suballiance* party party_num elected votes_h votes_j

replace votes_j=votes_h if votes_j==. // replace list votes with individual votes for lists with a single candidate

bysort party: egen cand=rank(-votes_h), unique // generate rank of candidate on party list

save example.dta, replace
}

*************************************
* D) Vote Change to get/lose a seat *
*************************************

use "example", clear

qui levelsof party_num, local(jlevels)
foreach j of local jlevels {
//di "New party number " "`j'" " with the following candidates:"
use "example", clear
qui levelsof cand if party_num ==`j', local(hlevels)
foreach h of local hlevels {

	
	use "example", clear
	qui levelsof ID, local(idh), if party_num ==`j' & cand==`h'
	//di  "New individual:" `idh'
	
	qui sum alliance_num if party_num ==`j'
	local l=r(max)
	qui sum suballiance_num if party_num ==`j'
	local s=r(max)
	
	
	qui sum votes_h if party_num ==`j' & cand==`h' 
	local Votes_orig_h=r(max)
	qui sum votes_j if party_num ==`j'
	local Votes_orig_j=r(max)
	
	* Define  level for start
	local Votes_h=1 // start value for votes
	local Votes_j= `Votes_orig_j' - `Votes_orig_h' + `Votes_h'
	
	local Rank_h=$n
	local No_Seats=0

		while `Rank_h'>`No_Seats'{ // (1) Find intervall for which seat change for candidate h happens
			use "example", clear
			local Votes_h = `Votes_h'*2 
			qui replace votes_h=`Votes_h' if party_num ==`j' & cand==`h'
			qui replace votes_j=`Votes_orig_j' - `Votes_orig_h' + `Votes_h' if party_num ==`j'
			qui egen votes_h_rank=rank(-votes_h) if party_num==`j' , unique
			qui sum votes_h_rank if party_num==`j' & cand==`h'
			local Rank_h=r(max)
			DHondt `l' `s' `j'
			qui sum seat if seat==1 & party_num==`j'
			local No_Seats=r(N)
			//di "Votes_h: " "`Votes_h'"
			//di "Votes_j: " "`Votes_orig_j' - `Votes_orig_h' + `Votes_h'"
		}
	
			local xlow = `Votes_h'/2 
			local xhigh = `Votes_h'
		
		while `xhigh'-`xlow'> $convcrit { // (2) Find point within intervall in which candidate h barely gets a seat
			use "example", clear
			local Votes_h = (`xlow'+`xhigh')/2 // note: go to the center of the last step increase 
			qui replace votes_h=`Votes_h' if party_num ==`j' & cand==`h'
			qui replace votes_j=`Votes_orig_j' - `Votes_orig_h' + `Votes_h' if party_num ==`j'
			qui egen votes_h_rank=rank(-votes_h) if party_num==`j' , unique
			qui sum votes_h_rank if party_num==`j' & cand==`h'
			local Rank_h=r(max)
			//di "`Rank_h'"
			DHondt `l' `s' `j'
			qui sum seat if seat==1 & party_num==`j'
			local No_Seats=r(N)
			//di "`No_Seats'"
			qui sum votes_j if party_num ==`j'
			//di "Votes_h: " "`Votes_h'"
			//di "Votes_j: " "`Votes_orig_j' - `Votes_orig_h' + `Votes_h'"
			if `Rank_h'>`No_Seats' local xlow=`Votes_h' 
			if `Rank_h'<=`No_Seats' local xhigh=`Votes_h' 
		}

		* Save
		
		
		preserve
		clear
		qui set obs 1
		qui gen votes_h_req=`Votes_h'
		qui gen number=`h'
		qui gen ID=`idh'
		qui gen party_num=`j'
		qui gen year=`yr'
		qui gen canton="`can'"
		qui gen votes=`Votes_orig_h'
		
		if `j' == 1 & `h'==1 & `yr'==1931 & "`can'"=="AG" {
		qui save "RV_Simuation_all.dta", replace
		}
		else{
		qui append using "RV_Simuation_all.dta"
		qui save "RV_Simuation_all.dta", replace
		}
		restore	

		} // end of h-loop
		} // end of j-loop

} // end of year loop
} // end of cantonal loop

	
erase "SeatsAlliance.dta"
erase "SeatsSuballiance.dta"
erase "example_input.dta"

timer off 1
timer list 1

*************************************
* E) Read-in results			    *
*************************************

//use "$path\1 Data\1_data\nationalraete_1931_2015.dta", clear

use "RV_Analytical_all.dta", clear
rename votemargin votemargin_ana
merge 1:1 ID year using "RV_Simuation_all.dta", gen(merge_new)
gen votemargin_sim=votes_h-votes_h_req
 
gen test=votemargin_ana-votemargin_sim
tab test

gen votemargin=votes_h-votes_h_req

replace votemargin = round(votemargin)
replace votes_h_req = round(votes_h_req)

keep ID id name first party_num votemargin elected party votes_h votes_h_req  votes_j  party alliance
order ID id name first party_num votemargin elected votes_h votes_h_req  votes_j party alliance

gsort party_num -votes_h
