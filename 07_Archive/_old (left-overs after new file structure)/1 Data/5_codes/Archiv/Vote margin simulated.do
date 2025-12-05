******************************************************
* Ex. D'Hondt running variable with Listenverbindung *
******************************************************

// Potential issues: 1. Start values for votes (local Votes_h=1 // start value for votes)
// 					 2. Rounding of votes in loop 

cd "D:\SchmidLu\Dropbox\Projekt Nationalräte\Running variable\ExamplesMS"
*cd "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\Running variable\ExamplesMS"

******************************************
* A) Set parameters							 *
******************************************

global n 4 // define total number of seats in district
global convcrit 0.00005
local can "SZ"
local yr 2011

******************************************
* B) Read-in Data 						 *
******************************************

use "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\nationalraete_1931_2015.dta", clear
*use "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\1 Data\1_data\nationalraete_1931_2015.dta", clear

drop if name==""
replace list=name if list==""
keep if canton=="`can'" & year==`yr'
rename votes votes_h
rename pvotes votes_j
rename list party

sort alliance suballiance party
egen party_num=group(party) // numeric version of party
egen alliance_num=group(alliance) // numeric version of alliance
egen suballiance_num=group(alliance suballiance) // numeric version of suballiance

egen help=group(party) if alliance_num==. // generate alliance_num for all parties that are not in an alliance
sum alliance_num
if r(max)!=. replace alliance_num=r(max)+help if alliance_num==.
if r(max)==. replace alliance_num=help if alliance_num==.
drop help

egen help=group(party) if suballiance_num==.  // generate suballiance_num for all parties that are not in an suballiance
sum suballiance_num
if r(max)!=. replace suballiance_num=r(max)+help if suballiance_num==.
if r(max)==. replace suballiance_num=alliance_num if suballiance_num==.
drop help

gsort alliance party -votes_h
bysort party: gen id=_n // indicator for alliance

keep ID id canton year name firstname alliance* suballiance* party party_num elected votes_h votes_j

bysort party: egen cand=rank(-votes_h), unique // generate rank of candidate on party list

save example.dta, replace

******************************************
* C) Program for Seat Allocation D'Hondt *
******************************************

clear all
capture erase "PartymarginDHondt.dta"

timer on 1

program DHondt

save "example_input.dta", replace

unique party_num if alliance_num==`1'
local Al_size=r(unique)

unique suballiance_num if alliance_num==`1'
local Sal_l_size=r(unique)

unique party_num if alliance_num==`1'
local P_l_size=r(unique)

unique party_num if suballiance_num==`1'
local Sal_size=r(unique)


* (i) Allocation of seats to alliances

collapse (first) votes_j alliance_num, by(party_num)
collapse (sum) votes_j (first) party_num, by(alliance_num)

replace party_num=`3' if alliance_num==`1'

forvalues k = 1(1)$n {
	generate quot`k' = votes_j/`k'
}

reshape long quot, i(alliance_num) j(number)

egen quot_rank=rank(-quot), unique

gen seat=0 
replace seat=1 if quot_rank<=$n

* (ii) Allocation of alliance seats to suballiances

if `Al_size'>1 & `Al_size'!=. {

keep if alliance_num ==`1'
collapse (sum) seats_l=seat, by(alliance_num)
save "SeatsAlliance.dta", replace

use "example_input.dta", clear
collapse (first) votes_j alliance_num suballiance_num, by(party_num)
collapse (sum) votes_j (first) alliance_num party_num, by(suballiance_num)
merge m:1 alliance_num using SeatsAlliance.dta
keep if _merge==3
		
forvalues k = 1(1)$n {
	generate quot`k' = votes_j/`k'
}

reshape long quot, i(suballiance_num) j(number) 
egen quot_rank=rank(-quot), unique 
	
gen seat=0 
replace seat=1 if quot_rank<=seats_l
keep suballiance_num party_num number votes_j seat alliance_num 

* (iii) Allocation of suballiance seats to parties

if `Sal_l_size'<`P_l_size'{

keep if suballiance_num ==`2'
collapse (sum) seats_s=seat, by(suballiance_num)
save "SeatsSuballiance.dta", replace

use "example_input.dta", clear
collapse (first) votes_j (first) suballiance_num, by(party_num)
merge m:1 suballiance_num using SeatsSuballiance.dta
keep if _merge==3
		
forvalues k = 1(1)$n {
	generate quot`k' = votes_j/`k'
}

reshape long quot, i(party_num) j(number)
	
egen quot_rank=rank(-quot), unique
	
gen seat=0 
replace seat=1 if quot_rank<=seats_s

keep party_num number votes_j seat  

}

if `Sal_size'>1 & `Sal_size'!=. {

use "example_input.dta", clear
collapse (first) votes_j (first) alliance_num, by(party_num)
merge m:1 alliance_num using SeatsAlliance.dta
keep if _merge==3
		
forvalues k = 1(1)$n {
	generate quot`k' = votes_j/`k'
}

reshape long quot, i(party_num) j(number)
	
egen quot_rank=rank(-quot), unique
	
gen seat=0 
replace seat=1 if quot_rank<=seats_l

keep party_num number votes_j seat  
}
}

else{
}
end

clear
cap log close
set more off
version 15

timer on 1

*************************************
* D) Vote Change to get/lose a seat *
*************************************

use "example", clear

levelsof party_num, local(jlevels)
foreach j of local jlevels {
use "example", clear
levelsof id if party_num ==`j', local(hlevels)
foreach h of local hlevels {

	use "example", clear

	qui sum alliance_num if party_num ==`j'
	local l=r(max)
	qui sum suballiance_num if party_num ==`j'
	local s=r(max)
	
	di "party:" "`j'"
	
	qui sum votes_h if party_num ==`j' & id==`h' 
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
			replace votes_h=`Votes_h' if party_num ==`j' & id==`h'
			replace votes_j=`Votes_orig_j' - `Votes_orig_h' + `Votes_h' if party_num ==`j'
			egen votes_h_rank=rank(-votes_h) if party_num==`j' , unique
			qui sum votes_h_rank if party_num==`j' & id==`h'
			local Rank_h=r(max)
			qui DHondt `l' `s' `j'
			qui sum seat if seat==1 & party_num==`j'
			local No_Seats=r(N)
			di "Votes_h: " "`Votes_h'"
		}
	
			local xlow = `Votes_h'/2 
			local xhigh = `Votes_h'
		
		while `xhigh'-`xlow'> $convcrit { // (2) Find point within intervall in which candidate h barely gets a seat
			use "example", clear
			local Votes_h = (`xlow'+`xhigh')/2 // note: go to the center of the last step increase 
			replace votes_h=`Votes_h' if party_num ==`j' & id==`h'
			replace votes_j=`Votes_orig_j' - `Votes_orig_h' + `Votes_h' if party_num ==`j'
			egen votes_h_rank=rank(-votes_h) if party_num==`j' , unique
			sum votes_h_rank if party_num==`j' & id==`h'
			local Rank_h=r(max)
			qui DHondt `l' `s' `j'
			qui sum seat if seat==1 & party_num==`j'
			local No_Seats=r(N)
			di "Votes_h: " "`Votes_h'"
			if `Rank_h'>`No_Seats' local xlow=`Votes_h' 
			if `Rank_h'<=`No_Seats' local xhigh=`Votes_h' 
		}

		* Save
		preserve
		clear
		set obs 1
		gen votes_h_req=`Votes_h'
		gen id=`h'
		gen party_num=`j'
		
		if `j' == 1 & `h'==1{
		save "PartymarginDHondt.dta", replace
		}
		else{
		append using "PartymarginDHondt.dta"
		save "PartymarginDHondt.dta", replace
		}
		restore	
	

}
}
	
erase "SeatsAlliance.dta"
erase "SeatsSuballiance.dta"
erase "example_input.dta"

timer off 1
timer list 1

*************************************
* E) Read-in results			    *
*************************************

use "example.dta", clear

merge 1:1 party_num id using "PartymarginDHondt.dta", gen(merge_new)

gen votemargin=votes_h-votes_h_req

replace votemargin = round(votemargin)
replace votes_h_req = round(votes_h_req)

keep ID id name first party_num votemargin elected party votes_h votes_h_req  votes_j  party alliance
order ID id name first party_num votemargin elected votes_h votes_h_req  votes_j party alliance

gsort party_num -votes_h
