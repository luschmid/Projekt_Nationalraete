******************************************************
* Ex. D'Hondt running variable with Listenverbindung *
******************************************************


// Potential issues: 1. Start values for votes (local Votes=1 // start value for votes)
// 					 2. Rounding of votes in loop 


******************************************
* A) Set parameters							 *
******************************************

cd "D:\SchmidLu\Dropbox\Projekt Nationalräte\Running variable\ExamplesMS"

global n 3 // define total number of seats in district
global convcrit 0.0005

******************************************
* B) Read-in Data 						 *
******************************************

clear
input str2 party votes alliance
P1 30 1
P1 20 1
P1 10 1
P2 25 2
P2 10 2
P2 5 2
P3 15 2
P3 10 2
P3 5 2
end

egen party_num=group(party) // indicator for alliance
egen alliance_num=group(alliance) // indicator for alliance
bysort alliance: egen sd_party=sd(party_num) // indicator for alliance
gen al_member= sd_party >0 & !missing(sd_party)

gsort alliance party -votes
bysort party: gen id=_n // indicator for alliance

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

* Collapse votes by alliance
collapse (sum) votes (first) party_num al_member, by(alliance_num)
replace party=`2' if alliance_num==`1'

* Generate quotient
forvalues k = 1(1)$n {
	generate quot`k' = votes/`k'
}

reshape long quot, i(alliance_num) j(number)

* Assign seats to alliance
egen quot_rank=rank(-quot), unique

gen seat=0 
replace seat=1 if quot_rank<=$n

* Only if the considered party is within an alliance
if `Al_size'>1 & `Al_size'!=. { 

sum alliance_num if party_num ==`2'
local Alliance_members=r(max)
keep if alliance_num==`Alliance_members'
collapse (sum) seat, by(alliance_num)
rename seat seats
save "SeatsAlliance.dta", replace

* With alliance: allocate Seats to Lists
use "example_input.dta", clear
collapse (sum) votes (first) alliance_num, by(party_num)
merge m:1 alliance_num using SeatsAlliance.dta
keep if _merge==3
	
* Generate quotient
forvalues k = 1(1)$n {
	generate quot`k' = votes/`k'
}

reshape long quot, i(party_num) j(number) 
	
* Assign seats to alliance
egen quot_rank=rank(-quot), unique
	
gen seat=0 
replace seat=1 if quot_rank<=seats

keep party_num number votes seat alliance_num 

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

	* Check if the considered party is in an alliance: Al_member=1 / Al_member=0 
	
	qui sum alliance_num if party_num ==`j'
	local l=r(max)

	
	* Define  level for start
	local Votes=1 // start value for votes

	local Rank_h=$n
	local No_Seats=0

		while `Rank_h'>`No_Seats'{ // (1) Find intervall for which seat change for candidate h happens
			use "example", clear
			local Votes = `Votes'*2 
			replace votes=`Votes' if party_num ==`j' & id==`h'
			egen votes_rank=rank(-votes) if party_num==`j' , unique
			qui sum votes_rank if party_num==`j' & id==`h'
			local Rank_h=r(max)
			qui DHondt `l' `j'
			qui sum seat if seat==1 & party_num==`j'
			local No_Seats=r(N)
			di "Votes: " "`Votes'"
		}
	
			local xlow = `Votes'/2 
			local xhigh = `Votes'
		
		while `xhigh'-`xlow'> $convcrit { // (2) Find point within intervall in which candidate h barely gets a seat
			use "example", clear
			local Votes = (`xlow'+`xhigh')/2 // note: go to the center of the last step increase 
			replace votes=`Votes' if party_num==`j' & id==`h'
			egen votes_rank=rank(-votes) if party_num==`j' , unique
			sum votes_rank if party_num==`j' & id==`h'
			local Rank_h=r(max)
			qui DHondt `l' `j'
			qui sum seat if seat==1 & party_num==`j'
			local No_Seats=r(N)
			di "Votes: " "`Votes'"
			if `Rank_h'>`No_Seats' local xlow=`Votes' 
			if `Rank_h'<=`No_Seats' local xhigh=`Votes' 
		}

		* Save
		preserve
		clear
		set obs 1
		gen votes_req=`Votes'
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
erase "example_input.dta"

timer off 1
timer list 1

use "example.dta", clear

merge 1:1 party_num id using "PartymarginDHondt.dta"

gen votemargin=votes-votes_req

replace votemargin = round(votemargin)
replace votes_req = round(votes_req)

keep party  alliance_num votes votes_req  votemargin

gsort party -votes
