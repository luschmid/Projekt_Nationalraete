******************************************************
* Ex. D'Hondt running variable with Listenverbindung *
******************************************************


// Potential issues: 1. Start values for votes (local Votes=1 // start value for votes)
// 					 2. Rounding of votes in loop 


******************************************
* A) Entry Data here					 *
******************************************

global n=3 // define total number of seats in district

clear
input str2 party votes alliance
P1 60 1
P2 40 2
P3 30 2
end


bysort alliance: gen count=_N // indicator for alliance
gen al_member=count>1	
drop count

save example.dta, replace


******************************************
* B) Program for Seat Allocation D'Hondt *
******************************************

clear all

timer on 1

program DHondt

save "Alliance.dta", replace

* Collapse votes by alliance
collapse (sum) votes (first) party al_member, by(alliance)

replace party="`3'" if alliance==`2'

* Generate quotient
forvalues k = 1(1)$n {
	generate quot`k' = votes/`k'
}

reshape long quot, i(alliance) j(number)

* Assign seats to alliance
egen quot_rank=rank(-quot), unique

gen seat=0 
replace seat=1 if quot_rank<=$n

* Only if the considered party is within an alliance
if `1'==1 {

preserve 
sum alliance if party =="`3'" 
local Alliance_members=r(max)
keep if alliance==`Alliance_members'
collapse (sum) seat, by(alliance)
rename seat seats
save "SeatsAlliance.dta", replace
restore


* With alliance: allocate Seats to Lists
use "Alliance.dta", clear
merge m:1 alliance using SeatsAlliance.dta
drop if _merge==1
drop _merge
	
		
* Generate quotient
forvalues k = 1(1)$n {
	generate quot`k' = votes/`k'
}

reshape long quot, i(party) j(number) string
destring number, replace
	
* Assign seats to alliance
bysort alliance: egen quot_rank=rank(-quot), unique
	
gen seat=0 
replace seat=1 if quot_rank<=seats

keep party number votes seat alliance

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
* C) Vote Change to get/lose a seat *
*************************************

use "example", clear

levelsof party, local(levels)
foreach j of local levels {

	use "example", clear

	* Check if the considered party is in an alliance: Al_member=1 / Al_member=0 
	
	sum al_member if party =="`j'" 
	local Al_member=r(max)
	sum alliance if party =="`j'" 
	local l=r(max)

	local i=1  //seat 
	
	* Define  level for start
	local Votes=1 // start value for votes

	while `i'<=$n{
	local Change=0
	local No_Seats=`i'-1

		while `Change'==0{
			use "example", clear
			local Votes = `Votes'*2 
			replace votes=`Votes' if party=="`j'"
			qui DHondt `Al_member' "`l'" "`j'"
			sum seat if seat==1 & party=="`j'"
			local Change=r(N)-`No_Seats'
		}
	
			local xlow = `Votes'/2 
			local xhigh = `Votes'
		
		while `xhigh'-`xlow'>1 {
			use "example", clear
			local Votes = round((`xlow'+`xhigh')/2) // note: go to the center of the last step increase 
			replace votes=`Votes' if party=="`j'"
			qui DHondt `Al_member' "`l'" "`j'"
			sum seat if seat==1 & party=="`j'"
			local Change=r(N)-`No_Seats'

			if `Change'==0 local xlow=`Votes' 
			if `Change'>0 local xhigh=`Votes' 
		}

		* Save
		preserve
		keep if party=="`j'" & number==`i'
		keep party number votes
		
		if `i'==1 & "`j'"=="P1"{
			save "PartymarginDHondt.dta", replace
		}
			
		else {
			append using "PartymarginDHondt.dta"
			sort party number
			save "PartymarginDHondt.dta", replace
		}	
		restore
		
		* Change indicator for seat
		local ++i
	}


}
	
erase "SeatsAlliance.dta"
erase "Alliance.dta"


timer off 1
timer list 1
