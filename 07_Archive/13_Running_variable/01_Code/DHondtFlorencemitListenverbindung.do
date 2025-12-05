******************************************************
* Ex. D'Hondt running variable with Listenverbindung *
******************************************************

******************************************
* A) Program for Seat Allocation D'Hondt *
******************************************

program DHondtLV

* Allocate seats to Listenverbindungen

save "Listenverbindungen.dta", replace

* Collapse votes by listenverbindung
collapse votesSimulated  votes (first) party, by(listenverbindung)

* Generate Bruchzahl
forvalues h = 1(1)$Seattotal {
	generate Bruchzahl`h' = votesSimulated/`h'
}

reshape long Bruchzahl, i(listenverbindung) j(number)


* Assign seats to lists
egen Bruchzahl_rank=rank(Bruchzahl), field

gen Seat=0 
replace Seat=1 if Bruchzahl_rank<=$Seattotal


* Only if the considered party is within a Listenverbindung
if $LV==1 {
destring listenverbindung, gen(lala) force
	
	
* Without Listenverbindung
preserve 
keep if lala != .
keep party number votesSimulated votes Seat
save "SeatsNoListenverbindung.dta", replace
restore
	
	
* With Listenverbindung
preserve 
keep if lala ==.
collapse (sum) Seat, by(listenverbindung)
rename Seat Seats
save "SeatsListenverbindung.dta", replace
restore


* With Listenverbindung: Allocate Seats to Lists
clear
use "Listenverbindungen.dta"
merge m:1 listenverbindung using SeatsListenverbindung.dta
drop if _merge==1
drop _merge
	
		
* Generate Bruchzahl
forvalues i = 1(1)$Seattotal {
	generate Bruchzahl`i' = votesSimulated/`i'
}

reshape long Bruchzahl, i(party) j(number) string
destring number, replace
	
* Assign seats to lists
bysort listenverbindung: egen Bruchzahl_rank=rank(Bruchzahl), field
	
gen Seat=0 
replace Seat=1 if Bruchzahl_rank<=Seats

keep party number votesSimulated votes Seat listenverbindung

* Add parties without Listenverbindung
append using SeatsNoListenverbindung.dta
}

else{
}
end





clear
cap log close
set more off
version 15

timer on 1

**********************
* B) Seat allocation *
**********************

use "Data/ExRunningVariableListenverbindung.dta", clear
// use "/Users/Florence/Dropbox/Projekt Nationalräte/Running variable/ExamplesMS/Data/ExRunningVariableListenverbindung.dta", clear
gen votesSimulated=votes

* Recode lists with no listenverbindung as a separate listenverbindung
egen Listenverbindungen_help=group(party) if listenverbindung=="" 
tostring Listenverbindungen_help, replace
replace listenverbindung=Listenverbindungen_help if listenverbindung==""
	
sum totseat
global Seattotal= r(max)

* In the first round we use the code for parties within a listenverbindung
global LV=1
qui DHondtLV
global LV=0

egen ActualSeats=sum(Seat), by(party)

destring party, replace
egen Partymax = max(party)
local Partymax = Partymax
drop Partymax

save "ActualAllocation", replace
* This is the actual seat allocation


*************************************
* C) Vote Change to get/lose a seat *
*************************************

clear
local j=1  //Party

while `j'<=`Partymax'{

	use "ActualAllocation"
	sum(ActualSeats) if party==`j'
	local ActualSeats=r(max)
	
	* Check if the considered party is in a listenverbindung: LV=1 yes / LV=0 no
	gen LV=0
	replace LV = 1 if party ==`j' & !missing(listenverbindung)
	sum(LV)
	global LV=r(max)


	* Get a Seat
	* Define xhigh start
	local i=`ActualSeats'+1  //first seat "to get"
	gen xhigh = 0.1*votes if party==`j' //add 10% of the votes
	sum (xhigh)
	local xhigh=r(max)
	drop xhigh


	* Check if xhigh start is already too high
	local xhigh = 2*`xhigh'
	local Check=0

	while `Check'==0{
		egen Listenverbindungen_help=group(party) if listenverbindung=="" 
		tostring Listenverbindungen_help, replace
		replace listenverbindung=Listenverbindungen_help if listenverbindung==""
		collapse votesSimulated  votes (first) listenverbindung, by(party)
		replace votesSimulated=votes+`xhigh' if party==`j'
		
		qui DHondtLV

		* Check
		gen Check=0
		replace Check=1 if Seat==0 & party==`j' & number==`i'
		sum Check
		local Check=r(max)	
		local xhigh = `xhigh'/2
	}
			
			
	while `i'<=$Seattotal{
	local Change=0

		while `Change'==0{
	
			local xhigh = 2*`xhigh'
	
			egen Listenverbindungen_help=group(party) if listenverbindung=="" 
			tostring Listenverbindungen_help, replace
			replace listenverbindung=Listenverbindungen_help if listenverbindung==""
			collapse votesSimulated  votes (first) listenverbindung, by(party)
			replace votesSimulated=votes+`xhigh' if party==`j'
			
			qui DHondtLV

			* Check
			gen Change=0
			replace Change=1 if Seat==1 & party==`j' & number==`i'
			sum Change
			local Change=r(max)	
			local xlow = `xhigh'/2
		}
	

		while `xhigh'-`xlow'>1 {
			local x = round((`xlow'+`xhigh')/2)
			
			egen Listenverbindungen_help=group(party) if listenverbindung=="" 
			tostring Listenverbindungen_help, replace
			replace listenverbindung=Listenverbindungen_help if listenverbindung==""
			collapse votesSimulated  votes (first) listenverbindung, by(party)
			replace votesSimulated=votes+`xhigh' if party==`j'
		
			qui DHondtLV

			* Check
			gen toohigh=1 if (Seat==1 & party==`j' & number==`i')
			gen toolow=1  if (Seat==0 & party==`j' & number==`i')
	
			egen toohighmax=max(toohigh)
			egen toolowmax=max(toolow)
			
			if toohighmax==1 local xhigh=`x' 
			if toolowmax==1 local xlow=`x'
		}

		* Save
		preserve
		keep if party==`j' & number==`i'
		keep party number votes votesSimulated
		gen Partymargin=votesSimulated-votes	
		
		if `i'==`ActualSeats'+1 & `j'==1 {
			save "PartymarginDHondtLVFlorence.dta", replace
		}
			
		else {
			append using "PartymarginDHondtLVFlorence.dta"
			sort party number
			save "PartymarginDHondtLVFlorence.dta", replace
		}	
		restore
		
		* Change indicator for seat
		local ++i
	}


	

if `ActualSeats'>0{	//only if the considered party got at least one seat
* Lose a Seat	
local i=`ActualSeats'
gen xhigh = 0.1*votes if party==`j' //subtract 10% of the votes
sum (xhigh)
local xhigh=r(max)
drop xhigh


* Check if xhigh start is already too high
local xhigh = 2*`xhigh'
local Check=0

	while `Check'==0{
		egen Listenverbindungen_help=group(party) if listenverbindung=="" 
		tostring Listenverbindungen_help, replace
		replace listenverbindung=Listenverbindungen_help if listenverbindung==""
		collapse votesSimulated  votes (first) listenverbindung, by(party)
		replace votesSimulated=votes-`xhigh' if party==`j'
			
		qui DHondtLV

		* Check
		gen Check=0
		replace Check=1 if Seat==1 & party==`j' & number==`i'
		sum Check
		local Check=r(max)	
		local xhigh = `xhigh'/2
	}


	while `i'>0{
	local Change=0

		while `Change'==0{
	
			local xhigh = 2*`xhigh'
		
			egen Listenverbindungen_help=group(party) if listenverbindung=="" 
			tostring Listenverbindungen_help, replace
			replace listenverbindung=Listenverbindungen_help if listenverbindung==""
			collapse votesSimulated  votes (first) listenverbindung, by(party)
			replace votesSimulated=votes-`xhigh' if party==`j'
			
			DHondtLV
			
			* Check
			gen Change=0
			replace Change=1 if Seat==0 & party==`j' & number==`i'
			sum Change
			local Change=r(max)	
			local xlow = `xhigh'/2
		}
	

		while `xhigh'-`xlow'>1 {
			local x = round((`xlow'+`xhigh')/2)
			
			egen Listenverbindungen_help=group(party) if listenverbindung=="" 
			tostring Listenverbindungen_help, replace
			replace listenverbindung=Listenverbindungen_help if listenverbindung==""
			collapse votesSimulated  votes (first) listenverbindung, by(party)
			replace votesSimulated=votes-`xhigh' if party==`j'
			
			qui DHondtLV

			* Check
			gen toohigh=1 if Seat==0 & party==`j' & number==`i'
			gen toolow=1  if Seat==1 & party==`j' & number==`i'
	
			egen toohighmax=max(toohigh)
			egen toolowmax=max(toolow)
			
			if toohighmax==1 local xhigh=`x' 
			if toolowmax==1 local xlow=`x'
		}

		* Save
		preserve
		keep if party==`j' & number==`i'
		keep party number votes votesSimulated
		gen Partymargin=votesSimulated-votes	
		
		append using "PartymarginDHondtLVFlorence.dta"
		sort party number
		save "PartymarginDHondtLVFlorence.dta", replace
		restore

		* Change indicator for seat
		local --i	
	}
}	
		* Change indicator for party
		local ++j
		clear
}
	
erase "ActualAllocation.dta"
erase "SeatsNoListenverbindung.dta"
erase "SeatsListenverbindung.dta"
erase "Listenverbindungen.dta"


timer off 1
timer list 1
