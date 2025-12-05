********************************
* Ex. D'Hondt running variable *
********************************

******************************************
* A) Program for Seat Allocation D'Hondt *
******************************************

program DHondt

*Calculate Bruchzahlen
forvalues h = 1(1) $Seattotal {
generate Bruchzahl`h' = votesSimulated/`h'
}
reshape long Bruchzahl, i(party) j(number)

* Assign seats to lists
egen Bruchzahl_rank=rank(Bruchzahl), field
gen Seat=0 
replace Seat=1 if Bruchzahl_rank<=$Seattotal

end


clear
cap log close
set more off
version 15

timer on 1

**********************
* B) Seat allocation *
**********************

use "Data/ExRunningVariable.dta", clear
gen votesSimulated=votes

sum totseat
global Seattotal= r(max)

qui DHondt

egen ActualSeats=sum(Seat), by(party)

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
		collapse votes votesSimulated, by (party)
		replace votesSimulated=votes+`xhigh' if party==`j'
			
		qui DHondt

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
	
			collapse votes votesSimulated, by (party)
			replace votesSimulated=votes+`xhigh' if party==`j'
			
			qui DHondt

			* Check
			gen Change=0
			replace Change=1 if Seat==1 & party==`j' & number==`i'
			sum Change
			local Change=r(max)	
			local xlow = `xhigh'/2
		}
	

		while `xhigh'-`xlow'>1 {
			local x = round((`xlow'+`xhigh')/2)
			
			collapse votes votesSimulated, by (party)
			replace votesSimulated=votes+`x' if party==`j'
		
			qui DHondt

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
			save "Results/PartymarginDHondtFlorence.dta", replace
		}
			
		else {
			append using "Results/PartymarginDHondtFlorence.dta"
			sort party number
			save "Results/PartymarginDHondtFlorence.dta", replace
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
		collapse votes votesSimulated, by (party)
		replace votesSimulated=votes+`xhigh' if party==`j'
			
		qui DHondt

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
		
			collapse votes votesSimulated, by (party)
			replace votesSimulated=votes-`xhigh' if party==`j'
			
			qui DHondt
			
			* Check
			gen Change=0
			replace Change=1 if Seat==0 & party==`j' & number==`i'
			sum Change
			local Change=r(max)	
			local xlow = `xhigh'/2
		}
	

		while `xhigh'-`xlow'>1 {
			local x = round((`xlow'+`xhigh')/2)
			
			collapse votes votesSimulated, by (party)
			replace votesSimulated=votes-`x' if party==`j'
			
			qui DHondt

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
		
		append using "Results/PartymarginDHondtFlorence.dta"
		sort party number
		save "Results/PartymarginDHondtFlorence.dta", replace
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

timer off 1
timer list 1
