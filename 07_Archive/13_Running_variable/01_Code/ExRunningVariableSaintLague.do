clear
cap log close
set more off
version 15

//cd "C:\Users\schelkem\Desktop\"
//cd "D:\SchmidLu\Dropbox\Projekt Nationalräte\Running variable\ExamplesMS"
cd "C:\Users\StempfeF\Dropbox\Projekt Nationalräte\Running variable\ExamplesMS"

********************************
* Ex. D'Hondt running variable *
********************************

* D'Hondt: seat allocation

use ExRunningVariable, clear

// Note that ExRunningVariable is the Wikipedia Example
// ExRunningVariable2 is a distribution of votes in which party A and B have the same remainder (0.5)

sum totseat
local imax = r(max)*2

forvalues i = 1(2) `imax' {
generate Bruchzahl`i' = votes/`i'
}

reshape long Bruchzahl, i(party) j(number)
destring number, replace
replace number=(number+1)/2

* (iv) Assign seats to lists

egen Bruchzahl_rank=rank(Bruchzahl), field

gen Seat=0 
replace Seat=1 if Bruchzahl_rank<=totseat

gsort -Bruchzahl

bysort number: egen SumVotes=sum(votes)

* Running varibale: how many more votes party j would have needed/lost to capture a specific seat?

gen partymargin_ji = .

sum party
local jmax = r(max)
local imax = 8

forvalues j = 1(1)`jmax' {
	gen votes_n`j' = Bruchzahl if party!=`j'   // take vector of bruchzahlen of other than your own party (votes-j)
	egen order_n`j' = rank(-votes_n`j')        // rank() assigns average if two have the same rank
	replace order_n`j' = int(order_n`j')       // int(rank()) gives the same absolute rank to all with same rank
	forvalues i = 1(1)`imax' {
		local k = `imax'-`i'+1
		gen help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'
		replace help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'-1 & help_votes_nj_rank_`j'`i'  == .  // if two have the same rank, take the previous
		replace help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'-2 & help_votes_nj_rank_`j'`i'  == .  // if three have the same rank, take the previous
		replace help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'-3 & help_votes_nj_rank_`j'`i'  == .  // if four have the same rank, take the previous
		egen votes_nj_rank_`j'`i' = min(help_votes_nj_rank_`j'`i')
		drop help_votes_nj_rank_`j'`i'
		replace partymargin_ji = votes - ((2 * `i')-1) * votes_nj_rank_`j'`i' if party == `j' & number == `i'
	}
}

replace SumVotes=SumVotes-partymargin_ji

keep party number votes partymargin_ji Bruchzahl_rank SumVotes
gsort party number


save "PartymarginSaintLague", replace

