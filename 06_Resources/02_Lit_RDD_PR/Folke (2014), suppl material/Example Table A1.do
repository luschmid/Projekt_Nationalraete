clear
cap log close
set more off
version 15

clear
set more off

*global path "D:\SchmidLu\Dropbox\Projekt Nationalräte"
*global path "C:\Dropbox\Projekt Nationalräte"
global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"

**************************************************************
*                                                            *
* Example from Folke (2013), Online Appendix Table A1, p. 19 *
*                                                            *
**************************************************************

* Modified Sainte-Laguë: Seat allocation

import exc using ///
	"$path\06_Resources\01_Lit_RDD_PR\Folke (2014), suppl material\Example Table A1.xlsx", ///
	cellra(A1:G10) first clear
g party = _n
g votes = voteshare * 100

sum seats
local totseat = r(sum)
local imax = r(sum)*2
local jmax = r(N)

forv i = 1(2)`imax' {
if `i' == 1 g quot`i' = votes/1.4
else g quot`i' = votes/`i'
}

reshape long quot, i(party) j(num)
replace num=(num+1)/2

egen quotrank=rank(quot), field

g seats_lss = 0 
replace seats_lss = 1 if quotrank <= `totseat'

g partymarg = .
local imax = `totseat' // SL: Previously it was "local imax = 8"

forv j = 1(1)`jmax' {
    g votes_n`j' = quot if party != `j' // Take vector of quotients of parties -j.
    egen order_n`j' = rank(-votes_n`j') // Assigns avrg for identical ranks.
    replace order_n`j' = int(order_n`j')
    forv i = 1(1)`imax' {               // SL: Why in steps of one and not two?
        local k = `imax'-`i'+1
        g help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'
        replace help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'-1 & help_votes_nj_rank_`j'`i'  == .  // if two have the same rank, take the previous
        replace help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'-2 & help_votes_nj_rank_`j'`i'  == .  // if three have the same rank, take the previous
        replace help_votes_nj_rank_`j'`i' = votes_n`j' if order_n`j' == `k'-3 & help_votes_nj_rank_`j'`i'  == .  // if four have the same rank, take the previous
        egen votes_nj_rank_`j'`i' = min(help_votes_nj_rank_`j'`i')
        drop help_votes_nj_rank_`j'`i'
        if `i' == 1 replace partymarg = votes - 1.4 * votes_nj_rank_`j'`i' if party == `j' & num == `i'
        if `i' > 1 replace partymarg = votes - ((2 * `i')-1) * votes_nj_rank_`j'`i' if party == `j' & num == `i'
	}
}

tsset party num

g dgain_i_lss = -1*partymarg if seats_lss == 0 & L.seats_lss == 1
g dlose_i_lss = partymarg if seats_lss == 1 & F.seats_lss == 0

collapse (first) voteshare seats_f (sum) seats_lss (first) dgain_i_f dgain_a_f ///
    (min) dgain_i_lss (first) dlose_i_f dlose_a_f (min) dlose_i_lss, ///
    by(partyname party)
replace dgain_i_lss = round(dgain_i_lss/100, 0.01)
replace dlose_i_lss = round(dlose_i_lss/100, 0.01)

sort party
drop party

