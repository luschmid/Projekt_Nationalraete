clear
cap log close
set more off
version 15
ssc install moreobs

cd "D:\SchmidLu\Dropbox\Projekt Nationalräte\Running variable\ExamplesMS"
*cd "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\Running variable\ExamplesMS"

********************************
* Ex. D'Hondt running variable *
********************************

* (i) Input data

local n=3 // define total number of seats in district

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

bysort party: egen cand=rank(-votes), unique // generate rank of candidate on party list

save example.dta, replace

* (ii) Allocation of seats to alliances

collapse (sum) votes (first) party, by(alliance)

forvalues k = 1(1)`n' {
generate quot`k' = votes/`k'
}

reshape long quot, i(alliance) j(seats_l)

gen pm_1=.

levelsof alliance, local(levels)
foreach l of local levels {
	preserve
	keep if alliance!=`l'
	drop seats_l alliance
	egen seats_l=rank(-quot),  unique // generate rank of all other quotients except alliance a
	rename quot quot_others
	keep if seats_l<=`n'
	replace seats_l=`n'-seats_l+1
	gen alliance=`l'
	save help_others.dta, replace
	restore
	merge 1:1 seats_l alliance using help_others.dta
	replace pm_1=votes-(seats_l*quot_others) if alliance==`l'
	drop quot_others _merge
}

drop votes quot
save example_pm1.dta, replace

* (iii) Allocation of seats within alliances

use example.dta, clear
collapse (sum) votes (first) alliance, by(party)

duplicates tag alliance , gen(dups)
drop if dups==0
drop dups

forvalues i = 1(1)`n' {
generate quot`i' = votes/`i'
}

reshape long quot, i(party) j(number)
save example_pm2_all.dta, replace
 
levelsof alliance, local(alevels)
foreach l of local alevels {
use example_pm2_all, clear
keep if alliance==`l'
	forvalues k = 1(1)`n'{
		gen pm_2`k'=.
		levelsof party, local(plevels)
		foreach j of local plevels {
			preserve
			keep if party!="`j'"
			drop number party
			egen number=rank(-quot),  unique // generate rank of all other quotients except party j 
			rename quot quot_others
			keep if number<=`k'
			replace number=`k'-number+1
			gen party="`j'"
			save help_others.dta, replace
			restore
			merge 1:1 number party using help_others.dta, update
			replace pm_2`k'=votes-(number*quot_others) if party=="`j'"
			drop quot_others _merge
		}
	}

reshape long pm_2, i(party number) j(seats_l)
drop votes quot

keep if number <= seats_l // drop all impossible combinations (#seats alliance <  #seats for party in alliance)

save example_pm2_`l'.dta, replace
}
erase example_pm2_all.dta

clear
foreach l of local alevels {
append using example_pm2_`l'.dta
erase example_pm2_`l'.dta
}

* (iv) Merge pm1 and pm2

merge m:1 seats_l alliance using example_pm1.dta, gen(merge)
replace number = seats_l if number == .
drop merge
erase example_pm1.dta

sort alliance party seats_l number
order party alliance number seats_l pm_1 pm_2

sort party number seats_l
egen partymarg=rowmin(pm_1 pm_2)
collapse (max) partymarg , by(alliance party number)
save example_partymarg.dta, replace

* (v) Generate vote margin

use example.dta, clear

levelsof party, local(plevels)
foreach j of local plevels {
	use example.dta, clear
	gen votesdiff=.
	keep if party=="`j'"
	save example_p`j'.dta, replace
	sum cand
	local last_cand=r(max)
	if `last_cand'>1{ // Case 1: More than one candidate on a party list
	levelsof cand, local(clevels)
	foreach h of local clevels {
		use example_p`j'.dta, clear
		keep if cand!=`h'
		drop cand
		egen number=rank(-votes),  unique // generate rank of all other quotients except party j 
		rename votes votes_others
		gen cand=`h'
		save help_others.dta, replace
		use example_p`j'.dta, clear
		joinby cand party using help_others.dta, update
		replace votesdiff=votes-votes_others if cand==`h' & party=="`j'"
		expand 2 in 1, gen(newobs) // Add observations for candidate margin for nth potential seat
		replace votesdiff = votes if newobs==1
		replace number = `last_cand' if newobs==1
		drop votes_others newobs
	if `h' == 1 save example_cm_`j'.dta, replace
	else {
		append using example_cm_`j'.dta
		save example_cm_`j'.dta, replace
	}
	}
	}
	else { // Case 2: Only one candidate on a party list (individual margin does not exist and only party margin is binding)
	use example_p`j'.dta
	gen number=1
	save example_cm_`j'.dta, replace
	}
	erase example_p`j'.dta
}
erase help_others.dta

clear
foreach j of local plevels {
append using example_cm_`j'.dta
erase example_cm_`j'.dta
}

merge m:1 alliance party number using example_partymarg.dta
drop if _merge==2 // drop observations for which no individual margin exists (cand <n)
erase example_partymarg.dta

egen votemargin = rowmin(votesdiff partymarg)
replace partymarg = . if cand!=number
collapse (min) partymarg (max) votemargin, by(alliance party cand votes)
sort alliance party cand
save example.dta, replace

