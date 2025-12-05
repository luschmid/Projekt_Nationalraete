
clear
cap log close
set more off
version 15

cd "D:\SchmidLu\Dropbox\Projekt Nationalräte\Running variable\ExamplesMS"

********************************
* Ex. D'Hondt running variable *
********************************

* (i) Input data

local n=3 // define total number of seats in district

clear
input str2 party votes alliance
P1 60 1
P2 40 2
P3 30 2
end

save example.dta, replace



* (ii) Allocation of seats to alliances

collapse  (sum) votes (first) party , by(alliance)

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

duplicates tag alliance , gen(dups)
drop if dups==0
drop dups

forvalues i = 1(1)`n' {
generate quot`i' = votes/`i'
}

reshape long quot, i(party) j(number)

forvalues k = 1(1)`n'{
gen pm_2`k'=.
levelsof party, local(levels)
foreach j of local levels {
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
merge 1:1 number party using help_others.dta
replace pm_2`k'=votes-(number*quot_others) if party=="`j'"
drop quot_others _merge
}
}

reshape long pm_2, i(party number) j(seats_l)

drop votes quot

drop if pm_2==. // drop all impossible combinations (#seats alliance <  #seats for party in alliance)

save example_pm2.dta, replace


* (iv) Merge pm1 and pm2

merge m:1 seats_l alliance using example_pm1.dta, gen(merge)
replace number = seats_l if number == .
drop merge

sort alliance party seats_l number
order party alliance number seats_l pm_1 pm_2

sort party number seats_l

egen pm_min=rowmin(pm_1 pm_2)

collapse (max) pm_min , by(alliance party number)

