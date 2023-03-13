

*********************************************************
* A) Read in seats data for loop over cantons and years *
*********************************************************


preserve
clear
import excel "$path\01_Raw_data\07_Cantons\Verteilung der Nationalratssitze auf die Kantone.xlsx", clear  first sheet("out")
reshape long Seats_ct, i(canton) j(year)
keep if year<2019 
save "$path\02_Processed_data\07_Cantons\Seats_perCtYear.dta", replace
restore

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
merge m:1 canton year using "$path\02_Processed_data\07_Cantons\Seats_perCtYear.dta", gen(merge_seats)
keep if votes!=. // exclude all elections with "stillen Wahlen"
collapse (mean) Seats_ct , by(canton year)
save "$path\02_Processed_data\07_Cantons\Seats_perCtYear.dta", replace

levelsof canton, local(cans) // 1. Loop over cantons
//local cans "VS" // temporary
foreach can of local cans{

use "$path\02_Processed_data\07_Cantons\Seats_perCtYear.dta", clear
keep if canton=="`can'" 

levelsof year, local(yrs) // 2. Loop over years (might be different for each canton depending on when stille Wahlen took place)
//local yrs 1967
foreach yr of local yrs{

use "$path\02_Processed_data\07_Cantons\Seats_perCtYear.dta", clear
sum Seats_ct if canton=="`can'" & year==`yr'
local n=r(max) // define total number of seats in district
di "`n'"
di "New year:" "`yr'"



******************************************
* B) Read-in Data 						 *
******************************************

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear

drop if name==""
replace list=name if list==""
drop party party_num
gen party=list
egen party_num=group(list canton year)

keep if canton=="`can'" & year==`yr'
rename pvotes votes_j

replace votes_j=votes_h if votes_j==. // replace list votes with individual votes for lists with a single candidate


replace alliance=party if alliance=="." | alliance==""
replace suballiance=party if suballianc=="." | suballiance==""

sort alliance suballiance party
egen alliance_num=group(alliance) // numeric version of alliance
egen suballiance_num=group(alliance suballiance) // numeric version of suballiance

gsort alliance party -votes_h
bysort party_num: gen id=_n // indicator for alliance

keep ID canton year name firstname alliance* suballiance* party party_num elected votes_h votes_j

bysort party_num: egen cand=rank(-votes_h), unique // generate rank of candidate on party list

save "$path\02_Processed_data\13_Running_variable\example.dta", replace



******************************************
* C) Calculate analytical margin		 *
******************************************

unique alliance_num
local Al_size=r(unique)

unique alliance_num suballiance_num
local Sal_size=r(unique)

unique party_num
local P_size=r(unique)

unique ID
local H_size=r(unique)

di "`Sal_size'"
di "`Al_size'"
di "`P_size'"
di "`H_size'"

* (i) Allocation of seats to alliances

collapse (first) votes_j alliance alliance_num suballiance_num,  by(party_num)
collapse (sum) votes_j (first) party_num suballiance_num, by(alliance_num)

forvalues k = 1(1)`n' {
generate quot`k' = votes_j/`k'
}

reshape long quot, i(alliance_num) j(seats_l)
label drop _all

gen pm_1=.

levelsof alliance_num, local(levels)
foreach l of local levels {
	preserve
	keep if alliance_num!=`l'
	drop seats_l alliance_num
	egen seats_l=rank(-quot),  unique // generate rank of all other quotients except alliance l
	rename quot quot_others
	keep if seats_l<=`n'
	replace seats_l=`n'-seats_l+1
	gen alliance_num=`l'
	save "$path\02_Processed_data\13_Running_variable\help_others.dta", replace
	restore
	merge 1:1 seats_l alliance_num using "$path\02_Processed_data\13_Running_variable\help_others.dta"
	replace pm_1=votes_j-(seats_l*quot_others) if alliance_num==`l'
	drop quot_others _merge
}

drop votes_j quot

save "$path\02_Processed_data\13_Running_variable\example_pm1.dta", replace



* (ii) Allocation of alliance seats to suballiances

if `Al_size'<`P_size'{

use "$path\02_Processed_data\13_Running_variable\example.dta", clear
collapse (first) votes_j alliance suballiance alliance_num suballiance_num, by(party_num)
collapse (sum) votes_j (first) alliance alliance_num party_num, by(suballiance_num)

duplicates tag alliance_num, gen(dups)
drop if dups==0
drop dups

forvalues i = 1(1)`n' {
generate quot`i' = votes_j/`i'
}

reshape long quot, i(suballiance_num) j(seats_s)
save "$path\02_Processed_data\13_Running_variable\example_pm2_all.dta", replace
 
levelsof alliance_num, local(alevels)
foreach l of local alevels {
use "$path\02_Processed_data\13_Running_variable\example_pm2_all.dta", clear
keep if alliance_num==`l'
	forvalues k = 1(1)`n'{
		gen pm_2`k'=.
		levelsof suballiance_num, local(salevels)
		foreach s of local salevels {
			preserve
			keep if suballiance_num!=`s'
			drop seats_s suballiance_num
			egen seats_s=rank(-quot), unique // generate rank of all other quotients except suballiance s 
			rename quot quot_others
			keep if seats_s<=`k'
			replace seats_s=`k'-seats_s+1
			gen suballiance_num=`s'
			save "$path\02_Processed_data\13_Running_variable\help_others.dta", replace
			restore
			merge 1:1 seats_s suballiance_num using "$path\02_Processed_data\13_Running_variable\help_others.dta", update
			replace pm_2`k'=votes_j-(seats_s*quot_others) if suballiance_num==`s'
			drop quot_others _merge
		}
	}

reshape long pm_2, i(suballiance_num seats_s) j(seats_l)
label drop _all
drop votes_j quot

keep if seats_s <= seats_l // drop all impossible combinations (#seats alliance <  #seats for suballiance)

save "$path\02_Processed_data\13_Running_variable\example_pm2_`l'.dta", replace
}
erase "$path\02_Processed_data\13_Running_variable\example_pm2_all.dta"

clear
foreach l of local alevels {
append using "$path\02_Processed_data\13_Running_variable\example_pm2_`l'.dta"
erase "$path\02_Processed_data\13_Running_variable\example_pm2_`l'.dta"
}

save "$path\02_Processed_data\13_Running_variable\example_pm2.dta", replace

}

* (iii) Allocation of suballiance seats to parties


if `Sal_size'<`P_size'{

use "$path\02_Processed_data\13_Running_variable\example.dta", clear
collapse (first) votes_j alliance_num suballiance_num, by(party_num)

duplicates tag suballiance_num, gen(dups)
drop if dups==0
drop dups



forvalues i = 1(1)`n' {
generate quot`i' = votes_j/`i'
}

reshape long quot, i(party_num) j(number)
label drop _all

save "$path\02_Processed_data\13_Running_variable\example_pm3_all.dta", replace

levelsof suballiance_num, local(slevels)
foreach s of local slevels {
use "$path\02_Processed_data\13_Running_variable\example_pm3_all", clear
keep if suballiance_num==`s'
	forvalues k = 1(1)`n'{
		gen pm_3`k'=.
		levelsof party_num, local(plevels)
		foreach j of local plevels {
			preserve
			keep if party_num!=`j'
			drop number party_num
			egen number=rank(-quot), unique // generate rank of all other quotients except party j 
			rename quot quot_others
			keep if number<=`k'
			replace number=`k'-number+1
			gen party_num=`j'
			save "$path\02_Processed_data\13_Running_variable\help_others.dta", replace
			restore
			merge 1:1 number party_num using "$path\02_Processed_data\13_Running_variable\help_others.dta", update
			replace pm_3`k'=votes_j-(number*quot_others) if party_num==`j'
			drop quot_others _merge
		}
	}
reshape long pm_3, i(party_num number) j(seats_s)
drop votes_j quot

keep if number <= seats_s // drop all impossible combinations (#seats suballiance <  #seats for party in suballiance)

save "$path\02_Processed_data\13_Running_variable\example_pm3_`s'.dta", replace
}
erase "$path\02_Processed_data\13_Running_variable\example_pm3_all.dta"

clear
foreach s of local slevels {
append using "$path\02_Processed_data\13_Running_variable\example_pm3_`s'.dta"
erase "$path\02_Processed_data\13_Running_variable\example_pm3_`s'.dta"
}

gen seats_l=.

save "$path\02_Processed_data\13_Running_variable\help.dta", replace
clear

forvalues l = 1(1)`n'{
append using "$path\02_Processed_data\13_Running_variable\help.dta"
replace seats_l=`l' if seats_l==.
}

order alliance_num suballiance_num party_num seats_l seats_s number
sort alliance_num suballiance_num party_num seats_l seats_s number
drop if seats_s>seats_l

save "$path\02_Processed_data\13_Running_variable\example_pm3.dta", replace
}

* (iv) Merge pm1, pm2, and pm3

if `Sal_size'<`P_size'{
use "$path\02_Processed_data\13_Running_variable\example_pm3.dta", clear
merge m:1 seats_s seats_l suballiance_num using "$path\02_Processed_data\13_Running_variable\example_pm2.dta", gen(merge1)  // note that all merge1==2 are for parties which are not in a suballiance
replace number= seats_s if number==. // replace number for all parties that are not in a suballiance
merge m:1 seats_l alliance_num using "$path\02_Processed_data\13_Running_variable\example_pm1.dta", gen(merge2)
replace number=seats_l if number==. // replace number for all parties that are not in a alliance
}
else if `Al_size'<`P_size'{
use "$path\02_Processed_data\13_Running_variable\example_pm2.dta", clear
merge m:1 seats_l alliance_num using "$path\02_Processed_data\13_Running_variable\example_pm1.dta", gen(merge2)
gen pm_3=.
gen number= seats_s
replace number=seats_l if number==. // replace number for all parties that are not in a alliance
}
else{
use "$path\02_Processed_data\13_Running_variable\example_pm1.dta"
gen pm_2=.
gen pm_3=.
gen number= seats_l
}


capture drop merge*
capture erase "$path\02_Processed_data\13_Running_variable\example_pm1.dta "
capture erase "$path\02_Processed_data\13_Running_variable\example_pm2.dta"
capture erase "$path\02_Processed_data\13_Running_variable\example_pm3.dta"
capture erase "$path\02_Processed_data\13_Running_variable\help.dta"

egen partymarg=rowmin(pm_1 pm_2 pm_3)
replace number=1 if number==. // create number for candidates who run with no other candidates
collapse (max) partymarg, by(alliance_num alliance suballiance_num suballiance party_num number)
save "$path\02_Processed_data\13_Running_variable\example_partymarg.dta", replace

* (iv) Generate vote margin

use "$path\02_Processed_data\13_Running_variable\example.dta", clear

levelsof party_num, local(plevels)
foreach j of local plevels {
	use "$path\02_Processed_data\13_Running_variable\example.dta", clear
	gen votesdiff=.
	keep if party_num==`j'
	save "$path\02_Processed_data\13_Running_variable\example_p`j'.dta", replace
	sum cand
	local last_cand=r(max)
	if `last_cand'>1{ // Case 1: More than one candidate on a party list
	levelsof cand, local(clevels)
	foreach h of local clevels {
		use "$path\02_Processed_data\13_Running_variable\example_p`j'.dta", clear
		keep if cand!=`h'
		drop cand
		egen number=rank(-votes_h),  unique // generate rank of all other quotients except party j 
		rename votes_h votes_others
		gen cand=`h'
		save "$path\02_Processed_data\13_Running_variable\help_others.dta", replace
		use "$path\02_Processed_data\13_Running_variable\example_p`j'.dta", clear
		joinby cand party_num using "$path\02_Processed_data\13_Running_variable\help_others.dta", update
		replace votesdiff=votes_h-votes_others if cand==`h' & party_num==`j'
		expand 2 in 1, gen(newobs) // Add observations for candidate margin for nth potential seat
		replace votesdiff = votes_h if newobs==1
		replace number = `last_cand' if newobs==1
		drop votes_others newobs
	if `h' == 1 save "$path\02_Processed_data\13_Running_variable\example_cm_`j'.dta", replace
	else {
		append using "$path\02_Processed_data\13_Running_variable\example_cm_`j'.dta"
		save "$path\02_Processed_data\13_Running_variable\example_cm_`j'.dta", replace
	}
	}
	}
	else { // Case 2: Only one candidate on a party list (individual margin does not exist and only party margin is binding)
	use "$path\02_Processed_data\13_Running_variable\example_p`j'.dta"
	gen number=1
	replace votesdiff=votes_h if votesdiff==.
	save "$path\02_Processed_data\13_Running_variable\example_cm_`j'.dta", replace
	}
	erase "$path\02_Processed_data\13_Running_variable\example_p`j'.dta"
}
capture erase "$path\02_Processed_data\13_Running_variable\help_others.dta"

clear
foreach j of local plevels {
append using "$path\02_Processed_data\13_Running_variable\example_cm_`j'.dta"
erase "$path\02_Processed_data\13_Running_variable\example_cm_`j'.dta"
}

merge m:1 party_num alliance_num suballiance_num number using "$path\02_Processed_data\13_Running_variable\example_partymarg.dta", gen(merge_suballiance)

drop if merge_suballiance==2 // drop observations for which no individual margin exists (cand <n)
erase "$path\02_Processed_data\13_Running_variable\example_partymarg.dta"

egen votemargin = rowmin(votesdiff partymarg)
// replace partymarg = . if cand!=number
collapse (min) partymarg votesdiff (max) votemargin (first) number , by(cand ID canton year name firstname alliance* suballiance* party_num party elected votes_h votes_j)

//order ID name first party_num votemargin elected alliance_num suballiance_num party_num  number votes_h votes_j   suballiance*
//gsort party_num -votes_h

//replace votemargin = ceil(votemargin) if votemargin<0 // round up to next integer for not elected candidates
//replace votemargin = floor(votemargin) if votemargin>=0 // round down to next integer for elected candidates

if `P_size'==1 & `H_size'==1  replace votemargin=votes_j // replace votemargin for cantons with one party and one candidate

keep  ID year name firstname votes_h votesdiff votemargin partymarg elected party alliance  party_num  
if `yr'==1931 & "`can'"=="AG" {
save "$path\02_Processed_data\13_Running_variable\RV_analytical_all.dta", replace
}
else{
append using "$path\02_Processed_data\13_Running_variable\RV_analytical_all.dta"
save "$path\02_Processed_data\13_Running_variable\RV_analytical_all.dta", replace
}

} // end of cantonal loop
} // end of year loop



erase "$path\02_Processed_data\07_Cantons\Seats_perCtYear.dta"
erase "$path\02_Processed_data\13_Running_variable\example.dta"
