******************************************************
* D'Hondt running variable with Listenverbindung *
******************************************************

clear
set more off

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "\\Srw-hpc5\e\Projekt Nationalräte"

******************************************
*  Set parameters							 *
******************************************


global control_seat_allocation 0 // 1 for check of seat allocation, 0 for random votes
global check_convergence 0 // 1 if we check convergence, 0 otherwise
global m 10000 // number of trials (only applied if check_convergence=0)
global M 20000 // number of bootstrap repetitions  (only applied if check_convergence=0)
global canton AG

do "$path\03_Code\13_Running_variable\DHondt_Program.do" 

set seed 05071975

******************************************
* A) Simulation program					 *
******************************************

capture program drop election_bootstrap
program define election_bootstrap // 1. canton 2. year 3. seats
global n `3'

* (i) Prepare election data for a particular canton and a particular year

quietly{
use "$path\02_Processed_data\13_Running_variable\election_bootstrap.dta", clear
keep if canton=="`1'" & year==`2'
replace votes_j=votes_h if votes_j==. // replace list votes with individual votes for lists with a single candidate

* (ii) Introduce "extra candidate" to account for empty party votes

bysort party_num: egen votes_h_sum=sum(votes_h)
gen votes_j_diff=votes_j-votes_h_sum
replace votes_j=votes_h_sum if votes_j_diff<0 // should be only for FDP in Zurich 1939
replace votes_j_diff=0 if votes_j_diff<0 // should be only for FDP in Zurich 1939

preserve
collapse (first) votes_j_diff canton year party alliance suballiance votes_j alliance_num suballiance_num, by(party_num)
rename votes_j_diff votes_h
gen ID="Extra Candidate"
save "$path\02_Processed_data\13_Running_variable\extra_candidates_$canton.dta", replace
restore

append using "$path\02_Processed_data\13_Running_variable\extra_candidates_$canton.dta"

bysort party_num: egen votes_h_sum_test=sum(votes_h)
gen votes_j_diff_test=votes_j-votes_h_sum_test
sum votes_j_diff_test
if `r(mean)'!=0 {
di "Total candidate votes not equal to party votes."
exit 
}
drop votes_h_sum_test votes_j_diff_test

sum(votes_h)
gen double votes_h_share=votes_h/`r(sum)'
gen id=_n

drop votes_j_diff votes_h_sum

gen m=$m 
gen M=$M 

save "$path\02_Processed_data\13_Running_variable\example_all_$canton.dta", replace


* (iii) Draw random votes, run D'Hondt procedure and merge it to initial dataset

if $control_seat_allocation==0{

forvalues i=1(1)$M{   // comment out for check 
*if mod(`i',1)==0 noisily di "New vote iteration:" "`i'" 
use "$path\02_Processed_data\13_Running_variable\example_all_$canton.dta", clear

mata:st_view(p=., ., "votes_h_share")
mata:st_view(m=., ., "m")
mata:st_view(M=., ., "M")
mata:m=m[1,1]
mata:M=M[1,1]
mata:votes_random=rdiscrete(m,1,p)

clear
set obs $m
getmata (votes_random)=votes_random

bysort votes_random: gen nvotes=_N
collapse (first) nvotes ,by(votes_random)
rename votes_random id
save "$path\02_Processed_data\13_Running_variable\random_votes_$canton.dta", replace

use "$path\02_Processed_data\13_Running_variable\example_all_$canton.dta",clear
merge 1:1 id using "$path\02_Processed_data\13_Running_variable\random_votes_$canton.dta"
count if _merge==2
if `r(N)'>0 {
di "Problem with using file"
exit 
}

drop _merge votes_h votes_j votes_h_share 
rename nvotes votes_h 
replace votes_h=0 if votes_h==.
bysort party_num: egen votes_j=sum(votes_h)
drop if ID=="Extra Candidate"
DHondt
rename number votes_h_rank
rename seat seat_`i'
if `i'==1 save "$path\02_Processed_data\13_Running_variable\temp_$canton.dta", replace
else{
merge 1:1 party_num votes_h_rank using  "$path\02_Processed_data\13_Running_variable\temp_$canton.dta", nogen
save  "$path\02_Processed_data\13_Running_variable\temp_$canton.dta", replace
}
}
gen canton="`1'"
gen year= `2' 
}

* (iv) Run D'Hondt procedure to check seat allocation

if $control_seat_allocation==1{
use "$path\02_Processed_data\13_Running_variable\example_all_$canton.dta", clear
drop votes_j
bysort party_num: egen votes_j=sum(votes_h)
drop if ID=="Extra Candidate"
 
DHondt
rename number votes_h_rank
gen canton="`1'"
gen year= `2'
}
}
end


capture program drop delete_data
program delete_data
cap erase "$path\02_Processed_data\13_Running_variable\example_all_$canton.dta"
cap erase "$path\02_Processed_data\13_Running_variable\extra_candidates_$canton.dta"
cap erase "$path\02_Processed_data\13_Running_variable\temp_$canton.dta"
cap erase "$path\02_Processed_data\13_Running_variable\SeatsSuballiance_$canton.dta"
cap erase "$path\02_Processed_data\13_Running_variable\SeatsAlliance_$canton.dta"
cap erase "$path\02_Processed_data\13_Running_variable\example_input_$canton.dta"
cap erase "$path\02_Processed_data\13_Running_variable\random_votes_$canton.dta"
cap erase "$path\02_Processed_data\13_Running_variable\candidate_seats_temp_$canton.dta"
end


************************************************************************
* B) Read-in data     												   *
************************************************************************


use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear

drop if name==""
replace list=name if list==""
rename votes votes_h
rename pvotes votes_j
rename list party

replace alliance=party if alliance=="." | alliance==""
replace suballiance=party if suballianc=="." | suballiance==""

sort alliance suballiance party
egen party_num=group(canton year party) // numeric version of party
egen alliance_num=group(canton year alliance) // numeric version of alliance
egen suballiance_num=group(canton year alliance suballiance) // numeric version of suballiance

gsort canton year alliance party -votes_h
bysort canton year party: gen id=_n // indicator for alliance

keep ID canton year name firstname alliance* suballiance* party party_num elected votes_h votes_j

save "$path\02_Processed_data\13_Running_variable\election_bootstrap.dta", replace


quietly{
clear
import excel "$path\01_Raw_data\07_Cantons\Verteilung der Nationalratssitze auf die Kantone.xlsx",  first sheet("out")
reshape long Seats_ct, i(canton) j(year)
keep if year<2019 
save "$path\02_Processed_data\07_Cantons\\Seats_perCtYear.dta", replace

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
merge m:1 canton year using "$path\02_Processed_data\07_Cantons\Seats_perCtYear.dta", gen(merge_seats)
keep if votes!=. // exclude all elections with "stillen Wahlen"
collapse (mean) Seats_ct , by(canton year)
save "$path\02_Processed_data\07_Cantons\Seats_perCtYear", replace
}


************************************************************************
* C) Check seat allocation over time   							   	   *
************************************************************************
/*
global control_seat_allocation 1 //1 for check of seat allocation, 0 for random votes

use "$path\02_Processed_data\13_Running_variable\election_bootstrap.dta", clear

qui levelsof canton, local(cans) // 1. Loop over cantons
local index 1
foreach can of local cans{
qui use "$path\02_Processed_data\07_Cantons\Seats_perCtYear", clear
qui keep if canton=="`can'" 
qui levelsof year, local(yrs) // 2. Loop over years
foreach yr of local yrs{
qui use "$path\02_Processed_data\07_Cantons\Seats_perCtYear", clear
qui sum Seats_ct if canton=="`can'"  & year==`yr'
di "`can'" " " "`yr'" " " "`r(mean)'"
election_bootstrap "`can'" `yr'	`r(mean)'
if `index'==1 qui save "$path\02_Processed_data\13_Running_variable\candidate_seats.dta", replace
else{
qui append using  "$path\02_Processed_data\13_Running_variable\candidate_seats.dta"
qui save  "$path\02_Processed_data\13_Running_variable\candidate_seats.dta", replace
}
local index=`index'+1
}
}


use "$path\02_Processed_data\13_Running_variable\election_bootstrap.dta", clear
bysort party_num: egen votes_h_rank=rank(-votes_h), unique  // generate rank of candidate on party list
merge 1:1 party_num votes_h_rank using "$path\02_Processed_data\13_Running_variable\candidate_seats.dta"
tab elected seats

* Result: looks fine except for TI 2011 in which year both candidates had the same #votes

*/

************************************************************************
* D) Run randomization of votes_h								   	   *
************************************************************************


use "$path\02_Processed_data\13_Running_variable\election_bootstrap.dta", clear

keep if canton=="$canton" // comment out for general loop

qui levelsof canton, local(cans) // 1. Loop over cantons
local index 1
foreach can of local cans{
qui use "$path\02_Processed_data\07_Cantons\Seats_perCtYear", clear
qui keep if canton=="`can'" 
qui levelsof year, local(yrs) // 2. Loop over years
foreach yr of local yrs{
qui use "$path\02_Processed_data\07_Cantons\Seats_perCtYear", clear
qui sum Seats_ct if canton=="`can'"  & year==`yr'
local seats_cant=`r(mean)'

*** (a) check convergence (see Kotakorpi et al. 2017, footnote 3)
if $check_convergence==1{
qui gen help=`seats_cant'*20 
qui sum help
global m `r(mean)' // number of trials
global M 2000 // number of bootstrap repetitions
local M_total=2000 // counter for number of M 
di "`can'" " " "`yr'" " " "`seats_cant'"
election_bootstrap "`can'" `yr'	`seats_cant'
qui egen seat_all=rowtotal(seat*)
gsort party_num -seat_all -votes_h elected firstname name // this is only to break ties in seat_all_rank the same way as we break ties in votes_h_rank
qui bysort party_num: egen seat_all_rank=rank(-seat_all), unique // generate rank of candidate on party list
keep party_num votes_h_rank  seat_all seat_all_rank
qui save  "$path\02_Processed_data\13_Running_variable\temp_all_$canton.dta", replace

preserve
use "$path\02_Processed_data\13_Running_variable\example_all_$canton.dta", clear
drop if ID=="Extra Candidate"
gsort party_num -votes_h elected firstname name // this is only to break ties in seat_all_rank the same way as we break ties in votes_h_rank
qui bysort party_num: egen votes_h_rank=rank(-votes_h), unique // generate rank of candidate on party list
qui merge 1:1 party_num votes_h_rank using "$path\02_Processed_data\13_Running_variable\temp_all_$canton.dta", nogen
qui gen rank_diff=abs(seat_all_rank-votes_h_rank)
qui sum rank_diff
restore

while `r(min)'> 0{ // do until convergence
local M_total=`M_total'+10000
use  "$path\02_Processed_data\13_Running_variable\temp_$canton.dta", replace
global M 10000 // number of bootstrap repetitions
election_bootstrap "`can'" `yr' `seats_cant'
qui egen seat_all_new=rowtotal(seat*)
merge 1:1 party_num votes_h_rank using "$path\02_Processed_data\13_Running_variable\temp_all_$canton.dta", nogen
qui replace seat_all=seat_all+seat_all_new
keep party_num votes_h_rank  seat_all
gsort party_num -seat_all -votes_h elected firstname name // this is only to break ties in seat_all_rank the same way as we break ties in votes_h_rank
qui bysort party_num: egen seat_all_rank=rank(-seat_all), unique // generate rank of candidate on party list
qui save  "$path\02_Processed_data\13_Running_variable\temp_all_$canton.dta", replace

preserve
use "$path\02_Processed_data\13_Running_variable\example_all_$canton.dta", clear
drop if ID=="Extra Candidate"
gsort party_num -votes_h -elected
qui bysort party_num: egen votes_h_rank=rank(-votes_h), unique // generate rank of candidate on party list
qui merge 1:1 party_num votes_h_rank using "$path\02_Processed_data\13_Running_variable\temp_all_$canton.dta", nogen
qui gen rank_diff=abs(seat_all_rank-votes_h_rank)
qui sum rank_diff
restore
}
gen p=seat_all/`M_total'
}
else{
}

*** (b) do not check convergence (fixed M)

if $check_convergence==0{
di "`can'" " " "`yr'" " " "`seats_cant'"
election_bootstrap "`can'" `yr'	`seats_cant'
qui egen seat_all=rowtotal(seat*)
gen p=seat_all/$M
keep party_num votes_h_rank p seat_all
}

*keep party_num votes_h_rank p

if `index'==1 qui save "$path\02_Processed_data\13_Running_variable\candidate_seats_$canton.dta", replace
else{
qui append using  "$path\02_Processed_data\13_Running_variable\candidate_seats_$canton.dta"
qui save  "$path\02_Processed_data\13_Running_variable\candidate_seats_$canton.dta", replace
}
local index=`index'+1
}
}

delete_data // delete all temporary files

/*

************************************************************************
* E) Append all cantonal files that were run separately			   	   *
************************************************************************

use "$path\02_Processed_data\13_Running_variable\election_bootstrap.dta", clear

capture erase "$path\02_Processed_data\13_Running_variable\candidate_seats_all.dta"
preserve
qui levelsof canton, local(cans) // 1. Loop over cantons
foreach can of local cans{
use "$path\02_Processed_data\13_Running_variable\candidate_seats_`can'.dta",clear
capture append using "$path\02_Processed_data\13_Running_variable\candidate_seats_all.dta"
sort party_num votes_h_rank
save "$path\02_Processed_data\13_Running_variable\candidate_seats_all.dta", replace
}
restore

bysort party_num: egen votes_h_rank =rank(-votes_h), unique  // generate rank of candidate on party list
sort party_num votes_h_rank
merge 1:1 party_num votes_h_rank using "$path\02_Processed_data\13_Running_variable\candidate_seats_all.dta"
keep ID year party_num seat_all seat_all_rank p
drop if ID==""
save "$path\02_Processed_data\13_Running_variable\candidate_seats_all.dta", replace


use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
merge 1:1 ID year using "$path\02_Processed_data\13_Running_variable\candidate_seats_all.dta"

*(i) Check whether we achieve convergence

gsort party_num -votes -elected name firstname
qui bysort party_num: egen votes_h_rank=rank(-votes), unique // generate rank of candidate on party list
gsort party_num -p -votes -elected name firstname
qui bysort party_num: egen p_rank=rank(-p), unique // generate rank of candidate on party list

tab canton year if p==.

gen rank_diff=votes_h_rank-p_rank
sum rank_diff
sum rank_diff if rank_diff!=0

qui bysort party_num: egen rank_diff_max=max(rank_diff)
sort party_num votes_h
br ID canton year name firstname elected votes votes_h_rank p p_rank if rank_diff_max!=0
br ID canton year name firstname elected votes votes_h_rank p p_rank if rank_diff!=0

*(ii) Rescale variable

qui bysort party_num: egen p_min_help=min(p) if elected==1
qui bysort party_num elected: egen p_max_help=max(p) if elected==0
qui bysort party_num: egen p_min=min(p_min_help) 
qui bysort party_num: egen p_max=min(p_max_help)
qui bysort party_num: gen number_elected_help=_N if elected==1
qui bysort party_num: gen number_notelected_help=_N if elected==0
qui bysort party_num: egen number_elected=min(number_elected_help) 
qui bysort party_num: egen number_notelected=min(number_notelected_help)
drop *_help
replace number_elected=0 if number_elected==.
replace number_notelected=0 if number_notelected==.

gen p_marginal=(p_min+p_max)/2
replace p_marginal=1 if number_elected==0
replace p_marginal=0 if number_notelected==0

gen p_transformed=p-p_marginal
hist p_transformed
bysort elected: sum p_transformed

br ID party_num canton year name firstname elected votes p p_marginal number_* if p_transformed<0.2 & elected==1
br ID party_num canton year name firstname elected votes p p_marginal number_* if p_transformed>0.2 & elected==0
br ID party_num canton year name firstname elected votes p p_marginal number_* if party_num==1989

cor votemargin p
cor votemargin_rel p

scatter votemargin p

sum elected if p==0
sum elected if p==1
sum elected if p!=1 & p!=0
tab elected if inrange(p,0.0000001,0.5)
tab elected if inrange(p,0.50000001,0.99999999999999)