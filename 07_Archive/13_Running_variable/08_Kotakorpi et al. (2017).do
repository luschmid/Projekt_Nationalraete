******************************************************
* D'Hondt running variable with Listenverbindung *
******************************************************

clear
set more off

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"

******************************************
*  Set parameters							 *
******************************************

timer on 1
global M 13 // number of bootstrap repetitions
global m 10000 // number of trials

do "$path\03_Code\13_Running_variable\DHondt_Program.do" 


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

*bysort party_num: egen cand=rank(-votes_h), unique // generate rank of candidate on party list


* (ii) Introduce "extra candidate" to account for empty party votes

bysort party_num: egen votes_h_sum=sum(votes_h)
gen votes_j_diff=votes_j-votes_h_sum
replace votes_j=votes_h_sum if votes_j_diff<0 // should be only for FDP in Zurich 1939
replace votes_j_diff=0 if votes_j_diff<0 // should be only for FDP in Zurich 1939

preserve
collapse (first) votes_j_diff canton year party alliance suballiance votes_j alliance_num suballiance_num, by(party_num)
rename votes_j_diff votes_h
gen ID="Extra Candidate"
save "$path\02_Processed_data\13_Running_variable\extra_candidates.dta", replace
restore

append using "$path\02_Processed_data\13_Running_variable\extra_candidates.dta"

bysort party_num: egen votes_h_sum_test=sum(votes_h)
gen votes_j_diff_test=votes_j-votes_h_sum_test
sum votes_j_diff_test
if `r(mean)'!=0 {
di "Total candidate votes not equal to party votes."
exit 
}

* (iii) Draw random votes and merge them to initial dataset

gen id=_n

sum(votes_h)
gen double votes_h_share=votes_h/`r(sum)'

drop votes_j_diff votes_h_sum
save "$path\02_Processed_data\13_Running_variable\example_all.dta", replace

gen m=$m 
gen M=$M 

mata:st_view(p=., ., "votes_h_share")
mata:st_view(m=., ., "m")
mata:st_view(M=., ., "M")
mata:m=m[1,1]
mata:M=M[1,1]
mata:votes_random=rdiscrete(m,M,p)

clear
set obs $m
getmata (votes_random*)=votes_random

forvalues i=1(1)$M{
preserve
bysort votes_random`i': gen nvotes`i'=_N
collapse (first) nvotes`i' ,by(votes_random`i')
rename votes_random id
save "$path\02_Processed_data\13_Running_variable\random_votes.dta", replace
use "$path\02_Processed_data\13_Running_variable\example_all.dta",clear
merge 1:1 id using "$path\02_Processed_data\13_Running_variable\random_votes.dta"
count if _merge==2
if `r(N)'>0 {
di "Problem with using file"
exit 
}
drop _merge
save "$path\02_Processed_data\13_Running_variable\example_all.dta", replace
restore
}

* (iv) Generate random party votes by aggregating random candidate votes

use "$path\02_Processed_data\13_Running_variable\example_all.dta", clear
*drop votes_j votes_h
drop votes_h_share
forvalues i=1(1)$M{
rename nvotes`i' votes_h_`i'
replace votes_h_`i'=0 if votes_h_`i'==.
bysort party_num: egen votes_j_`i'=sum(votes_h_`i')
}

drop if ID=="Extra Candidate"

save "$path\02_Processed_data\13_Running_variable\temp.dta", replace


* (v) Run D'Hondt procedure and merge it to initial dataset

forvalues i=1(1)$M{ // comment out for check 
use "$path\02_Processed_data\13_Running_variable\temp.dta", clear
replace votes_h=votes_h_`i' // comment out for check 
replace votes_j=votes_j_`i' // comment out for check 

qui levelsof party_num, local(jlevels)
foreach j of local jlevels {	
use "$path\02_Processed_data\13_Running_variable\temp.dta", clear
qui sum alliance_num if party_num ==`j'
local l=r(max)
di "`l'"
qui sum suballiance_num if party_num ==`j'
local s=r(max)
di "`s'"
DHondt `l' `s' `j'
keep if party_num==`j'
rename number votes_h_rank
gen canton="`1'"
gen year= `2'
keep canton year party_num votes_h_rank seat
rename seat seat_`i' // comment out for check 
local index 1
if `index'==1 save "$path\02_Processed_data\13_Running_variable\partyseats.dta", replace
else{
append using  "$path\02_Processed_data\13_Running_variable\partyseats.dta"
save  "$path\02_Processed_data\13_Running_variable\partyseats.dta", replace
local index=`index'+1
}
use "$path\02_Processed_data\13_Running_variable\election_bootstrap.dta", clear
merge 1:1 canton year party_num votes_h_rank using "$path\02_Processed_data\13_Running_variable\partyseats.dta", update
drop if _merge==2
drop _merge
save "$path\02_Processed_data\13_Running_variable\election_bootstrap.dta", replace
}
} // comment out for check 
erase "$path\02_Processed_data\13_Running_variable\temp.dta"
erase "$path\02_Processed_data\13_Running_variable\partyseats.dta"
erase "$path\02_Processed_data\13_Running_variable\example_all.dta"
erase "$path\02_Processed_data\13_Running_variable\random_votes.dta"
erase "$path\02_Processed_data\13_Running_variable\extra_candidates.dta"
}
end


************************************************************************
* B) Read-in data and runimulation program	over all cantons and years *
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

bysort canton year party_num: egen votes_h_rank=rank(-votes_h), unique


gsort canton year alliance party -votes_h
bysort canton year party: gen id=_n // indicator for alliance

keep ID canton year name firstname alliance* suballiance* party party_num elected votes_h votes_j votes_h_rank

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


qui levelsof canton, local(cans) // 1. Loop over cantons
foreach can of local cans{
qui use "$path\02_Processed_data\07_Cantons\Seats_perCtYear", clear
qui keep if canton=="`can'" 
qui levelsof year, local(yrs) // 2. Loop over years
foreach yr of local yrs{
qui use "$path\02_Processed_data\07_Cantons\Seats_perCtYear", clear
qui sum Seats_ct if canton=="`can'"  & year==`yr'
di "`can'" " " "`yr'" " " "`r(mean)'"
election_bootstrap "`can'" `yr'	`r(mean)'
}
}

/*
election_bootstrap AG	1931	12
election_bootstrap AG	1935	12
election_bootstrap AG	1939	12
election_bootstrap AG	1943	12
election_bootstrap AG	1947	12
election_bootstrap AG	1951	13
election_bootstrap AG	1955	13
election_bootstrap AG	1959	13
election_bootstrap AG	1963	13
election_bootstrap AG	1967	13
election_bootstrap AG	1971	14
election_bootstrap AG	1975	14
election_bootstrap AG	1979	14
election_bootstrap AG	1983	14
election_bootstrap AG	1987	14
election_bootstrap AG	1991	14
election_bootstrap AG	1995	15
election_bootstrap AG	1999	15
election_bootstrap AG	2003	15
election_bootstrap AG	2007	15
election_bootstrap AG	2011	15
election_bootstrap AG	2015	16
election_bootstrap AI	1931	1
election_bootstrap AI	1935	1
election_bootstrap AI	1939	1
election_bootstrap AI	1943	1
election_bootstrap AI	1947	1
election_bootstrap AI	1951	1
election_bootstrap AI	1955	1
election_bootstrap AI	1959	1
election_bootstrap AI	1963	1
election_bootstrap AI	1967	1
election_bootstrap AI	1971	1
election_bootstrap AI	1975	1
election_bootstrap AI	1979	1
election_bootstrap AI	1983	1
election_bootstrap AI	1987	1
election_bootstrap AI	1991	1
election_bootstrap AI	1995	1
election_bootstrap AI	1999	1
election_bootstrap AI	2003	1
election_bootstrap AI	2007	1
election_bootstrap AI	2011	1
election_bootstrap AI	2015	1
election_bootstrap AR	1931	2
election_bootstrap AR	1935	2
election_bootstrap AR	1951	2
election_bootstrap AR	1955	2
election_bootstrap AR	1971	2
election_bootstrap AR	1975	2
election_bootstrap AR	1983	2
election_bootstrap AR	1991	2
election_bootstrap AR	1995	2
election_bootstrap AR	1999	2
election_bootstrap AR	2003	1
election_bootstrap AR	2007	1
election_bootstrap AR	2011	1
election_bootstrap AR	2015	1
election_bootstrap BE	1931	31
election_bootstrap BE	1935	31
election_bootstrap BE	1939	31
election_bootstrap BE	1943	33
election_bootstrap BE	1947	33
election_bootstrap BE	1951	33
election_bootstrap BE	1955	33
election_bootstrap BE	1959	33
election_bootstrap BE	1963	33
election_bootstrap BE	1967	33
election_bootstrap BE	1971	31
election_bootstrap BE	1975	31
election_bootstrap BE	1979	29
election_bootstrap BE	1983	29
election_bootstrap BE	1987	29
election_bootstrap BE	1991	29
election_bootstrap BE	1995	27
election_bootstrap BE	1999	27
election_bootstrap BE	2003	26
election_bootstrap BE	2007	26
election_bootstrap BE	2011	26
election_bootstrap BE	2015	25
election_bootstrap BL	1931	4
election_bootstrap BL	1935	4
election_bootstrap BL	1939	4
election_bootstrap BL	1943	4
election_bootstrap BL	1947	4
election_bootstrap BL	1951	4
election_bootstrap BL	1955	4
election_bootstrap BL	1959	4
election_bootstrap BL	1963	5
election_bootstrap BL	1967	5
election_bootstrap BL	1971	7
election_bootstrap BL	1975	7
election_bootstrap BL	1979	7
election_bootstrap BL	1983	7
election_bootstrap BL	1987	7
election_bootstrap BL	1991	7
election_bootstrap BL	1995	7
election_bootstrap BL	1999	7
election_bootstrap BL	2003	7
election_bootstrap BL	2007	7
election_bootstrap BL	2011	7
election_bootstrap BL	2015	7
election_bootstrap BS	1931	7
election_bootstrap BS	1935	7
election_bootstrap BS	1939	7
election_bootstrap BS	1943	8
election_bootstrap BS	1947	8
election_bootstrap BS	1951	8
election_bootstrap BS	1955	8
election_bootstrap BS	1959	8
election_bootstrap BS	1963	8
election_bootstrap BS	1967	8
election_bootstrap BS	1971	7
election_bootstrap BS	1975	7
election_bootstrap BS	1979	7
election_bootstrap BS	1983	6
election_bootstrap BS	1987	6
election_bootstrap BS	1991	6
election_bootstrap BS	1995	6
election_bootstrap BS	1999	6
election_bootstrap BS	2003	5
election_bootstrap BS	2007	5
election_bootstrap BS	2011	5
election_bootstrap BS	2015	5
election_bootstrap FR	1931	7
election_bootstrap FR	1935	7
election_bootstrap FR	1939	7
election_bootstrap FR	1943	7
election_bootstrap FR	1947	7
election_bootstrap FR	1951	7
election_bootstrap FR	1955	7
election_bootstrap FR	1959	7
election_bootstrap FR	1963	6
election_bootstrap FR	1967	6
election_bootstrap FR	1971	6
election_bootstrap FR	1975	6
election_bootstrap FR	1979	6
election_bootstrap FR	1983	6
election_bootstrap FR	1987	6
election_bootstrap FR	1991	6
election_bootstrap FR	1995	6
election_bootstrap FR	1999	6
election_bootstrap FR	2003	7
election_bootstrap FR	2007	7
election_bootstrap FR	2011	7
election_bootstrap FR	2015	7
election_bootstrap GE	1931	8
election_bootstrap GE	1935	8
election_bootstrap GE	1939	8
election_bootstrap GE	1943	8
election_bootstrap GE	1947	8
election_bootstrap GE	1951	8
election_bootstrap GE	1955	8
election_bootstrap GE	1959	8
election_bootstrap GE	1963	10
election_bootstrap GE	1967	10
election_bootstrap GE	1971	11
election_bootstrap GE	1975	11
election_bootstrap GE	1979	11
election_bootstrap GE	1983	11
election_bootstrap GE	1987	11
election_bootstrap GE	1991	11
election_bootstrap GE	1995	11
election_bootstrap GE	1999	11
election_bootstrap GE	2003	11
election_bootstrap GE	2007	11
election_bootstrap GE	2011	11
election_bootstrap GE	2015	11
election_bootstrap GL	1931	2
election_bootstrap GL	1935	2
election_bootstrap GL	1939	2
election_bootstrap GL	1943	2
election_bootstrap GL	1947	2
election_bootstrap GL	1959	2
election_bootstrap GL	1971	1
election_bootstrap GL	1975	1
election_bootstrap GL	1979	1
election_bootstrap GL	1983	1
election_bootstrap GL	1987	1
election_bootstrap GL	1991	1
election_bootstrap GL	1995	1
election_bootstrap GL	1999	1
election_bootstrap GL	2003	1
election_bootstrap GL	2007	1
election_bootstrap GL	2011	1
election_bootstrap GL	2015	1
election_bootstrap GR	1931	6
election_bootstrap GR	1935	6
election_bootstrap GR	1939	6
election_bootstrap GR	1943	6
election_bootstrap GR	1947	6
election_bootstrap GR	1951	6
election_bootstrap GR	1955	6
election_bootstrap GR	1959	6
election_bootstrap GR	1963	5
election_bootstrap GR	1967	5
election_bootstrap GR	1971	5
election_bootstrap GR	1975	5
election_bootstrap GR	1979	5
election_bootstrap GR	1983	5
election_bootstrap GR	1987	5
election_bootstrap GR	1991	5
election_bootstrap GR	1995	5
election_bootstrap GR	1999	5
election_bootstrap GR	2003	5
election_bootstrap GR	2007	5
election_bootstrap GR	2011	5
election_bootstrap GR	2015	5
election_bootstrap JU	1979	2
election_bootstrap JU	1983	2
election_bootstrap JU	1987	2
election_bootstrap JU	1991	2
election_bootstrap JU	1995	2
election_bootstrap JU	1999	2
election_bootstrap JU	2003	2
election_bootstrap JU	2007	2
election_bootstrap JU	2011	2
election_bootstrap JU	2015	2
election_bootstrap LU	1931	9
election_bootstrap LU	1935	9
election_bootstrap LU	1943	9
election_bootstrap LU	1947	9
election_bootstrap LU	1951	9
election_bootstrap LU	1955	9
election_bootstrap LU	1959	9
election_bootstrap LU	1963	9
election_bootstrap LU	1967	9
election_bootstrap LU	1971	9
election_bootstrap LU	1975	9
election_bootstrap LU	1979	9
election_bootstrap LU	1983	9
election_bootstrap LU	1987	9
election_bootstrap LU	1991	9
election_bootstrap LU	1995	10
election_bootstrap LU	1999	10
election_bootstrap LU	2003	10
election_bootstrap LU	2007	10
election_bootstrap LU	2011	10
election_bootstrap LU	2015	10
election_bootstrap NE	1931	6
election_bootstrap NE	1935	6
election_bootstrap NE	1943	5
election_bootstrap NE	1947	5
election_bootstrap NE	1951	5
election_bootstrap NE	1955	5
election_bootstrap NE	1959	5
election_bootstrap NE	1963	5
election_bootstrap NE	1967	5
election_bootstrap NE	1971	5
election_bootstrap NE	1975	5
election_bootstrap NE	1979	5
election_bootstrap NE	1983	5
election_bootstrap NE	1987	5
election_bootstrap NE	1991	5
election_bootstrap NE	1995	5
election_bootstrap NE	1999	5
election_bootstrap NE	2003	5
election_bootstrap NE	2007	5
election_bootstrap NE	2011	5
election_bootstrap NE	2015	4
election_bootstrap NW	1931	1
election_bootstrap NW	1935	1
election_bootstrap NW	1939	1
election_bootstrap NW	1943	1
election_bootstrap NW	1947	1
election_bootstrap NW	1951	1
election_bootstrap NW	1955	1
election_bootstrap NW	1959	1
election_bootstrap NW	1963	1
election_bootstrap NW	1967	1
election_bootstrap NW	1971	1
election_bootstrap NW	1975	1
election_bootstrap NW	1979	1
election_bootstrap NW	1983	1
election_bootstrap NW	1987	1
election_bootstrap NW	1991	1
election_bootstrap NW	1995	1
election_bootstrap NW	1999	1
election_bootstrap NW	2003	1
election_bootstrap NW	2011	1
election_bootstrap NW	2015	1
election_bootstrap OW	1931	1
election_bootstrap OW	1935	1
election_bootstrap OW	1939	1
election_bootstrap OW	1943	1
election_bootstrap OW	1947	1
election_bootstrap OW	1951	1
election_bootstrap OW	1955	1
election_bootstrap OW	1959	1
election_bootstrap OW	1963	1
election_bootstrap OW	1967	1
election_bootstrap OW	1971	1
election_bootstrap OW	1975	1
election_bootstrap OW	1979	1
election_bootstrap OW	1983	1
election_bootstrap OW	1987	1
election_bootstrap OW	1991	1
election_bootstrap OW	1995	1
election_bootstrap OW	2003	1
election_bootstrap OW	2007	1
election_bootstrap OW	2011	1
election_bootstrap OW	2015	1
election_bootstrap SG	1931	13
election_bootstrap SG	1935	13
election_bootstrap SG	1939	13
election_bootstrap SG	1943	13
election_bootstrap SG	1947	13
election_bootstrap SG	1951	13
election_bootstrap SG	1955	13
election_bootstrap SG	1959	13
election_bootstrap SG	1963	13
election_bootstrap SG	1967	13
election_bootstrap SG	1971	12
election_bootstrap SG	1975	12
election_bootstrap SG	1979	12
election_bootstrap SG	1983	12
election_bootstrap SG	1987	12
election_bootstrap SG	1991	12
election_bootstrap SG	1995	12
election_bootstrap SG	1999	12
election_bootstrap SG	2003	12
election_bootstrap SG	2007	12
election_bootstrap SG	2011	12
election_bootstrap SG	2015	12
election_bootstrap SH	1931	2
election_bootstrap SH	1935	2
election_bootstrap SH	1939	2
election_bootstrap SH	1943	2
election_bootstrap SH	1947	2
election_bootstrap SH	1959	2
election_bootstrap SH	1963	2
election_bootstrap SH	1967	2
election_bootstrap SH	1971	2
election_bootstrap SH	1975	2
election_bootstrap SH	1979	2
election_bootstrap SH	1983	2
election_bootstrap SH	1987	2
election_bootstrap SH	1991	2
election_bootstrap SH	1995	2
election_bootstrap SH	1999	2
election_bootstrap SH	2003	2
election_bootstrap SH	2007	2
election_bootstrap SH	2011	2
election_bootstrap SH	2015	2
election_bootstrap SO	1931	7
election_bootstrap SO	1935	7
election_bootstrap SO	1943	7
election_bootstrap SO	1947	7
election_bootstrap SO	1951	7
election_bootstrap SO	1955	7
election_bootstrap SO	1959	7
election_bootstrap SO	1963	7
election_bootstrap SO	1967	7
election_bootstrap SO	1971	7
election_bootstrap SO	1975	7
election_bootstrap SO	1979	7
election_bootstrap SO	1983	7
election_bootstrap SO	1987	7
election_bootstrap SO	1991	7
election_bootstrap SO	1995	7
election_bootstrap SO	1999	7
election_bootstrap SO	2003	7
election_bootstrap SO	2007	7
election_bootstrap SO	2011	7
election_bootstrap SO	2015	6
election_bootstrap SZ	1931	3
election_bootstrap SZ	1935	3
election_bootstrap SZ	1943	3
election_bootstrap SZ	1947	3
election_bootstrap SZ	1951	3
election_bootstrap SZ	1955	3
election_bootstrap SZ	1959	3
election_bootstrap SZ	1963	3
election_bootstrap SZ	1971	3
election_bootstrap SZ	1975	3
election_bootstrap SZ	1979	3
election_bootstrap SZ	1983	3
election_bootstrap SZ	1987	3
election_bootstrap SZ	1991	3
election_bootstrap SZ	1995	3
election_bootstrap SZ	1999	3
election_bootstrap SZ	2003	4
election_bootstrap SZ	2007	4
election_bootstrap SZ	2011	4
election_bootstrap SZ	2015	4
election_bootstrap TG	1931	6
election_bootstrap TG	1935	6
election_bootstrap TG	1939	6
election_bootstrap TG	1943	6
election_bootstrap TG	1947	6
election_bootstrap TG	1951	6
election_bootstrap TG	1955	6
election_bootstrap TG	1959	6
election_bootstrap TG	1963	6
election_bootstrap TG	1967	6
election_bootstrap TG	1971	6
election_bootstrap TG	1975	6
election_bootstrap TG	1979	6
election_bootstrap TG	1983	6
election_bootstrap TG	1987	6
election_bootstrap TG	1991	6
election_bootstrap TG	1995	6
election_bootstrap TG	1999	6
election_bootstrap TG	2003	6
election_bootstrap TG	2007	6
election_bootstrap TG	2011	6
election_bootstrap TG	2015	6
election_bootstrap TI	1931	7
election_bootstrap TI	1935	7
election_bootstrap TI	1943	7
election_bootstrap TI	1947	7
election_bootstrap TI	1951	7
election_bootstrap TI	1955	7
election_bootstrap TI	1959	7
election_bootstrap TI	1963	7
election_bootstrap TI	1967	7
election_bootstrap TI	1971	8
election_bootstrap TI	1975	8
election_bootstrap TI	1979	8
election_bootstrap TI	1983	8
election_bootstrap TI	1987	8
election_bootstrap TI	1991	8
election_bootstrap TI	1995	8
election_bootstrap TI	1999	8
election_bootstrap TI	2003	8
election_bootstrap TI	2007	8
election_bootstrap TI	2011	8
election_bootstrap TI	2015	8
election_bootstrap UR	1931	1
election_bootstrap UR	1935	1
election_bootstrap UR	1939	1
election_bootstrap UR	1943	1
election_bootstrap UR	1947	1
election_bootstrap UR	1951	1
election_bootstrap UR	1955	1
election_bootstrap UR	1959	1
election_bootstrap UR	1963	1
election_bootstrap UR	1967	1
election_bootstrap UR	1971	1
election_bootstrap UR	1975	1
election_bootstrap UR	1979	1
election_bootstrap UR	1983	1
election_bootstrap UR	1987	1
election_bootstrap UR	1991	1
election_bootstrap UR	1995	1
election_bootstrap UR	1999	1
election_bootstrap UR	2003	1
election_bootstrap UR	2007	1
election_bootstrap UR	2011	1
election_bootstrap UR	2015	1
election_bootstrap VD	1931	15
election_bootstrap VD	1935	15
election_bootstrap VD	1943	16
election_bootstrap VD	1947	16
election_bootstrap VD	1951	16
election_bootstrap VD	1955	16
election_bootstrap VD	1959	16
election_bootstrap VD	1963	16
election_bootstrap VD	1967	16
election_bootstrap VD	1971	16
election_bootstrap VD	1975	16
election_bootstrap VD	1979	16
election_bootstrap VD	1983	17
election_bootstrap VD	1987	17
election_bootstrap VD	1991	17
election_bootstrap VD	1995	17
election_bootstrap VD	1999	17
election_bootstrap VD	2003	18
election_bootstrap VD	2007	18
election_bootstrap VD	2011	18
election_bootstrap VD	2015	18
election_bootstrap VS	1931	6
election_bootstrap VS	1935	6
election_bootstrap VS	1943	7
election_bootstrap VS	1947	7
election_bootstrap VS	1951	7
election_bootstrap VS	1955	7
election_bootstrap VS	1959	7
election_bootstrap VS	1963	7
election_bootstrap VS	1967	7
election_bootstrap VS	1971	7
election_bootstrap VS	1975	7
election_bootstrap VS	1979	7
election_bootstrap VS	1983	7
election_bootstrap VS	1987	7
election_bootstrap VS	1991	7
election_bootstrap VS	1995	7
election_bootstrap VS	1999	7
election_bootstrap VS	2003	7
election_bootstrap VS	2007	7
election_bootstrap VS	2011	7
election_bootstrap VS	2015	8
election_bootstrap ZG	1931	2
election_bootstrap ZG	1935	2
election_bootstrap ZG	1943	2
election_bootstrap ZG	1947	2
election_bootstrap ZG	1951	2
election_bootstrap ZG	1955	2
election_bootstrap ZG	1959	2
election_bootstrap ZG	1967	2
election_bootstrap ZG	1975	2
election_bootstrap ZG	1979	2
election_bootstrap ZG	1983	2
election_bootstrap ZG	1987	2
election_bootstrap ZG	1991	2
election_bootstrap ZG	1995	3
election_bootstrap ZG	1999	3
election_bootstrap ZG	2003	3
election_bootstrap ZG	2007	3
election_bootstrap ZG	2011	3
election_bootstrap ZG	2015	3
election_bootstrap ZH	1931	28
election_bootstrap ZH	1935	28
election_bootstrap ZH	1939	28
election_bootstrap ZH	1943	31
election_bootstrap ZH	1947	31
election_bootstrap ZH	1951	32
election_bootstrap ZH	1955	32
election_bootstrap ZH	1959	32
election_bootstrap ZH	1963	35
election_bootstrap ZH	1967	35
election_bootstrap ZH	1971	35
election_bootstrap ZH	1975	35
election_bootstrap ZH	1979	35
election_bootstrap ZH	1983	35
election_bootstrap ZH	1987	35
election_bootstrap ZH	1991	35
election_bootstrap ZH	1995	34
election_bootstrap ZH	1999	34
election_bootstrap ZH	2003	34
election_bootstrap ZH	2007	34
election_bootstrap ZH	2011	34
election_bootstrap ZH	2015	35

