global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
global serverpath "\\Srw-hpc5\e\Projekt Nationalräte"

/*
** 1.) Kotakorpi et al. (2017) procedure for CH in Stata: 
*       m=seat*20, M=2000, check of seat allocation: yes

use "$path\02_Processed_data\13_Running_variable\election_bootstrap.dta", clear

capture erase "$path\02_Processed_data\13_Running_variable\candidate_seats_all_Stata.dta"
preserve
qui levelsof canton, local(cans) // 1. Loop over cantons
foreach can of local cans{
use "$path\02_Processed_data\13_Running_variable\01_kotakorpietal2017_ch_Stata\candidate_seats_`can'.dta",clear
capture append using "$path\02_Processed_data\13_Running_variable\candidate_seats_all_Stata.dta"
sort party_num votes_h_rank
save "$path\02_Processed_data\13_Running_variable\candidate_seats_all_Stata.dta", replace
}
restore

bysort party_num: egen votes_h_rank =rank(-votes_h), unique  // generate rank of candidate on party list
sort party_num votes_h_rank
merge 1:1 party_num votes_h_rank using "$path\02_Processed_data\13_Running_variable\candidate_seats_all_Stata.dta"
keep ID year party_num seat_all seat_all_rank p
rename p p_1
drop if ID==""
save "$path\02_Processed_data\13_Running_variable\candidate_seats_all_Stata.dta", replace
*/


** 2.) Kotakorpi et al. (2017) procedure for CH in R: 
*       m=10000, M=20000, check of seat allocation: no

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
qui levelsof canton, local(cans) // 1. Loop over cantons

capture erase "$path\02_Processed_data\13_Running_variable\candidate_seats_all_R.dta"
foreach can of local cans{
import delimited "$path\02_Processed_data\13_Running_variable\02_kotakorpietal2017_ch_R\kotakorpietal2017_ch`can'.csv",clear
capture append using "$path\02_Processed_data\13_Running_variable\candidate_seats_all_R.dta"
save "$path\02_Processed_data\13_Running_variable\candidate_seats_all_R.dta", replace
}

gen p_2=elected_sum/m
rename id ID
save "$path\02_Processed_data\13_Running_variable\candidate_seats_all_R.dta", replace


** 3.) Kotakorpi et al. (2017) procedure for CH in R: 
*       m=10000, M=20000, check of seat allocation: yes



** 4.) Merge Stata and R data with original data

local vars "p_1 p_2"

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
merge 1:1 ID year using "$path\02_Processed_data\13_Running_variable\candidate_seats_all_Stata.dta", gen(merge_Stata)
merge 1:1 ID year using "$path\02_Processed_data\13_Running_variable\candidate_seats_all_R.dta", gen(merge_R)

*(i) Check whether we achieve convergence

gsort party_num -votes -elected name firstname
qui bysort party_num: egen votes_h_rank=rank(-votes), unique // generate rank of candidate on party list

foreach var of local vars{
gsort party_num -`var' -votes -elected name firstname
qui bysort party_num: egen `var'_rank=rank(-`var'), unique // generate rank of candidate on party list
gen `var'_rank_diff=votes_h_rank-`var'_rank
qui bysort party_num: egen `var'_rank_diff_max=max(`var'_rank_diff)
}



tab canton year if p_2==.

sum p_2_rank_diff if p_2_rank_diff!=0

bysort elected: sum p_2

bysort party_num: egen help1=min(p_2) if elected==1
bysort party_num: egen help2=max(p_2) if elected==0
bysort party_num: egen p_2_min=max(help1) 
bysort party_num: egen p_2_max=max(help2) 
gen p_2_diff=p_2_min-p_2_max
tab p_2_diff
unique party_num if p_2_diff<0

* Result: For 12 parties and 297 cases, the smallest p_2 of the elected candidates
* on a party list is smaller than the highest p_2 of the non-elected candidates on the same list. 

gsort party_num -votes 
br ID name firstname partyname elected votes p_2 if p_2_diff<0

sort party_num votes_h
br ID canton year name firstname elected votes votes_h_rank p_2 p_2_rank p_2_rank_diff if p_2_rank_diff_max!=0 & p_2_rank_diff_max!=.
br ID canton year name firstname elected votes votes_h_rank p_1 p_1_rank if p_2_rank_diff_max!=0 & p_2_rank_diff_max!=.

*(ii) Rescale variables

qui bysort party_num: gen number_elected_help=_N if elected==1
qui bysort party_num: gen number_notelected_help=_N if elected==0
qui bysort party_num: egen number_elected=min(number_elected_help) 
qui bysort party_num: egen number_notelected=min(number_notelected_help)
drop *_help
replace number_elected=0 if number_elected==.
replace number_notelected=0 if number_notelected==.


foreach var of local vars{
qui bysort party_num: egen `var'_min_help=min(`var') if elected==1
qui bysort party_num elected: egen `var'_max_help=max(`var') if elected==0
qui bysort party_num: egen `var'_min=min(`var'_min_help) 
qui bysort party_num: egen `var'_max=min(`var'_max_help)
gen `var'_marginal=(`var'_min+`var'_max)/2
replace `var'_marginal=1 if number_elected==0
replace `var'_marginal=0 if number_notelected==0
gen `var'_transformed=`var'-`var'_marginal
}

qui bysort party_num: egen votes_min_help=min(votes) if elected==1
qui bysort party_num elected: egen votes_max_help=max(votes) if elected==0
qui bysort party_num: egen votes_min=min(votes_min_help) 
qui bysort party_num: egen votes_max=min(votes_max_help)
gen votes_marginal=(votes_min+votes_max)/2
gen votes_diff=votes-votes_marginal

gsort party_num -votes
br ID name firstname votes elected votes_diff // seems correct
bysort canton year: egen votes_sum=sum(votes)

gen votes_transformed=votes_diff/eligible_cant

sum votes_diff if number_elected==0
gsort party_num -votes
br party_num ID elected votes votes_diff votes_transformed if number_elected==0
br party_num ID elected votes votes_diff votes_transformed if party_num==2


*(iv) Analysis 

hist p_2_transformed
bysort elected: sum p_2_transformed

sum p_2_transformed if elected==1 & p_2_transformed>=0
sum p_2_transformed if elected==1 & p_2_transformed<0
sum p_2_transformed if elected==0 & p_2_transformed>=0
sum p_2_transformed if elected==0 & p_2_transformed<0

sum p_2_transformed if elected==1 & p_2>=0
sum p_2_transformed if elected==1 & p_2<0
sum p_2_transformed if elected==0 & p_2>=0
sum p_2_transformed if elected==0 & p_2<0
br if elected==1 & p_2_transformed<0

gsort party_num -votes
br ID year name	firstname elected votes p_2* if party_num==4230


sum votes_diff if elected==1 & votes_diff>=0
sum votes_diff if elected==1 & votes_diff<0
sum votes_diff if elected==0 & votes_diff>=0
sum votes_diff if elected==0 & votes_diff<0


sum votes_transformed if elected==1 & votes_transformed>=0
sum votes_transformed if elected==1 & votes_transformed<0
sum votes_transformed if elected==0 & votes_transformed>=0
sum votes_transformed if elected==0 & votes_transformed<0


cor votemargin votemargin_rel p_2 p_2_transformed 
cor votemargin votemargin_rel p_2 p_2_transformed votes_diff votes_transformed


tab p_2_rank_diff
