global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
global serverpath "\\Srw-hpc5\e\Projekt Nationalräte"

local cans "AG BE SG VD ZH"
capture erase "$path\02_Processed_data\13_Running_variable\panu_temp.dta"
foreach can of local cans{
import delimited "$serverpath\02_Processed_data\13_Running_variable\temp_kotakorpietal2017_ch`can'.csv",clear
capture append using "$path\02_Processed_data\13_Running_variable\panu_temp.dta"
save "$path\02_Processed_data\13_Running_variable\panu_temp.dta", replace
}

rename id ID
save "$path\02_Processed_data\13_Running_variable\panu_temp.dta", replace


use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
merge 1:1 ID year using "$path\02_Processed_data\13_Running_variable\panu_temp.dta"


*(i) Check whether maximum votes of non-elected is always smaller than min of elected

bysort party_num: egen help1=min(p) if elected==1
bysort party_num: egen help2=max(p) if elected==0
bysort party_num: egen p_min=max(help1) 
bysort party_num: egen p_max=max(help2) 
gen p_diff=p_min-p_max
tab p_diff
unique party_num if p_diff<0


*(ii) Check whether we achieve convergence

gsort party_num -votes -elected name firstname
qui bysort party_num: egen votes_h_rank=rank(-votes), unique // generate rank of candidate on party list

local vars "p"
foreach var of local vars{
gsort party_num -`var' -votes -elected name firstname
qui bysort party_num: egen `var'_rank=rank(-`var'), unique // generate rank of candidate on party list
gen `var'_rank_diff=votes_h_rank-`var'_rank
qui bysort party_num: egen `var'_rank_diff_max=max(`var'_rank_diff)
}

gen p_rank_diff_indi=0
replace p_rank_diff_indi=1 if p_rank_diff!=0 & p_rank!=.
bysort party_num: egen p_rank_diff_indi_party=max(p_rank_diff_indi)

gsort party_num votes_h_rank
br canton year name firstname party_num votes m elected_sum p *rank if p_rank_diff_indi_party!=0 & m>=500000

gsort party_num -votes
br ID canton year name firstname party_num votes elected p *rank if canton=="ZH" & year==1931 & partyname=="76 - Revolutionäre Marxistische Liga / Kommunisten"