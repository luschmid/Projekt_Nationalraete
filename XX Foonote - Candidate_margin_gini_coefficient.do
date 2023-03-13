* 1. Set directories 

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

set more off

* 2. Read in data, generate running variable, election in next election,
* 	 party variables, and cantonal dummy variable

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear

* 3. Calculate candidate margin

bysort year canton list: egen temp1=min(votes) if elected==1
bysort year canton list: egen temp2=max(votes) if elected==0
bysort year canton list: egen votes_min=max(temp1) 
bysort year canton list: egen votes_max=max(temp2) 
drop temp1 temp2 

gen candidate_margin=.
replace candidate_margin=votes-votes_max if elected==1
replace candidate_margin=votes-votes_min if elected==0

bysort elected: sum candidate_margin
corr votemargin candidate_margin

gen candidate_margin_rel=candidate_margin/eligible_cant

* 3. Calculate Gini coefficient per group

bysort year canton list: gen no_candidates=_N
egen list_per_year=group(year canton list) if no_candidates>1
gen gini=.

sum list_per_year
forvalues i=1(1)`r(max)'{
if mod(`i',100)==0 display `i'
quietly ineqdeco votes if list_per_year==`i'
quietly replace gini=`r(gini)' if list_per_year==`i'
}

* 4. Save file

keep ID year candidate_margin candidate_margin_rel gini list_per_year no_candidates
save "$path\02_Processed_data\13_Running_variable\ch_gini.dta", replace

* 5. Correlation between Gini and difference between votemargin and candidate margin

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
merge 1:1 ID year using "$path\02_Processed_data\13_Running_variable\ch_gini.dta"

sort list_per_year votemargin_rel
br votes list_per_year gini if list_per_year<10 

gen votemargin_rel_abs=abs(votemargin_rel)
gen candidate_margin_rel_abs=abs(candidate_margin_rel)

gen diff_vm_cm=candidate_margin_rel_abs-votemargin_rel_abs

pwcorr gini diff_vm_cm, sig