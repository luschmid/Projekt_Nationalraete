use dta/AlternativeMargin2, clear
sort year SMD_district
save dta/AlternativeMargin2_sort, replace

use dta/Turnout_District_Level, clear 
sort year SMD_district
drop _merge
merge year SMD_district using dta/AlternativeMargin2_sort

drop if _merge==2 /* PR-districts not in sample */
drop if year>1927
replace margin=min_distance1 if year>1920  /* single-party unweighted */
*replace margin=min_distance if year>1920  /* multi-party unweighted */
*replace margin=min_distance1_iXv_i if year>1920  /* single-party weighted */
*replace margin=min_distance_iXv_i if year>1920  /* multi-party weighted */

gen marginBL=.
replace marginBL=margin if year<1920                
replace marginBL=min_distance1*magnitude if year>1920  /* Blais-Lago normalize single-party unweighted with district magnitude */

*************************************************************************************
*************************************************************************************
*************************************************************************************
eststo clear
reg turnout margin if year<1920, cluster(SMD_district)
eststo
reg turnout marginBL if year<1920, cluster(SMD_district)
eststo
reg turnout C1 if year<1920, cluster(SMD_district)
eststo
reg turnout C2 if year<1920, cluster(SMD_district)
eststo
reg turnout margin if year>1920, cluster(PR_district)
eststo
reg turnout marginBL if year>1920, cluster(PR_district)
eststo
reg turnout C1 if year>1920, cluster(PR_district)
eststo
reg turnout C2 if year>1920, cluster(PR_district)
eststo
esttab using tables/TableB1.tex, replace se b(%9.3f) se(%9.3f) r2 nostar mtitles("SMD" "SMD" "SMD" "SMD" "MMD" "MMD" "MMD" "MMD")  ///
coeflabels(_cons "Constant" margin "Traditional" marginBL "B-L" C1 "G-S 1" C2 "G-S 2") ///
prehead("\begin{tabular}{l*{@M}{rr}}" "\hline") posthead(\hline) postfoot("\hline" "\end{tabular}")  

*************************************************************************************
*************************************************************************************
*************************************************************************************
use dta/Panel.dta, clear
keep if country=="SWI"
*replace margin=min_distance if m>1 /* multi-party measure */

eststo clear
reg turnout margin if m==1, robust
eststo
reg turnout marginBL if m==1, robust
eststo
reg turnout C1 if m==1, robust
eststo
reg turnout C2 if m==1, robust
eststo
reg turnout margin if m>1, cluster(id)
eststo
reg turnout marginBL if m>1, cluster(id)
eststo
reg turnout C1 if m>1, cluster(id)
eststo
reg turnout C2 if m>1, cluster(id)
eststo
esttab using tables/TableB2.tex, replace se b(%9.3f) se(%9.3f) r2 nostar mtitles("SMD" "SMD" "SMD" "SMD" "MMD" "MMD" "MMD" "MMD")  ///
coeflabels(_cons "Constant" margin "Traditional" marginBL "B-L" C1 "G-S 1" C2 "G-S 2") ///
prehead("\begin{tabular}{l*{@M}{rr}}" "\hline") posthead(\hline) postfoot("\hline" "\end{tabular}")  

