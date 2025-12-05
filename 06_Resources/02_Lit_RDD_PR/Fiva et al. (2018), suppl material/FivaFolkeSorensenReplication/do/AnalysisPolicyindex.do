use ..\dta\MergedData, clear

sort knr ele
merge knr ele using "../dta/val med simuleringar DH4.dta"
drop _merge
bysort ele: sum dist_coal_margin
sort knr ele
merge knr ele using "../dta/val med simuleringar4.dta", update
bysort ele: sum dist_coal_margin

gen dist_coal_marginxMajLeft=dist_coal_margin*MajLeft

drop *DNA* /* REFERENCE CATEGORY */

keep if abs(dist_coal_margin)<.5  /* limiting sample to only 1122 obs. where we also have simulations for distance to threshold */

***************** T-TEST ************************
/*
g ZZ=1-MajLeft
eststo clear
estpost summarize DPTAX DPTAXres PTAX120sqm Fees childcare education elderlycare healthsocial culture transport centraladm other Total population* children young elderly unemployment by1911 rural ElectionPeriod0407 ElectionPeriod0811 SizeOfCouncil LVoteShareLEFT SeatShareLEFT if abs(dist_coal_margin)<.025 & MajLeft==1
eststo a
estpost summarize DPTAX DPTAXres PTAX120sqm Fees childcare education elderlycare healthsocial culture transport centraladm other Total population* children young elderly unemployment by1911 rural ElectionPeriod0407 ElectionPeriod0811 SizeOfCouncil LVoteShareLEFT SeatShareLEFT if abs(dist_coal_margin)<.025 & MajLeft==0
eststo b
estpost ttest DPTAX DPTAXres PTAX120sqm Fees childcare education elderlycare healthsocial culture transport centraladm other Total population* children young elderly unemployment by1911 rural ElectionPeriod0407 ElectionPeriod0811 SizeOfCouncil LVoteShareLEFT SeatShareLEFT if abs(dist_coal_margin)<.025 , by(ZZ)
eststo x
esttab a b x using ..\tables\Coalitions\DescrStatsMajLeft.tex, title(Summary Statistics for municipalities where less than 2.5 percent of the votes separate the left and right block\label{DescrStatsMajLeft}) style(tex) ///
stats(N) b(3) mtitles("Left-wing majority" "Right-wing majority" "Difference") collabels(Mean Estimate SD SE) starlevels(* 0.10 ** 0.05 *** 0.01) ///
cells("mean (fmt(3) pattern(1 1 0 0)) b (fmt(3) star pattern (0 0 1 1)) sd (fmt(3) par pattern(1 1 0 0)) se (fmt(3) par pattern(0 0 1 1))") ///
replace
*/

**************** ADDING POLYNOMIALS *************

gen first=dist_coal_margin
gen firstxMajLeft=dist_coal_margin*MajLeft

gen second=dist_coal_margin^2
gen secondxMajLeft=(dist_coal_margin^2)*MajLeft

gen third=dist_coal_margin^3
gen thirdxMajLeft=(dist_coal_margin^3)*MajLeft

gen fourth=dist_coal_margin^4
gen fourthxMajLeft=(dist_coal_margin^4)*MajLeft

gen MajRight=1-MajLeft



*************** REGRESSIONS *********************

foreach depvar in DPTAX DPTAXres Fees childcare education elderlycare healthsocial culture transport centraladm {
eststo clear
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* ,  cluster(knr) 
eststo
estadd local Method "IND"
estadd local VoteSHcontrol "No"
estadd local TimeFE "No"
estadd local Controls "No"
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* `depvar'Control,  cluster(knr) first 
eststo
estadd local Method "IND"
estadd local VoteSHcontrol "No"
estadd local TimeFE "No"
estadd local Controls "No"
reg `depvar' MajRight first* second* if abs(dist_coal_margin)<.5, cluster(knr)
eststo
estadd local Method "MAJ"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Second"
estadd local Bandwidth "50"
reg `depvar' MajRight first* second* third* if abs(dist_coal_margin)<.5, cluster(knr)
eststo
estadd local Method "MAJ"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Third"
estadd local Bandwidth "50"
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* MajRight first* second* if abs(dist_coal_margin)<.5, cluster(knr)
eststo
estadd local Method "BOTH"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Second"
estadd local Bandwidth "50"
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* MajRight first* second* third* if abs(dist_coal_margin)<.5, cluster(knr)
eststo
estadd local Method "BOTH"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Third"
estadd local Bandwidth "50"
esttab using ../tables/_`depvar'.tex, replace style(tex) ar2 se b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) drop(_cons *dum_w_* *first* *second* *third* *Control) ///
scalars("Method Method" "TimeFE TimeFE" "Controls DemContr" "Bandwidth Bandwidth" "Polynom Polynom") title(Political Representation and Fiscal Policy: `depvar' \label{Table`depvar'})
}


*************** REGRESSIONS w controls *********************

foreach depvar in DPTAX DPTAXres Fees childcare education elderlycare healthsocial culture transport centraladm {
eststo clear
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* population* children young elderly unemployment by1911 rural Ele*,  cluster(knr) 
eststo
estadd local Method "IND"
estadd local VoteSHcontrol "No"
estadd local TimeFE "No"
estadd local Controls "No"
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* `depvar'Control population* children young elderly unemployment by1911 rural Ele*,  cluster(knr) 
eststo
estadd local Method "IND"
estadd local VoteSHcontrol "No"
estadd local TimeFE "No"
estadd local Controls "No"
reg `depvar' MajRight first* second* population* children young elderly unemployment by1911 rural Ele* if abs(dist_coal_margin)<.5,  cluster(knr)
eststo
estadd local Method "MAJ"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Second"
estadd local Bandwidth "50"
reg `depvar' MajRight first* second* third* population* children young elderly unemployment by1911 rural Ele* if abs(dist_coal_margin)<.5, cluster(knr)
eststo
estadd local Method "MAJ"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Third"
estadd local Bandwidth "50"
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* MajRight first* second* population* children young elderly unemployment by1911 rural Ele* if abs(dist_coal_margin)<.5, cluster(knr)
eststo
estadd local Method "BOTH"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Second"
estadd local Bandwidth "50"
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* MajRight first* second* third* population* children young elderly unemployment by1911 rural Ele* if abs(dist_coal_margin)<.5, cluster(knr)
eststo
estadd local Method "BOTH"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Third"
estadd local Bandwidth "50"
esttab using ../tables/_`depvar'Controls.tex, replace style(tex) ar2 se b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) drop(_cons *dum_w_* *first* *second* *third* *Control population* children young elderly unemployment by1911 rural Ele*) ///
scalars("Method Method" "TimeFE TimeFE" "Controls DemContr" "Bandwidth Bandwidth" "Polynom Polynom") title(Political Representation and Fiscal Policy: `depvar' w. Controls \label{Table`depvar'})
}


*************** REGRESSIONS w controls + votesharecontrol *********************

foreach depvar in DPTAX DPTAXres Fees childcare education elderlycare healthsocial culture transport centraladm {
eststo clear
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* population* children young elderly unemployment by1911 rural Ele*,  cluster(knr) 
eststo
estadd local Method "IND"
estadd local VoteSHcontrol "No"
estadd local TimeFE "No"
estadd local Controls "No"
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* `depvar'Control population* children young elderly unemployment by1911 rural Ele*,  cluster(knr) 
eststo
estadd local Method "IND"
estadd local VoteSHcontrol "No"
estadd local TimeFE "No"
estadd local Controls "No"
reg `depvar' MajRight first* second* population* children young elderly unemployment by1911 rural Ele* if abs(dist_coal_margin)<.5,  cluster(knr)
eststo
estadd local Method "MAJ"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Second"
estadd local Bandwidth "50"
reg `depvar' MajRight first* second* third* population* children young elderly unemployment by1911 rural Ele* if abs(dist_coal_margin)<.5, cluster(knr)
eststo
estadd local Method "MAJ"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Third"
estadd local Bandwidth "50"
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* `depvar'Control MajRight first* second* population* children young elderly unemployment by1911 rural Ele* if abs(dist_coal_margin)<.5, cluster(knr)
eststo
estadd local Method "BOTH"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Second"
estadd local Bandwidth "50"
ivregress 2sls `depvar' (`depvar'Index=*treat_w_025*) *dum_w_025* `depvar'Control MajRight first* second* third* population* children young elderly unemployment by1911 rural Ele* if abs(dist_coal_margin)<.5, cluster(knr)
eststo
estadd local Method "BOTH"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Third"
estadd local Bandwidth "50"
esttab using ../tables/_`depvar'Controls_.tex, replace style(tex) ar2 se b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) drop(_cons *dum_w_* *first* *second* *third* *Control population* children young elderly unemployment by1911 rural Ele*) ///
scalars("Method Method" "TimeFE TimeFE" "Controls DemContr" "Bandwidth Bandwidth" "Polynom Polynom") title(Political Representation and Fiscal Policy: `depvar' w. Controls \label{Table`depvar'})
}
