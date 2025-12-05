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

**************** ADDING POLYNOMIALS *************

gen first=dist_coal_margin
gen firstxMajLeft=dist_coal_margin*MajLeft

gen second=dist_coal_margin^2
gen secondxMajLeft=(dist_coal_margin^2)*MajLeft

gen third=dist_coal_margin^3
gen thirdxMajLeft=(dist_coal_margin^3)*MajLeft

gen MajRight=1-MajLeft


*************** REGRESSIONS w PREDETERMINED*********************

*rdplot laglogpop dist_coal_margin, p(2) numbinl(5) numbinr(10)
*rdplot lagRural dist_coal_margin, p(2) numbinl(5) numbinr(10)


foreach depvar in laglogpop lagchildren lagyoung lagelderly lagunemployment by1911 lagRural {
eststo clear
ivregress 2sls `depvar' (LeftRightIndex=*treat_w_025*) *dum_w_025* LeftRightControl,  cluster(knr) 
eststo
estadd local Method "IND"
estadd local VoteSHcontrol "No"
estadd local TimeFE "No"
estadd local Controls "No"
reg `depvar' MajRight first* second* , cluster(knr)
eststo
estadd local Method "MAJ"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Second"
estadd local Bandwidth "50"
reg `depvar' MajRight first* second* third* , cluster(knr)
eststo
estadd local Method "MAJ"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Third"
estadd local Bandwidth "50"
ivregress 2sls `depvar' (LeftRightIndex=*treat_w_025*) *dum_w_025* MajRight first* second* LeftRightControl , cluster(knr)
eststo
estadd local Method "BOTH"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Second"
estadd local Bandwidth "50"
ivregress 2sls `depvar' (LeftRightIndex=*treat_w_025*) *dum_w_025* MajRight first* second* third* LeftRightControl , cluster(knr)
eststo
estadd local Method "BOTH"
estadd local TimeFE "No"
estadd local Controls "No"
estadd local Polynom "Third"
estadd local Bandwidth "50"
esttab using ../tables/Predetermined_`depvar'.tex, replace style(tex) ar2 se b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) drop(_cons *dum_w_* *first* *second* *third* *Control) ///
scalars("Method Method" "TimeFE TimeFE" "Controls DemContr" "Bandwidth Bandwidth" "Polynom Polynom") title(Political Representation and Fiscal Policy: `depvar' \label{Table`depvar'})
}
