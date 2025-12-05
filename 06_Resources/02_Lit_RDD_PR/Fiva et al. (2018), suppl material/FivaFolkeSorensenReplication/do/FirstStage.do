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

*************** FIRST STAGE *********************

/* STANDARDIZING THE CONTROL TO HAVE SD=1*/
replace LeftRightControl=LeftRightControl/.5510006

eststo clear
*reg LeftRightIndex *treat_w_025* *dum_w_025*,  cluster(knr) 
*eststo
reg LeftRightIndex *treat_w_025* *dum_w_025* LeftRightControl,  cluster(knr) 
test RV_treat_w_025 SV_treat_w_025 V_treat_w_025 SP_treat_w_025 KRF_treat_w_025 H_treat_w_025 FRP_treat_w_025 PP_treat_w_025 MDG_treat_w_025 JointL_treat_w_025 JOINTR_treat_w_025 Various_treat_w_025
eststo
reg LeftRightIndex *treat_w_025* *dum_w_025* MajRight first* second* LeftRightControl, cluster(knr)
test RV_treat_w_025 SV_treat_w_025 V_treat_w_025 SP_treat_w_025 KRF_treat_w_025 H_treat_w_025 FRP_treat_w_025 PP_treat_w_025 MDG_treat_w_025 JointL_treat_w_025 JOINTR_treat_w_025 Various_treat_w_025
eststo
reg LeftRightIndex *treat_w_025* *dum_w_025* MajRight first* second* third* LeftRightControl, cluster(knr)
test RV_treat_w_025 SV_treat_w_025 V_treat_w_025 SP_treat_w_025 KRF_treat_w_025 H_treat_w_025 FRP_treat_w_025 PP_treat_w_025 MDG_treat_w_025 JointL_treat_w_025 JOINTR_treat_w_025 Various_treat_w_025
eststo
esttab using ../tables/FirstStage.tex, replace style(tex) r2 se b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) drop(_cons *dum* MajRight first* second* third*) ///
title(First Stage \label{FirstStage})
