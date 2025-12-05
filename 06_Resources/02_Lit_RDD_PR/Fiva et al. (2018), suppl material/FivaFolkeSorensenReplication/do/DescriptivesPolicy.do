use ..\dta\MergedDataBeforeStandardizing, clear

sort knr ele
merge knr ele using "../dta/val med simuleringar DH4.dta"
drop _merge
bysort ele: sum dist_coal_margin
sort knr ele
merge knr ele using "../dta/val med simuleringar4.dta", update
bysort ele: sum dist_coal_margin

keep if abs(dist_coal_margin)<.5  /* limiting sample to only 1122 obs. where we also have simulations for distance to threshold */


/* Table 1: Descriptive statistics on fiscal policy outcomes  */

sutex DPTAX DPTAXres Fees childcare education elderlycare healthsocial culture transport centraladm other, minmax title(Descriptive Statistics Fiscal Policy \label{DescrStatsPolicyBeforeStandardizing}) ///
file(../tables/DescrStatsPolicyBeforeStandardizing.tex) replace

