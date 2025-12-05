use ..\dta\MergedData, clear

rename DPTAX Propertytax
rename Fees Usercharges
rename childcare Childcare
rename elderlycare Elderlycare
rename education Education
rename healthsocial Healthcare
rename culture Culture
rename transport Transport
rename centraladm Administration
rename other Other

foreach depvar in Propertytax DPTAXres PTAX120sqm Usercharges Childcare Education Elderlycare Healthcare Culture Transport Administration Other {
egen `depvar'mean=mean(`depvar')
replace `depvar'=`depvar'-`depvar'mean
}

sort knr ele
merge knr ele using "../dta/val med simuleringar DH4.dta"
drop _merge
bysort ele: sum dist_coal_margin
sort knr ele
merge knr ele using "../dta/val med simuleringar4.dta", update
bysort ele: sum dist_coal_margin
drop _merge

/* Local linear with IK bandwidth */

foreach depvar in DPTAXres Propertytax Usercharges Childcare Education Elderlycare Healthcare Culture Transport Administration Other {
rdrobust `depvar' dist_coal_margin, bwselect(IK)
}



/* Local linear with CCT bandwidth */

foreach depvar in DPTAXres Propertytax Usercharges Childcare Education Elderlycare Healthcare Culture Transport Administration Other {
rdrobust `depvar' dist_coal_margin, bwselect(CCT)
}



/* Local linear with varying bandwidth from min to max */

log using ../log/LOG_rdrobust.txt, text replace
forvalues i = 0.01(0.01)0.45 {
rdrobust DPTAXres dist_coal_margin, h(`i')
*estimate store bwidth`i'
}
log close


log using ../log/LOG_rdrobust2.txt, text replace
forvalues i = 0.01(0.01)0.45 {
rdrobust Propertytax dist_coal_margin, h(`i')
*estimate store bwidth`i'
}
log close




use ..\dta\rdrobust_bwidth.dta, clear

twoway (rbar CI_lower CI_upper bw if optimal==0, lwidth(vvthin) lpattern(solid) lcolor(gray) barwidth(.0001)) (scatter Coeff bw if optimal==0, msize(medium) mcolor(gray)) ///
(rbar CI_lower CI_upper bw if optimal==1, lwidth(medium) lpattern(solid) lcolor(orange) barwidth(.0001)) (scatter Coeff bw if optimal==1, msize(medium) mcolor(orange)) ///
, xline(0.166, lwidth(vvthin) lpattern(dash) lcolor(black)) xtitle(Bandwidth) title(Estimate and 95 pct CI) subtitle(Orange= IK optimal bw) legend(off)
graph export ..\figures\rdrobust_varying_bw.tif, replace
graph export ..\figures\rdrobust_varying_bw.eps, replace


twoway (rbar CI_lower CI_upper bw, lwidth(vvthin) lpattern(solid) lcolor(gray) barwidth(.0001)) (scatter Coeff bw, msize(medium) mcolor(gray)) ///
, xline(0.114, lwidth(vvthin) lpattern(dash) lcolor(gray)) xline(0.166, lwidth(vvthin) lpattern(dash) lcolor(gray)) xtitle(Bandwidth) ytitle(Estimated effect) subtitle("") legend(off) 
graph export ..\figures\rdrobust_varying_bw.tif, replace
graph export ..\figures\rdrobust_varying_bw.eps, replace


drop in 1
drop in 2
drop in 3
drop in 4
drop in 5
drop in 6
drop in 7
drop in 8
drop in 9
drop in 10
drop in 11
drop in 12
drop in 13
drop in 14
drop in 15
drop in 16
drop in 17
drop in 18
drop in 19
drop in 20
drop in 21
drop in 22


twoway (rbar CI_lower CI_upper bw, lwidth(vvthin) lpattern(solid) lcolor(gray) barwidth(.0001)) (scatter Coeff bw, msize(medium) mcolor(gray)) ///
, xline(0.114, lwidth(vvthin) lpattern(dash) lcolor(gray)) xline(0.166, lwidth(vvthin) lpattern(dash) lcolor(gray)) xtitle(Bandwidth) ytitle(Estimated effect) title(Residential property taxation) subtitle((Main dependent variable)) legend(off) graphregion(color(white))
graph export ..\figures\rdrobust_varying_bw.tif, replace
graph export ..\figures\rdrobust_varying_bw.eps, replace



use ..\dta\rdrobust_bwidth_DPTAX.dta, clear


drop in 1
drop in 2
drop in 3
drop in 4
drop in 5
drop in 6
drop in 7
drop in 8
drop in 9
drop in 10
drop in 11
drop in 12
drop in 13
drop in 14
drop in 15
drop in 16
drop in 17
drop in 18
drop in 19
drop in 20
drop in 21
drop in 22


twoway (rbar CI_lower CI_upper bw, lwidth(vvthin) lpattern(solid) lcolor(gray) barwidth(.0001)) (scatter Coeff bw, msize(medium) mcolor(gray)) ///
, xline(0.098, lwidth(vvthin) lpattern(dash) lcolor(gray)) xline(0.120, lwidth(vvthin) lpattern(dash) lcolor(gray)) xtitle(Bandwidth) ytitle(Estimated effect) title(Residential and commercial property taxation) subtitle((Alternative dependent variable)) legend(off) graphregion(color(white))
graph export ..\figures\rdrobust_varying_bw_DPTAX.tif, replace
graph export ..\figures\rdrobust_varying_bw_DPTAX.eps, replace

