*use ..\dta\MergedDataBeforeStandardizing  /* THIS IS FOR SEEING JUMP IN RAW DATA */
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

sort knr ele
merge knr ele using "../dta/val med simuleringar DH4.dta"
drop _merge
bysort ele: sum dist_coal_margin
sort knr ele
merge knr ele using "../dta/val med simuleringar4.dta", update
bysort ele: sum dist_coal_margin
drop _merge

keep if DPTAXres!=.    /* limiting sample to only 1122 obs. where we also have simulations for distance to threshold */
keep if abs(dist_coal_margin)<.5  /* limiting sample to only 1122 obs. where we also have simulations for distance to threshold */

foreach depvar in Propertytax DPTAXres PTAX120sqm Usercharges Childcare Education Elderlycare Healthcare Culture Transport Administration Other {
egen `depvar'mean=mean(`depvar')
replace `depvar'=`depvar'-`depvar'mean
}

/* creating bins */

xtile bins=dist_coal_margin if dist_coal_margin<0, nq(6)
xtile bins2=dist_coal_margin if dist_coal_margin>=0, nq(25)
replace bins=bins2+6 if missing(bins)
drop bins2

/*   BINSCATTER */ 

*foreach depvar in Propertytax DPTAXres PTAX120sqm Usercharges Childcare Education Elderlycare Healthcare Culture Transport Administration Other {
*binscatter `depvar' dist_coal_margin, rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(`depvar') xscale(range(-0.4 0.4)) xlabel(-0.4(0.1)0.4) xtick(-0.4(0.1)0.4) yscale(range(-1 1)) ylabel(-1(0.5)1) ymtick(-1(0.5)1)
*graph save ..\figures\gph\_`depvar'q.gph, replace
*}


binscatter DPTAXres dist_coal_margin, graphregion(color(ltbluishgray)) rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(Property tax) xscale(range(-0.4 0.4)) graphregion(fcolor(white)) xscale(off) xmtick(,nolabels) yscale(range(-1 1)) ylabel(-1(0.5)1) ymtick(-1(0.5)1)
*graph play BinscatterED
graph save ..\figures\gph\_DPTAXresq.gph, replace

binscatter Usercharges dist_coal_margin, graphregion(color(ltbluishgray)) rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(User charges) xscale(range(-0.4 0.4)) graphregion(fcolor(white)) xscale(off) xmtick(,nolabels) yscale(range(-1 1)) yscale(off)  ymtick(, nolabels) 
*graph play BinscatterED
graph save ..\figures\gph\_Userchargesq.gph, replace

binscatter Childcare dist_coal_margin, graphregion(color(ltbluishgray)) rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(Child care) xscale(range(-0.4 0.4)) graphregion(fcolor(white)) xscale(off) xmtick(,nolabels) yscale(range(-1 1)) yscale(off)  ymtick(, nolabels) 
*graph play BinscatterED
graph save ..\figures\gph\_Childcareq.gph, replace

binscatter Education dist_coal_margin, graphregion(color(ltbluishgray)) rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(Education) xscale(range(-0.4 0.4)) graphregion(fcolor(white)) xscale(off) xmtick(,nolabels) yscale(range(-1 1)) ylabel(-1(0.5)1) ymtick(-1(0.5)1)
*graph play BinscatterED
graph save ..\figures\gph\_Educationq.gph, replace

binscatter Elderlycare dist_coal_margin, graphregion(color(ltbluishgray)) rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(Elderly care) xscale(range(-0.4 0.4)) graphregion(fcolor(white)) xscale(off) xmtick(,nolabels) yscale(range(-1 1)) yscale(off)  ymtick(, nolabels) 
*graph play BinscatterED
graph save ..\figures\gph\_Elderlycareq.gph, replace

binscatter Healthcare dist_coal_margin, graphregion(color(ltbluishgray)) rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(Health care) xscale(range(-0.4 0.4)) graphregion(fcolor(white)) xscale(off) xmtick(,nolabels) yscale(range(-1 1)) yscale(off)  ymtick(, nolabels) 
*graph play BinscatterED
graph save ..\figures\gph\_Healthcareq.gph, replace

binscatter Culture dist_coal_margin, graphregion(color(ltbluishgray)) rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(Culture) xscale(range(-0.4 0.4)) graphregion(fcolor(white)) xlabel(-0.4(0.1)0.4) xtick(-0.4(0.1)0.4) yscale(range(-1 1)) ylabel(-1(0.5)1) ymtick(-1(0.5)1)
*graph play BinscatterED
graph save ..\figures\gph\_Cultureq.gph, replace

binscatter Transport dist_coal_margin, graphregion(color(ltbluishgray)) rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(Transport) xscale(range(-0.4 0.4)) graphregion(fcolor(white)) xlabel(-0.4(0.1)0.4) xtick(-0.4(0.1)0.4) yscale(range(-1 1)) yscale(off)  ymtick(, nolabels) 
*graph play BinscatterED
graph save ..\figures\gph\_Transportq.gph, replace

binscatter Administration dist_coal_margin, graphregion(color(ltbluishgray)) rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(Administration) xscale(range(-0.4 0.4)) graphregion(fcolor(white)) xlabel(-0.4(0.1)0.4) xtick(-0.4(0.1)0.4) yscale(range(-1 1)) yscale(off)  ymtick(, nolabels) 
*graph play BinscatterED
graph save ..\figures\gph\_Administrationq.gph, replace

cd ..\figures\gph\

graph combine _DPTAXresq.gph _Userchargesq.gph _childcareq.gph _educationq.gph _elderlycareq.gph _Healthcareq.gph _cultureq.gph _transportq.gph _Administrationq.gph , rows(3) cols(3) ycommon title("") graphregion(fcolor(white))
graph save ..\Figure3.gph, replace
graph export ..\Figure3.tif, replace
graph export ..\Figure3.eps, replace
*graph export ..\..\..\Parties\figuresFINAL\Figure3.eps, replace


cd ..\..\do\
