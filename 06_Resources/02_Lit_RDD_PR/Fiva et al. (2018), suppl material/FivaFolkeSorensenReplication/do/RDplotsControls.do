*use ..\dta\MergedDataBeforeStandardizing  /* THIS IS FOR SEEING JUMP IN RAW DATA */
use ..\dta\MergedData, clear

sort knr ele
merge knr ele using "../dta/val med simuleringar DH4.dta"
drop _merge
bysort ele: sum dist_coal_margin
sort knr ele
merge knr ele using "../dta/val med simuleringar4.dta", update
bysort ele: sum dist_coal_margin
drop _merge

keep if abs(dist_coal_margin)<.5  /* limiting sample to only 1122 obs. where we also have simulations for distance to threshold */

/* creating bins */

xtile bins=dist_coal_margin if dist_coal_margin<0, nq(6)
xtile bins2=dist_coal_margin if dist_coal_margin>=0, nq(25)
replace bins=bins2+6 if missing(bins)
drop bins2

/*   BINSCATTER -controls */ 

foreach depvar in laglogpop lagRural lagchildren lagyoung lagelderly lagunemployment by1911  {
egen `depvar'mean=mean(`depvar')
replace `depvar'=`depvar'-`depvar'mean
binscatter `depvar' dist_coal_margin, graphregion(color(white)) rd(0) linetype(qfit) xq(bins) xtitle("") ytitle("") subtitle(`depvar') xscale(range(-0.4 0.4)) xlabel(-0.4(0.1)0.4) xtick(-0.4(0.1)0.4) yscale(range(-1 1)) ylabel(-1(0.5)1) ymtick(-1(0.5)1)
*graph play BinscatterED
graph save ..\figures\gph\_`depvar'q.gph, replace
}


cd ..\figures\gph\

graph combine _laglogpopq.gph _lagRuralq.gph _lagchildrenq.gph _lagyoungq.gph _lagelderlyq.gph _lagunemploymentq.gph _by1911q.gph , rows(3) cols(3) ycommon title("") 
*graph play RDplot_Titles2
*graph play WhiteBackground
graph save ..\RDplot_Controls.gph, replace
graph export ..\RDplot_Controls.tif, replace
graph export ..\RDplot_Controls.eps, replace
*graph export ../../../Parties/figures/RDplot_Controls.eps, replace

cd ..\..\do\
