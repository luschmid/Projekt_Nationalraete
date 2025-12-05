use ..\dta\MergedData, clear

sort knr ele
merge knr ele using "../dta/val med simuleringar DH4.dta"
drop _merge
bysort ele: sum dist_coal_margin
sort knr ele
merge knr ele using "../dta/val med simuleringar4.dta", update
bysort ele: sum dist_coal_margin

gen dist_coal_marginxMajLeft=dist_coal_margin*MajLeft

keep if abs(dist_coal_margin)<.5  /* limiting sample to only 1122 obs. where we also have simulations for distance to threshold */


********************************************************************************************

foreach parti in RV SV DNA V SP KRF H FRP PP {
	gen margin_`parti'=-mindiff`parti'p
	replace margin_`parti'= mindiff`parti'n if  mindiff`parti'p > mindiff`parti'n
    }
*********
foreach parti in RV SV DNA V SP KRF H FRP PP {
gen bins`parti'=.
replace bins`parti'=1 if margin_`parti'>-0.010 & margin_`parti'<-0.009
replace bins`parti'=2 if margin_`parti'>-0.009 & margin_`parti'<-0.008
replace bins`parti'=3 if margin_`parti'>-0.008 & margin_`parti'<-0.007
replace bins`parti'=4 if margin_`parti'>-0.007 & margin_`parti'<-0.006
replace bins`parti'=5 if margin_`parti'>-0.006 & margin_`parti'<-0.005
replace bins`parti'=6 if margin_`parti'>-0.005 & margin_`parti'<-0.004
replace bins`parti'=7 if margin_`parti'>-0.004 & margin_`parti'<-0.003
replace bins`parti'=8 if margin_`parti'>-0.003 & margin_`parti'<-0.002
replace bins`parti'=9 if margin_`parti'>-0.002 & margin_`parti'<-0.001
replace bins`parti'=10 if margin_`parti'>-0.001 & margin_`parti'<0.000
replace bins`parti'=11 if margin_`parti'>0.000 & margin_`parti'<0.001
replace bins`parti'=12 if margin_`parti'>0.001 & margin_`parti'<0.002
replace bins`parti'=13 if margin_`parti'>0.002 & margin_`parti'<0.003
replace bins`parti'=14 if margin_`parti'>0.003 & margin_`parti'<0.004
replace bins`parti'=15 if margin_`parti'>0.004 & margin_`parti'<0.005
replace bins`parti'=16 if margin_`parti'>0.005 & margin_`parti'<0.006
replace bins`parti'=17 if margin_`parti'>0.006 & margin_`parti'<0.007
replace bins`parti'=18 if margin_`parti'>0.007 & margin_`parti'<0.008
replace bins`parti'=19 if margin_`parti'>0.008 & margin_`parti'<0.009
replace bins`parti'=20 if margin_`parti'>0.009 & margin_`parti'<0.010
}

*************** VOTE SHARES AROUND THRESHOLD FIGURES *****************

	cd ..\figures\gph
	foreach parti in SV DNA V SP KRF H FRP {
	binscatter LVoteShare`parti' margin_`parti' if abs(margin_`parti')<0.01, graphregion(color(white)) xq(bins`parti') rd(0) ytitle("") xtitle("") title(`parti') linetype(none)
	*graph play BinscatterED
	graph save _`parti'.gph, replace
	}

	graph combine _SV.gph _DNA.gph _V.gph _SP.gph _KRF.gph _H.gph _FRP.gph ,ycommon graphregion(color(white))
	graph export ..\VoteshareThresholdDistance.eps, replace

	
cd ..\..\do

*************** VOTE SHARES AROUND THRESHOLD FIGURES *****************

	cd ..\figures\gph
	foreach parti in SV DNA V SP KRF H FRP {
	binscatter LVoteShare`parti' margin_`parti' if abs(margin_`parti')<0.01, graphregion(color(white))xq(bins`parti') rd(0) ytitle("") xtitle("") title(`parti') linetype(lfit)
	*graph play BinscatterED
	graph save _`parti'.gph, replace
	}

	graph combine _SV.gph _DNA.gph _V.gph _SP.gph _KRF.gph _H.gph _FRP.gph ,ycommon graphregion(color(white))
	graph export ..\VoteshareThresholdDistance_Line.eps, replace

	
cd ..\..\do

