clear
clear matrix
set more off
set matsize 1000
use ..\dta\MergedDataBeforeCollapse

********************************************************************************************

foreach parti in RV SV DNA V SP KRF H FRP PP Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2{
	gen margin_`parti'=-mindiff`parti'p
	replace margin_`parti'= mindiff`parti'n if  mindiff`parti'p > mindiff`parti'n
    }

	
************ DENSITY FIGURES ***********************

foreach parti in SV DNA V SP KRF H FRP {
*foreach parti in RV SV DNA V SP KRF H FRP PP Indep1 Other1 JointL JointR1{
    twoway (hist margin_`parti' if LVoteShare`parti'!=0 & margin_`parti'>-0.01 & margin_`parti'<0, start(-0.01) width(0.001) freq bfcolor(gs10) blcolor(gs6)) ///
	(hist margin_`parti' if LVoteShare`parti'!=0 & margin_`parti'>0 & margin_`parti'<0.01, start(0.0) width(0.001) freq bfcolor(gs5) blcolor(gs6)), xtitle("") ytitle("") title(`parti') legend(off)
	graph save ..\figures\gph\_`parti'.gph, replace
	}
cd ..\figures\gph\

graph combine _SV.gph _DNA.gph _V.gph _SP.gph _KRF.gph _H.gph _FRP.gph , ycommon
*graph combine _RV.gph _SV.gph _DNA.gph _V.gph _SP.gph _KRF.gph _H.gph _FRP.gph _PP.gph _Indep1.gph _Other1.gph _JointL.gph _JointR1.gph
graph export ..\DensityTest.eps, replace
cd ..\..\do
