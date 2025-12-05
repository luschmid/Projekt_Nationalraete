*figure of incumbency effect at different lags
clear
use $DATAPATH/esmall
save $DATAPATH/tempdata,replace
clear
quietly capture erase $DATAPATH/incumbency_e.dta
set obs 1
gen temp=1
save $DATAPATH/incumbency_e, replace


* Prepare tempdata, create BINS by forcing variable (used only for RD graphs)
clear
use $DATAPATH/tempdata
cap drop y
gen y=nextElected
global DB=1 /* bin width */
egen bsbin=cut(x), at(-101($DB)101)
replace bsbin=bsbin+$DB/2
gen yxexist = (y!=. & x!=.)
bysort bsbin: egen N=count(yxexist)
drop yxexist
save $DATAPATH/tempdata,replace


rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
replace y=nextElected6
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly append using $DATAPATH/incumbency_e
save $DATAPATH/incumbency_e,replace


use $DATAPATH/tempdata
replace y=nextElected
rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
replace y=nextElected5
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly append using $DATAPATH/incumbency_e
save $DATAPATH/incumbency_e,replace

use $DATAPATH/tempdata
replace y=nextElected
rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
replace y=nextElected4
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly append using $DATAPATH/incumbency_e
save $DATAPATH/incumbency_e,replace

use $DATAPATH/tempdata
replace y=nextElected
rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
replace y=nextElected3
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly append using $DATAPATH/incumbency_e
save $DATAPATH/incumbency_e,replace

use $DATAPATH/tempdata
replace y=nextElected
rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
replace y=nextElected2
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly append using $DATAPATH/incumbency_e
save $DATAPATH/incumbency_e,replace

use $DATAPATH/tempdata
replace y=nextElected
rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth */
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
replace y=nextElected
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly append using $DATAPATH/incumbency_e
save $DATAPATH/incumbency_e,replace

use $DATAPATH/tempdata
replace y=nextElected
rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
replace y=lastElected
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly append using $DATAPATH/incumbency_e
save $DATAPATH/incumbency_e,replace

use $DATAPATH/tempdata
replace y=nextElected
rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
replace y=lastElected2
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly append using $DATAPATH/incumbency_e
save $DATAPATH/incumbency_e,replace

gen period=_n-3 if _n<3
replace period=_n-2 if _n>2
drop if temp==1
drop temp

compress


/*graph twoway scatter estimate min95 max95 period, msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) ///
	connect(l l l) lcolor(black gray gray) lpattern(solid dot dot) 
*/
graph twoway scatter estimate min95 max95 period if period<0, msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) ///
	connect(l l l) lcolor(black gray gray) lpattern(solid dot dot) ///
|| scatter estimate min95 max95 period if period>=0, msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) ///
	connect(l l l) lcolor(black gray gray) lpattern(solid dot dot) ///
|| , legend(off) xlabel(-2(1)6) b2title("", ring(4)) xtitle("Period") ytitle("") saving("$GRAPHPATH/IncumbencyEffect_Periods.gph",replace)


graph export "$OUTPATH/Figure_parliament_incumbency_periods.pdf",  fontface(Helvetica)  replace
graph export "$OUTPATH/Figure_parliament_incumbency_periods.eps",  fontface(Helvetica)  replace


clear
quietly capture erase $DATAPATH/incumbency_e.dta
use $DATAPATH/tempdata
eststo clear

*Incumbency effect 

quietly replace y=nextElected
quietly rdob_mod2 y x
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly rd y x, bw($band)  n(100) mbw(100)
eststo
estadd scalar Bandwidth = round($band,.001) 
global TITLE="Incumbency effect (Parliamentary elections)"
global XTITLE="Electoral closeness"
global YTITLE="Elected in subsequent election"
global B2TITLE=""
global GPHSTRING="$GRAPHPATH/IncumbencyEffect.gph"
global MARKERSIZE="vsmall"
MainRDgraph

esttab  using $OUTPATH/Incumbency.smcl, replace type label nonumber nodepvar se stats(N Bandwidth) mtitles("") ///
 title("Incumbency effect in parliamentary elections" ) 

* Bandwidth robustness
* Note: bwplot drops observations with missing $OUTCOME
clear
use $DATAPATH/tempdata
global oldOutcome = "$OUTCOME"
global OUTCOME = "nextElected"
bwplot nextElected x "Incumbency effect" "Bandwidth" " "
global OUTCOME=  "$oldOutcome"
global drop oldOutcome

* Graphs for incumbency advantage 
* global MARKERSIZE="vsmall"
* do incumbency_effect
