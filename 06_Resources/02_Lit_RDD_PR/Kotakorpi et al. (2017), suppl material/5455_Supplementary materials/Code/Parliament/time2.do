*figure of main effect at different lags

save $DATAPATH/tempdata.dta, replace
clear

quietly capture erase $DATAPATH/parmlags_e.dta
set obs 1
gen temp=1
save $DATAPATH/parmlags_e, replace

clear
use $DATAPATH/tempdata


clear 
use $DATAPATH/tempdata
quietly replace y=next1_3tulo
quietly rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly replace y=next21_23tulo
quietly parmby "rd y x if (alwaysIn5==1 | neverIn5==1), bw($band) n(100) mbw(100)" , norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_e
quietly save $DATAPATH/parmlags_e,replace

clear 
use $DATAPATH/tempdata
quietly replace y=next1_3tulo
quietly rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly replace y=next17_19tulo
quietly parmby "rd y x if (alwaysIn4==1 | neverIn4==1), bw($band) n(100) mbw(100)" , norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_e
quietly save $DATAPATH/parmlags_e,replace

clear 
use $DATAPATH/tempdata
quietly replace y=next1_3tulo
quietly rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly replace y=next13_15tulo
quietly parmby "rd y x if (alwaysIn3==1 | neverIn3==1), bw($band) n(100) mbw(100) ", norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_e
quietly save $DATAPATH/parmlags_e,replace

clear 
use $DATAPATH/tempdata
quietly replace y=next1_3tulo
quietly rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly replace y=next9_11tulo
quietly parmby "rd y x if (alwaysIn2==1 | neverIn2==1), bw($band) n(100) mbw(100) ", norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_e
quietly save $DATAPATH/parmlags_e,replace

clear 
use $DATAPATH/tempdata
quietly replace y=next1_3tulo
quietly rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly replace y=next5_7tulo
quietly parmby "rd y x if (elected==1 & nextElected==1) | (elected==0 & nextElected==0) , bw($band) n(100) mbw(100) ", norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_e
quietly save $DATAPATH/parmlags_e,replace

clear 
use $DATAPATH/tempdata
quietly replace y=next1_3tulo
quietly rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly replace y=next1_3tulo
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_e
quietly save $DATAPATH/parmlags_e,replace

clear 
use $DATAPATH/tempdata
quietly replace y=next1_3tulo
quietly rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly replace y=last1_3tulo
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_e
quietly save $DATAPATH/parmlags_e,replace

use $DATAPATH/tempdata
quietly replace y=next1_3tulo
quietly rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly replace y=last5_7tulo
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_e
quietly save $DATAPATH/parmlags_e,replace


gen period=_n
replace period=period-3
drop if temp==1
drop temp
compress


graph twoway scatter estimate min95 max95 period if period<0, msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) connect(l l l) lcolor(black gray gray) lpattern(solid dot dot) ///
|| scatter estimate min95 max95 period if period>=0, msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) connect(l l l) lcolor(black gray gray) lpattern(solid dot dot) ///
|| , legend(off) xlabel(-2(1)5) b2title("", ring(4)) xtitle("Electoral period") ytitle("") saving("$GRAPHPATH/EarningsPeriodsByStatus1.gph",replace)

graph export "$OUTPATH/Figure_parliament_earnings_periods_bystatus1.pdf", fontface(Helvetica)  replace
graph export "$OUTPATH/Figure_parliament_earnings_periods_bystatus1.eps", fontface(Helvetica)  replace

clear
quietly capture erase $DATAPATH/parmlags_e.dta
use $DATAPATH/tempdata
