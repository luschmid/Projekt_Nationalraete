*figure of main effect at different lags (electoral periods)

save $DATAPATH/tempdata.dta, replace
clear

quietly capture erase $DATAPATH/parmlags_e.dta
set obs 1
gen temp=1
quietly save $DATAPATH/parmlags_e,replace

clear
use $DATAPATH/tempdata
foreach var of varlist next25_27tulo next21_23tulo next17_19tulo next13_15tulo next9_11tulo next5_7tulo next1_3tulo last1_3tulo last5_7tulo {
clear 
use $DATAPATH/tempdata
quietly replace y=next1_3tulo
quietly quietly rdob_mod2 y x  
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly replace y=`var'
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_e
quietly save $DATAPATH/parmlags_e,replace
}

gen period=_n-2
drop if temp==1
drop temp

compress


graph twoway scatter estimate min95 max95 period if period<1,  ///
 msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) ///
 connect(l l l) lcolor(black gray gray) lpattern(solid dot dot) ///
|| scatter estimate min95 max95 period if period>=1,  xlabel(none -1 "-2" 0 "-1" 1 2 3 4 5 6 7) ///
msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) connect(l l l) ///
lcolor(black gray gray) lpattern(solid dot dot) ///
|| , legend(off) b2title("", ring(4)) xtitle("Electoral period") ytitle(" €/year") ///
 saving("$GRAPHPATH/EarningsPeriods.gph",replace)

graph export "$OUTPATH/Figure_parliament_earnings_periods.pdf", fontface(Helvetica)  replace
graph export "$OUTPATH/Figure_parliament_earnings_periods.eps", fontface(Helvetica)  replace



*Figure of main effect at different lags, comparing those elected at 0 and t, to those not elected at 0 and t

clear

quietly capture erase $DATAPATH/parmlags_e.dta
set obs 1
gen temp=1
quietly save $DATAPATH/parmlags_e,replace

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
replace y=next25_27tulo
quietly parmby "rd y x if (elected==1 & nextElected6==1) | (elected==0 & nextElected6==0), bw($band) n(100) mbw(100)" , norestore fast
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
quietly replace y=next21_23tulo
quietly parmby "rd y x if (elected==1 & nextElected5==1) | (elected==0 & nextElected5==0), bw($band) n(100) mbw(100)" , norestore fast
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
quietly parmby "rd y x if (elected==1 & nextElected4==1) | (elected==0 & nextElected4==0), bw($band) n(100) mbw(100)" , norestore fast
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
quietly parmby "rd y x if (elected==1 & nextElected3==1) | (elected==0 & nextElected3==0), bw($band) n(100) mbw(100) ", norestore fast
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
quietly parmby "rd y x if (elected==1 & nextElected2==1) | (elected==0 & nextElected2==0), bw($band) n(100) mbw(100) ", norestore fast
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
replace y=next1_3tulo
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
replace y=next1_3tulo
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
replace y=next1_3tulo
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


gen period=_n-2
drop if temp==1
drop temp
compress


graph twoway scatter estimate min95 max95 period if period>=1, ///
	msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) connect(l l l) ///
	lcolor(black gray gray) lpattern(solid dot dot) ///
|| , legend(off)  yscale(range(0 60000)) xlabel(1(1)7) b2title("", ring(4)) xtitle("Electoral period") ytitle(" €/year") ///
	saving("$GRAPHPATH/EarningsPeriodsByStatus.gph",replace)

graph export "$OUTPATH/Figure_parliament_earnings_periods_bystatus.pdf", fontface(Helvetica)  replace
graph export "$OUTPATH/Figure_parliament_earnings_periods_bystatus.eps", fontface(Helvetica)  replace

clear 
quietly capture erase $DATAPATH/parmlags_e.dta
use $DATAPATH/tempdata
