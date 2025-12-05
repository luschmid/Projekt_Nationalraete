*figure of main effect at different lags (electoral periods)

save $DATAPATH/tempdata.dta, replace
clear


* EARNINGS 

quietly capture erase $DATAPATH/parmlags_k.dta
set obs 1
gen temp=1
quietly save $DATAPATH/parmlags_k,replace

clear
use $DATAPATH/tempdata
foreach var of varlist next13_15tulo next9_11tulo next5_7tulo next1_3tulo last1_3tulo last5_7tulo {
clear 
use $DATAPATH/tempdata
quietly replace y=next1_3tulo
quietly quietly rdob_mod2 y x  
* global band=r(h_opt)
quietly replace y=`var'
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_k
quietly save $DATAPATH/parmlags_k,replace
}

gen period=_n-2
drop if temp==1
drop temp

compress


graph twoway scatter estimate min95 max95 period if period<1,  ///
 msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) ///
 connect(l l l) lcolor(black gray gray) lpattern(solid dot dot) ///
|| scatter estimate min95 max95 period if period>=1,  xlabel(none -1 "-2" 0 "-1" 1 2 3 4) ///
msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) connect(l l l) ///
lcolor(black gray gray) lpattern(solid dot dot) ///
|| , legend(off) b2title("", ring(4)) xtitle("Electoral period") ytitle(" €/year") ///
 saving("$GRAPHPATH/EarningsPeriods.gph",replace)

graph export "$OUTPATH/Figure_muni_earnings_periods.pdf", fontface(Helvetica)  replace
graph export "$OUTPATH/Figure_muni_earnings_periods.eps", fontface(Helvetica)  replace


clear
quietly capture erase $DATAPATH/parmlagk_e.dta


* CAPITAL INCOME 


quietly capture erase $DATAPATH/parmlags_k.dta
set obs 1
gen temp=1
quietly save $DATAPATH/parmlags_k,replace

clear
use $DATAPATH/tempdata
foreach var of varlist next13_15potulo next9_11potulo next5_7potulo next1_3potulo last1_3potulo last5_7potulo {
clear 
use $DATAPATH/tempdata
quietly replace y=next1_3potulo
quietly quietly rdob_mod2 y x  
* global band=r(h_opt)
quietly replace y=`var'
quietly parmby "rd y x, bw($band) n(100) mbw(100)", norestore fast
quietly gen nobs=e(N)
quietly append using $DATAPATH/parmlags_k
quietly save $DATAPATH/parmlags_k,replace
}

gen period=_n-2
drop if temp==1
drop temp

compress


graph twoway scatter estimate min95 max95 period if period<1,  ///
 msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) ///
 connect(l l l) lcolor(black gray gray) lpattern(solid dot dot) ///
|| scatter estimate min95 max95 period if period>=1,  xlabel(none -1 "-2" 0 "-1" 1 2 3 4) ///
msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) connect(l l l) ///
lcolor(black gray gray) lpattern(solid dot dot) ///
|| , legend(off) b2title("", ring(4)) xtitle("Electoral period") ytitle(" €/year") ///
 saving("$GRAPHPATH/EarningsPeriods.gph",replace)

graph export "$OUTPATH/Figure_muni_capincome_periods.pdf", fontface(Helvetica)  replace
graph export "$OUTPATH/Figure_muni_capincome_periods.eps", fontface(Helvetica)  replace

*
clear
use $DATAPATH/tempdata
