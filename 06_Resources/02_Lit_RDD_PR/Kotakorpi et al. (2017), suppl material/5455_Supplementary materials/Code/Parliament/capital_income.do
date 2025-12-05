*Impact on capital income

quietly sum next1_3potulo
quietly sum id if next1_3potulo==r(max)
global outlierID=r(max)

* RD by electoral period

*average capital income in (t+1)-(t+3)
replace y=next1_3potulo
replace y=. if id==$outlierID
rdob_mod2 y x /* calculate the optimal Imbens-Kalyanaraman bandwidth for rd estimation*/
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */

rd y x, bw($band) n(100) mbw(100)
display e(N)
global TITLE="Average capital income in e=1"
global GPHSTRING="$GRAPHPATH/next1_3potulo.gph"
global B2TITLE=""
MainRDgraph

*(t+5)-(t+7)

replace y=next5_7potulo
replace y=. if id==$outlierID
rd y x, bw($band) n(100) mbw(100)
display e(N)
global TITLE="Average capital income in e=2"
global GPHSTRING="$GRAPHPATH/next5_7potulo.gph"
global B2TITLE=""
MainRDgraph

*(t+9)-(t+11)

replace y=next9_11potulo
replace y=. if id==$outlierID
rd y x, bw($band) n(100) mbw(100)
display e(N)
global TITLE="Average capital income in e=3"
global GPHSTRING="$GRAPHPATH/next9_11potulo.gph"
global B2TITLE=""
MainRDgraph


replace y=next13_15potulo
replace y=. if id==$outlierID
rd y x, bw($band) n(100) mbw(100)
display e(N)
global TITLE="Average capital income in e=4"
global GPHSTRING="$GRAPHPATH/next13_15potulo.gph"
global B2TITLE=""
MainRDgraph

* average capital income after the election

replace y=future_avgCap
replace y=. if id==$outlierID
rd y x, bw($band) n(100) mbw(100)
display e(N)


* bwplot next1_3pot x "Capital income in e=1" Bandwidth" ""

* TIME SERIES graph
*figure of effect at different lags (electoral periods)

save $DATAPATH/tempdata.dta, replace

quietly replace y=next1_3potulo
quietly rdob_mod2 y x if id!=$outlierID
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */

clear
quietly capture erase $DATAPATH/parmlags_e.dta
set obs 1
gen temp=1
save $DATAPATH/parmlags_e,replace

clear 
use $DATAPATH/tempdata

foreach var of varlist next25_27potulo next21_23potulo next17_19potulo next13_15potulo next9_11potulo next5_7potulo next1_3potulo last1_3potulo last5_7potulo {
clear 
use $DATAPATH/tempdata
replace y=`var'
quietly parmby "rd y x if id!=$outlierID, bw($band) n(100) mbw(100)", norestore fast
gen nobs=e(N)
quietly append using $DATAPATH/parmlags_e
save $DATAPATH/parmlags_e,replace
}

gen period=_n-2
drop if temp==1
drop temp
compress

twoway scatter nobs period, connect(l) title("Observations") xtitle("Period") ///
 saving("$GRAPHPATH/lagNCapincomeSeries.gph",replace) 


graph twoway scatter estimate min95 max95 period if period<=0, ///
		msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) connect(l l l) ///
		lcolor(black gray gray) lpattern(solid dot dot) ///
|| scatter estimate min95 max95 period if period>0, ///
		msize(small tiny tiny) mcolor(black gray gray) mfcolor(white white white) ///
		connect(l l l) lcolor(black gray gray) lpattern(solid dot dot) ///
|| , legend(off) xlabel(none -1 "-2" 0 "-1" 1 2 3 4 5 6 7) b2title("", ring(4)) xtitle("Electoral period") ytitle("€/year") ///
   saving("$GRAPHPATH/CapIncomePeriods.gph",replace)

graph export "$OUTPATH/Figure parliament capincome periods.pdf", fontface(Helvetica)  replace
graph export "$OUTPATH/Figure parliament capincome periods.eps", fontface(Helvetica)  replace

clear
quietly capture erase $DATAPATH/parmlags_e.dta
use $DATAPATH/tempdata




