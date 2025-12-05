*save $DATAPATH/tempdata, replace
clear
use $DATAPATH/tempdata

* figures to examine outside earnings & earnings of marginal elected politician, pre & post 2000 

replace y=$OUTCOME
quietly rdob_mod2 y x /* calculate the optimal Imbens-Kalyanaraman bandwidth for rd estimation*/
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */

keep if year<1999
rd y x, bw($band) n(100) mbw(100) 
display e(N)
global TITLE="1991 & 1995 elections"
global GPHSTRING="$GRAPHPATH/next1_3pre2000.gph"
global B2TITLE=""
MainRDgraph

clear
use $DATAPATH/tempdata

replace y=$OUTCOME

keep if year>=1999


rd y x, bw($band) n(100) mbw(100) 
display e(N)
global TITLE="1999-2007 elections"
global GPHSTRING="$GRAPHPATH/next1_3post2000.gph"
global B2TITLE=""
keep if next1_3tulo<340000 /* remove outlier from figure */
MainRDgraph

graph combine $GRAPHPATH/next1_3pre2000.gph $GRAPHPATH/next1_3post2000.gph, ///
  b2title("Earnings in (t+1)-(t+3) before and after salary reform") xcommon ycommon
graph export "$OUTPATH/Salary_reform.pdf",  fontface(Helvetica) replace
graph export "$OUTPATH/Salary_reform.eps",  fontface(Helvetica) replace
*
