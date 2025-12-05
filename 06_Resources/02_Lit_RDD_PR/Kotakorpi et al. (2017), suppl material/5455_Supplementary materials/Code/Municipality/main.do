* Analyze data from MUNICIPAL elections
*

* PATH DEFINITIONS
global ROOTPATH E:/Temp/jleo/municipality /* root directory - customize this for local computer  */
*  
global DATAPATH $ROOTPATH/data   /*  Data directory.  Dropbox can be too slow - you may want to move data to a different drive */
global GRAPHPATH $ROOTPATH/gph  /* save Stata-graphs to this directory */
global OUTPATH $ROOTPATH/output  /* save output (log,tex,pdf etc) to this directory */

adopath + $ROOTPATH/ado/	 
capture mkdir $GRAPHPATH 	/* create subdirectory for Stata graphs */
capture mkdir $OUTPATH   	/* create subdirectory for output */

cd $ROOTPATH /* run do-files from current directory */

/* FILE DEPENDENCIES

This do-file calls and outputs, in this order: 
 make_ksmall.do 	-> ksmall.dta
 summary.do		-> Table I (part)
 earnings.do		-> Table VII
 altoutcomes.do	-> Table VIII
 covariates.do	-> Table IX
 balance.do		-> Figures
 subsamples.do	-> Table X
 incumbency.do	-> Table XI
 
ado-files:
  MainRDgraph.ado
  rdob_mod2.ado -- Calculates correct IK-bandwidth, by Miikka Rokkanen.
  DCdensity.ado -- McCrary-test, ado from http://emlab.berkeley.edu/~jmccrary/DCdensity/
 
DATA DEPENDENCIES

 kpanel.dta or ksmall.dta /* latter allows the skipping of make_ksmall.do */

*/

* LOAD DATA

/* starting from the main data set kpanel.dta. Can be slow!  */

/* clear
use $DATAPATH/kpanel
quietly do make_ksmall 
more
 */

/* ALTERNATIVE: start from a previously saved data created by make_ksmall.do  */
clear
use $DATAPATH/ksmall



* DEFINE the main OUTCOME variable 

global OUTCOME="next1_3tulo"
capture gen y = $OUTCOME
label var y "Average yearly income in e=1"

* the FORCING variable is x

* Create BINS by forcing variable (used only for graphs)
global DB=1 /* bin width */
egen bsbin=cut(x), at(-101($DB)101)
replace bsbin=bsbin+$DB/2
gen yxexist=(y!=. & x!=.)
bysort bsbin: egen N=count(yxexist)
drop yxexist

* DEFINE variables and constants

gen byte found = (avgtulo!=.)
label variable found "Register data found"

capture gen change_avg=future_avg-past_avg
capture gen int age2=age^2

char year[omit] 2008 /* Reference election year  */
global YEARS = "y96 y00 y04" 

set matsize 2000

set scheme s1color

* Examine the DISTRIBUTION of the forcing variable. Note large number of hopeless OTHER party candidates near x = -100

histogram x if x>=-99,width($DB) start(-100) freq xtitle("") ytitle("") title("Municipal elections") ///
 saving($GRAPHPATH/Histogram_muni.gph, replace)
count if x<-99
graph export "$OUTPATH/Histogram - Municipal.pdf", fontface(Helvetica) replace
graph export "$OUTPATH/Histogram - Municipal.eps", fontface(Helvetica) replace

DCdensity x, breakpoint(0) generate(Xj Yj r0 fhat se_fhat) nograph
quietly drop Xj Yj r0 fhat se_fhat
more

* MAIN ANALYSIS

do summary
  
do earnings


* Main RD Graph with fixed bandwidth (Figure 2)

global XMAX=50
global YTITLE="€/year"
global XTITLE="Electoral closeness"
global MARKERSIZE="small"

replace y= $OUTCOME
quietly rdob_mod2 y x
global band=r(h_opt)
*rd y x, bw($band) n(100) mbw(100)
global TITLE=""
global B2TITLE=""
global GPHSTRING="$GRAPHPATH/main.gph"
MainRDgraph
graph export "$OUTPATH/Main_RD_Graph_Municipality.pdf", fontface(Helvetica) replace
graph export "$OUTPATH/Main_RD_Graph_Municipality.eps", fontface(Helvetica) replace

* Bandwidth robustness

bwplot $OUTCOME x "" "Bandwidth" "Euros/year"
graph export "$OUTPATH/BandwidthRobustness - municipal.pdf", fontface(Helvetica)  replace
graph export "$OUTPATH/BandwidthRobustness - municipal.eps", fontface(Helvetica)  replace

*rdob_mod2 $OUTCOME x /* calculate the optimal Imbens-Kalyanaraman bandwidth for rd estimation*/
*global band=r(h_opt) /* save the optimal bandwidth*/
*quietly rd y x, bw($band) n(100) mbw(10(10)210) bdep ox
*graph export "$OUTPATH/BandwidthRobustness - municipal.pdf", as(pdf) fontface(Helvetica)  replace

do time

do altoutcomes
 
do covariates
 
do balance

do subsamples
 
do incumbency

do time

/*
use $DATAPATH/tempdata
capture replace y= $OUTCOME
do rd_poly
*/
