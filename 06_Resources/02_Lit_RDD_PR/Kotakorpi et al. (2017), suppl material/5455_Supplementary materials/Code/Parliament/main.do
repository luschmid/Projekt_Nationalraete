* Analyze data from PARLIAMENTARY elections
*

* PATH DEFINITIONS

global ROOTPATH "E:/temp/jleo/Parliament" /* root directory - customize this for local system  */
*
global DATAPATH $ROOTPATH/data  /* Data directory */
global GRAPHPATH $ROOTPATH/gph  /* save Stata-graphs here */
global OUTPATH $ROOTPATH/output /* save output (log,tex,pdf etc) here */

/* Marko varapolut */
*global ROOTPATH "C:/Dropbox/Elections/Parliament" 
*global DATAPATH "D:/Politicians/data"   

/* Kaisan polut:
global ROOTPATH "C:\Users\kailkot\Dropbox\Politicians (1)"
global DATAPATH "C:\Users\kailkot\Desktop\kaisan matskua\RD\Dbox"   
*/

adopath + $ROOTPATH/ado/
capture mkdir $GRAPHPATH 	/* create subdirectory for Stata graphs */
capture mkdir $OUTPATH   	/* create subdirectory for output */

cd $ROOTPATH /* run do-files from current directory */

/* FILE DEPENDENCIES

 mainP.do calls and outputs, in this order: 
 make_esmall.do -> esmall.dta
 summary.do		-> Table 1 (parliament columns)
 earnings.do	-> Table 2 
 time.do 		-> Figure 3
 time2.do 		-> Figure 6
 altoutcomes.do	-> Table A1
 capital_income.do -> Figure A3
 covariates.do	-> Table A2
 balance.do		-> Table A2, Figure A2, Figure 5
 bandwidths.do 	-> Figure 4
 subsamples.do	-> Table 3
 discussion.do 	-> Figure 7
 incumbency.do	-> Table 4 (parliament column)
 rd_poly.do 	-> Estimate in a footnote

main.do also directly outputs Figure 1 (top panel), Figure 2

 Some do-files produce also additional figures that are not included in the final manuscript. 
 
included ado-files:
  bwplot.ado 	-- 	plot estimate and confidence intervals as a function of bandwidth
  MainRDgraph.ado -- standard RD plot, but with circle size proportional to #obs in bin
  rdob_mod2.ado -- Calculates correct IK-bandwidth, by Miikka Rokkanen. (rd.ado has a known bug in calculating IK bandwidth)
  DCdensity.ado -- McCrary-test, from http://emlab.berkeley.edu/~jmccrary/DCdensity/ 
  
following ados called  (available via standard install): distinct, rd, outreg2.
 
 
DATA DEPENDENCIES

epanel.dta or esmall.dta /* latter allows the skipping of make_esmall.do */

*/



* LOAD DATA (Comment out one of the two alternatives)

/* ALTERNATIVE 1: Start from the main data set */
/*  
clear
use $DATAPATH/epanel
quietly do make_esmall 
more
*/

/* ALTERNATIVE 2: Start from a previously saved data created by make_esmall.do  */
clear
use $DATAPATH/esmall

* END of LOAD DATA


* DEFINE the main OUTCOME variable 

global OUTCOME="next1_3tulo"
capture gen y = $OUTCOME
label var y "Average yearly income in e=1"

* Note the FORCING variable is x

* Create BINS by forcing variable (used only for RD graphs)
global DB=1 /* bin width */
egen bsbin=cut(x), at(-101($DB)101)
replace bsbin=bsbin+$DB/2
gen yxexist = (y!=. & x!=.)
bysort bsbin: egen N=count(yxexist)
drop yxexist

* DEFINE variables and constants

quietly tab eldist, gen(Dist) /* Electoral district dummies */
char eldist[omit] 1 /* Helsinki is the reference district */
char year[omit] 2007 /* 2007 is the reference election year during the income data years */
global YEARS = "y03 y99 y95" /* pre-election earnings observed for 1995 election and later */

capture gen change_avg=future_avg-past_avg
capture gen int age2=age^2

gen byte found=(missing==0)
label variable found "Register data found"

gen byte prevelected=(neverelectedbefore==0)
label variable prevelected "Previously elected"

set scheme s1color


* Examine the DISTRIBUTION of the forcing variable. 

histogram x if x>=-99, width($DB) start(-100) freq xtitle("") ytitle("")  title("Parliamentary elections") ///
 saving($GRAPHPATH/Histogram_parliament.gph, replace)
count if x<-99
graph export "$OUTPATH/Histogram - Parliament.pdf", fontface(Helvetica) replace
graph export "$OUTPATH/Histogram - Parliament.eps", fontface(Helvetica) replace

DCdensity x, breakpoint(0) generate(Xj Yj r0 fhat se_fhat) nograph
quietly drop Xj Yj r0 fhat se_fhat



* SUMMARY statistics
do summary
distinct id

* Main Bandwidth choice
global forceband=. /* setting forceband=. results in using IK-optimal bandwidth. */

* MAIN RD estimates
do earnings


* Main RD Graph with fixed bandwidth (Figure 2)

global XMAX=50
global YTITLE="Euros/year"
global XTITLE="Electoral closeness"
global MARKERSIZE="vsmall"

replace y= $OUTCOME
quietly rdob_mod2 y x
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
global TITLE=""
global B2TITLE=""
global GPHSTRING="$GRAPHPATH/main.gph"
MainRDgraph
graph export "$OUTPATH/Main_RD_Graph.pdf",  fontface(Helvetica) replace
graph export "$OUTPATH/Main_RD_Graph.eps",  fontface(Helvetica) replace

* figure with estimates for different time horizons 
* and similar figure for those still in parliament at period e vs. those who are not 
do time
do time2

* alternative outcome variables
do altoutcomes

*closer look at capital income
do capital_income

* VALIDITY & ROBUSTNESS
do covariates 

* BALANCE of predetermined variables
do balance

* Bandwidth robustness
do bandwidths


* SUBSAMPLES: main specification after splitting the sample by covariates
capture replace y = $OUTCOME
do subsamples 

* Additional material for discussion and notes
do discussion

* Incumbency effect
do incumbency


* Polynomial specification
clear 
use $DATAPATH/esmall
capture gen y = $OUTCOME
do rd_poly

clear
window manage close graph
* END

/* This works only if Municipality/main.do has been executed previously: */
graph combine $GRAPHPATH/Histogram_parliament.gph ../Municipality/gph/Histogram_muni.gph, ///
  b2title("Electoral closeness") l2title("Frequency")
graph export "$OUTPATH/Combined_Histograms.pdf",  fontface(Helvetica) replace
graph export "$OUTPATH/Combined_Histograms.eps",  fontface(Helvetica) replace
