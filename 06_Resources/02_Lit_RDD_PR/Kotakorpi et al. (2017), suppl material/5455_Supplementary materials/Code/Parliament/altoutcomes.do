* RD for ALTERNATIVE OUTCOME VARIABLES (income measures)
preserve
clear
use $DATAPATH/esmall
quietly cap gen y=.
quietly cap gen diff3_3=next1_3tulo - last1_3tulo
quietly cap gen lognext1_3tulo=log(next1_3tulo)
quietly cap gen loglast1_3tulo=log(last1_3tulo)
quietly cap gen logdiff1_3tulo=lognext1_3tulo-loglast1_3tulo

/* remove outlier
quietly cap drop next1_3poB
quietly cap gen next1_3poB=next1_3potulo
quietly sum next1_3poB
quietly sum id if next1_3poB==r(max) 
replace next1_3poB=. if id==r(max)  
 */


eststo clear

local depvarlist diff3_3 lognext1_3tulo logdiff1_3tulo /* next1_3poB */

 
foreach depvar in `depvarlist' {
quietly replace y=`depvar'
quietly rdob_mod2 y x 
*global band=r(h_opt)
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
quietly rd y x, bw($band) n(100) mbw(100 50)
quietly eststo 
estadd scalar Bandwidth = round($band,.001) 
}

/* Display, save as smcl */
esttab  using $OUTPATH/AltOutcomes.smcl, /*page(tex)*/ replace type label nodepvar nonumber se  ///
 stats(N Bandwidth) ///
 mtitles("Difference" "Log level" "Log difference") ///
 title("Effect of getting elected: alternative outcome variables" ) 

restore
