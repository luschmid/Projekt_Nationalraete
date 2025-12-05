
 * RD for ALTERNATIVE OUTCOME VARIABLES (income measures)
 
replace y=$OUTCOME

cap gen diff3_3=next1_3tulo - last1_3tulo
cap gen lognext1_3tulo=log(next1_3tulo)
cap gen loglast1_3tulo=log(last1_3tulo)
cap gen logdiff1_3tulo=lognext1_3tulo-loglast1_3tulo

cap drop next1_3poB
cap gen next1_3poB=next1_3potulo
quietly sum next1_3poB
quietly sum id if next1_3poB==r(max) 
replace next1_3poB=. if id==r(max) /* remove outlier */

eststo clear
local depvarlist diff3_3 lognext1_3tulo logdiff1_3tulo next1_3poB
 
foreach depvar in `depvarlist' {
quietly replace y=`depvar'
quietly rdob_mod2 y x 
global band1=r(h_opt)
quietly rd y x, bw($band1) n(100) mbw(100)
eststo 
estadd scalar Bandwidth = round($band1,.001) 
}

/* Display, save as smcl */
esttab  using $OUTPATH/AltOutcomes.smcl, /*page(tex)*/ replace type label nodepvar nonumber se  ///
 stats(N Bandwidth) ///
 mtitles("Difference" "Log level" "Log difference" "Capital income") ///
 title("Effect of getting elected: alternative outcome variables" ) ///
 addnote("IK-optimal bandwidths")
