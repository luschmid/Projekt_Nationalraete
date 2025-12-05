* MAIN RD estimates
eststo clear

quietly replace y=$OUTCOME
quietly rdob_mod2 y x /* calculate the optimal Imbens-Kalyanaraman bandwidth for rd estimation*/
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */

quietly sum next1_3potulo
quietly sum id if next1_3potulo==r(max)
global outlierID=r(max)
capture gen next1_3potuloB=next1_3potulo
replace next1_3potuloB= . if id==$outlierID

foreach depvar in next1_3tulo next5_7tulo next9_11tulo future_avg next1_3potuloB {
quietly  replace y=`depvar'
quietly rd `depvar' x, bw($band)  n(100) mbw(100 50) 
quietly eststo
estadd scalar Bandwidth = $band 
quietly sum `depvar' if abs(x)<= $band
estadd scalar N_in_BW = r(N) 
}	

/* Display, save as smcl */
esttab  using $OUTPATH/Earnings.smcl, replace type label nodepvar nonumber se compress mtitles("1st" "2nd" "3rd" "Future" "CapInc") ///
 title(Effect of getting elected on earnings by electoral period)  /// 
 addnote("Average yearly earnings in 2011 Euros" "IK optimal bandwidth $ikband")
 
 /* Save as tex fragment 
esttab using $OUTPATH/Earnings.tex, page[(tex)]  replace  label nodepvar nonumber se compress mtitles("1st" "2nd" "3rd" "Future") */


/* Naive income regressions
reg y elected,robust
reg y elected $YEARS, robust
xi: areg y elected $YEARS age age2 female incumbent, robust absorb(eldist)
xi: areg y elected $YEARS age age2 female incumbent last1_3tulo, robust absorb(eldist)
*/
