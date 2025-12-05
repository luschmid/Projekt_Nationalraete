* MAIN RD estimates 
eststo clear

quietly replace y=$OUTCOME
quietly rdob_mod2 y x /* calculate the optimal Imbens-Kalyanaraman bandwidth for rd estimation*/
global band=r(h_opt) /* save the optimal bandwidth*/
*global band=5 /* force a bandwidth */

foreach depvar in next1_3tulo next5_7tulo next9_11tulo /*next13_15tulo*/ future_avg {
quietly  replace y=`depvar'
quietly rd `depvar' x, bw($band)  n(100) mbw(100)
eststo
quietly sum `depvar' if abs(x)<= $band
estadd scalar N_in_BW = r(N) 
}

/* Display, save as smcl */
esttab  using $OUTPATH/Earnings.smcl, replace type label nodepvar nonumber se compress mtitles("1st" "2nd" "3rd" /*"4th"*/ "Future") ///
 title(Effect of getting elected on earnings by electoral period) ///
  stats(N N_in_BW)  ///
 addnote("Average yearly earnings in 2011 Euros" "IK optimal bandwidth $band")

 
 /* Naive income regressions 
reg y elected,robust
reg y elected $YEARS, robust
xi: areg y elected $YEARS age age2 female incumbent, robust absorb(eldist)
xi: areg y elected $YEARS age age2 female incumbent last1_3tulo, robust absorb(eldist)
*/
