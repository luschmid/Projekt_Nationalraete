* RD by SUBSAMPLE

replace y=$OUTCOME
quietly rdob_mod2 y x /* calculate the optimal Imbens-Kalyanaraman bandwidth for rd estimation*/
global band=r(h_opt) /* save the optimal bandwidth */
*global band=10 /* force alternative bandwidth */

eststo clear

* Binary split loop
foreach splitvar in female incumbent {
*rdob_mod2 y x if `splitvar'
*global band=r(h_opt)
quietly rd y x if `splitvar', bw($band)  n(100) mbw(100)
eststo
estadd scalar Bandwidth = round($band,.001) 
quietly rd y x if !`splitvar', bw($band) n(100) mbw(100)
eststo
estadd scalar Bandwidth = round($band,.001) 
}
 
esttab using $OUTPATH/Subsamples_1.smcl, replace type label nonumber nodepvar se  ///
 stats(N Bandwidth)  ///
  title("Main effect by subsample  - split by binary variable" ) ///
 mtitles("Female" "Male" "Incumbent" "Non-incumbent") 

 eststo clear
 
 * Party split loop 
foreach splitvar in KESK KOK SDP OTHER {
*rdob_mod2 y x if `splitvar'
*global band=r(h_opt)
quietly rd y x if `splitvar', bw($band)  n(100) mbw(100)
eststo
estadd scalar Bandwidth = round($band,.001) 
}

esttab using $OUTPATH/Subsamples_2.smcl, replace type label nonumber nodepvar se  ///
 stats(N Bandwidth)  ///
  title("Main effect by subsample  - split by major party" ) ///
 mtitles("Centre" "NCP" "SDP" "Other") 
 
eststo clear 

* Median split loop:
foreach splitvar in nseats age past_avg year {
quietly sum `splitvar' if (`splitvar'!=. & y!=. & x!=.), detail
global CUTOFF = r(p50)
*rdob_mod2 y x if `splitvar' <= $CUTOFF
*global band=r(h_opt)
quietly rd y x if `splitvar' <= $CUTOFF, bw($band) n(100) mbw(100)
eststo
estadd scalar Bandwidth = round($band,.001) 
quietly sum `splitvar' if `splitvar'<= $CUTOFF & y!=. & x!=.
global CUTOFF0 = r(max)
estadd scalar Cutoff = round($CUTOFF0,.001) 
quietly rd y x if `splitvar' > $CUTOFF & `splitvar'!=., bw($band) n(100) mbw(100)
eststo
estadd scalar Bandwidth = round($band,.001) 
quietly sum `splitvar' if `splitvar'> $CUTOFF & y!=. & x!=.
global CUTOFF1 = r(min)
estadd scalar Cutoff = round($CUTOFF1,.001)
}

esttab using $OUTPATH/Subsamples_3.smcl, replace type label nonumber nodepvar se  ///
 stats(N Cutoff Bandwidth)  ///
 title("Main effect by subsample - split at median" ) ///
 mtitles("Small muni" "Large muni" "Young" "Old" "Low income" "High income" "Before 2000" "Since 2000") 
