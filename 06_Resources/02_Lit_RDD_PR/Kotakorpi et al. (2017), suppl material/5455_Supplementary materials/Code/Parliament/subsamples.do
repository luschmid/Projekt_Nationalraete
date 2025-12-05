* RD by SUBSAMPLES

quietly rdob_mod2 y x /* calculate the optimal Imbens-Kalyanaraman bandwidth for rd estimation*/
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
rd y x, bw($band) n(100) mbw(100)
display "N: " e(N) ", IK-optimal BW: " $ikband

eststo clear

* Binary split loop
foreach splitvar in female incumbent neverelected south {
*rdob_mod2 y x if `splitvar'
*global band=r(h_opt)
quietly rd y x if `splitvar', bw($band)  n(100) mbw(100 50)
eststo
estadd scalar Bandwidth = round($band,.001) 
quietly rd y x if !`splitvar', bw($band) n(100) mbw(100 50)
eststo
estadd scalar Bandwidth = round($band,.001) 
}
 
esttab using $OUTPATH/Subsamples_1.smcl, replace type label nonumber nodepvar se  ///
 stats(N Bandwidth)  ///
 title("Main effect by subsample  - split by binary variable" ) ///
 mtitles("Female" "Male" "Incumbent" "Non-incumbent" "Not previously elected" "Previously elected" "South" "North") 

 eststo clear
 
 * Party split loop 
foreach splitvar in KESK KOK SDP OTHER {
*rdob_mod2 y x if `splitvar'
*global band=r(h_opt)
quietly rd y x if `splitvar', bw($band)  n(100) mbw(100 50)
eststo
estadd scalar Bandwidth = round($band,.001) 
}

esttab using $OUTPATH/Subsamples_2.smcl, replace type label nonumber nodepvar se  ///
 stats(N Bandwidth)  ///
 title("Main effect by subsample  - split by major party" ) ///
 mtitles("Centre" "NCP" "SDP" "Other") 
 
eststo clear 

* Median split loop:
foreach splitvar in age past_avg year {
quietly sum `splitvar' if (`splitvar'!=. & y!=. & x!=.), detail
global CUTOFF = r(p50)
*rdob_mod2 y x if `splitvar' <= $CUTOFF
*global band=r(h_opt)
quietly rd y x if `splitvar' <= $CUTOFF, bw($band) n(100) mbw(100 50)
eststo
estadd scalar Bandwidth = round($band,.001) 
quietly sum `splitvar' if `splitvar'<= $CUTOFF & y!=. & x!=.
global CUTOFF0 = r(max)
estadd scalar Cutoff = round($CUTOFF0,.001) 
quietly rd y x if `splitvar' > $CUTOFF & `splitvar'!=., bw($band) n(100) mbw(100 50)
eststo
estadd scalar Bandwidth = round($band,.001) 
quietly sum `splitvar' if `splitvar'> $CUTOFF & y!=. & x!=.
global CUTOFF1 = r(min)
estadd scalar Cutoff = round($CUTOFF1,.001)
}

esttab using $OUTPATH/Subsamples_3.smcl, replace type label nonumber nodepvar se  ///
 stats(N Cutoff Bandwidth)  ///
 title("Main effect by subsample - split at median" ) ///
 mtitles("Young" "Old" "Low income" "High income" "Pre 2000" "Post 2000") 
