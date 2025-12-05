* Robustness to covariates

cap drop y
cap gen y = $OUTCOME
quietly rdob_mod2 y x
global ikband=r(h_opt) /* IK-optimal bandwidth*/
global band=$ikband
quietly if ($forceband!=.) global band=$forceband  /* forcing a bandwidth */
eststo clear

quietly rd y x, bw($band) covar($YEARS) mbw(100 50)
quietly eststo 
quietly rd y x, bw($band) covar($YEARS Dist*) mbw(100 50)
quietly eststo 
quietly rd y x, bw($band) covar($YEARS Dist* age age2 female incumbent) mbw(100 50)
quietly eststo 

/* Display, save as smcl */
esttab  using $OUTPATH/RD_Covariates.smcl, /*page(tex)*/ replace type label nodepvar nonumber se  ///
 title("Effect on future income of getting elected: alternative outcome variables" ) ///
  mtitles("Years" "+District" "+Indiv.") ///
 addnote("Bandwidth $band")
