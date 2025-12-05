* RD WITH COVARIATES

replace y=$OUTCOME
quietly rdob_mod2 y x
global band=r(h_opt)

eststo clear


quietly rd y x, bw($band) covar($YEARS) mbw(100)
eststo 
quietly rd y x, bw($band) covar($YEARS nseats) mbw(100)
eststo 
quietly rd y x, bw($band) covar($YEARS nseats age age2 female incumbent) mbw(100)
eststo 
*rd y x, bw($band) covar($YEARS Dist* age age2 female incumbent last1_3tulo) mbw(100)
*rd y x, bw($band) n(100) mbw(100)  x(past_avg last1_3tulo age female incumbent KOK KESK SDP OTHER votes south Dist*) 

/* Display, save as smcl */
esttab  using $OUTPATH/RD_Covariates.smcl, /*page(tex)*/ replace type label nodepvar nonumber se  ///
 title("Effect on future income of getting elected: alternative outcome variables" ) ///
  mtitles("Years" "+Muni. size" "+Indiv.") ///
 addnote("IK-optimal bandwidth $band") 
