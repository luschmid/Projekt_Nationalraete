
global TITLE="Polynomial specification"
global GPHSTRING="$GRAPHPATH/temp.gph"

foreach depvar in y {

local bswidth 1    /* width of x bins */
local polyorder 4  /* order of control polynomial, 2 or more */
local controls " "
local bsmin -50-`bswidth'/2 /*  dropping the no-hoper population */
 
quietly drop if x<`bsmin'
quietly drop if x>-`bsmin'
quietly drop if `depvar'==.
quietly drop if x==.

forvalues ii = 2/`polyorder' {
quietly gen b`ii'=x^`ii'
}
local polyvars "b2-b`polyorder'"

/* Taking out the year means 
quietly reg `depvar' i.year
quietly predict depvarhat
quietly replace `depvar' = `depvar'-depvarhat
*/

* Regressing winners
quietly reg `depvar' x `polyvars' `controls' if elected, robust

matrix beta1 = e(b)
matrix V1 = e(V)
scalar alpha1 = beta1[1,`polyorder'+1]
scalar var1 = V1[`polyorder'+1,`polyorder'+1]
quietly predict yhat1 if e(sample)
quietly predict sigma1 if e(sample), stdp

* Regressing losers
quietly reg `depvar' x `polyvars' `controls' if !elected, robust

matrix beta0 = e(b)
matrix V0 = e(V)
scalar alpha0 = beta0[1,`polyorder'+1]
scalar var0 = V0[`polyorder'+1,`polyorder'+1]
quietly predict yhat0 if e(sample)
quietly predict sigma0 if e(sample), stdp

* Combining fitted values
quietly gen yhat=yhat1
quietly replace yhat = yhat0 if !elected
quietly gen sigma=sigma1
quietly replace sigma=sigma0 if !elected
quietly gen yupper = yhat + 1.96*sigma
quietly gen ylower = yhat - 1.96*sigma

* RD estimate
scalar pihat = alpha1-alpha0
scalar pisigma = sqrt(var0+var1)

* Creating bins
quietly gen  binmedian = `bswidth'*round((x+`bswidth'/2)/`bswidth') - (`bswidth'/2)
label variable binmedia "Bin midpoint"
* tab binmedian elected

quietly xi: reg `depvar' i.binmedian, noconstant
quietly predict eybin if e(sample)

* Main figure
sort x
twoway scatter yhat yupper ylower x if elected, connect(l l l) lcolor(black gray gray) lpattern(solid shortdash shortdash) ///
   mcolor(none none none) msize(tiny) ///
   || scatter yhat yupper ylower x if !elected, connect(l l l) lcolor(black gray gray) lpattern(solid shortdash shortdash) ///
   mcolor(none none none) msize(tiny) ///
   || scatter eybin binmedian if eybin>0, msize(vsmall) mcolor(red) /// 
   ,xtick(-50(5)50) xlabel(-50(10)50) t1title("`: var label `depvar''") ///
   legend(off) saving("$GRAPHPATH/rdpoly.gph",replace)  /* nodraw */
   
display "4th order polynomial RD"
display "`depvar' impact: " pihat " (" pisigma ")   N: " e(N)   
}

drop `polyvars' _*
drop yhat1-eybin
quietly compress

*histogram binmedian, discrete freq saving("$GRAPHPATH/BinHistogram_x.gph",replace)
