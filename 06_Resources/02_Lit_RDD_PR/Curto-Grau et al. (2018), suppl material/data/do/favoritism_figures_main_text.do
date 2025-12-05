*******************************************
** REPLICATION FILES FOR
** Curto-Grau, Marta, Solé-Ollé, Albert, and Sorribas-Navarro, Pilar. 2018. "Does electoral competition curb party favoritism?" AEJ: Applied Economics
** Figures in the Main Text
*******************************************


cd "..."

clear all



use "db_main.dta", clear


su compreg1 if dca1==1
gen ecs1_bis = (compreg1 - r(mean)) if dca1==1

foreach i of numlist 2/15 {
su compreg1 if dca`i'==1
replace ecs1_bis = (compreg1 - r(mean)) if dca`i'==1
}


*******
** Figure 1 - Variation in seat margin
*******
 preserve
 drop ecs1_bis ecs1
 
 bysort codccaa t: keep if _n==1 //Keep one observation per regiona and year
 
 su compreg1 
 gen ecs1 = (compreg1 - r(mean))  
 


 
 su compreg1 if dca1==1
	gen ecs1_bis = (compreg1 - r(mean)) if dca1==1

	foreach i of numlist 2/15 {
	su compreg1 if dca`i'==1
	replace ecs1_bis = (compreg1 - r(mean)) if dca`i'==1
	}
	label var ecs1 "'between variation'"
	label var ecs1_bis "'within variation'"
	
 graph box ecs1_bis ecs1 , graphregion(color(white) fcolor(white))   ylabel(-20(5)32, grid glcol(gs14)) 
 
 restore
 
 

*******
** Figure 2 - Alignment vs. vote and seat margin
*******

preserve
//Generate bins with vote distance
	gen bin=1 if dist1>=0 & dist1<0.05
	replace bin=2 if dist1>=0.05 & dist1<0.1
	replace bin=3 if dist1>=0.1 & dist1<0.15
	replace bin=4 if dist1>=0.15 & dist1<0.2
	replace bin=5 if dist1>=0.2 & dist1<0.25
	replace bin=6 if dist1>=0.25 & dist1<0.3
	replace bin=7 if dist1>=0.3 & dist1<0.35
	replace bin=8 if dist1>=0.35 & dist1<0.4
	replace bin=9 if dist1>=0.4 & dist1<0.45
	replace bin=10 if dist1>=0.45 & dist1<0.5
	replace bin=11 if dist1>=0.5 & dist1<0.55
	replace bin=12 if dist1>=0.55 & dist1<0.6
	replace bin=13 if dist1>=0.60 & dist1<0.65
	replace bin=14 if dist1>=0.65 & dist1<0.70
	replace bin=21 if dist1<0 & dist1>=-0.05
	replace bin=22 if dist1<-0.05 & dist1>=-0.1
	replace bin=23 if dist1<-0.1 & dist1>=-0.15
	replace bin=24 if dist1<-0.15 & dist1>=-0.2
	replace bin=25 if dist1<-0.2 & dist1>=-0.25
	replace bin=26 if dist1<-0.25 & dist1>=-0.3
	replace bin=27 if dist1<-0.3 & dist1>=-0.35
	replace bin=28 if dist1<-0.35 & dist1>=-0.4
	replace bin=29 if dist1<-0.4 & dist1>=-0.45
	replace bin=30 if dist1<-0.45 & dist1>=-0.5
	replace bin=31 if dist1<-0.5 & dist1>=-0.55
	replace bin=32 if dist1<-0.55 & dist1>=-0.6
	replace bin=33 if dist1<-0.60 & dist1>=-0.65
	replace bin=34 if dist1<-0.65 & dist1>=-0.70

//Generate bins with seat distance
	gen bin2=1 if dif==1
	replace bin2=2 if dif==2
	replace bin2=3 if dif==3
	replace bin2=4 if dif==4
	replace bin2=5 if dif==5
	replace bin2=6 if dif==6
	replace bin2=7 if dif==7
	replace bin2=8 if dif==8
	replace bin2=9 if dif==9
	replace bin2=10 if dif==10
	replace bin2=11 if dif==11
	replace bin2=12 if dif==-1
	replace bin2=13 if dif==-2
	replace bin2=14 if dif==-3
	replace bin2=15 if dif==-4
	replace bin2=16 if dif==-5
	replace bin2=17 if dif==-6
	replace bin2=18 if dif==-7
	replace bin2=19 if dif==-8
	replace bin2=20 if dif==-9

// Average value in each bin
bysort bin: egen av_dist1=mean(dist1)
bysort bin2: egen av_dif=mean(dif)
bysort bin: egen av_alig_ab=mean(ab)
bysort bin2: egen av_alig_ab2=mean(ab)

keep if abs(dist1)<0.6 // Restrict bdw for better visualization

twoway (scatter av_alig_ab av_dist if dist1<0 & dist1>-0.6, msymbol(Oh) mcolor(black)) ///
 (scatter av_alig_ab av_dist if dist1>0 & dist1<0.6, msymbol(Oh) mcolor(black)) ///
 (lpolyci ab dist1 if  dist1<0, lpattern(shortdash) lcolor(black) bw(0.193) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  ab dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw(0.193) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly ab dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.193) degree(1) kernel(epa)   ) ///
 (lpoly  ab dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw(0.193) degree(1) kernel(epa)  ), ///
   xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(medlarge) margin(medsmall)) ///
  ytitle("Alignment Regional-Local", size(medlarge) margin(medsmall)) graphregion(fcolor(white))  legend(off)  ///
   ylabel(0(0.2)1, labsize(4) grid glcolor(gs14)) xlabel(-0.5 (0.5) 0.5, labsize(4)) yscale(range(0 1))  xscale(range(-0.6 0.6)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white))

twoway (scatter av_alig_ab2 av_dif if dif<0 & dif>=-9, msymbol(Oh) mcolor(black)) ///
 (scatter av_alig_ab2 av_dif if dif>0 & dif<=9, msymbol(Oh) mcolor(black)), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc seat margin", size(med) margin(medsmall)) ///
  ytitle("Alignment Regional-Local", size(med) margin(medsmall)) graphregion(fcolor(white))  legend(off)  ///
  ylabel(0(0.2)1, labsize(med)) xlabel(-9 (1) 9, labsize(med)) yscale(range(0 1))  xscale(range(-9 9)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white))

restore 




**********
** Figure 3 - Mccrary and histograms
**********
clear all

use "db_full.dta"

keep if dist1!=.

//McCrary
DCdensity dist1 , breakpoint(0) generate(Xj Yj r0 fhat se_fhat) 

//Histogram
twoway (histogram dist1 , start(-1) frequency width(0.1) fcolor(gs6)  lcolor(gs5)) (histogram dist1 , start(-1) frequency width(0.05) fcolor(gs9) lcolor(gs5)) (histogram dist1 , start(-1) frequency width(0.025) fcolor(gs12) lcolor(gs5)),  xline(0, lcolor(gray)) graphregion(color(white) lcolor(gs14)) xtitle("Regional incumbent's bloc vote margin", margin(medsmall)) ytitle("Frequency", margin(medsmall))





**********
** Figure 4 - Capital transfers vs. bloc margin
**********

use "db_main.dta", clear

preserve
//Generate bins with vote distance
	gen bin=1 if dist1>=0 & dist1<0.05
	replace bin=2 if dist1>=0.05 & dist1<0.1
	replace bin=3 if dist1>=0.1 & dist1<0.15
	replace bin=4 if dist1>=0.15 & dist1<0.2
	replace bin=5 if dist1>=0.2 & dist1<0.25
	replace bin=6 if dist1>=0.25 & dist1<0.3
	replace bin=7 if dist1>=0.3 & dist1<0.35
	replace bin=8 if dist1>=0.35 & dist1<0.4
	replace bin=9 if dist1>=0.4 & dist1<0.45
	replace bin=10 if dist1>=0.45 & dist1<0.5

	replace bin=21 if dist1<0 & dist1>=-0.05
	replace bin=22 if dist1<-0.05 & dist1>=-0.1
	replace bin=23 if dist1<-0.1 & dist1>=-0.15
	replace bin=24 if dist1<-0.15 & dist1>=-0.2
	replace bin=25 if dist1<-0.2 & dist1>=-0.25
	replace bin=26 if dist1<-0.25 & dist1>=-0.3
	replace bin=27 if dist1<-0.3 & dist1>=-0.35
	replace bin=28 if dist1<-0.35 & dist1>=-0.4
	replace bin=29 if dist1<-0.4 & dist1>=-0.45
	replace bin=30 if dist1<-0.45 & dist1>=-0.5


rdrobust tk dist1  , fuzzy(ab) vce(cluster codiine) kernel(epa) //Select optimal bandwidth

// Average value in each bin
bysort bin: egen av_dist1=mean(dist1)
bysort bin: egen av_alig_ab=mean(ab)
bysort bin: egen av_tk=mean(tk)

keep if abs(dist1)<0.5 // Restrict bdw for better visualization

twoway (scatter av_tk av_dist if dist1<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_tk av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci tk dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.213) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  tk dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw(0.213) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly tk dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.213) degree(1) kernel(epa)   ) ///
 (lpoly  tk dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw(0.213) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(medsmall)) ///
  ytitle("Regional transfers", size(vlarge) margin(medsmall)) graphregion(fcolor(white))  legend(off)  ///
  ylabel(80(20)180, labsize(vlarge) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(vlarge)) yscale(range(80 180))  xscale(range(-0.50 0.50)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white))


restore	



*****************
**  Figure 5 - Low vs High competition
*****************	
// Low competition
preserve

keep if above==1	//Margin of victory above average
//Generate bins with vote distance
	gen bin=1 if dist1>=0 & dist1<0.05
	replace bin=2 if dist1>=0.05 & dist1<0.1
	replace bin=3 if dist1>=0.1 & dist1<0.15
	replace bin=4 if dist1>=0.15 & dist1<0.2
	replace bin=5 if dist1>=0.2 & dist1<0.25
	replace bin=6 if dist1>=0.25 & dist1<0.3
	replace bin=7 if dist1>=0.3 & dist1<0.35
	replace bin=8 if dist1>=0.35 & dist1<0.4
	replace bin=9 if dist1>=0.4 & dist1<0.45
	replace bin=10 if dist1>=0.45 & dist1<0.5

	replace bin=21 if dist1<0 & dist1>=-0.05
	replace bin=22 if dist1<-0.05 & dist1>=-0.1
	replace bin=23 if dist1<-0.1 & dist1>=-0.15
	replace bin=24 if dist1<-0.15 & dist1>=-0.2
	replace bin=25 if dist1<-0.2 & dist1>=-0.25
	replace bin=26 if dist1<-0.25 & dist1>=-0.3
	replace bin=27 if dist1<-0.3 & dist1>=-0.35
	replace bin=28 if dist1<-0.35 & dist1>=-0.4
	replace bin=29 if dist1<-0.4 & dist1>=-0.45
	replace bin=30 if dist1<-0.45 & dist1>=-0.5


rdrobust tk dist1  , fuzzy(ab) vce(cluster codiine) kernel(epa)  //Select optimal bandwidth

// Average value in each bin
bysort bin: egen av_dist1=mean(dist1)
bysort bin: egen av_alig_ab=mean(ab)
bysort bin: egen av_tk=mean(tk)
bysort bin: egen av_ecs1=mean(ecs1)

keep if abs(dist1)<=0.50  // Restrict bdw for better visualization

twoway (scatter av_tk av_dist if dist1<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_tk av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci tk dist1 if  dist1<0, lpattern(shortdash) lcolor(black) bw(0.214) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  tk dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw(0.214) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly tk dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.214) degree(1) kernel(epa)  ) ///
 (lpoly  tk dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw(0.214) degree(1) kernel(epa) ), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(medsmall)) ///
  ytitle("Regional transfers", size(vlarge) margin(medsmall)   ) graphregion(fcolor(white))  legend(off)  ///
  ylabel(50(50)200, labsize(vlarge) grid glcolor(gs14)  ) xlabel(-0.50 (0.50) 0.50,  labsize(vlarge)) yscale(range(50 200) )  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) 

restore	 
				
				
// High competition
preserve

keep if above==0	//Margin of victory below average

	gen bin=1 if dist1>=0 & dist1<0.05
	replace bin=2 if dist1>=0.05 & dist1<0.1
	replace bin=3 if dist1>=0.1 & dist1<0.15
	replace bin=4 if dist1>=0.15 & dist1<0.2
	replace bin=5 if dist1>=0.2 & dist1<0.25
	replace bin=6 if dist1>=0.25 & dist1<0.3
	replace bin=7 if dist1>=0.3 & dist1<0.35
	replace bin=8 if dist1>=0.35 & dist1<0.4
	replace bin=9 if dist1>=0.4 & dist1<0.45
	replace bin=10 if dist1>=0.45 & dist1<0.5

	replace bin=21 if dist1<0 & dist1>=-0.05
	replace bin=22 if dist1<-0.05 & dist1>=-0.1
	replace bin=23 if dist1<-0.1 & dist1>=-0.15
	replace bin=24 if dist1<-0.15 & dist1>=-0.2
	replace bin=25 if dist1<-0.2 & dist1>=-0.25
	replace bin=26 if dist1<-0.25 & dist1>=-0.3
	replace bin=27 if dist1<-0.3 & dist1>=-0.35
	replace bin=28 if dist1<-0.35 & dist1>=-0.4
	replace bin=29 if dist1<-0.4 & dist1>=-0.45
	replace bin=30 if dist1<-0.45 & dist1>=-0.5

rdrobust tk dist1  , fuzzy(ab) vce(cluster codiine) kernel(epa)  //Select optimal bandwidth

// Average value in each bin 
bysort bin: egen av_dist1=mean(dist1)
bysort bin: egen av_alig_ab=mean(ab)
bysort bin: egen av_tk=mean(tk)
bysort bin: egen av_ecs1=mean(ecs1)

keep if abs(dist1)<=0.50 // Restrict bdw for better visualization

twoway (scatter av_tk av_dist if dist1<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_tk av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci tk dist1 if  dist1<0, lpattern(shortdash) lcolor(black) bw(0.228) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  tk dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw(0.228) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly tk dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.228) degree(1) kernel(epa)  ) ///
 (lpoly  tk dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw(0.228) degree(1) kernel(epa) ), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(medsmall)) ///
  ytitle("Regional transfers", size(vlarge) margin(medsmall)   ) graphregion(fcolor(white))  legend(off)  ///
  ylabel(50(50)200, labsize(vlarge) grid glcolor(gs14)  ) xlabel(-0.50 (0.50) 0.50,  labsize(vlarge)) yscale(range(50 200) )  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) 
restore	 

 
 
**********
** Figure 6 - Continuity of regional seat margin
**********
preserve
	gen bin=1 if dist1>=0 & dist1<0.05
	replace bin=2 if dist1>=0.05 & dist1<0.1
	replace bin=3 if dist1>=0.1 & dist1<0.15
	replace bin=4 if dist1>=0.15 & dist1<0.2
	replace bin=5 if dist1>=0.2 & dist1<0.25
	replace bin=6 if dist1>=0.25 & dist1<0.3
	replace bin=7 if dist1>=0.3 & dist1<0.35
	replace bin=8 if dist1>=0.35 & dist1<0.4
	replace bin=9 if dist1>=0.4 & dist1<0.45
	replace bin=10 if dist1>=0.45 & dist1<0.5

	replace bin=21 if dist1<0 & dist1>=-0.05
	replace bin=22 if dist1<-0.05 & dist1>=-0.1
	replace bin=23 if dist1<-0.1 & dist1>=-0.15
	replace bin=24 if dist1<-0.15 & dist1>=-0.2
	replace bin=25 if dist1<-0.2 & dist1>=-0.25
	replace bin=26 if dist1<-0.25 & dist1>=-0.3
	replace bin=27 if dist1<-0.3 & dist1>=-0.35
	replace bin=28 if dist1<-0.35 & dist1>=-0.4
	replace bin=29 if dist1<-0.4 & dist1>=-0.45
	replace bin=30 if dist1<-0.45 & dist1>=-0.5


// Average value in each bin
bysort bin: egen av_dist1=mean(dist1) 
bysort bin: egen av_ecs1=mean(ecs1)

rdrobust ecs1 dist1  , fuzzy(ab) vce(cluster codiine) kernel(epa) //Select optimal bandwidth

keep if abs(dist1)<=0.5 // Restrict bdw for better visualization

twoway (scatter av_ecs1 av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_ecs1  av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci ecs1  dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.1973) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  ecs1  dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.1973) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly ecs1  dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.1973) degree(1) kernel(epa)   ) ///
 (lpoly  ecs1  dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.1973) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(medsmall))  ///
  ytitle("Regional seat margin", size(vlarge) margin(medsmall)) graphregion(fcolor(white))  legend(off)  ///
  ylabel(-6 (2) 4.5, labsize(vlarge) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(vlarge))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white))

restore


	
**********
** Figure 7 - Marginal effects
**********
clear all

use "db_main.dta", clear

//Generate "within" regional competition variable
su compreg1 if dca1==1
gen ecs1_bis = (compreg1 - r(mean)) if dca1==1

foreach i of numlist 2/15 {
su compreg1 if dca`i'==1
replace ecs1_bis = (compreg1 - r(mean)) if dca`i'==1
}

ge aecs1_bis=ecs1_bis*ab 
ge decs1_bis=ecs1_bis*dab 

ge aecs1_bis_vsa=ecs1_bis*ab*dist1 
ge decs1_bisvsa=ecs1_bis*dab*dist1 
ge aecs1_bisvsa2=ecs1_bis*ab*dist2
ge decs1_bisvsa2=ecs1_bis*dab*dist2 

ge dist1_ecs1_bis=ecs1_bis*dist1
ge dist2_ecs1_bis=ecs1_bis*dist2

global controls2b "dist1_ecs1_bis dist2_ecs1_bis dist1 dist2 "
global regFE "i.codccaa"

foreach i of numlist 1/15{
su codiine if dca`i'==1  //To know how many observations in each region
}

ivregress 2sls tk  (dca_ab* aecs1_bis aecs1_bis_vsa aecs1_bisvsa2 dca_vsa* dca_2vsa* = dca_dab* decs1_bis decs1_bisvsa decs1_bisvsa2 dca_vda* dca_2vda*) $controls2b, vce(cluster codccaa)
predictnl marginal=_b[dca_ab1]*(1078/6050)+_b[dca_ab2]*(202/6050)+_b[dca_ab3]*(147/6050)+_b[dca_ab4]*(115/6050)+_b[dca_ab5]*(160/6050)+_b[dca_ab6]*(112/6050)+_b[dca_ab7]*(533/6050)+_b[dca_ab8]*(661/6050)+_b[dca_ab9]*(828/6050)+_b[dca_ab10]*(707/6050)+_b[dca_ab11]*(384/6050)+_b[dca_ab12]*(747/6050)+_b[dca_ab13]*(207/6050)+_b[dca_ab14]*(102/6050)+_b[dca_ab15]*(67/6050)+_b[aecs1_bis]*ecs1_bis  if e(sample), ci (conf1 conf2)
sort ecs1_bis
twoway (line marginal ecs1_bis if e(sample), yaxis(1) lc(black) lw(thin)) (line conf1 ecs1_bis if e(sample), yaxis(1) lc(black) lw(thin) lp(dash)) (line conf2 ecs1_bis if e(sample), yaxis(1) lc(black) lw(thin) lp(dash)) (kdensity ecs1_bis if e(sample), bwidth(4)  yaxis(2) lc(cranberry) lw(thin)),  yscale(range(-50 250)) yscale(range(0 0.15)) ylabel(-100(100)200, axis(1)) ylabel(0(0.01)0.06, axis(2)) xlabel(-12(4)12) xline(0, lp(dot) lc(black)) yline(0, lp(dot) lc(black)) legend(off) graphregion(fc(white)) ///
 ytitle("Marginal effect", margin(medsmall) axis(1) size(vlarge))  ytitle("Density", margin(medsmall) axis(2) size(vlarge))   xtitle("Regional seat margin", margin(medsmall) size(vlarge)) 


 ////////END
