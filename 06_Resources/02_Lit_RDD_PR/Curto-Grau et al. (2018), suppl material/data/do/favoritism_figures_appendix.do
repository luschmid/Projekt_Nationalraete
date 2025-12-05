*******************************************
** REPLICATION FILES FOR
** Curto-Grau, Marta, Solé-Ollé, Albert, and Sorribas-Navarro, Pilar. 2018. "Does electoral competition curb party favoritism?" AEJ: Applied Economics
** Figures in the Appendix
*******************************************


cd "..."

clear all


*************
*Figure A1. Covariates’ discontinuity tests
*************

clear all

use "db_main.dta", clear


*preserve
// Generate bins
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


* Average values in each bin
bysort bin: egen av_dist1=mean(dist1) 

replace pob=pob/1000

*Population has some outliers (very large cities) that should be excluded
graph box pob
su pob, d
replace pob=. if pob>r(p95)

foreach y in debt tipo vcp  pob density pob_mes6591 pob_5_1491 extr_noeu91 unempl income localp educ2 presscirc regre regde meanden regterm ecs1 { 

qui{
bysort bin: egen av_`y'=mean(`y')

rdrobust `y' dist1  , fuzzy(ab) vce(cluster codiine) kernel(epa)  //Select optimal bandwidth
ge b_`y' = e(h_r)
	}
}


keep if abs(dist1)<=0.5 //Restrict sample for better visualization

//Debt burden
twoway (scatter av_debt av_dist if dist1<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_debt av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci debt dist1 if  dist1<0, lpattern(shortdash) lcolor(black) bw(0.287) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  debt dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw(0.287) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly debt dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.287) degree(1) kernel(epa)   ) ///
 (lpoly  debt dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw(0.287) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
   ylabel(0.00 (0.03) 0.12, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white))  title("Debt burden", size(huge) color(black)) legend(off) 

//Property tax rate
twoway (scatter av_tipo av_dist if dist1<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_tipo av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci tipo dist1 if  dist1<0, lpattern(shortdash) lcolor(black) bw(0.2548) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  tipo dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw(0.2548) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly tipo dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.2548) degree(1) kernel(epa)   ) ///
 (lpoly  tipo dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw(0.2548) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
   ylabel(0.4 (0.1) 0.8, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("Property tax rate", size(huge) color(black)) legend(off)  
	
//Property value
twoway (scatter av_vcp av_dist if dist1<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_vcp av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci vcp dist1 if  dist1<0, lpattern(shortdash) lcolor(black) bw(0.2509) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  vcp dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw(0.2509) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly vcp dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.2509) degree(1) kernel(epa)   ) ///
 (lpoly  vcp dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw(0.2509) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
    ylabel(15(5) 30, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("Property value", size(huge) color(black)) legend(off) 



//Population (thousands)
twoway (scatter av_pob av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_pob av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci pob dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.2226) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  pob dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.2226) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly pob dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.2226) degree(1) kernel(epa)   ) ///
 (lpoly  pob dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.2226) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
   ylabel(2 (2) 13, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("Population (thousands)", size(huge) color(black)) legend(off) 

	
//Population density	
twoway (scatter av_density av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_density av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci density dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.1927) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  density dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.1927) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly density dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.1927) degree(1) kernel(epa)   ) ///
 (lpoly  density dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.1927) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
    ylabel(100 (200) 1200, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("Population density", size(huge) color(black)) legend(off) 


//% Old	
twoway (scatter av_pob_mes6591 av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_pob_mes6591 av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci pob_mes6591 dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.203) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  pob_mes6591 dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.203) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly pob_mes6591 dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.203) degree(1) kernel(epa)   ) ///
 (lpoly  pob_mes6591 dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.203) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
   ylabel(0.12 (0.02) 0.22, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("% Old", size(huge) color(black)) legend(off) 

//% Young
twoway (scatter av_pob_5_1491 av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_pob_5_1491 av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci pob_5_1491 dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.2140) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  pob_5_1491 dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.2140) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly pob_5_1491 dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.2140) degree(1) kernel(epa)   ) ///
 (lpoly  pob_5_1491 dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.2140) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
   ylabel(0.18 (0.02) 0.26, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("% Young", size(huge) color(black)) legend(off) 

//% Immigrant
twoway (scatter av_extr_noeu91 av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_extr_noeu91 av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci extr_noeu91 dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.2202) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  extr_noeu91 dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.2202) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly extr_noeu91 dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.2202) degree(1) kernel(epa)   ) ///
 (lpoly  extr_noeu91 dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.2202) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
   ylabel(0.00 (0.01) 0.035, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) title("% Immigrant", size(huge) color(black)) legend(off) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white))
	
//% Unemployed
twoway (scatter av_unempl av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_unempl av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci unempl dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.2624) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  unempl dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.2624) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly unempl dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.2624) degree(1) kernel(epa)   ) ///
 (lpoly  unempl dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.2624) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
   ylabel(5 (0.5) 7.2, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("% Unemployed", size(huge) color(black)) legend(off) 

//Income indicator
twoway (scatter av_income av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_income  av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci income  dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.232) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  income  dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.232) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly income  dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.232) degree(1) kernel(epa)   ) ///
 (lpoly  income  dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.232) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
   ylabel(0.8 (0.1) 1.1, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("Income indicator", size(huge) color(black)) legend(off) 
		
//% Educated
twoway (scatter av_educ2 av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_educ2  av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci educ2  dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.229) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  educ2  dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.229) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly educ2  dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.229) degree(1) kernel(epa)   ) ///
 (lpoly  educ2  dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.229) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
    ylabel(-0.5 (0.25) 1, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("% Educated", size(huge) color(black)) legend(off) 

//Press circulation
twoway (scatter av_presscirc av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_presscirc  av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci presscirc  dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.256) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  presscirc  dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.2560) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly presscirc  dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.2560) degree(1) kernel(epa)   ) ///
 (lpoly  presscirc  dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.2560) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  ///
   ylabel(-20 (10) 30, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("Press circulation", size(huge) color(black)) legend(off) 

//Regional revenues per capita
twoway (scatter av_regre av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_regre  av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci regre  dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.2173) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  regre  dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.2173) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly regre  dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.2173) degree(1) kernel(epa)   ) ///
 (lpoly  regre  dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.2173) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray)) xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(small)) ///
   ylabel(-300 (150) 300, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("Regional revenues per capita", size(huge) color(black)) legend(off) 
	
//Regional debt
twoway (scatter av_regde av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_regde  av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci regde  dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.204) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  regde  dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.204) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly regde  dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.204) degree(1) kernel(epa)   ) ///
 (lpoly  regde  dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.204) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray)) xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(small)) ///
    ylabel(-2 (1) 4.4, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("Regional debt", size(huge) color(black)) legend(off) 

//Municipal density
twoway (scatter av_meanden av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_meanden  av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci meanden  dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.2525) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  meanden  dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.2525) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly meanden  dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.2525) degree(1) kernel(epa)   ) ///
 (lpoly  meanden  dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.2525) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray)) xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(small))  ///
   ylabel(-150 (100) 290, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("Municipal density", size(huge) color(black)) legend(off) 

			
//Tenure in office		
twoway (scatter av_regterm av_dist if dist1<0  , msymbol(Oh) mcolor(black)) ///
 (scatter av_regterm  av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci regterm  dist1 if  dist1<0 , lpattern(shortdash) lcolor(black) bw(0.2396) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  regterm  dist1 if  dist1>0 , lpattern(shortdash)  lcolor(black) bw(0.2396) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly regterm  dist1 if  dist1<0 , lpattern(solid) lcolor(black) bw(0.2396) degree(1) kernel(epa)   ) ///
 (lpoly  regterm  dist1 if  dist1>0 , lpattern(solid)  lcolor(black) bw(0.2396) degree(1) kernel(epa)  ), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(small))  ///
   ylabel(0 (0.1) 0.75, labsize(6) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(6))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) title("Tenure in office", size(huge) color(black)) legend(off) 
	
				

	

*************
*Figure A2. Covariates’ discontinuity tests
*************
clear all

use "db_full.dta", clear

preserve

keep if dist_alt!=. & tk!=.

//Generate bins
gen bin=1 if dist_alt>=0 & dist_alt<0.05
replace bin=2 if dist_alt>=0.05 & dist_alt<0.1
replace bin=3 if dist_alt>=0.1 & dist_alt<0.15
replace bin=4 if dist_alt>=0.15 & dist_alt<0.2
replace bin=5 if dist_alt>=0.2 & dist_alt<0.25
replace bin=6 if dist_alt>=0.25 & dist_alt<0.3
replace bin=7 if dist_alt>=0.3 & dist_alt<0.35
replace bin=8 if dist_alt>=0.35 & dist_alt<0.4
replace bin=9 if dist_alt>=0.4 & dist_alt<0.45
replace bin=10 if dist_alt>=0.45 & dist_alt<0.5
replace bin=11 if dist_alt>=0.5 & dist_alt<0.55
replace bin=12 if dist_alt>=0.55 & dist_alt<0.6
replace bin=13 if dist_alt>=0.60 & dist_alt<0.65
replace bin=14 if dist_alt>=0.65 & dist_alt<0.70
replace bin=21 if dist_alt<0 & dist_alt>=-0.05
replace bin=22 if dist_alt<-0.05 & dist_alt>=-0.1
replace bin=23 if dist_alt<-0.1 & dist_alt>=-0.15
replace bin=24 if dist_alt<-0.15 & dist_alt>=-0.2
replace bin=25 if dist_alt<-0.2 & dist_alt>=-0.25
replace bin=26 if dist_alt<-0.25 & dist_alt>=-0.3
replace bin=27 if dist_alt<-0.3 & dist_alt>=-0.35
replace bin=28 if dist_alt<-0.35 & dist_alt>=-0.4
replace bin=29 if dist_alt<-0.4 & dist_alt>=-0.45
replace bin=30 if dist_alt<-0.45 & dist_alt>=-0.5
replace bin=31 if dist_alt<-0.5 & dist_alt>=-0.55
replace bin=32 if dist_alt<-0.55 & dist_alt>=-0.6
replace bin=33 if dist_alt<-0.60 & dist_alt>=-0.65
replace bin=34 if dist_alt<-0.65 & dist_alt>=-0.70


* Mean values in each bin
bysort bin: egen av_dist=mean(dist_alt)
bysort bin: egen av_alig_ab=mean(ab)
bysort bin: egen av_tk=mean(tk)

rdrobust tk dist_alt  , fuzzy(ab) vce(cluster codiine) kernel(epa) //Select optimal bandwith

keep if abs(dist_alt)<=0.50 //Restrict for better visualization


twoway (scatter av_tk av_dist if dist_alt<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_tk av_dist if dist_alt>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci tk dist_alt if  dist_alt<0, lpattern(shortdash) lcolor(black) bw(0.199) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  tk dist_alt if  dist_alt>0, lpattern(shortdash)  lcolor(black) bw(0.199) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly tk dist_alt if  dist_alt<0, lpattern(solid) lcolor(black) bw(0.199) degree(1) kernel(epa)  ) ///
 (lpoly  tk dist_alt if  dist_alt>0, lpattern(solid)  lcolor(black) bw(0.199) degree(1) kernel(epa) ), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(medlarge) margin(medsmall)) ///
  ytitle("Regional transfers", size(medlarge) margin(medsmall)) graphregion(fcolor(white))  legend(off)  ///
  ylabel(80(20)220, labsize(4) grid glcolor(gs14)) xlabel(-0.50 (0.50) 0.50, labsize(4)) yscale(range(80 220))  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) 

restore	 


*************
*Figure A3. See R file
*************
//Dataset to create plots
clear all

use "db_main.dta", clear

ge dist100 = dist1*100

matrix B = J(30, 4, .)
	
	quietly forvalues j=4(2)30 {
		rdrobust tk dist100 , fuzzy(ab) vce(cluster codiine) h(`j') kernel(epa)
		matrix B[`j', 1] =e(tau_cl)
		matrix B[`j', 2] =e(tau_cl)- 1.96 * e(se_tau_cl)
		matrix B[`j', 3] =e(tau_cl)+ 1.96 * e(se_tau_cl)
		matrix B[`j', 4] =`j'
		
	}
	
	preserve
	svmat B
	keep B1-B4
	keep if B1 != .
	saveold ".../for_LATEplot_CI_epa.dta", version(12) replace
	restore	

	
clear all

use "db_main.dta", clear

ge dist100 = dist1*100

matrix B = J(30, 4, .)
	
	quietly forvalues j=4(2)30 {
		rdrobust tk dist100 , fuzzy(ab) vce(cluster codiine) h(`j') kernel(tri)
		matrix B[`j', 1] =e(tau_cl)
		matrix B[`j', 2] =e(tau_cl)- 1.96 * e(se_tau_cl)
		matrix B[`j', 3] =e(tau_cl)+ 1.96 * e(se_tau_cl)
		matrix B[`j', 4] =`j'
		
	}
	
	preserve
	svmat B
	keep B1-B4
	keep if B1 != .
	saveold ".../for_LATEplot_CI_tri.dta", version(12) replace
	restore	



*************
*Figure A4. McCrary and Histogram of alternative forcing variable
*************
clear all

use "db_full.dta", clear

preserve
keep if dist_alt!=.

DCdensity dist_alt , breakpoint(0) generate(Xj Yj r0 fhat se_fhat)

twoway (histogram dist_alt , start(-1) frequency width(0.1) fcolor(gs6)  lcolor(gs5)) (histogram dist_alt , start(-1) frequency width(0.05) fcolor(gs9) lcolor(gs5)) (histogram dist_alt , start(-1) frequency width(0.025) fcolor(gs12) lcolor(gs5)),  xline(0, lcolor(gray)) graphregion(color(white) lcolor(gs14)) xtitle("Regional incumbent's bloc vote margin (alternative measure)", margin(medsmall)) ytitle("Frequency", margin(medsmall))


restore



*************
*Figure A5. Social spending
*************
clear all

use "db_main.dta", clear

drop if f34142_c6_c6==.

//Select optimal bandwidths
rdrobust f34142_c6_c6 dist1  if above==0 & left_auto==1, fuzzy(ab) vce(cluster codiine) kernel(epa) 
rdrobust f34142_c6_c6 dist1  if above==0 & left_auto==0, fuzzy(ab) vce(cluster codiine) kernel(epa) 
rdrobust f34142_c6_c6 dist1  if above==1 & left_auto==1, fuzzy(ab) vce(cluster codiine) kernel(epa) 
rdrobust f34142_c6_c6 dist1  if above==1 & left_auto==0, fuzzy(ab) vce(cluster codiine) kernel(epa) 


// Low competition & left
preserve

keep if above==1 & left_auto==1

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

* Mean values in each bin
bysort bin: egen av_dist1=mean(dist1)
bysort bin: egen av_alig_ab=mean(ab)
bysort bin: egen av_f34142_c6_c6=mean(f34142_c6_c6)

keep if abs(dist1)<=0.50

twoway (scatter av_f34142_c6_c6 av_dist if dist1<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_f34142_c6_c6 av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci f34142_c6_c6 dist1 if  dist1<0, lpattern(shortdash) lcolor(black) bw( 0.281) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  f34142_c6_c6 dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw( 0.281) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly f34142_c6_c6 dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.281) degree(1) kernel(epa)  ) ///
 (lpoly  f34142_c6_c6 dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw( 0.281) degree(1) kernel(epa) ), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(medsmall)) ///
  ytitle("% Social spending", size(vlarge) margin(medsmall)   ) graphregion(fcolor(white))  legend(off)  ///
  ylabel(5(5)20, labsize(vlarge) grid glcolor(gs14)  ) xlabel(-0.50 (0.50) 0.50,  labsize(vlarge)) yscale(range(5 20) )  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) 

restore	 	


// High competition & left
preserve

keep if above==0 & left_auto==1

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

* Mean values in each bin
bysort bin: egen av_dist1=mean(dist1)
bysort bin: egen av_alig_ab=mean(ab)
bysort bin: egen av_f34142_c6_c6=mean(f34142_c6_c6)
bysort bin: egen av_ecs1=mean(ecs1)

keep if abs(dist1)<=0.50

twoway (scatter av_f34142_c6_c6 av_dist if dist1<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_f34142_c6_c6 av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci f34142_c6_c6 dist1 if  dist1<0, lpattern(shortdash) lcolor(black) bw( 0.186) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  f34142_c6_c6 dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw( 0.186) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly f34142_c6_c6 dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.186) degree(1) kernel(epa)  ) ///
 (lpoly  f34142_c6_c6 dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw( 0.186) degree(1) kernel(epa) ), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(medsmall)) ///
  ytitle("% Social spending", size(vlarge) margin(medsmall)   ) graphregion(fcolor(white))  legend(off)  ///
  ylabel(0(5)20, labsize(vlarge) grid glcolor(gs14)  ) xlabel(-0.50 (0.50) 0.50,  labsize(vlarge)) yscale(range(0 15) )  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) 

restore	 	


// Low competition & right
preserve

keep if above==1 & left_auto==0

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

* Mean values in each bin
bysort bin: egen av_dist1=mean(dist1)
bysort bin: egen av_alig_ab=mean(ab)
bysort bin: egen av_f34142_c6_c6=mean(f34142_c6_c6)


keep if abs(dist1)<=0.50

twoway (scatter av_f34142_c6_c6 av_dist if dist1<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_f34142_c6_c6 av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci f34142_c6_c6 dist1 if  dist1<0, lpattern(shortdash) lcolor(black) bw( 0.151) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci  f34142_c6_c6 dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw(0.151) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly f34142_c6_c6 dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.151) degree(1) kernel(epa)  ) ///
 (lpoly  f34142_c6_c6 dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw( 0.151) degree(1) kernel(epa) ), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(medsmall)) ///
  ytitle("% Social spending", size(vlarge) margin(medsmall)   ) graphregion(fcolor(white))  legend(off)  ///
  ylabel(5(5)20, labsize(vlarge) grid glcolor(gs14)  ) xlabel(-0.50 (0.50) 0.50,  labsize(vlarge)) yscale(range(5 20) )  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) 

restore	 		


// High competition & right
preserve

keep if above==0 & left_auto==0

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

* Mean values in each bin
bysort bin: egen av_dist1=mean(dist1)
bysort bin: egen av_alig_ab=mean(ab)
bysort bin: egen av_f34142_c6_c6=mean(f34142_c6_c6)

keep if abs(dist1)<=0.50

twoway (scatter av_f34142_c6_c6 av_dist if dist1<0 , msymbol(Oh) mcolor(black)) ///
 (scatter av_f34142_c6_c6 av_dist if dist1>0 , msymbol(Oh) mcolor(black)) ///
 (lpolyci f34142_c6_c6 dist1 if  dist1<0, lpattern(shortdash) lcolor(black) bw( 0.209) degree(1) kernel(epa)  ciplot(rline) ) ///
 (lpolyci f34142_c6_c6 dist1 if  dist1>0, lpattern(shortdash)  lcolor(black) bw(0.209) degree(1) kernel(epa)  ciplot(rline)) ///
 (lpoly f34142_c6_c6 dist1 if  dist1<0, lpattern(solid) lcolor(black) bw(0.209) degree(1) kernel(epa)  ) ///
 (lpoly  f34142_c6_c6 dist1 if  dist1>0, lpattern(solid)  lcolor(black) bw( 0.209) degree(1) kernel(epa) ), ///
  xline(0, lcolor(gray))  xtitle("Regional incumbent's bloc vote margin", size(vlarge) margin(medsmall)) ///
  ytitle("% Social spending", size(vlarge) margin(medsmall)   ) graphregion(fcolor(white))  legend(off)  ///
  ylabel(5(5)20, labsize(vlarge) grid glcolor(gs14)  ) xlabel(-0.50 (0.50) 0.50,  labsize(vlarge)) yscale(range(5 20) )  xscale(range(-0.5 0.5)) ///
	graphregion(color(white) fcolor(white)) plotregion(color(white)) 

restore	 




