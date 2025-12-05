******************************************************************************************************************************************************
********************* EUROPE SURVEY ******************************************************************************************************************
******************************************************************************************************************************************************
use dta/620_CCS_Data_Wave1_v4.0.dta, clear 

****************************** District magnitude *************************************************************
replace t10=1 if t1==23 /* UK --> SMDs */
replace t6=1 if t1==8 /* one district in netherlands */
replace t6=1 if t1==23 /* set all UK obs. to be one */
replace t6=1 if t1==9 /* set all CAN obs. to be one */

gen M=t10*1 /* district magnitude */
egen min_M=min(M), by(t1 t3 t6)
egen avg_M=mean(M), by(t1 t3 t6)
egen max_M=max(M), by(t1 t3 t6)
drop if min_M!=max_M  /* 74 observations for Belgium */
gen obs=1
***************************************************************************************************************

************************************** Dropping the following countries ***************************************
drop if t10==.a /* missing district magnitude, e.g. Sweden, question not asked */
drop if t1==3 /* Germany, mixed system */
drop if t1==22 /* Malta ; district magnitude coded as 5.3? */
drop if t1==11 /* Estonia , no answers */
drop if t1==14 /* Austria , no answers */
drop if t1==18 /* Czech Republic doesn't answer b5a1 b5a2 */
drop if t1==1|t1==9|t1==20 /* drop countries outside Europe: Australia, Canada, and NZ */
***************************************************************************************************************

************************************** Dropping candidates running for upper house ****************************
drop if t9==3 /* "both" */
drop if t9==2 /* "lower house/Senate" */
***************************************************************************************************************

***************************************************************************************************************
***Candidates were asked how many hours they spent on various campaign activities per week during the last month before //
***of the election: (i) 0 hours, (ii) 1-5 hours (we code this as 3 hours), (iii) 5-10 hours (coded as 8 hours), //
***(iv) 10-15 hours (coded as 13 hours), (v) 15-20 hours (coded as 18 hours), and (vi) ``more than that'' (coded as 25 hours). 
*** Some countries collapse category (iv) and (v) (coded as 15 hours). 

foreach var in b5a1 b5a2 b5a3 b5a4 b5a5 b5a6 b5a7 b5a8 b5a9 b5a10 b5b1 b5b2 b5b3 b5b4 b5b5 b5b6 b5b7 {
gen h`var'=0
replace h`var'=3 if `var'==2  /* 1-5 hours */
replace h`var'=8 if `var'==3  /* 5-10 hours */
replace h`var'=13 if `var'==4  /* 10-15 hours */
replace h`var'=15 if `var'==4  & t1==2 /* switzerland 10-20 hours */
replace h`var'=15 if `var'==4  & t1==4 /* ireland 10-20 hours */
replace h`var'=15 if `var'==4  & t1==5 /* greece 10-20 hours */
replace h`var'=15 if `var'==4  & t1==6 /* finland 10-20 hours */
replace h`var'=15 if `var'==4  & t1==9 /* canada 10-20 hours */
replace h`var'=15 if `var'==4  & t1==10 /* portugal 10-20 hours */
replace h`var'=15 if `var'==4  & t1==12 /* iceland 11-20 hours */
replace h`var'=15 if `var'==4  & t1==13 /* hungary 11-20 hours */
replace h`var'=15 if `var'==4  & t1==17 /* romania 10-20 hours */
replace h`var'=15 if `var'==4  & t1==19 /* norway 11-20 hours */
replace h`var'=15 if `var'==4  & t1==21 /* italy 11-20 hours */
replace h`var'=15 if `var'==4  & t1==23 /* UK 11-20 hours */
replace h`var'=18 if `var'==5  /* 15-20 hours */
replace h`var'=25 if `var'==6  /* more than that */
replace h`var'=. if `var'==.a  /* question not asked */
replace h`var'=. if `var'==.b  /* question not answered */
}

gen canvassing=hb5a1
gen internet=hb5a5
gen rallies=hb5b3
gen nationalnews=hb5b5+hb5b7
replace nationalnews=hb5b7 if t1==16 /* Denmark doesn't ask about national newspapers (hb5b5) */
gen scalable=rallies+nationalnews+internet

save dta/ImportSurveyEur.dta, replace
