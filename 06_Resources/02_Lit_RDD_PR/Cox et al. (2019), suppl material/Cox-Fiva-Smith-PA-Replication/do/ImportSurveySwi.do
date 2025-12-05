******************************************************************************************************************************************************
********************* SWITZERLAND SURVEY *************************************************************************************************************
******************************************************************************************************************************************************

use "dta/CH71-03parties.dta" , clear
sort kt year
save "dta/CH71-03parties_sort.dta", replace
 
use "dta/495_Selects_CumulativeFile_Data_1971-2015_v1.0.dta" , clear
gen kt=sg2*1
sort kt year
merge kt year using "dta/CH71-03parties_sort.dta"

foreach var in ci1 ci2 ci3 ci4 ci5 ci6 ci7 ci8 ci9 ci10 ci11 ci12 ci13 ci14 {
gen `var'_use=.
replace `var'_use=1 if `var'>0    /* "use" also include 0.33 and 0.67 which means "rarely" and "frequently" for some questions */
replace `var'_use=0 if `var'==0   /* don't use */
replace `var'_use=0 if `var'==-1  /* don't know  / can't remember */
replace `var'_use=. if `var'==-2  /* na */
replace `var'_use=. if `var'==.  
}

foreach var in sc1 maritals educ year religion {
tab `var', gen(d_`var')
}

save dta/ImportSurveySwi.dta, replace
