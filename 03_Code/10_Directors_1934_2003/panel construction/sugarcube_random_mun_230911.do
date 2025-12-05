version 17
clear all
set more off
cap log close

*global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"
global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\" /* Office */
global dump "C:\Users\schmutzy\Desktop\Projekt Nationalräte\" /* Office Dumpster */
*global path "C:\Users\yanni\Dropbox\Projekt Nationalräte\" /* Home */
*global dump "C:\Users\yanni\Desktop\Projekt Nationalräte\" /* Home Dumpster */
 
********************************************************************************
*Creating subsample using random draws (in munlist_qmatch_`y')

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
		use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_temp.dta", clear
			merge m:1 gdenr using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\Random draw\munlist_qmatch_`y'.dta"
			keep if _merge == 3
			keep CID year cname gdename gdenr capital
		save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'.dta", replace
	}

********************************************************************************
*** Collapse + 10% + 20% on evaluated subsamples *******************************
********************************************************************************

*********************************
*Transform evaluated xlsx into dta for merging

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	foreach x in strict medium_2 {
import excel "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\subsample_qmatch_`y'_`x'_FP.xlsx", sheet("Sheet1") firstrow clear
rename FalsePositive FalsePositive_`x'
save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'_eval.dta", replace 
	}
}

*********************************
*Merging to get the exact same obs & Perfect match 

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'.dta", clear 

merge m:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_strict_eval.dta"
	keep if _merge == 3
	drop _merge
merge m:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_medium_2_eval.dta"
	keep if _merge == 3
	drop _merge UCID UCID_nb
	
bysort gdenr : strgroup cname, gen(UCID) thresh(0.00000000001) force
	bysort UCID : gen UCID_nb = _N
			
order CID year cname UCID UCID_nb gdename gdenr capital
save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_eval_1.dta", replace
}

* Collapse perfect match (UCID)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_eval_1.dta", clear

*Build unique CID_str for all UCID_1
gen CID_str = string(CID)
order CID CID_str

bysort UCID: replace CID_str = CID_str[_n-1] + "; " + CID_str if UCID_nb > 1
replace CID_str = " "+CID_str
replace CID_str = usubinstr(CID_str, " ; ", "",.)
by UCID: replace CID_str = CID_str[_N]

*Delete duplicates
by UCID:  gen dup = cond(_N==1,0,_n)
drop if dup > 1

save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_eval_2.dta", replace
}

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_eval_2.dta", clear 	
bysort gdenr : strgroup cname, gen(UCID_1) thresh(0.1) force	
	bysort UCID_1 : gen UCID_1_nb = _N			
bysort gdenr : strgroup cname, gen(UCID_2) thresh(0.2) force	
	bysort UCID_2 : gen UCID_2_nb = _N
save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_eval_3.dta", replace
}

*Putting everything together
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_eval_3.dta", replace
		
gen FP_10_new = sum(FalsePositive_strict)
gen FP_20_new = sum(FalsePositive_medium_2)
gen count_new = 1
collapse year (max) FP_10 FP_20 (sum) count_new	

save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_eval_4.dta", replace
}			

***Create dta with nb of FP and FN by year and threshold
use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_1934_eval_4.dta", clear 
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
append using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_eval_4.dta"		
}
merge 1:1 year using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\FN_FP_count_1934_2003.dta"
drop FN_strict FN_medium_2 FN_lenient p_FN_strict p_FN_medium_2 p_FN_lenient _merge FP_lenient p_FP_lenient p_FP_strict p_FP_medium_2
rename (FP_strict FP_medium_2 count) (FP_10_old FP_20_old count_old)
order year count_new

gen diff_temp = count_new - count_old	
gen diff = sum(diff_temp)
gen obs_total = sum(count_old)
drop diff_temp

save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\FP_count_1934_2003.dta", replace

********************************************************************************
*deleting temp files
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	foreach x in strict medium_2 {
	erase "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'_eval.dta"
	}
}

foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	erase "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_eval_1.dta"
}


********************************************************************************
********************************************************************************
********************************************************************************
***Create Excel with nb of FP and FN by year and threshold
**False Positive
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	foreach x in strict /*medium_1*/ medium_2 /*medium_3*/ lenient  {
		import excel "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\subsample_qmatch_`y'_`x'_FP.xlsx", sheet("Sheet1") firstrow clear
		gen FP_`x' = sum(FalsePositive)
		gen count = 1
		collapse year (max) FP_`x' (sum) count
		save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_`y'_`x'.dta", replace 
	}
}

*Everything in one file for comparison
*FP in one file for each year
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_`y'_strict.dta", clear		
	foreach x in /*medium_1*/ medium_2 /*medium_3*/ lenient  {
		merge 1:1 year using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_`y'_`x'.dta"
		drop _merge
		merge 1:1 year using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_`y'_`x'.dta"
		drop _merge
		save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_`y'.dta", replace 
	}
}

*FP all year all threshold
use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_1934.dta", clear 
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
append using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_`y'.dta"		
}
save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_1934_2003.dta", replace

	
**False Negative	
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	foreach x in strict /*medium_1*/ medium_2 /*medium_3*/ lenient  {
		import excel "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\subsample_qmatch_`y'_`x'_FN.xlsx", sheet("Sheet1") firstrow clear
		gen FN_`x' = sum(FalseNegative)
		gen count = 1
		collapse year (max) FN_`x' (sum) count
		save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_`y'_`x'.dta", replace 
	}
}
	
*Everything in one file for comparison
*FN in one file for each year
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_`y'_strict.dta", clear		
	foreach x in /*medium_1*/ medium_2 /*medium_3*/ lenient  {
		merge 1:1 year using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_`y'_`x'.dta"
		drop _merge
		merge 1:1 year using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_`y'_`x'.dta"
		drop _merge
		save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_`y'.dta", replace 
	}
}

*FP all year all threshold
use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_1934.dta", clear 
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
append using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_`y'.dta"		
}
save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_1934_2003.dta", replace

***Erase intermediary files
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
erase "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_`y'.dta"		
erase "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_`y'.dta"		
}

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	foreach x in strict /*medium_1*/ medium_2 /*medium_3*/ lenient  {
		erase "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_`y'_`x'.dta"
		erase "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_`y'_`x'.dta"
	}
}

***All errors in one file
use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_count_1934_2003.dta", clear
merge 1:1 year using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FP_count_1934_2003.dta"
drop _merge
order year count
gen p_FN_strict = (FN_strict/count)*100 
gen p_FN_medium_2 = (FN_medium_2/count)*100 
gen p_FN_lenient = (FN_lenient/count)*100
gen p_FP_strict = (FP_strict/count)*100 
gen p_FP_medium_2 = (FP_medium_2/count)*100 
gen p_FP_lenient = (FP_lenient/count)*100 
save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False negative\FN_FP_count_1934_2003.dta", replace
save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\False positive\FN_FP_count_1934_2003.dta", replace
save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\FN_FP_count_1934_2003.dta", replace