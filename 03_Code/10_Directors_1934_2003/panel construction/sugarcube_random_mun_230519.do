version 17
clear all
set more off
cap log close

*global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\" /* Office */
*global dump "C:\Users\schmutzy\Desktop\Projekt Nationalräte\" /* Office Dumpster */
global path "C:\Users\yanni\Dropbox\Projekt Nationalräte\" /* Home */
global dump "C:\Users\yanni\Desktop\Projekt Nationalräte\" /* Home Dumpster */
 
* Draw's already done, using the code below on 23rd of may *
/*Random draw of mun 
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_strict.dta", clear 
		quietly by gdenr: gen dup = cond(_N==1,0,_n)
		drop if dup == 0
		drop if dup<5
		drop if dup>5
		sample 15, count
		keep gdename gdenr
	save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\munlist_qmatch_`y'.dta", replace
}*/


/*Creating subsamples (matched within mun and capital)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	foreach x in strict /*medium_1*/ medium_2 /*medium_3*/ lenient {
		use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_`x'.dta", clear 
			drop _merge
			merge m:1 gdenr using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\munlist_qmatch_`y'.dta"
			keep if _merge == 3
			sort UCID
			keep CID UCID year cname gdename capital
		save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'.dta", replace
	}
}
*/

*Creating subsamples (matched within mun only)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	foreach x in strict /*medium_1*/ medium_2 /*medium_3*/ lenient {
		use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_`x'_mun.dta", clear 
			merge m:1 gdenr using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\munlist_qmatch_`y'.dta"
			keep if _merge == 3
			sort UCID
			keep CID UCID year cname gdename gdenr capital
		save "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'.dta", replace
	}
}

*Create Excel files to check false positive (FP)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	foreach x in strict /*medium_1*/ medium_2 /*medium_3*/ lenient {
		use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'.dta", clear
			bysort UCID : gen UCID_nb = _N
			order year cname gdename gdenr UCID UCID_nb capital CID
			sort gdenr UCID cname
		export excel using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'_FP.xlsx", replace firstrow(variables)
		putexcel set "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'_FP.xlsx", modify
		putexcel I1 = "False Positive"
	}
}
	
*Create Excel files to check false positive (FN)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	foreach x in strict /*medium_1*/ medium_2 /*medium_3*/ lenient {
		use "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'.dta", clear
			bysort UCID : gen UCID_nb = _N
			order year cname gdename UCID UCID_nb capital CID
			sort gdename cname UCID
		export excel using "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'_FN.xlsx", replace firstrow(variables)
		putexcel set "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'_FN.xlsx", modify
		putexcel H1 = "False Negative"
	}
}

*deleting temp files
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	foreach x in strict /*medium_1*/ medium_2 /*medium_3*/ lenient {
	erase "$dump\02_Processed_data\10_Directors_1934_2003\12_Random_Subsample\subsample_qmatch_`y'_`x'.dta"
	}
}

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