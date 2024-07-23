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
* Panel Construct - preparation										  		   *
********************************************************************************

*appending cleaned cross-sections
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_1934_match_3_clean.dta", clear
local agrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`a'_match_3_clean.dta"
}
*drop UCID_3 match_1 match_2 match_3 cname CID_str si_dummy_2 gdenr_2012 maxUCID_temp4 FP_group_2
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\complete_crosssections_1934_2003.dta", replace

*drop unecessary vars for matching and put into "characteristics.dta"
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`a'_match_3_clean.dta", clear
*keep UCID_3 CID year cname si_dummy_2 gdenr gdenr_2018 ctn cname_group CID_str 
order CID UCID_3 year cname gdenr gdenr_2018 si_dummy_2 cname_group CID_str
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_unique_id_`a'_match_3_clean_temp.dta", replace
}

********************************************************************************
* 1st step : Approach by pair of years to match over time 					   *
********************************************************************************

**********************************************************
*** construct pair of years (to match firms over time) ***

local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
    use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_unique_id_`a'_match_3_clean_temp.dta", clear
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_unique_id_`b'_match_3_clean_temp.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_id_`a'_`b'_match_3_clean.dta", replace 
}

**********************************************************
*** Generate unique ID over two years: perfect matches ***

local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002  
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
    use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_id_`a'_`b'_match_3_clean.dta", clear
		bysort gdenr_2018 : strgroup cname, gen(UCID_`a'_`b'_1) thresh(0.00000001) force
		order CID UCID_3 UCID_`a'_`b'
		sort UCID_`a'_`b' year
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_unique_id_`a'_`b'_match_1.dta", replace
}

********************************************************
*** Putting years together to form the Panel - Tests ***

* Put two pair of years together and if the "middle" year have been matched ///
* in each of them, it means the "1st" and "last" should be matched as well. ///
* We apply this logic sequentially to build the panel of firms exactly ///
* matched over time. One group = same name, same municipality, no holes in the ///
* serie (100% correct, same mun, same exact name, two "years" in a row).

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_unique_id_1934_1943_match_1.dta", clear
sort UCID_1934_1943_1 year
local agrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002  
local bgrp 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 
local cgrp 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 
local dgrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001  
foreach a of local agrp {    
    gettoken b bgrp : bgrp   
	gettoken c cgrp : cgrp
	gettoken d dgrp : dgrp
      
    *Merge datasets
		merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_unique_id_`a'_`b'_match_1.dta", generate(_merge`c')
		bysort UCID_`d'_`a'_1 (UCID_`a'_`b'_1) : replace UCID_`a'_`b'_1 = UCID_`a'_`b'_1[_n-1] if missing(UCID_`a'_`b'_1)
		sort UCID_`a'_`b'_1 year
	
	*Find firms disappearing between years, & were duplicates or unique
		quietly bysort UCID_`d'_`a'_1 : gen dup_2 = cond(_N==1,0,_n) if UCID_`a'_`b'_1 ==.
		sort UCID_`a'_`b'_1 UCID_`d'_`a'_1
	
	*Generate temp dataset for duplicate ones
		preserve
		drop if dup_2 ==. | dup_2 == 0
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_dup_1934_`b'_temp.dta", replace
		restore

	*Generate IDs for unique ones
		keep if dup_2 ==. | dup_2 == 0
		replace UCID_`a'_`b'_1 = UCID_`a'_`b'_1[_n-1] + 1 if missing(UCID_`a'_`b'_1)

	*Again, IDs but for duplicate ones
		append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_dup_1934_`b'_temp.dta"
		sort UCID_`a'_`b'_1 UCID_`d'_`a'_1
		egen temp = group(UCID_`d'_`a'_1) if missing(UCID_`a'_`b'_1)
		replace UCID_`a'_`b'_1 = UCID_`a'_`b'_1[_n-1] if missing(UCID_`a'_`b'_1)
		replace temp = 0 if missing(temp)
		gen UCID_`a'_`b'_temp = UCID_`a'_`b'_1 + temp
		drop UCID_`a'_`b'_1 temp dup_2
		rename UCID_`a'_`b'_temp UCID_`a'_`b'_1
		order CID UCID_3 UCID_`d'_`a'_1 UCID_`a'_`b'_1
}

order CID UCID_3
sort UCID_2002_2003 year
drop UCID_1934_1943 UCID_1943_1960 UCID_1960_1962 UCID_1962_1963 UCID_1963_1964 UCID_1964_1965 UCID_1965_1966 UCID_1966_1969 UCID_1969_1972 UCID_1972_1975 UCID_1975_1979 UCID_1979_1980 UCID_1980_1981 UCID_1981_1982 UCID_1982_1983 UCID_1983_1984 UCID_1984_1985 UCID_1985_1986 UCID_1986_1987 UCID_1987_1988 UCID_1988_1989 UCID_1989_1990 UCID_1990_1991 UCID_1991_1992 UCID_1992_1993 UCID_1993_1994 UCID_1994_1995 UCID_1995_1996 UCID_1996_1997 UCID_1997_1998 UCID_1998_1999 UCID_1999_2000 UCID_2000_2001 UCID_2001_2002 _merge1 _merge2 _merge3 _merge4 _merge5 _merge6 _merge7 _merge8 _merge9 _merge10 _merge11 _merge12 _merge13 _merge14 _merge15 _merge16 _merge17 _merge18 _merge19 _merge20 _merge21 _merge22 _merge23 _merge24 _merge25 _merge26 _merge27 _merge28 _merge29 _merge30 _merge31 _merge32 _merge33 _merge34

local agrp 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {     
		erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\comp_dup_1934_`a'_temp.dta"
}

rename UCID_2002_2003_1 UCID_panel_1

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003.dta", replace

********************************************************************************
* Separate first and last appearance of a firm within a UCID group
********************************************************************************

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003.dta", clear

sort UCID_panel_1 year
gen tag = 1
egen obs_total = total(tag), by(UCID_panel_1) //How many years in a row perfectly matched
quietly by UCID_panel_1:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_panel_1: egen maxyear = max(dup_panel)
bysort UCID_panel_1: egen minyear = min(dup_panel)

*Subsample of matched observations 
preserve
drop if dup_panel == 0
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_ok_full.dta", replace
restore

*Obs in bewteen first and last appearance (100% correct and not usefull to further match, same mun, same exact name, two "years" in a row)
preserve
drop if dup_panel == maxyear | dup_panel == minyear
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_middle_ok_1.dta", replace
restore

*First and last appareance in group in 1934 or 2003 (who are not single occurences) 
preserve 
keep if year == 1934 | year == 2003
keep if obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_first_end_ok_1.dta", replace
restore

*First and last appareance of a firm in each group
preserve
keep if dup_panel == maxyear | dup_panel == minyear
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_notok_1.dta", replace
restore

*Append two "ok" subsample
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_middle_ok_1.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_first_end_ok_1.dta"
sort UCID_panel_1 year
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_ok_1.dta", replace

**************************************************************************************
* Next matching step without "in-between" : use the different options in cname_group *
**************************************************************************************

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_notok_1.dta", replace

*Split in cross-sections
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp{
	preserve
		keep if year == `a'
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_`a'_notok_1_temp.dta", replace
	restore
}

* Build pair of years again
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
    use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_`a'_notok_1_temp.dta", clear
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_`b'_notok_1_temp.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_`a'_`b'_notok_1_temp.dta", replace 
}

*Split cname_group in multiple variables
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_`a'_`b'_notok_1_temp.dta", clear 
	
*Duplicates:
	bysort year UCID_panel_1: gen dup = cond(_N==1,0,_n) 	
	
	/*Delete potential duplicates coming from matching on gdenr :
	
	CID	UCID_3	UCID_panel_1	year	cname	gdenr	gdenr_2012	gdenr_2018
	17378	17667	815645	1960	millefiori sa	5002	5002	5002
	22636	17717	815645	1960	millefiori sa	5005	5005	5002

	=> They have similar UCID_panel_1 & year : wrong geolocalisation, decision to match on contemporary state of mun
	because two firms could exist with same name in two different municipalities before they merged. 
	It is even possible that they are registered as two distinct firms in the registry...
	
	=> We first need to correct the wrongly attributed municipalities
	*/
	
	*Corrections of wrongly attributed mun:
	
	*Luzern and Littau
	replace gdenr = 1060 if gdenr == 1061 & dup > 0
	replace gdename = "littau" if gdename == "luzern" & dup > 0
	
	*Lugano and multiple municipalities
	replace gdenr = 5172 if gdenr == 5192 & dup > 0 //Cassarate, Castagnola, Castagnola-Cassarate
	replace gdename = "castagnola" if gdename == "lugano" & dup > 0
	replace gdenr = 5158 if gdenr == 5192 & dup > 0 // Breganzona
	replace gdename = "breganzona" if gdename == "lugano" & dup > 0
	replace gdenr = 5147 if gdenr == 5192 & dup > 0 // Barbengo
	replace gdename = "barbengo" if gdename == "lugano" & dup > 0
	replace gdenr = 5168 if gdenr == 5192 & dup > 0 // Carabbia
	replace gdename = "carrabia" if gdename == "lugano" & dup > 0
	replace gdenr = 5184 if gdenr == 5192 & dup > 0 // Gandria
	replace gdename = "gandria" if gdename == "lugano" & dup > 0
	replace gdenr = 5211 if gdenr == 5192 & dup > 0 // Pazzalo
	replace gdename = "pazzalo" if gdename == "lugano" & dup > 0
	replace gdenr = 5215 if gdenr == 5192 & dup > 0 // Pregassona
	replace gdename = "pregassona" if gdename == "lugano" & dup > 0
	replace gdenr = 5234 if gdenr == 5192 & dup > 0 // Viganello
	replace gdename = "viganello" if gdename == "lugano" & dup > 0
	replace gdenr = 5209 if gdenr == 5192 & dup > 0 // Pambio noranco
	replace gdename = "pambio noranco" if gdename == "lugano" & dup > 0
	
	*Langanu LU
	replace gdenr = 1134 if gdenr == 1140 & dup > 0 //
	replace gdename = "langnau bei reiden" if gdename == "reiden" & dup > 0
	
	*Missing cname
	drop if missing(cname)
	
	
		*Count number of commas
		replace owners = usubinstr(owners,  ";", ",",.) if dup > 0
		gen comma = (length(trim(owners))-length(trim(subinstr(owners,",","",.))))
		sort UCID_panel_1 comma cname
		
		*No need to collapse cname, they have the same (perfect match)
		
		*Same with owners
		bysort UCID_panel_1 comma: replace owners = owners[_n-1] + "; " + owners
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		by UCID_panel_1: replace owners = owners[_N]

		*Same with function
		bysort UCID_panel_1 comma: replace function = function[_n-1] + "; " + function
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_panel_1: replace function = function[_N]

		*Same with signature
		bysort UCID_panel_1 comma: replace signature = signature[_n-1] + "; " + signature
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_panel_1: replace signature = signature[_N]

		*Same with page
		bysort UCID_panel_1 comma: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_panel_1: replace page = page[_N]

		*Delete duplicates
		drop if dup > 1
		sort UCID_panel_1 cname

	
	*Split cname_group into different variables
	split cname_group, p(";") g(cname_opt)
	order cname cname_opt*
		
		forval i = 1/2 {
			local j = `r(nvars)' + `i'
			gen cname_opt`j' = ""
			}
		
		forval v = 1/`r(nvars)' {
		replace cname_opt`v' = cname if cname_opt`v'==";"
		bysort cname_opt`v' gdenr_2018: gen dup`v' = cond(_N==1,0,_n)
		replace gdenr_2018 = 9999 if missing(gdenr_2018)
		egen UCID_panel_temp`v' = group(cname_opt`v' gdenr_2018)
		preserve
			keep if dup`v' > 0
			save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp`v'.dta", replace
		restore
			drop if dup`v' > 0
			local w = `v' + 1
			replace cname_opt`w' = cname_opt`v' if missing(cname_opt`w')
			drop dup`v'
		}
		
		*gen UCID_panel_temp = _n
		*save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_`a'_`b'_notok_2_temp.dta", replace		
}

*Append all temp subsamples by years
*1st to 4th round
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp1.dta", clear
forval i = 2/4{
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp`i'.dta"
	erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp`i'.dta"
}
	erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp1.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", replace
}

*add 5th for years where it exists
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 /*1992 1993 1994 1995*/ 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 /*1993 1994 1995 1996*/ 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", clear
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp5.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", replace
	erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp5.dta"
	}

*add 6th for years where it exists
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 /*1969 1972*/ 1975 1979 /*1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995*/ 1996 1997 1998 1999 2000 2001 /*2002*/
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 /*1972 1975*/ 1979 1980 /*1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996*/ 1997 1998 1999 2000 2001 2002 /*2003*/
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", clear
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp6.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", replace
	erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp6.dta"
	}
	
*add 7th for years where it exists
local agrp 1934 1943 1960 /*1962*/ 1963 1964 /*1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999*/ 2000 2001 /*2002*/
local bgrp 1943 1960 1962 /*1963*/ 1964 1965 /*1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000*/ 2001 2002 /*2003*/
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", clear
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp7.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", replace
	erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp7.dta"
	}
	
*add 8th for years where it exists
local agrp /*1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999*/ 2000 2001 /*2002*/
local bgrp /*1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000*/ 2001 2002 /*2003*/
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", clear
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp8.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", replace
	erase "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp8.dta"
	}

*Create unique group ID in those "cnameopt_`a'_`b'_ok_temp" subsamples
*Pair of years with UCID_panel_temp8
local agrp 2000 2001
local bgrp 2001 2002 
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", clear

	forval i = 1/8 {	
		replace UCID_panel_temp`i' =. if missing(dup`i')
	}
		egen maxUCID_panel_temp1=max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp1 + UCID_panel_temp2 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp2 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp2 + UCID_panel_temp3 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp3 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp3 + UCID_panel_temp4 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp4 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp4 + UCID_panel_temp5 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp5 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp5 + UCID_panel_temp6 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp6 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp6 + UCID_panel_temp7 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp7 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp7 + UCID_panel_temp7 if missing(UCID_panel_temp1)
	
	egen UCID_panel_2 = group(UCID_panel_temp1)
	
	drop UCID_panel_temp* cname_opt* dup* maxUCID_panel_temp*
	
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp_2.dta", replace
}


*Pair of years with UCID_panel_temp7
local agrp 1934 1943 1960 1963 1964 2000 2001
local bgrp 1943 1960 1962 1964 1965 2001 2002
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", clear

	forval i = 1/7 {	
		replace UCID_panel_temp`i' =. if missing(dup`i')
	}
		egen maxUCID_panel_temp1=max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp1 + UCID_panel_temp2 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp2 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp2 + UCID_panel_temp3 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp3 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp3 + UCID_panel_temp4 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp4 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp4 + UCID_panel_temp5 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp5 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp5 + UCID_panel_temp6 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp6 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp6 + UCID_panel_temp7 if missing(UCID_panel_temp1)
	
	egen UCID_panel_2 = group(UCID_panel_temp1)
	
	drop UCID_panel_temp* cname_opt* dup* maxUCID_panel_temp*
	
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp_2.dta", replace
}

*Pair of years with UCID_panel_temp6
local agrp 1962 1965 1966 1975 1979 1996 1997 1998 1999 2000 2001
local bgrp 1963 1966 1969 1979 1980 1997 1998 1999 2000 2001 2002
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", clear

	forval i = 1/6 {	
		replace UCID_panel_temp`i' =. if missing(dup`i')
	}
		egen maxUCID_panel_temp1=max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp1 + UCID_panel_temp2 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp2 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp2 + UCID_panel_temp3 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp3 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp3 + UCID_panel_temp4 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp4 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp4 + UCID_panel_temp5 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp5 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp5 + UCID_panel_temp6 if missing(UCID_panel_temp1)
		
	egen UCID_panel_2 = group(UCID_panel_temp1)
	
	drop UCID_panel_temp* cname_opt* dup* maxUCID_panel_temp*
		
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp_2.dta", replace
}

*Pair of years with UCID_panel_temp5
local agrp 1969 1972 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1999 2000 2001 2002
local bgrp 1972 1975 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", clear

	forval i = 1/5 {	
		replace UCID_panel_temp`i' =. if missing(dup`i')
	}
		egen maxUCID_panel_temp1=max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp1 + UCID_panel_temp2 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp2 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp2 + UCID_panel_temp3 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp3 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp3 + UCID_panel_temp4 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp4 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp4 + UCID_panel_temp5 if missing(UCID_panel_temp1)
		
	egen UCID_panel_2 = group(UCID_panel_temp1)
	
	drop UCID_panel_temp* cname_opt* dup* maxUCID_panel_temp* 
		
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp_2.dta", replace
}

*Pair of years with UCID_panel_temp4
local agrp 1992 1993 1994 1995
local bgrp 1993 1994 1995 1996
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp.dta", clear

	forval i = 1/4 {	
		replace UCID_panel_temp`i' =. if missing(dup`i')
	}
		egen maxUCID_panel_temp1=max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp1 + UCID_panel_temp2 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp2 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp2 + UCID_panel_temp3 if missing(UCID_panel_temp1)
		
		egen maxUCID_panel_temp3 = max(UCID_panel_temp1)
		replace UCID_panel_temp1 = maxUCID_panel_temp3 + UCID_panel_temp4 if missing(UCID_panel_temp1)
		
	egen UCID_panel_2 = group(UCID_panel_temp1)
	
	drop UCID_panel_temp* cname_opt* dup* maxUCID_panel_temp* 
		
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp_2.dta", replace
}
		
*Put that back into full data by merging on CID year & create yearly group ID 
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003.dta", clear
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp
		merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\cnameopt_`a'_`b'_ok_temp_2.dta", nogen
		rename UCID_panel_2 UCID_panel_2_`a'_`b'
		bysort UCID_panel_2_`a'_`b' : gen dup_UCID_panel_2_`a'_`b' = cond(_N==1,0,_n) if ! missing(UCID_panel_2_`a'_`b')
}

local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
local cgrp 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 
foreach a of local agrp {
	gettoken b bgrp : bgrp
	gettoken c cgrp : cgrp
		sort UCID_panel_2_`a'_`b' year
		by UCID_panel_2_`a'_`b': replace UCID_panel_2_`b'_`c' = UCID_panel_2_`b'_`c'[_N] if missing(UCID_panel_2_`b'_`c')
}
//check and correct loop, it works but stops only when files don't exist

local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {
	gettoken b bgrp : bgrp
		sort UCID_panel_1 year
		by UCID_panel_1: gen x = _n
		by UCID_panel_1: mipolate UCID_panel_2_`a'_`b' x, gen(y1) groupwise
		sort y1 year
		
		by y1: replace UCID_panel_1 = UCID_panel_1[_N] if ! missing(y1)
		drop UCID_panel_2_`a'_`b' dup_UCID_panel_2_`a'_`b' x 
		rename y1 UCID_panel_2_`a'_`b'
}

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_2.dta", replace
 
********************************************************************************
* Separate first and last appearance of a firm within a UCID group
* Sample with match over different ways of writting cname 
********************************************************************************

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_2.dta", clear

sort UCID_panel_1 year
drop tag obs_total maxyear minyear UCID_panel_2*

gen tag = 1
egen obs_total = total(tag), by(UCID_panel_1) //How many years in a row perfectly matched
quietly by UCID_panel_1:  gen dup_panel = cond(_N==1,0,_n)
bysort UCID_panel_1: egen maxyear = max(dup_panel)
bysort UCID_panel_1: egen minyear = min(dup_panel)

*Subsample of matched observations : example: two obs are matched together, three others, but they should be one single group of 5. We are searching for them in next steps.
preserve
drop if dup_panel == 0
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_ok_full_2.dta", replace
restore

*Obs in bewteen first and last appearance (100% correct and not usefull to further match, same mun, same exact name, two "years" in a row)
preserve
drop if dup_panel == maxyear | dup_panel == minyear
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_middle_ok_2.dta", replace
restore

*First and last appareance in group (not single occurences) 
preserve 
keep if year == 1934 | year == 2003
keep if obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_first_end_ok_2.dta", replace
restore

*First and last appareance of a firm in each group : unmatched subsample
preserve
keep if dup_panel == maxyear | dup_panel == minyear
drop if year == 2003 & obs_total > 1
drop if year == 1934 & obs_total > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_notok_2.dta", replace
restore

*Append two "ok" subsample
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_middle_ok_2.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_first_end_ok_2.dta"
sort UCID_panel_1 year
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_ok_2.dta", replace
 
********************************************************************************
* Fuzzy matches within municipalities over pairs of years time

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_notok_2.dta", clear

*Split in cross-sections
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp{
	preserve
		keep if year == `a'
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_`a'_notok_2_temp.dta", replace
	restore
}

* Build pair of years again
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
    use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_`a'_notok_2_temp.dta", clear
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_`b'_notok_2_temp.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_`a'_`b'_notok_2_temp.dta", replace 
}


***
*Directly across every periods within mun, run on server on 26th february ///
 and generated fuzzy_match dataset used below

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_notok_2.dta", clear
bysort gdenr_2018 : strgroup cname, gen(UCID_`a'_`b') thresh(0.1) force
sort UCID_1934_2003 year
save "$server\match_perfect_1934_2003_fuzzy_match_1.dta", replace
*Run on server 26th of february
***

*Resulting subsample of matched observations
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_fuzzy_match_1", clear

*drop unmatched
quietly by UCID_1934_2003:  gen dup = cond(_N==1,0,_n)
drop if dup == 0

*drop those who were matched in perfect matching and are not match with other (new) obs
sort UCID_panel_1 year
quietly by UCID_panel_1:  gen dup2 = cond(_N==1,0,_n)
by UCID_panel_1: egen maxdup2 = max(dup2)
bysort UCID_1934_2003:  egen maxdup = max(dup)
drop if maxdup2 == maxdup
sort UCID_1934_2003 year

*Create files to be evaluated by hand
gen periode =.
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
local bgrp 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
foreach a of local agrp { 
	gettoken b bgrp : bgrp
	replace periode = `b' if year == `a'
}

order CID UCID_panel_1 UCID_1934_2003 year periode cname gdenr_2018 dup_panel maxyear minyear owners capital function signature page
gen FalsePositive =.
splitsample, generate(split_var) cluster(UCID_1934_2003)
drop if missing(cname)
sort split_var UCID_1934_2003 year 

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_fuzzy_1934_2003_tocheck", replace

preserve
drop if split_var == 1
keep CID UCID_panel_1 UCID_1934_2003 year periode cname gdenr_2018 dup_panel maxyear minyear owners capital function signature page FalsePositive
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Evaluation\match_fuzzy_1934_2003_tocheck_ED.xlsx", replace firstrow(variables)
restore

preserve
drop if split_var == 2
keep CID UCID_panel_1 UCID_1934_2003 year periode cname gdenr_2018 dup_panel maxyear minyear owners capital function signature page FalsePositive
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Evaluation\match_fuzzy_1934_2003_tocheck_JU.xlsx", replace firstrow(variables)
restore

********************************************************************************
* Corrections following JU and ED evaluations

*save as dta
import excel "$dump\02_Processed_data\11_Matching_Eval\match_fuzzy_1934_2003_tocheck_ED.xlsx", sheet("Sheet1") firstrow clear
save "$dump\02_Processed_data\11_Matching_Eval\match_fuzzy_1934_2003_check_ED", replace

import excel "$dump\02_Processed_data\11_Matching_Eval\match_fuzzy_1934_2003_tocheck_JU.xlsx", sheet("Sheet1") firstrow clear
save "$dump\02_Processed_data\11_Matching_Eval\match_fuzzy_1934_2003_check_JU", replace

*append JU and ED eval
use "$dump\02_Processed_data\11_Matching_Eval\match_fuzzy_1934_2003_check_ED", clear
append using "$dump\02_Processed_data\11_Matching_Eval\match_fuzzy_1934_2003_check_JU"
save "$dump\02_Processed_data\11_Matching_Eval\match_fuzzy_1934_2003_check_ED_JU", replace

* "fuzzy_match_1" = JU + ED + some obs...
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_fuzzy_match_1", clear 
merge 1:1 CID year using "$dump\02_Processed_data\11_Matching_Eval\match_fuzzy_1934_2003_check_ED_JU", nogen

order CID UCID_3 UCID_panel_1 UCID_1934_2003 year periode FalsePositive
sort UCID_1934_2003 year

*add ok subsample (usefull here?)
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_ok_2.dta"
sort UCID_panel_1 year UCID_1934_2003

*Attribute UCID_1934_2003 to obs in the beginning (for 1934), middle and end (for 2003) of series
by UCID_panel_1: replace UCID_1934_2003 = UCID_1934_2003[_n-1] if missing(UCID_1934_2003)
forval x = 1/34{
by UCID_panel_1: replace UCID_1934_2003 = UCID_1934_2003[_n+1] if missing(UCID_1934_2003)
}

*Deal with those who are matched over the whole sample periods (1934 to 2003 have no UCID_1934_2003) & those who are only in one year (1934 or 2003) : attribute UCID_1934_2003 to obs 

egen maxUCID_1934_2003 = max(UCID_1934_2003)
replace UCID_1934_2003 = maxUCID_1934_2003 + UCID_panel_1 if missing(UCID_1934_2003)


********************************************************************************
*Deal with duplicates UCID_1934_2003 & year (Collapse?)
*tag them & exclude them first

bysort year UCID_1934_2003: gen dup = cond(_N==1,0,_n) 	
sort dup UCID_1934_2003 year 

sort UCID_1934_2003 year  UCID_panel_1  
by UCID_1934_2003: replace FalsePositive = FalsePositive[_n-1] if missing(FalsePositive)
forval x = 1/34{
by UCID_1934_2003: replace FalsePositive = FalsePositive[_n+1] if missing(FalsePositive)
}

*Subsample with groups where there is a FP case
preserve
keep if FalsePositive == 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\fuzzy10_falsepositive", replace
restore

*Subsample with ok matches
preserve
keep if FalsePositive == 0
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\fuzzy10_nofalsepositive", replace
restore

******************************
* START HERE*
******************************

replace UCID_1934_2003 = maxUCID_1934_2003 + UCID_panel_1 if FalsePositive == 1
bysort year UCID_1934_2003: gen dup2 = cond(_N==1,0,_n) 	

/*keep if dup2 > 0
keep if si_dummy_2 == 1
*/

*Witouht FP cases
drop if FalsePositive == 1

