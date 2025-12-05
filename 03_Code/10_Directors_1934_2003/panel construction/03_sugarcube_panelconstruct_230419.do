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
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003.dta", clear // this is were we reintroduce the correct matches from stop 1 (basically the "match_perfect_1934_2003_ok_1" file)
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

*Subsample of matched observations : example: two obs are matched together, three others, but they should be one single group of 5. We are searching for them in next steps. Note: just for informational purpose, not used later.
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
* Fuzzy matches within municipalities over all periods

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

/*use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_1934_2003_notok_2.dta", clear
bysort gdenr_2018 : strgroup cname, gen(UCID_1934_2003) thresh(0.1) force
sort UCID_1934_2003 year
save "$server\match_perfect_1934_2003_fuzzy_match_1.dta", replace
*Run on server 26th of february*/
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
*Deal with duplicates UCID_1934_2003 & year 
*tag them & exclude them first

*Identify duplicate UCID_1934_2003 in one year (twice same ID in a year)
bysort year UCID_1934_2003: gen dup = cond(_N==1,0,_n) 	
sort dup UCID_1934_2003 year 

*Identify group where ED and JU found FalsePositive
gen FalsePositive_orig = FalsePositive
rename FalsePositive FalsePositive_group
sort UCID_1934_2003 year UCID_panel_1  
by UCID_1934_2003: replace FalsePositive_group = FalsePositive_group[_n-1] if missing(FalsePositive_group)
forval x = 1/34{
by UCID_1934_2003: replace FalsePositive_group = FalsePositive_group[_n+1] if missing(FalsePositive_group)
}

order CID UCID_3 UCID_panel_1 UCID_1934_2003 year periode FalsePositive_group FalsePositive_orig

*Give back UCID_panel_1 groups if dup (going back to perfect match)
replace UCID_1934_2003 = maxUCID_1934_2003 + UCID_panel_1 if dup > 0
replace UCID_1934_2003 = maxUCID_1934_2003 + UCID_panel_1 if FalsePositive_group == 1

********************************************************************************
***Identify further duplicates (same ID in a year): those who were already duplicate before. This goes back to the cross-sections.

bysort year UCID_1934_2003: gen dup2 = cond(_N==1,0,_n) 	
sort dup2 UCID_1934_2003 year

*There are 1353 obs with same ID in one year, will deal with that later as they will come up anyway with another matching round
*Taged by dup2

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003.dta", replace

********************************************************************************
********************************************************************************
* Matching without geographical restrictions
********************************************************************************

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003.dta", clear

sort UCID_1934_2003 year

gen tag_2 = 1
egen obs_total_2 = total(tag_2), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel_2 = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxyear_2 = max(dup_panel_2)
bysort UCID_1934_2003: egen minyear_2 = min(dup_panel_2)

*Subsample of matched observations : example: two obs are matched together, three others, but they should be one single group of 5. We are searching for them in next steps. Note: just for informational purpose, not used later.
preserve
drop if dup_panel_2 == 0
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003_ok_full_2.dta", replace
restore

*Obs in bewteen first and last appearance (100% correct and not usefull to further match, same mun, same exact name, two "years" in a row)
preserve
drop if dup_panel_2 == maxyear_2 | dup_panel_2 == minyear_2
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003_middle_ok_2.dta", replace
restore

*First and last appareance in group starting in 1934 or ending in 2003 (not single occurences) 
preserve 
keep if year == 1934 | year == 2003
keep if obs_total_2 > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003_first_end_ok_2.dta", replace
restore

*First and last appareance of a firm in each group : unmatched subsample
preserve
keep if dup_panel_2 == maxyear_2 | dup_panel_2 == minyear_2
drop if year == 2003 & obs_total_2 > 1
drop if year == 1934 & obs_total_2 > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003_notok_2.dta", replace
restore

*Append two "ok" subsample
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003_middle_ok_2.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003_first_end_ok_2.dta"
sort UCID_1934_2003 year
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003_ok_2.dta", replace

********************************************************************************
* Perfect match across municipalities (pairs of years)
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003_notok_2.dta", clear
sort cname year

local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
local bgrp 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
foreach a of local agrp { 
	gettoken b bgrp : bgrp
	replace periode = `b' if year == `a'
}

drop if missing(cname)

*Split in cross-sections
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp{
	preserve
		keep if year == `a'
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_notok_2_temp.dta", replace
	restore
}

* Build pair of years again
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
    use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_notok_2_temp.dta", clear
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`b'_notok_2_temp.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_`b'_notok_2_temp.dta", replace 
}

*Match across mun, over pairs of year
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002  
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
    use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_`b'_notok_2_temp.dta", clear
		strgroup cname, gen(UCID_`a'_`b') thresh(0.000000000001) force
		order CID UCID_`a'_`b'
		sort UCID_`a'_`b' year
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_`a'_`b'.dta", replace
}

********************************************************************************
*Building panel back

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_1943.dta", clear
order CID UCID_1934_1943 

local agrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002  
local bgrp 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 
local cgrp 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 
local dgrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001  
foreach a of local agrp {    
    gettoken b bgrp : bgrp   
	gettoken c cgrp : cgrp
	gettoken d dgrp : dgrp

		merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_`a'_`b'.dta", generate(_merge`c')
		order CID UCID_`a'_`b' 
		sort UCID_`a'_`b' UCID_`d'_`a'
		bysort UCID_`d'_`a' (UCID_`a'_`b'): replace UCID_`a'_`b'  = UCID_`a'_`b'[_n-1] if missing(UCID_`a'_`b')
		sort UCID_`a'_`b' UCID_`d'_`a' year
		egen maxUCID_`a'_`b' = max(UCID_`a'_`b')
		replace UCID_`a'_`b' = maxUCID_`a'_`b'  + UCID_`d'_`a' if missing(UCID_`a'_`b')
}
	
drop UCID_1934_1943 UCID_1943_1960 UCID_1960_1962 UCID_1962_1963 UCID_1963_1964 UCID_1964_1965 UCID_1965_1966 UCID_1966_1969 UCID_1969_1972 UCID_1972_1975 UCID_1975_1979 UCID_1979_1980 UCID_1980_1981 UCID_1981_1982 UCID_1982_1983 UCID_1983_1984 UCID_1984_1985 UCID_1985_1986 UCID_1986_1987 UCID_1987_1988 UCID_1988_1989 UCID_1989_1990 UCID_1990_1991 UCID_1991_1992 UCID_1992_1993 UCID_1993_1994 UCID_1994_1995 UCID_1995_1996 UCID_1996_1997 UCID_1997_1998 UCID_1998_1999 UCID_1999_2000 UCID_2000_2001 UCID_2001_2002 _merge1 _merge2 _merge3 _merge4 _merge5 _merge6 _merge7 _merge8 _merge9 _merge10 _merge11 _merge12 _merge13 _merge14 _merge15 _merge16 _merge17 _merge18 _merge19 _merge20 _merge21 _merge22 _merge23 _merge24 _merge25 _merge26 _merge27 _merge28 _merge29 _merge30 _merge31 _merge32 _merge33 _merge34 maxUCID_1943_1960 maxUCID_1960_1962 maxUCID_1962_1963 maxUCID_1963_1964 maxUCID_1964_1965 maxUCID_1965_1966 maxUCID_1966_1969 maxUCID_1969_1972 maxUCID_1972_1975 maxUCID_1975_1979 maxUCID_1979_1980 maxUCID_1980_1981 maxUCID_1981_1982 maxUCID_1982_1983 maxUCID_1983_1984 maxUCID_1984_1985 maxUCID_1985_1986 maxUCID_1986_1987 maxUCID_1987_1988 maxUCID_1988_1989 maxUCID_1989_1990 maxUCID_1990_1991 maxUCID_1991_1992 maxUCID_1992_1993 maxUCID_1993_1994 maxUCID_1994_1995 maxUCID_1995_1996 maxUCID_1996_1997 maxUCID_1997_1998 maxUCID_1998_1999 maxUCID_1999_2000 maxUCID_2000_2001 maxUCID_2001_2002 maxUCID_2002_2003

rename UCID_2002_2003 UCID_panel_2
sort UCID_panel_2 year

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_all.dta", replace

********************************************************************************
* Match same year, same name, different municipalities
* Sometimes OCR (immobilien AG, 21 times)
* Sometimes due to cleaning of "Holding" (you donkey)
* Sometimes geocoding wrong (no mun, wrong mun)

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_all.dta", clear
sort UCID_panel_2 year
egen maxUCID_panel_2 = max(UCID_panel_2)
quietly by UCID_panel_2 year:  gen REAL_PROBLEM = cond(_N==1,0,_n) if UCID_panel_2 !=. 
replace UCID_panel_2 = UCID_panel_2 + UCID_1934_2003 + maxUCID_panel_2 if REAL_PROBLEM > 0
gen UCID_panel_2_duplicates = UCID_panel_2

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_all_temp.dta", replace

********************************************************************************
* Attribute same IDs

use "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_1934_2003.dta", clear 
sort CID year

quietly by CID year:  gen PROBLEM = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
drop if PROBLEM == 2
drop PROBLEM

merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_all_temp.dta"

gen UCID_1934_2003_temp = UCID_1934_2003
order CID UCID_3 UCID_1934_2003 UCID_1934_2003_temp UCID_panel_1 UCID_panel_2 year periode cname gdenr gdenr_2018

forval x = 1/8{
sort UCID_panel_2 year
by UCID_panel_2: replace UCID_1934_2003 = UCID_1934_2003[1] if UCID_panel_2 !=.
sort UCID_1934_2003_temp year
by UCID_1934_2003_temp: replace UCID_1934_2003 = UCID_1934_2003[1]
}

sort UCID_1934_2003 year
quietly by UCID_1934_2003 year:  gen REAL_PROBLEM_2 = cond(_N==1,0,_n) /* SAME ID IN SAME YEAR ?*/

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003.dta", replace

********************************************************************************
* Subsample of beginning and ending of series
********************************************************************************

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003.dta", clear /* Full sample after perfect match with UCID_1934_2003 for exact cname over pair of years w/o geo restrictions */
sort UCID_1934_2003 year

gen tag_3 = 1
egen obs_total_3 = total(tag_3), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel_3 = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxyear_3 = max(dup_panel_3)
bysort UCID_1934_2003: egen minyear_3 = min(dup_panel_3)

*Subsample of matched observations : example: two obs are matched together, three other together, but they should be one single group of 5. We are searching for them in next steps. Note: just for informational purpose, not used later.
preserve
drop if dup_panel_3 == 0
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003_ok_full_2.dta", replace
restore

*Obs in bewteen first and last appearance (100% correct and not usefull to further match, same mun, same exact name, two "years" in a row)
preserve
drop if dup_panel_3 == maxyear_3 | dup_panel_3 == minyear_3
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003_middle_ok_2.dta", replace
restore

*First and last appareance in group starting in 1934 or ending in 2003 (not single occurences) 
preserve 
keep if year == 1934 | year == 2003
keep if obs_total_3 > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003_first_end_ok_2.dta", replace
restore

*First and last appareance of a firm in each group : unmatched subsample
preserve
keep if dup_panel_3 == maxyear_3 | dup_panel_3 == minyear_3
drop if year == 2003 & obs_total_3 > 1
drop if year == 1934 & obs_total_3 > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003_notok_2.dta", replace
restore

*Append two "ok" subsample
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003_middle_ok_2.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003_first_end_ok_2.dta"
sort UCID_1934_2003 year
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003_ok_2.dta", replace

********************************************************************************
* 10% matches across municipalities 

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003_notok_2.dta", clear
sort cname year
drop if missing(cname) /* 6 Obs without any firm names*/

*Split in cross-sections
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp{
	preserve
		keep if year == `a'
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_notok_3_temp.dta", replace
	restore
}

* Build pair of years again
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
    use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_notok_3_temp.dta", clear
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`b'_notok_3_temp.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_`b'_notok_3_temp.dta", replace 
}

*Match across mun, over pairs of year
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002  
local bgrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
    use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_`b'_notok_3_temp.dta", clear
		strgroup cname, gen(UCID_`a'_`b') thresh(0.1) force
		order CID UCID_`a'_`b'
		sort UCID_`a'_`b' year
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_`a'_`b'_3.dta", replace
}

********************************************************************************
*Building panel back
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_1943_3.dta", clear
order CID UCID_1934_1943 

local agrp 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002  
local bgrp 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 
local cgrp 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 
local dgrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001  
foreach a of local agrp {    
    gettoken b bgrp : bgrp   
	gettoken c cgrp : cgrp
	gettoken d dgrp : dgrp

		merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_`a'_`b'_3.dta", generate(_merge`c')
		order CID UCID_`a'_`b' 
		sort UCID_`a'_`b' UCID_`d'_`a'
		bysort UCID_`d'_`a' (UCID_`a'_`b'): replace UCID_`a'_`b'  = UCID_`a'_`b'[_n-1] if missing(UCID_`a'_`b')
		sort UCID_`a'_`b' UCID_`d'_`a' year
		egen maxUCID_`a'_`b' = max(UCID_`a'_`b')
		replace UCID_`a'_`b' = maxUCID_`a'_`b'  + UCID_`d'_`a' if missing(UCID_`a'_`b')
}
	
drop UCID_1934_1943 UCID_1943_1960 UCID_1960_1962 UCID_1962_1963 UCID_1963_1964 UCID_1964_1965 UCID_1965_1966 UCID_1966_1969 UCID_1969_1972 UCID_1972_1975 UCID_1975_1979 UCID_1979_1980 UCID_1980_1981 UCID_1981_1982 UCID_1982_1983 UCID_1983_1984 UCID_1984_1985 UCID_1985_1986 UCID_1986_1987 UCID_1987_1988 UCID_1988_1989 UCID_1989_1990 UCID_1990_1991 UCID_1991_1992 UCID_1992_1993 UCID_1993_1994 UCID_1994_1995 UCID_1995_1996 UCID_1996_1997 UCID_1997_1998 UCID_1998_1999 UCID_1999_2000 UCID_2000_2001 UCID_2001_2002 _merge1 _merge2 _merge3 _merge4 _merge5 _merge6 _merge7 _merge8 _merge9 _merge10 _merge11 _merge12 _merge13 _merge14 _merge15 _merge16 _merge17 _merge18 _merge19 _merge20 _merge21 _merge22 _merge23 _merge24 _merge25 _merge26 _merge27 _merge28 _merge29 _merge30 _merge31 _merge32 _merge33 _merge34 maxUCID_1943_1960 maxUCID_1960_1962 maxUCID_1962_1963 maxUCID_1963_1964 maxUCID_1964_1965 maxUCID_1965_1966 maxUCID_1966_1969 maxUCID_1969_1972 maxUCID_1972_1975 maxUCID_1975_1979 maxUCID_1979_1980 maxUCID_1980_1981 maxUCID_1981_1982 maxUCID_1982_1983 maxUCID_1983_1984 maxUCID_1984_1985 maxUCID_1985_1986 maxUCID_1986_1987 maxUCID_1987_1988 maxUCID_1988_1989 maxUCID_1989_1990 maxUCID_1990_1991 maxUCID_1991_1992 maxUCID_1992_1993 maxUCID_1993_1994 maxUCID_1994_1995 maxUCID_1995_1996 maxUCID_1996_1997 maxUCID_1997_1998 maxUCID_1998_1999 maxUCID_1999_2000 maxUCID_2000_2001 maxUCID_2001_2002 maxUCID_2002_2003

rename UCID_2002_2003 UCID_panel_3
sort UCID_panel_3 year

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_all.dta", replace

********************************************************************************
*** Subsample of obs matched
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_all.dta", clear
	gen UCID_1934_2003_temp_2 = UCID_1934_2003 /* used later to attribute IDs over time */
	quietly by UCID_panel_3 year:  gen REAL_PROBLEM_3 = cond(_N==1,0,_n) if UCID_panel_3 !=.

	**droping unmatched obs
		sort UCID_panel_3 year
		quietly by UCID_panel_3:  gen matched = cond(_N==1,0,_n) /* Nb of obs within panel_2 groups */
		drop if matched == 0

	**droping those who have no additional match
		*Nb of obs matched in previous step (to compare if same)
		sort UCID_1934_2003 year
		quietly by UCID_1934_2003:  gen matched_UCID_1934_2003 = cond(_N==1,0,_n) /* Nb of obs within UCID_1934_2003 groups */
		sort UCID_panel_3 year

		*Compare the nb of obs in both groups
		by UCID_panel_3: egen max_matched = max(matched)
		by UCID_panel_3: egen max_matched_UCID_1934_2003 = max(matched_UCID_1934_2003)
		
		/* Doesn't work because two matched obs (max_matched == 2) could also have been matched to another obs which is not the new one (max_matched_UCID_1934_2003 == 2)
			Typically with beginning and end of a serie that are not matched together anymore
		*/
	
		*Drop if no new matches
		drop if max_matched == max_matched_UCID_1934_2003 & REAL_PROBLEM_3 == 0
		
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_all_temp.dta", replace

********************************************************************************
*** Files to evaluate by ED & JU

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_all_temp.dta", clear

drop UCID_3 UCID_1934_2003_temp UCID_panel_1 UCID_panel_2 FalsePositive_group FalsePositive_orig si_dummy_2 match_1 match_2 match_3 caddress_orig gdename_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract E_CNTR N_CNTR gdenr_2012 maxUCID_temp4 FP_group_2 comma tag maxUCID_1934_2003 dup dup2 tag_2  maxUCID_panel_2 REAL_PROBLEM UCID_panel_2_duplicates _merge REAL_PROBLEM_2 tag_3 UCID_1934_2003_temp_2 matched matched_UCID_1934_2003 max_matched max_matched_UCID_1934_2003 CID_str REAL_PROBLEM_3 obs_total dup_panel maxyear minyear obs_total_2 dup_panel_2 maxyear_2 minyear_2 maxyear_3 minyear_3 function signature page

order CID UCID_panel_3 UCID_1934_2003 year periode cname cname_orig 

local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
local bgrp 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
foreach a of local agrp { 
	gettoken b bgrp : bgrp
	replace periode = `b' if year == `a'
}
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_toeval.dta", replace

********************************************************************************
* Test to link owners ID with owners name
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_toeval.dta", clear

replace owners = usubinstr(owners, ";", ",",.)
replace owners = usubinstr(owners, "   ", "  ",.)
replace owners = usubinstr(owners, "  ", " ",.)
replace owners = usubinstr(owners, " ", "",.)
split owners, parse(,) generate(owner)

generate id = _n
reshape long owner, i(UCID_panel_3 id) j(test)

drop if missing(owner)
rename owner PID
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_toeval_long.dta", replace

*Cleaning of persons list' (to have single PID per year)
use PID year firstname lastname CID using "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", clear

quietly by PID year:  gen dup = cond(_N==1,0,_n)
drop if dup > 1

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\PID_year_list.dta", replace

*Merge
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_toeval_long.dta", clear
destring PID, replace
merge m:1 PID year using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\PID_year_list.dta"
keep if _merge == 3

g fullname = firstname + " " + lastname
drop firstname lastname dup _merge

*Reshape to have one unique CID with multiple var with owners
sort CID year
bysort CID year: replace fullname = fullname[_n-1] + "; " + fullname
replace fullname = " "+fullname
replace fullname = usubinstr(fullname, " ; ", "",.)
by CID year: replace fullname = fullname[_N]

by CID year:  gen dup = cond(_N==1,0,_n)
drop if dup > 1
sort UCID_panel_3 year
keep CID year cname fullname

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\CID_year_owners.dta", replace

*Adding full names to file to evaluate
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_toeval.dta", clear
merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\CID_year_owners.dta", nogen 
sort UCID_panel_3 year
order CID UCID_panel_3 UCID_1934_2003 year periode cname gdenr gdenr_2018 fullname
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_toeval_fullnames.dta", replace
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Evaluation\match_10_nogeo_toeval_fullnames_ED.xlsx", replace firstrow(variables)
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Evaluation\match_10_nogeo_toeval_fullnames_JU.xlsx", replace firstrow(variables)

*compare their evaluation
import excel "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_toeval_fullnames_ED_prep.xlsx", sheet("Sheet1") firstrow clear
rename group group_ED
save "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_toeval_fullnames_ED_prep", replace

import excel "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_toeval_fullnames_JU_prep.xlsx", sheet("Sheet1") firstrow clear
rename group group_JU
save "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_toeval_fullnames_JU_prep", replace

merge 1:1 CID year using "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_toeval_fullnames_ED_prep"
sort UCID_panel_3 year

gen group_JU_temp = 0 if missing(group_JU)
replace group_JU_temp = 1 if missing(group_JU_temp)
gen group_ED_temp = 0 if missing(group_ED)
replace group_ED_temp = 1 if missing(group_ED_temp)

gen match_ok = 1 if missing(group_JU) & missing(group_ED) // they agreed on wrong matching
replace match_ok = 1 if group_JU_temp == group_ED_temp // agreed on correct matching
replace match_ok = 0 if missing(match_ok) // disagreement

order CID UCID_panel_3 UCID_1934_2003 year group_JU group_JU_temp group_ED group_ED_temp match_ok
replace group_ED = 0 if missing(group_ED)
replace group_JU = 0 if missing(group_JU)

gen correct_match = .
order CID UCID_panel_3 UCID_1934_2003 year group_JU group_JU_temp group_ED group_ED_temp match_ok correct_match

replace correct_match = 0 if match_ok == 1 & group_ED_temp == 0
replace correct_match = 1 if match_ok == 1 & group_ED_temp == 1
drop _merge 

merge 1:1 CID year using "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_toeval_fullnames_YS_prep"
replace group_YS = 0 if missing(group_YS) & _merge == 3
replace correct_match = 0 if group_YS == 0
replace correct_match = 1 if group_YS > 0 & _merge == 3
replace group_YS = group_ED if missing(group_YS)

drop _merge

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_evaldone.dta", replace

/*/* subsample disagreeing 100% */
preserve
	gen group_notok = 1 if match_ok == 0
		by UCID_panel_3: replace group_notok = group_notok[_n-1] if missing(group_notok)
		forval x = 1/27{
		by UCID_panel_3: replace group_notok = group_notok[_n+1] if missing(group_notok)
		}
	keep if group_notok == 1
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_eval_3.dta", replace /*Decide on this again*/
restore
*/
/* Done by Yannick : 03 to 07 november */

********************************************************************************
/*/*subsample agreed on wrong matches*/
preserve
	gen group_notok = 0 if match_ok == 0
		by UCID_panel_3: replace group_notok = group_notok[_n-1] if missing(group_notok)
		forval x = 1/27{
		by UCID_panel_3: replace group_notok = group_notok[_n+1] if missing(group_notok)
		}
	drop if group_notok == 1
	bysort UCID_panel_3: keep if match_ok == 1 & group_ED_temp == 0
	gen wrong_match = 1
	drop _merge
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_eval_1.dta", replace
restore

/*subsample agreed on matches*/
preserve
	gen group_notok = 0 if match_ok == 0
		by UCID_panel_3: replace group_notok = group_notok[_n-1] if missing(group_notok)
		forval x = 1/27{
		by UCID_panel_3: replace group_notok = group_notok[_n+1] if missing(group_notok)
		}
	drop if group_notok == 1
	bysort UCID_panel_3: keep if match_ok == 1 & group_ED_temp == 1
	gen wrong_match = 0
	drop _merge
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_eval_2.dta", replace
restore

********************************************************************************
/* Deciding on 50/50 cases */

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_eval_3.dta", clear
drop group_JU_temp group_ED_temp S T _merge group_notok
gen group_YS = group_ED
order CID UCID_panel_3 UCID_1934_2003 year group_JU group_ED match_ok group_YS
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Evaluation\match_10_nogeo_toeval_fullnames_YS.xlsx", replace firstrow(variables)

*my decision on disagreement
import excel "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_toeval_fullnames_YS_prep.xlsx", sheet("Sheet1") firstrow clear
save "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_toeval_fullnames_YS_prep", replace
*/

********************************************************************************
/*What to do:
	- We have a subsample of beginning and ending of series where we applied the previous matching
	procedure. Thus, we need to attribute different IDs to those matched togehter or not here.
	- First, take the subsample and attribute the correct IDs. The unmatched ones must have different IDs. 
	Then the ones with the same group_YS must receive the same IDs. Different than other group_YS values.
	- Second, append the two subsamples and use the new UCID_panel_3 to correctly link series together.
	- Third, attribute the same UCID_1934_2003 to the complete serie (the two previously indep ones)
	*/

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_all.dta", clear	

sort CID year
quietly by CID year:  gen PROBLEM = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
drop if PROBLEM == 2
drop PROBLEM _merge

sort UCID_panel_3 year
drop UCID_1934_2003_temp
gen UCID_1934_2003_temp = UCID_1934_2003

merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_evaldone.dta"
sort UCID_panel_3 group_YS
gen UCID_panel_3_backup = UCID_panel_3	
order CID UCID_3 UCID_1934_2003 UCID_1934_2003_temp UCID_panel_1 UCID_panel_2 UCID_panel_3_backup UCID_panel_3 group_YS	
	
egen maxUCID_panel_3 = max(UCID_panel_3) /* 1 to 1391305 */
drop maxUCID_1934_2003
egen maxUCID_1934_2003 = max(UCID_1934_2003) /* 1 to 1367837 */	
	
replace group_YS = 0 if missing(group_YS)
replace UCID_panel_3 = UCID_1934_2003 + maxUCID_1934_2003 if group_YS == 0

forval i = 1/15 {
drop maxUCID_panel_3
egen maxUCID_panel_3 = max(UCID_panel_3) 
replace UCID_panel_3 = UCID_panel_3 + maxUCID_panel_3 if group_YS == `i'
}

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_all_corrected.dta", replace
	
********************************************************************************
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_perfect_nogeo_1934_2003_ok_2.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_all_corrected.dta"

sort CID year
quietly by CID year:  gen PROBLEM = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
drop if PROBLEM == 2
drop PROBLEM _merge
drop if missing(cname)

sort UCID_panel_3 year
drop UCID_1934_2003_temp
gen UCID_1934_2003_temp = UCID_1934_2003
order CID UCID_3 UCID_1934_2003 UCID_1934_2003_temp UCID_panel_1 UCID_panel_2 UCID_panel_3

replace group_YS = 0 if missing(group_YS)

sort UCID_panel_3 year
by UCID_panel_3: replace UCID_1934_2003 = UCID_1934_2003[1] if group_YS > 0
sort UCID_1934_2003_temp year
by UCID_1934_2003_temp: replace UCID_1934_2003 = UCID_1934_2003[1]
sort UCID_1934_2003 year

drop maxUCID_1934_2003
egen maxUCID_1934_2003 = max(UCID_1934_2003) 

********************************************************************************
/*Obs with georeferencing issues*/

/*What are the different issues:

1) Missing geolocalization now matched (attribute same gdenr and regroup)
2) Branch Offices (separation)
3) Holdings (regroup under the holding or firm (?))
4) Sociétés immobilières (regroup if same gdenr_2018)
5) Littau and Luzern (simply regroup (attention: do not regroup Luzern blindly, cases fall in 6))
6) Other cause (OCR, wrong geoloc, wrong evaluation) : discuss !

*/

quietly by UCID_1934_2003 year:  gen dupdupdup = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/

preserve
	keep if dupdupdup == 0
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_ok.dta", replace
restore 

drop if dupdupdup == 0		
drop if missing(cname)

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_not_ok.dta", replace

*** 1) Missing gdenr & gdenr_2012
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_not_ok.dta", clear

sort UCID_1934_2003 year gdenr_2012

g missing_gdenr = .
replace missing_gdenr = 1 if missing(gdenr)
g miss_gdenr_group =.
replace miss_gdenr_group = 1 if missing_gdenr == 1

by UCID_1934_2003: replace miss_gdenr_group = miss_gdenr_group[_n+1] if missing(miss_gdenr_group)
by UCID_1934_2003: replace miss_gdenr_group = miss_gdenr_group[_n-1] if missing(miss_gdenr_group)
replace miss_gdenr_group = 0 if missing(miss_gdenr_group)

preserve
	keep if miss_gdenr_group == 1
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_missing.dta", replace
restore

drop if miss_gdenr_group == 1

***2) Holding
g holding_dummy = 0

replace holding_dummy = 1 if regexm(cname_orig, "Holding")
replace holding_dummy = 1 if regexm(cname_orig, "holding")
sort UCID_1934_2003 year holding_dummy

g holding_group =.
replace holding_group = 1 if holding_dummy == 1
by UCID_1934_2003: replace holding_group = holding_group[_n+1] if missing(holding_group)
by UCID_1934_2003: replace holding_group = holding_group[_n-1] if missing(holding_group)
replace holding_group= 0 if missing(holding_group)

preserve
	keep if holding_group == 1
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_holding.dta", replace
restore

drop if holding_group == 1

*** 3) Succursales
replace branch_dummy = 1 if (regexm(cname_orig), "succursale") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(cname_orig), "branch") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(cname_orig), "zweigniederlassung") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(cname_orig), "filiale")

replace branch_dummy = 1 if (regexm(gdename_orig), "succursale") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(gdename_orig), "branch") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(gdename_orig), "zweigniederlassung") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(gdename_orig), "filiale")

replace branch_dummy = 1 if (regexm(cname_orig), "Succursale") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(cname_orig), "Branch") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(cname_orig), "Zweigniederlassung") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(cname_orig), "Filiale")

replace branch_dummy = 1 if (regexm(gdename_orig), "Succursale") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(gdename_orig), "Branch") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(gdename_orig), "Zweigniederlassung") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(gdename_orig), "Zweigniededassung") // indicating if succursale in the cname
replace branch_dummy = 1 if (regexm(gdename_orig), "Filiale")

g branch_group =.
replace branch_group = 1 if branch_dummy == 1

by UCID_1934_2003: replace branch_group = branch_group[_n+1] if missing(branch_group)
by UCID_1934_2003: replace branch_group = branch_group[_n-1] if missing(branch_group)
replace branch_group = 0 if missing(branch_group)

preserve
	keep if branch_group == 1
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_branch.dta", replace
restore

drop if branch_group == 1

*** 4) Sociétés immobilières (regroup if same gdenr_2018)
g si_dummy_group=.
replace si_dummy_group = 1 if si_dummy_2 == 1

by UCID_1934_2003: replace si_dummy_group = si_dummy_group[_n+1] if missing(si_dummy_group)
by UCID_1934_2003: replace si_dummy_group = si_dummy_group[_n-1] if missing(si_dummy_group)
replace si_dummy_group = 0 if missing(si_dummy_group)

preserve
	keep if si_dummy_group == 1
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_soc_immo.dta", replace
restore

drop if si_dummy_group == 1

*** 5) Littau and Luzern
*Somehow not merged previously although this is the same municipality in 2018.

g littau_lu = .
replace littau_lu = 1 if gdenr == 1060
g littau_lu_group =.
replace littau_lu_group = 1 if littau_lu == 1

by UCID_1934_2003: replace littau_lu_group = littau_lu_group[_n+1] if missing(littau_lu_group)
by UCID_1934_2003: replace littau_lu_group = littau_lu_group[_n-1] if missing(littau_lu_group)
replace littau_lu_group = 0 if missing(littau_lu_group)

preserve
	keep if littau_lu_group == 1
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_littau_lu.dta", replace
restore

drop if littau_lu_group == 1

***6) Other cases

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_prob.dta", replace

********************************************************************************
* Regrouping / Collapsing duplicates

*** 1) Missing geolocalisation
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_missing.dta", clear
sort UCID_1934_2003 year gdenr_2018

foreach y in gdenr gdenr_2012 gdenr_2018 {
	by UCID_1934_2003: replace `y' = `y'[_n+1] if missing(`y')
	by UCID_1934_2003: replace `y' = `y'[_n-1] if missing(`y')
}
		
	*Collapse cname & Rest of variables
	foreach y in owners function signature page CID_str cname_group {
		bysort UCID_1934_2003 year gdenr_2018: replace `y' = `y'[_n-1] + "; " + `y'
		replace `y' = " "+`y'
		replace `y' = usubinstr(`y', " ; ", "",.)
		bysort UCID_1934_2003 year gdenr_2018: replace `y' = `y'[_N]
		}
		
	*Delete duplicates
		quietly by UCID_1934_2003 year gdenr_2018:  gen dupdupdupdup = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
		drop if dupdupdupdup > 1
		sort UCID_1934_2003 year
	
	*Remaining dup in this subsample
	preserve
		keep if dupdupdupdup == 0
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_missing_dup.dta", replace
	restore

	*Cleaned subsample
	drop if dupdupdupdup == 0
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_missing_dedup.dta", replace

***2) Holding
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_holding.dta", clear
sort UCID_1934_2003 year holding_dummy

replace UCID_1934_2003 = maxUCID_1934_2003 + UCID_1934_2003 if holding_dummy == 1

	*Delete duplicates
	sort UCID_1934_2003 year gdenr_2018
	quietly by UCID_1934_2003 year gdenr_2018:  gen dupdupdupdup = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
	
	*Remaining dup in this sample
	preserve
		drop if dupdupdupdup == 0
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_holding_dup.dta", replace
	restore
	
	*Cleaned subsample
	keep if dupdupdupdup == 0
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_holding_dedup.dta", replace
	
*** 3) Succursales
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_branch.dta", clear

	*Collapse cname & Rest of variables
	foreach y in owners function signature page CID_str cname_group {
		bysort UCID_1934_2003 year gdenr_2018: replace `y' = `y'[_n-1] + "; " + `y'
		replace `y' = " "+`y'
		replace `y' = usubinstr(`y', " ; ", "",.)
		bysort UCID_1934_2003 year gdenr_2018: replace `y' = `y'[_N]
		}
		
	*Delete duplicates
		quietly by UCID_1934_2003 year gdenr_2018:  gen test = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
		drop if test > 1
		sort UCID_1934_2003 year
		
	*Separate UCID for branches coexisting
		bysort UCID_1934_2003 year:  gen dupdupdupdup = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
		replace UCID_1934_2003 = maxUCID_1934_2003 + gdenr_2018 if dupdupdupdup > 1
		
		bysort UCID_1934_2003 year gdenr_2018:  gen test2 = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
		
	*Remaining dup in this sample
	preserve
		drop if test2 == 0
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_branch_dup.dta", replace
	restore

	keep if test2 == 0
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_branch_dedup.dta", replace
		
		
*** 4) Sociétés immobilières (regroup if same gdenr_2018)
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_soc_immo.dta", clear
sort UCID_1934_2003 year gdenr_2018
		
	*Collapse cname & Rest of variables
	foreach y in owners function signature page CID_str cname_group {
		bysort UCID_1934_2003 year gdenr_2018: replace `y' = `y'[_n-1] + "; " + `y'
		replace `y' = " "+`y'
		replace `y' = usubinstr(`y', " ; ", "",.)
		bysort UCID_1934_2003 year gdenr_2018: replace `y' = `y'[_N]
		}
		
	*Delete duplicates
		quietly by UCID_1934_2003 year gdenr_2018:  gen test = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
		drop if test > 1
		sort UCID_1934_2003 year
		
	*Remaining dup in this subsample
		quietly by UCID_1934_2003 year:  gen test2 = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
	preserve
		drop if test2 == 0
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_soc_immo_dup.dta", replace
	restore

keep if test2 == 0
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_soc_immo_dedup.dta", replace

*** 5) Luzern and Littau
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_littau_lu.dta", clear
sort UCID_1934_2003 year gdenr_2018
		
	*Collapse cname & Rest of variables
	foreach y in owners function signature page CID_str cname_group {
		bysort UCID_1934_2003 year gdenr_2018: replace `y' = `y'[_n-1] + "; " + `y'
		replace `y' = " "+`y'
		replace `y' = usubinstr(`y', " ; ", "",.)
		bysort UCID_1934_2003 year gdenr_2018: replace `y' = `y'[_N]
		}
		
	*Delete duplicates
		quietly by UCID_1934_2003 year gdenr_2018:  gen test = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
		drop if test > 1
		sort UCID_1934_2003 year
		
	*Remaining dup in this subsample
		quietly by UCID_1934_2003 year:  gen test2 = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/
	preserve
		drop if test2 == 0
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_littau_lu_dup.dta", replace
	restore
	
	keep if test2 == 0
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_littau_lu_dedup.dta", replace

*** 6) Best of the rest
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_prob.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_littau_lu_dup.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_soc_immo_dup.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_branch_dup.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_holding_dup.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_missing_dup.dta"

g checked = .
order CID checked UCID_1934_2003 UCID_panel_3 year cname cname_orig gdenr gdenr_2012 gdenr_2018 gdename gdename_orig maxyear maxyear_2 maxyear_3 cname_group owners
sort UCID_1934_2003 year maxyear_3 gdenr_2018

replace gdename ="arbon" if UCID_1934_2003 == 241377
replace gdenr = 4401 if UCID_1934_2003 == 241377
replace gdenr_2012 = 4401 if UCID_1934_2003 == 241377
replace gdenr_2018 = 4401 if UCID_1934_2003 == 241377
replace checked = 1 if UCID_1934_2003 == 241377

replace gdename ="tesserete" if UCID_1934_2003 == 401654 
replace gdenr = 5226 if UCID_1934_2003 == 401654 
replace gdenr_2012 = 5226 if UCID_1934_2003 == 401654 
replace gdenr_2018 = 5226 if UCID_1934_2003 == 401654 
replace checked = 1 if UCID_1934_2003 ==401654

replace gdename ="bex" if UCID_1934_2003 == 413524 
replace gdenr = 5402 if UCID_1934_2003 == 413524 
replace gdenr_2012 = 5402 if UCID_1934_2003 == 413524 
replace gdenr_2018 = 5402 if UCID_1934_2003 == 413524 
replace checked = 1 if UCID_1934_2003 ==413524

replace gdename ="givisiez" if UCID_1934_2003 == 227992
replace gdenr = 2197 if UCID_1934_2003 == 227992
replace gdenr_2012 = 2197 if UCID_1934_2003 == 227992
replace gdenr_2018 = 2197 if UCID_1934_2003 == 227992
replace checked = 1 if UCID_1934_2003 ==227992

replace gdename ="luzern" if UCID_1934_2003 == 144587
replace gdenr = 1061 if UCID_1934_2003 == 144587
replace gdenr_2012 = 1061 if UCID_1934_2003 == 144587
replace gdenr_2018 = 1061 if UCID_1934_2003 == 144587
replace checked = 1 if UCID_1934_2003 ==144587

replace gdename ="horw" if UCID_1934_2003 == 135899
replace gdenr = 1058 if UCID_1934_2003 == 135899
replace gdenr_2012 = 1058 if UCID_1934_2003 == 135899
replace gdenr_2018 = 1058 if UCID_1934_2003 == 135899
replace checked = 1 if UCID_1934_2003 ==135899

replace gdename ="buchrain" if UCID_1934_2003 == 134908
replace gdenr = 1052 if UCID_1934_2003 == 134908
replace gdenr_2012 = 1052 if UCID_1934_2003 == 134908
replace gdenr_2018 = 1052 if UCID_1934_2003 == 134908
replace checked = 1 if UCID_1934_2003 ==134908

replace gdename ="baden" if UCID_1934_2003 == 133987
replace gdenr = 4021 if UCID_1934_2003 == 133987
replace gdenr_2012 = 4021 if UCID_1934_2003 == 133987
replace gdenr_2018 = 4021 if UCID_1934_2003 == 133987
replace checked = 1 if UCID_1934_2003 ==133987

replace gdename ="hasle lu" if UCID_1934_2003 == 132544
replace gdenr = 1005 if UCID_1934_2003 == 132544
replace gdenr_2012 = 1005 if UCID_1934_2003 == 132544
replace gdenr_2018 = 1005 if UCID_1934_2003 == 132544
replace checked = 1 if UCID_1934_2003 ==132544

replace gdename ="emmen" if UCID_1934_2003 == 133174
replace gdenr = 1024 if UCID_1934_2003 == 133174
replace gdenr_2012 = 1024 if UCID_1934_2003 == 133174
replace gdenr_2018 = 1024 if UCID_1934_2003 == 133174
replace checked = 1 if UCID_1934_2003 ==133174

replace gdename ="thoerigen" if UCID_1934_2003 == 132176
replace gdenr = 989 if UCID_1934_2003 == 132176
replace gdenr_2012 = 989 if UCID_1934_2003 == 132176
replace gdenr_2018 = 989 if UCID_1934_2003 == 132176
replace checked = 1 if UCID_1934_2003 ==132176

replace gdename ="weesen" if UCID_1934_2003 == 131334
replace gdenr = 3316 if UCID_1934_2003 == 131334
replace gdenr_2012 = 3316 if UCID_1934_2003 == 131334
replace gdenr_2018 = 3316 if UCID_1934_2003 == 131334
replace checked = 1 if UCID_1934_2003 ==131334

replace gdename ="biel (be)" if UCID_1934_2003 == 127777
replace gdenr = 371 if UCID_1934_2003 == 127777
replace gdenr_2012 = 371 if UCID_1934_2003 == 127777
replace gdenr_2018 = 371 if UCID_1934_2003 ==  127777
replace checked = 1 if UCID_1934_2003 ==127777

replace gdename ="saanen" if UCID_1934_2003 == 125989
replace gdenr = 843 if UCID_1934_2003 == 125989
replace gdenr_2012 = 843 if UCID_1934_2003 == 125989
replace gdenr_2018 = 843 if UCID_1934_2003 ==  125989
replace checked = 1 if UCID_1934_2003 ==125989

replace gdename ="saanen" if UCID_1934_2003 == 126412
replace gdenr = 843 if UCID_1934_2003 == 126412
replace gdenr_2012 = 843 if UCID_1934_2003 == 126412
replace gdenr_2018 = 843 if UCID_1934_2003 ==  125123
replace checked = 1 if UCID_1934_2003 ==125123

replace gdename ="oberwil im simmental" if UCID_1934_2003 == 125123
replace gdenr = 766 if UCID_1934_2003 == 125123
replace gdenr_2012 = 766 if UCID_1934_2003 == 125123
replace gdenr_2018 = 766 if UCID_1934_2003 == 125123
replace checked = 1 if UCID_1934_2003 ==125123

replace gdename ="nidau" if UCID_1934_2003 == 124403
replace gdenr = 743 if UCID_1934_2003 == 124403
replace gdenr_2012 = 743 if UCID_1934_2003 == 124403
replace gdenr_2018 = 743 if UCID_1934_2003 == 124403
replace checked = 1 if UCID_1934_2003 == 124403

replace gdename ="bruegg" if UCID_1934_2003 == 123678
replace gdenr = 733 if UCID_1934_2003 == 123678
replace gdenr_2012 = 733 if UCID_1934_2003 == 123678
replace gdenr_2018 = 733 if UCID_1934_2003 == 123678
replace checked = 1 if UCID_1934_2003 == 123678

replace gdename ="arth" if UCID_1934_2003 == 114138
replace gdenr = 1362 if UCID_1934_2003 == 114138
replace gdenr_2012 = 1362 if UCID_1934_2003 == 114138
replace gdenr_2018 = 1362 if UCID_1934_2003 == 114138
replace checked = 1 if UCID_1934_2003 == 114138

replace gdename ="ittigen" if UCID_1934_2003 == 109305
replace gdenr = 362 if UCID_1934_2003 == 109305
replace gdenr_2012 = 362 if UCID_1934_2003 == 109305
replace gdenr_2018 = 362 if UCID_1934_2003 == 109305
replace checked = 1 if UCID_1934_2003 == 109305
    
replace gdename ="horgen" if UCID_1934_2003 == 91180
replace gdenr = 133 if UCID_1934_2003 == 91180
replace gdenr_2012 = 133 if UCID_1934_2003 == 91180
replace gdenr_2018 = 295 if UCID_1934_2003 ==  91180
replace checked = 1 if UCID_1934_2003 == 91180

replace gdename ="elgg" if UCID_1934_2003 == 90959
replace gdenr = 217 if UCID_1934_2003 == 90959
replace gdenr_2012 = 217 if UCID_1934_2003 == 90959
replace gdenr_2018 = 294 if UCID_1934_2003 == 90959
replace checked = 1 if UCID_1934_2003 == 90959
 
replace gdename ="lotzwil" if UCID_1934_2003 == 94805
replace gdenr = 331 if UCID_1934_2003 == 94805
replace gdenr_2012 = 331 if UCID_1934_2003 == 94805
replace gdenr_2018 = 331 if UCID_1934_2003 == 94805
replace checked = 1 if UCID_1934_2003 == 94805
 
replace gdename ="murten" if UCID_1934_2003 == 102653
replace gdenr = 2275 if UCID_1934_2003 == 102653
replace gdenr_2012 = 2275 if UCID_1934_2003 == 102653
replace gdenr_2018 = 2275 if UCID_1934_2003 == 102653
replace checked = 1 if UCID_1934_2003 == 102653

replace gdename ="biel (be)" if UCID_1934_2003 == 104780
replace gdenr = 371 if UCID_1934_2003 == 104780
replace gdenr_2012 = 371 if UCID_1934_2003 == 104780
replace gdenr_2018 = 371 if UCID_1934_2003 == 104780
replace checked = 1 if UCID_1934_2003 == 104780

replace gdename ="ittigen" if UCID_1934_2003 == 109171
replace gdenr = 362 if UCID_1934_2003 == 109171
replace gdenr_2012 = 362 if UCID_1934_2003 == 109171
replace gdenr_2018 = 362 if UCID_1934_2003 == 109171 
replace checked = 1 if UCID_1934_2003 == 109171

replace gdename ="winterthur" if UCID_1934_2003 == 57822
replace gdenr = 230 if UCID_1934_2003 == 57822
replace gdenr_2012 = 230 if UCID_1934_2003 == 57822
replace gdenr_2018 = 230 if UCID_1934_2003 == 57822 
replace checked = 1 if UCID_1934_2003 == 57822

replace gdename ="interlaken" if UCID_1934_2003 == 59241
replace gdenr = 581 if UCID_1934_2003 == 59241
replace gdenr_2012 = 581 if UCID_1934_2003 == 59241
replace gdenr_2018 = 581 if UCID_1934_2003 == 59241
replace checked = 1 if UCID_1934_2003 == 59241

replace gdename ="udligenswil" if UCID_1934_2003 == 66512
replace gdenr = 1067 if UCID_1934_2003 == 66512
replace gdenr_2012 = 1067 if UCID_1934_2003 == 66512
replace gdenr_2018 = 1067 if UCID_1934_2003 == 66512
replace checked = 1 if UCID_1934_2003 == 66512
 
replace gdename ="zuerich" if UCID_1934_2003 == 80551
replace gdenr = 261 if UCID_1934_2003 == 80551 & year > 1989
replace gdenr = 253 if UCID_1934_2003 == 80551 & year < 1990
replace gdenr_2012 = 261 if UCID_1934_2003 == 80551
replace gdenr_2018 = 261 if UCID_1934_2003 ==  80551
replace checked = 1 if UCID_1934_2003 == 80551

replace gdename ="daellikon" if UCID_1934_2003 == 81511
replace gdenr = 84 if UCID_1934_2003 == 81511
replace gdenr_2012 = 84 if UCID_1934_2003 == 81511
replace gdenr_2018 = 84 if UCID_1934_2003 ==  81511
replace checked = 1 if UCID_1934_2003 == 81511

replace gdename ="dietikon" if UCID_1934_2003 == 37677
replace gdenr = 54 if UCID_1934_2003 == 37677
replace gdenr_2012 = 54 if UCID_1934_2003 == 37677
replace gdenr_2018 = 54 if UCID_1934_2003 ==  37677
replace checked = 1 if UCID_1934_2003 ==  37677

replace gdename ="hoeri" if UCID_1934_2003 == 4185
replace gdenr = 60 if UCID_1934_2003 == 4185
replace gdenr_2012 = 60 if UCID_1934_2003 == 4185
replace gdenr_2018 = 60 if UCID_1934_2003 == 4185
replace checked = 1 if UCID_1934_2003 == 4185 

replace gdename ="benken zh" if UCID_1934_2003 == 1303
replace gdenr = 22 if UCID_1934_2003 == 1303
replace gdenr_2012 = 22 if UCID_1934_2003 == 1303
replace gdenr_2018 = 22 if UCID_1934_2003 ==  1303
replace checked = 1 if UCID_1934_2003 == 1303

replace gdename ="affoltern am albis" if UCID_1934_2003 == 294
replace gdenr = 2 if UCID_1934_2003 == 294
replace gdenr_2012 = 2 if UCID_1934_2003 == 294
replace gdenr_2018 = 2 if UCID_1934_2003 ==  294
replace checked = 1 if UCID_1934_2003 == 294
  
replace gdename ="pully" if UCID_1934_2003 == 1362929
replace gdenr = 5590 if UCID_1934_2003 == 1362929
replace gdenr_2012 = 5590 if UCID_1934_2003 == 1362929
replace gdenr_2018 = 5590 if UCID_1934_2003 ==  1362929
replace checked = 1 if UCID_1934_2003 == 1362929

replace gdename ="lausanne" if UCID_1934_2003 == 1362744
replace gdenr = 5586 if UCID_1934_2003 == 1362744
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1362744
replace gdenr_2018 = 5586 if UCID_1934_2003 ==  1362744
replace checked = 1 if UCID_1934_2003 == 1362744

replace gdename ="montana" if UCID_1934_2003 == 1342562
replace gdenr = 6243 if UCID_1934_2003 == 1342562
replace gdenr_2012 = 6243 if UCID_1934_2003 == 1342562
replace gdenr_2018 = 6253 if UCID_1934_2003 ==  1342562
replace checked = 1 if UCID_1934_2003 == 1342562

replace gdename ="chermignon" if UCID_1934_2003 == 1342522
replace gdenr = 6234 if UCID_1934_2003 == 1342522
replace gdenr_2012 = 6234 if UCID_1934_2003 == 1342522
replace gdenr_2018 = 6253 if UCID_1934_2003 ==  1342522
replace checked = 1 if UCID_1934_2003 == 1342522

replace gdename ="zuerich" if UCID_1934_2003 == 1332022
replace gdenr = 261 if UCID_1934_2003 == 1332022 & year > 1989
replace gdenr = 253 if UCID_1934_2003 == 1332022 & year < 1990
replace gdenr_2012 = 261 if UCID_1934_2003 == 1332022
replace gdenr_2018 = 261 if UCID_1934_2003 ==  1332022
replace checked = 1 if UCID_1934_2003 == 1332022

replace gdename ="ligornetto" if UCID_1934_2003 == 1322859
replace gdenr = 5253 if UCID_1934_2003 == 1322859
replace gdenr_2012 = 5253 if UCID_1934_2003 == 1322859
replace gdenr_2018 = 5254 if UCID_1934_2003 ==  1322859
replace checked = 1 if UCID_1934_2003 == 1322859

replace gdename ="muehlebach" if UCID_1934_2003 == 1126663
replace gdenr = 6062 if UCID_1934_2003 == 1126663
replace gdenr_2012 = 6056 if UCID_1934_2003 == 1126663
replace gdenr_2018 = 6056 if UCID_1934_2003 ==  1126663
replace checked = 1 if UCID_1934_2003 == 1126663

replace gdename ="langnau bei reiden" if UCID_1934_2003 == 1089434
replace gdenr = 1134 if UCID_1934_2003 == 1089434
replace gdenr_2012 = 1140 if UCID_1934_2003 == 1089434
replace gdenr_2018 = 1140 if UCID_1934_2003 ==  1089434
replace checked = 1 if UCID_1934_2003 == 1089434

replace gdename ="zuerich" if UCID_1934_2003 == 935795
replace gdenr = 253 if UCID_1934_2003 == 935795 & year < 1989
replace gdenr = 261 if UCID_1934_2003 == 935795 & year > 1990
replace gdenr_2012 = 261 if UCID_1934_2003 == 935795
replace gdenr_2018 = 261 if UCID_1934_2003 ==  935795
replace checked = 1 if UCID_1934_2003 == 935795

replace gdename ="st antoenien" if UCID_1934_2003 == 902582
replace gdenr = 3893 if UCID_1934_2003 == 902582
replace gdenr_2012 = 3893 if UCID_1934_2003 == 902582
replace gdenr_2018 = 3891 if UCID_1934_2003 ==  902582
replace checked = 1 if UCID_1934_2003 == 902582

replace gdename ="muehlebach" if UCID_1934_2003 == 713533
replace gdenr = 6062 if UCID_1934_2003 == 713533
replace gdenr_2012 = 6056 if UCID_1934_2003 == 713533
replace gdenr_2018 = 6056 if UCID_1934_2003 == 713533
replace checked = 1 if UCID_1934_2003 == 713533
 
replace gdename ="sitterdorf" if UCID_1934_2003 == 680465
replace gdenr = 4518 if UCID_1934_2003 == 680465
replace gdenr_2012 = 4511 if UCID_1934_2003 == 680465
replace gdenr_2018 = 4511 if UCID_1934_2003 ==  680465
replace checked = 1 if UCID_1934_2003 == 680465

replace gdename ="prez vers siviriez" if UCID_1934_2003 == 628908
replace gdenr = 2094 if UCID_1934_2003 == 628908
replace gdenr_2012 = 2099 if UCID_1934_2003 == 628908
replace gdenr_2018 = 2099 if UCID_1934_2003 == 628908
replace checked = 1 if UCID_1934_2003 == 628908
 
replace gdename ="les friques" if UCID_1934_2003 == 628483
replace gdenr = 2021 if UCID_1934_2003 == 628483
replace gdenr_2012 = 2041 if UCID_1934_2003 == 628483
replace gdenr_2018 = 2041 if UCID_1934_2003 ==  628483
replace checked = 1 if UCID_1934_2003 == 628483

replace gdename ="wohlen ag" if UCID_1934_2003 == 549063
replace gdenr = 4082 if UCID_1934_2003 == 549063
replace gdenr_2012 = 4082 if UCID_1934_2003 == 549063
replace gdenr_2018 = 4082 if UCID_1934_2003 ==  549063
replace checked = 1 if UCID_1934_2003 == 549063
 
replace gdename ="wohlen ag" if UCID_1934_2003 == 546399
replace gdenr = 4082 if UCID_1934_2003 == 546399
replace gdenr_2012 = 4082 if UCID_1934_2003 == 546399
replace gdenr_2018 = 4082 if UCID_1934_2003 ==  546399  
replace checked = 1 if UCID_1934_2003 == 546399

replace gdename ="wohlen ag" if UCID_1934_2003 == 544208
replace gdenr = 4082 if UCID_1934_2003 == 544208
replace gdenr_2012 = 4082 if UCID_1934_2003 == 544208
replace gdenr_2018 = 4082 if UCID_1934_2003 ==  544208
replace checked = 1 if UCID_1934_2003 == 544208

replace gdename ="wohlen ag" if UCID_1934_2003 == 538542
replace gdenr = 4082 if UCID_1934_2003 == 538542
replace gdenr_2012 = 4082 if UCID_1934_2003 == 538542
replace gdenr_2018 = 4082 if UCID_1934_2003 ==  538542
replace checked = 1 if UCID_1934_2003 == 538542

*Tag wrong matches
g separate = .

replace separate = 1 if UCID_1934_2003 ==70484
replace separate = 1 if UCID_1934_2003 ==106199
replace separate = 1 if UCID_1934_2003 ==108515
replace separate = 1 if UCID_1934_2003 ==106637
replace separate = 1 if UCID_1934_2003 ==110216
replace separate = 1 if UCID_1934_2003 ==127074
replace separate = 1 if UCID_1934_2003 ==127754

replace checked = 1 if UCID_1934_2003 ==70484
replace checked = 1 if UCID_1934_2003 ==106199
replace checked = 1 if UCID_1934_2003 ==108515
replace checked = 1 if UCID_1934_2003 ==106637
replace checked = 1 if UCID_1934_2003 ==110216
replace checked = 1 if UCID_1934_2003 ==127074
replace checked = 1 if UCID_1934_2003 ==127754

*Correct, just collapse
replace separate = 0 if UCID_1934_2003 == 1341450
replace separate = 0 if UCID_1934_2003 == 1290070

replace checked = 1 if UCID_1934_2003 == 1341450
replace checked = 1 if UCID_1934_2003 == 1290070

*no idea

/*
export for evaluation
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Evaluation\duplicates_ED.xlsx", replace firstrow(variables)
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Evaluation\duplicates_JU.xlsx", replace firstrow(variables)
restore
*/

********************************************************************************
* Corrections following ED and Ju

*****
*Julie
*Corrections
replace gdename ="buelach" if UCID_1934_2003 == 2603
replace gdenr = 53 if UCID_1934_2003 == 2603
replace gdenr_2012 = 53 if UCID_1934_2003 == 2603
replace gdenr_2018 = 53 if UCID_1934_2003 ==  2603

replace gdename ="buelach" if UCID_1934_2003 == 2835
replace gdenr = 53 if UCID_1934_2003 == 2835
replace gdenr_2012 = 53 if UCID_1934_2003 == 2835
replace gdenr_2018 = 53 if UCID_1934_2003 ==  2835

replace gdename ="dietlikon" if UCID_1934_2003 == 3489
replace gdenr = 54 if UCID_1934_2003 == 3489
replace gdenr_2012 = 54 if UCID_1934_2003 == 3489
replace gdenr_2018 = 54 if UCID_1934_2003 ==  3489

replace gdename ="schuebelbach" if UCID_1934_2003 == 3611
replace gdenr = 1346 if UCID_1934_2003 == 3611
replace gdenr_2012 = 1346 if UCID_1934_2003 == 3611
replace gdenr_2018 = 1346 if UCID_1934_2003 == 3611

replace gdename ="hoeri" if UCID_1934_2003 == 4106
replace gdenr = 60 if UCID_1934_2003 == 4106
replace gdenr_2012 = 60 if UCID_1934_2003 == 4106
replace gdenr_2018 = 60 if UCID_1934_2003 == 4106

replace gdename ="kloten" if UCID_1934_2003 == 4440
replace gdenr = 62 if UCID_1934_2003 == 4440
replace gdenr_2012 = 62 if UCID_1934_2003 == 4440
replace gdenr_2018 = 62 if UCID_1934_2003 == 4440

replace gdename ="kloten" if UCID_1934_2003 == 4641
replace gdenr = 62 if UCID_1934_2003 == 4641
replace gdenr_2012 = 62 if UCID_1934_2003 == 4641
replace gdenr_2018 = 62 if UCID_1934_2003 == 4641

replace gdename ="kloten" if UCID_1934_2003 == 4828
replace gdenr = 62 if UCID_1934_2003 == 4828
replace gdenr_2012 = 62 if UCID_1934_2003 == 4828
replace gdenr_2018 = 62 if UCID_1934_2003 == 4828

replace gdename ="kloten" if UCID_1934_2003 == 6226
replace gdenr = 62 if UCID_1934_2003 == 6226
replace gdenr_2012 = 62 if UCID_1934_2003 == 6226
replace gdenr_2018 = 62 if UCID_1934_2003 == 6226

replace gdename ="opfikon" if UCID_1934_2003 == 6546
replace gdenr = 66 if UCID_1934_2003 ==6546
replace gdenr_2012 = 66 if UCID_1934_2003 ==6546
replace gdenr_2018 = 66 if UCID_1934_2003 ==6546

replace gdename ="dietikon" if UCID_1934_2003 == 8051
replace gdenr = 243 if UCID_1934_2003 == 8051
replace gdenr_2012 = 243 if UCID_1934_2003 == 8051
replace gdenr_2018 = 243 if UCID_1934_2003 == 8051

replace gdename ="dietlikon" if UCID_1934_2003 == 8175
replace gdenr = 54 if UCID_1934_2003 == 8175
replace gdenr_2012 = 54 if UCID_1934_2003 == 8175
replace gdenr_2018 = 54 if UCID_1934_2003 == 8175

replace gdename ="regensdorf" if UCID_1934_2003 == 9805
replace gdenr = 96 if UCID_1934_2003 == 9805
replace gdenr_2012 = 96 if UCID_1934_2003 == 9805
replace gdenr_2018 = 96 if UCID_1934_2003 == 9805

replace gdename ="schoefflisdorf" if UCID_1934_2003 == 12048
replace gdenr = 99 if UCID_1934_2003 == 12048
replace gdenr_2012 = 99 if UCID_1934_2003 == 12048
replace gdenr_2018 = 99 if UCID_1934_2003 == 12048

replace gdename ="gossau zh" if UCID_1934_2003 == 12975
replace gdenr = 115 if UCID_1934_2003 == 12975
replace gdenr_2012 = 115 if UCID_1934_2003 == 12975
replace gdenr_2018 = 115 if UCID_1934_2003 == 12975

replace gdename ="mezzovico vira" if UCID_1934_2003 == 15781
replace gdenr = 5199 if UCID_1934_2003 == 15781
replace gdenr_2012 = 5199 if UCID_1934_2003 == 15781
replace gdenr_2018 = 5199 if UCID_1934_2003 == 15781

replace gdename ="hallwil" if UCID_1934_2003 == 18311
replace gdenr = 4197 if UCID_1934_2003 == 18311
replace gdenr_2012 = 4197 if UCID_1934_2003 == 18311
replace gdenr_2018 = 4197 if UCID_1934_2003 == 18311

replace gdename ="herrliberg" if UCID_1934_2003 == 19660
replace gdenr = 152 if UCID_1934_2003 == 19660
replace gdenr_2012 = 152 if UCID_1934_2003 == 19660
replace gdenr_2018 = 152 if UCID_1934_2003 == 19660
 
replace gdename ="kuessnacht am rigi" if UCID_1934_2003 == 20311
replace gdenr = 1331 if UCID_1934_2003 == 20311
replace gdenr_2012 = 1331 if UCID_1934_2003 == 20311
replace gdenr_2018 = 1331 if UCID_1934_2003 == 20311

replace gdename ="zug" if UCID_1934_2003 == 20835
replace gdenr = 1711 if UCID_1934_2003 == 20835
replace gdenr_2012 = 1711 if UCID_1934_2003 == 20835
replace gdenr_2018 = 1711 if UCID_1934_2003 == 20835

replace gdename ="dielsdorf" if UCID_1934_2003 == 22564
replace gdenr = 86 if UCID_1934_2003 == 22564
replace gdenr_2012 = 86 if UCID_1934_2003 == 22564
replace gdenr_2018 = 86 if UCID_1934_2003 == 22564

replace gdename ="staefa" if UCID_1934_2003 == 22939
replace gdenr = 158 if UCID_1934_2003 == 22939
replace gdenr_2012 = 158 if UCID_1934_2003 == 22939
replace gdenr_2018 = 158 if UCID_1934_2003 == 22939

replace gdename ="uetikon" if UCID_1934_2003 == 23237
replace gdenr = 159 if UCID_1934_2003 == 23237
replace gdenr_2012 = 159 if UCID_1934_2003 == 23237
replace gdenr_2018 = 159 if UCID_1934_2003 == 23237

replace gdename ="schlieren" if UCID_1934_2003 == 24602
replace gdenr = 247 if UCID_1934_2003 == 24602
replace gdenr_2012 = 247 if UCID_1934_2003 == 24602
replace gdenr_2018 = 247 if UCID_1934_2003 == 24602
			
replace gdename ="zuerich" if UCID_1934_2003 == 24654
replace gdenr = 253 if UCID_1934_2003 == 24654 & year < 1990 
replace gdenr = 261 if UCID_1934_2003 == 24654 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 24654
replace gdenr_2018 = 261 if UCID_1934_2003 == 24654

replace gdename ="wila" if UCID_1934_2003 == 26581
replace gdenr = 181 if UCID_1934_2003 == 26581
replace gdenr_2012 = 181 if UCID_1934_2003 == 26581
replace gdenr_2018 = 181 if UCID_1934_2003 == 26581

replace gdename ="schwerzenbach" if UCID_1934_2003 == 28637
replace gdenr = 197 if UCID_1934_2003 == 28637
replace gdenr_2012 = 197 if UCID_1934_2003 == 28637
replace gdenr_2018 = 197 if UCID_1934_2003 == 28637

replace gdename ="duebendorf" if UCID_1934_2003 == 31704
replace gdenr = 191 if UCID_1934_2003 == 31704
replace gdenr_2012 = 191 if UCID_1934_2003 == 31704
replace gdenr_2018 = 191 if UCID_1934_2003 == 31704

replace gdename ="volketswil" if UCID_1934_2003 == 31751
replace gdenr = 199 if UCID_1934_2003 == 31751
replace gdenr_2012 = 199 if UCID_1934_2003 == 31751
replace gdenr_2018 = 199 if UCID_1934_2003 == 31751

replace gdename ="grossandelfingen" if UCID_1934_2003 == 33286
replace gdenr = 30 if UCID_1934_2003 == 33286
replace gdenr_2012 = 30 if UCID_1934_2003 == 33286
replace gdenr_2018 = 30 if UCID_1934_2003 == 33286
		
replace gdename ="zuerich" if UCID_1934_2003 == 33486
replace gdenr = 253 if UCID_1934_2003 == 33486 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 33486 & year > 1989 
replace gdenr_2012 = 261 if UCID_1934_2003 == 33486
replace gdenr_2018 = 261 if UCID_1934_2003 == 33486

replace gdename ="winterthur" if UCID_1934_2003 == 35184
replace gdenr = 230 if UCID_1934_2003 == 35184
replace gdenr_2012 = 230 if UCID_1934_2003 == 35184
replace gdenr_2018 = 230 if UCID_1934_2003 == 35184

replace gdename ="dietikon" if UCID_1934_2003 == 37440
replace gdenr = 243 if UCID_1934_2003 == 37440
replace gdenr_2012 = 243 if UCID_1934_2003 == 37440
replace gdenr_2018 = 243 if UCID_1934_2003 == 37440
				
replace gdename ="binningen" if UCID_1934_2003 == 37511
replace gdenr = 2806 if UCID_1934_2003 == 37511
replace gdenr_2012 = 2765 if UCID_1934_2003 == 37511
replace gdenr_2018 = 2765 if UCID_1934_2003 == 37511

replace gdename ="dietikon" if UCID_1934_2003 == 37677
replace gdenr = 243 if UCID_1934_2003 == 37677
replace gdenr_2012 = 243 if UCID_1934_2003 == 37677
replace gdenr_2018 = 243 if UCID_1934_2003 == 37677

replace gdename ="dierikon" if UCID_1934_2003 == 37706
replace gdenr = 1053 if UCID_1934_2003 == 37706
replace gdenr_2012 = 1053 if UCID_1934_2003 == 37706
replace gdenr_2018 = 1053 if UCID_1934_2003 == 37706
		
replace gdename ="dietlikon" if UCID_1934_2003 == 38842
replace gdenr = 54 if UCID_1934_2003 == 38842
replace gdenr_2012 = 54 if UCID_1934_2003 == 38842
replace gdenr_2018 = 54 if UCID_1934_2003 == 38842

replace gdename ="urdorf" if UCID_1934_2003 == 41294
replace gdenr = 250 if UCID_1934_2003 == 41294
replace gdenr_2012 = 250 if UCID_1934_2003 == 41294
replace gdenr_2018 = 250 if UCID_1934_2003 == 41294
	
replace gdename ="baden" if UCID_1934_2003 == 43676
replace gdenr = 4021 if UCID_1934_2003 == 43676
replace gdenr_2012 = 4021 if UCID_1934_2003 == 43676
replace gdenr_2018 = 4021 if UCID_1934_2003 == 43676

replace gdename ="baden" if UCID_1934_2003 == 44041
replace gdenr = 4021 if UCID_1934_2003 == 44041
replace gdenr_2012 = 4021 if UCID_1934_2003 == 44041
replace gdenr_2018 = 4021 if UCID_1934_2003 == 44041
	
replace gdename ="biel (be)" if UCID_1934_2003 == 47230
replace gdenr = 371 if UCID_1934_2003 == 47230
replace gdenr_2012 = 371 if UCID_1934_2003 == 47230
replace gdenr_2018 = 371 if UCID_1934_2003 == 47230

replace gdename ="urdorf" if UCID_1934_2003 == 48658
replace gdenr = 250 if UCID_1934_2003 == 48658
replace gdenr_2012 = 250 if UCID_1934_2003 == 48658
replace gdenr_2018 = 250 if UCID_1934_2003 == 48658
 	
replace gdename ="meilen" if UCID_1934_2003 == 51047
replace gdenr = 156 if UCID_1934_2003 == 51047
replace gdenr_2012 = 156 if UCID_1934_2003 == 51047
replace gdenr_2018 = 156 if UCID_1934_2003 == 51047

replace gdename ="zuerich" if UCID_1934_2003 == 51985
replace gdenr = 253 if UCID_1934_2003 == 51985 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 51985 & year > 1989 
replace gdenr_2012 = 261 if UCID_1934_2003 == 51985
replace gdenr_2018 = 261 if UCID_1934_2003 == 51985

replace gdename ="dietlikon" if UCID_1934_2003 == 52895
replace gdenr = 54 if UCID_1934_2003 == 52895
replace gdenr_2012 = 54 if UCID_1934_2003 == 52895
replace gdenr_2018 = 54 if UCID_1934_2003 == 52895

replace gdename ="zuerich" if UCID_1934_2003 == 56587
replace gdenr = 253 if UCID_1934_2003 == 56587 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 56587 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 56587
replace gdenr_2018 = 261 if UCID_1934_2003 == 56587

replace gdename ="wolhusen" if UCID_1934_2003 == 56858
replace gdenr = 1107 if UCID_1934_2003 == 56858
replace gdenr_2012 = 1107 if UCID_1934_2003 == 56858
replace gdenr_2018 = 1107 if UCID_1934_2003 == 56858

replace gdename ="zuerich" if UCID_1934_2003 == 57146
replace gdenr = 253 if UCID_1934_2003 == 57146 < 1990
replace gdenr = 261 if UCID_1934_2003 == 57146 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 57146
replace gdenr_2018 = 261 if UCID_1934_2003 == 57146

replace gdename ="interlaken" if UCID_1934_2003 == 59241
replace gdenr = 581 if UCID_1934_2003 == 59241
replace gdenr_2012 = 581 if UCID_1934_2003 == 59241
replace gdenr_2018 = 581 if UCID_1934_2003 == 59241

replace gdename ="kriegstetten" if UCID_1934_2003 == 67210
replace gdenr = 2525 if UCID_1934_2003 == 67210
replace gdenr_2012 = 2525 if UCID_1934_2003 == 67210
replace gdenr_2018 = 2525 if UCID_1934_2003 == 67210

replace gdename ="grenchen" if UCID_1934_2003 == 67294
replace gdenr = 2546 if UCID_1934_2003 == 67294
replace gdenr_2012 = 2546 if UCID_1934_2003 == 67294
replace gdenr_2018 = 2546 if UCID_1934_2003 == 67294

replace gdename ="zuerich" if UCID_1934_2003 == 69333
replace gdenr = 253 if UCID_1934_2003 == 69333 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 69333 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 69333
replace gdenr_2018 = 261 if UCID_1934_2003 == 69333

replace gdename ="geneve" if UCID_1934_2003 == 69795
replace gdenr = 6621 if UCID_1934_2003 == 69795
replace gdenr_2012 = 6621 if UCID_1934_2003 == 69795
replace gdenr_2018 = 6621 if UCID_1934_2003 == 69795

replace gdename ="zuerich" if UCID_1934_2003 == 73672
replace gdenr = 253 if UCID_1934_2003 == 73672 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 73672 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 73672
replace gdenr_2018 = 261 if UCID_1934_2003 == 73672 

replace gdename ="vilters" if UCID_1934_2003 == 77171
replace gdenr = 3297 if UCID_1934_2003 == 77171
replace gdenr_2012 = 3297 if UCID_1934_2003 == 77171
replace gdenr_2018 = 3297 if UCID_1934_2003 == 77171

replace gdename ="zuerich" if UCID_1934_2003 == 80551
replace gdenr = 253 if UCID_1934_2003 == 80551 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 80551 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 80551
replace gdenr_2018 = 261 if UCID_1934_2003 == 80551
	
replace gdename ="zug" if UCID_1934_2003 == 84230
replace gdenr = 1711 if UCID_1934_2003 == 84230
replace gdenr_2012 = 1711 if UCID_1934_2003 == 84230
replace gdenr_2018 = 1711 if UCID_1934_2003 == 84230
	
replace gdename ="ebikon" if UCID_1934_2003 == 84475
replace gdenr = 1054 if UCID_1934_2003 == 84475
replace gdenr_2012 = 1054 if UCID_1934_2003 == 84475
replace gdenr_2018 = 1054 if UCID_1934_2003 == 84475

replace gdename ="bern" if UCID_1934_2003 == 89341
replace gdenr = 351 if UCID_1934_2003 == 89341
replace gdenr_2012 = 351 if UCID_1934_2003 == 89341
replace gdenr_2018 = 351 if UCID_1934_2003 == 89341

replace gdename ="kuesnacht zh" if UCID_1934_2003 == 90040
replace gdenr = 154 if UCID_1934_2003 == 90040
replace gdenr_2012 = 154 if UCID_1934_2003 == 90040
replace gdenr_2018 = 154 if UCID_1934_2003 == 90040

replace gdename ="elgg" if UCID_1934_2003 == 90959
replace gdenr = 217 if UCID_1934_2003 == 90959
replace gdenr_2012 = 217 if UCID_1934_2003 == 90959
replace gdenr_2018 = 294 if UCID_1934_2003 == 90959

replace gdename ="horgen" if UCID_1934_2003 == 91180
replace gdenr = 133 if UCID_1934_2003 == 91180
replace gdenr_2012 = 133 if UCID_1934_2003 == 91180
replace gdenr_2018 = 295 if UCID_1934_2003 == 91180
				
replace gdename ="horgen" if UCID_1934_2003 == 91651
replace gdenr = 133 if UCID_1934_2003 == 91651
replace gdenr_2012 = 133 if UCID_1934_2003 == 91651
replace gdenr_2018 = 295 if UCID_1934_2003 == 91651
	
replace gdename ="pratteln" if UCID_1934_2003 == 92003
replace gdenr = 2831 if UCID_1934_2003 == 92003
replace gdenr_2012 = 2831 if UCID_1934_2003 == 92003
replace gdenr_2018 = 2831 if UCID_1934_2003 == 92003
				
replace gdename ="kyburg buchegg" if UCID_1934_2003 == 92261
replace gdenr = 2453 if UCID_1934_2003 == 92261
replace gdenr_2012 = 2453 if UCID_1934_2003 == 92261
replace gdenr_2018 = 2465 if UCID_1934_2003 == 92261
	
replace gdename ="baar" if UCID_1934_2003 == 92479
replace gdenr = 1701 if UCID_1934_2003 == 92479
replace gdenr_2012 = 1701 if UCID_1934_2003 == 92479
replace gdenr_2018 = 1701 if UCID_1934_2003 == 92479

replace gdename ="baar" if UCID_1934_2003 == 92527
replace gdenr = 1701 if UCID_1934_2003 == 92527
replace gdenr_2012 = 1701 if UCID_1934_2003 == 92527
replace gdenr_2018 = 1701 if UCID_1934_2003 == 92527

replace gdename ="st. gallen" if UCID_1934_2003 == 94550
replace gdenr = 3203 if UCID_1934_2003 == 94550
replace gdenr_2012 = 3203 if UCID_1934_2003 == 94550
replace gdenr_2018 = 3203 if UCID_1934_2003 == 94550
	
replace gdename ="lotzwil" if UCID_1934_2003 == 94805
replace gdenr = 331 if UCID_1934_2003 == 94805
replace gdenr_2012 = 331 if UCID_1934_2003 == 94805
replace gdenr_2018 = 331 if UCID_1934_2003 == 94805

replace gdename ="thunstetten" if UCID_1934_2003 == 95232
replace gdenr = 342 if UCID_1934_2003 == 95232
replace gdenr_2012 = 342 if UCID_1934_2003 == 95232
replace gdenr_2018 = 342 if UCID_1934_2003 == 95232
			
replace gdename ="zuerich" if UCID_1934_2003 == 95655
replace gdenr = 253 if UCID_1934_2003 == 95655 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 95655 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 95655
replace gdenr_2018 = 261 if UCID_1934_2003 == 95655

replace gdename ="bern" if UCID_1934_2003 == 98142
replace gdenr = 351 if UCID_1934_2003 == 98142
replace gdenr_2012 = 351 if UCID_1934_2003 == 98142
replace gdenr_2018 = 351 if UCID_1934_2003 == 98142

replace gdename ="ittigen" if UCID_1934_2003 == 98598
replace gdenr = 362 if UCID_1934_2003 == 98598
replace gdenr_2012 = 362 if UCID_1934_2003 == 98598
replace gdenr_2018 = 362 if UCID_1934_2003 ==  98598

replace gdename ="ittigen" if UCID_1934_2003 == 99076
replace gdenr = 362 if UCID_1934_2003 == 99076
replace gdenr_2012 = 362 if UCID_1934_2003 == 99076
replace gdenr_2018 = 362 if UCID_1934_2003 ==  99076
				
replace gdename ="zuerich" if UCID_1934_2003 == 99151
replace gdenr = 253 if UCID_1934_2003 == 99151
replace gdenr = 253 if UCID_1934_2003 == 99151 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 99151 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 99151
replace gdenr_2018 = 261 if UCID_1934_2003 == 99151

replace gdename ="schlieren" if UCID_1934_2003 == 99521
replace gdenr = 247 if UCID_1934_2003 == 99521
replace gdenr_2012 = 247 if UCID_1934_2003 == 99521
replace gdenr_2018 = 247 if UCID_1934_2003 == 99521

replace gdename ="locarno" if UCID_1934_2003 == 99703
replace gdenr = 5113 if UCID_1934_2003 == 99703
replace gdenr_2012 = 5113 if UCID_1934_2003 == 99703
replace gdenr_2018 = 5113 if UCID_1934_2003 == 99703

replace gdename ="bern" if UCID_1934_2003 == 100475
replace gdenr = 351 if UCID_1934_2003 == 100475
replace gdenr_2012 = 351 if UCID_1934_2003 == 100475
replace gdenr_2018 = 351 if UCID_1934_2003 == 100475

replace gdename ="bern" if UCID_1934_2003 == 100515
replace gdenr = 351 if UCID_1934_2003 == 100515
replace gdenr_2012 = 351 if UCID_1934_2003 == 100515
replace gdenr_2018 = 351 if UCID_1934_2003 == 100515

replace gdename ="bottighofen" if UCID_1934_2003 == 100627
replace gdenr = 4686 if UCID_1934_2003 == 100627 & year < 1994
replace gdenr = 4643 if UCID_1934_2003 == 100627 & year > 1993
replace gdenr_2012 = 4643 if UCID_1934_2003 == 100627
replace gdenr_2018 = 4643 if UCID_1934_2003 == 100627

replace gdename ="bern" if UCID_1934_2003 == 102006
replace gdenr = 351 if UCID_1934_2003 == 102006
replace gdenr_2012 = 351 if UCID_1934_2003 == 102006
replace gdenr_2018 = 351 if UCID_1934_2003 == 102006

replace gdename ="bern" if UCID_1934_2003 == 102297
replace gdenr = 351 if UCID_1934_2003 == 102297
replace gdenr_2012 = 351 if UCID_1934_2003 == 102297
replace gdenr_2018 = 351 if UCID_1934_2003 == 102297

replace gdename ="ittigen" if UCID_1934_2003 == 102503
replace gdenr = 362 if UCID_1934_2003 == 102503
replace gdenr_2012 = 362 if UCID_1934_2003 == 102503
replace gdenr_2018 = 362 if UCID_1934_2003 ==  102503
	
replace gdename ="murten" if UCID_1934_2003 == 102653
replace gdenr = 2275 if UCID_1934_2003 == 102653
replace gdenr_2012 = 2275 if UCID_1934_2003 == 102653
replace gdenr_2018 = 2275 if UCID_1934_2003 == 102653
	
replace gdename ="basel" if UCID_1934_2003 == 104405
replace gdenr = 2701 if UCID_1934_2003 == 104405
replace gdenr_2012 = 2701 if UCID_1934_2003 == 104405
replace gdenr_2018 = 2701 if UCID_1934_2003 == 104405

replace gdename ="biel (be)" if UCID_1934_2003 == 104780
replace gdenr = 371 if UCID_1934_2003 == 104780
replace gdenr_2012 = 371 if UCID_1934_2003 == 104780
replace gdenr_2018 = 371 if UCID_1934_2003 == 104780
	
replace gdename ="basel" if UCID_1934_2003 == 104864
replace gdenr = 2701 if UCID_1934_2003 == 104864
replace gdenr_2012 = 2701 if UCID_1934_2003 == 104864
replace gdenr_2018 = 2701 if UCID_1934_2003 == 104864

replace gdename ="ittigen" if UCID_1934_2003 == 105790
replace gdenr = 362 if UCID_1934_2003 == 105790
replace gdenr_2012 = 362 if UCID_1934_2003 == 105790
replace gdenr_2018 = 362 if UCID_1934_2003 == 105790 

replace gdename ="ittigen" if UCID_1934_2003 == 105798
replace gdenr = 362 if UCID_1934_2003 == 105798
replace gdenr_2012 = 362 if UCID_1934_2003 == 105798
replace gdenr_2018 = 362 if UCID_1934_2003 == 105798 

replace gdename ="ittigen" if UCID_1934_2003 == 105811
replace gdenr = 362 if UCID_1934_2003 == 105811
replace gdenr_2012 = 362 if UCID_1934_2003 == 105811
replace gdenr_2018 = 362 if UCID_1934_2003 == 105811
			
replace gdename ="bolligen" if UCID_1934_2003 == 105847
replace gdenr = 352 if UCID_1934_2003 == 105847
replace gdenr_2012 = 352 if UCID_1934_2003 == 105847
replace gdenr_2018 = 352 if UCID_1934_2003 == 105847

replace gdename ="ittigen" if UCID_1934_2003 == 105925
replace gdenr = 362 if UCID_1934_2003 == 105925
replace gdenr_2012 = 362 if UCID_1934_2003 == 105925
replace gdenr_2018 = 362 if UCID_1934_2003 == 105925

replace gdename ="ittigen" if UCID_1934_2003 == 105998
replace gdenr = 362 if UCID_1934_2003 == 105998
replace gdenr_2012 = 362 if UCID_1934_2003 == 105998
replace gdenr_2018 = 362 if UCID_1934_2003 == 105998

replace gdename ="koeniz" if UCID_1934_2003 == 106527
replace gdenr = 355 if UCID_1934_2003 == 106527
replace gdenr_2012 = 355 if UCID_1934_2003 == 106527
replace gdenr_2018 = 355 if UCID_1934_2003 == 106527

replace gdename ="ittigen" if UCID_1934_2003 == 109037
replace gdenr = 362 if UCID_1934_2003 == 109037
replace gdenr_2012 = 362 if UCID_1934_2003 == 109037
replace gdenr_2018 = 362 if UCID_1934_2003 == 109037 

replace gdename ="ittigen" if UCID_1934_2003 == 109046
replace gdenr = 362 if UCID_1934_2003 == 109046
replace gdenr_2012 = 362 if UCID_1934_2003 == 109046
replace gdenr_2018 = 362 if UCID_1934_2003 == 109046 

replace gdename ="ittigen" if UCID_1934_2003 == 109070
replace gdenr = 362 if UCID_1934_2003 == 109070
replace gdenr_2012 = 362 if UCID_1934_2003 == 109070
replace gdenr_2018 = 362 if UCID_1934_2003 == 109070

replace gdename ="ittigen" if UCID_1934_2003 == 109121
replace gdenr = 362 if UCID_1934_2003 == 109121
replace gdenr_2012 = 362 if UCID_1934_2003 == 109121
replace gdenr_2018 = 362 if UCID_1934_2003 == 109121

replace gdename ="ittigen" if UCID_1934_2003 == 109219
replace gdenr = 362 if UCID_1934_2003 == 109219
replace gdenr_2012 = 362 if UCID_1934_2003 == 109219
replace gdenr_2018 = 362 if UCID_1934_2003 ==  109219

replace gdename ="ittigen" if UCID_1934_2003 == 109261
replace gdenr = 362 if UCID_1934_2003 == 109261
replace gdenr_2012 = 362 if UCID_1934_2003 == 109261
replace gdenr_2018 = 362 if UCID_1934_2003 == 109261

replace gdename ="ittigen" if UCID_1934_2003 == 109275
replace gdenr = 362 if UCID_1934_2003 == 109275
replace gdenr_2012 = 362 if UCID_1934_2003 == 109275
replace gdenr_2018 = 362 if UCID_1934_2003 == 109275

replace gdename ="ittigen" if UCID_1934_2003 == 109323
replace gdenr = 362 if UCID_1934_2003 == 109323
replace gdenr_2012 = 362 if UCID_1934_2003 == 109323
replace gdenr_2018 = 362 if UCID_1934_2003 ==  109323

replace gdename ="ittigen" if UCID_1934_2003 == 109327
replace gdenr = 362 if UCID_1934_2003 == 109327
replace gdenr_2012 = 362 if UCID_1934_2003 == 109327
replace gdenr_2018 = 362 if UCID_1934_2003 == 109327

replace gdename ="ittigen" if UCID_1934_2003 == 109356
replace gdenr = 362 if UCID_1934_2003 == 109356
replace gdenr_2012 = 362 if UCID_1934_2003 == 109356
replace gdenr_2018 = 362 if UCID_1934_2003 == 109356

replace gdename ="ittigen" if UCID_1934_2003 == 109410
replace gdenr = 362 if UCID_1934_2003 == 109410
replace gdenr_2012 = 362 if UCID_1934_2003 == 109410
replace gdenr_2018 = 362 if UCID_1934_2003 == 109410

replace gdename ="ittigen" if UCID_1934_2003 == 109464
replace gdenr = 362 if UCID_1934_2003 == 109464
replace gdenr_2012 = 362 if UCID_1934_2003 == 109464
replace gdenr_2018 = 362 if UCID_1934_2003 == 109464

replace gdename ="ittigen" if UCID_1934_2003 == 109724
replace gdenr = 362 if UCID_1934_2003 == 109724
replace gdenr_2012 = 362 if UCID_1934_2003 == 109724
replace gdenr_2018 = 362 if UCID_1934_2003 == 109724

replace gdename ="ittigen" if UCID_1934_2003 == 109836
replace gdenr = 362 if UCID_1934_2003 == 109836
replace gdenr_2012 = 362 if UCID_1934_2003 == 109836
replace gdenr_2018 = 362 if UCID_1934_2003 == 109836

replace gdename ="ittigen" if UCID_1934_2003 == 109928
replace gdenr = 362 if UCID_1934_2003 == 109928
replace gdenr_2012 = 362 if UCID_1934_2003 == 109928
replace gdenr_2018 = 362 if UCID_1934_2003 == 109928

replace gdename ="biel (be)" if UCID_1934_2003 == 112254
replace gdenr = 371 if UCID_1934_2003 == 112254
replace gdenr_2012 = 371 if UCID_1934_2003 == 112254
replace gdenr_2018 = 371 if UCID_1934_2003 == 112254
				
replace gdename ="zuerich" if UCID_1934_2003 == 112980
replace gdenr = 253 if UCID_1934_2003 == 112980 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 112980 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 112980
replace gdenr_2018 = 261 if UCID_1934_2003 == 112980
	
replace gdename ="arth" if UCID_1934_2003 == 114138
replace gdenr = 1362 if UCID_1934_2003 == 114138
replace gdenr_2012 = 1362 if UCID_1934_2003 == 114138
replace gdenr_2018 = 1362 if UCID_1934_2003 == 114138
			
replace gdename ="waengi" if UCID_1934_2003 == 114786
replace gdenr = 4784 if UCID_1934_2003 == 114786 & year < 1969
replace gdenr = 4781 if UCID_1934_2003 == 114786 & year > 1968
replace gdenr_2012 = 4781 if UCID_1934_2003 == 114786
replace gdenr_2018 = 4781 if UCID_1934_2003 == 114786

replace gdename ="baeriswil" if UCID_1934_2003 == 114837
replace gdenr = 403 if UCID_1934_2003 == 114837
replace gdenr_2012 = 403 if UCID_1934_2003 == 114837
replace gdenr_2018 = 403 if UCID_1934_2003 == 114837

replace gdename ="herisau" if UCID_1934_2003 == 115838
replace gdenr = 3001 if UCID_1934_2003 == 115838
replace gdenr_2012 = 3001 if UCID_1934_2003 == 115838
replace gdenr_2018 = 3001 if UCID_1934_2003 == 115838
	
replace gdename ="gals" if UCID_1934_2003 == 117381
replace gdenr = 494 if UCID_1934_2003 == 117381
replace gdenr_2012 = 494 if UCID_1934_2003 == 117381
replace gdenr_2018 = 494 if UCID_1934_2003 == 117381
							
replace gdename ="scheunen" if UCID_1934_2003 == 117689
replace gdenr = 550 if UCID_1934_2003 == 117689
replace gdenr_2012 = 550 if UCID_1934_2003 == 117689
replace gdenr_2018 = 540 if UCID_1934_2003 == 117689
		
replace gdename ="jegenstorf" if UCID_1934_2003 == 117759
replace gdenr = 540 if UCID_1934_2003 == 117759
replace gdenr_2012 = 540 if UCID_1934_2003 == 117759
replace gdenr_2018 = 540 if UCID_1934_2003 == 117759
		
replace gdename ="bern" if UCID_1934_2003 == 117761
replace gdenr = 351 if UCID_1934_2003 == 117761
replace gdenr_2012 = 351 if UCID_1934_2003 == 117761
replace gdenr_2018 = 351 if UCID_1934_2003 == 117761
		
replace gdename ="urtenen" if UCID_1934_2003 == 118260
replace gdenr = 551 if UCID_1934_2003 == 118260
replace gdenr_2012 = 551 if UCID_1934_2003 == 118260
replace gdenr_2018 = 551 if UCID_1934_2003 == 118260
			
replace gdename ="nidau" if UCID_1934_2003 == 124403
replace gdenr = 743 if UCID_1934_2003 == 124403
replace gdenr_2012 = 743 if UCID_1934_2003 == 124403
replace gdenr_2018 = 743 if UCID_1934_2003 == 124403
				
replace gdename ="port" if UCID_1934_2003 == 124558
replace gdenr = 745 if UCID_1934_2003 == 124558
replace gdenr_2012 = 745 if UCID_1934_2003 == 124558
replace gdenr_2018 = 745 if UCID_1934_2003 == 124558
		
replace gdename ="oberwil im simmental" if UCID_1934_2003 == 125123
replace gdenr = 766 if UCID_1934_2003 == 125123
replace gdenr_2012 = 766 if UCID_1934_2003 == 125123
replace gdenr_2018 = 766 if UCID_1934_2003 ==125123
				
replace gdename ="leuk" if UCID_1934_2003 == 125803
replace gdenr = 6110 if UCID_1934_2003 == 125803
replace gdenr_2012 = 6110 if UCID_1934_2003 == 125803
replace gdenr_2018 = 6110 if UCID_1934_2003 == 125803

replace gdename ="saanen" if UCID_1934_2003 == 125989
replace gdenr = 843 if UCID_1934_2003 == 125989
replace gdenr_2012 = 843 if UCID_1934_2003 == 125989
replace gdenr_2018 = 843 if UCID_1934_2003 == 125989
				
replace gdename ="saanen" if UCID_1934_2003 == 126018
replace gdenr = 843 if UCID_1934_2003 == 126018
replace gdenr_2012 = 843 if UCID_1934_2003 == 126018
replace gdenr_2018 = 843 if UCID_1934_2003 == 126018

replace gdename ="saanen" if UCID_1934_2003 == 126110
replace gdenr = 843 if UCID_1934_2003 == 126110
replace gdenr_2012 = 843 if UCID_1934_2003 == 126110
replace gdenr_2018 = 843 if UCID_1934_2003 == 126110	
	
replace gdename ="saanen" if UCID_1934_2003 == 126465
replace gdenr = 843 if UCID_1934_2003 == 126465
replace gdenr_2012 = 843 if UCID_1934_2003 == 126465
replace gdenr_2018 = 843 if UCID_1934_2003 == 126465
				
replace gdename ="ruemlingen" if UCID_1934_2003 == 127299
replace gdenr = 2859 if UCID_1934_2003 == 127299
replace gdenr_2012 = 2859 if UCID_1934_2003 == 127299
replace gdenr_2018 = 2859 if UCID_1934_2003 == 127299
			
replace gdename ="toffen" if UCID_1934_2003 == 127353
replace gdenr = 884 if UCID_1934_2003 == 127353
replace gdenr_2012 = 884 if UCID_1934_2003 == 127353
replace gdenr_2018 = 884 if UCID_1934_2003 == 127353
		
replace gdename ="biel (be)" if UCID_1934_2003 == 127777
replace gdenr = 371 if UCID_1934_2003 == 127777
replace gdenr_2012 = 371 if UCID_1934_2003 == 127777
replace gdenr_2018 = 371 if UCID_1934_2003 == 127777
	
replace gdename ="affoltern im emmental" if UCID_1934_2003 == 130882
replace gdenr = 951 if UCID_1934_2003 == 130882
replace gdenr_2012 = 951 if UCID_1934_2003 == 130882
replace gdenr_2018 = 951 if UCID_1934_2003 == 130882
				
replace gdename ="weesen" if UCID_1934_2003 == 131334
replace gdenr = 3316 if UCID_1934_2003 == 131334
replace gdenr_2012 = 3316 if UCID_1934_2003 == 131334
replace gdenr_2018 = 3316 if UCID_1934_2003 == 131334
		
replace gdename ="thoerigen" if UCID_1934_2003 == 132176
replace gdenr = 989 if UCID_1934_2003 == 132176
replace gdenr_2012 = 989 if UCID_1934_2003 == 132176
replace gdenr_2018 = 989 if UCID_1934_2003 == 132176
			
replace gdename ="oensingen" if UCID_1934_2003 == 132299
replace gdenr = 2407 if UCID_1934_2003 == 132299
replace gdenr_2012 = 2407 if UCID_1934_2003 == 132299
replace gdenr_2018 = 2407 if UCID_1934_2003 == 132299
			
replace gdename ="sachseln" if UCID_1934_2003 == 132533
replace gdenr = 1406 if UCID_1934_2003 == 132533
replace gdenr_2012 = 1406 if UCID_1934_2003 == 132533
replace gdenr_2018 = 1406 if UCID_1934_2003 == 132533
			
replace gdename ="hasle lu" if UCID_1934_2003 == 132544
replace gdenr = 1005 if UCID_1934_2003 == 132544
replace gdenr_2012 = 1005 if UCID_1934_2003 == 132544
replace gdenr_2018 = 1005 if UCID_1934_2003 == 132544
			
replace gdename ="emmen" if UCID_1934_2003 == 133079
replace gdenr = 1024 if UCID_1934_2003 == 133079
replace gdenr_2012 = 1024 if UCID_1934_2003 == 133079
replace gdenr_2018 = 1024 if UCID_1934_2003 == 133079

replace gdename ="emmen" if UCID_1934_2003 == 133174
replace gdenr = 1024 if UCID_1934_2003 == 133174
replace gdenr_2012 = 1024 if UCID_1934_2003 == 133174
replace gdenr_2018 = 1024 if UCID_1934_2003 == 133174
	
replace gdename ="emmen" if UCID_1934_2003 == 133339
replace gdenr = 1024 if UCID_1934_2003 == 133339
replace gdenr_2012 = 1024 if UCID_1934_2003 == 133339
replace gdenr_2018 = 1024 if UCID_1934_2003 == 133339
	
replace gdename ="ernen" if UCID_1934_2003 == 133321
replace gdenr = 6056 if UCID_1934_2003 == 133321
replace gdenr_2012 = 6056 if UCID_1934_2003 == 133321
replace gdenr_2018 = 6056 if UCID_1934_2003 == 133321
						
replace gdename ="baden" if UCID_1934_2003 == 133987
replace gdenr = 4021 if UCID_1934_2003 == 133987
replace gdenr_2012 = 4021 if UCID_1934_2003 == 133987
replace gdenr_2018 = 4021 if UCID_1934_2003 == 133987

replace gdename ="buchrain" if UCID_1934_2003 == 134908
replace gdenr = 1052 if UCID_1934_2003 == 134908
replace gdenr_2012 = 1052 if UCID_1934_2003 == 134908
replace gdenr_2018 = 1052 if UCID_1934_2003 == 134908
				
replace gdename ="horw" if UCID_1934_2003 == 135507
replace gdenr = 1058 if UCID_1934_2003 == 135507
replace gdenr_2012 = 1058 if UCID_1934_2003 == 135507
replace gdenr_2018 = 1058 if UCID_1934_2003 == 135507

replace gdename ="flawil" if UCID_1934_2003 == 135884
replace gdenr = 3402 if UCID_1934_2003 == 135884
replace gdenr_2012 = 3402 if UCID_1934_2003 == 135884
replace gdenr_2018 = 3402 if UCID_1934_2003 == 135884
		
replace gdename ="horw" if UCID_1934_2003 == 135899
replace gdenr = 1058 if UCID_1934_2003 == 135899
replace gdenr_2012 = 1058 if UCID_1934_2003 == 135899
replace gdenr_2018 = 1058 if UCID_1934_2003 == 135899
			
replace gdename ="montana" if UCID_1934_2003 == 136426
replace gdenr = 6243 if UCID_1934_2003 == 136426
replace gdenr_2012 = 6243 if UCID_1934_2003 == 136426
replace gdenr_2018 = 6243 if UCID_1934_2003 == 136426
			
replace gdename ="kriens" if UCID_1934_2003 == 136809
replace gdenr = 1059 if UCID_1934_2003 == 136809
replace gdenr_2012 = 1059 if UCID_1934_2003 == 136809
replace gdenr_2018 = 1059 if UCID_1934_2003 == 136809
			
replace gdename ="stansstad" if UCID_1934_2003 == 137098
replace gdenr = 1510 if UCID_1934_2003 == 137098
replace gdenr_2012 = 1510 if UCID_1934_2003 == 137098
replace gdenr_2018 = 1510 if UCID_1934_2003 == 137098
		
replace gdenr = 253 if UCID_1934_2003 == 137771 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 137771 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 137771
replace gdenr_2018 = 261 if UCID_1934_2003 == 137771
				
replace gdename ="littau" if UCID_1934_2003 == 138400
replace gdenr = 1060 if UCID_1934_2003 == 138400
replace gdenr_2012 = 1061 if UCID_1934_2003 == 138400
replace gdenr_2018 = 1061 if UCID_1934_2003 == 138400
		
replace gdename ="littau" if UCID_1934_2003 == 139074
replace gdenr =1060  if UCID_1934_2003 == 139074
replace gdenr_2012 = 1061 if UCID_1934_2003 == 139074
replace gdenr_2018 = 1061 if UCID_1934_2003 == 139074
				
replace gdename ="littau" if UCID_1934_2003 == 144587
replace gdenr = 1060 if UCID_1934_2003 == 144587
replace gdenr_2012 = 1061 if UCID_1934_2003 == 144587
replace gdenr_2018 = 1061 if UCID_1934_2003 == 144587
	
replace gdename ="malters" if UCID_1934_2003 == 145124
replace gdenr = 1062 if UCID_1934_2003 == 145124
replace gdenr_2012 = 1062 if UCID_1934_2003 == 145124
replace gdenr_2018 = 1062 if UCID_1934_2003 == 145124
		
replace gdename ="risch" if UCID_1934_2003 == 145664
replace gdenr = 1707 if UCID_1934_2003 == 145664
replace gdenr_2012 = 1707 if UCID_1934_2003 == 145664
replace gdenr_2018 = 1707 if UCID_1934_2003 == 145664

replace gdename ="triengen" if UCID_1934_2003 == 148333
replace gdenr = 1104 if UCID_1934_2003 == 148333
replace gdenr_2012 = 1104 if UCID_1934_2003 == 148333
replace gdenr_2018 = 1104 if UCID_1934_2003 == 148333
					
replace gdename ="triengen" if UCID_1934_2003 == 148346
replace gdenr = 1104 if UCID_1934_2003 == 148346
replace gdenr_2012 = 1104 if UCID_1934_2003 == 148346
replace gdenr_2018 = 1104 if UCID_1934_2003 == 148346
			
replace gdename ="triengen" if UCID_1934_2003 == 148349
replace gdenr = 1104 if UCID_1934_2003 == 148349
replace gdenr_2012 = 1104 if UCID_1934_2003 == 148349
replace gdenr_2018 = 1104 if UCID_1934_2003 == 148349
			
replace gdename ="luthern" if UCID_1934_2003 == 148846
replace gdenr = 1135 if UCID_1934_2003 == 148846
replace gdenr_2012 = 1135 if UCID_1934_2003 == 148846
replace gdenr_2018 = 1135 if UCID_1934_2003 ==148846
		
replace gdename ="bern" if UCID_1934_2003 == 150553
replace gdenr = 351 if UCID_1934_2003 == 150553
replace gdenr_2012 = 351 if UCID_1934_2003 == 150553
replace gdenr_2018 = 351 if UCID_1934_2003 ==150553
	
replace gdename ="freienbach" if UCID_1934_2003 == 152037
replace gdenr = 1322 if UCID_1934_2003 == 152037
replace gdenr_2012 = 1322 if UCID_1934_2003 == 152037
replace gdenr_2018 = 1322 if UCID_1934_2003 == 152037
						
replace gdename ="freienbach" if UCID_1934_2003 == 152450
replace gdenr = 1322 if UCID_1934_2003 == 152450
replace gdenr_2012 = 1322 if UCID_1934_2003 == 152450
replace gdenr_2018 = 1322 if UCID_1934_2003 == 152450

replace gdename ="schuebelbach" if UCID_1934_2003 == 156482
replace gdenr = 1346 if UCID_1934_2003 == 156482
replace gdenr_2012 = 1346 if UCID_1934_2003 == 156482
replace gdenr_2018 = 1346 if UCID_1934_2003 ==156482

replace gdename ="tuggen" if UCID_1934_2003 == 156812
replace gdenr = 1347 if UCID_1934_2003 == 156812
replace gdenr_2012 = 1347 if UCID_1934_2003 == 156812
replace gdenr_2018 = 1347 if UCID_1934_2003 ==156812
			
replace gdename ="schwyz" if UCID_1934_2003 == 158552
replace gdenr = 1372 if UCID_1934_2003 == 158552
replace gdenr_2012 = 1372 if UCID_1934_2003 == 158552
replace gdenr_2018 = 1372 if UCID_1934_2003 ==158552
			
replace gdename ="fluehli" if UCID_1934_2003 == 160084
replace gdenr = 1004 if UCID_1934_2003 == 160084
replace gdenr_2012 = 1004 if UCID_1934_2003 == 160084
replace gdenr_2018 = 1004 if UCID_1934_2003 ==160084

replace gdename ="saanen" if UCID_1934_2003 == 160919
replace gdenr = 843 if UCID_1934_2003 == 160919
replace gdenr_2012 = 843 if UCID_1934_2003 == 160919
replace gdenr_2018 = 843 if UCID_1934_2003 ==160919
		
replace gdename ="bruegg" if UCID_1934_2003 == 161892
replace gdenr = 733 if UCID_1934_2003 == 161892
replace gdenr_2012 = 733 if UCID_1934_2003 == 161892
replace gdenr_2018 = 733 if UCID_1934_2003 ==161892
				
replace gdename ="stansstad" if UCID_1934_2003 == 166230
replace gdenr = 1510 if UCID_1934_2003 == 166230
replace gdenr_2012 = 1510 if UCID_1934_2003 == 166230
replace gdenr_2018 = 1510 if UCID_1934_2003 ==166230

replace gdename ="naefels" if UCID_1934_2003 == 166752
replace gdenr = 1619 if UCID_1934_2003 == 166752
replace gdenr_2012 = 1619 if UCID_1934_2003 ==166752
replace gdenr_2018 = 1630 if UCID_1934_2003 ==166752
		
replace gdename ="naefels" if UCID_1934_2003 == 166757
replace gdenr = 1619 if UCID_1934_2003 == 166757
replace gdenr_2012 = 1619 if UCID_1934_2003 == 166757
replace gdenr_2018 = 1630 if UCID_1934_2003 ==166757

replace gdename ="riedern" if UCID_1934_2003 == 168526
replace gdenr = 1625 if UCID_1934_2003 == 168526
replace gdenr_2012 = 1625 if UCID_1934_2003 == 168526
replace gdenr_2018 = 1632 if UCID_1934_2003 ==168526
		
replace gdename ="riedern" if UCID_1934_2003 == 170724
replace gdenr = 1625 if UCID_1934_2003 == 170724
replace gdenr_2012 = 1625 if UCID_1934_2003 == 170724
replace gdenr_2018 = 1632 if UCID_1934_2003 ==170724

replace gdename ="zuerich" if UCID_1934_2003 == 204818
replace gdenr = 253 if UCID_1934_2003 == 204818 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 204818 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 204818
replace gdenr_2018 = 261 if UCID_1934_2003 == 204818

replace gdename ="zuerich" if UCID_1934_2003 == 171491
replace gdenr = 253 if UCID_1934_2003 == 171491 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 171491 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 171491
replace gdenr_2018 = 261 if UCID_1934_2003 == 171491

replace gdename ="zuerich" if UCID_1934_2003 == 213554			
replace gdenr = 253 if UCID_1934_2003 == 213554 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 213554 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 213554
replace gdenr_2018 = 261 if UCID_1934_2003 == 213554
		
replace gdename ="riedern" if UCID_1934_2003 == 172277
replace gdenr = 1625 if UCID_1934_2003 == 172277
replace gdenr_2012 = 1625 if UCID_1934_2003 == 172277
replace gdenr_2018 = 1632 if UCID_1934_2003 ==172277
	
replace gdename ="zug" if UCID_1934_2003 == 172860
replace gdenr = 1711 if UCID_1934_2003 == 172860
replace gdenr_2012 = 1711 if UCID_1934_2003 == 172860
replace gdenr_2018 = 1711 if UCID_1934_2003 == 172860

replace gdename ="zug" if UCID_1934_2003 == 187028
replace gdenr = 1711 if UCID_1934_2003 == 187028
replace gdenr_2012 = 1711 if UCID_1934_2003 == 187028
replace gdenr_2018 = 1711 if UCID_1934_2003 ==187028
	
replace gdename ="zug" if UCID_1934_2003 == 207337
replace gdenr = 1711 if UCID_1934_2003 == 207337
replace gdenr_2012 = 1711 if UCID_1934_2003 == 207337
replace gdenr_2018 = 1711 if UCID_1934_2003 ==207337
				
replace gdename ="steinhausen" if UCID_1934_2003 == 197794
replace gdenr = 1708 if UCID_1934_2003 == 197794
replace gdenr_2012 = 1708 if UCID_1934_2003 == 197794
replace gdenr_2018 = 1708 if UCID_1934_2003 ==197794
			
replace gdename ="oberglatt" if UCID_1934_2003 == 208768
replace gdenr = 92 if UCID_1934_2003 == 208768
replace gdenr_2012 = 92 if UCID_1934_2003 == 208768
replace gdenr_2018 =92 if UCID_1934_2003 ==208768
						
replace gdename ="sion" if UCID_1934_2003 == 211887
replace gdenr = 6266 if UCID_1934_2003 == 211887
replace gdenr_2012 =6266  if UCID_1934_2003 == 211887
replace gdenr_2018 = 6266 if UCID_1934_2003 ==211887
						
replace gdename ="wilen bei wil" if UCID_1934_2003 == 212701
replace gdenr = 4752 if UCID_1934_2003 == 212701 & year < 1998
replace gdenr = 4786 if UCID_1934_2003 == 212701 & year > 1997
replace gdenr_2012 = 4786 if UCID_1934_2003 == 212701
replace gdenr_2018 = 4786 if UCID_1934_2003 ==212701
					
replace gdename ="manno" if UCID_1934_2003 == 216448
replace gdenr = 5194 if UCID_1934_2003 == 216448
replace gdenr_2012 = 5194 if UCID_1934_2003 == 216448
replace gdenr_2018 = 5194 if UCID_1934_2003 ==216448
			
replace gdename ="montreux" if UCID_1934_2003 == 217205
replace gdenr = 5886 if UCID_1934_2003 == 217205
replace gdenr_2012 = 5886 if UCID_1934_2003 == 217205
replace gdenr_2018 = 5886 if UCID_1934_2003 ==217205
		
replace gdename ="bulle" if UCID_1934_2003 == 218467
replace gdenr = 2125 if UCID_1934_2003 == 218467
replace gdenr_2012 = 2125 if UCID_1934_2003 == 218467
replace gdenr_2018 = 2125 if UCID_1934_2003 ==218467
				
replace gdename ="bulle" if UCID_1934_2003 == 218523
replace gdenr = 2125 if UCID_1934_2003 == 218523
replace gdenr_2012 = 2125 if UCID_1934_2003 == 218523
replace gdenr_2018 = 2125 if UCID_1934_2003 ==218523

replace gdename ="geroldswil" if UCID_1934_2003 == 218952
replace gdenr = 244 if UCID_1934_2003 == 218952
replace gdenr_2012 = 244 if UCID_1934_2003 == 218952
replace gdenr_2018 = 244 if UCID_1934_2003 ==218952
				
replace gdename ="charmey" if UCID_1934_2003 == 219381
replace gdenr = 2127 if UCID_1934_2003 == 219381
replace gdenr_2012 = 2127 if UCID_1934_2003 == 219381
replace gdenr_2018 = 2163 if UCID_1934_2003 == 219381
	
replace gdename ="fribourg" if UCID_1934_2003 == 221667
replace gdenr = 2196 if UCID_1934_2003 == 221667
replace gdenr_2012 = 2196 if UCID_1934_2003 == 221667
replace gdenr_2018 = 2196 if UCID_1934_2003 ==221667

replace gdename ="fribourg" if UCID_1934_2003 == 222008
replace gdenr = 2196 if UCID_1934_2003 == 222008
replace gdenr_2012 = 2196 if UCID_1934_2003 == 222008
replace gdenr_2018 = 2196 if UCID_1934_2003 ==222008

replace gdename ="fribourg" if UCID_1934_2003 == 232995
replace gdenr = 2196 if UCID_1934_2003 == 232995
replace gdenr_2012 = 2196 if UCID_1934_2003 == 232995
replace gdenr_2018 = 2196 if UCID_1934_2003 ==232995
	
replace gdename ="baar" if UCID_1934_2003 == 224926
replace gdenr = 1701 if UCID_1934_2003 == 224926
replace gdenr_2012 = 1701 if UCID_1934_2003 == 224926
replace gdenr_2018 = 1701 if UCID_1934_2003 ==224926

replace gdename ="bulle" if UCID_1934_2003 == 227820
replace gdenr = 2125 if UCID_1934_2003 == 227820
replace gdenr_2012 = 2125 if UCID_1934_2003 == 227820
replace gdenr_2018 = 2125 if UCID_1934_2003 ==227820

replace gdename ="givisiez" if UCID_1934_2003 == 227992
replace gdenr = 2197 if UCID_1934_2003 == 227992
replace gdenr_2012 = 2197 if UCID_1934_2003 == 227992
replace gdenr_2018 = 2197 if UCID_1934_2003 ==227992
	
replace gdename ="cressier fr" if UCID_1934_2003 == 230086
replace gdenr = 2257 if UCID_1934_2003 == 230086
replace gdenr_2012 = 2257 if UCID_1934_2003 == 230086
replace gdenr_2018 = 2257 if UCID_1934_2003 ==230086
	
replace gdename ="ollon" if UCID_1934_2003 == 233284
replace gdenr = 5409 if UCID_1934_2003 == 233284
replace gdenr_2012 = 5409 if UCID_1934_2003 == 233284
replace gdenr_2018 = 5409 if UCID_1934_2003 ==233284

replace gdename ="lalden" if UCID_1934_2003 == 233675
replace gdenr = 6286 if UCID_1934_2003 == 233675
replace gdenr_2012 = 6286 if UCID_1934_2003 == 233675
replace gdenr_2018 = 6286 if UCID_1934_2003 ==233675
	
replace gdename ="st antoni" if UCID_1934_2003 == 236106
replace gdenr = 2302 if UCID_1934_2003 == 236106
replace gdenr_2012 = 2302 if UCID_1934_2003 == 236106
replace gdenr_2018 = 2302 if UCID_1934_2003 ==236106
		
replace gdename ="oftringen" if UCID_1934_2003 == 237541
replace gdenr = 4280 if UCID_1934_2003 == 237541
replace gdenr_2012 = 4280 if UCID_1934_2003 == 237541
replace gdenr_2018 = 4280 if UCID_1934_2003 ==237541

replace gdename ="oensingen" if UCID_1934_2003 == 237755
replace gdenr = 2407 if UCID_1934_2003 == 237755
replace gdenr_2012 = 2407 if UCID_1934_2003 == 237755
replace gdenr_2018 = 2407 if UCID_1934_2003 ==237755

replace gdename ="oensingen" if UCID_1934_2003 == 237769
replace gdenr = 2407 if UCID_1934_2003 == 237769
replace gdenr_2012 = 2407 if UCID_1934_2003 == 237769
replace gdenr_2018 = 2407 if UCID_1934_2003 ==237769

replace gdename ="oensingen" if UCID_1934_2003 == 237785
replace gdenr = 2407 if UCID_1934_2003 == 237785
replace gdenr_2012 = 2407 if UCID_1934_2003 == 237785
replace gdenr_2018 = 2407 if UCID_1934_2003 ==237785

replace gdename ="oensingen" if UCID_1934_2003 == 237791
replace gdenr = 2407 if UCID_1934_2003 == 237791
replace gdenr_2012 = 2407 if UCID_1934_2003 == 237791
replace gdenr_2018 = 2407 if UCID_1934_2003 ==237791
	
replace gdename ="oensingen" if UCID_1934_2003 == 237798
replace gdenr = 2407 if UCID_1934_2003 == 237798
replace gdenr_2012 = 2407 if UCID_1934_2003 == 237798
replace gdenr_2018 = 2407 if UCID_1934_2003 ==237798
	
replace gdename ="breitenbach" if UCID_1934_2003 == 238016
replace gdenr = 2613 if UCID_1934_2003 == 238016
replace gdenr_2012 = 2613 if UCID_1934_2003 == 238016
replace gdenr_2018 = 2613 if UCID_1934_2003 ==238016
	
replace gdename ="hauenstein ifenthal" if UCID_1934_2003 == 238927
replace gdenr = 2491 if UCID_1934_2003 == 238927
replace gdenr_2012 = 2491 if UCID_1934_2003 == 238927
replace gdenr_2018 = 2491 if UCID_1934_2003 ==238927
	
replace gdename ="deitingen" if UCID_1934_2003 == 239712
replace gdenr = 2516 if UCID_1934_2003 == 239712
replace gdenr_2012 = 2516 if UCID_1934_2003 == 239712
replace gdenr_2018 = 2516 if UCID_1934_2003 ==239712
	
replace gdename ="bellach" if UCID_1934_2003 == 240871
replace gdenr = 2542 if UCID_1934_2003 == 240871
replace gdenr_2012 = 2542 if UCID_1934_2003 == 240871
replace gdenr_2018 = 2542 if UCID_1934_2003 ==240871
	
replace gdename ="bettlach" if UCID_1934_2003 == 240888
replace gdenr = 2543 if UCID_1934_2003 == 240888
replace gdenr_2012 = 2543 if UCID_1934_2003 == 240888
replace gdenr_2018 = 2543 if UCID_1934_2003 ==240888
	
replace gdename ="gunzgen" if UCID_1934_2003 == 241344
replace gdenr = 2578 if UCID_1934_2003 == 241344
replace gdenr_2012 = 2578 if UCID_1934_2003 == 241344
replace gdenr_2018 = 2578 if UCID_1934_2003 ==241344
	
replace gdename ="arbon" if UCID_1934_2003 == 241377
replace gdenr = 4401 if UCID_1934_2003 == 241377
replace gdenr_2012 = 4401 if UCID_1934_2003 == 241377
replace gdenr_2018 = 4401 if UCID_1934_2003 ==241377
	
replace gdename ="graenichen" if UCID_1934_2003 == 241411
replace gdenr = 4006 if UCID_1934_2003 == 241411
replace gdenr_2012 = 4006 if UCID_1934_2003 == 241411
replace gdenr_2018 = 4006 if UCID_1934_2003 ==241411

replace gdename ="daeniken" if UCID_1934_2003 == 242439
replace gdenr = 2572 if UCID_1934_2003 == 242439
replace gdenr_2012 = 2572 if UCID_1934_2003 == 242439
replace gdenr_2018 = 2572 if UCID_1934_2003 ==242439
	
replace gdename ="solothurn" if UCID_1934_2003 == 243311
replace gdenr = 2601 if UCID_1934_2003 == 243311
replace gdenr_2012 = 2601 if UCID_1934_2003 == 243311
replace gdenr_2018 = 2601 if UCID_1934_2003 ==243311
	
replace gdename ="beinwil am see" if UCID_1934_2003 == 246647
replace gdenr = 4131 if UCID_1934_2003 == 246647
replace gdenr_2012 =4131 if UCID_1934_2003 == 246647
replace gdenr_2018 =4131 if UCID_1934_2003 ==246647

replace gdename ="zuerich" if UCID_1934_2003 == 248947			
replace gdenr = 253 if UCID_1934_2003 == 248947 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 248947 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 248947
replace gdenr_2018 = 261 if UCID_1934_2003 == 248947

replace gdename ="zuerich" if UCID_1934_2003 == 290924			
replace gdenr = 253 if UCID_1934_2003 == 290924 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 290924 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 290924
replace gdenr_2018 = 261 if UCID_1934_2003 == 290924
				
replace gdename ="binningen" if UCID_1934_2003 == 253495
replace gdenr = 2806 if UCID_1934_2003 == 253495 & year < 1994
replace gdenr = 2765 if UCID_1934_2003 == 253495 & year > 1993
replace gdenr_2012 = 2765 if UCID_1934_2003 == 253495
replace gdenr_2018 = 2765 if UCID_1934_2003 ==253495
	
replace gdename ="basel" if UCID_1934_2003 == 253919
replace gdenr = 2701 if UCID_1934_2003 == 253919
replace gdenr_2012 = 2701 if UCID_1934_2003 == 253919
replace gdenr_2018 = 2701 if UCID_1934_2003 ==253919

replace gdename ="basel" if UCID_1934_2003 == 255294
replace gdenr = 2701 if UCID_1934_2003 == 255294
replace gdenr_2012 = 2701 if UCID_1934_2003 == 255294
replace gdenr_2018 = 2701 if UCID_1934_2003 ==255294

replace gdename ="basel" if UCID_1934_2003 == 255364
replace gdenr = 2701 if UCID_1934_2003 == 255364
replace gdenr_2012 = 2701 if UCID_1934_2003 == 255364
replace gdenr_2018 = 2701 if UCID_1934_2003 ==255364

replace gdename ="basel" if UCID_1934_2003 == 260917
replace gdenr = 2701 if UCID_1934_2003 == 260917
replace gdenr_2012 = 2701 if UCID_1934_2003 == 260917
replace gdenr_2018 = 2701 if UCID_1934_2003 ==260917

replace gdename ="basel" if UCID_1934_2003 == 261599
replace gdenr = 2701 if UCID_1934_2003 == 261599
replace gdenr_2012 = 2701 if UCID_1934_2003 == 261599
replace gdenr_2018 = 2701 if UCID_1934_2003 ==261599

replace gdename ="basel" if UCID_1934_2003 == 257248
replace gdenr = 2701 if UCID_1934_2003 == 257248
replace gdenr_2012 = 2701 if UCID_1934_2003 == 257248
replace gdenr_2018 = 2701 if UCID_1934_2003 ==257248
	
replace gdename ="mels" if UCID_1934_2003 == 256691
replace gdenr = 3293 if UCID_1934_2003 == 256691
replace gdenr_2012 = 3293 if UCID_1934_2003 == 256691
replace gdenr_2018 = 3293 if UCID_1934_2003 ==256691

replace gdename ="ittigen" if UCID_1934_2003 == 258384
replace gdenr = 362 if UCID_1934_2003 == 258384
replace gdenr_2012 = 362 if UCID_1934_2003 == 258384
replace gdenr_2018 = 362 if UCID_1934_2003 ==258384

replace gdename ="ittigen" if UCID_1934_2003 == 277995
replace gdenr = 362 if UCID_1934_2003 == 277995
replace gdenr_2012 = 362 if UCID_1934_2003 == 277995
replace gdenr_2018 = 362 if UCID_1934_2003 ==277995
	
replace gdename ="weggis" if UCID_1934_2003 == 263400
replace gdenr = 1069 if UCID_1934_2003 == 263400
replace gdenr_2012 = 1069 if UCID_1934_2003 == 263400
replace gdenr_2018 = 1069 if UCID_1934_2003 ==263400
	
replace gdename ="baettwil" if UCID_1934_2003 == 270696
replace gdenr = 2471 if UCID_1934_2003 == 270696
replace gdenr_2012 = 2471 if UCID_1934_2003 == 270696
replace gdenr_2018 = 2471 if UCID_1934_2003 ==270696
	
replace gdename ="metzerlen" if UCID_1934_2003 == 271120
replace gdenr = 2477 if UCID_1934_2003 == 271120
replace gdenr_2012 = 2477 if UCID_1934_2003 == 271120
replace gdenr_2018 = 2477 if UCID_1934_2003 ==271120

replace gdename ="muttenz" if UCID_1934_2003 == 271812
replace gdenr = 2811 if UCID_1934_2003 == 271812 & year < 1994
replace gdenr = 2770 if UCID_1934_2003 == 271812 & year > 1993
replace gdenr_2012 = 2770 if UCID_1934_2003 == 271812
replace gdenr_2018 = 2770 if UCID_1934_2003 ==271812

replace gdename ="wohlen ag" if UCID_1934_2003 == 274748
replace gdenr = 4082 if UCID_1934_2003 == 274748
replace gdenr_2012 = 4082 if UCID_1934_2003 == 274748
replace gdenr_2018 = 4082 if UCID_1934_2003 ==274748

replace gdename ="bubendorf" if UCID_1934_2003 == 275097
replace gdenr = 2823 if UCID_1934_2003 == 275097
replace gdenr_2012 = 2823 if UCID_1934_2003 == 275097
replace gdenr_2018 = 2823 if UCID_1934_2003 ==275097
	
replace gdename ="liestal" if UCID_1934_2003 == 276658
replace gdenr = 2829 if UCID_1934_2003 == 276658
replace gdenr_2012 = 2829 if UCID_1934_2003 == 276658
replace gdenr_2018 = 2829 if UCID_1934_2003 ==276658

replace gdename ="boeckten" if UCID_1934_2003 == 277653
replace gdenr = 2842 if UCID_1934_2003 == 277653
replace gdenr_2012 = 2842 if UCID_1934_2003 == 277653
replace gdenr_2018 = 2842 if UCID_1934_2003 ==277653

replace gdename ="boeckten" if UCID_1934_2003 == 277659
replace gdenr = 2842 if UCID_1934_2003 == 277659
replace gdenr_2012 = 2842 if UCID_1934_2003 == 277659
replace gdenr_2018 = 2842 if UCID_1934_2003 ==277659
	
replace gdename ="boeckten" if UCID_1934_2003 == 277673
replace gdenr = 2842 if UCID_1934_2003 == 277673
replace gdenr_2012 = 2842 if UCID_1934_2003 == 277673
replace gdenr_2018 = 2842 if UCID_1934_2003 ==277673
	
replace gdename ="diegten" if UCID_1934_2003 == 278731
replace gdenr = 2884 if UCID_1934_2003 == 278731
replace gdenr_2012 = 2884 if UCID_1934_2003 == 278731
replace gdenr_2018 = 2884 if UCID_1934_2003 ==278731

replace gdename ="balsthal" if UCID_1934_2003 == 279073
replace gdenr = 2422 if UCID_1934_2003 == 279073
replace gdenr_2012 = 2422 if UCID_1934_2003 == 279073
replace gdenr_2018 = 2422 if UCID_1934_2003 ==279073

replace gdename ="herisau" if UCID_1934_2003 == 283042
replace gdenr = 3001 if UCID_1934_2003 == 283042
replace gdenr_2012 = 3001 if UCID_1934_2003 == 283042
replace gdenr_2018 = 3001 if UCID_1934_2003 ==283042
	
replace gdename ="gais" if UCID_1934_2003 == 284517
replace gdenr = 3022 if UCID_1934_2003 == 284517
replace gdenr_2012 = 3022 if UCID_1934_2003 == 284517
replace gdenr_2018 = 3022 if UCID_1934_2003 ==284517

replace gdename ="walzenhausen" if UCID_1934_2003 == 286030
replace gdenr = 3037 if UCID_1934_2003 == 286030
replace gdenr_2012 = 3037 if UCID_1934_2003 == 286030
replace gdenr_2018 = 3037 if UCID_1934_2003 ==286030
				
replace gdename ="sirnach" if UCID_1934_2003 == 289027
replace gdenr = 4764 if UCID_1934_2003 == 289027 & year < 1997
replace gdenr = 4761 if UCID_1934_2003 == 289027 & year > 1996
replace gdenr_2012 = 4761 if UCID_1934_2003 == 289027
replace gdenr_2018 = 4761 if UCID_1934_2003 ==289027
	
replace gdename ="st. gallen" if UCID_1934_2003 == 290669
replace gdenr = 3203 if UCID_1934_2003 == 290669
replace gdenr_2012 = 3203 if UCID_1934_2003 == 290669
replace gdenr_2018 = 3203 if UCID_1934_2003 ==290669

replace gdename ="wittenbach" if UCID_1934_2003 == 294387
replace gdenr = 3204 if UCID_1934_2003 == 294387
replace gdenr_2012 = 3204 if UCID_1934_2003 == 294387
replace gdenr_2018 = 3204 if UCID_1934_2003 ==294387

replace gdename ="goldach" if UCID_1934_2003 == 294654
replace gdenr = 3213 if UCID_1934_2003 == 294654
replace gdenr_2012 = 3213 if UCID_1934_2003 == 294654
replace gdenr_2018 = 3213 if UCID_1934_2003 ==294654
	
replace gdename ="goldach" if UCID_1934_2003 == 294901
replace gdenr = 3213 if UCID_1934_2003 == 294901
replace gdenr_2012 = 3213 if UCID_1934_2003 == 294901
replace gdenr_2018 = 3213 if UCID_1934_2003 ==294901
	
replace gdename ="rorschacherberg" if UCID_1934_2003 == 295677
replace gdenr = 3216 if UCID_1934_2003 == 295677
replace gdenr_2012 = 3216 if UCID_1934_2003 == 295677
replace gdenr_2018 = 3216 if UCID_1934_2003 ==295677
	
replace gdename ="balgach" if UCID_1934_2003 == 296996
replace gdenr = 3232 if UCID_1934_2003 == 296996
replace gdenr_2012 = 3232 if UCID_1934_2003 == 296996
replace gdenr_2018 = 3232 if UCID_1934_2003 ==296996
	
replace gdename ="st. margrethen" if UCID_1934_2003 == 298328
replace gdenr = 3236 if UCID_1934_2003 == 298328
replace gdenr_2012 = 3236 if UCID_1934_2003 == 298328
replace gdenr_2018 = 3236 if UCID_1934_2003 ==298328
	
replace gdename ="widnau" if UCID_1934_2003 == 298910
replace gdenr = 3238 if UCID_1934_2003 == 298910
replace gdenr_2012 = 3238 if UCID_1934_2003 == 298910
replace gdenr_2018 = 3238 if UCID_1934_2003 ==298910

replace gdename ="buchs sg" if UCID_1934_2003 == 300988
replace gdenr = 3271 if UCID_1934_2003 == 300988
replace gdenr_2012 = 3271 if UCID_1934_2003 == 300988
replace gdenr_2018 = 3271 if UCID_1934_2003 ==300988
	
replace gdename ="vilters wangs" if UCID_1934_2003 == 303619
replace gdenr = 3297 if UCID_1934_2003 == 303619
replace gdenr_2012 = 3297 if UCID_1934_2003 == 303619
replace gdenr_2018 = 3297 if UCID_1934_2003 ==303619
	
replace gdename ="schmerikon" if UCID_1934_2003 == 304199
replace gdenr = 3338 if UCID_1934_2003 == 304199
replace gdenr_2012 = 3338 if UCID_1934_2003 == 304199
replace gdenr_2018 = 3338 if UCID_1934_2003 ==304199

replace gdename ="uznach" if UCID_1934_2003 == 304380
replace gdenr = 3339 if UCID_1934_2003 == 304380
replace gdenr_2012 = 3339 if UCID_1934_2003 == 304380
replace gdenr_2018 = 3339 if UCID_1934_2003 ==304380
				
replace gdename ="ernetschwil" if UCID_1934_2003 == 305741
replace gdenr = 3331 if UCID_1934_2003 == 305741
replace gdenr_2012 = 3331 if UCID_1934_2003 == 305741
replace gdenr_2018 = 3341 if UCID_1934_2003 ==305741

replace gdename ="wattwil" if UCID_1934_2003 == 306622
replace gdenr = 3337 if UCID_1934_2003 == 306622
replace gdenr_2012 = 3337 if UCID_1934_2003 == 306622
replace gdenr_2018 = 3379 if UCID_1934_2003 == 306622
				
replace gdename ="wil (sg)" if UCID_1934_2003 == 309286
replace gdenr = 3425 if UCID_1934_2003 == 309286
replace gdenr_2012 = 3425 if UCID_1934_2003 == 309286
replace gdenr_2018 = 3427 if UCID_1934_2003 ==309286

replace gdename ="wilen bei wil" if UCID_1934_2003 == 309327
replace gdenr = 4752 if UCID_1934_2003 == 309327 & year < 1998
replace gdenr = 4786 if UCID_1934_2003 == 309327 & year > 1997
replace gdenr_2012 = 4786 if UCID_1934_2003 == 309327
replace gdenr_2018 = 4786 if UCID_1934_2003 ==309327

replace gdename ="abtwil" if UCID_1934_2003 == 310426
replace gdenr = 4221 if UCID_1934_2003 == 310426
replace gdenr_2012 = 4221 if UCID_1934_2003 == 310426
replace gdenr_2018 = 4221 if UCID_1934_2003 ==310426
	
replace gdename ="gaiserwald" if UCID_1934_2003 == 310579
replace gdenr = 3442 if UCID_1934_2003 == 310579
replace gdenr_2012 = 3442 if UCID_1934_2003 == 310579
replace gdenr_2018 = 3442 if UCID_1934_2003 ==310579
	
replace gdename ="falera" if UCID_1934_2003 == 312260
replace gdenr = 3572 if UCID_1934_2003 == 312260
replace gdenr_2012 = 3572 if UCID_1934_2003 == 312260
replace gdenr_2018 = 3572 if UCID_1934_2003 ==312260
				
replace gdename ="ilanz" if UCID_1934_2003 == 312536
replace gdenr = 3574 if UCID_1934_2003 == 312536
replace gdenr_2012 = 3574 if UCID_1934_2003 == 312536
replace gdenr_2018 = 3619 if UCID_1934_2003 ==312536

replace gdename ="chur" if UCID_1934_2003 == 312741
replace gdenr = 3901 if UCID_1934_2003 == 312741
replace gdenr_2012 = 3901 if UCID_1934_2003 == 312741
replace gdenr_2018 = 3901 if UCID_1934_2003 ==312741
				
replace gdename ="mutten" if UCID_1934_2003 == 312869
replace gdenr = 3503 if UCID_1934_2003 == 312869
replace gdenr_2012 = 3503 if UCID_1934_2003 == 312869
replace gdenr_2018 = 3668 if UCID_1934_2003 ==312869

replace gdename ="muttenz" if UCID_1934_2003 == 312966
replace gdenr = 2811 if UCID_1934_2003 == 312966 & year < 1994
replace gdenr = 2770 if UCID_1934_2003 == 312966 & year > 1993
replace gdenr_2012 = 2770 if UCID_1934_2003 == 312966
replace gdenr_2018 = 2770 if UCID_1934_2003 ==312966
	
replace gdename ="scuol/schuls" if UCID_1934_2003 == 313802
replace gdenr = 3762 if UCID_1934_2003 == 313802
replace gdenr_2012 = 3762 if UCID_1934_2003 == 313802
replace gdenr_2018 = 3762 if UCID_1934_2003 ==313802

replace gdename ="saint maurice" if UCID_1934_2003 == 315192
replace gdenr = 6217 if UCID_1934_2003 == 315192
replace gdenr_2012 = 6217 if UCID_1934_2003 == 315192
replace gdenr_2018 = 6217 if UCID_1934_2003 ==315192
	
replace gdename ="jenaz" if UCID_1934_2003 == 320568
replace gdenr = 3863 if UCID_1934_2003 == 320568
replace gdenr_2012 = 3863 if UCID_1934_2003 == 320568
replace gdenr_2018 = 3863 if UCID_1934_2003 ==320568

replace gdename ="chur" if UCID_1934_2003 == 321086
replace gdenr = 3901 if UCID_1934_2003 == 321086
replace gdenr_2012 = 3901 if UCID_1934_2003 == 321086
replace gdenr_2018 = 3901 if UCID_1934_2003 ==321086

replace gdename ="chur" if UCID_1934_2003 == 324414
replace gdenr = 3901 if UCID_1934_2003 == 324414
replace gdenr_2012 = 3901 if UCID_1934_2003 == 324414
replace gdenr_2018 = 3901 if UCID_1934_2003 ==324414

replace gdename ="chur" if UCID_1934_2003 == 326257
replace gdenr = 3901 if UCID_1934_2003 == 326257
replace gdenr_2012 = 3901 if UCID_1934_2003 == 326257
replace gdenr_2018 = 3901 if UCID_1934_2003 ==326257

replace gdename ="malix" if UCID_1934_2003 == 327884
replace gdenr = 3912 if UCID_1934_2003 == 327884
replace gdenr_2012 = 3912 if UCID_1934_2003 == 327884
replace gdenr_2018 = 3911 if UCID_1934_2003 ==327884

replace gdename ="malix" if UCID_1934_2003 == 327899
replace gdenr = 3912 if UCID_1934_2003 == 327899
replace gdenr_2012 = 3912 if UCID_1934_2003 == 327899
replace gdenr_2018 = 3911 if UCID_1934_2003 ==327899

replace gdename ="haldenstein" if UCID_1934_2003 == 328393
replace gdenr = 3941 if UCID_1934_2003 == 328393
replace gdenr_2012 = 3941 if UCID_1934_2003 == 328393
replace gdenr_2018 = 3941 if UCID_1934_2003 ==328393
		
replace gdename ="suhr" if UCID_1934_2003 == 330785
replace gdenr = 4012 if UCID_1934_2003 == 330785
replace gdenr_2012 = 4012 if UCID_1934_2003 == 330785
replace gdenr_2018 = 4012 if UCID_1934_2003 ==330785

replace gdename ="buchs ag" if UCID_1934_2003 == 331173
replace gdenr = 4003 if UCID_1934_2003 == 331173
replace gdenr_2012 =4003 if UCID_1934_2003 == 331173
replace gdenr_2018 =4003 if UCID_1934_2003 ==331173

replace gdename ="graenichen" if UCID_1934_2003 == 332195
replace gdenr = 4006 if UCID_1934_2003 == 332195
replace gdenr_2012 = 4006 if UCID_1934_2003 == 332195
replace gdenr_2018 = 4006 if UCID_1934_2003 ==332195

*Tag wrong match
replace separate = 1 if UCID_1934_2003 ==11790
replace separate = 1 if UCID_1934_2003 ==13665
replace separate = 1 if UCID_1934_2003 ==20200
replace separate = 1 if UCID_1934_2003 ==23576
replace separate = 1 if UCID_1934_2003 ==25358
replace separate = 1 if UCID_1934_2003 ==25601
replace separate = 1 if UCID_1934_2003 ==26615
replace separate = 1 if UCID_1934_2003 ==27048
replace separate = 1 if UCID_1934_2003 ==30396
replace separate = 1 if UCID_1934_2003 ==30720
replace separate = 1 if UCID_1934_2003 ==32160
replace separate = 1 if UCID_1934_2003 ==34954
replace separate = 1 if UCID_1934_2003 ==44292
replace separate = 1 if UCID_1934_2003 ==47020
replace separate = 1 if UCID_1934_2003 ==50228
replace separate = 1 if UCID_1934_2003 ==50230
replace separate = 1 if UCID_1934_2003 ==58438
replace separate = 1 if UCID_1934_2003 ==62264
replace separate = 1 if UCID_1934_2003 ==63833
replace separate = 1 if UCID_1934_2003 ==65621
replace separate = 1 if UCID_1934_2003 ==69452
replace separate = 1 if UCID_1934_2003 ==70484
replace separate = 1 if UCID_1934_2003 ==75529
replace separate = 1 if UCID_1934_2003 ==76800
replace separate = 1 if UCID_1934_2003 ==79727
replace separate = 1 if UCID_1934_2003 ==89993
replace separate = 1 if UCID_1934_2003 ==101347
replace separate = 1 if UCID_1934_2003 ==103318
replace separate = 1 if UCID_1934_2003 ==104516
replace separate = 1 if UCID_1934_2003 ==105729
replace separate = 1 if UCID_1934_2003 ==106199
replace separate = 1 if UCID_1934_2003 ==106637
replace separate = 1 if UCID_1934_2003 ==108515
replace separate = 1 if UCID_1934_2003 ==110216
replace separate = 1 if UCID_1934_2003 ==112801
replace separate = 1 if UCID_1934_2003 ==115014

*Correct, just collapse
replace separate = 0 if UCID_1934_2003 == 27181

*****
*Emilie
*Corrections

replace gdename ="affoltern am albis" if UCID_1934_2003 == 294
replace gdenr = 2 if UCID_1934_2003 == 294
replace gdenr_2012 = 2 if UCID_1934_2003 == 294
replace gdenr_2018 = 2 if UCID_1934_2003 == 294

replace gdename ="winterthur" if UCID_1934_2003 == 915200
replace gdenr = 230 if UCID_1934_2003 == 915200
replace gdenr_2012 = 230 if UCID_1934_2003 == 915200
replace gdenr_2018 = 230 if UCID_1934_2003 == 915200

replace gdename ="zuerich" if UCID_1934_2003 == 1288091
replace gdenr = 253 if UCID_1934_2003 == 1288091 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 1288091 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 1288091
replace gdenr_2018 = 261 if UCID_1934_2003 == 1288091

replace gdename ="zuerich" if UCID_1934_2003 == 1147963
replace gdenr = 253 if UCID_1934_2003 == 1147963 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 1147963 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 1147963
replace gdenr_2018 = 261 if UCID_1934_2003 == 1147963

replace gdename ="zuerich" if UCID_1934_2003 == 1109141
replace gdenr = 253 if UCID_1934_2003 == 1109141 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 1109141 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 1109141
replace gdenr_2018 = 261 if UCID_1934_2003 == 1109141

replace gdename ="zuerich" if UCID_1934_2003 == 1015007
replace gdenr = 253 if UCID_1934_2003 == 1015007 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 1015007 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 1015007
replace gdenr_2018 = 261 if UCID_1934_2003 == 1015007

replace gdename ="zuerich" if UCID_1934_2003 == 966394
replace gdenr = 253 if UCID_1934_2003 == 966394 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 966394 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 966394
replace gdenr_2018 = 261 if UCID_1934_2003 == 966394

replace gdename ="zuerich" if UCID_1934_2003 == 905254
replace gdenr = 253 if UCID_1934_2003 == 905254 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 905254 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 905254
replace gdenr_2018 = 261 if UCID_1934_2003 == 905254

replace gdename ="zuerich" if UCID_1934_2003 == 580745580745
replace gdenr = 253 if UCID_1934_2003 == 580745 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 580745 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 580745
replace gdenr_2018 = 261 if UCID_1934_2003 == 580745

replace gdename ="zuerich" if UCID_1934_2003 == 571533
replace gdenr = 253 if UCID_1934_2003 == 571533 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 571533 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 571533
replace gdenr_2018 = 261 if UCID_1934_2003 == 571533

replace gdename ="bern" if UCID_1934_2003 == 915460
replace gdenr = 351 if UCID_1934_2003 == 915460
replace gdenr_2012 = 351 if UCID_1934_2003 == 915460
replace gdenr_2018 = 351 if UCID_1934_2003 == 915460

replace gdename ="interlaken" if UCID_1934_2003 == 540027
replace gdenr = 581 if UCID_1934_2003 == 540027
replace gdenr_2012 = 581 if UCID_1934_2003 == 540027
replace gdenr_2018 = 581 if UCID_1934_2003 == 540027

replace gdename ="littau" if UCID_1934_2003 == 1969713
replace gdenr = 1060 if UCID_1934_2003 == 1969713
replace gdenr_2012 = 1061 if UCID_1934_2003 == 1969713
replace gdenr_2018 = 1061 if UCID_1934_2003 == 1969713

replace gdename ="littau" if UCID_1934_2003 == 1969170
replace gdenr = 1060 if UCID_1934_2003 == 1969170
replace gdenr_2012 = 1061 if UCID_1934_2003 == 1969170
replace gdenr_2018 = 1061 if UCID_1934_2003 == 1969170

replace gdename ="littau" if UCID_1934_2003 == 1968947
replace gdenr = 1060 if UCID_1934_2003 == 1968947
replace gdenr_2012 = 1061 if UCID_1934_2003 == 1968947
replace gdenr_2018 = 1061 if UCID_1934_2003 == 1968947

replace gdename ="littau" if UCID_1934_2003 == 1968187
replace gdenr = 1060 if UCID_1934_2003 == 1968187
replace gdenr_2012 = 1061 if UCID_1934_2003 == 1968187
replace gdenr_2018 = 1061 if UCID_1934_2003 == 1968187

replace gdename ="littau" if UCID_1934_2003 == 601360
replace gdenr = 1060 if UCID_1934_2003 == 601360
replace gdenr_2012 = 1061 if UCID_1934_2003 == 601360
replace gdenr_2018 = 1061 if UCID_1934_2003 == 601360

replace gdename ="littau" if UCID_1934_2003 == 601147
replace gdenr = 1060 if UCID_1934_2003 == 601147
replace gdenr_2012 = 1061 if UCID_1934_2003 == 601147
replace gdenr_2018 = 1061 if UCID_1934_2003 == 601147

replace gdename ="littau" if UCID_1934_2003 == 600350
replace gdenr = 1060 if UCID_1934_2003 == 600350
replace gdenr_2012 = 1061 if UCID_1934_2003 == 600350
replace gdenr_2018 = 1061 if UCID_1934_2003 == 600350

replace gdename ="littau" if UCID_1934_2003 == 600037
replace gdenr = 1060 if UCID_1934_2003 == 600037
replace gdenr_2012 = 1061 if UCID_1934_2003 == 600037
replace gdenr_2018 = 1061 if UCID_1934_2003 == 600037

replace gdename ="weggis" if UCID_1934_2003 == 541946
replace gdenr = 1069 if UCID_1934_2003 == 541946
replace gdenr_2012 = 1069 if UCID_1934_2003 == 541946
replace gdenr_2018 = 1069 if UCID_1934_2003 == 541946
	
replace gdename ="glarus" if UCID_1934_2003 == 1070761
replace gdenr = 1609 if UCID_1934_2003 == 1070761
replace gdenr_2012 = 1632 if UCID_1934_2003 == 1070761
replace gdenr_2018 = 1632 if UCID_1934_2003 == 1070761

replace gdename ="zug" if UCID_1934_2003 == 1164946
replace gdenr = 1711 if UCID_1934_2003 == 1164946
replace gdenr_2012 = 1711 if UCID_1934_2003 == 1164946
replace gdenr_2018 = 1711 if UCID_1934_2003 == 1164946

replace gdename ="zug" if UCID_1934_2003 == 874950
replace gdenr = 1711 if UCID_1934_2003 == 874950
replace gdenr_2012 = 1711 if UCID_1934_2003 == 874950
replace gdenr_2018 = 1711 if UCID_1934_2003 == 874950

replace gdename ="zug" if UCID_1934_2003 == 815572
replace gdenr = 1711 if UCID_1934_2003 == 815572
replace gdenr_2012 = 1711 if UCID_1934_2003 == 815572
replace gdenr_2018 = 1711 if UCID_1934_2003 == 815572
	
replace gdename ="zug" if UCID_1934_2003 == 547002
replace gdenr = 1711 if UCID_1934_2003 == 547002
replace gdenr_2012 = 1711 if UCID_1934_2003 == 547002
replace gdenr_2018 = 1711 if UCID_1934_2003 == 547002
		
replace gdename ="zug" if UCID_1934_2003 == 546157
replace gdenr = 1711 if UCID_1934_2003 == 546157
replace gdenr_2012 = 1711 if UCID_1934_2003 == 546157
replace gdenr_2018 = 1711 if UCID_1934_2003 == 546157

replace gdename ="zug" if UCID_1934_2003 == 546060	
replace gdenr = 1711 if UCID_1934_2003 == 546060	
replace gdenr_2012 = 1711 if UCID_1934_2003 == 546060	
replace gdenr_2018 = 1711 if UCID_1934_2003 == 546060	
	
replace gdename ="zug" if UCID_1934_2003 == 508796
replace gdenr = 1711 if UCID_1934_2003 == 508796
replace gdenr_2012 = 1711 if UCID_1934_2003 == 508796
replace gdenr_2018 = 1711 if UCID_1934_2003 == 508796

replace gdename ="bulle" if UCID_1934_2003 == 480627
replace gdenr = 2125 if UCID_1934_2003 == 480627
replace gdenr_2012 = 2125 if UCID_1934_2003 == 480627
replace gdenr_2018 = 2125 if UCID_1934_2003 == 480627

replace gdename ="fribourg" if UCID_1934_2003 == 859503
replace gdenr = 2196 if UCID_1934_2003 == 859503
replace gdenr_2012 = 2196 if UCID_1934_2003 == 859503
replace gdenr_2018 = 2196 if UCID_1934_2003 == 859503

replace gdename ="fribourg" if UCID_1934_2003 == 784019
replace gdenr = 2196 if UCID_1934_2003 == 784019
replace gdenr_2012 = 2196 if UCID_1934_2003 == 784019
replace gdenr_2018 = 2196 if UCID_1934_2003 == 784019

replace gdename ="heitenried" if UCID_1934_2003 == 1093450
replace gdenr = 2296 if UCID_1934_2003 == 1093450
replace gdenr_2012 = 2296 if UCID_1934_2003 == 1093450
replace gdenr_2018 = 2296 if UCID_1934_2003 == 1093450
	
replace gdename ="basel" if UCID_1934_2003 == 1138058
replace gdenr = 2701 if UCID_1934_2003 == 1138058
replace gdenr_2012 = 2701 if UCID_1934_2003 == 1138058
replace gdenr_2018 = 2701 if UCID_1934_2003 == 1138058
	
replace gdename ="biel benken" if UCID_1934_2003 == 543157
replace gdenr = 2805 if UCID_1934_2003 == 543157 & year < 1994
replace gdenr = 2764 if UCID_1934_2003 == 543157 & year > 1993
replace gdenr_2012 = 2764 if UCID_1934_2003 == 543157
replace gdenr_2018 = 2764 if UCID_1934_2003 == 543157

replace gdename ="vilters wangs" if UCID_1934_2003 == 480721
replace gdenr = 3297 if UCID_1934_2003 == 480721
replace gdenr_2012 = 3297 if UCID_1934_2003 == 480721
replace gdenr_2018 = 3297 if UCID_1934_2003 == 480721
	
replace gdename ="vilters" if UCID_1934_2003 == 469797
replace gdenr = 3297 if UCID_1934_2003 == 469797
replace gdenr_2012 = 3297 if UCID_1934_2003 == 469797
replace gdenr_2018 = 3297 if UCID_1934_2003 == 469797

replace gdename ="rapperswil sg" if UCID_1934_2003 == 537974
replace gdenr = 3336 if UCID_1934_2003 == 537974
replace gdenr_2012 = 3340 if UCID_1934_2003 == 537974
replace gdenr_2018 = 3340 if UCID_1934_2003 == 537974

replace gdename ="st. moritz" if UCID_1934_2003 == 1097744
replace gdenr = 3787 if UCID_1934_2003 == 1097744
replace gdenr_2012 = 3787 if UCID_1934_2003 == 1097744
replace gdenr_2018 = 3787 if UCID_1934_2003 == 1097744

replace gdename ="st antoenien" if UCID_1934_2003 == 902581
replace gdenr = 3893 if UCID_1934_2003 == 902581
replace gdenr_2012 = 3893 if UCID_1934_2003 == 902581
replace gdenr_2018 = 3891 if UCID_1934_2003 == 902581
	
replace gdename ="chur" if UCID_1934_2003 == 480318
replace gdenr = 3901 if UCID_1934_2003 == 480318
replace gdenr_2012 = 3901 if UCID_1934_2003 == 480318
replace gdenr_2018 = 3901 if UCID_1934_2003 == 480318

replace gdename ="bremgarten (ag)" if UCID_1934_2003 == 1367428
replace gdenr = 4063 if UCID_1934_2003 == 1367428
replace gdenr_2012 = 4063 if UCID_1934_2003 == 1367428
replace gdenr_2018 = 4063 if UCID_1934_2003 == 1367428

replace gdename ="hermetschwil staffeln" if UCID_1934_2003 == 1155341
replace gdenr = 4069 if UCID_1934_2003 == 1155341
replace gdenr_2012 = 4069 if UCID_1934_2003 == 1155341
replace gdenr_2018 = 4063 if UCID_1934_2003 == 1155341

replace gdename ="hermetschwil staffeln" if UCID_1934_2003 == 1098803
replace gdenr = 4069 if UCID_1934_2003 == 1098803
replace gdenr_2012 = 4069 if UCID_1934_2003 == 1098803
replace gdenr_2018 = 4063 if UCID_1934_2003 == 1098803

replace gdename ="niederwil (ag)" if UCID_1934_2003 == 458439
replace gdenr = 4072 if UCID_1934_2003 == 458439
replace gdenr_2012 = 4072 if UCID_1934_2003 == 458439
replace gdenr_2018 = 4072 if UCID_1934_2003 == 458439

replace gdename ="wohlen ag" if UCID_1934_2003 == 550035
replace gdenr = 4082 if UCID_1934_2003 == 550035
replace gdenr_2012 = 4082 if UCID_1934_2003 == 550035
replace gdenr_2018 = 4082 if UCID_1934_2003 == 550035

replace gdename ="wohlen ag" if UCID_1934_2003 == 546829
replace gdenr = 4082 if UCID_1934_2003 == 546829
replace gdenr_2012 = 4082 if UCID_1934_2003 == 546829
replace gdenr_2018 = 4082 if UCID_1934_2003 == 546829

replace gdename ="wohlen ag" if UCID_1934_2003 == 544988
replace gdenr = 4082 if UCID_1934_2003 == 544988
replace gdenr_2012 = 4082 if UCID_1934_2003 == 544988
replace gdenr_2018 = 4082 if UCID_1934_2003 == 544988

replace gdename ="wohlen ag" if UCID_1934_2003 == 543693
replace gdenr = 4082 if UCID_1934_2003 == 543693
replace gdenr_2012 = 4082 if UCID_1934_2003 == 543693
replace gdenr_2018 = 4082 if UCID_1934_2003 == 543693

replace gdename ="wohlen ag" if UCID_1934_2003 == 542879
replace gdenr = 4082 if UCID_1934_2003 == 542879
replace gdenr_2012 = 4082 if UCID_1934_2003 == 542879
replace gdenr_2018 = 4082 if UCID_1934_2003 == 542879

replace gdename ="wohlen ag" if UCID_1934_2003 == 539907
replace gdenr = 4082 if UCID_1934_2003 == 539907
replace gdenr_2012 = 4082 if UCID_1934_2003 == 539907
replace gdenr_2018 = 4082 if UCID_1934_2003 == 539907
	
replace gdename ="birr" if UCID_1934_2003 == 458230
replace gdenr = 4092 if UCID_1934_2003 == 458230
replace gdenr_2012 = 4092 if UCID_1934_2003 == 458230
replace gdenr_2018 = 4092 if UCID_1934_2003 == 458230

replace gdename ="birr" if UCID_1934_2003 == 458227
replace gdenr = 4092 if UCID_1934_2003 == 458227
replace gdenr_2012 = 4092 if UCID_1934_2003 == 458227
replace gdenr_2018 = 4092 if UCID_1934_2003 == 458227
	
replace gdename ="lupfig" if UCID_1934_2003 == 488028
replace gdenr = 4104 if UCID_1934_2003 == 488028
replace gdenr_2012 = 4104 if UCID_1934_2003 == 488028
replace gdenr_2018 = 4104 if UCID_1934_2003 == 488028
	
replace gdename ="lenzburg" if UCID_1934_2003 == 962008
replace gdenr = 4201 if UCID_1934_2003 == 962008
replace gdenr_2012 = 4201 if UCID_1934_2003 == 962008
replace gdenr_2018 = 4201 if UCID_1934_2003 == 962008

replace gdename ="zihlschlacht" if UCID_1934_2003 == 1193596
replace gdenr = 4518 if UCID_1934_2003 == 1193596
replace gdenr_2012 = 4511 if UCID_1934_2003 == 1193596
replace gdenr_2018 = 4511 if UCID_1934_2003 == 1193596

replace gdename ="guendelhart hoerhausen" if UCID_1934_2003 == 1123290
replace gdenr = 4861 if UCID_1934_2003 == 1123290
replace gdenr_2012 = 4816 if UCID_1934_2003 == 1123290
replace gdenr_2018 = 4816 if UCID_1934_2003 == 1123290
				
replace gdename ="nussbaumen" if UCID_1934_2003 == 910298
replace gdenr = 4822 if UCID_1934_2003 == 910298
replace gdenr_2012 = 4821 if UCID_1934_2003 == 910298
replace gdenr_2018 = 4821 if UCID_1934_2003 == 910298
		
replace gdename ="giubiasco" if UCID_1934_2003 == 1366263
replace gdenr = 5005 if UCID_1934_2003 == 1366263
replace gdenr_2012 = 5005 if UCID_1934_2003 == 1366263
replace gdenr_2018 = 5002 if UCID_1934_2003 == 1366263
	
replace gdename ="sementina" if UCID_1934_2003 == 1215534
replace gdenr = 5019 if UCID_1934_2003 == 1215534
replace gdenr_2012 = 5019 if UCID_1934_2003 == 1215534
replace gdenr_2018 = 5002 if UCID_1934_2003 == 1215534
	
replace gdename ="bellinzona" if UCID_1934_2003 == 888379
replace gdenr = 5002 if UCID_1934_2003 == 888379
replace gdenr_2012 = 5002 if UCID_1934_2003 == 888379
replace gdenr_2018 = 5002 if UCID_1934_2003 == 888379

replace gdename ="carabbia" if UCID_1934_2003 == 1321091
replace gdenr = 5168 if UCID_1934_2003 == 1321091
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1321091
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1321091

replace gdename ="carabbia" if UCID_1934_2003 == 1350320
replace gdenr = 5168 if UCID_1934_2003 == 1350320
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1350320
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1350320

replace gdename ="castagnola" if UCID_1934_2003 == 1340515
replace gdenr = 5172 if UCID_1934_2003 == 1340515
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1340515
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1340515

replace gdename ="castagnola" if UCID_1934_2003 == 1332210
replace gdenr = 5172 if UCID_1934_2003 == 1332210
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1332210
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1332210

replace gdename ="castagnola" if UCID_1934_2003 == 1321135
replace gdenr = 5172 if UCID_1934_2003 == 1321135
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1321135
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1321135

replace gdename ="castagnola" if UCID_1934_2003 == 1123793
replace gdenr = 5172 if UCID_1934_2003 == 1123793
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1123793
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1123793

replace gdename ="castagnola" if UCID_1934_2003 == 1350421
replace gdenr = 5172 if UCID_1934_2003 == 1350421
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1350421
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1350421

replace gdename ="castagnola" if UCID_1934_2003 == 974891
replace gdenr = 5172 if UCID_1934_2003 == 974891
replace gdenr_2012 = 5192 if UCID_1934_2003 == 974891
replace gdenr_2018 = 5192 if UCID_1934_2003 == 974891

replace gdename ="viganello" if UCID_1934_2003 == 1320752
replace gdenr = 5234 if UCID_1934_2003 == 1320752
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1320752
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1320752

replace gdename ="viganello" if UCID_1934_2003 == 1100995
replace gdenr = 5234 if UCID_1934_2003 == 1100995
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1100995
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1100995

replace gdename ="viganello" if UCID_1934_2003 == 1303896
replace gdenr = 5234 if UCID_1934_2003 == 1303896
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1303896
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1303896

replace gdename ="viganello" if UCID_1934_2003 == 511064
replace gdenr = 5234 if UCID_1934_2003 == 511064
replace gdenr_2012 = 5192 if UCID_1934_2003 == 511064
replace gdenr_2018 = 5192 if UCID_1934_2003 == 511064

replace gdename ="breganzona" if UCID_1934_2003 == 1299007
replace gdenr = 5158 if UCID_1934_2003 == 1299007
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1299007
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1299007
	
replace gdename ="pambio noranco" if UCID_1934_2003 == 1225563
replace gdenr = 5209 if UCID_1934_2003 == 1225563
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1225563
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1225563

replace gdename ="pambio noranco" if UCID_1934_2003 == 1077437
replace gdenr = 5209 if UCID_1934_2003 == 1077437
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1077437
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1077437

replace gdename ="barbengo" if UCID_1934_2003 == 1215928
replace gdenr = 5147 if UCID_1934_2003 == 1215928
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1215928
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1215928

replace gdename ="lugano" if UCID_1934_2003 == 1215797
replace gdenr = 5192 if UCID_1934_2003 == 1215797
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1215797
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1215797

replace gdename ="lugano" if UCID_1934_2003 == 1142158
replace gdenr = 5192 if UCID_1934_2003 == 1142158
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1142158
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1142158

replace gdename ="lugano" if UCID_1934_2003 == 1100919
replace gdenr = 5192 if UCID_1934_2003 == 1100919
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1100919
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1100919

replace gdename ="lugano" if UCID_1934_2003 == 1043352
replace gdenr = 5192 if UCID_1934_2003 == 1043352
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1043352
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1043352

replace gdename ="lugano" if UCID_1934_2003 == 967763
replace gdenr = 5192 if UCID_1934_2003 == 967763
replace gdenr_2012 = 5192 if UCID_1934_2003 == 967763
replace gdenr_2018 = 5192 if UCID_1934_2003 == 967763

replace gdename ="lugano" if UCID_1934_2003 == 806018
replace gdenr = 5192 if UCID_1934_2003 == 806018
replace gdenr_2012 = 5192 if UCID_1934_2003 == 806018
replace gdenr_2018 = 5192 if UCID_1934_2003 == 806018

replace gdename ="lugano" if UCID_1934_2003 == 492304
replace gdenr = 5192 if UCID_1934_2003 == 492304
replace gdenr_2012 = 5192 if UCID_1934_2003 == 492304
replace gdenr_2018 = 5192 if UCID_1934_2003 == 492304

replace gdename ="lugano" if UCID_1934_2003 == 440777
replace gdenr = 5192 if UCID_1934_2003 == 440777
replace gdenr_2012 = 5192 if UCID_1934_2003 == 440777
replace gdenr_2018 = 5192 if UCID_1934_2003 == 440777

replace gdename ="pregassona" if UCID_1934_2003 == 1009325
replace gdenr = 5215 if UCID_1934_2003 == 1009325
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1009325
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1009325

replace gdename ="pregassona" if UCID_1934_2003 == 687832
replace gdenr = 5215 if UCID_1934_2003 == 687832
replace gdenr_2012 = 5192 if UCID_1934_2003 == 687832
replace gdenr_2018 = 5192 if UCID_1934_2003 == 687832
	
replace gdename ="manno" if UCID_1934_2003 == 871302
replace gdenr = 5194 if UCID_1934_2003 == 871302
replace gdenr_2012 = 5194 if UCID_1934_2003 == 871302
replace gdenr_2018 = 5194 if UCID_1934_2003 == 871302

replace gdename ="chiasso" if UCID_1934_2003 == 518320
replace gdenr = 5250 if UCID_1934_2003 == 518320
replace gdenr_2012 = 5250 if UCID_1934_2003 == 518320
replace gdenr_2018 = 5250 if UCID_1934_2003 == 518320

replace gdename ="bex" if UCID_1934_2003 == 458774
replace gdenr = 5402 if UCID_1934_2003 == 458774
replace gdenr_2012 = 5402 if UCID_1934_2003 == 458774
replace gdenr_2018 = 5402 if UCID_1934_2003 == 458774

replace gdename ="ollon" if UCID_1934_2003 == 546934
replace gdenr = 5409 if UCID_1934_2003 == 546934
replace gdenr_2012 = 5409 if UCID_1934_2003 == 546934
replace gdenr_2018 = 5409 if UCID_1934_2003 == 546934

replace gdename ="crissier" if UCID_1934_2003 == 493289
replace gdenr = 5583 if UCID_1934_2003 == 493289
replace gdenr_2012 = 5583 if UCID_1934_2003 == 493289
replace gdenr_2018 = 5583 if UCID_1934_2003 == 493289

replace gdename ="crissier" if UCID_1934_2003 == 444392
replace gdenr = 5583 if UCID_1934_2003 == 444392
replace gdenr_2012 = 5583 if UCID_1934_2003 == 444392
replace gdenr_2018 = 5583 if UCID_1934_2003 == 444392

replace gdename ="lausanne" if UCID_1934_2003 == 1366350
replace gdenr = 5586 if UCID_1934_2003 == 1366350
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1366350
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1366350

replace gdename ="lausanne" if UCID_1934_2003 == 1365214
replace gdenr = 5586 if UCID_1934_2003 == 1365214
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1365214
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1365214

replace gdename ="lausanne" if UCID_1934_2003 == 1365165
replace gdenr = 5586 if UCID_1934_2003 == 1365165
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1365165
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1365165

replace gdename ="lausanne" if UCID_1934_2003 == 1362783
replace gdenr = 5586 if UCID_1934_2003 == 1362783
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1362783
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1362783

replace gdename ="lausanne" if UCID_1934_2003 == 1359567
replace gdenr = 5586 if UCID_1934_2003 == 1359567
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1359567
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1359567

replace gdename ="lausanne" if UCID_1934_2003 == 1341637
replace gdenr = 5586 if UCID_1934_2003 == 1341637
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1341637
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1341637

replace gdename ="lausanne" if UCID_1934_2003 == 1341172
replace gdenr = 5586 if UCID_1934_2003 == 1341172
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1341172
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1341172

replace gdename ="lausanne" if UCID_1934_2003 == 1292656
replace gdenr = 5586 if UCID_1934_2003 == 1292656
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1292656
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1292656

replace gdename ="lausanne" if UCID_1934_2003 == 1254675
replace gdenr = 5586 if UCID_1934_2003 == 1254675
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1254675
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1254675

replace gdename ="lausanne" if UCID_1934_2003 == 1237153
replace gdenr = 5586 if UCID_1934_2003 == 1237153
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1237153
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1237153

replace gdename ="lausanne" if UCID_1934_2003 == 1143244
replace gdenr = 5586 if UCID_1934_2003 == 1143244
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1143244
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1143244

replace gdename ="lausanne" if UCID_1934_2003 == 1102394
replace gdenr = 5586 if UCID_1934_2003 == 1102394
replace gdenr_2012 = 5586 if UCID_1934_2003 == 1102394
replace gdenr_2018 = 5586 if UCID_1934_2003 == 1102394

replace gdename ="lausanne" if UCID_1934_2003 == 975647
replace gdenr = 5586 if UCID_1934_2003 == 975647
replace gdenr_2012 = 5586 if UCID_1934_2003 == 975647
replace gdenr_2018 = 5586 if UCID_1934_2003 == 975647

replace gdename ="lausanne" if UCID_1934_2003 == 965424
replace gdenr = 5586 if UCID_1934_2003 == 965424
replace gdenr_2012 = 5586 if UCID_1934_2003 == 965424
replace gdenr_2018 = 5586 if UCID_1934_2003 == 965424

replace gdename ="lausanne" if UCID_1934_2003 == 510162
replace gdenr = 5586 if UCID_1934_2003 == 510162
replace gdenr_2012 = 5586 if UCID_1934_2003 == 510162
replace gdenr_2018 = 5586 if UCID_1934_2003 == 510162
	
replace gdename ="lausanne" if UCID_1934_2003 == 436207
replace gdenr = 5586 if UCID_1934_2003 == 436207
replace gdenr_2012 = 5586 if UCID_1934_2003 == 436207
replace gdenr_2018 = 5586 if UCID_1934_2003 == 436207

replace gdename ="pully" if UCID_1934_2003 == 1341719
replace gdenr = 5590 if UCID_1934_2003 == 1341719
replace gdenr_2012 = 5590 if UCID_1934_2003 == 1341719
replace gdenr_2018 = 5590 if UCID_1934_2003 ==1341719

replace gdename ="pully" if UCID_1934_2003 == 1143373
replace gdenr = 5590 if UCID_1934_2003 == 1143373
replace gdenr_2012 = 5590 if UCID_1934_2003 == 1143373
replace gdenr_2018 = 5590 if UCID_1934_2003 == 1143373

replace gdename ="pully" if UCID_1934_2003 == 975673
replace gdenr = 5590 if UCID_1934_2003 == 975673
replace gdenr_2012 = 5590 if UCID_1934_2003 == 975673
replace gdenr_2018 = 5590 if UCID_1934_2003 ==975673

replace gdename ="pully" if UCID_1934_2003 == 871701
replace gdenr = 5590 if UCID_1934_2003 == 871701
replace gdenr_2012 = 5590 if UCID_1934_2003 == 871701
replace gdenr_2018 = 5590 if UCID_1934_2003 ==871701

replace gdename ="pully" if UCID_1934_2003 == 524931
replace gdenr = 5590 if UCID_1934_2003 == 524931
replace gdenr_2012 = 5590 if UCID_1934_2003 == 524931
replace gdenr_2018 = 5590 if UCID_1934_2003 ==524931

replace gdename ="pully" if UCID_1934_2003 == 494605
replace gdenr = 5590 if UCID_1934_2003 == 494605
replace gdenr_2012 = 5590 if UCID_1934_2003 == 494605
replace gdenr_2018 = 5590 if UCID_1934_2003 ==494605

replace gdename ="pully" if UCID_1934_2003 == 459660
replace gdenr = 5590 if UCID_1934_2003 == 459660
replace gdenr_2012 = 5590 if UCID_1934_2003 == 459660
replace gdenr_2018 = 5590 if UCID_1934_2003 ==459660

replace gdename ="pully" if UCID_1934_2003 == 459635
replace gdenr = 5590 if UCID_1934_2003 == 459635
replace gdenr_2012 = 5590 if UCID_1934_2003 == 459635
replace gdenr_2018 = 5590 if UCID_1934_2003 ==459635

replace gdename ="pully" if UCID_1934_2003 == 459624
replace gdenr = 5590 if UCID_1934_2003 == 459624
replace gdenr_2012 = 5590 if UCID_1934_2003 == 459624
replace gdenr_2018 = 5590 if UCID_1934_2003 ==459624

replace gdename ="pully" if UCID_1934_2003 == 459603
replace gdenr = 5590 if UCID_1934_2003 == 459603
replace gdenr_2012 = 5590 if UCID_1934_2003 == 459603
replace gdenr_2018 = 5590 if UCID_1934_2003 ==459603

replace gdename ="pully" if UCID_1934_2003 == 459568
replace gdenr = 5590 if UCID_1934_2003 == 459568
replace gdenr_2012 = 5590 if UCID_1934_2003 == 459568
replace gdenr_2018 = 5590 if UCID_1934_2003 ==459568

replace gdename ="pully" if UCID_1934_2003 == 459439
replace gdenr = 5590 if UCID_1934_2003 == 459439
replace gdenr_2012 = 5590 if UCID_1934_2003 == 459439
replace gdenr_2018 = 5590 if UCID_1934_2003 ==459439

replace gdename ="pully" if UCID_1934_2003 == 459428
replace gdenr = 5590 if UCID_1934_2003 == 459428
replace gdenr_2012 = 5590 if UCID_1934_2003 == 459428
replace gdenr_2018 = 5590 if UCID_1934_2003 ==459428

replace gdename ="pully" if UCID_1934_2003 == 459412
replace gdenr = 5590 if UCID_1934_2003 == 459412
replace gdenr_2012 = 5590 if UCID_1934_2003 == 459412
replace gdenr_2018 = 5590 if UCID_1934_2003 ==459412

replace gdename ="pully" if UCID_1934_2003 == 459356
replace gdenr = 5590 if UCID_1934_2003 == 459356
replace gdenr_2012 = 5590 if UCID_1934_2003 == 459356
replace gdenr_2018 = 5590 if UCID_1934_2003 ==459356

replace gdename ="pully" if UCID_1934_2003 == 459285
replace gdenr = 5590 if UCID_1934_2003 == 459285
replace gdenr_2012 = 5590 if UCID_1934_2003 == 459285
replace gdenr_2018 = 5590 if UCID_1934_2003 ==459285

replace gdename ="pully" if UCID_1934_2003 == 443066
replace gdenr = 5590 if UCID_1934_2003 == 443066
replace gdenr_2012 = 5590 if UCID_1934_2003 == 443066
replace gdenr_2018 = 5590 if UCID_1934_2003 ==443066

replace gdename ="pully" if UCID_1934_2003 == 437247
replace gdenr = 5590 if UCID_1934_2003 == 437247
replace gdenr_2012 = 5590 if UCID_1934_2003 == 437247
replace gdenr_2018 = 5590 if UCID_1934_2003 ==437247

replace gdename ="pully" if UCID_1934_2003 == 437155
replace gdenr = 5590 if UCID_1934_2003 == 437155
replace gdenr_2012 = 5590 if UCID_1934_2003 == 437155
replace gdenr_2018 = 5590 if UCID_1934_2003 ==437155

replace gdename ="pully" if UCID_1934_2003 == 437115
replace gdenr = 5590 if UCID_1934_2003 == 437115
replace gdenr_2012 = 5590 if UCID_1934_2003 == 437115
replace gdenr_2018 = 5590 if UCID_1934_2003 ==437115

replace gdename ="pully" if UCID_1934_2003 == 437072
replace gdenr = 5590 if UCID_1934_2003 == 437072
replace gdenr_2012 = 5590 if UCID_1934_2003 == 437072
replace gdenr_2018 = 5590 if UCID_1934_2003 ==437072

replace gdename ="pully" if UCID_1934_2003 == 437063
replace gdenr = 5590 if UCID_1934_2003 == 437063
replace gdenr_2012 = 5590 if UCID_1934_2003 == 437063
replace gdenr_2018 = 5590 if UCID_1934_2003 ==437063

replace gdename ="pully" if UCID_1934_2003 == 436980
replace gdenr = 5590 if UCID_1934_2003 == 436980
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436980
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436980

replace gdename ="pully" if UCID_1934_2003 == 436972
replace gdenr = 5590 if UCID_1934_2003 == 436972
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436972
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436972

replace gdename ="pully" if UCID_1934_2003 == 436965
replace gdenr = 5590 if UCID_1934_2003 == 436965
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436965
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436965

replace gdename ="pully" if UCID_1934_2003 == 436902
replace gdenr = 5590 if UCID_1934_2003 == 436902
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436902
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436902

replace gdename ="pully" if UCID_1934_2003 == 436757
replace gdenr = 5590 if UCID_1934_2003 == 436757
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436757
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436757

replace gdename ="pully" if UCID_1934_2003 == 436689
replace gdenr = 5590 if UCID_1934_2003 == 436689
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436689
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436689

replace gdename ="pully" if UCID_1934_2003 == 436625
replace gdenr = 5590 if UCID_1934_2003 == 436625
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436625
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436625

replace gdename ="pully" if UCID_1934_2003 == 436550
replace gdenr = 5590 if UCID_1934_2003 == 436550
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436550
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436550

replace gdename ="pully" if UCID_1934_2003 == 436545
replace gdenr = 5590 if UCID_1934_2003 == 436545
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436545
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436545

replace gdename ="pully" if UCID_1934_2003 == 436473
replace gdenr = 5590 if UCID_1934_2003 == 436473
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436473
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436473

replace gdename ="pully" if UCID_1934_2003 == 436427
replace gdenr = 5590 if UCID_1934_2003 == 436427
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436427
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436427

replace gdename ="pully" if UCID_1934_2003 == 436397
replace gdenr = 5590 if UCID_1934_2003 == 436397
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436397
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436397

replace gdename ="pully" if UCID_1934_2003 == 436250
replace gdenr = 5590 if UCID_1934_2003 == 436250
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436250
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436250

replace gdename ="pully" if UCID_1934_2003 ==436220 
replace gdenr = 5590 if UCID_1934_2003 == 436220
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436220
replace gdenr_2018 = 5590 if UCID_1934_2003 ==436220

replace gdename ="pully" if UCID_1934_2003 == 436182
replace gdenr = 5590 if UCID_1934_2003 == 436182
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436182
replace gdenr_2018 = 5590 if UCID_1934_2003 == 436182

replace gdename ="renens (vd)" if UCID_1934_2003 == 437518
replace gdenr = 5591 if UCID_1934_2003 == 437518
replace gdenr_2012 = 5591 if UCID_1934_2003 == 437518
replace gdenr_2018 = 5591 if UCID_1934_2003 == 437518
	
replace gdename ="chavannes pres renens" if UCID_1934_2003 == 439932
replace gdenr = 5627 if UCID_1934_2003 == 439932
replace gdenr_2012 = 5627 if UCID_1934_2003 == 439932
replace gdenr_2018 = 5627 if UCID_1934_2003 == 439932

replace gdename ="morges" if UCID_1934_2003 == 447174
replace gdenr = 5642 if UCID_1934_2003 == 447174
replace gdenr_2012 = 5642 if UCID_1934_2003 == 447174
replace gdenr_2018 = 5642 if UCID_1934_2003 == 447174
	
replace gdename ="moudon" if UCID_1934_2003 == 443467
replace gdenr = 5678 if UCID_1934_2003 == 443467
replace gdenr_2012 = 5678 if UCID_1934_2003 == 443467
replace gdenr_2018 = 5678 if UCID_1934_2003 == 443467
	
replace gdename ="crans pres celigny" if UCID_1934_2003 == 447021
replace gdenr = 5713 if UCID_1934_2003 == 447021
replace gdenr_2012 = 5713 if UCID_1934_2003 == 447021
replace gdenr_2018 = 5713 if UCID_1934_2003 == 447021

replace gdename ="crans pres celigny" if UCID_1934_2003 == 444302
replace gdenr = 5713 if UCID_1934_2003 == 444302
replace gdenr_2012 = 5713 if UCID_1934_2003 == 444302
replace gdenr_2018 = 5713 if UCID_1934_2003 == 444302

replace gdename ="crassier" if UCID_1934_2003 == 444430
replace gdenr = 5714 if UCID_1934_2003 == 444430
replace gdenr_2012 = 5714 if UCID_1934_2003 == 444430
replace gdenr_2018 = 5714 if UCID_1934_2003 == 444430

replace gdename ="nyon" if UCID_1934_2003 == 487642
replace gdenr = 5724 if UCID_1934_2003 == 487642
replace gdenr_2012 = 5724 if UCID_1934_2003 == 487642
replace gdenr_2018 = 5724 if UCID_1934_2003 == 487642

replace gdename ="nyon" if UCID_1934_2003 == 445466
replace gdenr = 5724 if UCID_1934_2003 == 445466
replace gdenr_2012 = 5724 if UCID_1934_2003 == 445466
replace gdenr_2018 = 5724 if UCID_1934_2003 ==  445466

replace gdename ="trelex" if UCID_1934_2003 == 465938
replace gdenr = 5730 if UCID_1934_2003 == 465938
replace gdenr_2012 = 5730 if UCID_1934_2003 == 465938
replace gdenr_2018 = 5730 if UCID_1934_2003 == 465938

replace gdename ="vuiteboeuf" if UCID_1934_2003 == 448193
replace gdenr = 5766 if UCID_1934_2003 == 448193
replace gdenr_2012 = 5766 if UCID_1934_2003 == 448193
replace gdenr_2018 = 5766 if UCID_1934_2003 == 448193	

replace gdename ="vuiteboeuf" if UCID_1934_2003 == 438511
replace gdenr = 5766 if UCID_1934_2003 == 438511
replace gdenr_2012 = 5766 if UCID_1934_2003 == 438511
replace gdenr_2018 = 5766 if UCID_1934_2003 == 438511

replace gdename ="villars tiercelin" if UCID_1934_2003 == 448204
replace gdenr = 5538 if UCID_1934_2003 == 448204
replace gdenr_2012 = 5804 if UCID_1934_2003 == 448204
replace gdenr_2018 = 5804 if UCID_1934_2003 == 448204

replace gdename ="corcelles pres payerne" if UCID_1934_2003 == 448066
replace gdenr = 5816 if UCID_1934_2003 == 448066
replace gdenr_2012 =5816  if UCID_1934_2003 == 448066
replace gdenr_2018 = 5816 if UCID_1934_2003 == 448066

replace gdename ="montreux" if UCID_1934_2003 == 542848
replace gdenr = 5886 if UCID_1934_2003 == 542848
replace gdenr_2012 = 5886 if UCID_1934_2003 == 542848
replace gdenr_2018 = 5886 if UCID_1934_2003 == 542848
	
replace gdename ="montreux" if UCID_1934_2003 == 450625
replace gdenr = 5886 if UCID_1934_2003 == 450625
replace gdenr_2012 = 5886 if UCID_1934_2003 == 450625
replace gdenr_2018 = 5886 if UCID_1934_2003 == 450625

replace gdename ="la tour de peilz" if UCID_1934_2003 == 975846
replace gdenr = 5889 if UCID_1934_2003 == 975846
replace gdenr_2012 = 5889 if UCID_1934_2003 == 975846
replace gdenr_2018 = 5889 if UCID_1934_2003 == 975846
	
replace gdename ="vevey" if UCID_1934_2003 == 453405
replace gdenr = 5890 if UCID_1934_2003 == 453405
replace gdenr_2012 = 5890 if UCID_1934_2003 == 453405
replace gdenr_2018 = 5890 if UCID_1934_2003 == 453405

replace gdename ="brig glis" if UCID_1934_2003 == 456647
replace gdenr = 6002 if UCID_1934_2003 == 456647
replace gdenr_2012 = 6002 if UCID_1934_2003 == 456647
replace gdenr_2018 = 6002 if UCID_1934_2003 ==456647

replace gdename ="brig glis" if UCID_1934_2003 == 456566
replace gdenr = 6002 if UCID_1934_2003 == 456566
replace gdenr_2012 = 6002 if UCID_1934_2003 == 456566
replace gdenr_2018 = 6002 if UCID_1934_2003 ==456566

replace gdename ="brig glis" if UCID_1934_2003 == 456544
replace gdenr = 6002 if UCID_1934_2003 == 456544
replace gdenr_2012 = 6002 if UCID_1934_2003 == 456544
replace gdenr_2018 = 6002 if UCID_1934_2003 == 456544
	
replace gdename ="naters" if UCID_1934_2003 == 456855
replace gdenr = 6007 if UCID_1934_2003 == 456855
replace gdenr_2012 = 6007 if UCID_1934_2003 == 456855
replace gdenr_2018 = 6007 if UCID_1934_2003 == 456855

replace gdename ="fully" if UCID_1934_2003 == 459400
replace gdenr = 6133 if UCID_1934_2003 == 459400
replace gdenr_2012 = 6133 if UCID_1934_2003 == 459400
replace gdenr_2018 = 6133 if UCID_1934_2003 == 459400

replace gdename ="fully" if UCID_1934_2003 == 459367
replace gdenr = 6133 if UCID_1934_2003 == 459367
replace gdenr_2012 = 6133 if UCID_1934_2003 == 459367
replace gdenr_2018 = 6133 if UCID_1934_2003 == 459367
	
replace gdename ="riddes" if UCID_1934_2003 == 460999
replace gdenr = 6139 if UCID_1934_2003 == 460999
replace gdenr_2012 =6139  if UCID_1934_2003 == 460999
replace gdenr_2018 = 6139 if UCID_1934_2003 == 460999
	
replace gdename ="raron" if UCID_1934_2003 == 462746
replace gdenr = 6199 if UCID_1934_2003 == 462746
replace gdenr_2012 = 6199 if UCID_1934_2003 == 462746
replace gdenr_2018 = 6199 if UCID_1934_2003 == 462746
	
replace gdename ="unterbaech" if UCID_1934_2003 == 462793
replace gdenr = 6201 if UCID_1934_2003 == 462793
replace gdenr_2012 = 6201 if UCID_1934_2003 == 462793
replace gdenr_2018 = 6201 if UCID_1934_2003 == 462793

replace gdename ="saint maurice" if UCID_1934_2003 == 463024
replace gdenr = 6217 if UCID_1934_2003 == 463024
replace gdenr_2012 = 6217 if UCID_1934_2003 == 463024
replace gdenr_2018 = 6217 if UCID_1934_2003 == 463024

replace gdename ="saint maurice" if UCID_1934_2003 == 463010
replace gdenr = 6217 if UCID_1934_2003 == 463010
replace gdenr_2012 = 6217 if UCID_1934_2003 == 463010
replace gdenr_2018 = 6217 if UCID_1934_2003 == 463010

replace gdename ="chalais" if UCID_1934_2003 == 463211
replace gdenr = 6232 if UCID_1934_2003 == 463211
replace gdenr_2012 = 6232 if UCID_1934_2003 == 463211
replace gdenr_2018 = 6232 if UCID_1934_2003 == 463211

replace gdename ="lens" if UCID_1934_2003 == 466099
replace gdenr = 6240 if UCID_1934_2003 == 466099
replace gdenr_2012 = 6240 if UCID_1934_2003 == 466099
replace gdenr_2018 = 6240 if UCID_1934_2003 == 466099

replace gdename ="lens" if UCID_1934_2003 == 465866
replace gdenr = 6240 if UCID_1934_2003 == 465866
replace gdenr_2012 = 6240 if UCID_1934_2003 == 465866
replace gdenr_2018 = 6240 if UCID_1934_2003 == 465866

replace gdename ="lens" if UCID_1934_2003 == 463622
replace gdenr = 6240 if UCID_1934_2003 == 463622
replace gdenr_2012 = 6240 if UCID_1934_2003 == 463622
replace gdenr_2018 = 6240 if UCID_1934_2003 == 463622

replace gdename ="lens" if UCID_1934_2003 == 463420
replace gdenr = 6240 if UCID_1934_2003 == 463420
replace gdenr_2012 = 6240 if UCID_1934_2003 == 463420
replace gdenr_2018 = 6240 if UCID_1934_2003 == 463420

replace gdename ="sierre" if UCID_1934_2003 == 464803
replace gdenr = 6248 if UCID_1934_2003 == 464803
replace gdenr_2012 = 6248 if UCID_1934_2003 == 464803
replace gdenr_2018 = 6248 if UCID_1934_2003 == 464803

replace gdename ="sierre" if UCID_1934_2003 == 464171
replace gdenr = 6248 if UCID_1934_2003 == 464171
replace gdenr_2012 = 6248 if UCID_1934_2003 == 464171
replace gdenr_2018 = 6248 if UCID_1934_2003 == 464171

replace gdename ="montana" if UCID_1934_2003 == 1342578
replace gdenr = 6243 if UCID_1934_2003 == 1342578
replace gdenr_2012 = 6243 if UCID_1934_2003 == 1342578
replace gdenr_2018 = 6253 if UCID_1934_2003 == 1342578

replace gdename ="montana" if UCID_1934_2003 == 1256521
replace gdenr = 6243 if UCID_1934_2003 == 1256521
replace gdenr_2012 = 6243 if UCID_1934_2003 == 1256521
replace gdenr_2018 = 6253 if UCID_1934_2003 == 1256521

replace gdename ="montana" if UCID_1934_2003 == 1144282
replace gdenr = 6243 if UCID_1934_2003 == 1144282
replace gdenr_2012 = 6243 if UCID_1934_2003 == 1144282
replace gdenr_2018 = 6253 if UCID_1934_2003 == 1144282

replace gdename ="montana" if UCID_1934_2003 == 466320
replace gdenr = 6243 if UCID_1934_2003 == 466320
replace gdenr_2012 = 6243 if UCID_1934_2003 == 466320
replace gdenr_2018 = 6253 if UCID_1934_2003 == 466320

replace gdename ="montana" if UCID_1934_2003 == 465647
replace gdenr = 6243 if UCID_1934_2003 == 465647
replace gdenr_2012 = 6243 if UCID_1934_2003 == 465647
replace gdenr_2018 = 6253 if UCID_1934_2003 == 465647

replace gdename ="montana" if UCID_1934_2003 == 465312
replace gdenr = 6243 if UCID_1934_2003 == 465312
replace gdenr_2012 = 6243 if UCID_1934_2003 == 465312
replace gdenr_2018 = 6253 if UCID_1934_2003 == 465312

replace gdename ="mollens vs" if UCID_1934_2003 == 465320
replace gdenr = 6242 if UCID_1934_2003 == 465320
replace gdenr_2012 = 6242 if UCID_1934_2003 == 465320
replace gdenr_2018 = 6253 if UCID_1934_2003 == 465320

replace gdename ="mollens vs" if UCID_1934_2003 == 465252
replace gdenr = 6242 if UCID_1934_2003 == 465252
replace gdenr_2012 = 6242 if UCID_1934_2003 == 465252
replace gdenr_2018 = 6253 if UCID_1934_2003 == 465252

replace gdename ="mollens vs" if UCID_1934_2003 == 460979
replace gdenr = 6242 if UCID_1934_2003 == 460979
replace gdenr_2012 = 6242 if UCID_1934_2003 == 460979
replace gdenr_2018 = 6253 if UCID_1934_2003 == 460979
		
replace gdename ="mollens vs" if UCID_1934_2003 == 460973
replace gdenr = 6242 if UCID_1934_2003 == 460973
replace gdenr_2012 = 6242 if UCID_1934_2003 == 460973
replace gdenr_2018 = 6253 if UCID_1934_2003 == 460973

replace gdename ="chermignon" if UCID_1934_2003 == 1227386
replace gdenr = 6234 if UCID_1934_2003 == 1227386
replace gdenr_2012 = 6234 if UCID_1934_2003 == 1227386
replace gdenr_2018 = 6253 if UCID_1934_2003 == 1227386
	
replace gdename ="chermignon" if UCID_1934_2003 == 1325470
replace gdenr = 6234 if UCID_1934_2003 == 1325470
replace gdenr_2012 = 6234 if UCID_1934_2003 == 1325470
replace gdenr_2018 = 6253 if UCID_1934_2003 == 1325470

replace gdename ="sion" if UCID_1934_2003 == 1144336
replace gdenr = 6266 if UCID_1934_2003 == 1144336
replace gdenr_2012 = 6266 if UCID_1934_2003 == 1144336
replace gdenr_2018 = 6266 if UCID_1934_2003 == 1144336

replace gdename ="sion" if UCID_1934_2003 == 466865
replace gdenr = 6266 if UCID_1934_2003 == 466865
replace gdenr_2012 = 6266 if UCID_1934_2003 == 466865
replace gdenr_2018 = 6266 if UCID_1934_2003 == 466865

replace gdename ="sion" if UCID_1934_2003 == 465928
replace gdenr = 6266 if UCID_1934_2003 == 465928
replace gdenr_2012 = 6266 if UCID_1934_2003 == 465928
replace gdenr_2018 = 6266 if UCID_1934_2003 == 465928
	
replace gdename ="randa" if UCID_1934_2003 == 541882
replace gdenr = 6287 if UCID_1934_2003 == 541882
replace gdenr_2012 = 6287 if UCID_1934_2003 == 541882
replace gdenr_2018 = 6287 if UCID_1934_2003 == 541882

replace gdename ="bole" if UCID_1934_2003 == 471206
replace gdenr = 6403 if UCID_1934_2003 == 471206
replace gdenr_2012 = 6403 if UCID_1934_2003 == 471206
replace gdenr_2018 = 6416 if UCID_1934_2003 == 471206

replace gdename ="bole" if UCID_1934_2003 == 471150
replace gdenr = 6403 if UCID_1934_2003 == 471150
replace gdenr_2012 = 6403 if UCID_1934_2003 == 471150
replace gdenr_2018 = 6416 if UCID_1934_2003 == 471150

replace gdename ="bole" if UCID_1934_2003 == 471096
replace gdenr = 6403 if UCID_1934_2003 == 471096
replace gdenr_2012 = 6403 if UCID_1934_2003 == 471096
replace gdenr_2018 = 6416 if UCID_1934_2003 == 471096

replace gdename ="la chaux de fonds" if UCID_1934_2003 == 473073
replace gdenr = 6421 if UCID_1934_2003 == 473073
replace gdenr_2012 = 6421 if UCID_1934_2003 == 473073
replace gdenr_2018 = 6421 if UCID_1934_2003 == 473073

replace gdename ="le locle" if UCID_1934_2003 == 474594
replace gdenr = 6436 if UCID_1934_2003 == 474594
replace gdenr_2012 = 6436 if UCID_1934_2003 == 474594
replace gdenr_2018 = 6436 if UCID_1934_2003 == 474594

replace gdename ="neuchatel" if UCID_1934_2003 == 479279
replace gdenr = 6458 if UCID_1934_2003 == 479279
replace gdenr_2012 = 6458 if UCID_1934_2003 == 479279
replace gdenr_2018 = 6458 if UCID_1934_2003 == 479279
	
replace gdename ="marin epagnier" if UCID_1934_2003 == 492419
replace gdenr = 6457 if UCID_1934_2003 == 492419
replace gdenr_2012 = 6461 if UCID_1934_2003 == 492419
replace gdenr_2018 = 6461 if UCID_1934_2003 == 492419

replace gdename ="bellevue" if UCID_1934_2003 == 481135
replace gdenr = 6606 if UCID_1934_2003 == 481135
replace gdenr_2012 = 6606 if UCID_1934_2003 == 481135
replace gdenr_2018 = 6606 if UCID_1934_2003 == 481135

replace gdename ="bernex" if UCID_1934_2003 == 484780
replace gdenr = 6607 if UCID_1934_2003 == 484780
replace gdenr_2012 = 6607 if UCID_1934_2003 == 484780
replace gdenr_2018 = 6607 if UCID_1934_2003 == 484780

replace gdename ="bernex" if UCID_1934_2003 == 481379
replace gdenr = 6607 if UCID_1934_2003 == 481379
replace gdenr_2012 = 6607 if UCID_1934_2003 == 481379
replace gdenr_2018 = 6607 if UCID_1934_2003 == 481379

replace gdename ="bernex" if UCID_1934_2003 == 481361
replace gdenr = 6607 if UCID_1934_2003 == 481361
replace gdenr_2012 = 6607 if UCID_1934_2003 == 481361
replace gdenr_2018 = 6607 if UCID_1934_2003 == 481361

replace gdename ="carouge ge" if UCID_1934_2003 == 483555
replace gdenr = 6608 if UCID_1934_2003 == 483555
replace gdenr_2012 = 6608 if UCID_1934_2003 == 483555
replace gdenr_2018 = 6608 if UCID_1934_2003 == 483555

replace gdename ="carouge ge" if UCID_1934_2003 == 481802
replace gdenr = 6608 if UCID_1934_2003 == 481802
replace gdenr_2012 = 6608 if UCID_1934_2003 == 481802
replace gdenr_2018 = 6608 if UCID_1934_2003 == 481802

replace gdename ="carouge ge" if UCID_1934_2003 == 481065
replace gdenr = 6608 if UCID_1934_2003 == 481065
replace gdenr_2012 = 6608 if UCID_1934_2003 == 481065
replace gdenr_2018 = 6608 if UCID_1934_2003 == 481065

replace gdename ="geneve" if UCID_1934_2003 == 1353939
replace gdenr = 6621 if UCID_1934_2003 == 1353939
replace gdenr_2012 = 6621 if UCID_1934_2003 == 1353939
replace gdenr_2018 = 6621 if UCID_1934_2003 == 1353939

replace gdename ="geneve" if UCID_1934_2003 == 1344831
replace gdenr = 6621 if UCID_1934_2003 == 1344831
replace gdenr_2012 = 6621 if UCID_1934_2003 == 1344831
replace gdenr_2018 = 6621 if UCID_1934_2003 == 1344831

replace gdename ="geneve" if UCID_1934_2003 == 1333250
replace gdenr = 6621 if UCID_1934_2003 == 1333250
replace gdenr_2012 = 6621 if UCID_1934_2003 == 1333250
replace gdenr_2018 = 6621 if UCID_1934_2003 == 1333250

replace gdename ="geneve" if UCID_1934_2003 == 1327967
replace gdenr = 6621 if UCID_1934_2003 == 1327967
replace gdenr_2012 = 6621 if UCID_1934_2003 == 1327967
replace gdenr_2018 = 6621 if UCID_1934_2003 == 1327967

replace gdename ="geneve" if UCID_1934_2003 == 1257914
replace gdenr = 6621 if UCID_1934_2003 == 1257914
replace gdenr_2012 = 6621 if UCID_1934_2003 == 1257914
replace gdenr_2018 = 6621 if UCID_1934_2003 == 1257914

replace gdename ="geneve" if UCID_1934_2003 == 947726
replace gdenr = 6621 if UCID_1934_2003 == 947726
replace gdenr_2012 = 6621 if UCID_1934_2003 == 947726
replace gdenr_2018 = 6621 if UCID_1934_2003 == 947726

replace gdename ="geneve" if UCID_1934_2003 == 947224
replace gdenr = 6621 if UCID_1934_2003 == 947224
replace gdenr_2012 = 6621 if UCID_1934_2003 == 947224
replace gdenr_2018 = 6621 if UCID_1934_2003 == 947224

replace gdename ="geneve" if UCID_1934_2003 == 947120
replace gdenr = 6621 if UCID_1934_2003 == 947120
replace gdenr_2012 = 6621 if UCID_1934_2003 == 947120
replace gdenr_2018 = 6621 if UCID_1934_2003 == 947120

replace gdename ="geneve" if UCID_1934_2003 == 877990
replace gdenr = 6621 if UCID_1934_2003 == 877990
replace gdenr_2012 = 6621 if UCID_1934_2003 == 877990
replace gdenr_2018 = 6621 if UCID_1934_2003 == 877990

replace gdename ="geneve" if UCID_1934_2003 == 733094
replace gdenr = 6621 if UCID_1934_2003 == 733094
replace gdenr_2012 = 6621 if UCID_1934_2003 == 733094
replace gdenr_2018 = 6621 if UCID_1934_2003 == 733094

replace gdename ="geneve" if UCID_1934_2003 == 532297
replace gdenr = 6621 if UCID_1934_2003 == 532297
replace gdenr_2012 = 6621 if UCID_1934_2003 == 532297
replace gdenr_2018 = 6621 if UCID_1934_2003 == 532297

replace gdename ="geneve" if UCID_1934_2003 == 520495
replace gdenr = 6621 if UCID_1934_2003 == 520495
replace gdenr_2012 = 6621 if UCID_1934_2003 == 520495
replace gdenr_2018 = 6621 if UCID_1934_2003 == 520495

replace gdename ="geneve" if UCID_1934_2003 == 518420
replace gdenr = 6621 if UCID_1934_2003 == 518420
replace gdenr_2012 = 6621 if UCID_1934_2003 == 518420
replace gdenr_2018 = 6621 if UCID_1934_2003 == 518420

replace gdename ="geneve" if UCID_1934_2003 == 517577
replace gdenr = 6621 if UCID_1934_2003 == 517577
replace gdenr_2012 = 6621 if UCID_1934_2003 == 517577
replace gdenr_2018 = 6621 if UCID_1934_2003 == 517577

replace gdename ="geneve" if UCID_1934_2003 == 502252
replace gdenr = 6621 if UCID_1934_2003 == 502252
replace gdenr_2012 = 6621 if UCID_1934_2003 == 502252
replace gdenr_2018 = 6621 if UCID_1934_2003 == 502252

replace gdename ="geneve" if UCID_1934_2003 == 501372
replace gdenr = 6621 if UCID_1934_2003 == 501372
replace gdenr_2012 = 6621 if UCID_1934_2003 == 501372
replace gdenr_2018 = 6621 if UCID_1934_2003 == 501372

replace gdename ="geneve" if UCID_1934_2003 == 501085
replace gdenr = 6621 if UCID_1934_2003 == 501085
replace gdenr_2012 = 6621 if UCID_1934_2003 == 501085
replace gdenr_2018 = 6621 if UCID_1934_2003 == 501085

replace gdename ="geneve" if UCID_1934_2003 == 497774
replace gdenr = 6621 if UCID_1934_2003 == 497774
replace gdenr_2012 = 6621 if UCID_1934_2003 == 497774
replace gdenr_2018 = 6621 if UCID_1934_2003 == 497774

replace gdename ="geneve" if UCID_1934_2003 == 496101
replace gdenr = 6621 if UCID_1934_2003 == 496101
replace gdenr_2012 = 6621 if UCID_1934_2003 == 496101
replace gdenr_2018 = 6621 if UCID_1934_2003 == 496101

replace gdename ="geneve" if UCID_1934_2003 == 495081
replace gdenr = 6621 if UCID_1934_2003 == 495081
replace gdenr_2012 = 6621 if UCID_1934_2003 == 495081
replace gdenr_2018 = 6621 if UCID_1934_2003 == 495081

replace gdename ="geneve" if UCID_1934_2003 == 491989
replace gdenr = 6621 if UCID_1934_2003 == 491989
replace gdenr_2012 = 6621 if UCID_1934_2003 == 491989
replace gdenr_2018 = 6621 if UCID_1934_2003 == 491989

replace gdename ="geneve" if UCID_1934_2003 == 489729
replace gdenr = 6621 if UCID_1934_2003 == 489729
replace gdenr_2012 = 6621 if UCID_1934_2003 == 489729
replace gdenr_2018 = 6621 if UCID_1934_2003 == 489729

replace gdename ="geneve" if UCID_1934_2003 == 489302
replace gdenr = 6621 if UCID_1934_2003 == 489302
replace gdenr_2012 = 6621 if UCID_1934_2003 == 489302
replace gdenr_2018 = 6621 if UCID_1934_2003 == 489302

replace gdename ="geneve" if UCID_1934_2003 == 486923
replace gdenr = 6621 if UCID_1934_2003 == 486923
replace gdenr_2012 = 6621 if UCID_1934_2003 == 486923
replace gdenr_2018 = 6621 if UCID_1934_2003 == 486923

replace gdename ="lancy" if UCID_1934_2003 == 494031
replace gdenr = 6628 if UCID_1934_2003 == 494031
replace gdenr_2012 = 6628 if UCID_1934_2003 == 494031
replace gdenr_2018 = 6628 if UCID_1934_2003 == 494031
	
replace gdename ="meyrin" if UCID_1934_2003 == 530244
replace gdenr = 6630 if UCID_1934_2003 == 530244
replace gdenr_2012 = 6630 if UCID_1934_2003 == 530244
replace gdenr_2018 = 6630 if UCID_1934_2003 == 530244

replace gdename ="onex" if UCID_1934_2003 == 530599
replace gdenr = 6631 if UCID_1934_2003 == 530599
replace gdenr_2012 = 6631 if UCID_1934_2003 == 530599
replace gdenr_2018 = 6631 if UCID_1934_2003 == 530599

replace gdename ="thonex" if UCID_1934_2003 == 532488
replace gdenr = 6640 if UCID_1934_2003 == 532488
replace gdenr_2012 = 6640 if UCID_1934_2003 == 532488
replace gdenr_2018 = 6640 if UCID_1934_2003 == 532488
	
replace gdename ="thonex" if UCID_1934_2003 == 532247
replace gdenr = 6640 if UCID_1934_2003 == 532247
replace gdenr_2012 = 6640 if UCID_1934_2003 == 532247
replace gdenr_2018 = 6640 if UCID_1934_2003 == 532247

replace gdename ="veyrier" if UCID_1934_2003 == 534191
replace gdenr = 6645 if UCID_1934_2003 == 534191
replace gdenr_2012 = 6645 if UCID_1934_2003 == 534191
replace gdenr_2018 = 6645 if UCID_1934_2003 == 534191
	
replace gdename ="veyrier" if UCID_1934_2003 == 503614
replace gdenr = 6645 if UCID_1934_2003 == 503614
replace gdenr_2012 = 6645 if UCID_1934_2003 == 503614
replace gdenr_2018 = 6645 if UCID_1934_2003 == 503614

replace gdename ="boecourt" if UCID_1934_2003 == 535721
replace gdenr = 6702 if UCID_1934_2003 == 535721
replace gdenr_2012 = 6702 if UCID_1934_2003 == 535721
replace gdenr_2018 = 6702 if UCID_1934_2003 == 535721

replace gdename ="delemont" if UCID_1934_2003 == 1360856
replace gdenr = 467 if UCID_1934_2003 == 1360856 & year < 1979
replace gdenr = 6711 if UCID_1934_2003 == 1360856 & year > 1978
replace gdenr_2012 = 6711 if UCID_1934_2003 == 1360856
replace gdenr_2018 = 6711 if UCID_1934_2003 == 1360856

replace gdename ="vicques" if UCID_1934_2003 == 469817
replace gdenr = 6727 if UCID_1934_2003 == 469817
replace gdenr_2012 = 6727 if UCID_1934_2003 == 469817
replace gdenr_2018 = 6730 if UCID_1934_2003 == 469817

replace gdename ="vicques" if UCID_1934_2003 == 469595
replace gdenr = 6727 if UCID_1934_2003 == 469595
replace gdenr_2012 = 6727 if UCID_1934_2003 == 469595
replace gdenr_2018 = 6730 if UCID_1934_2003 == 469595

replace gdename ="vicques" if UCID_1934_2003 == 536043
replace gdenr = 6727 if UCID_1934_2003 == 536043
replace gdenr_2012 = 6727 if UCID_1934_2003 == 536043
replace gdenr_2018 = 6730 if UCID_1934_2003 == 536043

replace gdename ="vermes" if UCID_1934_2003 == 533064
replace gdenr = 482 if UCID_1934_2003 == 1360856 & year < 1979
replace gdenr = 6726 if UCID_1934_2003 == 1360856 & year > 1978
replace gdenr_2012 = 6726 if UCID_1934_2003 == 533064
replace gdenr_2018 = 6730 if UCID_1934_2003 == 533064

replace gdename ="porrentruy" if UCID_1934_2003 == 537153
replace gdenr = 830 if UCID_1934_2003 == 537153 & year < 1979
replace gdenr = 6800 if UCID_1934_2003 == 1360856 & year > 1978
replace gdenr_2012 = 6800 if UCID_1934_2003 == 537153
replace gdenr_2018 = 6800 if UCID_1934_2003 == 537153

replace gdename ="buix" if UCID_1934_2003 == 536117
replace gdenr = 807 if UCID_1934_2003 == 536117 & year < 1979
replace gdenr = 6777 if UCID_1934_2003 == 1360856 & year > 1978
replace gdenr_2012 = 6807 if UCID_1934_2003 == 536117
replace gdenr_2018 = 6807 if UCID_1934_2003 == 536117

replace gdename ="pully" if UCID_1934_2003 == 436160
replace gdenr = 5590 if UCID_1934_2003 == 436160
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436160
replace gdenr_2018 = 5590 if UCID_1934_2003 == 436160
			
replace gdename ="pully" if UCID_1934_2003 == 436032
replace gdenr = 5590 if UCID_1934_2003 == 436032
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436032
replace gdenr_2018 = 5590 if UCID_1934_2003 == 436032
				
replace gdename ="pully" if UCID_1934_2003 == 436022
replace gdenr = 5590 if UCID_1934_2003 == 436022
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436022
replace gdenr_2018 = 5590 if UCID_1934_2003 == 436022
				
replace gdename ="pully" if UCID_1934_2003 == 436019
replace gdenr = 5590 if UCID_1934_2003 == 436019
replace gdenr_2012 = 5590 if UCID_1934_2003 == 436019
replace gdenr_2018 = 5590 if UCID_1934_2003 == 436019
				
replace gdename ="pully" if UCID_1934_2003 == 435990
replace gdenr = 5590 if UCID_1934_2003 == 435990
replace gdenr_2012 = 5590 if UCID_1934_2003 == 435990
replace gdenr_2018 = 5590 if UCID_1934_2003 == 435990
				
replace gdename ="pully" if UCID_1934_2003 == 435980
replace gdenr = 5590 if UCID_1934_2003 == 435980
replace gdenr_2012 = 5590 if UCID_1934_2003 == 435980
replace gdenr_2018 = 5590 if UCID_1934_2003 == 435980
			
replace gdename ="pully" if UCID_1934_2003 == 435898
replace gdenr = 5590 if UCID_1934_2003 == 435898
replace gdenr_2012 = 5590 if UCID_1934_2003 == 435898
replace gdenr_2018 = 5590 if UCID_1934_2003 == 435898
				
replace gdename ="pully" if UCID_1934_2003 == 435832
replace gdenr = 5590 if UCID_1934_2003 == 435832
replace gdenr_2012 = 5590 if UCID_1934_2003 == 435832
replace gdenr_2018 = 5590 if UCID_1934_2003 == 435832

replace gdename ="prilly" if UCID_1934_2003 == 435675
replace gdenr = 5589 if UCID_1934_2003 == 435675
replace gdenr_2012 = 5589 if UCID_1934_2003 == 435675
replace gdenr_2018 = 5589 if UCID_1934_2003 == 435675

replace gdename ="prilly" if UCID_1934_2003 == 435665
replace gdenr = 5589 if UCID_1934_2003 == 435665
replace gdenr_2012 = 5589 if UCID_1934_2003 == 435665
replace gdenr_2018 = 5589 if UCID_1934_2003 == 435665
							
replace gdename ="prilly" if UCID_1934_2003 == 435635
replace gdenr = 5589 if UCID_1934_2003 == 435635
replace gdenr_2012 = 5589 if UCID_1934_2003 == 435635
replace gdenr_2018 = 5589 if UCID_1934_2003 == 435635	
			
replace gdename ="prilly" if UCID_1934_2003 == 435622
replace gdenr = 5589 if UCID_1934_2003 == 435622
replace gdenr_2012 = 5589 if UCID_1934_2003 == 435622
replace gdenr_2018 = 5589 if UCID_1934_2003 ==435622 	
		
replace gdename ="prilly" if UCID_1934_2003 == 435619
replace gdenr = 5589 if UCID_1934_2003 == 435619
replace gdenr_2012 = 5589 if UCID_1934_2003 == 435619
replace gdenr_2018 = 5589 if UCID_1934_2003 == 435619	
			
replace gdename ="prilly" if UCID_1934_2003 == 435609
replace gdenr = 5589 if UCID_1934_2003 == 435609
replace gdenr_2012 = 5589 if UCID_1934_2003 == 435609
replace gdenr_2018 = 5589 if UCID_1934_2003 == 435609	
			
replace gdename ="prilly" if UCID_1934_2003 == 435592
replace gdenr = 5589 if UCID_1934_2003 == 435592
replace gdenr_2012 = 5589 if UCID_1934_2003 == 435592
replace gdenr_2018 = 5589 if UCID_1934_2003 == 435592
			
replace gdename ="prilly" if UCID_1934_2003 == 435528
replace gdenr = 5589 if UCID_1934_2003 == 435528
replace gdenr_2012 = 5589 if UCID_1934_2003 == 435528
replace gdenr_2018 = 5589 if UCID_1934_2003 == 435528
			
replace gdename ="prilly" if UCID_1934_2003 == 435498
replace gdenr = 5589 if UCID_1934_2003 == 435498
replace gdenr_2012 = 5589 if UCID_1934_2003 == 435498
replace gdenr_2018 = 5589 if UCID_1934_2003 == 435498

replace gdename ="prilly" if UCID_1934_2003 == 426422
replace gdenr = 5589 if UCID_1934_2003 == 426422
replace gdenr_2012 = 5589 if UCID_1934_2003 == 426422
replace gdenr_2018 = 5589 if UCID_1934_2003 == 426422

replace gdename ="prilly" if UCID_1934_2003 == 435470
replace gdenr = 5589 if UCID_1934_2003 == 435470
replace gdenr_2012 = 5589 if UCID_1934_2003 == 435470
replace gdenr_2018 = 5589 if UCID_1934_2003 == 435470

replace gdename ="prilly" if UCID_1934_2003 == 423919
replace gdenr = 5589 if UCID_1934_2003 == 423919
replace gdenr_2012 = 5589 if UCID_1934_2003 == 423919
replace gdenr_2018 = 5589 if UCID_1934_2003 == 423919				
		
replace gdename ="prilly" if UCID_1934_2003 == 423704
replace gdenr = 5589 if UCID_1934_2003 == 423704
replace gdenr_2012 = 5589 if UCID_1934_2003 == 423704
replace gdenr_2018 = 5589 if UCID_1934_2003 == 423704	

replace gdename ="pully" if UCID_1934_2003 == 434539
replace gdenr = 5590 if UCID_1934_2003 == 434539
replace gdenr_2012 = 5590 if UCID_1934_2003 == 434539
replace gdenr_2018 = 5590 if UCID_1934_2003 == 434539
			
replace gdename ="lausanne" if UCID_1934_2003 == 432973
replace gdenr = 5586 if UCID_1934_2003 == 432973
replace gdenr_2012 = 5586 if UCID_1934_2003 == 432973
replace gdenr_2018 = 5586 if UCID_1934_2003 == 432973
				
replace gdename ="pully" if UCID_1934_2003 == 432032
replace gdenr = 5590 if UCID_1934_2003 == 432032
replace gdenr_2012 = 5590 if UCID_1934_2003 == 432032
replace gdenr_2018 = 5590 if UCID_1934_2003 == 432032
				
replace gdename ="pully" if UCID_1934_2003 == 431243
replace gdenr = 5590 if UCID_1934_2003 == 431243
replace gdenr_2012 = 5590 if UCID_1934_2003 == 431243
replace gdenr_2018 = 5590 if UCID_1934_2003 == 431243
			
replace gdename ="pully" if UCID_1934_2003 == 431160
replace gdenr = 5590 if UCID_1934_2003 == 431160
replace gdenr_2012 = 5590 if UCID_1934_2003 ==431160 
replace gdenr_2018 = 5590 if UCID_1934_2003 == 431160

replace gdename ="pully" if UCID_1934_2003 == 430751
replace gdenr = 5590 if UCID_1934_2003 == 430751
replace gdenr_2012 = 5590 if UCID_1934_2003 == 430751
replace gdenr_2018 = 5590 if UCID_1934_2003 == 430751
						
replace gdename ="pully" if UCID_1934_2003 == 428295
replace gdenr = 5590 if UCID_1934_2003 == 428295
replace gdenr_2012 = 5590 if UCID_1934_2003 == 428295
replace gdenr_2018 = 5590 if UCID_1934_2003 == 428295
	
replace gdename ="pully" if UCID_1934_2003 == 427694
replace gdenr = 5590 if UCID_1934_2003 == 427694
replace gdenr_2012 = 5590 if UCID_1934_2003 == 427694
replace gdenr_2018 = 5590 if UCID_1934_2003 == 427694
				
replace gdename ="pully" if UCID_1934_2003 == 425873
replace gdenr = 5590 if UCID_1934_2003 == 425873
replace gdenr_2012 = 5590 if UCID_1934_2003 == 425873
replace gdenr_2018 = 5590 if UCID_1934_2003 == 425873

replace gdename ="pully" if UCID_1934_2003 == 425758
replace gdenr = 5590 if UCID_1934_2003 == 425758
replace gdenr_2012 = 5590 if UCID_1934_2003 == 425758
replace gdenr_2018 = 5590 if UCID_1934_2003 == 425758
			
replace gdename ="pully" if UCID_1934_2003 == 419847
replace gdenr = 5590 if UCID_1934_2003 == 419847
replace gdenr_2012 = 5590 if UCID_1934_2003 == 419847
replace gdenr_2018 = 5590 if UCID_1934_2003 == 419847
		
replace gdename ="pully" if UCID_1934_2003 == 419739
replace gdenr = 5590 if UCID_1934_2003 == 419739
replace gdenr_2012 = 5590 if UCID_1934_2003 == 419739
replace gdenr_2018 = 5590 if UCID_1934_2003 == 419739
		
replace gdename ="pully" if UCID_1934_2003 == 421305
replace gdenr = 5590 if UCID_1934_2003 == 421305
replace gdenr_2012 = 5590 if UCID_1934_2003 == 421305
replace gdenr_2018 = 5590 if UCID_1934_2003 == 421305

replace gdename ="renens (vd)" if UCID_1934_2003 == 429535
replace gdenr = 5591 if UCID_1934_2003 == 429535
replace gdenr_2012 = 5591 if UCID_1934_2003 == 429535
replace gdenr_2018 = 5591 if UCID_1934_2003 == 429535
			
replace gdename ="lausanne" if UCID_1934_2003 == 429215
replace gdenr = 5586 if UCID_1934_2003 == 429215
replace gdenr_2012 = 5586 if UCID_1934_2003 == 429215
replace gdenr_2018 = 5586 if UCID_1934_2003 == 429215

replace gdename ="lausanne" if UCID_1934_2003 == 427916
replace gdenr = 5586 if UCID_1934_2003 == 427916
replace gdenr_2012 = 5586 if UCID_1934_2003 == 427916
replace gdenr_2018 = 5586 if UCID_1934_2003 == 427916

replace gdename ="lausanne" if UCID_1934_2003 == 426776
replace gdenr = 5586 if UCID_1934_2003 == 426776
replace gdenr_2012 = 5586 if UCID_1934_2003 == 426776
replace gdenr_2018 = 5586 if UCID_1934_2003 == 426776
			
replace gdename ="lausanne" if UCID_1934_2003 == 421187
replace gdenr = 5586 if UCID_1934_2003 == 421187
replace gdenr_2012 = 5586 if UCID_1934_2003 == 421187
replace gdenr_2018 = 5586 if UCID_1934_2003 == 421187

replace gdename ="lausanne" if UCID_1934_2003 == 421024
replace gdenr = 5586 if UCID_1934_2003 == 421024
replace gdenr_2012 = 5586 if UCID_1934_2003 == 421024
replace gdenr_2018 = 5586 if UCID_1934_2003 == 421024
								
replace gdename ="collombey muraz" if UCID_1934_2003 == 429067
replace gdenr = 6152 if UCID_1934_2003 == 429067
replace gdenr_2012 = 6152 if UCID_1934_2003 == 429067
replace gdenr_2018 = 6152 if UCID_1934_2003 == 429067		
				
replace gdename ="geneve" if UCID_1934_2003 == 426016
replace gdenr = 6621 if UCID_1934_2003 == 426016
replace gdenr_2012 = 6621 if UCID_1934_2003 == 426016
replace gdenr_2018 = 6621 if UCID_1934_2003 == 426016
		
replace gdename ="geneve" if UCID_1934_2003 == 386005
replace gdenr = 6621 if UCID_1934_2003 == 386005
replace gdenr_2012 = 6621 if UCID_1934_2003 == 386005
replace gdenr_2018 = 6621 if UCID_1934_2003 == 386005
	
replace gdename ="geneve" if UCID_1934_2003 == 424045
replace gdenr = 6621 if UCID_1934_2003 == 424045
replace gdenr_2012 = 6621 if UCID_1934_2003 == 424045
replace gdenr_2018 = 6621 if UCID_1934_2003 == 424045 	
				
replace gdename ="crassier" if UCID_1934_2003 == 424940
replace gdenr = 5714 if UCID_1934_2003 == 424940
replace gdenr_2012 = 5714 if UCID_1934_2003 == 424940
replace gdenr_2018 = 5714 if UCID_1934_2003 == 424940
								
replace gdename ="cham" if UCID_1934_2003 == 424358
replace gdenr = 1702 if UCID_1934_2003 == 424358
replace gdenr_2012 = 1702 if UCID_1934_2003 == 424358
replace gdenr_2018 = 1702 if UCID_1934_2003 == 424358		
					
replace gdename ="paudex" if UCID_1934_2003 == 418660
replace gdenr = 5588 if UCID_1934_2003 == 418660
replace gdenr_2012 = 5588 if UCID_1934_2003 == 418660
replace gdenr_2018 = 5588 if UCID_1934_2003 == 418660
	
replace gdename ="crissier" if UCID_1934_2003 == 417821
replace gdenr = 5583 if UCID_1934_2003 == 417821
replace gdenr_2012 = 5583 if UCID_1934_2003 == 417821
replace gdenr_2018 = 5583 if UCID_1934_2003 == 417821
	
replace gdename ="crissier" if UCID_1934_2003 == 417746
replace gdenr = 5583 if UCID_1934_2003 == 417746
replace gdenr_2012 = 5583 if UCID_1934_2003 == 417746
replace gdenr_2018 = 5583 if UCID_1934_2003 == 417746 	
			
replace gdename ="yvorne" if UCID_1934_2003 == 415028
replace gdenr = 5415 if UCID_1934_2003 == 415028
replace gdenr_2012 = 5415 if UCID_1934_2003 ==415028 
replace gdenr_2018 = 5415 if UCID_1934_2003 ==415028 	
			
replace gdename ="montreux" if UCID_1934_2003 == 414755
replace gdenr = 5886 if UCID_1934_2003 == 414755
replace gdenr_2012 = 5886 if UCID_1934_2003 == 414755
replace gdenr_2018 = 5886 if UCID_1934_2003 == 414755
		
replace gdename ="ormont dessus" if UCID_1934_2003 == 414646
replace gdenr = 5411 if UCID_1934_2003 == 414646
replace gdenr_2012 = 5411 if UCID_1934_2003 == 414646
replace gdenr_2018 = 5411 if UCID_1934_2003 == 414646
					
replace gdename ="ollon" if UCID_1934_2003 == 414561
replace gdenr = 5409 if UCID_1934_2003 == 414561
replace gdenr_2012 = 5409 if UCID_1934_2003 == 414561
replace gdenr_2018 = 5409 if UCID_1934_2003 == 414561
	
replace gdename ="aigle" if UCID_1934_2003 == 413125
replace gdenr = 5401 if UCID_1934_2003 == 413125
replace gdenr_2012 = 5401 if UCID_1934_2003 == 413125
replace gdenr_2018 = 5401 if UCID_1934_2003 == 413125

replace gdename ="mendrisio" if UCID_1934_2003 == 410254
replace gdenr = 5254 if UCID_1934_2003 == 410254
replace gdenr_2012 = 5254 if UCID_1934_2003 == 410254
replace gdenr_2018 = 5254 if UCID_1934_2003 == 410254
					
replace gdename ="mendrisio" if UCID_1934_2003 == 410224
replace gdenr = 5254 if UCID_1934_2003 == 410224
replace gdenr_2012 = 5254 if UCID_1934_2003 == 410224
replace gdenr_2018 = 5254 if UCID_1934_2003 == 410224
					
replace gdename ="coldrerio" if UCID_1934_2003 == 409135
replace gdenr = 5251 if UCID_1934_2003 == 409135
replace gdenr_2012 = 5251 if UCID_1934_2003 == 409135
replace gdenr_2018 = 5251 if UCID_1934_2003 == 409135

replace gdename ="chiasso" if UCID_1934_2003 == 404817
replace gdenr = 5250 if UCID_1934_2003 == 404817
replace gdenr_2012 = 5250 if UCID_1934_2003 == 404817
replace gdenr_2018 = 5250 if UCID_1934_2003 == 404817
	
replace gdename ="chiasso" if UCID_1934_2003 == 406050
replace gdenr = 5250 if UCID_1934_2003 == 406050
replace gdenr_2012 = 5250 if UCID_1934_2003 == 406050
replace gdenr_2018 = 5250 if UCID_1934_2003 == 406050
		
replace gdename ="chur" if UCID_1934_2003 == 405233
replace gdenr = 3901 if UCID_1934_2003 == 405233
replace gdenr_2012 = 3901 if UCID_1934_2003 == 405233
replace gdenr_2018 = 3901 if UCID_1934_2003 == 405233
			
replace gdename ="castel san pietro" if UCID_1934_2003 == 403226
replace gdenr = 5249 if UCID_1934_2003 == 403226
replace gdenr_2012 = 5249 if UCID_1934_2003 == 403226
replace gdenr_2018 = 5249 if UCID_1934_2003 == 403226
				
replace gdename ="vezia" if UCID_1934_2003 == 402621
replace gdenr = 5231 if UCID_1934_2003 == 402621
replace gdenr_2012 = 5231 if UCID_1934_2003 == 402621
replace gdenr_2018 = 5231 if UCID_1934_2003 == 402621
								
replace gdename ="vezia" if UCID_1934_2003 == 401906
replace gdenr = 5231 if UCID_1934_2003 == 401906
replace gdenr_2012 = 5231 if UCID_1934_2003 == 401906
replace gdenr_2018 = 5231 if UCID_1934_2003 == 401906
				
replace gdename ="morcote" if UCID_1934_2003 == 399600
replace gdenr = 5203 if UCID_1934_2003 == 399600
replace gdenr_2012 = 5203  if UCID_1934_2003 == 399600
replace gdenr_2018 = 5203 if UCID_1934_2003 == 399600
				
replace gdename ="mesocco" if UCID_1934_2003 == 399560
replace gdenr = 3822 if UCID_1934_2003 == 399560
replace gdenr_2012 = 3822 if UCID_1934_2003 == 399560
replace gdenr_2018 = 3822 if UCID_1934_2003 == 399560
					
replace gdename ="mezzovico vira" if UCID_1934_2003 == 399398
replace gdenr = 5199 if UCID_1934_2003 == 399398
replace gdenr_2012 = 5199 if UCID_1934_2003 == 399398
replace gdenr_2018 = 5199 if UCID_1934_2003 == 399398
					
replace gdename ="bioggio" if UCID_1934_2003 == 395845
replace gdenr = 5151 if UCID_1934_2003 == 395845
replace gdenr_2012 = 5151 if UCID_1934_2003 == 395845
replace gdenr_2018 = 5151 if UCID_1934_2003 == 395845
					
replace gdename ="gandria" if UCID_1934_2003 == 392388
replace gdenr = 5184 if UCID_1934_2003 == 392388
replace gdenr_2012 = 5192 if UCID_1934_2003 == 392388
replace gdenr_2018 = 5192 if UCID_1934_2003 == 392388

replace gdename ="lugano" if UCID_1934_2003 == 392274
replace gdenr = 5192 if UCID_1934_2003 == 392274
replace gdenr_2012 = 5192 if UCID_1934_2003 == 392274
replace gdenr_2018 = 5192 if UCID_1934_2003 == 392274
			
replace gdename ="lugano" if UCID_1934_2003 == 382756
replace gdenr = 5192 if UCID_1934_2003 == 382756
replace gdenr_2012 = 5192 if UCID_1934_2003 == 382756
replace gdenr_2018 = 5192 if UCID_1934_2003 == 382756
			
replace gdename ="lugano" if UCID_1934_2003 == 381902
replace gdenr = 5192 if UCID_1934_2003 == 381902
replace gdenr_2012 = 5192 if UCID_1934_2003 == 381902
replace gdenr_2018 = 5192 if UCID_1934_2003 == 381902
				
replace gdename ="glarus" if UCID_1934_2003 == 384979
replace gdenr = 1609 if UCID_1934_2003 == 384979
replace gdenr_2012 = 1632 if UCID_1934_2003 == 384979
replace gdenr_2018 = 1632 if UCID_1934_2003 == 384979
							
replace gdename ="bioggio" if UCID_1934_2003 == 381106
replace gdenr = 5151 if UCID_1934_2003 == 381106
replace gdenr_2012 = 5151 if UCID_1934_2003 == 381106
replace gdenr_2018 = 5151 if UCID_1934_2003 == 381106
				
replace gdename ="vernate" if UCID_1934_2003 == 376349
replace gdenr = 5230 if UCID_1934_2003 == 376349
replace gdenr_2012 =5230  if UCID_1934_2003 == 376349
replace gdenr_2018 = 5230 if UCID_1934_2003 == 376349
				
replace gdename ="castagnola" if UCID_1934_2003 == 376348
replace gdenr = 5172 if UCID_1934_2003 == 376348
replace gdenr_2012 = 5192 if UCID_1934_2003 == 376348
replace gdenr_2018 =5192  if UCID_1934_2003 == 376348
					
replace gdename ="muzzano" if UCID_1934_2003 == 373608
replace gdenr = 5205 if UCID_1934_2003 == 373608
replace gdenr_2012 = 5205 if UCID_1934_2003 == 373608
replace gdenr_2018 =5205  if UCID_1934_2003 == 373608
				
replace gdename ="minusio" if UCID_1934_2003 == 372489
replace gdenr = 5118 if UCID_1934_2003 == 372489
replace gdenr_2012 = 5118 if UCID_1934_2003 == 372489
replace gdenr_2018 = 5118 if UCID_1934_2003 == 372489
				
replace gdename ="olivone" if UCID_1934_2003 == 371737
replace gdenr = 5043 if UCID_1934_2003 == 371737
replace gdenr_2012 = 5049 if UCID_1934_2003 == 371737
replace gdenr_2018 = 5049 if UCID_1934_2003 == 371737
					
replace gdename ="sant antonino" if UCID_1934_2003 == 369176
replace gdenr = 5017 if UCID_1934_2003 == 369176
replace gdenr_2012 = 5017 if UCID_1934_2003 == 369176
replace gdenr_2018 = 5017 if UCID_1934_2003 == 369176
					
replace gdename ="bellinzona" if UCID_1934_2003 == 368716
replace gdenr = 5002 if UCID_1934_2003 == 368716
replace gdenr_2012 = 5002 if UCID_1934_2003 == 368716
replace gdenr_2018 = 5002 if UCID_1934_2003 == 368716
		
replace gdename ="eschert" if UCID_1934_2003 == 364225
replace gdenr = 692 if UCID_1934_2003 == 364225
replace gdenr_2012 = 692 if UCID_1934_2003 == 364225
replace gdenr_2018 = 692 if UCID_1934_2003 == 364225
	
replace gdename ="wilen bei wil" if UCID_1934_2003 == 364130
replace gdenr = 4752 if UCID_1934_2003 == 364130
replace gdenr_2012 = 4786 if UCID_1934_2003 == 364130
replace gdenr_2018 = 4786 if UCID_1934_2003 == 364130
			
replace gdename ="wilen bei wil" if UCID_1934_2003 == 364090
replace gdenr = 4752 if UCID_1934_2003 == 364090
replace gdenr_2012 = 4786 if UCID_1934_2003 == 364090
replace gdenr_2018 = 4786 if UCID_1934_2003 == 364090
		
replace gdename ="muenchwilen tg" if UCID_1934_2003 == 363249
replace gdenr = 4746 if UCID_1934_2003 == 363249
replace gdenr_2012 = 4746 if UCID_1934_2003 == 363249
replace gdenr_2018 = 4746 if UCID_1934_2003 == 363249
			
replace gdename ="chur" if UCID_1934_2003 == 362792
replace gdenr = 3901 if UCID_1934_2003 == 362792
replace gdenr_2012 = 3901 if UCID_1934_2003 == 362792
replace gdenr_2018 = 3901 if UCID_1934_2003 == 362792

replace gdename ="frauenfeld" if UCID_1934_2003 == 359846
replace gdenr = 4566 if UCID_1934_2003 == 359846
replace gdenr_2012 = 4566 if UCID_1934_2003 == 359846
replace gdenr_2018 = 4566 if UCID_1934_2003 == 359846
		
replace gdename ="frauenfeld" if UCID_1934_2003 == 358846
replace gdenr = 4566 if UCID_1934_2003 == 358846
replace gdenr_2012 = 4566 if UCID_1934_2003 == 358846
replace gdenr_2018 = 4566 if UCID_1934_2003 == 358846

replace gdename ="frauenfeld" if UCID_1934_2003 == 362077
replace gdenr = 4566 if UCID_1934_2003 == 362077
replace gdenr_2012 = 4566 if UCID_1934_2003 == 362077
replace gdenr_2018 = 4566 if UCID_1934_2003 == 362077
					
replace gdename ="kreuzlingen" if UCID_1934_2003 == 361782
replace gdenr = 4671 if UCID_1934_2003 == 361782
replace gdenr_2012 = 4671 if UCID_1934_2003 == 361782
replace gdenr_2018 = 4671 if UCID_1934_2003 == 361782
		
replace gdename ="guntershausen bei aadorf" if UCID_1934_2003 == 358393
replace gdenr = 4554 if UCID_1934_2003 == 358393
replace gdenr_2012 = 4551 if UCID_1934_2003 == 358393
replace gdenr_2018 = 4551 if UCID_1934_2003 == 358393
					
replace gdename ="diessenhofen" if UCID_1934_2003 == 358331
replace gdenr = 4541 if UCID_1934_2003 == 358331
replace gdenr_2012 = 4545 if UCID_1934_2003 == 358331
replace gdenr_2018 = 4545 if UCID_1934_2003 == 358331

replace gdename ="doettingen" if UCID_1934_2003 == 354113
replace gdenr = 4304 if UCID_1934_2003 == 354113
replace gdenr_2012 = 4304 if UCID_1934_2003 == 354113
replace gdenr_2018 = 4304 if UCID_1934_2003 == 354113
			
replace gdename ="koelliken" if UCID_1934_2003 == 351777
replace gdenr = 4276 if UCID_1934_2003 == 351777
replace gdenr_2012 = 4276 if UCID_1934_2003 == 351777
replace gdenr_2018 = 4276 if UCID_1934_2003 == 351777

replace gdename ="koelliken" if UCID_1934_2003 == 351757
replace gdenr = 4276 if UCID_1934_2003 == 351757
replace gdenr_2012 = 4276 if UCID_1934_2003 == 351757
replace gdenr_2018 = 4276 if UCID_1934_2003 == 351757

replace gdename ="basel" if UCID_1934_2003 == 350938
replace gdenr = 2701 if UCID_1934_2003 == 350938
replace gdenr_2012 = 2701 if UCID_1934_2003 == 350938
replace gdenr_2018 = 2701 if UCID_1934_2003 ==350938
			
replace gdename ="taegerig" if UCID_1934_2003 == 350733
replace gdenr = 4077 if UCID_1934_2003 == 350733
replace gdenr_2012 = 4077 if UCID_1934_2003 == 350733
replace gdenr_2018 = 4077 if UCID_1934_2003 ==350733

replace gdename ="seon" if UCID_1934_2003 == 348996
replace gdenr = 4209 if UCID_1934_2003 == 348996
replace gdenr_2012 = 4209 if UCID_1934_2003 == 348996
replace gdenr_2018 = 4209 if UCID_1934_2003 ==348996
	
replace gdename ="triengen" if UCID_1934_2003 == 348794
replace gdenr = 1104 if UCID_1934_2003 == 348794
replace gdenr_2012 = 1104 if UCID_1934_2003 == 348794
replace gdenr_2018 = 1104 if UCID_1934_2003 ==348794
					
replace gdename ="windisch" if UCID_1934_2003 == 348733
replace gdenr = 4123 if UCID_1934_2003 == 348733
replace gdenr_2012 = 4123 if UCID_1934_2003 == 348733
replace gdenr_2018 = 4123 if UCID_1934_2003 ==348733

replace gdename ="graenichen" if UCID_1934_2003 == 347193
replace gdenr = 4006 if UCID_1934_2003 == 347193
replace gdenr_2012 = 4006 if UCID_1934_2003 == 347193
replace gdenr_2018 = 4006 if UCID_1934_2003 ==347193

replace gdename ="hallwil" if UCID_1934_2003 == 347051
replace gdenr = 4197 if UCID_1934_2003 == 347051
replace gdenr_2012 = 4197 if UCID_1934_2003 == 347051
replace gdenr_2018 = 4197 if UCID_1934_2003 ==347051
	
replace gdename ="hallwil" if UCID_1934_2003 == 347034
replace gdenr = 4197 if UCID_1934_2003 == 347034
replace gdenr_2012 = 4197 if UCID_1934_2003 == 347034
replace gdenr_2018 = 4197 if UCID_1934_2003 ==347034
	
replace gdename ="egliswil" if UCID_1934_2003 == 346846
replace gdenr = 4195 if UCID_1934_2003 == 346846
replace gdenr_2012 = 4195 if UCID_1934_2003 == 346846
replace gdenr_2018 = 4195 if UCID_1934_2003 ==346846
	
replace gdename ="oensingen" if UCID_1934_2003 == 346052
replace gdenr = 2407 if UCID_1934_2003 == 346052
replace gdenr_2012 = 2407 if UCID_1934_2003 == 346052
replace gdenr_2018 = 2407 if UCID_1934_2003 ==346052
		
replace gdename ="hallwil" if UCID_1934_2003 == 345582
replace gdenr = 4197 if UCID_1934_2003 == 345582
replace gdenr_2012 = 4197 if UCID_1934_2003 == 345582
replace gdenr_2018 = 4197 if UCID_1934_2003 ==345582

replace gdename ="naters" if UCID_1934_2003 == 344797
replace gdenr =6007  if UCID_1934_2003 == 344797
replace gdenr_2012 = 6007 if UCID_1934_2003 == 344797
replace gdenr_2018 = 6007 if UCID_1934_2003 ==344797

replace gdename ="brugg" if UCID_1934_2003 == 344457
replace gdenr = 4095 if UCID_1934_2003 == 344457
replace gdenr_2012 = 4095 if UCID_1934_2003 == 344457
replace gdenr_2018 = 4095 if UCID_1934_2003 ==344457
	
replace gdename ="effingen" if UCID_1934_2003 == 343696
replace gdenr = 4096 if UCID_1934_2003 == 343696
replace gdenr_2012 = 4096 if UCID_1934_2003 == 343696
replace gdenr_2018 = 4096 if UCID_1934_2003 ==343696
	
replace gdename ="birr" if UCID_1934_2003 == 342926
replace gdenr = 4092 if UCID_1934_2003 == 342926
replace gdenr_2012 = 4092 if UCID_1934_2003 == 342926
replace gdenr_2018 = 4092 if UCID_1934_2003 ==342926
	
replace gdename ="zufikon" if UCID_1934_2003 == 342798
replace gdenr = 4083 if UCID_1934_2003 == 342798
replace gdenr_2012 = 4083 if UCID_1934_2003 == 342798
replace gdenr_2018 = 4083 if UCID_1934_2003 ==342798

replace gdename ="wohlen ag" if UCID_1934_2003 == 342689
replace gdenr = 4082 if UCID_1934_2003 == 342689
replace gdenr_2012 = 4082 if UCID_1934_2003 == 342689
replace gdenr_2018 = 4082 if UCID_1934_2003 ==342689

replace gdename ="wohlen ag" if UCID_1934_2003 == 342686
replace gdenr = 4082 if UCID_1934_2003 == 342686
replace gdenr_2012 = 4082 if UCID_1934_2003 == 342686
replace gdenr_2018 = 4082 if UCID_1934_2003 ==342686

replace gdename ="wohlen ag" if UCID_1934_2003 == 342607
replace gdenr = 4082 if UCID_1934_2003 == 342607
replace gdenr_2012 = 4082 if UCID_1934_2003 == 342607
replace gdenr_2018 = 4082 if UCID_1934_2003 ==342607

replace gdename ="wohlen ag" if UCID_1934_2003 == 342597
replace gdenr = 4082 if UCID_1934_2003 == 342597
replace gdenr_2012 = 4082 if UCID_1934_2003 == 342597
replace gdenr_2018 = 4082 if UCID_1934_2003 ==342597

replace gdename ="wohlen ag" if UCID_1934_2003 == 342520
replace gdenr = 4082 if UCID_1934_2003 == 342520
replace gdenr_2012 = 4082 if UCID_1934_2003 == 342520
replace gdenr_2018 = 4082 if UCID_1934_2003 ==342520

replace gdename ="wohlen ag" if UCID_1934_2003 == 342475
replace gdenr = 4082 if UCID_1934_2003 == 342475
replace gdenr_2012 = 4082 if UCID_1934_2003 == 342475
replace gdenr_2018 = 4082 if UCID_1934_2003 ==342475

replace gdename ="wohlen ag" if UCID_1934_2003 == 342261
replace gdenr = 4082 if UCID_1934_2003 == 342261
replace gdenr_2012 = 4082 if UCID_1934_2003 == 342261
replace gdenr_2018 = 4082 if UCID_1934_2003 ==342261

replace gdename ="boeckten" if UCID_1934_2003 == 341978
replace gdenr = 2842 if UCID_1934_2003 == 341978
replace gdenr_2012 = 2842 if UCID_1934_2003 == 341978
replace gdenr_2018 = 2842 if UCID_1934_2003 ==341978

replace gdename ="wohlen ag" if UCID_1934_2003 == 341355
replace gdenr = 4082 if UCID_1934_2003 == 341355
replace gdenr_2012 = 4082 if UCID_1934_2003 == 341355
replace gdenr_2018 = 4082 if UCID_1934_2003 ==341355

replace gdename ="jonen" if UCID_1934_2003 == 341269
replace gdenr = 4071 if UCID_1934_2003 == 341269
replace gdenr_2012 = 4071 if UCID_1934_2003 == 341269
replace gdenr_2018 = 4071 if UCID_1934_2003 ==341269
	
replace gdename ="haegglingen" if UCID_1934_2003 == 341216
replace gdenr = 4068 if UCID_1934_2003 == 341216
replace gdenr_2012 = 4068 if UCID_1934_2003 == 341216
replace gdenr_2018 = 4068 if UCID_1934_2003 ==341216
				
replace gdename ="dietlikon" if UCID_1934_2003 == 340188
replace gdenr = 54 if UCID_1934_2003 == 340188
replace gdenr_2012 =54 if UCID_1934_2003 == 340188
replace gdenr_2018 =54 if UCID_1934_2003 ==340188
	
replace gdename ="turgi" if UCID_1934_2003 == 338546
replace gdenr = 4042 if UCID_1934_2003 == 338546
replace gdenr_2012 = 4042 if UCID_1934_2003 == 338546
replace gdenr_2018 = 4042 if UCID_1934_2003 ==338546

replace gdename ="zuerich" if UCID_1934_2003 == 338337
replace gdenr = 253 if UCID_1934_2003 == 338337 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 338337 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 338337
replace gdenr_2018 = 261 if UCID_1934_2003 == 338337
			
replace gdename ="zuerich" if UCID_1934_2003 == 338297
replace gdenr = 253 if UCID_1934_2003 == 338297 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 338297 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 338297
replace gdenr_2018 = 261 if UCID_1934_2003 == 338297
						
replace gdename ="ernetschwil" if UCID_1934_2003 == 337816
replace gdenr = 3331 if UCID_1934_2003 == 337816
replace gdenr_2012 = 3331 if UCID_1934_2003 == 337816
replace gdenr_2018 = 3341 if UCID_1934_2003 ==337816
	
replace gdename ="muhen" if UCID_1934_2003 == 336718
replace gdenr = 4009 if UCID_1934_2003 == 336718
replace gdenr_2012 = 4009 if UCID_1934_2003 == 336718
replace gdenr_2018 = 4009 if UCID_1934_2003 ==336718

replace gdename ="jonen" if UCID_1934_2003 == 336041
replace gdenr = 4071 if UCID_1934_2003 == 336041
replace gdenr_2012 = 4071 if UCID_1934_2003 == 336041
replace gdenr_2018 = 4071 if UCID_1934_2003 ==336041

replace gdename ="bern" if UCID_1934_2003 == 335749
replace gdenr = 351 if UCID_1934_2003 == 335749
replace gdenr_2012 = 351 if UCID_1934_2003 == 335749
replace gdenr_2018 = 351 if UCID_1934_2003 ==335749
	
replace gdename ="baden" if UCID_1934_2003 == 334886
replace gdenr = 4021 if UCID_1934_2003 == 334886
replace gdenr_2012 = 4021 if UCID_1934_2003 == 334886
replace gdenr_2018 = 4021 if UCID_1934_2003 ==334886
	
replace gdename ="baden" if UCID_1934_2003 == 334548
replace gdenr = 4021 if UCID_1934_2003 == 334548
replace gdenr_2012 = 4021 if UCID_1934_2003 == 334548
replace gdenr_2018 = 4021 if UCID_1934_2003 ==334548

replace gdename ="birr" if UCID_1934_2003 == 334054
replace gdenr = 4092 if UCID_1934_2003 == 334054
replace gdenr_2012 = 4092 if UCID_1934_2003 == 334054
replace gdenr_2018 = 4092 if UCID_1934_2003 ==334054
	
replace gdename ="ammerswil" if UCID_1934_2003 == 333688
replace gdenr = 4191 if UCID_1934_2003 == 333688
replace gdenr_2012 = 4191 if UCID_1934_2003 == 333688
replace gdenr_2018 = 4191 if UCID_1934_2003 ==333688
				
replace gdename ="suhr" if UCID_1934_2003 == 332925
replace gdenr = 4012 if UCID_1934_2003 == 332925
replace gdenr_2012 = 4012 if UCID_1934_2003 == 332925
replace gdenr_2018 = 4012 if UCID_1934_2003 ==332925

replace gdename ="oberentfelden" if UCID_1934_2003 == 332828
replace gdenr = 4010 if UCID_1934_2003 == 332828
replace gdenr_2012 = 4010 if UCID_1934_2003 == 332828
replace gdenr_2018 = 4010 if UCID_1934_2003 ==332828

replace gdename ="oberentfelden" if UCID_1934_2003 == 332757
replace gdenr = 4010 if UCID_1934_2003 == 332757
replace gdenr_2012 = 4010 if UCID_1934_2003 == 332757
replace gdenr_2018 = 4010 if UCID_1934_2003 ==332757

*Tag wrong match

replace separate = 1 if UCID_1934_2003 == 1374458
replace separate = 1 if UCID_1934_2003 == 1363888
replace separate = 1 if UCID_1934_2003 == 1340346
replace separate = 1 if UCID_1934_2003 == 1302059
replace separate = 1 if UCID_1934_2003 == 1300888
replace separate = 1 if UCID_1934_2003 == 1296371
replace separate = 1 if UCID_1934_2003 == 1279988
replace separate = 1 if UCID_1934_2003 == 1279678
replace separate = 1 if UCID_1934_2003 == 1276065
replace separate = 1 if UCID_1934_2003 == 1245761
replace separate = 1 if UCID_1934_2003 == 1245673
replace separate = 1 if UCID_1934_2003 == 1231289
replace separate = 1 if UCID_1934_2003 == 1185104
replace separate = 1 if UCID_1934_2003 == 1160195
replace separate = 1 if UCID_1934_2003 == 1129361
replace separate = 1 if UCID_1934_2003 == 1117970
replace separate = 1 if UCID_1934_2003 == 1117956
replace separate = 1 if UCID_1934_2003 == 1084238
replace separate = 1 if UCID_1934_2003 == 1080900
replace separate = 1 if UCID_1934_2003 == 1067339
replace separate = 1 if UCID_1934_2003 == 1067116
replace separate = 1 if UCID_1934_2003 == 1064118
replace separate = 1 if UCID_1934_2003 == 1046899
replace separate = 1 if UCID_1934_2003 == 1037557
replace separate = 1 if UCID_1934_2003 == 1032893
replace separate = 1 if UCID_1934_2003 == 1032253
replace separate = 1 if UCID_1934_2003 == 956888
replace separate = 1 if UCID_1934_2003 == 937924
replace separate = 1 if UCID_1934_2003 == 934493
replace separate = 1 if UCID_1934_2003 == 887051
replace separate = 1 if UCID_1934_2003 == 860938
replace separate = 1 if UCID_1934_2003 == 844514
replace separate = 1 if UCID_1934_2003 == 824749
replace separate = 1 if UCID_1934_2003 == 792963
replace separate = 1 if UCID_1934_2003 == 770916
replace separate = 1 if UCID_1934_2003 == 733038
replace separate = 1 if UCID_1934_2003 == 703900
replace separate = 1 if UCID_1934_2003 == 656820
replace separate = 1 if UCID_1934_2003 == 572516
replace separate = 1 if UCID_1934_2003 == 568348
replace separate = 1 if UCID_1934_2003 == 546896
replace separate = 1 if UCID_1934_2003 == 545383
replace separate = 1 if UCID_1934_2003 == 544049
replace separate = 1 if UCID_1934_2003 == 541328
replace separate = 1 if UCID_1934_2003 == 518537
replace separate = 1 if UCID_1934_2003 == 479766
replace separate = 1 if UCID_1934_2003 == 473655
replace separate = 1 if UCID_1934_2003 == 456306
replace separate = 1 if UCID_1934_2003 == 436213
replace separate = 1 if UCID_1934_2003 == 396869
replace separate = 1 if UCID_1934_2003 == 390799
replace separate = 1 if UCID_1934_2003 == 369869
replace separate = 1 if UCID_1934_2003 == 367681
replace separate = 1 if UCID_1934_2003 == 361986
replace separate = 1 if UCID_1934_2003 == 344994

*Correct, just collapse
replace separate = 0 if UCID_1934_2003 == 1969270
replace separate = 0 if UCID_1934_2003 == 988519
replace separate = 0 if UCID_1934_2003 == 977340

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_prob_solved.dta", replace
********************************************************************************
* Putting everything back together

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_ok.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_littau_lu_dedup.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_soc_immo_dedup.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_branch_dedup.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_holding_dedup.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_missing_dedup.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\duplicates_prob_solved.dta"
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_ok2.dta", replace

drop dupdupdup missing_gdenr miss_gdenr_group holding_dummy holding_group branch_group si_dummy_group littau_lu littau_lu_group test test2 dupdupdupdup
order CID UCID_1934_2003 UCID_panel_3 year cname cname_orig gdenr gdenr_2012 gdenr_2018 gdename gdename_orig maxyear maxyear_2 maxyear_3 cname_group owners
bysort UCID_1934_2003 year:  gen test2 = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/

*Corrections Yannick
replace gdename ="ittigen" if UCID_1934_2003 == 109094
replace gdenr = 362 if UCID_1934_2003 == 109094
replace gdenr_2012 = 362 if UCID_1934_2003 == 109094
replace gdenr_2018 = 362 if UCID_1934_2003 == 109094 

replace gdename ="breitenbach" if UCID_1934_2003 == 249441
replace gdenr = 2613 if UCID_1934_2003 == 249441
replace gdenr_2012 = 2613 if UCID_1934_2003 == 249441
replace gdenr_2018 = 2613 if UCID_1934_2003 == 249441 

replace gdename ="romont fr" if UCID_1934_2003 == 628869
replace gdenr = 2096 if UCID_1934_2003 == 628869
replace gdenr_2012 = 2096 if UCID_1934_2003 == 628869
replace gdenr_2018 = 2096 if UCID_1934_2003 == 628869 

replace gdename ="doettingen" if UCID_1934_2003 == 1394855
replace gdenr = 4304 if UCID_1934_2003 == 1394855
replace gdenr_2012 = 4304 if UCID_1934_2003 == 1394855
replace gdenr_2018 = 4304 if UCID_1934_2003 == 1394855 

replace gdename ="volketswil" if UCID_1934_2003 == 1413210
replace gdenr = 199 if UCID_1934_2003 == 1413210
replace gdenr_2012 = 199 if UCID_1934_2003 == 1413210
replace gdenr_2018 = 199 if UCID_1934_2003 == 1413210 

replace gdename ="basel" if UCID_1934_2003 == 1426989
replace gdenr = 2701 if UCID_1934_2003 == 1426989
replace gdenr_2012 = 2701 if UCID_1934_2003 == 1426989
replace gdenr_2018 = 2701 if UCID_1934_2003 == 1426989 

replace gdename ="altdorf ur" if UCID_1934_2003 == 1467843
replace gdenr = 1201 if UCID_1934_2003 == 1467843
replace gdenr_2012 = 1201 if UCID_1934_2003 == 1467843
replace gdenr_2018 = 1201 if UCID_1934_2003 == 1467843 

replace gdename ="saanen" if UCID_1934_2003 == 1494183
replace gdenr = 843 if UCID_1934_2003 == 1494183
replace gdenr_2012 = 843 if UCID_1934_2003 == 1494183
replace gdenr_2018 = 843 if UCID_1934_2003 == 1494183

replace gdename ="oberdiessbach" if UCID_1934_2003 == 1498103 & year > 1990
replace gdenr = 619 if UCID_1934_2003 == 1498103 & year > 1990
replace gdenr_2012 = 619 if UCID_1934_2003 == 1498103 & year > 1990
replace gdenr_2018 = 619 if UCID_1934_2003 == 1498103 & year > 1990

replace gdename ="thun" if UCID_1934_2003 == 1498103 & year < 1991
replace gdenr = 942 if UCID_1934_2003 == 1498103 & year < 1991
replace gdenr_2012 = 942 if UCID_1934_2003 == 1498103 & year < 1991
replace gdenr_2018 = 942 if UCID_1934_2003 == 1498103 & year < 1991     

replace gdename ="triengen" if UCID_1934_2003 == 1513244
replace gdenr = 1104 if UCID_1934_2003 == 1513244
replace gdenr_2012 = 1104 if UCID_1934_2003 == 1513244
replace gdenr_2018 = 1104 if UCID_1934_2003 == 1513244

replace gdename ="galgenen" if UCID_1934_2003 == 1524478
replace gdenr = 1342 if UCID_1934_2003 == 1524478
replace gdenr_2012 = 1342 if UCID_1934_2003 == 1524478
replace gdenr_2018 = 1342 if UCID_1934_2003 == 1524478

replace gdename ="zug" if UCID_1934_2003 == 1532585
replace gdenr = 1711 if UCID_1934_2003 == 1532585
replace gdenr_2012 = 1711 if UCID_1934_2003 == 1532585
replace gdenr_2018 = 1711 if UCID_1934_2003 == 1532585

replace gdename ="glarus" if UCID_1934_2003 == 1536279
replace gdenr = 1609 if UCID_1934_2003 == 1536279
replace gdenr_2012 = 1632 if UCID_1934_2003 == 1536279
replace gdenr_2018 = 1632 if UCID_1934_2003 == 1536279

replace gdename ="naters" if UCID_1934_2003 == 1543371 
replace gdenr = 6007 if UCID_1934_2003 == 1543371
replace gdenr_2012 = 6007 if UCID_1934_2003 == 1543371
replace gdenr_2018 = 6007 if UCID_1934_2003 == 1543371

replace gdename ="basel" if UCID_1934_2003 == 1631583
replace gdenr = 2701 if UCID_1934_2003 == 1631583
replace gdenr_2012 = 2701 if UCID_1934_2003 == 1631583
replace gdenr_2018 = 2701 if UCID_1934_2003 == 1631583  

replace gdename ="kaltbrunn" if UCID_1934_2003 == 1671765
replace gdenr = 3313 if UCID_1934_2003 == 1671765
replace gdenr_2012 = 3313 if UCID_1934_2003 == 1671765
replace gdenr_2018 = 3313 if UCID_1934_2003 == 1671765 

replace gdename ="zug" if UCID_1934_2003 == 1683980
replace gdenr = 1711 if UCID_1934_2003 == 1683980
replace gdenr_2012 = 1711 if UCID_1934_2003 == 1683980
replace gdenr_2018 = 1711 if UCID_1934_2003 == 1683980

replace gdename ="zug" if UCID_1934_2003 == 1704679
replace gdenr = 1711 if UCID_1934_2003 == 1704679
replace gdenr_2012 = 1711 if UCID_1934_2003 == 1704679
replace gdenr_2018 = 1711 if UCID_1934_2003 == 1704679

replace gdename ="romanshorn" if UCID_1934_2003 == 1724098
replace gdenr = 4436 if UCID_1934_2003 == 1724098
replace gdenr_2012 = 4436 if UCID_1934_2003 == 1724098
replace gdenr_2018 = 4436 if UCID_1934_2003 == 1724098

replace gdename ="lugano" if UCID_1934_2003 == 1744583
replace gdenr = 5192 if UCID_1934_2003 == 1744583
replace gdenr_2012 = 5192 if UCID_1934_2003 == 1744583
replace gdenr_2018 = 5192 if UCID_1934_2003 == 1744583 
 
replace gdename ="melide" if UCID_1934_2003 == 1778376
replace gdenr = 5198 if UCID_1934_2003 == 1778376
replace gdenr_2012 = 5198 if UCID_1934_2003 == 1778376
replace gdenr_2018 = 5198 if UCID_1934_2003 == 1778376  

replace gdename ="ayent" if UCID_1934_2003 == 1826357
replace gdenr = 6082 if UCID_1934_2003 == 1826357
replace gdenr_2012 = 6082 if UCID_1934_2003 == 1826357
replace gdenr_2018 = 6082 if UCID_1934_2003 == 1826357    

replace gdename ="pully" if UCID_1934_2003 == 1827539
replace gdenr = 5590 if UCID_1934_2003 == 1827539
replace gdenr_2012 = 5590 if UCID_1934_2003 == 1827539
replace gdenr_2018 = 5590 if UCID_1934_2003 == 1827539 

replace gdename ="zuerich" if UCID_1934_2003 == 2729058
replace gdenr = 253 if UCID_1934_2003 == 2729058 & year < 1990
replace gdenr = 261 if UCID_1934_2003 == 2729058 & year > 1989
replace gdenr_2012 = 261 if UCID_1934_2003 == 2729058
replace gdenr_2018 = 261 if UCID_1934_2003 == 2729058

replace gdename ="freienbach" if UCID_1934_2003 == 1393916
replace gdenr = 1322 if UCID_1934_2003 == 1393916
replace gdenr_2012 = 1322 if UCID_1934_2003 == 1393916
replace gdenr_2018 = 1322 if UCID_1934_2003 == 1393916

replace gdename ="zug" if UCID_1934_2003 == 1413321
replace gdenr = 1711 if UCID_1934_2003 == 1413321
replace gdenr_2012 = 1711 if UCID_1934_2003 == 1413321
replace gdenr_2018 = 1711 if UCID_1934_2003 == 1413321

replace gdename ="basel" if UCID_1934_2003 == 1426908
replace gdenr = 2701 if UCID_1934_2003 == 1426908
replace gdenr_2012 = 2701 if UCID_1934_2003 == 1426908
replace gdenr_2018 = 2701 if UCID_1934_2003 == 1426908        

replace gdename ="manno" if UCID_1934_2003 == 1765605
replace gdenr = 5194 if UCID_1934_2003 == 1765605
replace gdenr_2012 = 5194 if UCID_1934_2003 == 1765605
replace gdenr_2018 = 5194 if UCID_1934_2003 == 1765605 

replace gdename ="morges" if UCID_1934_2003 == 1809932
replace gdenr = 5642 if UCID_1934_2003 == 1809932
replace gdenr_2012 = 5642 if UCID_1934_2003 == 1809932
replace gdenr_2018 = 5642 if UCID_1934_2003 == 1809932

replace gdename ="zug" if UCID_1934_2003 == 1538702
replace gdenr = 1711 if UCID_1934_2003 == 1538702
replace gdenr_2012 = 1711 if UCID_1934_2003 == 1538702
replace gdenr_2018 = 1711 if UCID_1934_2003 == 1538702

replace gdename ="zug" if UCID_1934_2003 == 1561726
replace gdenr = 1711 if UCID_1934_2003 == 1561726
replace gdenr_2012 = 1711 if UCID_1934_2003 == 1561726
replace gdenr_2018 = 1711 if UCID_1934_2003 == 1561726

replace gdename ="vaz/obervaz" if UCID_1934_2003 == 1695633
replace gdenr = 3506 if UCID_1934_2003 == 1695633
replace gdenr_2012 = 3506 if UCID_1934_2003 == 1695633
replace gdenr_2018 = 3506 if UCID_1934_2003 == 1695633

replace gdename ="lausanne" if UCID_1934_2003 == 2256924
replace gdenr = 5586 if UCID_1934_2003 == 2256924
replace gdenr_2012 = 5586 if UCID_1934_2003 == 2256924
replace gdenr_2018 = 5586 if UCID_1934_2003 == 2256924

replace gdename ="nidau" if CID == 110861 & year == 1998
replace gdenr = 743 if CID == 110861 & year == 1998
replace gdenr_2012 = 743 if CID == 110861 & year == 1998
replace gdenr_2018 = 743 if CID == 110861 & year == 1998

replace gdename ="pully" if CID == 119699 & year == 1997
replace gdenr = 5590 if CID == 110861 & year == 1997
replace gdenr_2012 = 5590 if CID == 110861 & year == 1997
replace gdenr_2018 = 5590 if CID == 110861 & year == 1997

replace gdename ="lausanne" if CID == 41987 & year == 1992
replace gdenr = 5586 if CID == 41987 & year == 1992
replace gdenr_2012 = 5586 if CID == 41987 & year == 1992
replace gdenr_2018 = 5586 if CID == 41987 & year == 1992

replace gdename ="lausanne" if CID == 126341 & year == 1995
replace gdenr = 5586 if CID == 126341 & year == 1995
replace gdenr_2012 = 5586 if CID == 126341 & year == 1995
replace gdenr_2018 = 5586 if CID == 126341 & year == 1995 

replace gdename ="geneve" if CID == 50957 & year == 1972
replace gdenr = 6621 if CID == 50957 & year == 1972
replace gdenr_2012 = 6621 if CID == 50957 & year == 1972
replace gdenr_2018 = 6621 if CID == 50957 & year == 1972 

replace gdename ="buelach" if CID == 120547 & year == 1989
replace gdenr = 53 if CID == 120547 & year == 1989
replace gdenr_2012 = 53 if CID == 120547 & year == 1989
replace gdenr_2018 = 53 if CID == 120547 & year == 1989 

replace gdename ="dietlikon" if CID == 169163 & year == 1997
replace gdenr = 54 if CID == 169163 & year == 1997
replace gdenr_2012 = 54 if CID == 169163 & year == 1997
replace gdenr_2018 = 54 if CID == 169163 & year == 1997 

replace gdename ="bulle" if UCID_1934_2003 == 536588
replace gdenr = 2125 if UCID_1934_2003 == 536588
replace gdenr_2012 = 2125 if UCID_1934_2003 == 536588
replace gdenr_2018 = 2125 if UCID_1934_2003 == 536588

replace gdename ="vicques" if UCID_1934_2003 == 536044
replace gdenr = 6727 if UCID_1934_2003 == 536044
replace gdenr_2012 = 6727 if UCID_1934_2003 == 536044
replace gdenr_2018 = 6730 if UCID_1934_2003 == 536044

replace gdename ="vicques" if UCID_1934_2003 == 536074
replace gdenr = 6727 if UCID_1934_2003 == 536074
replace gdenr_2012 = 6727 if UCID_1934_2003 == 536074
replace gdenr_2018 = 6730 if UCID_1934_2003 == 536074

replace gdename ="onex" if UCID_1934_2003 == 534153
replace gdenr = 6631 if UCID_1934_2003 == 534153
replace gdenr_2012 = 6631 if UCID_1934_2003 == 534153
replace gdenr_2018 = 6631 if UCID_1934_2003 == 534153

*Tag wrong matches
replace separate = 1 if UCID_1934_2003 == 22564
replace separate = 1 if UCID_1934_2003 == 138699
replace separate = 1 if UCID_1934_2003 == 1777080

*Correct, just collapse
replace separate = 0 if UCID_1934_2003 == 86301
replace separate = 0 if UCID_1934_2003 == 459939
replace separate = 0 if UCID_1934_2003 == 546882
replace separate = 0 if UCID_1934_2003 == 600037
replace separate = 0 if UCID_1934_2003 == 601013
replace separate = 0 if UCID_1934_2003 == 601147
replace separate = 0 if UCID_1934_2003 == 976371
replace separate = 0 if UCID_1934_2003 == 1000996
replace separate = 0 if UCID_1934_2003 == 1017649
replace separate = 0 if UCID_1934_2003 == 1069813
replace separate = 0 if UCID_1934_2003 == 1293160
replace separate = 0 if UCID_1934_2003 == 1326652

****
egen maxUCID_1934_2003_2 = max(UCID_1934_2003)

* Separation if separate == 1 *
bysort UCID_1934_2003 gdenr: gen newid = _n if _n==1
replace newid = sum(newid)
sort UCID_1934_2003 year gdenr

*corrections
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 1 if UCID_1934_2003 == 11790 & cname_orig == "Baumgartner Training AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 2 if UCID_1934_2003 == 13665 & cname_orig == "Elfag Immobilien AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 3 if UCID_1934_2003 == 22564 & cname_orig == "Walder U. AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 4 if UCID_1934_2003 == 22564 & cname_orig == "Walder, U. AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 5 if UCID_1934_2003 == 23576 & newid == 21983
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 6 if UCID_1934_2003 == 25358 & newid == 23945
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 7 if UCID_1934_2003 == 25601 & newid == 24161
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 8 if UCID_1934_2003 == 26615 & newid == 25151
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 9 if UCID_1934_2003 == 27048 & newid == 25547
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 10 if UCID_1934_2003 == 30396 & newid == 28723
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 11 if UCID_1934_2003 == 30720 & newid == 29049
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 12 if UCID_1934_2003 == 32160 & newid == 30323
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 13 if newid == 33039
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 14 if cname_orig == "Solaret AG" & UCID_1934_2003 == 44292
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 15 if UCID_1934_2003 == 47020 & cname_orig == "Landhus AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 16 if UCID_1934_2003 == 50228 & newid == 51016
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 17 if UCID_1934_2003 == 50228 & newid == 51018
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 18 if UCID_1934_2003 == 50228 & newid == 51019
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 19 if UCID_1934_2003 == 50230 & cname_orig == "Cofitrac AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 19 if UCID_1934_2003 == 50230 & cname_orig == "Cofitrac AG."
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 19 if UCID_1934_2003 == 50230 & cname_orig == "CofitracAG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 20 if UCID_1934_2003 == 58438 & newid == 62096
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 21 if cname_orig == "S.A. Müller-Maschinen" & CID == 36834
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 22 if newid == 67248
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 23 if newid == 69334
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 24 if newid == 71727
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 25 if newid == 76626
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 26 if UCID_1934_2003 == 70484 & newid == 78004
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 27 if newid == 84727
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 28 if UCID_1934_2003 == 76800 & cname_orig == "Valinvest AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 29 if UCID_1934_2003 == 76800 & cname_orig == "Varinvest AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 30 if UCID_1934_2003 == 89993 & cname_orig == "Gedreide-Silo AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 31 if UCID_1934_2003 == 101347 & newid == 114900
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 32 if UCID_1934_2003 == 101347 & newid == 114901
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 33 if UCID_1934_2003 == 104516 & cname_orig == "Berimo B Immobilien AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 34 if UCID_1934_2003 == 106199 & newid == 119903
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 35 if UCID_1934_2003 == 106637 & cname_orig == "Gustav Renfer AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 36 if UCID_1934_2003 == 108515 & newid == 122096
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 37 if UCID_1934_2003 == 110216 & newid == 123786
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 38 if UCID_1934_2003 == 115014 & newid == 128608
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 39 if UCID_1934_2003 == 127074 & newid == 140252
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 40 if UCID_1934_2003 == 127754 & newid == 140900

replace UCID_1934_2003 = maxUCID_1934_2003_2 + 41 if UCID_1934_2003 == 138699 & cname_orig == "Weber Holding AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 42 if UCID_1934_2003 == 138699 & newid == 151573
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 43 if UCID_1934_2003 == 344994 & cname_orig == "AND-Bau AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 44 if UCID_1934_2003 == 361986 & newid == 373280
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 45 if UCID_1934_2003 == 367681 & newid == 379350
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 46 if UCID_1934_2003 == 369869 & cname_orig == "Fudalta SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 47 if UCID_1934_2003 == 390799 & cname_orig == "Sof-Tech SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 48 if UCID_1934_2003 == 390799 & cname_orig == "Sofiter S.A."
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 48 if UCID_1934_2003 == 390799 & cname_orig == "Sofiter SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 49 if UCID_1934_2003 == 390799 & cname_orig == "Sofinder SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 50 if UCID_1934_2003 == 396869 & newid == 409170
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 51 if UCID_1934_2003 == 436213 & cname_orig == "VRP & partners consulting SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 52 if UCID_1934_2003 == 456306 & cname_orig == "Valimag Immobilien AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 53 if UCID_1934_2003 == 456306 & cname_orig == "Arimag-Immobilien AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 53 if UCID_1934_2003 == 456306 & cname_orig == "Arimag-Immobilien AG."
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 54 if UCID_1934_2003 == 473655 & newid == 484699
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 55 if UCID_1934_2003 == 479766 & newid == 490863
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 56 if UCID_1934_2003 == 518537 & cname_orig == "TAG AVIATION HOLDING SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 57 if UCID_1934_2003 == 545383 & newid == 557865
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 58 if UCID_1934_2003 == 546896 & newid == 559744
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 59 if UCID_1934_2003 == 568348 & cname_orig == "Mentor Immobilien AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 60 if UCID_1934_2003 == 572516 & cname_orig == "IC Consult AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 61 if UCID_1934_2003 == 703900 & cname_orig == "Valoris SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 62 if UCID_1934_2003 == 733038 & cname_orig == "LDB FINANCE SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 63 if UCID_1934_2003 == 733038 & cname_orig == "PB Finance SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 64 if UCID_1934_2003 == 824749 & cname_orig == "SoninvestAG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 64 if UCID_1934_2003 == 824749 & cname_orig == "Soninvest AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 65 if UCID_1934_2003 == 860938 & cname_orig == "Sefinvest Holding SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 66 if UCID_1934_2003 == 860938 & cname_orig == "Cefinvest SA."
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 67 if UCID_1934_2003 == 887051 & cname_orig == "Gentras AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 68 if UCID_1934_2003 == 956888 & cname_orig == "Parfina SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 69 if UCID_1934_2003 == 956888 & cname_orig == "Partiga SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 70 if UCID_1934_2003 == 1032253 & cname_orig == "RR Capital Management AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 71 if UCID_1934_2003 == 1032253 & cname_orig == "M2 Capital Management AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 72 if UCID_1934_2003 == 1032893 & cname_orig == "HS Management Consulting AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 73 if UCID_1934_2003 == 1032893 & cname_orig == "RAL Management Consulting AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 74 if UCID_1934_2003 == 1037557 & cname_orig == "LS Consulting AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 75 if UCID_1934_2003 == 1046899 & cname_orig == "Progeco SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 76 if UCID_1934_2003 == 1064118 & cname_orig == "Mantrex SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 77 if UCID_1934_2003 == 1067116 & cname_orig == "Bütler & Partner AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 78 if UCID_1934_2003 == 1067339 & cname_orig == "CorrectaAG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 79 if UCID_1934_2003 == 1117970 & cname_orig == "WTT Fördertechnik AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 80 if UCID_1934_2003 == 1185104 & newid == 574703

replace UCID_1934_2003 = maxUCID_1934_2003_2 + 81 if UCID_1934_2003 == 1231289 & newid == 575598
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 82 if UCID_1934_2003 == 1245673 & cname_orig == "City-Immobilien AG, Luzern"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 82 if UCID_1934_2003 == 1245673 & cname_orig == "City-Immobilien AG. Luzern"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 83 if UCID_1934_2003 == 1245761 & cname_orig == "Limag Immobilien AG Luzern"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 84 if UCID_1934_2003 == 1276065 & cname_orig == "Artema Trasporti Internazionali S.A."
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 85 if UCID_1934_2003 == 1279678 & cname_orig == "Schaltag AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 86 if UCID_1934_2003 == 1279988 & cname_orig == "Primobau Immobilien AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 87 if UCID_1934_2003 == 1296371 & cname_orig == "Remonta AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 88 if UCID_1934_2003 == 1302059 & cname_orig == "Latinvest AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 89 if UCID_1934_2003 == 1340346 & cname_orig == "Casablanca SA"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 90 if UCID_1934_2003 == 1363888 & cname_orig == "liau-& Grundstücke AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 91 if UCID_1934_2003 == 1374458 & cname_orig == "Shearson Lehman Brothers Inc, Wilmington"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 91 if UCID_1934_2003 == 1374458 & cname_orig == "Lehman Brothers Inc, Wilmington"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 92 if UCID_1934_2003 == 1374458 & cname_orig == "American Express International, Inc, Wilmington"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 92 if UCID_1934_2003 == 1374458 & cname_orig == "American Express International. Inc, Wilmington"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 93 if UCID_1934_2003 == 1374458 & cname_orig == "Odyssey Advanced Financial Solutions SA, Luxembourg"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 94 if UCID_1934_2003 == 1374458 & cname_orig == "JCH PULP AND PAPER INd., Montréal, succursale de Genève, Genève SOCIETE GESTRADE ET FINANCE INC, Montreal, succursale de Lausanne"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 95 if UCID_1934_2003 == 1777080 & cname_orig == "Steminvest AG Holding"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 96 if UCID_1934_2003 == 2896 & cname_orig == "Clip Holding AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 97 if UCID_1934_2003 == 550151 & cname_orig == "CCS Compact Computer Systeme AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 98 if UCID_1934_2003 == 545079 & cname_orig == "Interautmation Holding AG"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 98 if UCID_1934_2003 == 545079 & cname_orig == "Interautmation Holding AC"
replace UCID_1934_2003 = maxUCID_1934_2003_2 + 99 if UCID_1934_2003 == 541583 & cname_orig == "Wyrsch AG"

* collapse if duplicates in terms of UCID_1934_2003 year gdenr_2018 (maybe w/o gdenr_2018)
sort UCID_1934_2003 year gdenr_2018
quietly by UCID_1934_2003 year gdenr_2018:  gen test3 = cond(_N==1,0,_n)
quietly by UCID_1934_2003 year:  gen test4 = cond(_N==1,0,_n) /*Exact duplicates for all vars!*/

drop comma
replace owners = usubinstr(owners,  ";", ",",.) 
gen nb_owners = (length(trim(owners))-length(trim(subinstr(owners,",","",.))))+1

sort UCID_1934_2003 year maxyear_3 maxyear nb_owners 
	foreach y in owners function signature page CID_str cname_group {
		by UCID_1934_2003 year: replace `y' = `y'[_n-1] + "; " + `y' if test3<test4
		replace `y' = " "+`y' if test3<test4
		replace `y' = usubinstr(`y', " ; ", "",.) if test3<test4
		by UCID_1934_2003 year: replace `y' = `y'[_N] if test3<test4
		}

sort UCID_1934_2003 year maxyear_3 maxyear nb_owners test4 	
by UCID_1934_2003 year: gen id2=(_n)
		
drop if id2 == 1 & test3<test4	

*manually kick ones left:

drop if CID == 29879 & year == 1975
drop if CID == 29876 & year == 1975
drop if CID == 1543 & year == 1969	
drop if CID == 71502 & year == 1975
drop if CID == 5126 & year == 1966
drop if CID == 4843 & year == 1965
drop if UCID_1934_2003 == 138763 & id2 == 2 
drop if UCID_1934_2003 == 245155 & id2 == 2 
drop if CID == 7133 & year == 1964
drop if CID == 56341 & year == 1965
drop if CID == 58224 & year == 1966
drop if CID == 51627 & year == 1972
drop if CID == 21285 & year == 1962
drop if CID == 39201 & year == 1960
drop if CID == 50337 & year == 1963
drop if CID == 53148 & year == 1964
drop if CID == 56094 & year == 1965
drop if CID == 57975 & year == 1966
drop if CID == 65215 & year == 1969
drop if CID == 28202 & year == 1963
drop if CID == 29914 & year == 1964
drop if CID == 31549 & year == 1965
drop if CID == 1784 & year == 1943
drop if CID == 10761 & year == 1963
drop if CID == 60770 & year == 1972
drop if CID == 53553 & year == 1966
drop if CID == 41990 & year == 1969
drop if CID == 122388 & year == 2000
drop if CID == 98231 & year == 1997
drop if CID == 95232 & year == 1980
		
* Save that file as ok3
bysort UCID_1934_2003 year: gen id3=(_n)

	foreach y in owners function signature page CID_str cname_group {
		bysort UCID_1934_2003 year gdenr_2018: replace `y' = `y'[_n-1] + "; " + `y'
		replace `y' = " "+`y'
		replace `y' = usubinstr(`y', " ; ", "",.)
		bysort UCID_1934_2003 year gdenr_2018: replace `y' = `y'[_N]
		}
	
drop if id3 == 2

drop maxUCID_temp4 FP_group_2 tag REAL_PROBLEM UCID_panel_2_duplicates REAL_PROBLEM_2 UCID_panel_3_backup group_YS group_JU group_JU_temp group_ED group_ED_temp match_ok correct_match S T checked separate test2 newid test3 test4 id2 id3 id4 dup dup2 tag_2 obs_total_2 dup_panel_2 minyear_2 maxUCID_panel_2 tag_3 obs_total_3 dup_panel_3 minyear_3 match_1 match_2 match_3 dup_panel minyear FalsePositive_group FalsePositive_orig obs_total maxUCID_panel_3 maxUCID_1934_2003 maxUCID_1934_2003_2

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_ok3.dta", replace

********************************************************************************
* Matching t and t+2 
********************************************************************************
********************************************************************************
* Subsample of beginning and ending of series
********************************************************************************

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_ok3.dta", clear 
sort UCID_1934_2003 year

gen tag_4 = 1
egen obs_total_4 = total(tag_4), by(UCID_1934_2003) //How many years in a row perfectly matched
quietly by UCID_1934_2003:  gen dup_panel_4 = cond(_N==1,0,_n)
bysort UCID_1934_2003: egen maxyear_4 = max(dup_panel_4)
bysort UCID_1934_2003: egen minyear_4 = min(dup_panel_4)

*Obs in bewteen first and last appearance (100% correct and not usefull to further match, same mun, same exact name, two "years" in a row)
preserve
drop if dup_panel_4 == maxyear_4 | dup_panel_4 == minyear_4
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_middle_ok_3.dta", replace
restore

*First and last appareance in group starting in 1934 or ending in 2003 (not single occurences) 
preserve 
keep if year == 1934 | year == 2003
keep if obs_total_4 > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_first_end_ok_3.dta", replace
restore

*First and last appareance of a firm in each group : unmatched subsample
preserve
keep if dup_panel_4 == maxyear_4 | dup_panel_4 == minyear_4
drop if year == 2003 & obs_total_4 > 1
drop if year == 1934 & obs_total_4 > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_notok_3.dta", replace
restore

*Append two "ok" subsample
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_middle_ok_3.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_first_end_ok_3.dta"
sort UCID_1934_2003 year
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_ok_last.dta", replace

********************************************************************************
* 10% matches across municipalities in year t and t+2 (to find firms that disappear one year and reappear)

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_notok_3.dta", clear
drop if missing(cname) /* 0 Obs without any firm name*/

*Split in cross-sections
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp{
	preserve
		keep if year == `a'
		save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_notok_4_temp.dta", replace
	restore
}

* Build pair of years in t and t+2
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001
local bgrp 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
    use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_notok_4_temp.dta", clear
	append using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`b'_notok_4_temp.dta"
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_`b'_notok_4_temp.dta", replace 
}

*Match across mun, over pairs of year
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001
local bgrp 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp      
    use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_geo_`a'_`b'_notok_4_temp.dta", clear
		strgroup cname, gen(UCID_`a'_`b') thresh(0.1) force
		order CID UCID_`a'_`b'
		sort UCID_`a'_`b' year
	save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_`a'_`b'_4.dta", replace
}


********************************************************************************
*Matches that should be evaluated
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001
local bgrp 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp  
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_`a'_`b'_4.dta", clear
bysort UCID_`a'_`b' :  gen n = cond(_N==1,0,_n)
drop if n == 0
bysort UCID_`a'_`b' year :  gen n2 = cond(_N==1,0,_n)
drop if n == n2
bysort UCID_`a'_`b' UCID_1934_2003: gen n3 = cond(_N==1,0,_n)
gen zero = 0 if n3 == 0
bysort UCID_`a'_`b': replace zero = zero[_n-1] if missing(zero)
forval i = 1/15 {
bysort UCID_`a'_`b': replace zero = zero[_n+1] if missing(zero)
}
keep if zero == 0
sort UCID_`a'_`b' year
save "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_`a'_`b'_toeval.dta", replace
}

********************************************************************************
* Test to link owners ID with owners name

local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001
local bgrp 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp
	
use "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_`a'_`b'_toeval.dta", clear

replace owners = usubinstr(owners, ";", ",",.)
replace owners = usubinstr(owners, "   ", "  ",.)
replace owners = usubinstr(owners, "  ", " ",.)
replace owners = usubinstr(owners, " ", "",.)
split owners, parse(,) generate(owner)

generate id = _n
reshape long owner, i(UCID_1934_2003 id) j(test)

drop if missing(owner)
rename owner PID
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_`a'_`b'_toeval_long.dta", replace

*Cleaning of persons list' (to have single PID per year)
use PID year firstname lastname CID using "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", clear

quietly by PID year:  gen dup = cond(_N==1,0,_n)
drop if dup > 1

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\PID_year_list.dta", replace

*Merge
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_`a'_`b'_toeval_long.dta", clear
destring PID, replace
merge m:1 PID year using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\PID_year_list.dta"
keep if _merge == 3
drop fullname
g fullname = firstname + " " + lastname
drop firstname lastname dup _merge

*Reshape to have one unique CID with multiple var with owners
sort CID year
bysort CID year: replace fullname = fullname[_n-1] + "; " + fullname
replace fullname = " "+fullname
replace fullname = usubinstr(fullname, " ; ", "",.)
by CID year: replace fullname = fullname[_N]

by CID year:  gen dup = cond(_N==1,0,_n)
drop if dup > 1
sort UCID_panel_3 year
keep CID year cname fullname

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\CID_year_owners_`a'_`b'_toeval_long.dta", replace
}

* Add full names of owners to obs to evaluate
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001
local bgrp 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {    
    gettoken b bgrp : bgrp  
use "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_`a'_`b'_toeval.dta", clear
drop fullname
merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\CID_year_owners_`a'_`b'_toeval_long.dta", nogen
sort UCID_`a'_`b' year
save "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_`a'_`b'_toeval_names.dta", replace
drop minyear_4	n	n2	n3	zero maxyear	maxyear_2	maxyear_3	cname_group	owners	UCID_3	UCID_1934_2003_temp	UCID_panel_1	UCID_panel_2	periode	si_dummy_2	CID_str	capital	function	signature	page	caddress_orig	geo_merge	branch_dummy	gdename_extract	capital_extract	cname_extract	ctn	E_CNTR	N_CNTR	nb_owners	tag_4	obs_total_4	UCID_panel_3

export excel using "$dump\02_Processed_data\11_Matching_Eval\match_10_nogeo_`a'_`b'_toeval_names.xlsx", replace firstrow(variables)
}


********************************************************************************
* Spotting gaps and constructing dataset with them
********************************************************************************

* Flagging creation of firm, relocation and deaths
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_ok3.dta", clear

*Create files to be evaluated by hand
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
local bgrp 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
foreach a of local agrp { 
	gettoken b bgrp : bgrp
	replace periode = `b' if year == `a'
}

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_ok3_temp.dta", replace
