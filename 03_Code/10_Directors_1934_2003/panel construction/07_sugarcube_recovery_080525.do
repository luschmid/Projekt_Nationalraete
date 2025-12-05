version 17
clear all
set more off
cap log close

global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\" /* Office */
global dump "C:\Users\schmutzy\Desktop\Projekt Nationalräte\" /* Office Dumpster */
global ED "C:\Users\schmutzy\Dropbox\Travaux_EmilieDousse\" /* Office ED */

******************************************
* Preparing clean dataset for Mark 


*********************************************
* Add complete original info

/*use  "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Key_CID_UCID.dta", clear

merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_complete.dta"

We still have 9752 missing observations compared to the original sugarcube digitalized raw data */

* Retrieve CIDs lost in deduplicate step
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_toeval.dta", clear
keep CID year CID_str
split CID_str, parse(;) generate(CID_recovered)
drop if missing(CID_recovered2)
rename CID_str CID_rec
drop CID_recovered*
save "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\30_Recovery\recovery_`y'.dta", replace
}

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale.dta", clear
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
merge 1:1 CID year using "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\30_Recovery\recovery_`y'.dta", keep(1 3)
replace CID_str = CID_str + ";" + CID_rec if CID_rec != ""
drop _merge
drop CID_rec
}

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale_2.dta", replace

********************************************************************************
* Final Panel - all CIDs 
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale_2.dta", clear

keep CID UCID_1934_2003 year periode CID_str
split CID_str, parse(;) generate(CID_grp)
drop CID_str
destring CID_grp*, replace force 

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale_2_parsed.dta", replace

********************************************************************************
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale_2_parsed.dta", clear

rename CID CID_orig
rename CID_grp1 CID

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
preserve
merge 1:1 CID year using "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\30_Recovery\recovery_`y'_2.dta", keep(1 3)

drop if missing(CID_rec)

save "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\30_Recovery\Missing_CID_`y'.dta", replace
restore
}

* Append all years
use "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\30_Recovery\Missing_CID_1934.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\30_Recovery\Missing_CID_`y'.dta"
}
drop CID_grp* CID _merge periode
rename CID_orig CID
save "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\30_Recovery\Missing_CID_1934_2003.dta", replace

********************************************************************************

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale_2.dta", clear
merge 1:1 CID year using "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\30_Recovery\Missing_CID_1934_2003.dta"

replace CID_str = CID_str + ";" + CID_rec if _merge == 3
drop _merge CID_rec

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale_3.dta", replace

****************************************
* All CIDs : run on server !
**************************************** 

/* Finale Panel
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale_3.dta", clear

keep CID UCID_1934_2003 year periode CID_str
split CID_str, parse(;) generate(CID_grp)
drop CID_str
destring CID_grp*, replace force 

generate id = _n
order id UCID_1934_2003
reshape long CID_grp, i(UCID_1934_2003 id) j(test)
drop if missing(CID_grp)
drop id test CID
rename CID_grp CID

bysort CID year:  gen dup = cond(_N==1,0,_n)
drop dup

sort UCID_1934_2003 year

save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Key_CID_UCID.dta", replace*/

*Run on server using cross-sections and reshape long: file is  "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Firms_Key_ID.dta"

********************************************************************************
use  "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Firms_Key_ID.dta", clear
drop dup
bysort CID year:  gen dup = cond(_N==1,0,_n)
drop if dup > 1 // I appended 1934 after using 1934 : the whole cross-section is a duplicate

merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_complete.dta"
********************************************************************************
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale_3.dta", clear

