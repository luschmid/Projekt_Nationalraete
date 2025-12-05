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
* Construct cross-sections of total panel

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Firm_Panel_1934_2003.dta", clear
merge m:m gdenr_2018 using "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018.dta", keep(1 3) nogen
rename language language_temp
g language = 1 if language_temp == "Deutsches Sprachgebiet"
replace language = 2 if language_temp == "Französisches Sprachgebiet"
replace language = 3 if language_temp == "Italienisches Sprachgebiet"
replace language = 4 if language_temp == "Rätoromanisches Sprachgebiet"
drop gem_cd gem_name 

local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp{
	preserve	
		keep if year == `a'
		save "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`a'.dta", replace
	restore
	drop if year == `a'
}


***************************************************
*** Generating list of all owners for each firm ***
***************************************************

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'.dta", clear

*** Generate list of owners with single firm IDs ***
replace owners = usubinstr(owners, ";", ",",.)
replace owners = usubinstr(owners, "   ", "  ",.)
replace owners = usubinstr(owners, "  ", " ",.)
replace owners = usubinstr(owners, " ", "",.)
split owners, parse(,) generate(owner)
keep UCID_1934_2003 owner*
drop owners
sort UCID_1934_2003 owner*
quietly by UCID_1934_2003 owner*: gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*Deal with multiple identical IDs
generate id = _n
order id UCID_1934_2003
reshape long owner, i(UCID_1934_2003 id) j(test)
drop id test
drop if owner==""
drop if UCID_1934_2003==.
byso UCID_1934_2003: gen id=(_n)
destring owner, replace
reshape wide owner, i(UCID_1934_2003) j(id)

save "$dump\02_Processed_data\10_Directors_1934_2003\11_Ownerlists\ownerlist_`y'.dta", replace
}

******************************************************************
*Nb of owners deduplicated by UCID_1934_2003 : to link with cross-sections
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'.dta", clear

*** Generate list of owners with single firm IDs ***
replace owners = usubinstr(owners, ";", ",",.)
replace owners = usubinstr(owners, "   ", "  ",.)
replace owners = usubinstr(owners, "  ", " ",.)
replace owners = usubinstr(owners, " ", "",.)
split owners, parse(,) generate(owner)
keep UCID_1934_2003 owner*
drop owners
sort UCID_1934_2003 owner*

*Deal with multiple identical IDs
generate id = _n
order id UCID_1934_2003
reshape long owner, i(UCID_1934_2003 id) j(test)
drop id test
drop if owner==""
drop if UCID_1934_2003==.
byso UCID_1934_2003 owner: gen dup = cond(_N==1,0,_n)
drop if dup>1
byso UCID_1934_2003: gen id=(_n)
	collapse (max) id, by(UCID_1934_2003)
	sort id
	rename id nb_owners

save "$dump\02_Processed_data\10_Directors_1934_2003\11_Ownerlists\nb_owner_`y'.dta", replace
}

***************************************************

***************************************************
*** Edgelist bipartite ***
***************************************************

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'.dta", clear

*** Create network of firms ***
replace owners = usubinstr(owners, ";", ",",.)
replace owners = usubinstr(owners, "   ", "  ",.)
replace owners = usubinstr(owners, "  ", " ",.)
replace owners = usubinstr(owners, " ", "",.)
split owners, parse(,) generate(owner)
keep UCID_1934_2003 owner*
drop owners
sort UCID_1934_2003 owner*
quietly by UCID_1934_2003 owner*: gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*Edge list
generate id = _n
order id UCID_1934_2003
reshape long owner, i(UCID_1934_2003 id) j(test)
drop id test
drop if owner==""
drop if UCID_1934_2003==.
rename owner id_sug 
destring id_sug, replace

tostring UCID_1934_2003, gen(UCIDstr)
tostring id_sug, gen(id_sugstr)
gen N_UCID = "C" + UCIDstr
gen N_id_sug = "P" + id_sugstr
order N_UCID N_id_sug

sort N_UCID N_id_sug
quietly by N_UCID N_id_sug: gen dup = cond(_N==1,0,_n)
drop if dup>1

keep N_UCID N_id_sug

save "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\edgeliste_`y'_bipartite.dta", replace
export excel "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\edgeliste_`y'_bipartite.xlsx", replace firstrow(variables)
export delimited using "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\edgeliste_`y'_bipartite.csv", delimiter(";") replace
}

***************************************************
*** List of persons (unique appearance) ***
***************************************************

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

*List of persons with unique entries
use "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person_CleanName-Gender-Geo.dta" , clear
keep if year == `y'
sort id_sug
quietly by id_sug: gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

tostring id_sug, gen(id_sugstr)
gen N_id_sug = "P" + id_sugstr

save "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_pers\pers_`y'_attribute_complete.dta", replace
export excel "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_pers\pers_`y'_attribute_complete.xlsx", replace firstrow(variables)
export delimited using "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\pers_`y'_attribute_complete.csv", delimiter(";") replace
}

***************************************************
*** List of firms (unique appearance) ***
***************************************************

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

*List of firms with unique entries & add N_UCID
use "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`y'.dta", clear

tostring UCID_1934_2003, gen(UCIDstr)
gen N_UCID = "C" + UCIDstr

save "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_`y'_attribute_complete.dta", replace
export excel "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_`y'_attribute_complete.xlsx", replace firstrow(variables)
export delimited using "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\comp_`y'_attribute_complete.csv", delimiter(";") replace
}


********************************************************
*** Building nodes attributes all nodes (bipartite)  ***
********************************************************

*** 1st persons attributes w/o duplicates***
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\edgeliste_`y'_bipartite.dta", clear

merge m:1 N_id_sug using "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_pers\pers_`y'_attribute_complete.dta"
keep if _merge == 3

quietly bysort N_id_sug: gen dup = cond(_N==1,0,_n)
drop if dup>1
drop N_UCID dup ID_dupl gde_options options4a options4b options_manual _merge id_sugstr namecorrection geo_merge countPLZ4 alternateaddress PLZ4
gen fullname = firstname + " " + lastname

save "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_pers\pers_attribute_`y'_deduplicates.dta", replace
export excel "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_pers\pers_attribute_`y'_deduplicates.xlsx", replace firstrow(variables)
export delimited using "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\pers_attribute_`y'_deduplicates.csv", delimiter(";") replace
}

*** 2nd companies attributes w/o duplicates***
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\edgeliste_`y'_bipartite.dta", clear

merge m:1 N_UCID using"$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_`y'_attribute_complete.dta"
keep if _merge == 3
drop UCIDstr caddress_orig gdename_orig _merge

quietly bysort N_UCID: gen dup2 = cond(_N==1,0,_n)
drop if dup2>1

save "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_attribute_`y'_deduplicates.dta", replace
export excel "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_attribute_`y'_deduplicates.xlsx", replace firstrow(variables)
export delimited using "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\comp_attribute_`y'_deduplicates.csv", delimiter(";") replace
}

*************************************************************
*** Projection of bipartite network into unimodal network ***
*************************************************************

/*
Done using Rstudio : "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\03_Code\10_Directors_1934_2003\R\..."

Files : "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\..."

*/

/*************************************************************
*** Building nodes attributes (Unimodal)                  ***
*************************************************************

*** Firms

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
import excel "$dump\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\firm_unimode_w_`y'.xlsx", sheet("firms") firstrow clear

drop A weight to
quietly bysort from: gen dup3 = cond(_N==1,0,_n)
drop if dup3>1
drop dup3

save "$dump\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\from_`y'.dta", replace

import excel "$dump\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\firm_unimode_w_`y'.xlsx", sheet("firms") firstrow clear

drop A weight from
quietly bysort to: gen dup3 = cond(_N==1,0,_n)
drop if dup3>1
drop dup3

append using "$dump\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\from_`y'.dta"
replace to = from if missing(to)
drop from
quietly bysort to: gen dup3 = cond(_N==1,0,_n)
drop if dup3>1
drop dup3 
rename to N_UCID

merge 1:1 N_UCID using"$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_`y'_attribute_complete.dta"
keep if _merge == 3

save "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_attribute_`y'_unimodal.dta", replace
export excel "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_attribute_`y'_unimodal.xlsx", replace firstrow(variables)
erase "$dump\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\from_`y'.dta"
}
*/

********************************************************************************
* Cleaning for Yinzuo
********************************************************************************

*Firms characteristics
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_`y'_attribute_complete.dta", clear 

order CID UCID_1934_2003 N_UCID 

save "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_`y'_attribute_complete_toshare.dta", replace
export delimited using "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\comp_`y'_attribute_complete_toshare.csv", delimiter(";") replace
}

*Directors characteristics
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_pers\pers_attribute_`y'_deduplicates.dta", clear 

preserve
keep N_id_sug PID year
save "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\pers_key_`y'.dta", replace
restore


drop titlejob firstname lastname firstname_orig lastname_orig year_sug fullname PID id_sug

save "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_pers\pers_`y'_attribute_complete_toshare.dta", replace
export delimited using "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\pers_`y'_attribute_complete_toshare.csv", delimiter(";") replace
}