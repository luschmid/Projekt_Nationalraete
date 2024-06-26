version 17
clear all
set more off
cap log close

*global path "C:\Users\yanni\Dropbox\Projekt Nationalräte\"
*global dump "C:\Users\yanni\Desktop\Projekt Nationalräte\"
global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte" /* Office */
global dump "C:\Users\schmutzy\Desktop\Projekt Nationalräte" /* Office */

***************************************************
*** Generating list of all owners for each firm ***
***************************************************

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", clear

*** Generate list of owners with single firm IDs ***
replace owners = usubinstr(owners, ";", ",",.)
replace owners = usubinstr(owners, "   ", "  ",.)
replace owners = usubinstr(owners, "  ", " ",.)
replace owners = usubinstr(owners, " ", "",.)
split owners, parse(,) generate(owner)
keep UCID_3 owner*
drop owners
sort UCID_3 owner*
quietly by UCID_3 owner*: gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*Deal with multiple identical IDs
generate id = _n
order id UCID_3
reshape long owner, i(UCID_3 id) j(test)
drop id test
drop if owner==""
drop if UCID_3==.
byso UCID_3: gen id=(_n)
destring owner, replace
reshape wide owner, i(UCID_3) j(id)

save "$dump\02_Processed_data\10_Directors_1934_2003\11_Ownerlists\ownerlist_`y'.dta", replace
}

******************************************************************
*Nb of owners deduplicated by UCID_3 : to link with cross-sections
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", clear

*** Generate list of owners with single firm IDs ***
replace owners = usubinstr(owners, ";", ",",.)
replace owners = usubinstr(owners, "   ", "  ",.)
replace owners = usubinstr(owners, "  ", " ",.)
replace owners = usubinstr(owners, " ", "",.)
split owners, parse(,) generate(owner)
keep UCID_3 owner*
drop owners
sort UCID_3 owner*

*Deal with multiple identical IDs
generate id = _n
order id UCID_3
reshape long owner, i(UCID_3 id) j(test)
drop id test
drop if owner==""
drop if UCID_3==.
byso UCID_3 owner: gen dup = cond(_N==1,0,_n)
drop if dup>1
byso UCID_3: gen id=(_n)
	collapse (max) id, by(UCID_3)
	sort id
	rename id nb_owners

save "$dump\02_Processed_data\10_Directors_1934_2003\11_Ownerlists\nb_owner_`y'.dta", replace
}

***************************************************

***************************************************
*** Edgelist bipartite ***
***************************************************

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", clear

*** Create network of firms ***
replace owners = usubinstr(owners, ";", ",",.)
replace owners = usubinstr(owners, "   ", "  ",.)
replace owners = usubinstr(owners, "  ", " ",.)
replace owners = usubinstr(owners, " ", "",.)
split owners, parse(,) generate(owner)
keep UCID_3 owner*
drop owners
sort UCID_3 owner*
quietly by UCID_3 owner*: gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*Edge list
generate id = _n
order id UCID_3
reshape long owner, i(UCID_3 id) j(test)
drop id test
drop if owner==""
drop if UCID_3==.
rename owner id_sug 
destring id_sug, replace

tostring UCID_3, gen(UCIDstr)
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
}

***************************************************
*** List of firms (unique appearance) ***
***************************************************

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

*List of firms with unique entries & add N_UCID
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", clear

tostring UCID_3, gen(UCIDstr)
gen N_UCID = "C" + UCIDstr

save "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_`y'_attribute_complete.dta", replace
export excel "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_`y'_attribute_complete.xlsx", replace firstrow(variables)
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
}

*** 2nd companies attributes w/o duplicates***
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\edgeliste_`y'_bipartite.dta", clear

merge m:1 N_UCID using"$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_`y'_attribute_complete.dta"
keep if _merge == 3
drop UCIDstr match_1 match_2 match_3 caddress_orig gdename_orig geo_merge _merge

quietly bysort N_UCID: gen dup = cond(_N==1,0,_n)
drop if dup>1

save "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_attribute_`y'_deduplicates.dta", replace
export excel "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_attribute_`y'_deduplicates.xlsx", replace firstrow(variables)
}

*************************************************************
*** Projection of bipartite network into unimodal network ***
*************************************************************

/*
Done using Rstudio : "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\03_Code\10_Directors_1934_2003\R\..."

Files : "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\..."

*/

*************************************************************
*** Building nodes attributes (Unimodal)                  ***
*************************************************************

*** Firms

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
import excel "$dump\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\firm_unimode_w_`y'.xlsx", sheet("firms") firstrow clear

drop A weight to
quietly bysort from: gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup 

save "$dump\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\from_`y'.dta", replace

import excel "$dump\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\firm_unimode_w_`y'.xlsx", sheet("firms") firstrow clear

drop A weight from
quietly bysort to: gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup 

append using "$dump\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\from_`y'.dta"
replace to = from if missing(to)
drop from
quietly bysort to: gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup 
rename to N_UCID

merge 1:1 N_UCID using"$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_`y'_attribute_complete.dta"
keep if _merge == 3

save "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_attribute_`y'_unimodal.dta", replace
export excel "$dump\02_Processed_data\10_Directors_1934_2003\10_Network_Bipartite\Node_attribute_comp\comp_attribute_`y'_unimodal.xlsx", replace firstrow(variables)
erase "$dump\02_Processed_data\10_Directors_1934_2003\14_Network_OneMode\from_`y'.dta"
}
