version 17
clear all
set more off
cap log close

*global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"
global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\" /* Office */
global dump "C:\Users\schmutzy\Desktop\Projekt Nationalräte\" /* Office Dumpster */
global ED "C:\Users\schmutzy\Dropbox\Travaux_EmilieDousse\" /* Office ED */
*global path "C:\Users\yanni\Dropbox\Projekt Nationalräte\" /* Home */
*global dump "C:\Users\yanni\Desktop\Projekt Nationalräte\" /* Home Dumpster */
*global ED "C:\Users\yanni\Dropbox\Travaux_EmilieDousse\" /* Home ED */

 
********************************************************************************
*** Append geocoded database from Sugarcube Data *******************************
********************************************************************************

use "$dump\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_complete.dta", clear //constructed in geocode dofile (near end)

********************************************************************************
*** String standardization ***
******************************

* Simple & obvious character corrections
foreach str in cname {
	replace `str' = usubinstr(`str', "ä", "ae",.)
	replace `str' = usubinstr(`str', "ö", "oe",.)
	replace `str' = usubinstr(`str', "ü", "ue",.)
	replace `str' = usubinstr(`str', "è", "e",.)
	replace `str' = usubinstr(`str', "é", "e",.)
	replace `str' = usubinstr(`str', "ê", "e",.)
	replace `str' = usubinstr(`str', "ë", "e",.)
	replace `str' = usubinstr(`str', "œ", "oe",.)
	replace `str' = usubinstr(`str', "à", "a",.)
	replace `str' = usubinstr(`str', "â", "a",.)
	replace `str' = usubinstr(`str', "û", "u",.)
	replace `str' = usubinstr(`str', "ù", "u",.)
	replace `str' = usubinstr(`str', "ô", "o",.)
	replace `str' = usubinstr(`str', "ò", "o",.)
	replace `str' = usubinstr(`str', "ó", "o",.)	
	replace `str' = usubinstr(`str', "î", "i",.)
	replace `str' = usubinstr(`str', "ì", "i",.)
	replace `str' = usubinstr(`str', "ï", "i",.)
	replace `str' = usubinstr(`str', "ç", "c",.)
	replace `str' = usubinstr(`str', "ß", "ss",.)
	replace `str' = usubinstr(`str', " und ", "&",.)
	replace `str' = usubinstr(`str', " et ", "&",.)
	replace `str' = usubinstr(`str', "& co.", "& co",.)
	replace `str' = usubinstr(`str', "et co.", "& co",.)
	replace `str' = usubinstr(`str', "und co.", "& co",.)
	replace `str' = usubinstr(`str', "& cie.", "& co",.)
	replace `str' = usubinstr(`str', "& cie", "& co",.)
	replace `str' = usubinstr(`str', "et cie", "& co",.)
	replace `str' = usubinstr(`str', "und cie", "& co",.)
	replace `str' = usubinstr(`str', "&cie", "& co",.)
	replace `str' = usubinstr(`str', "etcie", "& co",.)
	replace `str' = usubinstr(`str', "undcie", "& co",.)
	replace `str' = usubinstr(`str', "&fils", "& fils",.)
	replace `str' = usubinstr(`str', "etfils", "& fils",.)
	replace `str' = usubinstr(`str', "vormals", "vorm",.)
	replace `str' = usubinstr(`str', "$", "s",.)
	replace `str' = usubinstr(`str', "£", " ",.)
	replace `str' = usubinstr(`str', "^", " ",.)	
	replace `str' = usubinstr(`str', "~", " ",.)
	replace `str' = usubinstr(`str', "''", " ",.)
	replace `str' = usubinstr(`str', ",", " ",.)
	replace `str' = usubinstr(`str', "\", " ",.)
	replace `str' = usubinstr(`str', "/", " ",.)
	replace `str' = usubinstr(`str', "+", " ",.)
	replace `str' = usubinstr(`str', "%", " ",.)
	replace `str' = usubinstr(`str', "=", " ",.)
	replace `str' = usubinstr(`str', "[", " ",.)
	replace `str' = usubinstr(`str', "]", " ",.)
	replace `str' = usubinstr(`str', "{", " ",.)
	replace `str' = usubinstr(`str', "}", " ",.)
	replace `str' = usubinstr(`str', "^", " ",.)
	replace `str' = usubinstr(`str', ":", " ",.)
	replace `str' = usubinstr(`str', ";", " ",.)
	replace `str' = usubinstr(`str', "_", " ",.)
	replace `str' = usubinstr(`str', ",,", " ",.)
	replace `str' = usubinstr(`str', "*", " ",.)
	replace `str' = usubinstr(`str', "'", " ",.)
	replace `str' = usubinstr(`str', "succursale de", "",.)
	replace `str' = usubinstr(`str', "succursale", "",.)
	replace `str' = usubinstr(`str', "branch of", "",.)
	replace `str' = usubinstr(`str', "branch", "",.)
	replace `str' = usubinstr(`str', "zweigniederlassung", "",.)
	replace `str' = usubinstr(`str', "ag fuer", "ag",.)
	replace `str' = usubinstr(`str', "ag fur", "ag",.)
	replace `str' = usubinstr(`str', "ag.", "ag",.)
	replace `str' = usubinstr(`str', "aktiongesellschaft", "ag",.)
	replace `str' = usubinstr(`str', "aktiengesellschaft", "ag",.)
	replace `str' = usubinstr(`str', "aktiengeselischaft", "ag",.)
	replace `str' = usubinstr(`str', "(ag)", "ag",.)
	replace `str' = usubinstr(`str', "a.-g.", "ag",.)
	replace `str' = usubinstr(`str', "a.g.", "ag",.)
	replace `str' = usubinstr(`str', " sa ", " ",.)
	replace `str' = usubinstr(`str', "s.a.", "sa",.)
	replace `str' = usubinstr(`str', "(sa).", "sa",.)
	replace `str' = usubinstr(`str', "s. a. de l'immeuble", "s.i. ",.)
	replace `str' = usubinstr(`str', "s. a. immeuble", "s.i. ",.)
	replace `str' = usubinstr(`str', "s. a. des immeubles", "s.i. ",.)
	replace `str' = usubinstr(`str', "ste anonyme des immeubles", "s.i. ",.)
	replace `str' = usubinstr(`str', "ste anonyme", "sa",.)
	replace `str' = usubinstr(`str', "s. a.", "sa",.)
	replace `str' = usubinstr(`str', "s.a", "sa",.)
	replace `str' = usubinstr(`str', "sa.", "sa",.)
	replace `str' = usubinstr(`str', "schweiz.", "schweizerische",.)
	replace `str' = usubinstr(`str', "holding", "",.)
	replace `str' = usubinstr(`str', "vorm.", "vorm",.)
	replace `str' = usubinstr(`str', "compagnie", "",.)
	replace `str' = usubinstr(`str', "cie.", "& co",.)
	replace `str' = usubinstr(`str', "co.", "& co",.)
	replace `str' = usubinstr(`str', "-", " ",.)
	replace `str' = usubinstr(`str', "    ", " ",.)
	replace `str' = usubinstr(`str', "   ", " ",.)
	replace `str' = usubinstr(`str', "  ", " ",.)
	replace `str' = ustrtrim(`str')
	replace `str' = ustrlower(`str')
}

* Other corrections after looking at matching results

replace cname = "  "+cname+"  "

foreach str in cname {
	replace `str' = usubinstr(`str', "`", " ",.)
	replace `str' = usubinstr(`str', ">", " ",.)
	replace `str' = usubinstr(`str', "<", " ",.)
	replace `str' = usubinstr(`str', "«", " ",.)
	replace `str' = usubinstr(`str', "»", " ",.)
	replace `str' = usubinstr(`str', " u. dir ", " ",.)
	replace `str' = usubinstr(`str', " u. vr ", " ",.)
	replace `str' = usubinstr(`str', " et dir ", " ",.)
	replace `str' = usubinstr(`str', " dir ", " ",.)
	replace `str' = usubinstr(`str', " u. del ", " ",.)
	replace `str' = usubinstr(`str', " et del ", " ",.)
	replace `str' = usubinstr(`str', " u. ", " ",.)
	replace `str' = usubinstr(`str', " pr ", " ",.)
	replace `str' = usubinstr(`str', " a h ", " ",.)
	replace `str' = usubinstr(`str', " ag fuer ", "ag",.)
	replace `str' = usubinstr(`str', " fuer ", " ",.)
	replace `str' = usubinstr(`str', " fur ", " ",.)
	replace `str' = usubinstr(`str', " ges. fr ", " ges. ",.)
	replace `str' = usubinstr(`str', " a g. ", "ag",.)
	replace `str' = usubinstr(`str', " a. g ", "ag",.)
	replace `str' = usubinstr(`str', " a. g. ", "ag",.)
	replace `str' = usubinstr(`str', " s. a ", "sa",.)
	replace `str' = usubinstr(`str', "&", " & ",.)
	replace `str' = usubinstr(`str', " u ", " & ",.)
	replace `str' = usubinstr(`str', " and ", " & ",.)
	replace `str' = usubinstr(`str', " e ci ", " &co ",.)
	replace `str' = usubinstr(`str', " & co ", " &co",.)
	replace `str' = usubinstr(`str', " allgemeine ", " allg ",.)
	replace `str' = usubinstr(`str', " allgem ", " allg ",.)
	replace `str' = usubinstr(`str', " avenue ", " av ",.)
	replace `str' = usubinstr(`str', " ave ", " av ",.)
	replace `str' = usubinstr(`str', " gesellschaft ", " ges. ",.)
	replace `str' = usubinstr(`str', "    ", " ",.)
	replace `str' = usubinstr(`str', "   ", " ",.)
	replace `str' = usubinstr(`str', "  ", " ",.)
	replace `str' = ustrtrim(`str')
	replace `str' = ustrlower(`str')
}

********************************************************************************
*** Flagging special firms with extremely similar names (i.e. SI)    ***********
********************************************************************************
*Flag SI by exploiting structure of name (ending with single letter or number)

*Standardization of expression for S.I.
foreach str in cname {
	replace `str' = usubinstr(`str', "5.1.", "s.i. ",.)
	replace `str' = usubinstr(`str', "soc. immobiliere", "s.i. ",.)
	replace `str' = usubinstr(`str', "ste immobiliere", "s.i. ",.)
	replace `str' = usubinstr(`str', "sa de l'immeuble", "s.i. ",.)
	replace `str' = usubinstr(`str', "ste de l immeuble", "s.i. ",.)
	replace `str' = usubinstr(`str', "soc. de l'immeuble", "s.i. ",.)
	replace `str' = usubinstr(`str', "societe immob.", "s.i. ",.)
	replace `str' = usubinstr(`str', "soc. immob.", "s.i. ",.)
	replace `str' = usubinstr(`str', "ste immobiliere", "s.i. ",.)
	*replace `str' = usubinstr(`str', "s.i. ", "si ",.)
	replace `str' = usubinstr(`str', "s.a. de l'immeuble", "s.i. ",.)
	replace `str' = usubinstr(`str', "ste immeuble", "s.i. ",.)
	replace `str' = usubinstr(`str', "societe de l'immeuble", "s.i. ",.)	
	replace `str' = usubinstr(`str', "    ", " ",.)
	replace `str' = usubinstr(`str', "   ", " ",.)
	replace `str' = usubinstr(`str', "  ", " ",.)
	replace `str' = ustrtrim(`str')
	replace `str' = ustrlower(`str')
} 

*Strict way of flagging SI's
g cname_2 = cname +" " // adding of a space at the end to have " a " structure


*Flag SI by simply searching for "si " or "s.i. " 
g si_dummy_2 = 0
order CID year cname si_dummy_2
replace si_dummy_2 = 1 if regexm(cname_2, "s\.i\. ")
replace si_dummy_2 = 1 if regexm(cname_2, "s\. i\. ")
replace si_dummy_2 = 1 if regexm(cname_2, "s i\. ")
replace si_dummy_2 = 1 if regexm(cname_2, "s i ")
replace si_dummy_2 = 1 if regexm(cname_orig, "SI ")
replace si_dummy_2 = 1 if regexm(cname_orig, "S.1. ")
*Cases with only s and not si? Use cantons probably
replace si_dummy_2 = 1 if regexm(cname_orig, "^S ") & ctn == "ge"
replace si_dummy_2 = 1 if regexm(cname_orig, "^S ") & ctn == "vd"
replace si_dummy_2 = 1 if regexm(cname_orig, "^S ") & ctn == "fr"
replace si_dummy_2 = 1 if regexm(cname_orig, "^S ") & ctn == "vs"
replace si_dummy_2 = 1 if regexm(cname_orig, "^S ") & ctn == "ne"


* Company name cleaning
replace cname = usubinstr(cname, ".", "",.)
replace cname = usubinstr(cname, " no ", " ",.)
replace cname = usubinstr(cname, " nr ", " ",.)
replace cname = usubinstr(cname, "  ", " ",.)
replace cname = ustrtrim(cname)

foreach str in cname {
	replace `str' = usubinstr(`str', "societe ", "ste ",.)
	replace `str' = usubinstr(`str', "soc ", "ste ",.)
	replace `str' = usubinstr(`str', "!", "",.)
}

*save "$dump\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_complete.dta", replace 
*use "$dump\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_complete.dta", clear 

* separate years (to work with lighter files)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	preserve
		keep if year == `y'
		save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_temp.dta", replace
	restore	
}

********************************************************************************
* Matching and deduplicates "Flow" Approach
********************************************************************************

********************************************************************************
*** 1st matching : perfect matching	
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_temp.dta", clear
	*attributing company ID if cnames are exactly the same : 1st matching threshold
	*preserve
		sort gdenr cname 
		bysort gdenr : strgroup cname, gen(UCID_1) thresh(0.00000000000001) force
		order CID UCID_1 year
		duplicates tag UCID_1, gen(match_1)
		local i = 1000
		while `i' > 1 {
			local prev = `i' - 1
			replace match_1 = `i' if match_1 == `prev'
			local i = `prev'
		}
		save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_1.dta", replace
	*restore
}

*collapse multiple perfect matches in first round
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_1.dta", clear

*Build unique CID_str for all UCID_1
gen CID_str = string(CID)
order CID CID_str
sort UCID_1 cname

bysort UCID_1: replace CID_str = CID_str[_n-1] + "; " + CID_str if match_1 > 0
replace CID_str = " "+CID_str
replace CID_str = usubinstr(CID_str, " ; ", "",.)
by UCID_1: replace CID_str = CID_str[_N]

*Same with owners
bysort UCID_1: replace owners = owners[_n-1] + "; " + owners if match_1 > 0
replace owners = " "+owners
replace owners = usubinstr(owners, " ; ", "",.)
by UCID_1: replace owners = owners[_N]

*Same with function
bysort UCID_1: replace function = function[_n-1] + "; " + function if match_1 > 0
replace function = " "+function
replace function = usubinstr(function, " ; ", "",.)
by UCID_1: replace function = function[_N]

*Same with signature
bysort UCID_1: replace signature = signature[_n-1] + "; " + signature if match_1 > 0
replace signature = " "+signature
replace signature = usubinstr(signature, " ; ", "",.)
by UCID_1: replace signature = signature[_N]

*Same with page
bysort UCID_1: replace page = page[_n-1] + "; " + page if match_1 > 0
replace page = " "+page
replace page = usubinstr(page, " ; ", "",.)
by UCID_1: replace page = page[_N]

*Delete duplicates
by UCID_1:  gen dup = cond(_N==1,0,_n)
drop if dup > 1
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_1.dta", replace
}

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
preserve
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_1.dta", clear
	keep if match_1 > 0
	sort UCID_1 year
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_1_FP.dta", replace
	*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_temp.dta"
restore
}

********************************************************************************

********************************************************************************
*** 2nd matching: a) create subsample with SIs 
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_1.dta", clear
preserve
	keep if si_dummy_2 == 1
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_SI.dta", replace
restore

*2nd matching: b) rest (10%)		
preserve
	drop if si_dummy_2 == 1
	bysort gdenr: strgroup cname, gen(UCID_nonSI) thresh(0.1) force
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_nonSI.dta", replace
restore
}

*a) Match SIs according to procedure : 
*UCID_SI_a: same ID if exact same owners; 
*UCID_SI_b: same ID if same UCID_SI_a, and cname 20% different max; collapsed on UCID_SI_b
*UCID_SI_c: unique ID if cname 10% different max (UCID_SI_d) AND we extract the same single number or letter; collapsed as well

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_SI.dta", clear
	
	***Same owners within matched UCID_SI	
		replace owners = usubinstr(owners, "; ", ",",.)
		strgroup owners, gen(UCID_SI_a) thresh(0.00000000000001) force
		bysort UCID_SI_a: strgroup cname, gen(UCID_SI_b) thresh(0.2) force
			drop dup
			duplicates tag UCID_SI_b, gen(dup)
			sort UCID_SI_a UCID_SI_b cname
		
	*Match 10% to create UCID_SI_c	
		bysort gdenr: strgroup cname, gen(UCID_SI_d) thresh(0.10) force
		order CID year UCID_1 UCID_SI_a UCID_SI_b UCID_SI_d
		
	***collapse multiple perfect matches in first round
		*Build unique CID_str for all UCID_SI_b (perfect match on owners & 20% cname)
		order CID CID_str
		sort UCID_SI_b cname
		bysort UCID_SI_b: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_SI_b: replace CID_str = CID_str[_N]

		*Delete duplicates
		drop dup
		by UCID_SI_b:  gen dup = cond(_N==1,0,_n)
		drop if dup > 1
	
	drop cname_2
	g cname_2 = cname +" " // adding of a space at the end to have " a " structure
	gen SI_extract = ""
	foreach i in l a b c d e f g h i j k m n o p q r s t u v w x y z{
	replace SI_extract = "`i'" if (regexm(cname_2, " `i' "))
	}
	forvalues i = 1(1)999{
	replace SI_extract = "`i'" if (regexm(cname_2, " `i' "))
	}
	replace SI_extract = cname if missing(SI_extract)

	bysort UCID_SI_d: strgroup SI_extract, gen(UCID_SI_c) thresh(0.00000000000001) force
	order CID year UCID_1 UCID_SI_a UCID_SI_b SI_extract UCID_SI_c UCID_SI_d
	
		***collapse multiple perfect matches in first round
		*Build unique CID_str for all UCID_SI_c (perfect match on owners & 20% cname)
		order CID CID_str
		sort UCID_SI_c cname
		bysort UCID_SI_c: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_SI_c: replace CID_str = CID_str[_N]
		
		*Same with owners
		bysort UCID_SI_c: replace owners = owners[_n-1] + "; " + owners 
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		by UCID_SI_c: replace owners = owners[_N]

		*Same with function
		bysort UCID_SI_c: replace function = function[_n-1] + "; " + function 
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_SI_c: replace function = function[_N]

		*Same with signature
		bysort UCID_SI_c: replace signature = signature[_n-1] + "; " + signature
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_SI_c: replace signature = signature[_N]

		*Same with page
		bysort UCID_SI_c: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_SI_c: replace page = page[_N]

		*Delete duplicates
		drop dup
		by UCID_SI_c:  gen dup = cond(_N==1,0,_n)
		drop if dup > 1
	
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_SI_a.dta", replace
}

*Put collapsed SIs back in the match_2 dataset
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_nonSI.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_SI_a.dta" // Here we have two different UCID that we need to merge into 1...
order CID UCID_1 UCID_nonSI UCID_SI_c year cname

*Constructing UCID_2
gen UCID_2 = UCID_nonSI if si_dummy_2 == 0
egen maxUCID_nonSI=max(UCID_nonSI)
replace UCID_2 = UCID_SI_c + maxUCID_nonSI if si_dummy_2 == 1
drop UCID_nonSI UCID_SI_a UCID_SI_b UCID_SI_c 

*Tag newly matched obs (compared to match 1)
duplicates tag UCID_2, gen(match_2)
local i = 1000
		while `i' > 1 {
			local prev = `i' - 1
			replace match_2 = `i' if match_2 == `prev'
			local i = `prev'
		}
order CID UCID_1 UCID_2
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2.dta", replace
}

*Extract new matches tagged just before
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2.dta", clear

preserve
	gen match_2_temp = match_2 + match_1
	drop match_2
	rename match_2_temp match_2
	keep if match_2 > match_1
	sort UCID_2 year
	order CID year UCID_1 UCID_2 match_1 match_2 cname gdename owners capital function signature page 
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP.dta", replace
restore

preserve
	gen match_2_temp = match_2 + match_1
	drop match_2
	rename match_2_temp match_2
	drop if match_2 > match_1
	sort UCID_2 year
	order CID year UCID_1 UCID_2 match_1 match_2 cname gdename owners capital function signature page 
	drop UCID_2
	gen UCID_2_temp1 = _n
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_nonFP.dta", replace
restore
}

*Subsamples of big matches : suspiciously high number of obs matched together - require evaluation & subsample of rest
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP.dta", clear
	*subsamples of 3+ obs matched together
	preserve
		drop dup 
		duplicates tag UCID_2, gen(dup)
		replace dup = dup + 1 // change from nb of dup to nb of UCID_2 appearing
		drop if dup < 3 // drop dup = 1 and 2
		***Same owners within matched : we want to get rid of potential SIs not identified. So we identify firms with exact same owners and company names 20% different... 
		replace owners = usubinstr(owners, "; ", ",",.)
		strgroup owners, gen(UCID_SI_a) thresh(0.00000000000001) force
		bysort UCID_SI_a: strgroup cname, gen(UCID_SI_b) thresh(0.2) force
			drop dup 
			duplicates tag UCID_SI_b, gen(dup)
			sort UCID_SI_a UCID_SI_b cname
			order CID year UCID_1 UCID_SI_a UCID_SI_b 
		
		***collapse on UCID_SI_b : get rid of SIs with same owners
		*Build unique CID_str for all UCID_SI_b (perfect match on owners & 20% cname)
		order CID CID_str
		sort UCID_SI_b cname
		bysort UCID_SI_b: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_SI_b: replace CID_str = CID_str[_N]
		
		*Same with cname (create cname group)
		gen cname_group = cname
		bysort UCID_SI_b: replace cname_group = cname_group[_n-1] + "; " + cname_group 
		replace cname_group = " "+cname_group
		replace cname_group = usubinstr(cname_group, " ; ", "",.)
		replace cname_group = ustrtrim(cname_group)
		by UCID_SI_b: replace cname_group = cname_group[_N]

		*Same with function
		bysort UCID_SI_b: replace function = function[_n-1] + "; " + function 
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_SI_b: replace function = function[_N]

		*Same with signature
		bysort UCID_SI_b: replace signature = signature[_n-1] + "; " + signature 
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_SI_b: replace signature = signature[_N]

		*Same with page
		bysort UCID_SI_b: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_SI_b: replace page = page[_N]

		*Delete duplicates
		drop dup
		by UCID_SI_b:  gen dup = cond(_N==1,0,_n)
		drop if dup > 1
		sort UCID_2 cname
		save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_toeval.dta", replace
	restore

	*Subsample of rest : ok
	preserve	
		drop dup 
		duplicates tag UCID_2, gen(dup)
		replace dup = dup + 1
		drop if dup > 2 // drop dup 3, 4, 5, etc.
		
		gen cname_group = cname
		bysort UCID_2: replace cname_group = cname_group[_n-1] + "; " + cname_group 
		replace cname_group = " "+cname_group
		replace cname_group = usubinstr(cname_group, " ; ", "",.)
		replace cname_group = ustrtrim(cname_group)
		by UCID_2: replace cname_group = cname_group[_N]
		
		*Collapse on UCID_2 : those obs are assumed close enough
		order CID CID_str
		sort UCID_2 cname
		bysort UCID_2: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_2: replace CID_str = CID_str[_N]

		*Same with owners
		bysort UCID_2: replace owners = owners[_n-1] + "; " + owners 
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		by UCID_2: replace owners = owners[_N]

		*Same with function
		bysort UCID_2: replace function = function[_n-1] + "; " + function 
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_2: replace function = function[_N]

		*Same with signature
		bysort UCID_2: replace signature = signature[_n-1] + "; " + signature 
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_2: replace signature = signature[_N]

		*Same with page
		bysort UCID_2: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_2: replace page = page[_N]

		*Delete duplicates
		drop dup
		by UCID_2:  gen dup = cond(_N==1,0,_n)
		drop if dup > 1
		save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_ok.dta", replace
	restore
}

/*Export in excel sorted by UCID_2 cname

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
preserve
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_toeval.dta", clear
	sort UCID_2 cname	
	drop UCID_1 UCID_SI_a UCID_SI_b capital function signature page cname_orig caddress_orig gdename_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract E_CNTR N_CNTR cname_2 SI_extract maxUCID_nonSI dup CID_str match_1 match_2 si_dummy_2 gdenr
	export excel using "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_test.xlsx", replace firstrow(variables)
restore
}*/


/********************************************************************************
* Corrections of Emilie and Julie : identifying differences in coding FPs

*FP evaluations comparison: only single obs false positive or max one group false positive in a UCID_2 group
foreach y in 1934 1943 1960 /*1962 1963*/ 1964 1965 1966 1969 1972 1975 1979 /*1980 1981 1982 1983*/ 1984 /*1985*/ 1986 1987 /*1988 1989 */ 1990 1991 1992 1993 /*1994 1995 1996*/ 1997 1998 1999 /*2000 2001 2002 2003*/ {
import excel "$ED\FP check\Emilie\FP2_check_`y'", sheet("Sheet1") firstrow clear
rename (FP_single FP_group) (FP_single_ED FP_group_ED)
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_ED.dta", replace

import excel "$ED\FP check\Julie\FP2_check_`y'", sheet("Sheet1") firstrow clear
rename (FP_single FP_group) (FP_single_JU FP_group_JU)
merge 1:1 year CID using "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_ED.dta"
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_ED_JU.dta", replace

gen diff_check = 0 if FP_single_ED == FP_single_JU & FP_group_ED == FP_group_JU
replace diff_check = 1 if missing(diff_check)
levelsof UCID_2 if diff_check == 1, local(diff)
gen diff_check_2 = 0
foreach l in `diff' {
 	replace diff_check_2 = 1 if UCID_2 == `l' 
}
 
preserve
keep if diff_check_2 == 1
sort UCID_2 cname
drop CID _merge diff_check_2
gen FP_single = 0
gen FP_group = 0
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_diff.dta", replace
export excel using "$ED\FP check\Final\FP2_check_`y'.xlsx", replace firstrow(variables)
restore

preserve
drop if diff_check_2 == 1
gen FP_single = FP_single_ED*FP_single_JU
gen FP_group = FP_group_ED*FP_group_JU
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_same.dta", replace
restore
}

*FP evaluations comparison : with several groups that should be splitted in one UCID_2 group 
foreach y in /*1934 1943 1960*/ 1962 1963 /*1964 1965 1966 1969 1972 1975 1979*/ 1980 1981 1982 1983 /*1984*/ 1985 /*1986 1987*/ 1988 1989 /*1990 1991 1992 1993*/ 1994 1995 1996 /*1997 1998 1999*/ 2000 2001 2002 2003{
import excel "$ED\FP check\Emilie\FP2_check_`y'", sheet("Sheet1") firstrow clear
rename (FP_single FP_group FP_group_2) (FP_single_ED FP_group_ED FP_group_2_ED)
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_ED.dta", replace

import excel "$ED\FP check\Julie\FP2_check_`y'", sheet("Sheet1") firstrow clear
rename (FP_single FP_group FP_group_2) (FP_single_JU FP_group_JU FP_group_2_JU)
merge 1:1 year CID using "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_ED.dta"
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_ED_JU.dta", replace

gen diff_check = 0 if FP_single_ED == FP_single_JU & FP_group_ED == FP_group_JU & FP_group_2_ED == FP_group_2_JU
replace diff_check = 1 if missing(diff_check)
levelsof UCID_2 if diff_check == 1, local(diff)
gen diff_check_2 = 0
foreach l in `diff' {
 	replace diff_check_2 = 1 if UCID_2 == `l' 
 }
 
preserve
keep if diff_check_2 == 1
sort UCID_2 cname
drop CID _merge diff_check_2
gen FP_single = 0
gen FP_group = 0
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_diff.dta", replace
export excel using "$ED\FP check\Final\FP2_check_`y'.xlsx", replace firstrow(variables)
restore

preserve
drop if diff_check_2 == 1
gen FP_single = FP_single_ED*FP_single_JU //Can be 0 or 1, but they had the same evaluation so it generates only 0*0 or 1*1
gen FP_group = FP_group_ED*FP_group_JU //0*0 or 1*1 again
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_same.dta", replace
restore
}

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
erase "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_ED.dta"
erase "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP2_check_`y'_ED_JU.dta"
}*/


********************************************************************************
**Years without FP_group cases: in "...\02_Processed_data\10_Directors_1934_2003\17_FP_Check\Final\FP2_check_`y'_final.xlsx"!
foreach y in 1934 1943 /*1960 1962 1963 1964 1965 1966 1969 1972*/ 1975 /*1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000*/ 2001 2002 2003 {

***
*Collapse in FP_ok file: those were correct matches (subsample: we assumed ok)
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_ok.dta", clear
*Gen ID by obs
gen UCID_2_temp2 = _n
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_1.dta", replace

***
*FP_toeval
*Add corrections and post processing : used ED & JU xls, had a look at UCID_2 listed above in diff.dta, made decisions... Basically "FP_toeval" in xlsx.
import excel "$ED\FP check\Final\FP2_check_`y'_final.xlsx", sheet("Sheet1") firstrow clear

**Collapse those who were right : to_eval and were evaluated as correct by JU & ED
*Count number of commas to select "best name": gen comma = length(trim(var1))-length(trim(subinstr(nb_owners,",","",.))) where nb_owners is the count minus 1 (replace nb_owners = nb_owners + 1)
preserve
	keep if FP_single ==. & FP_group ==.
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2.dta", keep(3) nogen

	***collapse multiple perfect matches in first round
		*Build unique CID_str for all UCID_SI_b (perfect match on owners & 20% cname)
		*Count number of commas
		gen comma = (length(trim(owners))-length(trim(subinstr(owners,",","",.))))
		sort UCID_2 comma
		order CID CID_str comma
				
		bysort UCID_2: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_2: replace CID_str = CID_str[_N]
		
		*Same with cname (create cname group)
		gen cname_group = cname
		bysort UCID_2: replace cname_group = cname_group[_n-1] + "; " + cname_group 
		replace cname_group = " "+cname_group
		replace cname_group = usubinstr(cname_group, " ; ", "",.)
		replace cname_group = ustrtrim(cname_group)
		by UCID_2: replace cname_group = cname_group[_N]

		*Same with owners
		bysort UCID_2: replace owners = owners[_n-1] + "; " + owners
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		by UCID_2: replace owners = owners[_N]

		*Same with function
		bysort UCID_2: replace function = function[_n-1] + "; " + function
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_2: replace function = function[_N]

		*Same with signature
		bysort UCID_2: replace signature = signature[_n-1] + "; " + signature
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_2: replace signature = signature[_N]

		*Same with page
		bysort UCID_2: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_2: replace page = page[_N]

		*Delete duplicates
		drop dup
		by UCID_2:  gen dup = cond(_N==1,0,_n)
		bysort UCID_2 : keep if dup == dup[_N]
		sort UCID_2 cname
	drop UCID_2
	gen UCID_2_temp3 = _n	
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_2.dta", replace
restore

*Correction of FP_single who were wrong (give unique ID to those obs)
preserve
keep if FP_single==1 
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2.dta"
	keep if _merge ==3
	drop _merge
	gen UCID_2_temp4 = _n
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_3.dta", replace
restore

*Putting everything back together in "match_2_done"
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_nonFP.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_1.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_2.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_3.dta"
sort UCID_2_temp1 UCID_2_temp2 UCID_2_temp3 UCID_2_temp4
order UCID_2_temp1 UCID_2_temp2 UCID_2_temp3 UCID_2_temp4

* Create one single UCID_2

egen maxUCID_temp3=max(UCID_2_temp3)
replace UCID_2_temp3 = maxUCID_temp3 + UCID_2_temp4 if missing(UCID_2_temp3)
egen maxUCID_temp2=max(UCID_2_temp2)
replace UCID_2_temp2 = maxUCID_temp2 + UCID_2_temp3 if missing(UCID_2_temp2)
egen maxUCID_temp1=max(UCID_2_temp1)
replace UCID_2_temp1 = maxUCID_temp1 + UCID_2_temp2 if missing(UCID_2_temp1)

drop UCID_2_temp2 UCID_2_temp3 UCID_2_temp4 maxUCID_temp1 maxUCID_temp2 maxUCID_temp3 UCID_1 cname_2 dup SI_extract UCID_SI_d maxUCID_nonSI comma FP_single FP_group UCID_2
rename UCID_2_temp1 UCID_2 
sort gdenr_2018 cname
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_clean.dta", replace

erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_1.dta"
erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_2.dta"
erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_3.dta"
}


*****
*Years with FP_group cases (group_1 only): watch out which years have group dummy = 1 or not!

foreach y in /*1934 1943*/ 1960 1962 1963 1964 1965 1966 1969 1972 /*1975*/ 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 /*2001 2002 2003*/ {
*Collapse in FP_ok file
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_ok.dta", clear
	***collapse multiple perfect matches in first round
		*Build unique CID_str for all UCID_2 (ok subsample: we assumed ok)
		*Count number of commas
		gen comma = (length(trim(owners))-length(trim(subinstr(owners,",","",.))))
		sort UCID_2 comma
		order CID CID_str comma
				
		bysort UCID_2: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_2: replace CID_str = CID_str[_N]
		
		bysort UCID_2: replace cname_group = cname_group[_n-1] + "; " + cname_group 
		replace cname_group = " "+cname_group
		replace cname_group = usubinstr(cname_group, " ; ", "",.)
		replace cname_group = ustrtrim(cname_group)
		by UCID_2: replace cname_group = cname_group[_N]

		*Same with owners
		bysort UCID_2: replace owners = owners[_n-1] + "; " + owners
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		by UCID_2: replace owners = owners[_N]

		*Same with function
		bysort UCID_2: replace function = function[_n-1] + "; " + function
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_2: replace function = function[_N]

		*Same with signature
		bysort UCID_2: replace signature = signature[_n-1] + "; " + signature
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_2: replace signature = signature[_N]

		*Same with page
		bysort UCID_2: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_2: replace page = page[_N]

		*Delete duplicates
		drop dup
		by UCID_2:  gen dup = cond(_N==1,0,_n)
		drop if dup > 1
		sort UCID_2 cname
		drop UCID_2
		
		*Gen ID by obs
		gen UCID_2_temp2 = _n
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_1.dta", replace

*Add corrections and post processing
import excel "$ED\FP check\Final\FP2_check_`y'_final.xlsx", sheet("Sheet1") firstrow clear

**Collapse those who are right
*Count number of commas to select "best name": gen comma = length(trim(var1))-length(trim(subinstr(nb_owners,",","",.))) where nb_owners is the count minus 1 (replace nb_owners = nb_owners + 1)
preserve
	keep if FP_single ==. & FP_group ==.
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2.dta", keep(3) nogen
	
	***collapse multiple perfect matches in first round
		*Build unique CID_str for all UCID_SI_b (perfect match on owners & 20% cname)
		*Count number of commas
		gen comma = (length(trim(owners))-length(trim(subinstr(owners,",","",.))))
		sort UCID_2 comma
		order CID CID_str comma
				
		bysort UCID_2: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_2: replace CID_str = CID_str[_N]

		*Same with cname (create cname group)
		gen cname_group = cname
		bysort UCID_2: replace cname_group = cname_group[_n-1] + "; " + cname_group 
		replace cname_group = " "+cname_group
		replace cname_group = usubinstr(cname_group, " ; ", "",.)
		replace cname_group = ustrtrim(cname_group)
		by UCID_2: replace cname_group = cname_group[_N]
		
		*Same with owners
		bysort UCID_2: replace owners = owners[_n-1] + "; " + owners
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		by UCID_2: replace owners = owners[_N]

		*Same with function
		bysort UCID_2: replace function = function[_n-1] + "; " + function
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_2: replace function = function[_N]

		*Same with signature
		bysort UCID_2: replace signature = signature[_n-1] + "; " + signature
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_2: replace signature = signature[_N]

		*Same with page
		bysort UCID_2: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_2: replace page = page[_N]

		*Delete duplicates
		drop dup
		by UCID_2:  gen dup = cond(_N==1,0,_n)
		bysort UCID_2 : keep if dup == dup[_N]
		sort UCID_2 cname
	drop UCID_2
	gen UCID_2_temp3 = _n	
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_2.dta", replace
restore

*Correct FP_single who were wrong (give unique ID to those obs)
preserve
keep if FP_single==1 
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2.dta"
	keep if _merge ==3
	drop _merge
	gen UCID_2_temp4 = _n
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_3.dta", replace
restore

*Correct FP_group who were wrong (give unique ID to those group of obs) & collapse them
preserve
keep if FP_group==1 
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2.dta"
	keep if _merge ==3
	drop _merge
	
***collapse multiple perfect matches in first round
		*Build unique CID_str for all UCID_SI_b (perfect match on owners & 20% cname)
		*Count number of commas
		gen comma = (length(trim(owners))-length(trim(subinstr(owners,",","",.))))
		sort UCID_2 comma
		order CID CID_str comma
				
		bysort UCID_2: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_2: replace CID_str = CID_str[_N]
		
		*Same with cname (create cname group)
		gen cname_group = cname
		bysort UCID_2: replace cname_group = cname_group[_n-1] + "; " + cname_group 
		replace cname_group = " "+cname_group
		replace cname_group = usubinstr(cname_group, " ; ", "",.)
		replace cname_group = ustrtrim(cname_group)
		by UCID_2: replace cname_group = cname_group[_N]

		*Same with owners
		bysort UCID_2: replace owners = owners[_n-1] + "; " + owners
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		by UCID_2: replace owners = owners[_N]

		*Same with function
		bysort UCID_2: replace function = function[_n-1] + "; " + function
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_2: replace function = function[_N]

		*Same with signature
		bysort UCID_2: replace signature = signature[_n-1] + "; " + signature
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_2: replace signature = signature[_N]

		*Same with page
		bysort UCID_2: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_2: replace page = page[_N]

		*Delete duplicates
		drop dup
		by UCID_2:  gen dup = cond(_N==1,0,_n)
		bysort UCID_2 : keep if dup == dup[_N]
		sort UCID_2 cname
	drop UCID_2	
	
	gen UCID_2_temp5 = _n
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_4.dta", replace
restore

*Putting everything back together in "match_2_done"
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_nonFP.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_1.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_2.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_3.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_4.dta" 
sort UCID_2_temp1 UCID_2_temp2 UCID_2_temp3 UCID_2_temp4
order UCID_2_temp1 UCID_2_temp2 UCID_2_temp3 UCID_2_temp4

* Create one single UCID_2

egen maxUCID_temp4=max(UCID_2_temp4)
replace UCID_2_temp4 = maxUCID_temp4 + UCID_2_temp5 if missing(UCID_2_temp4)
egen maxUCID_temp3=max(UCID_2_temp3)
replace UCID_2_temp3 = maxUCID_temp3 + UCID_2_temp4 if missing(UCID_2_temp3)
egen maxUCID_temp2=max(UCID_2_temp2)
replace UCID_2_temp2 = maxUCID_temp2 + UCID_2_temp3 if missing(UCID_2_temp2)
egen maxUCID_temp1=max(UCID_2_temp1)
replace UCID_2_temp1 = maxUCID_temp1 + UCID_2_temp2 if missing(UCID_2_temp1)

drop UCID_2_temp2 UCID_2_temp3 UCID_2_temp4 UCID_2_temp5 maxUCID_temp1 maxUCID_temp2 maxUCID_temp3 UCID_1 cname_2 dup SI_extract UCID_SI_d maxUCID_nonSI comma FP_single FP_group UCID_2
rename UCID_2_temp1 UCID_2 
sort gdenr_2018 cname
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_clean.dta", replace

erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_1.dta"
erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_2.dta"
erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_3.dta"
erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_check_4.dta" 
}

********************************************************************************
*** 3rd matching: 20% threshold, after collapsing perfect and 10% matches in step match 1 & 2
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_clean.dta", clear

*Drop SI obs
preserve
	keep if si_dummy_2 == 1
	gen UCID_SI = _n
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_SI.dta", replace
restore

*Deduplicates w/o SIs
preserve
	drop if si_dummy_2 == 1
	bysort gdenr: strgroup cname, gen(UCID_nonSI) thresh(0.2) force
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_nonSI.dta", replace
restore

*Put the two subsamples together
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_nonSI.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_SI.dta" // Here we have two different UCID that we need to merge into 1...

*Constructing UCID_3
gen UCID_3 = UCID_nonSI if si_dummy_2 == 0
egen maxUCID_nonSI=max(UCID_nonSI)
replace UCID_3 = UCID_SI + maxUCID_nonSI if si_dummy_2 == 1
drop UCID_nonSI UCID_SI 

*Subsample of newly matched obs (compared to match 2)
duplicates tag UCID_3, gen(match_3)
local i = 500
		while `i' > 1 {
			local prev = `i' - 1
			replace match_1 = `i' if match_1 == `prev'
			local i = `prev'
		}
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3.dta", replace
}

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
*subsample of new matches compared to 10% matching 
preserve
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3.dta", clear
	gen match_3_temp = match_3 + match_2
	drop match_3
	rename match_3_temp match_3
	keep if match_3 > match_2
	sort UCID_3 year
	order CID year UCID_2 UCID_3 match_1 match_2 match_3 cname gdename owners capital function signature page 
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP.dta", replace
restore

*subsample of unmatched obs in 20% matching step
preserve
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3.dta", clear
	gen match_3_temp = match_3 + match_2
	drop match_3
	rename match_3_temp match_3
	drop if match_3 > match_2
	sort UCID_3 year
	order CID year UCID_2 UCID_3 match_1 match_2 match_3 cname gdename owners capital function signature page 
	gen UCID_3_temp1 = _n
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_nonFP.dta", replace
restore
}

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
preserve
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP.dta", clear
	sort UCID_3 cname	
	drop UCID_2 function signature page CID_str si_dummy_2 cname_orig caddress_orig gdename_orig geo_merge branch_dummy gdename_extract capital_extract cname_extract E_CNTR N_CNTR maxUCID_nonSI
	export excel using "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP3_check_`y'.xlsx", replace firstrow(variables)
restore
}

/* Corrections of Emilie and Julie : identifying differences in coding FPs

*FP evaluations comparison (w/o FP_group_2)
foreach y in 1934 1943 1960 1962 /*1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003*/ {
import excel "$ED\FP check 2\Emilie\FP3_check_`y'", sheet("Sheet1") firstrow clear
rename (FP_single FP_group) (FP_single_ED FP_group_ED)
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP3_check_`y'_ED.dta", replace

import excel "$ED\FP check 2\Julie\FP3_check_`y'", sheet("Sheet1") firstrow clear
rename (FP_single FP_group) (FP_single_JU FP_group_JU)
merge 1:1 year CID using "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP3_check_`y'_ED.dta"
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP3_check_`y'_ED_JU.dta", replace

gen diff_check = 0 if FP_single_ED == FP_single_JU & FP_group_ED == FP_group_JU
replace diff_check = 1 if missing(diff_check)
levelsof UCID_3 if diff_check == 1, local(diff)
gen diff_check_2 = 0
foreach l in `diff' {
 	replace diff_check_2 = 1 if UCID_3 == `l' 
 }
 
preserve
keep if diff_check_2 == 1
sort UCID_3 cname
drop match_2 match_3 _merge diff_check_2
gen FP_single = 0
gen FP_group = 0
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP3_check_`y'_diff.dta", replace
export excel using "$ED\FP check 2\Final\FP3_check_`y'.xlsx", replace firstrow(variables)
restore

preserve
drop if diff_check_2 == 1
save "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP3_check_`y'_same.dta", replace
restore
}

********************************************************************************
* Use final.xlsx and compare how many times I agree with ED or JU over the first 4 years: is it random or not?
* Prepare code to correct matched groups
********************************************************************************

* Comparisons JU, ED vs YS
foreach y in 1934 1943 1960 1962 /*1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003*/ {
import excel "$ED\FP check 2\Final\FP3_check_`y'_final", sheet("Sheet1") firstrow clear
replace FP_single_JU = 0 if missing(FP_single_JU)
replace FP_group_JU = 0 if missing(FP_group_JU)
replace FP_single_ED = 0 if missing(FP_single_ED)
replace FP_group_ED = 0 if missing(FP_group_ED)

gen JU_ok = 1 if FP_single_JU == FP_single & FP_group_JU == FP_group
gen ED_ok = 1 if FP_single_ED == FP_single & FP_group_ED == FP_group
gen nb_ju_ok = sum(JU_ok)
gen nb_ed_ok = sum(ED_ok)

save "$ED\FP check 2\Final\me_ed_ju_comparison_`y'.dta", replace
}*/


********************************************************************************
*Cleaning of duplicates 20% thresholds

* Files evaluated by JU, ED, and myself
foreach y in 1934 1943 1960 1962 {

***
*FP subsample 
*Add corrections and post processing : used ED & JU xls, had a look at UCID_2 listed above in diff.dta, made decisions...

*creating file with all evaluations incl my final decision 
import excel "$ED\FP check 2\Final\FP3_check_`y'_final.xlsx", sheet("Sheet1") firstrow clear
save "$ED\FP check 2\Final\FP3_check_`y'_final.dta", replace

use "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP3_check_`y'_ED_JU.dta", clear
drop _merge
merge 1:1 CID year using "$ED\FP check 2\Final\FP3_check_`y'_final.dta"

replace FP_single = FP_single_ED if FP_single_ED==FP_single_JU & missing(FP_single)
replace FP_group = FP_group_ED if FP_group_ED==FP_group_JU & missing(FP_group)
replace FP_single = 0 if missing(FP_single)
replace FP_group = 0 if missing(FP_group)

sort UCID_3 cname
drop FP_single_JU FP_group_JU FP_single_ED FP_group_ED diff_check _merge


**Collapse those who were right : FP and were evaluated as correct (FP_single and FP_group == 0)
*Count number of commas to select "best name": gen comma = length(trim(var1))-length(trim(subinstr(nb_owners,",","",.))) where nb_owners is the count minus 1 (replace nb_owners = nb_owners + 1)
preserve
	keep if FP_single ==0 & FP_group ==0
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP.dta"
	keep if _merge ==3
	drop _merge
	
	***collapse multiple perfect matches in first round
		*Build unique CID_str for all UCID_SI_b (perfect match on owners & 20% cname)
		*Count number of commas
		replace owners = usubinstr(owners,  ";", ",",.)
		gen comma = (length(trim(owners))-length(trim(subinstr(owners,",","",.))))
		sort UCID_3 comma cname
		order CID CID_str comma
				
		bysort UCID_3: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_3: replace CID_str = CID_str[_N]
		
		
		*Same with cname (create cname group)
		bysort UCID_3: replace cname_group = cname_group[_n-1] + "; " + cname_group 
		replace cname_group = " "+cname_group
		replace cname_group = usubinstr(cname_group, " ; ", "",.)
		replace cname_group = ustrtrim(cname_group)
		by UCID_3: replace cname_group = cname_group[_N]

		*Same with owners
		bysort UCID_3: replace owners = owners[_n-1] + "; " + owners
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		bysort UCID_3: replace owners = owners[_N]

		*Same with function
		bysort UCID_3: replace function = function[_n-1] + "; " + function
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_3: replace function = function[_N]

		*Same with signature
		bysort UCID_3: replace signature = signature[_n-1] + "; " + signature
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_3: replace signature = signature[_N]

		*Same with page
		bysort UCID_3: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_3: replace page = page[_N]

		*Delete duplicates
		sort UCID_3 comma cname
		by UCID_3:  gen dup = cond(_N==1,0,_n)
		order CID CID_str comma dup
		bysort UCID_3 : keep if dup == dup[_N]
	gen UCID_3_temp2 = _n // temp1 is in nonFP subsample	
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_1.dta", replace
restore

*Correct FP_single who were wrong (give unique ID to those obs)
preserve
keep if FP_single==1 
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3.dta"
	keep if _merge ==3
	drop _merge
	gen UCID_3_temp3 = _n
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_2.dta", replace
restore

*Correct FP_group who were wrong (give unique ID to those group of obs)
preserve
keep if FP_group==1 
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3.dta"
	keep if _merge ==3
	drop _merge
		
	***collapse multiple perfect matches in first round
		*Build unique CID_str for all UCID_SI_b (perfect match on owners & 20% cname)
		*Count number of commas
		replace owners = usubinstr(owners,  ";", ",",.)
		gen comma = (length(trim(owners))-length(trim(subinstr(owners,",","",.))))
		sort UCID_3 comma cname
		order CID CID_str comma
		
		*Same with cname (create cname group)
		bysort UCID_3: replace cname_group = cname_group[_n-1] + "; " + cname_group 
		replace cname_group = " "+cname_group
		replace cname_group = usubinstr(cname_group, " ; ", "",.)
		replace cname_group = ustrtrim(cname_group)
		by UCID_3: replace cname_group = cname_group[_N]
				
		bysort UCID_3: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_3: replace CID_str = CID_str[_N]

		*Same with owners
		bysort UCID_3: replace owners = owners[_n-1] + "; " + owners
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		bysort UCID_3: replace owners = owners[_N]

		*Same with function
		bysort UCID_3: replace function = function[_n-1] + "; " + function
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_3: replace function = function[_N]

		*Same with signature
		bysort UCID_3: replace signature = signature[_n-1] + "; " + signature
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_3: replace signature = signature[_N]

		*Same with page
		bysort UCID_3: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_3: replace page = page[_N]

		*Delete duplicates
		sort UCID_3 comma cname
		by UCID_3:  gen dup = cond(_N==1,0,_n)
		order CID CID_str comma dup
		bysort UCID_3 : keep if dup == dup[_N]
	gen UCID_3_temp4 = _n
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_3.dta", replace
restore

*Putting everything back together in "match_3_clean"
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_nonFP.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_1.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_2.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_3.dta"
sort UCID_3_temp1 UCID_3_temp2 UCID_3_temp3 UCID_3_temp4
order UCID_3_temp1 UCID_3_temp2 UCID_3_temp3 UCID_3_temp4

* Create one single UCID_3

egen maxUCID_temp3=max(UCID_3_temp3)
replace UCID_3_temp3 = maxUCID_temp3 + UCID_3_temp4 if missing(UCID_3_temp3)
egen maxUCID_temp2=max(UCID_3_temp2)
replace UCID_3_temp2 = maxUCID_temp2 + UCID_3_temp3 if missing(UCID_3_temp2)
egen maxUCID_temp1=max(UCID_3_temp1)
replace UCID_3_temp1 = maxUCID_temp1 + UCID_3_temp2 if missing(UCID_3_temp1)

drop UCID_3_temp2 UCID_3_temp3 UCID_3_temp4 maxUCID_temp1 maxUCID_temp2 maxUCID_temp3 UCID_2 dup maxUCID_nonSI comma FP_single FP_group UCID_3
rename UCID_3_temp1 UCID_3 
sort gdenr_2018 cname

*Cleaning of cname_group
replace cname_group = " "+cname_group
replace cname_group = usubinstr(cname_group, " ;", "",.)
replace cname_group = ustrtrim(cname_group)
replace cname_group = cname if missing(cname_group)

*Cleaning of duplicates coming from differences in sorting : both observations are 100% equivalent except UCID_3
bysort CID: gen CIDdup = cond(_N==1,0,_n)
drop if CIDdup > 1
drop CIDdup

save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean_temp.dta", replace
}

* Files evaluated only by JU or ED (not both)
foreach y in 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003{

***
*FP subsample only
*creating file with all evaluations incl my final decision 
import excel "$ED\FP check 2\Final\FP3_check_`y'_final.xlsx", sheet("Sheet1") firstrow clear
save "$ED\FP check 2\Final\FP3_check_`y'_final.dta", replace

use "$ED\FP check 2\Final\FP3_check_`y'_final.dta", clear
destring FP_single FP_group, replace
replace FP_single = 0 if missing(FP_single)
replace FP_group = 0 if missing(FP_group)

sort UCID_3 cname

**Collapse those who were right : FP and were evaluated as correct (FP_single and FP_group == 0)
*Count number of commas to select "best name": gen comma = length(trim(var1))-length(trim(subinstr(nb_owners,",","",.))) where nb_owners is the count minus 1 (replace nb_owners = nb_owners + 1)
preserve
	keep if FP_single ==0 & FP_group ==0
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP.dta"
	keep if _merge ==3
	drop _merge
	
	***collapse multiple perfect matches in first round
		*Build unique CID_str for all UCID_SI_b (perfect match on owners & 20% cname)
		*Count number of commas
		replace owners = usubinstr(owners,  ";", ",",.)
		gen comma = (length(trim(owners))-length(trim(subinstr(owners,",","",.))))
		sort UCID_3 comma cname
		order CID CID_str comma
				
		bysort UCID_3: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_3: replace CID_str = CID_str[_N]
		
		*Same with cname (create cname group)
		bysort UCID_3: replace cname_group = cname_group[_n-1] + "; " + cname_group 
		replace cname_group = " "+cname_group
		replace cname_group = usubinstr(cname_group, " ; ", "",.)
		replace cname_group = ustrtrim(cname_group)
		by UCID_3: replace cname_group = cname_group[_N]

		*Same with owners
		bysort UCID_3: replace owners = owners[_n-1] + "; " + owners
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		bysort UCID_3: replace owners = owners[_N]

		*Same with function
		bysort UCID_3: replace function = function[_n-1] + "; " + function
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_3: replace function = function[_N]

		*Same with signature
		bysort UCID_3: replace signature = signature[_n-1] + "; " + signature
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_3: replace signature = signature[_N]

		*Same with page
		bysort UCID_3: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_3: replace page = page[_N]

		*Delete duplicates
		sort UCID_3 comma cname
		by UCID_3:  gen dup = cond(_N==1,0,_n)
		order CID CID_str comma dup
		bysort UCID_3 : keep if dup == dup[_N]
	gen UCID_3_temp2 = _n // temp1 is in nonFP subsample	
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_1.dta", replace
restore

*Correct FP_single who were wrong (give unique ID to those obs)
preserve
keep if FP_single==1 
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3.dta"
	keep if _merge ==3
	drop _merge
	gen UCID_3_temp3 = _n
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_2.dta", replace
restore

*Correct FP_group who were wrong (give unique ID to those group of obs)
preserve
keep if FP_group==1 
	merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3.dta"
	keep if _merge ==3
	drop _merge
		
	***collapse multiple perfect matches in first round
		*Build unique CID_str for all UCID_SI_b (perfect match on owners & 20% cname)
		*Count number of commas
		replace owners = usubinstr(owners,  ";", ",",.)
		gen comma = (length(trim(owners))-length(trim(subinstr(owners,",","",.))))
		sort UCID_3 comma cname
		order CID CID_str comma
		
		*Same with cname (create cname group)
		bysort UCID_3: replace cname_group = cname_group[_n-1] + "; " + cname_group 
		replace cname_group = " "+cname_group
		replace cname_group = usubinstr(cname_group, " ; ", "",.)
		replace cname_group = ustrtrim(cname_group)
		by UCID_3: replace cname_group = cname_group[_N]
				
		bysort UCID_3: replace CID_str = CID_str[_n-1] + "; " + CID_str
		replace CID_str = " "+CID_str
		replace CID_str = usubinstr(CID_str, " ; ", "",.)
		by UCID_3: replace CID_str = CID_str[_N]

		*Same with owners
		bysort UCID_3: replace owners = owners[_n-1] + "; " + owners
		replace owners = " "+owners
		replace owners = usubinstr(owners, " ; ", "",.)
		bysort UCID_3: replace owners = owners[_N]

		*Same with function
		bysort UCID_3: replace function = function[_n-1] + "; " + function
		replace function = " "+function
		replace function = usubinstr(function, " ; ", "",.)
		by UCID_3: replace function = function[_N]

		*Same with signature
		bysort UCID_3: replace signature = signature[_n-1] + "; " + signature
		replace signature = " "+signature
		replace signature = usubinstr(signature, " ; ", "",.)
		by UCID_3: replace signature = signature[_N]

		*Same with page
		bysort UCID_3: replace page = page[_n-1] + "; " + page
		replace page = " "+page
		replace page = usubinstr(page, " ; ", "",.)
		by UCID_3: replace page = page[_N]

		*Delete duplicates
		sort UCID_3 comma cname
		by UCID_3:  gen dup = cond(_N==1,0,_n)
		order CID CID_str comma dup
		bysort UCID_3 : keep if dup == dup[_N]
	gen UCID_3_temp4 = _n
	save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_3.dta", replace
restore

*Putting everything back together in "match_3_clean"
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_nonFP.dta", clear
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_1.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_2.dta"
append using "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_3.dta"
sort UCID_3_temp1 UCID_3_temp2 UCID_3_temp3 UCID_3_temp4
order UCID_3_temp1 UCID_3_temp2 UCID_3_temp3 UCID_3_temp4

* Create one single UCID_3

egen maxUCID_temp3=max(UCID_3_temp3)
replace UCID_3_temp3 = maxUCID_temp3 + UCID_3_temp4 if missing(UCID_3_temp3)
egen maxUCID_temp2=max(UCID_3_temp2)
replace UCID_3_temp2 = maxUCID_temp2 + UCID_3_temp3 if missing(UCID_3_temp2)
egen maxUCID_temp1=max(UCID_3_temp1)
replace UCID_3_temp1 = maxUCID_temp1 + UCID_3_temp2 if missing(UCID_3_temp1)

drop UCID_3_temp2 UCID_3_temp3 UCID_3_temp4 maxUCID_temp1 maxUCID_temp2 maxUCID_temp3 UCID_2 dup maxUCID_nonSI comma FP_single FP_group UCID_3
rename UCID_3_temp1 UCID_3 
sort gdenr_2018 cname

*Cleaning of cname_group
replace cname_group = " "+cname_group
replace cname_group = usubinstr(cname_group, " ;", "",.)
replace cname_group = ustrtrim(cname_group)
replace cname_group = cname if missing(cname_group)

/*Cleaning of duplicates : both observations are 100% equivalent except UCID_3
bysort CID year: gen CIDdup = cond(_N==1,0,_n)
drop if CIDdup > 1
drop CIDdup*/

save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean_temp.dta", replace
}

********************************************************************************
* Add gdenr_2012 to contemporary state of mun
********************************************************************************
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003{
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean_temp.dta", clear
merge m:1 gdenr year using "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\08_Municipalities\Gde_1931_2018.dta", keep(1 3) nogen
drop if year != `y'
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", replace
}

*Final UCID_3 in cross-sections (to have proper 1 to N structure)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003{
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", replace
gen UCID_3_temp = _n
order UCID_3 UCID_3_temp
drop UCID_3
rename UCID_3_temp UCID_3
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", replace
}

********************************************************************************
* Once we have match 3 clean done and saved
* Delete every small intermediate steps
********************************************************************************
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003{
*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_1.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_1_FP.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_SI.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_nonSI.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_SI_a.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_nonFP.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_toeval.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_2_FP_ok.dta"
erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_1.dta"
erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_2.dta"
erase "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_FP_check_3.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP3_check_`y'_ED.dta"
*erase "$dump\02_Processed_data\10_Directors_1934_2003\17_FP_Check\FP3_check_`y'_ED_JU.dta"
}

********************************************************************************
*** Generate yearly tables+graphs  = comparison Sugarcube v Ragionenbuch     ***
********************************************************************************

*number of UCID by year and canton for each threshold

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", clear
		gen tag = 1
		egen ndistinct = total(tag), by(ctn)
		sort ctn 
		collapse (mean) ndistinct, by(ctn)
		rename ndistinct sugarcube
		gen year = `y'
		replace ctn = "n/a" if missing(ctn)
	save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count.dta", replace
}

* Construct single file with all years

use "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_count.dta", clear
foreach x in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`x'_count.dta"
}
save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_count_1934_2003_sug.dta", replace

* Add stock of SA+succursale as reported from the Ragionenbuch
import excel "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\Stock_company.xlsx", sheet("Sheet1") firstrow clear
keep ctn ctn_long year sa succursale
merge 1:1 ctn year using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_count_1934_2003_sug.dta"
drop if _merge == 1
drop _merge
gen ragionenbuch_full = sa + succursale
order ctn ctn_long year sa succursale ragionenbuch
drop succursale
rename sa ragionenbuch_sa

gen diff_ragio_sa = ragionenbuch_sa - sugarcube
gen diff_ragio_full = ragionenbuch_full - sugarcube

save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_count_1934_2003_sug_ragionenbuch.dta", replace

* Add HSSO data (Yearbook data)

import excel "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\S.05a.xlsx", sheet("Worksheet") cellrange(A3:AO70) firstrow clear
drop if missing(Jahr)
drop Zürich Basel Genève Bern Lau St Winter Luzern Biel LaChaux Städtevilles AN Année

foreach y in ZH BE LU UR SZ OW NW GL ZG FR SO BS BL SH AR AI SG GR AG TG TI VD VS NE GE JU CH {
preserve
	keep Jahr `y'
	gen ctn = "`y'"
	rename (Jahr `y') (year yearbook_sa)
	save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\yearbook_`y'.dta", replace
restore
}

use "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\yearbook_zh.dta", clear
foreach y in BE LU UR SZ OW NW GL ZG FR SO BS BL SH AR AI SG GR AG TG TI VD VS NE GE JU CH {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\yearbook_`y'.dta"
	replace ctn = ustrlower(ctn)
}

foreach y in ZH BE LU UR SZ OW NW GL ZG FR SO BS BL SH AR AI SG GR AG TG TI VD VS NE GE JU CH {
	erase "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\yearbook_`y'.dta"
}

merge 1:1 ctn year using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_count_1934_2003_sug_ragionenbuch.dta"
drop if _merge == 1
drop _merge

order ctn ctn_long year sugarcube ragionenbuch_sa ragionenbuch_full diff_ragio_sa diff_ragio_full yearbook_sa
gen diff_yearbook_sa = yearbook_sa - sugarcube
gen p_diff_ragio_sa = (diff_ragio_sa/ragionenbuch_sa)*100
gen p_diff_ragio_full = (diff_ragio_full/ragionenbuch_full)*100
gen p_diff_yearbook = (diff_yearbook_sa/yearbook_sa)*100

save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_count_1934_2003_sug_ragionenbuch_yearbook.dta", replace
export excel "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_count_1934_2003_sug_ragionenbuch_yearbook.xlsx", replace firstrow(variables)

*********************************************
*** Generate yearly stock of firms, avg owners by mun ***
*********************************************

*number of UCID by year and gdenr contemporary

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", clear
	
		*Uniformisation of ; and ,
		replace owners = usubinstr(owners, ";", ",",.)
		replace owners = usubinstr(owners, "  ", " ",.)
		replace owners = usubinstr(owners, " ", "",.)

		merge 1:1 UCID_3 using "$dump\02_Processed_data\10_Directors_1934_2003\11_Ownerlists\nb_owner_`y'.dta", keep(3) nogen 
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr)
		drop tag
		egen total_owners = total(nb_owners), by(gdenr)
		gen avg_owners = total_owners/ndistinct
		egen min_owners = min(nb_owners), by(gdenr)
		egen max_owners = max(nb_owners), by(gdenr)
		gen log_avg_owners = ln(avg_owners)
		
		summarize nb_owners, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_o = r(p`p') 
		}
		
		gen bottom50_dummy_o = 1 if nb_owners<=p50_thresh_o
		gen top50_dummy_o = 1 if nb_owners>p50_thresh_o
		gen top75_dummy_o = 1 if nb_owners>p75_thresh_o
		gen top90_dummy_o = 1 if nb_owners>p90_thresh_o
		gen top99_dummy_o = 1 if nb_owners>p99_thresh_o
		
		egen total_nb_bottom50_o = sum(bottom50_dummy_o), by(gdenr)
		egen total_nb_top50_o = sum(top50_dummy), by(gdenr)
		egen total_nb_top75_o = sum(top75_dummy), by(gdenr)
		egen total_nb_top90_o = sum(top90_dummy), by(gdenr)
		egen total_nb_top99_o = sum(top99_dummy), by(gdenr)
		
		gen total_nb_top50_to_75_o = total_nb_top50_o - total_nb_top75_o
		gen total_nb_top75_to_90_o = total_nb_top75_o - total_nb_top90_o
		gen total_nb_top90_to_99_o = total_nb_top90_o - total_nb_top99_o
		gen total_nb_top90_to_100_o = total_nb_top90_o
		gen total_nb_top99_to_100_o = total_nb_top99_o
		
		gen single_owner = 1 if nb_owners == 1
		gen multiple_owner = 1 if nb_owners > 1
		gen two_owners = 1 if nb_owners == 2
		gen three_owners = 1 if nb_owners == 3
		gen four_owners = 1 if nb_owners == 4
		gen five_owners = 1 if nb_owners == 5
		gen six_more_owners = 1 if nb_owners >= 6
		gen five_more_owners = 1 if nb_owners >= 5
		gen two_four_owners = 1 if nb_owners >= 2 & nb_owners < 5     
		
		replace single_owner = 0 if missing(single_owner)
		replace multiple_owner = 0 if missing(multiple_owner)
		replace two_owners = 0 if missing(two_owners)
		replace three_owners = 0 if missing(three_owners)
		replace four_owners = 0 if missing(four_owners)
		replace five_more_owners = 0 if missing(five_more_owners)
		replace five_owners = 0 if missing(single_owner)
		replace six_more_owners = 0 if missing(single_owner)
		replace two_four_owners = 0 if missing(two_four_owners)

		collapse (mean) nb_owners ndistinct total_owners avg_owners min_owners max_owners log_avg_owners p1_thresh_o p5_thresh_o p10_thresh_o p25_thresh_o p50_thresh_o p75_thresh_o p90_thresh_o p99_thresh_o bottom50_dummy_o top50_dummy_o top75_dummy_o top90_dummy_o top99_dummy_o total_nb_bottom50_o total_nb_top50_o total_nb_top75_o total_nb_top90_o total_nb_top99_o total_nb_top50_to_75_o total_nb_top75_to_90_o  total_nb_top90_to_100_o total_nb_top90_to_99_o total_nb_top99_to_100_o year (sum) single_owner multiple_owner two_owners three_owners four_owners five_owners two_four_owners five_more_owners six_more_owners, by(gdenr)
save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun.dta", replace
}

*number of UCID, nb of owners (and diff measures) by year and gdenr 2012

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", clear
		*Uniformisation of ; and ,
		replace owners = usubinstr(owners, ";", ",",.)
		replace owners = usubinstr(owners, "  ", " ",.)
		replace owners = usubinstr(owners, " ", "",.)

		merge 1:1 UCID_3 using "$dump\02_Processed_data\10_Directors_1934_2003\11_Ownerlists\nb_owner_`y'.dta", keep(3) nogen 
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr_2012)
		drop tag
		egen total_owners = total(nb_owners), by(gdenr_2012)
		gen avg_owners = total_owners/ndistinct
		egen min_owners = min(nb_owners), by(gdenr_2012)
		egen max_owners = max(nb_owners), by(gdenr_2012)
		gen log_avg_owners = ln(avg_owners)
		
		summarize nb_owners, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_o = r(p`p') 
		}
		
		gen bottom50_dummy_o = 1 if nb_owners<=p50_thresh_o
		gen top50_dummy_o = 1 if nb_owners>p50_thresh_o
		gen top75_dummy_o = 1 if nb_owners>p75_thresh_o
		gen top90_dummy_o = 1 if nb_owners>p90_thresh_o
		gen top99_dummy_o = 1 if nb_owners>p99_thresh_o
		
		egen total_nb_bottom50_o = sum(bottom50_dummy_o), by(gdenr_2012)
		egen total_nb_top50_o = sum(top50_dummy), by(gdenr_2012)
		egen total_nb_top75_o = sum(top75_dummy), by(gdenr_2012)
		egen total_nb_top90_o = sum(top90_dummy), by(gdenr_2012)
		egen total_nb_top99_o = sum(top99_dummy), by(gdenr_2012)
		
		gen total_nb_top50_to_75_o = total_nb_top50_o - total_nb_top75_o
		gen total_nb_top75_to_90_o = total_nb_top75_o - total_nb_top90_o
		gen total_nb_top90_to_99_o = total_nb_top90_o - total_nb_top99_o
		gen total_nb_top90_to_100_o = total_nb_top90_o
		gen total_nb_top99_to_100_o = total_nb_top99_o
		
		gen single_owner = 1 if nb_owners == 1
		gen multiple_owner = 1 if nb_owners > 1
		gen two_owners = 1 if nb_owners == 2
		gen three_owners = 1 if nb_owners == 3
		gen four_owners = 1 if nb_owners == 4
		gen five_owners = 1 if nb_owners == 5
		gen six_more_owners = 1 if nb_owners >= 6
		gen five_more_owners = 1 if nb_owners >= 5
		gen two_four_owners = 1 if nb_owners >= 2 & nb_owners < 5     
		
		replace single_owner = 0 if missing(single_owner)
		replace multiple_owner = 0 if missing(multiple_owner)
		replace two_owners = 0 if missing(two_owners)
		replace three_owners = 0 if missing(three_owners)
		replace four_owners = 0 if missing(four_owners)
		replace five_more_owners = 0 if missing(five_more_owners)
		replace five_owners = 0 if missing(single_owner)
		replace six_more_owners = 0 if missing(single_owner)
		replace two_four_owners = 0 if missing(two_four_owners)

		collapse (mean) nb_owners ndistinct total_owners avg_owners min_owners max_owners log_avg_owners p1_thresh_o p5_thresh_o p10_thresh_o p25_thresh_o p50_thresh_o p75_thresh_o p90_thresh_o p99_thresh_o bottom50_dummy_o top50_dummy_o top75_dummy_o top90_dummy_o top99_dummy_o total_nb_bottom50_o total_nb_top50_o total_nb_top75_o total_nb_top90_o total_nb_top99_o total_nb_top50_to_75_o total_nb_top75_to_90_o  total_nb_top90_to_100_o total_nb_top90_to_99_o total_nb_top99_to_100_o year (sum) single_owner multiple_owner two_owners three_owners four_owners five_owners two_four_owners five_more_owners six_more_owners, by(gdenr_2012)
		
	save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_12.dta", replace
}
	
* Append all years (gdenr_2012)
use "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_count_mun_12.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_12.dta"
}

*drop if missing(gdenr_2012)
save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_2003_count_owners_mun_12.dta", replace


*number of UCID, nb of owners (and diff measures) by year and gdenr 2018 

foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", clear
		*Uniformisation of ; and ,
		replace owners = usubinstr(owners, ";", ",",.)
		replace owners = usubinstr(owners, "  ", " ",.)
		replace owners = usubinstr(owners, " ", "",.)

		merge 1:1 UCID_3 using "$dump\02_Processed_data\10_Directors_1934_2003\11_Ownerlists\nb_owner_`y'.dta", keep(3) nogen 
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr_2018)
		drop tag
		egen total_owners = total(nb_owners), by(gdenr_2018)
		gen avg_owners = total_owners/ndistinct
		egen min_owners = min(nb_owners), by(gdenr_2018)
		egen max_owners = max(nb_owners), by(gdenr_2018)
		gen log_avg_owners = ln(avg_owners)
		
		summarize nb_owners, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_o = r(p`p') 
		}
		
		gen bottom50_dummy_o = 1 if nb_owners<=p50_thresh_o
		gen top50_dummy_o = 1 if nb_owners>p50_thresh_o
		gen top75_dummy_o = 1 if nb_owners>p75_thresh_o
		gen top90_dummy_o = 1 if nb_owners>p90_thresh_o
		gen top99_dummy_o = 1 if nb_owners>p99_thresh_o
		
		egen total_nb_bottom50_o = sum(bottom50_dummy_o), by(gdenr_2018)
		egen total_nb_top50_o = sum(top50_dummy), by(gdenr_2018)
		egen total_nb_top75_o = sum(top75_dummy), by(gdenr_2018)
		egen total_nb_top90_o = sum(top90_dummy), by(gdenr_2018)
		egen total_nb_top99_o = sum(top99_dummy), by(gdenr_2018)
		
		gen total_nb_top50_to_75_o = total_nb_top50_o - total_nb_top75_o
		gen total_nb_top75_to_90_o = total_nb_top75_o - total_nb_top90_o
		gen total_nb_top90_to_99_o = total_nb_top90_o - total_nb_top99_o
		gen total_nb_top90_to_100_o = total_nb_top90_o
		gen total_nb_top99_to_100_o = total_nb_top99_o
		
		gen single_owner = 1 if nb_owners == 1
		gen multiple_owner = 1 if nb_owners > 1
		gen two_owners = 1 if nb_owners == 2
		gen three_owners = 1 if nb_owners == 3
		gen four_owners = 1 if nb_owners == 4
		gen five_owners = 1 if nb_owners == 5
		gen six_more_owners = 1 if nb_owners >= 6
		gen five_more_owners = 1 if nb_owners >= 5
		gen two_four_owners = 1 if nb_owners >= 2 & nb_owners < 5     
		
		replace single_owner = 0 if missing(single_owner)
		replace multiple_owner = 0 if missing(multiple_owner)
		replace two_owners = 0 if missing(two_owners)
		replace three_owners = 0 if missing(three_owners)
		replace four_owners = 0 if missing(four_owners)
		replace five_more_owners = 0 if missing(five_more_owners)
		replace five_owners = 0 if missing(single_owner)
		replace six_more_owners = 0 if missing(single_owner)
		replace two_four_owners = 0 if missing(two_four_owners)

		collapse (mean) nb_owners ndistinct total_owners avg_owners min_owners max_owners log_avg_owners p1_thresh_o p5_thresh_o p10_thresh_o p25_thresh_o p50_thresh_o p75_thresh_o p90_thresh_o p99_thresh_o bottom50_dummy_o top50_dummy_o top75_dummy_o top90_dummy_o top99_dummy_o total_nb_bottom50_o total_nb_top50_o total_nb_top75_o total_nb_top90_o total_nb_top99_o total_nb_top50_to_75_o total_nb_top75_to_90_o  total_nb_top90_to_100_o total_nb_top90_to_99_o total_nb_top99_to_100_o year (sum) single_owner multiple_owner two_owners three_owners four_owners five_owners two_four_owners five_more_owners six_more_owners, by(gdenr_2018)
		
	save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_18.dta", replace
}

* Append all years (gdenr_2018)
use "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_count_mun_18.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_`y'_count_mun_18.dta"
}
*drop if missing(gdenr_2018)
save "$dump\02_Processed_data\10_Directors_1934_2003\13_Counts\comp_1934_2003_count_owners_mun_18.dta", replace

********************************************************************************
*** Clean nominal capital and generate measures at mun-level ***
****************************************************************

*** Corrections of capital
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {

use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\comp_unique_id_`y'_match_3_clean.dta", clear

rename capital capital_orig
gen capital = capital_orig
replace capital = ustrlower(capital)
replace capital = ustrregexra(capital,"[a-z]","")
replace capital = usubinstr(capital, ")", "",.)
replace capital = usubinstr(capital, "í", "",.)
replace capital = usubinstr(capital, "'", "",.)
replace capital = usubinstr(capital, "-", "",.)
replace capital = usubinstr(capital, ",", "",.)
replace capital = usubinstr(capital, "!", "",.)
replace capital = usubinstr(capital, "`", "",.)
replace capital = usubinstr(capital, "�", "",.)
replace capital = usubinstr(capital, "/", "",.)
replace capital = usubinstr(capital, "?", "",.)
replace capital = usubinstr(capital, ":", "",.)
replace capital = usubinstr(capital, "è", "",.)
replace capital = usubinstr(capital, "ó", "",.)
replace capital = usubinstr(capital, "ö", "",.)
replace capital = usubinstr(capital, "ü", "",.)
replace capital = usubinstr(capital, "ä", "",.)
replace capital = usubinstr(capital, "   ", " ",.)
replace capital = usubinstr(capital, "  ", " ",.)
replace capital = usubinstr(capital, " ", "",.)
replace capital = ustrtrim(capital)
destring capital, replace force 
replace capital = 0.5 if capital_orig == "0.S"
replace capital = 0.5 if capital_orig == "0.5)."
replace capital = 0.05 if capital_orig == "0:05"
replace capital = 0.3 if capital_orig == "0:3"
replace capital = 0.3 if capital_orig == "0.3) 5"
replace capital = 0.1 if capital_orig == "0-1"
replace capital = 0.1 if capital_orig == "0.'"
replace capital = 0.051 if capital_orig == "0.05l"
replace capital = 32.66 if capital_orig == "32.6602).5"
replace capital = 0.1 if capital_orig == "/0,1"
replace capital = 145.0 if capital_orig == "145:0"
replace capital = 0.1 if capital_orig == "0.1) 5"
replace capital = 0.98 if capital_orig == "O,98"
replace capital = 0.06 if capital_orig == "0.C6"
replace capital = . if capital_orig == "Liesta140,1"
replace capital = . if capital_orig == ",4,0"
replace capital = . if capital_orig == "41,"
replace capital = 0.1 if capital_orig == "0.i"
replace capital = 0.05 if capital_orig == "0..05"
replace capital = 0.05 if capital_orig == "0.0.5"
replace capital = 0.05 if capital_orig == "0.05) Pr.-"
replace capital = 0.1 if capital_orig == "0.1) Vr."
replace capital = 0.1 if capital_orig == "0.1."
replace capital = 0.15 if capital_orig == "0.15."
replace capital = 0.1 if capital_orig == "0.25)' t..."

preserve
	keep if missing(capital) | capital == .
	save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_missing_capital_mun.dta", replace
	keep UCID_3 CID year capital_orig cname_orig caddress_orig gdename_orig gdename_extract capital_extract cname_extract capital
restore

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_mun.dta", replace
}


*Missing capital for ED JU
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_missing_capital_mun.dta", clear
drop match_1 match_2 match_3 function signature	page CID_str si_dummy_2 geo_merge branch_dummy gdename_extract capital_extract	cname_extract gdenr gdenr_2018 ctn E_CNTR N_CNTR
export excel using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\Capital check\capital_check_`y'.xlsx", replace firstrow(variables)
}

*Include corrections they've done and add corrected ones
*transform xlsx into dta (for merging)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
import excel "C:\Users\schmutzy\Dropbox\Travaux_EmilieDousse\Capital Check\capital_check_`y'.xlsx", sheet("Sheet1") firstrow clear
keep UCID_3 CID year cname gdename owners capital_orig cname_orig capital capital_corrected
destring capital_corrected, replace

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_missing_capital_mun_corrected.dta", replace
}

*Put back into cross-sections
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_mun.dta", clear
merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_missing_capital_mun_corrected.dta", keep (1 3) nogen
replace capital = capital_corrected if missing(capital)
drop capital_corrected
save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_mun.dta", replace
}

*** Generate statistics on capital at municipal-level (gdenr_2012)
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_mun.dta", clear
		
		drop if missing(capital) | capital == .
		
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr_2012)
		drop tag
		egen total_capital = sum(capital), by(gdenr_2012)
		gen avg_capital = total_capital/ndistinct
		egen min_capital = min(capital), by(gdenr_2012)
		egen max_capital = max(capital), by(gdenr_2012)
		gen log_avg_capital = ln(avg_capital)
		
		summarize capital, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_c = r(p`p') 
		}
		
		gen bottom50_dummy_c = 1 if capital<=p50_thresh_c
		gen top50_dummy_c = 1 if capital>p50_thresh_c
		gen top75_dummy_c = 1 if capital>p75_thresh_c
		gen top90_dummy_c = 1 if capital>p90_thresh_c
		gen top99_dummy_c = 1 if capital>p99_thresh_c
		
		egen total_nb_bottom50_c = sum(bottom50_dummy_c), by(gdenr_2012)
		egen total_nb_top50_c = sum(top50_dummy_c), by(gdenr_2012)
		egen total_nb_top75_c = sum(top75_dummy_c), by(gdenr_2012)
		egen total_nb_top90_c = sum(top90_dummy_c), by(gdenr_2012)
		egen total_nb_top99_c = sum(top99_dummy_c), by(gdenr_2012)

		gen total_nb_top50_to_75_c = total_nb_top50_c - total_nb_top75_c
		gen total_nb_top75_to_90_c = total_nb_top75_c - total_nb_top90_c
		gen total_nb_top90_to_99_c = total_nb_top90_c - total_nb_top99_c
		gen total_nb_top90_to_100_c = total_nb_top90_c
		gen total_nb_top99_to_100_c = total_nb_top99_c
		
		collapse (mean) ndistinct total_capital avg_capital min_capital max_capital log_avg_capital	p1_thresh_c p5_thresh_c p10_thresh_c p25_thresh_c p50_thresh_c p75_thresh_c p90_thresh_c p99_thresh_c bottom50_dummy_c top50_dummy_c top75_dummy_c top90_dummy_c top99_dummy_c total_nb_bottom50_c total_nb_top50_c total_nb_top75_c total_nb_top90_c total_nb_top99_c total_nb_top50_to_75_c total_nb_top75_to_90_c total_nb_top90_to_100_c total_nb_top90_to_99_c total_nb_top99_to_100_c year, by(gdenr_2012)
		
save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_12.dta", replace
}

* Append all years (gdenr_2012)
use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_capital_wstat_mun_12.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_12.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_2003_capital_wstat_mun_12.dta", replace


*** Generate statistics on capital at municipal-level
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_mun.dta", clear
		
		drop if missing(capital) | capital == .
		
		gen tag = 1
		egen ndistinct = total(tag), by(gdenr_2018)
		drop tag
		egen total_capital = sum(capital), by(gdenr_2018)
		gen avg_capital = total_capital/ndistinct
		egen min_capital = min(capital), by(gdenr_2018)
		egen max_capital = max(capital), by(gdenr_2018)
		gen log_avg_capital = ln(avg_capital)
		
		summarize capital, detail
		
		local perc 1 5 10 25 50 75 90 99
		foreach p of local perc {
		gen p`p'_thresh_c = r(p`p') 
		}
		
		gen bottom50_dummy_c = 1 if capital<=p50_thresh_c
		gen top50_dummy_c = 1 if capital>p50_thresh_c
		gen top75_dummy_c = 1 if capital>p75_thresh_c
		gen top90_dummy_c = 1 if capital>p90_thresh_c
		gen top99_dummy_c = 1 if capital>p99_thresh_c
		
		egen total_nb_bottom50_c = sum(bottom50_dummy_c), by(gdenr_2018)
		egen total_nb_top50_c = sum(top50_dummy_c), by(gdenr_2018)
		egen total_nb_top75_c = sum(top75_dummy_c), by(gdenr_2018)
		egen total_nb_top90_c = sum(top90_dummy_c), by(gdenr_2018)
		egen total_nb_top99_c = sum(top99_dummy_c), by(gdenr_2018)

		gen total_nb_top50_to_75_c = total_nb_top50_c - total_nb_top75_c
		gen total_nb_top75_to_90_c = total_nb_top75_c - total_nb_top90_c
		gen total_nb_top90_to_99_c = total_nb_top90_c - total_nb_top99_c
		gen total_nb_top90_to_100_c = total_nb_top90_c
		gen total_nb_top99_to_100_c = total_nb_top99_c
		
		collapse (mean) ndistinct total_capital avg_capital min_capital max_capital log_avg_capital	p1_thresh_c p5_thresh_c p10_thresh_c p25_thresh_c p50_thresh_c p75_thresh_c p90_thresh_c p99_thresh_c bottom50_dummy_c top50_dummy_c top75_dummy_c top90_dummy_c top99_dummy_c total_nb_bottom50_c total_nb_top50_c total_nb_top75_c total_nb_top90_c total_nb_top99_c total_nb_top50_to_75_c total_nb_top75_to_90_c total_nb_top90_to_100_c total_nb_top90_to_99_c total_nb_top99_to_100_c year, by(gdenr_2018)
				
save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_18.dta", replace
}

* Append all years (gdenr_2018)
use "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_capital_wstat_mun_18.dta", clear
foreach y in 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
	append using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_capital_wstat_mun_18.dta"
}

save "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_1934_2003_capital_wstat_mun_18.dta", replace