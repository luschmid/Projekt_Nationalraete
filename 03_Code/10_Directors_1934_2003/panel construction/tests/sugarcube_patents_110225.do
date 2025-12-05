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

******************** 
*Regroup subsample into dataset of patents attributed in CH from 1980 to 2003
local agrp A1 A2 A3 A4 A5 A6 A7 R1 R2 R3 R4 R5 
foreach a of local agrp { 
	import delimited "C:\Users\schmutzy\Dropbox\Idea\Patent data\PATSTAT\P`a'.csv", clear
	save "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\Pat_`a'.dta", replace    
}

use "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\Pat_A1.dta", clear
local agrp A2 A3 A4 A5 A6 A7 R1 R2 R3 R4 R5 
foreach a of local agrp { 
	append using "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\Pat_`a'.dta"
}

*Standardization of the company name
gen name_temp = name
rename (name bulletin_year) (name_orig year)
order id name_orig name_temp

foreach str in name_temp {
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
	replace `str' = usubinstr(`str', "$", "s",.)
	replace `str' = usubinstr(`str', "£", " ",.)	
	replace `str' = usubinstr(`str', "^", " ",.)	
	replace `str' = usubinstr(`str', "'", " ",.)
	replace `str' = usubinstr(`str', "`", " ",.)
	replace `str' = usubinstr(`str', "´", " ",.)
	replace `str' = usubinstr(`str', "~", " ",.)	
	replace `str' = usubinstr(`str', ":", " ",.)
	replace `str' = usubinstr(`str', ";", " ",.)
	replace `str' = ustrlower(`str')
	replace `str' = ustrtrim(`str')
}

replace name_temp = usubinstr(name_temp, "    ", " ",.)
replace name_temp = usubinstr(name_temp, "   ", " ",.)
replace name_temp = usubinstr(name_temp, "  ", " ",.)

foreach str in name_temp {
	replace `str' = usubinstr(`str', " und ", "&",.)
	replace `str' = usubinstr(`str', " et ", "&",.)
	replace `str' = usubinstr(`str', "+ co", "& co",.)
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
	replace `str' = usubinstr(`str', "+", "&",.)
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
	replace `str' = usubinstr(`str', "a. g", "",.)
	replace `str' = usubinstr(`str', "s.a.", "sa",.)
	replace `str' = usubinstr(`str', "(sa).", "sa",.)
	replace `str' = usubinstr(`str', "s. a. de l'immeuble", "s.i. ",.)
	replace `str' = usubinstr(`str', "s. a. immeuble", "s.i. ",.)
	replace `str' = usubinstr(`str', "s. a. des immeubles", "s.i. ",.)
	replace `str' = usubinstr(`str', "ste anonyme des immeubles", "s.i. ",.)
	replace `str' = usubinstr(`str', "societe immobiliere", "s.i. ",.)
	replace `str' = usubinstr(`str', "ste anonyme", "sa",.)
	replace `str' = usubinstr(`str', "s. a.", "sa",.)
	replace `str' = usubinstr(`str', "s.a", "sa",.)
	replace `str' = usubinstr(`str', "sa.", "sa",.)
	replace `str' = usubinstr(`str', "schweiz.", "schweizerische",.)
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

replace name_temp = "  "+name_temp+"  "

foreach str in name_temp {
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
	replace `str' = usubinstr(`str', "akziegsellschaft", "ag",.)
	replace `str' = usubinstr(`str', "liebefeld gde", "",.)
	replace `str' = usubinstr(`str', ".", "",.)
	replace `str' = usubinstr(`str', "    ", " ",.)
	replace `str' = usubinstr(`str', "   ", " ",.)
	replace `str' = usubinstr(`str', "  ", " ",.)
	replace `str' = ustrtrim(`str')
	replace `str' = ustrlower(`str')
}


quietly bysort id year :  gen dup = cond(_N==1,0,_n)
sort name_temp year
drop if missing(name_temp)

save "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\CH_patents_1980_2003.dta", replace

*Prepare crosse-sections for matching 
use "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\CH_patents_1980_2003.dta", clear
local agrp 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {
preserve
	keep if year == `a'
	save "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\CH_patents_`a'.dta", replace
restore
}

use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\match_10_nogeo_1934_2003_ok3_temp8.dta", clear
xtset UCID_1934_2003 periode
local agrp 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {
preserve
	keep if year == `a'
	save "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`a'.dta", replace
restore
	drop if year == `a'
}

*Matching
local agrp 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {
use "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\CH_patents_`a'.dta", clear
matchit id name_temp using "$dump\02_Processed_data\10_Directors_1934_2003\26_Crosssections\Firms_`a'.dta", idu(UCID_1934_2003) txtu(cname_temp)
keep if similscore>.70
save "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\CH_patents_`a'_matched.dta"
}

***keep best 10 matches to analyze
local agrp 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {
use "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\CH_patents_`a'_matched.dta", clear

*Drop suboptimal matches in groups with perfect matches
bysort name_temp: gen newid = 1 if _n==1
replace newid = sum(newid)

gen perfect = 1 if similscore == 1
gen perfect_group = 1 if perfect == 1
forval i = 1/1200 {
bysort id (similscore):  replace perfect_group = perfect_group[_n+1] if missing(perfect_group)
}
	drop if perfect != perfect_group

*Keep only 10 best match options
	quietly bysort id (similscore) :  gen n = cond(_N==1,0,_n)
	quietly bysort id (n) :  gen n1 =_N
	keep if n > n1 - 10
	drop n n1	

*Extract perfect matches
preserve 
	keep if perfect == 1
	save "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\CH_patents_`a'_perfect_matched.dta", replace
restore
	drop if perfect == 1
	save "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\CH_patents_`a'_perfect_toeval.dta", replace
}


***Excel to evaluate matches
local agrp 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {
	use "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\CH_patents_`a'_perfect_toeval.dta", clear
	
		drop newid perfect perfect_group
		gen match = .
		sort name_temp id
		egen ID = group(id)
		bysort id UCID_1934_2003: gen dup = cond(_N==1,0,_n)
		drop if dup > 1
		
		sort name_temp id
		
	export excel using "$dump\02_Processed_data\10_Directors_1934_2003\27_Patents\CH_patents_`a'_toeval.xlsx", replace firstrow(variables)
}

/*
*** Join perfect + evaluated as correct
*Merge with initial patent data
*Extract those with no matches yet and rerun matching
local agrp 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
foreach a of local agrp {
	
}
*/