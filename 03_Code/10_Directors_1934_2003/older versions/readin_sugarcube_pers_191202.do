version 16
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"

******************************************
*** Import and Geo-Code Sugarcube Data ***
******************************************

** import from original csv files (saved as UTF-8 with BOM)
foreach val in 1934 1943 1960 1962 1963 1964 1966 1967 1970 1973 1976 1979_1980 1980_1981 1998_1999 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", delimiters(";") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta", replace
}
forv val = 1982(1)1998 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", delimiters(";") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta", replace
}
forv val = 1999(1)2003 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", delimiters(";") clear
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta", replace
}
use "$path\02_Processed_data\10_Directors_1934_2003\pers_1934_tmp.dta", clear
foreach val in 1943 1960 1962 1963 1964 1966 1967 1970 1973 1976 1979_1980 1980_1981 1998_1999 {
	append using "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta", force
	erase "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta"
}
erase "$path\02_Processed_data\10_Directors_1934_2003\pers_1934_tmp.dta"
forv val = 1982(1)2003 {
	append using "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta", force
	erase "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta"
}
save "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp1.dta", replace

* Normalize data *
use "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp1.dta", clear
rename ïid PID
label var PID "Person ID per year"
destring PID, replace
destring year, replace
save "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp2.dta", replace

/*
cd "$path\02_Processed_data\10_Directors_1934_2003\"
clear
unicode analyze *.dta
unicode encoding set Windows-1252
unicode translate *.dta, transutf8

unicode retranslate *.dta
*/

** Cleaning of special characters, etc. **
use "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp2.dta", clear
foreach str in titlejob firstname lastname address address2 alternateaddress city functions {
	replace `str' = usubinstr(`str', "Ã¤", "ae", .) // ä
	replace `str' = usubinstr(`str', "Ã", "ae", .) // Ä
	replace `str' = usubinstr(`str', "Ã ", "a", .) // à
	replace `str' = usubinstr(`str', "Ã¢", "a", .) // â
	replace `str' = usubinstr(`str', "Ã¨", "e", .) // è
	replace `str' = usubinstr(`str', "Ã©", "e", .) // é
	replace `str' = usubinstr(`str', "Ãª", "e", .) // ê
	replace `str' = usubinstr(`str', "Ã«", "e", .) // ë 
	replace `str' = usubinstr(`str', "Ã¶", "oe", .) // ö
	replace `str' = usubinstr(`str', "Ã", "oe", .) // Ö
	replace `str' = usubinstr(`str', "Ã²", "o", .) // ò
	replace `str' = usubinstr(`str', "Ã´", "o", .) // ô
	replace `str' = usubinstr(`str', "Ã¼", "ue", .) // ü
	replace `str' = usubinstr(`str', "Ã", "ue", .) // Ü
	replace `str' = usubinstr(`str', "Ã¹", "u", .) // ù
	replace `str' = usubinstr(`str', "Ã»", "u", .) // û
	replace `str' = usubinstr(`str', "Ã¯", "ï", .) // ï
	replace `str' = usubinstr(`str', "ã®", "î", .) // î
	replace `str' = usubinstr(`str', "Ã§", "c", .) // ç
	replace `str' = usubinstr(`str', "Ã", "ss", .) // ß
	replace `str' = usubinstr(`str', "ï¿½", "", .) // special character
	replace `str' = usubinstr(`str', "âº", "", .) // special character
	replace `str' = usubinstr(`str', "â", "", .) // special character
	replace `str' = usubinstr(`str', "â", "", .) // special character
	replace `str' = usubinstr(`str', "Âº", "", .) // special character
	replace `str' = usubinstr(`str', "Â", "", .) // special character
	replace `str' = usubinstr(`str', "Â", "", .) // special character
	replace `str' = usubinstr(`str', "Ã¿", "", .) // special character
	replace `str' = usubinstr(`str', "Â»", "", .) // special character «
	replace `str' = usubinstr(`str', "Â»", "", .) // special character »
	replace `str' = usubinstr(`str', "Â©", "", .) // special character
	replace `str' = usubinstr(`str', "Â°", "o", .) // °
	replace `str' = ustrlower(`str')
	replace `str' = ustrtrim(`str')
}

replace city = alternateaddress if (year>=1959 & year<=1966)  // in the 60ties the city was coded in "alternateadress" (Note that the year!=book publication year)
replace alternateaddress = "" if (year>=1959 & year<=1966)
rename postcode PLZ4 
destring PLZ4, replace
sort PLZ4 year
save "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp3.dta", replace

erase "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp1.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp2.dta"
*erase "$path\02_Processed_data\10_Directors_1934_2003\comp_1934-2003_tmp1.dta"
*erase "$path\02_Processed_data\10_Directors_1934_2003\comp_1934-2003_tmp2.dta"


** 1) Merge on PLZ: Discard information on WHEN PLZ was in use (available back until 1986)-- > take average coordinates within PLZ (over time) per GdeNr

use "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp3.dta", clear
bysort PID year: gen ID_dupl = _n  // duplicate entries w.r.t. PID year

gen gdename = city
replace gdename = alternateaddress if gdename == "" & alternateaddress != "" // (mostly) foreign location names are containted in alternateaddress

preserve
use "$path\02_Processed_data\08_Municipalities\uniquePLZ_(Avg)GdeNr2018-Geo.dta", clear
duplicates report PLZ4 
sort PLZ4 
save "$path\02_Processed_data\10_Directors_1934_2003\uPLZ_Geo.dta", replace
restore

preserve // save observation without municipality and PLZ information
keep if gdename == "" & PLZ4 == . // those that will remain unmerged
gen geo = -99
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo-99.dta", replace
restore

drop if gdename == "" & PLZ4 == .
replace gdename = usubinstr(gdename, "'", " ",.)
replace gdename = usubinstr(gdename, "`", " ",.)
replace gdename = usubinstr(gdename, "´", " ",.)
replace gdename = usubinstr(gdename, "+", " ",.)
replace gdename = usubinstr(gdename, "-", " ",.)
replace gdename = usubinstr(gdename, "$", "s",.)
replace gdename = usubinstr(gdename, "&", " ",.)
replace gdename = usubinstr(gdename, "?", " ",.)
replace gdename = usubinstr(gdename, "<", " ",.)
replace gdename = usubinstr(gdename, ">", " ",.)
replace gdename = usubinstr(gdename, "\", " ",.)
replace gdename = usubinstr(gdename , "(zuerich)" , "(zh)",.)
replace gdename = usubinstr(gdename , "(bern)" , "(be)",.)
replace gdename = usubinstr(gdename , "(berne)" , "(be)",.)
replace gdename = usubinstr(gdename , "(jura bernois)" , "(be)",.)
replace gdename = usubinstr(gdename , "(luzern)" , "(lu)",.)
replace gdename = usubinstr(gdename , "(lucerne)" , "(lu)",.)
replace gdename = usubinstr(gdename , "(uri)" , "(ur)",.)
replace gdename = usubinstr(gdename , "(schwyz)" , "(sz)",.)
replace gdename = usubinstr(gdename , "(obwalden)" , "(ow)",.)
replace gdename = usubinstr(gdename , "(nidwalden)" , "(nw)",.)
replace gdename = usubinstr(gdename , "(glarus)" , "(gl)",.)
replace gdename = usubinstr(gdename , "(zug)" , "(zg)",.)
replace gdename = usubinstr(gdename , "(fribourg)" , "(fr)",.)
replace gdename = usubinstr(gdename , "(freiburg)" , "(fr)",.)
replace gdename = usubinstr(gdename , "(solothurn)" , "(so)",.)
replace gdename = usubinstr(gdename , "(basel stadt)" , "(bs)",.)
replace gdename = usubinstr(gdename , "(basel land)" , "(bl)",.)
replace gdename = usubinstr(gdename , "(baselland)" , "(bl)",.)
replace gdename = usubinstr(gdename , "(basel landschaft)" , "(bl)",.)
replace gdename = usubinstr(gdename , "(schaffhausen)" , "(sh)",.)
replace gdename = usubinstr(gdename , "(appenzell a.r.)" , "(ar)",.)
replace gdename = usubinstr(gdename , "(appenzell a.rh.)" , "(ar)",.)
replace gdename = usubinstr(gdename , "(appenzell i.r.)" , "(ai)",.)
replace gdename = usubinstr(gdename , "(appenzell i.rh.)" , "(ai)",.)
replace gdename = usubinstr(gdename , "(appenzell a. r.)" , "(ar)",.)
replace gdename = usubinstr(gdename , "(appenzell a. rh.)" , "(ar)",.)
replace gdename = usubinstr(gdename , "(appenzell i. r.)" , "(ai)",.)
replace gdename = usubinstr(gdename , "(appenzell i. rh.)" , "(ai)",.)
replace gdename = usubinstr(gdename , "(st. gallen)" , "(sg)",.)
replace gdename = usubinstr(gdename , "(graubuenden)" , "(gr)",.)
replace gdename = usubinstr(gdename , "(aargau)" , "(ag)",.)
replace gdename = usubinstr(gdename , "(thurgau)" , "(tg)",.)
replace gdename = usubinstr(gdename , "(tessin)" , "(ti)",.)
replace gdename = usubinstr(gdename , "(ticino)" , "(ti)",.)
replace gdename = usubinstr(gdename , "(waadt)" , "(vd)",.)
replace gdename = usubinstr(gdename , "(vaud)" , "(vd)",.)
replace gdename = usubinstr(gdename , "(wallis)" , "(vs)",.)
replace gdename = usubinstr(gdename , "(valais)" , "(vs)",.)
replace gdename = usubinstr(gdename , "(neuenburg)" , "(ne)",.)
replace gdename = usubinstr(gdename , "(neuchatel)" , "(ne)",.)
replace gdename = usubinstr(gdename , "(genf)" , "(ge)",.)
replace gdename = usubinstr(gdename , "(geneve)" , "(ge)",.)
replace gdename = usubinstr(gdename , "(jura)" , "(ju)",.)
replace gdename = usubinstr(gdename, "    ", " ",.)
replace gdename = usubinstr(gdename, "   ", " ",.)
replace gdename = usubinstr(gdename, "  ", " ",.)
replace gdename = ustrtrim(gdename)
gen gdename_std1 = gdename


sort PLZ4
merge m:1 PLZ4 using "$path\02_Processed_data\10_Directors_1934_2003\uPLZ_Geo.dta"  // merge on PLZ

tab PLZ4 if PLZ4!=. & _merge==1
drop if _merge==2
rename _merge _merge_PLZ

preserve
keep if _merge_PLZ == 3
gen geo = 1
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo1.dta", replace  // save all merges from round 1 (merge on PLZ4)
restore

drop if _merge_PLZ == 3 // keep only those not merged on PLZ
drop countPLZ4 GdeNr_E_CNTR GdeNr_N_CNTR gdenr_2018
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo1_notmerged.dta", replace  // what's left to be merged


*** Merge geo-info on city names ***

// 1) merge on PLZ -- done
// 2) merge on city name & year
//    a) original name with minimum of standardization (as done above)
//    b) further standardize city names
// 3) merge on extracted PLZ (from string) 
// 4) merge on split city information (for those obs. with multiple info) & year
//    a) original name with minimum of standardization
//    b) further standardize city names excluding parentheses
//    c) delete canton info in parentheses completely, problem of duplicates in city names: loop over unique city observations
// 5) manual corrections 


** Prepare Geo-Datasets 
* Standardize municipality names in geo dataset **
use "$path\02_Processed_data\08_Municipalities\GdeNr2018_Geo_1931-2018.dta", clear
keep year ctn gdename gdenr gdenr_2018 E_CNTR N_CNTR
drop if year==. & gdename==""
duplicates report gdename year
duplicates drop gdename year, force // drop 4 observations
replace gdename = usubinstr(gdename, "Ä", "ae",.)
replace gdename = usubinstr(gdename, "Ö", "oe",.)
replace gdename = usubinstr(gdename, "Ü", "ue",.)
replace gdename = usubinstr(gdename, "ä", "ae",.)
replace gdename = usubinstr(gdename, "ö", "oe",.)
replace gdename = usubinstr(gdename, "ü", "ue",.)
replace gdename = usubinstr(gdename, "È", "e",.)
replace gdename = usubinstr(gdename, "É", "e",.)
replace gdename = usubinstr(gdename, "Ê", "e",.)
replace gdename = usubinstr(gdename, "è", "e",.)
replace gdename = usubinstr(gdename, "é", "e",.)
replace gdename = usubinstr(gdename, "ê", "e",.)
replace gdename = usubinstr(gdename, "À", "a",.)
replace gdename = usubinstr(gdename, "Â", "a",.)
replace gdename = usubinstr(gdename, "à", "a",.)
replace gdename = usubinstr(gdename, "â", "a",.)
replace gdename = usubinstr(gdename, "Û", "u",.)
replace gdename = usubinstr(gdename, "Û", "u",.)
replace gdename = usubinstr(gdename, "û", "u",.)
replace gdename = usubinstr(gdename, "ù", "u",.)
replace gdename = usubinstr(gdename, "Ô", "o",.)
replace gdename = usubinstr(gdename, "Ô", "o",.)
replace gdename = usubinstr(gdename, "ô", "o",.)
replace gdename = usubinstr(gdename, "ò", "o",.)
replace gdename = usubinstr(gdename, "Î", "i",.)
replace gdename = usubinstr(gdename, "Î", "i",.)
replace gdename = usubinstr(gdename, "î", "i",.)
replace gdename = usubinstr(gdename, "ì", "i",.)
replace gdename = usubinstr(gdename, "ç", "c",.)
replace gdename = usubinstr(gdename, "'", " ",.)
replace gdename = usubinstr(gdename, "`", " ",.)
replace gdename = usubinstr(gdename, "´", " ",.)
replace gdename = usubinstr(gdename, ",", " ",.)
replace gdename = usubinstr(gdename, "-", " ",.)
replace gdename = usubinstr(gdename, "    ", " ",.)
replace gdename = usubinstr(gdename, "   ", " ",.)
replace gdename = usubinstr(gdename, "  ", " ",.)
replace gdename = ustrlower(gdename)
replace gdename = ustrtrim(gdename)
replace ctn = ustrlower(ctn)
replace ctn = ustrtrim(ctn)
save "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp1.dta", replace // standardized gdenames

// delete ()
replace gdename = usubinstr(gdename, "(", " ",.)  
replace gdename = usubinstr(gdename, ")", " ",.)
replace gdename = usubinstr(gdename, "    ", " ",.)
replace gdename = usubinstr(gdename, "   ", " ",.)
replace gdename = usubinstr(gdename, "  ", " ",.)
replace gdename = ustrtrim(gdename)
save "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp2.dta", replace // standardized gdenames without ()

// delete "."
use "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp1.dta", clear
replace gdename = usubinstr(gdename, ".", " ",.)
replace gdename = usubinstr(gdename, "    ", " ",.)
replace gdename = usubinstr(gdename, "   ", " ",.)
replace gdename = usubinstr(gdename, "  ", " ",.)
replace gdename = ustrtrim(gdename)
save "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp3.dta", replace // standardized gdenames without "."

// delete both "." and ()
replace gdename = usubinstr(gdename, "(", " ",.)
replace gdename = usubinstr(gdename, ".", " ",.)
replace gdename = usubinstr(gdename, ")", " ",.)
replace gdename = usubinstr(gdename, "    ", " ",.)
replace gdename = usubinstr(gdename, "   ", " ",.)
replace gdename = usubinstr(gdename, "  ", " ",.)
replace gdename = ustrtrim(gdename)
save "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp4.dta", replace // standardized gdenames without "." and ()

// create one geo-dataset with different name possibilities
use "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp1.dta", clear 
forv a = 2(1)4 { 		
	append using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp`a'.dta"
	erase "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp`a'.dta"
}
erase "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp1.dta"
duplicates drop gdename year, force

sort gdename year
save "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta", replace  // dataset including gdenames with canton indication (case of potential duplicates)

// deal with duplicates when deleting canton-indication in parentheses
replace gdename = ustrregexra(gdename,"\([a-z]+\)","")  // delete canton-info in parentheses
replace gdename = usubinstr(gdename, "  ", "",.)
replace gdename = ustrtrim(gdename)

duplicates report gdename year
duplicates tag gdename year, gen(dupli)
bysort gdename year: gen ndup = _n if dupli!=0  // numbering of duplicates

forv m = 1(1)4 { // generate sub-datasets with duplicates versions
	preserve
	keep if ndup == `m'
	sort gdename year
	save "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_dupl`m'.dta", replace
	restore
}
// subset of gdename without duplicates
keep if dupli == 0
sort gdename year
save "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_dupl0.dta", replace  // dataset excluding duplicates 


* Create dataset with unique GdeNr_2018 for manual checks
use "$path\02_Processed_data\08_Municipalities\unique_Gdename_Gdenr2018-Geo.dta", clear
duplicates drop gdenr_2018, force
sort gdenr_2018
save "$path\02_Processed_data\10_Directors_1934_2003\unique_Gdenr2018-Geo.dta", replace



** String-merge with sugarcube data **

* 2a: merge on standardized city names
use "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo1_notmerged.dta", clear  // only take observations that did not merge on original PLZ4
sort gdename year
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge==2
rename _merge _merge_2a
rename E_CNTR E_CNTR_2a
rename N_CNTR N_CNTR_2a

preserve
keep if _merge_2a == 3
gen geo = 21
gen GdeNr_N_CNTR = N_CNTR_2a if geo == 21
gen GdeNr_E_CNTR = E_CNTR_2a if geo == 21
gen gdename_used = gdename_std1
keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
GdeNr_E_CNTR GdeNr_N_CNTR gdename_used gdenr_2018 geo
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo21.dta", replace   // save those who merged in 2a
restore

drop if _merge_2a==3
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo2a_notmerged.dta", replace  // save all that has not been successfully merged

* 2b: merge on further standardized city names w/o parenthesis 
replace gdename = usubinstr(gdename, "(", "",.)
replace gdename = usubinstr(gdename, ")", "",.)
replace gdename = usubinstr(gdename, "[", " ",.)
replace gdename = usubinstr(gdename, "]", " ",.)
replace gdename = usubinstr(gdename, ".", " ",.)
replace gdename = usubinstr(gdename, "/", " ",.)
replace gdename = usubinstr(gdename, ",", " ",.)
replace gdename = usubinstr(gdename, "    ", " ",.)
replace gdename = usubinstr(gdename, "   ", " ",.)
replace gdename = usubinstr(gdename, "  ", " ",.)
replace gdename = ustrregexra(gdename,"[0-9]+","")
replace gdename = ustrtrim(gdename)
gen gdename_std2 = gdename
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge==2
rename _merge _merge_2b
rename E_CNTR E_CNTR_2b
rename N_CNTR N_CNTR_2b

preserve
keep if _merge_2b == 3
gen geo = 22 
gen GdeNr_N_CNTR = N_CNTR_2b if geo == 22
gen GdeNr_E_CNTR = E_CNTR_2b if geo == 22
gen gdename_used = gdename_std2
keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
GdeNr_E_CNTR GdeNr_N_CNTR gdename_used gdenr_2018 geo
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo22.dta", replace
restore 

drop if _merge_2b==3
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo2b_notmerged.dta", replace


* 3: extract PLZ info from strings and merge on PLZ4
use "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo2b_notmerged.dta", clear
replace gdename = gdename_std1  // back to first standardization (including parenthesis with potential canton info)
replace gdename = ustrregexra(gdename,"3oll","3011")  // delete "3oll" in string (bern)
replace gdename = ustrregexra(gdename,"9oll","9011")  // delete "9oll" in string (st. gallen)
gen PLZ4_2 = regexs(1) if (regexm(gdename, "([0-9][0-9][0-9][0-9])")) // finds 4-digit numbers within string
destring PLZ4_2, replace
replace PLZ4 = PLZ4_2 if PLZ4 == . & PLZ4_2 != .
drop gdename gdenr_2018

sort PLZ4 
merge m:1 PLZ4 using "$path\02_Processed_data\10_Directors_1934_2003\uPLZ_Geo.dta"
drop if _merge == 2
rename _merge _merge_3

preserve  // check matches: check all entries on whether or not it is a PLZ4
keep if _merge_3 == 3
gen n = 1
merge m:1 gdenr_2018 using "$path\02_Processed_data\10_Directors_1934_2003\unique_Gdenr2018-Geo.dta"
keep if _merge == 3
collapse (sum) n, by(gdename_std1 gdename PLZ4)
save "$path\02_Processed_data\10_Directors_1934_2003\check_3.dta", replace
restore

// corrections 3: eliminate wrong matches
gen wrong_3 = .
replace wrong_3 = 1 if gdename_std1 == "spilistatt, 38181wald"
replace wrong_3 = 1 if gdename_std1 == "weston, mass. (usa) ciba geigy ag, basel (426, 1064)"
replace wrong_3 = 1 if gdename_std1 == "1711 cersalettes"
replace wrong_3 = 1 if gdename_std1 == "72061918"
replace wrong_3 = 1 if gdename_std1 == "ag bern, case bern manuela, (0.3)3112 allmendingen bet"
replace wrong_3 = 1 if gdename_std1 == "ag. 3933 staldenried staldenried(1.014) vp"
replace wrong_3 = 1 if gdename_std1 == "baeraustr. 43, .3550 langnau i.e, schweizerische kaeseunion ag, bern (11.563)"
replace wrong_3 = 1 if gdename_std1 == "bahnhofstr. 2540 grenchen, hotel restaurant bahnhof"
replace wrong_3 = 1 if gdename_std1 == "bern cie du chemin de fer montreux oberland bernois, montreux (3, 54169)"
replace wrong_3 = 1 if gdename_std1 == "du hotelier, glacier de 3962 la plaine montanaf rte sa, montana (2.65)"
replace wrong_3 = 1 if gdename_std1 == "eichstr. 24 comicro ag. 8045 zuerich"
replace wrong_3 = 1 if gdename_std1 == "gilbert, henri geneve oscar, rte (0.05)de 1296 cara, coppet1249 presinge"
replace wrong_3 = 1 if gdename_std1 == "grand communal, chemin du 16, 1245 collonge bellerive"
replace wrong_3 = 1 if gdename_std1 == "grand hotel regina. 3818 grindelwald"
replace wrong_3 = 1 if gdename_std1 == "herren8964 ( )"
replace wrong_3 = 1 if gdename_std1 == "hotel astoria, 7050 arosa"
replace wrong_3 = 1 if gdename_std1 == "hotel christiania. 3906 saas fee"
replace wrong_3 = 1 if gdename_std1 == "hotel krone. 3770 zweisimmen"
replace wrong_3 = 1 if gdename_std1 == "hotel krone. 9410 heiden"
replace wrong_3 = 1 if gdename_std1 == "hotel kuessnachterhof, 6403 kuessnacht"
replace wrong_3 = 1 if gdename_std1 == "hotel panorama, 3656aeschien"
replace wrong_3 = 1 if gdename_std1 == "hotel post. 7563 samnaun"
replace wrong_3 = 1 if gdename_std1 == "hotel primavera, 3901.saas grund"
replace wrong_3 = 1 if gdename_std1 == "hotel valbella jnn, 7077 valbella"
replace wrong_3 = 1 if gdename_std1 == "hotel viralago. 6574 vira (gambarogno)"
replace wrong_3 = 1 if gdename_std1 == "in firma s. denoth soehne, i7551 ramosch"
replace wrong_3 = 1 if gdename_std1 == "k clinique des grangettes kublystrasse sa, 2, chene bougeries 9016 st. gallen(1.0)"
replace wrong_3 = 1 if gdename_std1 == "kleinweld ag, 11, ennetbuergen 6330 cham(0.2)"
replace wrong_3 = 1 if gdename_std1 == "morges, prematex s.a., morges (0, 25055)"
replace wrong_3 = 1 if gdename_std1 == "obermoos. 3416 affoltern i.e, spar und leihkasse sumiswald. sumiswald (3.2)"
replace wrong_3 = 1 if gdename_std1 == "pery 2603 c/o nisus sa"
replace wrong_3 = 1 if gdename_std1 == "poststr. 26. 6300 zug grand garage acklin ag.. zug (0.4)"
replace wrong_3 = 1 if gdename_std1 == "r. general dufour 153, 2502151e1"
replace wrong_3 = 1 if gdename_std1 == "r. jardiniere 71, c accurist sa. 2300 la chaux de fonds"
replace wrong_3 = 1 if gdename_std1 == "reiterstr. 1. 4054 basel, c/o reiter haus ag, c/o reiter haus ag"
replace wrong_3 = 1 if gdename_std1 == "rustica unt. s.a, , rlnderweldstr, 19, muralto (0.17)4632 trlmbach"
replace wrong_3 = 1 if gdename_std1 == "s.i. parc du bude e, geneve (0, 1377)"
replace wrong_3 = 1 if gdename_std1 == "strassen, promo sport 2000 sa, luxembourg, lugano"
replace wrong_3 = 1 if gdename_std1 == "suzanne, dulliken 4657 (0.5)dulliken"
replace wrong_3 = 1 if gdename_std1 == "tubulaires rte sauverny jed 223, sa, carouge 1290 versoix ge"
replace wrong_3 = 1 if gdename_std1 == "ulrich, holz 3422 8965 ag, und alchenfluehberikonbern baustoffe (0, 15)ag., urdorf (0.95) vp"
replace wrong_3 = 1 if gdename_std1 == "zuerichstr. inter 135, ag, kuesnacht (0.5) 8700 kuesnacht"
replace wrong_3 = 1 if gdename_std1 == ", 8524iselisberg"
replace wrong_3 = 1 if gdename_std1 == "tigl31.10.8525, merkli jun., hansmwettingen"

replace PLZ4 = . if wrong_3 == 1
replace GdeNr_E_CNTR = . if wrong_3 == 1
replace GdeNr_N_CNTR = . if wrong_3 == 1
replace gdenr_2018 = . if wrong_3 == 1

// corrections 3: change to correct PLZ4
replace wrong_3 = 2 if gdename_std1 == "unterdorf 1405. 9042 speicher"
replace wrong_3 = 2 if gdename_std1 == "gyrixweg 95000 aarau"
replace wrong_3 = 2 if gdename_std1 == "hirzbodenweg 10/54052 basel"
replace wrong_3 = 2 if gdename_std1 == "hochspueelstr. 1683. 5722 graenichen"
replace wrong_3 = 2 if gdename_std1 == "hoehenweg 393965 berikon"
replace wrong_3 = 2 if gdename_std1 == "hofgasse 4, 41144 ãrlesheim"
replace wrong_3 = 2 if gdename_std1 == "sonnenbergweg 1115704 egliswii"
replace wrong_3 = 2 if gdename_std1 == "staatsstr. 1002. 9464 ruethi (rheintal)"
replace wrong_3 = 2 if gdename_std1 == "steinacherstr. 4, 8306lllnau effretikon"
replace wrong_3 = 2 if gdename_std1 == "steinenbergstr. 7532 feldbrunnen st. niklaus"
replace wrong_3 = 2 if gdename_std1 == "stettbachstr. 1290, b051 zuerich"
replace wrong_3 = 2 if gdename_std1 == "unterdorfstr. 1248. 9444 diepoldsau"
replace wrong_3 = 2 if gdename_std1 == "welerwies 1260, 9113.degersheim"

replace PLZ4 = 9042 if gdename_std1 == "unterdorf 1405. 9042 speicher"
replace PLZ4 = 5000 if gdename_std1 == "gyrixweg 95000 aarau"
replace PLZ4 = 4052 if gdename_std1 == "hirzbodenweg 10/54052 basel"
replace PLZ4 = 5722 if gdename_std1 == "hochspueelstr. 1683. 5722 graenichen"
replace PLZ4 = 3965 if gdename_std1 == "hoehenweg 393965 berikon"
replace PLZ4 = 4144 if gdename_std1 == "hofgasse 4, 41144 ãrlesheim"
replace PLZ4 = 5704 if gdename_std1 == "sonnenbergweg 1115704 egliswii"
replace PLZ4 = 9464 if gdename_std1 == "staatsstr. 1002. 9464 ruethi (rheintal)"
replace PLZ4 = 8308 if gdename_std1 == "steinacherstr. 4, 8306lllnau effretikon"
replace PLZ4 = 4532 if gdename_std1 == "steinenbergstr. 7532 feldbrunnen st. niklaus"
replace PLZ4 = 8051 if gdename_std1 == "stettbachstr. 1290, b051 zuerich"
replace PLZ4 = 9444 if gdename_std1 == "unterdorfstr. 1248. 9444 diepoldsau"
replace PLZ4 = 9113 if gdename_std1 == "welerwies 1260, 9113.degersheim"

preserve  // preserve the cases without necessary corrections
drop if wrong_3 == 2
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo3_tmp.dta", replace
restore

keep if wrong_3 == 2
drop countPLZ4 GdeNr_E_CNTR GdeNr_N_CNTR gdenr_2018
merge m:1 PLZ4 using "$path\02_Processed_data\10_Directors_1934_2003\uPLZ_Geo.dta"
drop if _merge == 2
drop _merge
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo3_tmp.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo3_tmp.dta"

preserve
keep if _merge_3 == 3
gen geo = 3
keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
countPLZ4 GdeNr_E_CNTR GdeNr_N_CNTR gdenr_2018 geo
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo3.dta", replace
restore

drop if _merge_3 == 3
drop countPLZ4 GdeNr_E_CNTR GdeNr_N_CNTR gdenr_2018
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo3_notmerged.dta", replace


* 4a/b/c: Split gdename to find cities in longer strings
use "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo3_notmerged.dta", clear
gen gdename = gdename_std1  // back to first standardization (including parenthesis with potential canton info)
gen nomcap = regexs(0) if regexm(gdename, "\([0-9]+\.[0-9]+\)") // extracts information in parenthesis (the info on nominal capital is represented that way)
replace gdename = ustrregexra(gdename,"([0-9][0-9][0-9][0-9])","")  // deletes 4-digit numbers coded as PLZ4_2
replace gdename = ustrregexra(gdename,"\([0-9]+\.[0-9]+\)","")  // deletes nominal capital information
replace gdename = usubinstr(gdename, "    ", " ",.)
replace gdename = usubinstr(gdename, "   ", " ",.)
replace gdename = usubinstr(gdename, "  ", " ",.)
replace gdename = ustrtrim(gdename)
split gdename, parse(,)
*local maxvar = r(nvars)  // number of variables created

drop gdename

duplicates report PID year
duplicates report PID year ID_dupl
reshape long gdename, i(PID year ID_dupl) j(n_obs_4) 

replace gdename = ustrtrim(gdename)
drop if gdename == ""
gen gdename_std4 = gdename  // preserve 4rd round standarization

save "$path\02_Processed_data\10_Directors_1934_2003\reshaped_4rd-round.dta", replace

use "$path\02_Processed_data\10_Directors_1934_2003\reshaped_4rd-round.dta", clear

// 4a: municipality names, split cells
sort gdename year
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"  // no duplicates (gdename includes canton-info)
drop if _merge==2
rename _merge _merge_4a
rename gdenr_2018 gdenr_2018_4a
rename E_CNTR E_CNTR_4a
rename N_CNTR N_CNTR_4a
gen gdename_4a = gdename

// 4b: municipality names without parentheses, split cells
replace gdename = usubinstr(gdename, "(", " ",.)
replace gdename = usubinstr(gdename, ")", " ",.)
replace gdename = usubinstr(gdename, "[", " ",.)
replace gdename = usubinstr(gdename, "]", " ",.)
replace gdename = usubinstr(gdename, "    ", " ",.)
replace gdename = usubinstr(gdename, "   ", " ",.)
replace gdename = usubinstr(gdename, "  ", " ",.)
replace gdename = ustrtrim(gdename)

sort gdename year
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"  // no duplicates (gdename includes canton-info)
drop if _merge==2
rename _merge _merge_4b
rename gdenr_2018 gdenr_2018_4b
rename E_CNTR E_CNTR_4b
rename N_CNTR N_CNTR_4b
rename gdename gdename_4b 

// 4c: municipality names without cantonal info in () and other special characters, split cells (loop through duplicates-municipal-geo datasests)
gen gdename = gdename_std4
replace gdename = ustrregexra(gdename,"\([a-z]+\)","")  // delete canton-info in parentheses
replace gdename = usubinstr(gdename, "(", " ",.)
replace gdename = usubinstr(gdename, ")", " ",.)
replace gdename = usubinstr(gdename, "/", " ",.)
replace gdename = usubinstr(gdename, "-", " ",.)
replace gdename = usubinstr(gdename, ".", " ",.)
replace gdename = usubinstr(gdename, "*", " ",.)

gen gdename_new = "" // delete all special strings, including numbers, etc.
gen length = length(gdename) 
su length, meanonly 
forval i = 1/`r(max)' { 
     local char substr(gdename, `i', 1) 
     local OK inrange(`char', "a", "z") | inlist(`char', " ")
     replace gdename_new = gdename_new + `char' if `OK' 
}
replace gdename = gdename_new
replace gdename = usubinstr(gdename, "    ", " ",.)
replace gdename = usubinstr(gdename, "   ", " ",.)
replace gdename = usubinstr(gdename, "  ", " ",.)
replace gdename = ustrtrim(gdename)
drop gdename_new

// merge on gdenames without canton info, on sequence of duplicate geo-datasets
sort gdename year
forv n = 0(1)4 {
	merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_dupl`n'.dta" // sub-datasets with duplicates
	drop if _merge==2
	rename _merge _merge_4c`n'
	rename gdenr_2018 gdenr_2018_4c`n' 
	rename E_CNTR E_CNTR_4c`n'
	rename N_CNTR N_CNTR_4c`n'
	rename ctn ctn_`n'
}

rename gdename gdename_4c
gen gdename = gdename_std4

preserve
keep if (_merge_4a == 3 | _merge_4b == 3 | ///
_merge_4c1 == 3 | _merge_4c2 == 3 | _merge_4c3 == 3 | _merge_4c4 == 3)

gen geo = 41 if _merge_4a == 3 // Step 4: use information after splitting
replace geo = 42 if geo==. & _merge_4b == 3 // Step 4: use information after splitting w\o parentheses
replace geo = 43 if geo==. & (_merge_4c1 == 3 | _merge_4c2 == 3 | _merge_4c3 == 3 | _merge_4c4 == 3) // Step 4: use information after splitting w/o cantonal info in parentheses
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo4abc.dta", replace
restore

drop if (_merge_4a == 3 | _merge_4b == 3 | ///
_merge_4c1 == 3 | _merge_4c2 == 3 | _merge_4c3 == 3 | _merge_4c4 == 3)
keep if n_obs_4 == 1  // keep only 1 obs. per PID-year-ID_dupl: in step 4, we reshape after split and create duplicate entries 
drop gdename 
gen geo = -9
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo4_notmerged.dta", replace


*** Check and Corrections ***
** Procedures
* Round 4 matches:
*     4a: Correct 1033 observations (hand-checked)
*     4b: only a few
*     4c: identify matches from this round in data, extract them, reshape long (more than one potential match) including all information, add corrections
* Missing matches with PLZ: check all cases by hand (3075 obs)
* Missing matches without PLZ: several merges & by hand
* Append all


* 4a corrections
use "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo4abc.dta", clear
keep if geo == 41 

// produce list with unique gde-info for manual checks
preserve
gen n = -1
collapse (sum) n, by (gdename_std1 gdename gdename_4a)
order n 
sort n
save "$path\02_Processed_data\10_Directors_1934_2003\check_4a.dta", replace
restore

gen wrong4a = .
replace wrong4a = 1 if gdename_std1 == "0 unterflurbetankungsanlage flughafen zuerich., zuerich" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "05)export., uster (0.5)" &  gdename_4a == "uster"
replace wrong4a = 1 if gdename_std1 == "accra, nassar group s.a., lausanne (0.5)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "adliswil, ladoco ag, zug (0.05)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "afim s.a., givrins (0.12)" &  gdename_4a == "givrins"
replace wrong4a = 1 if gdename_std1 == "ag fuer erstellung von arbeiterwohnungen, zuerich" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "ag fuer sand u. kiesverwertung nidau, nidau (0.4)" &  gdename_4a == "nidau"
replace wrong4a = 1 if gdename_std1 == "ag fuer sand u. kiesvewertung nidau, nidau (0.4)" &  gdename_4a == "nidau"
replace wrong4a = 1 if gdename_std1 == "ag heimtierbedarf en gros josti hunde katzenartikel, dietikon (0.3)" &  gdename_4a == "dietikon"
replace wrong4a = 1 if gdename_std1 == "ag, illnau effretikon (0.05)" &  gdename_4a == "illnau effretikon"
replace wrong4a = 1 if gdename_std1 == "ag., grindelwald (0.3)" &  gdename_4a == "grindelwald"
replace wrong4a = 1 if gdename_std1 == "ag., lyss (0.26)" &  gdename_4a == "lyss"
replace wrong4a = 1 if gdename_std1 == "altishofen gebr. graf ag, dagmersellen (0.12)" &  gdename_4a == "dagmersellen"
replace wrong4a = 1 if gdename_std1 == "amsterdam, succursale di lugano, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "arledhei, lonza ag, gampel (75.0)" &  gdename_4a == "gampel"
replace wrong4a = 1 if gdename_std1 == "arleshe im cusi ag, basel (2.2)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "arosa ag hotel bahnhof arosa, arosa (0.105)" &  gdename_4a == "arosa"
replace wrong4a = 1 if gdename_std1 == "ascot gb, alexander stenhouse (international) sa., geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "asphaltgeschaeft, zuerich" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "aspi, seedorf be" &  gdename_4a == "seedorf be"
replace wrong4a = 1 if gdename_std1 == "bachmann franz musica ag, luzern (0.05)" &  gdename_4a == "luzern"
replace wrong4a = 1 if gdename_std1 == "banque de depots, geneve (10.0)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "banque pour le commerce suisse israelien s.a., geneve (30.0)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "basel neptun, transport und schiffahrts ag, basel (2.4)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "basel silur ag, basel (0.05)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "basel, cobir ag, birsfelden (0.1)" &  gdename_4a == "birsfelden"
replace wrong4a = 1 if gdename_std1 == "basel, dr. r. stupp ag, rheinfelden (0.1)" &  gdename_4a == "rheinfelden"
replace wrong4a = 1 if gdename_std1 == "basel, j.f. mueller co. ag, therwil (1.0)" &  gdename_4a == "therwil"
replace wrong4a = 1 if gdename_std1 == "basel, schnyder cie ag, birsfelden (0.05)" &  gdename_4a == "birsfelden"
replace wrong4a = 1 if gdename_std1 == "bazenheid., kirchberg sg (0.05)" &  gdename_4a == "kirchberg sg"
replace wrong4a = 1 if gdename_std1 == "beirut (libanon) , unisusa sa, geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "bellevue, usa" &  gdename_4a == "bellevue"
replace wrong4a = 1 if gdename_std1 == "bellevue, wash usa" &  gdename_4a == "bellevue"
replace wrong4a = 1 if gdename_std1 == "bellevue, washington (usa)" &  gdename_4a == "bellevue"
replace wrong4a = 1 if gdename_std1 == "bellevue, washington state" &  gdename_4a == "bellevue"
replace wrong4a = 1 if gdename_std1 == "bellevue, washington usa" &  gdename_4a == "bellevue"
replace wrong4a = 1 if gdename_std1 == "bellexue, prodorite s.a., geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "belmare sa, chiasso (0.05)" &  gdename_4a == "chiasso"
replace wrong4a = 1 if gdename_std1 == "beratungs und verwaltungs ag., horw (0.2)" &  gdename_4a == "horw"
replace wrong4a = 1 if gdename_std1 == "bergamo romedi co. ag, veltlinerweine, madulain (0.12)" &  gdename_4a == "madulain"
replace wrong4a = 1 if gdename_std1 == "berger et plate s.a., geneve (0.1)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "berlin , sullmesa ag, glarus (0.5)" &  gdename_4a == "glarus"
replace wrong4a = 1 if gdename_std1 == "bern , silva ag, basel (0.05)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "bern , wilhelm und schaerz ag fuer oberflaechentechnik, zug (0.5)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "bern. kaufhaus zum erker ag, bern (0.17)" &  gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "bernische kraftwerke ag. beteiligungs ges., bern (20.0)" &  gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "biel ferba ag, bern (0.15)" &  gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "bill, eska aarau ag, aarau (0.05)" &  gdename_4a == "aarau"
replace wrong4a = 1 if gdename_std1 == "binningen, ag zum krokus, basel (0.05)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "binz gotthard immobilien ag, andermatt (0.21)" &  gdename_4a == "andermatt"
replace wrong4a = 1 if gdename_std1 == "birra bellinzona s.a., bellinzona (0.6)" &  gdename_4a == "bellinzona"
replace wrong4a = 1 if gdename_std1 == "bonn a. rh, hch. bertrams ag, basel (0.5)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "brotzer heierli ag, chur (0.05)" &  gdename_4a == "chur"
replace wrong4a = 1 if gdename_std1 == "bruessel , , sparco ag, basel (0.1)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "brugg k. ruetschi ag pumpenbau brugg, brugg (2.0)" &  gdename_4a == "brugg"
replace wrong4a = 1 if gdename_std1 == "brugg, buchdruckerei werde r ag, windisch (0.18)" &  gdename_4a == "windisch"
replace wrong4a = 1 if gdename_std1 == "bruxelles , protector s.a., lucens (0.8)" &  gdename_4a == "lucens"
replace wrong4a = 1 if gdename_std1 == "bruxelles b, bank fuer internationalen zahlungsausgleich, basel" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "bruxelles b, bank fuer internationalen zahlungsausgleich, basel (1500.0)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "bruxelles b, bankfuer internationalen zahlungsausgleich, basel" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "bubendorf bella lui s.a., montana (0.05)" &  gdename_4a == "montana"
replace wrong4a = 1 if gdename_std1 == "buergstadt a. m, sand ag neuheim, neuheim (1.5)" &  gdename_4a == "neuheim"
replace wrong4a = 1 if gdename_std1 == "buergstadt a. m, sand ag, schwyz (1.5)" &  gdename_4a == "schwyz"
replace wrong4a = 1 if gdename_std1 == "bueroeinrichtung und planung ag, zug" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "bueron arnold et cie ag, bueron (0.12)" &  gdename_4a == "bueron"
replace wrong4a = 1 if gdename_std1 == "burgstadt a. m, sand ag, schwyz (1.5)" &  gdename_4a == "schwyz"
replace wrong4a = 1 if gdename_std1 == "buss ag, basel (3.6)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "c cotram sa, collombey muraz" &  gdename_4a == "collombey muraz"
replace wrong4a = 1 if gdename_std1 == "c o merrill, lynch, p., f., r. rhone 65, geneve" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "c o roner ag, bern, bennenboden" &  gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o hewaplastic ag, luzern" &  gdename_4a == "luzern"
replace wrong4a = 1 if gdename_std1 == "ca dei galli, besazio" &  gdename_4a == "besazio"
replace wrong4a = 1 if gdename_std1 == "campbell s soups s.a., geneve (2.0)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "cayman ind, asg arlesheimer sporthallen a g, arlesheim (1.3)" &  gdename_4a == "arlesheim"
replace wrong4a = 1 if gdename_std1 == "cayman ind, asg arlesheimer sporthallen ag, arlesheim (0.26)" &  gdename_4a == "arlesheim"
replace wrong4a = 1 if gdename_std1 == "cayman ind, asg arlesheimer sporthallen ag., arlesheim (0.26)" &  gdename_4a == "arlesheim"
replace wrong4a = 1 if gdename_std1 == "cayman ind, asg arlesheimer sporthallen aktiengesellschaft, arlesheim" &  gdename_4a == "arlesheim"
replace wrong4a = 1 if gdename_std1 == "cayman ind, curlingzentruin arlesheim ag, arlesheim (0.26)" &  gdename_4a == "arlesheim"
replace wrong4a = 1 if gdename_std1 == "cayman ind, curlingzentrum arlesheim ag, arlesheim (0.26)" &  gdename_4a == "arlesheim"
replace wrong4a = 1 if gdename_std1 == "charpentes metalliques sa, sion" &  gdename_4a == "sion"
replace wrong4a = 1 if gdename_std1 == "collonge bellerive, garage desjacques s.a., geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "cologny commodities trading company s.a., geneve (0.15)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "compagnie d assurances, lausanne (60.0)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "compagnie d assurances, lausanne (7.0)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "conte di sane, geneve" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "conte di sarre, geneve" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "corsier la genevoise, cie d assurances sur la vie, geneve (12.0)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "crema s.i. silver s.a., lugano (0.06)" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "crocifisso di savosi radio citta t.v. s.a., vezia (0.05)" &  gdename_4a == "vezia"
replace wrong4a = 1 if gdename_std1 == "curacao n.v, tcpl finance nv, curacao, zug branch, zug" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "curacao na, addax trading management nv, curacao, succursale de geneve, geneve" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "dillenburg (d) ., erismann ag, neunkirch (1.0)" &  gdename_4a == "neunkirch"
replace wrong4a = 1 if gdename_std1 == "dorf, st. etienne (f)" &  gdename_4a == "dorf"
replace wrong4a = 1 if gdename_std1 == "dorf, ursenbach" &  gdename_4a == "dorf"
replace wrong4a = 1 if gdename_std1 == "dr. med, a. w. graf ag weberei, illnau (0.2)" &  gdename_4a == "illnau"
replace wrong4a = 1 if gdename_std1 == "dr., ag, duggingen (0.05)" &  gdename_4a == "duggingen"
replace wrong4a = 1 if gdename_std1 == "du hotelier, glacier de 3962 la plaine montanaf rte sa, montana (2.65)" &  gdename_4a == "montana"
replace wrong4a = 1 if gdename_std1 == "duesseldorf , maschinenfabrik horgen ag in horgen, horgen (0.1)" &  gdename_4a == "horgen"
replace wrong4a = 1 if gdename_std1 == "ebnat, schuhhaus brunner ag, wattwil (0.15)" &  gdename_4a == "wattwil"
replace wrong4a = 1 if gdename_std1 == "eduard bolli ag, schaffhausen (0.8)" &  gdename_4a == "schaffhausen"
replace wrong4a = 1 if gdename_std1 == "elan ag., bern (0.15)" &  gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "eldima ag, zuerich (0.6)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "elektrodol ag, ebikon (0.05)" &  gdename_4a == "ebikon"
replace wrong4a = 1 if gdename_std1 == "etablissement pre fleuti sauvabelin s.a., lausanne (0.3)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "fafnir roulements international s.a., fribourg" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "fafnir roulements international s.a., fribourg (8.05)" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "ferlor ag, zug (0.3)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "fib finanzierungs immobilien und beteiligungs ag, herzogenbuchsee (0.1)" &  gdename_4a == "herzogenbuchsee"
replace wrong4a = 1 if gdename_std1 == "fils, pully, s.i. de rouvenoz, lausanne (0.05)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "flawil , stahlbau und montagen ag altwil sg, gaiserwald" &  gdename_4a == "gaiserwald"
replace wrong4a = 1 if gdename_std1 == "fluids ag, zollikon (0.1)" &  gdename_4a == "zollikon"
replace wrong4a = 1 if gdename_std1 == "forest row, bp benzin petroleum ag, zuerich (20.0)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frackfurt a m, phoenix armaturen ag, basel (0.1)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt a m, electromation ag, basel (0.75)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt a m, electromation ag, basel (5.0)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt a m, electromation ag, basel (5.45)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt a m, phoenix armaturen ag, basel (0.1)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, ameropa immobilien investment ag, zuerich (1.0)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, brehm s beteiligungs ag brebag, zuerich (0.5)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, contrasmog ag, zuerich (0.05)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, dr. leistritz schoop ag, baden (0.25)" &  gdename_4a == "baden"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, farbwerke hoechst investment ag, zuerich (2.0)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, hoechst pharma ag, zuerich (0.5)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, intermit ag fuer internationale miteigentumswerte, zollikon (1.0)" &  gdename_4a == "zollikon"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, internationale genossenschaftsbank ag, basel (10.0)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, panagrar ag, glarus (0.5)" &  gdename_4a == "glarus"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, phoenix atmaturcn ag, basel (0.1)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, selmi ag, zuerich (0.13)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frankfurt a. m, sperry rand international corporation, lausanne" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, â«interfondsâ» internationale investmenttrust gesellschaft, basel (1.0)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, â«interfondsâ» internationale investmenttrust gesetschaft, basel (1.0)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, brehm s beteiligungs ag brebag, zuerich (0.5)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, cioccolata pretty sa, caslano (0.05)" &  gdename_4a == "caslano"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, cioccolata pretty sa, lugano (0.05)" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, cioccolata titlis s.a., caslano (0.65)" &  gdename_4a == "caslano"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, control data holding ag, luzern (1.3)" &  gdename_4a == "luzern"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, creda s.a. etablissement de credit, meyrin (0.05)" &  gdename_4a == "meyrin"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, dr leistritz schoop ag, baden (0.25)" &  gdename_4a == "baden"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, elbag elektronik beteiligungs ag, zuerich (0.05)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, farbwerke hoechst investment ag, zuerich (2.0)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, maedler ag, zuerich (0.1)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, maedler holding ag, zug (0.1)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, s.i. les verjus, geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, schmetzer ag, basel (0.2)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, schoyca holding s.a., lugano (0.05)" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, selbst waschautomaten betrieb ag, geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, tenigerbad ag, wuerenlos (0.33)" &  gdename_4a == "wuerenlos"
replace wrong4a = 1 if gdename_std1 == "frankfurt a.m, vacanza holding ag, basel (0.05)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "frankfurt thermalbad und hotel betriebs ag bad vals, vals" &  gdename_4a == "vals"
replace wrong4a = 1 if gdename_std1 == "frauenfeld, espa immobilien ag, frauenfeld (0.5)" &  gdename_4a == "frauenfeld"
replace wrong4a = 1 if gdename_std1 == "freiburg 1.br, eden ag, chur (0.225)" &  gdename_4a == "chur"
replace wrong4a = 1 if gdename_std1 == "freiburg i. b, elektrizitaetswerk rheinau ag, rheinau (20.0)" &  gdename_4a == "rheinau"
replace wrong4a = 1 if gdename_std1 == "freiburg i. br, eden ag 1 ierrenwueschefabrik, chur (0.225)" &  gdename_4a == "chur"
replace wrong4a = 1 if gdename_std1 == "freiburg i. br, eden ag herrenwaeschefabrik, chur (0.225)" &  gdename_4a == "chur"
replace wrong4a = 1 if gdename_std1 == "freiburg i. br, radio und fernseh entwicklungsges. ag, zug (0.05)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "freiburg i. br, tulanum ag, zuerich (0.3)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "freiburg i.b, elektrizitaetswerk rheinau ag, rheinau (20.0)" &  gdename_4a == "rheinau"
replace wrong4a = 1 if gdename_std1 == "freiburg i.b, herbert ammann co ag, riehen (0.1)" &  gdename_4a == "riehen"
replace wrong4a = 1 if gdename_std1 == "freiburg i.b, isobale ag, basel (0.1)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "freiburg i.b, multi contact ag, basel (0.42)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "funchal, succursale di lugano, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "garrett international s.a., geneve (0.16)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "gebr. kuenzli ag kunstverlag, zuerich (0.5)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "genese, atlas tea company, geneve (0.2)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "geneve , brown fintube s.a., fribourg (0.15)" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "geneve , cinema d echallens s.a., echallens (0.05)" &  gdename_4a == "echallens"
replace wrong4a = 1 if gdename_std1 == "geneve , ludo s.a., carabbia (0.05)" &  gdename_4a == "carabbia"
replace wrong4a = 1 if gdename_std1 == "geneve., gallet et co s.a., la chaux de fonds (0.2)" &  gdename_4a == "la chaux de fonds"
replace wrong4a = 1 if gdename_std1 == "gepla s.a. des metaux plaques, geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "gerald dambach societe anonyme., tolochenaz (0.2)" &  gdename_4a == "tolochenaz"
replace wrong4a = 1 if gdename_std1 == "giessen a.l, cordial ag., zuerich (0.1)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "gladbeck i.w, ems gelsenberg ag, domat/ems (15.0)" &  gdename_4a == "domat/ems"
replace wrong4a = 1 if gdename_std1 == "gladbeck i.w, ems gelsenberg ag, domat/ems (25.0)" &  gdename_4a == "domat/ems"
replace wrong4a = 1 if gdename_std1 == "goldbach adrema ag, zuerich (0.2)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "gordola , , s.a. aerocentro ticinese, locarno (0.05)" &  gdename_4a == "locarno"
replace wrong4a = 1 if gdename_std1 == "grafinma commerce sa, geneve (0.1)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "guldenmann ad, basel (0.2)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "h.c. fehr blockfioetenbau ag, zuerich (0.1)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "h.c. fehr blockfloetenbau ag, zuerich (0.1)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "haegglingen, m. geissmann co. ag, habsburg (0.1)" &  gdename_4a == "habsburg"
replace wrong4a = 1 if gdename_std1 == "hamilton., stelco sa, geneve (0.5)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "hans scherrer erben ag, opfikon (3.15)" &  gdename_4a == "opfikon"
replace wrong4a = 1 if gdename_std1 == "hans u. bon ag, zuerich (0.2)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "harrison n.y, federal pacific electric overseas s.a., zug (4.0)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "haslum (n) noralu a.g., duebendorf (0.1)" &  gdename_4a == "duebendorf"
replace wrong4a = 1 if gdename_std1 == "haury co. ag, st gallen (0.55)" &  gdename_4a == "st gallen"
replace wrong4a = 1 if gdename_std1 == "hausen stocker bau ag, brugg (0.2)" &  gdename_4a == "brugg"
replace wrong4a = 1 if gdename_std1 == "heberlein co ag, wattwil (3.0)" &  gdename_4a == "wattwil"
replace wrong4a = 1 if gdename_std1 == "heid en eskapol ag, heiden (0.05)" &  gdename_4a == "heiden"
replace wrong4a = 1 if gdename_std1 == "helsinki/fins, walter kadler ag, zuerich (0.05)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "henau, investment handels ag, zollikon (0.05)" &  gdename_4a == "zollikon"
replace wrong4a = 1 if gdename_std1 == "herfeld ag, stein am rhein (0.25)" &  gdename_4a == "stein am rhein"
replace wrong4a = 1 if gdename_std1 == "hessigkof en uhrenfabrik wega ag, grenchen (0.2)" &  gdename_4a == "grenchen"
replace wrong4a = 1 if gdename_std1 == "hoershoim daen, nielsen, otto, emballage ag, oetwil am see (0.05)" &  gdename_4a == "oetwil am see"
replace wrong4a = 1 if gdename_std1 == "hoersholm daen, nielsen, otto, emballage ag, oetwil am see (0.05)" &  gdename_4a == "oetwil am see"
replace wrong4a = 1 if gdename_std1 == "hoersholm daen, nielsen, otto, emballage ag., oetwil am see (0.05)" &  gdename_4a == "oetwil am see"
replace wrong4a = 1 if gdename_std1 == "hoersholm daen, nielsen, otto. emballage ag., oetwil am see (0.05)" &  gdename_4a == "oetwil am see"
replace wrong4a = 1 if gdename_std1 == "hohenrain, 8564waldi" &  gdename_4a == "hohenrain"
replace wrong4a = 1 if gdename_std1 == "holz und chaletbau ag, leissigen" &  gdename_4a == "leissigen"
replace wrong4a = 1 if gdename_std1 == "homburg, saar rfa" &  gdename_4a == "homburg"
replace wrong4a = 1 if gdename_std1 == "hong kong hk, progin poissons sa, neuchatel (0.05)" &  gdename_4a == "neuchatel"
replace wrong4a = 1 if gdename_std1 == "horlogere sa., la chaux de fonds (0.051)" &  gdename_4a == "la chaux de fonds"
replace wrong4a = 1 if gdename_std1 == "horw, temp ag, zug (0.2)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "howald lux schneider partner sa, luxembourg, zweigniederlassung zuerich, zuerich" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "ile de sark guernesey, clearton ltd, dublino, succursale di lugano, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "imfela ag, altstaetten (0.1)" &  gdename_4a == "altstaetten"
replace wrong4a = 1 if gdename_std1 == "irnfela ag, altstaetten (0.1)" &  gdename_4a == "altstaetten"
replace wrong4a = 1 if gdename_std1 == "k ericsson ag, wangen bruettisellen (0.1)" &  gdename_4a == "wangen bruettisellen"
replace wrong4a = 1 if gdename_std1 == "k. s.a. geneverenens, renens (0.05) miche henri albert, geneve" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "kbaktrade ltd, kuesnacht zh (0.1)" &  gdename_4a == "kuesnacht zh"
replace wrong4a = 1 if gdename_std1 == "kbush ag, zuerich (0.05)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "kenwood schumpf ag, baar (0.3)" &  gdename_4a == "baar"
replace wrong4a = 1 if gdename_std1 == "kieswerk frauholz, steinen" &  gdename_4a == "steinen"
replace wrong4a = 1 if gdename_std1 == "kiew uk ukraine international airlines, kiew, zweigniederlassung kloten, kloten" &  gdename_4a == "kloten"
replace wrong4a = 1 if gdename_std1 == "koganei/jap, nippon yakin kogyo sa, geneve, geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "kumag ag maschinenfabrik, zuerich (0.3)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "kurrag ag maschinenfabrik, zuerich" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "la chaux de fonds j. wyss s.a., neuchatel (0.05)" &  gdename_4a == "neuchatel"
replace wrong4a = 1 if gdename_std1 == "lanark gb viminex sa, geneve (2.0)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "langen a m, phoenix armaturen ag, basel (0.1)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "lantac ag, zug (0.05)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "lausanne, en plan s.a. s.i., morges (0.05)" &  gdename_4a == "morges"
replace wrong4a = 1 if gdename_std1 == "lausanne, fabrique de machines et d articles plastiques s.a., vevey (0.05)" &  gdename_4a == "vevey"
replace wrong4a = 1 if gdename_std1 == "lausanne, s.i. de la place ronjat, vevey (0.05)" &  gdename_4a == "vevey"
replace wrong4a = 1 if gdename_std1 == "le franc montagnard, saignelegier (0.112)" &  gdename_4a == "saignelegier"
replace wrong4a = 1 if gdename_std1 == "lettre b, geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "lettre f, geneve" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "liestal, anfos immobilien ag, basel (0.05)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "litemaster ag, c o erwin c.kasp, zug" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "lofty detfaux, fribourg" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "london , southern musik ag, volketswil (0.05)" &  gdename_4a == "volketswil"
replace wrong4a = 1 if gdename_std1 == "london gb, marks spencer (basel) ag, basel (0.05)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "london, kingsfire productions ag, chur (0.05)" &  gdename_4a == "chur"
replace wrong4a = 1 if gdename_std1 == "london, succursale de geneve, geneve" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "londra succursale di lugano, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "londra, succursale di lugano, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "londres gb k 1wm independent wealth management sa, geneve" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "londres gb, handelsfinanz midland bank., geneve (40.0)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "lubag ag, luzern (0.05)" &  gdename_4a == "luzern"
replace wrong4a = 1 if gdename_std1 == "lugano, brossa ag, st. gallen (0.05)" &  gdename_4a == "st. gallen"
replace wrong4a = 1 if gdename_std1 == "lugano, lnpraegnierwerke brittnau/wikon ag, brittnau" &  gdename_4a == "brittnau"
replace wrong4a = 1 if gdename_std1 == "lugano, rohimpag roh materialien ag, basel (0.05)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "luxembourg, kiwi ii management company sa, luxembourg, lugano branch, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "luzerg, gebr. macchi ag automobile, dietlikon (0.06)" &  gdename_4a == "dietlikon"
replace wrong4a = 1 if gdename_std1 == "magnetic elektromotoren ag, liestal (0.3)" &  gdename_4a == "liestal"
replace wrong4a = 1 if gdename_std1 == "magnetic elektromotoren ag, liestal (0.5)" &  gdename_4a == "liestal"
replace wrong4a = 1 if gdename_std1 == "maidstone gb, fidelity investments international, hildenborough, zurich branch, zuerich" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "maria enzersdorf a, oberflaechentechnik klus ag, balsthal (0.15)" &  gdename_4a == "balsthal"
replace wrong4a = 1 if gdename_std1 == "marseille, aerodrome regional de montreux s.a., rennaz (1.5)" &  gdename_4a == "rennaz"
replace wrong4a = 1 if gdename_std1 == "martigny bourg , walter fr. moser, patent service s.a., geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "maschinenfabrik ag, kuesnacht zh (0.39)" &  gdename_4a == "kuesnacht zh"
replace wrong4a = 1 if gdename_std1 == "matters. , kuehlhaus luzern, luzern (0.575)" &  gdename_4a == "luzern"
replace wrong4a = 1 if gdename_std1 == "mauren fl, succursale di lugano, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "medidenta ag, st. gallen (0.6)" &  gdename_4a == "st. gallen"
replace wrong4a = 1 if gdename_std1 == "meggen, cheddite plastic s.a., liestal (2.0)" &  gdename_4a == "liestal"
replace wrong4a = 1 if gdename_std1 == "meggen, progrebras ag, basel (0.6)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "meilen, familag herisau, holding ag, herisau (2.0)" &  gdename_4a == "herisau"
replace wrong4a = 1 if gdename_std1 == "meli de longavia ag, lugano (0.05)" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "melle b , frigosuisse ag, basel (1.3)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "metzger spina, churwalden" &  gdename_4a == "churwalden"
replace wrong4a = 1 if gdename_std1 == "midland u.s.a, dow chemical ag, zuerich (420.0)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "mies , ste des eaux commugny mies, commugny (0.06)" &  gdename_4a == "commugny"
replace wrong4a = 1 if gdename_std1 == "milan, alfa holding company s.a., geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "milano , lucilla finanziaria sa., bellinzona (0.05)" &  gdename_4a == "bellinzona"
replace wrong4a = 1 if gdename_std1 == "minitax ag, zuerich (0.1)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "moebel neuhof ag, winterthur (0.15)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "montreal , chemcell sa, zug (0.4)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "montreux, pres" &  gdename_4a == "montreux"
replace wrong4a = 1 if gdename_std1 == "moulins rod s.a., orbe (0.95)" &  gdename_4a == "orbe"
replace wrong4a = 1 if gdename_std1 == "muell schob ) partner ag, zuerich (0.3)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "muralto , tipografia alla motta s.a., locarno (0.3)" &  gdename_4a == "locarno"
replace wrong4a = 1 if gdename_std1 == "muttenz tacoffor ag, basel (0.05)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "n. pedolin s erben ag, chur (0.4)" &  gdename_4a == "chur"
replace wrong4a = 1 if gdename_std1 == "neuff en esge electric s.a., glarus (0.05)" &  gdename_4a == "glarus"
replace wrong4a = 1 if gdename_std1 == "neuheim, 627513allwil" &  gdename_4a == "neuheim"
replace wrong4a = 1 if gdename_std1 == "neuhhausen a.r, lega ag, bern (0.09)" &  gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "neuilly sur seine , saint gobain international, fribourg (230.0)" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "new trend sa, fribourg (0.1)" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "new york u.s.a, ag fuer verwaltung von investment trusts (intrag), zuerich (3.0)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "new york u.s.a, intrag ag verwaltung von investmenttrusts, zuerich" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "new york, succursale di lugano, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "normatrol ag, oberrieden (0.1)" &  gdename_4a == "oberrieden"
replace wrong4a = 1 if gdename_std1 == "novastampa sa, barbengo (0.05)" &  gdename_4a == "barbengo"
replace wrong4a = 1 if gdename_std1 == "o holding ag., zug (4.0)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "oberkirch, d" &  gdename_4a == "oberkirch"
replace wrong4a = 1 if gdename_std1 == "oberstdorf , immo ag neugasse st. gallen, st. gallen (0.53)" &  gdename_4a == "st. gallen"
replace wrong4a = 1 if gdename_std1 == "obstverwertung beromuenster, beromuenster (0.15)" &  gdename_4a == "beromuenster"
replace wrong4a = 1 if gdename_std1 == "oceline, zuerich" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "one), s.i. ekeko, geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "opfikon westa ag liegenschaften , herrliberg (0.05)" &  gdename_4a == "herrliberg"
replace wrong4a = 1 if gdename_std1 == "ordrup. , hefefabriken ag, olten (1.0)" &  gdename_4a == "olten"
replace wrong4a = 1 if gdename_std1 == "orgelbau maag ag, zuerich (0.12)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "oslo, incinerator ag, zug (0.05)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "ostenda b, cis distribution ltd, dublino, succursale di lugano, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "p. stazione, domus rustica sa, muralto" &  gdename_4a == "muralto"
replace wrong4a = 1 if gdename_std1 == "palo alto, friden international sa, fribourg (4.3)" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "panama, lugano branch, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "panissieres f, sefar ag, rueschlikon (4.8)" &  gdename_4a == "rueschlikon"
replace wrong4a = 1 if gdename_std1 == "panissieres f, sefar ag, rueschlikon (7.5)" &  gdename_4a == "rueschlikon"
replace wrong4a = 1 if gdename_std1 == "panissieres f, sefar holding ag, rueschlikon (7.5)" &  gdename_4a == "rueschlikon"
replace wrong4a = 1 if gdename_std1 == "panissleres f, sefar ag, rueschlikon (4.8)" &  gdename_4a == "rueschlikon"
replace wrong4a = 1 if gdename_std1 == "papierfabrik albert ziegler ag, grellingen (4.5)" &  gdename_4a == "grellingen"
replace wrong4a = 1 if gdename_std1 == "paris banque de paris et des pays bas (suisse) sa, geneve (85.0)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "paris glantex s.a., gland (0.5)" &  gdename_4a == "gland"
replace wrong4a = 1 if gdename_std1 == "paris resist s.a., fabrique de ressorts, neuchatel (0.05)" &  gdename_4a == "neuchatel"
replace wrong4a = 1 if gdename_std1 == "paris/f, ste continentale de gestion financieres sa, lausanne" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "passavant iselin cie. ag, allschwil(1.29)" &  gdename_4a == "allschwil"
replace wrong4a = 1 if gdename_std1 == "peabody usa, 1, electrographic sa, lausanne (0.05)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "perrin, forster crosby inc., philadelphie, meyrinmeyrin branch, meyrin" &  gdename_4a == "meyrin"
replace wrong4a = 1 if gdename_std1 == "philips ag, zuerich (20.0)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "piaseczno pl, azyxsa, fribourg (0.1)" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "pini associati, ingegneri consulenti sa, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "plan les ouates , labiol s.a., geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "playa s.juan, campello, alicante e" &  gdename_4a == "campello"
replace wrong4a = 1 if gdename_std1 == "prag cz , czechlmpex ag, zuerich (0.1)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "pratteln, post augst" &  gdename_4a == "pratteln"
replace wrong4a = 1 if gdename_std1 == "preverenges. s.i. valency champrilly a, lausanne (0.05)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "publicitas sa, suisse de publicite, lausanne (12.0)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "pull, j. p. schmidt ag holzgeschaeft, chur (0.2)" &  gdename_4a == "chur"
replace wrong4a = 1 if gdename_std1 == "pulli, j. p. schmidt ag holzgeschaeft, chur (0.2)" &  gdename_4a == "chur"
replace wrong4a = 1 if gdename_std1 == "pulls, s.i. les tilleuls s.a., yverdon (0.05)" &  gdename_4a == "yverdon"
replace wrong4a = 1 if gdename_std1 == "pulls, savary s.a., lausanne (0.05)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "racine, racine und vickers armstrong ag, zug (1.25)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "rau ag, st. gallen (0.15)" &  gdename_4a == "st. gallen"
replace wrong4a = 1 if gdename_std1 == "raytheon ag, zug (0.5)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "rovio, g. tacchella s.a., arogno (0.1)" &  gdename_4a == "arogno"
replace wrong4a = 1 if gdename_std1 == "s.a., genevelettre b, geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "s.a., lausanne (0.05)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "s.diego cal, sulzer energieconsulting ag, winterthur (0.05)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "s.diego cal, sulzer energieconsulting ag., winterthur (0.5)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "s.diego cal, sulzer general atomic waermetauscher ag, winterthur (0.05)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "s.dlego cal, sulzer energieconsulting ag, winterthur (0.05)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "s.dlego cal, sulzer energieconsulting ag, winterthur (0.5)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "s.dlego cal, sulzer energieconsulting ag., winterthur (0.5)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "s.dtego cal, sulzer energieconsulting ag, winterthur (0.500)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "s.i. levant bellevue i, lausanne (0.06)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "sa fribourg, fribourg (0.15)" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "sa, geneve (50.0)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "sa., gordola (0.4)" &  gdename_4a == "gordola"
replace wrong4a = 1 if gdename_std1 == "saas i.pr, klosters madrisa bergbahnen ag, klosters (3.0)" &  gdename_4a == "klosters"
replace wrong4a = 1 if gdename_std1 == "san diego cal, pgmm schweiz ag, winterthur (0.5)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "san diego cal, sulzer energieconsulting ag, winterthur (0.5)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "sarpa ag, zug (0.05)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "schagen nl inter pharm holding ag, kreuzlingen (0.05)" &  gdename_4a == "kreuzlingen"
replace wrong4a = 1 if gdename_std1 == "schwanenhaus immobilien ag, goldach (0.099)" &  gdename_4a == "goldach"
replace wrong4a = 1 if gdename_std1 == "schweizerische suedostbahn, waedenswil (2.25)" &  gdename_4a == "waedenswil"
replace wrong4a = 1 if gdename_std1 == "seiler co. ag, gelterkinden (0.75)" &  gdename_4a == "gelterkinden"
replace wrong4a = 1 if gdename_std1 == "selzach scierie boillat s.a., les breuleux (0.15)" &  gdename_4a == "les breuleux"
replace wrong4a = 1 if gdename_std1 == "senningen l zamasport distribution sa, lussemburgo, succursale di chiasso, chiasso" &  gdename_4a == "chiasso"
replace wrong4a = 1 if gdename_std1 == "settelcn ag, basel (0.3)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "settelen ag, basel (0.3)" &  gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "sibromo ag, bern (0.1)" &  gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "sidermetal ltd, douglas, succursale di lugano, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "sils i. d, satzelektronik ag, chur (0.06)" &  gdename_4a == "chur"
replace wrong4a = 1 if gdename_std1 == "sion jean vicarino meyer s.a., fribourg (0.1)" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "sperry rand international corporation, lausanne" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "sporthaus, grindelwald" &  gdename_4a == "grindelwald"
replace wrong4a = 1 if gdename_std1 == "st cloud f, cits sa, geneve (0.2)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "st jean d aulps f, mcm partners, lyon, succursale de geneve, geneve" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "st martin, ile jersey (gb)" &  gdename_4a == "st martin"
replace wrong4a = 1 if gdename_std1 == "st martin, jersey gb" &  gdename_4a == "st martin"
replace wrong4a = 1 if gdename_std1 == "st martin, lie jersey gb" &  gdename_4a == "st martin"
replace wrong4a = 1 if gdename_std1 == "st. gallisch appenzellische kraftwerke ag, st. gallen" &  gdename_4a == "st. gallen"
replace wrong4a = 1 if gdename_std1 == "st. georgen strasse 191b, 9oll st. kreditallenanstal grabs, grabs (5.0)" &  gdename_4a == "grabs"
replace wrong4a = 1 if gdename_std1 == "st. martin, jersey gb" &  gdename_4a == "st. martin"
replace wrong4a = 1 if gdename_std1 == "st. martin, uk" &  gdename_4a == "st. martin"
replace wrong4a = 1 if gdename_std1 == "standford calif, ullberg s.a., neuchatel (0.05)" &  gdename_4a == "neuchatel"
replace wrong4a = 1 if gdename_std1 == "stanztechnik, bettlach" &  gdename_4a == "bettlach"
replace wrong4a = 1 if gdename_std1 == "steinen, kreis loerrach d" &  gdename_4a == "steinen"
replace wrong4a = 1 if gdename_std1 == "stockholm, hallberg adhesive products s.a., lausanne (0.3)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "strasbourg (f), cofrastra ag, bern (0.1)" &  gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "strassen, promo sport 2000 sa, luxembourg, lugano" &  gdename_4a == "lugano"
replace wrong4a = 1 if gdename_std1 == "stuttgart ., prolibro ag, zug (0.06)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "succursale de lutry, lutry" &  gdename_4a == "lutry"
replace wrong4a = 1 if gdename_std1 == "succursale de sembrancher, sembrancher" &  gdename_4a == "sembrancher"
replace wrong4a = 1 if gdename_std1 == "sur i.o, kieswerk mulegns ag, mulegns (0.14)" &  gdename_4a == "mulegns"
replace wrong4a = 1 if gdename_std1 == "sutton, reservations et investissements lausanne s.a., lausanne (0.05)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "taipei (taiwan. r.o.c.), china airlines, ltd, taipei, taiwan, r.o.c., zweigniederlassung glattbrugg, opfikon" &  gdename_4a == "opfikon"
replace wrong4a = 1 if gdename_std1 == "tehag techno hydraulik ag, zuerich (0.1)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "teheran ir micronet systems sa, lausanne (0.1)" &  gdename_4a == "lausanne"
replace wrong4a = 1 if gdename_std1 == "teppichhaus j. oetterli ag, luzern (0.15)" &  gdename_4a == "luzern"
replace wrong4a = 1 if gdename_std1 == "teuf en allwell express food inc., st gallen (0.1)" &  gdename_4a == "st gallen"
replace wrong4a = 1 if gdename_std1 == "thoerishaus charles kaufmann fribourg s.a., fribourg (0.175)" &  gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "thun, baustoffe bern ag, bern (0.4)" &  gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "tokyo j abb asea brown boveri ltd, zuerich (2768.0)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "torre, sta. bancaria ticinese, bellinzona (1.0)" &  gdename_4a == "bellinzona"
replace wrong4a = 1 if gdename_std1 == "transtour invest, meyrin (0.1)" &  gdename_4a == "meyrin"
replace wrong4a = 1 if gdename_std1 == "trend solutions ag, moeriken wildegg" &  gdename_4a == "moeriken wildegg"
replace wrong4a = 1 if gdename_std1 == "tretorn ag, zuerich (0.1)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "treuhand und verwaltungs ag., chur, (0.1)" &  gdename_4a == "chur"
replace wrong4a = 1 if gdename_std1 == "turin , cotra autotransport ag, rothrist (1.0)" &  gdename_4a == "rothrist"
replace wrong4a = 1 if gdename_std1 == "u.s.a, automatic musical instruments s.a., geneve (0.25)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "u.s.a, hewlett packard s.a., geneve (0.2)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "u.s.a, new britain machine company international, geneve" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "ueberlingen a. b, spaeth ag, arbon (0.05)" &  gdename_4a == "arbon"
replace wrong4a = 1 if gdename_std1 == "uitikon , a. a, rapid motormaeher ag, dietikon (4.0)" &  gdename_4a == "dietikon"
replace wrong4a = 1 if gdename_std1 == "usa, zweigniederlassung duebendorf, duebendorf" &  gdename_4a == "duebendorf"
replace wrong4a = 1 if gdename_std1 == "utimaco s.a., chene bourg (0.2)" &  gdename_4a == "chene bourg"
replace wrong4a = 1 if gdename_std1 == "utimaco s.a., geneve (0.2)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "vaduz (fl), succursale di chiasso, chiasso" &  gdename_4a == "chiasso"
replace wrong4a = 1 if gdename_std1 == "vaduz, handels und warenfinanzierungs ag, zuerich (0.05)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "vaduz, succursale di chiasso, chiasso" &  gdename_4a == "chiasso"
replace wrong4a = 1 if gdename_std1 == "vallechiara sa, capolago (0.1)" &  gdename_4a == "capolago"
replace wrong4a = 1 if gdename_std1 == "ventura/calif, oelfeld bohrgestaenge dienst ag, zug (0.99)" &  gdename_4a == "zug"
replace wrong4a = 1 if gdename_std1 == "verkehrsbetriebe des zuercher oberlandes, grueningen" &  gdename_4a == "grueningen"
replace wrong4a = 1 if gdename_std1 == "verlagsanstalt benzinger co. ag, einsiedeln (1.53)" &  gdename_4a == "einsiedeln"
replace wrong4a = 1 if gdename_std1 == "versailles f, le previsa sa, geneve (0.5)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "vesinet f, alexander stenhouse (international) sa, geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "vesinet f, alexander stenhouse (international) sa., geneve (0.05)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "veyrier centre d information et de public relations â« cipr â», geneve (0.1)" &  gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "vezia zuerichag, zuerich (0.075)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "vouvry , atelier de constructions mecaniques et metalliques de **" &  gdename_4a == "vouvry"
replace wrong4a = 1 if gdename_std1 == "waldstatt nufer co ag zwirnerei garnhandel, urnaesch" &  gdename_4a == "urnaesch"
replace wrong4a = 1 if gdename_std1 == "waschanstalt winterthur ag (wawag), winterthur (0.1)" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "wassenaar holl, shell (switzerland), zuerich (60.0)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "werkzeuge maschinen ag, buelach" &  gdename_4a == "buelach"
replace wrong4a = 1 if gdename_std1 == "william pemberton harry, hamburgad, zuerich (0.28)" &  gdename_4a == "zuerich"
replace wrong4a = 1 if gdename_std1 == "winterthur, calorifer ag, elgg (0.4)" &  gdename_4a == "elgg"
replace wrong4a = 1 if gdename_std1 == "winterthurâ» lebensversicherungs gesellschaft, winterthur" &  gdename_4a == "winterthur"
replace wrong4a = 1 if gdename_std1 == "worb, verzinkerei wattenwil ag, wattenwil (0.05)" &  gdename_4a == "wattenwil"
replace wrong4a = 1 if gdename_std1 == "zolli o, ditona ag, obergoesgen (0.5)" &  gdename_4a == "obergoesgen"
replace wrong4a = 1 if gdename_std1 == "zollikon, gluehlampenfabrik gloria ag, aarau (0.124)" &  gdename_4a == "aarau"
replace wrong4a = 1 if gdename_std1 == "zuerich, tegula ag, niederurnen (5.0)" &  gdename_4a == "niederurnen"
replace wrong4a = 1 if gdename_std1 == "aarbergg. 29 / schw. coiffeurv, fachsch., 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "ador, quai gustave 26, 1207 geneve" & gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "ador, quai gustave 56, 1207geneve" & gdename_4a == "geneve"
replace wrong4a = 1 if gdename_std1 == "amtshausg. 1 drahtseilbahn marzili stadt, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "amtshausg. 1/drahtseilbahn marzili stadt, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "amtshausg. 1/drahtseilbahn marzill stadt, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "amtshausg. 1/drahtseilbahn marzlll stadt, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "bankg. 8, c/o brunschwyler ag, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c o airtour suisse sa, bahnhofstr. 10, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c o herag, spitalgasse 32, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c o lonestar ag, spitalgasse 37, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c o mardlicka, st. georgenstr., 9oll st. gallen" & gdename_4a == "st. gallen"
replace wrong4a = 1 if gdename_std1 == "c. o immob. ges. de la gare, neueng. 25, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o bauspiegel verlag, neuengasse 39, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o buchh. faehre ag, maulbeerstr. 10, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o calzaferi sa, 7165 breil/brigels" & gdename_4a == "breil/brigels"
replace wrong4a = 1 if gdename_std1 == "c/o finaz ag, hirzbodenpark 10, 4052 basel" & gdename_4a == "basel"
replace wrong4a = 1 if gdename_std1 == "c/o herag, spitalgasse 32, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o immob. ges. de la gare, neueng. 25, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o kanalgasse 3 ag, bahnhofplatz 5, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o lenz, waisenhauspiatz 14, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o lenz, waisenhausplatz 14, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o lonestar ag, spitalgasse 37, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o ruprecht, christoffelgasse 7, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o schaerer ag, schauplatzgasse 26, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o vibo ag, rathausgasse 18/8 9, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "c/o w. p. moor ag, gutenbergstr. 14, 3oll bern" & gdename_4a == "bern"
replace wrong4a = 1 if gdename_std1 == "lausanne fleurier watch co s.a., fleurier (0.5)" & gdename_4a == "fleurier"
replace wrong4a = 1 if gdename_std1 == "lausanne, au richelieu s.a., yverdon (0.05)" & gdename_4a == "yverdon"
replace wrong4a = 1 if gdename_std1 == "lausanne. la motte s.a., fribourg (0.1)" & gdename_4a == "fribourg"
replace wrong4a = 1 if gdename_std1 == "soral, route de 19a, 1232confignon" & gdename_4a == "soral"

// corrections 3a
replace gdename = "affoltern am albis" if gdename_std1 == "affoltern a. a, autropa ag., zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "affoltern am albis" if gdename_std1 == "affoltern a. a, kyburz kieswerk obfelden ag, obfelden (0.06)" & gdename_4a == "obfelden"
replace gdename = "affoltern am albis" if gdename_std1 == "affoltern a.a, erich schwegler ag, ottenbach (0.15)" & gdename_4a == "ottenbach"
replace gdename = "affoltern am albis" if gdename_std1 == "affoltern a.a, frego ottenbach ag, ottenbach (0.1)" & gdename_4a == "ottenbach"
replace gdename = "affoltern am albis" if gdename_std1 == "affoltern a.a, kyburz kieswerk obfelden ag, obfelden (0.06)" & gdename_4a == "obfelden"
replace gdename = "affoltern am albis" if gdename_std1 == "affoltern a.a, leuthard soehne ag, merenschwand (0.3)" & gdename_4a == "merenschwand"
replace gdename = "affoltern am albis" if gdename_std1 == "affoltern a.a, tecalit ag, doettingen (0.05)" & gdename_4a == "doettingen"
replace gdename = "auvernier" if gdename_std1 == "auveruier, sportswear s.a., neuchatel (0.05)" & gdename_4a == "neuchatel"
replace gdename = "langnau im emmental" if gdename_std1 == "baeraustr. 43, .3550 langnau i.e, schweizerische kaeseunion ag, bern (11.563)" & gdename_4a == "bern"
replace gdename = "beinwil am see" if gdename_std1 == "beinwil a. s, rogelin watch ag, grenchen (0.05)" & gdename_4a == "grenchen"
replace gdename = "beinwil am see" if gdename_std1 == "beinwil a.s, cerapor ag, unterkulm (0.825)" & gdename_4a == "unterkulm"
replace gdename = "beinwil am see" if gdename_std1 == "beinwil a.s, gamma radiatoren ag, gontenschwil (0.75)" & gdename_4a == "gontenschwil"
replace gdename = "beinwil am see" if gdename_std1 == "beinwil a.s, sabruma ag, zuerich (0.12)" & gdename_4a == "zuerich"
replace gdename = "biel/bienne" if gdename_std1 == "biel. carrosseriewerke ag biel nidau, nidau (0.08)" & gdename_4a == "nidau"
replace gdename = "iseltwald" if gdename_std1 == "blatten, 3807iseltwaid" & gdename_4a == "blatten"
replace gdename = "bremgarten bei bern" if gdename_std1 == "brislachstr. 5, 422 , breitenbach" & gdename_4a == "breitenbach"
replace gdename = "bueren an der aare" if gdename_std1 == "bueren a. a, wollweberei rothrist ag, rothrist (0.15)" & gdename_4a == "rothrist"
replace gdename = "corsier sur vevey" if gdename_std1 == "corsier/vevey, wuton ag fuer feinmechanik und elektrotechnik, bern (0.2)" & gdename_4a == "bern"
replace gdename = "davos" if gdename_std1 == "davos platz ag luftseilbahn parsenn weissl9uhgipfel, langwies (0.4)" & gdename_4a == "langwies"
replace gdename = "haeutligen" if gdename_std1 == "dorf, 351c haeutligen" & gdename_4a == "dorf"
replace gdename = "adelboden" if gdename_std1 == "dorf, 3715adelboden" & gdename_4a == "dorf"
replace gdename = "zell (lu)" if gdename_std1 == "dorf, 61.44 zell lu" & gdename_4a == "dorf"
replace gdename = "zell (lu)" if gdename_std1 == "dorf, 614.4 zell lu" & gdename_4a == "dorf"
replace gdename = "zell (lu)" if gdename_std1 == "dorf, 614a zell lu" & gdename_4a == "dorf"
replace gdename = "hospental" if gdename_std1 == "egg, ss493 hospental" & gdename_4a == "egg"
replace gdename = "nyon" if gdename_std1 == "eysins, chemin d 42.1260 nyon" & gdename_4a == "eysins"
replace gdename = "geneve" if gdename_std1 == "geneve ., cominter s.a., fribourg (0.6)" & gdename_4a == "fribourg"
replace gdename = "geneve" if gdename_std1 == "geneve arts menagers s.a. (amsa), chene bougeries (0.5)" & gdename_4a == "chene bougeries"
replace gdename = "langnau am albis" if gdename_std1 == "gontenbach langnau a. a, proconsult ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "hausen am albis" if gdename_std1 == "hausen a. a, aschmann scheller ag buchdruckerei zur froschau, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "hausen am albis" if gdename_std1 == "hausen a. a, aschmann schuller ag buchdruckerei zur froschau, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "hausen am albis" if gdename_std1 == "hausen a. a, immobilien ag zur froschau, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "hausen am albis" if gdename_std1 == "hausen a. a, loring ag, mettmenstetten (0.1)" & gdename_4a == "mettmenstetten"
replace gdename = "kuessnacht (sz)" if gdename_std1 == "kiissnacht a.r, prohydro ag, zollikon (0.2)" & gdename_4a == "zollikon"
replace gdename = "kilchberg (zh)" if gdename_std1 == "kilchber9, philips ag, zuerich (20.0)" & gdename_4a == "zuerich"
replace gdename = "kuesnacht (zh)" if gdename_std1 == "kuesnacht, , baumwoll zwirnerei mols ag, quarten (0.2)" & gdename_4a == "quarten"
replace gdename = "kuesnacht (zh)" if gdename_std1 == "kuesnacht, j. ouboter ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "kuessnacht am rigi" if gdename_std1 == "kuessnacht a r, prohydro ag, zollikon (0.2)" & gdename_4a == "zollikon"
replace gdename = "kuessnacht am rigi" if gdename_std1 == "kuessnacht a. r, orso immobilien ag, aarau (0.01)" & gdename_4a == "aarau"
replace gdename = "kuessnacht am rigi" if gdename_std1 == "kuessnacht a. r, trafina ag, basel (0.1)" & gdename_4a == "basel"
replace gdename = "kuessnacht am rigi" if gdename_std1 == "kuessnacht a. r, villa sagirain, kuessnacht am rigi, schwyz (0.1)" & gdename_4a == "schwyz"
replace gdename = "kuessnacht am rigi" if gdename_std1 == "kuessnacht a.r, emu ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "kuessnacht am rigi" if gdename_std1 == "kuessnacht a.r, emu ag, zuerich (1.0)" & gdename_4a == "zuerich"
replace gdename = "diesse" if gdename_std1 == "lamboing, rte de, .2517 diesse" & gdename_4a == "lamboing"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a. a, atlantis parkhotel ag, zuerich (1.0)" & gdename_4a == "zuerich"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a. a, barverkaufs ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a. a, fivat ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a. a, melitta ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a. a, notag, e. nebel orientteppich transithandel ag, zuerich (1.0)" & gdename_4a == "zuerich"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a. a, texsana reinigung hirschmatt luzern ag, luzern" & gdename_4a == "luzern"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a. a, trevin ag, treuhand inkasso, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a.a, profib ag, hinwil (0.25)" & gdename_4a == "hinwil"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a.a, profib ag, zuerich (0.25)" & gdename_4a == "zuerich"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a.a, reflo ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "langnau am albis" if gdename_std1 == "langnau a.a, skilift ag tschappina, tschappina (0.25)" & gdename_4a == "tschappina"
replace gdename = "langnau im emmental" if gdename_std1 == "langnau i. e, pferdemetzgerei horisberger ag, aarau (0.05)" & gdename_4a == "aarau"
replace gdename = "langnau im emmental" if gdename_std1 == "langnau i. e, tuileries de corbieres s.a., corbieres (0.24)" & gdename_4a == "corbieres"
replace gdename = "langnau im emmental" if gdename_std1 == "langnau i.e, holdingges. pilatus ag, burgdorf (0.306)" & gdename_4a == "burgdorf"
replace gdename = "lausanne" if gdename_std1 == "lausanne mela mecanisation laminoire s.a., bex (0.05)" & gdename_4a == "bex"
replace gdename = "lausanne" if gdename_std1 == "lausanne., s.i. monts et soleil s.a. verbier, bagnes (0.05)" & gdename_4a == "bagnes"
replace gdename = "chexbres" if gdename_std1 == "lignieres, t605 chexbres" & gdename_4a == "lignieres"
replace gdename = "muri bei bern" if gdename_std1 == "mari b. b, kaiser co ag, bern (1.2)" & gdename_4a == "bern"
replace gdename = "meggen" if gdename_std1 == "meggen bauunternehmung fritschi ag, luzern (0.05)" & gdename_4a == "luzern"
replace gdename = "muri bei bern" if gdename_std1 == "mud b. b, eika ag, bern (1.0)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muni b. b, j. hirter co ag, bern (0.8)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muni b. b, schoenholzer ag, bern (0.06)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b b, berner handelsbank ag, bern (0.1)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, bauges. nubag ag, bern (0.05)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, bauges. nuhag ag, bern (0.05)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, eika ag, bern (1.0)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, f. anker cie ag, bern (0.32)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, f. anker cie. ag, bern (0.32)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, financiere de participations internationales s.a., bern (1.35)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, garage egghoelzli ag, bern (0.1)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, imlo immobilien ag, bern (0.4)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, j. hinter co ag, bern (0.8)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, j. hirter co ag, bern (0.8)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, j. hitter co ag, bern (0.8)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, jucker cie ag, bern (0.07)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, kaiser co ag, bern (1.2)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, privatklinik engeried, bern (0.16)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, privatklinik engeried, bern (1.06)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, rud. kull ag, bern (0.15)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, schoenholer ag, bern (0.06)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. b, schoenholzer ag, bern (0.06)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. e, privatklinik engeried, bern (1.06)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b. r, privatklinik engeried, bern (1.06)" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri b., bern" & gdename_4a == "bern"
replace gdename = "muri bei bern" if gdename_std1 == "muri, bern" & gdename_4a == "bern"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a r, kieswerk siggenthal ag, untersiggenthal (1.0)" & gdename_4a == "untersiggenthal"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a r, rovo claude ag, zuerich (0.275)" & gdename_4a == "zuerich"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. r, aiag engineering ag, chippis (0.5)" & gdename_4a == "chippis"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. r, gebrueder bruehlmann ag, schaffhausen (0.8)" & gdename_4a == "schaffhausen"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. r, ges. der schaffhauser kaffeehallen, schaffhausen (0.2)" & gdename_4a == "schaffhausen"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. r, kieswerk siggenthal ag, untersiggenthal (1.0)" & gdename_4a == "untersiggenthal"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. r, mag engineering ag, chippis (0.5)" & gdename_4a == "chippis"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. r, preluwag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. r, rhenum metall ag, schaffhausen (0.05)" & gdename_4a == "schaffhausen"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. r, schaffhauser strickmaschinenfabrik, schaffhausen (1.0)" & gdename_4a == "schaffhausen"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. r, wohnbauges. der georg fischer ag, schaffhausen" & gdename_4a == "schaffhausen"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. r, wohnbauges. der georg fisher ag, schaffhausen (0.5)" & gdename_4a == "schaffhausen"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a. rh, buerox ag, bueren an der aare (1.4)" & gdename_4a == "bueren an der aare"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a.r, bruetsch leu ag. hoch und tiefbau, schaffhausen" & gdename_4a == "schaffhausen"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a.r, comestibles g. barras s.a., montana (0.08)" & gdename_4a == "montana"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a.r, e. stuber ag, bern (0.1)" & gdename_4a == "bern"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a.r, lega ag, bern (0.09)" & gdename_4a == "bern"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a.r, thenum metall ag, schaffhausen (0.05)" & gdename_4a == "schaffhausen"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a.r, villa la foresta ag, basel (0.05)" & gdename_4a == "basel"
replace gdename = "neuhausen am rheinfall" if gdename_std1 == "neuhausen a.r, villa la forests ag, basel (0.05)" & gdename_4a == "basel"
replace gdename = "oetwil am see" if gdename_std1 == "oetwil a.s, skilift barga parsenn ag , langwies (1.0)" & gdename_4a == "langwies"
replace gdename = "oetwil am see" if gdename_std1 == "oetwil a.s, skilift parsenn fondel ag, schiers (0.1)" & gdename_4a == "schiers"
replace gdename = "oetwil am see" if gdename_std1 == "oetwil a.s, suvretta bau ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "ostermundigen" if gdename_std1 == "ostermundigen, sonorfilm ag, bern (0.21)" & gdename_4a == "bern"
replace gdename = "ilanz" if gdename_std1 == "schnaus, 71301lanz" & gdename_4a == "schnaus"
replace gdename = "uetikon am see" if gdename_std1 == "steinmaur, 8707.uetikon am see" & gdename_4a == "steinmaur"
replace gdename = "torricella taverne" if gdename_std1 == "taverne, dipharm ag, zug (0.05)" & gdename_4a == "zug"
replace gdename = "teufen (ar)" if gdename_std1 == "teuf en city grund ag, st gallen (0.05)" & gdename_4a == "st gallen"
replace gdename = "teufen (ar)" if gdename_std1 == "teuf en durofer ag, st. gallen (0.25)" & gdename_4a == "st. gallen"
replace gdename = "teufen (ar)" if gdename_std1 == "teuf en pharmag ag, herisau (0.05)" & gdename_4a == "herisau"
replace gdename = "teufen (ar)" if gdename_std1 == "teuf en werner hachen ag, st. gallen (0.07)" & gdename_4a == "st. gallen"
replace gdename = "thonex" if gdename_std1 == "thonex4, brissonneau et lotz international s.a., geneve (0.5)" & gdename_4a == "geneve"
replace gdename = "uetikon am see" if gdename_std1 == "uctikon a.a, s.i. plainpalais centre a, geneve (0.05)" & gdename_4a == "geneve"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a s, chemap ag, maennedorf (0.05)" & gdename_4a == "maennedorf"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a. s, hydroma ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a. s, memmel co. ag, basel (0.5)" & gdename_4a == "basel"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a. s, r.holliger co ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a. s, r.holliger co ag, zuerich (0.5)" & gdename_4a == "zuerich"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a.a, s.i. plainpalais centre a, geneve (0.05)" & gdename_4a == "geneve"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a.s, bomina ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a.s, gewerbebank maennedorf, maennedorf (1.2)" & gdename_4a == "maennedorf"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a.s, inrescor internationale forschungs ges. (ag), schwerzenbach (3.0)" & gdename_4a == "schwerzenbach"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a.s, johann weibel ag, mechanische ziegelei, eschlikon" & gdename_4a == "eschlikon"
replace gdename = "uetikon am see" if gdename_std1 == "uetikon a.s, memmel co. ag, basel (0.5)" & gdename_4a == "basel"
replace gdename = "uitikon" if gdename_std1 == "ui ikon a.a, scale , ia ag, chur (0.05)" & gdename_4a == "chur"
replace gdename = "uitikon" if gdename_std1 == "uitikon a a, aircraft parts ag., zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, a. blum co. ag, zuerich (1.0)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, aircraft parts ag., zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, alex landau ag, zuerich (0.3)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, aventa ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, boutique robat ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, carpentier ag, zuerich (0.5)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, gebr. kuoni ag, zuerich (0.5)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, gebr. kuoni. ag, zuerich (0.5)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, haus delphin ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, ibag ag, uitikon (0.05)" & gdename_4a == "uitikon"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, neue buecher ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, progress lederwaren und sportartikelfabrik ag, zuerich (0.08)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, rapid baumaschinen ag, dietikon (0.3)" & gdename_4a == "dietikon"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, rapid motormaeher ag, zuerich (3.0)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, robat ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, schneider co ag, dietikon (0.1)" & gdename_4a == "dietikon"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, segna ag, zuerich (0.25)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, siemens asia investments ag, zuerich (15.0)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, siemens europa beteiligungen ag, zuerich (50.0)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, stala immobilien ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a. a, willy grob ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, agia ag fuer industrie automation, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, aschmann scheller ag buchdruckerei zur froschau, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, baechtold ag. immobilien und treuhand gesellschaft, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, daetwyler optik ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, ernst mittelholzer ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, hamaco ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, kurt huber ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, scaletta ag, chur (0.3)" & gdename_4a == "chur"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, semperit ag, zuerich (0.25)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, systembau ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, systembau ag, zuerich (0.25)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, systembau ag, zuerich (0.65)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, verwaltungsges. fuer investment trusts (vit), zuerich" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, vibro technik ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.a, woh bauelement ag, zuerich (0.1)" & gdename_4a == "zuerich"
replace gdename = "uitikon" if gdename_std1 == "uitikon a.s, si tro sicherheitstresor ag, zuerich (0.35)" & gdename_4a == "zuerich"
replace gdename = "blonay" if gdename_std1 == "vevey, route de 65b. 1807 blonay" & gdename_4a == "vevey"
replace gdename = "sumiswald" if gdename_std1 == "wasen i. e, daniel sinzig ag, thun (0.16)" & gdename_4a == "thun"
replace gdename = "sumiswald" if gdename_std1 == "wasen i. e, metallveredelun8s und spritzwerk ag, thun(0.2)" & gdename_4a == "thun"
replace gdename = "sumiswald" if gdename_std1 == "wasen i. e, metallveredelungs und spritzwerk ag, thun (0.2)" & gdename_4a == "thun"
replace gdename = "sumiswald" if gdename_std1 == "wasen i. e, metallveredelungs und spritzwerk ag, thun(0.2)" & gdename_4a == "thun"
replace gdename = "wettswil" if gdename_std1 == "wettswil a. a, airtechnik ag, zuerich (0.05)" & gdename_4a == "zuerich"
replace gdename = "wettswil" if gdename_std1 == "wettswil a.a, laubscher spiegel ag, zuerich (0.25)" & gdename_4a == "zuerich"
replace gdename = "waedenswil" if gdename_std1 == "widen, 8820waedenswll" & gdename_4a == "widen"
replace gdename = "waedenswil" if gdename_std1 == "widen, 882owaedenswll" & gdename_4a == "widen"
replace gdename = "innertkirchen" if gdename_std1 == "winkel, 38621nnertkirchen" & gdename_4a == "winkel"
replace gdename = "ittigen" if gdename_std1 == "worblauf en ziegelei tiefenau ag, bern (0.4)" & gdename_4a == "bern"
replace gdename = "zollikofen" if gdename_std1 == "zollikof en heag ag, bern (0.236)" & gdename_4a == "bern"
replace gdename = "zuerich" if gdename_std1 == "zuerich ., cartonnagefabrik waedenswil ag, waedenswil (0.15)" & gdename_4a == "waedenswil"
replace gdename = "zuerich" if gdename_std1 == "zuerich. balimex ag, glarus (0.1)" & gdename_4a == "glarus"
replace gdename = "zuerich" if gdename_std1 == "zuerich. drahtseilbahn interlaken heimwehfluh ag, interlaken" & gdename_4a == "interlaken"
replace gdename = "zug" if gdename_std1 == "zup, abnox ag, cham (0.25)" & gdename_4a == "cham"

// drop irrelevant duplicates (due to reshape and merge)
replace gdename = "" if wrong4a == 1
bysort PID year ID_dupl gdename_std1 gdename: gen nr_double = _n  // identify duplicates after correction of wrong assignments in 4a
drop if nr_double > 1
duplicates tag PID year ID_dupl gdename_std1, gen(help_dupl4a)  // identify  cases with multiple original entries. If only correct assignment: keep; if only wrong assignment: keep empty gdename
drop if help_dupl4a > 0 & wrong4a == 1  // drop multiple empty gdenames per ID-year-ID_dupl observation
drop help_dupl4a nr_double

bysort PID year ID_dupl gdename_std1: gen options4a = _n   // some observations contain two potentially correct gdenames.

sort gdename
merge m:1 gdename using "$path\02_Processed_data\08_Municipalities\unique_Gdename_Gdenr2018-Geo.dta"
drop if _merge==2
rename _merge _merge4a2
gen gdename_used = gdename
gen GdeNr_E_CNTR = E_CNTR
gen GdeNr_N_CNTR = N_CNTR
replace geo = -9 if wrong4a == 1

keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
GdeNr_E_CNTR GdeNr_N_CNTR gdename_used gdenr_2018 geo options4a
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo41.dta", replace

* 4b corrections
use "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo4abc.dta", clear
keep if geo == 42  // only use observations from round 3b

gen gdename_used = gdename_4b
gen gdenr_2018 = gdenr_2018_4b
gen GdeNr_E_CNTR = E_CNTR_4b
gen GdeNr_N_CNTR = N_CNTR_4b

keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
GdeNr_E_CNTR GdeNr_N_CNTR gdename_used gdenr_2018 geo
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo42.dta", replace


* 4c corrections
use "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo4abc.dta", clear
keep if geo == 43  // only use observations from round 4c

drop E_CNTR_2a N_CNTR_2a _merge_2a E_CNTR_2b N_CNTR_2b _merge_2b ///
_merge_3 ///
gdename_4a gdename_4b E_CNTR_4a N_CNTR_4a _merge_4a gdename_4a ///
gdenr_2018_4b E_CNTR_4b N_CNTR_4b _merge_4b gdename_4b gdename 

reshape long ctn_ _merge_4c gdenr_2018_4c E_CNTR_4c N_CNTR_4c, i(PID year ID_dupl n_obs_4) j(n_obs_4c)
drop if _merge_4c == 1

// produce list with unique gde-info for manual checks
preserve
gen n = -1
collapse (sum) n, by (gdename_std1 gdenr_2018_4c gdename_4c ctn)
order n 
sort n
save "$path\02_Processed_data\10_Directors_1934_2003\check_4c.dta", replace
restore

// drop wrong matches (by taking gdenr_2018 into account, there is no risk of dropping original observation entirely)
drop if gdename_std1 == "accum ag, gossau (0.4)" & gdenr_2018_4c == 3443
drop if gdename_std1 == "accum ag, gossau (0.4)" & gdenr_2018_4c == 115
drop if gdename_std1 == "aesch (basel)" & gdenr_2018_4c == 1021
drop if gdename_std1 == "aesch (bld)" & gdenr_2018_4c == 1021
drop if gdename_std1 == "apotheloz s.a., colombier (0.15)" & gdenr_2018_4c == 6416
drop if gdename_std1 == "apotheloz s.a., colombier (0.15)" & gdenr_2018_4c == 5634
drop if gdename_std1 == "aubrigstr. 15, rueti" & gdenr_2018_4c == 1631
drop if gdename_std1 == "augartenstr. 36, andwil" & gdenr_2018_4c == 4901
drop if gdename_std1 == "av. tir fedederal 94, ecublens" & gdenr_2018_4c == 2072
drop if gdename_std1 == "bad hamburg v.d.h. decsmholdingund finanz ag, buchs, buchs" & gdenr_2018_4c == 3271
drop if gdename_std1 == "bad hamburg v.d.h. decsmholdingund finanz ag, buchs, buchs" & gdenr_2018_4c == 83
drop if gdename_std1 == "bad hamburg v.d.h. decsmholdingund finanz ag, buchs, buchs" & gdenr_2018_4c == 4003
drop if gdename_std1 == "bad hamburg v.d.h. decsmholdingund finanz ag, buchs, buchs" & gdenr_2018_4c == 1125
drop if gdename_std1 == "baerenfelsweg 14, aesch" & gdenr_2018_4c == 1021
drop if gdename_std1 == "bahnhofstr. 21, biel" & gdenr_2018_4c == 6077
drop if gdename_std1 == "berg (d)" & gdenr_2018_4c == 3211
drop if gdename_std1 == "berg (d)" & gdenr_2018_4c == 4891
drop if gdename_std1 == "berg , starnberger see d" & gdenr_2018_4c == 4891
drop if gdename_std1 == "berg , starnberger see d" & gdenr_2018_4c == 3211
drop if gdename_std1 == "berg, 6390.engelberg" & gdenr_2018_4c == 4891
drop if gdename_std1 == "berg, 6390.engelberg" & gdenr_2018_4c == 3211
drop if gdename_std1 == "berg, starnberger see d" & gdenr_2018_4c == 3211
drop if gdename_std1 == "berg, starnberger see d" & gdenr_2018_4c == 4891
drop if gdename_std1 == "biel ., schnyder liechti ag, biel (0.2)" & gdenr_2018_4c == 6077
drop if gdename_std1 == "biel ., schnyder liechti ag, biel (0.2)" & gdenr_2018_4c == 2764
drop if gdename_std1 == "biel, imaco ag uhrenfabrik, biel (0.4)" & gdenr_2018_4c == 6077
drop if gdename_std1 == "biel, imaco ag uhrenfabrik, biel (0.4)" & gdenr_2018_4c == 2764
drop if gdename_std1 == "blockmetall ag, buchs (0.2)" & gdenr_2018_4c == 83
drop if gdename_std1 == "blockmetall ag, buchs (0.2)" & gdenr_2018_4c == 3271
drop if gdename_std1 == "blockmetall ag, buchs (0.2)" & gdenr_2018_4c == 1125
drop if gdename_std1 == "blockmetall ag, buchs (0.2)" & gdenr_2018_4c == 4003
drop if gdename_std1 == "blockmetall ag, buchs (0.4)" & gdenr_2018_4c == 1125
drop if gdename_std1 == "blockmetall ag, buchs (0.4)" & gdenr_2018_4c == 83
drop if gdename_std1 == "blockmetall ag, buchs (0.4)" & gdenr_2018_4c == 3271
drop if gdename_std1 == "blockmetall ag, buchs (0.4)" & gdenr_2018_4c == 4003
drop if gdename_std1 == "campo (bienio)" & gdenr_2018_4c == 5307
drop if gdename_std1 == "campo (blenlo)" & gdenr_2018_4c == 5307
drop if gdename_std1 == "ch. champ colomb 8, ecublens" & gdenr_2018_4c == 2072
drop if gdename_std1 == "ch. de l ormet 14, ecublens" & gdenr_2018_4c == 2072
drop if gdename_std1 == "ch. du stand 15b, ecublens" & gdenr_2018_4c == 2072
drop if gdename_std1 == "ch. tortiquet 15, villeneuve" & gdenr_2018_4c == 2044
drop if gdename_std1 == "cherventino, gerra" & gdenr_2018_4c == 5138
drop if gdename_std1 == "der ey, kirchberg" & gdenr_2018_4c == 3392
drop if gdename_std1 == "dorfstr. 36, kilchberg" & gdenr_2018_4c == 2851
drop if gdename_std1 == "engleiten, bad ischl, schauberger biotechnik ag, wetzikon (0.15)" & gdenr_2018_4c == 4611
drop if gdename_std1 == "engleiten, bad ischl, schauberger biotechnik ag, wetzikon (0.15)" & gdenr_2018_4c == 121
drop if gdename_std1 == "ettingerstr. 47, reinach" & gdenr_2018_4c == 4141
drop if gdename_std1 == "eujiyoshida (j), matile sa, biel (0.3)" & gdenr_2018_4c == 371
drop if gdename_std1 == "eujiyoshida (j), matile sa, biel (0.3)" & gdenr_2018_4c == 6077
drop if gdename_std1 == "fleischbachstr. 9, reinach" & gdenr_2018_4c == 4141
drop if gdename_std1 == "forel (lavauz)" & gdenr_2018_4c == 2054
drop if gdename_std1 == "forel (uvaux)" & gdenr_2018_4c == 2054
drop if gdename_std1 == "forel (vd)" & gdenr_2018_4c == 2054
drop if gdename_std1 == "fribourg sipra s.a., rossens (0.05)" & gdenr_2018_4c == 2236
drop if gdename_std1 == "fribourg sipra s.a., rossens (0.05)" & gdenr_2018_4c == 5830
drop if gdename_std1 == "friesenheid, schmitten" & gdenr_2018_4c == 3514
drop if gdename_std1 == "gartenstr. 13, 8636, wald" & gdenr_2018_4c == 3036
drop if gdename_std1 == "glion s.1. le bluet s.a. villeneuve vd., villeneuve (0.05)" & gdenr_2018_4c == 5414
drop if gdename_std1 == "glion s.1. le bluet s.a. villeneuve vd., villeneuve (0.05)" & gdenr_2018_4c == 2044
drop if gdename_std1 == "heikendorf (d) ., cdc st. sulpice car dry clean sa, saint sulpice (0.05)" & gdenr_2018_4c == 5648
drop if gdename_std1 == "heikendorf (d) ., cdc st. sulpice car dry clean sa, saint sulpice (0.05)" & gdenr_2018_4c == 6512
drop if gdename_std1 == "hoeheweg 24, biel" & gdenr_2018_4c == 6077
drop if gdename_std1 == "ifangstr. 17, weiningen" & gdenr_2018_4c == 4621
drop if gdename_std1 == "imobu ag, biel (0.1)" & gdenr_2018_4c == 2764
drop if gdename_std1 == "imobu ag, biel (0.1)" & gdenr_2018_4c == 6077
drop if gdename_std1 == "imobu ag, biel (0.1)" & gdenr_2018_4c == 371
drop if gdename_std1 == "jegensdorferstr. 3b, zuzwil" & gdenr_2018_4c == 3426
drop if gdename_std1 == "kanalweg 4, kirchberg" & gdenr_2018_4c == 3392
drop if gdename_std1 == "kappel (sg)" & gdenr_2018_4c == 2580
drop if gdename_std1 == "le pommeret, le paquier" & gdenr_2018_4c == 6487
drop if gdename_std1 == "les cerisiers, 1754, onnens" & gdenr_2018_4c == 5565
drop if gdename_std1 == "marbach (d)" & gdenr_2018_4c == 1010
drop if gdename_std1 == "marbach (d)" & gdenr_2018_4c == 3253
drop if gdename_std1 == "neue jonastr. 81, rapperswil" & gdenr_2018_4c == 310
drop if gdename_std1 == "oberdorf (basel)" & gdenr_2018_4c == 1508
drop if gdename_std1 == "oberdorf (basel)" & gdenr_2018_4c == 2553
drop if gdename_std1 == "oberdorf (bld)" & gdenr_2018_4c == 2553
drop if gdename_std1 == "oberdorf (bld)" & gdenr_2018_4c == 1508
drop if gdename_std1 == "oberfeldring 18, arni" & gdenr_2018_4c == 602
drop if gdename_std1 == "oberwil (basel)" & gdenr_2018_4c == 4074
drop if gdename_std1 == "oberwil (basel)" & gdenr_2018_4c == 4571
drop if gdename_std1 == "oberwil (bid)" & gdenr_2018_4c == 4074
drop if gdename_std1 == "oberwil (bid)" & gdenr_2018_4c == 4571
drop if gdename_std1 == "r. de bassenges 2, ecublens" & gdenr_2018_4c == 2072
drop if gdename_std1 == "r. fortifications 8, villeneuve" & gdenr_2018_4c == 2044
drop if gdename_std1 == "rapperswil (zh)" & gdenr_2018_4c == 310
drop if gdename_std1 == "reinach (basel)" & gdenr_2018_4c == 4141
drop if gdename_std1 == "reinach (bid)" & gdenr_2018_4c == 4141
drop if gdename_std1 == "rickenbach (d)" & gdenr_2018_4c == 1097
drop if gdename_std1 == "rickenbach (d)" & gdenr_2018_4c == 2582
drop if gdename_std1 == "rickenbach (d)" & gdenr_2018_4c == 225
drop if gdename_std1 == "rickenbach (d)" & gdenr_2018_4c == 2857
drop if gdename_std1 == "riedweg 10, biel" & gdenr_2018_4c == 6077
drop if gdename_std1 == "rigi schiffahrtsgesellschaft ag, kilchberg (0.6)" & gdenr_2018_4c == 135
drop if gdename_std1 == "rigi schiffahrtsgesellschaft ag, kilchberg (0.6)" & gdenr_2018_4c == 2851
drop if gdename_std1 == "rolliweg 15, lengnau" & gdenr_2018_4c == 4312
drop if gdename_std1 == "rossweg, schmitten" & gdenr_2018_4c == 2305
drop if gdename_std1 == "rte arrufens 8, romont" & gdenr_2018_4c == 442
drop if gdename_std1 == "rte arruffens 57, romont" & gdenr_2018_4c == 442
drop if gdename_std1 == "schmitten (albula)" & gdenr_2018_4c == 2305
drop if gdename_std1 == "schmitten (d)" & gdenr_2018_4c == 3514
drop if gdename_std1 == "schmitten (d)" & gdenr_2018_4c == 2305
drop if gdename_std1 == "schnegg tools ag, biel (0.1)" & gdenr_2018_4c == 371
drop if gdename_std1 == "schnegg tools ag, biel (0.1)" & gdenr_2018_4c == 6077
drop if gdename_std1 == "schuetzengasse, altdorf" & gdenr_2018_4c == 2920
drop if gdename_std1 == "schulweg 21, kirchberg" & gdenr_2018_4c == 3392
drop if gdename_std1 == "seestr. 129, kilchberg" & gdenr_2018_4c == 2851
drop if gdename_std1 == "silvretta liegenschaften ag, kilchberg (0.1)" & gdenr_2018_4c == 2851
drop if gdename_std1 == "silvretta liegenschaften ag, kilchberg (0.1)" & gdenr_2018_4c == 135
drop if gdename_std1 == "sonnenberg 224, holderbank" & gdenr_2018_4c == 4199
drop if gdename_std1 == "storchenweg 6, oberwil" & gdenr_2018_4c == 4571
drop if gdename_std1 == "studen, biel" & gdenr_2018_4c == 2764
drop if gdename_std1 == "studen, biel" & gdenr_2018_4c == 6077
drop if gdename_std1 == "wangen, druessel ag kohlen heizoele, wangen (0.075)" & gdenr_2018_4c == 200
drop if gdename_std1 == "wangen, druessel ag kohlen heizoele, wangen (0.075)" & gdenr_2018_4c == 1349
drop if gdename_std1 == "weiningen (baden)" & gdenr_2018_4c == 4621
drop if gdename_std1 == "weissensteinstr. 146, oberdorf" & gdenr_2018_4c == 2892
drop if gdename_std1 == "weissensteinstr. 146, oberdorf" & gdenr_2018_4c == 1508
drop if gdename_std1 == "wetzikon, hypoterra ag, wetzikon (0.5)" & gdenr_2018_4c == 4611
drop if gdename_std1 == "zell (a)" & gdenr_2018_4c == 1150
drop if gdename_std1 == "zell (a)" & gdenr_2018_4c == 231
drop if gdename_std1 == "zell (d)" & gdenr_2018_4c == 231
drop if gdename_std1 == "zell (d)" & gdenr_2018_4c == 1150
drop if gdename_std1 == "zuzwil, protametal ag, zuzwil (0.25)" & gdenr_2018_4c == 3426

// corrections (correction of gdenr_2018)
gen wrong4c =.
replace wrong4c = 2 if gdename_std1 == "aesch (zh)" & gdenr_2018_4c == 2761
replace wrong4c = 2 if gdename_std1 == "aesch (zh)" & gdenr_2018_4c == 1021
replace wrong4c = 2 if gdename_std1 == "chapelle, route 40, gran lancy" & gdenr_2018_4c == 2009
replace wrong4c = 2 if gdename_std1 == "chapelle, route 40, gran lancy" & gdenr_2018_4c == 2066
replace wrong4c = 2 if gdename_std1 == "granges (vd)" & gdenr_2018_4c == 6248
replace wrong4c = 2 if gdename_std1 == "granges (vd)" & gdenr_2018_4c == 2328
replace wrong4c = 2 if gdename_std1 == "jun., stein, am rhein" & gdenr_2018_4c == 4260
replace wrong4c = 2 if gdename_std1 == "jun., stein, am rhein" & gdenr_2018_4c == 3360
replace wrong4c = 2 if gdename_std1 == "jun., stein, am rhein" & gdenr_2018_4c == 3005
replace wrong4c = 2 if gdename_std1 == "langnau i. e, sesselbahn und skilift marbachegg ag, marbach (0.25)" & gdenr_2018_4c == 3253
replace wrong4c = 2 if gdename_std1 == "langnau i. e, sesselbahn und skilift marbachegg ag, marbach (0.25)" & gdenr_2018_4c == 1010
replace wrong4c = 2 if gdename_std1 == "montet (vully)" & gdenr_2018_4c == 2050
replace wrong4c = 2 if gdename_std1 == "montet (vully)" & gdenr_2018_4c == 2089
replace wrong4c = 2 if gdename_std1 == "rickenbach (sg)" & gdenr_2018_4c == 2857
replace wrong4c = 2 if gdename_std1 == "rickenbach (sg)" & gdenr_2018_4c == 1097
replace wrong4c = 2 if gdename_std1 == "rickenbach (sg)" & gdenr_2018_4c == 2582
replace wrong4c = 2 if gdename_std1 == "rickenbach (sg)" & gdenr_2018_4c == 225
replace wrong4c = 2 if gdename_std1 == "rickenbach (tg)" & gdenr_2018_4c == 1097
replace wrong4c = 2 if gdename_std1 == "rickenbach (tg)" & gdenr_2018_4c == 225
replace wrong4c = 2 if gdename_std1 == "rickenbach (tg)" & gdenr_2018_4c == 2857
replace wrong4c = 2 if gdename_std1 == "rickenbach (tg)" & gdenr_2018_4c == 2582
replace wrong4c = 2 if gdename_std1 == "rickenbachstrasse 155, rickenbach" & gdenr_2018_4c == 2857
replace wrong4c = 2 if gdename_std1 == "rickenbachstrasse 155, rickenbach" & gdenr_2018_4c == 1097
replace wrong4c = 2 if gdename_std1 == "rickenbachstrasse 155, rickenbach" & gdenr_2018_4c == 2582
replace wrong4c = 2 if gdename_std1 == "rickenbachstrasse 155, rickenbach" & gdenr_2018_4c == 225
replace wrong4c = 2 if gdename_std1 == "schwyz , ag projekta ingenieurbuero, altdorf (0.08)" & gdenr_2018_4c == 2920
replace wrong4c = 2 if gdename_std1 == "schwyz , ag projekta ingenieurbuero, altdorf (0.08)" & gdenr_2018_4c == 1201
replace wrong4c = 2 if gdename_std1 == "sytirain, rickenbach" & gdenr_2018_4c == 2857
replace wrong4c = 2 if gdename_std1 == "sytirain, rickenbach" & gdenr_2018_4c == 1097
replace wrong4c = 2 if gdename_std1 == "sytirain, rickenbach" & gdenr_2018_4c == 225
replace wrong4c = 2 if gdename_std1 == "sytirain, rickenbach" & gdenr_2018_4c == 2582
replace wrong4c = 2 if gdename_std1 == "sytiraln, rickenbach" & gdenr_2018_4c == 225
replace wrong4c = 2 if gdename_std1 == "sytiraln, rickenbach" & gdenr_2018_4c == 1097
replace wrong4c = 2 if gdename_std1 == "sytiraln, rickenbach" & gdenr_2018_4c == 2582
replace wrong4c = 2 if gdename_std1 == "sytiraln, rickenbach" & gdenr_2018_4c == 2857
replace wrong4c = 2 if gdename_std1 == "wangen a. a, ersparniskasse des amtbezirks ag, wangen (1.0)" & gdenr_2018_4c == 1349
replace wrong4c = 2 if gdename_std1 == "wangen a. a, ersparniskasse des amtbezirks ag, wangen (1.0)" & gdenr_2018_4c == 200
replace wrong4c = 2 if gdename_std1 == "wangen a. a, esparniskasse des amtsbezirks ag, wangen (1.0)" & gdenr_2018_4c == 200
replace wrong4c = 2 if gdename_std1 == "wangen a. a, esparniskasse des amtsbezirks ag, wangen (1.0)" & gdenr_2018_4c == 1349
replace wrong4c = 2 if gdename_std1 == "wangen a. a, oberaareauische automobilkurse ag, wangen" & gdenr_2018_4c == 200
replace wrong4c = 2 if gdename_std1 == "wangen a. a, oberaareauische automobilkurse ag, wangen" & gdenr_2018_4c == 1349
replace wrong4c = 2 if gdename_std1 == "wangen a. a, oberaargauische automobikurse ag, wangen (0.105)" & gdenr_2018_4c == 200
replace wrong4c = 2 if gdename_std1 == "wangen a. a, oberaargauische automobikurse ag, wangen (0.105)" & gdenr_2018_4c == 1349
replace wrong4c = 2 if gdename_std1 == "wangen a. a, oberaargauische automobilkurse ag, wangen" & gdenr_2018_4c == 1349
replace wrong4c = 2 if gdename_std1 == "wangen a. a, oberaargauische automobilkurse ag, wangen" & gdenr_2018_4c == 200
replace wrong4c = 2 if gdename_std1 == "wangen a. a, obrecht co ag, wangen (0.8)" & gdenr_2018_4c == 200
replace wrong4c = 2 if gdename_std1 == "wangen a. a, obrecht co ag, wangen (0.8)" & gdenr_2018_4c == 1349
replace wrong4c = 2 if gdename_std1 == "wangen a. a, r. schweizer co ag, wangen (0.3)" & gdenr_2018_4c == 200
replace wrong4c = 2 if gdename_std1 == "wangen a. a, r. schweizer co ag, wangen (0.3)" & gdenr_2018_4c == 1349
replace wrong4c = 2 if gdename_std1 == "wangen a.a, r. schweizer co ag, wangen (0.5)" & gdenr_2018_4c == 200
replace wrong4c = 2 if gdename_std1 == "wangen a.a, r. schweizer co ag, wangen (0.5)" & gdenr_2018_4c == 1349
replace wrong4c = 2 if gdename_std1 == "wannen a. a, r. schweizer co ag, wangen (0.3)" & gdenr_2018_4c == 1349
replace wrong4c = 2 if gdename_std1 == "wannen a. a, r. schweizer co ag, wangen (0.3)" & gdenr_2018_4c == 200
replace wrong4c = 2 if gdename_std1 == "zuerich. arbo tuerenfabrik ag, buchs (0.06)" & gdenr_2018_4c == 4003
replace wrong4c = 2 if gdename_std1 == "zuerich. arbo tuerenfabrik ag, buchs (0.06)" & gdenr_2018_4c == 3271
replace wrong4c = 2 if gdename_std1 == "zuerich. arbo tuerenfabrik ag, buchs (0.06)" & gdenr_2018_4c == 83
replace wrong4c = 2 if gdename_std1 == "zuerich. arbo tuerenfabrik ag, buchs (0.06)" & gdenr_2018_4c == 1125
replace wrong4c = 2 if gdename_std1 == "zumikon rueegsegger ag, gossau (0.11)" & gdenr_2018_4c == 115
replace wrong4c = 2 if gdename_std1 == "zumikon rueegsegger ag, gossau (0.11)" & gdenr_2018_4c == 3443

replace gdenr_2018_4c = 241 if gdename_std1 == "aesch (zh)" & gdenr_2018_4c == 2761
replace gdenr_2018_4c = 241 if gdename_std1 == "aesch (zh)" & gdenr_2018_4c == 1021
replace gdenr_2018_4c = 6628 if gdename_std1 == "chapelle, route 40, gran lancy" & gdenr_2018_4c == 2009
replace gdenr_2018_4c = 6628 if gdename_std1 == "chapelle, route 40, gran lancy" & gdenr_2018_4c == 2066
replace gdenr_2018_4c = 5831 if gdename_std1 == "granges (vd)" & gdenr_2018_4c == 6248
replace gdenr_2018_4c = 5831 if gdename_std1 == "granges (vd)" & gdenr_2018_4c == 2328
replace gdenr_2018_4c = 2964 if gdename_std1 == "jun., stein, am rhein" & gdenr_2018_4c == 4260
replace gdenr_2018_4c = 2964 if gdename_std1 == "jun., stein, am rhein" & gdenr_2018_4c == 3360
replace gdenr_2018_4c = 2964 if gdename_std1 == "jun., stein, am rhein" & gdenr_2018_4c == 3005
replace gdenr_2018_4c = 902 if gdename_std1 == "langnau i. e, sesselbahn und skilift marbachegg ag, marbach (0.25)" & gdenr_2018_4c == 3253
replace gdenr_2018_4c = 902 if gdename_std1 == "langnau i. e, sesselbahn und skilift marbachegg ag, marbach (0.25)" & gdenr_2018_4c == 1010
replace gdenr_2018_4c = 5456 if gdename_std1 == "montet (vully)" & gdenr_2018_4c == 2050
replace gdenr_2018_4c = 5456 if gdename_std1 == "montet (vully)" & gdenr_2018_4c == 2089
replace gdenr_2018_4c = 4751 if gdename_std1 == "rickenbach (sg)" & gdenr_2018_4c == 2857
replace gdenr_2018_4c = 4751 if gdename_std1 == "rickenbach (sg)" & gdenr_2018_4c == 1097
replace gdenr_2018_4c = 4751 if gdename_std1 == "rickenbach (sg)" & gdenr_2018_4c == 2582
replace gdenr_2018_4c = 4751 if gdename_std1 == "rickenbach (sg)" & gdenr_2018_4c == 225
replace gdenr_2018_4c = 4751 if gdename_std1 == "rickenbach (tg)" & gdenr_2018_4c == 1097
replace gdenr_2018_4c = 4751 if gdename_std1 == "rickenbach (tg)" & gdenr_2018_4c == 225
replace gdenr_2018_4c = 4751 if gdename_std1 == "rickenbach (tg)" & gdenr_2018_4c == 2857
replace gdenr_2018_4c = 4751 if gdename_std1 == "rickenbach (tg)" & gdenr_2018_4c == 2582
replace gdenr_2018_4c = 1372 if gdename_std1 == "rickenbachstrasse 155, rickenbach" & gdenr_2018_4c == 2857
replace gdenr_2018_4c = 1372 if gdename_std1 == "rickenbachstrasse 155, rickenbach" & gdenr_2018_4c == 1097
replace gdenr_2018_4c = 1372 if gdename_std1 == "rickenbachstrasse 155, rickenbach" & gdenr_2018_4c == 2582
replace gdenr_2018_4c = 1372 if gdename_std1 == "rickenbachstrasse 155, rickenbach" & gdenr_2018_4c == 225
replace gdenr_2018_4c = 1372 if gdename_std1 == "schwyz , ag projekta ingenieurbuero, altdorf (0.08)" & gdenr_2018_4c == 2920
replace gdenr_2018_4c = 1372 if gdename_std1 == "schwyz , ag projekta ingenieurbuero, altdorf (0.08)" & gdenr_2018_4c == 1201
replace gdenr_2018_4c = 1372 if gdename_std1 == "sytirain, rickenbach" & gdenr_2018_4c == 2857
replace gdenr_2018_4c = 1372 if gdename_std1 == "sytirain, rickenbach" & gdenr_2018_4c == 1097
replace gdenr_2018_4c = 1372 if gdename_std1 == "sytirain, rickenbach" & gdenr_2018_4c == 225
replace gdenr_2018_4c = 1372 if gdename_std1 == "sytirain, rickenbach" & gdenr_2018_4c == 2582
replace gdenr_2018_4c = 1372 if gdename_std1 == "sytiraln, rickenbach" & gdenr_2018_4c == 225
replace gdenr_2018_4c = 1372 if gdename_std1 == "sytiraln, rickenbach" & gdenr_2018_4c == 1097
replace gdenr_2018_4c = 1372 if gdename_std1 == "sytiraln, rickenbach" & gdenr_2018_4c == 2582
replace gdenr_2018_4c = 1372 if gdename_std1 == "sytiraln, rickenbach" & gdenr_2018_4c == 2857
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, ersparniskasse des amtbezirks ag, wangen (1.0)" & gdenr_2018_4c == 1349
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, ersparniskasse des amtbezirks ag, wangen (1.0)" & gdenr_2018_4c == 200
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, esparniskasse des amtsbezirks ag, wangen (1.0)" & gdenr_2018_4c == 200
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, esparniskasse des amtsbezirks ag, wangen (1.0)" & gdenr_2018_4c == 1349
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, oberaareauische automobilkurse ag, wangen" & gdenr_2018_4c == 200
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, oberaareauische automobilkurse ag, wangen" & gdenr_2018_4c == 1349
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, oberaargauische automobikurse ag, wangen (0.105)" & gdenr_2018_4c == 200
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, oberaargauische automobikurse ag, wangen (0.105)" & gdenr_2018_4c == 1349
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, oberaargauische automobilkurse ag, wangen" & gdenr_2018_4c == 1349
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, oberaargauische automobilkurse ag, wangen" & gdenr_2018_4c == 200
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, obrecht co ag, wangen (0.8)" & gdenr_2018_4c == 200
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, obrecht co ag, wangen (0.8)" & gdenr_2018_4c == 1349
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, r. schweizer co ag, wangen (0.3)" & gdenr_2018_4c == 200
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a. a, r. schweizer co ag, wangen (0.3)" & gdenr_2018_4c == 1349
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a.a, r. schweizer co ag, wangen (0.5)" & gdenr_2018_4c == 200
replace gdenr_2018_4c = 992 if gdename_std1 == "wangen a.a, r. schweizer co ag, wangen (0.5)" & gdenr_2018_4c == 1349
replace gdenr_2018_4c = 992 if gdename_std1 == "wannen a. a, r. schweizer co ag, wangen (0.3)" & gdenr_2018_4c == 1349
replace gdenr_2018_4c = 992 if gdename_std1 == "wannen a. a, r. schweizer co ag, wangen (0.3)" & gdenr_2018_4c == 200
replace gdenr_2018_4c = 261 if gdename_std1 == "zuerich. arbo tuerenfabrik ag, buchs (0.06)" & gdenr_2018_4c == 4003
replace gdenr_2018_4c = 261 if gdename_std1 == "zuerich. arbo tuerenfabrik ag, buchs (0.06)" & gdenr_2018_4c == 3271
replace gdenr_2018_4c = 261 if gdename_std1 == "zuerich. arbo tuerenfabrik ag, buchs (0.06)" & gdenr_2018_4c == 83
replace gdenr_2018_4c = 261 if gdename_std1 == "zuerich. arbo tuerenfabrik ag, buchs (0.06)" & gdenr_2018_4c == 1125
replace gdenr_2018_4c = 160 if gdename_std1 == "zumikon rueegsegger ag, gossau (0.11)" & gdenr_2018_4c == 115
replace gdenr_2018_4c = 160 if gdename_std1 == "zumikon rueegsegger ag, gossau (0.11)" & gdenr_2018_4c == 3443

preserve // take geo-data and produce dataset with unique gdenr_2018_4a
use "$path\02_Processed_data\08_Municipalities\unique_Gdename_Gdenr2018-Geo.dta", clear
duplicates drop gdenr_2018, force
sort gdenr_2018
save "$path\02_Processed_data\10_Directors_1934_2003\unique_Gdenr2018-Geo.dta", replace
restore

rename gdenr_2018_4c gdenr_2018
sort gdenr_2018
merge m:1 gdenr_2018 using "$path\02_Processed_data\10_Directors_1934_2003\unique_Gdenr2018-Geo.dta"  
drop if _merge == 2
rename gdenr_2018 gdenr_2018_4c

replace gdename_4c = gdename if wrong4c == 2
replace E_CNTR_4c = E_CNTR if wrong4c == 2
replace N_CNTR_4c = N_CNTR if wrong4c == 2
drop gdename E_CNTR N_CNTR

duplicates report PID year ID_dupl gdenr_2018_4c  // CHECK!!!! (eliminated "n_obs_4")
duplicates drop PID year ID_dupl gdenr_2018_4c, force 
bysort PID year ID_dupl gdename_std1: gen options4c = _n   // number of options per PID-year-ID_dupl bundle

gen GdeNr_N_CNTR = N_CNTR_4c
gen GdeNr_E_CNTR = E_CNTR_4c
gen gdename_used = gdename_4c
gen gdenr_2018 = gdenr_2018_4c
gen countPLZ4 = .

// produce list with unique gde-info for manual checks
preserve
gen n = -1
collapse (sum) n, by (gdename_std1 gdenr_2018_4c gdename_4c)
order n 
sort n
save "$path\02_Processed_data\10_Directors_1934_2003\check_4c2.dta", replace
restore

keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
countPLZ4 GdeNr_E_CNTR GdeNr_N_CNTR gdename_used gdenr_2018 geo options4c

save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo43.dta", replace  


* Manual corrections 
//Missing matches with PLZ
use "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo4_notmerged.dta", clear
gen gdename = gdename_std1
replace PLZ4 = 3097 if PLZ4 == 3028 & gdename == "bem"
replace PLZ4 = 3097 if PLZ4 == 3028 & gdename == "liebefeld"
replace PLZ4 = 3038 if PLZ4 == 3028 & gdename == "oberlindach"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "s"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "s pegel bei bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "solegel"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiege bei bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegei bei bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel ."
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel / be"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel b. be"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel b. bem"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel b. bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel b. koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel b.b"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel b.bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel be"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei bem"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei bem/koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei berg"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei bern /koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei bern/ koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei bern/koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei bern/koenlz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei bern/koniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bei n"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bel bem"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bel bem/koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bel bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bel bern /koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bel bern /koenlz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bel bern/koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bel bern/koenlz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bet bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel bet bern/koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel koenitz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel koenlz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel/be"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel/bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel/koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegel/koenlz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spiegelb. bern"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "spirgrl bei koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "splegel koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "splegel koenlz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "splegel/koeniz"
replace PLZ4 = 3095 if PLZ4 == 3028 & gdename == "splegelb. bern"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == ".mutschellen"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "berikon mutschellen"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "berlkon mutschellen"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "muetschellen"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mulschellen"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutscheelen"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutscheilen"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "mutscheilen/zufikon"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutscheiten"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschelen"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschelfen"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschelien"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschelien/berikon"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschellem"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschellen"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschellen (ag)"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "mutschellen (zufikon)"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "mutschellen (zutikon)"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschellen ag"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschellen berikon"
replace PLZ4 = 8964 if PLZ4 == 8968 & gdename == "mutschellen friedlisberg"
replace PLZ4 = 8964 if PLZ4 == 8968 & gdename == "mutschellen rudolfstetten"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "mutschellen zufikon"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "mutschellen zuflkon"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschellen/"
replace PLZ4 = 8964 if PLZ4 == 8968 & gdename == "mutschellen/ rudolfstetten"
replace PLZ4 = 8964 if PLZ4 == 8968 & gdename == "mutschellen/ rudolfstetten friedlisberg"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschellen/berikon"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschellen/berlkon"
replace PLZ4 = 8964 if PLZ4 == 8968 & gdename == "mutschellen/rudolfstetten friedlisberg"
replace PLZ4 = 8967 if PLZ4 == 8968 & gdename == "mutschellen/widen"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "mutschellen/zufikon"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "mutschellen/zuflkon"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutscheller"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschellien"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutschelten"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "mutschelten/zufikon"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mutscheuen"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "mãºtschellen"
replace PLZ4 = 8966 if PLZ4 == 8968 & gdename == "oberwll lleli"
replace PLZ4 = 8964 if PLZ4 == 8968 & gdename == "rudotfstetten friedlisberg"
replace PLZ4 = 6966 if PLZ4 == 8968 & gdename == "villa luganoâ»"
replace PLZ4 = 8965 if PLZ4 == 8968 & gdename == "widen/mutschellen"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "zufikon mutschellen"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "zufikon/ mutschellen"
replace PLZ4 = 5621 if PLZ4 == 8968 & gdename == "zuflkon"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "marsken wildegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "merã¬ken wildegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "milken (ag)"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moe riken wildegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriiken wildegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken (ag)"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken ag"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken wadegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken waldegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken widegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken wiidegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken wldegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken wlidegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken wlldegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken wã¬ldegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moeriken wïldegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moerikon"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moerken (ag)"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moerlken"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moerlken (ag)"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moerlken ag"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moerlken wadegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moerlken waldegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moerlken wildegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moerlken wlldegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moriken (ag)"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moriken ag"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "moriken wildegg"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "morlken (ag)"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "mã³riken (ag)"
replace PLZ4 = 5415 if PLZ4 == 5115 & gdename == "obersiggental"
replace PLZ4 = 5103 if PLZ4 == 5115 & gdename == "wildegg"
replace PLZ4 = 1095 if PLZ4 == 1602 & gdename == "bossiere lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "croix"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "croix lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "croix sur lutry"
replace PLZ4 = 1095 if PLZ4 == 1602 & gdename == "escherin"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la chroix"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la chrolx"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la cr"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croise s lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix (luini)"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix (luny)"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix (lutry)"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix / lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix s lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix s/ lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix s/lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix sl lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix sur lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix(/lurty)"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix/lurty"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la croix/lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la crolx lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la crolx sur lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la crolx/lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "la crolxllutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "lacroix s/lutry"
replace PLZ4 = 1090 if PLZ4 == 1602 & gdename == "le croix sur lutry"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "aalenwinden"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "ailenwinden"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "allenwinden"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "allenwinden baar"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "allenwinden/baar"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "allenwinden/baer"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "allenwlnden"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "allenwlnden baar"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "allenwlnden saar"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "altenwinden"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "altenwinden baar"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "altenwinden saar"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "altenwinden/baar"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "altenwinden/baer"
replace PLZ4 = 6319 if PLZ4 == 6311 & gdename == "asienwinden"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "ediibach"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "edlibach"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "edlibach (menzingen)"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "edlibach menzingen"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "edllbach"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "finstersee"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "finstersee /menzingen"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "finstersee menzingen"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "finstersee/menzingen"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "flnstersee"
replace PLZ4 = 6313 if PLZ4 == 6311 & gdename == "flnstersee/menzingen"
replace PLZ4 = 6315 if PLZ4 == 6311 & gdename == "morgarten"
replace PLZ4 = 6315 if PLZ4 == 6311 & gdename == "morgarten oberaegeri"
replace PLZ4 = 6315 if PLZ4 == 6311 & gdename == "morgarten oberaegerl"
replace PLZ4 = 6314 if PLZ4 == 6311 & gdename == "neuaegeri"
replace PLZ4 = 6314 if PLZ4 == 6311 & gdename == "pffistersee"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "(cogne"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == ".groene"
replace PLZ4 = 3951 if PLZ4 == 3941 & gdename == "agam"
replace PLZ4 = 3951 if PLZ4 == 3941 & gdename == "agern"
replace PLZ4 = 3955 if PLZ4 == 3941 & gdename == "aibinen"
replace PLZ4 = 3955 if PLZ4 == 3941 & gdename == "alblnen"
replace PLZ4 = 3963 if PLZ4 == 3941 & gdename == "champzabe"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "cogne"
replace PLZ4 = 3963 if PLZ4 == 3941 & gdename == "crans"
replace PLZ4 = 3963 if PLZ4 == 3941 & gdename == "crans lens"
replace PLZ4 = 3963 if PLZ4 == 3941 & gdename == "crans s"
replace PLZ4 = 3963 if PLZ4 == 3941 & gdename == "crans s/lens"
replace PLZ4 = 3963 if PLZ4 == 3941 & gdename == "crans sur lens"
replace PLZ4 = 3963 if PLZ4 == 3941 & gdename == "crans/lens"
replace PLZ4 = 3943 if PLZ4 == 3941 & gdename == "eilscholl"
replace PLZ4 = 3943 if PLZ4 == 3941 & gdename == "eischola"
replace PLZ4 = 3943 if PLZ4 == 3941 & gdename == "elscholl"
replace PLZ4 = 3943 if PLZ4 == 3941 & gdename == "etscholl"
replace PLZ4 = 1978 if PLZ4 == 3941 & gdename == "fianthey"
replace PLZ4 = 1978 if PLZ4 == 3941 & gdename == "flantey"
replace PLZ4 = 1978 if PLZ4 == 3941 & gdename == "flanthey"
replace PLZ4 = 1978 if PLZ4 == 3941 & gdename == "flanthey lens"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "grane"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "grbne"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "grbno"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "greene"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "grene"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "grime"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "grmne"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "groene"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "grã³ne"
replace PLZ4 = 3956 if PLZ4 == 3941 & gdename == "gullet"
replace PLZ4 = 3956 if PLZ4 == 3941 & gdename == "gutten"
replace PLZ4 = 3956 if PLZ4 == 3941 & gdename == "gutten vs"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "lcogne"
replace PLZ4 = 3979 if PLZ4 == 3941 & gdename == "lcogre"
replace PLZ4 = 3960 if PLZ4 == 3941 & gdename == "no"
replace PLZ4 = 3960 if PLZ4 == 3941 & gdename == "noaes/sierre"
replace PLZ4 = 3960 if PLZ4 == 3941 & gdename == "nobs"
replace PLZ4 = 3960 if PLZ4 == 3941 & gdename == "noees"
replace PLZ4 = 3960 if PLZ4 == 3941 & gdename == "noes"
replace PLZ4 = 3960 if PLZ4 == 3941 & gdename == "noes sierre"
replace PLZ4 = 3960 if PLZ4 == 3941 & gdename == "noes/sierre"
replace PLZ4 = 3960 if PLZ4 == 3941 & gdename == "noes/slerre"
replace PLZ4 = 3960 if PLZ4 == 3941 & gdename == "nooes/sierre"
replace PLZ4 = 3942 if PLZ4 == 3941 & gdename == "rarogne"
replace PLZ4 = 3944 if PLZ4 == 3941 & gdename == "unterbach"
replace PLZ4 = 3944 if PLZ4 == 3941 & gdename == "unterbaech vs"
replace PLZ4 = 3944 if PLZ4 == 3941 & gdename == "unterbaechn"
replace PLZ4 = 8132 if PLZ4 == 8128 & gdename == "hinteregg"
replace PLZ4 = 8132 if PLZ4 == 8128 & gdename == "hinteregg bei zuerich"
replace PLZ4 = 8132 if PLZ4 == 8128 & gdename == "hinteregg egg"
replace PLZ4 = 8132 if PLZ4 == 8128 & gdename == "hinteregg zh"
replace PLZ4 = 8132 if PLZ4 == 8128 & gdename == "hinteregg/egg bei zuerich"
replace PLZ4 = 8132 if PLZ4 == 8128 & gdename == "hinteregg/egg bel zuerich"
replace PLZ4 = 8132 if PLZ4 == 8128 & gdename == "hinteregg/egg zh"
replace PLZ4 = 8132 if PLZ4 == 8128 & gdename == "hinterwegg"
replace PLZ4 = 8132 if PLZ4 == 8128 & gdename == "hlnteregg"
replace PLZ4 = 8132 if PLZ4 == 8128 & gdename == "hlnteregg zh"
replace PLZ4 = 8126 if PLZ4 == 8128 & gdename == "zumlkon"
replace PLZ4 = 1454 if PLZ4 == 1451 & gdename == "auberson"
replace PLZ4 = 1454 if PLZ4 == 1451 & gdename == "l auberson"
replace PLZ4 = 1454 if PLZ4 == 1451 & gdename == "l auberson /ste croix"
replace PLZ4 = 1450 if PLZ4 == 1451 & gdename == "la bagne vd"
replace PLZ4 = 1450 if PLZ4 == 1451 & gdename == "la sagne pres ste croix"
replace PLZ4 = 1450 if PLZ4 == 1451 & gdename == "la sagne vd"
replace PLZ4 = 1450 if PLZ4 == 1451 & gdename == "la segne vd"
replace PLZ4 = 1452 if PLZ4 == 1451 & gdename == "lea rassea"
replace PLZ4 = 1452 if PLZ4 == 1451 & gdename == "les basses"
replace PLZ4 = 1453 if PLZ4 == 1451 & gdename == "les cluds"
replace PLZ4 = 1453 if PLZ4 == 1451 & gdename == "les cluds riere bullet"
replace PLZ4 = 1452 if PLZ4 == 1451 & gdename == "les r a sse s"
replace PLZ4 = 1452 if PLZ4 == 1451 & gdename == "les rasses"
replace PLZ4 = 1452 if PLZ4 == 1451 & gdename == "les rassies"
replace PLZ4 = 6939 if PLZ4 == 6911 & gdename == "aroslo"
replace PLZ4 = 6930 if PLZ4 == 6911 & gdename == "badano"
replace PLZ4 = 6917 if PLZ4 == 6911 & gdename == "barbengo/lugano"
replace PLZ4 = 6937 if PLZ4 == 6911 & gdename == "brano"
replace PLZ4 = 6827 if PLZ4 == 6911 & gdename == "brasino arsizio"
replace PLZ4 = 6979 if PLZ4 == 6911 & gdename == "bre"
replace PLZ4 = 6979 if PLZ4 == 6911 & gdename == "bre di lugano"
replace PLZ4 = 6979 if PLZ4 == 6911 & gdename == "bre s/lugano"
replace PLZ4 = 6979 if PLZ4 == 6911 & gdename == "bre sopra lugano"
replace PLZ4 = 6979 if PLZ4 == 6911 & gdename == "brit dl lugano"
replace PLZ4 = 6827 if PLZ4 == 6911 & gdename == "brusino"
replace PLZ4 = 6827 if PLZ4 == 6911 & gdename == "brusino arsizlo"
replace PLZ4 = 6827 if PLZ4 == 6911 & gdename == "brusino arslzlo"
replace PLZ4 = 6827 if PLZ4 == 6911 & gdename == "brusio arizo"
replace PLZ4 = 6827 if PLZ4 == 6911 & gdename == "bruslno arsizio"
replace PLZ4 = 6827 if PLZ4 == 6911 & gdename == "brusã¬no arsizio"
replace PLZ4 = 6827 if PLZ4 == 6911 & gdename == "busino arsizio"
replace PLZ4 = 6913 if PLZ4 == 6911 & gdename == "carabbio"
replace PLZ4 = 6913 if PLZ4 == 6911 & gdename == "carabbla"
replace PLZ4 = 6913 if PLZ4 == 6911 & gdename == "carabblia"
replace PLZ4 = 6926 if PLZ4 == 6911 & gdename == "carabietta montagnoia"
replace PLZ4 = 6926 if PLZ4 == 6911 & gdename == "carabietta montagnola"
replace PLZ4 = 6926 if PLZ4 == 6911 & gdename == "carabletta"
replace PLZ4 = 6926 if PLZ4 == 6911 & gdename == "carabletta montagnola"
replace PLZ4 = 6913 if PLZ4 == 6911 & gdename == "carrabia"
replace PLZ4 = 6949 if PLZ4 == 6911 & gdename == "ccmano"
replace PLZ4 = 6949 if PLZ4 == 6911 & gdename == "comano /ti"
replace PLZ4 = 6949 if PLZ4 == 6911 & gdename == "coniano"
replace PLZ4 = 6949 if PLZ4 == 6911 & gdename == "cumano"
replace PLZ4 = 6938 if PLZ4 == 6911 & gdename == "fescoggla"
replace PLZ4 = 6926 if PLZ4 == 6911 & gdename == "gentillno"
replace PLZ4 = 6916 if PLZ4 == 6911 & gdename == "grancla"
replace PLZ4 = 6916 if PLZ4 == 6911 & gdename == "grand e"
replace PLZ4 = 6916 if PLZ4 == 6911 & gdename == "grande"
replace PLZ4 = 6916 if PLZ4 == 6911 & gdename == "granola"
replace PLZ4 = 6916 if PLZ4 == 6911 & gdename == "grenela"
replace PLZ4 = 6929 if PLZ4 == 6911 & gdename == "grumo di gravesano"
replace PLZ4 = 6929 if PLZ4 == 6911 & gdename == "grumo dl gravesano"
replace PLZ4 = 6929 if PLZ4 == 6911 & gdename == "lugano gravesano"
replace PLZ4 = 6928 if PLZ4 == 6911 & gdename == "manna"
replace PLZ4 = 6926 if PLZ4 == 6911 & gdename == "montagnola carabietta"
replace PLZ4 = 6926 if PLZ4 == 6911 & gdename == "montagnola carabletta"
replace PLZ4 = 6900 if PLZ4 == 6911 & gdename == "monte bre lugano"
replace PLZ4 = 6900 if PLZ4 == 6911 & gdename == "monte bre/lugano"
replace PLZ4 = 6915 if PLZ4 == 6911 & gdename == "noracnco"
replace PLZ4 = 6915 if PLZ4 == 6911 & gdename == "noranco"
replace PLZ4 = 6915 if PLZ4 == 6911 & gdename == "pambio"
replace PLZ4 = 6915 if PLZ4 == 6911 & gdename == "pambio norano"
replace PLZ4 = 6915 if PLZ4 == 6911 & gdename == "pamblo"
replace PLZ4 = 6915 if PLZ4 == 6911 & gdename == "pamblo noranco"
replace PLZ4 = 6915 if PLZ4 == 6911 & gdename == "pampio noranco"
replace PLZ4 = 6930 if PLZ4 == 6911 & gdename == "sedano"
replace PLZ4 = 6938 if PLZ4 == 6911 & gdename == "vezlo"
replace PLZ4 = 6921 if PLZ4 == 6911 & gdename == "vice morcote"
replace PLZ4 = 6921 if PLZ4 == 6911 & gdename == "vico miorcote"
replace PLZ4 = 6921 if PLZ4 == 6911 & gdename == "vida morcote"
replace PLZ4 = 6921 if PLZ4 == 6911 & gdename == "vlco morcote"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "8t leonard"
replace PLZ4 = 1966 if PLZ4 == 3958 & gdename == "signese st leonhard"
replace PLZ4 = 1966 if PLZ4 == 3958 & gdename == "slignese st leonhard"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "st leonard"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "st leonard uvrier"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "st leonard/sion"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "st leonard/slon"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "st leonhard"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "st. leonard"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "st. leonhard"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvrier"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvrier gjslon"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvrier s/sion"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvrier sion"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvrier st leonard"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvrier st leonhard"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvrier/sion"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvrier/st leonard"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvrier/st. leonard"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvriersion"
replace PLZ4 = 1958 if PLZ4 == 3958 & gdename == "uvrïer/st leonard"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipf 0berfrick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipf oberf rick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipf oberflick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipf oberfrack"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipf oberfriek"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipf oberfrlck"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipf oberirlck"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipf obertrick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipf.obertrick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipt oberfrack"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gipt oberfrick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "glpf bberfrick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "glpf oberfrack"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "glpf oberfrick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "glpf oberfrlck"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "glpf oberfrlek"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "glpf oberirick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "glpf obertrick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "gã¬pf oberfrick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "oberf rick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "oberfrick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "oberfrlck"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "obertrick"
replace PLZ4 = 5073 if PLZ4 == 5264 & gdename == "obertrickâ"
replace PLZ4 = 7265 if PLZ4 == 7299 & gdename == "davos laret"
replace PLZ4 = 7235 if PLZ4 == 7299 & gdename == "fideris dorf"
replace PLZ4 = 7235 if PLZ4 == 7299 & gdename == "fideris station"
replace PLZ4 = 7235 if PLZ4 == 7299 & gdename == "fiderls dorf"
replace PLZ4 = 7235 if PLZ4 == 7299 & gdename == "fiderls station"
replace PLZ4 = 7235 if PLZ4 == 7299 & gdename == "flderis dorf"
replace PLZ4 = 7235 if PLZ4 == 7299 & gdename == "flderla"
replace PLZ4 = 7235 if PLZ4 == 7299 & gdename == "flderls station"
replace PLZ4 = 7232 if PLZ4 == 7299 & gdename == "fuma"
replace PLZ4 = 7232 if PLZ4 == 7299 & gdename == "fuma station"
replace PLZ4 = 7232 if PLZ4 == 7299 & gdename == "fume"
replace PLZ4 = 7232 if PLZ4 == 7299 & gdename == "furna dorf"
replace PLZ4 = 7232 if PLZ4 == 7299 & gdename == "furna dort"
replace PLZ4 = 7232 if PLZ4 == 7299 & gdename == "furna station"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saas 1. praettigau"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saas i. p"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saas i. praettigau"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saas i. prattigau"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saas i.p"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saas im praettigau"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saas praedigau"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saas praedlgau"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saas praettigau"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saes"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saes i. praettigau"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "saes im praettigau"
replace PLZ4 = 7493 if PLZ4 == 7299 & gdename == "schmitten seewis"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "seas"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "seas i. praettigau"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "seas im praettigau"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewcs"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewcs dorf"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewcs dort"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewcs station"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewis"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewis dorf"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewis i. pr."
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewis i. praettigaue"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewis im praettigeu"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewis station"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewls"
replace PLZ4 = 7212 if PLZ4 == 7299 & gdename == "seewls dorf"
replace PLZ4 = 4411 if PLZ4 == 7299 & gdename == "seitisberg"
replace PLZ4 = 7250 if PLZ4 == 7299 & gdename == "serneus"
replace PLZ4 = 7260 if PLZ4 == 7299 & gdename == "sertig"
replace PLZ4 = 7265 if PLZ4 == 7299 & gdename == "wolfgang"
replace PLZ4 = 7265 if PLZ4 == 7299 & gdename == "wolfgang bei davos"
replace PLZ4 = 7265 if PLZ4 == 7299 & gdename == "wolfgang bel davos"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "altstaetten/luechingen"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luchingen"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luchingen (sg)"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luchingen sg"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luechangen"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luechingen"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luechingen (sg)"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luechingen/ altstaetten"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luechingen/ altstaetten sg"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luechingen/altstaetten"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luechingen/altstaetten sg"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luechlngen"
replace PLZ4 = 9450 if PLZ4 == 9438 & gdename == "luechlngen (sg)"
replace PLZ4 = 6810 if PLZ4 == 6849 & gdename == "leone"
replace PLZ4 = 6810 if PLZ4 == 6849 & gdename == "lsone"
replace PLZ4 = 6809 if PLZ4 == 6849 & gdename == "medegiia"
replace PLZ4 = 6809 if PLZ4 == 6849 & gdename == "medeglia manz isolazioni sa manne"
replace PLZ4 = 6809 if PLZ4 == 6849 & gdename == "medeglla"
replace PLZ4 = 6805 if PLZ4 == 6849 & gdename == "meuovico"
replace PLZ4 = 6805 if PLZ4 == 6849 & gdename == "mezzovicco"
replace PLZ4 = 6805 if PLZ4 == 6849 & gdename == "mezzovico"
replace PLZ4 = 6805 if PLZ4 == 6849 & gdename == "mezzovlco"
replace PLZ4 = 6822 if PLZ4 == 6849 & gdename == "pugerna"
replace PLZ4 = 6821 if PLZ4 == 6849 & gdename == "rovia"
replace PLZ4 = 6821 if PLZ4 == 6849 & gdename == "roviio"
replace PLZ4 = 6821 if PLZ4 == 6849 & gdename == "rovio ti"
replace PLZ4 = 6821 if PLZ4 == 6849 & gdename == "rovioi"
replace PLZ4 = 6821 if PLZ4 == 6849 & gdename == "rovlo"
replace PLZ4 = 6806 if PLZ4 == 6849 & gdename == "sigirino ti"
replace PLZ4 = 6806 if PLZ4 == 6849 & gdename == "siglrino"
replace PLZ4 = 6806 if PLZ4 == 6849 & gdename == "slgirino"
replace PLZ4 = 6806 if PLZ4 == 6849 & gdename == "slglrino"
replace PLZ4 = 1123 if PLZ4 == 1111 & gdename == "aciens"
replace PLZ4 = 1121 if PLZ4 == 1111 & gdename == "brembiens"
replace PLZ4 = 1121 if PLZ4 == 1111 & gdename == "bremlens"
replace PLZ4 = 1121 if PLZ4 == 1111 & gdename == "brernblens"
replace PLZ4 = 1127 if PLZ4 == 1111 & gdename == "ciarmont"
replace PLZ4 = 1127 if PLZ4 == 1111 & gdename == "clermont"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "colombier s/morges"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "colombier sur marges"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "colombier sur morges"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "colombler sur morges"
replace PLZ4 = 1116 if PLZ4 == 1111 & gdename == "cottens vaud"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "echiens"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "echlchens"
replace PLZ4 = 1134 if PLZ4 == 1111 & gdename == "en chantemerte"
replace PLZ4 = 1124 if PLZ4 == 1111 & gdename == "goiiion"
replace PLZ4 = 1124 if PLZ4 == 1111 & gdename == "goilion"
replace PLZ4 = 1124 if PLZ4 == 1111 & gdename == "golfion"
replace PLZ4 = 1124 if PLZ4 == 1111 & gdename == "goliion"
replace PLZ4 = 1124 if PLZ4 == 1111 & gdename == "golllon"
replace PLZ4 = 1117 if PLZ4 == 1111 & gdename == "grancy vd"
replace PLZ4 = 1132 if PLZ4 == 1111 & gdename == "luffy (vd)"
replace PLZ4 = 1132 if PLZ4 == 1111 & gdename == "luily"
replace PLZ4 = 1132 if PLZ4 == 1111 & gdename == "luliy"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "luliy sur morges"
replace PLZ4 = 1132 if PLZ4 == 1111 & gdename == "lullt vd"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lully morges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lully s marges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lully s morges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lully s. morges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lully s/morges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lully sur mordes"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lully sur morges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lully/morges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lusey sur morges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lussy s morges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lussy s. morges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lussy s/morges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lussy sur moorges"
replace PLZ4 = 1167 if PLZ4 == 1111 & gdename == "lussy sur morgen"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "monnaz sur morges"
replace PLZ4 = 1128 if PLZ4 == 1111 & gdename == "reveroile"
replace PLZ4 = 1128 if PLZ4 == 1111 & gdename == "reverolie"
replace PLZ4 = 1122 if PLZ4 == 1111 & gdename == "romane) s morges"
replace PLZ4 = 1032 if PLZ4 == 1111 & gdename == "romanei sur lausanne"
replace PLZ4 = 1122 if PLZ4 == 1111 & gdename == "romanei sur morges"
replace PLZ4 = 1122 if PLZ4 == 1111 & gdename == "romanel s morges"
replace PLZ4 = 1304 if PLZ4 == 1111 & gdename == "senarciens"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "st saphorin"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "st saphorin (merges)"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "st saphorin (morges)"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "st saphorin morges"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "st saphorin s morges"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "st saphorin s/morges"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "st saphorin sur morges"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "st saphorin sur morses"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "st saphorã­n sur morges"
replace PLZ4 = 1112 if PLZ4 == 1111 & gdename == "st. saphorin sur morges"
replace PLZ4 = 1131 if PLZ4 == 1111 & gdename == "toiochena"
replace PLZ4 = 1131 if PLZ4 == 1111 & gdename == "tolochena"
replace PLZ4 = 1115 if PLZ4 == 1111 & gdename == "vuillerens"
replace PLZ4 = 1115 if PLZ4 == 1111 & gdename == "vuliierens"
replace PLZ4 = 8625 if PLZ4 == 8611 & gdename == "bertschikon (gossau zh)"
replace PLZ4 = 8625 if PLZ4 == 8611 & gdename == "bertschikon b. gossau zh"
replace PLZ4 = 8625 if PLZ4 == 8611 & gdename == "bertschikon b. uster"
replace PLZ4 = 8625 if PLZ4 == 8611 & gdename == "bertschikon gossau"
replace PLZ4 = 8625 if PLZ4 == 8611 & gdename == "bertschikon zh"
replace PLZ4 = 8625 if PLZ4 == 8611 & gdename == "bertschlkon"
replace PLZ4 = 8625 if PLZ4 == 8611 & gdename == "bertschlkon (gossau zh)"
replace PLZ4 = 8625 if PLZ4 == 8611 & gdename == "bertschlkon gossau"
replace PLZ4 = 8625 if PLZ4 == 8611 & gdename == "bertschlkon zh"
replace PLZ4 = 8610 if PLZ4 == 8611 & gdename == "freudwil"
replace PLZ4 = 8610 if PLZ4 == 8611 & gdename == "freudwll"
replace PLZ4 = 8625 if PLZ4 == 8611 & gdename == "gossau bertschikon"
replace PLZ4 = 8610 if PLZ4 == 8611 & gdename == "riedikon"
replace PLZ4 = 8614 if PLZ4 == 8611 & gdename == "sulzbach"
replace PLZ4 = 8614 if PLZ4 == 8611 & gdename == "sulzbach b. uster"
replace PLZ4 = 8614 if PLZ4 == 8611 & gdename == "sulzbach bei uster"
replace PLZ4 = 8615 if PLZ4 == 8611 & gdename == "wermatsw l"
replace PLZ4 = 8615 if PLZ4 == 8611 & gdename == "wermatswii uster"
replace PLZ4 = 8615 if PLZ4 == 8611 & gdename == "wermatswil"
replace PLZ4 = 8615 if PLZ4 == 8611 & gdename == "wermatswil uster"
replace PLZ4 = 8615 if PLZ4 == 8611 & gdename == "wermatswil/uster"
replace PLZ4 = 8615 if PLZ4 == 8611 & gdename == "wermatswll"
replace PLZ4 = 8615 if PLZ4 == 8611 & gdename == "wermatswll uster"
replace PLZ4 = 8615 if PLZ4 == 8611 & gdename == "wermetswil"
replace PLZ4 = 7416 if PLZ4 == 7499 & gdename == "almens gr"
replace PLZ4 = 7492 if PLZ4 == 7499 & gdename == "alvaneu bad"
replace PLZ4 = 7492 if PLZ4 == 7499 & gdename == "alvaneu dorf"
replace PLZ4 = 7492 if PLZ4 == 7499 & gdename == "alveneu"
replace PLZ4 = 7492 if PLZ4 == 7499 & gdename == "alveneu bad"
replace PLZ4 = 7416 if PLZ4 == 7499 & gdename == "atmens"
replace PLZ4 = 7408 if PLZ4 == 7499 & gdename == "cazls"
replace PLZ4 = 7404 if PLZ4 == 7499 & gdename == "feidis veulden"
replace PLZ4 = 7404 if PLZ4 == 7499 & gdename == "feldas veulden"
replace PLZ4 = 7404 if PLZ4 == 7499 & gdename == "feldas/veulden"
replace PLZ4 = 7404 if PLZ4 == 7499 & gdename == "feldis"
replace PLZ4 = 7404 if PLZ4 == 7499 & gdename == "feldis veulden"
replace PLZ4 = 7404 if PLZ4 == 7499 & gdename == "feldisneulden"
replace PLZ4 = 7404 if PLZ4 == 7499 & gdename == "feldisrveulden"
replace PLZ4 = 7414 if PLZ4 == 7499 & gdename == "fuerstenaubruck"
replace PLZ4 = 7414 if PLZ4 == 7499 & gdename == "furstenaubruck"
replace PLZ4 = 7403 if PLZ4 == 7499 & gdename == "rhaezuns"
replace PLZ4 = 7403 if PLZ4 == 7499 & gdename == "rhazuens"
replace PLZ4 = 7403 if PLZ4 == 7499 & gdename == "rhazuns"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "s is i. domleschg"
replace PLZ4 = 7493 if PLZ4 == 7499 & gdename == "schmitten/albula"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sii i.d"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sil i"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sil i.d"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils 1. domleschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils 1.domieschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils domieschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils domleschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils i. d"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils i. domieschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils i. domleschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils i.d"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils i.domieschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils i.domleschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils im domieschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils im domleechg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils l.domieschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sils ã­. d"
replace PLZ4 = 7514 if PLZ4 == 7499 & gdename == "sils/e"
replace PLZ4 = 7514 if PLZ4 == 7499 & gdename == "sils/sega maria"
replace PLZ4 = 7514 if PLZ4 == 7499 & gdename == "sils/segl maria"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "slls i. domieschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "slls i.domieschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "slls i.domleschg"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "stuls"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "sã­ls i.domleschg"
replace PLZ4 = 7418 if PLZ4 == 7499 & gdename == "tumegi/tomils"
replace PLZ4 = 7418 if PLZ4 == 7499 & gdename == "tumegl"
replace PLZ4 = 7418 if PLZ4 == 7499 & gdename == "tumegl.tomils"
replace PLZ4 = 7418 if PLZ4 == 7499 & gdename == "tumegl/tamils"
replace PLZ4 = 7418 if PLZ4 == 7499 & gdename == "tumegl/tomiis"
replace PLZ4 = 7418 if PLZ4 == 7499 & gdename == "tumegl/tomlls"
replace PLZ4 = 7418 if PLZ4 == 7499 & gdename == "tumeglitamils"
replace PLZ4 = 7411 if PLZ4 == 7499 & gdename == "us im domlesahg"
replace PLZ4 = 1946 if PLZ4 == 1931 & gdename == "bourg st pierre"
replace PLZ4 = 1946 if PLZ4 == 1931 & gdename == "bourg st plerre"
replace PLZ4 = 1946 if PLZ4 == 1931 & gdename == "bourg st. pierre"
replace PLZ4 = 1932 if PLZ4 == 1931 & gdename == "bovemier"
replace PLZ4 = 1932 if PLZ4 == 1931 & gdename == "bovemler"
replace PLZ4 = 1932 if PLZ4 == 1931 & gdename == "bovernler"
replace PLZ4 = 1934 if PLZ4 == 1931 & gdename == "champsec"
replace PLZ4 = 1937 if PLZ4 == 1931 & gdename == "la fouiy"
replace PLZ4 = 1937 if PLZ4 == 1931 & gdename == "la fouly"
replace PLZ4 = 1941 if PLZ4 == 1931 & gdename == "leveron"
replace PLZ4 = 1945 if PLZ4 == 1931 & gdename == "liddes vs"
replace PLZ4 = 1937 if PLZ4 == 1931 & gdename == "praz de fort"
replace PLZ4 = 1934 if PLZ4 == 1931 & gdename == "sarreyer"
replace PLZ4 = 1934 if PLZ4 == 1931 & gdename == "v rsegeres"
replace PLZ4 = 1934 if PLZ4 == 1931 & gdename == "vergenes bagnes"
replace PLZ4 = 1934 if PLZ4 == 1931 & gdename == "versegere"
replace PLZ4 = 1934 if PLZ4 == 1931 & gdename == "versegeres"
replace PLZ4 = 1934 if PLZ4 == 1931 & gdename == "versegeres / prarreyer"
replace PLZ4 = 1934 if PLZ4 == 1931 & gdename == "versegeres i prarreyer"
replace PLZ4 = 1934 if PLZ4 == 1931 & gdename == "versegeres prarreyer"
replace PLZ4 = 1941 if PLZ4 == 1931 & gdename == "voileges"
replace PLZ4 = 1941 if PLZ4 == 1931 & gdename == "volieges"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "mtn st pierre"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "viilaz st pierre"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "viliaz st pierre"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "viliaz st. pierre"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "villan st pierre"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "villaz lussy"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "villaz st pierre"
replace PLZ4 = 1694 if PLZ4 == 1758 & gdename == "villaz st pierre/chavannes sous orsonnens"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "villaz st plerre"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "villaz st. pierre"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "vlllaz st pierre"
replace PLZ4 = 1690 if PLZ4 == 1758 & gdename == "vlllaz st plerre"
replace PLZ4 = 7019 if PLZ4 == 7099 & gdename == "fidaz"
replace PLZ4 = 7019 if PLZ4 == 7099 & gdename == "fidaz films"
replace PLZ4 = 7019 if PLZ4 == 7099 & gdename == "fidaz flims"
replace PLZ4 = 7019 if PLZ4 == 7099 & gdename == "fidaz/flims"
replace PLZ4 = 7019 if PLZ4 == 7099 & gdename == "fldaz films"
replace PLZ4 = 7019 if PLZ4 == 7099 & gdename == "fldaz flims"
replace PLZ4 = 7050 if PLZ4 == 7099 & gdename == "langwien"
replace PLZ4 = 7050 if PLZ4 == 7099 & gdename == "langwies gr"
replace PLZ4 = 7050 if PLZ4 == 7099 & gdename == "langwles"
replace PLZ4 = 7083 if PLZ4 == 7099 & gdename == "lantsch"
replace PLZ4 = 7083 if PLZ4 == 7099 & gdename == "lantsch lenz"
replace PLZ4 = 7083 if PLZ4 == 7099 & gdename == "lantsch%lenz"
replace PLZ4 = 7083 if PLZ4 == 7099 & gdename == "lantsch. lenz"
replace PLZ4 = 7083 if PLZ4 == 7099 & gdename == "lantsch/l.crnz"
replace PLZ4 = 7083 if PLZ4 == 7099 & gdename == "lantsch/lerz"
replace PLZ4 = 7083 if PLZ4 == 7099 & gdename == "lantschilenz"
replace PLZ4 = 7050 if PLZ4 == 7099 & gdename == "litzirueti"
replace PLZ4 = 7074 if PLZ4 == 7099 & gdename == "malfix"
replace PLZ4 = 7082 if PLZ4 == 7099 & gdename == "muldain vaz"
replace PLZ4 = 7082 if PLZ4 == 7099 & gdename == "obervaz"
replace PLZ4 = 7050 if PLZ4 == 7099 & gdename == "paglg"
replace PLZ4 = 7050 if PLZ4 == 7099 & gdename == "pelst"
replace PLZ4 = 7015 if PLZ4 == 7099 & gdename == "reichenau"
replace PLZ4 = 7050 if PLZ4 == 7099 & gdename == "sapuen"
replace PLZ4 = 7014 if PLZ4 == 7099 & gdename == "trier mulin"
replace PLZ4 = 7014 if PLZ4 == 7099 & gdename == "triin"
replace PLZ4 = 7014 if PLZ4 == 7099 & gdename == "trin digg"
replace PLZ4 = 7014 if PLZ4 == 7099 & gdename == "trin mulin"
replace PLZ4 = 7014 if PLZ4 == 7099 & gdename == "trin mulln"
replace PLZ4 = 7014 if PLZ4 == 7099 & gdename == "trin mutin"
replace PLZ4 = 7082 if PLZ4 == 7099 & gdename == "vaz"
replace PLZ4 = 7082 if PLZ4 == 7099 & gdename == "vaz obervaz"
replace PLZ4 = 7082 if PLZ4 == 7099 & gdename == "vaziobervaz"
replace PLZ4 = 7082 if PLZ4 == 7099 & gdename == "zorten"
replace PLZ4 = 7082 if PLZ4 == 7099 & gdename == "zorten gde. vaz ober"
replace PLZ4 = 7082 if PLZ4 == 7099 & gdename == "zorten gde. vaz obervaz"
replace PLZ4 = 8248 if PLZ4 == 8448 & gdename == "laufen uhwlesen"
replace PLZ4 = 8248 if PLZ4 == 8448 & gdename == "lauten uhwlesen"
replace PLZ4 = 8248 if PLZ4 == 8448 & gdename == "uhwieaen"
replace PLZ4 = 8248 if PLZ4 == 8448 & gdename == "uhwiesen"
replace PLZ4 = 8248 if PLZ4 == 8448 & gdename == "uhwiesen laufen"
replace PLZ4 = 8248 if PLZ4 == 8448 & gdename == "uhwlesen"
replace PLZ4 = 4229 if PLZ4 == 4249 & gdename == "beinwal (so)"
replace PLZ4 = 4229 if PLZ4 == 4249 & gdename == "beinwal so"
replace PLZ4 = 4229 if PLZ4 == 4249 & gdename == "beinwll"
replace PLZ4 = 4229 if PLZ4 == 4249 & gdename == "belnwil (so)"
replace PLZ4 = 4223 if PLZ4 == 4249 & gdename == "blauen be"
replace PLZ4 = 4232 if PLZ4 == 4249 & gdename == "fahren"
replace PLZ4 = 4247 if PLZ4 == 4249 & gdename == "grinde)"
replace PLZ4 = 4247 if PLZ4 == 4249 & gdename == "grindel so"
replace PLZ4 = 4204 if PLZ4 == 4249 & gdename == "hammelried"
replace PLZ4 = 4204 if PLZ4 == 4249 & gdename == "hammelried so"
replace PLZ4 = 4204 if PLZ4 == 4249 & gdename == "hilmmeiried"
replace PLZ4 = 4204 if PLZ4 == 4249 & gdename == "himmeiried"
replace PLZ4 = 4204 if PLZ4 == 4249 & gdename == "himmelrfed"
replace PLZ4 = 4204 if PLZ4 == 4249 & gdename == "himmelried so"
replace PLZ4 = 4204 if PLZ4 == 4249 & gdename == "himmelrled"
replace PLZ4 = 4204 if PLZ4 == 4249 & gdename == "himmelrled so"
replace PLZ4 = 4204 if PLZ4 == 4249 & gdename == "hlmmeiried"
replace PLZ4 = 4204 if PLZ4 == 4249 & gdename == "hlmmelried"
replace PLZ4 = 4233 if PLZ4 == 4249 & gdename == "meitingen"
replace PLZ4 = 4233 if PLZ4 == 4249 & gdename == "mettingen"
replace PLZ4 = 4208 if PLZ4 == 4249 & gdename == "nenziingen"
replace PLZ4 = 4229 if PLZ4 == 4249 & gdename == "oberbeinwil"
replace PLZ4 = 4229 if PLZ4 == 4249 & gdename == "unterbeinwil"
replace PLZ4 = 4229 if PLZ4 == 4249 & gdename == "unterbelnwil"
replace PLZ4 = 4246 if PLZ4 == 4249 & gdename == "wahiern"
replace PLZ4 = 4246 if PLZ4 == 4249 & gdename == "wahlen b."
replace PLZ4 = 4246 if PLZ4 == 4249 & gdename == "wahlen b. laufen"
replace PLZ4 = 4246 if PLZ4 == 4249 & gdename == "wahlen bei"
replace PLZ4 = 4246 if PLZ4 == 4249 & gdename == "wahlen bei laufen"
replace PLZ4 = 4246 if PLZ4 == 4249 & gdename == "wahlen bel laufen"
replace PLZ4 = 4234 if PLZ4 == 4249 & gdename == "zuliwil"
replace PLZ4 = 4234 if PLZ4 == 4249 & gdename == "zuliwll"
replace PLZ4 = 4234 if PLZ4 == 4249 & gdename == "zutiwii"
replace PLZ4 = 9423 if PLZ4 == 9499 & gdename == "altenrhein"
replace PLZ4 = 9467 if PLZ4 == 9499 & gdename == "fruemsen"
replace PLZ4 = 9467 if PLZ4 == 9499 & gdename == "fruemsen sennwald"
replace PLZ4 = 9469 if PLZ4 == 9499 & gdename == "haag"
replace PLZ4 = 9469 if PLZ4 == 9499 & gdename == "haag (rheintal)"
replace PLZ4 = 9469 if PLZ4 == 9499 & gdename == "haag (sennwaid)"
replace PLZ4 = 9469 if PLZ4 == 9499 & gdename == "haag (sennwald)"
replace PLZ4 = 9469 if PLZ4 == 9499 & gdename == "haag (sg)"
replace PLZ4 = 9469 if PLZ4 == 9499 & gdename == "haag sennwald"
replace PLZ4 = 9469 if PLZ4 == 9499 & gdename == "haag sg"
replace PLZ4 = 9468 if PLZ4 == 9499 & gdename == "sax"
replace PLZ4 = 9427 if PLZ4 == 9499 & gdename == "zeig"
replace PLZ4 = 9427 if PLZ4 == 9499 & gdename == "zelg"
replace PLZ4 = 8269 if PLZ4 == 8557 & gdename == "fruthwilen"
replace PLZ4 = 8269 if PLZ4 == 8557 & gdename == "fruthwllen"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "hattenhausen/ lipperswii"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "hattenhausen/ lipperswil"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "hattenhausen/ upperswil"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "hattenhausen/lipperswil"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "hattenhausen/lipperswll"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "hefenhausen/ lipperswil"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "hefenhausen/ lipperswll"
replace PLZ4 = 8556 if PLZ4 == 8557 & gdename == "ilihart"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "lamperswil"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "lipperewil"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "lipperswi)"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "lipperswil"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "lipperswll"
replace PLZ4 = 8556 if PLZ4 == 8557 & gdename == "lllhart"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "llpperswii"
replace PLZ4 = 8564 if PLZ4 == 8557 & gdename == "llpperswll"
replace PLZ4 = 8558 if PLZ4 == 8557 & gdename == "raperswilen tg"
replace PLZ4 = 8558 if PLZ4 == 8557 & gdename == "raperswllen"
replace PLZ4 = 8558 if PLZ4 == 8557 & gdename == "reperswilen"
replace PLZ4 = 6963 if PLZ4 == 8944 & gdename == "cureglla"
replace PLZ4 = 6340 if PLZ4 == 8944 & gdename == "s hlbrugg dorf"
replace PLZ4 = 6340 if PLZ4 == 8944 & gdename == "sihibrugg"
replace PLZ4 = 6340 if PLZ4 == 8944 & gdename == "sihibrugg dorf"
replace PLZ4 = 6340 if PLZ4 == 8944 & gdename == "sihibrugg dort"
replace PLZ4 = 6340 if PLZ4 == 8944 & gdename == "sihlbrugg"
replace PLZ4 = 6340 if PLZ4 == 8944 & gdename == "sihlbrugg dorf"
replace PLZ4 = 6340 if PLZ4 == 8944 & gdename == "sihlbrugg dorf/baar"
replace PLZ4 = 6340 if PLZ4 == 8944 & gdename == "sihlbrugg dort"
replace PLZ4 = 6340 if PLZ4 == 8944 & gdename == "slhlbrugg"
replace PLZ4 = 6340 if PLZ4 == 8944 & gdename == "slhlbrugg dorf"
replace PLZ4 = 1288 if PLZ4 == 1249 & gdename == "aire"
replace PLZ4 = 1285 if PLZ4 == 1249 & gdename == "athenay"
replace PLZ4 = 1285 if PLZ4 == 1249 & gdename == "athenaz"
replace PLZ4 = 1237 if PLZ4 == 1249 & gdename == "avuily"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "bossy"
replace PLZ4 = 1244 if PLZ4 == 1249 & gdename == "chevrier"
replace PLZ4 = 1244 if PLZ4 == 1249 & gdename == "chevrier choulex"
replace PLZ4 = 1735 if PLZ4 == 1249 & gdename == "chevrilles"
replace PLZ4 = 1244 if PLZ4 == 1249 & gdename == "chouiex"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "codex bossy"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "coilex"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "coilex bossy"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "coliex"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "coliex bossy"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "coljex"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "colles"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "collex"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "collex bos y"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "collex bossy ge"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "collez"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "collez bossy"
replace PLZ4 = 1283 if PLZ4 == 1249 & gdename == "daradagny"
replace PLZ4 = 1283 if PLZ4 == 1249 & gdename == "la plaine"
replace PLZ4 = 1243 if PLZ4 == 1249 & gdename == "presigne"
replace PLZ4 = 1243 if PLZ4 == 1249 & gdename == "presinge/meinier"
replace PLZ4 = 1243 if PLZ4 == 1249 & gdename == "presinoe"
replace PLZ4 = 1241 if PLZ4 == 1249 & gdename == "pupiinge"
replace PLZ4 = 1241 if PLZ4 == 1249 & gdename == "puplinges"
replace PLZ4 = 1281 if PLZ4 == 1249 & gdename == "rusaln"
replace PLZ4 = 1285 if PLZ4 == 1249 & gdename == "sezegnin"
replace PLZ4 = 1285 if PLZ4 == 1249 & gdename == "sezegnin/ge"
replace PLZ4 = 1285 if PLZ4 == 1249 & gdename == "sezegnln"
replace PLZ4 = 1286 if PLZ4 == 1249 & gdename == "sorel"
replace PLZ4 = 1239 if PLZ4 == 1249 & gdename == "vireloup collex ge"
replace PLZ4 = 1321 if PLZ4 == 1349 & gdename == "amex"
replace PLZ4 = 1321 if PLZ4 == 1349 & gdename == "amex sur orbe"
replace PLZ4 = 1321 if PLZ4 == 1349 & gdename == "annex s/orbe"
replace PLZ4 = 1321 if PLZ4 == 1349 & gdename == "arnex"
replace PLZ4 = 1321 if PLZ4 == 1349 & gdename == "arnex s/orbe"
replace PLZ4 = 1329 if PLZ4 == 1349 & gdename == "bretonnleres"
replace PLZ4 = 1148 if PLZ4 == 1349 & gdename == "coudre l isle"
replace PLZ4 = 1322 if PLZ4 == 1349 & gdename == "croy/romainmdtier"
replace PLZ4 = 1322 if PLZ4 == 1349 & gdename == "croy/romainmã³tier"
replace PLZ4 = 1306 if PLZ4 == 1349 & gdename == "daiilens"
replace PLZ4 = 1306 if PLZ4 == 1349 & gdename == "dailiens"
replace PLZ4 = 1306 if PLZ4 == 1349 & gdename == "dalliens"
replace PLZ4 = 1306 if PLZ4 == 1349 & gdename == "dallions"
replace PLZ4 = 1306 if PLZ4 == 1349 & gdename == "dalllens"
replace PLZ4 = 1312 if PLZ4 == 1349 & gdename == "eciepens"
replace PLZ4 = 1326 if PLZ4 == 1349 & gdename == "jurions"
replace PLZ4 = 1326 if PLZ4 == 1349 & gdename == "jurlens"
replace PLZ4 = 1148 if PLZ4 == 1349 & gdename == "la coudre"
replace PLZ4 = 1148 if PLZ4 == 1349 & gdename == "la coudre (l jsle)"
replace PLZ4 = 1148 if PLZ4 == 1349 & gdename == "la coudre l isle"
replace PLZ4 = 1148 if PLZ4 == 1349 & gdename == "la coudre l jsle"
replace PLZ4 = 1148 if PLZ4 == 1349 & gdename == "la coudre/l isie"
replace PLZ4 = 1148 if PLZ4 == 1349 & gdename == "moiry vd"
replace PLZ4 = 1148 if PLZ4 == 1349 & gdename == "molry"
replace PLZ4 = 1303 if PLZ4 == 1349 & gdename == "penthas vd"
replace PLZ4 = 1318 if PLZ4 == 1349 & gdename == "pomgaples"
replace PLZ4 = 1323 if PLZ4 == 1349 & gdename == "romainmbtier"
replace PLZ4 = 1323 if PLZ4 == 1349 & gdename == "romainmoetier"
replace PLZ4 = 1323 if PLZ4 == 1349 & gdename == "romainmoetler"
replace PLZ4 = 1323 if PLZ4 == 1349 & gdename == "romainmotier"
replace PLZ4 = 1323 if PLZ4 == 1349 & gdename == "romainmotler"
replace PLZ4 = 1323 if PLZ4 == 1349 & gdename == "romainmotter"
replace PLZ4 = 1323 if PLZ4 == 1349 & gdename == "romainmã³tier"
replace PLZ4 = 1323 if PLZ4 == 1349 & gdename == "romainmã³tler"
replace PLZ4 = 1323 if PLZ4 == 1349 & gdename == "romalnmotler"
replace PLZ4 = 1323 if PLZ4 == 1349 & gdename == "romalnmã³tier"
replace PLZ4 = 1148 if PLZ4 == 1349 & gdename == "st. denis"
replace PLZ4 = 1325 if PLZ4 == 1349 & gdename == "vauiion"
replace PLZ4 = 1325 if PLZ4 == 1349 & gdename == "vaullon"
replace PLZ4 = 1325 if PLZ4 == 1349 & gdename == "vaution"
replace PLZ4 = 1302 if PLZ4 == 1349 & gdename == "vuffoens la ville"
replace PLZ4 = 4657 if PLZ4 == 4857 & gdename == "dulllken"
replace PLZ4 = 4853 if PLZ4 == 4857 & gdename == "murgenthal/riken"
replace PLZ4 = 4853 if PLZ4 == 4857 & gdename == "riken"
replace PLZ4 = 4853 if PLZ4 == 4857 & gdename == "riken (ag)"
replace PLZ4 = 4853 if PLZ4 == 4857 & gdename == "riken ag"
replace PLZ4 = 4853 if PLZ4 == 4857 & gdename == "riken bei murgenthal"
replace PLZ4 = 4853 if PLZ4 == 4857 & gdename == "riken murgenthal"
replace PLZ4 = 4853 if PLZ4 == 4857 & gdename == "riken/ murgenthal"
replace PLZ4 = 4853 if PLZ4 == 4857 & gdename == "riken/murgenthal"
replace PLZ4 = 4853 if PLZ4 == 4857 & gdename == "rlken"
replace PLZ4 = 7504 if PLZ4 == 7749 & gdename == "bernina suof"
replace PLZ4 = 7748 if PLZ4 == 7749 & gdename == "campascio"
replace PLZ4 = 7748 if PLZ4 == 7749 & gdename == "campascio brusio"
replace PLZ4 = 7748 if PLZ4 == 7749 & gdename == "campasclo"
replace PLZ4 = 7748 if PLZ4 == 7749 & gdename == "campasclo brusio"
replace PLZ4 = 7745 if PLZ4 == 7749 & gdename == "li curt"
replace PLZ4 = 7710 if PLZ4 == 7749 & gdename == "ospizio bernina"
replace PLZ4 = 7745 if PLZ4 == 7749 & gdename == "prada"
replace PLZ4 = 7745 if PLZ4 == 7749 & gdename == "prada annunziata"
replace PLZ4 = 7741 if PLZ4 == 7749 & gdename == "s. carlo"
replace PLZ4 = 7741 if PLZ4 == 7749 & gdename == "s. carlo (poschiavo)"
replace PLZ4 = 7745 if PLZ4 == 7749 & gdename == "s.antonio"
replace PLZ4 = 7745 if PLZ4 == 7749 & gdename == "s.antonio/poschiavo"
replace PLZ4 = 7741 if PLZ4 == 7749 & gdename == "s.carlo"
replace PLZ4 = 7741 if PLZ4 == 7749 & gdename == "s.carlo (poschiavo)"
replace PLZ4 = 1792 if PLZ4 == 1781 & gdename == "grossguschelmuth"
replace PLZ4 = 1789 if PLZ4 == 1781 & gdename == "lugnorre"
replace PLZ4 = 1789 if PLZ4 == 1781 & gdename == "lugnorre fr"
replace PLZ4 = 1789 if PLZ4 == 1781 & gdename == "lugnorre haut vully"
replace PLZ4 = 1789 if PLZ4 == 1781 & gdename == "lugnorre vuily le haut"
replace PLZ4 = 1789 if PLZ4 == 1781 & gdename == "lugnorre vully le haut"
replace PLZ4 = 1789 if PLZ4 == 1781 & gdename == "lugnorre/haut vuliy"
replace PLZ4 = 1789 if PLZ4 == 1781 & gdename == "lugnorre/haut vully"
replace PLZ4 = 1787 if PLZ4 == 1781 & gdename == "lugnorre/motier"
replace PLZ4 = 1797 if PLZ4 == 1781 & gdename == "muenchenwiier"
replace PLZ4 = 1797 if PLZ4 == 1781 & gdename == "muenchenwller"
replace PLZ4 = 1797 if PLZ4 == 1781 & gdename == "muenchwiler"
replace PLZ4 = 1786 if PLZ4 == 1781 & gdename == "nieder wistenlach"
replace PLZ4 = 1786 if PLZ4 == 1781 & gdename == "nieder wlstenlach"
replace PLZ4 = 1786 if PLZ4 == 1781 & gdename == "nleder wistenlach"
replace PLZ4 = 1788 if PLZ4 == 1781 & gdename == "praz"
replace PLZ4 = 1788 if PLZ4 == 1781 & gdename == "praz (vuily)"
replace PLZ4 = 1788 if PLZ4 == 1781 & gdename == "praz (vully)"
replace PLZ4 = 1788 if PLZ4 == 1781 & gdename == "praz vully"
replace PLZ4 = 1794 if PLZ4 == 1781 & gdename == "saivenach"
replace PLZ4 = 4102 if PLZ4 == 4113 & gdename == "biningen"
replace PLZ4 = 4114 if PLZ4 == 4113 & gdename == "fiueh"
replace PLZ4 = 4114 if PLZ4 == 4113 & gdename == "fiueh so"
replace PLZ4 = 4114 if PLZ4 == 4113 & gdename == "flueh"
replace PLZ4 = 4114 if PLZ4 == 4113 & gdename == "flueh so"
replace PLZ4 = 4114 if PLZ4 == 4113 & gdename == "flueh/so"
replace PLZ4 = 4114 if PLZ4 == 4113 & gdename == "fluh"
replace PLZ4 = 4114 if PLZ4 == 4113 & gdename == "fluh so"
replace PLZ4 = 7187 if PLZ4 == 7181 & gdename == "camischolas"
replace PLZ4 = 7184 if PLZ4 == 7181 & gdename == "curaglia"
replace PLZ4 = 7184 if PLZ4 == 7181 & gdename == "curaglla"
replace PLZ4 = 7189 if PLZ4 == 7181 & gdename == "rueras"
replace PLZ4 = 7189 if PLZ4 == 7181 & gdename == "rueras tujetsch"
replace PLZ4 = 7189 if PLZ4 == 7181 & gdename == "rueras/tujetsch"
replace PLZ4 = 7189 if PLZ4 == 7181 & gdename == "ruerasftujetsch"
replace PLZ4 = 7017 if PLZ4 == 7181 & gdename == "segnes"
replace PLZ4 = 7189 if PLZ4 == 7181 & gdename == "selva"
replace PLZ4 = 7189 if PLZ4 == 7181 & gdename == "selva (tavetsch)"
replace PLZ4 = 7536 if PLZ4 == 7181 & gdename == "st. maria"
replace PLZ4 = 7146 if PLZ4 == 7131 & gdename == "(gels vattiz"
replace PLZ4 = 7133 if PLZ4 == 7131 & gdename == "affeier"
replace PLZ4 = 7133 if PLZ4 == 7131 & gdename == "affeler"
replace PLZ4 = 7153 if PLZ4 == 7131 & gdename == "faiera"
replace PLZ4 = 7153 if PLZ4 == 7131 & gdename == "fellers"
replace PLZ4 = 7137 if PLZ4 == 7131 & gdename == "fiond"
replace PLZ4 = 7146 if PLZ4 == 7131 & gdename == "igels"
replace PLZ4 = 7146 if PLZ4 == 7131 & gdename == "igels vattiz"
replace PLZ4 = 7031 if PLZ4 == 7131 & gdename == "laax b. ilanz"
replace PLZ4 = 7031 if PLZ4 == 7131 & gdename == "laax b. jlanz"
replace PLZ4 = 7155 if PLZ4 == 7131 & gdename == "ladlr"
replace PLZ4 = 7148 if PLZ4 == 7131 & gdename == "lumbreln"
replace PLZ4 = 6476 if PLZ4 == 7131 & gdename == "luren"
replace PLZ4 = 7134 if PLZ4 == 7131 & gdename == "obersaken"
replace PLZ4 = 7134 if PLZ4 == 7131 & gdename == "obersaxen tusen"
replace PLZ4 = 7134 if PLZ4 == 7131 & gdename == "obesaxen"
replace PLZ4 = 7152 if PLZ4 == 7131 & gdename == "sagoon"
replace PLZ4 = 7151 if PLZ4 == 7131 & gdename == "schieuis"
replace PLZ4 = 7151 if PLZ4 == 7131 & gdename == "schiuein"
replace PLZ4 = 7151 if PLZ4 == 7131 & gdename == "schleuis"
replace PLZ4 = 7151 if PLZ4 == 7131 & gdename == "schleuls"
replace PLZ4 = 7151 if PLZ4 == 7131 & gdename == "schleunis"
replace PLZ4 = 7151 if PLZ4 == 7131 & gdename == "schlueln"
replace PLZ4 = 7127 if PLZ4 == 7131 & gdename == "sergein"
replace PLZ4 = 7127 if PLZ4 == 7131 & gdename == "sergeln"
replace PLZ4 = 7127 if PLZ4 == 7131 & gdename == "sevgeln"
replace PLZ4 = 7152 if PLZ4 == 7131 & gdename == "sgogn"
replace PLZ4 = 7138 if PLZ4 == 7131 & gdename == "surcuoim"
replace PLZ4 = 7114 if PLZ4 == 7131 & gdename == "uors"
replace PLZ4 = 7138 if PLZ4 == 7131 & gdename == "vaiata obersaxen"
replace PLZ4 = 7138 if PLZ4 == 7131 & gdename == "valata obersaxen"
replace PLZ4 = 7138 if PLZ4 == 7131 & gdename == "valataobersaxen"
replace PLZ4 = 7146 if PLZ4 == 7131 & gdename == "vattiz"
replace PLZ4 = 7147 if PLZ4 == 7131 & gdename == "vigens"
replace PLZ4 = 7144 if PLZ4 == 7131 & gdename == "villa"
replace PLZ4 = 7144 if PLZ4 == 7131 & gdename == "villa gr"
replace PLZ4 = 6713 if PLZ4 == 6712 & gdename == "maivaglia chiesa"
replace PLZ4 = 6713 if PLZ4 == 6712 & gdename == "malvaglia chiesa"
replace PLZ4 = 6713 if PLZ4 == 6712 & gdename == "malvaglia chiesu"
replace PLZ4 = 6713 if PLZ4 == 6712 & gdename == "malvaglia chlesa"
replace PLZ4 = 6713 if PLZ4 == 6712 & gdename == "malvaglla chiesa"
replace PLZ4 = 6710 if PLZ4 == 6712 & gdename == "pontirone"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "beilerive"
replace PLZ4 = 1589 if PLZ4 == 1581 & gdename == "chaberey"
replace PLZ4 = 1589 if PLZ4 == 1581 & gdename == "chambery"
replace PLZ4 = 1584 if PLZ4 == 1581 & gdename == "les friques"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "saiavaux"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "salavaux"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "salavaux bellerive"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "salavaux vd"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "vailamand dessous"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "valiamand"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "valiamand dessous"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "valiamand dessus"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "vallamand dessous"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "vallamand dessu3"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "vallamand dessus"
replace PLZ4 = 1585 if PLZ4 == 1581 & gdename == "vallamand.dessous"
replace PLZ4 = 4127 if PLZ4 == 4122 & gdename == "blirsfelden"
replace PLZ4 = 4114 if PLZ4 == 4122 & gdename == "flueh"
replace PLZ4 = 4123 if PLZ4 == 4122 & gdename == "neu allschwil"
replace PLZ4 = 4123 if PLZ4 == 4122 & gdename == "neuailschwil"
replace PLZ4 = 4123 if PLZ4 == 4122 & gdename == "neualischwil"
replace PLZ4 = 4123 if PLZ4 == 4122 & gdename == "neualischwll"
replace PLZ4 = 4123 if PLZ4 == 4122 & gdename == "neuallschwii"
replace PLZ4 = 4123 if PLZ4 == 4122 & gdename == "neuallschwil"
replace PLZ4 = 3536 if PLZ4 == 3549 & gdename == "aeschau"
replace PLZ4 = 3556 if PLZ4 == 3549 & gdename == "fankhaus (trub)"
replace PLZ4 = 3556 if PLZ4 == 3549 & gdename == "fankhaus trub"
replace PLZ4 = 3556 if PLZ4 == 3549 & gdename == "fankhaus(trub)"
replace PLZ4 = 3550 if PLZ4 == 3549 & gdename == "gohl"
replace PLZ4 = 3510 if PLZ4 == 3549 & gdename == "gysenstein"
replace PLZ4 = 3510 if PLZ4 == 3549 & gdename == "gysenstein konolfingen"
replace PLZ4 = 3510 if PLZ4 == 3549 & gdename == "gysensteln"
replace PLZ4 = 3531 if PLZ4 == 3549 & gdename == "oberthai"
replace PLZ4 = 3111 if PLZ4 == 3549 & gdename == "taegertschl"
replace PLZ4 = 3111 if PLZ4 == 3549 & gdename == "tagertschi"
replace PLZ4 = 3083 if PLZ4 == 3549 & gdename == "triimsteln"
replace PLZ4 = 3083 if PLZ4 == 3549 & gdename == "trimatein"
replace PLZ4 = 3083 if PLZ4 == 3549 & gdename == "trimstcin"
replace PLZ4 = 3083 if PLZ4 == 3549 & gdename == "trimstein"
replace PLZ4 = 3083 if PLZ4 == 3549 & gdename == "trimsteln"
replace PLZ4 = 3083 if PLZ4 == 3549 & gdename == "trlmstein"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "wolf halden zeig"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "wolfhalden zeig"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "wolfhalden zelg"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "wollhalden zeig"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zeig"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zeig (wolf halden)"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zeig (wolfhaiden)"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zeig (wolfhalden)"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zeig (wollhalden)"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zeig wolfhalden"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zeig/woifhaiden"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zelg"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zelg (woifhalden)"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zelg (wolf halden)"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zelg (wolfhaiden)"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zelg (wolfhalden)"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zelg (wollhalden)"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zelg wolfhaiden"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zelg wolfhalden"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zelg/wolfhalden"
replace PLZ4 = 9427 if PLZ4 == 9429 & gdename == "zolg wolfhalden"
replace PLZ4 = 8880 if PLZ4 == 8891 & gdename == "berschis"
replace PLZ4 = 8880 if PLZ4 == 8891 & gdename == "berschis (walenstadt)"
replace PLZ4 = 8880 if PLZ4 == 8891 & gdename == "berschis walenstadt"
replace PLZ4 = 8880 if PLZ4 == 8891 & gdename == "berschis/walenstadt"
replace PLZ4 = 8880 if PLZ4 == 8891 & gdename == "berschls"
replace PLZ4 = 8890 if PLZ4 == 8891 & gdename == "grossberg fiums"
replace PLZ4 = 8890 if PLZ4 == 8891 & gdename == "grossberg flums"
replace PLZ4 = 8890 if PLZ4 == 8891 & gdename == "tannenbodenalp"
replace PLZ4 = 8890 if PLZ4 == 8891 & gdename == "tannenbodenalp (flums)"
replace PLZ4 = 8890 if PLZ4 == 8891 & gdename == "tannenheim"
replace PLZ4 = 8890 if PLZ4 == 8891 & gdename == "tannenheim fiums"
replace PLZ4 = 8890 if PLZ4 == 8891 & gdename == "tannenheim flume"
replace PLZ4 = 8890 if PLZ4 == 8891 & gdename == "tannenheim flums"
replace PLZ4 = 8890 if PLZ4 == 8891 & gdename == "tannenheim plums"
replace PLZ4 = 8890 if PLZ4 == 8891 & gdename == "tannenhelm flums"
replace PLZ4 = 1610 if PLZ4 == 1599 & gdename == "chaetillens"
replace PLZ4 = 1610 if PLZ4 == 1599 & gdename == "chatiilens"
replace PLZ4 = 1610 if PLZ4 == 1599 & gdename == "chatiliens"
replace PLZ4 = 1610 if PLZ4 == 1599 & gdename == "chatlllens"
replace PLZ4 = 1024 if PLZ4 == 1599 & gdename == "ecubiens"
replace PLZ4 = 1525 if PLZ4 == 1599 & gdename == "les treize cantons"
replace PLZ4 = 1524 if PLZ4 == 1599 & gdename == "mamand"
replace PLZ4 = 1524 if PLZ4 == 1599 & gdename == "marmand"
replace PLZ4 = 1524 if PLZ4 == 1599 & gdename == "marriand"
replace PLZ4 = 1607 if PLZ4 == 1599 & gdename == "palesieux vd"
replace PLZ4 = 1607 if PLZ4 == 1599 & gdename == "palezieux vd"
replace PLZ4 = 1607 if PLZ4 == 1599 & gdename == "palezieux village"
replace PLZ4 = 1607 if PLZ4 == 1599 & gdename == "palezleux vd"
replace PLZ4 = 1607 if PLZ4 == 1599 & gdename == "palezleux village"
replace PLZ4 = 1525 if PLZ4 == 1599 & gdename == "seigneuz"
replace PLZ4 = 1525 if PLZ4 == 1599 & gdename == "soigneux"
replace PLZ4 = 1607 if PLZ4 == 1599 & gdename == "tavernes"
replace PLZ4 = 1525 if PLZ4 == 1599 & gdename == "treize cantcns"
replace PLZ4 = 1525 if PLZ4 == 1599 & gdename == "treize canton"
replace PLZ4 = 1525 if PLZ4 == 1599 & gdename == "treize cantons"
replace PLZ4 = 1610 if PLZ4 == 1599 & gdename == "vuibroye chaetillens"
replace PLZ4 = 1610 if PLZ4 == 1599 & gdename == "vuibroye chatiilens"
replace PLZ4 = 1610 if PLZ4 == 1599 & gdename == "vuibroye chatiliens"
replace PLZ4 = 1610 if PLZ4 == 1599 & gdename == "vuibroye chatillens"
replace PLZ4 = 1610 if PLZ4 == 1599 & gdename == "vulbroye"
replace PLZ4 = 1610 if PLZ4 == 1599 & gdename == "vulbroye chatlllens"
replace PLZ4 = 1856 if PLZ4 == 1861 & gdename == "corbeyraer"
replace PLZ4 = 1856 if PLZ4 == 1861 & gdename == "corbeyrler"
replace PLZ4 = 1866 if PLZ4 == 1861 & gdename == "forclaz"
replace PLZ4 = 1884 if PLZ4 == 1861 & gdename == "huemoz"
replace PLZ4 = 1884 if PLZ4 == 1861 & gdename == "huemoz sur ollon"
replace PLZ4 = 1884 if PLZ4 == 1861 & gdename == "huemoz sur otton"
replace PLZ4 = 1862 if PLZ4 == 1861 & gdename == "la comballaz"
replace PLZ4 = 1866 if PLZ4 == 1861 & gdename == "la forclaz"
replace PLZ4 = 1866 if PLZ4 == 1861 & gdename == "la forclaz/evolene"
replace PLZ4 = 1862 if PLZ4 == 1861 & gdename == "les mosses"
replace PLZ4 = 1862 if PLZ4 == 1861 & gdename == "les mosses (ormont dessous)"
replace PLZ4 = 1862 if PLZ4 == 1861 & gdename == "mosses/ ormont dessous"
replace PLZ4 = 1862 if PLZ4 == 1861 & gdename == "mosses/ormont dessous"
replace PLZ4 = 1867 if PLZ4 == 1861 & gdename == "panex"
replace PLZ4 = 1867 if PLZ4 == 1861 & gdename == "panex ollon"
replace PLZ4 = 1867 if PLZ4 == 1861 & gdename == "panex sur ollon"
replace PLZ4 = 1867 if PLZ4 == 1861 & gdename == "panex sur olson"
replace PLZ4 = 1867 if PLZ4 == 1861 & gdename == "panez ollon"
replace PLZ4 = 6874 if PLZ4 == 6851 & gdename == "cortegiia"
replace PLZ4 = 6874 if PLZ4 == 6851 & gdename == "corteglia"
replace PLZ4 = 6874 if PLZ4 == 6851 & gdename == "corteglia die castel san pietro"
replace PLZ4 = 6874 if PLZ4 == 6851 & gdename == "corteglia/ castel s. pietro"
replace PLZ4 = 6874 if PLZ4 == 6851 & gdename == "corteglia/castel s. pietro"
replace PLZ4 = 6874 if PLZ4 == 6851 & gdename == "corteglla"
replace PLZ4 = 6850 if PLZ4 == 6851 & gdename == "mendrisio sortite"
replace PLZ4 = 8725 if PLZ4 == 8731 & gdename == "emetschwil"
replace PLZ4 = 8725 if PLZ4 == 8731 & gdename == "ernetschwll"
replace PLZ4 = 8725 if PLZ4 == 8731 & gdename == "gebertangen"
replace PLZ4 = 8725 if PLZ4 == 8731 & gdename == "gebertingen"
replace PLZ4 = 8737 if PLZ4 == 8731 & gdename == "gomiswald"
replace PLZ4 = 8853 if PLZ4 == 8731 & gdename == "lachen sz"
replace PLZ4 = 8726 if PLZ4 == 8731 & gdename == "ricken"
replace PLZ4 = 8726 if PLZ4 == 8731 & gdename == "ricken (sg)"
replace PLZ4 = 8739 if PLZ4 == 8731 & gdename == "rieden sg"
replace PLZ4 = 8739 if PLZ4 == 8731 & gdename == "rleden"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetiiberg"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetliberg"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetliburg"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetliburg b.gommiswald"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetliburg bei gommiswald"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetliburg bel gommiswald"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetliburg gommiswaid"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetliburg gommiswald"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetllberg"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetllburg"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetllburg gommiswald"
replace PLZ4 = 8738 if PLZ4 == 8731 & gdename == "uetlã¬burg"
replace PLZ4 = 1070 if PLZ4 == 1604 & gdename == "pueldoux gare"
replace PLZ4 = 1070 if PLZ4 == 1604 & gdename == "puidoux gare"
replace PLZ4 = 1070 if PLZ4 == 1604 & gdename == "puidoux village"
replace PLZ4 = 1070 if PLZ4 == 1604 & gdename == "puldoux"
replace PLZ4 = 1070 if PLZ4 == 1604 & gdename == "puldoux gare"
replace PLZ4 = 1070 if PLZ4 == 1604 & gdename == "puldoux village"
replace PLZ4 = 9413 if PLZ4 == 9432 & gdename == "bueriswilen"
replace PLZ4 = 9413 if PLZ4 == 9432 & gdename == "bueriswilen oberegg"
replace PLZ4 = 9413 if PLZ4 == 9432 & gdename == "bueriswilen/oberegg"
replace PLZ4 = 9413 if PLZ4 == 9432 & gdename == "bueriswilenioberegg"
replace PLZ4 = 9413 if PLZ4 == 9432 & gdename == "bueriswllen oberegg"
replace PLZ4 = 9413 if PLZ4 == 9432 & gdename == "buerlswilen/oberegg"
replace PLZ4 = 9428 if PLZ4 == 9432 & gdename == "platz"
replace PLZ4 = 9428 if PLZ4 == 9432 & gdename == "platz (ar)"
replace PLZ4 = 9428 if PLZ4 == 9432 & gdename == "platz (waizenhausen)"
replace PLZ4 = 9428 if PLZ4 == 9432 & gdename == "platz (walzenhausen)"
replace PLZ4 = 9428 if PLZ4 == 9432 & gdename == "platz waizenhausen"
replace PLZ4 = 9428 if PLZ4 == 9432 & gdename == "platz walzenhausen"
replace PLZ4 = 9428 if PLZ4 == 9432 & gdename == "waizenhausen platz"
replace PLZ4 = 9428 if PLZ4 == 9432 & gdename == "walzenhausen piatz"
replace PLZ4 = 9428 if PLZ4 == 9432 & gdename == "walzenhausen platz"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wailiswil bel niederbipp"
replace PLZ4 = 3377 if PLZ4 == 4705 & gdename == "waliiswil"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "waliiswil bei niederbipp"
replace PLZ4 = 3377 if PLZ4 == 4705 & gdename == "walliswal"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "walliswal bei niederbipp"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "walliswal bel niederbipp"
replace PLZ4 = 3377 if PLZ4 == 4705 & gdename == "walliswii"
replace PLZ4 = 3377 if PLZ4 == 4705 & gdename == "walliswil"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "walliswil b. niederbipp"
replace PLZ4 = 3377 if PLZ4 == 4705 & gdename == "walliswil b. wangen"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "walliswil bel niederbipp"
replace PLZ4 = 3377 if PLZ4 == 4705 & gdename == "walliswll"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "walltswil bei niederbipp"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a aare"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a. a"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a. a."
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a. aare"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a. albas"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a. albis"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a. d"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a. d. aare"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a.a"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a.aare"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a.d. aare"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a.d./a."
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a.d.a"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen a.d.aare"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangen en der aare"
replace PLZ4 = 3380 if PLZ4 == 4705 & gdename == "wangena.a"
replace PLZ4 = 5524 if PLZ4 == 5523 & gdename == "nesseinbach"
replace PLZ4 = 5524 if PLZ4 == 5523 & gdename == "nesseinbach ag"
replace PLZ4 = 5524 if PLZ4 == 5523 & gdename == "nesselnbach"
replace PLZ4 = 5524 if PLZ4 == 5523 & gdename == "nesselnbach /niederwil ag"
replace PLZ4 = 5524 if PLZ4 == 5523 & gdename == "nesselnbach ag"
replace PLZ4 = 5524 if PLZ4 == 5523 & gdename == "nesselnbach/niederwal ag"
replace PLZ4 = 5524 if PLZ4 == 5523 & gdename == "nesselnbach/niederwil"
replace PLZ4 = 5524 if PLZ4 == 5523 & gdename == "nesselnbach/niederwil ag"
replace PLZ4 = 5524 if PLZ4 == 5523 & gdename == "nesselnbach/niederwll"
replace PLZ4 = 5524 if PLZ4 == 5523 & gdename == "niederwal ag"
replace PLZ4 = 1219 if PLZ4 == 1210 & gdename == "chaetelaine"
replace PLZ4 = 1219 if PLZ4 == 1210 & gdename == "chaetelaine vernier"
replace PLZ4 = 1219 if PLZ4 == 1210 & gdename == "chatelaine"
replace PLZ4 = 1219 if PLZ4 == 1210 & gdename == "chatelaine vernier"
replace PLZ4 = 1218 if PLZ4 == 1210 & gdename == "grand saconnex"
replace PLZ4 = 1219 if PLZ4 == 1210 & gdename == "la chatelaine"
replace PLZ4 = 2053 if PLZ4 == 1531 & gdename == "cerniez"
replace PLZ4 = 1553 if PLZ4 == 1531 & gdename == "chaetonnaye"
replace PLZ4 = 1553 if PLZ4 == 1531 & gdename == "chatennaye"
replace PLZ4 = 1545 if PLZ4 == 1531 & gdename == "chevraux"
replace PLZ4 = 1536 if PLZ4 == 1531 & gdename == "combremont ie petit"
replace PLZ4 = 1536 if PLZ4 == 1531 & gdename == "combremont le petlt"
replace PLZ4 = 1682 if PLZ4 == 1531 & gdename == "domplerre"
replace PLZ4 = 1532 if PLZ4 == 1531 & gdename == "etigny"
replace PLZ4 = 1532 if PLZ4 == 1531 & gdename == "fetlgny"
replace PLZ4 = 1544 if PLZ4 == 1531 & gdename == "gietterens"
replace PLZ4 = 1544 if PLZ4 == 1531 & gdename == "gletternens"
replace PLZ4 = 1533 if PLZ4 == 1531 & gdename == "memeres"
replace PLZ4 = 1533 if PLZ4 == 1531 & gdename == "menleres"
replace PLZ4 = 1542 if PLZ4 == 1531 & gdename == "rueyres le pres"
replace PLZ4 = 1534 if PLZ4 == 1531 & gdename == "sasse)"
replace PLZ4 = 1534 if PLZ4 == 1531 & gdename == "sassei"
replace PLZ4 = 1534 if PLZ4 == 1531 & gdename == "sassel vd"
replace PLZ4 = 1534 if PLZ4 == 1531 & gdename == "sasses"
replace PLZ4 = 1534 if PLZ4 == 1531 & gdename == "sessel"
replace PLZ4 = 1534 if PLZ4 == 1531 & gdename == "sessel vd"
replace PLZ4 = 1551 if PLZ4 == 1531 & gdename == "vers .chez perrin"
replace PLZ4 = 1551 if PLZ4 == 1531 & gdename == "vers chez perrin"
replace PLZ4 = 1551 if PLZ4 == 1531 & gdename == "vers chez perrin/payeme"
replace PLZ4 = 1551 if PLZ4 == 1531 & gdename == "vers chez perrin/payerne"
replace PLZ4 = 1551 if PLZ4 == 1531 & gdename == "verschez payern"
replace PLZ4 = 1551 if PLZ4 == 1531 & gdename == "verschez perrin/payeme"
replace PLZ4 = 1551 if PLZ4 == 1531 & gdename == "verschez perrin/payerne"
replace PLZ4 = 1555 if PLZ4 == 1531 & gdename == "viliarzel"
replace PLZ4 = 1555 if PLZ4 == 1531 & gdename == "villarzell"
replace PLZ4 = 1555 if PLZ4 == 1531 & gdename == "vlllarzel"
replace PLZ4 = 1555 if PLZ4 == 1531 & gdename == "vlllarzell"
replace PLZ4 = 7203 if PLZ4 == 7209 & gdename == "trimmas dorf"
replace PLZ4 = 7203 if PLZ4 == 7209 & gdename == "trimmas dort"
replace PLZ4 = 7203 if PLZ4 == 7209 & gdename == "trimmfis dorf"
replace PLZ4 = 7203 if PLZ4 == 7209 & gdename == "trimmis dorf"
replace PLZ4 = 7203 if PLZ4 == 7209 & gdename == "trimmis dorl"
replace PLZ4 = 7203 if PLZ4 == 7209 & gdename == "trimmis dort"
replace PLZ4 = 7203 if PLZ4 == 7209 & gdename == "trlmmis dorf"
replace PLZ4 = 6436 if PLZ4 == 6437 & gdename == "hinterthai"
replace PLZ4 = 6436 if PLZ4 == 6437 & gdename == "hinterthal"
replace PLZ4 = 6436 if PLZ4 == 6437 & gdename == "hinterthal /muotathal"
replace PLZ4 = 6436 if PLZ4 == 6437 & gdename == "hinterthal muotathai"
replace PLZ4 = 6436 if PLZ4 == 6437 & gdename == "hinterthal muotathal"
replace PLZ4 = 6436 if PLZ4 == 6437 & gdename == "hinterthal/muotathal"
replace PLZ4 = 6436 if PLZ4 == 6437 & gdename == "hlnterthal"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "oetlikon"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "warenlos"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "werenlos"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "wueenlos"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "wuer"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "wuerenios"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "wuerenloa"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "wuerenlos ag"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "wuerenlos/ag"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "wuerentas"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "wuerentos"
replace PLZ4 = 5436 if PLZ4 == 8116 & gdename == "wurenlos"
replace PLZ4 = 1096 if PLZ4 == 1603 & gdename == "araen"
replace PLZ4 = 1096 if PLZ4 == 1603 & gdename == "aran"
replace PLZ4 = 1091 if PLZ4 == 1603 & gdename == "aran grandvaux"
replace PLZ4 = 1096 if PLZ4 == 1603 & gdename == "aran sur villette"
replace PLZ4 = 1096 if PLZ4 == 1603 & gdename == "aran villette"
replace PLZ4 = 1096 if PLZ4 == 1603 & gdename == "aren villette"
replace PLZ4 = 1096 if PLZ4 == 1603 & gdename == "arlin"
replace PLZ4 = 1096 if PLZ4 == 1603 & gdename == "aron villette"
replace PLZ4 = 1091 if PLZ4 == 1603 & gdename == "chenaux"
replace PLZ4 = 1091 if PLZ4 == 1603 & gdename == "chenaux grandvaux"
replace PLZ4 = 1091 if PLZ4 == 1603 & gdename == "granduaux"
replace PLZ4 = 1091 if PLZ4 == 1603 & gdename == "grandvaux vd"
replace PLZ4 = 1091 if PLZ4 == 1603 & gdename == "grandveau"
replace PLZ4 = 1091 if PLZ4 == 1603 & gdename == "grardvaux"
replace PLZ4 = 1052 if PLZ4 == 1603 & gdename == "mont sur lausanne"
replace PLZ4 = 8514 if PLZ4 == 8515 & gdename == "amilkon"
replace PLZ4 = 8514 if PLZ4 == 8515 & gdename == "amlikon"
replace PLZ4 = 8514 if PLZ4 == 8515 & gdename == "amlikon tg"
replace PLZ4 = 8514 if PLZ4 == 8515 & gdename == "amllkon"
replace PLZ4 = 8514 if PLZ4 == 8515 & gdename == "amlã¬kon"
replace PLZ4 = 9556 if PLZ4 == 8515 & gdename == "maltbach"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "muestalr"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "mustair"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "mãºstair"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "st. maria"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "st. maria 1. m."
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "st. maria i. m"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta maria"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta maria i.m"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta maria im muenstertal"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta maria im munstertal"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta. maria"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta. maria i. m"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta. maria i. m."
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta. maria im muenstertai"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta. maria im muenstertal"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta. maria im munstertal"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta. marla"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta. marla i. m"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "sta. marla im muenstertal"
replace PLZ4 = 7537 if PLZ4 == 7531 & gdename == "ste. maria im muenstertal"
replace PLZ4 = 7532 if PLZ4 == 7531 & gdename == "tschier"
replace PLZ4 = 7532 if PLZ4 == 7531 & gdename == "tschiery"
replace PLZ4 = 6280 if PLZ4 == 6282 & gdename == "urswe"
replace PLZ4 = 6280 if PLZ4 == 6282 & gdename == "urswii/hochdorf"
replace PLZ4 = 6280 if PLZ4 == 6282 & gdename == "urswil"
replace PLZ4 = 6280 if PLZ4 == 6282 & gdename == "urswil/hochdorf"
replace PLZ4 = 6280 if PLZ4 == 6282 & gdename == "urswll"
replace PLZ4 = 8370 if PLZ4 == 9572 & gdename == "busawll"
replace PLZ4 = 8370 if PLZ4 == 9572 & gdename == "busswal"
replace PLZ4 = 8370 if PLZ4 == 9572 & gdename == "busswal tg"
replace PLZ4 = 8370 if PLZ4 == 9572 & gdename == "busswal tg/sirnach"
replace PLZ4 = 8370 if PLZ4 == 9572 & gdename == "busswii tg"
replace PLZ4 = 8370 if PLZ4 == 9572 & gdename == "busswil tg"
replace PLZ4 = 8370 if PLZ4 == 9572 & gdename == "busswil tg/egnach"
replace PLZ4 = 8370 if PLZ4 == 9572 & gdename == "busswil tg/sirnach"
replace PLZ4 = 8370 if PLZ4 == 9572 & gdename == "busswll tg"
replace PLZ4 = 8370 if PLZ4 == 9572 & gdename == "buswal tg"
replace PLZ4 = 8335 if PLZ4 == 8336 & gdename == "duerstelen"
replace PLZ4 = 8335 if PLZ4 == 8336 & gdename == "duersteren"
replace PLZ4 = 8335 if PLZ4 == 8336 & gdename == "duersteten"
replace PLZ4 = 8335 if PLZ4 == 8336 & gdename == "ober hittnau"
replace PLZ4 = 8335 if PLZ4 == 8336 & gdename == "oberhittnau"
replace PLZ4 = 8335 if PLZ4 == 8336 & gdename == "oberhlttnau"
replace PLZ4 = 8330 if PLZ4 == 8336 & gdename == "pfaeffikon zh"
replace PLZ4 = 1820 if PLZ4 == 1842 & gdename == "montreux territet"
replace PLZ4 = 1820 if PLZ4 == 1842 & gdename == "territet"
replace PLZ4 = 1820 if PLZ4 == 1842 & gdename == "territtet"
replace PLZ4 = 1820 if PLZ4 == 1842 & gdename == "terrltet"
replace PLZ4 = 5046 if PLZ4 == 5047 & gdename == "walde"
replace PLZ4 = 5046 if PLZ4 == 5047 & gdename == "walde (ag)"
replace PLZ4 = 6696 if PLZ4 == 6671 & gdename == ".fusio mogno"
replace PLZ4 = 6677 if PLZ4 == 6671 & gdename == "auegno"
replace PLZ4 = 6677 if PLZ4 == 6671 & gdename == "aurlgeno"
replace PLZ4 = 6685 if PLZ4 == 6671 & gdename == "bosco"
replace PLZ4 = 6685 if PLZ4 == 6671 & gdename == "bosco gurin"
replace PLZ4 = 6696 if PLZ4 == 6671 & gdename == "brogllo"
replace PLZ4 = 6678 if PLZ4 == 6671 & gdename == "cogllo"
replace PLZ4 = 6696 if PLZ4 == 6671 & gdename == "fusio magna"
replace PLZ4 = 6696 if PLZ4 == 6671 & gdename == "fusio mogno"
replace PLZ4 = 6678 if PLZ4 == 6671 & gdename == "giumagliopus"
replace PLZ4 = 6678 if PLZ4 == 6671 & gdename == "giumagllo"
replace PLZ4 = 6678 if PLZ4 == 6671 & gdename == "glnmaglio"
replace PLZ4 = 6682 if PLZ4 == 6671 & gdename == "linesclo"
replace PLZ4 = 6695 if PLZ4 == 6671 & gdename == "peccla"
replace PLZ4 = 6674 if PLZ4 == 6671 & gdename == "riveo"
replace PLZ4 = 6674 if PLZ4 == 6671 & gdename == "riveo di someo"
replace PLZ4 = 6674 if PLZ4 == 6671 & gdename == "rlveo"
replace PLZ4 = 6677 if PLZ4 == 6671 & gdename == "ronchina aurigeno"
replace PLZ4 = 1407 if PLZ4 == 1411 & gdename == "donneioye"
replace PLZ4 = 1417 if PLZ4 == 1411 & gdename == "essertides sur yverdon"
replace PLZ4 = 1417 if PLZ4 == 1411 & gdename == "essertines yverdon"
replace PLZ4 = 1417 if PLZ4 == 1411 & gdename == "essertlnes sur yverdon"
replace PLZ4 = 1420 if PLZ4 == 1411 & gdename == "flez"
replace PLZ4 = 1422 if PLZ4 == 1411 & gdename == "fontaines/grandson"
replace PLZ4 = 1423 if PLZ4 == 1411 & gdename == "fontanezler"
replace PLZ4 = 1429 if PLZ4 == 1411 & gdename == "gier"
replace PLZ4 = 1429 if PLZ4 == 1411 & gdename == "giez vd"
replace PLZ4 = 1429 if PLZ4 == 1411 & gdename == "glez"
replace PLZ4 = 1429 if PLZ4 == 1411 & gdename == "glez vd"
replace PLZ4 = 1421 if PLZ4 == 1411 & gdename == "grandevent vd"
replace PLZ4 = 1422 if PLZ4 == 1411 & gdename == "les tuileries"
replace PLZ4 = 1422 if PLZ4 == 1411 & gdename == "les tuileries de grandson"
replace PLZ4 = 1422 if PLZ4 == 1411 & gdename == "les tuilerles de grandson"
replace PLZ4 = 1422 if PLZ4 == 1411 & gdename == "les tulleries de grandson"
replace PLZ4 = 1422 if PLZ4 == 1411 & gdename == "les tullerles/grandson"
replace PLZ4 = 1416 if PLZ4 == 1411 & gdename == "pallly"
replace PLZ4 = 1408 if PLZ4 == 1411 & gdename == "prahlns"
replace PLZ4 = 1422 if PLZ4 == 1411 & gdename == "tuileries de grandson"
replace PLZ4 = 1423 if PLZ4 == 1411 & gdename == "viliars burquin"
replace PLZ4 = 1423 if PLZ4 == 1411 & gdename == "villars bourquin"
replace PLZ4 = 1423 if PLZ4 == 1411 & gdename == "villars burquln"
replace PLZ4 = 1418 if PLZ4 == 1411 & gdename == "vuarrengei"
replace PLZ4 = 1418 if PLZ4 == 1411 & gdename == "vuarrengel"
replace PLZ4 = 1431 if PLZ4 == 1411 & gdename == "vugelies la mothe"
replace PLZ4 = 3617 if PLZ4 == 3611 & gdename == "fahmi bel thun"
replace PLZ4 = 3617 if PLZ4 == 3611 & gdename == "fahml"
replace PLZ4 = 3617 if PLZ4 == 3611 & gdename == "fahrei b. thun"
replace PLZ4 = 3617 if PLZ4 == 3611 & gdename == "fahrni b. thun"
replace PLZ4 = 3617 if PLZ4 == 3611 & gdename == "fahrni bei thun"
replace PLZ4 = 3617 if PLZ4 == 3611 & gdename == "fahrni thun"
replace PLZ4 = 3617 if PLZ4 == 3611 & gdename == "fahrnl bel thun"
replace PLZ4 = 3636 if PLZ4 == 3611 & gdename == "forst b. laengenbuehl"
replace PLZ4 = 3631 if PLZ4 == 3611 & gdename == "hoefen bei thun"
replace PLZ4 = 3631 if PLZ4 == 3611 & gdename == "hoefen bel thun"
replace PLZ4 = 3631 if PLZ4 == 3611 & gdename == "hoefen/thun"
replace PLZ4 = 3632 if PLZ4 == 3611 & gdename == "nlederstocken"
replace PLZ4 = 3638 if PLZ4 == 3611 & gdename == "pohiern"
replace PLZ4 = 3618 if PLZ4 == 3611 & gdename == "suederen"
replace PLZ4 = 3635 if PLZ4 == 3611 & gdename == "uebeschl"
replace PLZ4 = 3635 if PLZ4 == 3611 & gdename == "ueheschi"
replace PLZ4 = 3618 if PLZ4 == 3611 & gdename == "wachseidorn"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "althaeusern"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "arisf au"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "aristau birri"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "arlstau birri"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "artisau"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "artstau"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "barri aristau"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "birri"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "birri aristau"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "birrl"
replace PLZ4 = 5628 if PLZ4 == 5649 & gdename == "birrl artstau"
replace PLZ4 = 5626 if PLZ4 == 5649 & gdename == "hermetschwii staffeln"
replace PLZ4 = 5626 if PLZ4 == 5649 & gdename == "hermetschwii stallein"
replace PLZ4 = 5626 if PLZ4 == 5649 & gdename == "hermetschwil"
replace PLZ4 = 5626 if PLZ4 == 5649 & gdename == "hermetschwil staff ein"
replace PLZ4 = 5626 if PLZ4 == 5649 & gdename == "hermetschwil staffein"
replace PLZ4 = 5626 if PLZ4 == 5649 & gdename == "hermetschwil stallein"
replace PLZ4 = 5626 if PLZ4 == 5649 & gdename == "hermetschwll staffeln"
replace PLZ4 = 5608 if PLZ4 == 5649 & gdename == "steffen"
replace PLZ4 = 5608 if PLZ4 == 5649 & gdename == "steffen (ag)"
replace PLZ4 = 5608 if PLZ4 == 5649 & gdename == "steffen ag"
replace PLZ4 = 5608 if PLZ4 == 5649 & gdename == "stehen"
replace PLZ4 = 5608 if PLZ4 == 5649 & gdename == "stetterr ag"
replace PLZ4 = 7315 if PLZ4 == 7311 & gdename == "vaettis"
replace PLZ4 = 7315 if PLZ4 == 7311 & gdename == "vaettis/pfaefers"
replace PLZ4 = 7315 if PLZ4 == 7311 & gdename == "vaettls"
replace PLZ4 = 7315 if PLZ4 == 7311 & gdename == "valens"
replace PLZ4 = 7315 if PLZ4 == 7311 & gdename == "vattis"
replace PLZ4 = 3664 if PLZ4 == 3134 & gdename == "burgfisteln"
replace PLZ4 = 3664 if PLZ4 == 3134 & gdename == "burgfisteln dorf"
replace PLZ4 = 3664 if PLZ4 == 3134 & gdename == "burgistein dorf"
replace PLZ4 = 3664 if PLZ4 == 3134 & gdename == "burgistein dort"
replace PLZ4 = 3664 if PLZ4 == 3134 & gdename == "burgistein station"
replace PLZ4 = 3664 if PLZ4 == 3134 & gdename == "burgisteln dorf"
replace PLZ4 = 3664 if PLZ4 == 3134 & gdename == "burglstein"
replace PLZ4 = 3664 if PLZ4 == 3134 & gdename == "burglsteln station"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "kaeppelimatt"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "ostergau"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "ostergau/wi llisau"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "ostergau/wiliisau land"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "rohrmatt"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "rohrmatt wallisau land"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "rohrmatt willisau land"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "rohrmatt wlllisau land"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "wallisau land"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "wallisau ostergau"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "wiliisau ostergau"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "willisau ostergau"
replace PLZ4 = 6130 if PLZ4 == 6131 & gdename == "wlllisau ostergau"
replace PLZ4 = 4615 if PLZ4 == 4699 & gdename == "allerheiligenberg"
replace PLZ4 = 4615 if PLZ4 == 4699 & gdename == "allerheiligenberg so"
replace PLZ4 = 4615 if PLZ4 == 4699 & gdename == "allerhelligenberg"
replace PLZ4 = 4615 if PLZ4 == 4699 & gdename == "allerhelllgenberg"
replace PLZ4 = 4615 if PLZ4 == 4699 & gdename == "alterheiligenberg"
replace PLZ4 = 4633 if PLZ4 == 4699 & gdename == "hauenstein"
replace PLZ4 = 4633 if PLZ4 == 4699 & gdename == "hauenstein ifenthai"
replace PLZ4 = 4633 if PLZ4 == 4699 & gdename == "hauenstein lfenthal"
replace PLZ4 = 4633 if PLZ4 == 4699 & gdename == "hauensteln"
replace PLZ4 = 4633 if PLZ4 == 4699 & gdename == "hauensteln ifenthal"
replace PLZ4 = 4633 if PLZ4 == 4699 & gdename == "nauenstein"
replace PLZ4 = 4633 if PLZ4 == 4699 & gdename == "nauenstein ifenthal"
replace PLZ4 = 4634 if PLZ4 == 4699 & gdename == "wasen"
replace PLZ4 = 9315 if PLZ4 == 9307 & gdename == "winden"
replace PLZ4 = 9315 if PLZ4 == 9307 & gdename == "winden egnach"
replace PLZ4 = 9315 if PLZ4 == 9307 & gdename == "winden tg"
replace PLZ4 = 9315 if PLZ4 == 9307 & gdename == "winden/egnach"
replace PLZ4 = 3672 if PLZ4 == 3516 & gdename == "aeachlen b. oberdlessbach"
replace PLZ4 = 3672 if PLZ4 == 3516 & gdename == "aeschien"
replace PLZ4 = 3672 if PLZ4 == 3516 & gdename == "aeschien b. oberdiessbach"
replace PLZ4 = 3655 if PLZ4 == 3516 & gdename == "aeschien bei oberdiessbach/sigriswil"
replace PLZ4 = 3672 if PLZ4 == 3516 & gdename == "aeschlen b. oberdiessbach"
replace PLZ4 = 3672 if PLZ4 == 3516 & gdename == "aeschlen bei oberdiesbach"
replace PLZ4 = 3672 if PLZ4 == 3516 & gdename == "aeschlen bei oberdiessbach"
replace PLZ4 = 3655 if PLZ4 == 3516 & gdename == "aeschlen bei oberdiessbach/ sigriswil"
replace PLZ4 = 3655 if PLZ4 == 3516 & gdename == "aeschlen bei oberdiessbach/sigriswil"
replace PLZ4 = 3672 if PLZ4 == 3516 & gdename == "aeschlen bel oberdiessbach"
replace PLZ4 = 3655 if PLZ4 == 3516 & gdename == "aeschlen bel oberdiessbach/ sigriswfl"
replace PLZ4 = 3672 if PLZ4 == 3516 & gdename == "aeschlen kuekenvertriebs ag"
replace PLZ4 = 3674 if PLZ4 == 3516 & gdename == "bielken bei oberdiessbach"
replace PLZ4 = 3674 if PLZ4 == 3516 & gdename == "bleiken"
replace PLZ4 = 3674 if PLZ4 == 3516 & gdename == "bleiken bel oberdiessbach"
replace PLZ4 = 3673 if PLZ4 == 3516 & gdename == "jassbach"
replace PLZ4 = 3673 if PLZ4 == 3516 & gdename == "jassbach zur linden"
replace PLZ4 = 3615 if PLZ4 == 3516 & gdename == "wangelen"
replace PLZ4 = 2336 if PLZ4 == 2311 & gdename == "boechet"
replace PLZ4 = 2300 if PLZ4 == 2311 & gdename == "chaux de fonds"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "emibois be"
replace PLZ4 = 2333 if PLZ4 == 2311 & gdename == "la chaux d abel"
replace PLZ4 = 2314 if PLZ4 == 2311 & gdename == "la corbatiere"
replace PLZ4 = 2314 if PLZ4 == 2311 & gdename == "la torbatiere"
replace PLZ4 = 2336 if PLZ4 == 2311 & gdename == "le boechet"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "les emibois"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "les emibois muriaux"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "les emibols"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "les emibols murlaux"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "les emlbois"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "les emlbois muraux"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "les emlbois muriaux"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "les emlbols"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "les emã¬bols"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "murlaux"
replace PLZ4 = 2338 if PLZ4 == 2311 & gdename == "murtaux"
replace PLZ4 = 3155 if PLZ4 == 3151 & gdename == "heigisried rohrbach"
replace PLZ4 = 3088 if PLZ4 == 3151 & gdename == "heigisried rueeggisberg"
replace PLZ4 = 3155 if PLZ4 == 3151 & gdename == "helgisried rohrbach"
replace PLZ4 = 3088 if PLZ4 == 3151 & gdename == "helgisried rueeggisberg"
replace PLZ4 = 3088 if PLZ4 == 3151 & gdename == "helgisried rueegglsberg"
replace PLZ4 = 3158 if PLZ4 == 3151 & gdename == "hirschmatt"
replace PLZ4 = 3158 if PLZ4 == 3151 & gdename == "riffenmatt"
replace PLZ4 = 3154 if PLZ4 == 3151 & gdename == "rueschegg graben"
replace PLZ4 = 1183 if PLZ4 == 1181 & gdename == "bursitis"
replace PLZ4 = 1183 if PLZ4 == 1181 & gdename == "burslns"
replace PLZ4 = 1186 if PLZ4 == 1181 & gdename == "essertine sur rolle"
replace PLZ4 = 1186 if PLZ4 == 1181 & gdename == "essertines/rolle"
replace PLZ4 = 1186 if PLZ4 == 1181 & gdename == "essertlne sur rolle"
replace PLZ4 = 1186 if PLZ4 == 1181 & gdename == "essertlnes sur rolle"
replace PLZ4 = 1182 if PLZ4 == 1181 & gdename == "gaily"
replace PLZ4 = 1182 if PLZ4 == 1181 & gdename == "giily"
replace PLZ4 = 1182 if PLZ4 == 1181 & gdename == "giity"
replace PLZ4 = 1182 if PLZ4 == 1181 & gdename == "glily"
replace PLZ4 = 1182 if PLZ4 == 1181 & gdename == "gllly"
replace PLZ4 = 1184 if PLZ4 == 1181 & gdename == "lulns"
replace PLZ4 = 1185 if PLZ4 == 1181 & gdename == "mont rolle"
replace PLZ4 = 1185 if PLZ4 == 1181 & gdename == "mont s rolle"
replace PLZ4 = 1185 if PLZ4 == 1181 & gdename == "mont s. rolle"
replace PLZ4 = 1185 if PLZ4 == 1181 & gdename == "mont s/rolie"
replace PLZ4 = 1185 if PLZ4 == 1181 & gdename == "mont s/rolle"
replace PLZ4 = 1185 if PLZ4 == 1181 & gdename == "mont sun rolle"
replace PLZ4 = 1185 if PLZ4 == 1181 & gdename == "mont sur roller"
replace PLZ4 = 1185 if PLZ4 == 1181 & gdename == "mont sur roule"
replace PLZ4 = 1185 if PLZ4 == 1181 & gdename == "monts rolle"
replace PLZ4 = 1187 if PLZ4 == 1181 & gdename == "st oyens"
replace PLZ4 = 1180 if PLZ4 == 1181 & gdename == "tartegniri"
replace PLZ4 = 1184 if PLZ4 == 1181 & gdename == "vinzei"
replace PLZ4 = 1184 if PLZ4 == 1181 & gdename == "vinzel vd"
replace PLZ4 = 1184 if PLZ4 == 1181 & gdename == "vlnzel"
replace PLZ4 = 1184 if PLZ4 == 1181 & gdename == "wins"
replace PLZ4 = 2352 if PLZ4 == 2726 & gdename == "cemievillers"
replace PLZ4 = 2352 if PLZ4 == 2726 & gdename == "cernieviliers"
replace PLZ4 = 2352 if PLZ4 == 2726 & gdename == "cernievillers"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "les"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "les ceriatez"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "les cerlatez/saignelegier"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "les cerlatez/salgnelegier"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "saigneiegier"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "saignelegier ju"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "saignelegiers"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "saignelegiier"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "saignelegler"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "saignelegter"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "saignetegier"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "saine legier"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "sainelegier"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "salgnelegier"
replace PLZ4 = 2350 if PLZ4 == 2726 & gdename == "salgnelegler"
replace PLZ4 = 5466 if PLZ4 == 8434 & gdename == "kaiserstuhi"
replace PLZ4 = 5466 if PLZ4 == 8434 & gdename == "kaiserstuhl ag"
replace PLZ4 = 5466 if PLZ4 == 8434 & gdename == "kalserstuhl"
replace PLZ4 = 5466 if PLZ4 == 8434 & gdename == "kalserstuhl (ag)"
replace PLZ4 = 5466 if PLZ4 == 8434 & gdename == "kalserstuhl ag"
replace PLZ4 = 4614 if PLZ4 == 6352 & gdename == "hagendorn"
replace PLZ4 = 6353 if PLZ4 == 6352 & gdename == "hertenetein"
replace PLZ4 = 6353 if PLZ4 == 6352 & gdename == "hertenstein"
replace PLZ4 = 6353 if PLZ4 == 6352 & gdename == "hertenstein weggis"
replace PLZ4 = 6353 if PLZ4 == 6352 & gdename == "hertensteln"
replace PLZ4 = 6353 if PLZ4 == 6352 & gdename == "hertensteln weggis"
replace PLZ4 = 2345 if PLZ4 == 2724 & gdename == "breuieux"
replace PLZ4 = 2345 if PLZ4 == 2724 & gdename == "breuleux"
replace PLZ4 = 2345 if PLZ4 == 2724 & gdename == "cerneux veusil"
replace PLZ4 = 2345 if PLZ4 == 2724 & gdename == "cerneux veusll"
replace PLZ4 = 2345 if PLZ4 == 2724 & gdename == "le breuleux"
replace PLZ4 = 2345 if PLZ4 == 2724 & gdename == "le cemeux veusil"
replace PLZ4 = 2345 if PLZ4 == 2724 & gdename == "le cerneux veusil"
replace PLZ4 = 2345 if PLZ4 == 2724 & gdename == "le roselet"
replace PLZ4 = 2345 if PLZ4 == 2724 & gdename == "lea breuleux"
replace PLZ4 = 2345 if PLZ4 == 2724 & gdename == "les breuieux"
replace PLZ4 = 6048 if PLZ4 == 6046 & gdename == "niklausen"
replace PLZ4 = 6048 if PLZ4 == 6046 & gdename == "st. niklausen"
replace PLZ4 = 6048 if PLZ4 == 6046 & gdename == "st. niklausen (lu)"
replace PLZ4 = 6048 if PLZ4 == 6046 & gdename == "st. niklausen horw"
replace PLZ4 = 6048 if PLZ4 == 6046 & gdename == "st. niklausen/lu"
replace PLZ4 = 6048 if PLZ4 == 6046 & gdename == "st. nlklausen"
replace PLZ4 = 6048 if PLZ4 == 6046 & gdename == "st.niklausen"
replace PLZ4 = 3046 if PLZ4 == 6046 & gdename == "wahlendorf"
replace PLZ4 = 4525 if PLZ4 == 4511 & gdename == "balm bel guensberg"
replace PLZ4 = 4525 if PLZ4 == 4511 & gdename == "balm guensberg"
replace PLZ4 = 4525 if PLZ4 == 4511 & gdename == "balmberg"
replace PLZ4 = 4539 if PLZ4 == 4511 & gdename == "famem"
replace PLZ4 = 4539 if PLZ4 == 4511 & gdename == "famern"
replace PLZ4 = 4539 if PLZ4 == 4511 & gdename == "farnem"
replace PLZ4 = 4539 if PLZ4 == 4511 & gdename == "fernern"
replace PLZ4 = 4558 if PLZ4 == 4511 & gdename == "hersiwll"
replace PLZ4 = 4558 if PLZ4 == 4511 & gdename == "herslwil"
replace PLZ4 = 4557 if PLZ4 == 4511 & gdename == "horriwii"
replace PLZ4 = 4557 if PLZ4 == 4511 & gdename == "horriwil so"
replace PLZ4 = 4557 if PLZ4 == 4511 & gdename == "horriwll"
replace PLZ4 = 4557 if PLZ4 == 4511 & gdename == "horriwll so"
replace PLZ4 = 4557 if PLZ4 == 4511 & gdename == "horrlwil"
replace PLZ4 = 4557 if PLZ4 == 4511 & gdename == "horrlwll"
replace PLZ4 = 4557 if PLZ4 == 4511 & gdename == "horrwil"
replace PLZ4 = 4535 if PLZ4 == 4511 & gdename == "huberadorf"
replace PLZ4 = 4535 if PLZ4 == 4511 & gdename == "hubersdorf so"
replace PLZ4 = 4535 if PLZ4 == 4511 & gdename == "kammersohr"
replace PLZ4 = 4535 if PLZ4 == 4511 & gdename == "kammersrohr hubersdorf"
replace PLZ4 = 4523 if PLZ4 == 4511 & gdename == "niederwal"
replace PLZ4 = 4523 if PLZ4 == 4511 & gdename == "niederwal so"
replace PLZ4 = 4539 if PLZ4 == 4511 & gdename == "rumasberg"
replace PLZ4 = 4539 if PLZ4 == 4511 & gdename == "rumlaberg"
replace PLZ4 = 4539 if PLZ4 == 4511 & gdename == "rumlsberg"
replace PLZ4 = 4558 if PLZ4 == 4511 & gdename == "winistorf"
replace PLZ4 = 4558 if PLZ4 == 4511 & gdename == "wlnistorf"
replace PLZ4 = 4556 if PLZ4 == 3361 & gdename == "aeschl"
replace PLZ4 = 4556 if PLZ4 == 3361 & gdename == "aeschl so"
replace PLZ4 = 4556 if PLZ4 == 3361 & gdename == "bocken"
replace PLZ4 = 4556 if PLZ4 == 3361 & gdename == "boiken"
replace PLZ4 = 4556 if PLZ4 == 3361 & gdename == "bolken ebersol"
replace PLZ4 = 4556 if PLZ4 == 3361 & gdename == "burgaeschi"
replace PLZ4 = 4556 if PLZ4 == 3361 & gdename == "burgaeschl"
replace PLZ4 = 3429 if PLZ4 == 3361 & gdename == "heilsau"
replace PLZ4 = 3429 if PLZ4 == 3361 & gdename == "helisau"
replace PLZ4 = 3373 if PLZ4 == 3361 & gdename == "helmenhausen"
replace PLZ4 = 3429 if PLZ4 == 3361 & gdename == "hoerstetten"
replace PLZ4 = 3373 if PLZ4 == 3361 & gdename == "roethenbach"
replace PLZ4 = 3373 if PLZ4 == 3361 & gdename == "roethenbach b.herzogenbuch"
replace PLZ4 = 3373 if PLZ4 == 3361 & gdename == "roethenbach bel herzogenbuchsee"
replace PLZ4 = 3373 if PLZ4 == 3361 & gdename == "roethenbach bet herzogenbuchsee"
replace PLZ4 = 3373 if PLZ4 == 3361 & gdename == "roethenbach i.e"
replace PLZ4 = 3373 if PLZ4 == 3361 & gdename == "rothenbach bei herzogenbuchsee"
replace PLZ4 = 3374 if PLZ4 == 3361 & gdename == "wangenrled"
replace PLZ4 = 3374 if PLZ4 == 3361 & gdename == "wangeriried"
replace PLZ4 = 3372 if PLZ4 == 3361 & gdename == "wanzenwil"
replace PLZ4 = 3372 if PLZ4 == 3361 & gdename == "wanzenwll"
replace PLZ4 = 1965 if PLZ4 == 1956 & gdename == "st germain"
replace PLZ4 = 1955 if PLZ4 == 1956 & gdename == "st pierce de clages"
replace PLZ4 = 1955 if PLZ4 == 1956 & gdename == "st pierre de cl ages"
replace PLZ4 = 1955 if PLZ4 == 1956 & gdename == "st pierre de clages"
replace PLZ4 = 1955 if PLZ4 == 1956 & gdename == "st pierre de clages/chamoson"
replace PLZ4 = 1955 if PLZ4 == 1956 & gdename == "st. pierre de ciages"
replace PLZ4 = 1955 if PLZ4 == 1956 & gdename == "st. pierre de clages"
replace PLZ4 = 6260 if PLZ4 == 6261 & gdename == "reidermoos"
replace PLZ4 = 6260 if PLZ4 == 6261 & gdename == "reidermoos/reiden"
replace PLZ4 = 6260 if PLZ4 == 6261 & gdename == "reldermoos"
replace PLZ4 = 5058 if PLZ4 == 6261 & gdename == "wiliberg hintermoos"
replace PLZ4 = 5058 if PLZ4 == 6261 & gdename == "wlliberg"
replace PLZ4 = 5058 if PLZ4 == 6261 & gdename == "wlliberg hintermoos"
replace PLZ4 = 8254 if PLZ4 == 8251 & gdename == "basedingen"
replace PLZ4 = 8254 if PLZ4 == 8251 & gdename == "bassdingen"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "neu paradies unterschlatt"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "neu paradies/unterschlaff"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "neu paradies/unterschlatt"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "paradies/schlatt tg"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schiatt bei diessenhofen"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schlaff b. diessenhofen"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schlaff bei diessenhofen"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schlaft"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schlaft b. diessenhofen"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schlatt b. diessenhofen"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schlatt b. dlessenhofen"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schlatt bei diessenhofen"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schlatt diesenhofen"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schlatt dlesenhofen"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "schlatt tg"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "unterschiatt"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "unterschlaff tg"
replace PLZ4 = 8252 if PLZ4 == 8251 & gdename == "unterschlatt tg"
replace PLZ4 = 9633 if PLZ4 == 9128 & gdename == "bachli (hemberg)"
replace PLZ4 = 9633 if PLZ4 == 9128 & gdename == "baechii"
replace PLZ4 = 9633 if PLZ4 == 9128 & gdename == "baechli"
replace PLZ4 = 9633 if PLZ4 == 9128 & gdename == "baechli (hemberg)"
replace PLZ4 = 9633 if PLZ4 == 9128 & gdename == "baechll"
replace PLZ4 = 9633 if PLZ4 == 9128 & gdename == "baechll (hemberg)"
replace PLZ4 = 9126 if PLZ4 == 9128 & gdename == "necker"
replace PLZ4 = 1619 if PLZ4 == 1622 & gdename == "les paccots"
replace PLZ4 = 4573 if PLZ4 == 4572 & gdename == "ammannseg"
replace PLZ4 = 4573 if PLZ4 == 4572 & gdename == "ammannsegg"
replace PLZ4 = 4573 if PLZ4 == 4572 & gdename == "ammannsegg/kriegstetten"
replace PLZ4 = 4573 if PLZ4 == 4572 & gdename == "ammansegg"
replace PLZ4 = 4573 if PLZ4 == 4572 & gdename == "ammansegg (so)"
replace PLZ4 = 4573 if PLZ4 == 4572 & gdename == "arrmannsegg"
replace PLZ4 = 7325 if PLZ4 == 7321 & gdename == "schwendi (weiestannental)"
replace PLZ4 = 7325 if PLZ4 == 7321 & gdename == "schwendi (weisstannental)"
replace PLZ4 = 7325 if PLZ4 == 7321 & gdename == "schwendl (welsstannental)"
replace PLZ4 = 7326 if PLZ4 == 7321 & gdename == "weisstannen"
replace PLZ4 = 6719 if PLZ4 == 6711 & gdename == "aquilia"
replace PLZ4 = 6720 if PLZ4 == 6711 & gdename == "campa/blenio"
replace PLZ4 = 6720 if PLZ4 == 6711 & gdename == "campo bienio"
replace PLZ4 = 6720 if PLZ4 == 6711 & gdename == "campo blenlo"
replace PLZ4 = 6720 if PLZ4 == 6711 & gdename == "campo/blenlo"
replace PLZ4 = 6722 if PLZ4 == 6711 & gdename == "corzoneso piano"
replace PLZ4 = 6722 if PLZ4 == 6711 & gdename == "corzoneso plano"
replace PLZ4 = 6722 if PLZ4 == 6711 & gdename == "corzoresco"
replace PLZ4 = 6724 if PLZ4 == 6711 & gdename == "largarlo"
replace PLZ4 = 6716 if PLZ4 == 6711 & gdename == "loittigna"
replace PLZ4 = 6716 if PLZ4 == 6711 & gdename == "lottlgna"
replace PLZ4 = 6721 if PLZ4 == 6711 & gdename == "ludlam)"
replace PLZ4 = 6721 if PLZ4 == 6711 & gdename == "ludlano"
replace PLZ4 = 6723 if PLZ4 == 6711 & gdename == "maroita"
replace PLZ4 = 6723 if PLZ4 == 6711 & gdename == "marotta"
replace PLZ4 = 6720 if PLZ4 == 6711 & gdename == "motto"
replace PLZ4 = 6720 if PLZ4 == 6711 & gdename == "motto (bienio)"
replace PLZ4 = 6720 if PLZ4 == 6711 & gdename == "motto (blenio)"
replace PLZ4 = 6720 if PLZ4 == 6711 & gdename == "motto (blenio)/dongio"
replace PLZ4 = 6724 if PLZ4 == 6711 & gdename == "ponte valentino"
replace PLZ4 = 6723 if PLZ4 == 6711 & gdename == "pruglasco"
replace PLZ4 = 5333 if PLZ4 == 8439 & gdename == "baidingen"
replace PLZ4 = 5333 if PLZ4 == 8439 & gdename == "baldiingen"
replace PLZ4 = 5333 if PLZ4 == 8439 & gdename == "betdingen"
replace PLZ4 = 5334 if PLZ4 == 8439 & gdename == "boeblkon"
replace PLZ4 = 5463 if PLZ4 == 8439 & gdename == "mellstorf/wislikofen"
replace PLZ4 = 5463 if PLZ4 == 8439 & gdename == "mellstorf/wisllkofen"
replace PLZ4 = 5464 if PLZ4 == 8439 & gdename == "ruemikon ag"
replace PLZ4 = 5464 if PLZ4 == 8439 & gdename == "ruemlkon"
replace PLZ4 = 5464 if PLZ4 == 8439 & gdename == "ruemlkon ag"
replace PLZ4 = 5464 if PLZ4 == 8439 & gdename == "rumikon"
replace PLZ4 = 5464 if PLZ4 == 8439 & gdename == "rumikon ag"
replace PLZ4 = 5462 if PLZ4 == 8439 & gdename == "sigiistorf"
replace PLZ4 = 5462 if PLZ4 == 8439 & gdename == "sigilstorf"
replace PLZ4 = 5462 if PLZ4 == 8439 & gdename == "siglisdorf"
replace PLZ4 = 5462 if PLZ4 == 8439 & gdename == "siglistori"
replace PLZ4 = 5462 if PLZ4 == 8439 & gdename == "siglistort"
replace PLZ4 = 5462 if PLZ4 == 8439 & gdename == "sigllstorf"
replace PLZ4 = 5462 if PLZ4 == 8439 & gdename == "sigã¬istorf"
replace PLZ4 = 5462 if PLZ4 == 8439 & gdename == "slglisdorf"
replace PLZ4 = 5462 if PLZ4 == 8439 & gdename == "slglistorf"
replace PLZ4 = 5462 if PLZ4 == 8439 & gdename == "sïgliistorf"
replace PLZ4 = 5463 if PLZ4 == 8439 & gdename == "wislikof en"
replace PLZ4 = 5463 if PLZ4 == 8439 & gdename == "wisllkofen"
replace PLZ4 = 5463 if PLZ4 == 8439 & gdename == "wlsllkofen"
replace PLZ4 = 1955 if PLZ4 == 1916 & gdename == "st pierre d.c"
replace PLZ4 = 1955 if PLZ4 == 1916 & gdename == "st pierre de clages"
replace PLZ4 = 1955 if PLZ4 == 1916 & gdename == "st plerre de clages"
replace PLZ4 = 1955 if PLZ4 == 1916 & gdename == "st. pier de clages"
replace PLZ4 = 1955 if PLZ4 == 1916 & gdename == "st. pierre"
replace PLZ4 = 1955 if PLZ4 == 1916 & gdename == "st. pierre de clages"
replace PLZ4 = 1955 if PLZ4 == 1916 & gdename == "st. pierrre de clages"
replace PLZ4 = 1955 if PLZ4 == 1916 & gdename == "st. plerre de clages"
replace PLZ4 = 1955 if PLZ4 == 1916 & gdename == "st. plerrre de clages"
replace PLZ4 = 5316 if PLZ4 == 4354 & gdename == "felsenau"
replace PLZ4 = 5316 if PLZ4 == 4354 & gdename == "felsenau (ag)"
replace PLZ4 = 5324 if PLZ4 == 4354 & gdename == "fuli reuenthal"
replace PLZ4 = 5324 if PLZ4 == 4354 & gdename == "full"
replace PLZ4 = 5324 if PLZ4 == 4354 & gdename == "full reuenthai"
replace PLZ4 = 5324 if PLZ4 == 4354 & gdename == "reuenthal"
replace PLZ4 = 5037 if PLZ4 == 5038 & gdename == "ober muhen"
replace PLZ4 = 5037 if PLZ4 == 5038 & gdename == "obermuhen"
replace PLZ4 = 1660 if PLZ4 == 1834 & gdename == "les moulins"
replace PLZ4 = 1660 if PLZ4 == 1834 & gdename == "les moulins/ chateau d oex"
replace PLZ4 = 1660 if PLZ4 == 1834 & gdename == "les moulins/chateau d oex"
replace PLZ4 = 9562 if PLZ4 == 9563 & gdename == "buch bei maerwil"
replace PLZ4 = 9565 if PLZ4 == 9563 & gdename == "oppikon"
replace PLZ4 = 9565 if PLZ4 == 9563 & gdename == "opplkon"
replace PLZ4 = 9565 if PLZ4 == 9563 & gdename == "schmidshof"
replace PLZ4 = 9565 if PLZ4 == 9563 & gdename == "schmidshof buch"
replace PLZ4 = 9565 if PLZ4 == 9563 & gdename == "schmidshof/bussnang"
replace PLZ4 = 9565 if PLZ4 == 9563 & gdename == "schmldshof"
replace PLZ4 = 1607 if PLZ4 == 1501 & gdename == "paiezieux village"
replace PLZ4 = 1607 if PLZ4 == 1501 & gdename == "palezieux village"
replace PLZ4 = 1607 if PLZ4 == 1501 & gdename == "palezleux village"
replace PLZ4 = 7223 if PLZ4 == 7221 & gdename == "buchen"
replace PLZ4 = 7223 if PLZ4 == 7221 & gdename == "buchen i. p"
replace PLZ4 = 7223 if PLZ4 == 7221 & gdename == "buchen i.p"
replace PLZ4 = 7223 if PLZ4 == 7221 & gdename == "buchen staad"
replace PLZ4 = 7223 if PLZ4 == 7221 & gdename == "buchen stead"
replace PLZ4 = 7220 if PLZ4 == 7221 & gdename == "schuders"
replace PLZ4 = 7550 if PLZ4 == 7221 & gdename == "schuls"
replace PLZ4 = 7550 if PLZ4 == 7221 & gdename == "scuol"
replace PLZ4 = 7550 if PLZ4 == 7221 & gdename == "scuols/schuls"
replace PLZ4 = 8259 if PLZ4 == 8256 & gdename == "etzwflen"
replace PLZ4 = 8259 if PLZ4 == 8256 & gdename == "etzwilen"
replace PLZ4 = 8259 if PLZ4 == 8256 & gdename == "etzwilen tg"
replace PLZ4 = 8259 if PLZ4 == 8256 & gdename == "etzwilen/kaltenbach"
replace PLZ4 = 8259 if PLZ4 == 8256 & gdename == "etzwlien"
replace PLZ4 = 8259 if PLZ4 == 8256 & gdename == "etzwllen"
replace PLZ4 = 8259 if PLZ4 == 8256 & gdename == "etzwllen/kaltenbach"
replace PLZ4 = 8259 if PLZ4 == 8256 & gdename == "rheinklingen"
replace PLZ4 = 8532 if PLZ4 == 8534 & gdename == "weiningen tg"
replace PLZ4 = 8532 if PLZ4 == 8534 & gdename == "welningen"
replace PLZ4 = 8532 if PLZ4 == 8534 & gdename == "welningen tg"
replace PLZ4 = 8532 if PLZ4 == 8534 & gdename == "wetningen"
replace PLZ4 = 8532 if PLZ4 == 8534 & gdename == "wetningen tg"
replace PLZ4 = 5074 if PLZ4 == 5268 & gdename == "elken"
replace PLZ4 = 1880 if PLZ4 == 1881 & gdename == "frenieres"
replace PLZ4 = 1880 if PLZ4 == 1881 & gdename == "les devens"
replace PLZ4 = 1880 if PLZ4 == 1881 & gdename == "les pians sur bex"
replace PLZ4 = 1880 if PLZ4 == 1881 & gdename == "les plans"
replace PLZ4 = 1880 if PLZ4 == 1881 & gdename == "les plans sur bex"
replace PLZ4 = 1880 if PLZ4 == 1881 & gdename == "pians sur bex"
replace PLZ4 = 1880 if PLZ4 == 1881 & gdename == "plans sur bex"
replace PLZ4 = 1880 if PLZ4 == 1881 & gdename == "posses bex"
replace PLZ4 = 2149 if PLZ4 == 2093 & gdename == "brot dessus"
replace PLZ4 = 3454 if PLZ4 == 3451 & gdename == "griesbach sumiswald"
replace PLZ4 = 3416 if PLZ4 == 3451 & gdename == "haeusernmoos"
replace PLZ4 = 3416 if PLZ4 == 3451 & gdename == "haeusernmoos i. e."
replace PLZ4 = 3464 if PLZ4 == 3451 & gdename == "schmidigen"
replace PLZ4 = 3464 if PLZ4 == 3451 & gdename == "schmidigen muehleweg"
replace PLZ4 = 3464 if PLZ4 == 3451 & gdename == "schmidigen muehlweg"
replace PLZ4 = 3464 if PLZ4 == 3451 & gdename == "schmidigen muhleweg"
replace PLZ4 = 3464 if PLZ4 == 3451 & gdename == "schmidlgen muehleweg"
replace PLZ4 = 3464 if PLZ4 == 3451 & gdename == "schmldigen muehleweg"
replace PLZ4 = 3462 if PLZ4 == 3451 & gdename == "weber i e"
replace PLZ4 = 3462 if PLZ4 == 3451 & gdename == "weier be"
replace PLZ4 = 3462 if PLZ4 == 3451 & gdename == "weier i.e"
replace PLZ4 = 3462 if PLZ4 == 3451 & gdename == "weier im emmental"
replace PLZ4 = 3462 if PLZ4 == 3451 & gdename == "weler im emmental"
replace PLZ4 = 8914 if PLZ4 == 8920 & gdename == "aeugstertal"
replace PLZ4 = 8914 if PLZ4 == 8920 & gdename == "aeugstertal/aeugst"
replace PLZ4 = 8914 if PLZ4 == 8920 & gdename == "aeugstertal/aeugst am albis"
replace PLZ4 = 1884 if PLZ4 == 1883 & gdename == "arveyes"
replace PLZ4 = 1884 if PLZ4 == 1883 & gdename == "arveyes ollon"
replace PLZ4 = 1884 if PLZ4 == 1883 & gdename == "arveyes villars"
replace PLZ4 = 1884 if PLZ4 == 1883 & gdename == "arveyes/ollon"
replace PLZ4 = 2825 if PLZ4 == 2801 & gdename == "ceurchapolx"
replace PLZ4 = 2843 if PLZ4 == 2801 & gdename == "chaetillon"
replace PLZ4 = 2843 if PLZ4 == 2801 & gdename == "chatilion"
replace PLZ4 = 2823 if PLZ4 == 2801 & gdename == "courcellon"
replace PLZ4 = 2823 if PLZ4 == 2801 & gdename == "courcelon"
replace PLZ4 = 2825 if PLZ4 == 2801 & gdename == "courchapolx"
replace PLZ4 = 2802 if PLZ4 == 2801 & gdename == "devetier"
replace PLZ4 = 2827 if PLZ4 == 2801 & gdename == "la scheulten"
replace PLZ4 = 2827 if PLZ4 == 2801 & gdename == "merveller"
replace PLZ4 = 2828 if PLZ4 == 2801 & gdename == "montseveller"
replace PLZ4 = 2812 if PLZ4 == 2801 & gdename == "moveller"
replace PLZ4 = 2842 if PLZ4 == 2801 & gdename == "rossemalson"
replace PLZ4 = 2842 if PLZ4 == 2801 & gdename == "rossmaisen"
replace PLZ4 = 2842 if PLZ4 == 2801 & gdename == "rossmaison"
replace PLZ4 = 2829 if PLZ4 == 2801 & gdename == "vernies"
replace PLZ4 = 6373 if PLZ4 == 6366 & gdename == "buergenstock"
replace PLZ4 = 6373 if PLZ4 == 6366 & gdename == "buergenstock/ ennetbuergen"
replace PLZ4 = 6373 if PLZ4 == 6366 & gdename == "buergenstock/ennetbuergen"
replace PLZ4 = 6373 if PLZ4 == 6366 & gdename == "buergenstock/stansstad"
replace PLZ4 = 6373 if PLZ4 == 6366 & gdename == "burgenstock"
replace PLZ4 = 8404 if PLZ4 == 8473 & gdename == "reutlingen"
replace PLZ4 = 8404 if PLZ4 == 8473 & gdename == "reutlingen (winterthur)"
replace PLZ4 = 8404 if PLZ4 == 8473 & gdename == "reutlingen winterthur"
replace PLZ4 = 3927 if PLZ4 == 3921 & gdename == "embal"
replace PLZ4 = 3927 if PLZ4 == 3921 & gdename == "herbriggen"
replace PLZ4 = 3927 if PLZ4 == 3921 & gdename == "herbrlggen"
replace PLZ4 = 3928 if PLZ4 == 3921 & gdename == "rande"
replace PLZ4 = 3928 if PLZ4 == 3921 & gdename == "renda"
replace PLZ4 = 3924 if PLZ4 == 3921 & gdename == "st. nlklaus"
replace PLZ4 = 3933 if PLZ4 == 3921 & gdename == "staddenried"
replace PLZ4 = 3933 if PLZ4 == 3921 & gdename == "staidenried"
replace PLZ4 = 3933 if PLZ4 == 3921 & gdename == "staldenrled"
replace PLZ4 = 3929 if PLZ4 == 3921 & gdename == "taesch vs"
replace PLZ4 = 3929 if PLZ4 == 3921 & gdename == "tisch"
replace PLZ4 = 4543 if PLZ4 == 4707 & gdename == "deftingen"
replace PLZ4 = 4543 if PLZ4 == 4707 & gdename == "deiingen"
replace PLZ4 = 4543 if PLZ4 == 4707 & gdename == "delitingen"
replace PLZ4 = 4543 if PLZ4 == 4707 & gdename == "dellingen"
replace PLZ4 = 4543 if PLZ4 == 4707 & gdename == "deltingen"
replace PLZ4 = 4543 if PLZ4 == 4707 & gdename == "dettingen"
replace PLZ4 = 1164 if PLZ4 == 3211 & gdename == "buchiilon"
replace PLZ4 = 3215 if PLZ4 == 3211 & gdename == "buechsien"
replace PLZ4 = 3213 if PLZ4 == 3211 & gdename == "kieinboesingen"
replace PLZ4 = 3213 if PLZ4 == 3211 & gdename == "klelnboesingen"
replace PLZ4 = 3213 if PLZ4 == 3211 & gdename == "liebislorf"
replace PLZ4 = 3213 if PLZ4 == 3211 & gdename == "liebistorf gruenenburg"
replace PLZ4 = 3213 if PLZ4 == 3211 & gdename == "lieblstorf gruenenburg"
replace PLZ4 = 3213 if PLZ4 == 3211 & gdename == "llebistorf"
replace PLZ4 = 3213 if PLZ4 == 3211 & gdename == "llebistorf gruenenburg"
replace PLZ4 = 3215 if PLZ4 == 3211 & gdename == "lutrigen"
replace PLZ4 = 3215 if PLZ4 == 3211 & gdename == "lutrlgen"
replace PLZ4 = 3215 if PLZ4 == 3211 & gdename == "luttigen3452307"
replace PLZ4 = 3185 if PLZ4 == 3211 & gdename == "ried bei berg (schmitten)"
replace PLZ4 = 3216 if PLZ4 == 3211 & gdename == "ried bel karzers"
replace PLZ4 = 3216 if PLZ4 == 3211 & gdename == "ried bel kerzers"
replace PLZ4 = 3216 if PLZ4 == 3211 & gdename == "ried/kerzers"
replace PLZ4 = 3216 if PLZ4 == 3211 & gdename == "rled bel kerzers"
replace PLZ4 = 3182 if PLZ4 == 3211 & gdename == "uebistorf gruenenburg"
replace PLZ4 = 3214 if PLZ4 == 3211 & gdename == "uimiz"
replace PLZ4 = 3214 if PLZ4 == 3211 & gdename == "ulmiz sg"
replace PLZ4 = 3214 if PLZ4 == 3211 & gdename == "ulmlz"
replace PLZ4 = 7278 if PLZ4 == 7275 & gdename == "davos monstein"
replace PLZ4 = 7276 if PLZ4 == 7275 & gdename == "frauenkirch"
replace PLZ4 = 7276 if PLZ4 == 7275 & gdename == "frauenkirch/davos"
replace PLZ4 = 7277 if PLZ4 == 7275 & gdename == "glaris davos"
replace PLZ4 = 7278 if PLZ4 == 7275 & gdename == "monstein"
replace PLZ4 = 1742 if PLZ4 == 1751 & gdename == "autlgny"
replace PLZ4 = 1741 if PLZ4 == 1751 & gdename == "cottans"
replace PLZ4 = 1741 if PLZ4 == 1751 & gdename == "cottene fr"
replace PLZ4 = 1741 if PLZ4 == 1751 & gdename == "cottons"
replace PLZ4 = 1749 if PLZ4 == 1751 & gdename == "middes fr"
replace PLZ4 = 1740 if PLZ4 == 1751 & gdename == "neyrãºz"
replace PLZ4 = 1746 if PLZ4 == 1751 & gdename == "prez"
replace PLZ4 = 1746 if PLZ4 == 1751 & gdename == "prez vers noreez"
replace PLZ4 = 1748 if PLZ4 == 1751 & gdename == "tomy le grand"
replace PLZ4 = 1748 if PLZ4 == 1751 & gdename == "torny ie grand"
replace PLZ4 = 1695 if PLZ4 == 1751 & gdename == "viilareel le gibloux"
replace PLZ4 = 1695 if PLZ4 == 1751 & gdename == "viilarsel le giboux"
replace PLZ4 = 1690 if PLZ4 == 1751 & gdename == "viliarimboud"
replace PLZ4 = 1695 if PLZ4 == 1751 & gdename == "viliarsel le gibloux"
replace PLZ4 = 1690 if PLZ4 == 1751 & gdename == "villarlmboud"
replace PLZ4 = 1752 if PLZ4 == 1751 & gdename == "villars sur glaene"
replace PLZ4 = 1752 if PLZ4 == 1751 & gdename == "villars sur glene"
replace PLZ4 = 1695 if PLZ4 == 1751 & gdename == "villarsel"
replace PLZ4 = 1695 if PLZ4 == 1751 & gdename == "villarsel le gibioux"
replace PLZ4 = 1695 if PLZ4 == 1751 & gdename == "villarsel le giboux"
replace PLZ4 = 1695 if PLZ4 == 1751 & gdename == "villarsel le glboux"
replace PLZ4 = 1690 if PLZ4 == 1751 & gdename == "villaz st pierre"
replace PLZ4 = 1690 if PLZ4 == 1751 & gdename == "vlllarimboud"
replace PLZ4 = 1690 if PLZ4 == 1751 & gdename == "vlllarlmboud"
replace PLZ4 = 1695 if PLZ4 == 1751 & gdename == "vlllarsel le gibloux"
replace PLZ4 = 3262 if PLZ4 == 3055 & gdename == "suberg"
replace PLZ4 = 3262 if PLZ4 == 3055 & gdename == "surberg"
replace PLZ4 = 3306 if PLZ4 == 3349 & gdename == "etzeikofen"
replace PLZ4 = 3309 if PLZ4 == 3349 & gdename == "kemenried"
replace PLZ4 = 3309 if PLZ4 == 3349 & gdename == "kemenrled"
replace PLZ4 = 3309 if PLZ4 == 3349 & gdename == "kernenrled"
replace PLZ4 = 3315 if PLZ4 == 3349 & gdename == "kraeiligen"
replace PLZ4 = 3315 if PLZ4 == 3349 & gdename == "kraeillgen"
replace PLZ4 = 3315 if PLZ4 == 3349 & gdename == "kraelligen"
replace PLZ4 = 3315 if PLZ4 == 3349 & gdename == "kraellligen"
replace PLZ4 = 3315 if PLZ4 == 3349 & gdename == "krailigen"
replace PLZ4 = 3305 if PLZ4 == 3349 & gdename == "lffwil"
replace PLZ4 = 3303 if PLZ4 == 3349 & gdename == "zuzwll"
replace PLZ4 = 3303 if PLZ4 == 3349 & gdename == "zuzwll (be)"
replace PLZ4 = 3303 if PLZ4 == 3349 & gdename == "zuzwll be"
replace PLZ4 = 6110 if PLZ4 == 6116 & gdename == "wolhi"
replace PLZ4 = 6110 if PLZ4 == 6116 & gdename == "wolhusen markfjwerthenstein"
replace PLZ4 = 6110 if PLZ4 == 6116 & gdename == "wolhusen markt"
replace PLZ4 = 6110 if PLZ4 == 6116 & gdename == "wolhusen markt/ werthenstein"
replace PLZ4 = 6110 if PLZ4 == 6116 & gdename == "wolhusen markt/werthenstein"
replace PLZ4 = 6110 if PLZ4 == 6116 & gdename == "wolhusen markt/werthensteln"
replace PLZ4 = 6110 if PLZ4 == 6116 & gdename == "wolhusen/ werthenstein"
replace PLZ4 = 1172 if PLZ4 == 1171 & gdename == "bougv villars"
replace PLZ4 = 1172 if PLZ4 == 1171 & gdename == "bougy"
replace PLZ4 = 1172 if PLZ4 == 1171 & gdename == "bougy viliars"
replace PLZ4 = 1172 if PLZ4 == 1171 & gdename == "bougy villard"
replace PLZ4 = 1173 if PLZ4 == 1171 & gdename == "fechy dessus"
replace PLZ4 = 1174 if PLZ4 == 1171 & gdename == "plzy"
replace PLZ4 = 1176 if PLZ4 == 1171 & gdename == "st livres"
replace PLZ4 = 1176 if PLZ4 == 1171 & gdename == "st. livres"
replace PLZ4 = 8139 if PLZ4 == 3136 & gdename == "gattikon"
replace PLZ4 = 3662 if PLZ4 == 3136 & gdename == "saftigen"
replace PLZ4 = 3662 if PLZ4 == 3136 & gdename == "seftlgen"
replace PLZ4 = 3662 if PLZ4 == 3136 & gdename == "seifigen"
replace PLZ4 = 3662 if PLZ4 == 3136 & gdename == "seltigen"
replace PLZ4 = 3475 if PLZ4 == 3354 & gdename == "hermiewil"
replace PLZ4 = 3475 if PLZ4 == 3354 & gdename == "hermiswll"
replace PLZ4 = 3475 if PLZ4 == 3354 & gdename == "hermlswil"
replace PLZ4 = 3475 if PLZ4 == 3354 & gdename == "riedtwal"
replace PLZ4 = 3475 if PLZ4 == 3354 & gdename == "riedtwii"
replace PLZ4 = 3475 if PLZ4 == 3354 & gdename == "riedtwil"
replace PLZ4 = 3475 if PLZ4 == 3354 & gdename == "rietwil"
replace PLZ4 = 3475 if PLZ4 == 3354 & gdename == "rledtwil"
replace PLZ4 = 4922 if PLZ4 == 3357 & gdename == "buetzberg"
replace PLZ4 = 3672 if PLZ4 == 3515 & gdename == "oberdaessbach"
replace PLZ4 = 3672 if PLZ4 == 3515 & gdename == "oberdiesabach"
replace PLZ4 = 3672 if PLZ4 == 3515 & gdename == "oberdiesbach"
replace PLZ4 = 3672 if PLZ4 == 3515 & gdename == "oberdiessbagh"
replace PLZ4 = 3672 if PLZ4 == 3515 & gdename == "oberdlesabach"
replace PLZ4 = 3672 if PLZ4 == 3515 & gdename == "oberdlessbach"
replace PLZ4 = 3672 if PLZ4 == 3515 & gdename == "oberdliessbach"
replace PLZ4 = 3672 if PLZ4 == 3515 & gdename == "oberdã­essbach"
replace PLZ4 = 3672 if PLZ4 == 3515 & gdename == "obergiessbach"
replace PLZ4 = 3766 if PLZ4 == 3767 & gdename == "boitigen/weissenbach"
replace PLZ4 = 3766 if PLZ4 == 3767 & gdename == "boltigen/weissenbach"
replace PLZ4 = 3766 if PLZ4 == 3767 & gdename == "weissenbach"
replace PLZ4 = 3766 if PLZ4 == 3767 & gdename == "weissenbach boltigen"
replace PLZ4 = 3766 if PLZ4 == 3767 & gdename == "weissenbach/boltigen"
replace PLZ4 = 3766 if PLZ4 == 3767 & gdename == "welissenbach/boltigen"
replace PLZ4 = 3766 if PLZ4 == 3767 & gdename == "welssenbach"
replace PLZ4 = 2340 if PLZ4 == 2725 & gdename == "le noir mont"
replace PLZ4 = 2340 if PLZ4 == 2725 & gdename == "le nolrmont"
replace PLZ4 = 2340 if PLZ4 == 2725 & gdename == "noirmont"
replace PLZ4 = 2340 if PLZ4 == 2725 & gdename == "noirmont le"
replace PLZ4 = 2340 if PLZ4 == 2725 & gdename == "nolrmont"
replace PLZ4 = 8758 if PLZ4 == 8875 & gdename == "obstaiden"
replace PLZ4 = 8758 if PLZ4 == 8875 & gdename == "obstalgen"
replace PLZ4 = 1610 if PLZ4 == 1502 & gdename == "chaetiliens"
replace PLZ4 = 1610 if PLZ4 == 1502 & gdename == "chaetillens"
replace PLZ4 = 1610 if PLZ4 == 1502 & gdename == "chaetlllens"
replace PLZ4 = 1610 if PLZ4 == 1502 & gdename == "chatiilens"
replace PLZ4 = 1610 if PLZ4 == 1502 & gdename == "chatiliens"
replace PLZ4 = 1610 if PLZ4 == 1502 & gdename == "chatitlens"
replace PLZ4 = 1610 if PLZ4 == 1502 & gdename == "chatlllens"
replace PLZ4 = 1610 if PLZ4 == 1502 & gdename == "chetillens"
replace PLZ4 = 1610 if PLZ4 == 1502 & gdename == "vulbroye"
replace PLZ4 = 2416 if PLZ4 == 2401 & gdename == "fretes"
replace PLZ4 = 2405 if PLZ4 == 2401 & gdename == "la chaux du millieu"
replace PLZ4 = 2314 if PLZ4 == 2401 & gdename == "la combe boudry"
replace PLZ4 = 2405 if PLZ4 == 2401 & gdename == "le cachot"
replace PLZ4 = 2405 if PLZ4 == 2401 & gdename == "le cachot/la chaux du milieu"
replace PLZ4 = 2416 if PLZ4 == 2401 & gdename == "les fretes"
replace PLZ4 = 3555 if PLZ4 == 6199 & gdename == "kroeschenbrunnen"
replace PLZ4 = 3555 if PLZ4 == 6199 & gdename == "kroeschenbrunnen/trub"
replace PLZ4 = 3421 if PLZ4 == 3327 & gdename == "lissach"
replace PLZ4 = 3421 if PLZ4 == 3327 & gdename == "llssach"
replace PLZ4 = 3421 if PLZ4 == 3327 & gdename == "rueti b. lyssach"
replace PLZ4 = 3421 if PLZ4 == 3327 & gdename == "rueti bel lyssach"
replace PLZ4 = 3421 if PLZ4 == 3327 & gdename == "rueti lyssach"
replace PLZ4 = 3421 if PLZ4 == 3327 & gdename == "ruetl lyssach"
replace PLZ4 = 3421 if PLZ4 == 3327 & gdename == "ruti bei lyssach"
replace PLZ4 = 1071 if PLZ4 == 1605 & gdename == "cherbres"
replace PLZ4 = 1071 if PLZ4 == 1605 & gdename == "cheschrez"
replace PLZ4 = 1071 if PLZ4 == 1605 & gdename == "chexbrea"
replace PLZ4 = 1071 if PLZ4 == 1605 & gdename == "chexbree"
replace PLZ4 = 1071 if PLZ4 == 1605 & gdename == "chexbres vd"
replace PLZ4 = 1071 if PLZ4 == 1605 & gdename == "cremieres"
replace PLZ4 = 1071 if PLZ4 == 1605 & gdename == "cremieres chexbres"
replace PLZ4 = 1071 if PLZ4 == 1605 & gdename == "cremieres sur chexbres"
replace PLZ4 = 1071 if PLZ4 == 1605 & gdename == "lignieres sur chexbres"
replace PLZ4 = 2071 if PLZ4 == 2018 & gdename == "perreux"
replace PLZ4 = 6074 if PLZ4 == 6075 & gdename == "grossteil"
replace PLZ4 = 6074 if PLZ4 == 6075 & gdename == "grossteil/giswil"
replace PLZ4 = 6074 if PLZ4 == 6075 & gdename == "grosstell"
replace PLZ4 = 6162 if PLZ4 == 6165 & gdename == "finsterwaid"
replace PLZ4 = 6162 if PLZ4 == 6165 & gdename == "finsterwald"
replace PLZ4 = 6162 if PLZ4 == 6165 & gdename == "finsterwald b. entlebuch"
replace PLZ4 = 6162 if PLZ4 == 6165 & gdename == "finsterwald bei entiebuch"
replace PLZ4 = 6162 if PLZ4 == 6165 & gdename == "finsterwald bei entlebuch"
replace PLZ4 = 6162 if PLZ4 == 6165 & gdename == "finsterwald bel entlebuch"
replace PLZ4 = 6162 if PLZ4 == 6165 & gdename == "finsterwald/ entlebuch"
replace PLZ4 = 6162 if PLZ4 == 6165 & gdename == "finsterwald/entlebuch"
replace PLZ4 = 7606 if PLZ4 == 7649 & gdename == "bonde"
replace PLZ4 = 7605 if PLZ4 == 7649 & gdename == "casaccia"
replace PLZ4 = 7605 if PLZ4 == 7649 & gdename == "casaccla"
replace PLZ4 = 7606 if PLZ4 == 7649 & gdename == "promontogno"
replace PLZ4 = 7516 if PLZ4 == 7649 & gdename == "stampa maloja"
replace PLZ4 = 1353 if PLZ4 == 1351 & gdename == "boffiens"
replace PLZ4 = 1356 if PLZ4 == 1351 & gdename == "la russilie/les clees"
replace PLZ4 = 1356 if PLZ4 == 1351 & gdename == "la russille vd"
replace PLZ4 = 1356 if PLZ4 == 1351 & gdename == "la russille/les clees"
replace PLZ4 = 1355 if PLZ4 == 1351 & gdename == "le"
replace PLZ4 = 1355 if PLZ4 == 1351 & gdename == "le vail loud l abergement"
replace PLZ4 = 1355 if PLZ4 == 1351 & gdename == "le vailloud l abergement"
replace PLZ4 = 1356 if PLZ4 == 1351 & gdename == "les ciees"
replace PLZ4 = 1357 if PLZ4 == 1351 & gdename == "ligneroile"
replace PLZ4 = 1357 if PLZ4 == 1351 & gdename == "lignerolie"
replace PLZ4 = 1354 if PLZ4 == 1351 & gdename == "orbe montcherand"
replace PLZ4 = 1439 if PLZ4 == 1351 & gdename == "rances vd"
replace PLZ4 = 1358 if PLZ4 == 1351 & gdename == "vaieyres sous rances"
replace PLZ4 = 1358 if PLZ4 == 1351 & gdename == "valeyers sous rances"
replace PLZ4 = 1358 if PLZ4 == 1351 & gdename == "valeyres sous pances"
replace PLZ4 = 1358 if PLZ4 == 1351 & gdename == "valeyres/rances"
replace PLZ4 = 2830 if PLZ4 == 2764 & gdename == "ccurrendlin"
replace PLZ4 = 2830 if PLZ4 == 2764 & gdename == "courrendiin"
replace PLZ4 = 2830 if PLZ4 == 2764 & gdename == "courrendiln"
replace PLZ4 = 2830 if PLZ4 == 2764 & gdename == "courrendli"
replace PLZ4 = 2830 if PLZ4 == 2764 & gdename == "courrendlln"
replace PLZ4 = 2830 if PLZ4 == 2764 & gdename == "courrer dlin"
replace PLZ4 = 8524 if PLZ4 == 8533 & gdename == "buch b. frauenfeld"
replace PLZ4 = 8524 if PLZ4 == 8533 & gdename == "buch bei frauenfeld"
replace PLZ4 = 8524 if PLZ4 == 8533 & gdename == "buch bel frauenfeld"
replace PLZ4 = 8608 if PLZ4 == 8533 & gdename == "wolfhausen"
replace PLZ4 = 3863 if PLZ4 == 3861 & gdename == "nessental"
replace PLZ4 = 4628 if PLZ4 == 4855 & gdename == "waltwil"
replace PLZ4 = 4628 if PLZ4 == 4855 & gdename == "woifwil"
replace PLZ4 = 4628 if PLZ4 == 4855 & gdename == "woifwil be"
replace PLZ4 = 4628 if PLZ4 == 4855 & gdename == "woitwil"
replace PLZ4 = 4628 if PLZ4 == 4855 & gdename == "wolf wil"
replace PLZ4 = 4628 if PLZ4 == 4855 & gdename == "wolfwii"
replace PLZ4 = 4628 if PLZ4 == 4855 & gdename == "wolfwil be"
replace PLZ4 = 4628 if PLZ4 == 4855 & gdename == "wolfwll"
replace PLZ4 = 4628 if PLZ4 == 4855 & gdename == "wolfwll be"
replace PLZ4 = 4628 if PLZ4 == 4855 & gdename == "woltwil"
replace PLZ4 = 6078 if PLZ4 == 6077 & gdename == "buerglen ow"
replace PLZ4 = 6078 if PLZ4 == 6077 & gdename == "buerglen ow/lungern"
replace PLZ4 = 6078 if PLZ4 == 6077 & gdename == "buergten ow"
replace PLZ4 = 6078 if PLZ4 == 6077 & gdename == "buergten ow/lungern"
replace PLZ4 = 6078 if PLZ4 == 6077 & gdename == "kaiserstuhl ow"
replace PLZ4 = 6078 if PLZ4 == 6077 & gdename == "kaiserstuhl ow/lungern"
replace PLZ4 = 6078 if PLZ4 == 6077 & gdename == "kalserstuhl"
replace PLZ4 = 6078 if PLZ4 == 6077 & gdename == "kalserstuhl ow"
replace PLZ4 = 6162 if PLZ4 == 6164 & gdename == "rengg"
replace PLZ4 = 6162 if PLZ4 == 6164 & gdename == "rengg entiebuch"
replace PLZ4 = 6162 if PLZ4 == 6164 & gdename == "rengg entlebuch"
replace PLZ4 = 6410 if PLZ4 == 6411 & gdename == "rigi"
replace PLZ4 = 6410 if PLZ4 == 6411 & gdename == "rigi staffel"
replace PLZ4 = 6410 if PLZ4 == 6411 & gdename == "rigi staffel/arth"
replace PLZ4 = 6410 if PLZ4 == 6411 & gdename == "rigl staffel"
replace PLZ4 = 9223 if PLZ4 == 9221 & gdename == "schweizershoiz"
replace PLZ4 = 9225 if PLZ4 == 9221 & gdename == "st. pelagiberg"
replace PLZ4 = 9225 if PLZ4 == 9221 & gdename == "st. pelagiberg gotteshaus"
replace PLZ4 = 9225 if PLZ4 == 9221 & gdename == "st. pelaglberg gotteshaus"
replace PLZ4 = 9225 if PLZ4 == 9221 & gdename == "walen gottshaus"
replace PLZ4 = 9225 if PLZ4 == 9221 & gdename == "wilen gottshaus"
replace PLZ4 = 9225 if PLZ4 == 9221 & gdename == "wlien gottshaus"
replace PLZ4 = 1936 if PLZ4 == 1935 & gdename == "verbier"
replace PLZ4 = 1936 if PLZ4 == 1935 & gdename == "verbier bagnes"
replace PLZ4 = 1936 if PLZ4 == 1935 & gdename == "verbler"
replace PLZ4 = 2904 if PLZ4 == 2901 & gdename == "bressancourt"
replace PLZ4 = 2947 if PLZ4 == 2901 & gdename == "charmoille ju"
replace PLZ4 = 2947 if PLZ4 == 2901 & gdename == "charmollle"
replace PLZ4 = 2908 if PLZ4 == 2901 & gdename == "grandfontalne"
replace PLZ4 = 2946 if PLZ4 == 2901 & gdename == "mlecourt"
replace PLZ4 = 2912 if PLZ4 == 2901 & gdename == "reciere"
replace PLZ4 = 5083 if PLZ4 == 4349 & gdename == "ittersthal"
replace PLZ4 = 5083 if PLZ4 == 4349 & gdename == "lttenthal"
replace PLZ4 = 5274 if PLZ4 == 4349 & gdename == "meltau"
replace PLZ4 = 5274 if PLZ4 == 4349 & gdename == "menthe"
replace PLZ4 = 5085 if PLZ4 == 4349 & gdename == "rheinsulz"
replace PLZ4 = 5085 if PLZ4 == 4349 & gdename == "sulz b. laufenburg"
replace PLZ4 = 5085 if PLZ4 == 4349 & gdename == "sulz bei laufenburg"
replace PLZ4 = 5085 if PLZ4 == 4349 & gdename == "sulz bel laufenburg"
replace PLZ4 = 5276 if PLZ4 == 4349 & gdename == "wii"
replace PLZ4 = 5276 if PLZ4 == 4349 & gdename == "wll ag"
replace PLZ4 = 7527 if PLZ4 == 7549 & gdename == "brail"
replace PLZ4 = 7526 if PLZ4 == 7549 & gdename == "cinuos chel"
replace PLZ4 = 7545 if PLZ4 == 7549 & gdename == "guards"
replace PLZ4 = 7522 if PLZ4 == 7549 & gdename == "la punt"
replace PLZ4 = 7522 if PLZ4 == 7549 & gdename == "la punt chamues ch*"
replace PLZ4 = 7522 if PLZ4 == 7549 & gdename == "la punt charnues ch"
replace PLZ4 = 7523 if PLZ4 == 7549 & gdename == "madulaln"
replace PLZ4 = 8564 if PLZ4 == 8563 & gdename == "engwilen"
replace PLZ4 = 8564 if PLZ4 == 8563 & gdename == "engwllen"
replace PLZ4 = 8564 if PLZ4 == 8563 & gdename == "sonterswil"
replace PLZ4 = 1000 if PLZ4 == 1075 & gdename == "chalet a gobet"
replace PLZ4 = 1000 if PLZ4 == 1075 & gdename == "chalet ae gobet"
replace PLZ4 = 1000 if PLZ4 == 1075 & gdename == "le chalet a gobet"
replace PLZ4 = 1000 if PLZ4 == 1075 & gdename == "le chalet ae gobet"
replace PLZ4 = 1000 if PLZ4 == 1075 & gdename == "lp chalet a gobet"
replace PLZ4 = 3473 if PLZ4 == 3399 & gdename == "aichenstorf"
replace PLZ4 = 3473 if PLZ4 == 3399 & gdename == "alchensdorf"
replace PLZ4 = 3473 if PLZ4 == 3399 & gdename == "alchenstorf be"
replace PLZ4 = 3476 if PLZ4 == 3399 & gdename == "oschwand"
replace PLZ4 = 3474 if PLZ4 == 3399 & gdename == "rueedisbach"
replace PLZ4 = 3377 if PLZ4 == 4706 & gdename == "wailiswil bei wangen"
replace PLZ4 = 3377 if PLZ4 == 4706 & gdename == "walliswal bei wangen"
replace PLZ4 = 3377 if PLZ4 == 4706 & gdename == "walliswal bel wangen"
replace PLZ4 = 3377 if PLZ4 == 4706 & gdename == "walliswil"
replace PLZ4 = 3377 if PLZ4 == 4706 & gdename == "walliswil b. wangen"
replace PLZ4 = 3377 if PLZ4 == 4706 & gdename == "walliswil bel wangen"
replace PLZ4 = 3377 if PLZ4 == 4706 & gdename == "walliswll bei wangen"
replace PLZ4 = 3377 if PLZ4 == 4706 & gdename == "walllawil bel wangen"
replace PLZ4 = 3377 if PLZ4 == 4706 & gdename == "walllswil bel wangen"
replace PLZ4 = 3377 if PLZ4 == 4706 & gdename == "walllswlil bei wangen"
replace PLZ4 = 3629 if PLZ4 == 3117 & gdename == "klesen"
replace PLZ4 = 3629 if PLZ4 == 3117 & gdename == "oppiigen"
replace PLZ4 = 3629 if PLZ4 == 3117 & gdename == "oppiingen"
replace PLZ4 = 3629 if PLZ4 == 3117 & gdename == "oppilgen"
replace PLZ4 = 3629 if PLZ4 == 3117 & gdename == "oppl gen"
replace PLZ4 = 3629 if PLZ4 == 3117 & gdename == "oppligen be"
replace PLZ4 = 3629 if PLZ4 == 3117 & gdename == "opplingen"
replace PLZ4 = 3629 if PLZ4 == 3117 & gdename == "oppllgen"
replace PLZ4 = 6122 if PLZ4 == 6124 & gdename == "twenenegg"
replace PLZ4 = 6122 if PLZ4 == 6124 & gdename == "twerenegg"
replace PLZ4 = 8136 if PLZ4 == 8138 & gdename == "gattlkon"
replace PLZ4 = 3661 if PLZ4 == 8138 & gdename == "uetendort"
replace PLZ4 = 8143 if PLZ4 == 8138 & gdename == "uetliberg"
replace PLZ4 = 8143 if PLZ4 == 8138 & gdename == "uettiberg"
replace PLZ4 = 8915 if PLZ4 == 8943 & gdename == "sihibrugg"
replace PLZ4 = 8915 if PLZ4 == 8943 & gdename == "sihlbrugg"
replace PLZ4 = 8915 if PLZ4 == 8943 & gdename == "sihlbrugg station"
replace PLZ4 = 8915 if PLZ4 == 8943 & gdename == "slhibrugg"
replace PLZ4 = 8915 if PLZ4 == 8943 & gdename == "slhlbrugg"
replace PLZ4 = 2744 if PLZ4 == 2741 & gdename == "balprahon"
replace PLZ4 = 2744 if PLZ4 == 2741 & gdename == "beiprahon"
replace PLZ4 = 2744 if PLZ4 == 2741 & gdename == "belpranon"
replace PLZ4 = 2743 if PLZ4 == 2741 & gdename == "escihert"
replace PLZ4 = 2742 if PLZ4 == 2741 & gdename == "perefitte"
replace PLZ4 = 2742 if PLZ4 == 2741 & gdename == "perrefltte"
replace PLZ4 = 2742 if PLZ4 == 2741 & gdename == "perretitte"
replace PLZ4 = 2742 if PLZ4 == 2741 & gdename == "peyrefitte"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "wdlflinswll"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woelfiinswil"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woelflanswil"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woelflinswii"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woelflinswl"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woelflinswll"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woelfllinswil"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woelfllnswil"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woelfllnswll"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woeltllnswil"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "wolflinswil"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woltiinswil"
replace PLZ4 = 5063 if PLZ4 == 5266 & gdename == "woltlinswil"
replace PLZ4 = 6122 if PLZ4 == 6111 & gdename == "fontannen bei wolhusen"
replace PLZ4 = 6122 if PLZ4 == 6111 & gdename == "fontannen bel wolhusen"
replace PLZ4 = 6114 if PLZ4 == 6111 & gdename == "steinhuserberg"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gassau"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gassau sg"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gaussa sg"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gcssau sg"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gessau (sg)"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "goaaau sg"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "goasau (sg)"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "goss au sg"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gossaeu (sg)"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gossau sa"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gossau sg e"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gossaue"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gossaus"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gossen sg"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gowssau"
replace PLZ4 = 9200 if PLZ4 == 9202 & gdename == "gussau sg"
replace PLZ4 = 9428 if PLZ4 == 9412 & gdename == "lachen ar"
replace PLZ4 = 9428 if PLZ4 == 9412 & gdename == "lachen walzenhausen"
replace PLZ4 = 9428 if PLZ4 == 9412 & gdename == "lachen/ar"
replace PLZ4 = 9428 if PLZ4 == 9412 & gdename == "lachen/walzenhausen"
replace PLZ4 = 1563 if PLZ4 == 1557 & gdename == "dampierre vd"
replace PLZ4 = 1563 if PLZ4 == 1557 & gdename == "domplerre vd"
replace PLZ4 = 1563 if PLZ4 == 1557 & gdename == "granges de dompierre"
replace PLZ4 = 1563 if PLZ4 == 1557 & gdename == "granges de domplerre"
replace PLZ4 = 1563 if PLZ4 == 1557 & gdename == "granges de=dompierre"
replace PLZ4 = 2353 if PLZ4 == 2727 & gdename == "cernievillers"
replace PLZ4 = 2353 if PLZ4 == 2727 & gdename == "cernlevlllers"
replace PLZ4 = 2353 if PLZ4 == 2727 & gdename == "lea pommerats"
replace PLZ4 = 2353 if PLZ4 == 2727 & gdename == "les pommerais"
replace PLZ4 = 2353 if PLZ4 == 2727 & gdename == "les pommerets"
replace PLZ4 = 2353 if PLZ4 == 2727 & gdename == "pommerats"
replace PLZ4 = 3375 if PLZ4 == 4555 & gdename == "inkwil roethenbach"
replace PLZ4 = 3375 if PLZ4 == 4555 & gdename == "inkwll"
replace PLZ4 = 3375 if PLZ4 == 4555 & gdename == "lnkwil"
replace PLZ4 = 3375 if PLZ4 == 4555 & gdename == "lnkwil (roethenbach)"
replace PLZ4 = 3375 if PLZ4 == 4555 & gdename == "lnkwll (roethenbach)"
replace PLZ4 = 6170 if PLZ4 == 6172 & gdename == "kiusstalden /schuepfheim"
replace PLZ4 = 6170 if PLZ4 == 6172 & gdename == "klusstaiden"
replace PLZ4 = 6170 if PLZ4 == 6172 & gdename == "klusstalden"
replace PLZ4 = 6170 if PLZ4 == 6172 & gdename == "klusstalden /schuepf heim"
replace PLZ4 = 6170 if PLZ4 == 6172 & gdename == "klusstalden /schuepfheim"
replace PLZ4 = 6170 if PLZ4 == 6172 & gdename == "ktusatalden"
replace PLZ4 = 8243 if PLZ4 == 8211 & gdename == "astdorf"
replace PLZ4 = 8239 if PLZ4 == 8211 & gdename == "dorflingen"
replace PLZ4 = 8234 if PLZ4 == 8211 & gdename == "steffen"
replace PLZ4 = 8234 if PLZ4 == 8211 & gdename == "steffen (sh)"
replace PLZ4 = 8234 if PLZ4 == 8211 & gdename == "steffen sh"
replace PLZ4 = 8234 if PLZ4 == 8211 & gdename == "stetten 8h"
replace PLZ4 = 9470 if PLZ4 == 9474 & gdename == "raefis"
replace PLZ4 = 9470 if PLZ4 == 9474 & gdename == "raefis buchs"
replace PLZ4 = 9470 if PLZ4 == 9474 & gdename == "raefls"
replace PLZ4 = 4717 if PLZ4 == 4711 & gdename == "ramiswil"
replace PLZ4 = 8187 if PLZ4 == 8433 & gdename == "walach"
replace PLZ4 = 8187 if PLZ4 == 8433 & gdename == "welach"
replace PLZ4 = 8187 if PLZ4 == 8433 & gdename == "welsch"
replace PLZ4 = 8815 if PLZ4 == 8811 & gdename == "horgenberg"
replace PLZ4 = 8815 if PLZ4 == 8811 & gdename == "horgerberg"
replace PLZ4 = 8825 if PLZ4 == 8821 & gdename == "hueften"
replace PLZ4 = 8824 if PLZ4 == 8821 & gdename == "schoenenherg"
replace PLZ4 = 8824 if PLZ4 == 8821 & gdename == "schonenberg"
replace PLZ4 = 8824 if PLZ4 == 8821 & gdename == "schonenberg (zh)"
replace PLZ4 = 8824 if PLZ4 == 8821 & gdename == "schã³nenberg (zh)"
replace PLZ4 = 8824 if PLZ4 == 8821 & gdename == "tanne schoenenberg"
replace PLZ4 = 8824 if PLZ4 == 8821 & gdename == "tanne waedenswil"
replace PLZ4 = 9514 if PLZ4 == 9516 & gdename == "hagenbuch am nollen"
replace PLZ4 = 9514 if PLZ4 == 9516 & gdename == "hagenbuch tg"
replace PLZ4 = 8577 if PLZ4 == 9516 & gdename == "toos"
replace PLZ4 = 8577 if PLZ4 == 9516 & gdename == "toos/schoenholzerswilen"
replace PLZ4 = 1885 if PLZ4 == 1855 & gdename == "chesieres ollon"
replace PLZ4 = 1885 if PLZ4 == 1855 & gdename == "chesieres olson"
replace PLZ4 = 1867 if PLZ4 == 1855 & gdename == "st triphon"
replace PLZ4 = 1867 if PLZ4 == 1855 & gdename == "st trlphon"
replace PLZ4 = 3671 if PLZ4 == 3526 & gdename == "brenziikofen"
replace PLZ4 = 3671 if PLZ4 == 3526 & gdename == "brenzlkofen"
replace PLZ4 = 3256 if PLZ4 == 3526 & gdename == "dieterswil"
replace PLZ4 = 3671 if PLZ4 == 3526 & gdename == "herbiigen"
replace PLZ4 = 3671 if PLZ4 == 3526 & gdename == "herbilgen"
replace PLZ4 = 3671 if PLZ4 == 3526 & gdename == "herbligen be"
replace PLZ4 = 3671 if PLZ4 == 3526 & gdename == "herbligen brenzikofen"
replace PLZ4 = 3671 if PLZ4 == 3526 & gdename == "herblingen"
replace PLZ4 = 3671 if PLZ4 == 3526 & gdename == "herbllgen"
replace PLZ4 = 8965 if PLZ4 == 8986 & gdename == "mutschellen"
replace PLZ4 = 5621 if PLZ4 == 8986 & gdename == "mutschellen/zufikon"
replace PLZ4 = 5082 if PLZ4 == 4336 & gdename == "kaisten ag"
replace PLZ4 = 5082 if PLZ4 == 4336 & gdename == "kalsten"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "farel"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "farel (lavaux)"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "forai (lavauz)"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "forei (lavaux)"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "forel laveaux"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "forer lavaux"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "foret (lavaux)"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "foret (lavauz)"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "foret (uvaux)"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "foret/lavaux"
replace PLZ4 = 1072 if PLZ4 == 1606 & gdename == "lavaux"
replace PLZ4 = 1656 if PLZ4 == 1655 & gdename == "im fang"
replace PLZ4 = 3419 if PLZ4 == 3411 & gdename == "biembach"
replace PLZ4 = 3419 if PLZ4 == 3411 & gdename == "blembach"
replace PLZ4 = 3417 if PLZ4 == 3411 & gdename == "fiueegsau"
replace PLZ4 = 3417 if PLZ4 == 3411 & gdename == "rueegaau"
replace PLZ4 = 3417 if PLZ4 == 3411 & gdename == "ruegsau"
replace PLZ4 = 8192 if PLZ4 == 8432 & gdename == "giattfelden/zweidien"
replace PLZ4 = 8192 if PLZ4 == 8432 & gdename == "glattfelden/zweidlen"
replace PLZ4 = 8192 if PLZ4 == 8432 & gdename == "glattfelden/zweldlen"
replace PLZ4 = 8192 if PLZ4 == 8432 & gdename == "glattfeldern/zweidlen"
replace PLZ4 = 6432 if PLZ4 == 8432 & gdename == "rickenbach b. sch"
replace PLZ4 = 8192 if PLZ4 == 8432 & gdename == "zweidien"
replace PLZ4 = 8192 if PLZ4 == 8432 & gdename == "zweidlen"
replace PLZ4 = 8192 if PLZ4 == 8432 & gdename == "zweldlen"
replace PLZ4 = 8484 if PLZ4 == 8485 & gdename == "theilingen"
replace PLZ4 = 8484 if PLZ4 == 8485 & gdename == "theilingen/weisslingen"
replace PLZ4 = 8484 if PLZ4 == 8485 & gdename == "thellingen"
replace PLZ4 = 8512 if PLZ4 == 8513 & gdename == "lustdorf"
replace PLZ4 = 8512 if PLZ4 == 8513 & gdename == "wetzikon tg"
replace PLZ4 = 8512 if PLZ4 == 8513 & gdename == "wetzlkon tg"
replace PLZ4 = 8706 if PLZ4 == 8705 & gdename == "feldmeilen"
replace PLZ4 = 2054 if PLZ4 == 2055 & gdename == "st martin ne"
replace PLZ4 = 2054 if PLZ4 == 2055 & gdename == "st martin ne/ chezard st martin"
replace PLZ4 = 2054 if PLZ4 == 2055 & gdename == "st. martin ne"
replace PLZ4 = 3437 if PLZ4 == 3431 & gdename == "randflueh"
replace PLZ4 = 3437 if PLZ4 == 3431 & gdename == "ranflueh"
replace PLZ4 = 3437 if PLZ4 == 3431 & gdename == "ranfluh"
replace PLZ4 = 3437 if PLZ4 == 3431 & gdename == "ranflãºh"
replace PLZ4 = 3433 if PLZ4 == 3431 & gdename == "schwanden i.e"
replace PLZ4 = 5325 if PLZ4 == 4353 & gdename == "lelbatadt"
replace PLZ4 = 5325 if PLZ4 == 4353 & gdename == "lelbstadt"
replace PLZ4 = 5053 if PLZ4 == 5052 & gdename == "wittwil"
replace PLZ4 = 5053 if PLZ4 == 5052 & gdename == "wittwil staffelbach"
replace PLZ4 = 5053 if PLZ4 == 5052 & gdename == "wittwilstaffelbach"
replace PLZ4 = 5053 if PLZ4 == 5052 & gdename == "wittwll"
replace PLZ4 = 6837 if PLZ4 == 6831 & gdename == "bruzeila"
replace PLZ4 = 6837 if PLZ4 == 6831 & gdename == "bruzelia"
replace PLZ4 = 6838 if PLZ4 == 6831 & gdename == "cabblo"
replace PLZ4 = 6837 if PLZ4 == 6831 & gdename == "canegglo"
replace PLZ4 = 6838 if PLZ4 == 6831 & gdename == "gabbio"
replace PLZ4 = 6838 if PLZ4 == 6831 & gdename == "mugglo"
replace PLZ4 = 6839 if PLZ4 == 6831 & gdename == "segno"
replace PLZ4 = 5330 if PLZ4 == 8437 & gdename == "bad zurzach"
replace PLZ4 = 5330 if PLZ4 == 8437 & gdename == "zuerzach"
replace PLZ4 = 5330 if PLZ4 == 8437 & gdename == "zur.ach"
replace PLZ4 = 5330 if PLZ4 == 8437 & gdename == "zurzach ag"
replace PLZ4 = 5330 if PLZ4 == 8437 & gdename == "zurzadh"
replace PLZ4 = 5330 if PLZ4 == 8437 & gdename == "zurzech"
replace PLZ4 = 3628 if PLZ4 == 3118 & gdename == "uttingen"
replace PLZ4 = 3628 if PLZ4 == 3118 & gdename == "uttlgen"
replace PLZ4 = 8514 if PLZ4 == 8531 & gdename == "amilkon"
replace PLZ4 = 8514 if PLZ4 == 8531 & gdename == "amllkon"
replace PLZ4 = 8512 if PLZ4 == 8531 & gdename == "aufhofen thundorf"
replace PLZ4 = 8514 if PLZ4 == 8531 & gdename == "leutmerken"
replace PLZ4 = 3145 if PLZ4 == 3 & gdename == "145 niederscherli koeniz"
replace PLZ4 = 4000 if PLZ4 == 3 & gdename == "basei"
replace PLZ4 = 3073 if PLZ4 == 3 & gdename == "guemligen"
replace PLZ4 = 1510 if PLZ4 == 1087 & gdename == "bressonnaz"
replace PLZ4 = 1510 if PLZ4 == 1087 & gdename == "bressonnaz/moudon"
replace PLZ4 = 1510 if PLZ4 == 1087 & gdename == "bressonnez"
replace PLZ4 = 6156 if PLZ4 == 6157 & gdename == "luthern bad"
replace PLZ4 = 6156 if PLZ4 == 6157 & gdename == "luthern bad/luthern"
replace PLZ4 = 1042 if PLZ4 == 1049 & gdename == "assena"
replace PLZ4 = 1042 if PLZ4 == 1049 & gdename == "bioiey orjulaz"
replace PLZ4 = 1042 if PLZ4 == 1049 & gdename == "bloley magnoux"
replace PLZ4 = 1042 if PLZ4 == 1049 & gdename == "bloley orjulaz"
replace PLZ4 = 1035 if PLZ4 == 1049 & gdename == "boumens"
replace PLZ4 = 1044 if PLZ4 == 1049 & gdename == "fey vd"
replace PLZ4 = 1036 if PLZ4 == 1049 & gdename == "suilens"
replace PLZ4 = 1036 if PLZ4 == 1049 & gdename == "suliens"
replace PLZ4 = 1372 if PLZ4 == 1399 & gdename == "bavols"
replace PLZ4 = 1372 if PLZ4 == 1399 & gdename == "bavons"
replace PLZ4 = 1372 if PLZ4 == 1399 & gdename == "bevels"
replace PLZ4 = 1374 if PLZ4 == 1399 & gdename == "corcelles chavornay"
replace PLZ4 = 1374 if PLZ4 == 1399 & gdename == "corcelles sur chavomay"
replace PLZ4 = 1374 if PLZ4 == 1399 & gdename == "corcelles/chavornay"
replace PLZ4 = 1376 if PLZ4 == 1399 & gdename == "goumoen la ville"
replace PLZ4 = 2854 if PLZ4 == 2862 & gdename == "berlincourt"
replace PLZ4 = 3076 if PLZ4 == 3135 & gdename == "walltenwll"
replace PLZ4 = 3076 if PLZ4 == 3135 & gdename == "watlenwil"
replace PLZ4 = 3076 if PLZ4 == 3135 & gdename == "wattenwii"
replace PLZ4 = 3076 if PLZ4 == 3135 & gdename == "wattenwil (wort))"
replace PLZ4 = 3076 if PLZ4 == 3135 & gdename == "wattenwll"
replace PLZ4 = 3076 if PLZ4 == 3135 & gdename == "wattenwlã¬l"
replace PLZ4 = 3076 if PLZ4 == 3135 & gdename == "wattwile"
replace PLZ4 = 3550 if PLZ4 == 3352 & gdename == "baerau i.e"
replace PLZ4 = 3472 if PLZ4 == 3352 & gdename == "ruinendingen"
replace PLZ4 = 3472 if PLZ4 == 3352 & gdename == "wynigen be"
replace PLZ4 = 3472 if PLZ4 == 3352 & gdename == "wynlgen"
replace PLZ4 = 3674 if PLZ4 == 3518 & gdename == "bleiken"
replace PLZ4 = 3674 if PLZ4 == 3518 & gdename == "bleiken b. oberdiessbach"
replace PLZ4 = 3674 if PLZ4 == 3518 & gdename == "bleiken bei oberdlessbach"
replace PLZ4 = 3674 if PLZ4 == 3518 & gdename == "bleiken bel oberdiessbach"
replace PLZ4 = 3674 if PLZ4 == 3518 & gdename == "bleiken bel oberdlessbach"
replace PLZ4 = 3674 if PLZ4 == 3518 & gdename == "blelken bei oberdiessbach"
replace PLZ4 = 3674 if PLZ4 == 3518 & gdename == "blelken bel oberdiessbach"
replace PLZ4 = 3615 if PLZ4 == 3519 & gdename == "wangelen b. oberdiessbach/buchholterberg"
replace PLZ4 = 3615 if PLZ4 == 3519 & gdename == "wangelen bei oberdiessbach"
replace PLZ4 = 3615 if PLZ4 == 3519 & gdename == "wangelen bei oberdiessbach/ buchholterberg"
replace PLZ4 = 3615 if PLZ4 == 3519 & gdename == "wangelen bei oberdiessbach/buchholterberg"
replace PLZ4 = 7605 if PLZ4 == 7949 & gdename == "casaccia"
replace PLZ4 = 3415 if PLZ4 == 3514 & gdename == "schaffhausen im emmental"
replace PLZ4 = 3415 if PLZ4 == 3514 & gdename == "schafhausen be"
replace PLZ4 = 3415 if PLZ4 == 3514 & gdename == "schafhausen i.e"
replace PLZ4 = 3415 if PLZ4 == 3514 & gdename == "schafhausen im emmental/hasle bei burgdorf"
replace PLZ4 = 3415 if PLZ4 == 3514 & gdename == "schafhausen im emmental/haste bel burgdorf"
replace PLZ4 = 3415 if PLZ4 == 3514 & gdename == "schafhausen untergomerkinden"
replace PLZ4 = 3415 if PLZ4 == 3514 & gdename == "schoefhausen be"
replace PLZ4 = 3415 if PLZ4 == 3514 & gdename == "schoenhausen be"
replace PLZ4 = 5085 if PLZ4 == 4338 & gdename == "rheinsulz"
replace PLZ4 = 4944 if PLZ4 == 4931 & gdename == "answil"
replace PLZ4 = 4944 if PLZ4 == 4931 & gdename == "auswll"
replace PLZ4 = 4932 if PLZ4 == 4931 & gdename == "gutenberg"
replace PLZ4 = 4932 if PLZ4 == 4931 & gdename == "gutenburg be"
replace PLZ4 = 4932 if PLZ4 == 4931 & gdename == "gutenburg/thunstetten"
replace PLZ4 = 4935 if PLZ4 == 4931 & gdename == "leimiswii"
replace PLZ4 = 4935 if PLZ4 == 4931 & gdename == "leimiswll"
replace PLZ4 = 4935 if PLZ4 == 4931 & gdename == "lelmiswil"
replace PLZ4 = 4942 if PLZ4 == 4931 & gdename == "walterswii (be)"
replace PLZ4 = 2950 if PLZ4 == 2892 & gdename == "courgena"
replace PLZ4 = 2950 if PLZ4 == 2892 & gdename == "courgenax"
replace PLZ4 = 2950 if PLZ4 == 2892 & gdename == "courgnay"
replace PLZ4 = 2950 if PLZ4 == 2892 & gdename == "courtemautruy"
replace PLZ4 = 4542 if PLZ4 == 4708 & gdename == "luterbach so"
replace PLZ4 = 4542 if PLZ4 == 4708 & gdename == "luterbach/so"
replace PLZ4 = 4542 if PLZ4 == 4708 & gdename == "luterbech"
replace PLZ4 = 8708 if PLZ4 == 4708 & gdename == "mannedorf"
replace PLZ4 = 6055 if PLZ4 == 6099 & gdename == "piatus kulm"
replace PLZ4 = 6055 if PLZ4 == 6099 & gdename == "pilatus kulm"
replace PLZ4 = 6055 if PLZ4 == 6099 & gdename == "platus kulm"
replace PLZ4 = 6055 if PLZ4 == 6099 & gdename == "pllatus kulm"
replace PLZ4 = 7082 if PLZ4 == 7999 & gdename == "vaz obervaz"
replace PLZ4 = 1134 if PLZ4 == 1 & gdename == ".134 vufflens le chaeteau"
replace PLZ4 = 1162 if PLZ4 == 1 & gdename == "162 st prex"
replace PLZ4 = 1212 if PLZ4 == 1 & gdename == "212 grand lancy"
replace PLZ4 = 1222 if PLZ4 == 1 & gdename == "222 vesenaz"
replace PLZ4 = 1255 if PLZ4 == 1 & gdename == "255 v eyrier"
replace PLZ4 = 1618 if PLZ4 == 1 & gdename == "618 chatel st denis"
replace PLZ4 = 1292 if PLZ4 == 1 & gdename == "avenue de tourney. 1292 chambesy"
replace PLZ4 = 3932 if PLZ4 == 1 & gdename == "visperterm"
replace PLZ4 = 3207 if PLZ4 == 3249 & gdename == "goiaten"
replace PLZ4 = 3207 if PLZ4 == 3249 & gdename == "gurlern"
replace PLZ4 = 3206 if PLZ4 == 3249 & gdename == "jerisberg ferenbaim"
replace PLZ4 = 3206 if PLZ4 == 3249 & gdename == "jerisberg ferenbalm"
replace PLZ4 = 3234 if PLZ4 == 3249 & gdename == "vineiz"
replace PLZ4 = 3234 if PLZ4 == 3249 & gdename == "vinetz"
replace PLZ4 = 3207 if PLZ4 == 3249 & gdename == "wilerottigen"
replace PLZ4 = 3267 if PLZ4 == 3258 & gdename == "baggwii seedorf be"
replace PLZ4 = 3267 if PLZ4 == 3258 & gdename == "baggwil"
replace PLZ4 = 3267 if PLZ4 == 3258 & gdename == "baggwil seedorf"
replace PLZ4 = 3267 if PLZ4 == 3258 & gdename == "baggwil seedorf be"
replace PLZ4 = 3267 if PLZ4 == 3258 & gdename == "baggwll seedorf be"
replace PLZ4 = 3267 if PLZ4 == 3258 & gdename == "frienisberg"
replace PLZ4 = 5076 if PLZ4 == 5254 & gdename == "boezen ag"
replace PLZ4 = 5076 if PLZ4 == 5254 & gdename == "boezen/"
replace PLZ4 = 5076 if PLZ4 == 5254 & gdename == "bozen"
replace PLZ4 = 5405 if PLZ4 == 5513 & gdename == "daettwil"
replace PLZ4 = 5405 if PLZ4 == 5513 & gdename == "daettwil (ag)"
replace PLZ4 = 5405 if PLZ4 == 5513 & gdename == "daettwil baden"
replace PLZ4 = 5405 if PLZ4 == 5513 & gdename == "daettwil/baden"
replace PLZ4 = 5405 if PLZ4 == 5513 & gdename == "dattwil"
replace PLZ4 = 5405 if PLZ4 == 5513 & gdename == "dattwil baden"
replace PLZ4 = 6313 if PLZ4 == 8313 & gdename == "edllbach"
replace PLZ4 = 8310 if PLZ4 == 8313 & gdename == "ottikon b. kemptthal"
replace PLZ4 = 8310 if PLZ4 == 8313 & gdename == "ottikon bei kempttal"
replace PLZ4 = 8310 if PLZ4 == 8313 & gdename == "ottikon bei kemptthal"
replace PLZ4 = 8310 if PLZ4 == 8313 & gdename == "ottlkon bel kemptthal"
replace PLZ4 = 3627 if PLZ4 == 3527 & gdename == "halmberg"
replace PLZ4 = 3627 if PLZ4 == 3527 & gdename == "heimberg be"
replace PLZ4 = 3627 if PLZ4 == 3527 & gdename == "helmberg"
replace PLZ4 = 3627 if PLZ4 == 3527 & gdename == "holmberg"
replace PLZ4 = 6474 if PLZ4 == 6499 & gdename == "intschi"
replace PLZ4 = 6482 if PLZ4 == 6499 & gdename == "intschi gurtnellen"
replace PLZ4 = 6474 if PLZ4 == 6499 & gdename == "intschl"
replace PLZ4 = 6482 if PLZ4 == 6499 & gdename == "intschl gurtnellen"
replace PLZ4 = 6424 if PLZ4 == 6499 & gdename == "lauen"
replace PLZ4 = 6474 if PLZ4 == 6499 & gdename == "lntschi"
replace PLZ4 = 8213 if PLZ4 == 7105 & gdename == "neunkiroh"
replace PLZ4 = 7107 if PLZ4 == 7105 & gdename == "safien platz"
replace PLZ4 = 7107 if PLZ4 == 7105 & gdename == "saften platz"
replace PLZ4 = 1608 if PLZ4 == 1672 & gdename == "chapelle sur oron"
replace PLZ4 = 1610 if PLZ4 == 1672 & gdename == "oran la ville"
replace PLZ4 = 1610 if PLZ4 == 1672 & gdename == "oron"
replace PLZ4 = 1610 if PLZ4 == 1672 & gdename == "oron ia ville"
replace PLZ4 = 1610 if PLZ4 == 1672 & gdename == "oron la viiie"
replace PLZ4 = 1610 if PLZ4 == 1672 & gdename == "own la ville"
replace PLZ4 = 3612 if PLZ4 == 3528 & gdename == "steffisburg station"
replace PLZ4 = 3612 if PLZ4 == 3528 & gdename == "stefflsburg"
replace PLZ4 = 3267 if PLZ4 == 6267 & gdename == "baggwil"
replace PLZ4 = 3267 if PLZ4 == 6267 & gdename == "baggwil/seedorf be"
replace PLZ4 = 6300 if PLZ4 == 6316 & gdename == "zugerberg"
replace PLZ4 = 6300 if PLZ4 == 6316 & gdename == "zugerberg/zug"
replace PLZ4 = 6383 if PLZ4 == 6385 & gdename == "niederrickenbach"
replace PLZ4 = 6383 if PLZ4 == 6385 & gdename == "niederrickenbach nw"
replace PLZ4 = 8570 if PLZ4 == 8571 & gdename == "weerswilen"
replace PLZ4 = 5406 if PLZ4 == 5508 & gdename == "ruetihof"
replace PLZ4 = 5406 if PLZ4 == 5508 & gdename == "ruetihof b. baden"
replace PLZ4 = 5406 if PLZ4 == 5508 & gdename == "ruetihof baden"
replace PLZ4 = 5406 if PLZ4 == 5508 & gdename == "ruetihof/baden"
replace PLZ4 = 5742 if PLZ4 == 5747 & gdename == "kolliken"
replace PLZ4 = 4665 if PLZ4 == 5747 & gdename == "kuengoidingen"
replace PLZ4 = 4665 if PLZ4 == 5747 & gdename == "kuengoldingen"
replace PLZ4 = 6203 if PLZ4 == 6230 & gdename == "sempach station"
replace PLZ4 = 6659 if PLZ4 == 6651 & gdename == "camedo/borgnone"
replace PLZ4 = 6655 if PLZ4 == 6651 & gdename == "rasa"
replace PLZ4 = 6827 if PLZ4 == 6923 & gdename == "brusino"
replace PLZ4 = 6827 if PLZ4 == 6923 & gdename == "brusino arsizlo"
replace PLZ4 = 6827 if PLZ4 == 6923 & gdename == "brusio arizio"
replace PLZ4 = 6827 if PLZ4 == 6923 & gdename == "brusio arsizio"
replace PLZ4 = 6827 if PLZ4 == 6923 & gdename == "bruslno arsizio"
replace PLZ4 = 6883 if PLZ4 == 8682 & gdename == "brusata"
replace PLZ4 = 6382 if PLZ4 == 8682 & gdename == "bueren/oberdorf nw"
replace PLZ4 = 9565 if PLZ4 == 9518 & gdename == "rothenhausen"
replace PLZ4 = 1884 if PLZ4 == 1857 & gdename == "huemoz"
replace PLZ4 = 1884 if PLZ4 == 1857 & gdename == "huemoz villars"
replace PLZ4 = 1884 if PLZ4 == 1857 & gdename == "huemoz/ollon"
replace PLZ4 = 3084 if PLZ4 == 3884 & gdename == "wabern"
replace PLZ4 = 2500 if PLZ4 == 5202 & gdename == "biel/bienne"
replace PLZ4 = 5070 if PLZ4 == 5262 & gdename == "frack"
replace PLZ4 = 5070 if PLZ4 == 5262 & gdename == "frick ag"
replace PLZ4 = 4074 if PLZ4 == 5262 & gdename == "oberwil ag"
replace PLZ4 = 4074 if PLZ4 == 5262 & gdename == "oberwll ag"
replace PLZ4 = 8545 if PLZ4 == 5635 & gdename == "rickenbach ag"
replace PLZ4 = 8545 if PLZ4 == 5635 & gdename == "rlckenbach"
replace PLZ4 = 7563 if PLZ4 == 7363 & gdename == "sammnaun"
replace PLZ4 = 8468 if PLZ4 == 8469 & gdename == "guntalingen"
replace PLZ4 = 1231 if PLZ4 == 2131 & gdename == "conches"
replace PLZ4 = 4571 if PLZ4 == 4575 & gdename == "ichertswil"
replace PLZ4 = 4571 if PLZ4 == 4575 & gdename == "ichertswll"
replace PLZ4 = 4571 if PLZ4 == 4575 & gdename == "lchertswil"
replace PLZ4 = 6343 if PLZ4 == 6334 & gdename == "rotkreuz"
replace PLZ4 = 6377 if PLZ4 == 6446 & gdename == "seeiisberg"
replace PLZ4 = 6377 if PLZ4 == 6446 & gdename == "seeilsberg"
replace PLZ4 = 6377 if PLZ4 == 6446 & gdename == "seellsberg"
replace PLZ4 = 6598 if PLZ4 == 6698 & gdename == "tenero"
replace PLZ4 = 3150 if PLZ4 == 8150 & gdename == "schwarzen burg wahiern"
replace PLZ4 = 3150 if PLZ4 == 8150 & gdename == "schwarzen burg wahlern"
replace PLZ4 = 3150 if PLZ4 == 8150 & gdename == "schwarzenburg wahlern"
replace PLZ4 = 8162 if PLZ4 == 8163 & gdename == "obersteinmaur"
replace PLZ4 = 8126 if PLZ4 == 8163 & gdename == "zumlkon"
replace PLZ4 = 6648 if PLZ4 == 8848 & gdename == "mlnuslo"
replace PLZ4 = 8840 if PLZ4 == 8848 & gdename == "trachslau"
replace PLZ4 = 8840 if PLZ4 == 8848 & gdename == "trachslau/einsiedeln"
replace PLZ4 = 8840 if PLZ4 == 8848 & gdename == "trachstau/einsiedeln"
replace PLZ4 = 1000 if PLZ4 == 1067 & gdename == "vers chez les blanc"
replace PLZ4 = 1820 if PLZ4 == 1843 & gdename == "grandchamp"
replace PLZ4 = 1820 if PLZ4 == 1843 & gdename == "veytaux chillon"
replace PLZ4 = 1820 if PLZ4 == 1843 & gdename == "veytaux/chillon"
replace PLZ4 = 2500 if PLZ4 == 2509 & gdename == "biel/bienne"
replace PLZ4 = 4955 if PLZ4 == 4918 & gdename == "gondiswii"
replace PLZ4 = 4955 if PLZ4 == 4918 & gdename == "gondiswll"
replace PLZ4 = 4955 if PLZ4 == 4918 & gdename == "gondlawil"
replace PLZ4 = 4955 if PLZ4 == 4918 & gdename == "gondlswil"
replace PLZ4 = 4955 if PLZ4 == 4918 & gdename == "gondlswll"
replace PLZ4 = 6260 if PLZ4 == 6266 & gdename == "reidermoos"
replace PLZ4 = 6260 if PLZ4 == 6266 & gdename == "reidermoos/beiden"
replace PLZ4 = 6260 if PLZ4 == 6266 & gdename == "reidermoos/reiden"
replace PLZ4 = 8266 if PLZ4 == 6266 & gdename == "steckbom"
replace PLZ4 = 7513 if PLZ4 == 7713 & gdename == "silvaplana/surlei"
replace PLZ4 = 7513 if PLZ4 == 7713 & gdename == "slivaplana/surlei"
replace PLZ4 = 5273 if PLZ4 == 8361 & gdename == "oberhofen zh"
replace PLZ4 = 5273 if PLZ4 == 8361 & gdename == "oberhufen zh"
replace PLZ4 = 7017 if PLZ4 == 7 & gdename == "017 flims dorf"
replace PLZ4 = 7077 if PLZ4 == 7 & gdename == "077 valbella/vaz/obervaz"
replace PLZ4 = 7270 if PLZ4 == 7 & gdename == "270 davos platz"
replace PLZ4 = 7235 if PLZ4 == 7 & gdename == "299 fi feris station"
replace PLZ4 = 1222 if PLZ4 == 7 & gdename == "ch.des epines. 1222 vesenaz"
replace PLZ4 = 9249 if PLZ4 == 249 & gdename == "algetshausen"
replace PLZ4 = 1514 if PLZ4 == 1511 & gdename == "bussy moudon"
replace PLZ4 = 2318 if PLZ4 == 2092 & gdename == "petits ponts"
replace PLZ4 = 2336 if PLZ4 == 2337 & gdename == "le boechet"
replace PLZ4 = 3097 if PLZ4 == 3079 & gdename == "liebefeid"
replace PLZ4 = 3097 if PLZ4 == 3079 & gdename == "liebefeld"
replace PLZ4 = 3510 if PLZ4 == 3511 & gdename == "konoifingen gruenegg"
replace PLZ4 = 3510 if PLZ4 == 3511 & gdename == "konolfingen gruenegg"
replace PLZ4 = 3963 if PLZ4 == 3980 & gdename == "corin montana"
replace PLZ4 = 3960 if PLZ4 == 3980 & gdename == "sierra"
replace PLZ4 = 3960 if PLZ4 == 3980 & gdename == "slerre"
replace PLZ4 = 5272 if PLZ4 == 4346 & gdename == "gansingen ag"
replace PLZ4 = 5272 if PLZ4 == 4346 & gdename == "gansirngen"
replace PLZ4 = 5272 if PLZ4 == 4346 & gdename == "gensingen"
replace PLZ4 = 5079 if PLZ4 == 5256 & gdename == "zechen"
replace PLZ4 = 5079 if PLZ4 == 5256 & gdename == "zelhen"
replace PLZ4 = 5626 if PLZ4 == 5449 & gdename == "staffeln"
replace PLZ4 = 5637 if PLZ4 == 5638 & gdename == "geitwil"
replace PLZ4 = 5637 if PLZ4 == 5638 & gdename == "winterschwil"
replace PLZ4 = 5643 if PLZ4 == 5648 & gdename == "alikon"
replace PLZ4 = 5643 if PLZ4 == 5648 & gdename == "allkon"
replace PLZ4 = 2072 if PLZ4 == 7072 & gdename == "st bialse"
replace PLZ4 = 2072 if PLZ4 == 7072 & gdename == "st blaise"
replace PLZ4 = 5467 if PLZ4 == 8435 & gdename == "fis4bach"
replace PLZ4 = 5467 if PLZ4 == 8435 & gdename == "fisibach/ag"
replace PLZ4 = 5467 if PLZ4 == 8435 & gdename == "fislbach"
replace PLZ4 = 5467 if PLZ4 == 8435 & gdename == "flalbach"
replace PLZ4 = 5467 if PLZ4 == 8435 & gdename == "flslbach"
replace PLZ4 = 5467 if PLZ4 == 8435 & gdename == "halbach"
replace PLZ4 = 8586 if PLZ4 == 8686 & gdename == "riedt erlen"
replace PLZ4 = 1957 if PLZ4 == 1917 & gdename == "anion"
replace PLZ4 = 1957 if PLZ4 == 1917 & gdename == "ardon vd"
replace PLZ4 = 1957 if PLZ4 == 1917 & gdename == "ardon vs"
replace PLZ4 = 1957 if PLZ4 == 1917 & gdename == "ardow"
replace PLZ4 = 2830 if PLZ4 == 2763 & gdename == "choindez"
replace PLZ4 = 2832 if PLZ4 == 2763 & gdename == "rebeuveller"
replace PLZ4 = 2885 if PLZ4 == 2851 & gdename == "epauviliers"
replace PLZ4 = 2885 if PLZ4 == 2851 & gdename == "epauvlllers"
replace PLZ4 = 2885 if PLZ4 == 2851 & gdename == "les rangiers"
replace PLZ4 = 3075 if PLZ4 == 3057 & gdename == "ruefenacht"
replace PLZ4 = 3075 if PLZ4 == 3057 & gdename == "ruefenacht be"
replace PLZ4 = 3326 if PLZ4 == 3316 & gdename == "hub bei krauchthai"
replace PLZ4 = 3326 if PLZ4 == 3316 & gdename == "hub bei krauchthal"
replace PLZ4 = 4629 if PLZ4 == 4854 & gdename == "fueenbach"
replace PLZ4 = 4629 if PLZ4 == 4854 & gdename == "fuienbach"
replace PLZ4 = 6020 if PLZ4 == 6620 & gdename == "emmenbruecke"
replace PLZ4 = 6436 if PLZ4 == 8438 & gdename == "muotathel"
replace PLZ4 = 5323 if PLZ4 == 8438 & gdename == "riethelm"
replace PLZ4 = 5323 if PLZ4 == 8438 & gdename == "rletheim"
replace PLZ4 = 5323 if PLZ4 == 8438 & gdename == "rlethelm"
replace PLZ4 = 8645 if PLZ4 == 8654 & gdename == "jona sg"
replace PLZ4 = 8854 if PLZ4 == 8654 & gdename == "siebnen"
replace PLZ4 = 8775 if PLZ4 == 8776 & gdename == "haerzingen"
replace PLZ4 = 8775 if PLZ4 == 8776 & gdename == "hatzingen"
replace PLZ4 = 9604 if PLZ4 == 9235 & gdename == "liitisburg"
replace PLZ4 = 9604 if PLZ4 == 9235 & gdename == "lirtisburg"
replace PLZ4 = 9604 if PLZ4 == 9235 & gdename == "luetlsburg"
replace PLZ4 = 9215 if PLZ4 == 9235 & gdename == "schoenenberg a.thur"
replace PLZ4 = 9315 if PLZ4 == 9375 & gdename == "neukiirch egnach"
replace PLZ4 = 9315 if PLZ4 == 9375 & gdename == "neukirch egnach"
replace PLZ4 = 8704 if PLZ4 == 87604 & gdename == "herrliberg zh"
replace PLZ4 = 8704 if PLZ4 == 87604 & gdename == "herrllberg zh"
replace PLZ4 = 1410 if PLZ4 == 1065 & gdename == "thlerrens"
replace PLZ4 = 1162 if PLZ4 == 1102 & gdename == "st prei"
replace PLZ4 = 1162 if PLZ4 == 1102 & gdename == "st prex"
replace PLZ4 = 1880 if PLZ4 == 1887 & gdename == "fenalet sur bex"
replace PLZ4 = 1955 if PLZ4 == 1915 & gdename == "chamosen"
replace PLZ4 = 1815 if PLZ4 == 1915 & gdename == "clarens"
replace PLZ4 = 1815 if PLZ4 == 1915 & gdename == "clarens montreux"
replace PLZ4 = 3098 if PLZ4 == 3093 & gdename == "koeniz/kschliern"
replace PLZ4 = 3098 if PLZ4 == 3093 & gdename == "koeniz/schiiern"
replace PLZ4 = 3098 if PLZ4 == 3093 & gdename == "koeniz/schliern"
replace PLZ4 = 3663 if PLZ4 == 3137 & gdename == "gurzeien"
replace PLZ4 = 3298 if PLZ4 == 3290 & gdename == "oberwil bueren"
replace PLZ4 = 3298 if PLZ4 == 3290 & gdename == "oberwll bueren"
replace PLZ4 = 4114 if PLZ4 == 4149 & gdename == "hofsteffen"
replace PLZ4 = 4114 if PLZ4 == 4149 & gdename == "hofsteffen so"
replace PLZ4 = 4116 if PLZ4 == 4149 & gdename == "mariastein"
replace PLZ4 = 4116 if PLZ4 == 4149 & gdename == "meterlen"
replace PLZ4 = 5326 if PLZ4 == 4352 & gdename == "schwaderioch"
replace PLZ4 = 5064 if PLZ4 == 5265 & gdename == "wi .tnau"
replace PLZ4 = 5064 if PLZ4 == 5265 & gdename == "wlttnau"
replace PLZ4 = 6331 if PLZ4 == 6339 & gdename == "huenenberg cham"
replace PLZ4 = 6744 if PLZ4 == 6799 & gdename == "personlco"
replace PLZ4 = 6772 if PLZ4 == 6799 & gdename == "prata leventina"
replace PLZ4 = 6772 if PLZ4 == 6799 & gdename == "prato levantina"
replace PLZ4 = 7250 if PLZ4 == 7251 & gdename == "monbiei"
replace PLZ4 = 7250 if PLZ4 == 7251 & gdename == "monbiel"
replace PLZ4 = 7017 if PLZ4 == 7817 & gdename == "flims dorf"
replace PLZ4 = 8135 if PLZ4 == 8137 & gdename == "langnau a. albas"
replace PLZ4 = 8135 if PLZ4 == 8137 & gdename == "langnau a. albis"
replace PLZ4 = 8280 if PLZ4 == 8220 & gdename == "kreuzlangen"
replace PLZ4 = 8280 if PLZ4 == 8220 & gdename == "krezlingen"
replace PLZ4 = 8646 if PLZ4 == 8446 & gdename == "wagen"
replace PLZ4 = 9545 if PLZ4 == 9544 & gdename == "rosental"
replace PLZ4 = 9545 if PLZ4 == 9544 & gdename == "rosental waengi"
replace PLZ4 = 2300 if PLZ4 == 2320 & gdename == "oa chaux de fonds"
replace PLZ4 = 3673 if PLZ4 == 3517 & gdename == "linden b. oberdiessbach"
replace PLZ4 = 3673 if PLZ4 == 3517 & gdename == "linden be"
replace PLZ4 = 3673 if PLZ4 == 3517 & gdename == "linden bei oberdiessbach"
replace PLZ4 = 3421 if PLZ4 == 3630 & gdename == "rueti lyssach"
replace PLZ4 = 5080 if PLZ4 == 4335 & gdename == "laufen ourg"
replace PLZ4 = 5080 if PLZ4 == 4335 & gdename == "laufenbur9"
replace PLZ4 = 5080 if PLZ4 == 4335 & gdename == "laufer..burg"
replace PLZ4 = 5080 if PLZ4 == 4335 & gdename == "lautenburg"
replace PLZ4 = 5277 if PLZ4 == 4348 & gdename == "hottwil e wag wernli ag, granichen del"
replace PLZ4 = 5277 if PLZ4 == 4348 & gdename == "hottwll"
replace PLZ4 = 5600 if PLZ4 == 5460 & gdename == "lenzburg stauffen"
replace PLZ4 = 5642 if PLZ4 == 5641 & gdename == "muehiau unterhuenenberg"
replace PLZ4 = 5642 if PLZ4 == 5641 & gdename == "muehlau unterhuenenberg"
replace PLZ4 = 6808 if PLZ4 == 6847 & gdename == "taverne torricella"
replace PLZ4 = 6989 if PLZ4 == 6961 & gdename == "purasca"
replace PLZ4 = 7554 if PLZ4 == 7555 & gdename == "post crusch"
replace PLZ4 = 9450 if PLZ4 == 9439 & gdename == "luechingen"
replace PLZ4 = 9428 if PLZ4 == 9482 & gdename == "waizenhausen"
replace PLZ4 = 9428 if PLZ4 == 9482 & gdename == "walzen hausen"
replace PLZ4 = 1254 if PLZ4 == 4 & gdename == "ch. de la monaisse. 1254 jussy"
replace PLZ4 = 6414 if PLZ4 == 4 & gdename == "oberarth/arth"
replace PLZ4 = 1202 if PLZ4 == 19 & gdename == "petit saconnex"
replace PLZ4 = 1264 if PLZ4 == 1364 & gdename == "st cergue"
replace PLZ4 = 4144 if PLZ4 == 1444 & gdename == "ariesheim"
replace PLZ4 = 4144 if PLZ4 == 1444 & gdename == "arlesheirn"
replace PLZ4 = 4144 if PLZ4 == 1444 & gdename == "arleshelm"
replace PLZ4 = 2000 if PLZ4 == 2200 & gdename == "neuchaetel"
replace PLZ4 = 2000 if PLZ4 == 2200 & gdename == "neuchatei"
replace PLZ4 = 2300 if PLZ4 == 2309 & gdename == "les bulles"
replace PLZ4 = 2520 if PLZ4 == 2529 & gdename == "la neuveville be"
replace PLZ4 = 2363 if PLZ4 == 2728 & gdename == "cernlevlllers"
replace PLZ4 = 2353 if PLZ4 == 2728 & gdename == "goumots"
replace PLZ4 = 3124 if PLZ4 == 3199 & gdename == "beipberg"
replace PLZ4 = 3126 if PLZ4 == 3199 & gdename == "geiterfingen"
replace PLZ4 = 3126 if PLZ4 == 3199 & gdename == "geitertingen"
replace PLZ4 = 3262 if PLZ4 == 3221 & gdename == "suberg"
replace PLZ4 = 3250 if PLZ4 == 3260 & gdename == "lyse"
replace PLZ4 = 3960 if PLZ4 == 3260 & gdename == "slerre"
replace PLZ4 = 3110 if PLZ4 == 3310 & gdename == "muensingen be"
replace PLZ4 = 3465 if PLZ4 == 3458 & gdename == "durrenroth"
replace PLZ4 = 3507 if PLZ4 == 3570 & gdename == "bigien"
replace PLZ4 = 3507 if PLZ4 == 3570 & gdename == "blgien"
replace PLZ4 = 3858 if PLZ4 == 3585 & gdename == "hofstetten"
replace PLZ4 = 4938 if PLZ4 == 4038 & gdename == "rohrbach b. huttwil"
replace PLZ4 = 4310 if PLZ4 == 4130 & gdename == "rhelnfelden"
replace PLZ4 = 5246 if PLZ4 == 5118 & gdename == "scherz ag"
replace PLZ4 = 5077 if PLZ4 == 5255 & gdename == "elf ingen"
replace PLZ4 = 5077 if PLZ4 == 5255 & gdename == "elflangen"
replace PLZ4 = 5453 if PLZ4 == 5433 & gdename == "busslingen"
replace PLZ4 = 6600 if PLZ4 == 6607 & gdename == "solduno"
replace PLZ4 = 3962 if PLZ4 == 7229 & gdename == "schiere"
replace PLZ4 = 3962 if PLZ4 == 7229 & gdename == "schlers"
replace PLZ4 = 3084 if PLZ4 == 8084 & gdename == "wabern"
replace PLZ4 = 8344 if PLZ4 == 8346 & gdename == "neuthai"
replace PLZ4 = 8344 if PLZ4 == 8346 & gdename == "neuthal"
replace PLZ4 = 8330 if PLZ4 == 8380 & gdename == "pfaeffikon zh"
replace PLZ4 = 8707 if PLZ4 == 8701 & gdename == "uetikon a. see"
replace PLZ4 = 8707 if PLZ4 == 8701 & gdename == "uetlkon a. see"
replace PLZ4 = 8926 if PLZ4 == 8923 & gdename == "kappel a. albis"
replace PLZ4 = 8957 if PLZ4 == 8958 & gdename == "spreltenbach"
replace PLZ4 = 9304 if PLZ4 == 9394 & gdename == "bernhardzeii"
replace PLZ4 = 9304 if PLZ4 == 9394 & gdename == "bernhardzell"
replace PLZ4 = 9500 if PLZ4 == 9600 & gdename == "wii"
replace PLZ4 = 9500 if PLZ4 == 9600 & gdename == "wii sg"
replace PLZ4 = 9500 if PLZ4 == 9600 & gdename == "wit sg"
replace PLZ4 = 9642 if PLZ4 == 9644 & gdename == "wintersberg"
replace PLZ4 = 9642 if PLZ4 == 9644 & gdename == "wintersberg sg"
replace PLZ4 = 8362 if PLZ4 == 25587 & gdename == "balterswil tg"
replace PLZ4 = 8127 if PLZ4 == 81127 & gdename == "forch zh"
replace PLZ4 = 1245 if PLZ4 == 2 & gdename == "vesenaz"
replace PLZ4 = 1117 if PLZ4 == 11 & gdename == "i 1 grancy"
replace PLZ4 = 1920 if PLZ4 == 11 & gdename == "rue de l hopital. 1920 martigny"
replace PLZ4 = 1802 if PLZ4 == 12 & gdename == "rte. de lavaux. 1802 corseaux"
replace PLZ4 = 4124 if PLZ4 == 24 & gdename == "schoenenbuch bl"
replace PLZ4 = 3044 if PLZ4 == 304 & gdename == "lnnerberg"
replace PLZ4 = 3043 if PLZ4 == 304 & gdename == "uettligen"
replace PLZ4 = 1233 if PLZ4 == 1388 & gdename == "sezenove"
replace PLZ4 = 1422 if PLZ4 == 1392 & gdename == "3randson corcelettes"
replace PLZ4 = 1422 if PLZ4 == 1392 & gdename == "grandson corceiettes"
replace PLZ4 = 1803 if PLZ4 == 1810 & gdename == "le mont pellerin"
replace PLZ4 = 4104 if PLZ4 == 2104 & gdename == "oberwal"
replace PLZ4 = 4104 if PLZ4 == 2104 & gdename == "oberwii"
replace PLZ4 = 2735 if PLZ4 == 2135 & gdename == "malleray bevilard"
replace PLZ4 = 2564 if PLZ4 == 2156 & gdename == "4ellmund"
replace PLZ4 = 2564 if PLZ4 == 2156 & gdename == "belimund"
replace PLZ4 = 2718 if PLZ4 == 2719 & gdename == "fornet dessus"
replace PLZ4 = 2718 if PLZ4 == 2719 & gdename == "fornet dessus/lajoux"
replace PLZ4 = 2353 if PLZ4 == 2749 & gdename == "pommerats"
replace PLZ4 = 2733 if PLZ4 == 2749 & gdename == "poutenet"
replace PLZ4 = 7272 if PLZ4 == 2772 & gdename == "clavadel"
replace PLZ4 = 3979 if PLZ4 == 2941 & gdename == "groene"
replace PLZ4 = 3000 if PLZ4 == 3021 & gdename == "berne"
replace PLZ4 = 3065 if PLZ4 == 3064 & gdename == "roerswil bolligen"
replace PLZ4 = 3084 if PLZ4 == 3064 & gdename == "wabern"
replace PLZ4 = 3661 if PLZ4 == 3138 & gdename == "iietendorf"
replace PLZ4 = 3661 if PLZ4 == 3138 & gdename == "uetendort"
replace PLZ4 = 3184 if PLZ4 == 3164 & gdename == "wuennewil"
replace PLZ4 = 3365 if PLZ4 == 3265 & gdename == "grasswil"
replace PLZ4 = 4414 if PLZ4 == 4141 & gdename == "fuelllnsdorf"
replace PLZ4 = 4414 if PLZ4 == 4141 & gdename == "fullinsdorf"
replace PLZ4 = 5276 if PLZ4 == 4347 & gdename == "wi ag"
replace PLZ4 = 5276 if PLZ4 == 4347 & gdename == "wii ag"
replace PLZ4 = 4452 if PLZ4 == 4403 & gdename == "itigen"
replace PLZ4 = 4452 if PLZ4 == 4403 & gdename == "ltingen"
replace PLZ4 = 4632 if PLZ4 == 4832 & gdename == "ti imbach"
replace PLZ4 = 4632 if PLZ4 == 4832 & gdename == "trlmbach"
replace PLZ4 = 5075 if PLZ4 == 5257 & gdename == "homussen"
replace PLZ4 = 5062 if PLZ4 == 5267 & gdename == "oberhot"
replace PLZ4 = 5273 if PLZ4 == 5267 & gdename == "oherhofen"
replace PLZ4 = 5074 if PLZ4 == 5286 & gdename == "eikerl"
replace PLZ4 = 5074 if PLZ4 == 5286 & gdename == "elken"
replace PLZ4 = 4813 if PLZ4 == 5743 & gdename == "uekheim"
replace PLZ4 = 4813 if PLZ4 == 5743 & gdename == "uerkhelm"
replace PLZ4 = 5630 if PLZ4 == 5830 & gdename == "murl"
replace PLZ4 = 5630 if PLZ4 == 5830 & gdename == "murl (ag)"
replace PLZ4 = 6952 if PLZ4 == 5952 & gdename == "canobblo"
replace PLZ4 = 6600 if PLZ4 == 6001 & gdename == "locarno/muralto"
replace PLZ4 = 8636 if PLZ4 == 6134 & gdename == "huebeli"
replace PLZ4 = 8636 if PLZ4 == 6134 & gdename == "huebell"
replace PLZ4 = 6166 if PLZ4 == 6168 & gdename == "heiligkreuz lu"
replace PLZ4 = 6402 if PLZ4 == 6492 & gdename == "merlischachen"
replace PLZ4 = 6454 if PLZ4 == 6554 & gdename == "floeleri ur"
replace PLZ4 = 6454 if PLZ4 == 6554 & gdename == "flueeien ur"
replace PLZ4 = 9620 if PLZ4 == 6920 & gdename == "lichtenstelg"
replace PLZ4 = 9620 if PLZ4 == 6920 & gdename == "liechtenstein"
replace PLZ4 = 8127 if PLZ4 == 8129 & gdename == "forch"
replace PLZ4 = 6170 if PLZ4 == 8170 & gdename == "schoepfhelm"
replace PLZ4 = 6170 if PLZ4 == 8170 & gdename == "schuepfhelm"
replace PLZ4 = 8310 if PLZ4 == 8319 & gdename == "kemptal lindau"
replace PLZ4 = 8566 if PLZ4 == 8516 & gdename == "neuwilen/lippoldswilen"
replace PLZ4 = 8566 if PLZ4 == 8516 & gdename == "neuwilen/llppoldswilen"
replace PLZ4 = 8121 if PLZ4 == 8721 & gdename == "benglen faellanden"
replace PLZ4 = 8739 if PLZ4 == 8760 & gdename == "mellen"
replace PLZ4 = 8739 if PLZ4 == 8760 & gdename == "rieden gl"
replace PLZ4 = 8840 if PLZ4 == 8814 & gdename == "einsiedein"
replace PLZ4 = 8840 if PLZ4 == 8814 & gdename == "lamons"
replace PLZ4 = 8840 if PLZ4 == 8838 & gdename == "bennau"
replace PLZ4 = 8038 if PLZ4 == 8838 & gdename == "zurich"
replace PLZ4 = 6900 if PLZ4 == 8900 & gdename == "lugano/massagno"
replace PLZ4 = 6900 if PLZ4 == 8900 & gdename == "paradlso"
replace PLZ4 = 8926 if PLZ4 == 8928 & gdename == "kappel am albfis"
replace PLZ4 = 8926 if PLZ4 == 8928 & gdename == "uerzllkon"
replace PLZ4 = 6981 if PLZ4 == 8981 & gdename == "crogllo"
replace PLZ4 = 8964 if PLZ4 == 8984 & gdename == "rudolfstetten"
replace PLZ4 = 8890 if PLZ4 == 8997 & gdename == "flumserberg"
replace PLZ4 = 9107 if PLZ4 == 9106 & gdename == "zuerchersmuehle"
replace PLZ4 = 9444 if PLZ4 == 9440 & gdename == "dlepoldsau"
replace PLZ4 = 7165 if PLZ4 == 5 & gdename == "brell/brigels"
replace PLZ4 = 1246 if PLZ4 == 22 & gdename == "ch. des usses. 1246 corsier"
replace PLZ4 = 1258 if PLZ4 == 25 & gdename == "ch du village. 125e perly"
replace PLZ4 = 8127 if PLZ4 == 27 & gdename == "forch"
replace PLZ4 = 1201 if PLZ4 == 48 & gdename == "r.de lausanne41201 geneve"
replace PLZ4 = 9053 if PLZ4 == 52 & gdename == "hiederteufen"
replace PLZ4 = 6962 if PLZ4 == 62 & gdename == "viganello ti"
replace PLZ4 = 1214 if PLZ4 == 65 & gdename == "rte de peney. 1214 vernier"
replace PLZ4 = 8121 if PLZ4 == 121 & gdename == "benglen"
replace PLZ4 = 1245 if PLZ4 == 122 & gdename == "vesenaz"
replace PLZ4 = 1261 if PLZ4 == 126 & gdename == "1261 burtigriy vd"
replace PLZ4 = 1257 if PLZ4 == 127 & gdename == "croix de rozon"
replace PLZ4 = 8152 if PLZ4 == 152 & gdename == "glattbrugg"
replace PLZ4 = 1752 if PLZ4 == 175 & gdename == "l/illars sur glane"
replace PLZ4 = 1892 if PLZ4 == 189 & gdename == "lavey village"
replace PLZ4 = 2520 if PLZ4 == 214 & gdename == "chavannes (be)"
replace PLZ4 = 1218 if PLZ4 == 218 & gdename == "grand saconnex"
replace PLZ4 = 3074 if PLZ4 == 274 & gdename == "muri b. bern"
replace PLZ4 = 3322 if PLZ4 == 322 & gdename == "urtenen schoenbuehl"
replace PLZ4 = 4334 if PLZ4 == 334 & gdename == "sisseln ag"
replace PLZ4 = 8354 if PLZ4 == 353 & gdename == "wenzikon"
replace PLZ4 = 8560 if PLZ4 == 562 & gdename == "maerstetten dorf"
replace PLZ4 = 8570 if PLZ4 == 570 & gdename == "welnfelden"
replace PLZ4 = 8702 if PLZ4 == 702 & gdename == "zotlikssn"
replace PLZ4 = 8704 if PLZ4 == 704 & gdename == "herrtlberg"
replace PLZ4 = 7250 if PLZ4 == 720 & gdename == "klosters"
replace PLZ4 = 8102 if PLZ4 == 810 & gdename == "l.)berengstringen"
replace PLZ4 = 8632 if PLZ4 == 830 & gdename == "tann rueti"
replace PLZ4 = 8133 if PLZ4 == 833 & gdename == "esslingen"
replace PLZ4 = 8555 if PLZ4 == 855 & gdename == "muellheim dorf"
replace PLZ4 = 4853 if PLZ4 == 857 & gdename == "riken/ murgenthal"
replace PLZ4 = 1400 if PLZ4 == 1100 & gdename == "yverdon"
replace PLZ4 = 1245 if PLZ4 == 1221 & gdename == "vesenaz/"
replace PLZ4 = 1023 if PLZ4 == 1238 & gdename == "crissler"
replace PLZ4 = 1025 if PLZ4 == 1250 & gdename == "st.sulpice"
replace PLZ4 = 1212 if PLZ4 == 1282 & gdename == "grand lancy"
replace PLZ4 = 1832 if PLZ4 == 1332 & gdename == "chamby"
replace PLZ4 = 1820 if PLZ4 == 1620 & gdename == "le chaetelard (vd)"
replace PLZ4 = 1815 if PLZ4 == 1850 & gdename == "clarens"
replace PLZ4 = 2954 if PLZ4 == 1851 & gdename == "les rangiers"
replace PLZ4 = 1950 if PLZ4 == 1952 & gdename == "slon"
replace PLZ4 = 3628 if PLZ4 == 2118 & gdename == "uttlgen"
replace PLZ4 = 2406 if PLZ4 == 2128 & gdename == "le brouillet"
replace PLZ4 = 2300 if PLZ4 == 2230 & gdename == "4a chaux de fonds"
replace PLZ4 = 2557 if PLZ4 == 2357 & gdename == "studen be"
replace PLZ4 = 4204 if PLZ4 == 2449 & gdename == "hammelried"
replace PLZ4 = 2500 if PLZ4 == 2507 & gdename == "biel/bienne"
replace PLZ4 = 2572 if PLZ4 == 2527 & gdename == "moeriger"
replace PLZ4 = 2557 if PLZ4 == 2567 & gdename == "studen b. bruegg"
replace PLZ4 = 2560 if PLZ4 == 2580 & gdename == "nldau"
replace PLZ4 = 2533 if PLZ4 == 2633 & gdename == "evllard"
replace PLZ4 = 2545 if PLZ4 == 2645 & gdename == "setzach"
replace PLZ4 = 7302 if PLZ4 == 2707 & gdename == "landquart fabriken"
replace PLZ4 = 2832 if PLZ4 == 2765 & gdename == "rebeuveller"
replace PLZ4 = 3825 if PLZ4 == 2835 & gdename == "murren"
replace PLZ4 = 3930 if PLZ4 == 2930 & gdename == "viege"
replace PLZ4 = 3963 if PLZ4 == 2963 & gdename == "crans sur sierre"
replace PLZ4 = 2072 if PLZ4 == 2972 & gdename == "st. blaise"
replace PLZ4 = 3065 if PLZ4 == 3060 & gdename == "ferenberg bolligen"
replace PLZ4 = 6332 if PLZ4 == 3131 & gdename == "hagendorn"
replace PLZ4 = 8142 if PLZ4 == 3142 & gdename == "waldegg b. uitikon"
replace PLZ4 = 3178 if PLZ4 == 3187 & gdename == "besingen"
replace PLZ4 = 3664 if PLZ4 == 3194 & gdename == "burgisteln"
replace PLZ4 = 3270 if PLZ4 == 3288 & gdename == "lobsagen/aarberg"
replace PLZ4 = 3303 if PLZ4 == 3304 & gdename == "zuzwll be"
replace PLZ4 = 3322 if PLZ4 == 3320 & gdename == "schoenbuehl urtenen"
replace PLZ4 = 8353 if PLZ4 == 3353 & gdename == "elgq"
replace PLZ4 = 4556 if PLZ4 == 3364 & gdename == "stelnhof"
replace PLZ4 = 3362 if PLZ4 == 3382 & gdename == "nlederoenz"
replace PLZ4 = 3097 if PLZ4 == 3397 & gdename == "liebefeld"
replace PLZ4 = 8442 if PLZ4 == 3442 & gdename == "hettlingerr"
replace PLZ4 = 8483 if PLZ4 == 3483 & gdename == "kollbrunn"
replace PLZ4 = 3773 if PLZ4 == 3500 & gdename == "matten"
replace PLZ4 = 3625 if PLZ4 == 3525 & gdename == "heiligenschwendl"
replace PLZ4 = 8546 if PLZ4 == 3546 & gdename == "lslikon"
replace PLZ4 = 3654 if PLZ4 == 3554 & gdename == "gunten"
replace PLZ4 = 3754 if PLZ4 == 3573 & gdename == "oey"
replace PLZ4 = 8580 if PLZ4 == 3580 & gdename == "niedersommeri"
replace PLZ4 = 8645 if PLZ4 == 3640 & gdename == "rapperswii sg"
replace PLZ4 = 3782 if PLZ4 == 3732 & gdename == "lauenen bei gstaad"
replace PLZ4 = 3780 if PLZ4 == 3786 & gdename == "gstaad"
replace PLZ4 = 3778 if PLZ4 == 3788 & gdename == "schoenried"
replace PLZ4 = 3608 if PLZ4 == 3808 & gdename == "thun allmendingen"
replace PLZ4 = 3617 if PLZ4 == 3817 & gdename == "fahrni b. thun"
replace PLZ4 = 3074 if PLZ4 == 3874 & gdename == "muri b. bern"
replace PLZ4 = 3000 if PLZ4 == 3897 & gdename == "bem"
replace PLZ4 = 4123 if PLZ4 == 4163 & gdename == "neu allschwil"
replace PLZ4 = 4233 if PLZ4 == 4198 & gdename == "neffingen"
replace PLZ4 = 5274 if PLZ4 == 4344 & gdename == "meltau"
replace PLZ4 = 6373 if PLZ4 == 4366 & gdename == "buergenstock"
replace PLZ4 = 4303 if PLZ4 == 4394 & gdename == "kalseraugst"
replace PLZ4 = 4434 if PLZ4 == 4454 & gdename == "holstein"
replace PLZ4 = 4448 if PLZ4 == 4488 & gdename == "laufelfingen"
replace PLZ4 = 4588 if PLZ4 == 4589 & gdename == "oberramsem"
replace PLZ4 = 4632 if PLZ4 == 4631 & gdename == "trlmbach"
replace PLZ4 = 4633 if PLZ4 == 4695 & gdename == "hauensteln"
replace PLZ4 = 4616 if PLZ4 == 4816 & gdename == "kappet so"
replace PLZ4 = 4625 if PLZ4 == 4825 & gdename == "oberbuchstten"
replace PLZ4 = 5022 if PLZ4 == 5021 & gdename == "flombach"
replace PLZ4 = 6055 if PLZ4 == 5055 & gdename == "alpnach dorf"
replace PLZ4 = 8193 if PLZ4 == 5193 & gdename == "eglã¬sau"
replace PLZ4 = 5078 if PLZ4 == 5253 & gdename == "effinaen"
replace PLZ4 = 9423 if PLZ4 == 5499 & gdename == "altenrhein"
replace PLZ4 = 6710 if PLZ4 == 5710 & gdename == "blasco"
replace PLZ4 = 5605 if PLZ4 == 5805 & gdename == "dottlkon"
replace PLZ4 = 5616 if PLZ4 == 5818 & gdename == "melsterschwanden"
replace PLZ4 = 8807 if PLZ4 == 5832 & gdename == "wilen freienbach sz"
replace PLZ4 = 9556 if PLZ4 == 5901 & gdename == "buch bel maerwil"
replace PLZ4 = 8905 if PLZ4 == 5905 & gdename == "arni islisberg"
replace PLZ4 = 6925 if PLZ4 == 5925 & gdename == "gentiilno"
replace PLZ4 = 5430 if PLZ4 == 5930 & gdename == "wettlngen"
replace PLZ4 = 6020 if PLZ4 == 6029 & gdename == "emmenbruecke"
replace PLZ4 = 6078 if PLZ4 == 6087 & gdename == "lungem"
replace PLZ4 = 8127 if PLZ4 == 6127 & gdename == "forch"
replace PLZ4 = 6146 if PLZ4 == 6148 & gdename == "grossdletwil"
replace PLZ4 = 5210 if PLZ4 == 6200 & gdename == "windlsch"
replace PLZ4 = 8308 if PLZ4 == 6308 & gdename == "illnau"
replace PLZ4 = 6410 if PLZ4 == 6358 & gdename == "rigl"
replace PLZ4 = 8903 if PLZ4 == 6413 & gdename == "birmenedorf"
replace PLZ4 = 6410 if PLZ4 == 6420 & gdename == "goldau"
replace PLZ4 = 6438 if PLZ4 == 6435 & gdename == "ibach"
replace PLZ4 = 6410 if PLZ4 == 6470 & gdename == "goldau"
replace PLZ4 = 6598 if PLZ4 == 6588 & gdename == "tenero"
replace PLZ4 = 8617 if PLZ4 == 6617 & gdename == "moenchaltdorf"
replace PLZ4 = 8632 if PLZ4 == 6630 & gdename == "tann rueti"
replace PLZ4 = 8704 if PLZ4 == 6704 & gdename == "herrllberg"
replace PLZ4 = 8706 if PLZ4 == 6706 & gdename == "mellen"
replace PLZ4 = 6760 if PLZ4 == 6782 & gdename == "faldo"
replace PLZ4 = 6644 if PLZ4 == 6844 & gdename == "orseline"
replace PLZ4 = 6976 if PLZ4 == 6909 & gdename == "lugano cassarate"
replace PLZ4 = 6986 if PLZ4 == 6985 & gdename == "cuvio"
replace PLZ4 = 7442 if PLZ4 == 7042 & gdename == "ciugin"
replace PLZ4 = 7260 if PLZ4 == 7280 & gdename == "devos dorf"
replace PLZ4 = 7186 if PLZ4 == 7486 & gdename == "segnas"
replace PLZ4 = 7505 if PLZ4 == 7507 & gdename == "celerina"
replace PLZ4 = 7505 if PLZ4 == 7585 & gdename == "celerina/ schlarigna"
replace PLZ4 = 7536 if PLZ4 == 7631 & gdename == "santa marla im muenstertal"
replace PLZ4 = 7606 if PLZ4 == 7849 & gdename == "promontogno"
replace PLZ4 = 8000 if PLZ4 == 8007 & gdename == "zurich"
replace PLZ4 = 3073 if PLZ4 == 8073 & gdename == "guemligen"
replace PLZ4 = 6132 if PLZ4 == 8131 & gdename == "rohrmatt wllllsau land"
replace PLZ4 = 8162 if PLZ4 == 8161 & gdename == "niedersteinmaur"
replace PLZ4 = 6260 if PLZ4 == 8281 & gdename == "reidermoos"
replace PLZ4 = 8800 if PLZ4 == 8300 & gdename == "thatwit"
replace PLZ4 = 8309 if PLZ4 == 8369 & gdename == "birchwll/ nuerensdorf"
replace PLZ4 = 8758 if PLZ4 == 8375 & gdename == "obstaiden"
replace PLZ4 = 8362 if PLZ4 == 8382 & gdename == "balterswll"
replace PLZ4 = 8586 if PLZ4 == 8386 & gdename == "riedt bei erlen"
replace PLZ4 = 6390 if PLZ4 == 8390 & gdename == "engeiberg"
replace PLZ4 = 5332 if PLZ4 == 8436 & gdename == "reklngen"
replace PLZ4 = 6440 if PLZ4 == 8440 & gdename == "brunnen"
replace PLZ4 = 8624 if PLZ4 == 8642 & gdename == "gruet"
replace PLZ4 = 6644 if PLZ4 == 8644 & gdename == "oreellna"
replace PLZ4 = 6648 if PLZ4 == 8648 & gdename == "minuslo"
replace PLZ4 = 8952 if PLZ4 == 8652 & gdename == "schlieren zh"
replace PLZ4 = 6653 if PLZ4 == 8653 & gdename == "versclo"
replace PLZ4 = 8618 if PLZ4 == 8661 & gdename == "oetwil a/see"
replace PLZ4 = 6760 if PLZ4 == 8780 & gdename == "faldo"
replace PLZ4 = 8617 if PLZ4 == 8817 & gdename == "moenchaitorf"
replace PLZ4 = 8620 if PLZ4 == 8823 & gdename == "wetzlkon"
replace PLZ4 = 6828 if PLZ4 == 8828 & gdename == "belerna"
replace PLZ4 = 8757 if PLZ4 == 8876 & gdename == "fllzbach"
replace PLZ4 = 8135 if PLZ4 == 8935 & gdename == "langnau"
replace PLZ4 = 8966 if PLZ4 == 8936 & gdename == "lieii ag"
replace PLZ4 = 6977 if PLZ4 == 8977 & gdename == "ruvigliana"
replace PLZ4 = 8965 if PLZ4 == 8988 & gdename == "mutschellen"
replace PLZ4 = 9062 if PLZ4 == 9082 & gdename == "lustmuehle"
replace PLZ4 = 9244 if PLZ4 == 9224 & gdename == "niederuzwil"
replace PLZ4 = 9230 if PLZ4 == 9290 & gdename == "flawll"
replace PLZ4 = 8309 if PLZ4 == 9309 & gdename == "oberwii b. nuerensdorf"
replace PLZ4 = 9313 if PLZ4 == 9316 & gdename == "mouien"
replace PLZ4 = 9314 if PLZ4 == 9324 & gdename == "steinebrunn egnach"
replace PLZ4 = 9032 if PLZ4 == 9332 & gdename == "engel burg"
replace PLZ4 = 9052 if PLZ4 == 9352 & gdename == "nlederteufen"
replace PLZ4 = 9410 if PLZ4 == 9416 & gdename == "helden"
replace PLZ4 = 9621 if PLZ4 == 9421 & gdename == "oberheifenschwil"
replace PLZ4 = 9450 if PLZ4 == 9433 & gdename == "luechingen"
replace PLZ4 = 9477 if PLZ4 == 9447 & gdename == "truebbach"
replace PLZ4 = 9630 if PLZ4 == 9530 & gdename == "wattwit"
replace PLZ4 = 9524 if PLZ4 == 9624 & gdename == "zuzwll sg"
replace PLZ4 = 9545 if PLZ4 == 9645 & gdename == "waengl"
replace PLZ4 = 9562 if PLZ4 == 9682 & gdename == "maerwil"
replace PLZ4 = 9604 if PLZ4 == 9804 & gdename == "luetleburg"
replace PLZ4 = 9606 if PLZ4 == 9806 & gdename == "buetechwll"
replace PLZ4 = 9630 if PLZ4 == 9830 & gdename == "wattwll"
replace PLZ4 = 9055 if PLZ4 == 9955 & gdename == "buehlen"
replace PLZ4 = 1212 if PLZ4 == 11212 & gdename == "grand lancy"
replace PLZ4 = 1273 if PLZ4 == 11261 & gdename == "le muids"
replace PLZ4 = 1226 if PLZ4 == 12261 & gdename == "hã³nex"
replace PLZ4 = 1233 if PLZ4 == 12323 & gdename == "bernez"
replace PLZ4 = 8307 if PLZ4 == 13307 & gdename == "effretikon"
replace PLZ4 = 2300 if PLZ4 == 23004 & gdename == "chaux de fonds"
replace PLZ4 = 3464 if PLZ4 == 33451 & gdename == "schmidigen"
replace PLZ4 = 4000 if PLZ4 == 40000 & gdename == "basei"
replace PLZ4 = 4450 if PLZ4 == 44150 & gdename == "slssach"
replace PLZ4 = 6010 if PLZ4 == 60101 & gdename == "riens"
replace PLZ4 = 6808 if PLZ4 == 61307 & gdename == "taverne"
replace PLZ4 = 8360 if PLZ4 == 63660 & gdename == "eschlikon wallenwil"
replace PLZ4 = 7206 if PLZ4 == 72061 & gdename == "gis"
replace PLZ4 = 8000 if PLZ4 == 80313 & gdename == "zurich"
replace PLZ4 = 8413 if PLZ4 == 84112 & gdename == "aesch/neftenbach"
replace PLZ4 = 5463 if PLZ4 == 84398 & gdename == "mellsdorf"
replace PLZ4 = 8547 if PLZ4 == 85474 & gdename == "oberopferstofen"
replace PLZ4 = 8707 if PLZ4 == 88707 & gdename == "uetikon a. see"
replace PLZ4 = 8127 if PLZ4 == 89127 & gdename == "forch/aesch"
replace PLZ4 = 9500 if PLZ4 == 95010 & gdename == "wii (sg)"
replace PLZ4 = 9621 if PLZ4 == 96210 & gdename == "berhelfenschwi i"
replace PLZ4 = 9306 if PLZ4 == 9307 & gdename == "freihof"
replace PLZ4 = 6078 if PLZ4 == 6099 & gdename == "buergten"
replace PLZ4 = 1231 if PLZ4 == 1 & gdename == "conches"
replace PLZ4 = 5634 if PLZ4 == 5635 & gdename == "hagnau"
replace PLZ4 = 3067 if PLZ4 == 1067 & gdename == "boll"
replace PLZ4 = 1209 if PLZ4 == 19 & gdename == "av. mervelet"
replace PLZ4 = 8492 if PLZ4 == 6492 & gdename == "wile"
replace PLZ4 = 6340 if PLZ4 == 5340 & gdename == "saar"
replace PLZ4 = 8594 if PLZ4 == 8504 & gdename == "goettingen"

replace gdename = "koeniz" if gdename == "bem" & PLZ4 == 3097
replace gdename = "koeniz" if gdename == "liebefeld" & PLZ4 == 3097
replace gdename = "kirchlindach" if gdename == "oberlindach" & PLZ4 == 3038
replace gdename = "koeniz" if gdename == "s" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "s pegel bei bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "solegel" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiege bei bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegei bei bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel ." & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel / be" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel b. be" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel b. bem" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel b. bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel b. koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel b.b" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel b.bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel be" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei bem" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei bem/koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei berg" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei bern /koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei bern/ koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei bern/koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei bern/koenlz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei bern/koniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bei n" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bel bem" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bel bem/koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bel bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bel bern /koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bel bern /koenlz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bel bern/koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bel bern/koenlz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bet bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel bet bern/koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel koenitz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel koenlz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel/be" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel/bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel/koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegel/koenlz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spiegelb. bern" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "spirgrl bei koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "splegel koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "splegel koenlz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "splegel/koeniz" & PLZ4 == 3095
replace gdename = "koeniz" if gdename == "splegelb. bern" & PLZ4 == 3095
replace gdename = "berikon" if gdename == ".mutschellen" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "berikon mutschellen" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "berlkon mutschellen" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "muetschellen" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mulschellen" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutscheelen" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutscheilen" & PLZ4 == 8965
replace gdename = "zufikon" if gdename == "mutscheilen/zufikon" & PLZ4 == 5621
replace gdename = "berikon" if gdename == "mutscheiten" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschelen" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschelfen" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschelien" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschelien/berikon" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschellem" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschellen" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschellen (ag)" & PLZ4 == 8965
replace gdename = "zufikon" if gdename == "mutschellen (zufikon)" & PLZ4 == 5621
replace gdename = "zufikon" if gdename == "mutschellen (zutikon)" & PLZ4 == 5621
replace gdename = "berikon" if gdename == "mutschellen ag" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschellen berikon" & PLZ4 == 8965
replace gdename = "rudolfstetten" if gdename == "mutschellen friedlisberg" & PLZ4 == 8964
replace gdename = "rudolfstetten" if gdename == "mutschellen rudolfstetten" & PLZ4 == 8964
replace gdename = "zufikon" if gdename == "mutschellen zufikon" & PLZ4 == 5621
replace gdename = "zufikon" if gdename == "mutschellen zuflkon" & PLZ4 == 5621
replace gdename = "berikon" if gdename == "mutschellen/" & PLZ4 == 8965
replace gdename = "rudolfstetten" if gdename == "mutschellen/ rudolfstetten" & PLZ4 == 8964
replace gdename = "rudolfstetten" if gdename == "mutschellen/ rudolfstetten friedlisberg" & PLZ4 == 8964
replace gdename = "berikon" if gdename == "mutschellen/berikon" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschellen/berlkon" & PLZ4 == 8965
replace gdename = "rudolfstetten" if gdename == "mutschellen/rudolfstetten friedlisberg" & PLZ4 == 8964
replace gdename = "widen" if gdename == "mutschellen/widen" & PLZ4 == 8967
replace gdename = "zufikon" if gdename == "mutschellen/zufikon" & PLZ4 == 5621
replace gdename = "zufikon" if gdename == "mutschellen/zuflkon" & PLZ4 == 5621
replace gdename = "berikon" if gdename == "mutscheller" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschellien" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mutschelten" & PLZ4 == 8965
replace gdename = "zufikon" if gdename == "mutschelten/zufikon" & PLZ4 == 5621
replace gdename = "berikon" if gdename == "mutscheuen" & PLZ4 == 8965
replace gdename = "berikon" if gdename == "mãºtschellen" & PLZ4 == 8965
replace gdename = "oberwil lieli" if gdename == "oberwll lleli" & PLZ4 == 8966
replace gdename = "rudolfstetten" if gdename == "rudotfstetten friedlisberg" & PLZ4 == 8964
replace gdename = "lugano" if gdename == "villa luganoâ»" & PLZ4 == 6966
replace gdename = "berikon" if gdename == "widen/mutschellen" & PLZ4 == 8965
replace gdename = "zufikon" if gdename == "zufikon mutschellen" & PLZ4 == 5621
replace gdename = "zufikon" if gdename == "zufikon/ mutschellen" & PLZ4 == 5621
replace gdename = "zufikon" if gdename == "zuflkon" & PLZ4 == 5621
replace gdename = "moeriken wildegg" if gdename == "marsken wildegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "merã¬ken wildegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "milken (ag)" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moe riken wildegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriiken wildegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken (ag)" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken ag" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken wadegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken waldegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken widegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken wiidegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken wldegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken wlidegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken wlldegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken wã¬ldegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moeriken wïldegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moerikon" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moerken (ag)" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moerlken" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moerlken (ag)" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moerlken ag" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moerlken wadegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moerlken waldegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moerlken wildegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moerlken wlldegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moriken (ag)" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moriken ag" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "moriken wildegg" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "morlken (ag)" & PLZ4 == 5103
replace gdename = "moeriken wildegg" if gdename == "mã³riken (ag)" & PLZ4 == 5103
replace gdename = "obersiggenthal" if gdename == "obersiggental" & PLZ4 == 5415
replace gdename = "moeriken wildegg" if gdename == "wildegg" & PLZ4 == 5103
replace gdename = "lutry" if gdename == "bossiere lutry" & PLZ4 == 1095
replace gdename = "lutry" if gdename == "croix" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "croix lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "croix sur lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "escherin" & PLZ4 == 1095
replace gdename = "lutry" if gdename == "la" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la chroix" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la chrolx" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la cr" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croise s lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix (luini)" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix (luny)" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix (lutry)" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix / lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix s lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix s/ lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix s/lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix sl lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix sur lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix(/lurty)" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix/lurty" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la croix/lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la crolx lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la crolx sur lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la crolx/lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "la crolxllutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "lacroix s/lutry" & PLZ4 == 1090
replace gdename = "lutry" if gdename == "le croix sur lutry" & PLZ4 == 1090
replace gdename = "baar" if gdename == "aalenwinden" & PLZ4 == 6319
replace gdename = "baar" if gdename == "ailenwinden" & PLZ4 == 6319
replace gdename = "baar" if gdename == "allenwinden" & PLZ4 == 6319
replace gdename = "baar" if gdename == "allenwinden baar" & PLZ4 == 6319
replace gdename = "baar" if gdename == "allenwinden/baar" & PLZ4 == 6319
replace gdename = "baar" if gdename == "allenwinden/baer" & PLZ4 == 6319
replace gdename = "baar" if gdename == "allenwlnden" & PLZ4 == 6319
replace gdename = "baar" if gdename == "allenwlnden baar" & PLZ4 == 6319
replace gdename = "baar" if gdename == "allenwlnden saar" & PLZ4 == 6319
replace gdename = "baar" if gdename == "altenwinden" & PLZ4 == 6319
replace gdename = "baar" if gdename == "altenwinden baar" & PLZ4 == 6319
replace gdename = "baar" if gdename == "altenwinden saar" & PLZ4 == 6319
replace gdename = "baar" if gdename == "altenwinden/baar" & PLZ4 == 6319
replace gdename = "baar" if gdename == "altenwinden/baer" & PLZ4 == 6319
replace gdename = "baar" if gdename == "asienwinden" & PLZ4 == 6319
replace gdename = "menzingen" if gdename == "ediibach" & PLZ4 == 6313
replace gdename = "menzingen" if gdename == "edlibach" & PLZ4 == 6313
replace gdename = "menzingen" if gdename == "edlibach (menzingen)" & PLZ4 == 6313
replace gdename = "menzingen" if gdename == "edlibach menzingen" & PLZ4 == 6313
replace gdename = "menzingen" if gdename == "edllbach" & PLZ4 == 6313
replace gdename = "menzingen" if gdename == "finstersee" & PLZ4 == 6313
replace gdename = "menzingen" if gdename == "finstersee /menzingen" & PLZ4 == 6313
replace gdename = "menzingen" if gdename == "finstersee menzingen" & PLZ4 == 6313
replace gdename = "menzingen" if gdename == "finstersee/menzingen" & PLZ4 == 6313
replace gdename = "menzingen" if gdename == "flnstersee" & PLZ4 == 6313
replace gdename = "menzingen" if gdename == "flnstersee/menzingen" & PLZ4 == 6313
replace gdename = "oberaegeri" if gdename == "morgarten" & PLZ4 == 6315
replace gdename = "oberaegeri" if gdename == "morgarten oberaegeri" & PLZ4 == 6315
replace gdename = "oberaegeri" if gdename == "morgarten oberaegerl" & PLZ4 == 6315
replace gdename = "unteraegeri" if gdename == "neuaegeri" & PLZ4 == 6314
replace gdename = "unteraegeri" if gdename == "pffistersee" & PLZ4 == 6314
replace gdename = "grone" if gdename == "(cogne" & PLZ4 == 3979
replace gdename = "grone" if gdename == ".groene" & PLZ4 == 3979
replace gdename = "agarn" if gdename == "agam" & PLZ4 == 3951
replace gdename = "agarn" if gdename == "agern" & PLZ4 == 3951
replace gdename = "albinen" if gdename == "aibinen" & PLZ4 == 3955
replace gdename = "albinen" if gdename == "alblnen" & PLZ4 == 3955
replace gdename = "crans montana" if gdename == "champzabe" & PLZ4 == 3963
replace gdename = "grone" if gdename == "cogne" & PLZ4 == 3979
replace gdename = "crans montana" if gdename == "crans" & PLZ4 == 3963
replace gdename = "crans montana" if gdename == "crans lens" & PLZ4 == 3963
replace gdename = "crans montana" if gdename == "crans s" & PLZ4 == 3963
replace gdename = "crans montana" if gdename == "crans s/lens" & PLZ4 == 3963
replace gdename = "crans montana" if gdename == "crans sur lens" & PLZ4 == 3963
replace gdename = "crans montana" if gdename == "crans/lens" & PLZ4 == 3963
replace gdename = "eischoll" if gdename == "eilscholl" & PLZ4 == 3943
replace gdename = "eischoll" if gdename == "eischola" & PLZ4 == 3943
replace gdename = "eischoll" if gdename == "elscholl" & PLZ4 == 3943
replace gdename = "eischoll" if gdename == "etscholl" & PLZ4 == 3943
replace gdename = "lens" if gdename == "fianthey" & PLZ4 == 1978
replace gdename = "lens" if gdename == "flantey" & PLZ4 == 1978
replace gdename = "lens" if gdename == "flanthey" & PLZ4 == 1978
replace gdename = "lens" if gdename == "flanthey lens" & PLZ4 == 1978
replace gdename = "grone" if gdename == "grane" & PLZ4 == 3979
replace gdename = "grone" if gdename == "grbne" & PLZ4 == 3979
replace gdename = "grone" if gdename == "grbno" & PLZ4 == 3979
replace gdename = "grone" if gdename == "greene" & PLZ4 == 3979
replace gdename = "grone" if gdename == "grene" & PLZ4 == 3979
replace gdename = "grone" if gdename == "grime" & PLZ4 == 3979
replace gdename = "grone" if gdename == "grmne" & PLZ4 == 3979
replace gdename = "grone" if gdename == "groene" & PLZ4 == 3979
replace gdename = "grone" if gdename == "grã³ne" & PLZ4 == 3979
replace gdename = "guttet feschel" if gdename == "gullet" & PLZ4 == 3956
replace gdename = "guttet feschel" if gdename == "gutten" & PLZ4 == 3956
replace gdename = "guttet feschel" if gdename == "gutten vs" & PLZ4 == 3956
replace gdename = "grone" if gdename == "lcogne" & PLZ4 == 3979
replace gdename = "grone" if gdename == "lcogre" & PLZ4 == 3979
replace gdename = "sierre" if gdename == "no" & PLZ4 == 3960
replace gdename = "sierre" if gdename == "noaes/sierre" & PLZ4 == 3960
replace gdename = "sierre" if gdename == "nobs" & PLZ4 == 3960
replace gdename = "sierre" if gdename == "noees" & PLZ4 == 3960
replace gdename = "sierre" if gdename == "noes" & PLZ4 == 3960
replace gdename = "sierre" if gdename == "noes sierre" & PLZ4 == 3960
replace gdename = "sierre" if gdename == "noes/sierre" & PLZ4 == 3960
replace gdename = "sierre" if gdename == "noes/slerre" & PLZ4 == 3960
replace gdename = "sierre" if gdename == "nooes/sierre" & PLZ4 == 3960
replace gdename = "rarogne" if gdename == "rarogne" & PLZ4 == 3942
replace gdename = "unterbaech" if gdename == "unterbach" & PLZ4 == 3944
replace gdename = "unterbaech" if gdename == "unterbaech vs" & PLZ4 == 3944
replace gdename = "unterbaech" if gdename == "unterbaechn" & PLZ4 == 3944
replace gdename = "egg" if gdename == "hinteregg" & PLZ4 == 8132
replace gdename = "egg" if gdename == "hinteregg bei zuerich" & PLZ4 == 8132
replace gdename = "egg" if gdename == "hinteregg egg" & PLZ4 == 8132
replace gdename = "egg" if gdename == "hinteregg zh" & PLZ4 == 8132
replace gdename = "egg" if gdename == "hinteregg/egg bei zuerich" & PLZ4 == 8132
replace gdename = "egg" if gdename == "hinteregg/egg bel zuerich" & PLZ4 == 8132
replace gdename = "egg" if gdename == "hinteregg/egg zh" & PLZ4 == 8132
replace gdename = "egg" if gdename == "hinterwegg" & PLZ4 == 8132
replace gdename = "egg" if gdename == "hlnteregg" & PLZ4 == 8132
replace gdename = "egg" if gdename == "hlnteregg zh" & PLZ4 == 8132
replace gdename = "zumikon" if gdename == "zumlkon" & PLZ4 == 8126
replace gdename = "sainte croix" if gdename == "auberson" & PLZ4 == 1454
replace gdename = "sainte croix" if gdename == "l auberson" & PLZ4 == 1454
replace gdename = "sainte croix" if gdename == "l auberson /ste croix" & PLZ4 == 1454
replace gdename = "sainte croix" if gdename == "la bagne vd" & PLZ4 == 1450
replace gdename = "sainte croix" if gdename == "la sagne pres ste croix" & PLZ4 == 1450
replace gdename = "sainte croix" if gdename == "la sagne vd" & PLZ4 == 1450
replace gdename = "sainte croix" if gdename == "la segne vd" & PLZ4 == 1450
replace gdename = "bullet" if gdename == "lea rassea" & PLZ4 == 1452
replace gdename = "bullet" if gdename == "les basses" & PLZ4 == 1452
replace gdename = "bullet" if gdename == "les cluds" & PLZ4 == 1453
replace gdename = "bullet" if gdename == "les cluds riere bullet" & PLZ4 == 1453
replace gdename = "bullet" if gdename == "les r a sse s" & PLZ4 == 1452
replace gdename = "bullet" if gdename == "les rasses" & PLZ4 == 1452
replace gdename = "bullet" if gdename == "les rassies" & PLZ4 == 1452
replace gdename = "alto malcantone" if gdename == "aroslo" & PLZ4 == 6939
replace gdename = "bedano" if gdename == "badano" & PLZ4 == 6930
replace gdename = "lugano" if gdename == "barbengo/lugano" & PLZ4 == 6917
replace gdename = "alto malcantone" if gdename == "brano" & PLZ4 == 6937
replace gdename = "brusino arsizio" if gdename == "brasino arsizio" & PLZ4 == 6827
replace gdename = "lugano" if gdename == "bre" & PLZ4 == 6979
replace gdename = "lugano" if gdename == "bre di lugano" & PLZ4 == 6979
replace gdename = "lugano" if gdename == "bre s/lugano" & PLZ4 == 6979
replace gdename = "lugano" if gdename == "bre sopra lugano" & PLZ4 == 6979
replace gdename = "lugano" if gdename == "brit dl lugano" & PLZ4 == 6979
replace gdename = "brusino arsizio" if gdename == "brusino" & PLZ4 == 6827
replace gdename = "brusino arsizio" if gdename == "brusino arsizlo" & PLZ4 == 6827
replace gdename = "brusino arsizio" if gdename == "brusino arslzlo" & PLZ4 == 6827
replace gdename = "brusino arsizio" if gdename == "brusio arizo" & PLZ4 == 6827
replace gdename = "brusino arsizio" if gdename == "bruslno arsizio" & PLZ4 == 6827
replace gdename = "brusino arsizio" if gdename == "brusã¬no arsizio" & PLZ4 == 6827
replace gdename = "brusino arsizio" if gdename == "busino arsizio" & PLZ4 == 6827
replace gdename = "lugano" if gdename == "carabbio" & PLZ4 == 6913
replace gdename = "lugano" if gdename == "carabbla" & PLZ4 == 6913
replace gdename = "lugano" if gdename == "carabblia" & PLZ4 == 6913
replace gdename = "collina d'oro" if gdename == "carabietta montagnoia" & PLZ4 == 6926
replace gdename = "collina d'oro" if gdename == "carabietta montagnola" & PLZ4 == 6926
replace gdename = "collina d'oro" if gdename == "carabletta" & PLZ4 == 6926
replace gdename = "collina d'oro" if gdename == "carabletta montagnola" & PLZ4 == 6926
replace gdename = "lugano" if gdename == "carrabia" & PLZ4 == 6913
replace gdename = "comano" if gdename == "ccmano" & PLZ4 == 6949
replace gdename = "comano" if gdename == "comano /ti" & PLZ4 == 6949
replace gdename = "comano" if gdename == "coniano" & PLZ4 == 6949
replace gdename = "comano" if gdename == "cumano" & PLZ4 == 6949
replace gdename = "alto malcantone" if gdename == "fescoggla" & PLZ4 == 6938
replace gdename = "collina d'oro" if gdename == "gentillno" & PLZ4 == 6926
replace gdename = "grancia" if gdename == "grancla" & PLZ4 == 6916
replace gdename = "grancia" if gdename == "grand e" & PLZ4 == 6916
replace gdename = "grancia" if gdename == "grande" & PLZ4 == 6916
replace gdename = "grancia" if gdename == "granola" & PLZ4 == 6916
replace gdename = "grancia" if gdename == "grenela" & PLZ4 == 6916
replace gdename = "gravesano" if gdename == "grumo di gravesano" & PLZ4 == 6929
replace gdename = "gravesano" if gdename == "grumo dl gravesano" & PLZ4 == 6929
replace gdename = "gravesano" if gdename == "lugano gravesano" & PLZ4 == 6929
replace gdename = "manno" if gdename == "manna" & PLZ4 == 6928
replace gdename = "collina d'oro" if gdename == "montagnola carabietta" & PLZ4 == 6926
replace gdename = "collina d'oro" if gdename == "montagnola carabletta" & PLZ4 == 6926
replace gdename = "lugano" if gdename == "monte bre lugano" & PLZ4 == 6900
replace gdename = "lugano" if gdename == "monte bre/lugano" & PLZ4 == 6900
replace gdename = "lugano" if gdename == "noracnco" & PLZ4 == 6915
replace gdename = "lugano" if gdename == "noranco" & PLZ4 == 6915
replace gdename = "lugano" if gdename == "pambio" & PLZ4 == 6915
replace gdename = "lugano" if gdename == "pambio norano" & PLZ4 == 6915
replace gdename = "lugano" if gdename == "pamblo" & PLZ4 == 6915
replace gdename = "lugano" if gdename == "pamblo noranco" & PLZ4 == 6915
replace gdename = "lugano" if gdename == "pampio noranco" & PLZ4 == 6915
replace gdename = "bedano" if gdename == "sedano" & PLZ4 == 6930
replace gdename = "alto malcantone" if gdename == "vezlo" & PLZ4 == 6938
replace gdename = "vico morcote" if gdename == "vice morcote" & PLZ4 == 6921
replace gdename = "vico morcote" if gdename == "vico miorcote" & PLZ4 == 6921
replace gdename = "vico morcote" if gdename == "vida morcote" & PLZ4 == 6921
replace gdename = "vico morcote" if gdename == "vlco morcote" & PLZ4 == 6921
replace gdename = "saint leonard" if gdename == "8t leonard" & PLZ4 == 1958
replace gdename = "ayent" if gdename == "signese st leonhard" & PLZ4 == 1966
replace gdename = "ayent" if gdename == "slignese st leonhard" & PLZ4 == 1966
replace gdename = "saint leonard" if gdename == "st leonard" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "st leonard uvrier" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "st leonard/sion" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "st leonard/slon" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "st leonhard" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "st. leonard" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "st. leonhard" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvrier" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvrier gjslon" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvrier s/sion" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvrier sion" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvrier st leonard" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvrier st leonhard" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvrier/sion" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvrier/st leonard" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvrier/st. leonard" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvriersion" & PLZ4 == 1958
replace gdename = "saint leonard" if gdename == "uvrïer/st leonard" & PLZ4 == 1958
replace gdename = "gipf oberfrick" if gdename == "gipf 0berfrick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gipf oberf rick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gipf oberflick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gipf oberfrack" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gipf oberfriek" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gipf oberfrlck" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gipf oberirlck" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gipf obertrick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gipf.obertrick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gipt oberfrack" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gipt oberfrick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "glpf bberfrick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "glpf oberfrack" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "glpf oberfrick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "glpf oberfrlck" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "glpf oberfrlek" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "glpf oberirick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "glpf obertrick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "gã¬pf oberfrick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "oberf rick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "oberfrick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "oberfrlck" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "obertrick" & PLZ4 == 5073
replace gdename = "gipf oberfrick" if gdename == "obertrickâ" & PLZ4 == 5073
replace gdename = "davos" if gdename == "davos laret" & PLZ4 == 7265
replace gdename = "fideris" if gdename == "fideris dorf" & PLZ4 == 7235
replace gdename = "fideris" if gdename == "fideris station" & PLZ4 == 7235
replace gdename = "fideris" if gdename == "fiderls dorf" & PLZ4 == 7235
replace gdename = "fideris" if gdename == "fiderls station" & PLZ4 == 7235
replace gdename = "fideris" if gdename == "flderis dorf" & PLZ4 == 7235
replace gdename = "fideris" if gdename == "flderla" & PLZ4 == 7235
replace gdename = "fideris" if gdename == "flderls station" & PLZ4 == 7235
replace gdename = "furna" if gdename == "fuma" & PLZ4 == 7232
replace gdename = "furna" if gdename == "fuma station" & PLZ4 == 7232
replace gdename = "furna" if gdename == "fume" & PLZ4 == 7232
replace gdename = "furna" if gdename == "furna dorf" & PLZ4 == 7232
replace gdename = "furna" if gdename == "furna dort" & PLZ4 == 7232
replace gdename = "furna" if gdename == "furna station" & PLZ4 == 7232
replace gdename = "klosters serneus" if gdename == "saas 1. praettigau" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saas i. p" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saas i. praettigau" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saas i. prattigau" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saas i.p" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saas im praettigau" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saas praedigau" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saas praedlgau" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saas praettigau" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saes" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saes i. praettigau" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "saes im praettigau" & PLZ4 == 7250
replace gdename = "schmitten (gr)" if gdename == "schmitten seewis" & PLZ4 == 7493
replace gdename = "klosters serneus" if gdename == "seas" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "seas i. praettigau" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "seas im praettigau" & PLZ4 == 7250
replace gdename = "seewis im praettigau" if gdename == "seewcs" & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewcs dorf" & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewcs dort" & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewcs station" & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewis" & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewis dorf" & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewis i. pr." & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewis i. praettigaue" & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewis im praettigeu" & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewis station" & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewls" & PLZ4 == 7212
replace gdename = "seewis im praettigau" if gdename == "seewls dorf" & PLZ4 == 7212
replace gdename = "seltisberg" if gdename == "seitisberg" & PLZ4 == 4411
replace gdename = "klosters serneus" if gdename == "serneus" & PLZ4 == 7250
replace gdename = "davos" if gdename == "sertig" & PLZ4 == 7260
replace gdename = "davos" if gdename == "wolfgang" & PLZ4 == 7265
replace gdename = "davos" if gdename == "wolfgang bei davos" & PLZ4 == 7265
replace gdename = "davos" if gdename == "wolfgang bel davos" & PLZ4 == 7265
replace gdename = "altstaetten" if gdename == "altstaetten/luechingen" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luchingen" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luchingen (sg)" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luchingen sg" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luechangen" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luechingen" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luechingen (sg)" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luechingen/ altstaetten" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luechingen/ altstaetten sg" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luechingen/altstaetten" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luechingen/altstaetten sg" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luechlngen" & PLZ4 == 9450
replace gdename = "altstaetten" if gdename == "luechlngen (sg)" & PLZ4 == 9450
replace gdename = "isone" if gdename == "leone" & PLZ4 == 6810
replace gdename = "isone" if gdename == "lsone" & PLZ4 == 6810
replace gdename = "monteceneri" if gdename == "medegiia" & PLZ4 == 6809
replace gdename = "monteceneri" if gdename == "medeglia manz isolazioni sa manne" & PLZ4 == 6809
replace gdename = "monteceneri" if gdename == "medeglla" & PLZ4 == 6809
replace gdename = "mezzovico vira" if gdename == "meuovico" & PLZ4 == 6805
replace gdename = "mezzovico vira" if gdename == "mezzovicco" & PLZ4 == 6805
replace gdename = "mezzovico vira" if gdename == "mezzovico" & PLZ4 == 6805
replace gdename = "mezzovico vira" if gdename == "mezzovlco" & PLZ4 == 6805
replace gdename = "arogno" if gdename == "pugerna" & PLZ4 == 6822
replace gdename = "rovio" if gdename == "rovia" & PLZ4 == 6821
replace gdename = "rovio" if gdename == "roviio" & PLZ4 == 6821
replace gdename = "rovio" if gdename == "rovio ti" & PLZ4 == 6821
replace gdename = "rovio" if gdename == "rovioi" & PLZ4 == 6821
replace gdename = "rovio" if gdename == "rovlo" & PLZ4 == 6821
replace gdename = "monteceneri" if gdename == "sigirino ti" & PLZ4 == 6806
replace gdename = "monteceneri" if gdename == "siglrino" & PLZ4 == 6806
replace gdename = "monteceneri" if gdename == "slgirino" & PLZ4 == 6806
replace gdename = "monteceneri" if gdename == "slglrino" & PLZ4 == 6806
replace gdename = "aclens" if gdename == "aciens" & PLZ4 == 1123
replace gdename = "bremblens" if gdename == "brembiens" & PLZ4 == 1121
replace gdename = "bremblens" if gdename == "bremlens" & PLZ4 == 1121
replace gdename = "bremblens" if gdename == "brernblens" & PLZ4 == 1121
replace gdename = "clarmont" if gdename == "ciarmont" & PLZ4 == 1127
replace gdename = "clarmont" if gdename == "clermont" & PLZ4 == 1127
replace gdename = "echichens" if gdename == "colombier s/morges" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "colombier sur marges" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "colombier sur morges" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "colombler sur morges" & PLZ4 == 1112
replace gdename = "cottens (vd)" if gdename == "cottens vaud" & PLZ4 == 1116
replace gdename = "echichens" if gdename == "echiens" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "echlchens" & PLZ4 == 1112
replace gdename = "chigny" if gdename == "en chantemerte" & PLZ4 == 1134
replace gdename = "gollion" if gdename == "goiiion" & PLZ4 == 1124
replace gdename = "gollion" if gdename == "goilion" & PLZ4 == 1124
replace gdename = "gollion" if gdename == "golfion" & PLZ4 == 1124
replace gdename = "gollion" if gdename == "goliion" & PLZ4 == 1124
replace gdename = "gollion" if gdename == "golllon" & PLZ4 == 1124
replace gdename = "grancy" if gdename == "grancy vd" & PLZ4 == 1117
replace gdename = "lully (vd)" if gdename == "luffy (vd)" & PLZ4 == 1132
replace gdename = "lully (vd)" if gdename == "luily" & PLZ4 == 1132
replace gdename = "lully (vd)" if gdename == "luliy" & PLZ4 == 1132
replace gdename = "lully (vd)" if gdename == "luliy sur morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lullt vd" & PLZ4 == 1132
replace gdename = "lully vd" if gdename == "lully morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lully s marges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lully s morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lully s. morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lully s/morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lully sur mordes" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lully sur morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lully/morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lusey sur morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lussy s morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lussy s. morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lussy s/morges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lussy sur moorges" & PLZ4 == 1167
replace gdename = "lully (vd)" if gdename == "lussy sur morgen" & PLZ4 == 1167
replace gdename = "echichens" if gdename == "monnaz sur morges" & PLZ4 == 1112
replace gdename = "reverolle" if gdename == "reveroile" & PLZ4 == 1128
replace gdename = "reverolle" if gdename == "reverolie" & PLZ4 == 1128
replace gdename = "romanel sur morges" if gdename == "romane) s morges" & PLZ4 == 1122
replace gdename = "romanel sur lausanne" if gdename == "romanei sur lausanne" & PLZ4 == 1032
replace gdename = "romanel sur morges" if gdename == "romanei sur morges" & PLZ4 == 1122
replace gdename = "romanel sur morges" if gdename == "romanel s morges" & PLZ4 == 1122
replace gdename = "senarclens" if gdename == "senarciens" & PLZ4 == 1304
replace gdename = "echichens" if gdename == "st saphorin" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "st saphorin (merges)" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "st saphorin (morges)" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "st saphorin morges" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "st saphorin s morges" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "st saphorin s/morges" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "st saphorin sur morges" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "st saphorin sur morses" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "st saphorã­n sur morges" & PLZ4 == 1112
replace gdename = "echichens" if gdename == "st. saphorin sur morges" & PLZ4 == 1112
replace gdename = "tolochenaz" if gdename == "toiochena" & PLZ4 == 1131
replace gdename = "tolochenaz" if gdename == "tolochena" & PLZ4 == 1131
replace gdename = "vullierens" if gdename == "vuillerens" & PLZ4 == 1115
replace gdename = "vullierens" if gdename == "vuliierens" & PLZ4 == 1115
replace gdename = "gossau" if gdename == "bertschikon (gossau zh)" & PLZ4 == 8625
replace gdename = "gossau" if gdename == "bertschikon b. gossau zh" & PLZ4 == 8625
replace gdename = "gossau" if gdename == "bertschikon b. uster" & PLZ4 == 8625
replace gdename = "gossau" if gdename == "bertschikon gossau" & PLZ4 == 8625
replace gdename = "gossau" if gdename == "bertschikon zh" & PLZ4 == 8625
replace gdename = "gossau" if gdename == "bertschlkon" & PLZ4 == 8625
replace gdename = "gossau" if gdename == "bertschlkon (gossau zh)" & PLZ4 == 8625
replace gdename = "gossau" if gdename == "bertschlkon gossau" & PLZ4 == 8625
replace gdename = "gossau" if gdename == "bertschlkon zh" & PLZ4 == 8625
replace gdename = "uster" if gdename == "freudwil" & PLZ4 == 8610
replace gdename = "uster" if gdename == "freudwll" & PLZ4 == 8610
replace gdename = "gossau" if gdename == "gossau bertschikon" & PLZ4 == 8625
replace gdename = "uster" if gdename == "riedikon" & PLZ4 == 8610
replace gdename = "uster" if gdename == "sulzbach" & PLZ4 == 8614
replace gdename = "uster" if gdename == "sulzbach b. uster" & PLZ4 == 8614
replace gdename = "uster" if gdename == "sulzbach bei uster" & PLZ4 == 8614
replace gdename = "uster" if gdename == "wermatsw l" & PLZ4 == 8615
replace gdename = "uster" if gdename == "wermatswii uster" & PLZ4 == 8615
replace gdename = "uster" if gdename == "wermatswil" & PLZ4 == 8615
replace gdename = "uster" if gdename == "wermatswil uster" & PLZ4 == 8615
replace gdename = "uster" if gdename == "wermatswil/uster" & PLZ4 == 8615
replace gdename = "uster" if gdename == "wermatswll" & PLZ4 == 8615
replace gdename = "uster" if gdename == "wermatswll uster" & PLZ4 == 8615
replace gdename = "uster" if gdename == "wermetswil" & PLZ4 == 8615
replace gdename = "domleschg" if gdename == "almens gr" & PLZ4 == 7416
replace gdename = "albula" if gdename == "alvaneu bad" & PLZ4 == 7492
replace gdename = "albula" if gdename == "alvaneu dorf" & PLZ4 == 7492
replace gdename = "albula" if gdename == "alveneu" & PLZ4 == 7492
replace gdename = "albula" if gdename == "alveneu bad" & PLZ4 == 7492
replace gdename = "domleschg" if gdename == "atmens" & PLZ4 == 7416
replace gdename = "cazis" if gdename == "cazls" & PLZ4 == 7408
replace gdename = "domleschg" if gdename == "feidis veulden" & PLZ4 == 7404
replace gdename = "domleschg" if gdename == "feldas veulden" & PLZ4 == 7404
replace gdename = "domleschg" if gdename == "feldas/veulden" & PLZ4 == 7404
replace gdename = "domleschg" if gdename == "feldis" & PLZ4 == 7404
replace gdename = "domleschg" if gdename == "feldis veulden" & PLZ4 == 7404
replace gdename = "domleschg" if gdename == "feldisneulden" & PLZ4 == 7404
replace gdename = "domleschg" if gdename == "feldisrveulden" & PLZ4 == 7404
replace gdename = "fuerstenau" if gdename == "fuerstenaubruck" & PLZ4 == 7414
replace gdename = "fuerstenau" if gdename == "furstenaubruck" & PLZ4 == 7414
replace gdename = "rhaezuens" if gdename == "rhaezuns" & PLZ4 == 7403
replace gdename = "rhaezuens" if gdename == "rhazuens" & PLZ4 == 7403
replace gdename = "rhaezuens" if gdename == "rhazuns" & PLZ4 == 7403
replace gdename = "sils im domleschg" if gdename == "s is i. domleschg" & PLZ4 == 7411
replace gdename = "schmitten (gr)" if gdename == "schmitten/albula" & PLZ4 == 7493
replace gdename = "sils im domleschg" if gdename == "sii i.d" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sil i" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sil i.d" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils 1. domleschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils 1.domieschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils domieschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils domleschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils i. d" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils i. domieschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils i. domleschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils i.d" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils i.domieschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils i.domleschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils im domieschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils im domleechg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils l.domieschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sils ã­. d" & PLZ4 == 7411
replace gdename = "sils im engadin" if gdename == "sils/e" & PLZ4 == 7514
replace gdename = "sils im engadin" if gdename == "sils/sega maria" & PLZ4 == 7514
replace gdename = "sils im engadin" if gdename == "sils/segl maria" & PLZ4 == 7514
replace gdename = "sils im domleschg" if gdename == "slls i. domieschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "slls i.domieschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "slls i.domleschg" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "stuls" & PLZ4 == 7411
replace gdename = "sils im domleschg" if gdename == "sã­ls i.domleschg" & PLZ4 == 7411
replace gdename = "domleschg" if gdename == "tumegi/tomils" & PLZ4 == 7418
replace gdename = "domleschg" if gdename == "tumegl" & PLZ4 == 7418
replace gdename = "domleschg" if gdename == "tumegl.tomils" & PLZ4 == 7418
replace gdename = "domleschg" if gdename == "tumegl/tamils" & PLZ4 == 7418
replace gdename = "domleschg" if gdename == "tumegl/tomiis" & PLZ4 == 7418
replace gdename = "domleschg" if gdename == "tumegl/tomlls" & PLZ4 == 7418
replace gdename = "domleschg" if gdename == "tumeglitamils" & PLZ4 == 7418
replace gdename = "sils im domleschg" if gdename == "us im domlesahg" & PLZ4 == 7411
replace gdename = "bourg saint pierre" if gdename == "bourg st pierre" & PLZ4 == 1946
replace gdename = "bourg saint pierre" if gdename == "bourg st plerre" & PLZ4 == 1946
replace gdename = "bourg saint pierre" if gdename == "bourg st. pierre" & PLZ4 == 1946
replace gdename = "bovernier" if gdename == "bovemier" & PLZ4 == 1932
replace gdename = "bovernier" if gdename == "bovemler" & PLZ4 == 1932
replace gdename = "bovernier" if gdename == "bovernler" & PLZ4 == 1932
replace gdename = "bagnes" if gdename == "champsec" & PLZ4 == 1934
replace gdename = "orsieres" if gdename == "la fouiy" & PLZ4 == 1937
replace gdename = "orsieres" if gdename == "la fouly" & PLZ4 == 1937
replace gdename = "volleges" if gdename == "leveron" & PLZ4 == 1941
replace gdename = "liddes" if gdename == "liddes vs" & PLZ4 == 1945
replace gdename = "orsieres" if gdename == "praz de fort" & PLZ4 == 1937
replace gdename = "bagnes" if gdename == "sarreyer" & PLZ4 == 1934
replace gdename = "bagnes" if gdename == "v rsegeres" & PLZ4 == 1934
replace gdename = "bagnes" if gdename == "vergenes bagnes" & PLZ4 == 1934
replace gdename = "bagnes" if gdename == "versegere" & PLZ4 == 1934
replace gdename = "bagnes" if gdename == "versegeres" & PLZ4 == 1934
replace gdename = "bagnes" if gdename == "versegeres / prarreyer" & PLZ4 == 1934
replace gdename = "bagnes" if gdename == "versegeres i prarreyer" & PLZ4 == 1934
replace gdename = "bagnes" if gdename == "versegeres prarreyer" & PLZ4 == 1934
replace gdename = "volleges" if gdename == "voileges" & PLZ4 == 1941
replace gdename = "volleges" if gdename == "volieges" & PLZ4 == 1941
replace gdename = "villaz saint pierre" if gdename == "mtn st pierre" & PLZ4 == 1690
replace gdename = "villaz saint pierre" if gdename == "viilaz st pierre" & PLZ4 == 1690
replace gdename = "villaz saint pierre" if gdename == "viliaz st pierre" & PLZ4 == 1690
replace gdename = "villaz saint pierre" if gdename == "viliaz st. pierre" & PLZ4 == 1690
replace gdename = "villaz saint pierre" if gdename == "villan st pierre" & PLZ4 == 1690
replace gdename = "villaz saint pierre" if gdename == "villaz lussy" & PLZ4 == 1690
replace gdename = "villaz saint pierre" if gdename == "villaz st pierre" & PLZ4 == 1690
replace gdename = "villorsonnens" if gdename == "villaz st pierre/chavannes sous orsonnens" & PLZ4 == 1694
replace gdename = "villaz saint pierre" if gdename == "villaz st plerre" & PLZ4 == 1690
replace gdename = "villaz saint pierre" if gdename == "villaz st. pierre" & PLZ4 == 1690
replace gdename = "villaz saint pierre" if gdename == "vlllaz st pierre" & PLZ4 == 1690
replace gdename = "villaz saint pierre" if gdename == "vlllaz st plerre" & PLZ4 == 1690
replace gdename = "flims" if gdename == "fidaz" & PLZ4 == 7019
replace gdename = "flims" if gdename == "fidaz films" & PLZ4 == 7019
replace gdename = "flims" if gdename == "fidaz flims" & PLZ4 == 7019
replace gdename = "flims" if gdename == "fidaz/flims" & PLZ4 == 7019
replace gdename = "flims" if gdename == "fldaz films" & PLZ4 == 7019
replace gdename = "flims" if gdename == "fldaz flims" & PLZ4 == 7019
replace gdename = "arosa" if gdename == "langwien" & PLZ4 == 7050
replace gdename = "arosa" if gdename == "langwies gr" & PLZ4 == 7050
replace gdename = "arosa" if gdename == "langwles" & PLZ4 == 7050
replace gdename = "lantsch" if gdename == "lantsch" & PLZ4 == 7083
replace gdename = "lantsch" if gdename == "lantsch lenz" & PLZ4 == 7083
replace gdename = "lantsch" if gdename == "lantsch%lenz" & PLZ4 == 7083
replace gdename = "lantsch" if gdename == "lantsch. lenz" & PLZ4 == 7083
replace gdename = "lantsch" if gdename == "lantsch/l.crnz" & PLZ4 == 7083
replace gdename = "lantsch" if gdename == "lantsch/lerz" & PLZ4 == 7083
replace gdename = "lantsch" if gdename == "lantschilenz" & PLZ4 == 7083
replace gdename = "arosa" if gdename == "litzirueti" & PLZ4 == 7050
replace gdename = "churwalden" if gdename == "malfix" & PLZ4 == 7074
replace gdename = "vaz" if gdename == "muldain vaz" & PLZ4 == 7082
replace gdename = "vaz" if gdename == "obervaz" & PLZ4 == 7082
replace gdename = "arosa" if gdename == "paglg" & PLZ4 == 7050
replace gdename = "arosa" if gdename == "pelst" & PLZ4 == 7050
replace gdename = "tamins" if gdename == "reichenau" & PLZ4 == 7015
replace gdename = "arosa" if gdename == "sapuen" & PLZ4 == 7050
replace gdename = "trin" if gdename == "trier mulin" & PLZ4 == 7014
replace gdename = "trin" if gdename == "triin" & PLZ4 == 7014
replace gdename = "trin" if gdename == "trin digg" & PLZ4 == 7014
replace gdename = "trin" if gdename == "trin mulin" & PLZ4 == 7014
replace gdename = "trin" if gdename == "trin mulln" & PLZ4 == 7014
replace gdename = "trin" if gdename == "trin mutin" & PLZ4 == 7014
replace gdename = "vaz" if gdename == "vaz" & PLZ4 == 7082
replace gdename = "vaz" if gdename == "vaz obervaz" & PLZ4 == 7082
replace gdename = "vaz" if gdename == "vaziobervaz" & PLZ4 == 7082
replace gdename = "vaz" if gdename == "zorten" & PLZ4 == 7082
replace gdename = "vaz" if gdename == "zorten gde. vaz ober" & PLZ4 == 7082
replace gdename = "vaz" if gdename == "zorten gde. vaz obervaz" & PLZ4 == 7082
replace gdename = "laufen uhwiesen" if gdename == "laufen uhwlesen" & PLZ4 == 8248
replace gdename = "laufen uhwiesen" if gdename == "lauten uhwlesen" & PLZ4 == 8248
replace gdename = "laufen uhwiesen" if gdename == "uhwieaen" & PLZ4 == 8248
replace gdename = "laufen uhwiesen" if gdename == "uhwiesen" & PLZ4 == 8248
replace gdename = "laufen uhwiesen" if gdename == "uhwiesen laufen" & PLZ4 == 8248
replace gdename = "laufen uhwiesen" if gdename == "uhwlesen" & PLZ4 == 8248
replace gdename = "beinwil " if gdename == "beinwal (so)" & PLZ4 == 4229
replace gdename = "beinwil " if gdename == "beinwal so" & PLZ4 == 4229
replace gdename = "beinwil " if gdename == "beinwll" & PLZ4 == 4229
replace gdename = "beinwil " if gdename == "belnwil (so)" & PLZ4 == 4229
replace gdename = "blauen" if gdename == "blauen be" & PLZ4 == 4223
replace gdename = "fehren" if gdename == "fahren" & PLZ4 == 4232
replace gdename = "grindel" if gdename == "grinde)" & PLZ4 == 4247
replace gdename = "grindel" if gdename == "grindel so" & PLZ4 == 4247
replace gdename = "himmelried" if gdename == "hammelried" & PLZ4 == 4204
replace gdename = "himmelried" if gdename == "hammelried so" & PLZ4 == 4204
replace gdename = "himmelried" if gdename == "hilmmeiried" & PLZ4 == 4204
replace gdename = "himmelried" if gdename == "himmeiried" & PLZ4 == 4204
replace gdename = "himmelried" if gdename == "himmelrfed" & PLZ4 == 4204
replace gdename = "himmelried" if gdename == "himmelried so" & PLZ4 == 4204
replace gdename = "himmelried" if gdename == "himmelrled" & PLZ4 == 4204
replace gdename = "himmelried" if gdename == "himmelrled so" & PLZ4 == 4204
replace gdename = "himmelried" if gdename == "hlmmeiried" & PLZ4 == 4204
replace gdename = "himmelried" if gdename == "hlmmelried" & PLZ4 == 4204
replace gdename = "meltingen" if gdename == "meitingen" & PLZ4 == 4233
replace gdename = "meltingen" if gdename == "mettingen" & PLZ4 == 4233
replace gdename = "nunningen" if gdename == "nenziingen" & PLZ4 == 4208
replace gdename = "beinwil " if gdename == "oberbeinwil" & PLZ4 == 4229
replace gdename = "beinwil " if gdename == "unterbeinwil" & PLZ4 == 4229
replace gdename = "beinwil " if gdename == "unterbelnwil" & PLZ4 == 4229
replace gdename = "wahlen" if gdename == "wahiern" & PLZ4 == 4246
replace gdename = "wahlen" if gdename == "wahlen b." & PLZ4 == 4246
replace gdename = "wahlen" if gdename == "wahlen b. laufen" & PLZ4 == 4246
replace gdename = "wahlen" if gdename == "wahlen bei" & PLZ4 == 4246
replace gdename = "wahlen" if gdename == "wahlen bei laufen" & PLZ4 == 4246
replace gdename = "wahlen" if gdename == "wahlen bel laufen" & PLZ4 == 4246
replace gdename = "zullwil" if gdename == "zuliwil" & PLZ4 == 4234
replace gdename = "zullwil" if gdename == "zuliwll" & PLZ4 == 4234
replace gdename = "zullwil" if gdename == "zutiwii" & PLZ4 == 4234
replace gdename = "thal" if gdename == "altenrhein" & PLZ4 == 9423
replace gdename = "sennwald" if gdename == "fruemsen" & PLZ4 == 9467
replace gdename = "sennwald" if gdename == "fruemsen sennwald" & PLZ4 == 9467
replace gdename = "sennwald" if gdename == "haag" & PLZ4 == 9469
replace gdename = "sennwald" if gdename == "haag (rheintal)" & PLZ4 == 9469
replace gdename = "sennwald" if gdename == "haag (sennwaid)" & PLZ4 == 9469
replace gdename = "sennwald" if gdename == "haag (sennwald)" & PLZ4 == 9469
replace gdename = "sennwald" if gdename == "haag (sg)" & PLZ4 == 9469
replace gdename = "sennwald" if gdename == "haag sennwald" & PLZ4 == 9469
replace gdename = "sennwald" if gdename == "haag sg" & PLZ4 == 9469
replace gdename = "sennwald" if gdename == "sax" & PLZ4 == 9468
replace gdename = "wolfhalden" if gdename == "zeig" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zelg" & PLZ4 == 9427
replace gdename = "salenstein" if gdename == "fruthwilen" & PLZ4 == 8269
replace gdename = "salenstein" if gdename == "fruthwllen" & PLZ4 == 8269
replace gdename = "waeldi" if gdename == "hattenhausen/ lipperswii" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "hattenhausen/ lipperswil" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "hattenhausen/ upperswil" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "hattenhausen/lipperswil" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "hattenhausen/lipperswll" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "hefenhausen/ lipperswil" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "hefenhausen/ lipperswll" & PLZ4 == 8564
replace gdename = "wigoltingen" if gdename == "ilihart" & PLZ4 == 8556
replace gdename = "waeldi" if gdename == "lamperswil" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "lipperewil" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "lipperswi)" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "lipperswil" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "lipperswll" & PLZ4 == 8564
replace gdename = "wigoltingen" if gdename == "lllhart" & PLZ4 == 8556
replace gdename = "waeldi" if gdename == "llpperswii" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "llpperswll" & PLZ4 == 8564
replace gdename = "raperswilen" if gdename == "raperswilen tg" & PLZ4 == 8558
replace gdename = "raperswilen" if gdename == "raperswllen" & PLZ4 == 8558
replace gdename = "raperswilen" if gdename == "reperswilen" & PLZ4 == 8558
replace gdename = "lugano" if gdename == "cureglla" & PLZ4 == 6963
replace gdename = "baar" if gdename == "s hlbrugg dorf" & PLZ4 == 6340
replace gdename = "baar" if gdename == "sihibrugg" & PLZ4 == 6340
replace gdename = "baar" if gdename == "sihibrugg dorf" & PLZ4 == 6340
replace gdename = "baar" if gdename == "sihibrugg dort" & PLZ4 == 6340
replace gdename = "baar" if gdename == "sihlbrugg" & PLZ4 == 6340
replace gdename = "baar" if gdename == "sihlbrugg dorf" & PLZ4 == 6340
replace gdename = "baar" if gdename == "sihlbrugg dorf/baar" & PLZ4 == 6340
replace gdename = "baar" if gdename == "sihlbrugg dort" & PLZ4 == 6340
replace gdename = "baar" if gdename == "slhlbrugg" & PLZ4 == 6340
replace gdename = "baar" if gdename == "slhlbrugg dorf" & PLZ4 == 6340
replace gdename = "aire la ville" if gdename == "aire" & PLZ4 == 1288
replace gdename = "avusy" if gdename == "athenay" & PLZ4 == 1285
replace gdename = "avusy" if gdename == "athenaz" & PLZ4 == 1285
replace gdename = "avully" if gdename == "avuily" & PLZ4 == 1237
replace gdename = "collex bossy" if gdename == "bossy" & PLZ4 == 1239
replace gdename = "choulex" if gdename == "chevrier" & PLZ4 == 1244
replace gdename = "choulex" if gdename == "chevrier choulex" & PLZ4 == 1244
replace gdename = "giffers" if gdename == "chevrilles" & PLZ4 == 1735
replace gdename = "choulex" if gdename == "chouiex" & PLZ4 == 1244
replace gdename = "collex bossy" if gdename == "codex bossy" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "coilex" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "coilex bossy" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "coliex" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "coliex bossy" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "coljex" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "colles" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "collex" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "collex bos y" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "collex bossy ge" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "collez" & PLZ4 == 1239
replace gdename = "collex bossy" if gdename == "collez bossy" & PLZ4 == 1239
replace gdename = "dardagny" if gdename == "daradagny" & PLZ4 == 1283
replace gdename = "dardagny" if gdename == "la plaine" & PLZ4 == 1283
replace gdename = "presinge" if gdename == "presigne" & PLZ4 == 1243
replace gdename = "presinge" if gdename == "presinge/meinier" & PLZ4 == 1243
replace gdename = "presinge" if gdename == "presinoe" & PLZ4 == 1243
replace gdename = "puplinge" if gdename == "pupiinge" & PLZ4 == 1241
replace gdename = "puplinge" if gdename == "puplinges" & PLZ4 == 1241
replace gdename = "russin" if gdename == "rusaln" & PLZ4 == 1281
replace gdename = "avusy" if gdename == "sezegnin" & PLZ4 == 1285
replace gdename = "avusy" if gdename == "sezegnin/ge" & PLZ4 == 1285
replace gdename = "avusy" if gdename == "sezegnln" & PLZ4 == 1285
replace gdename = "soral" if gdename == "sorel" & PLZ4 == 1286
replace gdename = "collex bossy" if gdename == "vireloup collex ge" & PLZ4 == 1239
replace gdename = "arnex sur orbe" if gdename == "amex" & PLZ4 == 1321
replace gdename = "arnex sur orbe" if gdename == "amex sur orbe" & PLZ4 == 1321
replace gdename = "arnex sur orbe" if gdename == "annex s/orbe" & PLZ4 == 1321
replace gdename = "arnex sur orbe" if gdename == "arnex" & PLZ4 == 1321
replace gdename = "arnex sur orbe" if gdename == "arnex s/orbe" & PLZ4 == 1321
replace gdename = "bretonnieres" if gdename == "bretonnleres" & PLZ4 == 1329
replace gdename = "l'isle" if gdename == "coudre l isle" & PLZ4 == 1148
replace gdename = "croy" if gdename == "croy/romainmdtier" & PLZ4 == 1322
replace gdename = "croy" if gdename == "croy/romainmã³tier" & PLZ4 == 1322
replace gdename = "daillens" if gdename == "daiilens" & PLZ4 == 1306
replace gdename = "daillens" if gdename == "dailiens" & PLZ4 == 1306
replace gdename = "daillens" if gdename == "dalliens" & PLZ4 == 1306
replace gdename = "daillens" if gdename == "dallions" & PLZ4 == 1306
replace gdename = "daillens" if gdename == "dalllens" & PLZ4 == 1306
replace gdename = "eclepens" if gdename == "eciepens" & PLZ4 == 1312
replace gdename = "juriens" if gdename == "jurions" & PLZ4 == 1326
replace gdename = "juriens" if gdename == "jurlens" & PLZ4 == 1326
replace gdename = "l'isle" if gdename == "la coudre" & PLZ4 == 1148
replace gdename = "l'isle" if gdename == "la coudre (l jsle)" & PLZ4 == 1148
replace gdename = "l'isle" if gdename == "la coudre l isle" & PLZ4 == 1148
replace gdename = "l'isle" if gdename == "la coudre l jsle" & PLZ4 == 1148
replace gdename = "l'isle" if gdename == "la coudre/l isie" & PLZ4 == 1148
replace gdename = "moiry" if gdename == "moiry vd" & PLZ4 == 1148
replace gdename = "moiry" if gdename == "molry" & PLZ4 == 1148
replace gdename = "penthaz" if gdename == "penthas vd" & PLZ4 == 1303
replace gdename = "pompaples" if gdename == "pomgaples" & PLZ4 == 1318
replace gdename = "romainmotier envy" if gdename == "romainmbtier" & PLZ4 == 1323
replace gdename = "romainmotier envy" if gdename == "romainmoetier" & PLZ4 == 1323
replace gdename = "romainmotier envy" if gdename == "romainmoetler" & PLZ4 == 1323
replace gdename = "romainmotier envy" if gdename == "romainmotier" & PLZ4 == 1323
replace gdename = "romainmotier envy" if gdename == "romainmotler" & PLZ4 == 1323
replace gdename = "romainmotier envy" if gdename == "romainmotter" & PLZ4 == 1323
replace gdename = "romainmotier envy" if gdename == "romainmã³tier" & PLZ4 == 1323
replace gdename = "romainmotier envy" if gdename == "romainmã³tler" & PLZ4 == 1323
replace gdename = "romainmotier envy" if gdename == "romalnmotler" & PLZ4 == 1323
replace gdename = "romainmotier envy" if gdename == "romalnmã³tier" & PLZ4 == 1323
replace gdename = "chavannes le veyron" if gdename == "st. denis" & PLZ4 == 1148
replace gdename = "vaulion" if gdename == "vauiion" & PLZ4 == 1325
replace gdename = "vaulion" if gdename == "vaullon" & PLZ4 == 1325
replace gdename = "vaulion" if gdename == "vaution" & PLZ4 == 1325
replace gdename = "vufflens la ville" if gdename == "vuffoens la ville" & PLZ4 == 1302
replace gdename = "dulliken" if gdename == "dulllken" & PLZ4 == 4657
replace gdename = "murgenthal" if gdename == "murgenthal/riken" & PLZ4 == 4853
replace gdename = "murgenthal" if gdename == "riken" & PLZ4 == 4853
replace gdename = "murgenthal" if gdename == "riken (ag)" & PLZ4 == 4853
replace gdename = "murgenthal" if gdename == "riken ag" & PLZ4 == 4853
replace gdename = "murgenthal" if gdename == "riken bei murgenthal" & PLZ4 == 4853
replace gdename = "murgenthal" if gdename == "riken murgenthal" & PLZ4 == 4853
replace gdename = "murgenthal" if gdename == "riken/ murgenthal" & PLZ4 == 4853
replace gdename = "murgenthal" if gdename == "riken/murgenthal" & PLZ4 == 4853
replace gdename = "murgenthal" if gdename == "rlken" & PLZ4 == 4853
replace gdename = "pontresina" if gdename == "bernina suof" & PLZ4 == 7504
replace gdename = "brusio" if gdename == "campascio" & PLZ4 == 7748
replace gdename = "brusio" if gdename == "campascio brusio" & PLZ4 == 7748
replace gdename = "brusio" if gdename == "campasclo" & PLZ4 == 7748
replace gdename = "brusio" if gdename == "campasclo brusio" & PLZ4 == 7748
replace gdename = "poschiavo" if gdename == "li curt" & PLZ4 == 7745
replace gdename = "poschiavo" if gdename == "ospizio bernina" & PLZ4 == 7710
replace gdename = "poschiavo" if gdename == "prada" & PLZ4 == 7745
replace gdename = "poschiavo" if gdename == "prada annunziata" & PLZ4 == 7745
replace gdename = "poschiavo" if gdename == "s. carlo" & PLZ4 == 7741
replace gdename = "poschiavo" if gdename == "s. carlo (poschiavo)" & PLZ4 == 7741
replace gdename = "poschiavo" if gdename == "s.antonio" & PLZ4 == 7745
replace gdename = "poschiavo" if gdename == "s.antonio/poschiavo" & PLZ4 == 7745
replace gdename = "poschiavo" if gdename == "s.carlo" & PLZ4 == 7741
replace gdename = "poschiavo" if gdename == "s.carlo (poschiavo)" & PLZ4 == 7741
replace gdename = "gurmels" if gdename == "grossguschelmuth" & PLZ4 == 1792
replace gdename = "mont vully" if gdename == "lugnorre" & PLZ4 == 1789
replace gdename = "mont vully" if gdename == "lugnorre fr" & PLZ4 == 1789
replace gdename = "mont vully" if gdename == "lugnorre haut vully" & PLZ4 == 1789
replace gdename = "mont vully" if gdename == "lugnorre vuily le haut" & PLZ4 == 1789
replace gdename = "mont vully" if gdename == "lugnorre vully le haut" & PLZ4 == 1789
replace gdename = "mont vully" if gdename == "lugnorre/haut vuliy" & PLZ4 == 1789
replace gdename = "mont vully" if gdename == "lugnorre/haut vully" & PLZ4 == 1789
replace gdename = "mont vully" if gdename == "lugnorre/motier" & PLZ4 == 1787
replace gdename = "muenchenwiler" if gdename == "muenchenwiier" & PLZ4 == 1797
replace gdename = "muenchenwiler" if gdename == "muenchenwller" & PLZ4 == 1797
replace gdename = "muenchenwiler" if gdename == "muenchwiler" & PLZ4 == 1797
replace gdename = "mont vully" if gdename == "nieder wistenlach" & PLZ4 == 1786
replace gdename = "mont vully" if gdename == "nieder wlstenlach" & PLZ4 == 1786
replace gdename = "mont vully" if gdename == "nleder wistenlach" & PLZ4 == 1786
replace gdename = "mont vully" if gdename == "praz" & PLZ4 == 1788
replace gdename = "mont vully" if gdename == "praz (vuily)" & PLZ4 == 1788
replace gdename = "mont vully" if gdename == "praz (vully)" & PLZ4 == 1788
replace gdename = "mont vully" if gdename == "praz vully" & PLZ4 == 1788
replace gdename = "morat" if gdename == "saivenach" & PLZ4 == 1794
replace gdename = "binningen" if gdename == "biningen" & PLZ4 == 4102
replace gdename = "hofstetten flueh" if gdename == "fiueh" & PLZ4 == 4114
replace gdename = "hofstetten flueh" if gdename == "fiueh so" & PLZ4 == 4114
replace gdename = "hofstetten flueh" if gdename == "flueh" & PLZ4 == 4114
replace gdename = "hofstetten flueh" if gdename == "flueh so" & PLZ4 == 4114
replace gdename = "hofstetten flueh" if gdename == "flueh/so" & PLZ4 == 4114
replace gdename = "hofstetten flueh" if gdename == "fluh" & PLZ4 == 4114
replace gdename = "hofstetten flueh" if gdename == "fluh so" & PLZ4 == 4114
replace gdename = "tujetsch" if gdename == "camischolas" & PLZ4 == 7187
replace gdename = "medel" if gdename == "curaglia" & PLZ4 == 7184
replace gdename = "medel" if gdename == "curaglla" & PLZ4 == 7184
replace gdename = "tujetsch" if gdename == "rueras" & PLZ4 == 7189
replace gdename = "tujetsch" if gdename == "rueras tujetsch" & PLZ4 == 7189
replace gdename = "tujetsch" if gdename == "rueras/tujetsch" & PLZ4 == 7189
replace gdename = "tujetsch" if gdename == "ruerasftujetsch" & PLZ4 == 7189
replace gdename = "flims" if gdename == "segnes" & PLZ4 == 7017
replace gdename = "tujetsch" if gdename == "selva" & PLZ4 == 7189
replace gdename = "tujetsch" if gdename == "selva (tavetsch)" & PLZ4 == 7189
replace gdename = "val muestair" if gdename == "st. maria" & PLZ4 == 7536
replace gdename = "lumnezia" if gdename == "(gels vattiz" & PLZ4 == 7146
replace gdename = "obersaxen mundaun" if gdename == "affeier" & PLZ4 == 7133
replace gdename = "obersaxen mundaun" if gdename == "affeler" & PLZ4 == 7133
replace gdename = "falera" if gdename == "faiera" & PLZ4 == 7153
replace gdename = "falera" if gdename == "fellers" & PLZ4 == 7153
replace gdename = "obersaxen mundaun" if gdename == "fiond" & PLZ4 == 7137
replace gdename = "lumnezia" if gdename == "igels" & PLZ4 == 7146
replace gdename = "lumnezia" if gdename == "igels vattiz" & PLZ4 == 7146
replace gdename = "laax" if gdename == "laax b. ilanz" & PLZ4 == 7031
replace gdename = "laax" if gdename == "laax b. jlanz" & PLZ4 == 7031
replace gdename = "ilanz" if gdename == "ladlr" & PLZ4 == 7155
replace gdename = "lumnezia" if gdename == "lumbreln" & PLZ4 == 7148
replace gdename = "silenen" if gdename == "luren" & PLZ4 == 6476
replace gdename = "obersaxen mundaun" if gdename == "obersaken" & PLZ4 == 7134
replace gdename = "obersaxen mundaun" if gdename == "obersaxen tusen" & PLZ4 == 7134
replace gdename = "obersaxen mundaun" if gdename == "obesaxen" & PLZ4 == 7134
replace gdename = "sagogn" if gdename == "sagoon" & PLZ4 == 7152
replace gdename = "schluein" if gdename == "schieuis" & PLZ4 == 7151
replace gdename = "schluein" if gdename == "schiuein" & PLZ4 == 7151
replace gdename = "schluein" if gdename == "schleuis" & PLZ4 == 7151
replace gdename = "schluein" if gdename == "schleuls" & PLZ4 == 7151
replace gdename = "schluein" if gdename == "schleunis" & PLZ4 == 7151
replace gdename = "schluein" if gdename == "schlueln" & PLZ4 == 7151
replace gdename = "ilanz" if gdename == "sergein" & PLZ4 == 7127
replace gdename = "ilanz" if gdename == "sergeln" & PLZ4 == 7127
replace gdename = "ilanz" if gdename == "sevgeln" & PLZ4 == 7127
replace gdename = "sagogn" if gdename == "sgogn" & PLZ4 == 7152
replace gdename = "obersaxen mundaun" if gdename == "surcuoim" & PLZ4 == 7138
replace gdename = "uors" if gdename == "uors" & PLZ4 == 7114
replace gdename = "obersaxen mundaun" if gdename == "vaiata obersaxen" & PLZ4 == 7138
replace gdename = "obersaxen mundaun" if gdename == "valata obersaxen" & PLZ4 == 7138
replace gdename = "obersaxen mundaun" if gdename == "valataobersaxen" & PLZ4 == 7138
replace gdename = "lumnezia" if gdename == "vattiz" & PLZ4 == 7146
replace gdename = "lumnezia" if gdename == "vigens" & PLZ4 == 7147
replace gdename = "lumnezia" if gdename == "villa" & PLZ4 == 7144
replace gdename = "lumnezia" if gdename == "villa gr" & PLZ4 == 7144
replace gdename = "serravalle" if gdename == "maivaglia chiesa" & PLZ4 == 6713
replace gdename = "serravalle" if gdename == "malvaglia chiesa" & PLZ4 == 6713
replace gdename = "serravalle" if gdename == "malvaglia chiesu" & PLZ4 == 6713
replace gdename = "serravalle" if gdename == "malvaglia chlesa" & PLZ4 == 6713
replace gdename = "serravalle" if gdename == "malvaglla chiesa" & PLZ4 == 6713
replace gdename = "biasca" if gdename == "pontirone" & PLZ4 == 6710
replace gdename = "vully les lacs" if gdename == "beilerive" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "chaberey" & PLZ4 == 1589
replace gdename = "vully les lacs" if gdename == "chambery" & PLZ4 == 1589
replace gdename = "saint aubin (fr)" if gdename == "les friques" & PLZ4 == 1584
replace gdename = "vully les lacs" if gdename == "saiavaux" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "salavaux" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "salavaux bellerive" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "salavaux vd" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "vailamand dessous" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "valiamand" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "valiamand dessous" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "valiamand dessus" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "vallamand dessous" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "vallamand dessu3" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "vallamand dessus" & PLZ4 == 1585
replace gdename = "vully les lacs" if gdename == "vallamand.dessous" & PLZ4 == 1585
replace gdename = "birsfelden" if gdename == "blirsfelden" & PLZ4 == 4127
replace gdename = "hofstetten flueh" if gdename == "flueh" & PLZ4 == 4114
replace gdename = "allschwil" if gdename == "neu allschwil" & PLZ4 == 4123
replace gdename = "allschwil" if gdename == "neuailschwil" & PLZ4 == 4123
replace gdename = "allschwil" if gdename == "neualischwil" & PLZ4 == 4123
replace gdename = "allschwil" if gdename == "neualischwll" & PLZ4 == 4123
replace gdename = "allschwil" if gdename == "neuallschwii" & PLZ4 == 4123
replace gdename = "allschwil" if gdename == "neuallschwil" & PLZ4 == 4123
replace gdename = "eggiwil" if gdename == "aeschau" & PLZ4 == 3536
replace gdename = "trub" if gdename == "fankhaus (trub)" & PLZ4 == 3556
replace gdename = "trub" if gdename == "fankhaus trub" & PLZ4 == 3556
replace gdename = "trub" if gdename == "fankhaus(trub)" & PLZ4 == 3556
replace gdename = "langnau (be)" if gdename == "gohl" & PLZ4 == 3550
replace gdename = "konolfingen" if gdename == "gysenstein" & PLZ4 == 3510
replace gdename = "konolfingen" if gdename == "gysenstein konolfingen" & PLZ4 == 3510
replace gdename = "konolfingen" if gdename == "gysensteln" & PLZ4 == 3510
replace gdename = "oberthal" if gdename == "oberthai" & PLZ4 == 3531
replace gdename = "muensingen" if gdename == "taegertschl" & PLZ4 == 3111
replace gdename = "muensingen" if gdename == "tagertschi" & PLZ4 == 3111
replace gdename = "muensingen" if gdename == "triimsteln" & PLZ4 == 3083
replace gdename = "muensingen" if gdename == "trimatein" & PLZ4 == 3083
replace gdename = "muensingen" if gdename == "trimstcin" & PLZ4 == 3083
replace gdename = "muensingen" if gdename == "trimstein" & PLZ4 == 3083
replace gdename = "muensingen" if gdename == "trimsteln" & PLZ4 == 3083
replace gdename = "muensingen" if gdename == "trlmstein" & PLZ4 == 3083
replace gdename = "wolfhalden" if gdename == "wolf halden zeig" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "wolfhalden zeig" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "wolfhalden zelg" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "wollhalden zeig" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zeig" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zeig (wolf halden)" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zeig (wolfhaiden)" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zeig (wolfhalden)" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zeig (wollhalden)" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zeig wolfhalden" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zeig/woifhaiden" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zelg" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zelg (woifhalden)" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zelg (wolf halden)" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zelg (wolfhaiden)" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zelg (wolfhalden)" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zelg (wollhalden)" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zelg wolfhaiden" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zelg wolfhalden" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zelg/wolfhalden" & PLZ4 == 9427
replace gdename = "wolfhalden" if gdename == "zolg wolfhalden" & PLZ4 == 9427
replace gdename = "walenstadt" if gdename == "berschis" & PLZ4 == 8880
replace gdename = "walenstadt" if gdename == "berschis (walenstadt)" & PLZ4 == 8880
replace gdename = "walenstadt" if gdename == "berschis walenstadt" & PLZ4 == 8880
replace gdename = "walenstadt" if gdename == "berschis/walenstadt" & PLZ4 == 8880
replace gdename = "walenstadt" if gdename == "berschls" & PLZ4 == 8880
replace gdename = "flums" if gdename == "grossberg fiums" & PLZ4 == 8890
replace gdename = "flums" if gdename == "grossberg flums" & PLZ4 == 8890
replace gdename = "flums" if gdename == "tannenbodenalp" & PLZ4 == 8890
replace gdename = "flums" if gdename == "tannenbodenalp (flums)" & PLZ4 == 8890
replace gdename = "flums" if gdename == "tannenheim" & PLZ4 == 8890
replace gdename = "flums" if gdename == "tannenheim fiums" & PLZ4 == 8890
replace gdename = "flums" if gdename == "tannenheim flume" & PLZ4 == 8890
replace gdename = "flums" if gdename == "tannenheim flums" & PLZ4 == 8890
replace gdename = "flums" if gdename == "tannenheim plums" & PLZ4 == 8890
replace gdename = "flums" if gdename == "tannenhelm flums" & PLZ4 == 8890
replace gdename = "oron" if gdename == "chaetillens" & PLZ4 == 1610
replace gdename = "oron" if gdename == "chatiilens" & PLZ4 == 1610
replace gdename = "oron" if gdename == "chatiliens" & PLZ4 == 1610
replace gdename = "oron" if gdename == "chatlllens" & PLZ4 == 1610
replace gdename = "ecublens" if gdename == "ecubiens" & PLZ4 == 1024
replace gdename = "valbroye" if gdename == "les treize cantons" & PLZ4 == 1525
replace gdename = "valbroye" if gdename == "mamand" & PLZ4 == 1524
replace gdename = "valbroye" if gdename == "marmand" & PLZ4 == 1524
replace gdename = "valbroye" if gdename == "marriand" & PLZ4 == 1524
replace gdename = "oron" if gdename == "palesieux vd" & PLZ4 == 1607
replace gdename = "oron" if gdename == "palezieux vd" & PLZ4 == 1607
replace gdename = "oron" if gdename == "palezieux village" & PLZ4 == 1607
replace gdename = "oron" if gdename == "palezleux vd" & PLZ4 == 1607
replace gdename = "oron" if gdename == "palezleux village" & PLZ4 == 1607
replace gdename = "valbroye" if gdename == "seigneuz" & PLZ4 == 1525
replace gdename = "valbroye" if gdename == "soigneux" & PLZ4 == 1525
replace gdename = "oron" if gdename == "tavernes" & PLZ4 == 1607
replace gdename = "valbroye" if gdename == "treize cantcns" & PLZ4 == 1525
replace gdename = "valbroye" if gdename == "treize canton" & PLZ4 == 1525
replace gdename = "valbroye" if gdename == "treize cantons" & PLZ4 == 1525
replace gdename = "oron" if gdename == "vuibroye chaetillens" & PLZ4 == 1610
replace gdename = "oron" if gdename == "vuibroye chatiilens" & PLZ4 == 1610
replace gdename = "oron" if gdename == "vuibroye chatiliens" & PLZ4 == 1610
replace gdename = "oron" if gdename == "vuibroye chatillens" & PLZ4 == 1610
replace gdename = "oron" if gdename == "vulbroye" & PLZ4 == 1610
replace gdename = "oron" if gdename == "vulbroye chatlllens" & PLZ4 == 1610
replace gdename = "corbeyrier" if gdename == "corbeyraer" & PLZ4 == 1856
replace gdename = "corbeyrier" if gdename == "corbeyrler" & PLZ4 == 1856
replace gdename = "ormont dessous" if gdename == "forclaz" & PLZ4 == 1866
replace gdename = "ollon" if gdename == "huemoz" & PLZ4 == 1884
replace gdename = "ollon" if gdename == "huemoz sur ollon" & PLZ4 == 1884
replace gdename = "ollon" if gdename == "huemoz sur otton" & PLZ4 == 1884
replace gdename = "ormont dessous" if gdename == "la comballaz" & PLZ4 == 1862
replace gdename = "ormont dessous" if gdename == "la forclaz" & PLZ4 == 1866
replace gdename = "ormont dessous" if gdename == "la forclaz/evolene" & PLZ4 == 1866
replace gdename = "ormont dessous" if gdename == "les mosses" & PLZ4 == 1862
replace gdename = "ormont dessous" if gdename == "les mosses (ormont dessous)" & PLZ4 == 1862
replace gdename = "ormont dessous" if gdename == "mosses/ ormont dessous" & PLZ4 == 1862
replace gdename = "ormont dessous" if gdename == "mosses/ormont dessous" & PLZ4 == 1862
replace gdename = "ollon" if gdename == "panex" & PLZ4 == 1867
replace gdename = "ollon" if gdename == "panex ollon" & PLZ4 == 1867
replace gdename = "ollon" if gdename == "panex sur ollon" & PLZ4 == 1867
replace gdename = "ollon" if gdename == "panex sur olson" & PLZ4 == 1867
replace gdename = "ollon" if gdename == "panez ollon" & PLZ4 == 1867
replace gdename = "castel san pietro" if gdename == "cortegiia" & PLZ4 == 6874
replace gdename = "castel san pietro" if gdename == "corteglia" & PLZ4 == 6874
replace gdename = "castel san pietro" if gdename == "corteglia die castel san pietro" & PLZ4 == 6874
replace gdename = "castel san pietro" if gdename == "corteglia/ castel s. pietro" & PLZ4 == 6874
replace gdename = "castel san pietro" if gdename == "corteglia/castel s. pietro" & PLZ4 == 6874
replace gdename = "castel san pietro" if gdename == "corteglla" & PLZ4 == 6874
replace gdename = "mendrisio" if gdename == "mendrisio sortite" & PLZ4 == 6850
replace gdename = "gommiswald" if gdename == "emetschwil" & PLZ4 == 8725
replace gdename = "gommiswald" if gdename == "ernetschwll" & PLZ4 == 8725
replace gdename = "gommiswald" if gdename == "gebertangen" & PLZ4 == 8725
replace gdename = "gommiswald" if gdename == "gebertingen" & PLZ4 == 8725
replace gdename = "gommiswald" if gdename == "gomiswald" & PLZ4 == 8737
replace gdename = "lachen (sz)" if gdename == "lachen sz" & PLZ4 == 8853
replace gdename = "gommiswald" if gdename == "ricken" & PLZ4 == 8726
replace gdename = "gommiswald" if gdename == "ricken (sg)" & PLZ4 == 8726
replace gdename = "gommiswald" if gdename == "rieden sg" & PLZ4 == 8739
replace gdename = "gommiswald" if gdename == "rleden" & PLZ4 == 8739
replace gdename = "gommiswald" if gdename == "uetiiberg" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetliberg" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetliburg" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetliburg b.gommiswald" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetliburg bei gommiswald" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetliburg bel gommiswald" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetliburg gommiswaid" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetliburg gommiswald" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetllberg" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetllburg" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetllburg gommiswald" & PLZ4 == 8738
replace gdename = "gommiswald" if gdename == "uetlã¬burg" & PLZ4 == 8738
replace gdename = "puidoux" if gdename == "pueldoux gare" & PLZ4 == 1070
replace gdename = "puidoux" if gdename == "puidoux gare" & PLZ4 == 1070
replace gdename = "puidoux" if gdename == "puidoux village" & PLZ4 == 1070
replace gdename = "puidoux" if gdename == "puldoux" & PLZ4 == 1070
replace gdename = "puidoux" if gdename == "puldoux gare" & PLZ4 == 1070
replace gdename = "puidoux" if gdename == "puldoux village" & PLZ4 == 1070
replace gdename = "oberegg" if gdename == "bueriswilen" & PLZ4 == 9413
replace gdename = "oberegg" if gdename == "bueriswilen oberegg" & PLZ4 == 9413
replace gdename = "oberegg" if gdename == "bueriswilen/oberegg" & PLZ4 == 9413
replace gdename = "oberegg" if gdename == "bueriswilenioberegg" & PLZ4 == 9413
replace gdename = "oberegg" if gdename == "bueriswllen oberegg" & PLZ4 == 9413
replace gdename = "oberegg" if gdename == "buerlswilen/oberegg" & PLZ4 == 9413
replace gdename = "walzenhausen" if gdename == "platz" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "platz (ar)" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "platz (waizenhausen)" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "platz (walzenhausen)" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "platz waizenhausen" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "platz walzenhausen" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "waizenhausen platz" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "walzenhausen piatz" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "walzenhausen platz" & PLZ4 == 9428
replace gdename = "walliswil bei niederbipp" if gdename == "wailiswil bel niederbipp" & PLZ4 == 3380
replace gdename = "walliswil bei wangen" if gdename == "waliiswil" & PLZ4 == 3377
replace gdename = "walliswil bei niederbipp" if gdename == "waliiswil bei niederbipp" & PLZ4 == 3380
replace gdename = "walliswil bei wangen" if gdename == "walliswal" & PLZ4 == 3377
replace gdename = "walliswil bei niederbipp" if gdename == "walliswal bei niederbipp" & PLZ4 == 3380
replace gdename = "walliswil bei niederbipp" if gdename == "walliswal bel niederbipp" & PLZ4 == 3380
replace gdename = "walliswil bei wangen" if gdename == "walliswii" & PLZ4 == 3377
replace gdename = "walliswil bei wangen" if gdename == "walliswil" & PLZ4 == 3377
replace gdename = "walliswil bei niederbipp" if gdename == "walliswil b. niederbipp" & PLZ4 == 3380
replace gdename = "walliswil bei wangen" if gdename == "walliswil b. wangen" & PLZ4 == 3377
replace gdename = "walliswil bei niederbipp" if gdename == "walliswil bel niederbipp" & PLZ4 == 3380
replace gdename = "walliswil bei wangen" if gdename == "walliswll" & PLZ4 == 3377
replace gdename = "walliswil bei niederbipp" if gdename == "walltswil bei niederbipp" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a aare" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a. a" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a. a." & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a. aare" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a. albas" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a. albis" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a. d" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a. d. aare" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a.a" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a.aare" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a.d. aare" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a.d./a." & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a.d.a" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen a.d.aare" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangen en der aare" & PLZ4 == 3380
replace gdename = "wangen an der aare" if gdename == "wangena.a" & PLZ4 == 3380
replace gdename = "niederwil" if gdename == "nesseinbach" & PLZ4 == 5524
replace gdename = "niederwil" if gdename == "nesseinbach ag" & PLZ4 == 5524
replace gdename = "niederwil" if gdename == "nesselnbach" & PLZ4 == 5524
replace gdename = "niederwil" if gdename == "nesselnbach /niederwil ag" & PLZ4 == 5524
replace gdename = "niederwil" if gdename == "nesselnbach ag" & PLZ4 == 5524
replace gdename = "niederwil" if gdename == "nesselnbach/niederwal ag" & PLZ4 == 5524
replace gdename = "niederwil" if gdename == "nesselnbach/niederwil" & PLZ4 == 5524
replace gdename = "niederwil" if gdename == "nesselnbach/niederwil ag" & PLZ4 == 5524
replace gdename = "niederwil" if gdename == "nesselnbach/niederwll" & PLZ4 == 5524
replace gdename = "niederwil" if gdename == "niederwal ag" & PLZ4 == 5524
replace gdename = "vernier" if gdename == "chaetelaine" & PLZ4 == 1219
replace gdename = "vernier" if gdename == "chaetelaine vernier" & PLZ4 == 1219
replace gdename = "vernier" if gdename == "chatelaine" & PLZ4 == 1219
replace gdename = "vernier" if gdename == "chatelaine vernier" & PLZ4 == 1219
replace gdename = "le grand saconnex" if gdename == "grand saconnex" & PLZ4 == 1218
replace gdename = "vernier" if gdename == "la chatelaine" & PLZ4 == 1219
replace gdename = "val de ruz" if gdename == "cerniez" & PLZ4 == 2053
replace gdename = "chatonnaye" if gdename == "chaetonnaye" & PLZ4 == 1553
replace gdename = "chatonnaye" if gdename == "chatennaye" & PLZ4 == 1553
replace gdename = "chevroux" if gdename == "chevraux" & PLZ4 == 1545
replace gdename = "valbroye" if gdename == "combremont ie petit" & PLZ4 == 1536
replace gdename = "valbroye" if gdename == "combremont le petlt" & PLZ4 == 1536
replace gdename = "dompierre" if gdename == "domplerre" & PLZ4 == 1682
replace gdename = "fetigny" if gdename == "etigny" & PLZ4 == 1532
replace gdename = "fetigny" if gdename == "fetlgny" & PLZ4 == 1532
replace gdename = "gletterens" if gdename == "gietterens" & PLZ4 == 1544
replace gdename = "gletterens" if gdename == "gletternens" & PLZ4 == 1544
replace gdename = "menieres" if gdename == "memeres" & PLZ4 == 1533
replace gdename = "menieres" if gdename == "menleres" & PLZ4 == 1533
replace gdename = "estavayer" if gdename == "rueyres le pres" & PLZ4 == 1542
replace gdename = "valbroye" if gdename == "sasse)" & PLZ4 == 1534
replace gdename = "valbroye" if gdename == "sassei" & PLZ4 == 1534
replace gdename = "valbroye" if gdename == "sassel vd" & PLZ4 == 1534
replace gdename = "valbroye" if gdename == "sasses" & PLZ4 == 1534
replace gdename = "valbroye" if gdename == "sessel" & PLZ4 == 1534
replace gdename = "valbroye" if gdename == "sessel vd" & PLZ4 == 1534
replace gdename = "payerne" if gdename == "vers .chez perrin" & PLZ4 == 1551
replace gdename = "payerne" if gdename == "vers chez perrin" & PLZ4 == 1551
replace gdename = "payerne" if gdename == "vers chez perrin/payeme" & PLZ4 == 1551
replace gdename = "payerne" if gdename == "vers chez perrin/payerne" & PLZ4 == 1551
replace gdename = "payerne" if gdename == "verschez payern" & PLZ4 == 1551
replace gdename = "payerne" if gdename == "verschez perrin/payeme" & PLZ4 == 1551
replace gdename = "payerne" if gdename == "verschez perrin/payerne" & PLZ4 == 1551
replace gdename = "villarzel" if gdename == "viliarzel" & PLZ4 == 1555
replace gdename = "villarzel" if gdename == "villarzell" & PLZ4 == 1555
replace gdename = "villarzel" if gdename == "vlllarzel" & PLZ4 == 1555
replace gdename = "villarzel" if gdename == "vlllarzell" & PLZ4 == 1555
replace gdename = "trimmis" if gdename == "trimmas dorf" & PLZ4 == 7203
replace gdename = "trimmis" if gdename == "trimmas dort" & PLZ4 == 7203
replace gdename = "trimmis" if gdename == "trimmfis dorf" & PLZ4 == 7203
replace gdename = "trimmis" if gdename == "trimmis dorf" & PLZ4 == 7203
replace gdename = "trimmis" if gdename == "trimmis dorl" & PLZ4 == 7203
replace gdename = "trimmis" if gdename == "trimmis dort" & PLZ4 == 7203
replace gdename = "trimmis" if gdename == "trlmmis dorf" & PLZ4 == 7203
replace gdename = "muotathal" if gdename == "hinterthai" & PLZ4 == 6436
replace gdename = "muotathal" if gdename == "hinterthal" & PLZ4 == 6436
replace gdename = "muotathal" if gdename == "hinterthal /muotathal" & PLZ4 == 6436
replace gdename = "muotathal" if gdename == "hinterthal muotathai" & PLZ4 == 6436
replace gdename = "muotathal" if gdename == "hinterthal muotathal" & PLZ4 == 6436
replace gdename = "muotathal" if gdename == "hinterthal/muotathal" & PLZ4 == 6436
replace gdename = "muotathal" if gdename == "hlnterthal" & PLZ4 == 6436
replace gdename = "wuerenlos" if gdename == "oetlikon" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "warenlos" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "werenlos" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "wueenlos" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "wuer" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "wuerenios" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "wuerenloa" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "wuerenlos ag" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "wuerenlos/ag" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "wuerentas" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "wuerentos" & PLZ4 == 5436
replace gdename = "wuerenlos" if gdename == "wurenlos" & PLZ4 == 5436
replace gdename = "bourg en lavaux" if gdename == "araen" & PLZ4 == 1096
replace gdename = "bourg en lavaux" if gdename == "aran" & PLZ4 == 1096
replace gdename = "bourg en lavaux" if gdename == "aran grandvaux" & PLZ4 == 1091
replace gdename = "bourg en lavaux" if gdename == "aran sur villette" & PLZ4 == 1096
replace gdename = "bourg en lavaux" if gdename == "aran villette" & PLZ4 == 1096
replace gdename = "bourg en lavaux" if gdename == "aren villette" & PLZ4 == 1096
replace gdename = "bourg en lavaux" if gdename == "arlin" & PLZ4 == 1096
replace gdename = "bourg en lavaux" if gdename == "aron villette" & PLZ4 == 1096
replace gdename = "bourg en lavaux" if gdename == "chenaux" & PLZ4 == 1091
replace gdename = "bourg en lavaux" if gdename == "chenaux grandvaux" & PLZ4 == 1091
replace gdename = "bourg en lavaux" if gdename == "granduaux" & PLZ4 == 1091
replace gdename = "bourg en lavaux" if gdename == "grandvaux vd" & PLZ4 == 1091
replace gdename = "bourg en lavaux" if gdename == "grandveau" & PLZ4 == 1091
replace gdename = "bourg en lavaux" if gdename == "grardvaux" & PLZ4 == 1091
replace gdename = "le mont sur lausanne" if gdename == "mont sur lausanne" & PLZ4 == 1052
replace gdename = "amlikon bissegg" if gdename == "amilkon" & PLZ4 == 8514
replace gdename = "amlikon bissegg" if gdename == "amlikon" & PLZ4 == 8514
replace gdename = "amlikon bissegg" if gdename == "amlikon tg" & PLZ4 == 8514
replace gdename = "amlikon bissegg" if gdename == "amllkon" & PLZ4 == 8514
replace gdename = "amlikon bissegg" if gdename == "amlã¬kon" & PLZ4 == 8514
replace gdename = "affeltragen" if gdename == "maltbach" & PLZ4 == 9556
replace gdename = "val muestair" if gdename == "muestalr" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "mustair" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "mãºstair" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "st. maria" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "st. maria 1. m." & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "st. maria i. m" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta maria" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta maria i.m" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta maria im muenstertal" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta maria im munstertal" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta. maria" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta. maria i. m" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta. maria i. m." & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta. maria im muenstertai" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta. maria im muenstertal" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta. maria im munstertal" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta. marla" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta. marla i. m" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "sta. marla im muenstertal" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "ste. maria im muenstertal" & PLZ4 == 7537
replace gdename = "val muestair" if gdename == "tschier" & PLZ4 == 7532
replace gdename = "val muestair" if gdename == "tschiery" & PLZ4 == 7532
replace gdename = "hochdorf" if gdename == "urswe" & PLZ4 == 6280
replace gdename = "hochdorf" if gdename == "urswii/hochdorf" & PLZ4 == 6280
replace gdename = "hochdorf" if gdename == "urswil" & PLZ4 == 6280
replace gdename = "hochdorf" if gdename == "urswil/hochdorf" & PLZ4 == 6280
replace gdename = "hochdorf" if gdename == "urswll" & PLZ4 == 6280
replace gdename = "sirnach" if gdename == "busawll" & PLZ4 == 8370
replace gdename = "sirnach" if gdename == "busswal" & PLZ4 == 8370
replace gdename = "sirnach" if gdename == "busswal tg" & PLZ4 == 8370
replace gdename = "sirnach" if gdename == "busswal tg/sirnach" & PLZ4 == 8370
replace gdename = "sirnach" if gdename == "busswii tg" & PLZ4 == 8370
replace gdename = "sirnach" if gdename == "busswil tg" & PLZ4 == 8370
replace gdename = "sirnach" if gdename == "busswil tg/egnach" & PLZ4 == 8370
replace gdename = "sirnach" if gdename == "busswil tg/sirnach" & PLZ4 == 8370
replace gdename = "sirnach" if gdename == "busswll tg" & PLZ4 == 8370
replace gdename = "sirnach" if gdename == "buswal tg" & PLZ4 == 8370
replace gdename = "hittnau" if gdename == "duerstelen" & PLZ4 == 8335
replace gdename = "hittnau" if gdename == "duersteren" & PLZ4 == 8335
replace gdename = "hittnau" if gdename == "duersteten" & PLZ4 == 8335
replace gdename = "hittnau" if gdename == "ober hittnau" & PLZ4 == 8335
replace gdename = "hittnau" if gdename == "oberhittnau" & PLZ4 == 8335
replace gdename = "hittnau" if gdename == "oberhlttnau" & PLZ4 == 8335
replace gdename = "pfaeffikon" if gdename == "pfaeffikon zh" & PLZ4 == 8330
replace gdename = "montreux" if gdename == "montreux territet" & PLZ4 == 1820
replace gdename = "montreux" if gdename == "territet" & PLZ4 == 1820
replace gdename = "montreux" if gdename == "territtet" & PLZ4 == 1820
replace gdename = "montreux" if gdename == "terrltet" & PLZ4 == 1820
replace gdename = "schmiedrued" if gdename == "walde" & PLZ4 == 5046
replace gdename = "schmiedrued" if gdename == "walde (ag)" & PLZ4 == 5046
replace gdename = "lavizzara" if gdename == ".fusio mogno" & PLZ4 == 6696
replace gdename = "maggia" if gdename == "auegno" & PLZ4 == 6677
replace gdename = "maggia" if gdename == "aurlgeno" & PLZ4 == 6677
replace gdename = "bosco/gurin" if gdename == "bosco" & PLZ4 == 6685
replace gdename = "bosco/gurin" if gdename == "bosco gurin" & PLZ4 == 6685
replace gdename = "lavizzara" if gdename == "brogllo" & PLZ4 == 6696
replace gdename = "maggia" if gdename == "cogllo" & PLZ4 == 6678
replace gdename = "lavizzara" if gdename == "fusio magna" & PLZ4 == 6696
replace gdename = "lavizzara" if gdename == "fusio mogno" & PLZ4 == 6696
replace gdename = "maggia" if gdename == "giumagliopus" & PLZ4 == 6678
replace gdename = "maggia" if gdename == "giumagllo" & PLZ4 == 6678
replace gdename = "maggia" if gdename == "glnmaglio" & PLZ4 == 6678
replace gdename = "linescio" if gdename == "linesclo" & PLZ4 == 6682
replace gdename = "lavizzara" if gdename == "peccla" & PLZ4 == 6695
replace gdename = "maggia" if gdename == "riveo" & PLZ4 == 6674
replace gdename = "maggia" if gdename == "riveo di someo" & PLZ4 == 6674
replace gdename = "maggia" if gdename == "rlveo" & PLZ4 == 6674
replace gdename = "maggia" if gdename == "ronchina aurigeno" & PLZ4 == 6677
replace gdename = "donneloye" if gdename == "donneioye" & PLZ4 == 1407
replace gdename = "essertines sur yverdon" if gdename == "essertides sur yverdon" & PLZ4 == 1417
replace gdename = "essertines sur yverdon" if gdename == "essertines yverdon" & PLZ4 == 1417
replace gdename = "essertines sur yverdon" if gdename == "essertlnes sur yverdon" & PLZ4 == 1417
replace gdename = "fiez" if gdename == "flez" & PLZ4 == 1420
replace gdename = "grandson" if gdename == "fontaines/grandson" & PLZ4 == 1422
replace gdename = "tevenon" if gdename == "fontanezler" & PLZ4 == 1423
replace gdename = "giez" if gdename == "gier" & PLZ4 == 1429
replace gdename = "giez" if gdename == "giez vd" & PLZ4 == 1429
replace gdename = "giez" if gdename == "glez" & PLZ4 == 1429
replace gdename = "giez" if gdename == "glez vd" & PLZ4 == 1429
replace gdename = "grandevent" if gdename == "grandevent vd" & PLZ4 == 1421
replace gdename = "grandson" if gdename == "les tuileries" & PLZ4 == 1422
replace gdename = "grandson" if gdename == "les tuileries de grandson" & PLZ4 == 1422
replace gdename = "grandson" if gdename == "les tuilerles de grandson" & PLZ4 == 1422
replace gdename = "grandson" if gdename == "les tulleries de grandson" & PLZ4 == 1422
replace gdename = "grandson" if gdename == "les tullerles/grandson" & PLZ4 == 1422
replace gdename = "pailly" if gdename == "pallly" & PLZ4 == 1416
replace gdename = "prahins" if gdename == "prahlns" & PLZ4 == 1408
replace gdename = "grandson" if gdename == "tuileries de grandson" & PLZ4 == 1422
replace gdename = "tevenon" if gdename == "viliars burquin" & PLZ4 == 1423
replace gdename = "tevenon" if gdename == "villars bourquin" & PLZ4 == 1423
replace gdename = "tevenon" if gdename == "villars burquln" & PLZ4 == 1423
replace gdename = "vuarrens" if gdename == "vuarrengei" & PLZ4 == 1418
replace gdename = "vuarrens" if gdename == "vuarrengel" & PLZ4 == 1418
replace gdename = "vugelles la mothe" if gdename == "vugelies la mothe" & PLZ4 == 1431
replace gdename = "fahrni" if gdename == "fahmi bel thun" & PLZ4 == 3617
replace gdename = "fahrni" if gdename == "fahml" & PLZ4 == 3617
replace gdename = "fahrni" if gdename == "fahrei b. thun" & PLZ4 == 3617
replace gdename = "fahrni" if gdename == "fahrni b. thun" & PLZ4 == 3617
replace gdename = "fahrni" if gdename == "fahrni bei thun" & PLZ4 == 3617
replace gdename = "fahrni" if gdename == "fahrni thun" & PLZ4 == 3617
replace gdename = "fahrni" if gdename == "fahrnl bel thun" & PLZ4 == 3617
replace gdename = "forst laengenbuehl" if gdename == "forst b. laengenbuehl" & PLZ4 == 3636
replace gdename = "stocken hoefen" if gdename == "hoefen bei thun" & PLZ4 == 3631
replace gdename = "stocken hoefen" if gdename == "hoefen bel thun" & PLZ4 == 3631
replace gdename = "stocken hoefen" if gdename == "hoefen/thun" & PLZ4 == 3631
replace gdename = "stocken hoefen" if gdename == "nlederstocken" & PLZ4 == 3632
replace gdename = "pohlern" if gdename == "pohiern" & PLZ4 == 3638
replace gdename = "wachseldorn" if gdename == "suederen" & PLZ4 == 3618
replace gdename = "uebeschi" if gdename == "uebeschl" & PLZ4 == 3635
replace gdename = "uebeschi" if gdename == "ueheschi" & PLZ4 == 3635
replace gdename = "wachseldorn" if gdename == "wachseidorn" & PLZ4 == 3618
replace gdename = "aristau" if gdename == "althaeusern" & PLZ4 == 5628
replace gdename = "aristau" if gdename == "arisf au" & PLZ4 == 5628
replace gdename = "aristau" if gdename == "aristau birri" & PLZ4 == 5628
replace gdename = "aristau" if gdename == "arlstau birri" & PLZ4 == 5628
replace gdename = "aristau" if gdename == "artisau" & PLZ4 == 5628
replace gdename = "aristau" if gdename == "artstau" & PLZ4 == 5628
replace gdename = "aristau" if gdename == "barri aristau" & PLZ4 == 5628
replace gdename = "aristau" if gdename == "birri" & PLZ4 == 5628
replace gdename = "aristau" if gdename == "birri aristau" & PLZ4 == 5628
replace gdename = "aristau" if gdename == "birrl" & PLZ4 == 5628
replace gdename = "aristau" if gdename == "birrl artstau" & PLZ4 == 5628
replace gdename = "bremgarten" if gdename == "hermetschwii staffeln" & PLZ4 == 5626
replace gdename = "bremgarten" if gdename == "hermetschwii stallein" & PLZ4 == 5626
replace gdename = "bremgarten" if gdename == "hermetschwil" & PLZ4 == 5626
replace gdename = "bremgarten" if gdename == "hermetschwil staff ein" & PLZ4 == 5626
replace gdename = "bremgarten" if gdename == "hermetschwil staffein" & PLZ4 == 5626
replace gdename = "bremgarten" if gdename == "hermetschwil stallein" & PLZ4 == 5626
replace gdename = "bremgarten" if gdename == "hermetschwll staffeln" & PLZ4 == 5626
replace gdename = "stetten" if gdename == "steffen" & PLZ4 == 5608
replace gdename = "stetten" if gdename == "steffen (ag)" & PLZ4 == 5608
replace gdename = "stetten" if gdename == "steffen ag" & PLZ4 == 5608
replace gdename = "stetten" if gdename == "stehen" & PLZ4 == 5608
replace gdename = "stetten" if gdename == "stetterr ag" & PLZ4 == 5608
replace gdename = "pfaefers" if gdename == "vaettis" & PLZ4 == 7315
replace gdename = "pfaefers" if gdename == "vaettis/pfaefers" & PLZ4 == 7315
replace gdename = "pfaefers" if gdename == "vaettls" & PLZ4 == 7315
replace gdename = "pfaefers" if gdename == "valens" & PLZ4 == 7315
replace gdename = "pfaefers" if gdename == "vattis" & PLZ4 == 7315
replace gdename = "burgistein" if gdename == "burgfisteln" & PLZ4 == 3664
replace gdename = "burgistein" if gdename == "burgfisteln dorf" & PLZ4 == 3664
replace gdename = "burgistein" if gdename == "burgistein dorf" & PLZ4 == 3664
replace gdename = "burgistein" if gdename == "burgistein dort" & PLZ4 == 3664
replace gdename = "burgistein" if gdename == "burgistein station" & PLZ4 == 3664
replace gdename = "burgistein" if gdename == "burgisteln dorf" & PLZ4 == 3664
replace gdename = "burgistein" if gdename == "burglstein" & PLZ4 == 3664
replace gdename = "burgistein" if gdename == "burglsteln station" & PLZ4 == 3664
replace gdename = "willisau" if gdename == "kaeppelimatt" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "ostergau" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "ostergau/wi llisau" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "ostergau/wiliisau land" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "rohrmatt" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "rohrmatt wallisau land" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "rohrmatt willisau land" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "rohrmatt wlllisau land" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "wallisau land" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "wallisau ostergau" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "wiliisau ostergau" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "willisau ostergau" & PLZ4 == 6130
replace gdename = "willisau" if gdename == "wlllisau ostergau" & PLZ4 == 6130
replace gdename = "haegendorf" if gdename == "allerheiligenberg" & PLZ4 == 4615
replace gdename = "haegendorf" if gdename == "allerheiligenberg so" & PLZ4 == 4615
replace gdename = "haegendorf" if gdename == "allerhelligenberg" & PLZ4 == 4615
replace gdename = "haegendorf" if gdename == "allerhelllgenberg" & PLZ4 == 4615
replace gdename = "haegendorf" if gdename == "alterheiligenberg" & PLZ4 == 4615
replace gdename = "hauenstein ifenthal" if gdename == "hauenstein" & PLZ4 == 4633
replace gdename = "hauenstein ifenthal" if gdename == "hauenstein ifenthai" & PLZ4 == 4633
replace gdename = "hauenstein ifenthal" if gdename == "hauenstein lfenthal" & PLZ4 == 4633
replace gdename = "hauenstein ifenthal" if gdename == "hauensteln" & PLZ4 == 4633
replace gdename = "hauenstein ifenthal" if gdename == "hauensteln ifenthal" & PLZ4 == 4633
replace gdename = "hauenstein ifenthal" if gdename == "nauenstein" & PLZ4 == 4633
replace gdename = "hauenstein ifenthal" if gdename == "nauenstein ifenthal" & PLZ4 == 4633
replace gdename = "wisen" if gdename == "wasen" & PLZ4 == 4634
replace gdename = "egnach" if gdename == "winden" & PLZ4 == 9315
replace gdename = "egnach" if gdename == "winden egnach" & PLZ4 == 9315
replace gdename = "egnach" if gdename == "winden tg" & PLZ4 == 9315
replace gdename = "egnach" if gdename == "winden/egnach" & PLZ4 == 9315
replace gdename = "oberdiessbach" if gdename == "aeachlen b. oberdlessbach" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "aeschien" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "aeschien b. oberdiessbach" & PLZ4 == 3672
replace gdename = "sigriswil" if gdename == "aeschien bei oberdiessbach/sigriswil" & PLZ4 == 3655
replace gdename = "oberdiessbach" if gdename == "aeschlen b. oberdiessbach" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "aeschlen bei oberdiesbach" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "aeschlen bei oberdiessbach" & PLZ4 == 3672
replace gdename = "sigriswil" if gdename == "aeschlen bei oberdiessbach/ sigriswil" & PLZ4 == 3655
replace gdename = "sigriswil" if gdename == "aeschlen bei oberdiessbach/sigriswil" & PLZ4 == 3655
replace gdename = "oberdiessbach" if gdename == "aeschlen bel oberdiessbach" & PLZ4 == 3672
replace gdename = "sigriswil" if gdename == "aeschlen bel oberdiessbach/ sigriswfl" & PLZ4 == 3655
replace gdename = "oberdiessbach" if gdename == "aeschlen kuekenvertriebs ag" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "bielken bei oberdiessbach" & PLZ4 == 3674
replace gdename = "oberdiessbach" if gdename == "bleiken" & PLZ4 == 3674
replace gdename = "oberdiessbach" if gdename == "bleiken bel oberdiessbach" & PLZ4 == 3674
replace gdename = "linden" if gdename == "jassbach" & PLZ4 == 3673
replace gdename = "linden" if gdename == "jassbach zur linden" & PLZ4 == 3673
replace gdename = "buchholterberg" if gdename == "wangelen" & PLZ4 == 3615
replace gdename = "les bois" if gdename == "boechet" & PLZ4 == 2336
replace gdename = "la chaux de fonds" if gdename == "chaux de fonds" & PLZ4 == 2300
replace gdename = "muriaux" if gdename == "emibois be" & PLZ4 == 2338
replace gdename = "la ferriere" if gdename == "la chaux d abel" & PLZ4 == 2333
replace gdename = "la sagne" if gdename == "la corbatiere" & PLZ4 == 2314
replace gdename = "la sagne" if gdename == "la torbatiere" & PLZ4 == 2314
replace gdename = "les bois" if gdename == "le boechet" & PLZ4 == 2336
replace gdename = "muriaux" if gdename == "les emibois" & PLZ4 == 2338
replace gdename = "muriaux" if gdename == "les emibois muriaux" & PLZ4 == 2338
replace gdename = "muriaux" if gdename == "les emibols" & PLZ4 == 2338
replace gdename = "muriaux" if gdename == "les emibols murlaux" & PLZ4 == 2338
replace gdename = "muriaux" if gdename == "les emlbois" & PLZ4 == 2338
replace gdename = "muriaux" if gdename == "les emlbois muraux" & PLZ4 == 2338
replace gdename = "muriaux" if gdename == "les emlbois muriaux" & PLZ4 == 2338
replace gdename = "muriaux" if gdename == "les emlbols" & PLZ4 == 2338
replace gdename = "muriaux" if gdename == "les emã¬bols" & PLZ4 == 2338
replace gdename = "muriaux" if gdename == "murlaux" & PLZ4 == 2338
replace gdename = "muriaux" if gdename == "murtaux" & PLZ4 == 2338
replace gdename = "rueeggisberg" if gdename == "heigisried rohrbach" & PLZ4 == 3155
replace gdename = "rueeggisberg" if gdename == "heigisried rueeggisberg" & PLZ4 == 3088
replace gdename = "rueeggisberg" if gdename == "helgisried rohrbach" & PLZ4 == 3155
replace gdename = "rueeggisberg" if gdename == "helgisried rueeggisberg" & PLZ4 == 3088
replace gdename = "rueeggisberg" if gdename == "helgisried rueegglsberg" & PLZ4 == 3088
replace gdename = "guggisberg" if gdename == "hirschmatt" & PLZ4 == 3158
replace gdename = "guggisberg" if gdename == "riffenmatt" & PLZ4 == 3158
replace gdename = "rueschegg" if gdename == "rueschegg graben" & PLZ4 == 3154
replace gdename = "bursins" if gdename == "bursitis" & PLZ4 == 1183
replace gdename = "bursins" if gdename == "burslns" & PLZ4 == 1183
replace gdename = "essertines sur rolle" if gdename == "essertine sur rolle" & PLZ4 == 1186
replace gdename = "essertines sur rolle" if gdename == "essertines/rolle" & PLZ4 == 1186
replace gdename = "essertines sur rolle" if gdename == "essertlne sur rolle" & PLZ4 == 1186
replace gdename = "essertines sur rolle" if gdename == "essertlnes sur rolle" & PLZ4 == 1186
replace gdename = "gilly" if gdename == "gaily" & PLZ4 == 1182
replace gdename = "gilly" if gdename == "giily" & PLZ4 == 1182
replace gdename = "gilly" if gdename == "giity" & PLZ4 == 1182
replace gdename = "gilly" if gdename == "glily" & PLZ4 == 1182
replace gdename = "gilly" if gdename == "gllly" & PLZ4 == 1182
replace gdename = "luins" if gdename == "lulns" & PLZ4 == 1184
replace gdename = "mont sur rolle" if gdename == "mont rolle" & PLZ4 == 1185
replace gdename = "mont sur rolle" if gdename == "mont s rolle" & PLZ4 == 1185
replace gdename = "mont sur rolle" if gdename == "mont s. rolle" & PLZ4 == 1185
replace gdename = "mont sur rolle" if gdename == "mont s/rolie" & PLZ4 == 1185
replace gdename = "mont sur rolle" if gdename == "mont s/rolle" & PLZ4 == 1185
replace gdename = "mont sur rolle" if gdename == "mont sun rolle" & PLZ4 == 1185
replace gdename = "mont sur rolle" if gdename == "mont sur roller" & PLZ4 == 1185
replace gdename = "mont sur rolle" if gdename == "mont sur roule" & PLZ4 == 1185
replace gdename = "mont sur rolle" if gdename == "monts rolle" & PLZ4 == 1185
replace gdename = "saint oyens" if gdename == "st oyens" & PLZ4 == 1187
replace gdename = "tartegnin" if gdename == "tartegniri" & PLZ4 == 1180
replace gdename = "vinzel" if gdename == "vinzei" & PLZ4 == 1184
replace gdename = "vinzel" if gdename == "vinzel vd" & PLZ4 == 1184
replace gdename = "vinzel" if gdename == "vlnzel" & PLZ4 == 1184
replace gdename = "luins" if gdename == "wins" & PLZ4 == 1184
replace gdename = "saignelegier" if gdename == "cemievillers" & PLZ4 == 2352
replace gdename = "saignelegier" if gdename == "cernieviliers" & PLZ4 == 2352
replace gdename = "saignelegier" if gdename == "cernievillers" & PLZ4 == 2352
replace gdename = "saignelegier" if gdename == "les" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "les ceriatez" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "les cerlatez/saignelegier" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "les cerlatez/salgnelegier" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "saigneiegier" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "saignelegier ju" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "saignelegiers" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "saignelegiier" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "saignelegler" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "saignelegter" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "saignetegier" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "saine legier" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "sainelegier" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "salgnelegier" & PLZ4 == 2350
replace gdename = "saignelegier" if gdename == "salgnelegler" & PLZ4 == 2350
replace gdename = "kaiserstuhl" if gdename == "kaiserstuhi" & PLZ4 == 5466
replace gdename = "kaiserstuhl" if gdename == "kaiserstuhl ag" & PLZ4 == 5466
replace gdename = "kaiserstuhl" if gdename == "kalserstuhl" & PLZ4 == 5466
replace gdename = "kaiserstuhl" if gdename == "kalserstuhl (ag)" & PLZ4 == 5466
replace gdename = "kaiserstuhl" if gdename == "kalserstuhl ag" & PLZ4 == 5466
replace gdename = "haegendorf" if gdename == "hagendorn" & PLZ4 == 4614
replace gdename = "weggis" if gdename == "hertenetein" & PLZ4 == 6353
replace gdename = "weggis" if gdename == "hertenstein" & PLZ4 == 6353
replace gdename = "weggis" if gdename == "hertenstein weggis" & PLZ4 == 6353
replace gdename = "weggis" if gdename == "hertensteln" & PLZ4 == 6353
replace gdename = "weggis" if gdename == "hertensteln weggis" & PLZ4 == 6353
replace gdename = "les breuleux" if gdename == "breuieux" & PLZ4 == 2345
replace gdename = "les breuleux" if gdename == "breuleux" & PLZ4 == 2345
replace gdename = "muriaux" if gdename == "cerneux veusil" & PLZ4 == 2345
replace gdename = "muriaux" if gdename == "cerneux veusll" & PLZ4 == 2345
replace gdename = "les breuleux" if gdename == "le breuleux" & PLZ4 == 2345
replace gdename = "muriaux" if gdename == "le cemeux veusil" & PLZ4 == 2345
replace gdename = "muriaux" if gdename == "le cerneux veusil" & PLZ4 == 2345
replace gdename = "muriaux" if gdename == "le roselet" & PLZ4 == 2345
replace gdename = "les breuleux" if gdename == "lea breuleux" & PLZ4 == 2345
replace gdename = "les breuleux" if gdename == "les breuieux" & PLZ4 == 2345
replace gdename = "horw" if gdename == "niklausen" & PLZ4 == 6048
replace gdename = "horw" if gdename == "st. niklausen" & PLZ4 == 6048
replace gdename = "horw" if gdename == "st. niklausen (lu)" & PLZ4 == 6048
replace gdename = "horw" if gdename == "st. niklausen horw" & PLZ4 == 6048
replace gdename = "horw" if gdename == "st. niklausen/lu" & PLZ4 == 6048
replace gdename = "horw" if gdename == "st. nlklausen" & PLZ4 == 6048
replace gdename = "horw" if gdename == "st.niklausen" & PLZ4 == 6048
replace gdename = "meikirch" if gdename == "wahlendorf" & PLZ4 == 3046
replace gdename = "balm bei guensberg" if gdename == "balm bel guensberg" & PLZ4 == 4525
replace gdename = "balm bei guensberg" if gdename == "balm guensberg" & PLZ4 == 4525
replace gdename = "balm bei guensberg" if gdename == "balmberg" & PLZ4 == 4525
replace gdename = "farnern" if gdename == "famem" & PLZ4 == 4539
replace gdename = "farnern" if gdename == "famern" & PLZ4 == 4539
replace gdename = "farnern" if gdename == "farnem" & PLZ4 == 4539
replace gdename = "farnern" if gdename == "fernern" & PLZ4 == 4539
replace gdename = "drei hoefe" if gdename == "hersiwll" & PLZ4 == 4558
replace gdename = "drei hoefe" if gdename == "herslwil" & PLZ4 == 4558
replace gdename = "horriwil" if gdename == "horriwii" & PLZ4 == 4557
replace gdename = "horriwil" if gdename == "horriwil so" & PLZ4 == 4557
replace gdename = "horriwil" if gdename == "horriwll" & PLZ4 == 4557
replace gdename = "horriwil" if gdename == "horriwll so" & PLZ4 == 4557
replace gdename = "horriwil" if gdename == "horrlwil" & PLZ4 == 4557
replace gdename = "horriwil" if gdename == "horrlwll" & PLZ4 == 4557
replace gdename = "horriwil" if gdename == "horrwil" & PLZ4 == 4557
replace gdename = "hubersdorf" if gdename == "huberadorf" & PLZ4 == 4535
replace gdename = "hubersdorf" if gdename == "hubersdorf so" & PLZ4 == 4535
replace gdename = "kammersrohr" if gdename == "kammersohr" & PLZ4 == 4535
replace gdename = "kammersrohr" if gdename == "kammersrohr hubersdorf" & PLZ4 == 4535
replace gdename = "riedholz" if gdename == "niederwal" & PLZ4 == 4523
replace gdename = "riedholz" if gdename == "niederwal so" & PLZ4 == 4523
replace gdename = "rumisberg" if gdename == "rumasberg" & PLZ4 == 4539
replace gdename = "rumisberg" if gdename == "rumlaberg" & PLZ4 == 4539
replace gdename = "rumisberg" if gdename == "rumlsberg" & PLZ4 == 4539
replace gdename = "drei hoefe" if gdename == "winistorf" & PLZ4 == 4558
replace gdename = "drei hoefe" if gdename == "wlnistorf" & PLZ4 == 4558
replace gdename = "aeschi (so)" if gdename == "aeschl" & PLZ4 == 4556
replace gdename = "aeschi (so)" if gdename == "aeschl so" & PLZ4 == 4556
replace gdename = "bolken" if gdename == "bocken" & PLZ4 == 4556
replace gdename = "bolken" if gdename == "boiken" & PLZ4 == 4556
replace gdename = "bolken" if gdename == "bolken ebersol" & PLZ4 == 4556
replace gdename = "aeschi (so)" if gdename == "burgaeschi" & PLZ4 == 4556
replace gdename = "aeschi (so)" if gdename == "burgaeschl" & PLZ4 == 4556
replace gdename = "hellsau" if gdename == "heilsau" & PLZ4 == 3429
replace gdename = "hellsau" if gdename == "helisau" & PLZ4 == 3429
replace gdename = "heimenhausen" if gdename == "helmenhausen" & PLZ4 == 3373
replace gdename = "hoechstetten" if gdename == "hoerstetten" & PLZ4 == 3429
replace gdename = "heimenhausen" if gdename == "roethenbach" & PLZ4 == 3373
replace gdename = "heimenhausen" if gdename == "roethenbach b.herzogenbuch" & PLZ4 == 3373
replace gdename = "heimenhausen" if gdename == "roethenbach bel herzogenbuchsee" & PLZ4 == 3373
replace gdename = "heimenhausen" if gdename == "roethenbach bet herzogenbuchsee" & PLZ4 == 3373
replace gdename = "heimenhausen" if gdename == "roethenbach i.e" & PLZ4 == 3373
replace gdename = "heimenhausen" if gdename == "rothenbach bei herzogenbuchsee" & PLZ4 == 3373
replace gdename = "wangenried" if gdename == "wangenrled" & PLZ4 == 3374
replace gdename = "wangenried" if gdename == "wangeriried" & PLZ4 == 3374
replace gdename = "heimenhausen" if gdename == "wanzenwil" & PLZ4 == 3372
replace gdename = "heimenhausen" if gdename == "wanzenwll" & PLZ4 == 3372
replace gdename = "saviese" if gdename == "st germain" & PLZ4 == 1965
replace gdename = "chamoson" if gdename == "st pierce de clages" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st pierre de cl ages" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st pierre de clages" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st pierre de clages/chamoson" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st. pierre de ciages" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st. pierre de clages" & PLZ4 == 1955
replace gdename = "reiden" if gdename == "reidermoos" & PLZ4 == 6260
replace gdename = "reiden" if gdename == "reidermoos/reiden" & PLZ4 == 6260
replace gdename = "reiden" if gdename == "reldermoos" & PLZ4 == 6260
replace gdename = "wiliberg" if gdename == "wiliberg hintermoos" & PLZ4 == 5058
replace gdename = "wiliberg" if gdename == "wlliberg" & PLZ4 == 5058
replace gdename = "wiliberg" if gdename == "wlliberg hintermoos" & PLZ4 == 5058
replace gdename = "basadingen schlattingen" if gdename == "basedingen" & PLZ4 == 8254
replace gdename = "basadingen schlattingen" if gdename == "bassdingen" & PLZ4 == 8254
replace gdename = "schlatt" if gdename == "neu paradies unterschlatt" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "neu paradies/unterschlaff" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "neu paradies/unterschlatt" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "paradies/schlatt tg" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schiatt bei diessenhofen" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schlaff b. diessenhofen" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schlaff bei diessenhofen" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schlaft" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schlaft b. diessenhofen" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schlatt b. diessenhofen" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schlatt b. dlessenhofen" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schlatt bei diessenhofen" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schlatt diesenhofen" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schlatt dlesenhofen" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "schlatt tg" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "unterschiatt" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "unterschlaff tg" & PLZ4 == 8252
replace gdename = "schlatt" if gdename == "unterschlatt tg" & PLZ4 == 8252
replace gdename = "hemberg" if gdename == "bachli (hemberg)" & PLZ4 == 9633
replace gdename = "hemberg" if gdename == "baechii" & PLZ4 == 9633
replace gdename = "hemberg" if gdename == "baechli" & PLZ4 == 9633
replace gdename = "hemberg" if gdename == "baechli (hemberg)" & PLZ4 == 9633
replace gdename = "hemberg" if gdename == "baechll" & PLZ4 == 9633
replace gdename = "hemberg" if gdename == "baechll (hemberg)" & PLZ4 == 9633
replace gdename = "neckertal" if gdename == "necker" & PLZ4 == 9126
replace gdename = "chatel saint denis" if gdename == "les paccots" & PLZ4 == 1619
replace gdename = "lohn ammannsegg" if gdename == "ammannseg" & PLZ4 == 4573
replace gdename = "lohn ammannsegg" if gdename == "ammannsegg" & PLZ4 == 4573
replace gdename = "lohn ammannsegg" if gdename == "ammannsegg/kriegstetten" & PLZ4 == 4573
replace gdename = "lohn ammannsegg" if gdename == "ammansegg" & PLZ4 == 4573
replace gdename = "lohn ammannsegg" if gdename == "ammansegg (so)" & PLZ4 == 4573
replace gdename = "lohn ammannsegg" if gdename == "arrmannsegg" & PLZ4 == 4573
replace gdename = "mels" if gdename == "schwendi (weiestannental)" & PLZ4 == 7325
replace gdename = "mels" if gdename == "schwendi (weisstannental)" & PLZ4 == 7325
replace gdename = "mels" if gdename == "schwendl (welsstannental)" & PLZ4 == 7325
replace gdename = "mels" if gdename == "weisstannen" & PLZ4 == 7326
replace gdename = "blenio" if gdename == "aquilia" & PLZ4 == 6719
replace gdename = "blenio" if gdename == "campa/blenio" & PLZ4 == 6720
replace gdename = "blenio" if gdename == "campo bienio" & PLZ4 == 6720
replace gdename = "blenio" if gdename == "campo blenlo" & PLZ4 == 6720
replace gdename = "blenio" if gdename == "campo/blenlo" & PLZ4 == 6720
replace gdename = "acquarossa" if gdename == "corzoneso piano" & PLZ4 == 6722
replace gdename = "acquarossa" if gdename == "corzoneso plano" & PLZ4 == 6722
replace gdename = "acquarossa" if gdename == "corzoresco" & PLZ4 == 6722
replace gdename = "acquarossa" if gdename == "largarlo" & PLZ4 == 6724
replace gdename = "acquarossa" if gdename == "loittigna" & PLZ4 == 6716
replace gdename = "acquarossa" if gdename == "lottlgna" & PLZ4 == 6716
replace gdename = "serravalle" if gdename == "ludlam)" & PLZ4 == 6721
replace gdename = "serravalle" if gdename == "ludlano" & PLZ4 == 6721
replace gdename = "acquarossa" if gdename == "maroita" & PLZ4 == 6723
replace gdename = "acquarossa" if gdename == "marotta" & PLZ4 == 6723
replace gdename = "blenio" if gdename == "motto" & PLZ4 == 6720
replace gdename = "blenio" if gdename == "motto (bienio)" & PLZ4 == 6720
replace gdename = "blenio" if gdename == "motto (blenio)" & PLZ4 == 6720
replace gdename = "blenio" if gdename == "motto (blenio)/dongio" & PLZ4 == 6720
replace gdename = "acquarossa" if gdename == "ponte valentino" & PLZ4 == 6724
replace gdename = "acquarossa" if gdename == "pruglasco" & PLZ4 == 6723
replace gdename = "baldingen" if gdename == "baidingen" & PLZ4 == 5333
replace gdename = "baldingen" if gdename == "baldiingen" & PLZ4 == 5333
replace gdename = "baldingen" if gdename == "betdingen" & PLZ4 == 5333
replace gdename = "boebikon" if gdename == "boeblkon" & PLZ4 == 5334
replace gdename = "wislikofen" if gdename == "mellstorf/wislikofen" & PLZ4 == 5463
replace gdename = "wislikofen" if gdename == "mellstorf/wisllkofen" & PLZ4 == 5463
replace gdename = "ruemikon" if gdename == "ruemikon ag" & PLZ4 == 5464
replace gdename = "ruemikon" if gdename == "ruemlkon" & PLZ4 == 5464
replace gdename = "ruemikon" if gdename == "ruemlkon ag" & PLZ4 == 5464
replace gdename = "ruemikon" if gdename == "rumikon" & PLZ4 == 5464
replace gdename = "ruemikon" if gdename == "rumikon ag" & PLZ4 == 5464
replace gdename = "siglistorf" if gdename == "sigiistorf" & PLZ4 == 5462
replace gdename = "siglistorf" if gdename == "sigilstorf" & PLZ4 == 5462
replace gdename = "siglistorf" if gdename == "siglisdorf" & PLZ4 == 5462
replace gdename = "siglistorf" if gdename == "siglistori" & PLZ4 == 5462
replace gdename = "siglistorf" if gdename == "siglistort" & PLZ4 == 5462
replace gdename = "siglistorf" if gdename == "sigllstorf" & PLZ4 == 5462
replace gdename = "siglistorf" if gdename == "sigã¬istorf" & PLZ4 == 5462
replace gdename = "siglistorf" if gdename == "slglisdorf" & PLZ4 == 5462
replace gdename = "siglistorf" if gdename == "slglistorf" & PLZ4 == 5462
replace gdename = "siglistorf" if gdename == "sïgliistorf" & PLZ4 == 5462
replace gdename = "wislikofen" if gdename == "wislikof en" & PLZ4 == 5463
replace gdename = "wislikofen" if gdename == "wisllkofen" & PLZ4 == 5463
replace gdename = "wislikofen" if gdename == "wlsllkofen" & PLZ4 == 5463
replace gdename = "chamoson" if gdename == "st pierre d.c" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st pierre de clages" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st plerre de clages" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st. pier de clages" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st. pierre" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st. pierre de clages" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st. pierrre de clages" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st. plerre de clages" & PLZ4 == 1955
replace gdename = "chamoson" if gdename == "st. plerrre de clages" & PLZ4 == 1955
replace gdename = "leuggern" if gdename == "felsenau" & PLZ4 == 5316
replace gdename = "leuggern" if gdename == "felsenau (ag)" & PLZ4 == 5316
replace gdename = "full reuenthal" if gdename == "fuli reuenthal" & PLZ4 == 5324
replace gdename = "full reuenthal" if gdename == "full" & PLZ4 == 5324
replace gdename = "full reuenthal" if gdename == "full reuenthai" & PLZ4 == 5324
replace gdename = "full reuenthal" if gdename == "reuenthal" & PLZ4 == 5324
replace gdename = "muhen" if gdename == "ober muhen" & PLZ4 == 5037
replace gdename = "muhen" if gdename == "obermuhen" & PLZ4 == 5037
replace gdename = "chateau d oex" if gdename == "les moulins" & PLZ4 == 1660
replace gdename = "chateau d oex" if gdename == "les moulins/ chateau d oex" & PLZ4 == 1660
replace gdename = "chateau d oex" if gdename == "les moulins/chateau d oex" & PLZ4 == 1660
replace gdename = "affeltrangen" if gdename == "buch bei maerwil" & PLZ4 == 9562
replace gdename = "bussnang" if gdename == "oppikon" & PLZ4 == 9565
replace gdename = "bussnang" if gdename == "opplkon" & PLZ4 == 9565
replace gdename = "bussnang" if gdename == "schmidshof" & PLZ4 == 9565
replace gdename = "bussnang" if gdename == "schmidshof buch" & PLZ4 == 9565
replace gdename = "bussnang" if gdename == "schmidshof/bussnang" & PLZ4 == 9565
replace gdename = "bussnang" if gdename == "schmldshof" & PLZ4 == 9565
replace gdename = "oron" if gdename == "paiezieux village" & PLZ4 == 1607
replace gdename = "oron" if gdename == "palezieux village" & PLZ4 == 1607
replace gdename = "oron" if gdename == "palezleux village" & PLZ4 == 1607
replace gdename = "luzein" if gdename == "buchen" & PLZ4 == 7223
replace gdename = "luzein" if gdename == "buchen i. p" & PLZ4 == 7223
replace gdename = "luzein" if gdename == "buchen i.p" & PLZ4 == 7223
replace gdename = "luzein" if gdename == "buchen staad" & PLZ4 == 7223
replace gdename = "luzein" if gdename == "buchen stead" & PLZ4 == 7223
replace gdename = "schiers" if gdename == "schuders" & PLZ4 == 7220
replace gdename = "scuol" if gdename == "schuls" & PLZ4 == 7550
replace gdename = "scuol" if gdename == "scuol" & PLZ4 == 7550
replace gdename = "scuol" if gdename == "scuols/schuls" & PLZ4 == 7550
replace gdename = "wagenhausen" if gdename == "etzwflen" & PLZ4 == 8259
replace gdename = "wagenhausen" if gdename == "etzwilen" & PLZ4 == 8259
replace gdename = "wagenhausen" if gdename == "etzwilen tg" & PLZ4 == 8259
replace gdename = "wagenhausen" if gdename == "etzwilen/kaltenbach" & PLZ4 == 8259
replace gdename = "wagenhausen" if gdename == "etzwlien" & PLZ4 == 8259
replace gdename = "wagenhausen" if gdename == "etzwllen" & PLZ4 == 8259
replace gdename = "wagenhausen" if gdename == "etzwllen/kaltenbach" & PLZ4 == 8259
replace gdename = "wagenhausen" if gdename == "rheinklingen" & PLZ4 == 8259
replace gdename = "warth weiningen" if gdename == "weiningen tg" & PLZ4 == 8532
replace gdename = "warth weiningen" if gdename == "welningen" & PLZ4 == 8532
replace gdename = "warth weiningen" if gdename == "welningen tg" & PLZ4 == 8532
replace gdename = "warth weiningen" if gdename == "wetningen" & PLZ4 == 8532
replace gdename = "warth weiningen" if gdename == "wetningen tg" & PLZ4 == 8532
replace gdename = "eiken" if gdename == "elken" & PLZ4 == 5074
replace gdename = "bex" if gdename == "frenieres" & PLZ4 == 1880
replace gdename = "bex" if gdename == "les devens" & PLZ4 == 1880
replace gdename = "bex" if gdename == "les pians sur bex" & PLZ4 == 1880
replace gdename = "bex" if gdename == "les plans" & PLZ4 == 1880
replace gdename = "bex" if gdename == "les plans sur bex" & PLZ4 == 1880
replace gdename = "bex" if gdename == "pians sur bex" & PLZ4 == 1880
replace gdename = "bex" if gdename == "plans sur bex" & PLZ4 == 1880
replace gdename = "bex" if gdename == "posses bex" & PLZ4 == 1880
replace gdename = "rochefort" if gdename == "brot dessus" & PLZ4 == 2149
replace gdename = "sumiswald" if gdename == "griesbach sumiswald" & PLZ4 == 3454
replace gdename = "affoltern im emmental" if gdename == "haeusernmoos" & PLZ4 == 3416
replace gdename = "affoltern im emmental" if gdename == "haeusernmoos i. e." & PLZ4 == 3416
replace gdename = "walterswil" if gdename == "schmidigen" & PLZ4 == 3464
replace gdename = "walterswil" if gdename == "schmidigen muehleweg" & PLZ4 == 3464
replace gdename = "walterswil" if gdename == "schmidigen muehlweg" & PLZ4 == 3464
replace gdename = "walterswil" if gdename == "schmidigen muhleweg" & PLZ4 == 3464
replace gdename = "walterswil" if gdename == "schmidlgen muehleweg" & PLZ4 == 3464
replace gdename = "walterswil" if gdename == "schmldigen muehleweg" & PLZ4 == 3464
replace gdename = "affoltern im emmental" if gdename == "weber i e" & PLZ4 == 3462
replace gdename = "affoltern im emmental" if gdename == "weier be" & PLZ4 == 3462
replace gdename = "affoltern im emmental" if gdename == "weier i.e" & PLZ4 == 3462
replace gdename = "affoltern im emmental" if gdename == "weier im emmental" & PLZ4 == 3462
replace gdename = "affoltern im emmental" if gdename == "weler im emmental" & PLZ4 == 3462
replace gdename = "aeugst am albis" if gdename == "aeugstertal" & PLZ4 == 8914
replace gdename = "aeugst am albis" if gdename == "aeugstertal/aeugst" & PLZ4 == 8914
replace gdename = "aeugst am albis" if gdename == "aeugstertal/aeugst am albis" & PLZ4 == 8914
replace gdename = "ollon" if gdename == "arveyes" & PLZ4 == 1884
replace gdename = "ollon" if gdename == "arveyes ollon" & PLZ4 == 1884
replace gdename = "ollon" if gdename == "arveyes villars" & PLZ4 == 1884
replace gdename = "ollon" if gdename == "arveyes/ollon" & PLZ4 == 1884
replace gdename = "courchapoix" if gdename == "ceurchapolx" & PLZ4 == 2825
replace gdename = "chatillon" if gdename == "chaetillon" & PLZ4 == 2843
replace gdename = "chatillon" if gdename == "chatilion" & PLZ4 == 2843
replace gdename = "courroux" if gdename == "courcellon" & PLZ4 == 2823
replace gdename = "courroux" if gdename == "courcelon" & PLZ4 == 2823
replace gdename = "courchapoix" if gdename == "courchapolx" & PLZ4 == 2825
replace gdename = "develier" if gdename == "devetier" & PLZ4 == 2802
replace gdename = "la scheulte" if gdename == "la scheulten" & PLZ4 == 2827
replace gdename = "mervelier" if gdename == "merveller" & PLZ4 == 2827
replace gdename = "val terbi" if gdename == "montseveller" & PLZ4 == 2828
replace gdename = "movelier" if gdename == "moveller" & PLZ4 == 2812
replace gdename = "rossemaison" if gdename == "rossemalson" & PLZ4 == 2842
replace gdename = "rossemaison" if gdename == "rossmaisen" & PLZ4 == 2842
replace gdename = "rossemaison" if gdename == "rossmaison" & PLZ4 == 2842
replace gdename = "val terbi" if gdename == "vernies" & PLZ4 == 2829
replace gdename = "ennetbuergen" if gdename == "buergenstock" & PLZ4 == 6373
replace gdename = "ennetbuergen" if gdename == "buergenstock/ ennetbuergen" & PLZ4 == 6373
replace gdename = "ennetbuergen" if gdename == "buergenstock/ennetbuergen" & PLZ4 == 6373
replace gdename = "ennetbuergen" if gdename == "buergenstock/stansstad" & PLZ4 == 6373
replace gdename = "ennetbuergen" if gdename == "burgenstock" & PLZ4 == 6373
replace gdename = "winterthur" if gdename == "reutlingen" & PLZ4 == 8404
replace gdename = "winterthur" if gdename == "reutlingen (winterthur)" & PLZ4 == 8404
replace gdename = "winterthur" if gdename == "reutlingen winterthur" & PLZ4 == 8404
replace gdename = "st niklaus" if gdename == "embal" & PLZ4 == 3927
replace gdename = "st niklaus" if gdename == "herbriggen" & PLZ4 == 3927
replace gdename = "st niklaus" if gdename == "herbrlggen" & PLZ4 == 3927
replace gdename = "randa" if gdename == "rande" & PLZ4 == 3928
replace gdename = "randa" if gdename == "renda" & PLZ4 == 3928
replace gdename = "st niklaus" if gdename == "st. nlklaus" & PLZ4 == 3924
replace gdename = "staldenried" if gdename == "staddenried" & PLZ4 == 3933
replace gdename = "staldenried" if gdename == "staidenried" & PLZ4 == 3933
replace gdename = "staldenried" if gdename == "staldenrled" & PLZ4 == 3933
replace gdename = "taesch" if gdename == "taesch vs" & PLZ4 == 3929
replace gdename = "taesch" if gdename == "tisch" & PLZ4 == 3929
replace gdename = "deitingen" if gdename == "deftingen" & PLZ4 == 4543
replace gdename = "deitingen" if gdename == "deiingen" & PLZ4 == 4543
replace gdename = "deitingen" if gdename == "delitingen" & PLZ4 == 4543
replace gdename = "deitingen" if gdename == "dellingen" & PLZ4 == 4543
replace gdename = "deitingen" if gdename == "deltingen" & PLZ4 == 4543
replace gdename = "deitingen" if gdename == "dettingen" & PLZ4 == 4543
replace gdename = "buchillon" if gdename == "buchiilon" & PLZ4 == 1164
replace gdename = "morat" if gdename == "buechsien" & PLZ4 == 3215
replace gdename = "kleinboesingen" if gdename == "kieinboesingen" & PLZ4 == 3213
replace gdename = "kleinboesingen" if gdename == "klelnboesingen" & PLZ4 == 3213
replace gdename = "gurmels" if gdename == "liebislorf" & PLZ4 == 3213
replace gdename = "gurmels" if gdename == "liebistorf gruenenburg" & PLZ4 == 3213
replace gdename = "gurmels" if gdename == "lieblstorf gruenenburg" & PLZ4 == 3213
replace gdename = "gurmels" if gdename == "llebistorf" & PLZ4 == 3213
replace gdename = "gurmels" if gdename == "llebistorf gruenenburg" & PLZ4 == 3213
replace gdename = "morat" if gdename == "lutrigen" & PLZ4 == 3215
replace gdename = "morat" if gdename == "lutrlgen" & PLZ4 == 3215
replace gdename = "morat" if gdename == "luttigen3452307" & PLZ4 == 3215
replace gdename = "schmitten" if gdename == "ried bei berg (schmitten)" & PLZ4 == 3185
replace gdename = "ried bei kerzers" if gdename == "ried bel karzers" & PLZ4 == 3216
replace gdename = "ried bei kerzers" if gdename == "ried bel kerzers" & PLZ4 == 3216
replace gdename = "ried bei kerzers" if gdename == "ried/kerzers" & PLZ4 == 3216
replace gdename = "ried bei kerzers" if gdename == "rled bel kerzers" & PLZ4 == 3216
replace gdename = "ueberstorf" if gdename == "uebistorf gruenenburg" & PLZ4 == 3182
replace gdename = "ulmiz" if gdename == "uimiz" & PLZ4 == 3214
replace gdename = "ulmiz" if gdename == "ulmiz sg" & PLZ4 == 3214
replace gdename = "ulmiz" if gdename == "ulmlz" & PLZ4 == 3214
replace gdename = "davos" if gdename == "davos monstein" & PLZ4 == 7278
replace gdename = "davos" if gdename == "frauenkirch" & PLZ4 == 7276
replace gdename = "davos" if gdename == "frauenkirch/davos" & PLZ4 == 7276
replace gdename = "davos" if gdename == "glaris davos" & PLZ4 == 7277
replace gdename = "davos" if gdename == "monstein" & PLZ4 == 7278
replace gdename = "autigny" if gdename == "autlgny" & PLZ4 == 1742
replace gdename = "cottens (fr)" if gdename == "cottans" & PLZ4 == 1741
replace gdename = "cottens (fr)" if gdename == "cottene fr" & PLZ4 == 1741
replace gdename = "cottens (fr)" if gdename == "cottons" & PLZ4 == 1741
replace gdename = "torny" if gdename == "middes fr" & PLZ4 == 1749
replace gdename = "neyruz" if gdename == "neyrãºz" & PLZ4 == 1740
replace gdename = "prez vers noreaz" if gdename == "prez" & PLZ4 == 1746
replace gdename = "prez vers noreaz" if gdename == "prez vers noreez" & PLZ4 == 1746
replace gdename = "torny" if gdename == "tomy le grand" & PLZ4 == 1748
replace gdename = "torny" if gdename == "torny ie grand" & PLZ4 == 1748
replace gdename = "gibloux" if gdename == "viilareel le gibloux" & PLZ4 == 1695
replace gdename = "gibloux" if gdename == "viilarsel le giboux" & PLZ4 == 1695
replace gdename = "la folliaz" if gdename == "viliarimboud" & PLZ4 == 1690
replace gdename = "gibloux" if gdename == "viliarsel le gibloux" & PLZ4 == 1695
replace gdename = "la folliaz" if gdename == "villarlmboud" & PLZ4 == 1690
replace gdename = "villars sur glane" if gdename == "villars sur glaene" & PLZ4 == 1752
replace gdename = "villars sur glane" if gdename == "villars sur glene" & PLZ4 == 1752
replace gdename = "gibloux" if gdename == "villarsel" & PLZ4 == 1695
replace gdename = "gibloux" if gdename == "villarsel le gibioux" & PLZ4 == 1695
replace gdename = "gibloux" if gdename == "villarsel le giboux" & PLZ4 == 1695
replace gdename = "gibloux" if gdename == "villarsel le glboux" & PLZ4 == 1695
replace gdename = "villaz saint pierre" if gdename == "villaz st pierre" & PLZ4 == 1690
replace gdename = "la folliaz" if gdename == "vlllarimboud" & PLZ4 == 1690
replace gdename = "la folliaz" if gdename == "vlllarlmboud" & PLZ4 == 1690
replace gdename = "gibloux" if gdename == "vlllarsel le gibloux" & PLZ4 == 1695
replace gdename = "grossaffoltern" if gdename == "suberg" & PLZ4 == 3262
replace gdename = "grossaffoltern" if gdename == "surberg" & PLZ4 == 3262
replace gdename = "fraubrunnen" if gdename == "etzeikofen" & PLZ4 == 3306
replace gdename = "kernenried" if gdename == "kemenried" & PLZ4 == 3309
replace gdename = "kernenried" if gdename == "kemenrled" & PLZ4 == 3309
replace gdename = "kernenried" if gdename == "kernenrled" & PLZ4 == 3309
replace gdename = "baetterkinden" if gdename == "kraeiligen" & PLZ4 == 3315
replace gdename = "baetterkinden" if gdename == "kraeillgen" & PLZ4 == 3315
replace gdename = "baetterkinden" if gdename == "kraelligen" & PLZ4 == 3315
replace gdename = "baetterkinden" if gdename == "kraellligen" & PLZ4 == 3315
replace gdename = "baetterkinden" if gdename == "krailigen" & PLZ4 == 3315
replace gdename = "iffwil" if gdename == "lffwil" & PLZ4 == 3305
replace gdename = "zuzwil" if gdename == "zuzwll" & PLZ4 == 3303
replace gdename = "zuzwil" if gdename == "zuzwll (be)" & PLZ4 == 3303
replace gdename = "zuzwil" if gdename == "zuzwll be" & PLZ4 == 3303
replace gdename = "wolhusen" if gdename == "wolhi" & PLZ4 == 6110
replace gdename = "wolhusen" if gdename == "wolhusen markfjwerthenstein" & PLZ4 == 6110
replace gdename = "wolhusen" if gdename == "wolhusen markt" & PLZ4 == 6110
replace gdename = "wolhusen" if gdename == "wolhusen markt/ werthenstein" & PLZ4 == 6110
replace gdename = "wolhusen" if gdename == "wolhusen markt/werthenstein" & PLZ4 == 6110
replace gdename = "wolhusen" if gdename == "wolhusen markt/werthensteln" & PLZ4 == 6110
replace gdename = "wolhusen" if gdename == "wolhusen/ werthenstein" & PLZ4 == 6110
replace gdename = "bougy villars" if gdename == "bougv villars" & PLZ4 == 1172
replace gdename = "bougy villars" if gdename == "bougy" & PLZ4 == 1172
replace gdename = "bougy villars" if gdename == "bougy viliars" & PLZ4 == 1172
replace gdename = "bougy villars" if gdename == "bougy villard" & PLZ4 == 1172
replace gdename = "fechy" if gdename == "fechy dessus" & PLZ4 == 1173
replace gdename = "aubonne" if gdename == "plzy" & PLZ4 == 1174
replace gdename = "saint livres" if gdename == "st livres" & PLZ4 == 1176
replace gdename = "saint livres" if gdename == "st. livres" & PLZ4 == 1176
replace gdename = "thalwil" if gdename == "gattikon" & PLZ4 == 8139
replace gdename = "seftigen" if gdename == "saftigen" & PLZ4 == 3662
replace gdename = "seftigen" if gdename == "seftlgen" & PLZ4 == 3662
replace gdename = "seftigen" if gdename == "seifigen" & PLZ4 == 3662
replace gdename = "seftigen" if gdename == "seltigen" & PLZ4 == 3662
replace gdename = "seeberg" if gdename == "hermiewil" & PLZ4 == 3475
replace gdename = "seeberg" if gdename == "hermiswll" & PLZ4 == 3475
replace gdename = "seeberg" if gdename == "hermlswil" & PLZ4 == 3475
replace gdename = "seeberg" if gdename == "riedtwal" & PLZ4 == 3475
replace gdename = "seeberg" if gdename == "riedtwii" & PLZ4 == 3475
replace gdename = "seeberg" if gdename == "riedtwil" & PLZ4 == 3475
replace gdename = "seeberg" if gdename == "rietwil" & PLZ4 == 3475
replace gdename = "seeberg" if gdename == "rledtwil" & PLZ4 == 3475
replace gdename = "thunstetten" if gdename == "buetzberg" & PLZ4 == 4922
replace gdename = "oberdiessbach" if gdename == "oberdaessbach" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "oberdiesabach" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "oberdiesbach" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "oberdiessbagh" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "oberdlesabach" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "oberdlessbach" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "oberdliessbach" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "oberdã­essbach" & PLZ4 == 3672
replace gdename = "oberdiessbach" if gdename == "obergiessbach" & PLZ4 == 3672
replace gdename = "boltigen" if gdename == "boitigen/weissenbach" & PLZ4 == 3766
replace gdename = "boltigen" if gdename == "boltigen/weissenbach" & PLZ4 == 3766
replace gdename = "boltigen" if gdename == "weissenbach" & PLZ4 == 3766
replace gdename = "boltigen" if gdename == "weissenbach boltigen" & PLZ4 == 3766
replace gdename = "boltigen" if gdename == "weissenbach/boltigen" & PLZ4 == 3766
replace gdename = "boltigen" if gdename == "welissenbach/boltigen" & PLZ4 == 3766
replace gdename = "boltigen" if gdename == "welssenbach" & PLZ4 == 3766
replace gdename = "le noirmont" if gdename == "le noir mont" & PLZ4 == 2340
replace gdename = "le noirmont" if gdename == "le nolrmont" & PLZ4 == 2340
replace gdename = "le noirmont" if gdename == "noirmont" & PLZ4 == 2340
replace gdename = "le noirmont" if gdename == "noirmont le" & PLZ4 == 2340
replace gdename = "le noirmont" if gdename == "nolrmont" & PLZ4 == 2340
replace gdename = "glarus nord" if gdename == "obstaiden" & PLZ4 == 8758
replace gdename = "glarus nord" if gdename == "obstalgen" & PLZ4 == 8758
replace gdename = "ollon" if gdename == "chaetiliens" & PLZ4 == 1610
replace gdename = "ollon" if gdename == "chaetillens" & PLZ4 == 1610
replace gdename = "ollon" if gdename == "chaetlllens" & PLZ4 == 1610
replace gdename = "ollon" if gdename == "chatiilens" & PLZ4 == 1610
replace gdename = "ollon" if gdename == "chatiliens" & PLZ4 == 1610
replace gdename = "ollon" if gdename == "chatitlens" & PLZ4 == 1610
replace gdename = "ollon" if gdename == "chatlllens" & PLZ4 == 1610
replace gdename = "ollon" if gdename == "chetillens" & PLZ4 == 1610
replace gdename = "ollon" if gdename == "vulbroye" & PLZ4 == 1610
replace gdename = "les brenets" if gdename == "fretes" & PLZ4 == 2416
replace gdename = "la chaux du milieu" if gdename == "la chaux du millieu" & PLZ4 == 2405
replace gdename = "la sagne" if gdename == "la combe boudry" & PLZ4 == 2314
replace gdename = "la chaux du milieu" if gdename == "le cachot" & PLZ4 == 2405
replace gdename = "la chaux du milieu" if gdename == "le cachot/la chaux du milieu" & PLZ4 == 2405
replace gdename = "les brenets" if gdename == "les fretes" & PLZ4 == 2416
replace gdename = "trub" if gdename == "kroeschenbrunnen" & PLZ4 == 3555
replace gdename = "trub" if gdename == "kroeschenbrunnen/trub" & PLZ4 == 3555
replace gdename = "lyssach" if gdename == "lissach" & PLZ4 == 3421
replace gdename = "lyssach" if gdename == "llssach" & PLZ4 == 3421
replace gdename = "rueti bei lyssach" if gdename == "rueti b. lyssach" & PLZ4 == 3421
replace gdename = "rueti bei lyssach" if gdename == "rueti bel lyssach" & PLZ4 == 3421
replace gdename = "rueti bei lyssach" if gdename == "rueti lyssach" & PLZ4 == 3421
replace gdename = "rueti bei lyssach" if gdename == "ruetl lyssach" & PLZ4 == 3421
replace gdename = "rueti bei lyssach" if gdename == "ruti bei lyssach" & PLZ4 == 3421
replace gdename = "chexbres" if gdename == "cherbres" & PLZ4 == 1071
replace gdename = "chexbres" if gdename == "cheschrez" & PLZ4 == 1071
replace gdename = "chexbres" if gdename == "chexbrea" & PLZ4 == 1071
replace gdename = "chexbres" if gdename == "chexbree" & PLZ4 == 1071
replace gdename = "chexbres" if gdename == "chexbres vd" & PLZ4 == 1071
replace gdename = "puidoux" if gdename == "cremieres" & PLZ4 == 1071
replace gdename = "puidoux" if gdename == "cremieres chexbres" & PLZ4 == 1071
replace gdename = "puidoux" if gdename == "cremieres sur chexbres" & PLZ4 == 1071
replace gdename = "puidoux" if gdename == "lignieres sur chexbres" & PLZ4 == 1071
replace gdename = "boudry" if gdename == "perreux" & PLZ4 == 2071
replace gdename = "giswil" if gdename == "grossteil" & PLZ4 == 6074
replace gdename = "giswil" if gdename == "grossteil/giswil" & PLZ4 == 6074
replace gdename = "giswil" if gdename == "grosstell" & PLZ4 == 6074
replace gdename = "entlebuch" if gdename == "finsterwaid" & PLZ4 == 6162
replace gdename = "entlebuch" if gdename == "finsterwald" & PLZ4 == 6162
replace gdename = "entlebuch" if gdename == "finsterwald b. entlebuch" & PLZ4 == 6162
replace gdename = "entlebuch" if gdename == "finsterwald bei entiebuch" & PLZ4 == 6162
replace gdename = "entlebuch" if gdename == "finsterwald bei entlebuch" & PLZ4 == 6162
replace gdename = "entlebuch" if gdename == "finsterwald bel entlebuch" & PLZ4 == 6162
replace gdename = "entlebuch" if gdename == "finsterwald/ entlebuch" & PLZ4 == 6162
replace gdename = "entlebuch" if gdename == "finsterwald/entlebuch" & PLZ4 == 6162
replace gdename = "bregaglia" if gdename == "bonde" & PLZ4 == 7606
replace gdename = "bregaglia" if gdename == "casaccia" & PLZ4 == 7605
replace gdename = "bregaglia" if gdename == "casaccla" & PLZ4 == 7605
replace gdename = "bregaglia" if gdename == "promontogno" & PLZ4 == 7606
replace gdename = "bregaglia" if gdename == "stampa maloja" & PLZ4 == 7516
replace gdename = "bofflens" if gdename == "boffiens" & PLZ4 == 1353
replace gdename = "les clees" if gdename == "la russilie/les clees" & PLZ4 == 1356
replace gdename = "les clees" if gdename == "la russille vd" & PLZ4 == 1356
replace gdename = "les clees" if gdename == "la russille/les clees" & PLZ4 == 1356
replace gdename = "l'abergement " if gdename == "le" & PLZ4 == 1355
replace gdename = "l'abergement " if gdename == "le vail loud l abergement" & PLZ4 == 1355
replace gdename = "l'abergement " if gdename == "le vailloud l abergement" & PLZ4 == 1355
replace gdename = "les clees" if gdename == "les ciees" & PLZ4 == 1356
replace gdename = "lignerolle" if gdename == "ligneroile" & PLZ4 == 1357
replace gdename = "lignerolle" if gdename == "lignerolie" & PLZ4 == 1357
replace gdename = "montcherand" if gdename == "orbe montcherand" & PLZ4 == 1354
replace gdename = "rances" if gdename == "rances vd" & PLZ4 == 1439
replace gdename = "valeyres sous rances" if gdename == "vaieyres sous rances" & PLZ4 == 1358
replace gdename = "valeyres sous rances" if gdename == "valeyers sous rances" & PLZ4 == 1358
replace gdename = "valeyres sous rances" if gdename == "valeyres sous pances" & PLZ4 == 1358
replace gdename = "valeyres sous rances" if gdename == "valeyres/rances" & PLZ4 == 1358
replace gdename = "courrendlin" if gdename == "ccurrendlin" & PLZ4 == 2830
replace gdename = "courrendlin" if gdename == "courrendiin" & PLZ4 == 2830
replace gdename = "courrendlin" if gdename == "courrendiln" & PLZ4 == 2830
replace gdename = "courrendlin" if gdename == "courrendli" & PLZ4 == 2830
replace gdename = "courrendlin" if gdename == "courrendlln" & PLZ4 == 2830
replace gdename = "courrendlin" if gdename == "courrer dlin" & PLZ4 == 2830
replace gdename = "uesslingen buch" if gdename == "buch b. frauenfeld" & PLZ4 == 8524
replace gdename = "uesslingen buch" if gdename == "buch bei frauenfeld" & PLZ4 == 8524
replace gdename = "uesslingen buch" if gdename == "buch bel frauenfeld" & PLZ4 == 8524
replace gdename = "bubikon" if gdename == "wolfhausen" & PLZ4 == 8608
replace gdename = "innertkirchen" if gdename == "nessental" & PLZ4 == 3863
replace gdename = "wolfwil" if gdename == "waltwil" & PLZ4 == 4628
replace gdename = "wolfwil" if gdename == "woifwil" & PLZ4 == 4628
replace gdename = "wolfwil" if gdename == "woifwil be" & PLZ4 == 4628
replace gdename = "wolfwil" if gdename == "woitwil" & PLZ4 == 4628
replace gdename = "wolfwil" if gdename == "wolf wil" & PLZ4 == 4628
replace gdename = "wolfwil" if gdename == "wolfwii" & PLZ4 == 4628
replace gdename = "wolfwil" if gdename == "wolfwil be" & PLZ4 == 4628
replace gdename = "wolfwil" if gdename == "wolfwll" & PLZ4 == 4628
replace gdename = "wolfwil" if gdename == "wolfwll be" & PLZ4 == 4628
replace gdename = "wolfwil" if gdename == "woltwil" & PLZ4 == 4628
replace gdename = "lungern" if gdename == "buerglen ow" & PLZ4 == 6078
replace gdename = "lungern" if gdename == "buerglen ow/lungern" & PLZ4 == 6078
replace gdename = "lungern" if gdename == "buergten ow" & PLZ4 == 6078
replace gdename = "lungern" if gdename == "buergten ow/lungern" & PLZ4 == 6078
replace gdename = "lungern" if gdename == "kaiserstuhl ow" & PLZ4 == 6078
replace gdename = "lungern" if gdename == "kaiserstuhl ow/lungern" & PLZ4 == 6078
replace gdename = "lungern" if gdename == "kalserstuhl" & PLZ4 == 6078
replace gdename = "lungern" if gdename == "kalserstuhl ow" & PLZ4 == 6078
replace gdename = "entlebuch" if gdename == "rengg" & PLZ4 == 6162
replace gdename = "entlebuch" if gdename == "rengg entiebuch" & PLZ4 == 6162
replace gdename = "entlebuch" if gdename == "rengg entlebuch" & PLZ4 == 6162
replace gdename = "arth" if gdename == "rigi" & PLZ4 == 6410
replace gdename = "arth" if gdename == "rigi staffel" & PLZ4 == 6410
replace gdename = "arth" if gdename == "rigi staffel/arth" & PLZ4 == 6410
replace gdename = "arth" if gdename == "rigl staffel" & PLZ4 == 6410
replace gdename = "bischofszell" if gdename == "schweizershoiz" & PLZ4 == 9223
replace gdename = "hauptwil gottshaus" if gdename == "st. pelagiberg" & PLZ4 == 9225
replace gdename = "hauptwil gottshaus" if gdename == "st. pelagiberg gotteshaus" & PLZ4 == 9225
replace gdename = "hauptwil gottshaus" if gdename == "st. pelaglberg gotteshaus" & PLZ4 == 9225
replace gdename = "hauptwil gottshaus" if gdename == "walen gottshaus" & PLZ4 == 9225
replace gdename = "hauptwil gottshaus" if gdename == "wilen gottshaus" & PLZ4 == 9225
replace gdename = "hauptwil gottshaus" if gdename == "wlien gottshaus" & PLZ4 == 9225
replace gdename = "bagnes" if gdename == "verbier" & PLZ4 == 1936
replace gdename = "bagnes" if gdename == "verbier bagnes" & PLZ4 == 1936
replace gdename = "bagnes" if gdename == "verbler" & PLZ4 == 1936
replace gdename = "fontenais" if gdename == "bressancourt" & PLZ4 == 2904
replace gdename = "la baroche" if gdename == "charmoille ju" & PLZ4 == 2947
replace gdename = "la baroche" if gdename == "charmollle" & PLZ4 == 2947
replace gdename = "grandfontaine" if gdename == "grandfontalne" & PLZ4 == 2908
replace gdename = "la baroche" if gdename == "mlecourt" & PLZ4 == 2946
replace gdename = "haute ajoie" if gdename == "reciere" & PLZ4 == 2912
replace gdename = "kaisten" if gdename == "ittersthal" & PLZ4 == 5083
replace gdename = "kaisten" if gdename == "lttenthal" & PLZ4 == 5083
replace gdename = "mettauertal" if gdename == "meltau" & PLZ4 == 5274
replace gdename = "mettauertal" if gdename == "menthe" & PLZ4 == 5274
replace gdename = "laufenburg" if gdename == "rheinsulz" & PLZ4 == 5085
replace gdename = "laufenburg" if gdename == "sulz b. laufenburg" & PLZ4 == 5085
replace gdename = "laufenburg" if gdename == "sulz bei laufenburg" & PLZ4 == 5085
replace gdename = "laufenburg" if gdename == "sulz bel laufenburg" & PLZ4 == 5085
replace gdename = "mettauertal" if gdename == "wii" & PLZ4 == 5276
replace gdename = "mettauertal" if gdename == "wll ag" & PLZ4 == 5276
replace gdename = "zernez" if gdename == "brail" & PLZ4 == 7527
replace gdename = "s chanf" if gdename == "cinuos chel" & PLZ4 == 7526
replace gdename = "scuol" if gdename == "guards" & PLZ4 == 7545
replace gdename = "la punt chamues ch" if gdename == "la punt" & PLZ4 == 7522
replace gdename = "la punt chamues ch" if gdename == "la punt chamues ch*" & PLZ4 == 7522
replace gdename = "la punt chamues ch" if gdename == "la punt charnues ch" & PLZ4 == 7522
replace gdename = "madulain" if gdename == "madulaln" & PLZ4 == 7523
replace gdename = "waeldi" if gdename == "engwilen" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "engwllen" & PLZ4 == 8564
replace gdename = "waeldi" if gdename == "sonterswil" & PLZ4 == 8564
replace gdename = "lausanne" if gdename == "chalet a gobet" & PLZ4 == 1000
replace gdename = "lausanne" if gdename == "chalet ae gobet" & PLZ4 == 1000
replace gdename = "lausanne" if gdename == "le chalet a gobet" & PLZ4 == 1000
replace gdename = "lausanne" if gdename == "le chalet ae gobet" & PLZ4 == 1000
replace gdename = "lausanne" if gdename == "lp chalet a gobet" & PLZ4 == 1000
replace gdename = "alchenstorf" if gdename == "aichenstorf" & PLZ4 == 3473
replace gdename = "alchenstorf" if gdename == "alchensdorf" & PLZ4 == 3473
replace gdename = "alchenstorf" if gdename == "alchenstorf be" & PLZ4 == 3473
replace gdename = "ochlenberg" if gdename == "oschwand" & PLZ4 == 3476
replace gdename = "wynigen" if gdename == "rueedisbach" & PLZ4 == 3474
replace gdename = "walliswil bei wangen" if gdename == "wailiswil bei wangen" & PLZ4 == 3377
replace gdename = "walliswil bei wangen" if gdename == "walliswal bei wangen" & PLZ4 == 3377
replace gdename = "walliswil bei wangen" if gdename == "walliswal bel wangen" & PLZ4 == 3377
replace gdename = "walliswil bei wangen" if gdename == "walliswil" & PLZ4 == 3377
replace gdename = "walliswil bei wangen" if gdename == "walliswil b. wangen" & PLZ4 == 3377
replace gdename = "walliswil bei wangen" if gdename == "walliswil bel wangen" & PLZ4 == 3377
replace gdename = "walliswil bei wangen" if gdename == "walliswll bei wangen" & PLZ4 == 3377
replace gdename = "walliswil bei wangen" if gdename == "walllawil bel wangen" & PLZ4 == 3377
replace gdename = "walliswil bei wangen" if gdename == "walllswil bel wangen" & PLZ4 == 3377
replace gdename = "walliswil bei wangen" if gdename == "walllswlil bei wangen" & PLZ4 == 3377
replace gdename = "kiesen" if gdename == "klesen" & PLZ4 == 3629
replace gdename = "oppligen" if gdename == "oppiigen" & PLZ4 == 3629
replace gdename = "oppligen" if gdename == "oppiingen" & PLZ4 == 3629
replace gdename = "oppligen" if gdename == "oppilgen" & PLZ4 == 3629
replace gdename = "oppligen" if gdename == "oppl gen" & PLZ4 == 3629
replace gdename = "oppligen" if gdename == "oppligen be" & PLZ4 == 3629
replace gdename = "oppligen" if gdename == "opplingen" & PLZ4 == 3629
replace gdename = "oppligen" if gdename == "oppllgen" & PLZ4 == 3629
replace gdename = "menznau" if gdename == "twenenegg" & PLZ4 == 6122
replace gdename = "menznau" if gdename == "twerenegg" & PLZ4 == 6122
replace gdename = "thalwil" if gdename == "gattlkon" & PLZ4 == 8136
replace gdename = "uetendorf" if gdename == "uetendort" & PLZ4 == 3661
replace gdename = "zuerich" if gdename == "uetliberg" & PLZ4 == 8143
replace gdename = "zuerich" if gdename == "uettiberg" & PLZ4 == 8143
replace gdename = "hausen am albis" if gdename == "sihibrugg" & PLZ4 == 8915
replace gdename = "hausen am albis" if gdename == "sihlbrugg" & PLZ4 == 8915
replace gdename = "hausen am albis" if gdename == "sihlbrugg station" & PLZ4 == 8915
replace gdename = "hausen am albis" if gdename == "slhibrugg" & PLZ4 == 8915
replace gdename = "hausen am albis" if gdename == "slhlbrugg" & PLZ4 == 8915
replace gdename = "belphraon" if gdename == "balprahon" & PLZ4 == 2744
replace gdename = "belphraon" if gdename == "beiprahon" & PLZ4 == 2744
replace gdename = "belphraon" if gdename == "belpranon" & PLZ4 == 2744
replace gdename = "eschert" if gdename == "escihert" & PLZ4 == 2743
replace gdename = "perrefitte" if gdename == "perefitte" & PLZ4 == 2742
replace gdename = "perrefitte" if gdename == "perrefltte" & PLZ4 == 2742
replace gdename = "perrefitte" if gdename == "perretitte" & PLZ4 == 2742
replace gdename = "perrefitte" if gdename == "peyrefitte" & PLZ4 == 2742
replace gdename = "woelflinswil" if gdename == "wdlflinswll" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woelfiinswil" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woelflanswil" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woelflinswii" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woelflinswl" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woelflinswll" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woelfllinswil" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woelfllnswil" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woelfllnswll" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woeltllnswil" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "wolflinswil" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woltiinswil" & PLZ4 == 5063
replace gdename = "woelflinswil" if gdename == "woltlinswil" & PLZ4 == 5063
replace gdename = "menznau" if gdename == "fontannen bei wolhusen" & PLZ4 == 6122
replace gdename = "menznau" if gdename == "fontannen bel wolhusen" & PLZ4 == 6122
replace gdename = "wolhusen" if gdename == "steinhuserberg" & PLZ4 == 6114
replace gdename = "gossau (sg)" if gdename == "gassau" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gassau sg" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gaussa sg" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gcssau sg" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gessau (sg)" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "goaaau sg" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "goasau (sg)" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "goss au sg" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gossaeu (sg)" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gossau sa" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gossau sg e" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gossaue" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gossaus" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gossen sg" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gowssau" & PLZ4 == 9200
replace gdename = "gossau (sg)" if gdename == "gussau sg" & PLZ4 == 9200
replace gdename = "walzenhausen" if gdename == "lachen ar" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "lachen walzenhausen" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "lachen/ar" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "lachen/walzenhausen" & PLZ4 == 9428
replace gdename = "belmont broye" if gdename == "dampierre vd" & PLZ4 == 1563
replace gdename = "belmont broye" if gdename == "domplerre vd" & PLZ4 == 1563
replace gdename = "belmont broye" if gdename == "granges de dompierre" & PLZ4 == 1563
replace gdename = "belmont broye" if gdename == "granges de domplerre" & PLZ4 == 1563
replace gdename = "belmont broye" if gdename == "granges de=dompierre" & PLZ4 == 1563
replace gdename = "saignelegier" if gdename == "cernievillers" & PLZ4 == 2353
replace gdename = "saignelegier" if gdename == "cernlevlllers" & PLZ4 == 2353
replace gdename = "saignelegier" if gdename == "lea pommerats" & PLZ4 == 2353
replace gdename = "saignelegier" if gdename == "les pommerais" & PLZ4 == 2353
replace gdename = "saignelegier" if gdename == "les pommerets" & PLZ4 == 2353
replace gdename = "saignelegier" if gdename == "pommerats" & PLZ4 == 2353
replace gdename = "inkwil" if gdename == "inkwil roethenbach" & PLZ4 == 3375
replace gdename = "inkwil" if gdename == "inkwll" & PLZ4 == 3375
replace gdename = "inkwil" if gdename == "lnkwil" & PLZ4 == 3375
replace gdename = "inkwil" if gdename == "lnkwil (roethenbach)" & PLZ4 == 3375
replace gdename = "inkwil" if gdename == "lnkwll (roethenbach)" & PLZ4 == 3375
replace gdename = "schuepfheim" if gdename == "kiusstalden /schuepfheim" & PLZ4 == 6170
replace gdename = "schuepfheim" if gdename == "klusstaiden" & PLZ4 == 6170
replace gdename = "schuepfheim" if gdename == "klusstalden" & PLZ4 == 6170
replace gdename = "schuepfheim" if gdename == "klusstalden /schuepf heim" & PLZ4 == 6170
replace gdename = "schuepfheim" if gdename == "klusstalden /schuepfheim" & PLZ4 == 6170
replace gdename = "schuepfheim" if gdename == "ktusatalden" & PLZ4 == 6170
replace gdename = "thayngen" if gdename == "astdorf" & PLZ4 == 8243
replace gdename = "doerflingen" if gdename == "dorflingen" & PLZ4 == 8239
replace gdename = "stetten" if gdename == "steffen" & PLZ4 == 8234
replace gdename = "stetten" if gdename == "steffen (sh)" & PLZ4 == 8234
replace gdename = "stetten" if gdename == "steffen sh" & PLZ4 == 8234
replace gdename = "stetten" if gdename == "stetten 8h" & PLZ4 == 8234
replace gdename = "buchs" if gdename == "raefis" & PLZ4 == 9470
replace gdename = "buchs" if gdename == "raefis buchs" & PLZ4 == 9470
replace gdename = "buchs" if gdename == "raefls" & PLZ4 == 9470
replace gdename = "muemliswil ramiswil" if gdename == "ramiswil" & PLZ4 == 4717
replace gdename = "weiach" if gdename == "walach" & PLZ4 == 8187
replace gdename = "weiach" if gdename == "welach" & PLZ4 == 8187
replace gdename = "weiach" if gdename == "welsch" & PLZ4 == 8187
replace gdename = "horgen" if gdename == "horgenberg" & PLZ4 == 8815
replace gdename = "horgen" if gdename == "horgerberg" & PLZ4 == 8815
replace gdename = "waedenswil" if gdename == "hueften" & PLZ4 == 8825
replace gdename = "waedenswil" if gdename == "schoenenherg" & PLZ4 == 8824
replace gdename = "waedenswil" if gdename == "schonenberg" & PLZ4 == 8824
replace gdename = "waedenswil" if gdename == "schonenberg (zh)" & PLZ4 == 8824
replace gdename = "waedenswil" if gdename == "schã³nenberg (zh)" & PLZ4 == 8824
replace gdename = "waedenswil" if gdename == "tanne schoenenberg" & PLZ4 == 8824
replace gdename = "waedenswil" if gdename == "tanne waedenswil" & PLZ4 == 8824
replace gdename = "wuppenau" if gdename == "hagenbuch am nollen" & PLZ4 == 9514
replace gdename = "wuppenau" if gdename == "hagenbuch tg" & PLZ4 == 9514
replace gdename = "schoenholzerswilen" if gdename == "toos" & PLZ4 == 8577
replace gdename = "schoenholzerswilen" if gdename == "toos/schoenholzerswilen" & PLZ4 == 8577
replace gdename = "ollon" if gdename == "chesieres ollon" & PLZ4 == 1885
replace gdename = "ollon" if gdename == "chesieres olson" & PLZ4 == 1885
replace gdename = "ollon" if gdename == "st triphon" & PLZ4 == 1867
replace gdename = "ollon" if gdename == "st trlphon" & PLZ4 == 1867
replace gdename = "brenzikofen" if gdename == "brenziikofen" & PLZ4 == 3671
replace gdename = "brenzikofen" if gdename == "brenzlkofen" & PLZ4 == 3671
replace gdename = "rapperswil" if gdename == "dieterswil" & PLZ4 == 3256
replace gdename = "herbligen" if gdename == "herbiigen" & PLZ4 == 3671
replace gdename = "herbligen" if gdename == "herbilgen" & PLZ4 == 3671
replace gdename = "herbligen" if gdename == "herbligen be" & PLZ4 == 3671
replace gdename = "herbligen" if gdename == "herbligen brenzikofen" & PLZ4 == 3671
replace gdename = "herbligen" if gdename == "herblingen" & PLZ4 == 3671
replace gdename = "herbligen" if gdename == "herbllgen" & PLZ4 == 3671
replace gdename = "berikon" if gdename == "mutschellen" & PLZ4 == 8965
replace gdename = "zufikon" if gdename == "mutschellen/zufikon" & PLZ4 == 5621
replace gdename = "kaisten" if gdename == "kaisten ag" & PLZ4 == 5082
replace gdename = "kaisten" if gdename == "kalsten" & PLZ4 == 5082
replace gdename = "forel (lavaux)" if gdename == "farel" & PLZ4 == 1072
replace gdename = "forel (lavaux)" if gdename == "farel (lavaux)" & PLZ4 == 1072
replace gdename = "forel (lavaux)" if gdename == "forai (lavauz)" & PLZ4 == 1072
replace gdename = "forel (lavaux)" if gdename == "forei (lavaux)" & PLZ4 == 1072
replace gdename = "forel (lavaux)" if gdename == "forel laveaux" & PLZ4 == 1072
replace gdename = "forel (lavaux)" if gdename == "forer lavaux" & PLZ4 == 1072
replace gdename = "forel (lavaux)" if gdename == "foret (lavaux)" & PLZ4 == 1072
replace gdename = "forel (lavaux)" if gdename == "foret (lavauz)" & PLZ4 == 1072
replace gdename = "forel (lavaux)" if gdename == "foret (uvaux)" & PLZ4 == 1072
replace gdename = "forel (lavaux)" if gdename == "foret/lavaux" & PLZ4 == 1072
replace gdename = "forel (lavaux)" if gdename == "lavaux" & PLZ4 == 1072
replace gdename = "jaun" if gdename == "im fang" & PLZ4 == 1656
replace gdename = "hasle" if gdename == "biembach" & PLZ4 == 3419
replace gdename = "hasle" if gdename == "blembach" & PLZ4 == 3419
replace gdename = "rueegsau" if gdename == "fiueegsau" & PLZ4 == 3417
replace gdename = "rueegsau" if gdename == "rueegaau" & PLZ4 == 3417
replace gdename = "rueegsau" if gdename == "ruegsau" & PLZ4 == 3417
replace gdename = "glattfelden" if gdename == "giattfelden/zweidien" & PLZ4 == 8192
replace gdename = "glattfelden" if gdename == "glattfelden/zweidlen" & PLZ4 == 8192
replace gdename = "glattfelden" if gdename == "glattfelden/zweldlen" & PLZ4 == 8192
replace gdename = "glattfelden" if gdename == "glattfeldern/zweidlen" & PLZ4 == 8192
replace gdename = "schwytz" if gdename == "rickenbach b. sch" & PLZ4 == 6432
replace gdename = "glattfelden" if gdename == "zweidien" & PLZ4 == 8192
replace gdename = "glattfelden" if gdename == "zweidlen" & PLZ4 == 8192
replace gdename = "glattfelden" if gdename == "zweldlen" & PLZ4 == 8192
replace gdename = "weisslingen" if gdename == "theilingen" & PLZ4 == 8484
replace gdename = "weisslingen" if gdename == "theilingen/weisslingen" & PLZ4 == 8484
replace gdename = "weisslingen" if gdename == "thellingen" & PLZ4 == 8484
replace gdename = "thundorf" if gdename == "lustdorf" & PLZ4 == 8512
replace gdename = "thundorf" if gdename == "wetzikon tg" & PLZ4 == 8512
replace gdename = "thundorf" if gdename == "wetzlkon tg" & PLZ4 == 8512
replace gdename = "meilen" if gdename == "feldmeilen" & PLZ4 == 8706
replace gdename = "val de ruz" if gdename == "st martin ne" & PLZ4 == 2054
replace gdename = "val de ruz" if gdename == "st martin ne/ chezard st martin" & PLZ4 == 2054
replace gdename = "val de ruz" if gdename == "st. martin ne" & PLZ4 == 2054
replace gdename = "ruederswil" if gdename == "randflueh" & PLZ4 == 3437
replace gdename = "ruederswil" if gdename == "ranflueh" & PLZ4 == 3437
replace gdename = "ruederswil" if gdename == "ranfluh" & PLZ4 == 3437
replace gdename = "ruederswil" if gdename == "ranflãºh" & PLZ4 == 3437
replace gdename = "ruederswil" if gdename == "schwanden i.e" & PLZ4 == 3433
replace gdename = "leibstadt" if gdename == "lelbatadt" & PLZ4 == 5325
replace gdename = "leibstadt" if gdename == "lelbstadt" & PLZ4 == 5325
replace gdename = "staffelbach" if gdename == "wittwil" & PLZ4 == 5053
replace gdename = "staffelbach" if gdename == "wittwil staffelbach" & PLZ4 == 5053
replace gdename = "staffelbach" if gdename == "wittwilstaffelbach" & PLZ4 == 5053
replace gdename = "staffelbach" if gdename == "wittwll" & PLZ4 == 5053
replace gdename = "breggia" if gdename == "bruzeila" & PLZ4 == 6837
replace gdename = "breggia" if gdename == "bruzelia" & PLZ4 == 6837
replace gdename = "breggia" if gdename == "cabblo" & PLZ4 == 6838
replace gdename = "breggia" if gdename == "canegglo" & PLZ4 == 6837
replace gdename = "breggia" if gdename == "gabbio" & PLZ4 == 6838
replace gdename = "breggia" if gdename == "mugglo" & PLZ4 == 6838
replace gdename = "breggia" if gdename == "segno" & PLZ4 == 6839
replace gdename = "bad zurzach" if gdename == "bad zurzach" & PLZ4 == 5330
replace gdename = "bad zurzach" if gdename == "zuerzach" & PLZ4 == 5330
replace gdename = "bad zurzach" if gdename == "zur.ach" & PLZ4 == 5330
replace gdename = "bad zurzach" if gdename == "zurzach ag" & PLZ4 == 5330
replace gdename = "bad zurzach" if gdename == "zurzadh" & PLZ4 == 5330
replace gdename = "bad zurzach" if gdename == "zurzech" & PLZ4 == 5330
replace gdename = "uttigen" if gdename == "uttingen" & PLZ4 == 3628
replace gdename = "uttigen" if gdename == "uttlgen" & PLZ4 == 3628
replace gdename = "amlikon bissegg" if gdename == "amilkon" & PLZ4 == 8514
replace gdename = "amlikon bissegg" if gdename == "amllkon" & PLZ4 == 8514
replace gdename = "thundorf" if gdename == "aufhofen thundorf" & PLZ4 == 8512
replace gdename = "amlikon bissegg" if gdename == "leutmerken" & PLZ4 == 8514
replace gdename = "koeniz" if gdename == "145 niederscherli koeniz" & PLZ4 == 3145
replace gdename = "basel" if gdename == "basei" & PLZ4 == 4000
replace gdename = "muri bei bern" if gdename == "guemligen" & PLZ4 == 3073
replace gdename = "moudon" if gdename == "bressonnaz" & PLZ4 == 1510
replace gdename = "moudon" if gdename == "bressonnaz/moudon" & PLZ4 == 1510
replace gdename = "moudon" if gdename == "bressonnez" & PLZ4 == 1510
replace gdename = "luthern" if gdename == "luthern bad" & PLZ4 == 6156
replace gdename = "luthern" if gdename == "luthern bad/luthern" & PLZ4 == 6156
replace gdename = "assens" if gdename == "assena" & PLZ4 == 1042
replace gdename = "bioley orjulaz" if gdename == "bioiey orjulaz" & PLZ4 == 1042
replace gdename = "bioley orjulaz" if gdename == "bloley magnoux" & PLZ4 == 1042
replace gdename = "bioley orjulaz" if gdename == "bloley orjulaz" & PLZ4 == 1042
replace gdename = "bournens" if gdename == "boumens" & PLZ4 == 1035
replace gdename = "fey" if gdename == "fey vd" & PLZ4 == 1044
replace gdename = "sullens" if gdename == "suilens" & PLZ4 == 1036
replace gdename = "sullens" if gdename == "suliens" & PLZ4 == 1036
replace gdename = "bavois" if gdename == "bavols" & PLZ4 == 1372
replace gdename = "bavois" if gdename == "bavons" & PLZ4 == 1372
replace gdename = "bavois" if gdename == "bevels" & PLZ4 == 1372
replace gdename = "chavornay" if gdename == "corcelles chavornay" & PLZ4 == 1374
replace gdename = "chavornay" if gdename == "corcelles sur chavomay" & PLZ4 == 1374
replace gdename = "chavornay" if gdename == "corcelles/chavornay" & PLZ4 == 1374
replace gdename = "goumoens" if gdename == "goumoen la ville" & PLZ4 == 1376
replace gdename = "haute sorne" if gdename == "berlincourt" & PLZ4 == 2854
replace gdename = "worb" if gdename == "walltenwll" & PLZ4 == 3076
replace gdename = "worb" if gdename == "watlenwil" & PLZ4 == 3076
replace gdename = "worb" if gdename == "wattenwii" & PLZ4 == 3076
replace gdename = "worb" if gdename == "wattenwil (wort))" & PLZ4 == 3076
replace gdename = "worb" if gdename == "wattenwll" & PLZ4 == 3076
replace gdename = "worb" if gdename == "wattenwlã¬l" & PLZ4 == 3076
replace gdename = "worb" if gdename == "wattwile" & PLZ4 == 3076
replace gdename = "langnau (be)" if gdename == "baerau i.e" & PLZ4 == 3550
replace gdename = "rumendingen" if gdename == "ruinendingen" & PLZ4 == 3472
replace gdename = "wynigen" if gdename == "wynigen be" & PLZ4 == 3472
replace gdename = "wynigen" if gdename == "wynlgen" & PLZ4 == 3472
replace gdename = "oberdiessbach" if gdename == "bleiken" & PLZ4 == 3674
replace gdename = "oberdiessbach" if gdename == "bleiken b. oberdiessbach" & PLZ4 == 3674
replace gdename = "oberdiessbach" if gdename == "bleiken bei oberdlessbach" & PLZ4 == 3674
replace gdename = "oberdiessbach" if gdename == "bleiken bel oberdiessbach" & PLZ4 == 3674
replace gdename = "oberdiessbach" if gdename == "bleiken bel oberdlessbach" & PLZ4 == 3674
replace gdename = "oberdiessbach" if gdename == "blelken bei oberdiessbach" & PLZ4 == 3674
replace gdename = "oberdiessbach" if gdename == "blelken bel oberdiessbach" & PLZ4 == 3674
replace gdename = "buchholterberg" if gdename == "wangelen b. oberdiessbach/buchholterberg" & PLZ4 == 3615
replace gdename = "buchholterberg" if gdename == "wangelen bei oberdiessbach" & PLZ4 == 3615
replace gdename = "buchholterberg" if gdename == "wangelen bei oberdiessbach/ buchholterberg" & PLZ4 == 3615
replace gdename = "buchholterberg" if gdename == "wangelen bei oberdiessbach/buchholterberg" & PLZ4 == 3615
replace gdename = "bregaglia" if gdename == "casaccia" & PLZ4 == 7605
replace gdename = "hasle bei burgdorf" if gdename == "schaffhausen im emmental" & PLZ4 == 3415
replace gdename = "hasle bei burgdorf" if gdename == "schafhausen be" & PLZ4 == 3415
replace gdename = "hasle bei burgdorf" if gdename == "schafhausen i.e" & PLZ4 == 3415
replace gdename = "hasle bei burgdorf" if gdename == "schafhausen im emmental/hasle bei burgdorf" & PLZ4 == 3415
replace gdename = "hasle bei burgdorf" if gdename == "schafhausen im emmental/haste bel burgdorf" & PLZ4 == 3415
replace gdename = "hasle bei burgdorf" if gdename == "schafhausen untergomerkinden" & PLZ4 == 3415
replace gdename = "hasle bei burgdorf" if gdename == "schoefhausen be" & PLZ4 == 3415
replace gdename = "hasle bei burgdorf" if gdename == "schoenhausen be" & PLZ4 == 3415
replace gdename = "laufenburg" if gdename == "rheinsulz" & PLZ4 == 5085
replace gdename = "auswil" if gdename == "answil" & PLZ4 == 4944
replace gdename = "auswil" if gdename == "auswll" & PLZ4 == 4944
replace gdename = "madiswil" if gdename == "gutenberg" & PLZ4 == 4932
replace gdename = "madiswil" if gdename == "gutenburg be" & PLZ4 == 4932
replace gdename = "madiswil" if gdename == "gutenburg/thunstetten" & PLZ4 == 4932
replace gdename = "madiswil" if gdename == "leimiswii" & PLZ4 == 4935
replace gdename = "madiswil" if gdename == "leimiswll" & PLZ4 == 4935
replace gdename = "madiswil" if gdename == "lelmiswil" & PLZ4 == 4935
replace gdename = "walterswil" if gdename == "walterswii (be)" & PLZ4 == 4942
replace gdename = "courgenay" if gdename == "courgena" & PLZ4 == 2950
replace gdename = "courgenay" if gdename == "courgenax" & PLZ4 == 2950
replace gdename = "courgenay" if gdename == "courgnay" & PLZ4 == 2950
replace gdename = "courgenay" if gdename == "courtemautruy" & PLZ4 == 2950
replace gdename = "luterbach" if gdename == "luterbach so" & PLZ4 == 4542
replace gdename = "luterbach" if gdename == "luterbach/so" & PLZ4 == 4542
replace gdename = "luterbach" if gdename == "luterbech" & PLZ4 == 4542
replace gdename = "maennedorf" if gdename == "mannedorf" & PLZ4 == 8708
replace gdename = "alpnach" if gdename == "piatus kulm" & PLZ4 == 6055
replace gdename = "alpnach" if gdename == "pilatus kulm" & PLZ4 == 6055
replace gdename = "alpnach" if gdename == "platus kulm" & PLZ4 == 6055
replace gdename = "alpnach" if gdename == "pllatus kulm" & PLZ4 == 6055
replace gdename = "vaz/obervaz" if gdename == "vaz obervaz" & PLZ4 == 7082
replace gdename = "vufflens le chateau" if gdename == ".134 vufflens le chaeteau" & PLZ4 == 1134
replace gdename = "saint prex" if gdename == "162 st prex" & PLZ4 == 1162
replace gdename = "lancy" if gdename == "212 grand lancy" & PLZ4 == 1212
replace gdename = "collonge bellerive" if gdename == "222 vesenaz" & PLZ4 == 1222
replace gdename = "veyrier" if gdename == "255 v eyrier" & PLZ4 == 1255
replace gdename = "chatel saint denis" if gdename == "618 chatel st denis" & PLZ4 == 1618
replace gdename = "pregny chambesy" if gdename == "avenue de tourney. 1292 chambesy" & PLZ4 == 1292
replace gdename = "visperterminen" if gdename == "visperterm" & PLZ4 == 3932
replace gdename = "kallnach" if gdename == "goiaten" & PLZ4 == 3207
replace gdename = "kallnach" if gdename == "gurlern" & PLZ4 == 3207
replace gdename = "ferenbalm" if gdename == "jerisberg ferenbaim" & PLZ4 == 3206
replace gdename = "ferenbalm" if gdename == "jerisberg ferenbalm" & PLZ4 == 3206
replace gdename = "vinelz" if gdename == "vineiz" & PLZ4 == 3234
replace gdename = "vinelz" if gdename == "vinetz" & PLZ4 == 3234
replace gdename = "wileroltigen" if gdename == "wilerottigen" & PLZ4 == 3207
replace gdename = "seedorf" if gdename == "baggwii seedorf be" & PLZ4 == 3267
replace gdename = "seedorf" if gdename == "baggwil" & PLZ4 == 3267
replace gdename = "seedorf" if gdename == "baggwil seedorf" & PLZ4 == 3267
replace gdename = "seedorf" if gdename == "baggwil seedorf be" & PLZ4 == 3267
replace gdename = "seedorf" if gdename == "baggwll seedorf be" & PLZ4 == 3267
replace gdename = "seedorf" if gdename == "frienisberg" & PLZ4 == 3267
replace gdename = "boezen" if gdename == "boezen ag" & PLZ4 == 5076
replace gdename = "boezen" if gdename == "boezen\" & PLZ4 == 5076
replace gdename = "boezen" if gdename == "bozen" & PLZ4 == 5076
replace gdename = "baden" if gdename == "daettwil" & PLZ4 == 5405
replace gdename = "baden" if gdename == "daettwil (ag)" & PLZ4 == 5405
replace gdename = "baden" if gdename == "daettwil baden" & PLZ4 == 5405
replace gdename = "baden" if gdename == "daettwil/baden" & PLZ4 == 5405
replace gdename = "baden" if gdename == "dattwil" & PLZ4 == 5405
replace gdename = "baden" if gdename == "dattwil baden" & PLZ4 == 5405
replace gdename = "edlibach" if gdename == "edllbach" & PLZ4 == 6313
replace gdename = "lindau" if gdename == "ottikon b. kemptthal" & PLZ4 == 8310
replace gdename = "lindau" if gdename == "ottikon bei kempttal" & PLZ4 == 8310
replace gdename = "lindau" if gdename == "ottikon bei kemptthal" & PLZ4 == 8310
replace gdename = "lindau" if gdename == "ottlkon bel kemptthal" & PLZ4 == 8310
replace gdename = "heimberg" if gdename == "halmberg" & PLZ4 == 3627
replace gdename = "heimberg" if gdename == "heimberg be" & PLZ4 == 3627
replace gdename = "heimberg" if gdename == "helmberg" & PLZ4 == 3627
replace gdename = "heimberg" if gdename == "holmberg" & PLZ4 == 3627
replace gdename = "gurtnellen" if gdename == "intschi" & PLZ4 == 6474
replace gdename = "gurtnellen" if gdename == "intschi gurtnellen" & PLZ4 == 6482
replace gdename = "gurtnellen" if gdename == "intschl" & PLZ4 == 6474
replace gdename = "gurtnellen" if gdename == "intschl gurtnellen" & PLZ4 == 6482
replace gdename = "lauerz" if gdename == "lauen" & PLZ4 == 6424
replace gdename = "gurtnellen" if gdename == "lntschi" & PLZ4 == 6474
replace gdename = "neunkirch" if gdename == "neunkiroh" & PLZ4 == 8213
replace gdename = "safiental" if gdename == "safien platz" & PLZ4 == 7107
replace gdename = "safiental" if gdename == "saften platz" & PLZ4 == 7107
replace gdename = "chapelle" if gdename == "chapelle sur oron" & PLZ4 == 1608
replace gdename = "oron" if gdename == "oran la ville" & PLZ4 == 1610
replace gdename = "oron" if gdename == "oron" & PLZ4 == 1610
replace gdename = "oron" if gdename == "oron ia ville" & PLZ4 == 1610
replace gdename = "oron" if gdename == "oron la viiie" & PLZ4 == 1610
replace gdename = "oron" if gdename == "own la ville" & PLZ4 == 1610
replace gdename = "steffisburg" if gdename == "steffisburg station" & PLZ4 == 3612
replace gdename = "steffisburg" if gdename == "stefflsburg" & PLZ4 == 3612
replace gdename = "seedorf" if gdename == "baggwil" & PLZ4 == 3267
replace gdename = "seedorf" if gdename == "baggwil/seedorf be" & PLZ4 == 3267
replace gdename = "zug" if gdename == "zugerberg" & PLZ4 == 6300
replace gdename = "zug" if gdename == "zugerberg/zug" & PLZ4 == 6300
replace gdename = "oberdorf" if gdename == "niederrickenbach" & PLZ4 == 6383
replace gdename = "oberdorf" if gdename == "niederrickenbach nw" & PLZ4 == 6383
replace gdename = "weinfelden" if gdename == "weerswilen" & PLZ4 == 8570
replace gdename = "baden" if gdename == "ruetihof" & PLZ4 == 5406
replace gdename = "baden" if gdename == "ruetihof b. baden" & PLZ4 == 5406
replace gdename = "baden" if gdename == "ruetihof baden" & PLZ4 == 5406
replace gdename = "baden" if gdename == "ruetihof/baden" & PLZ4 == 5406
replace gdename = "koelliken" if gdename == "kolliken" & PLZ4 == 5742
replace gdename = "oftringen" if gdename == "kuengoidingen" & PLZ4 == 4665
replace gdename = "oftringen" if gdename == "kuengoldingen" & PLZ4 == 4665
replace gdename = "neunkirch" if gdename == "sempach station" & PLZ4 == 6203
replace gdename = "centovalli" if gdename == "camedo/borgnone" & PLZ4 == 6659
replace gdename = "centovalli" if gdename == "rasa" & PLZ4 == 6655
replace gdename = "brusino arsizio" if gdename == "brusino" & PLZ4 == 6827
replace gdename = "brusino arsizio" if gdename == "brusino arsizlo" & PLZ4 == 6827
replace gdename = "brusino arsizio" if gdename == "brusio arizio" & PLZ4 == 6827
replace gdename = "brusino arsizio" if gdename == "brusio arsizio" & PLZ4 == 6827
replace gdename = "brusino arsizio" if gdename == "bruslno arsizio" & PLZ4 == 6827
replace gdename = "novazzano" if gdename == "brusata" & PLZ4 == 6883
replace gdename = "oberdorf" if gdename == "bueren/oberdorf nw" & PLZ4 == 6382
replace gdename = "bussnang" if gdename == "rothenhausen" & PLZ4 == 9565
replace gdename = "ollon" if gdename == "huemoz" & PLZ4 == 1884
replace gdename = "ollon" if gdename == "huemoz villars" & PLZ4 == 1884
replace gdename = "ollon" if gdename == "huemoz/ollon" & PLZ4 == 1884
replace gdename = "koeniz" if gdename == "wabern" & PLZ4 == 3084
replace gdename = "bienne" if gdename == "biel/bienne" & PLZ4 == 2500
replace gdename = "frick" if gdename == "frack" & PLZ4 == 5070
replace gdename = "frick" if gdename == "frick ag" & PLZ4 == 5070
replace gdename = "oberwil lieli" if gdename == "oberwil ag" & PLZ4 == 4074
replace gdename = "oberwil lieli" if gdename == "oberwll ag" & PLZ4 == 4074
replace gdename = "rickenbach" if gdename == "rickenbach ag" & PLZ4 == 8545
replace gdename = "rickenbach" if gdename == "rlckenbach" & PLZ4 == 8545
replace gdename = "samnaun" if gdename == "sammnaun" & PLZ4 == 7563
replace gdename = "stammheim" if gdename == "guntalingen" & PLZ4 == 8468
replace gdename = "chene bougeries" if gdename == "conches" & PLZ4 == 1231
replace gdename = "lueterkofen ichertswil" if gdename == "ichertswil" & PLZ4 == 4571
replace gdename = "lueterkofen ichertswil" if gdename == "ichertswll" & PLZ4 == 4571
replace gdename = "lueterkofen ichertswil" if gdename == "lchertswil" & PLZ4 == 4571
replace gdename = "risch" if gdename == "rotkreuz" & PLZ4 == 6343
replace gdename = "seelisberg" if gdename == "seeiisberg" & PLZ4 == 6377
replace gdename = "seelisberg" if gdename == "seeilsberg" & PLZ4 == 6377
replace gdename = "seelisberg" if gdename == "seellsberg" & PLZ4 == 6377
replace gdename = "tenero contra" if gdename == "tenero" & PLZ4 == 6598
replace gdename = "schwarzenburg" if gdename == "schwarzen burg wahiern" & PLZ4 == 3150
replace gdename = "schwarzenburg" if gdename == "schwarzen burg wahlern" & PLZ4 == 3150
replace gdename = "schwarzenburg" if gdename == "schwarzenburg wahlern" & PLZ4 == 3150
replace gdename = "steinmaur" if gdename == "obersteinmaur" & PLZ4 == 8162
replace gdename = "zumikon" if gdename == "zumlkon" & PLZ4 == 8126
replace gdename = "minusio" if gdename == "mlnuslo" & PLZ4 == 6648
replace gdename = "einsiedeln" if gdename == "trachslau" & PLZ4 == 8840
replace gdename = "einsiedeln" if gdename == "trachslau/einsiedeln" & PLZ4 == 8840
replace gdename = "einsiedeln" if gdename == "trachstau/einsiedeln" & PLZ4 == 8840
replace gdename = "lausanne" if gdename == "vers chez les blanc" & PLZ4 == 1000
replace gdename = "veytaux" if gdename == "grandchamp" & PLZ4 == 1820
replace gdename = "veytaux" if gdename == "veytaux chillon" & PLZ4 == 1820
replace gdename = "veytaux" if gdename == "veytaux/chillon" & PLZ4 == 1820
replace gdename = "bienne" if gdename == "biel/bienne" & PLZ4 == 2500
replace gdename = "gondiswil" if gdename == "gondiswii" & PLZ4 == 4955
replace gdename = "gondiswil" if gdename == "gondiswll" & PLZ4 == 4955
replace gdename = "gondiswil" if gdename == "gondlawil" & PLZ4 == 4955
replace gdename = "gondiswil" if gdename == "gondlswil" & PLZ4 == 4955
replace gdename = "gondiswil" if gdename == "gondlswll" & PLZ4 == 4955
replace gdename = "reiden" if gdename == "reidermoos" & PLZ4 == 6260
replace gdename = "reiden" if gdename == "reidermoos/beiden" & PLZ4 == 6260
replace gdename = "reiden" if gdename == "reidermoos/reiden" & PLZ4 == 6260
replace gdename = "steckborn" if gdename == "steckbom" & PLZ4 == 8266
replace gdename = "silvaplana" if gdename == "silvaplana/surlei" & PLZ4 == 7513
replace gdename = "silvaplana" if gdename == "slivaplana/surlei" & PLZ4 == 7513
replace gdename = "mettauertal" if gdename == "oberhofen zh" & PLZ4 == 5273
replace gdename = "mettauertal" if gdename == "oberhufen zh" & PLZ4 == 5273
replace gdename = "flims" if gdename == "017 flims dorf" & PLZ4 == 7017
replace gdename = "vaz/obervaz" if gdename == "077 valbella/vaz/obervaz" & PLZ4 == 7077
replace gdename = "davos" if gdename == "270 davos platz" & PLZ4 == 7270
replace gdename = "fideris" if gdename == "299 fi feris station" & PLZ4 == 7235
replace gdename = "collonge bellerive" if gdename == "ch.des epines. 1222 vesenaz" & PLZ4 == 1222
replace gdename = "uzwil" if gdename == "algetshausen" & PLZ4 == 9249
replace gdename = "bussy sur moudon" if gdename == "bussy moudon" & PLZ4 == 1514
replace gdename = "brot plamboz" if gdename == "petits ponts" & PLZ4 == 2318
replace gdename = "les bois" if gdename == "le boechet" & PLZ4 == 2336
replace gdename = "koeniz" if gdename == "liebefeid" & PLZ4 == 3097
replace gdename = "koeniz" if gdename == "liebefeld" & PLZ4 == 3097
replace gdename = "konolfingen" if gdename == "konoifingen gruenegg" & PLZ4 == 3510
replace gdename = "konolfingen" if gdename == "konolfingen gruenegg" & PLZ4 == 3510
replace gdename = "crans montana" if gdename == "corin montana" & PLZ4 == 3963
replace gdename = "sierre" if gdename == "sierra" & PLZ4 == 3960
replace gdename = "sierre" if gdename == "slerre" & PLZ4 == 3960
replace gdename = "gansingen" if gdename == "gansingen ag" & PLZ4 == 5272
replace gdename = "gansingen" if gdename == "gansirngen" & PLZ4 == 5272
replace gdename = "gansingen" if gdename == "gensingen" & PLZ4 == 5272
replace gdename = "zeihen" if gdename == "zechen" & PLZ4 == 5079
replace gdename = "zeihen" if gdename == "zelhen" & PLZ4 == 5079
replace gdename = "bremgarten" if gdename == "staffeln" & PLZ4 == 5626
replace gdename = "geltwil" if gdename == "geitwil" & PLZ4 == 5637
replace gdename = "geltwil" if gdename == "winterschwil" & PLZ4 == 5637
replace gdename = "alikon" if gdename == "alikon" & PLZ4 == 5643
replace gdename = "alikon" if gdename == "allkon" & PLZ4 == 5643
replace gdename = "saint blaise" if gdename == "st bialse" & PLZ4 == 2072
replace gdename = "saint blaise" if gdename == "st blaise" & PLZ4 == 2072
replace gdename = "fisibach" if gdename == "fis4bach" & PLZ4 == 5467
replace gdename = "fisibach" if gdename == "fisibach/ag" & PLZ4 == 5467
replace gdename = "fisibach" if gdename == "fislbach" & PLZ4 == 5467
replace gdename = "fisibach" if gdename == "flalbach" & PLZ4 == 5467
replace gdename = "fisibach" if gdename == "flslbach" & PLZ4 == 5467
replace gdename = "fisibach" if gdename == "halbach" & PLZ4 == 5467
replace gdename = "erlen" if gdename == "riedt erlen" & PLZ4 == 8586
replace gdename = "ardon" if gdename == "anion" & PLZ4 == 1957
replace gdename = "ardon" if gdename == "ardon vd" & PLZ4 == 1957
replace gdename = "ardon" if gdename == "ardon vs" & PLZ4 == 1957
replace gdename = "ardon" if gdename == "ardow" & PLZ4 == 1957
replace gdename = "courrendlin" if gdename == "choindez" & PLZ4 == 2830
replace gdename = "courrendlin" if gdename == "rebeuveller" & PLZ4 == 2832
replace gdename = "clos du doubs" if gdename == "epauviliers" & PLZ4 == 2885
replace gdename = "clos du doubs" if gdename == "epauvlllers" & PLZ4 == 2885
replace gdename = "clos du doubs" if gdename == "les rangiers" & PLZ4 == 2885
replace gdename = "worb" if gdename == "ruefenacht" & PLZ4 == 3075
replace gdename = "worb" if gdename == "ruefenacht be" & PLZ4 == 3075
replace gdename = "krauchthal" if gdename == "hub bei krauchthai" & PLZ4 == 3326
replace gdename = "krauchthal" if gdename == "hub bei krauchthal" & PLZ4 == 3326
replace gdename = "fulenbach" if gdename == "fueenbach" & PLZ4 == 4629
replace gdename = "fulenbach" if gdename == "fuienbach" & PLZ4 == 4629
replace gdename = "emmen" if gdename == "emmenbruecke" & PLZ4 == 6020
replace gdename = "muotathal" if gdename == "muotathel" & PLZ4 == 6436
replace gdename = "rietheim" if gdename == "riethelm" & PLZ4 == 5323
replace gdename = "rietheim" if gdename == "rletheim" & PLZ4 == 5323
replace gdename = "rietheim" if gdename == "rlethelm" & PLZ4 == 5323
replace gdename = "rapperswil jona" if gdename == "jona sg" & PLZ4 == 8645
replace gdename = "schuebelbach" if gdename == "siebnen" & PLZ4 == 8854
replace gdename = "glarus sued" if gdename == "haerzingen" & PLZ4 == 8775
replace gdename = "glarus sued" if gdename == "hatzingen" & PLZ4 == 8775
replace gdename = "luestisburg" if gdename == "liitisburg" & PLZ4 == 9604
replace gdename = "luestisburg" if gdename == "lirtisburg" & PLZ4 == 9604
replace gdename = "luestisburg" if gdename == "luetlsburg" & PLZ4 == 9604
replace gdename = "kradolf schoenenberg" if gdename == "schoenenberg a.thur" & PLZ4 == 9215
replace gdename = "egnach" if gdename == "neukiirch egnach" & PLZ4 == 9315
replace gdename = "egnach" if gdename == "neukirch egnach" & PLZ4 == 9315
replace gdename = "herrliberg" if gdename == "herrliberg zh" & PLZ4 == 8704
replace gdename = "herrliberg" if gdename == "herrllberg zh" & PLZ4 == 8704
replace gdename = "montanaire" if gdename == "thlerrens" & PLZ4 == 1410
replace gdename = "saint prex" if gdename == "st prei" & PLZ4 == 1162
replace gdename = "saint prex" if gdename == "st prex" & PLZ4 == 1162
replace gdename = "bex" if gdename == "fenalet sur bex" & PLZ4 == 1880
replace gdename = "chamoson" if gdename == "chamosen" & PLZ4 == 1955
replace gdename = "montreux" if gdename == "clarens" & PLZ4 == 1815
replace gdename = "montreux" if gdename == "clarens montreux" & PLZ4 == 1815
replace gdename = "koeniz" if gdename == "koeniz/kschliern" & PLZ4 == 3098
replace gdename = "koeniz" if gdename == "koeniz/schiiern" & PLZ4 == 3098
replace gdename = "koeniz" if gdename == "koeniz/schliern" & PLZ4 == 3098
replace gdename = "gurzelen" if gdename == "gurzeien" & PLZ4 == 3663
replace gdename = "oberwil bei bueren" if gdename == "oberwil bueren" & PLZ4 == 3298
replace gdename = "oberwil bei bueren" if gdename == "oberwll bueren" & PLZ4 == 3298
replace gdename = "hofstetten flueh" if gdename == "hofsteffen" & PLZ4 == 4114
replace gdename = "hofstetten flueh" if gdename == "hofsteffen so" & PLZ4 == 4114
replace gdename = "metzerlen mariastein" if gdename == "mariastein" & PLZ4 == 4116
replace gdename = "metzerlen mariastein" if gdename == "meterlen" & PLZ4 == 4116
replace gdename = "schwaderloch" if gdename == "schwaderioch" & PLZ4 == 5326
replace gdename = "wittnau" if gdename == "wi .tnau" & PLZ4 == 5064
replace gdename = "wittnau" if gdename == "wlttnau" & PLZ4 == 5064
replace gdename = "huenenberg" if gdename == "huenenberg cham" & PLZ4 == 6331
replace gdename = "personico" if gdename == "personlco" & PLZ4 == 6744
replace gdename = "prato" if gdename == "prata leventina" & PLZ4 == 6772
replace gdename = "prato" if gdename == "prato levantina" & PLZ4 == 6772
replace gdename = "klosters serneus" if gdename == "monbiei" & PLZ4 == 7250
replace gdename = "klosters serneus" if gdename == "monbiel" & PLZ4 == 7250
replace gdename = "flims" if gdename == "flims dorf" & PLZ4 == 7017
replace gdename = "langnau am albis" if gdename == "langnau a. albas" & PLZ4 == 8135
replace gdename = "langnau am albis" if gdename == "langnau a. albis" & PLZ4 == 8135
replace gdename = "kreuzlingen" if gdename == "kreuzlangen" & PLZ4 == 8280
replace gdename = "kreuzlingen" if gdename == "krezlingen" & PLZ4 == 8280
replace gdename = "rapperswil jona" if gdename == "wagen" & PLZ4 == 8646
replace gdename = "waengi" if gdename == "rosental" & PLZ4 == 9545
replace gdename = "waengi" if gdename == "rosental waengi" & PLZ4 == 9545
replace gdename = "la chaux de fonds" if gdename == "oa chaux de fonds" & PLZ4 == 2300
replace gdename = "linden" if gdename == "linden b. oberdiessbach" & PLZ4 == 3673
replace gdename = "linden" if gdename == "linden be" & PLZ4 == 3673
replace gdename = "linden" if gdename == "linden bei oberdiessbach" & PLZ4 == 3673
replace gdename = "rueti bei lyssach" if gdename == "rueti lyssach" & PLZ4 == 3421
replace gdename = "laufenburg" if gdename == "laufen ourg" & PLZ4 == 5080
replace gdename = "laufenburg" if gdename == "laufenbur9" & PLZ4 == 5080
replace gdename = "laufenburg" if gdename == "laufer..burg" & PLZ4 == 5080
replace gdename = "laufenburg" if gdename == "lautenburg" & PLZ4 == 5080
replace gdename = "mettauertal" if gdename == "hottwil e wag wernli ag, granichen del" & PLZ4 == 5277
replace gdename = "mettauertal" if gdename == "hottwll" & PLZ4 == 5277
replace gdename = "lenzburg" if gdename == "lenzburg stauffen" & PLZ4 == 5600
replace gdename = "muehlau" if gdename == "muehiau unterhuenenberg" & PLZ4 == 5642
replace gdename = "muehlau" if gdename == "muehlau unterhuenenberg" & PLZ4 == 5642
replace gdename = "torricella taverne" if gdename == "taverne torricella" & PLZ4 == 6808
replace gdename = "croglio" if gdename == "purasca" & PLZ4 == 6989
replace gdename = "scuol" if gdename == "post crusch" & PLZ4 == 7554
replace gdename = "altstaetten" if gdename == "luechingen" & PLZ4 == 9450
replace gdename = "walzenhausen" if gdename == "waizenhausen" & PLZ4 == 9428
replace gdename = "walzenhausen" if gdename == "walzen hausen" & PLZ4 == 9428
replace gdename = "jussy" if gdename == "ch. de la monaisse. 1254 jussy" & PLZ4 == 1254
replace gdename = "arth" if gdename == "oberarth/arth" & PLZ4 == 6414
replace gdename = "geneve" if gdename == "petit saconnex" & PLZ4 == 1202
replace gdename = "saint cergue" if gdename == "st cergue" & PLZ4 == 1264
replace gdename = "arlesheim" if gdename == "ariesheim" & PLZ4 == 4144
replace gdename = "arlesheim" if gdename == "arlesheirn" & PLZ4 == 4144
replace gdename = "arlesheim" if gdename == "arleshelm" & PLZ4 == 4144
replace gdename = "neuchatel" if gdename == "neuchaetel" & PLZ4 == 2000
replace gdename = "neuchatel" if gdename == "neuchatei" & PLZ4 == 2000
replace gdename = "la chaux de fonds" if gdename == "les bulles" & PLZ4 == 2300
replace gdename = "la neuveville" if gdename == "la neuveville be" & PLZ4 == 2520
replace gdename = "les enfers" if gdename == "cernlevlllers" & PLZ4 == 2363
replace gdename = "saignelegier" if gdename == "goumots" & PLZ4 == 2353
replace gdename = "belp" if gdename == "beipberg" & PLZ4 == 3124
replace gdename = "kirchdorf" if gdename == "geiterfingen" & PLZ4 == 3126
replace gdename = "kirchdorf" if gdename == "geitertingen" & PLZ4 == 3126
replace gdename = "grossaffoltern" if gdename == "suberg" & PLZ4 == 3262
replace gdename = "lyss" if gdename == "lyse" & PLZ4 == 3250
replace gdename = "sierre" if gdename == "slerre" & PLZ4 == 3960
replace gdename = "muensingen" if gdename == "muensingen be" & PLZ4 == 3110
replace gdename = "duerrenroth" if gdename == "durrenroth" & PLZ4 == 3465
replace gdename = "biglen" if gdename == "bigien" & PLZ4 == 3507
replace gdename = "biglen" if gdename == "blgien" & PLZ4 == 3507
replace gdename = "hofstetten bei brienz" if gdename == "hofstetten" & PLZ4 == 3858
replace gdename = "rohrbach" if gdename == "rohrbach b. huttwil" & PLZ4 == 4938
replace gdename = "rheinfelden" if gdename == "rhelnfelden" & PLZ4 == 4310
replace gdename = "lupfig" if gdename == "scherz ag" & PLZ4 == 5246
replace gdename = "elfingen" if gdename == "elf ingen" & PLZ4 == 5077
replace gdename = "elfingen" if gdename == "elflangen" & PLZ4 == 5077
replace gdename = "remetschwil" if gdename == "busslingen" & PLZ4 == 5453
replace gdename = "locarno" if gdename == "solduno" & PLZ4 == 6600
replace gdename = "schiers" if gdename == "schiere" & PLZ4 == 3962
replace gdename = "schiers" if gdename == "schlers" & PLZ4 == 3962
replace gdename = "koeniz" if gdename == "wabern" & PLZ4 == 3084
replace gdename = "baeretswil" if gdename == "neuthai" & PLZ4 == 8344
replace gdename = "baeretswil" if gdename == "neuthal" & PLZ4 == 8344
replace gdename = "pfaeffikon" if gdename == "pfaeffikon zh" & PLZ4 == 8330
replace gdename = "uetikon am see" if gdename == "uetikon a. see" & PLZ4 == 8707
replace gdename = "uetikon am see" if gdename == "uetlkon a. see" & PLZ4 == 8707
replace gdename = "kappel am albis" if gdename == "kappel a. albis" & PLZ4 == 8926
replace gdename = "spreitenbach" if gdename == "spreltenbach" & PLZ4 == 8957
replace gdename = "waldkirch" if gdename == "bernhardzeii" & PLZ4 == 9304
replace gdename = "waldkirch" if gdename == "bernhardzell" & PLZ4 == 9304
replace gdename = "wil sg" if gdename == "wii" & PLZ4 == 9500
replace gdename = "wil sg" if gdename == "wii sg" & PLZ4 == 9500
replace gdename = "wil sg" if gdename == "wit sg" & PLZ4 == 9500
replace gdename = "ebnat kappel" if gdename == "wintersberg" & PLZ4 == 9642
replace gdename = "ebnat kappel" if gdename == "wintersberg sg" & PLZ4 == 9642
replace gdename = "bichelsee balterswil" if gdename == "balterswil tg" & PLZ4 == 8362
replace gdename = "kuesnacht" if gdename == "forch zh" & PLZ4 == 8127
replace gdename = "vesenaz" if gdename == "vesenaz" & PLZ4 == 1245
replace gdename = "grancy" if gdename == "i 1 grancy" & PLZ4 == 1117
replace gdename = "martigny" if gdename == "rue de l hopital. 1920 martigny" & PLZ4 == 1920
replace gdename = "corseaux" if gdename == "rte. de lavaux. 1802 corseaux" & PLZ4 == 1802
replace gdename = "schoenenbuch" if gdename == "schoenenbuch bl" & PLZ4 == 4124
replace gdename = "wohlen bei bern" if gdename == "lnnerberg" & PLZ4 == 3044
replace gdename = "wohlen bei bern" if gdename == "uettligen" & PLZ4 == 3043
replace gdename = "bernex" if gdename == "sezenove" & PLZ4 == 1233
replace gdename = "grandson" if gdename == "3randson corcelettes" & PLZ4 == 1422
replace gdename = "grandson" if gdename == "grandson corceiettes" & PLZ4 == 1422
replace gdename = "chardonne" if gdename == "le mont pellerin" & PLZ4 == 1803
replace gdename = "oberwil" if gdename == "oberwal" & PLZ4 == 4104
replace gdename = "oberwil" if gdename == "oberwii" & PLZ4 == 4104
replace gdename = "valbirse" if gdename == "malleray bevilard" & PLZ4 == 2735
replace gdename = "bellmund" if gdename == "4ellmund" & PLZ4 == 2564
replace gdename = "bellmund" if gdename == "belimund" & PLZ4 == 2564
replace gdename = "lajoux" if gdename == "fornet dessus" & PLZ4 == 2718
replace gdename = "lajoux" if gdename == "fornet dessus/lajoux" & PLZ4 == 2718
replace gdename = "saignelegier" if gdename == "pommerats" & PLZ4 == 2353
replace gdename = "valbirse" if gdename == "poutenet" & PLZ4 == 2733
replace gdename = "davos" if gdename == "clavadel" & PLZ4 == 7272
replace gdename = "grone" if gdename == "groene" & PLZ4 == 3979
replace gdename = "bern" if gdename == "berne" & PLZ4 == 3000
replace gdename = "bolligen" if gdename == "roerswil bolligen" & PLZ4 == 3065
replace gdename = "koeniz" if gdename == "wabern" & PLZ4 == 3084
replace gdename = "uetendorf" if gdename == "iietendorf" & PLZ4 == 3661
replace gdename = "uetendorf" if gdename == "uetendort" & PLZ4 == 3661
replace gdename = "wuennewil flamatt" if gdename == "wuennewil" & PLZ4 == 3184
replace gdename = "seeberg" if gdename == "grasswil" & PLZ4 == 3365
replace gdename = "fuellinsdorf" if gdename == "fuelllnsdorf" & PLZ4 == 4414
replace gdename = "fuellinsdorf" if gdename == "fullinsdorf" & PLZ4 == 4414
replace gdename = "mettauertal" if gdename == "wi ag" & PLZ4 == 5276
replace gdename = "mettauertal" if gdename == "wii ag" & PLZ4 == 5276
replace gdename = "itingen" if gdename == "itigen" & PLZ4 == 4452
replace gdename = "itingen" if gdename == "ltingen" & PLZ4 == 4452
replace gdename = "trimbach" if gdename == "ti imbach" & PLZ4 == 4632
replace gdename = "trimbach" if gdename == "trlmbach" & PLZ4 == 4632
replace gdename = "hornussen" if gdename == "homussen" & PLZ4 == 5075
replace gdename = "oberhof" if gdename == "oberhot" & PLZ4 == 5062
replace gdename = "oberhofen" if gdename == "oherhofen" & PLZ4 == 5273
replace gdename = "eiken" if gdename == "eikerl" & PLZ4 == 5074
replace gdename = "eiken" if gdename == "elken" & PLZ4 == 5074
replace gdename = "uerkheim" if gdename == "uekheim" & PLZ4 == 4813
replace gdename = "uerkheim" if gdename == "uerkhelm" & PLZ4 == 4813
replace gdename = "muri (ag)" if gdename == "murl" & PLZ4 == 5630
replace gdename = "muri (ag)" if gdename == "murl (ag)" & PLZ4 == 5630
replace gdename = "canobbio" if gdename == "canobblo" & PLZ4 == 6952
replace gdename = "locarno" if gdename == "locarno/muralto" & PLZ4 == 6600
replace gdename = "wald (zh)" if gdename == "huebeli" & PLZ4 == 8636
replace gdename = "wald (zh)" if gdename == "huebell" & PLZ4 == 8636
replace gdename = "hasle" if gdename == "heiligkreuz lu" & PLZ4 == 6166
replace gdename = "kuessnacht" if gdename == "merlischachen" & PLZ4 == 6402
replace gdename = "flueelen" if gdename == "floeleri ur" & PLZ4 == 6454
replace gdename = "flueelen" if gdename == "flueeien ur" & PLZ4 == 6454
replace gdename = "lichtensteig" if gdename == "lichtenstelg" & PLZ4 == 9620
replace gdename = "lichtensteig" if gdename == "liechtenstein" & PLZ4 == 9620
replace gdename = "kuesnacht" if gdename == "forch" & PLZ4 == 8127
replace gdename = "schuepfheim" if gdename == "schoepfhelm" & PLZ4 == 6170
replace gdename = "schuepfheim" if gdename == "schuepfhelm" & PLZ4 == 6170
replace gdename = "lindau" if gdename == "kemptal lindau" & PLZ4 == 8310
replace gdename = "kemmental" if gdename == "neuwilen/lippoldswilen" & PLZ4 == 8566
replace gdename = "kemmental" if gdename == "neuwilen/llppoldswilen" & PLZ4 == 8566
replace gdename = "faellanden" if gdename == "benglen faellanden" & PLZ4 == 8121
replace gdename = "gommiswald" if gdename == "mellen" & PLZ4 == 8739
replace gdename = "gommiswald" if gdename == "rieden gl" & PLZ4 == 8739
replace gdename = "einsiedeln" if gdename == "einsiedein" & PLZ4 == 8840
replace gdename = "einsiedeln" if gdename == "lamons" & PLZ4 == 8840
replace gdename = "einsiedeln" if gdename == "bennau" & PLZ4 == 8840
replace gdename = "zuerich" if gdename == "zurich" & PLZ4 == 8038
replace gdename = "massagno" if gdename == "lugano/massagno" & PLZ4 == 6900
replace gdename = "paradiso" if gdename == "paradlso" & PLZ4 == 6900
replace gdename = "kappel am albis" if gdename == "kappel am albfis" & PLZ4 == 8926
replace gdename = "kappel am albis" if gdename == "uerzllkon" & PLZ4 == 8926
replace gdename = "croglio" if gdename == "crogllo" & PLZ4 == 6981
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolfstetten" & PLZ4 == 8964
replace gdename = "flums" if gdename == "flumserberg" & PLZ4 == 8890
replace gdename = "urnaesch" if gdename == "zuerchersmuehle" & PLZ4 == 9107
replace gdename = "diepoldsau" if gdename == "dlepoldsau" & PLZ4 == 9444
replace gdename = "brigels" if gdename == "brell/brigels" & PLZ4 == 7165
replace gdename = "corsier" if gdename == "ch. des usses. 1246 corsier" & PLZ4 == 1246
replace gdename = "perly certoux" if gdename == "ch du village. 125e perly" & PLZ4 == 1258
replace gdename = "kuesnacht" if gdename == "forch" & PLZ4 == 8127
replace gdename = "geneve" if gdename == "r.de lausanne41201 geneve" & PLZ4 == 1201
replace gdename = "teufen" if gdename == "hiederteufen" & PLZ4 == 9053
replace gdename = "lugano" if gdename == "viganello ti" & PLZ4 == 6962
replace gdename = "vernier" if gdename == "rte de peney. 1214 vernier" & PLZ4 == 1214
replace gdename = "faellanden" if gdename == "benglen" & PLZ4 == 8121
replace gdename = "vesenaz" if gdename == "vesenaz" & PLZ4 == 1245
replace gdename = "burtigny" if gdename == "1261 burtigriy vd" & PLZ4 == 1261
replace gdename = "bardonnex" if gdename == "croix de rozon" & PLZ4 == 1257
replace gdename = "opfikon" if gdename == "glattbrugg" & PLZ4 == 8152
replace gdename = "villars sur glane" if gdename == "l/illars sur glane" & PLZ4 == 1752
replace gdename = "lavey morcles" if gdename == "lavey village" & PLZ4 == 1892
replace gdename = "la neuveville" if gdename == "chavannes (be)" & PLZ4 == 2520
replace gdename = "le grand saconnex" if gdename == "grand saconnex" & PLZ4 == 1218
replace gdename = "muri bei bern" if gdename == "muri b. bern" & PLZ4 == 3074
replace gdename = "urtenen schoenbuehl" if gdename == "urtenen schoenbuehl" & PLZ4 == 3322
replace gdename = "sisseln" if gdename == "sisseln ag" & PLZ4 == 4334
replace gdename = "elgg" if gdename == "wenzikon" & PLZ4 == 8354
replace gdename = "maerstetten" if gdename == "maerstetten dorf" & PLZ4 == 8560
replace gdename = "weinfelden" if gdename == "welnfelden" & PLZ4 == 8570
replace gdename = "zollikon" if gdename == "zotlikssn" & PLZ4 == 8702
replace gdename = "herrliberg" if gdename == "herrtlberg" & PLZ4 == 8704
replace gdename = "klosters serneus" if gdename == "klosters" & PLZ4 == 7250
replace gdename = "oberengstringen" if gdename == "l.)berengstringen" & PLZ4 == 8102
replace gdename = "rueti" if gdename == "tann rueti" & PLZ4 == 8632
replace gdename = "esslingen" if gdename == "esslingen" & PLZ4 == 8133
replace gdename = "muellheim" if gdename == "muellheim dorf" & PLZ4 == 8555
replace gdename = "murgenthal" if gdename == "riken/ murgenthal" & PLZ4 == 4853
replace gdename = "yverdon les bains" if gdename == "yverdon" & PLZ4 == 1400
replace gdename = "vesenaz" if gdename == "vesenaz/" & PLZ4 == 1245
replace gdename = "crissier" if gdename == "crissler" & PLZ4 == 1023
replace gdename = "saint sulpice (vd)" if gdename == "st.sulpice" & PLZ4 == 1025
replace gdename = "lancy" if gdename == "grand lancy" & PLZ4 == 1212
replace gdename = "montreux" if gdename == "chamby" & PLZ4 == 1832
replace gdename = "montreux" if gdename == "le chaetelard (vd)" & PLZ4 == 1820
replace gdename = "montreux" if gdename == "clarens" & PLZ4 == 1815
replace gdename = "la baroche" if gdename == "les rangiers" & PLZ4 == 2954
replace gdename = "sion" if gdename == "slon" & PLZ4 == 1950
replace gdename = "uttigen" if gdename == "uttlgen" & PLZ4 == 3628
replace gdename = "la brevine" if gdename == "le brouillet" & PLZ4 == 2406
replace gdename = "la chaux de fonds" if gdename == "4a chaux de fonds" & PLZ4 == 2300
replace gdename = "studen" if gdename == "studen be" & PLZ4 == 2557
replace gdename = "himmelried" if gdename == "hammelried" & PLZ4 == 4204
replace gdename = "bienne" if gdename == "biel/bienne" & PLZ4 == 2500
replace gdename = "moerigen" if gdename == "moeriger" & PLZ4 == 2572
replace gdename = "studen" if gdename == "studen b. bruegg" & PLZ4 == 2557
replace gdename = "nidau" if gdename == "nldau" & PLZ4 == 2560
replace gdename = "evilard" if gdename == "evllard" & PLZ4 == 2533
replace gdename = "selzach" if gdename == "setzach" & PLZ4 == 2545
replace gdename = "landquart" if gdename == "landquart fabriken" & PLZ4 == 7302
replace gdename = "rebeuvelier" if gdename == "rebeuveller" & PLZ4 == 2832
replace gdename = "lauterbrunnen" if gdename == "murren" & PLZ4 == 3825
replace gdename = "viege" if gdename == "viege" & PLZ4 == 3930
replace gdename = "crans montana" if gdename == "crans sur sierre" & PLZ4 == 3963
replace gdename = "saint blaise" if gdename == "st. blaise" & PLZ4 == 2072
replace gdename = "bolligen" if gdename == "ferenberg bolligen" & PLZ4 == 3065
replace gdename = "cham" if gdename == "hagendorn" & PLZ4 == 6332
replace gdename = "uitikon" if gdename == "waldegg b. uitikon" & PLZ4 == 8142
replace gdename = "boesingen" if gdename == "besingen" & PLZ4 == 3178
replace gdename = "burgistein" if gdename == "burgisteln" & PLZ4 == 3664
replace gdename = "aarberg" if gdename == "lobsagen/aarberg" & PLZ4 == 3270
replace gdename = "zuzwil (be)" if gdename == "zuzwll be" & PLZ4 == 3303
replace gdename = "urtenen schoenbuehl" if gdename == "schoenbuehl urtenen" & PLZ4 == 3322
replace gdename = "elgg" if gdename == "elgq" & PLZ4 == 8353
replace gdename = "aeschi" if gdename == "stelnhof" & PLZ4 == 4556
replace gdename = "niederoenz" if gdename == "nlederoenz" & PLZ4 == 3362
replace gdename = "liebefeld" if gdename == "liebefeld" & PLZ4 == 3097
replace gdename = "hettlingen" if gdename == "hettlingerr" & PLZ4 == 8442
replace gdename = "zell" if gdename == "kollbrunn" & PLZ4 == 8483
replace gdename = "sankt stephan" if gdename == "matten" & PLZ4 == 3773
replace gdename = "heiligenschwendi" if gdename == "heiligenschwendl" & PLZ4 == 3625
replace gdename = "gachnang" if gdename == "lslikon" & PLZ4 == 8546
replace gdename = "sigriswil" if gdename == "gunten" & PLZ4 == 3654
replace gdename = "diemtigen" if gdename == "oey" & PLZ4 == 3754
replace gdename = "sommeri" if gdename == "niedersommeri" & PLZ4 == 8580
replace gdename = "rapperswil jona" if gdename == "rapperswii sg" & PLZ4 == 8645
replace gdename = "lauenen" if gdename == "lauenen bei gstaad" & PLZ4 == 3782
replace gdename = "saanen" if gdename == "gstaad" & PLZ4 == 3780
replace gdename = "saanen" if gdename == "schoenried" & PLZ4 == 3778
replace gdename = "thun" if gdename == "thun allmendingen" & PLZ4 == 3608
replace gdename = "fahrni" if gdename == "fahrni b. thun" & PLZ4 == 3617
replace gdename = "muri bei bern" if gdename == "muri b. bern" & PLZ4 == 3074
replace gdename = "bern" if gdename == "bem" & PLZ4 == 3000
replace gdename = "allschwil" if gdename == "neu allschwil" & PLZ4 == 4123
replace gdename = "meltingen" if gdename == "neffingen" & PLZ4 == 4233
replace gdename = "mettauertal" if gdename == "meltau" & PLZ4 == 5274
replace gdename = "ennetbuergen" if gdename == "buergenstock" & PLZ4 == 6373
replace gdename = "kaiseraugst" if gdename == "kalseraugst" & PLZ4 == 4303
replace gdename = "hoelstein" if gdename == "holstein" & PLZ4 == 4434
replace gdename = "laeufelfingen" if gdename == "laufelfingen" & PLZ4 == 4448
replace gdename = "messen" if gdename == "oberramsem" & PLZ4 == 4588
replace gdename = "trimbach" if gdename == "trlmbach" & PLZ4 == 4632
replace gdename = "hauenstein ifenthal" if gdename == "hauensteln" & PLZ4 == 4633
replace gdename = "kappel" if gdename == "kappet so" & PLZ4 == 4616
replace gdename = "oberbuchsiten" if gdename == "oberbuchstten" & PLZ4 == 4625
replace gdename = "kuettigen" if gdename == "flombach" & PLZ4 == 5022
replace gdename = "alpnach" if gdename == "alpnach dorf" & PLZ4 == 6055
replace gdename = "eglisau" if gdename == "eglã¬sau" & PLZ4 == 8193
replace gdename = "effingen" if gdename == "effinaen" & PLZ4 == 5078
replace gdename = "altenrhein" if gdename == "altenrhein" & PLZ4 == 9423
replace gdename = "biasca" if gdename == "blasco" & PLZ4 == 6710
replace gdename = "dottikon" if gdename == "dottlkon" & PLZ4 == 5605
replace gdename = "meisterschwanden" if gdename == "melsterschwanden" & PLZ4 == 5616
replace gdename = "freienbach" if gdename == "wilen freienbach sz" & PLZ4 == 8807
replace gdename = "affeltrangen" if gdename == "buch bel maerwil" & PLZ4 == 9556
replace gdename = "arni" if gdename == "arni islisberg" & PLZ4 == 8905
replace gdename = "collina d'oro" if gdename == "gentiilno" & PLZ4 == 6925
replace gdename = "wettingen" if gdename == "wettlngen" & PLZ4 == 5430
replace gdename = "emmen" if gdename == "emmenbruecke" & PLZ4 == 6020
replace gdename = "lungern" if gdename == "lungem" & PLZ4 == 6078
replace gdename = "kuesnacht" if gdename == "forch" & PLZ4 == 8127
replace gdename = "grossdietwil" if gdename == "grossdletwil" & PLZ4 == 6146
replace gdename = "windisch" if gdename == "windlsch" & PLZ4 == 5210
replace gdename = "illnau effretikon" if gdename == "illnau" & PLZ4 == 8308
replace gdename = "arth" if gdename == "rigl" & PLZ4 == 6410
replace gdename = "birmensdorf" if gdename == "birmenedorf" & PLZ4 == 8903
replace gdename = "arth" if gdename == "goldau" & PLZ4 == 6410
replace gdename = "ibach" if gdename == "ibach" & PLZ4 == 6438
replace gdename = "arth" if gdename == "goldau" & PLZ4 == 6410
replace gdename = "tenero contra" if gdename == "tenero" & PLZ4 == 6598
replace gdename = "moenchaltdorf" if gdename == "moenchaltdorf" & PLZ4 == 8617
replace gdename = "rueti" if gdename == "tann rueti" & PLZ4 == 8632
replace gdename = "herrliberg" if gdename == "herrllberg" & PLZ4 == 8704
replace gdename = "meilen" if gdename == "mellen" & PLZ4 == 8706
replace gdename = "faido" if gdename == "faldo" & PLZ4 == 6760
replace gdename = "orselina" if gdename == "orseline" & PLZ4 == 6644
replace gdename = "lugano" if gdename == "lugano cassarate" & PLZ4 == 6976
replace gdename = "curio" if gdename == "cuvio" & PLZ4 == 6986
replace gdename = "clugin" if gdename == "ciugin" & PLZ4 == 7442
replace gdename = "davos" if gdename == "devos dorf" & PLZ4 == 7260
replace gdename = "disentis/muster" if gdename == "segnas" & PLZ4 == 7186
replace gdename = "celerina/schlarigna" if gdename == "celerina" & PLZ4 == 7505
replace gdename = "celerina/schlarigna" if gdename == "celerina/ schlarigna" & PLZ4 == 7505
replace gdename = "val muestair" if gdename == "santa marla im muenstertal" & PLZ4 == 7536
replace gdename = "bregaglia" if gdename == "promontogno" & PLZ4 == 7606
replace gdename = "zuerich" if gdename == "zurich" & PLZ4 == 8000
replace gdename = "muri bei bern" if gdename == "guemligen" & PLZ4 == 3073
replace gdename = "willisau" if gdename == "rohrmatt wllllsau land" & PLZ4 == 6132
replace gdename = "steinmaur" if gdename == "niedersteinmaur" & PLZ4 == 8162
replace gdename = "reiden" if gdename == "reidermoos" & PLZ4 == 6260
replace gdename = "thalwil" if gdename == "thatwit" & PLZ4 == 8800
replace gdename = "nuerensdorf" if gdename == "birchwll/ nuerensdorf" & PLZ4 == 8309
replace gdename = "glaris nord" if gdename == "obstaiden" & PLZ4 == 8758
replace gdename = "bichelsee balterswil" if gdename == "balterswll" & PLZ4 == 8362
replace gdename = "erlen" if gdename == "riedt bei erlen" & PLZ4 == 8586
replace gdename = "engelberg" if gdename == "engeiberg" & PLZ4 == 6390
replace gdename = "rekingen" if gdename == "reklngen" & PLZ4 == 5332
replace gdename = "ingenbohl" if gdename == "brunnen" & PLZ4 == 6440
replace gdename = "gossau (zh)" if gdename == "gruet" & PLZ4 == 8624
replace gdename = "orselina" if gdename == "oreellna" & PLZ4 == 6644
replace gdename = "minusio" if gdename == "minuslo" & PLZ4 == 6648
replace gdename = "schlieren" if gdename == "schlieren zh" & PLZ4 == 8952
replace gdename = "terre di pedemonte" if gdename == "versclo" & PLZ4 == 6653
replace gdename = "oetwil am see" if gdename == "oetwil a/see" & PLZ4 == 8618
replace gdename = "faido" if gdename == "faldo" & PLZ4 == 6760
replace gdename = "moenchaltdorf" if gdename == "moenchaitorf" & PLZ4 == 8617
replace gdename = "wetzikon" if gdename == "wetzlkon" & PLZ4 == 8620
replace gdename = "balerna" if gdename == "belerna" & PLZ4 == 6828
replace gdename = "filzbach" if gdename == "fllzbach" & PLZ4 == 8757
replace gdename = "langnau am albis" if gdename == "langnau" & PLZ4 == 8135
replace gdename = "oberwil lieli" if gdename == "lieii ag" & PLZ4 == 8966
replace gdename = "lugano" if gdename == "ruvigliana" & PLZ4 == 6977
replace gdename = "berikon" if gdename == "mutschellen" & PLZ4 == 8965
replace gdename = "teufen" if gdename == "lustmuehle" & PLZ4 == 9062
replace gdename = "uzwil" if gdename == "niederuzwil" & PLZ4 == 9244
replace gdename = "flawil" if gdename == "flawll" & PLZ4 == 9230
replace gdename = "nuerensdorf" if gdename == "oberwii b. nuerensdorf" & PLZ4 == 8309
replace gdename = "muolen" if gdename == "mouien" & PLZ4 == 9313
replace gdename = "egnach" if gdename == "steinebrunn egnach" & PLZ4 == 9314
replace gdename = "gaiserwald" if gdename == "engel burg" & PLZ4 == 9032
replace gdename = "teufen" if gdename == "nlederteufen" & PLZ4 == 9052
replace gdename = "heiden" if gdename == "helden" & PLZ4 == 9410
replace gdename = "oberhelfenschwil" if gdename == "oberheifenschwil" & PLZ4 == 9621
replace gdename = "altstaetten" if gdename == "luechingen" & PLZ4 == 9450
replace gdename = "wartau" if gdename == "truebbach" & PLZ4 == 9477
replace gdename = "wattwil" if gdename == "wattwit" & PLZ4 == 9630
replace gdename = "zuzwil (sg)" if gdename == "zuzwll sg" & PLZ4 == 9524
replace gdename = "waengi" if gdename == "waengl" & PLZ4 == 9545
replace gdename = "affeltrangen" if gdename == "maerwil" & PLZ4 == 9562
replace gdename = "luestisburg" if gdename == "luetleburg" & PLZ4 == 9604
replace gdename = "buetschwil ganterschwil" if gdename == "buetechwll" & PLZ4 == 9606
replace gdename = "wattwil" if gdename == "wattwll" & PLZ4 == 9630
replace gdename = "buehler" if gdename == "buehlen" & PLZ4 == 9055
replace gdename = "lancy" if gdename == "grand lancy" & PLZ4 == 1212
replace gdename = "arzier le muids" if gdename == "le muids" & PLZ4 == 1273
replace gdename = "thonex" if gdename == "hã³nex" & PLZ4 == 1226
replace gdename = "bernex" if gdename == "bernez" & PLZ4 == 1233
replace gdename = "illnau effretikon" if gdename == "effretikon" & PLZ4 == 8307
replace gdename = "la chaux de fonds" if gdename == "chaux de fonds" & PLZ4 == 2300
replace gdename = "walterswil" if gdename == "schmidigen" & PLZ4 == 3464
replace gdename = "basel" if gdename == "basei" & PLZ4 == 4000
replace gdename = "sissach" if gdename == "slssach" & PLZ4 == 4450
replace gdename = "kriens" if gdename == "riens" & PLZ4 == 6010
replace gdename = "torricella taverne" if gdename == "taverne" & PLZ4 == 6808
replace gdename = "eschlikon" if gdename == "eschlikon wallenwil" & PLZ4 == 8360
replace gdename = "landquart" if gdename == "gis" & PLZ4 == 7206
replace gdename = "zuerich" if gdename == "zurich" & PLZ4 == 8000
replace gdename = "neftenbach" if gdename == "aesch/neftenbach" & PLZ4 == 8413
replace gdename = "wislikofen" if gdename == "mellsdorf" & PLZ4 == 5463
replace gdename = "gachnang" if gdename == "oberopferstofen" & PLZ4 == 8547
replace gdename = "uetikon am see" if gdename == "uetikon a. see" & PLZ4 == 8707
replace gdename = "kuesnacht" if gdename == "forch/aesch" & PLZ4 == 8127
replace gdename = "wil (sg)" if gdename == "wii (sg)" & PLZ4 == 9500
replace gdename = "oberhelfenschwil" if gdename == "berhelfenschwi i" & PLZ4 == 9621
replace gdename = "roggwil" if gdename == "freihof" & PLZ4 == 9306
replace gdename = "lungern" if gdename == "buergten" & PLZ4 == 6078
replace gdename = "chene bougeries" if gdename == "conches" & PLZ4 == 1231
replace gdename = "merenschwand" if gdename == "hagnau" & PLZ4 == 5634
replace gdename = "vechigen" if gdename == "boll" & PLZ4 == 3067
replace gdename = "geneve" if gdename == "av. mervelet" & PLZ4 == 1209
replace gdename = "wila" if gdename == "wile" & PLZ4 == 8492
replace gdename = "baar" if gdename == "saar" & PLZ4 == 6340
replace gdename = "guettingen" if gdename == "goettingen" & PLZ4 == 8594


// futher corrections (without PLZ4)
gen n_gde = 1
replace n_gde = 3 if gdename == "affoltern"
replace n_gde = 3 if gdename == "hausen"
replace n_gde = 2 if gdename == "oetwil"
replace n_gde = 3 if gdename == "ried"
replace n_gde = 2 if gdename == "sils"
replace n_gde = 2 if gdename == "st. sulpice"
replace n_gde = 2 if gdename == "st sulpice"
replace n_gde = 2 if gdename == "vuisternens"
replace n_gde = 2 if gdename == "wasen"
replace n_gde = 2 if gdename == "munchwilen"
replace n_gde = 2 if gdename == "kilchherg"
replace n_gde = 2 if gdename == "le coudre"
replace n_gde = 2 if gdename == "romanel"
replace n_gde = 2 if gdename == "killchberg"
replace n_gde = 8 if gdename == "chavannes"
replace n_gde = 6 if gdename == "wilen"
replace n_gde = 3 if gdename == "langnau"
replace n_gde = 2 if gdename == "cheseaux"
replace n_gde = 2 if gdename == "morbio"
replace n_gde = 2 if gdename == "oulens"
replace n_gde = 2 if gdename == "pfaelfikon"
replace n_gde = 2 if gdename == "roethenbach"
replace n_gde = 2 if gdename == "st niklausen"
replace n_gde = 2 if gdename == "st saphorin"
replace n_gde = 2 if gdename == "st. niklausen"
replace n_gde = 2 if gdename == "taverne"
replace n_gde = 2 if gdename == "uttwil romanshorn"
replace n_gde = 2 if gdename == "wasen"
replace n_gde = 2 if gdename == "zermatt et sion"
replace n_gde = 6 if gdename == "corcelles"
replace n_gde = 3 if gdename == "langnau"
replace n_gde = 3 if gdename == "oberhofen"
replace n_gde = 3 if gdename == "burg"
replace n_gde = 2 if gdename == "erlenbach"
replace n_gde = 2 if gdename == "lussy"
replace n_gde = 2 if gdename == "schwanden"
replace n_gde = 2 if gdename == "roveredo"
replace n_gde = 2 if gdename == "hergiswil"
replace n_gde = 2 if gdename == "hasle"
replace n_gde = 2 if gdename == "corsier"
replace n_gde = 2 if gdename == "wiler"
replace n_gde = 2 if gdename == "bueren"
replace n_gde = 2 if gdename == "bellerive"
replace n_gde = 2 if gdename == "wohlen"
replace n_gde = 2 if gdename == "muri"
replace n_gde = 4 if gdename == "hofstetten"
replace n_gde = 4 if gdename == "buchs"
replace n_gde = 3 if gdename == "busswil"
replace n_gde = 3 if gdename == "beinwil"
replace n_gde = 2 if gdename == "gossau"
replace n_gde = 2 if gdename == "oberried"
replace n_gde = 2 if gdename == "stetten"
replace n_gde = 2 if gdename == "gossau"
replace n_gde = 2 if gdename == "seedorf"

replace gdename = "wisen (so)" if gdename == "wiesen"
replace gdename = "biel/bienne" if gdename == "bienne"
replace gdename = "le locle" if gdename == "le lode"
replace gdename = "muri bei bern" if gdename == "muri/bern"
replace gdename = "davos" if gdename == "davos platz"
replace gdename = "sainte croix" if gdename == "ste croix"
replace gdename = "saint imier" if gdename == "st imier"
replace gdename = "le grand saconnex" if gdename == "le grand-saconnex"
replace gdename = "le chenit" if gdename == "le brassus"
replace gdename = "le chenit" if gdename == "le sentier"
replace gdename = "berneck" if gdename == "heerbrugg"
replace gdename = "vaz/obervaz" if gdename == "lenzerheide"
replace gdename = "pregny chambesy" if gdename == "pregny"
replace gdename = "kuessnacht (sz)" if gdename == "kuessnacht"
replace gdename = "kuesnacht (zh)" if gdename == "kiisnacht"
replace gdename = "evilard" if gdename == "leubringen"
replace gdename = "saint blaise" if gdename == "saint-blaise"
replace gdename = "saint maurice" if gdename == "st maurice"
replace gdename = "lancy" if gdename == "petit lancy"
replace gdename = "rueschlikon" if gdename == "ruschlikon"
replace gdename = "saint imier" if gdename == "st. imier"
replace gdename = "kirchberg (sg)" if gdename == "bazenheid"
replace gdename = "muenchenstein" if gdename == "munchenstein"
replace gdename = "saint prex" if gdename == "saint-prex"
replace gdename = "teufen" if gdename == "niederteufen"
replace gdename = "crans montana" if gdename == "crans-montana"
replace gdename = "collonge bellerive" if gdename == "collonges bellerive"
replace gdename = "neuhausen am rheinfall" if gdename == "neuhausen"
replace gdename = "zollikon" if gdename == "zollikerberg"
replace gdename = "st. stephan" if gdename == "sankt stephan"
replace gdename = "vandoeuvres" if gdename == "vandceuvres"
replace gdename = "saint aubin (fr)" if gdename == "st aubin"
replace gdename = "pregny chambesy" if gdename == "chambesy"
replace gdename = "waedenswil" if gdename == "wadenswil"
replace gdename = "le mont sur lausanne" if gdename == "mont/lausanne"
replace gdename = "davos" if gdename == "davos dorf"
replace gdename = "sarnen" if gdename == "samen"
replace gdename = "lauterbrunnen" if gdename == "wengen"
replace gdename = "lutry" if gdename == "la conversion"
replace gdename = "la tour de peilz" if gdename == "la tour de peitz"
replace gdename = "hilterfingen" if gdename == "huenibach"
replace gdename = "beinwil am see" if gdename == "beinwil a. see"
replace gdename = "koeniz" if gdename == "liebefeld koeniz"
replace gdename = "tujetsch" if gdename == "sedrun"
replace gdename = "collonge bellerive" if gdename == "vesenaz"
replace gdename = "doettingen" if gdename == "dottingen"
replace gdename = "wuennewil flamatt" if gdename == "flamatt"
replace gdename = "fully" if gdename == "frilly"
replace gdename = "vandoeuvres" if gdename == "vandåuvres"
replace gdename = "disentis/muster" if gdename == "disentis"
replace gdename = "muensingen" if gdename == "munsingen"
replace gdename = "littau" if gdename == "reussbuehl"
replace gdename = "chatel saint denis" if gdename == "chatel st denis"
replace gdename = "rueschlikon" if gdename == "riischlikon"
replace gdename = "la chaux de fonds" if gdename == "la-chaux-de-fonds"
replace gdename = "horw" if gdename == "kastanienbaum"
replace gdename = "la neuveville" if gdename == "neuveville"
replace gdename = "bern" if gdename == "bern buemplitz"
replace gdename = "arbedo castione" if gdename == "arbedo"
replace gdename = "muemliswil ramiswil" if gdename == "muemliswil"
replace gdename = "tramelan" if gdename == "tramelan dessus"
replace gdename = "corcelles cormondreche" if gdename == "cormondreche"
replace gdename = "samedan" if gdename == "samaden"
replace gdename = "staefa" if gdename == "stafa"
replace gdename = "staefa" if gdename == "uerikon"
replace gdename = "freienbach" if gdename == "baech"
replace gdename = "le landeron" if gdename == "landeron"
replace gdename = "le locle" if gdename == "le lock"
replace gdename = "crans montana" if gdename == "montana vermala"
replace gdename = "saint maurice" if gdename == "st. maurice"
replace gdename = "granges pres marnand" if gdename == "granges marnand"
replace gdename = "murten" if gdename == "morat"
replace gdename = "ittigen" if gdename == "worblaufen"
replace gdename = "illnau effretikon" if gdename == "illnau-effretikon"
replace gdename = "feldbrunnen st. niklaus" if gdename == "feldbrunnen"
replace gdename = "langnau im emmental" if gdename == "langnau i. e."
replace gdename = "feusisberg" if gdename == "schindellegi"
replace gdename = "ruederswil" if gdename == "zollbrueck"
replace gdename = "delemont" if gdename == "delsberg"
replace gdename = "reichenbach im kandertal" if gdename == "reichenbach"
replace gdename = "altstaetten" if gdename == "alstaetten"
replace gdename = "muri bei bern" if gdename == "guemligen muri"
replace gdename = "sainte croix" if gdename == "ste-croix"
replace gdename = "marin epagnier" if gdename == "marin"
replace gdename = "vilters wangs" if gdename == "wangs"
replace gdename = "crans montana" if gdename == "crans/chermignon"
replace gdename = "diessbach bei bueren" if gdename == "diessbach"
replace gdename = "kuessnacht (zh)" if gdename == "kuesnacht/zuerich"
replace gdename = "collonge bellerive" if gdename == "la capite"
replace gdename = "l abbaye" if gdename == "le pont"
replace gdename = "saint legier la chiesaz" if gdename == "st legier"
replace gdename = "sainte croix" if gdename == "ste-croix"
replace gdename = "bern" if gdename == "bern buempliz"
replace gdename = "bussigny" if gdename == "bussigny s/morges"
replace gdename = "chezard saint martin" if gdename == "chezard"
replace gdename = "nuglar st. pantaleon" if gdename == "nuglar"
replace gdename = "le chenit" if gdename == "orient"
replace gdename = "villars sur glane" if gdename == "villars/glane"
replace gdename = "ollon" if gdename == "villars/ollon"
replace gdename = "muri bei bern" if gdename == "gumligen"
replace gdename = "muenchenstein" if gdename == "neuewelt muenchenstein"
replace gdename = "salgesch" if gdename == "salquenen"
replace gdename = "walenstadt" if gdename == "wallenstadt"
replace gdename = "crans montana" if gdename == "crans/sierre"
replace gdename = "zuerich" if gdename == "dr., zurich"
replace gdename = "kuesnacht (zh)" if gdename == "goldbach"
replace gdename = "montreux" if gdename == "montreux clarens"
replace gdename = "muri bei bern" if gdename == "muri bern"
replace gdename = "geneve" if gdename == "genf"
replace gdename = "schwyz" if gdename == "ibach"
replace gdename = "ittigen" if gdename == "papiermuehle"
replace gdename = "saint blaise" if gdename == "saint-blaise"
replace gdename = "saint ursanne" if gdename == "st ursanne"
replace gdename = "belmont sur lausanne" if gdename == "belmont/lausanne"
replace gdename = "lugano" if gdename == "cassarate"
replace gdename = "davos" if gdename == "davos glaris"
replace gdename = "la tour de peilz" if gdename == "la tour de peilt"
replace gdename = "le mont sur lausanne" if gdename == "mont s/lausanne"
replace gdename = "crans montana" if gdename == "montana/vermala"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolfstetten-friedlisberg"
replace gdename = "saint leonard" if gdename == "st-leonard"
replace gdename = "muenchwilen (tg)" if gdename == "st. margarethen"
replace gdename = "saint legier la chiesaz" if gdename == "st legier la chiesaz"
replace gdename = "vandoeuvres" if gdename == "vandteuvres"
replace gdename = "cossonay" if gdename == "cossonay gare"
replace gdename = "davos" if gdename == "davos platz"
replace gdename = "gruyeres" if gdename == "epagny"
replace gdename = "hospental" if gdename == "hospenthal"
replace gdename = "jouxtens mezery" if gdename == "jouxtens"
replace gdename = "kuesnacht (zh)" if gdename == "kusnacht"
replace gdename = "koeniz" if gdename == "niederwangen"
replace gdename = "saint ursanne" if gdename == "st. ursanne"
replace gdename = "uitikon am albis" if gdename == "uitikon a. albis"
replace gdename = "corsier sur vevey" if gdename == "corsier/vevey"
replace gdename = "kuesnacht (zh)" if gdename == "goldbach kuesnacht"
replace gdename = "herrliberg" if gdename == "heerliberg"
replace gdename = "kuesnacht (zh)" if gdename == "kosnacht"
replace gdename = "le mont sur lausanne" if gdename == "le-mont-sur-lausanne"
replace gdename = "neukirch an der thur" if gdename == "neukirch"
replace gdename = "arth" if gdename == "oberarth"
replace gdename = "zuerich" if gdename == "ziirich"
replace gdename = "nendaz" if gdename == "aproz"
replace gdename = "chene bougeries" if gdename == "chenes bougeries"
replace gdename = "montreux" if gdename == "chernex"
replace gdename = "chene bougeries" if gdename == "chene-bougeries"
replace gdename = "sarnen" if gdename == "kaegiswil"
replace gdename = "koeniz" if gdename == "niederscherli"
replace gdename = "ostermundigen" if gdename == "ostermundingen"
replace gdename = "st. gallen" if gdename == "sankt gallen"
replace gdename = "stein am rhein" if gdename == "stein a. rhein"
replace gdename = "la tour de peilz" if gdename == "tour de peilz"
replace gdename = "villars sur glane" if gdename == "villars s/glane"
replace gdename = "gorgier" if gdename == "chez le bart"
replace gdename = "wohlen bei bern" if gdename == "hinterkappelen"
replace gdename = "kriens" if gdename == "krienz"
replace gdename = "malters" if gdename == "matters"
replace gdename = "muenchenstein" if gdename == "neuewelt"
replace gdename = "lancy" if gdename == "petit lancy (ge)"
replace gdename = "lancy" if gdename == "petit laney"
replace gdename = "rueegsau" if gdename == "rueegsauschachen"
replace gdename = "rueegsau" if gdename == "ruegsauschachen"
replace gdename = "sonceboz sombeval" if gdename == "sonceboz"
replace gdename = "saint prex" if gdename == "st. prex"
replace gdename = "saint legier la chiesaz" if gdename == "st legier la chiesaz"
replace gdename = "tafers" if gdename == "tavel"
replace gdename = "villaz saint pierre" if gdename == "villaz-st-pierre"
replace gdename = "boudry" if gdename == "areuse"
replace gdename = "attalens" if gdename == "attalens/tatroz"
replace gdename = "boettstein" if gdename == "bottstein"
replace gdename = "chateau d oex" if gdename == "chateau d åx"
replace gdename = "crans montana" if gdename == "crans-montana"
replace gdename = "duebenedorf" if gdename == "dubendorf"
replace gdename = "wigoltingen" if gdename == "hasli wigoltingen"
replace gdename = "kilchberg (zh)" if gdename == "kilchberg/zuerich"
replace gdename = "landeron combes" if gdename == "le landeron combes"
replace gdename = "niederhelfenschwil" if gdename == "lenggenwil"
replace gdename = "montilliez" if gdename == "montilier"
replace gdename = "vaz/obervaz" if gdename == "vaz"
replace gdename = "saint cergue" if gdename == "saint-cergue"
replace gdename = "torricella taverne" if gdename == "torricella-taverne"
replace gdename = "koeniz" if gdename == "thoerishaus"
replace gdename = "visp" if gdename == "viege"
replace gdename = "bern" if gdename == "bern/buempliz"
replace gdename = "pregny chambesy" if gdename == "chambesy pregny"
replace gdename = "montagny (fr)" if gdename == "cousset"
replace gdename = "crans montana" if gdename == "crans s/sierre"
replace gdename = "crans montana" if gdename == "crans-montana"
replace gdename = "faellanden" if gdename == "fallanden"
replace gdename = "koeniz" if gdename == "gasel"
replace gdename = "duedingen" if gdename == "guin"
replace gdename = "hilterfingen" if gdename == "hunibach"
replace gdename = "evilard" if gdename == "magglingen"
replace gdename = "magliaso" if gdename == "magliasco"
replace gdename = "chardonne" if gdename == "mont pelerin"
replace gdename = "zuerich" if gdename == "oerlikon"
replace gdename = "schaffhausen" if gdename == "schaff hausen"
replace gdename = "le chenit" if gdename == "sentier"
replace gdename = "saint legier la chiesaz" if gdename == "st. legier la chiesaz"
replace gdename = "duernten" if gdename == "tann"
replace gdename = "quarten" if gdename == "unterterzen"
replace gdename = "zollikon" if gdename == "zolliken"
replace gdename = "kandergrund" if gdename == "blausee"
replace gdename = "neuenegg" if gdename == "brueggelbach"
replace gdename = "bern" if gdename == "buempliz"
replace gdename = "bussigny" if gdename == "bussigny s/ morges"
replace gdename = "meinier" if gdename == "carre d amont"
replace gdename = "chigny" if gdename == "chigny s/morges"
replace gdename = "coldrerio" if gdename == "colderio"
replace gdename = "muri bei bern" if gdename == "giimligen"
replace gdename = "giffers" if gdename == "guglera"
replace gdename = "freienbach" if gdename == "hurden"
replace gdename = "kuessnacht (zh)" if gdename == "kilsnacht"
replace gdename = "kuessnacht (sz)" if gdename == "kussnacht"
replace gdename = "bagnes" if gdename == "le chable"
replace gdename = "crans montana" if gdename == "montana station"
replace gdename = "ruemlang" if gdename == "rumlang"
replace gdename = "simplon" if gdename == "simplon dorf"
replace gdename = "sonceboz sombeval" if gdename == "sombeval"
replace gdename = "saint legier la chiesaz" if gdename == "st. legier"
replace gdename = "egnach" if gdename == "steinebrunn"
replace gdename = "saint imier" if gdename == "st lmier"
replace gdename = "taegerwilen" if gdename == "tagerwilen"
replace gdename = "glarus nord" if gdename == "ziegelbruecke"
replace gdename = "zollikon" if gdename == "zollikon (zch.)"
replace gdename = "carouge (ge)" if gdename == "carouge geneve"
replace gdename = "ollon" if gdename == "chesieres"
replace gdename = "collombey muraz" if gdename == "collombey"
replace gdename = "crans pres celigny" if gdename == "crans/celigny"
replace gdename = "crans montana" if gdename == "crans chermignon"
replace gdename = "grenchen" if gdename == "grenchen (sol.)"
replace gdename = "kuessnacht (sz)" if gdename == "immensee"
replace gdename = "lausanne" if gdename == "lau sanne"
replace gdename = "faido" if gdename == "lavorgo"
replace gdename = "le chatelard" if gdename == "les avants"
replace gdename = "mont sur rolle" if gdename == "mont-sur-rolle"
replace gdename = "collombey muraz" if gdename == "muraz"
replace gdename = "oetwil am see" if gdename == "oetwil a. see"
replace gdename = "sumvitg" if gdename == "rabius"
replace gdename = "worb" if gdename == "ruefenacht worb"
replace gdename = "sainte croix" if gdename == "ste croix (vd)"
replace gdename = "torricella taverne" if gdename == "torricella-taverne"
replace gdename = "tenero contra" if gdename == "tenero-contra"
replace gdename = "veyras" if gdename == "veyraz/sierre"
replace gdename = "waltensburg/vuorz" if gdename == "waltensburg"
replace gdename = "schattenhalb" if gdename == "willigen"
replace gdename = "adligenswil" if gdename == "adligenschwil"
replace gdename = "riedholz" if gdename == "attisholz"
replace gdename = "basserdorf" if gdename == "basserdorf"
replace gdename = "port valais" if gdename == "bouveret"
replace gdename = "buelach" if gdename == "bulach"
replace gdename = "bussigny" if gdename == "bussigny sur morges"
replace gdename = "schuebelbach" if gdename == "buttikon"
replace gdename = "castagnola" if gdename == "castagnolla"
replace gdename = "chatel saint denis" if gdename == "chatel st denis"
replace gdename = "crans montana" if gdename == "crans-montana"
replace gdename = "ormont dessus" if gdename == "diablerets"
replace gdename = "sumiswald" if gdename == "gruenen/sumiswald"
replace gdename = "muri bei bern" if gdename == "gumlingen"
replace gdename = "koeniz" if gdename == "liebefeld/koeniz"
replace gdename = "lueterkofen" if gdename == "luterkofen"
replace gdename = "luetzelflueh" if gdename == "lutzelflueh"
replace gdename = "crans montana" if gdename == "montana station"
replace gdename = "neckertal" if gdename == "nassen"
replace gdename = "neuchatel" if gdename == "neuenburg"
replace gdename = "medels im rheinwald" if gdename == "rheinwald"
replace gdename = "schoenenwerd" if gdename == "schonenwerd"
replace gdename = "la grande beroche" if gdename == "st aubin sauges"
replace gdename = "sainte croix" if gdename == "ste croix"
replace gdename = "saint george" if gdename == "st george"
replace gdename = "thonex" if gdename == "thoenex"
replace gdename = "trun" if gdename == "truns"
replace gdename = "vandoeuvres" if gdename == "vandreuvres"
replace gdename = "collonge bellerive" if gdename == "vesenaz (ge)"
replace gdename = "schwende" if gdename == "weissbad"
replace gdename = "wartau" if gdename == "weite"
replace gdename = "brig" if gdename == "brigue"
replace gdename = "chavannes pres renens" if gdename == "chavannes renens"
replace gdename = "chevenez" if gdename == "chevenaz"
replace gdename = "cologny" if gdename == "cologny geneve"
replace gdename = "egg" if gdename == "esslingen"
replace gdename = "kuesnacht (zh)" if gdename == "kuesnacht"
replace gdename = "le grand saconnex" if gdename == "le grand-saconnex"
replace gdename = "hermetschwil staffeln" if gdename == "bremgarten"
replace gdename = "wetzikon (zh)" if gdename == "kempten wetzikon"
replace gdename = "boettstein" if gdename == "kleindottingen"
replace gdename = "kuesnacht (zh)" if gdename == "kuesnacbt"
replace gdename = "kuesnacht (zh)" if gdename == "kuesnacht (zch.)"
replace gdename = "kuessnacht am rigi" if gdename == "kuessnacht a. rigi"
replace gdename = "le chenit" if gdename == "l orient"
replace gdename = "evilard" if gdename == "macolin"
replace gdename = "montreux" if gdename == "montreux les planches"
replace gdename = "moerschwil" if gdename == "morschwil"
replace gdename = "lauterbrunnen" if gdename == "muerren"
replace gdename = "pratteln" if gdename == "prattein"
replace gdename = "schaffhausen" if gdename == "schauffhausen"
replace gdename = "montreux" if gdename == "territet montreux"
replace gdename = "collonge bellerive" if gdename == "vesenaz geneve"
replace gdename = "zug" if gdename == "zug oberwil"
replace gdename = "bre aldesago" if gdename == "aldesago"
replace gdename = "berguen/bravuogn" if gdename == "berguen"
replace gdename = "bern" if gdename == "bern bumpliz"
replace gdename = "pregny chambesy" if gdename == "chambesy (ge)"
replace gdename = "chezard saint martin" if gdename == "chezard st martin"
replace gdename = "chigny" if gdename == "chigny/morges"
replace gdename = "corsier sur vevey" if gdename == "corsier s/vevey"
replace gdename = "davesco soragno" if gdename == "davesco"
replace gdename = "emmebruecke" if gdename == "emmenbriicke"
replace gdename = "muttenz" if gdename == "freidorf muttenz"
replace gdename = "les geneveys sur coffrane" if gdename == "geneveys sur coffrane"
replace gdename = "emmen" if gdename == "gerliswil"
replace gdename = "opfikon" if gdename == "glattbrugg opfikon"
replace gdename = "reichenbach im kandertal" if gdename == "kandertal"
replace gdename = "feuerthalen" if gdename == "langwiesen"
replace gdename = "montreux chatelard" if gdename == "le chatelard montreux"
replace gdename = "ormont dessus" if gdename == "les diablerets"
replace gdename = "leuk" if gdename == "leuk susten"
replace gdename = "lonay" if gdename == "lonay/morges"
replace gdename = "bregaglia" if gdename == "maloja"
replace gdename = "montagny pres yverdon" if gdename == "montagny/yverdon"
replace gdename = "nesslau" if gdename == "neu st johann"
replace gdename = "wald (zh)" if gdename == "neuthal wald"
replace gdename = "pfaeffikon" if gdename == "pfaffikon"
replace gdename = "veyrier" if gdename == "pinchat veyrier"
replace gdename = "bad ragaz" if gdename == "ragaz"
replace gdename = "zell (zh)" if gdename == "rikon"
replace gdename = "schuepfheim" if gdename == "schupfheim"
replace gdename = "sils im engadin/segl" if gdename == "sils maria"
replace gdename = "sion" if gdename == "sitten"
replace gdename = "st. moritz" if gdename == "st mortiz"
replace gdename = "st. gallen" if gdename == "st, gallen"
replace gdename = "taeuffelen" if gdename == "tauffelen"
replace gdename = "treytorrens (payerne)" if gdename == "treytorrents"
replace gdename = "duedingen" if gdename == "villars les joncs"
replace gdename = "zofingen" if gdename == "zofmgen"
replace gdename = "bre aldesago" if gdename == "aldesago di bre"
replace gdename = "alpnach" if gdename == "alpnachstad"
replace gdename = "quinto" if gdename == "ambri di quinto"
replace gdename = "wartau" if gdename == "azmoos"
replace gdename = "baden" if gdename == "baden (aarg.)"
replace gdename = "beatenberg" if gdename == "beatenbucht"
replace gdename = "buchs (ag)" if gdename == "buchs bei aarau"
replace gdename = "brusino arsizio" if gdename == "busto arsizio"
replace gdename = "genthod" if gdename == "creux de genthod"
replace gdename = "genthod" if gdename == "creux de genthod"
replace gdename = "deisswil bei muenchenbuchsee" if gdename == "deisswil"
replace gdename = "basel" if gdename == "dr. basel"
replace gdename = "marthalen" if gdename == "ellikon"
replace gdename = "emmen" if gdename == "emmenbruecke emmen"
replace gdename = "ennetbuergen" if gdename == "ennetburgen"
replace gdename = "erlenbach (zh)" if gdename == "erlenhach"
replace gdename = "flims" if gdename == "films"
replace gdename = "geneve" if gdename == "ganeve"
replace gdename = "muri bei bern" if gdename == "guemligen (be)"
replace gdename = "reichenbach im kandertal" if gdename == "kiental"
replace gdename = "boettstein" if gdename == "kleindoettingen"
replace gdename = "zell (zh)" if gdename == "kollbrunn zell"
replace gdename = "kuettigen" if gdename == "kuettigen rombach"
replace gdename = "langnau im emmental" if gdename == "langnau i/e."
replace gdename = "le mont sur lausanne" if gdename == "le mont s/lausanne"
replace gdename = "lengnau (be)" if gdename == "lengnau bei biel"
replace gdename = "vaz/obervaz" if gdename == "lenzerheid"
replace gdename = "leuk" if gdename == "loeche"
replace gdename = "teufen (ar)" if gdename == "lustmuehle teufen"
replace gdename = "rochefort" if gdename == "montezillon rochefort"
replace gdename = "uster" if gdename == "naenikon uster"
replace gdename = "neuhausen am rheinfall" if gdename == "neuhausen a. rheinfall"
replace gdename = "otringen" if gdename == "of tringen"
replace gdename = "olten" if gdename == "often"
replace gdename = "buchrain" if gdename == "perlen"
replace gdename = "port" if gdename == "port bei nidau"
replace gdename = "porrentruy" if gdename == "pruntrut"
replace gdename = "riehen" if gdename == "riehen/basel"
replace gdename = "risch" if gdename == "risch-rotkreuz"
replace gdename = "urtenen schoenbuehl" if gdename == "schoenbuehl"
replace gdename = "capriasca" if gdename == "tesserette"
replace gdename = "uetikon am see" if gdename == "uetikon a see"
replace gdename = "zumikon" if gdename == "umikon"
replace gdename = "vandoeuvres" if gdename == "vandmuvres"
replace gdename = "les verrieres" if gdename == "verrieres"
replace gdename = "visp" if gdename == "vieques"
replace gdename = "ollon" if gdename == "villars s/ollon"
replace gdename = "daerstetten" if gdename == "weissenburg"
replace gdename = "wil (sg)" if gdename == "wil (st. g.)"
replace gdename = "niederhelfenschwil" if gdename == "zuckenriet"
replace gdename = "allschwil" if gdename == "alischwil"
replace gdename = "quinto" if gdename == "ambri"
replace gdename = "arbedo castione" if gdename == "arbedo gastione"
replace gdename = "arlesheim" if gdename == "arlesbeim"
replace gdename = "gossau (sg)" if gdename == "arnegg gossau"
replace gdename = "arth" if gdename == "arth goldau"
replace gdename = "wiesendangen" if gdename == "attikon"
replace gdename = "basel" if gdename == "bale"
replace gdename = "bellinzona" if gdename == "bellinzona e lugano"
replace gdename = "bellinzona" if gdename == "bellinzone"
replace gdename = "bern" if gdename == "bern rossfeld"
replace gdename = "lugano" if gdename == "biogno"
replace gdename = "malters" if gdename == "blatten malters"
replace gdename = "bourg saint pierre" if gdename == "bourg-st-pierre"
replace gdename = "le chenit" if gdename == "brassus"
replace gdename = "bremgarten bei bern" if gdename == "bremgarten b. bern"
replace gdename = "ruete" if gdename == "bruelisau"
replace gdename = "buchs (ag)" if gdename == "buchs b. aarau"
replace gdename = "bussigny pres lausanne" if gdename == "bussigny/lausanne"
replace gdename = "lugano" if gdename == "casserate"
replace gdename = "davesco soragno" if gdename == "castel do davesco"
replace gdename = "le cerneux pequignot" if gdename == "cerneux pequignot"
replace gdename = "chateau d oex" if gdename == "chaeteau d oex"
replace gdename = "rochefort" if gdename == "chambrelien rochefort"
replace gdename = "orsieres" if gdename == "champex"
replace gdename = "gals" if gdename == "chateau de thielle"
replace gdename = "chateau d oex" if gdename == "chateau d tex"
replace gdename = "chene bougeries" if gdename == "chenebougeries"
replace gdename = "lucens" if gdename == "chesalles s/moudon"
replace gdename = "chezard saint martin" if gdename == "chezard st martin"
replace gdename = "chur" if gdename == "chur zuerich"
replace gdename = "commugny" if gdename == "commungy"
replace gdename = "corseaux" if gdename == "corseaux sur vevey"
replace gdename = "courrendlin" if gdename == "courrendin"
replace gdename = "crans pres celigny" if gdename == "crans (vd)"
replace gdename = "crans montana" if gdename == "crans s. sierre"
replace gdename = "lutry" if gdename == "croix s/lutry"
replace gdename = "stettlen" if gdename == "deisswil stettlen"
replace gdename = "meilen" if gdename == "dr. feldmeilen"
replace gdename = "muri bei bern" if gdename == "dr., muri/bern"
replace gdename = "bern" if gdename == "dr., wabern"
replace gdename = "egerkingen" if gdename == "egerkinden"
replace gdename = "flawil" if gdename == "egg flawil"
replace gdename = "flawil" if gdename == "egg flawil"
replace gdename = "einsiedeln" if gdename == "euthal"
replace gdename = "littau" if gdename == "fluhmuehle"
replace gdename = "kirchberg (sg)" if gdename == "gaehwil"
replace gdename = "tafers" if gdename == "galtern"
replace gdename = "geneve" if gdename == "gen, ve"
replace gdename = "genthod" if gdename == "genthod bellevue"
replace gdename = "opfikon" if gdename == "glattburg"
replace gdename = "grosshoechstetten" if gdename == "grosshochstetten"
replace gdename = "volketswil" if gdename == "hegnau"
replace gdename = "heremence" if gdename == "heremance"
replace gdename = "altstaetten" if gdename == "hinterforst"
replace gdename = "daenikon" if gdename == "huenikon"
replace gdename = "bauen" if gdename == "isleten"
replace gdename = "ittigen" if gdename == "ittigen bolligen"
replace gdename = "horw" if gdename == "kastanienbaum horw"
replace gdename = "kilchberg (zh)" if gdename == "kilchberg (zch.)"
replace gdename = "kuesnacht (zh)" if gdename == "kues nacht (zh)"
replace gdename = "kuettigen" if gdename == "kuttigen"
replace gdename = "saint martin (vs)" if gdename == "la luette"
replace gdename = "leuk" if gdename == "la souste de loeche"
replace gdename = "langnau im emmental" if gdename == "langnau (be)"
replace gdename = "langnau bei reiden" if gdename == "langnau b. reiden"
replace gdename = "lausanne" if gdename == "lausannne"
replace gdename = "waldkirch" if gdename == "lehn waldkirch"
replace gdename = "lenzburg" if gdename == "lenzbourg"
replace gdename = "arzier le muids" if gdename == "les muids"
replace gdename = "luetzelflueh" if gdename == "luetzelffueh"
replace gdename = "lugano" if gdename == "lugano und oberengstringen"
replace gdename = "teufen (ar)" if gdename == "teufen"
replace gdename = "luetisburg" if gdename == "lutisberg"
replace gdename = "madiswil" if gdename == "madliswil"
replace gdename = "maroggia" if gdename == "maroggia melano"
replace gdename = "eschenbach (lu)" if gdename == "mettlen inwil"
replace gdename = "muenchenstein" if gdename == "miinchenstein"
replace gdename = "troistorrents" if gdename == "morgins"
replace gdename = "chateau d oex" if gdename == "moulins pres chateau d oex"
replace gdename = "mosnang" if gdename == "muehlrueti mosnang"
replace gdename = "muri bei bern" if gdename == "muri.b.bern"
replace gdename = "niederwil (ag)" if gdename == "niederwil"
replace gdename = "neyruz (fr)" if gdename == "neyrus"
replace gdename = "koeniz" if gdename == "niederwangen koeniz"
replace gdename = "uzwil" if gdename == "niederzuwil"
replace gdename = "gossau (sg)" if gdename == "ober arnegg"
replace gdename = "quarten" if gdename == "oberterzen"
replace gdename = "zug" if gdename == "oberwil zug"
replace gdename = "oetwil an der limmat" if gdename == "oetwil a. limmat"
replace gdename = "ormont dessus" if gdename == "ormond dessus"
replace gdename = "val de ruz" if gdename == "paquier"
replace gdename = "perly certoux" if gdename == "perly"
replace gdename = "lancy" if gdename == "petit lancy"
replace gdename = "geneve" if gdename == "petit saconnex (ge)"
replace gdename = "veyrier" if gdename == "petit veyrier"
replace gdename = "pfaeffikon" if gdename == "pf aeffikon"
replace gdename = "faellanden" if gdename == "pfaffhausen"
replace gdename = "quinto" if gdename == "piotta di quinto"
replace gdename = "sion" if gdename == "pont de la morge/sion"
replace gdename = "mont vully" if gdename == "mont-vully"
replace gdename = "prilly" if gdename == "prilly sur lausanne"
replace gdename = "lancy" if gdename == "pt lancy"
replace gdename = "zell (zh)" if gdename == "raemismuehle zell"
replace gdename = "raron" if gdename == "rarogne"
replace gdename = "regensberg" if gdename == "regenberg"
replace gdename = "riehen" if gdename == "riehen basel"
replace gdename = "zell (zh)" if gdename == "rikon zell"
replace gdename = "romanshorn" if gdename == "romanshom"
replace gdename = "sennwald" if gdename == "salaz"
replace gdename = "sennwald" if gdename == "salez"
replace gdename = "lutry" if gdename == "savuit sur lutry"
replace gdename = "s chanf" if gdename == "scanfs"
replace gdename = "schaffhausen" if gdename == "schffhausen"
replace gdename = "seengen" if gdename == "seegen"
replace gdename = "seftigen" if gdename == "seftingen"
replace gdename = "lully (fr)" if gdename == "seiry sur estavayer"
replace gdename = "signy avenex" if gdename == "signy"
replace gdename = "saint aubin (fr)" if gdename == "st. aubin"
replace gdename = "stansstad" if gdename == "stanstad"
replace gdename = "oberrohrdorf" if gdename == "staretschwil"
replace gdename = "wetzikon" if gdename == "stegen wetzikon"
replace gdename = "fischenthal" if gdename == "steg fischenthal"
replace gdename = "turbenthal" if gdename == "steinenbach turbenthal"
replace gdename = "saint legier la chiesaz" if gdename == "st legier (vd)"
replace gdename = "sutz lattrigen" if gdename == "sutz"
replace gdename = "aadorf" if gdename == "taennikon"
replace gdename = "taeuffelen" if gdename == "taeufellen"
replace gdename = "turbenthal" if gdename == "thurbenthal"
replace gdename = "tramelan" if gdename == "tramelan dessous"
replace gdename = "wohlen bei bern" if gdename == "uettlingen wohlen"
replace gdename = "uetikon am see" if gdename == "uitikon a. see"
replace gdename = "unteraegeri" if gdename == "unterageri"
replace gdename = "uetikon am see" if gdename == "utikon am see"
replace gdename = "vechigen" if gdename == "utzigen"
replace gdename = "vernamiege" if gdename == "varnamiege"
replace gdename = "vaux sur morges" if gdename == "vaux/morges"
replace gdename = "vevey" if gdename == "vevey et nyon"
replace gdename = "worb" if gdename == "vielbringen"
replace gdename = "lugano" if gdename == "vignello"
replace gdename = "villars sous yens" if gdename == "villars s/yens"
replace gdename = "villarsel sur marly" if gdename == "villarsel s/ marly"
replace gdename = "vira (gambarogno)" if gdename == "vira cambarogno"
replace gdename = "mont vully" if gdename == "vully"
replace gdename = "scuol" if gdename == "vulpera tarasp"
replace gdename = "saint aubin (fr)" if gdename == "vve, st aubin"
replace gdename = "koeniz" if gdename == "wabern (be)"
replace gdename = "koeniz" if gdename == "wabern koeniz"
replace gdename = "wallisellen" if gdename == "walisellen"
replace gdename = "wallisellen" if gdename == "walliselen"
replace gdename = "bergdietikon" if gdename == "wiesenthal bergdietikon"
replace gdename = "schattenhalb" if gdename == "willigen schattenhalb"
replace gdename = "wynau" if gdename == "winau"
replace gdename = "wolfhalden" if gdename == "wollhalden"
replace gdename = "wuerenlos" if gdename == "wuerenloos"
replace gdename = "wynau" if gdename == "wynau bern"
replace gdename = "zuerich" if gdename == "zuerich hoengg"
replace gdename = "arbedo castione" if gdename == "arbede"
replace gdename = "ascona" if gdename == "ascona moscia"
replace gdename = "kloten" if gdename == "augwil kloten"
replace gdename = "nendaz" if gdename == "baar/nendaz"
replace gdename = "basel" if gdename == "base]"
replace gdename = "belmont sur lausanne" if gdename == "belmont s lausanne"
replace gdename = "belp" if gdename == "belp/bern"
replace gdename = "bueren an der aare" if gdename == "bueren a. a."
replace gdename = "burgdorf" if gdename == "buergdorf"
replace gdename = "gurmels" if gdename == "cormondes"
replace gdename = "duebendorf" if gdename == "diibendorf"
replace gdename = "le locle" if gdename == "dr., le lode"
replace gdename = "lausanne" if gdename == "lausanne et genese"
replace gdename = "altstaetten" if gdename == "altstetten"
replace gdename = "silenen" if gdename == "amsteg"
replace gdename = "airolo" if gdename == "ariolo"
replace gdename = "avenches" if gdename == "avanches"
replace gdename = "balsthal" if gdename == "balstahl"
replace gdename = "nendaz" if gdename == "basse nendaz"
replace gdename = "bellevue" if gdename == "bellevue geneve"
replace gdename = "bern" if gdename == "bern buemplitz"
replace gdename = "bern" if gdename == "bern/biimpliz"
replace gdename = "beromuenster" if gdename == "beromunster"
replace gdename = "affoltern am albis" if gdename == "affoltern a. a."
replace gdename = "altstaetten" if gdename == "alststatten"
replace gdename = "alvaneu" if gdename == "alvenau bad"
replace gdename = "anieres" if gdename == "ani eres"
replace gdename = "arlesheim" if gdename == "arleshei m"
replace gdename = "arlesheim" if gdename == "arlesheim (bld.)"
replace gdename = "baeretswil" if gdename == "baeretzwil"
replace gdename = "beinwil am see" if gdename == "beinwil a s"
replace gdename = "beinwil am see" if gdename == "beinwill am see"
replace gdename = "beromuenster" if gdename == "beromuester"
replace gdename = "bettlach" if gdename == "bettfach"
replace gdename = "evilard" if gdename == "(leubringen)"
replace gdename = "adelboden" if gdename == "adelhoden"
replace gdename = "adliswil" if gdename == "adliswii"
replace gdename = "affeltrangen" if gdename == "affeltragen"
replace gdename = "agettes" if gdename == "agettes (vs)"
replace gdename = "ruedtligen alchenflueh" if gdename == "alchenflueh"
replace gdename = "allschwil" if gdename == "allscbwil"
replace gdename = "allschwil" if gdename == "allschwill"
replace gdename = "alpnach" if gdename == "alnnach"
replace gdename = "altstaetten" if gdename == "altstetten (zh)"
replace gdename = "alvaneu" if gdename == "albula"
replace gdename = "bremgarten (ag)" if gdename == "anglikon (ag)"
replace gdename = "anieres" if gdename == "anieres geneve"
replace gdename = "arlesheim" if gdename == "arles heim"
replace gdename = "arlesheim" if gdename == "arlesheim bl"
replace gdename = "arlesheim" if gdename == "arleslteim"
replace gdename = "arth" if gdename == "artb"
replace gdename = "riedholz" if gdename == "attisholz (so)"
replace gdename = "au (sg)" if gdename == "au (st. g.)"
replace gdename = "kloten" if gdename == "augwil/b. kloten"
replace gdename = "waedenswil" if gdename == "au waedenswil"
replace gdename = "bad ragaz" if gdename == "bad ragan"
replace gdename = "bad ragaz" if gdename == "bad ragaz st. gallen"
replace gdename = "freienbach" if gdename == "baech (sz)"
replace gdename = "freienbach" if gdename == "baech freienbach"
replace gdename = "baetterkinden" if gdename == "baetterkindeu"
replace gdename = "balsthal" if gdename == "balstal (so)"
replace gdename = "basel" if gdename == "basal"
replace gdename = "basel" if gdename == "base)"
replace gdename = "baetterkinden" if gdename == "batterkinden"
replace gdename = "beinwil am see" if gdename == "beinwil/see"
replace gdename = "bellinzona" if gdename == "bellineona"
replace gdename = "neuchatel" if gdename == "sjeuchatel"
replace gdename = "saint maurice" if gdename == "st maurice"
replace gdename = "geneve" if gdename == "(, eneve"
replace gdename = "genege" if gdename == "( eneva"
replace gdename = "rapperswil" if gdename == ")tapperswil"
replace gdename = "bellinzona" if gdename == ")teilinzona"
replace gdename = "lugano" if gdename == ", iiugano"
replace gdename = "geneve" if gdename == ". genthod geneve"
replace gdename = "landquart" if gdename == ". landquart"
replace gdename = "pregny chambesy" if gdename == ".pregny"
replace gdename = "carouge (ge)" if gdename == "acacias"
replace gdename = "carouge (ge)" if gdename == "acacias (ge)"
replace gdename = "aegerten" if gdename == "aergerten"
replace gdename = "aesch bei birmensdorf" if gdename == "aesch b/birmensdorf"
replace gdename = "aesch (bl)" if gdename == "aesch/arlesheim"
replace gdename = "st. ursen" if gdename == "aeschlenberg"
replace gdename = "haettenschwil" if gdename == "aettenschwil"
replace gdename = "affoltern im emmental" if gdename == "affoltern (be)"
replace gdename = "cugnasco" if gdename == "agarone"
replace gdename = "agno" if gdename == "agno/lugano"
replace gdename = "muzzano" if gdename == "agnuzzo (lugano)"
replace gdename = "vernier" if gdename == "aire/vernier"
replace gdename = "allschwil" if gdename == "allachwil"
replace gdename = "allschwil" if gdename == "allschwil (basald.)"
replace gdename = "allschwil" if gdename == "allschwil (bild.)"
replace gdename = "allschwil" if gdename == "allschyvil"
replace gdename = "altstaetten" if gdename == "alstiitten"
replace gdename = "altstaetten" if gdename == "altstetten/zuerich"
replace gdename = "altstaetten" if gdename == "altstetten zuerich"
replace gdename = "altstaetten" if gdename == "altstiitten"
replace gdename = "altstaetten" if gdename == "alttaetten"
replace gdename = "alvaneu" if gdename == "alvaneu dorf (gr)"
replace gdename = "nyon" if gdename == "amex s/ nyon"
replace gdename = "nyon" if gdename == "amex sur nyon"
replace gdename = "andelfingen" if gdename == "andelfingen (zh)"
replace gdename = "bremgarten (ag)" if gdename == "anglikon"
replace gdename = "anieres" if gdename == "anieres (suiffe )"
replace gdename = "nendaz" if gdename == "apnoz"
replace gdename = "boudry" if gdename == "areuse boudry"
replace gdename = "arlesheim" if gdename == "arlesheiin"
replace gdename = "arlesheim" if gdename == "arlesheim (basel land)"
replace gdename = "arlesheim" if gdename == "arlesheim/bl"
replace gdename = "arlesheim" if gdename == "arlesheinl"
replace gdename = "arlesheim" if gdename == "arleshiim"
replace gdename = "arni (ag)" if gdename == "arni isliberg"
replace gdename = "arni (ag)" if gdename == "arni lsliberg"
replace gdename = "arbon" if gdename == "arobn"
replace gdename = "arth" if gdename == "arth schwyz"
replace gdename = "bex" if gdename == "arveyes s. bex"
replace gdename = "waedenswil" if gdename == "au bei waedenswil (zh)"
replace gdename = "aubonne" if gdename == "auhonne"
replace gdename = "avry sur matran" if gdename == "avran sur matran"
replace gdename = "avry devant pont" if gdename == "avrey devant pont"
replace gdename = "avully" if gdename == "avully geneve"
replace gdename = "birmensdorf" if gdename == "b, rmensdorf"
replace gdename = "ried bei kerzers" if gdename == "b. kerzers"
replace gdename = "weinfelden" if gdename == "bachtobel weinfelden"
replace gdename = "wollerau" if gdename == "baech wollerau"
replace gdename = "freienbach" if gdename == "baech freienbach (sz)"
replace gdename = "baetterkinden" if gdename == "baettwil (sol.)"
replace gdename = "balsthal" if gdename == "balstal"
replace gdename = "basel" if gdename == "base/"
replace gdename = "basel" if gdename == "bas el"
replace gdename = "basel" if gdename == "ba s el"
replace gdename = "basel" if gdename == "basel (und paris)"
replace gdename = "basel" if gdename == "basel basel"
replace gdename = "basel" if gdename == "basel und burgdorf"
replace gdename = "basel" if gdename == "basel dornach"
replace gdename = "basel" if gdename == "basle"
replace gdename = "basserdorf" if gdename == "bass t rsdorf"
replace gdename = "kirchberg (sg)" if gdename == "bazenheid (sg)"
replace gdename = "beatenberg" if gdename == "beatenbucht (be)"
replace gdename = "beckenried" if gdename == "bechenried"
replace gdename = "beinwil am see" if gdename == "beinwil a. see (ag)"
replace gdename = "beinwil am see" if gdename == "beinwil a/see (ag)"
replace gdename = "beinwil am see" if gdename == "beinwil ans see"
replace gdename = "beinwil am see" if gdename == "beinwil bei muri (ag)"
replace gdename = "beinwil am see" if gdename == "beinwil muri (ag)"
replace gdename = "belp" if gdename == "beip/bern"
replace gdename = "bellinzona" if gdename == "bel li nzona"
replace gdename = "bellach" if gdename == "belach"
replace gdename = "collonge bellerive" if gdename == "bellerive geneve"
replace gdename = "bellevue" if gdename == "bellesue"
replace gdename = "bellevue" if gdename == "bellevue (sui te)"
replace gdename = "bellevue" if gdename == "bellevue pres geneve"
replace gdename = "bellevue" if gdename == "bellevue/geneve"
replace gdename = "bellinzona" if gdename == "belli nzona"
replace gdename = "bellinzona" if gdename == "bellinzoba"
replace gdename = "belmont sur lausanne" if gdename == "belmont s. lausanne"
replace gdename = "belfaux" if gdename == "beltaux"
replace gdename = "bern" if gdename == "bem/buempliz"
replace gdename = "bern" if gdename == "ber n"
replace gdename = "bergdietikon" if gdename == "berg dietikon"
replace gdename = "berguen/bravuogn" if gdename == "berguen (gr)"
replace gdename = "bern" if gdename == "bernâ"
replace gdename = "bern" if gdename == "bern baemplitz"
replace gdename = "bern" if gdename == "bern bethlehem"
replace gdename = "bern" if gdename == "bern bumplitz"
replace gdename = "berneck" if gdename == "berneck (st. g.)"
replace gdename = "bern" if gdename == "bern steinhoelzli"
replace gdename = "lausanne" if gdename == "bethusy lausanne"
replace gdename = "bettlach" if gdename == "bettlach (sol.)"
replace gdename = "zuerich" if gdename == "dr., zuerich"
replace gdename = "basel" if gdename == "dr., basel"

// further round of manual coding (Geo_coding_gdename2 -- > iteratioin of "stillmissing")
replace gdename = "l abbaye" if gdename == "abbaye"
replace gdename = "affoltern am albis" if gdename == "affoltern a. a, a. h. vollenweider, zum merkur ag, affoltern"
replace gdename = "affoltern am albis" if gdename == "affoltern a.a, ernst tanner ag, affoltern a.a."
replace gdename = "thal" if gdename == "altenrhein"
replace gdename = "augst" if gdename == "angst"
replace gdename = "sainte croix" if gdename == "auberson"
replace gdename = "bassersdorf" if gdename == "basserdorf"
replace gdename = "bern" if gdename == "berne"
replace gdename = "waldkirch" if gdename == "bernhardzell"
replace gdename = "bettwiesen" if gdename == "bettweisen"
replace gdename = "nendaz" if gdename == "beuson nendaz"
replace gdename = "nendaz" if gdename == "beuson/nendaz"
replace gdename = "valbirse" if gdename == "bevillard"
replace gdename = "bex" if gdename == "bex et casablanca"
replace gdename = "crans montana" if gdename == "bluche/randogne"
replace gdename = "lutry" if gdename == "bossieres/lutry"
replace gdename = "bourg saint pierre" if gdename == "bourg st pierre"
replace gdename = "bre aldesago" if gdename == "bre"
replace gdename = "bremblens" if gdename == "bremblens s/ morges"
replace gdename = "ingenbohl" if gdename == "brunnen"
replace gdename = "brusino arsizio" if gdename == "brusino"
replace gdename = "buchs (ag)" if gdename == "buchs b. a, teleferique de zabona s.a., crans s/sierre"
replace gdename = "bueren an der aare" if gdename == "bueren a.d.a, grossmetzgerei bigler ag, bueren a.d.a."
replace gdename = "thunstetten" if gdename == "buetzberg"
replace gdename = "duedingen" if gdename == "bundtels"
replace gdename = "bussigny" if gdename == "bussigny/morges"
replace gdename = "tujetsch" if gdename == "camischolas"
replace gdename = "brusio" if gdename == "campascio"
replace gdename = "brusio" if gdename == "campocologno"
replace gdename = "chene bougeries" if gdename == "cbene bougeries"
replace gdename = "celerina/schlarigna" if gdename == "celerina"
replace gdename = "montreux" if gdename == "chailly"
replace gdename = "montreux chatelard" if gdename == "chamby"
replace gdename = "bagnes" if gdename == "champsec"
replace gdename = "vernier" if gdename == "chatelaine"
replace gdename = "vernier" if gdename == "chatelaine/vernier"
replace gdename = "la chaux de fonds" if gdename == "chaux de fonds"
replace gdename = "ollon" if gdename == "chesieres s/ollon"
replace gdename = "ollon" if gdename == "chesieres/ollon"
replace gdename = "giffers" if gdename == "chevrilles"
replace gdename = "chezard saint martin" if gdename == "chezard st marin"
replace gdename = "chene bourg" if gdename == "chine bourg"
replace gdename = "montreux" if gdename == "clarens"
replace gdename = "montreux" if gdename == "clarens montreux"
replace gdename = "bardonnex" if gdename == "compesieres"
replace gdename = "corcelles le jorat" if gdename == "corcelles/jorat"
replace gdename = "corcelles le jorat" if gdename == "corcelles/payerne"
replace gdename = "corseaux" if gdename == "corseaux s/vevey"
replace gdename = "corzoneso" if gdename == "corzonesco"
replace gdename = "corzoneso" if gdename == "corzoneso piano"
replace gdename = "valcolla" if gdename == "cozza di valcolla"
replace gdename = "crans montana" if gdename == "crans"
replace gdename = "crans montana" if gdename == "crans lens"
replace gdename = "crans montana" if gdename == "crans sur sierre"
replace gdename = "crans montana" if gdename == "crans/lens"
replace gdename = "crans pres celigny" if gdename == "crans/nyon"
replace gdename = "bardonnex" if gdename == "croix de rozon"
replace gdename = "davos" if gdename == "davos wolfgang"
replace gdename = "delemont" if gdename == "dielsberg"
replace gdename = "rapperswil (be)" if gdename == "dieterswil"
replace gdename = "hilterfingen" if gdename == "dr. huenibach"
replace gdename = "mendrisio" if gdename == "dr. mendrisio"
replace gdename = "hausen am albis" if gdename == "ebertswil"
replace gdename = "ollon" if gdename == "ecovets"
replace gdename = "illnau effretikon" if gdename == "effretikon"
replace gdename = "emmen" if gdename == "emmebruecke"
replace gdename = "emmen" if gdename == "emmenbruecke"
replace gdename = "feldis/veulden" if gdename == "feldis"
replace gdename = "meilen" if gdename == "feldmeilen"
replace gdename = "bern" if gdename == "felsenau"
replace gdename = "flims" if gdename == "flims dorf"
replace gdename = "maur" if gdename == "forch"
replace gdename = "maur" if gdename == "forch maur"
replace gdename = "lancy" if gdename == "gd lancy"
replace gdename = "val de ruz" if gdename == "geneveys/coffrane"
replace gdename = "gerra (gambarogno)" if gdename == "gerra gamborogno"
replace gdename = "giez" if gdename == "giez/grandson"
replace gdename = "geneve" if gdename == "ginevra"
replace gdename = "opfikon" if gdename == "glattbrugg"
replace gdename = "arth" if gdename == "goldau"
replace gdename = "thun" if gdename == "goldiwil"
replace gdename = "langnau am albis" if gdename == "gontenbach"
replace gdename = "lancy" if gdename == "grand lancy"
replace gdename = "lancy" if gdename == "grand lancy (ge)"
replace gdename = "le grand saconnex" if gdename == "grand saconnex"
replace gdename = "le grand saconnex" if gdename == "grand saconnex"
replace gdename = "le grand saconnex" if gdename == "grand saconnex (ge)"
replace gdename = "riviera" if gdename == "gresciano"
replace gdename = "saanen" if gdename == "gstaad"
replace gdename = "muri bei bern" if gdename == "guemligen"
replace gdename = "muehleberg" if gdename == "guemmenen"
replace gdename = "sigriswil" if gdename == "gunten"
replace gdename = "merenschwand" if gdename == "hagnau merenschwand"
replace gdename = "hausen am albis" if gdename == "hausen a. a, aktienbaeckerei hausen a. a., hausen a. a."
replace gdename = "hausen am albis" if gdename == "hausen a. a, magnet ag, hausen a. a."
replace gdename = "nendaz" if gdename == "haute nendaz"
replace gdename = "val de ruz" if gdename == "hauts geneveys"
replace gdename = "buchholterberg" if gdename == "heimenschwand"
replace gdename = "waengi" if gdename == "heiterschen"
replace gdename = "henniez" if gdename == "henniez les bains"
replace gdename = "bremgarten (ag)" if gdename == "hermetschwil"
replace gdename = "weggis" if gdename == "hertenstein"
replace gdename = "zuerich" if gdename == "hoengg (zh)"
replace gdename = "hoelstein" if gdename == "holstein"
replace gdename = "warth weiningen" if gdename == "ittingen"
replace gdename = "kaiseraugst" if gdename == "kaisersaugst"
replace gdename = "mauensee" if gdename == "kaltbach"
replace gdename = "koelliken" if gdename == "koellikon"
replace gdename = "koelliken" if gdename == "kolliken"
replace gdename = "koeniz" if gdename == "koniz"
replace gdename = "kreuzlingen" if gdename == "kreuzungen"
replace gdename = "wittenbach" if gdename == "kroenbuehl"
replace gdename = "wittenbach" if gdename == "kronbuehl"
replace gdename = "kuesnacht (zh)" if gdename == "ktisnacht"
replace gdename = "kuessnacht (sz)" if gdename == "kuessnacht a. r, buchdruckerei v. kreienbuehl soehne ag, kuessnacht"
replace gdename = "kuessnacht (sz)" if gdename == "kuessnacht a. r, distillerie raeber ag, kuessnacht a. r."
replace gdename = "kuessnacht (sz)" if gdename == "kuessnacht a. r, heinrich saredi, ing. hoch und tiefbau ag, kuessnacht a. r."
replace gdename = "sainte croix" if gdename == "l auberson"
replace gdename = "cologny" if gdename == "la belotte"
replace gdename = "lutry" if gdename == "la conversion/lutry"
replace gdename = "lutry" if gdename == "la croix s/lutry"
replace gdename = "saint cergue" if gdename == "la cure"
replace gdename = "langendorf" if gdename == "langerdorf"
replace gdename = "berneck" if gdename == "langmoos"
replace gdename = "langnau am albis" if gdename == "langnau a. a, spinnerei langnau, langnau a. a."
replace gdename = "langnau im emmental" if gdename == "langnau i. e, bigler ag, langnau i. e."
replace gdename = "langnau im emmental" if gdename == "langnau i. e, bosshart co ag, langnau i. e."
replace gdename = "langnau im emmental" if gdename == "langnau i. e, gebr. wuethrich ag, langnau i. e."
replace gdename = "langnau im emmental" if gdename == "langnau i. e, gebr. wuthrich ag, langnau i. e."
replace gdename = "langnau im emmental" if gdename == "langnau i. e, immobilien ag, langnau i. e."
replace gdename = "langnau im emmental" if gdename == "langnau i. e, maschinenfabrik liechti co ag, langnau i. e."
replace gdename = "langnau im emmental" if gdename == "langnau i. e, tuchfabrik zuercher cie ag, langnau i. e."
replace gdename = "langnau im emmental" if gdename == "langnau i. e, voegeli moser ag, langnau i. e."
replace gdename = "langnau im emmental" if gdename == "langnau i.e, jakob kipfer ag, langengrund"
replace gdename = "l abbaye" if gdename == "les bioux"
replace gdename = "bullet" if gdename == "les cluds/bullet"
replace gdename = "montreux" if gdename == "les planches"
replace gdename = "koeniz" if gdename == "liebefeld"
replace gdename = "lichtensteig" if gdename == "liechtensteig"
replace gdename = "luzern" if gdename == "liftau"
replace gdename = "lonay" if gdename == "lonay s/morges"
replace gdename = "altstaetten" if gdename == "luechingen"
replace gdename = "luetisburg" if gdename == "luetisberg"
replace gdename = "lugano" if gdename == "lugano cassarate"
replace gdename = "lugano" if gdename == "lugano paradiso"
replace gdename = "teufen (ar)" if gdename == "lustmuehle"
replace gdename = "malters" if gdename == "maltera"
replace gdename = "maennedorf" if gdename == "mannedorf"
replace gdename = "massagno" if gdename == "massgno"
replace gdename = "matten bei interlaken" if gdename == "matten"
replace gdename = "mellingen" if gdename == "melligen"
replace gdename = "sigriswil" if gdename == "merligen"
replace gdename = "mezzovico vira" if gdename == "mezzovico"
replace gdename = "moeriken wildegg" if gdename == "moeriken"
replace gdename = "zuerich" if gdename == "mombasa zurich"
replace gdename = "muentschemier" if gdename == "monsmier"
replace gdename = "le mont sur lausanne" if gdename == "mont lausanne"
replace gdename = "mont sur rolle" if gdename == "mont s/rolle"
replace gdename = "le mont sur lausanne" if gdename == "mont sur lausanne"
replace gdename = "mont sur rolle" if gdename == "mont/rolle"
replace gdename = "chatel sur montsalvens" if gdename == "montsalvens/bulle"
replace gdename = "morges" if gdename == "morgen"
replace gdename = "moeriken wildegg" if gdename == "moriken wildegg"
replace gdename = "muenchenstein" if gdename == "muechenstein"
replace gdename = "arzier le muids" if gdename == "muids/nyon"
replace gdename = "aarwangen" if gdename == "mumenthal"
replace gdename = "muri bei bern" if gdename == "muni b. bern"
replace gdename = "sierre" if gdename == "muraz/sierre"
replace gdename = "quarten" if gdename == "murg"
replace gdename = "muri bei bern" if gdename == "muri b. bern"
replace gdename = "lauterbrunnen" if gdename == "murren"
replace gdename = "niederwil (ag)" if gdename == "nesselnbach"
replace gdename = "allschwil" if gdename == "neu allschwil"
replace gdename = "allschwil" if gdename == "neuallschwil"
replace gdename = "neuchatel" if gdename == "neuchaetel"
replace gdename = "nidau" if gdename == "nideau"
replace gdename = "uzwil" if gdename == "niederuzwil"
replace gdename = "pambio noranco" if gdename == "noranco"
replace gdename = "seeberg" if gdename == "obergrasswil"
replace gdename = "wynau" if gdename == "obermurgenthal wynau"
replace gdename = "obersiggenthal" if gdename == "obersiggental"
replace gdename = "vaz/obervaz" if gdename == "obervaz"
replace gdename = "oberwil lieli" if gdename == "oberwil/bremgarten"
replace gdename = "oensingen" if gdename == "oensigen"
replace gdename = "oetwil am see" if gdename == "oetwil a.s, gadola ag, oetwil a.s."
replace gdename = "oppligen" if gdename == "opplingen"
replace gdename = "geneve" if gdename == "oslo geneve"
replace gdename = "poschiavo" if gdename == "ospizio bernina"
replace gdename = "oftringen" if gdename == "otringen"
replace gdename = "pailly" if gdename == "pally"
replace gdename = "bellerive (vd)" if gdename == "paris bellerive"
replace gdename = "geneve" if gdename == "paris et geneve"
replace gdename = "prangins" if gdename == "paris et prangins"
replace gdename = "veyrier" if gdename == "pinchat"
replace gdename = "chiasso" if gdename == "ponte chiasso"
replace gdename = "ronco sopra ascona" if gdename == "porto ronco"
replace gdename = "lugano" if gdename == "povro breganzona"
replace gdename = "pregny chambesy" if gdename == "pregny (ge)"
replace gdename = "bondo" if gdename == "promontognio"
replace gdename = "pully" if gdename == "puffy"
replace gdename = "pully" if gdename == "pullt"
replace gdename = "pully" if gdename == "putty"
replace gdename = "chalais" if gdename == "rechy"
replace gdename = "richenthal" if gdename == "richental"
replace gdename = "nesslau krummenau" if gdename == "rietbad"
replace gdename = "zuerich" if gdename == "riverside zuerich"
replace gdename = "ostermundigen" if gdename == "roerswil"
replace gdename = "roggenburg" if gdename == "roggenbourg"
replace gdename = "romanel sur lausanne" if gdename == "romanel s/lausanne"
replace gdename = "romanshorn" if gdename == "romanhorn"
replace gdename = "kuettigen" if gdename == "rombach"
replace gdename = "ronco sopra ascona" if gdename == "ronco"
replace gdename = "risch" if gdename == "rotkreuz"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolfstetten"
replace gdename = "ruedtligen alchenflueh" if gdename == "ruedtligen (be)"
replace gdename = "worb" if gdename == "ruefenacht"
replace gdename = "worb" if gdename == "ruefenacht/worb"
replace gdename = "tujetsch" if gdename == "rueras"
replace gdename = "worb" if gdename == "rufenacht"
replace gdename = "rubigen" if gdename == "ruhigen"
replace gdename = "lugano" if gdename == "ruvigliana"
replace gdename = "gruyeres" if gdename == "saussivue"
replace gdename = "savosa" if gdename == "savosa paese"
replace gdename = "lutry" if gdename == "savuit"
replace gdename = "masein" if gdename == "schauenstein masein"
replace gdename = "saanen" if gdename == "schoenried"
replace gdename = "schoetz" if gdename == "schuetz"
replace gdename = "scuol" if gdename == "schuls"
replace gdename = "bischofszell" if gdename == "schweizerholz"
replace gdename = "le lieu" if gdename == "sechey"
replace gdename = "sempach" if gdename == "sempach station"
replace gdename = "sempach" if gdename == "sempach station"
replace gdename = "selzach" if gdename == "setzach"
replace gdename = "schuebelbach" if gdename == "siebnen"
replace gdename = "baar" if gdename == "sihlbrugg"
replace gdename = "zweisimmen" if gdename == "simmenthal"
replace gdename = "vechigen" if gdename == "sinneringen"
replace gdename = "solothurn" if gdename == "soluthurn"
replace gdename = "koeniz" if gdename == "spiegel"
replace gdename = "koeniz" if gdename == "spiegel koeniz"
replace gdename = "koeniz" if gdename == "spiegel/koeniz"
replace gdename = "saint blaise" if gdename == "st blaise"
replace gdename = "saint cergue" if gdename == "st cergue"
replace gdename = "saint george" if gdename == "st georges"
replace gdename = "saint legier la chiesaz" if gdename == "st legier/la chiesaz"
replace gdename = "saint leonard" if gdename == "st leonard"
replace gdename = "saint prex" if gdename == "st prex"
replace gdename = "saint blaise" if gdename == "st. blaise"
replace gdename = "st. gallen" if gdename == "st. gallen tablat"
replace gdename = "saint gingolph" if gdename == "st. gingolph"
replace gdename = "saint leonard" if gdename == "st. leonard"
replace gdename = "staefa" if gdename == "staef a"
replace gdename = "oberrohrdorf" if gdename == "starretschwil oberrohrdorf"
replace gdename = "staufen" if gdename == "stauten"
replace gdename = "suhr" if gdename == "suhr (aarg.)"
replace gdename = "flums" if gdename == "tannenbodenalp"
replace gdename = "torricella taverne" if gdename == "taverne torricella"
replace gdename = "tenero contra" if gdename == "tenero"
replace gdename = "montreux" if gdename == "territet"
replace gdename = "freienstein teufen" if gdename == "teufen freienstein"
replace gdename = "thonex" if gdename == "thbnex"
replace gdename = "termen" if gdename == "thermen (vs)"
replace gdename = "thun" if gdename == "thun gwatt"
replace gdename = "thun" if gdename == "thunersee"
replace gdename = "triengen" if gdename == "tirengen"
replace gdename = "wartau" if gdename == "truebbach"
replace gdename = "wartau" if gdename == "truebbach/wartau"
replace gdename = "twann tuescherz" if gdename == "tuescherz"
replace gdename = "grandson" if gdename == "tuileries de grandson"
replace gdename = "uetikon am see" if gdename == "uetikon a. see"
replace gdename = "wohlen bei bern" if gdename == "uettligen"
replace gdename = "uitikon" if gdename == "uitikon a a"
replace gdename = "uitikon" if gdename == "uitikon a. a, immobilien ag konradstrasse zuerich"
replace gdename = "uitikon" if gdename == "uitikon a. a, voehringer ag, uitikon a. a."
replace gdename = "uitikon" if gdename == "uitikon am albis"
replace gdename = "uitikon" if gdename == "uitikon waldegg"
replace gdename = "konolfingen" if gdename == "ursellen"
replace gdename = "sion" if gdename == "uvrier sion"
replace gdename = "pfaefers" if gdename == "vaettis"
replace gdename = "vaz/obervaz" if gdename == "valbella"
replace gdename = "vandoeuvres" if gdename == "vanda uvres"
replace gdename = "vandoeuvres" if gdename == "vandaeuvres"
replace gdename = "vandoeuvres" if gdename == "vandeuvres"
replace gdename = "bagnes" if gdename == "verbier"
replace gdename = "chalais" if gdename == "vercorin"
replace gdename = "crans montana" if gdename == "vermalla"
replace gdename = "montreux" if gdename == "vernex montreux"
replace gdename = "sainte croix" if gdename == "vers chez jaccard"
replace gdename = "lausanne" if gdename == "vers chez les blanc"
replace gdename = "lausanne" if gdename == "vers chez les blancs"
replace gdename = "lausanne" if gdename == "vers chez les blancs"
replace gdename = "villars sainte croix" if gdename == "villars ste croix"
replace gdename = "ollon" if gdename == "villars sur ollon"
replace gdename = "villaz saint pierre" if gdename == "villaz st pierre"
replace gdename = "vuarrens" if gdename == "vuarrengel"
replace gdename = "scuol" if gdename == "vulpera"
replace gdename = "koeniz" if gdename == "wabern"
replace gdename = "koeniz" if gdename == "wabern b. koeniz"
replace gdename = "wangen an der aare" if gdename == "wangen a. a."
replace gdename = "wangen bei olten" if gdename == "wangen b. olten"
replace gdename = "vilters wangs" if gdename == "wangs vilters"
replace gdename = "sumiswald" if gdename == "wasen i. e, rudolf ruch ag, wasen i. e."
replace gdename = "sumiswald" if gdename == "wasen i. e, wiedmer soehne ag, wasen i. e."
replace gdename = "regensdorf" if gdename == "watt regensdorf"
replace gdename = "affoltern im emmental" if gdename == "weier i. e, ernst schaerlig ag grossmetzgerei, weier"
replace gdename = "weingarten" if gdename == "weingarten triltschen"
replace gdename = "grabs" if gdename == "werdenberg"
replace gdename = "wettswil am albis" if gdename == "wettswil a. a, mueller barbieri ag, wettswil a. a."
replace gdename = "wetzikon (zh)" if gdename == "wetzikon"
replace gdename = "wiedlisbach" if gdename == "wiedisbach"
replace gdename = "moeriken wildegg" if gdename == "wildegg"
replace gdename = "moeriken wildegg" if gdename == "wildegg moerikon"
replace gdename = "winterthur" if gdename == "winter thur"
replace gdename = "wolhusen" if gdename == "wohlhusen"
replace gdename = "bubikon" if gdename == "wolfhausen"
replace gdename = "bubikon" if gdename == "wolfhausen bubikon"
replace gdename = "wolfwil" if gdename == "wolfil"
replace gdename = "zuerich" if gdename == "zilrich"
replace gdename = "zuerich" if gdename == "zprich"
replace gdename = "zuerich" if gdename == "zuerieh"
replace gdename = "zuerich" if gdename == "zurich"
replace gdename = "carouge" if gdename == "carouge"
replace gdename = "renens (vd)" if gdename == "renens"
replace gdename = "neyruz (fr)" if gdename == "neyruz"
replace gdename = "waedenswil" if gdename == "schoenenberg"
replace gdename = "birmensdorf (zh)" if gdename == "birmensdorf"
replace gdename = "rekingen (ag)" if gdename == "rekingen"
replace gdename = "illnau effretikon" if gdename == "8308lllnau effretikon"
replace gdename = "langnau am albis" if gdename == "langnau a.a"
replace gdename = "oberriet (sg)" if gdename == "oberriet"
replace gdename = "illnau effretikon" if gdename == "8308iiinau effretikon"
replace gdename = "illnau effretikon" if gdename == "8308iilnau effretikon"
replace gdename = "emmen" if gdename == "emmenbruecke/emmen"
replace gdename = "st. Gallen" if gdename == "9oll st.gallen"
replace gdename = "basel" if gdename == "4oll basel"
replace gdename = "zuerich" if gdename == "dr."
replace gdename = "bremgarten bei bern" if gdename == "bremgarten be"
replace gdename = "wettswil am albis" if gdename == "wettswil a.a"
replace gdename = "zuerich" if gdename == "zollikon (zh)"
replace gdename = "basel" if gdename == "dr."
replace gdename = "ipsach" if gdename == "25631psach"
replace gdename = "carouge (ge)" if gdename == "carouge"
replace gdename = "zuerich" if gdename == "8oll zuerich"
replace gdename = "wiesen (gr)" if gdename == "wiesen"
replace gdename = "illnau effretikon" if gdename == "8308iiinau"
replace gdename = "lausanne" if gdename == "1oll lausanne"
replace gdename = "affoltern am albis" if gdename == "affoltern a.a"
replace gdename = "horw" if gdename == "kastanienbaum/horw"
replace gdename = "ringgenberg (be)" if gdename == "ringgenberg"
replace gdename = "steffisburg" if gdename == "steffisburg dorf"
replace gdename = "le mont sur lausanne" if gdename == "1052"
replace gdename = "oetwil an der limmat" if gdename == "oetwil a.d.l"
replace gdename = "illnau effretikon" if gdename == "8308iilnau"
replace gdename = "ingenbohl" if gdename == "brunnen/ingenbohl"
replace gdename = "ittigen" if gdename == "30631ttigen"
replace gdename = "littau" if gdename == "reussbuehl/littau"
replace gdename = "lachen" if gdename == "lachen sz"
replace gdename = "luzern" if gdename == "lucerne"
replace gdename = "thalwil" if gdename == "8800thalwil"
replace gdename = "uetikon am see" if gdename == "uetikon a.s"
replace gdename = "mesocco" if gdename == "san bernardino/mesocco"
replace gdename = "kuessnacht am rigi" if gdename == "kuessnacht a.r"
replace gdename = "hausen am albis" if gdename == "hausen a.a"
replace gdename = "bern" if gdename == "dr., bern"
replace gdename = "freienbach" if gdename == "pfaeffikon sz"
replace gdename = "chene bougeries" if gdename == "chene bougeries (ge)"
replace gdename = "teufenthal (ag)" if gdename == "teufenthal"
replace gdename = "bardonnex" if gdename == "1257la croix de rozon"
replace gdename = "birmensdorf (zh)" if gdename == "birmensdorf"
replace gdename = "reckingen (vs)" if gdename == "reckingen"
replace gdename = "ried bei moerel" if gdename == "ried bei morel"
replace gdename = "oetwil am see" if gdename == "oetwil a.s"
replace gdename = "grandson" if gdename == "les tuileries/grandson"
replace gdename = "muri bei bern" if gdename == "muri (be)"
replace gdename = "grandson" if gdename == "les tuileries de grandson"
replace gdename = "rueschegg" if gdename == "heubach rueschegg"
replace gdename = "kuesnacht (zh)" if gdename == "kuessnacht (zh)"
replace gdename = "zuerich" if gdename == "zurigo"
replace gdename = "bern" if gdename == "3oll bem"
replace gdename = "renens (vd)" if gdename == "renens"
replace gdename = "langnau im emmental" if gdename == "langnau i.e"
replace gdename = "zuerich" if gdename == "dr., zuerich"
replace gdename = "basel" if gdename == "dr., basel"
replace gdename = "st gallen" if gdename == "dr., st. gallen"

replace gdename = "boettstein" if gdename == "boettstein kleindoettingen"
replace gdename = "grub (ar)" if gdename == "grub"
replace gdename = "zug" if gdename == "oberwil/zug"
replace gdename = "st gallen" if gdename == "st gall"
replace gdename = "veltheim (ag)" if gdename == "veltheim"
replace gdename = "kriens" if gdename == "amlehnstrasse 38, 6oll kriens"
replace gdename = "basel" if gdename == "postfach, 4oll basel"
replace gdename = "gaiserwald" if gdename == "abtwil/gaiserwald"
replace gdename = "risch" if gdename == "rotkreuz/risch"
replace gdename = "lauterbrunnen" if gdename == "wengen/lauterbrunnen"
replace gdename = "lugano" if gdename == "dr., lugano"
replace gdename = "solothurn" if gdename == "dr., solothurn"
replace gdename = "la chaux de fonds" if gdename == "fils, la chaux de fonds"
replace gdename = "motiers (ne)" if gdename == "motiers"
replace gdename = "gommiswald" if gdename == "rieden (sg)"
replace gdename = "treytorrens (payerne)" if gdename == "treytorrens"
replace gdename = "aesch bei birmensdorf" if gdename == "aesch b. birmensdorf"
replace gdename = "mont vully" if gdename == "nant/bas vully f"
replace gdename = "koeniz" if gdename == "schliern"
replace gdename = "hofen" if gdename == "hofen (sh)"
replace gdename = "ruethi (rheintal)" if gdename == "ruethi"
replace gdename = "freienbach" if gdename == "wilen/freienbach"
replace gdename = "horw" if gdename == "st. niklausen/horw"
replace gdename = "langnau am albis" if gdename == "langnau a. a"
replace gdename = "uetikon am see" if gdename == "uetikon a. s"
replace gdename = "ittigen" if gdename == "neuhausweg 5, ittigen bolligen"
replace gdename = "lausanne" if gdename == "rue du bugnon 22, 1oll lausanne"
replace gdename = "oberwil (ag)" if gdename == "oberwil ag"
replace gdename = "wartau" if gdename == "weite/wartau"
replace gdename = "ayent" if gdename == "luc/ayent"
replace gdename = "birmenstorf (ag)" if gdename == "birmenstorf"
replace gdename = "duebendorf" if gdename == "duebenedorf"
replace gdename = "kaiserstuhl" if gdename == "kaiserstuhl (ag)"
replace gdename = "pfaeffikon" if gdename == "pfaeffikon zh"
replace gdename = "bissegg" if gdename == "eutenberg/bissegg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolfstette friedlisberg"
replace gdename = "rickenbach bei wil" if gdename == "rickenbach b. wil"
replace gdename = "koeniz" if gdename == "niederscherli koeniz"
replace gdename = "lutry" if gdename == "bossieres sur lutry"
replace gdename = "versoix" if gdename == "versoix (ge)"
replace gdename = "saint sulpice (vd)" if gdename == "st sulpice vd"
replace gdename = "wohlen bei bern" if gdename == "wohlen be"
replace gdename = "schwyz" if gdename == "seewen/schwyz"
replace gdename = "feusisberg" if gdename == "schindellegi/feusisberg"
replace gdename = "wisen (so)" if gdename == "wisen"
replace gdename = "thun" if gdename == "gwatt"
replace gdename = "lausanne" if gdename == "losanna"
replace gdename = "brugg" if gdename == "brugg (ag)"
replace gdename = "kuesnacht (zh)" if gdename == "dr., kuesnacht"
replace gdename = "lausanne" if gdename == "dr., lausanne"
replace gdename = "langenthal" if gdename == "langenthal (be)"
replace gdename = "fribourg" if gdename == "vve, fribourg"
replace gdename = "borex" if gdename == "borex vd"
replace gdename = "faellanden" if gdename == "faellenden"
replace gdename = "vaz/obervaz" if gdename == "valbella/vaz/obervaz"
replace gdename = "flums" if gdename == "flumserberg"
replace gdename = "fribourg" if gdename == "dr., fribourg"
replace gdename = "monthey" if gdename == "monthey (vs)"
replace gdename = "bern" if gdename == "phil., bern"
replace gdename = "pully" if gdename == "pully (vd)"
replace gdename = "reute (ar)" if gdename == "reute"
replace gdename = "sierre" if gdename == "siders"
replace gdename = "lausanne" if gdename == "r. du bugnon 8, 1oll lausanne"
replace gdename = "nussbaumen" if gdename == "obersiggenthal nussbaumen"
replace gdename = "ormont dessus" if gdename == "les diablerets/ormont dessus"
replace gdename = "speicher" if gdename == "speicherschwendi/speicher"
replace gdename = "zuerich" if gdename == "eschenweg 7, 8oll zuerich"
replace gdename = "l abbaye" if gdename == "les bioux vd"
replace gdename = "schwende" if gdename == "blueemlisalp, wasserauen/schwende"
replace gdename = "bex" if gdename == "les posses/bex"
replace gdename = "muri bei bern" if gdename == "guemligen muri be"
replace gdename = "werthenstein" if gdename == "schachen/werthenstein"
replace gdename = "trachselwald" if gdename == "heimisbach trachselwald"
replace gdename = "montreux" if gdename == "rue du lac 2, clarens montreu"
replace gdename = "corcelles cormondreche" if gdename == "corcelles ne"
replace gdename = "menzingen" if gdename == "edlibach/menzingen"
replace gdename = "gebenstorf" if gdename == "gebensdorf"
replace gdename = "aesch (zh)" if gdename == "aesch zh"
replace gdename = "kirchlindach" if gdename == "herrenschwanden/kirchlindach"
replace gdename = "ballaigues" if gdename == "ballaigues (vd)"
replace gdename = "biglen" if gdename == "biglen (be)"
replace gdename = "oberstammheim" if gdename == "stammheim"
replace gdename = "oberburg" if gdename == "oberburg (be)"
replace gdename = "riehen" if gdename == "riehen (basel)"


preserve
drop if n_gde >1
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo90a.dta", replace
restore 

// a few ambiguous cases (more than one municipality possible)
drop if n_gde == 1

gen gdenameopt1 = ""
gen gdenameopt2 = ""
gen gdenameopt3 = ""
gen gdenameopt4 = ""
gen gdenameopt5 = ""
gen gdenameopt6 = ""
gen gdenameopt7 = ""
gen gdenameopt8 = ""

replace gdenameopt1 = "saint sulpice (ne)" if gdename == "st sulpice"
replace gdenameopt2 = "saint sulpice (vd)" if gdename == "st sulpice"
replace gdenameopt1 = "saint sulpice (ne)" if gdename == "st. sulpice"
replace gdenameopt2 = "saint sulpice (vd)" if gdename == "st. sulpice"
replace gdenameopt1 = "sils im engadin/segl" if gdename == "sils"
replace gdenameopt2 = "sils im domleschg" if gdename == "sils"
replace gdenameopt1 = "ried bei kerzers" if gdename == "ried"
replace gdenameopt2 = "ried bei moerel" if gdename == "ried"
replace gdenameopt3 = "ried bei brig" if gdename == "ried"
replace gdenameopt1 = "sumiswald" if gdename == "wasen"
replace gdenameopt2 = "schleinikon" if gdename == "wasen"
replace gdenameopt1 = "affoltern am albis" if gdename == "affoltern"
replace gdenameopt2 = "affoltern im emmental" if gdename == "affoltern"
replace gdenameopt3 = "zuerich" if gdename == "affoltern"
replace gdenameopt1 = "hausen (ag)" if gdename == "hausen"
replace gdenameopt2 = "hausen am albis" if gdename == "hausen"
replace gdenameopt3 = "hausen bei brugg" if gdename == "hausen"
replace gdenameopt1 = "vuisternens devant romont" if gdename == "vuisternens"
replace gdenameopt2 = "vuisternens en ogoz" if gdename == "vuisternens"
replace gdenameopt1 = "oetwil an der limmat" if gdename == "oetwil"
replace gdenameopt2 = "oetwil am see" if gdename == "oetwil"
replace gdenameopt1 = "muenchwilen (ag)" if gdename == "munchwilen"
replace gdenameopt2 = "muenchwilen (tg)" if gdename == "munchwilen"
replace gdenameopt1 = "kilchberg (zh)" if gdename == "kilchherg"
replace gdenameopt2 = "kilchberg (bl)" if gdename == "kilchherg"
replace gdenameopt1 = "neuchatel" if gdename == "le coudre"
replace gdenameopt2 = "l isle" if gdename == "le coudre"
replace gdenameopt1 = "romanel sur lausanne" if gdename == "romanel"
replace gdenameopt2 = "romanel sur morges" if gdename == "romanel"
replace gdenameopt1 = "uitikon" if gdename == "waldegg"
replace gdenameopt2 = "rickenbach (bl)" if gdename == "waldegg"
replace gdenameopt1 = "kilchberg (zh)" if gdename == "killchberg"
replace gdenameopt2 = "kilchberg (bl)" if gdename == "killchberg"
replace gdenameopt1 = "chavannes de bogis" if gdename == "chavannes"
replace gdenameopt2 = "chavannes des bois" if gdename == "chavannes"
replace gdenameopt3 = "chavannes le chene" if gdename == "chavannes"
replace gdenameopt4 = "chavannes le veyron" if gdename == "chavannes"
replace gdenameopt5 = "chavannes les forts" if gdename == "chavannes"
replace gdenameopt6 = "chavannes pres renens" if gdename == "chavannes"
replace gdenameopt7 = "chavannes sous orsonnens" if gdename == "chavannes"
replace gdenameopt8 = "chavannes sur moudon" if gdename == "chavannes"
replace gdenameopt1 = "wilen (tg)" if gdename == "wilen"
replace gdenameopt2 = "wilen bei neunforn" if gdename == "wilen"
replace gdenameopt3 = "hauptwil gottshaus" if gdename == "wilen"
replace gdenameopt4 = "herisau" if gdename == "wilen"
replace gdenameopt5 = "sarnen" if gdename == "wilen"
replace gdenameopt1 = "langnau am albis" if gdename == "langnau"
replace gdenameopt2 = "langnau bei reiden" if gdename == "langnau"
replace gdenameopt3 = "langnau im emmental" if gdename == "langnau"
replace gdenameopt1 = "cheseaux noreaz" if gdename == "cheseaux"
replace gdenameopt2 = "cheseaux sur lausanne" if gdename == "cheseaux"
replace gdenameopt1 = "morbio inferiore" if gdename == "morbio"
replace gdenameopt2 = "morbio superiore" if gdename == "morbio"
replace gdenameopt1 = "oulens sur lucens" if gdename == "oulens"
replace gdenameopt2 = "oulens sous echallens" if gdename == "oulens"
replace gdenameopt1 = "pfaeffikon" if gdename == "pfaelfikon"
replace gdenameopt2 = "freienbach" if gdename == "pfaelfikon"
replace gdenameopt1 = "roethenbach bei herzogenbuchsee" if gdename == "roethenbach"
replace gdenameopt2 = "roethenbach im emmental" if gdename == "roethenbach"
replace gdenameopt1 = "horw" if gdename == "st niklausen"
replace gdenameopt2 = "kerns" if gdename == "st niklausen"
replace gdenameopt1 = "saint saphorin (lavaux)" if gdename == "st saphorin"
replace gdenameopt2 = "saint saphorin sur morges" if gdename == "st saphorin"
replace gdenameopt1 = "horw" if gdename == "st. niklausen"
replace gdenameopt2 = "kerns" if gdename == "st. niklausen"
replace gdenameopt1 = "torricella taverne" if gdename == "taverne"
replace gdenameopt2 = "les tavernes" if gdename == "taverne"
replace gdenameopt1 = "uttwil" if gdename == "uttwil romanshorn"
replace gdenameopt2 = "romanshorn" if gdename == "uttwil romanshorn"
replace gdenameopt1 = "sumiswald" if gdename == "wasen"
replace gdenameopt2 = "schleinikon" if gdename == "wasen"
replace gdenameopt1 = "zermatt" if gdename == "zermatt et sion"
replace gdenameopt2 = "sion" if gdename == "zermatt et sion"
replace gdenameopt1 = "erlenbach (zh)" if gdename == "erlenbach"
replace gdenameopt2 = "erlenbach im simmental" if gdename == "erlenbach"
replace gdenameopt1 = "muri (ag)" if gdename == "muri"
replace gdenameopt2 = "muri bei bern" if gdename == "muri"
replace gdenameopt1 = "lussy (fr)" if gdename == "lussy"
replace gdenameopt2 = "lussy sur morges" if gdename == "lussy"
replace gdenameopt1 = "schwanden (gl)" if gdename == "schwanden"
replace gdenameopt2 = "schwanden bei brienz" if gdename == "schwanden"
replace gdenameopt1 = "roveredo (gr)" if gdename == "roveredo"
replace gdenameopt2 = "roveredo capriasca" if gdename == "roveredo"
replace gdenameopt1 = "hergiswil (nw)" if gdename == "hergiswil"
replace gdenameopt2 = "hergiswil bei willisau" if gdename == "hergiswil"
replace gdenameopt1 = "hasle (lu)" if gdename == "hasle"
replace gdenameopt2 = "hasle bei burgdorf" if gdename == "hasle"
replace gdenameopt1 = "corsier (ge)" if gdename == "corsier"
replace gdenameopt2 = "corsier sur vevey" if gdename == "corsier"
replace gdenameopt1 = "wiler (loetschen)" if gdename == "wiler"
replace gdenameopt2 = "wiler bei utzenstorf" if gdename == "wiler"
replace gdenameopt1 = "bueren (so)" if gdename == "bueren"
replace gdenameopt2 = "bueren an der aare" if gdename == "bueren"
replace gdenameopt1 = "bellerive (vd)" if gdename == "bellerive"
replace gdenameopt2 = "collonge bellerive" if gdename == "bellerive"
replace gdenameopt1 = "burg (ag)" if gdename == "burg"
replace gdenameopt2 = "burg bei murten" if gdename == "burg"
replace gdenameopt3 = "burg im leimental" if gdename == "burg"
replace gdenameopt1 = "oberhofen (ag)" if gdename == "oberhofen"
replace gdenameopt2 = "oberhofen am thunersee" if gdename == "oberhofen"
replace gdenameopt3 = "oberhofen bei kreuzlingen" if gdename == "oberhofen"
replace gdenameopt1 = "corcelles (be)" if gdename == "corcelles"
replace gdenameopt2 = "corcelles cormondreche" if gdename == "corcelles"
replace gdenameopt3 = "corcelles le jorat" if gdename == "corcelles"
replace gdenameopt4 = "corcelles pres concise" if gdename == "corcelles"
replace gdenameopt5 = "corcelles pres payerne" if gdename == "corcelles"
replace gdenameopt6 = "corcelles sur chavornay" if gdename == "corcelles"
replace gdenameopt1 = "wohlen (ag)" if gdename == "wohlen"
replace gdenameopt2 = "wohlen bei bern" if gdename == "wohlen"
replace gdenameopt1 = "gossau (sg)" if gdename == "gossau"
replace gdenameopt1 = "oberried (fr)" if gdename == "oberried"
replace gdenameopt1 = "stetten (ag)" if gdename == "stetten"
replace gdenameopt1 = "seedorf (be)" if gdename == "seedorf"
replace gdenameopt2 = "gossau (zh)" if gdename == "gossau"
replace gdenameopt2 = "oberried am brienzersee" if gdename == "oberried"
replace gdenameopt2 = "stetten (sh)" if gdename == "stetten"
replace gdenameopt2 = "seedorf (ur)" if gdename == "seedorf"
replace gdenameopt1 = "busswil (tg)" if gdename == "busswil"
replace gdenameopt1 = "beinwil (freiamt)" if gdename == "beinwil"
replace gdenameopt2 = "busswil bei bueren" if gdename == "busswil"
replace gdenameopt2 = "beinwil (so)" if gdename == "beinwil"
replace gdenameopt3 = "busswil bei melchnau" if gdename == "busswil"
replace gdenameopt3 = "beinwil am see" if gdename == "beinwil"
replace gdenameopt1 = "hofstetten (so)" if gdename == "hofstetten"
replace gdenameopt1 = "buchs (ag)" if gdename == "buchs"
replace gdenameopt2 = "hofstetten bei brienz" if gdename == "hofstetten"
replace gdenameopt2 = "buchs (lu)" if gdename == "buchs"
replace gdenameopt3 = "hofstetten bei elgg" if gdename == "hofstetten"
replace gdenameopt3 = "buchs (sg)" if gdename == "buchs"
replace gdenameopt4 = "hofstetten flueh" if gdename == "hofstetten"
replace gdenameopt4 = "buchs (zh)" if gdename == "buchs"

reshape long gdenameopt, i(PID year ID_dupl) j(options_manual)
drop if gdenameopt == ""
replace gdename = gdenameopt 
drop gdenameopt

save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo90b.dta", replace
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo90a.dta"
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo90.dta", replace

use "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo90.dta", clear


// merge on PLZ4
sort PLZ4
merge m:1 PLZ4 using "$path\02_Processed_data\10_Directors_1934_2003\uPLZ_Geo.dta"
drop if _merge == 2

preserve
keep if _merge == 3
replace geo = 91
gen gdename_used = gdename
drop _merge
keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
countPLZ4 GdeNr_E_CNTR GdeNr_N_CNTR gdenr_2018 gdename_used geo options_manual
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo91.dta", replace
restore

// merge on gdename & year for remaining non-merges
drop if _merge == 3
sort gdename year 
drop _merge countPLZ4 GdeNr_E_CNTR GdeNr_N_CNTR gdenr_2018
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2

preserve
keep if _merge == 3
replace geo = 92
rename E_CNTR GdeNr_E_CNTR 
rename N_CNTR GdeNr_N_CNTR 
gen gdename_used = gdename
drop _merge
keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
GdeNr_E_CNTR GdeNr_N_CNTR gdenr_2018 gdename_used geo options_manual
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo92.dta", replace
restore

// merge on gdename without year (!) for remaining non-merges
drop if _merge == 3
drop _merge E_CNTR N_CNTR gdenr_2018
sort gdename
merge m:1 gdename using "$path\02_Processed_data\08_Municipalities\unique_Gdename_Gdenr2018-Geo.dta"
drop if _merge == 2

preserve
keep if _merge == 3
replace geo = 93
rename E_CNTR GdeNr_E_CNTR 
rename N_CNTR GdeNr_N_CNTR 
gen gdename_used = gdename
drop _merge
keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
GdeNr_E_CNTR GdeNr_N_CNTR gdenr_2018 gdename_used geo options_manual
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo93.dta", replace
restore

keep if _merge == 1
drop _merge E_CNTR N_CNTR

// foreign cities (from hand coding)

gen nCH = .
replace nCH = 1 if PLZ4 == 9490 & gdename == "liechtenstein"
replace nCH = 1 if PLZ4 == 9490 & gdename == "vadu fl"
replace nCH = 1 if PLZ4 == 9490 & gdename == "vaduz"
replace nCH = 1 if PLZ4 == 9490 & gdename == "vaduz (fl"
replace nCH = 1 if PLZ4 == 9490 & gdename == "vaduz (fl)"
replace nCH = 1 if PLZ4 == 9490 & gdename == "vaduz ."
replace nCH = 1 if PLZ4 == 9490 & gdename == "vaduz fl"
replace nCH = 1 if PLZ4 == 9490 & gdename == "vaduz. parc des sports chateau d oex sa"
replace nCH = 1 if PLZ4 == 9490 & gdename == "vaduz/fl"
replace nCH = 1 if PLZ4 == 9494 & gdename == "planken"
replace nCH = 1 if PLZ4 == 9494 & gdename == "planken fl"
replace nCH = 1 if PLZ4 == 9494 & gdename == "schaan"
replace nCH = 1 if PLZ4 == 9494 & gdename == "schaan (fl)"
replace nCH = 1 if PLZ4 == 9494 & gdename == "schaan fl"
replace nCH = 1 if PLZ4 == 9494 & gdename == "schaan licht"
replace nCH = 1 if PLZ4 == 9494 & gdename == "schaan/fl"
replace nCH = 1 if PLZ4 == 9494 & gdename == "schean"
replace nCH = 1 if PLZ4 == 9494 & gdename == "scheer) (fl)"
replace nCH = 1 if PLZ4 == 9494 & gdename == "schoen"
replace nCH = 1 if PLZ4 == 9494 & gdename == "schoen fl"
replace nCH = 1 if PLZ4 == 9494 & gdename == "schwan (fl)"
replace nCH = 1 if PLZ4 == 9494 & gdename == "schwan fl"
replace nCH = 1 if PLZ4 == 9494 & gdename == "sclgn fl"
replace nCH = 1 if PLZ4 == 9494 & gdename == "vaduz"
replace nCH = 1 if PLZ4 == 9495 & gdename == "teesen"
replace nCH = 1 if PLZ4 == 9495 & gdename == "thesen fl"
replace nCH = 1 if PLZ4 == 9495 & gdename == "triesen"
replace nCH = 1 if PLZ4 == 9495 & gdename == "triesen (fl)"
replace nCH = 1 if PLZ4 == 9495 & gdename == "triesen fl"
replace nCH = 1 if PLZ4 == 9495 & gdename == "triesen/fi"
replace nCH = 1 if PLZ4 == 9495 & gdename == "triesen/fl"
replace nCH = 1 if PLZ4 == 9495 & gdename == "triessn"
replace nCH = 1 if PLZ4 == 9495 & gdename == "trilesen"
replace nCH = 1 if PLZ4 == 9495 & gdename == "trlesen"
replace nCH = 1 if PLZ4 == 9495 & gdename == "trlesen (fl)"
replace nCH = 1 if PLZ4 == 9495 & gdename == "trlesen fl"
replace nCH = 1 if PLZ4 == 9495 & gdename == "trtesen fl"
replace nCH = 1 if PLZ4 == 9496 & gdename == "bafzers fl"
replace nCH = 1 if PLZ4 == 9496 & gdename == "baizers"
replace nCH = 1 if PLZ4 == 9496 & gdename == "baizers (fl)"
replace nCH = 1 if PLZ4 == 9496 & gdename == "baizers fl"
replace nCH = 1 if PLZ4 == 9496 & gdename == "balms"
replace nCH = 1 if PLZ4 == 9496 & gdename == "balzens"
replace nCH = 1 if PLZ4 == 9496 & gdename == "balzens fl"
replace nCH = 1 if PLZ4 == 9496 & gdename == "balzera fl"
replace nCH = 1 if PLZ4 == 9496 & gdename == "balzers"
replace nCH = 1 if PLZ4 == 9496 & gdename == "balzers (fl)"
replace nCH = 1 if PLZ4 == 9496 & gdename == "balzers fl"
replace nCH = 1 if PLZ4 == 9496 & gdename == "triesen"
replace nCH = 1 if PLZ4 == 9496 & gdename == "trleaen"
replace nCH = 1 if PLZ4 == 3941 & gdename == "varan"
replace nCH = 1 if PLZ4 == 6911 & gdename == "campione"
replace nCH = 1 if PLZ4 == 6911 & gdename == "campione d italia"
replace nCH = 1 if PLZ4 == 6911 & gdename == "campione d italia (i)"
replace nCH = 1 if PLZ4 == 6911 & gdename == "camplone"
replace nCH = 1 if PLZ4 == 6911 & gdename == "camplone d italla"
replace nCH = 1 if PLZ4 == 6849 & gdename == "sons"
replace nCH = 1 if PLZ4 == 9491 & gdename == "bendern"
replace nCH = 1 if PLZ4 == 9491 & gdename == "bendern (fl)"
replace nCH = 1 if PLZ4 == 9491 & gdename == "bendern fl"
replace nCH = 1 if PLZ4 == 9491 & gdename == "bendern/fl"
replace nCH = 1 if PLZ4 == 9491 & gdename == "gampf rein"
replace nCH = 1 if PLZ4 == 9491 & gdename == "gampfrein"
replace nCH = 1 if PLZ4 == 9491 & gdename == "gamprin"
replace nCH = 1 if PLZ4 == 9491 & gdename == "gamprin bendern (fl)"
replace nCH = 1 if PLZ4 == 9491 & gdename == "gamprin fl"
replace nCH = 1 if PLZ4 == 9491 & gdename == "nendeln (fl)"
replace nCH = 1 if PLZ4 == 9491 & gdename == "ruggeil (fl)"
replace nCH = 1 if PLZ4 == 9491 & gdename == "ruggell"
replace nCH = 1 if PLZ4 == 9491 & gdename == "ruggell (fl)"
replace nCH = 1 if PLZ4 == 9491 & gdename == "ruggell fl"
replace nCH = 1 if PLZ4 == 9491 & gdename == "ruggell/fl"
replace nCH = 1 if PLZ4 == 9491 & gdename == "schaanwaid (fl)"
replace nCH = 1 if PLZ4 == 9491 & gdename == "schaanwald"
replace nCH = 1 if PLZ4 == 9491 & gdename == "schaanwald (fl)"
replace nCH = 1 if PLZ4 == 9491 & gdename == "schaanwald fl"
replace nCH = 1 if PLZ4 == 9491 & gdename == "schellenberg"
replace nCH = 1 if PLZ4 == 9491 & gdename == "schellenberg (fl)"
replace nCH = 1 if PLZ4 == 9491 & gdename == "schellenberg/fl"
replace nCH = 1 if PLZ4 == 9492 & gdename == "eschen"
replace nCH = 1 if PLZ4 == 9492 & gdename == "eschen (fl)"
replace nCH = 1 if PLZ4 == 9492 & gdename == "eschen fl"
replace nCH = 1 if PLZ4 == 9492 & gdename == "eschen/fl"
replace nCH = 1 if PLZ4 == 4249 & gdename == "hingen"
replace nCH = 1 if PLZ4 == 9493 & gdename == "mauren fl"
replace nCH = 1 if PLZ4 == 9493 & gdename == "schaanwald"
replace nCH = 1 if PLZ4 == 1249 & gdename == "les cocandes"
replace nCH = 1 if PLZ4 == 4857 & gdename == "aiken"
replace nCH = 1 if PLZ4 == 9497 & gdename == "triesen (fl)"
replace nCH = 1 if PLZ4 == 9497 & gdename == "triesenberg"
replace nCH = 1 if PLZ4 == 9497 & gdename == "triesenberg (fl)"
replace nCH = 1 if PLZ4 == 9497 & gdename == "triesenberg fl"
replace nCH = 1 if PLZ4 == 9497 & gdename == "triesenberg/fl"
replace nCH = 1 if PLZ4 == 9497 & gdename == "trleeenberg fl"
replace nCH = 1 if PLZ4 == 9497 & gdename == "trlesenberg fl"
replace nCH = 1 if PLZ4 == 7131 & gdename == "morlesen"
replace nCH = 1 if PLZ4 == 1210 & gdename == "ferney"
replace nCH = 1 if PLZ4 == 9307 & gdename == "freihof"
replace nCH = 1 if PLZ4 == 8439 & gdename == "boston usa"
replace nCH = 1 if PLZ4 == 4354 & gdename == "datteln/ deutschland"
replace nCH = 1 if PLZ4 == 4354 & gdename == "datteln/deutschland"
replace nCH = 1 if PLZ4 == 9488 & gdename == "schellenberg"
replace nCH = 1 if PLZ4 == 9488 & gdename == "schellenberg (fl)"
replace nCH = 1 if PLZ4 == 9488 & gdename == "schellenberg fl"
replace nCH = 1 if PLZ4 == 9487 & gdename == "gamprin"
replace nCH = 1 if PLZ4 == 9487 & gdename == "gamprin (fl)"
replace nCH = 1 if PLZ4 == 9487 & gdename == "gamprin bendem"
replace nCH = 1 if PLZ4 == 9487 & gdename == "gamprin bendern"
replace nCH = 1 if PLZ4 == 9487 & gdename == "gamprin fl"
replace nCH = 1 if PLZ4 == 9487 & gdename == "garnprin fl"
replace nCH = 1 if PLZ4 == 9487 & gdename == "gemprin"
replace nCH = 1 if PLZ4 == 8238 & gdename == "buesingen (d)"
replace nCH = 1 if PLZ4 == 8238 & gdename == "busingen (d)"
replace nCH = 1 if PLZ4 == 9485 & gdename == "nendein (fl)"
replace nCH = 1 if PLZ4 == 9485 & gdename == "nendeln"
replace nCH = 1 if PLZ4 == 9485 & gdename == "nendeln (fl)"
replace nCH = 1 if PLZ4 == 1655 & gdename == "la villette fr"
replace nCH = 1 if PLZ4 == 9486 & gdename == "schaanwald"
replace nCH = 1 if PLZ4 == 3 & gdename == "les charbonnieres"
replace nCH = 1 if PLZ4 == 6099 & gdename == "buergten"
replace nCH = 1 if PLZ4 == 1 & gdename == ""
replace nCH = 1 if PLZ4 == 1 & gdename == "conches"
replace nCH = 1 if PLZ4 == 1 & gdename == "vr"
replace nCH = 1 if PLZ4 == 5254 & gdename == "when"
replace nCH = 1 if PLZ4 == 5254 & gdename == "wizen"
replace nCH = 1 if PLZ4 == 9498 & gdename == "planken"
replace nCH = 1 if PLZ4 == 1067 & gdename == "boll"
replace nCH = 1 if PLZ4 == 7 & gdename == ""
replace nCH = 1 if PLZ4 == 7 & gdename == "boulevard de suisse"
replace nCH = 1 if PLZ4 == 2490 & gdename == "vaduz (fl)"
replace nCH = 1 if PLZ4 == 9482 & gdename == "eschen (fl)"
replace nCH = 1 if PLZ4 == 4 & gdename == ""
replace nCH = 1 if PLZ4 == 19 & gdename == "av. mervelet"
replace nCH = 1 if PLZ4 == 58040 & gdename == "punta ala (gr) toscana i"
replace nCH = 1 if PLZ4 == 58040 & gdename == "punta ala (gr) toscana italien"
replace nCH = 1 if PLZ4 == 0 & gdename == ""
replace nCH = 1 if PLZ4 == 2 & gdename == ""
replace nCH = 1 if PLZ4 == 12 & gdename == "bis rte. de veyrier"
replace nCH = 1 if PLZ4 == 24 & gdename == "les acacias"
replace nCH = 1 if PLZ4 == 8490 & gdename == "vaduz"
replace nCH = 1 if PLZ4 == 8981 & gdename == "sense"
replace nCH = 1 if PLZ4 == 9440 & gdename == "vaduz"
replace nCH = 1 if PLZ4 == 41844 & gdename == "wegberg d"
replace nCH = 1 if PLZ4 == 8 & gdename == "muenchen d"
replace nCH = 1 if PLZ4 == 9 & gdename == ""
replace nCH = 1 if PLZ4 == 13 & gdename == "palettes"
replace nCH = 1 if PLZ4 == 23 & gdename == ""
replace nCH = 1 if PLZ4 == 41 & gdename == ""
replace nCH = 1 if PLZ4 == 67 & gdename == "dino"
replace nCH = 1 if PLZ4 == 74 & gdename == "victoria street"
replace nCH = 1 if PLZ4 == 485 & gdename == "mast"
replace nCH = 1 if PLZ4 == 757 & gdename == "third avenue 21 st"
replace nCH = 1 if PLZ4 == 1790 & gdename == "ariborgo"
replace nCH = 1 if PLZ4 == 5340 & gdename == "saar"
replace nCH = 1 if PLZ4 == 5487 & gdename == "halbach"
replace nCH = 1 if PLZ4 == 8504 & gdename == "goettingen"
replace nCH = 1 if PLZ4 == 17025 & gdename == "loano (i)"
replace nCH = 1 if PLZ4 == 68220 & gdename == "buschwiller f"
replace nCH = 1 if PLZ4 == 99490 & gdename == "vaduz"

replace nCH = 1 if gdename == "aiken"
replace nCH = 1 if gdename == "ariborgo"
replace nCH = 1 if gdename == "av. mervelet"
replace nCH = 1 if gdename == "bafzers fl"
replace nCH = 1 if gdename == "baizers"
replace nCH = 1 if gdename == "baizers (fl)"
replace nCH = 1 if gdename == "baizers fl"
replace nCH = 1 if gdename == "balms"
replace nCH = 1 if gdename == "balzens"
replace nCH = 1 if gdename == "balzens fl"
replace nCH = 1 if gdename == "balzera fl"
replace nCH = 1 if gdename == "balzers"
replace nCH = 1 if gdename == "balzers (fl)"
replace nCH = 1 if gdename == "balzers fl"
replace nCH = 1 if gdename == "bendern"
replace nCH = 1 if gdename == "bendern (fl)"
replace nCH = 1 if gdename == "bendern fl"
replace nCH = 1 if gdename == "bendern/fl"
replace nCH = 1 if gdename == "boll"
replace nCH = 1 if gdename == "boston usa"
replace nCH = 1 if gdename == "boulevard de suisse"
replace nCH = 1 if gdename == "buergten"
replace nCH = 1 if gdename == "buesingen (d)"
replace nCH = 1 if gdename == "buschwiller f"
replace nCH = 1 if gdename == "busingen (d)"
replace nCH = 1 if gdename == "campione"
replace nCH = 1 if gdename == "campione d italia"
replace nCH = 1 if gdename == "campione d italia (i)"
replace nCH = 1 if gdename == "camplone"
replace nCH = 1 if gdename == "camplone d italla"
replace nCH = 1 if gdename == "chemin de l ancienne"
replace nCH = 1 if gdename == "conches"
replace nCH = 1 if gdename == "datteln/ deutschland"
replace nCH = 1 if gdename == "datteln/deutschland"
replace nCH = 1 if gdename == "dino"
replace nCH = 1 if gdename == "eschen"
replace nCH = 1 if gdename == "eschen (fl)"
replace nCH = 1 if gdename == "eschen fl"
replace nCH = 1 if gdename == "eschen/fl"
replace nCH = 1 if gdename == "ferney"
replace nCH = 1 if gdename == "freihof"
replace nCH = 1 if gdename == "gampf rein"
replace nCH = 1 if gdename == "gampfrein"
replace nCH = 1 if gdename == "gamprin"
replace nCH = 1 if gdename == "gamprin (fl)"
replace nCH = 1 if gdename == "gamprin bendem"
replace nCH = 1 if gdename == "gamprin bendern"
replace nCH = 1 if gdename == "gamprin bendern (fl)"
replace nCH = 1 if gdename == "gamprin fl"
replace nCH = 1 if gdename == "garnprin fl"
replace nCH = 1 if gdename == "gemprin"
replace nCH = 1 if gdename == "goettingen"
replace nCH = 1 if gdename == "hagnau"
replace nCH = 1 if gdename == "halbach"
replace nCH = 1 if gdename == "hingen"
replace nCH = 1 if gdename == "la villette fr"
replace nCH = 1 if gdename == "les acacias"
replace nCH = 1 if gdename == "les charbonnieres"
replace nCH = 1 if gdename == "les cocandes"
replace nCH = 1 if gdename == "letendorf"
replace nCH = 1 if gdename == "liechtenstein"
replace nCH = 1 if gdename == "loano (i)"
replace nCH = 1 if gdename == "mast"
replace nCH = 1 if gdename == "mauren fl"
replace nCH = 1 if gdename == "morlesen"
replace nCH = 1 if gdename == "muenchen d"
replace nCH = 1 if gdename == "milano"
replace nCH = 1 if gdename == "nendein (fl)"
replace nCH = 1 if gdename == "nendeln"
replace nCH = 1 if gdename == "nendeln (fl)"
replace nCH = 1 if gdename == "palettes"
replace nCH = 1 if gdename == "pasoux"
replace nCH = 1 if gdename == "planken"
replace nCH = 1 if gdename == "planken fl"
replace nCH = 1 if gdename == "ruggeil (fl)"
replace nCH = 1 if gdename == "ruggell"
replace nCH = 1 if gdename == "ruggell (fl)"
replace nCH = 1 if gdename == "ruggell fl"
replace nCH = 1 if gdename == "ruggell/fl"
replace nCH = 1 if gdename == "saar"
replace nCH = 1 if gdename == "schaan"
replace nCH = 1 if gdename == "schaan (fl)"
replace nCH = 1 if gdename == "schaan fl"
replace nCH = 1 if gdename == "schaan licht"
replace nCH = 1 if gdename == "schaan/fl"
replace nCH = 1 if gdename == "schaanwaid (fl)"
replace nCH = 1 if gdename == "schaanwald"
replace nCH = 1 if gdename == "schaanwald (fl)"
replace nCH = 1 if gdename == "schaanwald fl"
replace nCH = 1 if gdename == "schean"
replace nCH = 1 if gdename == "scheer) (fl)"
replace nCH = 1 if gdename == "schellenberg"
replace nCH = 1 if gdename == "schellenberg (fl)"
replace nCH = 1 if gdename == "schellenberg fl"
replace nCH = 1 if gdename == "schellenberg/fl"
replace nCH = 1 if gdename == "schoen"
replace nCH = 1 if gdename == "schoen fl"
replace nCH = 1 if gdename == "schwan (fl)"
replace nCH = 1 if gdename == "schwan fl"
replace nCH = 1 if gdename == "sclgn fl"
replace nCH = 1 if gdename == "sense"
replace nCH = 1 if gdename == "teesen"
replace nCH = 1 if gdename == "thesen fl"
replace nCH = 1 if gdename == "triesen"
replace nCH = 1 if gdename == "triesen (fl)"
replace nCH = 1 if gdename == "triesen fl"
replace nCH = 1 if gdename == "triesen/fi"
replace nCH = 1 if gdename == "triesen/fl"
replace nCH = 1 if gdename == "triesenberg"
replace nCH = 1 if gdename == "triesenberg (fl)"
replace nCH = 1 if gdename == "triesenberg fl"
replace nCH = 1 if gdename == "triesenberg/fl"
replace nCH = 1 if gdename == "triessn"
replace nCH = 1 if gdename == "trilesen"
replace nCH = 1 if gdename == "trleaen"
replace nCH = 1 if gdename == "trleeenberg fl"
replace nCH = 1 if gdename == "trlesen"
replace nCH = 1 if gdename == "trlesen (fl)"
replace nCH = 1 if gdename == "trlesen fl"
replace nCH = 1 if gdename == "trlesenberg fl"
replace nCH = 1 if gdename == "trtesen fl"
replace nCH = 1 if gdename == "vadu fl"
replace nCH = 1 if gdename == "vaduz"
replace nCH = 1 if gdename == "vaduz (fl"
replace nCH = 1 if gdename == "vaduz (fl)"
replace nCH = 1 if gdename == "vaduz ."
replace nCH = 1 if gdename == "vaduz fl"
replace nCH = 1 if gdename == "vaduz/fl"
replace nCH = 1 if gdename == "varan"
replace nCH = 1 if gdename == "victoria street"
replace nCH = 1 if gdename == "wegberg d"
replace nCH = 1 if gdename == "when"
replace nCH = 1 if gdename == "wile"
replace nCH = 1 if gdename == "wizen"

replace nCH = 1 if gdename == "accra"
replace nCH = 1 if gdename == "ada"
replace nCH = 1 if gdename == "adelaide"
replace nCH = 1 if gdename == "akron"
replace nCH = 1 if gdename == "alexandria"
replace nCH = 1 if gdename == "alhambra"
replace nCH = 1 if gdename == "alpena"
replace nCH = 1 if gdename == "amman"
replace nCH = 1 if gdename == "amsterdam"
replace nCH = 1 if gdename == "ancona"
replace nCH = 1 if gdename == "andorra"
replace nCH = 1 if gdename == "ankara"
replace nCH = 1 if gdename == "annecy"
replace nCH = 1 if gdename == "antwerpen"
replace nCH = 1 if gdename == "aosta"
replace nCH = 1 if gdename == "arcadia"
replace nCH = 1 if gdename == "arcata"
replace nCH = 1 if gdename == "ardmore"
replace nCH = 1 if gdename == "arnhem"
replace nCH = 1 if gdename == "arnold"
replace nCH = 1 if gdename == "arras"
replace nCH = 1 if gdename == "athens"
replace nCH = 1 if gdename == "atherton"
replace nCH = 1 if gdename == "atlanta"
replace nCH = 1 if gdename == "atlantic city"
replace nCH = 1 if gdename == "augsburg"
replace nCH = 1 if gdename == "aurora"
replace nCH = 1 if gdename == "baltimore"
replace nCH = 1 if gdename == "bangkok"
replace nCH = 1 if gdename == "barcelona"
replace nCH = 1 if gdename == "bari"
replace nCH = 1 if gdename == "barranquilla"
replace nCH = 1 if gdename == "barrington"
replace nCH = 1 if gdename == "batavia"
replace nCH = 1 if gdename == "beaumont"
replace nCH = 1 if gdename == "beirut"
replace nCH = 1 if gdename == "belleville"
replace nCH = 1 if gdename == "belmont"
replace nCH = 1 if gdename == "benton harbor"
replace nCH = 1 if gdename == "bergamo"
replace nCH = 1 if gdename == "bergen"
replace nCH = 1 if gdename == "berkeley"
replace nCH = 1 if gdename == "berlin"
replace nCH = 1 if gdename == "berthoud"
replace nCH = 1 if gdename == "besancon"
replace nCH = 1 if gdename == "beverly hills"
replace nCH = 1 if gdename == "bielefeld"
replace nCH = 1 if gdename == "birmingham"
replace nCH = 1 if gdename == "bled"
replace nCH = 1 if gdename == "bloomfield"
replace nCH = 1 if gdename == "bologna"
replace nCH = 1 if gdename == "bolzano"
replace nCH = 1 if gdename == "bonn"
replace nCH = 1 if gdename == "bordeaux"
replace nCH = 1 if gdename == "boston"
replace nCH = 1 if gdename == "bound brook"
replace nCH = 1 if gdename == "bournemouth"
replace nCH = 1 if gdename == "bracknell"
replace nCH = 1 if gdename == "bradford"
replace nCH = 1 if gdename == "braunschweig"
replace nCH = 1 if gdename == "bregenz"
replace nCH = 1 if gdename == "bremen"
replace nCH = 1 if gdename == "brentwood"
replace nCH = 1 if gdename == "briarcliff manor"
replace nCH = 1 if gdename == "bridgeport"
replace nCH = 1 if gdename == "bromley"
replace nCH = 1 if gdename == "bronxville"
replace nCH = 1 if gdename == "brookfield"
replace nCH = 1 if gdename == "brooklyn"
replace nCH = 1 if gdename == "budapest"
replace nCH = 1 if gdename == "buenos aires"
replace nCH = 1 if gdename == "burlingame"
replace nCH = 1 if gdename == "butler"
replace nCH = 1 if gdename == "cagliari"
replace nCH = 1 if gdename == "cairo"
replace nCH = 1 if gdename == "cambridge"
replace nCH = 1 if gdename == "camden"
replace nCH = 1 if gdename == "canton"
replace nCH = 1 if gdename == "cape town"
replace nCH = 1 if gdename == "caracas"
replace nCH = 1 if gdename == "carmel"
replace nCH = 1 if gdename == "casablanca"
replace nCH = 1 if gdename == "cedarhurst"
replace nCH = 1 if gdename == "centreville"
replace nCH = 1 if gdename == "chatham"
replace nCH = 1 if gdename == "chester"
replace nCH = 1 if gdename == "chicago"
replace nCH = 1 if gdename == "chillicothe"
replace nCH = 1 if gdename == "cincinnati"
replace nCH = 1 if gdename == "clairton"
replace nCH = 1 if gdename == "clearwater"
replace nCH = 1 if gdename == "clermont ferrand"
replace nCH = 1 if gdename == "cleveland"
replace nCH = 1 if gdename == "cleveland heights"
replace nCH = 1 if gdename == "coburg"
replace nCH = 1 if gdename == "cody"
replace nCH = 1 if gdename == "cold spring"
replace nCH = 1 if gdename == "cologne"
replace nCH = 1 if gdename == "columbus"
replace nCH = 1 if gdename == "como"
replace nCH = 1 if gdename == "concord"
replace nCH = 1 if gdename == "corning"
replace nCH = 1 if gdename == "cuyahoga falls"
replace nCH = 1 if gdename == "dallas"
replace nCH = 1 if gdename == "danbury"
replace nCH = 1 if gdename == "danville"
replace nCH = 1 if gdename == "darien"
replace nCH = 1 if gdename == "dayton"
replace nCH = 1 if gdename == "decatur"
replace nCH = 1 if gdename == "dekalb"
replace nCH = 1 if gdename == "demarest"
replace nCH = 1 if gdename == "denver"
replace nCH = 1 if gdename == "detroit"
replace nCH = 1 if gdename == "dieppe"
replace nCH = 1 if gdename == "dijon"
replace nCH = 1 if gdename == "dortmund"
replace nCH = 1 if gdename == "dresden"
replace nCH = 1 if gdename == "dublin"
replace nCH = 1 if gdename == "duisburg"
replace nCH = 1 if gdename == "dunbar"
replace nCH = 1 if gdename == "durban"
replace nCH = 1 if gdename == "dusseldorf"
replace nCH = 1 if gdename == "east aurora"
replace nCH = 1 if gdename == "easton"
replace nCH = 1 if gdename == "edgerton"
replace nCH = 1 if gdename == "edinburgh"
replace nCH = 1 if gdename == "eindhoven"
replace nCH = 1 if gdename == "el monte"
replace nCH = 1 if gdename == "eldorado"
replace nCH = 1 if gdename == "elizabeth"
replace nCH = 1 if gdename == "elkhart"
replace nCH = 1 if gdename == "elm grove"
replace nCH = 1 if gdename == "enfield"
replace nCH = 1 if gdename == "esbjerg"
replace nCH = 1 if gdename == "eschen"
replace nCH = 1 if gdename == "esperanza"
replace nCH = 1 if gdename == "essen"
replace nCH = 1 if gdename == "evanston"
replace nCH = 1 if gdename == "evansville"
replace nCH = 1 if gdename == "fairfield"
replace nCH = 1 if gdename == "fairport"
replace nCH = 1 if gdename == "findlay"
replace nCH = 1 if gdename == "fitchburg"
replace nCH = 1 if gdename == "florida"
replace nCH = 1 if gdename == "fort wayne"
replace nCH = 1 if gdename == "frankfurt"
replace nCH = 1 if gdename == "fraser"
replace nCH = 1 if gdename == "freiburg"
replace nCH = 1 if gdename == "fremont"
replace nCH = 1 if gdename == "fresno"
replace nCH = 1 if gdename == "fullerton"
replace nCH = 1 if gdename == "furth"
replace nCH = 1 if gdename == "geelong"
replace nCH = 1 if gdename == "geneva"
replace nCH = 1 if gdename == "genoa"
replace nCH = 1 if gdename == "giessen"
replace nCH = 1 if gdename == "glasgow"
replace nCH = 1 if gdename == "glen ellyn"
replace nCH = 1 if gdename == "glendale"
replace nCH = 1 if gdename == "glendora"
replace nCH = 1 if gdename == "gottingen"
replace nCH = 1 if gdename == "grand rapids"
replace nCH = 1 if gdename == "graz"
replace nCH = 1 if gdename == "great neck"
replace nCH = 1 if gdename == "greensboro"
replace nCH = 1 if gdename == "greenville"
replace nCH = 1 if gdename == "grenoble"
replace nCH = 1 if gdename == "grosse pointe woods"
replace nCH = 1 if gdename == "haifa"
replace nCH = 1 if gdename == "halifax"
replace nCH = 1 if gdename == "hamburg"
replace nCH = 1 if gdename == "hamilton"
replace nCH = 1 if gdename == "hannover"
replace nCH = 1 if gdename == "hanover"
replace nCH = 1 if gdename == "harrison"
replace nCH = 1 if gdename == "hartford"
replace nCH = 1 if gdename == "heidelberg"
replace nCH = 1 if gdename == "helsingborg"
replace nCH = 1 if gdename == "helsinki"
replace nCH = 1 if gdename == "hendersonville"
replace nCH = 1 if gdename == "herndon"
replace nCH = 1 if gdename == "herrin"
replace nCH = 1 if gdename == "hesperia"
replace nCH = 1 if gdename == "highland park"
replace nCH = 1 if gdename == "hillside"
replace nCH = 1 if gdename == "hof"
replace nCH = 1 if gdename == "hollywood"
replace nCH = 1 if gdename == "homewood"
replace nCH = 1 if gdename == "hong kong"
replace nCH = 1 if gdename == "houston"
replace nCH = 1 if gdename == "hove"
replace nCH = 1 if gdename == "huntingdon"
replace nCH = 1 if gdename == "huntington"
replace nCH = 1 if gdename == "ilford"
replace nCH = 1 if gdename == "indianapolis"
replace nCH = 1 if gdename == "inglewood"
replace nCH = 1 if gdename == "innsbruck"
replace nCH = 1 if gdename == "istanbul"
replace nCH = 1 if gdename == "jackson"
replace nCH = 1 if gdename == "jerusalem"
replace nCH = 1 if gdename == "johannesburg"
replace nCH = 1 if gdename == "kabul"
replace nCH = 1 if gdename == "kalmar"
replace nCH = 1 if gdename == "kanpur"
replace nCH = 1 if gdename == "kansas city"
replace nCH = 1 if gdename == "karachi"
replace nCH = 1 if gdename == "karlsruhe"
replace nCH = 1 if gdename == "kassel"
replace nCH = 1 if gdename == "kenilworth"
replace nCH = 1 if gdename == "kensington"
replace nCH = 1 if gdename == "khartoum"
replace nCH = 1 if gdename == "kiel"
replace nCH = 1 if gdename == "kingsport"
replace nCH = 1 if gdename == "kingston"
replace nCH = 1 if gdename == "kobe"
replace nCH = 1 if gdename == "kobenhavn"
replace nCH = 1 if gdename == "kokomo"
replace nCH = 1 if gdename == "kumasi"
replace nCH = 1 if gdename == "la crosse"
replace nCH = 1 if gdename == "la paz"
replace nCH = 1 if gdename == "lagos"
replace nCH = 1 if gdename == "lake forest"
replace nCH = 1 if gdename == "lancaster"
replace nCH = 1 if gdename == "larchmont"
replace nCH = 1 if gdename == "latrobe"
replace nCH = 1 if gdename == "le havre"
replace nCH = 1 if gdename == "leipzig"
replace nCH = 1 if gdename == "leominster"
replace nCH = 1 if gdename == "lewiston"
replace nCH = 1 if gdename == "libertyville"
replace nCH = 1 if gdename == "liege"
replace nCH = 1 if gdename == "lima"
replace nCH = 1 if gdename == "limoges"
replace nCH = 1 if gdename == "linz"
replace nCH = 1 if gdename == "lisbon"
replace nCH = 1 if gdename == "littleton"
replace nCH = 1 if gdename == "liverpool"
replace nCH = 1 if gdename == "livingston"
replace nCH = 1 if gdename == "livorno"
replace nCH = 1 if gdename == "london"
replace nCH = 1 if gdename == "long beach"
replace nCH = 1 if gdename == "los altos"
replace nCH = 1 if gdename == "los angeles"
replace nCH = 1 if gdename == "louisville"
replace nCH = 1 if gdename == "luanda"
replace nCH = 1 if gdename == "lubeck"
replace nCH = 1 if gdename == "lulea"
replace nCH = 1 if gdename == "luxembourg"
replace nCH = 1 if gdename == "lynn"
replace nCH = 1 if gdename == "lyon"
replace nCH = 1 if gdename == "madison"
replace nCH = 1 if gdename == "madrid"
replace nCH = 1 if gdename == "maidenhead"
replace nCH = 1 if gdename == "mainz"
replace nCH = 1 if gdename == "mamaroneck"
replace nCH = 1 if gdename == "managua"
replace nCH = 1 if gdename == "manchester"
replace nCH = 1 if gdename == "manila"
replace nCH = 1 if gdename == "mannheim"
replace nCH = 1 if gdename == "maplewood"
replace nCH = 1 if gdename == "marseille"
replace nCH = 1 if gdename == "medan"
replace nCH = 1 if gdename == "meknes"
replace nCH = 1 if gdename == "melbourne"
replace nCH = 1 if gdename == "melrose park"
replace nCH = 1 if gdename == "menlo park"
replace nCH = 1 if gdename == "messina"
replace nCH = 1 if gdename == "mexico"
replace nCH = 1 if gdename == "mexico city"
replace nCH = 1 if gdename == "michigan city"
replace nCH = 1 if gdename == "middletown"
replace nCH = 1 if gdename == "midland"
replace nCH = 1 if gdename == "milan"
replace nCH = 1 if gdename == "milford"
replace nCH = 1 if gdename == "milton"
replace nCH = 1 if gdename == "milwaukee"
replace nCH = 1 if gdename == "minden"
replace nCH = 1 if gdename == "minneapolis"
replace nCH = 1 if gdename == "mishawaka"
replace nCH = 1 if gdename == "mission"
replace nCH = 1 if gdename == "modena"
replace nCH = 1 if gdename == "moline"
replace nCH = 1 if gdename == "mombasa"
replace nCH = 1 if gdename == "monaco"
replace nCH = 1 if gdename == "monroe"
replace nCH = 1 if gdename == "monrovia"
replace nCH = 1 if gdename == "montclair"
replace nCH = 1 if gdename == "montevideo"
replace nCH = 1 if gdename == "montreal"
replace nCH = 1 if gdename == "morris plains"
replace nCH = 1 if gdename == "morristown"
replace nCH = 1 if gdename == "mount kisco"
replace nCH = 1 if gdename == "mount vernon"
replace nCH = 1 if gdename == "mulhouse"
replace nCH = 1 if gdename == "munich"
replace nCH = 1 if gdename == "muscatine"
replace nCH = 1 if gdename == "muskegon"
replace nCH = 1 if gdename == "nairobi"
replace nCH = 1 if gdename == "nancy"
replace nCH = 1 if gdename == "napoleon"
replace nCH = 1 if gdename == "nashville"
replace nCH = 1 if gdename == "nassau"
replace nCH = 1 if gdename == "new bedford"
replace nCH = 1 if gdename == "new britain"
replace nCH = 1 if gdename == "new brunswick"
replace nCH = 1 if gdename == "new delhi"
replace nCH = 1 if gdename == "new haven"
replace nCH = 1 if gdename == "new hope"
replace nCH = 1 if gdename == "new orleans"
replace nCH = 1 if gdename == "new rochelle"
replace nCH = 1 if gdename == "new york"
replace nCH = 1 if gdename == "newark"
replace nCH = 1 if gdename == "newburgh"
replace nCH = 1 if gdename == "newbury"
replace nCH = 1 if gdename == "newport"
replace nCH = 1 if gdename == "niagara falls"
replace nCH = 1 if gdename == "nice"
replace nCH = 1 if gdename == "norfolk"
replace nCH = 1 if gdename == "norristown"
replace nCH = 1 if gdename == "norrkoping"
replace nCH = 1 if gdename == "north chicago"
replace nCH = 1 if gdename == "north tonawanda"
replace nCH = 1 if gdename == "northbrook"
replace nCH = 1 if gdename == "northwood"
replace nCH = 1 if gdename == "norwalk"
replace nCH = 1 if gdename == "norwich"
replace nCH = 1 if gdename == "nottingham"
replace nCH = 1 if gdename == "novara"
replace nCH = 1 if gdename == "oakland"
replace nCH = 1 if gdename == "oldenburg"
replace nCH = 1 if gdename == "ontario"
replace nCH = 1 if gdename == "oradell"
replace nCH = 1 if gdename == "oran"
replace nCH = 1 if gdename == "orange"
replace nCH = 1 if gdename == "orinda"
replace nCH = 1 if gdename == "orleans"
replace nCH = 1 if gdename == "oslo"
replace nCH = 1 if gdename == "osnabruck"
replace nCH = 1 if gdename == "ossining"
replace nCH = 1 if gdename == "ottumwa"
replace nCH = 1 if gdename == "palatine"
replace nCH = 1 if gdename == "palermo"
replace nCH = 1 if gdename == "palm springs"
replace nCH = 1 if gdename == "palo alto"
replace nCH = 1 if gdename == "paradise"
replace nCH = 1 if gdename == "paris"
replace nCH = 1 if gdename == "parma"
replace nCH = 1 if gdename == "pasadena"
replace nCH = 1 if gdename == "pawtucket"
replace nCH = 1 if gdename == "peabody"
replace nCH = 1 if gdename == "pelham"
replace nCH = 1 if gdename == "pelham manor"
replace nCH = 1 if gdename == "peoria"
replace nCH = 1 if gdename == "peterborough"
replace nCH = 1 if gdename == "philadelphia"
replace nCH = 1 if gdename == "phoenix"
replace nCH = 1 if gdename == "piedmont"
replace nCH = 1 if gdename == "pittsburg"
replace nCH = 1 if gdename == "pittsburgh"
replace nCH = 1 if gdename == "plymouth"
replace nCH = 1 if gdename == "port chester"
replace nCH = 1 if gdename == "portland"
replace nCH = 1 if gdename == "porto"
replace nCH = 1 if gdename == "powell"
replace nCH = 1 if gdename == "prague"
replace nCH = 1 if gdename == "pretoria"
replace nCH = 1 if gdename == "princeton"
replace nCH = 1 if gdename == "rabat"
replace nCH = 1 if gdename == "racine"
replace nCH = 1 if gdename == "rahway"
replace nCH = 1 if gdename == "ravenna"
replace nCH = 1 if gdename == "reading"
replace nCH = 1 if gdename == "regensburg"
replace nCH = 1 if gdename == "riberalta"
replace nCH = 1 if gdename == "richmond"
replace nCH = 1 if gdename == "ridgefield"
replace nCH = 1 if gdename == "ridgewood"
replace nCH = 1 if gdename == "riga"
replace nCH = 1 if gdename == "rio de janeiro"
replace nCH = 1 if gdename == "river forest"
replace nCH = 1 if gdename == "riverside"
replace nCH = 1 if gdename == "rochester"
replace nCH = 1 if gdename == "rock island"
replace nCH = 1 if gdename == "rockford"
replace nCH = 1 if gdename == "roma"
replace nCH = 1 if gdename == "rome"
replace nCH = 1 if gdename == "roseland"
replace nCH = 1 if gdename == "rotterdam"
replace nCH = 1 if gdename == "rouen"
replace nCH = 1 if gdename == "rye"
replace nCH = 1 if gdename == "saint louis"
replace nCH = 1 if gdename == "salem"
replace nCH = 1 if gdename == "salt lake city"
replace nCH = 1 if gdename == "salzburg"
replace nCH = 1 if gdename == "san bernardino"
replace nCH = 1 if gdename == "san carlos"
replace nCH = 1 if gdename == "san diego"
replace nCH = 1 if gdename == "san francisco"
replace nCH = 1 if gdename == "san sebastian"
replace nCH = 1 if gdename == "santa ana"
replace nCH = 1 if gdename == "santa monica"
replace nCH = 1 if gdename == "santos"
replace nCH = 1 if gdename == "sao paulo"
replace nCH = 1 if gdename == "sarnia"
replace nCH = 1 if gdename == "scarsdale"
replace nCH = 1 if gdename == "schaan"
replace nCH = 1 if gdename == "schiller park"
replace nCH = 1 if gdename == "sea cliff"
replace nCH = 1 if gdename == "seattle"
replace nCH = 1 if gdename == "seoul"
replace nCH = 1 if gdename == "shaker heights"
replace nCH = 1 if gdename == "shanghai"
replace nCH = 1 if gdename == "sheffield"
replace nCH = 1 if gdename == "sidney"
replace nCH = 1 if gdename == "slough"
replace nCH = 1 if gdename == "sofia"
replace nCH = 1 if gdename == "solihull"
replace nCH = 1 if gdename == "south bend"
replace nCH = 1 if gdename == "south pasadena"
replace nCH = 1 if gdename == "sparta"
replace nCH = 1 if gdename == "springfield"
replace nCH = 1 if gdename == "stamford"
replace nCH = 1 if gdename == "stanton"
replace nCH = 1 if gdename == "stockholm"
replace nCH = 1 if gdename == "strasbourg"
replace nCH = 1 if gdename == "strasburg"
replace nCH = 1 if gdename == "stuttgart"
replace nCH = 1 if gdename == "summit"
replace nCH = 1 if gdename == "sundsvall"
replace nCH = 1 if gdename == "sutton"
replace nCH = 1 if gdename == "sylvania"
replace nCH = 1 if gdename == "syracuse"
replace nCH = 1 if gdename == "tampere"
replace nCH = 1 if gdename == "tamworth"
replace nCH = 1 if gdename == "tokyo"
replace nCH = 1 if gdename == "toledo"
replace nCH = 1 if gdename == "toronto"
replace nCH = 1 if gdename == "torrance"
replace nCH = 1 if gdename == "toulouse"
replace nCH = 1 if gdename == "trenton"
replace nCH = 1 if gdename == "triesen"
replace nCH = 1 if gdename == "triesenberg"
replace nCH = 1 if gdename == "trieste"
replace nCH = 1 if gdename == "tripoli"
replace nCH = 1 if gdename == "troy"
replace nCH = 1 if gdename == "troyes"
replace nCH = 1 if gdename == "tulsa"
replace nCH = 1 if gdename == "tunis"
replace nCH = 1 if gdename == "turin"
replace nCH = 1 if gdename == "ulm"
replace nCH = 1 if gdename == "union"
replace nCH = 1 if gdename == "upland"
replace nCH = 1 if gdename == "utrecht"
replace nCH = 1 if gdename == "vaduz"
replace nCH = 1 if gdename == "valencia"
replace nCH = 1 if gdename == "valparaiso"
replace nCH = 1 if gdename == "vancouver"
replace nCH = 1 if gdename == "verona"
replace nCH = 1 if gdename == "versailles"
replace nCH = 1 if gdename == "vienna"
replace nCH = 1 if gdename == "vincennes"
replace nCH = 1 if gdename == "virginia"
replace nCH = 1 if gdename == "wakefield"
replace nCH = 1 if gdename == "waltham"
replace nCH = 1 if gdename == "warren"
replace nCH = 1 if gdename == "warwick"
replace nCH = 1 if gdename == "washington"
replace nCH = 1 if gdename == "waterloo"
replace nCH = 1 if gdename == "waukegan"
replace nCH = 1 if gdename == "wayne"
replace nCH = 1 if gdename == "waynesboro"
replace nCH = 1 if gdename == "wembley"
replace nCH = 1 if gdename == "westbury"
replace nCH = 1 if gdename == "westfield"
replace nCH = 1 if gdename == "weston"
replace nCH = 1 if gdename == "westport"
replace nCH = 1 if gdename == "white plains"
replace nCH = 1 if gdename == "wichita"
replace nCH = 1 if gdename == "wiesbaden"
replace nCH = 1 if gdename == "willoughby"
replace nCH = 1 if gdename == "wilmette"
replace nCH = 1 if gdename == "wilmington"
replace nCH = 1 if gdename == "windsor"
replace nCH = 1 if gdename == "winnetka"
replace nCH = 1 if gdename == "winston salem"
replace nCH = 1 if gdename == "woodside"
replace nCH = 1 if gdename == "woodstock"
replace nCH = 1 if gdename == "wuppertal"
replace nCH = 1 if gdename == "youngstown"
replace nCH = 1 if gdename == "zagreb"
replace nCH = 1 if gdename == "zeeland"

replace nCH = 1 if gdename == "bruxelles"
replace nCH = 1 if gdename == "wien"
replace nCH = 1 if gdename == "muenchen"
replace nCH = 1 if gdename == "torino"
replace nCH = 1 if gdename == "laney"
replace nCH = 1 if gdename == "usa"
replace nCH = 1 if gdename == "duesseldorf"
replace nCH = 1 if gdename == "beyrouth"
replace nCH = 1 if gdename == "genova"
replace nCH = 1 if gdename == "kjoebenhavn"
replace nCH = 1 if gdename == "koeln"
replace nCH = 1 if gdename == "den haag"
replace nCH = 1 if gdename == "konstanz"
replace nCH = 1 if gdename == "sao paolo"
replace nCH = 1 if gdename == "goeteborg"
replace nCH = 1 if gdename == "varese"
replace nCH = 1 if gdename == "buenos ayres"
replace nCH = 1 if gdename == "tel aviv"
replace nCH = 1 if gdename == "nuernberg"
replace nCH = 1 if gdename == "mailand"
replace nCH = 1 if gdename == "neuilly/seine"
replace nCH = 1 if gdename == "greenwich"
replace nCH = 1 if gdename == "tokio"
replace nCH = 1 if gdename == "teheran"
replace nCH = 1 if gdename == "firenze"
replace nCH = 1 if gdename == "strassburg"
replace nCH = 1 if gdename == "bad homburg"
replace nCH = 1 if gdename == "bellegarde"
replace nCH = 1 if gdename == "monte carlo"
replace nCH = 1 if gdename == "annemasse"
replace nCH = 1 if gdename == "bad godesberg"
replace nCH = 1 if gdename == "bombay"
replace nCH = 1 if gdename == "dornbirn"
replace nCH = 1 if gdename == "malmoe"
replace nCH = 1 if gdename == "longeau"
replace nCH = 1 if gdename == "lahr"
replace nCH = 1 if gdename == "le vesinet"
replace nCH = 1 if gdename == "neuilly s/seine"
replace nCH = 1 if gdename == "schaetz"
replace nCH = 1 if gdename == "baden baden"
replace nCH = 1 if gdename == "frankfurt a. main"
replace nCH = 1 if gdename == "tanger"
replace nCH = 1 if gdename == "uccle"
replace nCH = 1 if gdename == "cannes"
replace nCH = 1 if gdename == "charbonnieres"
replace nCH = 1 if gdename == "riva s. vitale"
replace nCH = 1 if gdename == "bruessel"
replace nCH = 1 if gdename == "monte carlo"
replace nCH = 1 if gdename == "bloemendaal"
replace nCH = 1 if gdename == "napoli"
replace nCH = 1 if gdename == "neuilly sur seine"
replace nCH = 1 if gdename == "new canaan"
replace nCH = 1 if gdename == "heilbronn"
replace nCH = 1 if gdename == "hilversum"
replace nCH = 1 if gdename == "leopoldville"
replace nCH = 1 if gdename == "neuilly/paris"
replace nCH = 1 if gdename == "st. etienne"
replace nCH = 1 if gdename == "anvers"
replace nCH = 1 if gdename == "s gravenhage"
replace nCH = 1 if gdename == "short hills"
replace nCH = 1 if gdename == "wittlaer"
replace nCH = 1 if gdename == "alger"
replace nCH = 1 if gdename == "antwerp"
replace nCH = 1 if gdename == "enschede"
replace nCH = 1 if gdename == "londres"
replace nCH = 1 if gdename == "aerdenhout"
replace nCH = 1 if gdename == "cagialo"
replace nCH = 1 if gdename == "henan"
replace nCH = 1 if gdename == "island gb"
replace nCH = 1 if gdename == "kronberg"
replace nCH = 1 if gdename == "loerrach"
replace nCH = 1 if gdename == "neuss"
replace nCH = 1 if gdename == "oberhausen"
replace nCH = 1 if gdename == "rieben"
replace nCH = 1 if gdename == "rom"
replace nCH = 1 if gdename == "verviers"
replace nCH = 1 if gdename == "wassenaar"
replace nCH = 1 if gdename == "beaconsfield"
replace nCH = 1 if gdename == "beverley hills"
replace nCH = 1 if gdename == "carpi"
replace nCH = 1 if gdename == "celle"
replace nCH = 1 if gdename == "courtemaã®che"
replace nCH = 1 if gdename == "darmstadt"
replace nCH = 1 if gdename == "dietfurt"
replace nCH = 1 if gdename == "grand laney"
replace nCH = 1 if gdename == "luxemburg"
replace nCH = 1 if gdename == "rheinberg"
replace nCH = 1 if gdename == "sevenoaks"
replace nCH = 1 if gdename == "siegen"
replace nCH = 1 if gdename == "st cloud"
replace nCH = 1 if gdename == "vigevano"
replace nCH = 1 if gdename == "berlin dahlem"
replace nCH = 1 if gdename == "chateauneuf"
replace nCH = 1 if gdename == "cremona"
replace nCH = 1 if gdename == "croix/roubaix"
replace nCH = 1 if gdename == "domodossola"
replace nCH = 1 if gdename == "falkenstein"
replace nCH = 1 if gdename == "frankfurt a/main"
replace nCH = 1 if gdename == "furtwangen"
replace nCH = 1 if gdename == "haeusern"
replace nCH = 1 if gdename == "ivrea"
replace nCH = 1 if gdename == "kempten"
replace nCH = 1 if gdename == "pacific palisades"
replace nCH = 1 if gdename == "recklinghausen"
replace nCH = 1 if gdename == "roslyn"
replace nCH = 1 if gdename == "saarbruecken"
replace nCH = 1 if gdename == "springe"
replace nCH = 1 if gdename == "st etienne"
replace nCH = 1 if gdename == "waldshut"
replace nCH = 1 if gdename == "almelo"
replace nCH = 1 if gdename == "amersfoort"
replace nCH = 1 if gdename == "bellagio"
replace nCH = 1 if gdename == "bollate"
replace nCH = 1 if gdename == "boulogne billancourt"
replace nCH = 1 if gdename == "buederich"
replace nCH = 1 if gdename == "cobham"
replace nCH = 1 if gdename == "crainhem"
replace nCH = 1 if gdename == "detmold"
replace nCH = 1 if gdename == "dole"
replace nCH = 1 if gdename == "elisabethville"
replace nCH = 1 if gdename == "hallstahammar"
replace nCH = 1 if gdename == "hambourg"
replace nCH = 1 if gdename == "linz a.d. donau"
replace nCH = 1 if gdename == "neuffen"
replace nCH = 1 if gdename == "neuilly"
replace nCH = 1 if gdename == "princetown"
replace nCH = 1 if gdename == "rothenbach"
replace nCH = 1 if gdename == "schrofen"
replace nCH = 1 if gdename == "belgrad"
replace nCH = 1 if gdename == "breda"
replace nCH = 1 if gdename == "cernobbio"
replace nCH = 1 if gdename == "chamalieres"
replace nCH = 1 if gdename == "drize"
replace nCH = 1 if gdename == "frankfurt am main"
replace nCH = 1 if gdename == "isle"
replace nCH = 1 if gdename == "kjobenhavn"
replace nCH = 1 if gdename == "limpsfield"
replace nCH = 1 if gdename == "neuilly s/ seine"
replace nCH = 1 if gdename == "nijkerk"
replace nCH = 1 if gdename == "padova"
replace nCH = 1 if gdename == "ramat gan"
replace nCH = 1 if gdename == "reconvillier"
replace nCH = 1 if gdename == "sceaux"
replace nCH = 1 if gdename == "st. louis"
replace nCH = 1 if gdename == "tollen"
replace nCH = 1 if gdename == "wuppertal barmen"
replace nCH = 1 if gdename == "(russland)"
replace nCH = 1 if gdename == "amersfoot"
replace nCH = 1 if gdename == "amphion"
replace nCH = 1 if gdename == "beograd"
replace nCH = 1 if gdename == "bioux"
replace nCH = 1 if gdename == "bochum"
replace nCH = 1 if gdename == "boll vechigen"
replace nCH = 1 if gdename == "boulogne s/seine"
replace nCH = 1 if gdename == "britannique"
replace nCH = 1 if gdename == "dorridge"
replace nCH = 1 if gdename == "essex gb"
replace nCH = 1 if gdename == "francaises"
replace nCH = 1 if gdename == "gelsenkirchen"
replace nCH = 1 if gdename == "lewittown"
replace nCH = 1 if gdename == "luebeck"
replace nCH = 1 if gdename == "meudon"
replace nCH = 1 if gdename == "montmorency"
replace nCH = 1 if gdename == "pforzheim"
replace nCH = 1 if gdename == "remscheid"
replace nCH = 1 if gdename == "san sebastien"
replace nCH = 1 if gdename == "schuh"
replace nCH = 1 if gdename == "st mande"
replace nCH = 1 if gdename == "stocksund"
replace nCH = 1 if gdename == "st paul"
replace nCH = 1 if gdename == "struemp"
replace nCH = 1 if gdename == "ville d avray"
replace nCH = 1 if gdename == "villeurbanne"
replace nCH = 1 if gdename == "wilhelmshaven"
replace nCH = 1 if gdename == "(deutschland)"
replace nCH = 1 if gdename == "bad wiessee"
replace nCH = 1 if gdename == "cernusco"
replace nCH = 1 if gdename == "collonges sous saleve"
replace nCH = 1 if gdename == "dilbeek"
replace nCH = 1 if gdename == "esher"
replace nCH = 1 if gdename == "gerrards cross"
replace nCH = 1 if gdename == "goldern"
replace nCH = 1 if gdename == "guildford"
replace nCH = 1 if gdename == "heidenheim"
replace nCH = 1 if gdename == "helsingfors"
replace nCH = 1 if gdename == "high wycombe"
replace nCH = 1 if gdename == "hilden"
replace nCH = 1 if gdename == "krefeld"
replace nCH = 1 if gdename == "levittown"
replace nCH = 1 if gdename == "lisboa"
replace nCH = 1 if gdename == "lister"
replace nCH = 1 if gdename == "loiano"
replace nCH = 1 if gdename == "ludwigsburg"
replace nCH = 1 if gdename == "manilla"
replace nCH = 1 if gdename == "munchen"
replace nCH = 1 if gdename == "neuenstadt"
replace nCH = 1 if gdename == "new jersey"
replace nCH = 1 if gdename == "newyork"
replace nCH = 1 if gdename == "paterno"
replace nCH = 1 if gdename == "rehme"
replace nCH = 1 if gdename == "st paul"
replace nCH = 1 if gdename == "st germain en laye"
replace nCH = 1 if gdename == "thonon"
replace nCH = 1 if gdename == "thonon les bains"
replace nCH = 1 if gdename == "venezia"
replace nCH = 1 if gdename == "vesinet"
replace nCH = 1 if gdename == "waiblingen"
replace nCH = 1 if gdename == "wuppertal barmen"
replace nCH = 1 if gdename == "yorkschire gb"
replace nCH = 1 if gdename == "zoagli"
replace nCH = 1 if gdename == "(sri lanka)"
replace nCH = 1 if gdename == "addis abeba"
replace nCH = 1 if gdename == "ahlen"
replace nCH = 1 if gdename == "aix en provence"
replace nCH = 1 if gdename == "alkmaar"
replace nCH = 1 if gdename == "ammann"
replace nCH = 1 if gdename == "angleur les liege"
replace nCH = 1 if gdename == "ansbach"
replace nCH = 1 if gdename == "apeldoorn"
replace nCH = 1 if gdename == "ashtead"
replace nCH = 1 if gdename == "auenhofen"
replace nCH = 1 if gdename == "baveno"
replace nCH = 1 if gdename == "bierth"
replace nCH = 1 if gdename == "bottrop"
replace nCH = 1 if gdename == "bourgival"
replace nCH = 1 if gdename == "brescia"
replace nCH = 1 if gdename == "bueckenburg"
replace nCH = 1 if gdename == "buehlertal"
replace nCH = 1 if gdename == "calcutta"
replace nCH = 1 if gdename == "capetown"
replace nCH = 1 if gdename == "cascais"
replace nCH = 1 if gdename == "castelletto di brenzone"
replace nCH = 1 if gdename == "charlottenlund"
replace nCH = 1 if gdename == "chervilles"
replace nCH = 1 if gdename == "chiavenna"
replace nCH = 1 if gdename == "connecticut"
replace nCH = 1 if gdename == "courbevoie"
replace nCH = 1 if gdename == "couvin"
replace nCH = 1 if gdename == "dammarie les lys"
replace nCH = 1 if gdename == "deutschland"
replace nCH = 1 if gdename == "diursholm"
replace nCH = 1 if gdename == "djursholm"
replace nCH = 1 if gdename == "domblans"
replace nCH = 1 if gdename == "dover plains"
replace nCH = 1 if gdename == "dr., arlaching"
replace nCH = 1 if gdename == "dunfermline"
replace nCH = 1 if gdename == "eastbourne"
replace nCH = 1 if gdename == "ebersbach"
replace nCH = 1 if gdename == "edgbaston"
replace nCH = 1 if gdename == "einbeck"
replace nCH = 1 if gdename == "ekeren antwerpen"
replace nCH = 1 if gdename == "epping"
replace nCH = 1 if gdename == "erzenholz"
replace nCH = 1 if gdename == "essen ruhr"
replace nCH = 1 if gdename == "esslingen am neckar"
replace nCH = 1 if gdename == "fairlight"
replace nCH = 1 if gdename == "fallbrook"
replace nCH = 1 if gdename == "far hills"
replace nCH = 1 if gdename == "felstead"
replace nCH = 1 if gdename == "forest brussel"
replace nCH = 1 if gdename == "frankenstein"
replace nCH = 1 if gdename == "frankfurt/main"
replace nCH = 1 if gdename == "fuerth"
replace nCH = 1 if gdename == "gentbrugge"
replace nCH = 1 if gdename == "genua"
replace nCH = 1 if gdename == "gerardmer"
replace nCH = 1 if gdename == "graz kroisbach"
replace nCH = 1 if gdename == "grez doiceau"
replace nCH = 1 if gdename == "gruenwald"
replace nCH = 1 if gdename == "gruenwald/muenchen"
replace nCH = 1 if gdename == "gustavsberg"
replace nCH = 1 if gdename == "haddenham"
replace nCH = 1 if gdename == "hameln"
replace nCH = 1 if gdename == "haste"
replace nCH = 1 if gdename == "haven ct"
replace nCH = 1 if gdename == "heepen"
replace nCH = 1 if gdename == "heilbron boecklingen"
replace nCH = 1 if gdename == "heilbronn boeckingen"
replace nCH = 1 if gdename == "herne"
replace nCH = 1 if gdename == "herrlingen"
replace nCH = 1 if gdename == "hovas"
replace nCH = 1 if gdename == "huizen"
replace nCH = 1 if gdename == "hutington"
replace nCH = 1 if gdename == "kaessel"
replace nCH = 1 if gdename == "kapfenberg"
replace nCH = 1 if gdename == "kerteminde"
replace nCH = 1 if gdename == "kierspe"
replace nCH = 1 if gdename == "kopenhagen"
replace nCH = 1 if gdename == "kraainem"
replace nCH = 1 if gdename == "kumasi (west afrika)"
replace nCH = 1 if gdename == "la louviere"
replace nCH = 1 if gdename == "la mulatiere"
replace nCH = 1 if gdename == "landeck"
replace nCH = 1 if gdename == "langenbruecke"
replace nCH = 1 if gdename == "lavagna"
replace nCH = 1 if gdename == "loana"
replace nCH = 1 if gdename == "lobberich"
replace nCH = 1 if gdename == "lomazzo"
replace nCH = 1 if gdename == "long island"
replace nCH = 1 if gdename == "lund"
replace nCH = 1 if gdename == "luzenac"
replace nCH = 1 if gdename == "manz"
replace nCH = 1 if gdename == "mariastein metzerlen"
replace nCH = 1 if gdename == "mettlach"
replace nCH = 1 if gdename == "miex"
replace nCH = 1 if gdename == "missiones"
replace nCH = 1 if gdename == "montecarasso"
replace nCH = 1 if gdename == "montendrey"
replace nCH = 1 if gdename == "moulins"
replace nCH = 1 if gdename == "muelheim/ruhr"
replace nCH = 1 if gdename == "neuilly s. seine (france)"
replace nCH = 1 if gdename == "neuss am rhein"
replace nCH = 1 if gdename == "neviges"
replace nCH = 1 if gdename == "new york city"
replace nCH = 1 if gdename == "opladen"
replace nCH = 1 if gdename == "osteno"
replace nCH = 1 if gdename == "paris (france)"
replace nCH = 1 if gdename == "paris neuilly"
replace nCH = 1 if gdename == "perles"
replace nCH = 1 if gdename == "perstorp"
replace nCH = 1 if gdename == "prince, garches"
replace nCH = 1 if gdename == "purmarend"
replace nCH = 1 if gdename == "radolfzell"
replace nCH = 1 if gdename == "ratingen"
replace nCH = 1 if gdename == "reichenbach/fils"
replace nCH = 1 if gdename == "rezzato"
replace nCH = 1 if gdename == "rhode st genise"
replace nCH = 1 if gdename == "rodenkirchen"
replace nCH = 1 if gdename == "rolling hills"
replace nCH = 1 if gdename == "rotzloch"
replace nCH = 1 if gdename == "roubaix"
replace nCH = 1 if gdename == "rueil malmaison"
replace nCH = 1 if gdename == "salonica"
replace nCH = 1 if gdename == "san felin de guixols"
replace nCH = 1 if gdename == "san abbondio"
replace nCH = 1 if gdename == "sassuolo"
replace nCH = 1 if gdename == "schottland"
replace nCH = 1 if gdename == "selb ploessberg"
replace nCH = 1 if gdename == "septmoncel"
replace nCH = 1 if gdename == "sevres"
replace nCH = 1 if gdename == "sewickley heights"
replace nCH = 1 if gdename == "seyssel"
replace nCH = 1 if gdename == "simone/vacallo"
replace nCH = 1 if gdename == "sondrio"
replace nCH = 1 if gdename == "sopo"
replace nCH = 1 if gdename == "southbury"
replace nCH = 1 if gdename == "st antonio"
replace nCH = 1 if gdename == "st. georgen"
replace nCH = 1 if gdename == "st. ingbert"
replace nCH = 1 if gdename == "starnberg a, see"
replace nCH = 1 if gdename == "staufberg"
replace nCH = 1 if gdename == "steyr"
replace nCH = 1 if gdename == "stresa"
replace nCH = 1 if gdename == "st roch par sallanches"
replace nCH = 1 if gdename == "st sebastien s/loire"
replace nCH = 1 if gdename == "sturbridge"
replace nCH = 1 if gdename == "sulhamstead"
replace nCH = 1 if gdename == "surabaia"
replace nCH = 1 if gdename == "sutton coldfield"
replace nCH = 1 if gdename == "tangers"
replace nCH = 1 if gdename == "telfs"
replace nCH = 1 if gdename == "thiais"
replace nCH = 1 if gdename == "thionville"
replace nCH = 1 if gdename == "tientsin"
replace nCH = 1 if gdename == "tivon bei haifa"
replace nCH = 1 if gdename == "torno"
replace nCH = 1 if gdename == "tourcoing"
replace nCH = 1 if gdename == "treiton"
replace nCH = 1 if gdename == "tuebingen"
replace nCH = 1 if gdename == "tuscherz alfermee"
replace nCH = 1 if gdename == "uccles"
replace nCH = 1 if gdename == "varzi"
replace nCH = 1 if gdename == "villingen"
replace nCH = 1 if gdename == "villy"
replace nCH = 1 if gdename == "viotho hollweisen"
replace nCH = 1 if gdename == "vogelzang"
replace nCH = 1 if gdename == "waldhaus"
replace nCH = 1 if gdename == "wauthier braine"
replace nCH = 1 if gdename == "weil der stadt"
replace nCH = 1 if gdename == "welwyn garden city"
replace nCH = 1 if gdename == "west hartford"
replace nCH = 1 if gdename == "westcliff on sea"
replace nCH = 1 if gdename == "withefish bay"
replace nCH = 1 if gdename == "wyhlen"
replace nCH = 1 if gdename == "(suedafrika)"
replace nCH = 1 if gdename == "anet"
replace nCH = 1 if gdename == "ashridge"
replace nCH = 1 if gdename == "athenae"
replace nCH = 1 if gdename == "baconsfield"
replace nCH = 1 if gdename == "baden/wien"
replace nCH = 1 if gdename == "beech creek"
replace nCH = 1 if gdename == "den hag"
replace nCH = 1 if gdename == "denville"
replace nCH = 1 if gdename == "dr. ing., wuppertal"
replace nCH = 1 if gdename == "fairsfield"
replace nCH = 1 if gdename == "(buckinghamshire. gb)"
replace nCH = 1 if gdename == "(massachusetts usa)"
replace nCH = 1 if gdename == "(west sussex gb)"
replace nCH = 1 if gdename == ". new york"
replace nCH = 1 if gdename == "afeka"
replace nCH = 1 if gdename == "amstelveen"
replace nCH = 1 if gdename == "angmering on sea"
replace nCH = 1 if gdename == "ann harbor"
replace nCH = 1 if gdename == "atvidaberg"
replace nCH = 1 if gdename == "bad wildungen"
replace nCH = 1 if gdename == "berlin w"
replace nCH = 1 if gdename == "(luxembourg)"
replace nCH = 1 if gdename == "(michigan usa)"
replace nCH = 1 if gdename == "aix les bains"
replace nCH = 1 if gdename == "alep"
replace nCH = 1 if gdename == "amsterdam/usa"
replace nCH = 1 if gdename == "appelhueslen"
replace nCH = 1 if gdename == "asiago"
replace nCH = 1 if gdename == "beiruth"
replace nCH = 1 if gdename == "bergneustadt"
replace nCH = 1 if gdename == "berlin charlottenburg"
replace nCH = 1 if gdename == "berlin dahlem"
replace nCH = 1 if gdename == "bernei"
replace nCH = 1 if gdename == "(monaco)"
replace nCH = 1 if gdename == "(virginia usa)"
replace nCH = 1 if gdename == "aix les bains"
replace nCH = 1 if gdename == "ameres"
replace nCH = 1 if gdename == "amersham"
replace nCH = 1 if gdename == "ammenhausen"
replace nCH = 1 if gdename == "anglo normandes"
replace nCH = 1 if gdename == "arare"
replace nCH = 1 if gdename == "argenteau"
replace nCH = 1 if gdename == "argenteuil"
replace nCH = 1 if gdename == "arhedo"
replace nCH = 1 if gdename == "athenes"
replace nCH = 1 if gdename == "azcoita"
replace nCH = 1 if gdename == "bad aussee"
replace nCH = 1 if gdename == "bad kissingen"
replace nCH = 1 if gdename == "bad salzuften"
replace nCH = 1 if gdename == "bad kreuznach"
replace nCH = 1 if gdename == "bedford hills"
replace nCH = 1 if gdename == "beirut (libanon)"
replace nCH = 1 if gdename == "belfort"
replace nCH = 1 if gdename == "berlin grunewald"
replace nCH = 1 if gdename == "(mets)"
replace nCH = 1 if gdename == "(pennsylvania usa)"
replace nCH = 1 if gdename == "(portugal)"
replace nCH = 1 if gdename == "(sussex. gb)"
replace nCH = 1 if gdename == ". milano"
replace nCH = 1 if gdename == ".kjobenhavn"
replace nCH = 1 if gdename == ".toulouse"
replace nCH = 1 if gdename == "ik, bloemendaal"
replace nCH = 1 if gdename == ".oeteborg"
replace nCH = 1 if gdename == ", london"
replace nCH = 1 if gdename == "aevdenhout"
replace nCH = 1 if gdename == "ailmendingen"
replace nCH = 1 if gdename == "aix lest bains"
replace nCH = 1 if gdename == "alexandrie"
replace nCH = 1 if gdename == "alexandrien"
replace nCH = 1 if gdename == "alfortville"
replace nCH = 1 if gdename == "alimendingen"
replace nCH = 1 if gdename == "allison park"
replace nCH = 1 if gdename == "altusried"
replace nCH = 1 if gdename == "alwoodly"
replace nCH = 1 if gdename == "amhearst"
replace nCH = 1 if gdename == "amphion (france)"
replace nCH = 1 if gdename == "amphion (haute savoie france)"
replace nCH = 1 if gdename == "amphion (hte savoie)"
replace nCH = 1 if gdename == "amstetten"
replace nCH = 1 if gdename == "andorre"
replace nCH = 1 if gdename == "angenstein"
replace nCH = 1 if gdename == "annecy (france)"
replace nCH = 1 if gdename == "annemass e"
replace nCH = 1 if gdename == "anteerpen"
replace nCH = 1 if gdename == "antibes"
replace nCH = 1 if gdename == "antillen"
replace nCH = 1 if gdename == "antwerpen (belgien)"
replace nCH = 1 if gdename == "apenrade (daenemark)"
replace nCH = 1 if gdename == "arabie saoudite"
replace nCH = 1 if gdename == "aran s. villette"
replace nCH = 1 if gdename == "arcachon"
replace nCH = 1 if gdename == "arma"
replace nCH = 1 if gdename == "ascot"
replace nCH = 1 if gdename == "asnieres"
replace nCH = 1 if gdename == "athen"
replace nCH = 1 if gdename == "audenshaw (grossbritannien)"
replace nCH = 1 if gdename == "aureil (france)"
replace nCH = 1 if gdename == "b uxelles"
replace nCH = 1 if gdename == "baarn"
replace nCH = 1 if gdename == "babel"
replace nCH = 1 if gdename == "backnang"
replace nCH = 1 if gdename == "bad hombure"
replace nCH = 1 if gdename == "bad honnef"
replace nCH = 1 if gdename == "bad kreuznach"
replace nCH = 1 if gdename == "bad pyrntont"
replace nCH = 1 if gdename == "bad raeaz"
replace nCH = 1 if gdename == "bad raga."
replace nCH = 1 if gdename == "bad salzgitter"
replace nCH = 1 if gdename == "bad soden"
replace nCH = 1 if gdename == "badhoevedorp"
replace nCH = 1 if gdename == "bad homburg"
replace nCH = 1 if gdename == "badmatte lengnau (be)"
replace nCH = 1 if gdename == "bassecourt (j. b.)"
replace nCH = 1 if gdename == "bassecourt. (j. b.)"
replace nCH = 1 if gdename == "baveno (italien)"
replace nCH = 1 if gdename == "beaune (france)"
replace nCH = 1 if gdename == "bellgarde (france)"
replace nCH = 1 if gdename == "belmonte"
replace nCH = 1 if gdename == "bergisch gladbach"
replace nCH = 1 if gdename == "berkeley (californien)"
replace nCH = 1 if gdename == "berlin west"
replace nCH = 1 if gdename == "berlin gruenau"
replace nCH = 1 if gdename == "berlin gruenewald"
replace nCH = 1 if gdename == "berlin lichterfelde"
replace nCH = 1 if gdename == "berlin llchterfelde"
replace nCH = 1 if gdename == "berlinyen"
replace nCH = 1 if gdename == "berndorf (deutschland)"
replace nCH = 1 if gdename == "bernhausen"

// further round of manual coding (Geo_coding_gdename2 -- > last iteratioin of "stillmissing")
replace nCH = 1 if gdename == "111, grand rapids"
replace nCH = 1 if gdename == "alassio"
replace nCH = 1 if gdename == "almhult"
replace nCH = 1 if gdename == "ambilly"
replace nCH = 1 if gdename == "ancourt/dieppe"
replace nCH = 1 if gdename == "aprilia"
replace nCH = 1 if gdename == "arolsen"
replace nCH = 1 if gdename == "arona"
replace nCH = 1 if gdename == "arpajon"
replace nCH = 1 if gdename == "arveyres/ollon"
replace nCH = 1 if gdename == "aubemas"
replace nCH = 1 if gdename == "baerchen"
replace nCH = 1 if gdename == "bahia"
replace nCH = 1 if gdename == "bambach/furth"
replace nCH = 1 if gdename == "barberier"
replace nCH = 1 if gdename == "baron, forest les bruxelles (belg.)"
replace nCH = 1 if gdename == "bedburg"
replace nCH = 1 if gdename == "bensberg"
replace nCH = 1 if gdename == "bewdly"
replace nCH = 1 if gdename == "billund"
replace nCH = 1 if gdename == "blomfield"
replace nCH = 1 if gdename == "bloomfield hills"
replace nCH = 1 if gdename == "bofors"
replace nCH = 1 if gdename == "bondy"
replace nCH = 1 if gdename == "bonneval"
replace nCH = 1 if gdename == "bourgne"
replace nCH = 1 if gdename == "bretigny s. orge"
replace nCH = 1 if gdename == "bridgenorth"
replace nCH = 1 if gdename == "brignon"
replace nCH = 1 if gdename == "broadview"
replace nCH = 1 if gdename == "brookline"
replace nCH = 1 if gdename == "bryon"
replace nCH = 1 if gdename == "caderhurst"
replace nCH = 1 if gdename == "canelli"
replace nCH = 1 if gdename == "caracas como"
replace nCH = 1 if gdename == "casate nuovo"
replace nCH = 1 if gdename == "chambery"
replace nCH = 1 if gdename == "charlton kings"
replace nCH = 1 if gdename == "choisy le roi"
replace nCH = 1 if gdename == "codogno"
replace nCH = 1 if gdename == "croissy/seine"
replace nCH = 1 if gdename == "danderyd"
replace nCH = 1 if gdename == "dinslaken"
replace nCH = 1 if gdename == "driebergen"
replace nCH = 1 if gdename == "dundrum"
replace nCH = 1 if gdename == "ebingen"
replace nCH = 1 if gdename == "ecully"
replace nCH = 1 if gdename == "emirate"
replace nCH = 1 if gdename == "emmenau"
replace nCH = 1 if gdename == "emsdetten"
replace nCH = 1 if gdename == "england"
replace nCH = 1 if gdename == "erkrath/duesseldorf"
replace nCH = 1 if gdename == "essen werden"
replace nCH = 1 if gdename == "etang la ville"
replace nCH = 1 if gdename == "fetcham"
replace nCH = 1 if gdename == "france"
replace nCH = 1 if gdename == "frankfurt a m"
replace nCH = 1 if gdename == "frankfurt s/main"
replace nCH = 1 if gdename == "frechen bei koeln"
replace nCH = 1 if gdename == "friedrichshafen"
replace nCH = 1 if gdename == "gaillard"
replace nCH = 1 if gdename == "gand"
replace nCH = 1 if gdename == "gandia"
replace nCH = 1 if gdename == "garmisch partenkirchen"
replace nCH = 1 if gdename == "gastonbury"
replace nCH = 1 if gdename == "gates mill"
replace nCH = 1 if gdename == "gehringen"
replace nCH = 1 if gdename == "geislingen"
replace nCH = 1 if gdename == "genese"
replace nCH = 1 if gdename == "gerrard cross"
replace nCH = 1 if gdename == "gevelsberg"
replace nCH = 1 if gdename == "ghent"
replace nCH = 1 if gdename == "giuncarico"
replace nCH = 1 if gdename == "gottmadingen"
replace nCH = 1 if gdename == "graf, muenchen"
replace nCH = 1 if gdename == "great witley"
replace nCH = 1 if gdename == "grindsted"
replace nCH = 1 if gdename == "grosse point"
replace nCH = 1 if gdename == "gutenberg"
replace nCH = 1 if gdename == "hamburg rissen"
replace nCH = 1 if gdename == "hannover langenhagen"
replace nCH = 1 if gdename == "hassloc pfalz"
replace nCH = 1 if gdename == "headcorn"
replace nCH = 1 if gdename == "heiligenhaus"
replace nCH = 1 if gdename == "hellenthal eifel"
replace nCH = 1 if gdename == "hem"
replace nCH = 1 if gdename == "herdecke"
replace nCH = 1 if gdename == "hoesel/duesseldorf"
replace nCH = 1 if gdename == "hoexter"
replace nCH = 1 if gdename == "hoogboom kapellen"
replace nCH = 1 if gdename == "hueckelhoven"
replace nCH = 1 if gdename == "hyogo ken"
replace nCH = 1 if gdename == "issy les moulinaux"
replace nCH = 1 if gdename == "ixelles"
replace nCH = 1 if gdename == "kautzen"
replace nCH = 1 if gdename == "knivsta"
replace nCH = 1 if gdename == "koeln marienburg"
replace nCH = 1 if gdename == "kopparberg"
replace nCH = 1 if gdename == "koweit"
replace nCH = 1 if gdename == "krefeld traar ascona"
replace nCH = 1 if gdename == "la clusette et mijoux"
replace nCH = 1 if gdename == "la havana"
replace nCH = 1 if gdename == "laag keppel"
replace nCH = 1 if gdename == "lamberhurst"
replace nCH = 1 if gdename == "laren"
replace nCH = 1 if gdename == "lasaraaz"
replace nCH = 1 if gdename == "lerthe"
replace nCH = 1 if gdename == "lindenfels"
replace nCH = 1 if gdename == "loncin"
replace nCH = 1 if gdename == "ludvika"
replace nCH = 1 if gdename == "luino"
replace nCH = 1 if gdename == "lustanau"
replace nCH = 1 if gdename == "maasland"
replace nCH = 1 if gdename == "macon"
replace nCH = 1 if gdename == "maresfield"
replace nCH = 1 if gdename == "mari"
replace nCH = 1 if gdename == "masevaux"
replace nCH = 1 if gdename == "meinherzhagen"
replace nCH = 1 if gdename == "meise"
replace nCH = 1 if gdename == "mont s/marchienne"
replace nCH = 1 if gdename == "monte video"
replace nCH = 1 if gdename == "mount royal"
replace nCH = 1 if gdename == "muenster westfalen"
replace nCH = 1 if gdename == "murcie"
replace nCH = 1 if gdename == "nancy et chelles"
replace nCH = 1 if gdename == "neu ulm"
replace nCH = 1 if gdename == "new malden"
replace nCH = 1 if gdename == "newdigate"
replace nCH = 1 if gdename == "ninove"
replace nCH = 1 if gdename == "nordhorn"
replace nCH = 1 if gdename == "novedrate"
replace nCH = 1 if gdename == "novozzano"
replace nCH = 1 if gdename == "oakington"
replace nCH = 1 if gdename == "oberstdorf"
replace nCH = 1 if gdename == "old greenwich"
replace nCH = 1 if gdename == "olympia fields"
replace nCH = 1 if gdename == "ordrup"
replace nCH = 1 if gdename == "osnabrueck"
replace nCH = 1 if gdename == "partinico"
replace nCH = 1 if gdename == "pau"
replace nCH = 1 if gdename == "penhurst"
replace nCH = 1 if gdename == "pieter willem, amsterdam"
replace nCH = 1 if gdename == "porlezza"
replace nCH = 1 if gdename == "port elisabeth"
replace nCH = 1 if gdename == "postal"
replace nCH = 1 if gdename == "prestbury"
replace nCH = 1 if gdename == "radenthein"
replace nCH = 1 if gdename == "radnor"
replace nCH = 1 if gdename == "ratingen duesseldorf"
replace nCH = 1 if gdename == "redwood city"
replace nCH = 1 if gdename == "reggio emilia italia"
replace nCH = 1 if gdename == "reutlingen"
replace nCH = 1 if gdename == "rhinebeck"
replace nCH = 1 if gdename == "rochford"
replace nCH = 1 if gdename == "saronno"
replace nCH = 1 if gdename == "scants"
replace nCH = 1 if gdename == "schoepfheim"
replace nCH = 1 if gdename == "seeburg"
replace nCH = 1 if gdename == "selva"
replace nCH = 1 if gdename == "sevres (france)"
replace nCH = 1 if gdename == "shalford"
replace nCH = 1 if gdename == "siantar"
replace nCH = 1 if gdename == "sierra"
replace nCH = 1 if gdename == "sinn"
replace nCH = 1 if gdename == "skytop"
replace nCH = 1 if gdename == "sogogn"
replace nCH = 1 if gdename == "solligen"
replace nCH = 1 if gdename == "southborough"
replace nCH = 1 if gdename == "st john s"
replace nCH = 1 if gdename == "st maur des fosses"
replace nCH = 1 if gdename == "st petersburg"
replace nCH = 1 if gdename == "st. andreasberg"
replace nCH = 1 if gdename == "st. thomas"
replace nCH = 1 if gdename == "stellenbosch"
replace nCH = 1 if gdename == "suor, parma"
replace nCH = 1 if gdename == "surbiton"
replace nCH = 1 if gdename == "torno/como"
replace nCH = 1 if gdename == "torremolinos"
replace nCH = 1 if gdename == "trelleborg"
replace nCH = 1 if gdename == "triest"
replace nCH = 1 if gdename == "trivero"
replace nCH = 1 if gdename == "tullard drumbo"
replace nCH = 1 if gdename == "valbrona"
replace nCH = 1 if gdename == "valence"
replace nCH = 1 if gdename == "velp"
replace nCH = 1 if gdename == "vexin"
replace nCH = 1 if gdename == "vierlingsbeek"
replace nCH = 1 if gdename == "vossem"
replace nCH = 1 if gdename == "wasquehal"
replace nCH = 1 if gdename == "wavre"
replace nCH = 1 if gdename == "wels"
replace nCH = 1 if gdename == "wentworth"
replace nCH = 1 if gdename == "wessling oberbayern"
replace nCH = 1 if gdename == "weybridge"
replace nCH = 1 if gdename == "wiedenbrueck"
replace nCH = 1 if gdename == "wiel"
replace nCH = 1 if gdename == "wilmongton"
replace nCH = 1 if gdename == "wolowe saint pierre"
replace nCH = 1 if gdename == "wuppertal/barmen"
replace nCH = 1 if gdename == "wupperthal"
replace nCH = 1 if gdename == "yeovil"
replace nCH = 1 if gdename == "zandvoort"

replace nCH = 1 if gdename == "paris f"
replace nCH = 1 if gdename == "london gb"
replace nCH = 1 if gdename == "milano i"
replace nCH = 1 if gdename == "new york usa"
replace nCH = 1 if gdename == "paris (f)"
replace nCH = 1 if gdename == "wien a"
replace nCH = 1 if gdename == "hamburg d"
replace nCH = 1 if gdename == "london (gb)"
replace nCH = 1 if gdename == "stuttgart d"
replace nCH = 1 if gdename == "duesseldorf d"
replace nCH = 1 if gdename == "bruxelles b"
replace nCH = 1 if gdename == "madrid e"
replace nCH = 1 if gdename == "koeln d"
replace nCH = 1 if gdename == "neuilly sur seine f"
replace nCH = 1 if gdename == "monaco mc"
replace nCH = 1 if gdename == "berlin d"
replace nCH = 1 if gdename == "monte carlo mc"
replace nCH = 1 if gdename == "rom i"
replace nCH = 1 if gdename == "londres gb"
replace nCH = 1 if gdename == "mailand i"
replace nCH = 1 if gdename == "amsterdam nl"
replace nCH = 1 if gdename == "muenchen (d)"
replace nCH = 1 if gdename == "mailand (i)"
replace nCH = 1 if gdename == "luxembourg l"
replace nCH = 1 if gdename == "frankfurt am main d"
replace nCH = 1 if gdename == "new york (usa)"
replace nCH = 1 if gdename == "roma i"
replace nCH = 1 if gdename == "istanbul tr"
replace nCH = 1 if gdename == "tokyo j"
replace nCH = 1 if gdename == "surrey gb"
replace nCH = 1 if gdename == "stockholm s"
replace nCH = 1 if gdename == "new jersey usa"
replace nCH = 1 if gdename == "wiesbaden d"
replace nCH = 1 if gdename == "konstanz d"
replace nCH = 1 if gdename == "hamburg (d)"
replace nCH = 1 if gdename == "tokio j"
replace nCH = 1 if gdename == "essen d"
replace nCH = 1 if gdename == "hong kong hk"
replace nCH = 1 if gdename == "wien (a)"
replace nCH = 1 if gdename == "bruessel b"
replace nCH = 1 if gdename == "barcelona e"
replace nCH = 1 if gdename == "torino i"
replace nCH = 1 if gdename == "oslo n"
replace nCH = 1 if gdename == "lyon f"
replace nCH = 1 if gdename == "loerrach d"
replace nCH = 1 if gdename == "turin i"
replace nCH = 1 if gdename == "como i"
replace nCH = 1 if gdename == "athen gr"
replace nCH = 1 if gdename == "helsinki sf"
replace nCH = 1 if gdename == "nuernberg d"
replace nCH = 1 if gdename == "texas usa"
replace nCH = 1 if gdename == "moskau rus"
replace nCH = 1 if gdename == "rom (i)"
replace nCH = 1 if gdename == "luxemburg l"
replace nCH = 1 if gdename == "hannover d"
replace nCH = 1 if gdename == "wassenaar nl"
replace nCH = 1 if gdename == "toronto cdn"
replace nCH = 1 if gdename == "los angeles usa"
replace nCH = 1 if gdename == "stuttgart (d)"
replace nCH = 1 if gdename == "rheinfelden d"
replace nCH = 1 if gdename == "karlsruhe d"
replace nCH = 1 if gdename == "connecticut usa"
replace nCH = 1 if gdename == "rotterdam nl"
replace nCH = 1 if gdename == "dornbirn a"
replace nCH = 1 if gdename == "bremen d"
replace nCH = 1 if gdename == "illinois usa"
replace nCH = 1 if gdename == "bruessel (b)"
replace nCH = 1 if gdename == "milano (i)"
replace nCH = 1 if gdename == "california usa"
replace nCH = 1 if gdename == "duesseldorf (d)"
replace nCH = 1 if gdename == "frankfurt d"
replace nCH = 1 if gdename == "tel aviv il"
replace nCH = 1 if gdename == "wuppertal d"
replace nCH = 1 if gdename == "panama pa"
replace nCH = 1 if gdename == "versailles f"
replace nCH = 1 if gdename == "annemasse f"
replace nCH = 1 if gdename == "dortmund d"
replace nCH = 1 if gdename == "montreal cdn"
replace nCH = 1 if gdename == "bad homburg d"
replace nCH = 1 if gdename == "panama"
replace nCH = 1 if gdename == "koeln (d)"
replace nCH = 1 if gdename == "pforzheim d"
replace nCH = 1 if gdename == "varese i"
replace nCH = 1 if gdename == "strasbourg f"
replace nCH = 1 if gdename == "madrid (e)"
replace nCH = 1 if gdename == "heidelberg d"
replace nCH = 1 if gdename == "bonn d"
replace nCH = 1 if gdename == "kalifornien usa"
replace nCH = 1 if gdename == "bergamo i"
replace nCH = 1 if gdename == "linz a"
replace nCH = 1 if gdename == "goeteborg s"
replace nCH = 1 if gdename == "gaillard f"
replace nCH = 1 if gdename == "mannheim d"
replace nCH = 1 if gdename == "chicago usa"
replace nCH = 1 if gdename == "darmstadt d"
replace nCH = 1 if gdename == "divonne les bains f"
replace nCH = 1 if gdename == "paris/f"
replace nCH = 1 if gdename == "stockholm (s)"
replace nCH = 1 if gdename == "neuilly f"
replace nCH = 1 if gdename == "kopenhagen dk"
replace nCH = 1 if gdename == "rome i"
replace nCH = 1 if gdename == "hong kong cn"
replace nCH = 1 if gdename == "bielefeld d"
replace nCH = 1 if gdename == "kent gb"
replace nCH = 1 if gdename == "bombay ind"
replace nCH = 1 if gdename == "baden baden d"
replace nCH = 1 if gdename == "lissabon p"
replace nCH = 1 if gdename == "graz a"
replace nCH = 1 if gdename == "genova i"
replace nCH = 1 if gdename == "monte carlo (mc)"
replace nCH = 1 if gdename == "meerbusch d"
replace nCH = 1 if gdename == "malmoe s"
replace nCH = 1 if gdename == "beyrouth rl"
replace nCH = 1 if gdename == "le vesinet f"
replace nCH = 1 if gdename == "thonon les bains f"
replace nCH = 1 if gdename == "tokio (j)"
replace nCH = 1 if gdename == "budapest h"
replace nCH = 1 if gdename == "augsburg d"
replace nCH = 1 if gdename == "reutlingen d"
replace nCH = 1 if gdename == "weil am rhein d"
replace nCH = 1 if gdename == "barcelone e"
replace nCH = 1 if gdename == "boulogne billancourt f"
replace nCH = 1 if gdename == "athenes gr"
replace nCH = 1 if gdename == "surrey (gb)"
replace nCH = 1 if gdename == "mainz d"
replace nCH = 1 if gdename == "milan i"
replace nCH = 1 if gdename == "kronberg d"
replace nCH = 1 if gdename == "bregenz a"
replace nCH = 1 if gdename == "mulhouse f"
replace nCH = 1 if gdename == "florida usa"
replace nCH = 1 if gdename == "krefeld d"
replace nCH = 1 if gdename == "istanbul (tr)"
replace nCH = 1 if gdename == "bussum nl"
replace nCH = 1 if gdename == "den haag nl"
replace nCH = 1 if gdename == "hilversum nl"
replace nCH = 1 if gdename == "hongkong hk"
replace nCH = 1 if gdename == "antwerpen b"
replace nCH = 1 if gdename == "pennsylvania usa"
replace nCH = 1 if gdename == "lustenau a"
replace nCH = 1 if gdename == "brasschaat b"
replace nCH = 1 if gdename == "ontario cdn"
replace nCH = 1 if gdename == "singapur sgp"
replace nCH = 1 if gdename == "freiburg im breisgau d"
replace nCH = 1 if gdename == "monza i"
replace nCH = 1 if gdename == "waterloo b"
replace nCH = 1 if gdename == "moscou rus"
replace nCH = 1 if gdename == "marseille f"
replace nCH = 1 if gdename == "oslo (n)"
replace nCH = 1 if gdename == "koenigstein d"
replace nCH = 1 if gdename == "collonges sous saleve f"
replace nCH = 1 if gdename == "ferney voltaire f"
replace nCH = 1 if gdename == "bangkok t"
replace nCH = 1 if gdename == "besancon f"
replace nCH = 1 if gdename == "ratingen d"
replace nCH = 1 if gdename == "bologna i"
replace nCH = 1 if gdename == "bochum d"
replace nCH = 1 if gdename == "rueil malmaison f"
replace nCH = 1 if gdename == "moskau udssr"
replace nCH = 1 if gdename == "muenster d"
replace nCH = 1 if gdename == "kairo et"
replace nCH = 1 if gdename == "holte dk"
replace nCH = 1 if gdename == "athen (gr)"
replace nCH = 1 if gdename == "salzburg a"
replace nCH = 1 if gdename == "riyadh sa"
replace nCH = 1 if gdename == "aerdenhout nl"
replace nCH = 1 if gdename == "taeby s"
replace nCH = 1 if gdename == "sao paulo br"
replace nCH = 1 if gdename == "padova i"
replace nCH = 1 if gdename == "weinheim d"
replace nCH = 1 if gdename == "jeddah sa"
replace nCH = 1 if gdename == "jersey gb"
replace nCH = 1 if gdename == "dublin irl"
replace nCH = 1 if gdename == "duisburg d"
replace nCH = 1 if gdename == "st julien en genevois f"
replace nCH = 1 if gdename == "dallas usa"
replace nCH = 1 if gdename == "blaricum nl"
replace nCH = 1 if gdename == "laren nl"
replace nCH = 1 if gdename == "michigan usa"
replace nCH = 1 if gdename == "napoli i"
replace nCH = 1 if gdename == "washington usa"
replace nCH = 1 if gdename == "london/gb"
replace nCH = 1 if gdename == "heilbronn d"
replace nCH = 1 if gdename == "ambilly f"
replace nCH = 1 if gdename == "maryland usa"
replace nCH = 1 if gdename == "berlin (d)"
replace nCH = 1 if gdename == "espoo sf"
replace nCH = 1 if gdename == "grenoble f"
replace nCH = 1 if gdename == "essen (d)"
replace nCH = 1 if gdename == "saarbruecken d"
replace nCH = 1 if gdename == "conn. usa"
replace nCH = 1 if gdename == "oberursel d"
replace nCH = 1 if gdename == "buenos aires ra"
replace nCH = 1 if gdename == "turin (i)"
replace nCH = 1 if gdename == "bergisch gladbach d"
replace nCH = 1 if gdename == "aachen d"
replace nCH = 1 if gdename == "muri be"
replace nCH = 1 if gdename == "amsterdam (nl)"
replace nCH = 1 if gdename == "beirut (rl)"
replace nCH = 1 if gdename == "kassel d"
replace nCH = 1 if gdename == "monaco (mc)"
replace nCH = 1 if gdename == "gruenwald d"
replace nCH = 1 if gdename == "firenze i"
replace nCH = 1 if gdename == "frankfurt (d)"
replace nCH = 1 if gdename == "neuilly sur seine (f)"
replace nCH = 1 if gdename == "sofia bg"
replace nCH = 1 if gdename == "ulm d"
replace nCH = 1 if gdename == "brescia i"
replace nCH = 1 if gdename == "ohio usa"
replace nCH = 1 if gdename == "jerusalem il"
replace nCH = 1 if gdename == "bruxelles (b)"
replace nCH = 1 if gdename == "berkshire gb"
replace nCH = 1 if gdename == "hongkong"
replace nCH = 1 if gdename == "overijse b"
replace nCH = 1 if gdename == "courbevoie f"
replace nCH = 1 if gdename == "genua (i)"
replace nCH = 1 if gdename == "djursholm s"
replace nCH = 1 if gdename == "verona i"
replace nCH = 1 if gdename == "middlesex gb"
replace nCH = 1 if gdename == "anvers b"
replace nCH = 1 if gdename == "bremen (d)"
replace nCH = 1 if gdename == "ludwigsburg d"
replace nCH = 1 if gdename == "rio de janeiro br"
replace nCH = 1 if gdename == "bromma s"
replace nCH = 1 if gdename == "sindelfingen d"
replace nCH = 1 if gdename == "boulogne f"
replace nCH = 1 if gdename == "treviso i"
replace nCH = 1 if gdename == "wiesbaden (d)"
replace nCH = 1 if gdename == "innsbruck a"
replace nCH = 1 if gdename == "ausland"
replace nCH = 1 if gdename == "johannesburg za"
replace nCH = 1 if gdename == "weybridge gb"
replace nCH = 1 if gdename == "greenwich usa"
replace nCH = 1 if gdename == "archamps f"
replace nCH = 1 if gdename == "tuebingen d"
replace nCH = 1 if gdename == "genua i"
replace nCH = 1 if gdename == "garches f"
replace nCH = 1 if gdename == "annecy f"
replace nCH = 1 if gdename == "cheshire gb"
replace nCH = 1 if gdename == "beirut rl"
replace nCH = 1 if gdename == "neuilly (f)"
replace nCH = 1 if gdename == "meudon f"
replace nCH = 1 if gdename == "huizen nl"
replace nCH = 1 if gdename == "tokyo (j)"
replace nCH = 1 if gdename == "lasne b"
replace nCH = 1 if gdename == "frankfurt a.m. d"
replace nCH = 1 if gdename == "lidingoe s"
replace nCH = 1 if gdename == "fuerth d"
replace nCH = 1 if gdename == "cernobbio i"
replace nCH = 1 if gdename == "teheran ir"
replace nCH = 1 if gdename == "freiburg d"
replace nCH = 1 if gdename == "st cloud f"
replace nCH = 1 if gdename == "mallorca e"
replace nCH = 1 if gdename == "nice f"
replace nCH = 1 if gdename == "atlanta usa"
replace nCH = 1 if gdename == "leverkusen d"
replace nCH = 1 if gdename == "wassenaar (nl)"
replace nCH = 1 if gdename == "vicenza i"
replace nCH = 1 if gdename == "kiel d"
replace nCH = 1 if gdename == "grenzach wyhlen d"
replace nCH = 1 if gdename == "konstanz (d)"
replace nCH = 1 if gdename == "frankfurt a. m. d"
replace nCH = 1 if gdename == "sydney aus"
replace nCH = 1 if gdename == "miami usa"
replace nCH = 1 if gdename == "los angeles (usa)"
replace nCH = 1 if gdename == "massachusetts usa"
replace nCH = 1 if gdename == "vaucresson f"
replace nCH = 1 if gdename == "naarden nl"
replace nCH = 1 if gdename == "moskau su"
replace nCH = 1 if gdename == "vetraz monthoux f"
replace nCH = 1 if gdename == "goeppingen d"
replace nCH = 1 if gdename == "alicante e"
replace nCH = 1 if gdename == "mannheim (d)"
replace nCH = 1 if gdename == "chatou f"
replace nCH = 1 if gdename == "guernsey gb"
replace nCH = 1 if gdename == "saint cloud f"
replace nCH = 1 if gdename == "frankfurt a.m"
replace nCH = 1 if gdename == "neuss d"
replace nCH = 1 if gdename == "klagenfurt a"
replace nCH = 1 if gdename == "livorno i"
replace nCH = 1 if gdename == "solingen d"
replace nCH = 1 if gdename == "frankfurt am main (d)"
replace nCH = 1 if gdename == "sark (channel islands)"
replace nCH = 1 if gdename == "casablanca ma"
replace nCH = 1 if gdename == "san francisco usa"
replace nCH = 1 if gdename == "bordeaux f"
replace nCH = 1 if gdename == "rotterdam (nl)"
replace nCH = 1 if gdename == "feldkirch a"
replace nCH = 1 if gdename == "remscheid d"
replace nCH = 1 if gdename == "friedrichshafen d"
replace nCH = 1 if gdename == "campione i"
replace nCH = 1 if gdename == "cal. usa"
replace nCH = 1 if gdename == "yokohama j"
replace nCH = 1 if gdename == "philadelphia usa"
replace nCH = 1 if gdename == "belgrad yu"
replace nCH = 1 if gdename == "new jersey (usa)"
replace nCH = 1 if gdename == "charlottenlund dk"
replace nCH = 1 if gdename == "helsinki (sf)"
replace nCH = 1 if gdename == "singen d"
replace nCH = 1 if gdename == "wuerzburg d"
replace nCH = 1 if gdename == "bilthoven nl"
replace nCH = 1 if gdename == "novara i"
replace nCH = 1 if gdename == "hegenheim f"
replace nCH = 1 if gdename == "hannover (d)"
replace nCH = 1 if gdename == "sciez f"
replace nCH = 1 if gdename == "luino i"
replace nCH = 1 if gdename == "amstelveen nl"
replace nCH = 1 if gdename == "waldshut tiengen d"
replace nCH = 1 if gdename == "campione d italia i"
replace nCH = 1 if gdename == "gex f"
replace nCH = 1 if gdename == "virginia usa"
replace nCH = 1 if gdename == "bucks gb"
replace nCH = 1 if gdename == "kraainem b"
replace nCH = 1 if gdename == "heidelberg (d)"
replace nCH = 1 if gdename == "luebeck d"
replace nCH = 1 if gdename == "freiburg i. breisgau d"
replace nCH = 1 if gdename == "nuernberg (d)"
replace nCH = 1 if gdename == "vught nl"
replace nCH = 1 if gdename == "ratingen (d)"
replace nCH = 1 if gdename == "ettlingen d"
replace nCH = 1 if gdename == "ville d avray f"
replace nCH = 1 if gdename == "douvaine f"
replace nCH = 1 if gdename == "mass. usa"
replace nCH = 1 if gdename == "marbella e"
replace nCH = 1 if gdename == "recklinghausen d"
replace nCH = 1 if gdename == "hoechst a"
replace nCH = 1 if gdename == "ischgl a"
replace nCH = 1 if gdename == "hagen d"
replace nCH = 1 if gdename == "barcelona (e)"
replace nCH = 1 if gdename == "york gb"
replace nCH = 1 if gdename == "prevessin moens f"
replace nCH = 1 if gdename == "wermelskirchen d"
replace nCH = 1 if gdename == "starnberg d"
replace nCH = 1 if gdename == "trieste i"
replace nCH = 1 if gdename == "goeteborg (s)"
replace nCH = 1 if gdename == "ravensburg d"
replace nCH = 1 if gdename == "leonberg d"
replace nCH = 1 if gdename == "curacao na"
replace nCH = 1 if gdename == "zagreb yu"
replace nCH = 1 if gdename == "bad saeckingen d"
replace nCH = 1 if gdename == "veigy foncenex f"
replace nCH = 1 if gdename == "cambridge gb"
replace nCH = 1 if gdename == "hampshire gb"
replace nCH = 1 if gdename == "italien"
replace nCH = 1 if gdename == "zeist nl"
replace nCH = 1 if gdename == "uccle b"
replace nCH = 1 if gdename == "erlangen d"
replace nCH = 1 if gdename == "sark gb"
replace nCH = 1 if gdename == "lisbonne p"
replace nCH = 1 if gdename == "vancouver cdn"
replace nCH = 1 if gdename == "moskau gus"
replace nCH = 1 if gdename == "new canaan usa"
replace nCH = 1 if gdename == "buenos aires (ra)"
replace nCH = 1 if gdename == "porto p"
replace nCH = 1 if gdename == "new york/usa"
replace nCH = 1 if gdename == "gauting d"
replace nCH = 1 if gdename == "kairo (et)"
replace nCH = 1 if gdename == "parma i"
replace nCH = 1 if gdename == "lund s"
replace nCH = 1 if gdename == "lindau d"
replace nCH = 1 if gdename == "houston usa"
replace nCH = 1 if gdename == "hard a"
replace nCH = 1 if gdename == "lyon (f)"
replace nCH = 1 if gdename == "isle of man"
replace nCH = 1 if gdename == "illinois (usa)"
replace nCH = 1 if gdename == "thoiry f"
replace nCH = 1 if gdename == "luxemburg (l)"
replace nCH = 1 if gdename == "kiew uk"
replace nCH = 1 if gdename == "westport usa"
replace nCH = 1 if gdename == "milano/i"
replace nCH = 1 if gdename == "hellerup dk"
replace nCH = 1 if gdename == "evian les bains f"
replace nCH = 1 if gdename == "soest nl"
replace nCH = 1 if gdename == "bad soden d"
replace nCH = 1 if gdename == "tel aviv (il)"
replace nCH = 1 if gdename == "louveciennes f"
replace nCH = 1 if gdename == "st germain en laye f"
replace nCH = 1 if gdename == "florenz i"
replace nCH = 1 if gdename == "dubai"
replace nCH = 1 if gdename == "hong kong (hk)"
replace nCH = 1 if gdename == "efringen kirchen d"
replace nCH = 1 if gdename == "como (i)"
replace nCH = 1 if gdename == "warschau pl"
replace nCH = 1 if gdename == "erkrath d"
replace nCH = 1 if gdename == "st louis f"
replace nCH = 1 if gdename == "teheran (ir)"
replace nCH = 1 if gdename == "san donato milanese i"
replace nCH = 1 if gdename == "milano 1"
replace nCH = 1 if gdename == "strassen l"
replace nCH = 1 if gdename == "schopfheim d"
replace nCH = 1 if gdename == "ville la grand f"
replace nCH = 1 if gdename == "detmold d"
replace nCH = 1 if gdename == "egg zh"
replace nCH = 1 if gdename == "oxford gb"
replace nCH = 1 if gdename == "sollentuna s"
replace nCH = 1 if gdename == "lauterach a"
replace nCH = 1 if gdename == "eindhoven nl"
replace nCH = 1 if gdename == "paderborn d"
replace nCH = 1 if gdename == "edinburgh gb"
replace nCH = 1 if gdename == "sevenoaks gb"
replace nCH = 1 if gdename == "alger dz"
replace nCH = 1 if gdename == "caracas yv"
replace nCH = 1 if gdename == "douglas"
replace nCH = 1 if gdename == "liege b"
replace nCH = 1 if gdename == "dakar sn"
replace nCH = 1 if gdename == "helsingborg s"
replace nCH = 1 if gdename == "leiden nl"
replace nCH = 1 if gdename == "dueren d"
replace nCH = 1 if gdename == "strassburg f"
replace nCH = 1 if gdename == "birmingham gb"
replace nCH = 1 if gdename == "troisdorf d"
replace nCH = 1 if gdename == "gb"
replace nCH = 1 if gdename == "osaka j"
replace nCH = 1 if gdename == "limassol cy"
replace nCH = 1 if gdename == "ecully f"
replace nCH = 1 if gdename == "quebec cdn"
replace nCH = 1 if gdename == "georgia usa"
replace nCH = 1 if gdename == "marseille (f)"
replace nCH = 1 if gdename == "bondues f"
replace nCH = 1 if gdename == "venezuela"
replace nCH = 1 if gdename == "toronto (cdn)"
replace nCH = 1 if gdename == "new york ny usa"
replace nCH = 1 if gdename == "calif. usa"
replace nCH = 1 if gdename == "reichenau d"
replace nCH = 1 if gdename == "kuwait"
replace nCH = 1 if gdename == "luedenscheid d"
replace nCH = 1 if gdename == "helsinki fin"
replace nCH = 1 if gdename == "linz (a)"
replace nCH = 1 if gdename == "evian f"
replace nCH = 1 if gdename == "kempten d"
replace nCH = 1 if gdename == "rhode st genese b"
replace nCH = 1 if gdename == "tripoli lar"
replace nCH = 1 if gdename == "uppsala s"
replace nCH = 1 if gdename == "muehlheim d"
replace nCH = 1 if gdename == "kaarst d"
replace nCH = 1 if gdename == "lauchringen d"
replace nCH = 1 if gdename == "mettmann d"
replace nCH = 1 if gdename == "antwerpen (b)"
replace nCH = 1 if gdename == "vesinet f"
replace nCH = 1 if gdename == "grenzach d"
replace nCH = 1 if gdename == "leymen f"
replace nCH = 1 if gdename == "aix en provence f"
replace nCH = 1 if gdename == "wisconsin usa"
replace nCH = 1 if gdename == "kuwait kwt"
replace nCH = 1 if gdename == "annecy le vieux f"
replace nCH = 1 if gdename == "osnabrueck d"
replace nCH = 1 if gdename == "villette"
replace nCH = 1 if gdename == "hertfordshire gb"
replace nCH = 1 if gdename == "lagos wan"
replace nCH = 1 if gdename == "lottstetten d"
replace nCH = 1 if gdename == "radolfzell d"
replace nCH = 1 if gdename == "buesingen d"
replace nCH = 1 if gdename == "pittsburgh usa"
replace nCH = 1 if gdename == "schramberg d"
replace nCH = 1 if gdename == "arizona usa"
replace nCH = 1 if gdename == "genes i"
replace nCH = 1 if gdename == "windsor gb"
replace nCH = 1 if gdename == "kopenhagen (dk)"
replace nCH = 1 if gdename == "campione (i)"
replace nCH = 1 if gdename == "taipei (taiwan. r.o.c.)"
replace nCH = 1 if gdename == "miami"
replace nCH = 1 if gdename == "friedberg d"
replace nCH = 1 if gdename == "buckinghamshire gb"
replace nCH = 1 if gdename == "safat kwt"
replace nCH = 1 if gdename == "nassau bs"
replace nCH = 1 if gdename == "chicago (usa)"
replace nCH = 1 if gdename == "steinen d"
replace nCH = 1 if gdename == "reims f"
replace nCH = 1 if gdename == "kanagawa j"
replace nCH = 1 if gdename == "suresnes f"
replace nCH = 1 if gdename == "johannesburg (za)"
replace nCH = 1 if gdename == "bielefeld (d)"
replace nCH = 1 if gdename == "versailles (f)"
replace nCH = 1 if gdename == "gottmadingen d"
replace nCH = 1 if gdename == "krefeld (d)"
replace nCH = 1 if gdename == "laufenburg d"
replace nCH = 1 if gdename == "stamford usa"
replace nCH = 1 if gdename == "javea e"
replace nCH = 1 if gdename == "worms d"
replace nCH = 1 if gdename == "ohio (usa)"
replace nCH = 1 if gdename == "beaconsfield gb"
replace nCH = 1 if gdename == "biberach d"
replace nCH = 1 if gdename == "indiana usa"
replace nCH = 1 if gdename == "bertrange l"
replace nCH = 1 if gdename == "karlsruhe (d)"
replace nCH = 1 if gdename == "pulheim d"
replace nCH = 1 if gdename == "catania i"
replace nCH = 1 if gdename == "driebergen nl"
replace nCH = 1 if gdename == "california (usa)"
replace nCH = 1 if gdename == "houston (usa)"
replace nCH = 1 if gdename == "aalen d"
replace nCH = 1 if gdename == "oberkirch d"
replace nCH = 1 if gdename == "esslingen d"
replace nCH = 1 if gdename == "moenchengladbach d"
replace nCH = 1 if gdename == "dortmund (d)"
replace nCH = 1 if gdename == "sceaux f"
replace nCH = 1 if gdename == "moskau"
replace nCH = 1 if gdename == "belgrad (yu)"
replace nCH = 1 if gdename == "oss nl"
replace nCH = 1 if gdename == "rankweil a"
replace nCH = 1 if gdename == "singapore sgp"
replace nCH = 1 if gdename == "jestetten d"
replace nCH = 1 if gdename == "essex (gb)"
replace nCH = 1 if gdename == "connecticut (usa)"
replace nCH = 1 if gdename == "salzburg (a)"
replace nCH = 1 if gdename == "prag cz"
replace nCH = 1 if gdename == "tervuren b"
replace nCH = 1 if gdename == "montmorency f"
replace nCH = 1 if gdename == "modena i"
replace nCH = 1 if gdename == "riad (sa)"
replace nCH = 1 if gdename == "lima pe"
replace nCH = 1 if gdename == "harpenden gb"
replace nCH = 1 if gdename == "kent (gb)"
replace nCH = 1 if gdename == "sevres f"
replace nCH = 1 if gdename == "le vesinet (f)"
replace nCH = 1 if gdename == "gentofte dk"
replace nCH = 1 if gdename == "lagos (wan)"
replace nCH = 1 if gdename == "texas (usa)"
replace nCH = 1 if gdename == "haan d"
replace nCH = 1 if gdename == "nogent sur marne f"
replace nCH = 1 if gdename == "michigan (usa)"
replace nCH = 1 if gdename == "rixensart b"
replace nCH = 1 if gdename == "san diego usa"
replace nCH = 1 if gdename == "bloemendaal nl"
replace nCH = 1 if gdename == "colorado usa"
replace nCH = 1 if gdename == "kapellen b"
replace nCH = 1 if gdename == "wutoeschingen d"
replace nCH = 1 if gdename == "cannes f"
replace nCH = 1 if gdename == "dreieich d"
replace nCH = 1 if gdename == "ascot gb"
replace nCH = 1 if gdename == "neu isenburg d"
replace nCH = 1 if gdename == "singapur"
replace nCH = 1 if gdename == "viry f"
replace nCH = 1 if gdename == "channel"
replace nCH = 1 if gdename == "den haag (nl)"
replace nCH = 1 if gdename == "monte carlo (monaco)"
replace nCH = 1 if gdename == "bochum (d)"
replace nCH = 1 if gdename == "wayne usa"
replace nCH = 1 if gdename == "gb"
replace nCH = 1 if gdename == "amman hkj"
replace nCH = 1 if gdename == "emsdetten d"
replace nCH = 1 if gdename == "bermuda"
replace nCH = 1 if gdename == "danderyd s"
replace nCH = 1 if gdename == "dornbirn (a)"
replace nCH = 1 if gdename == "dijon f"
replace nCH = 1 if gdename == "jeddah (sa)"
replace nCH = 1 if gdename == "sark"
replace nCH = 1 if gdename == "cremona i"
replace nCH = 1 if gdename == "koblenz d"
replace nCH = 1 if gdename == "goettingen d"
replace nCH = 1 if gdename == "bottrop d"
replace nCH = 1 if gdename == "pretoria za"
replace nCH = 1 if gdename == "kuweit (kwt)"
replace nCH = 1 if gdename == "regensburg d"
replace nCH = 1 if gdename == "neustadt d"
replace nCH = 1 if gdename == "west sussex gb"
replace nCH = 1 if gdename == "frankfurt a. main d"
replace nCH = 1 if gdename == "erba i"
replace nCH = 1 if gdename == "croix f"
replace nCH = 1 if gdename == "cranves sales f"
replace nCH = 1 if gdename == "bad homburg (d)"
replace nCH = 1 if gdename == "aschaffenburg d"
replace nCH = 1 if gdename == "waldshut d"
replace nCH = 1 if gdename == "muenchen dl"
replace nCH = 1 if gdename == "leinfelden echterdingen d"
replace nCH = 1 if gdename == "bregenz (a)"
replace nCH = 1 if gdename == "ludwigshafen d"
replace nCH = 1 if gdename == "jeddah ksa"
replace nCH = 1 if gdename == "hilversum (nl)"
replace nCH = 1 if gdename == "schilde b"
replace nCH = 1 if gdename == "ludwigshafen (d)"
replace nCH = 1 if gdename == "taunusstein d"
replace nCH = 1 if gdename == "antony f"
replace nCH = 1 if gdename == "billund dk"
replace nCH = 1 if gdename == "loerrach (d)"
replace nCH = 1 if gdename == "mouvaux f"
replace nCH = 1 if gdename == "maidenhead gb"
replace nCH = 1 if gdename == "kairo"
replace nCH = 1 if gdename == "n.j. usa"
replace nCH = 1 if gdename == "bensheim d"
replace nCH = 1 if gdename == "reignier f"
replace nCH = 1 if gdename == "nottingham gb"
replace nCH = 1 if gdename == "neu isenburg (d)"
replace nCH = 1 if gdename == "woking"
replace nCH = 1 if gdename == "inzlingen d"
replace nCH = 1 if gdename == "ma usa"
replace nCH = 1 if gdename == "beverly hills usa"
replace nCH = 1 if gdename == "linkoeping s"
replace nCH = 1 if gdename == "schoten b"
replace nCH = 1 if gdename == "atherton usa"
replace nCH = 1 if gdename == "porlezza i"
replace nCH = 1 if gdename == "fontainebleau f"
replace nCH = 1 if gdename == "ornex f"
replace nCH = 1 if gdename == "tuttlingen d"
replace nCH = 1 if gdename == "braunschweig d"
replace nCH = 1 if gdename == "velbert d"
replace nCH = 1 if gdename == "freiburg (d)"
replace nCH = 1 if gdename == "philadelphia (usa)"
replace nCH = 1 if gdename == "lille f"
replace nCH = 1 if gdename == "northwood gb"
replace nCH = 1 if gdename == "aoste i"
replace nCH = 1 if gdename == "wichita usa"
replace nCH = 1 if gdename == "oldenburg d"
replace nCH = 1 if gdename == "bridel l"
replace nCH = 1 if gdename == "muellheim d"
replace nCH = 1 if gdename == "lissabon (p)"
replace nCH = 1 if gdename == "sandton za"
replace nCH = 1 if gdename == "bilbao e"
replace nCH = 1 if gdename == "koenigstein im taunus d"
replace nCH = 1 if gdename == "le pecq f"
replace nCH = 1 if gdename == "muelheim d"
replace nCH = 1 if gdename == "pa usa"
replace nCH = 1 if gdename == "keerbergen b"
replace nCH = 1 if gdename == "maisons laffitte f"
replace nCH = 1 if gdename == "pforzheim (d)"
replace nCH = 1 if gdename == "tourcoing f"
replace nCH = 1 if gdename == "goetzis a"
replace nCH = 1 if gdename == "graefelfing d"
replace nCH = 1 if gdename == "annemasse (f)"
replace nCH = 1 if gdename == "bonn (d)"
replace nCH = 1 if gdename == "morteau f"
replace nCH = 1 if gdename == "regensburg (d)"
replace nCH = 1 if gdename == "heusenstamm d"
replace nCH = 1 if gdename == "n.y. usa"
replace nCH = 1 if gdename == "montreal (cdn)"
replace nCH = 1 if gdename == "nassau bahamas"
replace nCH = 1 if gdename == "buehl d"
replace nCH = 1 if gdename == "palma de mallorca e"
replace nCH = 1 if gdename == "rungsted kyst dk"
replace nCH = 1 if gdename == "nice (f)"
replace nCH = 1 if gdename == "mauren (fl)"
replace nCH = 1 if gdename == "braine l alleud b"
replace nCH = 1 if gdename == "toronto can"
replace nCH = 1 if gdename == "nagano j"
replace nCH = 1 if gdename == "heemstede nl"
replace nCH = 1 if gdename == "offenburg d"
replace nCH = 1 if gdename == "tutzing d"
replace nCH = 1 if gdename == "schwanden"
replace nCH = 1 if gdename == "calgary cdn"
replace nCH = 1 if gdename == "bristol gb"
replace nCH = 1 if gdename == "liphook gb"
replace nCH = 1 if gdename == "minneapolis usa"
replace nCH = 1 if gdename == "hockenheim d"
replace nCH = 1 if gdename == "leonding a"
replace nCH = 1 if gdename == "gruenwald (d)"
replace nCH = 1 if gdename == "st didier au mont d or f"
replace nCH = 1 if gdename == "princeton usa"
replace nCH = 1 if gdename == "montevideo rou"
replace nCH = 1 if gdename == "milwaukee usa"
replace nCH = 1 if gdename == "san jose cr"
replace nCH = 1 if gdename == "bossey f"
replace nCH = 1 if gdename == "meerbusch (d)"
replace nCH = 1 if gdename == "lahti sf"
replace nCH = 1 if gdename == "seoul rok"
replace nCH = 1 if gdename == "langen d"
replace nCH = 1 if gdename == "moedling a"
replace nCH = 1 if gdename == "apeldoorn nl"
replace nCH = 1 if gdename == "ulm (d)"
replace nCH = 1 if gdename == "aarhus dk"
replace nCH = 1 if gdename == "singapour sgp"
replace nCH = 1 if gdename == "grenoble (f)"
replace nCH = 1 if gdename == "vincennes f"
replace nCH = 1 if gdename == "riedisheim f"
replace nCH = 1 if gdename == "greenwich (usa)"
replace nCH = 1 if gdename == "kehl d"
replace nCH = 1 if gdename == "saint julien en genevois f"
replace nCH = 1 if gdename == "langenfeld d"
replace nCH = 1 if gdename == "amstelveen (nl)"
replace nCH = 1 if gdename == "minden d"
replace nCH = 1 if gdename == "minnesota usa"
replace nCH = 1 if gdename == "lake forest usa"
replace nCH = 1 if gdename == "levallois perret f"
replace nCH = 1 if gdename == "ottobrunn d"
replace nCH = 1 if gdename == "darien usa"
replace nCH = 1 if gdename == "ca"
replace nCH = 1 if gdename == "knittlingen d"
replace nCH = 1 if gdename == "clermont ferrand f"
replace nCH = 1 if gdename == "dietzenbach d"
replace nCH = 1 if gdename == "nuertingen d"
replace nCH = 1 if gdename == "roma (i)"
replace nCH = 1 if gdename == "weston usa"
replace nCH = 1 if gdename == "walldorf d"
replace nCH = 1 if gdename == "palo alto usa"
replace nCH = 1 if gdename == "piacenza i"
replace nCH = 1 if gdename == "8645jona"
replace nCH = 1 if gdename == "saeroe s"
replace nCH = 1 if gdename == "kandern d"
replace nCH = 1 if gdename == "luxembourg lux"
replace nCH = 1 if gdename == "richmond gb"
replace nCH = 1 if gdename == "koenigswinter d"
replace nCH = 1 if gdename == "moers d"
replace nCH = 1 if gdename == "hanau d"
replace nCH = 1 if gdename == "manchester gb"
replace nCH = 1 if gdename == "ohain b"
replace nCH = 1 if gdename == "binzen d"
replace nCH = 1 if gdename == "triest i"
replace nCH = 1 if gdename == "moskau russland"
replace nCH = 1 if gdename == "fillinges f"
replace nCH = 1 if gdename == "tampere sf"
replace nCH = 1 if gdename == "gerrards cross gb"
replace nCH = 1 if gdename == "abidjan ci"
replace nCH = 1 if gdename == "utrecht nl"
replace nCH = 1 if gdename == "iserlohn d"
replace nCH = 1 if gdename == "cannobio i"
replace nCH = 1 if gdename == "karben d"
replace nCH = 1 if gdename == "chens sur leman f"
replace nCH = 1 if gdename == "strassburg (f)"
replace nCH = 1 if gdename == "biella i"
replace nCH = 1 if gdename == "skopje mk"
replace nCH = 1 if gdename == "korntal muenchingen d"
replace nCH = 1 if gdename == "weil am rhein (d)"
replace nCH = 1 if gdename == "ankara tr"
replace nCH = 1 if gdename == "moskau (russland)"
replace nCH = 1 if gdename == "postfach"
replace nCH = 1 if gdename == "kronberg (d)"
replace nCH = 1 if gdename == "congleton gb"
replace nCH = 1 if gdename == "melbourne aus"
replace nCH = 1 if gdename == "offenbach d"
replace nCH = 1 if gdename == "dubai uae"
replace nCH = 1 if gdename == "casale monferrato i"
replace nCH = 1 if gdename == "edegem b"
replace nCH = 1 if gdename == "meudon (f)"
replace nCH = 1 if gdename == "uden nl"
replace nCH = 1 if gdename == "breda nl"
replace nCH = 1 if gdename == "villingen schwenningen d"
replace nCH = 1 if gdename == "mumbai ind"
replace nCH = 1 if gdename == "frankfurt/main d"
replace nCH = 1 if gdename == "gerlingen d"
replace nCH = 1 if gdename == "segrate i"
replace nCH = 1 if gdename == "kansas usa"
replace nCH = 1 if gdename == "californie usa"
replace nCH = 1 if gdename == "wuppertal (d)"
replace nCH = 1 if gdename == "london uk"
replace nCH = 1 if gdename == "la"
replace nCH = 1 if gdename == "karlstad s"
replace nCH = 1 if gdename == "ruesselsheim d"
replace nCH = 1 if gdename == "sassuolo i"
replace nCH = 1 if gdename == "biberach an der riss d"
replace nCH = 1 if gdename == "lund (s)"
replace nCH = 1 if gdename == "taipei rc"
replace nCH = 1 if gdename == "baarn nl"
replace nCH = 1 if gdename == "purley gb"
replace nCH = 1 if gdename == "missouri usa"
replace nCH = 1 if gdename == "maintal d"
replace nCH = 1 if gdename == "vedbaek dk"
replace nCH = 1 if gdename == "singen (d)"
replace nCH = 1 if gdename == "murg d"
replace nCH = 1 if gdename == "rixheim f"
replace nCH = 1 if gdename == "n.j. (usa)"
replace nCH = 1 if gdename == "wiehl d"
replace nCH = 1 if gdename == "vaesteras s"
replace nCH = 1 if gdename == "baldham d"
replace nCH = 1 if gdename == "saratoga"
replace nCH = 1 if gdename == "bratislava sk"
replace nCH = 1 if gdename == "veyrier du lac f"
replace nCH = 1 if gdename == "karachi pak"
replace nCH = 1 if gdename == "nancy f"
replace nCH = 1 if gdename == "meylan f"
replace nCH = 1 if gdename == "ca usa"
replace nCH = 1 if gdename == "haslemere gb"
replace nCH = 1 if gdename == "ridgefield usa"
replace nCH = 1 if gdename == "holte (dk)"
replace nCH = 1 if gdename == "cascais p"
replace nCH = 1 if gdename == "tripoli lb"
replace nCH = 1 if gdename == "la celle st cloud f"
replace nCH = 1 if gdename == "hamburg dl"
replace nCH = 1 if gdename == "sao paulo (br)"
replace nCH = 1 if gdename == "gouda nl"
replace nCH = 1 if gdename == "pontarlier f"
replace nCH = 1 if gdename == "emmendingen d"
replace nCH = 1 if gdename == "kuessaberg d"
replace nCH = 1 if gdename == "arese i"
replace nCH = 1 if gdename == "berkshire (gb)"
replace nCH = 1 if gdename == "wehr d"
replace nCH = 1 if gdename == "conn. (usa)"
replace nCH = 1 if gdename == "waterloo (b)"
replace nCH = 1 if gdename == "walferdange l"
replace nCH = 1 if gdename == "oegstgeest nl"
replace nCH = 1 if gdename == "bad vilbel d"
replace nCH = 1 if gdename == "oerebro s"
replace nCH = 1 if gdename == "wilrijk b"
replace nCH = 1 if gdename == "phoenix usa"
replace nCH = 1 if gdename == "pino torinese i"
replace nCH = 1 if gdename == "dallas (usa)"
replace nCH = 1 if gdename == "colombo cl"
replace nCH = 1 if gdename == "aerdenhout (nl)"
replace nCH = 1 if gdename == "versonnex f"
replace nCH = 1 if gdename == "bombay (ind)"
replace nCH = 1 if gdename == "wavre b"
replace nCH = 1 if gdename == "riyad sa"
replace nCH = 1 if gdename == "herts gb"
replace nCH = 1 if gdename == "muenchen/d"
replace nCH = 1 if gdename == "birkenfeld d"
replace nCH = 1 if gdename == "zagreb hr"
replace nCH = 1 if gdename == "lima peru"
replace nCH = 1 if gdename == "caracas (yv)"
replace nCH = 1 if gdename == "berchem b"
replace nCH = 1 if gdename == "wentzwiller f"
replace nCH = 1 if gdename == "muenster"
replace nCH = 1 if gdename == "schwenningen d"
replace nCH = 1 if gdename == "perrignier f"
replace nCH = 1 if gdename == "mass. (usa)"
replace nCH = 1 if gdename == "hamilton (bermudes)"
replace nCH = 1 if gdename == "frankreich"
replace nCH = 1 if gdename == "antibes f"
replace nCH = 1 if gdename == "lohmar d"
replace nCH = 1 if gdename == "vaexjoe s"
replace nCH = 1 if gdename == "koenigstein (d)"
replace nCH = 1 if gdename == "sigmaringen d"
replace nCH = 1 if gdename == "johannesburg (south africa)"
replace nCH = 1 if gdename == "enschede nl"
replace nCH = 1 if gdename == "frankfurt a. m. (d)"
replace nCH = 1 if gdename == "panama city pa"
replace nCH = 1 if gdename == "ingelfingen d"
replace nCH = 1 if gdename == "delle f"
replace nCH = 1 if gdename == "woking gb"
replace nCH = 1 if gdename == "massachusetts (usa)"
replace nCH = 1 if gdename == "nuenen nl"
replace nCH = 1 if gdename == "torino (i)"
replace nCH = 1 if gdename == "cal. (usa)"
replace nCH = 1 if gdename == "colombes f"
replace nCH = 1 if gdename == "friedrichsdorf d"
replace nCH = 1 if gdename == "saint germain en laye f"
replace nCH = 1 if gdename == "blaricum (nl)"
replace nCH = 1 if gdename == "carimate i"
replace nCH = 1 if gdename == "lindenberg d"
replace nCH = 1 if gdename == "bludenz a"
replace nCH = 1 if gdename == "doha"
replace nCH = 1 if gdename == "huddinge s"
replace nCH = 1 if gdename == "trier d"
replace nCH = 1 if gdename == "valencia e"
replace nCH = 1 if gdename == "darmstadt (d)"
replace nCH = 1 if gdename == "roesrath d"
replace nCH = 1 if gdename == "pullach d"
replace nCH = 1 if gdename == "weinstadt d"
replace nCH = 1 if gdename == "marburg d"
replace nCH = 1 if gdename == "kuwait (kwt)"
replace nCH = 1 if gdename == "idstein d"
replace nCH = 1 if gdename == "wattens a"
replace nCH = 1 if gdename == "californien usa"
replace nCH = 1 if gdename == "wokingham gb"
replace nCH = 1 if gdename == "fort worth"
replace nCH = 1 if gdename == "houston texas usa"
replace nCH = 1 if gdename == "landshut d"
replace nCH = 1 if gdename == "frechen d"
replace nCH = 1 if gdename == "kelkheim d"
replace nCH = 1 if gdename == "new hampshire usa"
replace nCH = 1 if gdename == "colombo (sri lanka)"
replace nCH = 1 if gdename == "beyrouth (rl)"
replace nCH = 1 if gdename == "englewood"
replace nCH = 1 if gdename == "kirchheim teck d"
replace nCH = 1 if gdename == "genes"
replace nCH = 1 if gdename == "filderstadt d"
replace nCH = 1 if gdename == "vancouver (cdn)"
replace nCH = 1 if gdename == "espelkamp d"
replace nCH = 1 if gdename == "munich d"
replace nCH = 1 if gdename == "montecarlo mc"
replace nCH = 1 if gdename == "bergen n"
replace nCH = 1 if gdename == "leeds gb"
replace nCH = 1 if gdename == "rio de janeiro (br)"
replace nCH = 1 if gdename == "st nom la breteche f"
replace nCH = 1 if gdename == "northampton gb"
replace nCH = 1 if gdename == "unna d"
replace nCH = 1 if gdename == "winchester gb"
replace nCH = 1 if gdename == "mierlo nl"
replace nCH = 1 if gdename == "klampenborg dk"
replace nCH = 1 if gdename == "cruseilles f"
replace nCH = 1 if gdename == "gladbeck d"
replace nCH = 1 if gdename == "burgwedel d"
replace nCH = 1 if gdename == "croissy sur seine f"
replace nCH = 1 if gdename == "haarlem nl"
replace nCH = 1 if gdename == "san francisco (usa)"
replace nCH = 1 if gdename == "le chesnay f"
replace nCH = 1 if gdename == "boeblingen d"
replace nCH = 1 if gdename == "sussex gb"
replace nCH = 1 if gdename == "le caire et"
replace nCH = 1 if gdename == "hofheim d"
replace nCH = 1 if gdename == "neapel i"
replace nCH = 1 if gdename == "sao paolo br"
replace nCH = 1 if gdename == "spartanburg usa"
replace nCH = 1 if gdename == "arcore i"
replace nCH = 1 if gdename == "vienne f"
replace nCH = 1 if gdename == "sausheim f"
replace nCH = 1 if gdename == "greenwich connecticut usa"
replace nCH = 1 if gdename == "baden baden (d)"
replace nCH = 1 if gdename == "arnhem nl"
replace nCH = 1 if gdename == "veigy f"
replace nCH = 1 if gdename == "koweit kwt"
replace nCH = 1 if gdename == "hilden d"
replace nCH = 1 if gdename == "bad krozingen d"
replace nCH = 1 if gdename == "bromma (s)"
replace nCH = 1 if gdename == "ostfildern d"
replace nCH = 1 if gdename == "montgomery usa"
replace nCH = 1 if gdename == "atlanta georgia usa"
replace nCH = 1 if gdename == "guildford gb"
replace nCH = 1 if gdename == "malland i"
replace nCH = 1 if gdename == "reggio emilia i"
replace nCH = 1 if gdename == "wolfsburg d"
replace nCH = 1 if gdename == "jaerfaella s"
replace nCH = 1 if gdename == "reutlingen (d)"
replace nCH = 1 if gdename == "florenz (i)"
replace nCH = 1 if gdename == "melle d"
replace nCH = 1 if gdename == "mainz (d)"
replace nCH = 1 if gdename == "kalmthout b"
replace nCH = 1 if gdename == "schorndorf d"
replace nCH = 1 if gdename == "bologne i"
replace nCH = 1 if gdename == "n.y. (usa)"
replace nCH = 1 if gdename == "ivrea i"
replace nCH = 1 if gdename == "bailly f"
replace nCH = 1 if gdename == "montpellier f"
replace nCH = 1 if gdename == "rheine d"
replace nCH = 1 if gdename == "nj usa"
replace nCH = 1 if gdename == "weingarten d"
replace nCH = 1 if gdename == "puteaux f"
replace nCH = 1 if gdename == "steinsel l"
replace nCH = 1 if gdename == "niefern oeschelbronn d"
replace nCH = 1 if gdename == "lampertheim d"
replace nCH = 1 if gdename == "basiglio i"
replace nCH = 1 if gdename == "hattersheim d"
replace nCH = 1 if gdename == "le havre f"
replace nCH = 1 if gdename == "sussex (gb)"
replace nCH = 1 if gdename == "alabama usa"
replace nCH = 1 if gdename == "ferney voltaire (f)"
replace nCH = 1 if gdename == "neu anspach d"
replace nCH = 1 if gdename == "reinbek d"
replace nCH = 1 if gdename == "damas syr"
replace nCH = 1 if gdename == "waldshut (d)"
replace nCH = 1 if gdename == "st cergues f"
replace nCH = 1 if gdename == "menden d"
replace nCH = 1 if gdename == "obertshausen d"
replace nCH = 1 if gdename == "audincourt f"
replace nCH = 1 if gdename == "bordeaux (f)"
replace nCH = 1 if gdename == "arnsberg d"
replace nCH = 1 if gdename == "viersen d"
replace nCH = 1 if gdename == "feldkirch (a)"
replace nCH = 1 if gdename == "wasserburg d"
replace nCH = 1 if gdename == "ny usa"
replace nCH = 1 if gdename == "waldkirch d"
replace nCH = 1 if gdename == "cobham gb"
replace nCH = 1 if gdename == "esher gb"
replace nCH = 1 if gdename == "issy les moulineaux f"
replace nCH = 1 if gdename == "laren (nl)"
replace nCH = 1 if gdename == "giessen d"
replace nCH = 1 if gdename == "amstetten a"
replace nCH = 1 if gdename == "calgary"
replace nCH = 1 if gdename == "reeuwijk nl"
replace nCH = 1 if gdename == "moelndal s"
replace nCH = 1 if gdename == "villach a"
replace nCH = 1 if gdename == "monheim d"
replace nCH = 1 if gdename == "metzingen d"
replace nCH = 1 if gdename == "plano"
replace nCH = 1 if gdename == "strassbourg f"
replace nCH = 1 if gdename == "monrovia lb"
replace nCH = 1 if gdename == "wolfurt a"
replace nCH = 1 if gdename == "iffezheim d"
replace nCH = 1 if gdename == "fabriano i"
replace nCH = 1 if gdename == "meersburg d"
replace nCH = 1 if gdename == "henley on thames gb"
replace nCH = 1 if gdename == "aix les bains f"
replace nCH = 1 if gdename == "lagos nigeria"
replace nCH = 1 if gdename == "clamart f"
replace nCH = 1 if gdename == "monza (i)"
replace nCH = 1 if gdename == "lidingoe (s)"
replace nCH = 1 if gdename == "graz (a)"
replace nCH = 1 if gdename == "vicenza (i)"
replace nCH = 1 if gdename == "ibiza e"
replace nCH = 1 if gdename == "weilheim d"
replace nCH = 1 if gdename == "pennsylvania (usa)"
replace nCH = 1 if gdename == "echenevex f"
replace nCH = 1 if gdename == "thessaloniki gr"
replace nCH = 1 if gdename == "raedersdorf f"
replace nCH = 1 if gdename == "bad nauheim d"
replace nCH = 1 if gdename == "bogota co"
replace nCH = 1 if gdename == "annemasse/f"
replace nCH = 1 if gdename == "meran i"
replace nCH = 1 if gdename == "kalifornien (usa)"
replace nCH = 1 if gdename == "monza"
replace nCH = 1 if gdename == "oxted gb"
replace nCH = 1 if gdename == "bocholt d"
replace nCH = 1 if gdename == "charlotte usa"
replace nCH = 1 if gdename == "collonges sous saleve (f)"
replace nCH = 1 if gdename == "sierentz f"
replace nCH = 1 if gdename == "tunbridge wells gb"
replace nCH = 1 if gdename == "ennepetal d"
replace nCH = 1 if gdename == "fontenay sous bois f"
replace nCH = 1 if gdename == "walton on thames gb"
replace nCH = 1 if gdename == "kolding dk"
replace nCH = 1 if gdename == "stuttgart dl"
replace nCH = 1 if gdename == "schwetzingen d"
replace nCH = 1 if gdename == "grez doiceau b"
replace nCH = 1 if gdename == "bons en chablais f"
replace nCH = 1 if gdename == "denver usa"
replace nCH = 1 if gdename == "san sebastian e"
replace nCH = 1 if gdename == "herrenberg d"
replace nCH = 1 if gdename == "taeby (s)"
replace nCH = 1 if gdename == "kobe j"
replace nCH = 1 if gdename == "innsbruck (a)"
replace nCH = 1 if gdename == "kristiansand n"
replace nCH = 1 if gdename == "zaventem b"
replace nCH = 1 if gdename == "venezia i"
replace nCH = 1 if gdename == "nijmegen nl"
replace nCH = 1 if gdename == "zagreb kroatien"
replace nCH = 1 if gdename == "klagenfurt (a)"
replace nCH = 1 if gdename == "geislingen d"
replace nCH = 1 if gdename == "idar oberstein d"
replace nCH = 1 if gdename == "albstadt d"
replace nCH = 1 if gdename == "buende d"
replace nCH = 1 if gdename == "hongkong (hk)"
replace nCH = 1 if gdename == "korschenbroich d"
replace nCH = 1 if gdename == "eimeldingen d"
replace nCH = 1 if gdename == "beckenham gb"
replace nCH = 1 if gdename == "schwaebisch gmuend d"
replace nCH = 1 if gdename == "zeist (nl)"
replace nCH = 1 if gdename == "woerden nl"
replace nCH = 1 if gdename == "moenchengladbach (d)"
replace nCH = 1 if gdename == "heiloo nl"
replace nCH = 1 if gdename == "warmond nl"
replace nCH = 1 if gdename == "saint cloud (f)"
replace nCH = 1 if gdename == "andover"
replace nCH = 1 if gdename == "rambouillet f"
replace nCH = 1 if gdename == "tettnang d"
replace nCH = 1 if gdename == "ljubljana yu"
replace nCH = 1 if gdename == "novara (i)"
replace nCH = 1 if gdename == "lexington"
replace nCH = 1 if gdename == "guetersloh d"
replace nCH = 1 if gdename == "coventry gb"
replace nCH = 1 if gdename == "miltenberg a.m. d"
replace nCH = 1 if gdename == "bennekom nl"
replace nCH = 1 if gdename == "saitama j"
replace nCH = 1 if gdename == "eschborn d"
replace nCH = 1 if gdename == "s. donato milanese i"
replace nCH = 1 if gdename == "york"
replace nCH = 1 if gdename == "ambilly (f)"
replace nCH = 1 if gdename == "marly le roi f"
replace nCH = 1 if gdename == "oosterhout nl"
replace nCH = 1 if gdename == "hemel hempstead"
replace nCH = 1 if gdename == "bethersden ashford/kent gb"
replace nCH = 1 if gdename == "bruehl d"
replace nCH = 1 if gdename == "mexico mex"
replace nCH = 1 if gdename == "cape town za"
replace nCH = 1 if gdename == "tokyo japan"
replace nCH = 1 if gdename == "trento i"
replace nCH = 1 if gdename == "douala"
replace nCH = 1 if gdename == "kulmbach d"
replace nCH = 1 if gdename == "fellbach d"
replace nCH = 1 if gdename == "chambourcy f"
replace nCH = 1 if gdename == "rottenburg d"
replace nCH = 1 if gdename == "straubenhardt d"
replace nCH = 1 if gdename == "berkeley usa"
replace nCH = 1 if gdename == "vejle dk"
replace nCH = 1 if gdename == "goeppingen (d)"
replace nCH = 1 if gdename == "bergamo (i)"
replace nCH = 1 if gdename == "rhode island usa"
replace nCH = 1 if gdename == "paderborn (d)"
replace nCH = 1 if gdename == "annecy (f)"
replace nCH = 1 if gdename == "oelde d"
replace nCH = 1 if gdename == "mailand/i"
replace nCH = 1 if gdename == "meise b"
replace nCH = 1 if gdename == "collonge sous saleve f"
replace nCH = 1 if gdename == "feldkirchen d"
replace nCH = 1 if gdename == "herblay f"
replace nCH = 1 if gdename == "erlenbach"
replace nCH = 1 if gdename == "lauterach (a)"
replace nCH = 1 if gdename == "maulburg d"
replace nCH = 1 if gdename == "virum dk"
replace nCH = 1 if gdename == "short hills usa"
replace nCH = 1 if gdename == "simferopol ukraine"
replace nCH = 1 if gdename == "murnau d"
replace nCH = 1 if gdename == "champagnole f"
replace nCH = 1 if gdename == "kamerun"
replace nCH = 1 if gdename == "bad duerkheim d"
replace nCH = 1 if gdename == "donaueschingen d"
replace nCH = 1 if gdename == "weil a. rh. d"
replace nCH = 1 if gdename == "hamburg/d"
replace nCH = 1 if gdename == "st. peter port"
replace nCH = 1 if gdename == "augsburg (d)"
replace nCH = 1 if gdename == "sandton rsa"
replace nCH = 1 if gdename == "fulda (d)"
replace nCH = 1 if gdename == "fla usa"
replace nCH = 1 if gdename == "busto arsizio i"
replace nCH = 1 if gdename == "glashuetten d"
replace nCH = 1 if gdename == "hagenthal le bas f"
replace nCH = 1 if gdename == "muenster (d)"
replace nCH = 1 if gdename == "denia e"
replace nCH = 1 if gdename == "seeheim jugenheim d"
replace nCH = 1 if gdename == "mte carlo"
replace nCH = 1 if gdename == "breisach d"
replace nCH = 1 if gdename == "cluses f"
replace nCH = 1 if gdename == "tucson"
replace nCH = 1 if gdename == "quebec ca"
replace nCH = 1 if gdename == "nevada usa"
replace nCH = 1 if gdename == "kassel (d)"
replace nCH = 1 if gdename == "djibouti"
replace nCH = 1 if gdename == "tuttlingen (d)"
replace nCH = 1 if gdename == "beograd yu"
replace nCH = 1 if gdename == "gutach (d)"
replace nCH = 1 if gdename == "zandvoort nl"
replace nCH = 1 if gdename == "holzkirchen d"
replace nCH = 1 if gdename == "bethesda"
replace nCH = 1 if gdename == "aubagne f"
replace nCH = 1 if gdename == "abu dhabi"
replace nCH = 1 if gdename == "kirn d"
replace nCH = 1 if gdename == "wimbledon gb"
replace nCH = 1 if gdename == "pfullingen d"
replace nCH = 1 if gdename == "gent b"
replace nCH = 1 if gdename == "berkhamsted gb"
replace nCH = 1 if gdename == "alkmaar nl"
replace nCH = 1 if gdename == "cannes (f)"
replace nCH = 1 if gdename == "espoo fin"
replace nCH = 1 if gdename == "bietigheim bissingen d"
replace nCH = 1 if gdename == "uccle (b)"
replace nCH = 1 if gdename == "stockach d"
replace nCH = 1 if gdename == "tel aviv israel"
replace nCH = 1 if gdename == "mechelen b"
replace nCH = 1 if gdename == "brookfield usa"
replace nCH = 1 if gdename == "kuala lumpur mal"
replace nCH = 1 if gdename == "quai"
replace nCH = 1 if gdename == "esch sur alzette l"
replace nCH = 1 if gdename == "varese (i)"
replace nCH = 1 if gdename == "worms (d)"
replace nCH = 1 if gdename == "bennebroek nl"
replace nCH = 1 if gdename == "danville usa"
replace nCH = 1 if gdename == "lochau a"
replace nCH = 1 if gdename == "ponte tresa i"
replace nCH = 1 if gdename == "bussum (nl)"
replace nCH = 1 if gdename == "djeddah sa"
replace nCH = 1 if gdename == "heidenheim d"
replace nCH = 1 if gdename == "basingstoke gb"
replace nCH = 1 if gdename == "chichester gb"
replace nCH = 1 if gdename == "fort worth usa"
replace nCH = 1 if gdename == "guernsey"
replace nCH = 1 if gdename == "gehrden d"
replace nCH = 1 if gdename == "andorre and"
replace nCH = 1 if gdename == "dormagen d"
replace nCH = 1 if gdename == "boca raton usa"
replace nCH = 1 if gdename == "warwickshire gb"
replace nCH = 1 if gdename == "riad sa"
replace nCH = 1 if gdename == "gailingen d"
replace nCH = 1 if gdename == "hoquiam"
replace nCH = 1 if gdename == "gummersbach d"
replace nCH = 1 if gdename == "pekin rc"
replace nCH = 1 if gdename == "rastatt d"
replace nCH = 1 if gdename == "edina mn usa"
replace nCH = 1 if gdename == "viuz en sallaz f"
replace nCH = 1 if gdename == "weil a. rhein d"
replace nCH = 1 if gdename == "sintra p"
replace nCH = 1 if gdename == "boston (usa)"
replace nCH = 1 if gdename == "la souterraine f"
replace nCH = 1 if gdename == "singapore"
replace nCH = 1 if gdename == "limai b"
replace nCH = 1 if gdename == "noisy le roi f"
replace nCH = 1 if gdename == "cairo et"
replace nCH = 1 if gdename == "blotzheim f"
replace nCH = 1 if gdename == "watford gb"
replace nCH = 1 if gdename == "oakville"
replace nCH = 1 if gdename == "kingston gb"
replace nCH = 1 if gdename == "rheinfelden (d)"
replace nCH = 1 if gdename == "bernardsville"
replace nCH = 1 if gdename == "funchal p"
replace nCH = 1 if gdename == "rickmansworth gb"
replace nCH = 1 if gdename == "offenbach am main d"
replace nCH = 1 if gdename == "surrey"
replace nCH = 1 if gdename == "tokio japan"
replace nCH = 1 if gdename == "faloppio i"
replace nCH = 1 if gdename == "hesingue f"
replace nCH = 1 if gdename == "roskilde dk"
replace nCH = 1 if gdename == "flein d"
replace nCH = 1 if gdename == "bucarest ro"
replace nCH = 1 if gdename == "frankfurt a.m. (d)"
replace nCH = 1 if gdename == "rancho santa fe usa"
replace nCH = 1 if gdename == "la jolla"
replace nCH = 1 if gdename == "hambourg d"
replace nCH = 1 if gdename == "bourg la reine f"
replace nCH = 1 if gdename == "barcelone"
replace nCH = 1 if gdename == "freiburg i. br. d"
replace nCH = 1 if gdename == "santiago rch"
replace nCH = 1 if gdename == "damaskus syr"
replace nCH = 1 if gdename == "novokouznetzk rus"
replace nCH = 1 if gdename == "valleiry f"
replace nCH = 1 if gdename == "frederiksberg dk"
replace nCH = 1 if gdename == "cincinnati usa"
replace nCH = 1 if gdename == "seeshaupt d"
replace nCH = 1 if gdename == "metz f"
replace nCH = 1 if gdename == "fairfield usa"
replace nCH = 1 if gdename == "cascais (p)"
replace nCH = 1 if gdename == "porto (p)"
replace nCH = 1 if gdename == "rom/i"
replace nCH = 1 if gdename == "gallarate i"
replace nCH = 1 if gdename == "grasse f"
replace nCH = 1 if gdename == "ostwald f"
replace nCH = 1 if gdename == "marblehead usa"
replace nCH = 1 if gdename == "wangen im allgaeu d"
replace nCH = 1 if gdename == "obanos e"
replace nCH = 1 if gdename == "jerusalem (il)"
replace nCH = 1 if gdename == "saulgau d"
replace nCH = 1 if gdename == "new delhi ind"
replace nCH = 1 if gdename == "klosterneuburg a"
replace nCH = 1 if gdename == "vaernamo s"
replace nCH = 1 if gdename == "selb d"
replace nCH = 1 if gdename == "malland (i)"
replace nCH = 1 if gdename == "hard (a)"
replace nCH = 1 if gdename == "moerfelden walldorf d"
replace nCH = 1 if gdename == "chavenay f"
replace nCH = 1 if gdename == "gibraltar gb"
replace nCH = 1 if gdename == "lecco i"
replace nCH = 1 if gdename == "wiltshire gb"
replace nCH = 1 if gdename == "salt lake city usa"
replace nCH = 1 if gdename == "eichenau d"
replace nCH = 1 if gdename == "wilmington usa"
replace nCH = 1 if gdename == "madrid/e"
replace nCH = 1 if gdename == "washington (usa)"
replace nCH = 1 if gdename == "bretten d"
replace nCH = 1 if gdename == "bologna (i)"
replace nCH = 1 if gdename == "rickenbach d"
replace nCH = 1 if gdename == "novate milanese i"
replace nCH = 1 if gdename == "neuhausen a.r"
replace nCH = 1 if gdename == "venlo nl"
replace nCH = 1 if gdename == "brezzo di bedero i"
replace nCH = 1 if gdename == "tunis tn"
replace nCH = 1 if gdename == "kleve d"
replace nCH = 1 if gdename == "thonon f"
replace nCH = 1 if gdename == "oxon gb"
replace nCH = 1 if gdename == "saltsjoebaden s"
replace nCH = 1 if gdename == "dogern d"
replace nCH = 1 if gdename == "huy b"
replace nCH = 1 if gdename == "moscou russie"
replace nCH = 1 if gdename == "hauts de seine f"
replace nCH = 1 if gdename == "aidlingen d"
replace nCH = 1 if gdename == "wattens (a)"
replace nCH = 1 if gdename == "tokio (japan)"
replace nCH = 1 if gdename == "hamilton bm"
replace nCH = 1 if gdename == "bjerringbro dk"
replace nCH = 1 if gdename == "saverne f"
replace nCH = 1 if gdename == "schutterwald d"
replace nCH = 1 if gdename == "koenigstein/taunus d"
replace nCH = 1 if gdename == "saint nom la breteche f"
replace nCH = 1 if gdename == "berg d"
replace nCH = 1 if gdename == "dachau d"
replace nCH = 1 if gdename == "bergisch gladbach (d)"
replace nCH = 1 if gdename == "strasbourg (f)"
replace nCH = 1 if gdename == "wilhelmshaven d"
replace nCH = 1 if gdename == "norderstedt d"
replace nCH = 1 if gdename == "taufkirchen d"
replace nCH = 1 if gdename == "epernay f"
replace nCH = 1 if gdename == "caire et"
replace nCH = 1 if gdename == "tilburg nl"
replace nCH = 1 if gdename == "st. poelten a"
replace nCH = 1 if gdename == "landeck a"
replace nCH = 1 if gdename == "peterborough gb"
replace nCH = 1 if gdename == "delle (f)"
replace nCH = 1 if gdename == "kentucky usa"
replace nCH = 1 if gdename == "russland"
replace nCH = 1 if gdename == "arona i"
replace nCH = 1 if gdename == "fieberbrunn a"
replace nCH = 1 if gdename == "weiterstadt d"
replace nCH = 1 if gdename == "holzminden d"
replace nCH = 1 if gdename == "washington d.c. usa"
replace nCH = 1 if gdename == "richmond usa"
replace nCH = 1 if gdename == "cunardo i"
replace nCH = 1 if gdename == "herrsching d"
replace nCH = 1 if gdename == "frastanz a"
replace nCH = 1 if gdename == "hohenems a"
replace nCH = 1 if gdename == "peron f"
replace nCH = 1 if gdename == "abidjan cl"
replace nCH = 1 if gdename == "hameln d"
replace nCH = 1 if gdename == "weiden d"
replace nCH = 1 if gdename == "boulogne (f)"
replace nCH = 1 if gdename == "markt einersheim d"
replace nCH = 1 if gdename == "croydon gb"
replace nCH = 1 if gdename == "hemer d"
replace nCH = 1 if gdename == "monte carlo monaco"
replace nCH = 1 if gdename == "chiba j"
replace nCH = 1 if gdename == "la celle saint cloud f"
replace nCH = 1 if gdename == "groningen nl"
replace nCH = 1 if gdename == "peschiera borromeo i"
replace nCH = 1 if gdename == "brunn am gebirge a"
replace nCH = 1 if gdename == "marchirolo i"
replace nCH = 1 if gdename == "wavre/thielle wavre"
replace nCH = 1 if gdename == "djursholm (s)"
replace nCH = 1 if gdename == "solingen (d)"
replace nCH = 1 if gdename == "santo domingo dom"
replace nCH = 1 if gdename == "montreal (cnd)"
replace nCH = 1 if gdename == "korntal d"
replace nCH = 1 if gdename == "newport beach usa"
replace nCH = 1 if gdename == "choisy f"
replace nCH = 1 if gdename == "newbury gb"
replace nCH = 1 if gdename == "rhode saint genese b"
replace nCH = 1 if gdename == "horsholm dk"
replace nCH = 1 if gdename == "stocksund s"
replace nCH = 1 if gdename == "reno"
replace nCH = 1 if gdename == "prag cs"
replace nCH = 1 if gdename == "caracas (venezuela)"
replace nCH = 1 if gdename == "riyadh (sa)"
replace nCH = 1 if gdename == "willich d"
replace nCH = 1 if gdename == "waldbronn d"
replace nCH = 1 if gdename == "huningue f"
replace nCH = 1 if gdename == "st genis pouilly f"
replace nCH = 1 if gdename == "villers le lac f"
replace nCH = 1 if gdename == "syosset usa"
replace nCH = 1 if gdename == "djakarta ri"
replace nCH = 1 if gdename == "wangen d"
replace nCH = 1 if gdename == "bursa tr"
replace nCH = 1 if gdename == "bonneville f"
replace nCH = 1 if gdename == "herdecke d"
replace nCH = 1 if gdename == "bassano del grappa i"
replace nCH = 1 if gdename == "memmingen d"
replace nCH = 1 if gdename == "hohentengen d"
replace nCH = 1 if gdename == "luxenburg l"
replace nCH = 1 if gdename == "anaheim usa"
replace nCH = 1 if gdename == "wemmel b"
replace nCH = 1 if gdename == "eurasburg d"
replace nCH = 1 if gdename == "garbsen d"
replace nCH = 1 if gdename == "rouen f"
replace nCH = 1 if gdename == "ahlen d"
replace nCH = 1 if gdename == "leipzig d"
replace nCH = 1 if gdename == "rinteln d"
replace nCH = 1 if gdename == "le chesnay (f)"
replace nCH = 1 if gdename == "soresina i"
replace nCH = 1 if gdename == "duisburg (d)"
replace nCH = 1 if gdename == "billund (dk)"
replace nCH = 1 if gdename == "eltville d"
replace nCH = 1 if gdename == "vers f"
replace nCH = 1 if gdename == "merate i"
replace nCH = 1 if gdename == "tokyo (japan)"
replace nCH = 1 if gdename == "chatenay malabry f"
replace nCH = 1 if gdename == "moskau ru"
replace nCH = 1 if gdename == "son nl"
replace nCH = 1 if gdename == "trelleborg s"
replace nCH = 1 if gdename == "etterbeek b"
replace nCH = 1 if gdename == "sunningdale gb"
replace nCH = 1 if gdename == "east horsley"
replace nCH = 1 if gdename == "bad salzuflen d"
replace nCH = 1 if gdename == "florida (usa)"
replace nCH = 1 if gdename == "rheinberg d"
replace nCH = 1 if gdename == "chagrin falls usa"
replace nCH = 1 if gdename == "jersey"
replace nCH = 1 if gdename == "burlington"
replace nCH = 1 if gdename == "chilly mazarin f"
replace nCH = 1 if gdename == "capiago intimiano i"
replace nCH = 1 if gdename == "melsungen d"
replace nCH = 1 if gdename == "lodi i"
replace nCH = 1 if gdename == "heiligenhaus d"
replace nCH = 1 if gdename == "florenz"
replace nCH = 1 if gdename == "amberg d"
replace nCH = 1 if gdename == "waalwijk nl"
replace nCH = 1 if gdename == "rockville usa"
replace nCH = 1 if gdename == "troyes f"
replace nCH = 1 if gdename == "holland"
replace nCH = 1 if gdename == "coye la foret f"
replace nCH = 1 if gdename == "senden d"
replace nCH = 1 if gdename == "muelheim an der ruhr d"
replace nCH = 1 if gdename == "(usa)"
replace nCH = 1 if gdename == "frankfurt/m. d"
replace nCH = 1 if gdename == "messery f"
replace nCH = 1 if gdename == "feldafing d"
replace nCH = 1 if gdename == "mississauga cdn"
replace nCH = 1 if gdename == "grilly f"
replace nCH = 1 if gdename == "charlottenlund (dk)"
replace nCH = 1 if gdename == "vaterstetten d"
replace nCH = 1 if gdename == "walton on thames"
replace nCH = 1 if gdename == "lindingoe s"
replace nCH = 1 if gdename == "twickenham gb"
replace nCH = 1 if gdename == "voorschoten nl"
replace nCH = 1 if gdename == "noisy le grand f"
replace nCH = 1 if gdename == "jar n"
replace nCH = 1 if gdename == "emmerich d"
replace nCH = 1 if gdename == "sindelfingen (d)"
replace nCH = 1 if gdename == "waldenbuch d"
replace nCH = 1 if gdename == "dueren (d)"
replace nCH = 1 if gdename == "akarp s"
replace nCH = 1 if gdename == "velp nl"
replace nCH = 1 if gdename == "kissing d"
replace nCH = 1 if gdename == "aichtal d"
replace nCH = 1 if gdename == "andilly f"
replace nCH = 1 if gdename == "weinheim (d)"
replace nCH = 1 if gdename == "tortola bvi"
replace nCH = 1 if gdename == "colmar f"
replace nCH = 1 if gdename == "terni i"
replace nCH = 1 if gdename == "pavia i"
replace nCH = 1 if gdename == "ontario can"
replace nCH = 1 if gdename == "jakarta ri"
replace nCH = 1 if gdename == "woluwe st pierre b"
replace nCH = 1 if gdename == "wlen a"
replace nCH = 1 if gdename == "stolberg d"
replace nCH = 1 if gdename == "riyadh ksa"
replace nCH = 1 if gdename == "durach d"
replace nCH = 1 if gdename == "accra gh"
replace nCH = 1 if gdename == "darien (usa)"
replace nCH = 1 if gdename == "malta"
replace nCH = 1 if gdename == "bietigheim d"
replace nCH = 1 if gdename == "vaucresson (f)"
replace nCH = 1 if gdename == "oehningen d"
replace nCH = 1 if gdename == "bartenheim f"
replace nCH = 1 if gdename == "toronto (kanada)"
replace nCH = 1 if gdename == "montegranaro i"
replace nCH = 1 if gdename == "caracas w"
replace nCH = 1 if gdename == "alberta cdn"
replace nCH = 1 if gdename == "koetz d"
replace nCH = 1 if gdename == "sprockhoevel d"
replace nCH = 1 if gdename == "gouvieux f"
replace nCH = 1 if gdename == "hoechst (a)"
replace nCH = 1 if gdename == "lindau (d)"
replace nCH = 1 if gdename == "overveen nl"
replace nCH = 1 if gdename == "harrogate gb"
replace nCH = 1 if gdename == "ljungby s"
replace nCH = 1 if gdename == "ingatestone gb"
replace nCH = 1 if gdename == "epinal (f)"
replace nCH = 1 if gdename == "3063ittigen"
replace nCH = 1 if gdename == "epping gb"
replace nCH = 1 if gdename == "(gb)"
replace nCH = 1 if gdename == "cleveland usa"
replace nCH = 1 if gdename == "wien/a"
replace nCH = 1 if gdename == "oberursel (d)"
replace nCH = 1 if gdename == "prato i"
replace nCH = 1 if gdename == "untere"
replace nCH = 1 if gdename == "mailand 1"
replace nCH = 1 if gdename == "ostuni i"
replace nCH = 1 if gdename == "maroc"
replace nCH = 1 if gdename == "espoo (sf)"
replace nCH = 1 if gdename == "leatherhead gb"
replace nCH = 1 if gdename == "xanten d"
replace nCH = 1 if gdename == "la roche sur foron f"
replace nCH = 1 if gdename == "curacao"
replace nCH = 1 if gdename == "winston salem usa"
replace nCH = 1 if gdename == "vallentuna s"
replace nCH = 1 if gdename == "culoz f"
replace nCH = 1 if gdename == "jersey (gb)"
replace nCH = 1 if gdename == "paguera"
replace nCH = 1 if gdename == "grosshansdorf d"
replace nCH = 1 if gdename == "cheltenham gb"
replace nCH = 1 if gdename == "seligenstadt d"
replace nCH = 1 if gdename == "gelsenkirchen d"
replace nCH = 1 if gdename == "lahr d"
replace nCH = 1 if gdename == "dordrecht nl"
replace nCH = 1 if gdename == "eastbourne gb"
replace nCH = 1 if gdename == "koenigsberg d"
replace nCH = 1 if gdename == "kingston upon thames gb"
replace nCH = 1 if gdename == "schwanengasse 8"
replace nCH = 1 if gdename == "floda s"
replace nCH = 1 if gdename == "wetzlar d"
replace nCH = 1 if gdename == "la havane c"
replace nCH = 1 if gdename == "montgomery"
replace nCH = 1 if gdename == "wachtendonk d"
replace nCH = 1 if gdename == "marbella (e)"
replace nCH = 1 if gdename == "casablanca (ma)"
replace nCH = 1 if gdename == "reading gb"
replace nCH = 1 if gdename == "schweinfurt d"
replace nCH = 1 if gdename == "bad lippspringe d"
replace nCH = 1 if gdename == "memphis (usa)"
replace nCH = 1 if gdename == "waiblingen d"
replace nCH = 1 if gdename == "klettgau d"
replace nCH = 1 if gdename == "pescara i"
replace nCH = 1 if gdename == "barcelone (e)"
replace nCH = 1 if gdename == "tempe usa"
replace nCH = 1 if gdename == "st. louis f"
replace nCH = 1 if gdename == "algier dz"
replace nCH = 1 if gdename == "damaskus (syr)"
replace nCH = 1 if gdename == "castricum nl"
replace nCH = 1 if gdename == "leverkusen"
replace nCH = 1 if gdename == "perchtoldsdorf a"
replace nCH = 1 if gdename == "stanwick gb"
replace nCH = 1 if gdename == "leitershofen d"
replace nCH = 1 if gdename == "suessen d"
replace nCH = 1 if gdename == "bornheim d"
replace nCH = 1 if gdename == "herzogenaurach d"
replace nCH = 1 if gdename == "cattolica i"
replace nCH = 1 if gdename == "oisterwijk nl"
replace nCH = 1 if gdename == "bruxelles/b"
replace nCH = 1 if gdename == "essen dl"
replace nCH = 1 if gdename == "stuttgart/d"
replace nCH = 1 if gdename == "gif sur yvette f"
replace nCH = 1 if gdename == "chambery f"
replace nCH = 1 if gdename == "luebeck (d)"
replace nCH = 1 if gdename == "rennes f"
replace nCH = 1 if gdename == "helmond nl"
replace nCH = 1 if gdename == "ukraine"
replace nCH = 1 if gdename == "kuweit kwt"
replace nCH = 1 if gdename == "lippstadt d"
replace nCH = 1 if gdename == "voorthuizen nl"
replace nCH = 1 if gdename == "montevideo u"
replace nCH = 1 if gdename == "villeurbanne f"
replace nCH = 1 if gdename == "waldshut brd"
replace nCH = 1 if gdename == "herbrechtingen d"
replace nCH = 1 if gdename == "dresden d"
replace nCH = 1 if gdename == "los altos hills"
replace nCH = 1 if gdename == "maintal (d)"
replace nCH = 1 if gdename == "cdn"

replace nCH = 1 if gdename == "houston, texas usa"
replace nCH = 1 if gdename == "greenwich, connecticut usa"
replace nCH = 1 if gdename == "dallas, texas usa"
replace nCH = 1 if gdename == "caracas, venezuela"
replace nCH = 1 if gdename == "chicago, illinois usa"
replace nCH = 1 if gdename == "esher, surrey gb"
replace nCH = 1 if gdename == "douglas, isle of man"
replace nCH = 1 if gdename == "cobham, surrey gb"
replace nCH = 1 if gdename == "princeton, new jersey usa"
replace nCH = 1 if gdename == "miami, florida usa"
replace nCH = 1 if gdename == "northwood, middlesex gb"
replace nCH = 1 if gdename == "greenwich, conn. usa"
replace nCH = 1 if gdename == "new york, usa"
replace nCH = 1 if gdename == "richmond, virginia usa"
replace nCH = 1 if gdename == "fort worth, texas usa"
replace nCH = 1 if gdename == "atlanta, georgia usa"
replace nCH = 1 if gdename == "pittsburgh, pennsylvania usa"
replace nCH = 1 if gdename == "darien, connecticut usa"
replace nCH = 1 if gdename == "palo alto, california usa"
replace nCH = 1 if gdename == "sevenoaks, kent gb"
replace nCH = 1 if gdename == "woking, surrey gb"
replace nCH = 1 if gdename == "douala, kamerun"
replace nCH = 1 if gdename == "york, pennsylvania usa"
replace nCH = 1 if gdename == "midland, michigan usa"
replace nCH = 1 if gdename == "tucson, arizona usa"
replace nCH = 1 if gdename == "la jolla, kalifornien usa"
replace nCH = 1 if gdename == "racine, wisconsin usa"
replace nCH = 1 if gdename == "new canaan, connecticut usa"
replace nCH = 1 if gdename == "andover, hampshire gb"
replace nCH = 1 if gdename == "short hills, new jersey usa"
replace nCH = 1 if gdename == "los angeles, kalifornien usa"
replace nCH = 1 if gdename == "cincinnati, ohio usa"
replace nCH = 1 if gdename == "hamilton, bermuda"
replace nCH = 1 if gdename == "oakville, ontario cdn"
replace nCH = 1 if gdename == "rochester, new york usa"
replace nCH = 1 if gdename == "guildford, surrey gb"
replace nCH = 1 if gdename == "rabat, maroc"
replace nCH = 1 if gdename == "villars"
replace nCH = 1 if gdename == "brasschaat (b)"
replace nCH = 1 if gdename == "mydrecht nl"
replace nCH = 1 if gdename == "boras s"
replace nCH = 1 if gdename == "weil der stadt d"
replace nCH = 1 if gdename == "louisville usa"
replace nCH = 1 if gdename == "gauting (d)"
replace nCH = 1 if gdename == "unterschleissheim d"
replace nCH = 1 if gdename == "loehne d"
replace nCH = 1 if gdename == "verl d"
replace nCH = 1 if gdename == "maastricht nl"
replace nCH = 1 if gdename == "lumezzane i"
replace nCH = 1 if gdename == "trosa s"
replace nCH = 1 if gdename == "radolfzell (d)"
replace nCH = 1 if gdename == "eindhoven (nl)"
replace nCH = 1 if gdename == "dubai vae"
replace nCH = 1 if gdename == "sondrio i"
replace nCH = 1 if gdename == "summit usa"
replace nCH = 1 if gdename == "bad zwischenahn d"
replace nCH = 1 if gdename == "khartoum sud"
replace nCH = 1 if gdename == "tripolis (lt)"
replace nCH = 1 if gdename == "pittsburgh (usa)"
replace nCH = 1 if gdename == "burlington, ma usa"
replace nCH = 1 if gdename == "andernach d"
replace nCH = 1 if gdename == "geretsried d"
replace nCH = 1 if gdename == "cannobbio i"
replace nCH = 1 if gdename == "busseto i"
replace nCH = 1 if gdename == "calgary (cdn)"
replace nCH = 1 if gdename == "neapel (i)"
replace nCH = 1 if gdename == "leoben a"
replace nCH = 1 if gdename == "vallendar d"
replace nCH = 1 if gdename == "oldenburg (d)"
replace nCH = 1 if gdename == "bayreuth d"
replace nCH = 1 if gdename == "nassau (bahamas)"
replace nCH = 1 if gdename == "waldenburg d"
replace nCH = 1 if gdename == "francfort sur le main d"
replace nCH = 1 if gdename == "lich d"
replace nCH = 1 if gdename == "iggesund s"
replace nCH = 1 if gdename == "hillerod dk"
replace nCH = 1 if gdename == "aix en provence (f)"
replace nCH = 1 if gdename == "saint ulrich f"
replace nCH = 1 if gdename == "tegernsee d"
replace nCH = 1 if gdename == "triest (i)"
replace nCH = 1 if gdename == "leibnitz a"
replace nCH = 1 if gdename == "memphis usa"
replace nCH = 1 if gdename == "hilchenbach d"
replace nCH = 1 if gdename == "wangen (d)"
replace nCH = 1 if gdename == "waldstetten d"
replace nCH = 1 if gdename == "viroflay f"
replace nCH = 1 if gdename == "pasadena usa"
replace nCH = 1 if gdename == "new barnet gb"
replace nCH = 1 if gdename == "wimbledon/london gb"
replace nCH = 1 if gdename == "vienne a"
replace nCH = 1 if gdename == "fourqueux f"
replace nCH = 1 if gdename == "instanbul tr"
replace nCH = 1 if gdename == "altenkirchen d"
replace nCH = 1 if gdename == "muehlem d"
replace nCH = 1 if gdename == "plano, texas usa"
replace nCH = 1 if gdename == "karachi pk"
replace nCH = 1 if gdename == "beaconsfield (gb)"
replace nCH = 1 if gdename == "prague cz"
replace nCH = 1 if gdename == "boston, mass. usa"
replace nCH = 1 if gdename == "hanau (d)"
replace nCH = 1 if gdename == "fontenay aux roses f"
replace nCH = 1 if gdename == "pamplona e"
replace nCH = 1 if gdename == "hillsborough usa"
replace nCH = 1 if gdename == "marcq en baroeul f"
replace nCH = 1 if gdename == "naarden (nl)"
replace nCH = 1 if gdename == "calpe e"
replace nCH = 1 if gdename == "uithoorn nl"
replace nCH = 1 if gdename == "oberhausen d"
replace nCH = 1 if gdename == "ybbsitz a"
replace nCH = 1 if gdename == "hirschberg d"
replace nCH = 1 if gdename == "loisin f"
replace nCH = 1 if gdename == "neunkirchen d"
replace nCH = 1 if gdename == "wels a"
replace nCH = 1 if gdename == "chevreuse f"
replace nCH = 1 if gdename == "oviedo e"
replace nCH = 1 if gdename == "prevessin f"
replace nCH = 1 if gdename == "esslingen am neckar d"
replace nCH = 1 if gdename == "riyad (sa)"
replace nCH = 1 if gdename == "neapel"
replace nCH = 1 if gdename == "kinshasa zre"
replace nCH = 1 if gdename == "guernsey (channel islands)"
replace nCH = 1 if gdename == "massongy f"
replace nCH = 1 if gdename == "wain (d)"
replace nCH = 1 if gdename == "paguera, mallorca e"
replace nCH = 1 if gdename == "waalre nl"
replace nCH = 1 if gdename == "neydens f"
replace nCH = 1 if gdename == "klaus a"
replace nCH = 1 if gdename == "portland oregon usa"
replace nCH = 1 if gdename == "endhoven nl"
replace nCH = 1 if gdename == "riga lettland"
replace nCH = 1 if gdename == "conegliano i"
replace nCH = 1 if gdename == "palafrugell e"
replace nCH = 1 if gdename == "denzlingen d"
replace nCH = 1 if gdename == "sarrebourg f"
replace nCH = 1 if gdename == "maulberg d"
replace nCH = 1 if gdename == "hilleroed dk"
replace nCH = 1 if gdename == "krommenie nl"
replace nCH = 1 if gdename == "camberley gb"
replace nCH = 1 if gdename == "garmisch partenkirchen d"
replace nCH = 1 if gdename == "monte carlo, mc"
replace nCH = 1 if gdename == "leutkirch d"
replace nCH = 1 if gdename == "accra (gh)"
replace nCH = 1 if gdename == "la garenne colombes f"
replace nCH = 1 if gdename == "arlon b"
replace nCH = 1 if gdename == "maarssen nl"
replace nCH = 1 if gdename == "aichwald (d)"
replace nCH = 1 if gdename == "nairobi kenya"
replace nCH = 1 if gdename == "chislehurst gb"
replace nCH = 1 if gdename == "badenweiler d"
replace nCH = 1 if gdename == "izmir tr"
replace nCH = 1 if gdename == "hasselt b"
replace nCH = 1 if gdename == "dudelange l"
replace nCH = 1 if gdename == "deurne nl"
replace nCH = 1 if gdename == "moskau (rus)"
replace nCH = 1 if gdename == "wangen zh"
replace nCH = 1 if gdename == "fulda d"
replace nCH = 1 if gdename == "weiler a"
replace nCH = 1 if gdename == "muehlacker d"
replace nCH = 1 if gdename == "aeugst a.a"
replace nCH = 1 if gdename == "cecina i"
replace nCH = 1 if gdename == "gaiberg d"
replace nCH = 1 if gdename == "frankfurt a/m d"
replace nCH = 1 if gdename == "abbevillers f"
replace nCH = 1 if gdename == "bad bellingen d"
replace nCH = 1 if gdename == "warstein d"
replace nCH = 1 if gdename == "huizen (nl)"
replace nCH = 1 if gdename == "fulmer gb"
replace nCH = 1 if gdename == "durban za"
replace nCH = 1 if gdename == "euskirchen d"
replace nCH = 1 if gdename == "aosta i"
replace nCH = 1 if gdename == "goldbach d"
replace nCH = 1 if gdename == "st julien en genevois (f)"
replace nCH = 1 if gdename == "hattingen d"
replace nCH = 1 if gdename == "marbella/malaga e"
replace nCH = 1 if gdename == "high wycombe gb"
replace nCH = 1 if gdename == "tour de scay f"
replace nCH = 1 if gdename == "toulouse f"
replace nCH = 1 if gdename == "scarsdale usa"
replace nCH = 1 if gdename == "bernardsville, new jersey usa"
replace nCH = 1 if gdename == "haan (d)"
replace nCH = 1 if gdename == "wolfenbuettel d"
replace nCH = 1 if gdename == "lake forest, illinois usa"
replace nCH = 1 if gdename == "kauniainen fin"
replace nCH = 1 if gdename == "newcastle upon tyne gb"
replace nCH = 1 if gdename == "rielasingen d"
replace nCH = 1 if gdename == "sidney aus"
replace nCH = 1 if gdename == "kaiserslautern d"
replace nCH = 1 if gdename == "dillenburg d"
replace nCH = 1 if gdename == "charbonnieres les bains f"
replace nCH = 1 if gdename == "mareil marly f"
replace nCH = 1 if gdename == "schwendi d"
replace nCH = 1 if gdename == "weil d"
replace nCH = 1 if gdename == "pinner gb"
replace nCH = 1 if gdename == "padborg dk"
replace nCH = 1 if gdename == "pallanza i"
replace nCH = 1 if gdename == "friedrichshafen (d)"
replace nCH = 1 if gdename == "tirano (i)"
replace nCH = 1 if gdename == "limburgerhof d"
replace nCH = 1 if gdename == "seattle usa"
replace nCH = 1 if gdename == "ede nl"
replace nCH = 1 if gdename == "ochsenhausen d"
replace nCH = 1 if gdename == "sark (channel islands) gb"
replace nCH = 1 if gdename == "busto garolfo i"
replace nCH = 1 if gdename == "harken drive (cdn)"
replace nCH = 1 if gdename == "gaillard (f)"
replace nCH = 1 if gdename == "moscou r"
replace nCH = 1 if gdename == "roermond nl"
replace nCH = 1 if gdename == "detroit usa"
replace nCH = 1 if gdename == "neubiberg d"
replace nCH = 1 if gdename == "pringy f"
replace nCH = 1 if gdename == "chaumont f"
replace nCH = 1 if gdename == "hallein a"
replace nCH = 1 if gdename == "bagdad irq"
replace nCH = 1 if gdename == "sandviken s"
replace nCH = 1 if gdename == "winnenden d"
replace nCH = 1 if gdename == "ashtead gb"
replace nCH = 1 if gdename == "formigine i"
replace nCH = 1 if gdename == "neuilly s/seine f"
replace nCH = 1 if gdename == "donzdorf d"
replace nCH = 1 if gdename == "indianapolis usa"
replace nCH = 1 if gdename == "glendale usa"
replace nCH = 1 if gdename == "knokke b"
replace nCH = 1 if gdename == "nassau/bahamas"
replace nCH = 1 if gdename == "fort lauderdale usa"
replace nCH = 1 if gdename == "dornstetten d"
replace nCH = 1 if gdename == "charleroi b"
replace nCH = 1 if gdename == "delden nl"
replace nCH = 1 if gdename == "pirmasens d"
replace nCH = 1 if gdename == "frohnleiten a"
replace nCH = 1 if gdename == "biarritz f"
replace nCH = 1 if gdename == "corcelles/c. cormondreche"
replace nCH = 1 if gdename == "hemel hempstead, herts gb"
replace nCH = 1 if gdename == "cessy f"
replace nCH = 1 if gdename == "saint genis pouilly f"
replace nCH = 1 if gdename == "schwanau d"
replace nCH = 1 if gdename == "kaufering d"
replace nCH = 1 if gdename == "montecarlo"
replace nCH = 1 if gdename == "herrsching am ammersee d"
replace nCH = 1 if gdename == "rimini i"
replace nCH = 1 if gdename == "haiti rh"
replace nCH = 1 if gdename == "butzbach (d)"
replace nCH = 1 if gdename == "merdingen d"
replace nCH = 1 if gdename == "mestre i"
replace nCH = 1 if gdename == "kamp lintfort d"
replace nCH = 1 if gdename == "turnhout b"
replace nCH = 1 if gdename == "laufen d"
replace nCH = 1 if gdename == "solna s"
replace nCH = 1 if gdename == "redhill, surrey gb"
replace nCH = 1 if gdename == "offenburg (d)"
replace nCH = 1 if gdename == "albinea i"
replace nCH = 1 if gdename == "stadtallendorf d"
replace nCH = 1 if gdename == "kirchzarten d"
replace nCH = 1 if gdename == "brasschat b"
replace nCH = 1 if gdename == "fuerth (d)"
replace nCH = 1 if gdename == "lamorlaye f"
replace nCH = 1 if gdename == "san rafael usa"
replace nCH = 1 if gdename == "michelbach le bas f"
replace nCH = 1 if gdename == "alassio i"
replace nCH = 1 if gdename == "allmendingen d"
replace nCH = 1 if gdename == "hunt valley usa"
replace nCH = 1 if gdename == "lake forest (usa)"
replace nCH = 1 if gdename == "schwaebisch hall d"
replace nCH = 1 if gdename == "toronto (cnd)"
replace nCH = 1 if gdename == "walton on thames, surrey gb"
replace nCH = 1 if gdename == "koweït kwt"
replace nCH = 1 if gdename == "hove b"
replace nCH = 1 if gdename == "desio i"
replace nCH = 1 if gdename == "altach a"
replace nCH = 1 if gdename == "reigate surrey gb"
replace nCH = 1 if gdename == "aachen (d)"
replace nCH = 1 if gdename == "roedermark d"
replace nCH = 1 if gdename == "st jean d aulps f"
replace nCH = 1 if gdename == "turin 1"
replace nCH = 1 if gdename == "zaandam nl"
replace nCH = 1 if gdename == "belgien"
replace nCH = 1 if gdename == "kirchheim d"
replace nCH = 1 if gdename == "camporosso i"
replace nCH = 1 if gdename == "gorle i"
replace nCH = 1 if gdename == "stanmore gb"
replace nCH = 1 if gdename == "kontich b"
replace nCH = 1 if gdename == "willich (d)"
replace nCH = 1 if gdename == "monaco mco"
replace nCH = 1 if gdename == "richmond surrey gb"
replace nCH = 1 if gdename == "farnham gb"
replace nCH = 1 if gdename == "cambridge usa"
replace nCH = 1 if gdename == "new rochelle usa"
replace nCH = 1 if gdename == "bethesda, maryland usa"
replace nCH = 1 if gdename == "patalavaca/gran canaria"
replace nCH = 1 if gdename == "hoellviken s"
replace nCH = 1 if gdename == "alella e"
replace nCH = 1 if gdename == "saint louis f"
replace nCH = 1 if gdename == "bogota (co)"
replace nCH = 1 if gdename == "ditzingen d"
replace nCH = 1 if gdename == "neutraubling d"
replace nCH = 1 if gdename == "ravensburg (d)"
replace nCH = 1 if gdename == "tours f"
replace nCH = 1 if gdename == "wolfach (d)"
replace nCH = 1 if gdename == "hinsdale, illinois usa"
replace nCH = 1 if gdename == "ascot (gb)"
replace nCH = 1 if gdename == "dornhan d"
replace nCH = 1 if gdename == "nernier f"
replace nCH = 1 if gdename == "kristianstad s"
replace nCH = 1 if gdename == "vallauris (f)"
replace nCH = 1 if gdename == "belfort f"
replace nCH = 1 if gdename == "roettenbach d"
replace nCH = 1 if gdename == "wehrheim d"
replace nCH = 1 if gdename == "padova (i)"
replace nCH = 1 if gdename == "rantzwiller f"
replace nCH = 1 if gdename == "maidenhead (gb)"
replace nCH = 1 if gdename == "abu dhabi uae"
replace nCH = 1 if gdename == "rom 1"
replace nCH = 1 if gdename == "sandhausen d"
replace nCH = 1 if gdename == "stockport gb"
replace nCH = 1 if gdename == "valkenswaard nl"
replace nCH = 1 if gdename == "bamberg d"
replace nCH = 1 if gdename == "wettolsheim f"
replace nCH = 1 if gdename == "alfdorf d"
replace nCH = 1 if gdename == "imola i"
replace nCH = 1 if gdename == "bad schoenborn d"
replace nCH = 1 if gdename == "draguignan f"
replace nCH = 1 if gdename == "vught (nl)"
replace nCH = 1 if gdename == "miami florida usa"
replace nCH = 1 if gdename == "zell im wiesental d"
replace nCH = 1 if gdename == "silkeborg dk"
replace nCH = 1 if gdename == "lima (pe)"
replace nCH = 1 if gdename == "larchmont usa"
replace nCH = 1 if gdename == "vogelenzang nl"
replace nCH = 1 if gdename == "legnano i"
replace nCH = 1 if gdename == "reno, nevada usa"
replace nCH = 1 if gdename == "curacao nl"
replace nCH = 1 if gdename == "zagreb (yu)"
replace nCH = 1 if gdename == "norwalk (usa)"
replace nCH = 1 if gdename == "bangui rca"
replace nCH = 1 if gdename == "bedburg d"
replace nCH = 1 if gdename == "dubai (et)"
replace nCH = 1 if gdename == "colchester gb"
replace nCH = 1 if gdename == "ridgewood usa"
replace nCH = 1 if gdename == "aalen (d)"
replace nCH = 1 if gdename == "east horsley (gb)"
replace nCH = 1 if gdename == "offenbach"
replace nCH = 1 if gdename == "maria enzersdorf a"
replace nCH = 1 if gdename == "la varenne saint hilaire f"
replace nCH = 1 if gdename == "unternach d"
replace nCH = 1 if gdename == "obernburg d"
replace nCH = 1 if gdename == "buedelsdorf d"
replace nCH = 1 if gdename == "minden (d)"
replace nCH = 1 if gdename == "lorrez le bocage preaux f"
replace nCH = 1 if gdename == "sherman oaks usa"
replace nCH = 1 if gdename == "st dizier f"
replace nCH = 1 if gdename == "trivero i"
replace nCH = 1 if gdename == "rehau d"
replace nCH = 1 if gdename == "bristol (gb)"
replace nCH = 1 if gdename == "muenchberg d"
replace nCH = 1 if gdename == "staufen d"
replace nCH = 1 if gdename == "schwieberdingen d"
replace nCH = 1 if gdename == "conshohocken usa"
replace nCH = 1 if gdename == "st. louis usa"
replace nCH = 1 if gdename == "gee cross gb"
replace nCH = 1 if gdename == "saltsjoebaden (s)"
replace nCH = 1 if gdename == "deggendorf d"
replace nCH = 1 if gdename == "tremblay f"
replace nCH = 1 if gdename == "stamford (usa)"
replace nCH = 1 if gdename == "mombasa eak"
replace nCH = 1 if gdename == "diedorf d"
replace nCH = 1 if gdename == "norwalk usa"
replace nCH = 1 if gdename == "raffles lane usa"
replace nCH = 1 if gdename == "saratoga, california usa"
replace nCH = 1 if gdename == "duingt f"
replace nCH = 1 if gdename == "bad kreuznach (d)"
replace nCH = 1 if gdename == "allonzier la caille f"
replace nCH = 1 if gdename == "skopje yu"
replace nCH = 1 if gdename == "duesseldorf dl"
replace nCH = 1 if gdename == "greenwich ct usa"
replace nCH = 1 if gdename == "wickede d"
replace nCH = 1 if gdename == "menlo park usa"
replace nCH = 1 if gdename == "wichita, kansas usa"
replace nCH = 1 if gdename == "althuette d"
replace nCH = 1 if gdename == "siegen d"
replace nCH = 1 if gdename == "seynod f"
replace nCH = 1 if gdename == "ann arbor, michigan usa"
replace nCH = 1 if gdename == "lyngby dk"
replace nCH = 1 if gdename == "cabrera del mar e"
replace nCH = 1 if gdename == "havant gb"
replace nCH = 1 if gdename == "charlotte nc usa"
replace nCH = 1 if gdename == "buchholz d"
replace nCH = 1 if gdename == "offenbach (d)"
replace nCH = 1 if gdename == "neuss (d)"
replace nCH = 1 if gdename == "marnaz f"
replace nCH = 1 if gdename == "coburg d"
replace nCH = 1 if gdename == "maule f"
replace nCH = 1 if gdename == "vanves f"
replace nCH = 1 if gdename == "frankfurt a. main (d)"
replace nCH = 1 if gdename == "oeiras p"
replace nCH = 1 if gdename == "genemuiden nl"
replace nCH = 1 if gdename == "betzdorf d"
replace nCH = 1 if gdename == "aumuehle d"
replace nCH = 1 if gdename == "reigate gb"
replace nCH = 1 if gdename == "mailand (1)"
replace nCH = 1 if gdename == "schopp d"
replace nCH = 1 if gdename == "chesnay f"
replace nCH = 1 if gdename == "estoril p"
replace nCH = 1 if gdename == "bereldange l"
replace nCH = 1 if gdename == "yokohama japan"
replace nCH = 1 if gdename == "marner l"
replace nCH = 1 if gdename == "putten nl"
replace nCH = 1 if gdename == "austin, texas usa"
replace nCH = 1 if gdename == "gibraltar"
replace nCH = 1 if gdename == "bammental d"
replace nCH = 1 if gdename == "epsom gb"
replace nCH = 1 if gdename == "villeneuve loubet f"
replace nCH = 1 if gdename == "jerusalem (israel)"
replace nCH = 1 if gdename == "bad pyrmont d"
replace nCH = 1 if gdename == "hongkong cn"
replace nCH = 1 if gdename == "einbeck d"
replace nCH = 1 if gdename == "haegersten s"
replace nCH = 1 if gdename == "gaevle s"
replace nCH = 1 if gdename == "piraeus gr"
replace nCH = 1 if gdename == "boulder usa"
replace nCH = 1 if gdename == "saint paul de vence f"
replace nCH = 1 if gdename == "caracas venezuela"
replace nCH = 1 if gdename == "johannesburg rsa"
replace nCH = 1 if gdename == "morristown usa"
replace nCH = 1 if gdename == "cognac f"
replace nCH = 1 if gdename == "waldshut gurtweil d"
replace nCH = 1 if gdename == "st cloud (f)"
replace nCH = 1 if gdename == "hyogo j"
replace nCH = 1 if gdename == "gaertringen d"
replace nCH = 1 if gdename == "pressbaum a"
replace nCH = 1 if gdename == "quay house (guernsey) gb"
replace nCH = 1 if gdename == "pittsford usa"
replace nCH = 1 if gdename == "cantu i"
replace nCH = 1 if gdename == "hirschau d"
replace nCH = 1 if gdename == "hassocks gb"
replace nCH = 1 if gdename == "dakar"
replace nCH = 1 if gdename == "verviers b"
replace nCH = 1 if gdename == "sulzbach d"
replace nCH = 1 if gdename == "attenschwiller f"
replace nCH = 1 if gdename == "nesoya n"
replace nCH = 1 if gdename == "basilea"
replace nCH = 1 if gdename == "glion sur montreux"
replace nCH = 1 if gdename == "st. clement, jersey gb"
replace nCH = 1 if gdename == "lincoln usa"
replace nCH = 1 if gdename == "uxbridge usa"
replace nCH = 1 if gdename == "naples/florida (usa)"
replace nCH = 1 if gdename == "bjaerred s"
replace nCH = 1 if gdename == "ettelbruck l"
replace nCH = 1 if gdename == "st germain en laye (f)"
replace nCH = 1 if gdename == "horsens dk"
replace nCH = 1 if gdename == "hattem nl"
replace nCH = 1 if gdename == "eggingen d"
replace nCH = 1 if gdename == "somerset west rsa"
replace nCH = 1 if gdename == "wellesley usa"
replace nCH = 1 if gdename == "bougival f"
replace nCH = 1 if gdename == "herzlia il"
replace nCH = 1 if gdename == "morristown (usa)"
replace nCH = 1 if gdename == "feigeres f"
replace nCH = 1 if gdename == "ross ca usa"
replace nCH = 1 if gdename == "challex/ain f"
replace nCH = 1 if gdename == "chamonix f"
replace nCH = 1 if gdename == "luedinghausen d"
replace nCH = 1 if gdename == "westport (usa)"
replace nCH = 1 if gdename == "saratoga usa"
replace nCH = 1 if gdename == "orinda, cal. usa"
replace nCH = 1 if gdename == "walddorf d"
replace nCH = 1 if gdename == "ile de sark guernesey"
replace nCH = 1 if gdename == "nabern teck d"
replace nCH = 1 if gdename == "calonge e"
replace nCH = 1 if gdename == "los altos usa"
replace nCH = 1 if gdename == "ehringshausen d"
replace nCH = 1 if gdename == "st. cloud f"
replace nCH = 1 if gdename == "ninove b"
replace nCH = 1 if gdename == "gross umstadt d"
replace nCH = 1 if gdename == "marin/m. epagnier"
replace nCH = 1 if gdename == "aeschi"
replace nCH = 1 if gdename == "grenchen (so)"
replace nCH = 1 if gdename == "l'isle"
replace nCH = 1 if gdename == "menziken (ag)"
replace nCH = 1 if gdename == "veltheim"
replace nCH = 1 if gdename == "moosburg d"
replace nCH = 1 if gdename == "osnabrueck (d)"
replace nCH = 1 if gdename == "valenza i"
replace nCH = 1 if gdename == "kirchanschoering d"
replace nCH = 1 if gdename == "besancon (f)"
replace nCH = 1 if gdename == "trevise i"
replace nCH = 1 if gdename == "duesseldorf/d"
replace nCH = 1 if gdename == "orgeval f"
replace nCH = 1 if gdename == "alcobendas e"
replace nCH = 1 if gdename == "verbania i"
replace nCH = 1 if gdename == "kierspe d"
replace nCH = 1 if gdename == "ramat gan il"
replace nCH = 1 if gdename == "wolfratshausen d"
replace nCH = 1 if gdename == "jestetten (d)"
replace nCH = 1 if gdename == "badhoevedorp nl"
replace nCH = 1 if gdename == "brooklyn, new york usa"
replace nCH = 1 if gdename == "stauten d"
replace nCH = 1 if gdename == "haiterbach d"
replace nCH = 1 if gdename == "rosenheim d"
replace nCH = 1 if gdename == "dublin (irl)"
replace nCH = 1 if gdename == "bagshot gb"
replace nCH = 1 if gdename == "sergy f"
replace nCH = 1 if gdename == "kanagawa ken j"
replace nCH = 1 if gdename == "hofstetten so"
replace nCH = 1 if gdename == "edingen d"
replace nCH = 1 if gdename == "savyon isr"
replace nCH = 1 if gdename == "frankfurt/m d"
replace nCH = 1 if gdename == "overath d"
replace nCH = 1 if gdename == "tunbridge wells, kent gb"
replace nCH = 1 if gdename == "ob."
replace nCH = 1 if gdename == "renaix b"
replace nCH = 1 if gdename == "de pinte b"
replace nCH = 1 if gdename == "bueren a.d.a"
replace nCH = 1 if gdename == "dietzenbach (d)"
replace nCH = 1 if gdename == "verneuil sur igneraie f"
replace nCH = 1 if gdename == "langnau a. albis"
replace nCH = 1 if gdename == "lebanon usa"
replace nCH = 1 if gdename == "homburg d"
replace nCH = 1 if gdename == "allensbach d"
replace nCH = 1 if gdename == "rijswijk nl"
replace nCH = 1 if gdename == "denver, colorado usa"
replace nCH = 1 if gdename == "putzbrunn d"
replace nCH = 1 if gdename == "luvinate i"
replace nCH = 1 if gdename == "ruhstorf d"
replace nCH = 1 if gdename == "san jose usa"
replace nCH = 1 if gdename == "canelli d asti i"
replace nCH = 1 if gdename == "saltsjoe duvnaes s"
replace nCH = 1 if gdename == "kortenhoef nl"
replace nCH = 1 if gdename == "au d"
replace nCH = 1 if gdename == "nizza (f)"
replace nCH = 1 if gdename == "aurora, ontario cdn"
replace nCH = 1 if gdename == "lommel b"
replace nCH = 1 if gdename == "baerum n"
replace nCH = 1 if gdename == "santa cruz usa"
replace nCH = 1 if gdename == "leonberg (d)"
replace nCH = 1 if gdename == "bad duerrheim (d)"
replace nCH = 1 if gdename == "konstanz dl"
replace nCH = 1 if gdename == "veldhoven nl"
replace nCH = 1 if gdename == "dilbeek b"
replace nCH = 1 if gdename == "phoenix (usa)"
replace nCH = 1 if gdename == "kluisbergen b"
replace nCH = 1 if gdename == "sark ci"
replace nCH = 1 if gdename == "muelheim/ruhr d"
replace nCH = 1 if gdename == "penzberg d"
replace nCH = 1 if gdename == "skanoer s"
replace nCH = 1 if gdename == "hoberge uerentrup d"
replace nCH = 1 if gdename == "kuweit city (kwt)"
replace nCH = 1 if gdename == "anvers (b)"
replace nCH = 1 if gdename == "buggingen d"
replace nCH = 1 if gdename == "ottawa cdn"
replace nCH = 1 if gdename == "etang la ville (f)"
replace nCH = 1 if gdename == "limoges f"
replace nCH = 1 if gdename == "ottawa cnd"
replace nCH = 1 if gdename == "st. germain en laye f"
replace nCH = 1 if gdename == "st tropez (f)"
replace nCH = 1 if gdename == "vulbens (f)"
replace nCH = 1 if gdename == "bad oldesloe d"
replace nCH = 1 if gdename == "neuilly s seine f"
replace nCH = 1 if gdename == "selters d"
replace nCH = 1 if gdename == "mexico city mex"
replace nCH = 1 if gdename == "west hartford, conn. usa"
replace nCH = 1 if gdename == "kamakura city japan"
replace nCH = 1 if gdename == "budapest (ungarn)"
replace nCH = 1 if gdename == "reykjavik is"
replace nCH = 1 if gdename == "carrolton usa"
replace nCH = 1 if gdename == "gross gerau d"
replace nCH = 1 if gdename == "sevilla e"
replace nCH = 1 if gdename == "hinterkappelen/wohlen be"
replace nCH = 1 if gdename == "st eulalia e"
replace nCH = 1 if gdename == "hechendorf d"
replace nCH = 1 if gdename == "icking d"
replace nCH = 1 if gdename == "moessingen d"
replace nCH = 1 if gdename == "mougins f"
replace nCH = 1 if gdename == "abidjan (ci)"
replace nCH = 1 if gdename == "parma (i)"
replace nCH = 1 if gdename == "hinterbruehl a"
replace nCH = 1 if gdename == "hamm d"
replace nCH = 1 if gdename == "hamilton gb"
replace nCH = 1 if gdename == "new dehli ind"
replace nCH = 1 if gdename == "lesmo i"
replace nCH = 1 if gdename == "linkebeek b"
replace nCH = 1 if gdename == "san donato (i)"
replace nCH = 1 if gdename == "fredensborg dk"
replace nCH = 1 if gdename == "asker n"
replace nCH = 1 if gdename == "weert nl"
replace nCH = 1 if gdename == "bucarest r"
replace nCH = 1 if gdename == "cologny (ge)"
replace nCH = 1 if gdename == "zofingen (ag)"
replace nCH = 1 if gdename == "vauvert f"
replace nCH = 1 if gdename == "buckenhof d"
replace nCH = 1 if gdename == "steinen hofen d"
replace nCH = 1 if gdename == "tempe/arizona usa"
replace nCH = 1 if gdename == "watermael boitsfort b"
replace nCH = 1 if gdename == "santa rosa, cal. usa"
replace nCH = 1 if gdename == "barrington rhode island usa"
replace nCH = 1 if gdename == "salchendorf d"
replace nCH = 1 if gdename == "ramat gan (il)"
replace nCH = 1 if gdename == "schiedam nl"
replace nCH = 1 if gdename == "rosenfeld d"
replace nCH = 1 if gdename == "oberndorf/neckar d"
replace nCH = 1 if gdename == "abdullah al salem kwt"
replace nCH = 1 if gdename == "malmoe (s)"
replace nCH = 1 if gdename == "huntingdon gb"
replace nCH = 1 if gdename == "neukirchen vluyn d"
replace nCH = 1 if gdename == "overijse (b)"
replace nCH = 1 if gdename == "leersum nl"
replace nCH = 1 if gdename == "lanzo d intelvi i"
replace nCH = 1 if gdename == "neu ulm/donau d"
replace nCH = 1 if gdename == "bron f"
replace nCH = 1 if gdename == "maccagno i"
replace nCH = 1 if gdename == "heilbronn (d)"
replace nCH = 1 if gdename == "girona e"
replace nCH = 1 if gdename == "penshurst, kent gb"
replace nCH = 1 if gdename == "atlanta, georgie usa"
replace nCH = 1 if gdename == "gevelsberg (d)"
replace nCH = 1 if gdename == "markdorf d"
replace nCH = 1 if gdename == "vence f"
replace nCH = 1 if gdename == "lissabon"
replace nCH = 1 if gdename == "offenbach/bieber d"
replace nCH = 1 if gdename == "rochetrejoux f"
replace nCH = 1 if gdename == "oberderdingen d"
replace nCH = 1 if gdename == "gengenbach bergach d"
replace nCH = 1 if gdename == "waldkraiburg d"
replace nCH = 1 if gdename == "carate urio i"
replace nCH = 1 if gdename == "gibraltar gib"
replace nCH = 1 if gdename == "bellevue usa"
replace nCH = 1 if gdename == "halver d"
replace nCH = 1 if gdename == "pasadena, california usa"
replace nCH = 1 if gdename == "walluf d"
replace nCH = 1 if gdename == "suresne f"
replace nCH = 1 if gdename == "fontaine la mallet f"
replace nCH = 1 if gdename == "griessen d"
replace nCH = 1 if gdename == "lisse nl"
replace nCH = 1 if gdename == "olm l"
replace nCH = 1 if gdename == "neuwiller f"
replace nCH = 1 if gdename == "wollbach d"
replace nCH = 1 if gdename == "heemstede (nl)"
replace nCH = 1 if gdename == "st paul alberta can"
replace nCH = 1 if gdename == "karachi (pak)"
replace nCH = 1 if gdename == "arsenigo i"
replace nCH = 1 if gdename == "ehningen d"
replace nCH = 1 if gdename == "st. charles usa"
replace nCH = 1 if gdename == "march d"
replace nCH = 1 if gdename == "nassau (bs)"
replace nCH = 1 if gdename == "verrieres le buisson (f)"
replace nCH = 1 if gdename == "boston massachusetts usa"
replace nCH = 1 if gdename == "hiroschima j"
replace nCH = 1 if gdename == "bangkok (t)"
replace nCH = 1 if gdename == "gummersbach (d)"
replace nCH = 1 if gdename == "austin usa"
replace nCH = 1 if gdename == "huerth d"
replace nCH = 1 if gdename == "roma 1"
replace nCH = 1 if gdename == "eitorf d"
replace nCH = 1 if gdename == "noyarey f"
replace nCH = 1 if gdename == "bardolino i"
replace nCH = 1 if gdename == "bendestorf (d)"
replace nCH = 1 if gdename == "badevel f"
replace nCH = 1 if gdename == "granby cdn"
replace nCH = 1 if gdename == "datteln d"
replace nCH = 1 if gdename == "hohentengen (d)"
replace nCH = 1 if gdename == "winnenden rfg"
replace nCH = 1 if gdename == "orpesa (e)"
replace nCH = 1 if gdename == "naarden/nl"
replace nCH = 1 if gdename == "wentorf d"
replace nCH = 1 if gdename == "st. georgen (d)"
replace nCH = 1 if gdename == "drachten nl"
replace nCH = 1 if gdename == "ammerthal d"
replace nCH = 1 if gdename == "loeffingen d"
replace nCH = 1 if gdename == "bad kreuznach d"
replace nCH = 1 if gdename == "hamilton bermuda"
replace nCH = 1 if gdename == "chagrin falls (usa)"
replace nCH = 1 if gdename == "newton usa"
replace nCH = 1 if gdename == "neufraunhofen d"
replace nCH = 1 if gdename == "st lambert cdn"
replace nCH = 1 if gdename == "orlando usa"
replace nCH = 1 if gdename == "koenigstein"
replace nCH = 1 if gdename == "arbent f"
replace nCH = 1 if gdename == "krumbach d"
replace nCH = 1 if gdename == "s gravenhage nl"
replace nCH = 1 if gdename == "krailling/muenchen d"
replace nCH = 1 if gdename == "radolfszell d"
replace nCH = 1 if gdename == "santpoort nl"
replace nCH = 1 if gdename == "shaker heights usa"
replace nCH = 1 if gdename == "seoul (korea)"
replace nCH = 1 if gdename == "andover gb"
replace nCH = 1 if gdename == "caluire et cuire f"
replace nCH = 1 if gdename == "san donato milanese (i)"
replace nCH = 1 if gdename == "louisville (usa)"
replace nCH = 1 if gdename == "schwalbach am taunus d"
replace nCH = 1 if gdename == "vougeot f"
replace nCH = 1 if gdename == "neckarsteinach d"
replace nCH = 1 if gdename == "montrouge f"
replace nCH = 1 if gdename == "siegburg d"
replace nCH = 1 if gdename == "boege f"
replace nCH = 1 if gdename == "schorndorf weiler d"
replace nCH = 1 if gdename == "cologne d"
replace nCH = 1 if gdename == "lohmar/koeln d"
replace nCH = 1 if gdename == "malland"
replace nCH = 1 if gdename == "eckental d"
replace nCH = 1 if gdename == "bad homburg von der hoehe d"
replace nCH = 1 if gdename == "insel reichenau brd"
replace nCH = 1 if gdename == "borken dillich d"
replace nCH = 1 if gdename == "cincinnati ohio usa"
replace nCH = 1 if gdename == "summit/nj usa"
replace nCH = 1 if gdename == "hosle n"
replace nCH = 1 if gdename == "gelsenkirchen buer d"
replace nCH = 1 if gdename == "gent nl"
replace nCH = 1 if gdename == "lake bluff, illinois usa"
replace nCH = 1 if gdename == "esher surrey gb"
replace nCH = 1 if gdename == "ankara (tr)"
replace nCH = 1 if gdename == "lainate i"
replace nCH = 1 if gdename == "chiavenna i"
replace nCH = 1 if gdename == "mamaroneck ny usa"
replace nCH = 1 if gdename == "spiesen d"
replace nCH = 1 if gdename == "woking/surrey gb"
replace nCH = 1 if gdename == "kitzbuehel a"
replace nCH = 1 if gdename == "kressbronn (d)"
replace nCH = 1 if gdename == "frederikshavn dk"
replace nCH = 1 if gdename == "kronstorf a"
replace nCH = 1 if gdename == "voorburg nl"
replace nCH = 1 if gdename == "rutesheim d"
replace nCH = 1 if gdename == "nenzing a"
replace nCH = 1 if gdename == "ville d avray (f)"
replace nCH = 1 if gdename == "ettenheim/altdorf d"
replace nCH = 1 if gdename == "meckenheim d"
replace nCH = 1 if gdename == "zutphen nl"
replace nCH = 1 if gdename == "groebenzell d"
replace nCH = 1 if gdename == "new canaan, conn. usa"
replace nCH = 1 if gdename == "manama/bahrain"
replace nCH = 1 if gdename == "aabenraa dk"
replace nCH = 1 if gdename == "chartres f"
replace nCH = 1 if gdename == "wiler f"
replace nCH = 1 if gdename == "blaustein d"
replace nCH = 1 if gdename == "puebla (mex)"
replace nCH = 1 if gdename == "megeve f"
replace nCH = 1 if gdename == "kuala lumpur"
replace nCH = 1 if gdename == "berks gb"
replace nCH = 1 if gdename == "chambrelien/rochefort"
replace nCH = 1 if gdename == "ammann jor"
replace nCH = 1 if gdename == "lacq/basses pyrenees f"
replace nCH = 1 if gdename == "asarum s"
replace nCH = 1 if gdename == "usingen d"
replace nCH = 1 if gdename == "moosach altenburg d"
replace nCH = 1 if gdename == "bad soden (d)"
replace nCH = 1 if gdename == "teddington gb"
replace nCH = 1 if gdename == "stockdorf (d)"
replace nCH = 1 if gdename == "hilzingen d"
replace nCH = 1 if gdename == "valdagno (i)"
replace nCH = 1 if gdename == "los altos, california usa"
replace nCH = 1 if gdename == "palaiseau f"
replace nCH = 1 if gdename == "siegen (d)"
replace nCH = 1 if gdename == "rheinfelden brd"
replace nCH = 1 if gdename == "frielendorf d"
replace nCH = 1 if gdename == "bergneustadt d"
replace nCH = 1 if gdename == "schotten (d)"
replace nCH = 1 if gdename == "arequipa perou"
replace nCH = 1 if gdename == "los angeles, usa"
replace nCH = 1 if gdename == "lisbonne"
replace nCH = 1 if gdename == "neuilly sur seine/f"
replace nCH = 1 if gdename == "zaltbommel nl"
replace nCH = 1 if gdename == "la chaux"
replace nCH = 1 if gdename == "nl"
replace nCH = 1 if gdename == "bruchkoebel d"
replace nCH = 1 if gdename == "leusden nl"
replace nCH = 1 if gdename == "sonntag a"
replace nCH = 1 if gdename == "ryad (sa)"
replace nCH = 1 if gdename == "taipei taiwan"
replace nCH = 1 if gdename == "waltenheim f"
replace nCH = 1 if gdename == "rosengarten (d)"
replace nCH = 1 if gdename == "henley on thames (gb)"
replace nCH = 1 if gdename == "nordborg dk"
replace nCH = 1 if gdename == "zell unterentersbach d"
replace nCH = 1 if gdename == "west horsley, surrey gb"
replace nCH = 1 if gdename == "woodside (california), usa"
replace nCH = 1 if gdename == "miami beach usa"
replace nCH = 1 if gdename == "villingen d"
replace nCH = 1 if gdename == "sterzing i"
replace nCH = 1 if gdename == "glinde d"
replace nCH = 1 if gdename == "saint die f"
replace nCH = 1 if gdename == "san pedro de alcantara e"
replace nCH = 1 if gdename == "margencel f"
replace nCH = 1 if gdename == "wanfried/h. d"
replace nCH = 1 if gdename == "bremen brd"
replace nCH = 1 if gdename == "salt lake city, utha usa"
replace nCH = 1 if gdename == "sevenoaks, kent (gb)"
replace nCH = 1 if gdename == "ravels (b)"
replace nCH = 1 if gdename == "beyrouth lbn"
replace nCH = 1 if gdename == "gaensbach d"
replace nCH = 1 if gdename == "munkebo dk"
replace nCH = 1 if gdename == "pfullendorf d"
replace nCH = 1 if gdename == "schwerte d"
replace nCH = 1 if gdename == "seeheim d"
replace nCH = 1 if gdename == "chereng f"
replace nCH = 1 if gdename == "fujisawa j"
replace nCH = 1 if gdename == "bad ems d"
replace nCH = 1 if gdename == "kiffis f"
replace nCH = 1 if gdename == "sausalito usa"
replace nCH = 1 if gdename == "castiglione delle stivere i"
replace nCH = 1 if gdename == "erkrath (d)"
replace nCH = 1 if gdename == "betchworth"
replace nCH = 1 if gdename == "st mande f"
replace nCH = 1 if gdename == "guenzburg d"
replace nCH = 1 if gdename == "muellheim an der ruhr d"
replace nCH = 1 if gdename == "hoquiam, washington usa"
replace nCH = 1 if gdename == "eislingen d"
replace nCH = 1 if gdename == "hamilton (bermuda)"
replace nCH = 1 if gdename == "new canaan/connect. usa"
replace nCH = 1 if gdename == "le perreux f"
replace nCH = 1 if gdename == "pregny f"
replace nCH = 1 if gdename == "wedemark d"
replace nCH = 1 if gdename == "mettmann (d)"
replace nCH = 1 if gdename == "passau d"
replace nCH = 1 if gdename == "vaellingby s"
replace nCH = 1 if gdename == "erps kwerps b"
replace nCH = 1 if gdename == "victoria aus"
replace nCH = 1 if gdename == "ribe dk"
replace nCH = 1 if gdename == "winterbach bei stuttgart d"
replace nCH = 1 if gdename == "pfaffenhofen d"
replace nCH = 1 if gdename == "boeblingen"
replace nCH = 1 if gdename == "nendeln fl"
replace nCH = 1 if gdename == "winterthur i"
replace nCH = 1 if gdename == "porto cristo e"
replace nCH = 1 if gdename == "kronbuehl/wittenbach"
replace nCH = 1 if gdename == "muespach f"
replace nCH = 1 if gdename == "buenos aires arg"
replace nCH = 1 if gdename == "denver (usa)"
replace nCH = 1 if gdename == "grosshesslohe muenchen d"
replace nCH = 1 if gdename == "lauterbach d"
replace nCH = 1 if gdename == "ipswich gb"
replace nCH = 1 if gdename == "modena (i)"
replace nCH = 1 if gdename == "kernen d"
replace nCH = 1 if gdename == "muehlheim/ruhr d"
replace nCH = 1 if gdename == "schorndorf (d)"
replace nCH = 1 if gdename == "altenstadt d"
replace nCH = 1 if gdename == "neuhausen a. rhf"
replace nCH = 1 if gdename == "gerpinnes b"
replace nCH = 1 if gdename == "hangenbieten f"
replace nCH = 1 if gdename == "reggio calabria i"
replace nCH = 1 if gdename == "saarbruecken (d)"
replace nCH = 1 if gdename == "oberkirch (d)"
replace nCH = 1 if gdename == "igersheim d"
replace nCH = 1 if gdename == "broadhambury gb"
replace nCH = 1 if gdename == "ludwigshafen"
replace nCH = 1 if gdename == "au bei freiburg d"
replace nCH = 1 if gdename == "hamburg sassel d"
replace nCH = 1 if gdename == "muecke d"
replace nCH = 1 if gdename == "dubai sv"
replace nCH = 1 if gdename == "gebhardshain d"
replace nCH = 1 if gdename == "waasmunster b"
replace nCH = 1 if gdename == "chaudfontaine b"
replace nCH = 1 if gdename == "frangy f"
replace nCH = 1 if gdename == "fagersta s"
replace nCH = 1 if gdename == "gustav heinemann ufer 92, koeln d"
replace nCH = 1 if gdename == "bilthoven"
replace nCH = 1 if gdename == "northridge, cal. usa"
replace nCH = 1 if gdename == "scottsdale usa"
replace nCH = 1 if gdename == "apcn c62, porto christo novo, mallorca e"
replace nCH = 1 if gdename == "neviges d"
replace nCH = 1 if gdename == "gross bieberau d"
replace nCH = 1 if gdename == "voerendaal, limburg nl"
replace nCH = 1 if gdename == "lexington, n.c. (usa)"
replace nCH = 1 if gdename == "umtata rsr"
replace nCH = 1 if gdename == "burgbrohl d"
replace nCH = 1 if gdename == "cleveland ohio usa"
replace nCH = 1 if gdename == "saeffle s"
replace nCH = 1 if gdename == "boernsen d"
replace nCH = 1 if gdename == "greenville usa"
replace nCH = 1 if gdename == "ryad sa"
replace nCH = 1 if gdename == "feldafging d"
replace nCH = 1 if gdename == "weissenhorn d"
replace nCH = 1 if gdename == "purley/surrey gb"
replace nCH = 1 if gdename == "zwijndrecht nl"
replace nCH = 1 if gdename == "hannover dl"
replace nCH = 1 if gdename == "limal b"
replace nCH = 1 if gdename == "odenthal steinhaus d"
replace nCH = 1 if gdename == "luxembourg b"
replace nCH = 1 if gdename == "wilnsdorf d"
replace nCH = 1 if gdename == "mastholte d"
replace nCH = 1 if gdename == "hildesheim d"
replace nCH = 1 if gdename == "grimstad (n)"
replace nCH = 1 if gdename == "germering d"
replace nCH = 1 if gdename == "ohain (b)"
replace nCH = 1 if gdename == "yardley usa"
replace nCH = 1 if gdename == "orchard park, new york usa"
replace nCH = 1 if gdename == "wezembeek oppem b"
replace nCH = 1 if gdename == "malacca, malaysia"
replace nCH = 1 if gdename == "spijkenisse nl"
replace nCH = 1 if gdename == "wyhlen d"
replace nCH = 1 if gdename == "maiche f"
replace nCH = 1 if gdename == "bridge of allan gb"
replace nCH = 1 if gdename == "tahoma usa"
replace nCH = 1 if gdename == "autrey les gray f"
replace nCH = 1 if gdename == "stanmore, middlesex gb"
replace nCH = 1 if gdename == "lesigny f"
replace nCH = 1 if gdename == "danbury usa"
replace nCH = 1 if gdename == "heroldsberg d"
replace nCH = 1 if gdename == "stafford gb"
replace nCH = 1 if gdename == "wissous f"
replace nCH = 1 if gdename == "uhldingen d"
replace nCH = 1 if gdename == "reichelsheim d"
replace nCH = 1 if gdename == "ferrara i"
replace nCH = 1 if gdename == "well am rhein d"
replace nCH = 1 if gdename == "crozet f"
replace nCH = 1 if gdename == "versonnex/ain f"
replace nCH = 1 if gdename == "zwijndrecht b"
replace nCH = 1 if gdename == "esery f"
replace nCH = 1 if gdename == "prevessin (f)"
replace nCH = 1 if gdename == "cobham (gb)"
replace nCH = 1 if gdename == "nancy (f)"
replace nCH = 1 if gdename == "landau/pfalz d"
replace nCH = 1 if gdename == "rietberg d"
replace nCH = 1 if gdename == "kehl (d)"
replace nCH = 1 if gdename == "besozzo i"
replace nCH = 1 if gdename == "yerres f"
replace nCH = 1 if gdename == "oehringen d"
replace nCH = 1 if gdename == "lingen d"
replace nCH = 1 if gdename == "cachan f"
replace nCH = 1 if gdename == "birkerod dk"
replace nCH = 1 if gdename == "biot f"
replace nCH = 1 if gdename == "saskatoon cdn"
replace nCH = 1 if gdename == "east molesey gb"
replace nCH = 1 if gdename == "altea e"
replace nCH = 1 if gdename == "hamilton bermudas"
replace nCH = 1 if gdename == "coesfeld d"
replace nCH = 1 if gdename == "seattle, washington usa"
replace nCH = 1 if gdename == "iggesund (s)"
replace nCH = 1 if gdename == "monte carlo, monaco"
replace nCH = 1 if gdename == "pfinztal d"
replace nCH = 1 if gdename == "malakoff f"
replace nCH = 1 if gdename == "obrigheim d"
replace nCH = 1 if gdename == "east horsley, surrey gb"
replace nCH = 1 if gdename == "forest row gb"
replace nCH = 1 if gdename == "schoenau d"
replace nCH = 1 if gdename == "almere nl"
replace nCH = 1 if gdename == "nassau usa"
replace nCH = 1 if gdename == "piraeus (gr)"
replace nCH = 1 if gdename == "kirschweiler d"
replace nCH = 1 if gdename == "haren nl"
replace nCH = 1 if gdename == "ashbourne/derbyshire gb"
replace nCH = 1 if gdename == "uggiate trevano i"
replace nCH = 1 if gdename == "algier (dz)"
replace nCH = 1 if gdename == "guatemala"
replace nCH = 1 if gdename == "sassenheim nl"
replace nCH = 1 if gdename == "bremervoerde d"
replace nCH = 1 if gdename == "roeselare b"
replace nCH = 1 if gdename == "hof, saale d"
replace nCH = 1 if gdename == "remagen oberwinter d"
replace nCH = 1 if gdename == "herlev dk"
replace nCH = 1 if gdename == "gatteo i"
replace nCH = 1 if gdename == "rio de janeiro (brasilien)"
replace nCH = 1 if gdename == "san francisco, cal., usa"
replace nCH = 1 if gdename == "djeddah (sa)"
replace nCH = 1 if gdename == "setauket n.y. usa"
replace nCH = 1 if gdename == "les terrasses 9, aix les bains f"
replace nCH = 1 if gdename == "moskau, russland"
replace nCH = 1 if gdename == "st. erhard/knutwil"
replace nCH = 1 if gdename == "nervesa della battaglia i"
replace nCH = 1 if gdename == "schallbach d"
replace nCH = 1 if gdename == "oerlinghausen (d)"
replace nCH = 1 if gdename == "kuenzelsau/wuertt. brd"
replace nCH = 1 if gdename == "miesbach (d)"
replace nCH = 1 if gdename == "herent b"
replace nCH = 1 if gdename == "helsinki (finnland)"
replace nCH = 1 if gdename == "sasbach bei achern d"
replace nCH = 1 if gdename == "herne d"
replace nCH = 1 if gdename == "kleinneusiedel a"
replace nCH = 1 if gdename == "sofia (bg)"
replace nCH = 1 if gdename == "westport, conn. usa"
replace nCH = 1 if gdename == "barendrecht nl"
replace nCH = 1 if gdename == "istanbul tk"
replace nCH = 1 if gdename == "palagano i"
replace nCH = 1 if gdename == "korb d"
replace nCH = 1 if gdename == "kuenzelsau d"
replace nCH = 1 if gdename == "cassina rizzardi i"
replace nCH = 1 if gdename == "turin/i"
replace nCH = 1 if gdename == "summit, new jersey usa"
replace nCH = 1 if gdename == "agoura, calif. usa"
replace nCH = 1 if gdename == "mclean usa"
replace nCH = 1 if gdename == "montesson f"
replace nCH = 1 if gdename == "londres/gb"
replace nCH = 1 if gdename == "letchworth gb"
replace nCH = 1 if gdename == "scrofiano i"
replace nCH = 1 if gdename == "arenzano i"
replace nCH = 1 if gdename == "buedingen d"
replace nCH = 1 if gdename == "dallas, texas (usa)"
replace nCH = 1 if gdename == "norwich gb"
replace nCH = 1 if gdename == "santa monica ca usa"
replace nCH = 1 if gdename == "uitikon a.a"
replace nCH = 1 if gdename == "tervuren (b)"
replace nCH = 1 if gdename == "bad bentheim d"
replace nCH = 1 if gdename == "roedental d"
replace nCH = 1 if gdename == "freising d"
replace nCH = 1 if gdename == "voorthuizen (nl)"
replace nCH = 1 if gdename == "tirano i"
replace nCH = 1 if gdename == "marl d"
replace nCH = 1 if gdename == "louveciennes (f)"
replace nCH = 1 if gdename == "newhall/californie usa"
replace nCH = 1 if gdename == "salles chimay b"
replace nCH = 1 if gdename == "annemasse/hte savoie, f"
replace nCH = 1 if gdename == "fairport, new york usa"
replace nCH = 1 if gdename == "buchen d"
replace nCH = 1 if gdename == "york p.a. usa"
replace nCH = 1 if gdename == "sudbury/mass. usa"
replace nCH = 1 if gdename == "st galmier f"
replace nCH = 1 if gdename == "weston (usa)"
replace nCH = 1 if gdename == "molln a"
replace nCH = 1 if gdename == "la haye nl"
replace nCH = 1 if gdename == "dobel d"
replace nCH = 1 if gdename == "irvine usa"
replace nCH = 1 if gdename == "denver colorado usa"
replace nCH = 1 if gdename == "sint niklaas b"
replace nCH = 1 if gdename == "st aubin/st a. sauges"
replace nCH = 1 if gdename == "carson city usa"
replace nCH = 1 if gdename == "corenc f"
replace nCH = 1 if gdename == "recklingshausen d"
replace nCH = 1 if gdename == "premia mar e"
replace nCH = 1 if gdename == "neustadt an der weinstrasse d"
replace nCH = 1 if gdename == "bournemouth gb"
replace nCH = 1 if gdename == "stockdorf d"
replace nCH = 1 if gdename == "aelmhult s"
replace nCH = 1 if gdename == "ewersbach d"
replace nCH = 1 if gdename == "loerrach dl"
replace nCH = 1 if gdename == "oberschleissheim d"
replace nCH = 1 if gdename == "murg haenner d"
replace nCH = 1 if gdename == "shipley gb"
replace nCH = 1 if gdename == "silver spring, maryland usa"
replace nCH = 1 if gdename == "cluny f"
replace nCH = 1 if gdename == "manila rp"
replace nCH = 1 if gdename == "ruedersheim d"
replace nCH = 1 if gdename == "atlanta ga usa"
replace nCH = 1 if gdename == "darien, connecticut (usa)"
replace nCH = 1 if gdename == "esslingen a. n. d"
replace nCH = 1 if gdename == "schopfheim (d)"
replace nCH = 1 if gdename == "chenevert merignac f"
replace nCH = 1 if gdename == "auderghem b"
replace nCH = 1 if gdename == "stockholm/s"
replace nCH = 1 if gdename == "lebanon n.j. (usa)"
replace nCH = 1 if gdename == "casate nuovo i"
replace nCH = 1 if gdename == "ternat b"
replace nCH = 1 if gdename == "la hulpe b"
replace nCH = 1 if gdename == "beijing cn"
replace nCH = 1 if gdename == "althausen d"
replace nCH = 1 if gdename == "bailly (yvelines) f"
replace nCH = 1 if gdename == "vaglia i"
replace nCH = 1 if gdename == "durban (rsa)"
replace nCH = 1 if gdename == "thueringen a"
replace nCH = 1 if gdename == "schondorf am ammersee d"
replace nCH = 1 if gdename == "aschau d"
replace nCH = 1 if gdename == "calpe, alicante e"
replace nCH = 1 if gdename == "plettenberg ohle d"
replace nCH = 1 if gdename == "haifa il"
replace nCH = 1 if gdename == "malta m"
replace nCH = 1 if gdename == "riethoven nl"
replace nCH = 1 if gdename == "neuhausen a.rhf"
replace nCH = 1 if gdename == "la molina pe"
replace nCH = 1 if gdename == "ruechersdorf d"
replace nCH = 1 if gdename == "kreuzlangen"
replace nCH = 1 if gdename == "leeds (gb)"
replace nCH = 1 if gdename == "rheinmuenster d"
replace nCH = 1 if gdename == "geispitzen f"
replace nCH = 1 if gdename == "neerpelt b"
replace nCH = 1 if gdename == "barzano i"
replace nCH = 1 if gdename == "los altos hills usa"
replace nCH = 1 if gdename == "ladenburg d"
replace nCH = 1 if gdename == "bad honnef d"
replace nCH = 1 if gdename == "plano texas usa"
replace nCH = 1 if gdename == "damaskus sr"
replace nCH = 1 if gdename == "piedmont, california usa"
replace nCH = 1 if gdename == "nassau d"
replace nCH = 1 if gdename == "chigny sur morges"
replace nCH = 1 if gdename == "schliengen d"
replace nCH = 1 if gdename == "wilnsdorf obersdorf d"
replace nCH = 1 if gdename == "gaienhofen d"
replace nCH = 1 if gdename == "casciago i"
replace nCH = 1 if gdename == "lahnstein d"
replace nCH = 1 if gdename == "mozzo i"
replace nCH = 1 if gdename == "essex, fells, new jersey usa"
replace nCH = 1 if gdename == "glattbach d"
replace nCH = 1 if gdename == "laupheim brd"
replace nCH = 1 if gdename == "noordwijk nl"
replace nCH = 1 if gdename == "fischbach d"
replace nCH = 1 if gdename == "st egreve f"
replace nCH = 1 if gdename == "koerich l"
replace nCH = 1 if gdename == "braunschweig (d)"
replace nCH = 1 if gdename == "senlis f"
replace nCH = 1 if gdename == "sydney (aus)"
replace nCH = 1 if gdename == "beyrut rl"
replace nCH = 1 if gdename == "niedernhall d"
replace nCH = 1 if gdename == "orleans f"
replace nCH = 1 if gdename == "stadlern (d)"
replace nCH = 1 if gdename == "kirjala sf"
replace nCH = 1 if gdename == "tuebingen (d)"
replace nCH = 1 if gdename == "pittsburg usa"
replace nCH = 1 if gdename == "leicester (gb)"
replace nCH = 1 if gdename == "svedala s"
replace nCH = 1 if gdename == "sorengo i"
replace nCH = 1 if gdename == "angleur b"
replace nCH = 1 if gdename == "crowthorne, berkshire gb"
replace nCH = 1 if gdename == "kungsangen s"
replace nCH = 1 if gdename == "minneapolis/minnesota usa"
replace nCH = 1 if gdename == "teheran (iran)"
replace nCH = 1 if gdename == "salt lake city (usa)"
replace nCH = 1 if gdename == "much d"
replace nCH = 1 if gdename == "alsemberg b"
replace nCH = 1 if gdename == "alpbach a"
replace nCH = 1 if gdename == "dayton usa"
replace nCH = 1 if gdename == "lanark gb"
replace nCH = 1 if gdename == "bacolod city pi"
replace nCH = 1 if gdename == "washington dc(usa)"
replace nCH = 1 if gdename == "chably, quebec cdn"
replace nCH = 1 if gdename == "monrovia, liberia"
replace nCH = 1 if gdename == "schoten anvers"
replace nCH = 1 if gdename == "l orient/le chenit"
replace nCH = 1 if gdename == "scarsdale, new york usa"
replace nCH = 1 if gdename == "westerhorn d"
replace nCH = 1 if gdename == "espergaerde dk"
replace nCH = 1 if gdename == "brasilia br"
replace nCH = 1 if gdename == "virginia water gb"
replace nCH = 1 if gdename == "bronxville usa"
replace nCH = 1 if gdename == "drolshagen d"
replace nCH = 1 if gdename == "croissy sur seine (f)"
replace nCH = 1 if gdename == "epinay f"
replace nCH = 1 if gdename == "new orleans usa"
replace nCH = 1 if gdename == "christ church, barbados"
replace nCH = 1 if gdename == "eglharting d"
replace nCH = 1 if gdename == "mount kisco usa"
replace nCH = 1 if gdename == "basilio i"
replace nCH = 1 if gdename == "chaville f"
replace nCH = 1 if gdename == "palma de mallorca"
replace nCH = 1 if gdename == "chantilly f"
replace nCH = 1 if gdename == "makati rp"
replace nCH = 1 if gdename == "kampala uganda"
replace nCH = 1 if gdename == "stockheim d"
replace nCH = 1 if gdename == "leersum (nl)"
replace nCH = 1 if gdename == "reinbeck d"
replace nCH = 1 if gdename == "courbevoie (f)"
replace nCH = 1 if gdename == "pordenone i"
replace nCH = 1 if gdename == "tamm d"
replace nCH = 1 if gdename == "angouleme f"
replace nCH = 1 if gdename == "knokke heist b"
replace nCH = 1 if gdename == "amersfoort nl"
replace nCH = 1 if gdename == "neuilly/f"
replace nCH = 1 if gdename == "norwegen"
replace nCH = 1 if gdename == "singen hausen d"
replace nCH = 1 if gdename == "ingelheim d"
replace nCH = 1 if gdename == "sprendlingen d"
replace nCH = 1 if gdename == "zwiesel d"
replace nCH = 1 if gdename == "toijale sf"
replace nCH = 1 if gdename == "kungsbacka s"
replace nCH = 1 if gdename == "oxfordshire gb"
replace nCH = 1 if gdename == "namur b"
replace nCH = 1 if gdename == "broomall, penns. usa"
replace nCH = 1 if gdename == "san pancrazio i"
replace nCH = 1 if gdename == "kirchheim/teck d"
replace nCH = 1 if gdename == "winston salem (usa)"
replace nCH = 1 if gdename == "heppenheim d"
replace nCH = 1 if gdename == "fairfield, connecticut usa"
replace nCH = 1 if gdename == "kembs f"
replace nCH = 1 if gdename == "brenthonne f"
replace nCH = 1 if gdename == "elmshorn d"
replace nCH = 1 if gdename == "kaltenbach a"
replace nCH = 1 if gdename == "thonon les bain f"
replace nCH = 1 if gdename == "schweinfurt"
replace nCH = 1 if gdename == "pasadena (usa)"
replace nCH = 1 if gdename == "zoug"
replace nCH = 1 if gdename == "krasnodar gus"
replace nCH = 1 if gdename == "sax/sennwald"
replace nCH = 1 if gdename == "kirchheim teck (d)"
replace nCH = 1 if gdename == "tiana e"
replace nCH = 1 if gdename == "hong kong gb"
replace nCH = 1 if gdename == "olbia i"
replace nCH = 1 if gdename == "margaux (f)"
replace nCH = 1 if gdename == "kenilworth gb"
replace nCH = 1 if gdename == "redwitz mannsgereuth d"
replace nCH = 1 if gdename == "gorxheimertal d"
replace nCH = 1 if gdename == "heidenrod d"
replace nCH = 1 if gdename == "sark (gb)"
replace nCH = 1 if gdename == "herford d"
replace nCH = 1 if gdename == "neuhausen d"
replace nCH = 1 if gdename == "sosva (dom)"
replace nCH = 1 if gdename == "vestenbergsgreuth d"
replace nCH = 1 if gdename == "johannesbg"
replace nCH = 1 if gdename == "pers jussy f"
replace nCH = 1 if gdename == "rheinbach d"
replace nCH = 1 if gdename == "kokkedal dk"
replace nCH = 1 if gdename == "muehlheim (d)"
replace nCH = 1 if gdename == "auburn (usa)"
replace nCH = 1 if gdename == "trossingen d"
replace nCH = 1 if gdename == "altdorf d"
replace nCH = 1 if gdename == "attendorn d"
replace nCH = 1 if gdename == "muecke, nieder ohmen d"
replace nCH = 1 if gdename == "erbanno i"
replace nCH = 1 if gdename == "baederich"
replace nCH = 1 if gdename == "haibach d"
replace nCH = 1 if gdename == "moscou udssr"
replace nCH = 1 if gdename == "serignan f"
replace nCH = 1 if gdename == "vogt bei ravensburg brd"
replace nCH = 1 if gdename == "meaux f"
replace nCH = 1 if gdename == "baele"
replace nCH = 1 if gdename == "abidjan (cote d ivoire)"
replace nCH = 1 if gdename == "hounslow gb"
replace nCH = 1 if gdename == "alayor menorca e"
replace nCH = 1 if gdename == "peyrens f"
replace nCH = 1 if gdename == "chevychase maryland usa"
replace nCH = 1 if gdename == "deisenhofen d"
replace nCH = 1 if gdename == "beverley hills usa"
replace nCH = 1 if gdename == "salzkotten d"
replace nCH = 1 if gdename == "montgomery, alabama usa"
replace nCH = 1 if gdename == "nieuw vennep nl"
replace nCH = 1 if gdename == "tadworth gb"
replace nCH = 1 if gdename == "princeton new jersey usa"
replace nCH = 1 if gdename == "klosterneuburg (a)"
replace nCH = 1 if gdename == "krailling d"
replace nCH = 1 if gdename == "willingen d"
replace nCH = 1 if gdename == "gmunden a"
replace nCH = 1 if gdename == "muelheim (d)"
replace nCH = 1 if gdename == "purley (gb)"
replace nCH = 1 if gdename == "mijdrecht nl"
replace nCH = 1 if gdename == "racine usa"
replace nCH = 1 if gdename == "schagen nl"
replace nCH = 1 if gdename == "morteau (f)"
replace nCH = 1 if gdename == "wilrijk (b)"
replace nCH = 1 if gdename == "hasselt (b)"
replace nCH = 1 if gdename == "isernhagen d"
replace nCH = 1 if gdename == "sandefjord n"
replace nCH = 1 if gdename == "cantello i"
replace nCH = 1 if gdename == "warwick gb"
replace nCH = 1 if gdename == "sandefjord (n)"
replace nCH = 1 if gdename == "merignac f"
replace nCH = 1 if gdename == "sjovagen taby s"
replace nCH = 1 if gdename == "doorn nl"
replace nCH = 1 if gdename == "germignaga"
replace nCH = 1 if gdename == "servance f"
replace nCH = 1 if gdename == "louisville, kentucky usa"
replace nCH = 1 if gdename == "bad wiessee d"
replace nCH = 1 if gdename == "saint cloud"
replace nCH = 1 if gdename == "raanana il"
replace nCH = 1 if gdename == "trebbo i"
replace nCH = 1 if gdename == "le plessis robinson f"
replace nCH = 1 if gdename == "freeport bahamas"
replace nCH = 1 if gdename == "castelletto di brenzone i"
replace nCH = 1 if gdename == "marietta georgia usa"
replace nCH = 1 if gdename == "niederhausen d"
replace nCH = 1 if gdename == "klagenfurt"
replace nCH = 1 if gdename == "reno usa"
replace nCH = 1 if gdename == "calgary, alberta cdn"
replace nCH = 1 if gdename == "bad zwischenhahn/ekern d"
replace nCH = 1 if gdename == "guatemala city gca"
replace nCH = 1 if gdename == "rothenbuch d"
replace nCH = 1 if gdename == "nidderau d"
replace nCH = 1 if gdename == "umhlanga rocks za"
replace nCH = 1 if gdename == "leipheim d"
replace nCH = 1 if gdename == "kauniainen sf"
replace nCH = 1 if gdename == "ravenna i"
replace nCH = 1 if gdename == "gummerbach d"
replace nCH = 1 if gdename == "nova milanese i"
replace nCH = 1 if gdename == "affalterbach d"
replace nCH = 1 if gdename == "togliatti rus"
replace nCH = 1 if gdename == "holland, michigan usa"
replace nCH = 1 if gdename == "wellington (nz)"
replace nCH = 1 if gdename == "albstadt (d)"
replace nCH = 1 if gdename == "oldham gb"
replace nCH = 1 if gdename == "kellmuenz d"
replace nCH = 1 if gdename == "codroipo i"
replace nCH = 1 if gdename == "millbrook lane usa"
replace nCH = 1 if gdename == "giesshuebl (d)"
replace nCH = 1 if gdename == "reichenbach d"
replace nCH = 1 if gdename == "imst a"
replace nCH = 1 if gdename == "quatre bornes mus"
replace nCH = 1 if gdename == "jeddah uae"
replace nCH = 1 if gdename == "mantagnat f"
replace nCH = 1 if gdename == "geradstetten d"
replace nCH = 1 if gdename == "gerlingen (d)"
replace nCH = 1 if gdename == "pfullendorf (d)"
replace nCH = 1 if gdename == "koblach a"
replace nCH = 1 if gdename == "cugnaux f"
replace nCH = 1 if gdename == "fujisawa city japan"
replace nCH = 1 if gdename == "muskegon usa"
replace nCH = 1 if gdename == "son servera e"
replace nCH = 1 if gdename == "ingolstadt d"
replace nCH = 1 if gdename == "rheda wiedenbrueck d"
replace nCH = 1 if gdename == "vignola i"
replace nCH = 1 if gdename == "brugherio i"
replace nCH = 1 if gdename == "kaarina sf"
replace nCH = 1 if gdename == "arnhem (nl)"
replace nCH = 1 if gdename == "mexico (mex)"
replace nCH = 1 if gdename == "groesbeek (nl)"
replace nCH = 1 if gdename == "oakland usa"
replace nCH = 1 if gdename == "abidjan/cote d ivoire"
replace nCH = 1 if gdename == "south carolina usa"
replace nCH = 1 if gdename == "brislington bristol gb"
replace nCH = 1 if gdename == "witten d"
replace nCH = 1 if gdename == "born nl"
replace nCH = 1 if gdename == "kraainem (b)"
replace nCH = 1 if gdename == "lier n"
replace nCH = 1 if gdename == "engelwood, n.j. usa"
replace nCH = 1 if gdename == "rosenau f"
replace nCH = 1 if gdename == "st cannat/ bouches du rhone f"
replace nCH = 1 if gdename == "ghirla i"
replace nCH = 1 if gdename == "pont de la morge"
replace nCH = 1 if gdename == "st peter port, guernsey gb"
replace nCH = 1 if gdename == "eschwege brd"
replace nCH = 1 if gdename == "el cerrito usa"
replace nCH = 1 if gdename == "fritzens a"
replace nCH = 1 if gdename == "leidschendam nl"
replace nCH = 1 if gdename == "naperville, illinois usa"
replace nCH = 1 if gdename == "esbo sf"
replace nCH = 1 if gdename == "eschweiler d"
replace nCH = 1 if gdename == "taichung rc"
replace nCH = 1 if gdename == "dreieich buchschlag d"
replace nCH = 1 if gdename == "thorberg krauchthal"
replace nCH = 1 if gdename == "kanada"
replace nCH = 1 if gdename == "hildenborough gb"
replace nCH = 1 if gdename == "lugrin f"
replace nCH = 1 if gdename == "lama kara (tg)"
replace nCH = 1 if gdename == "wetzlar (d)"
replace nCH = 1 if gdename == "waidhofen an der ybbs a"
replace nCH = 1 if gdename == "indianapolis (usa)"
replace nCH = 1 if gdename == "breukelen nl"
replace nCH = 1 if gdename == "betchworth gb"
replace nCH = 1 if gdename == "wezembeek oppem (b)"
replace nCH = 1 if gdename == "parabiago i"
replace nCH = 1 if gdename == "putte nl"
replace nCH = 1 if gdename == "virsbo s"
replace nCH = 1 if gdename == "porto valtravaglia i"
replace nCH = 1 if gdename == "neusaess d"
replace nCH = 1 if gdename == "hilpoltstein d"
replace nCH = 1 if gdename == "mainz kostheim d"
replace nCH = 1 if gdename == "hovas s"
replace nCH = 1 if gdename == "barsbuettel d"
replace nCH = 1 if gdename == "tauberbischofsheim d"
replace nCH = 1 if gdename == "windsor (gb)"
replace nCH = 1 if gdename == "randers dk"
replace nCH = 1 if gdename == "manila (philippinen)"
replace nCH = 1 if gdename == "ruesselsheim (d)"
replace nCH = 1 if gdename == "hoehr grenzhausen"
replace nCH = 1 if gdename == "achim d"
replace nCH = 1 if gdename == "little neck, n.y. usa"
replace nCH = 1 if gdename == "oltingue f"
replace nCH = 1 if gdename == "kowloon hk"
replace nCH = 1 if gdename == "oss (nl)"
replace nCH = 1 if gdename == "saint ismier f"
replace nCH = 1 if gdename == "mettlach d"
replace nCH = 1 if gdename == "moscou"
replace nCH = 1 if gdename == "hohentengen herdern d"
replace nCH = 1 if gdename == "lasnes chapelle st lambert b"
replace nCH = 1 if gdename == "greven d"
replace nCH = 1 if gdename == "viken s"
replace nCH = 1 if gdename == "offenbach/main d"
replace nCH = 1 if gdename == "hoellviksnaes (s)"
replace nCH = 1 if gdename == "biberach brd"
replace nCH = 1 if gdename == "lengo d"
replace nCH = 1 if gdename == "joao pessoa br"
replace nCH = 1 if gdename == "s. donato i"
replace nCH = 1 if gdename == "long grove, illinois usa"
replace nCH = 1 if gdename == "boston ma usa"
replace nCH = 1 if gdename == "warstein (d)"
replace nCH = 1 if gdename == "hellerup (dk)"
replace nCH = 1 if gdename == "ile de man gb"
replace nCH = 1 if gdename == "milton keynes gb"
replace nCH = 1 if gdename == "sarcelles f"
replace nCH = 1 if gdename == "grossgmain a"
replace nCH = 1 if gdename == "villa guardia i"
replace nCH = 1 if gdename == "almunecar e"
replace nCH = 1 if gdename == "nieul f"
replace nCH = 1 if gdename == "aichach d"
replace nCH = 1 if gdename == "lutterbach f"
replace nCH = 1 if gdename == "san diego (usa)"
replace nCH = 1 if gdename == "lossburg d"
replace nCH = 1 if gdename == "ranst b"
replace nCH = 1 if gdename == "istambul (tr)"
replace nCH = 1 if gdename == "kanagawa japan"
replace nCH = 1 if gdename == "badaro rl"
replace nCH = 1 if gdename == "nahariya israel"
replace nCH = 1 if gdename == "ghiffa i"
replace nCH = 1 if gdename == "toyonaka j"
replace nCH = 1 if gdename == "miami (usa)"
replace nCH = 1 if gdename == "north hampton usa"
replace nCH = 1 if gdename == "rheinstetten d"
replace nCH = 1 if gdename == "troisvierges l"
replace nCH = 1 if gdename == "fukuoka j"
replace nCH = 1 if gdename == "torrington, conn. usa"
replace nCH = 1 if gdename == "ciudad py"
replace nCH = 1 if gdename == "rye, new jersey usa"
replace nCH = 1 if gdename == "ickenham gb"
replace nCH = 1 if gdename == "hainburg d"
replace nCH = 1 if gdename == "reutte (a)"
replace nCH = 1 if gdename == "ratingen lintorf d"
replace nCH = 1 if gdename == "oslo/n"
replace nCH = 1 if gdename == "agrate brianza i"
replace nCH = 1 if gdename == "punta paitilla pa"
replace nCH = 1 if gdename == "himeji city j"
replace nCH = 1 if gdename == "neuenbuerg d"
replace nCH = 1 if gdename == "woking surrey gb"
replace nCH = 1 if gdename == "singapur sg"
replace nCH = 1 if gdename == "la rochelle f"
replace nCH = 1 if gdename == "pennington, new jersey usa"
replace nCH = 1 if gdename == "donnbronn d"
replace nCH = 1 if gdename == "huntington usa"
replace nCH = 1 if gdename == "erftstadt d"
replace nCH = 1 if gdename == "poquott usa"
replace nCH = 1 if gdename == "partille s"
replace nCH = 1 if gdename == "criquetot l esneval f"
replace nCH = 1 if gdename == "mersch l"
replace nCH = 1 if gdename == "bromley, kent gb"
replace nCH = 1 if gdename == "wedel d"
replace nCH = 1 if gdename == "bethersden ashford gb"
replace nCH = 1 if gdename == "sundern melschede d"
replace nCH = 1 if gdename == "weybridge, surrey gb"
replace nCH = 1 if gdename == "rottweil d"
replace nCH = 1 if gdename == "denville n.y. usa"
replace nCH = 1 if gdename == "san mateo usa"
replace nCH = 1 if gdename == "muscatine iowa/usa"
replace nCH = 1 if gdename == "howald l"
replace nCH = 1 if gdename == "locmariquer f"
replace nCH = 1 if gdename == "kennelbach a"
replace nCH = 1 if gdename == "perth australien"
replace nCH = 1 if gdename == "leende nl"
replace nCH = 1 if gdename == "wentorf b. hamburg d"
replace nCH = 1 if gdename == "calonge/gerona e"
replace nCH = 1 if gdename == "benajarafe/malaga e"
replace nCH = 1 if gdename == "matadepera e"
replace nCH = 1 if gdename == "sora i"
replace nCH = 1 if gdename == "istambul tr"
replace nCH = 1 if gdename == "tampere (sf)"
replace nCH = 1 if gdename == "beaufays b"
replace nCH = 1 if gdename == "madison, wisconsin usa"
replace nCH = 1 if gdename == "bettlach f"
replace nCH = 1 if gdename == "bloomfield hills usa"
replace nCH = 1 if gdename == "hessisch oldendorf d"
replace nCH = 1 if gdename == "cagliari i"
replace nCH = 1 if gdename == "maison laffitte f"
replace nCH = 1 if gdename == "bad muender d"
replace nCH = 1 if gdename == "st. peter port (guernsey)"
replace nCH = 1 if gdename == "kiel (d)"
replace nCH = 1 if gdename == "oisterwijk (nl)"
replace nCH = 1 if gdename == "senningerberg l"
replace nCH = 1 if gdename == "new canaan (usa)"
replace nCH = 1 if gdename == "hartley witney gb"
replace nCH = 1 if gdename == "new caanan usa"
replace nCH = 1 if gdename == "essen b"
replace nCH = 1 if gdename == "morbier f"
replace nCH = 1 if gdename == "ludwigsburg (d)"
replace nCH = 1 if gdename == "pau f"
replace nCH = 1 if gdename == "glendale/californie, usa"
replace nCH = 1 if gdename == "neunkirchen/saar d"
replace nCH = 1 if gdename == "eschborn (d)"
replace nCH = 1 if gdename == "wiesensteig d"
replace nCH = 1 if gdename == "remagen d"
replace nCH = 1 if gdename == "feucherolles f"
replace nCH = 1 if gdename == "ponte di brenta i"
replace nCH = 1 if gdename == "salvador, bahia br"
replace nCH = 1 if gdename == "granges so"
replace nCH = 1 if gdename == "whitinsville usa"
replace nCH = 1 if gdename == "wiesloch d"
replace nCH = 1 if gdename == "newark, new jersey usa"
replace nCH = 1 if gdename == "orba alicante e"
replace nCH = 1 if gdename == "saint etienne f"
replace nCH = 1 if gdename == "montreal ca"
replace nCH = 1 if gdename == "crainheim b"
replace nCH = 1 if gdename == "vaexjoe (s)"
replace nCH = 1 if gdename == "binzen b. loerrach d"
replace nCH = 1 if gdename == "stuttgard (d)"
replace nCH = 1 if gdename == "fayet f"
replace nCH = 1 if gdename == "roubaix f"
replace nCH = 1 if gdename == "assmannshardt d"
replace nCH = 1 if gdename == "domburg nl"
replace nCH = 1 if gdename == "seefeld (d)"
replace nCH = 1 if gdename == "broumana rl"
replace nCH = 1 if gdename == "chartronges f"
replace nCH = 1 if gdename == "dinslaken d"
replace nCH = 1 if gdename == "los angeles/usa"
replace nCH = 1 if gdename == "washington d.c. (usa)"
replace nCH = 1 if gdename == "hettenhausen d"
replace nCH = 1 if gdename == "ridgefield (connecticut usa)"
replace nCH = 1 if gdename == "trezzano sul naviglio, milano"
replace nCH = 1 if gdename == "schenefeld d"
replace nCH = 1 if gdename == "ebertsheim (d)"
replace nCH = 1 if gdename == "schwenningen (d)"
replace nCH = 1 if gdename == "beckenham, kent gb"
replace nCH = 1 if gdename == "ridgefield connecticut usa"
replace nCH = 1 if gdename == "wirral gb"
replace nCH = 1 if gdename == "chenex f"
replace nCH = 1 if gdename == "assesse b"
replace nCH = 1 if gdename == "toceno i"
replace nCH = 1 if gdename == "hanau hohe tanne d"
replace nCH = 1 if gdename == "van nuys usa"
replace nCH = 1 if gdename == "poecking d"
replace nCH = 1 if gdename == "lanaken b"
replace nCH = 1 if gdename == "la gaubretiere f"
replace nCH = 1 if gdename == "haemeenlinna sf"
replace nCH = 1 if gdename == "bad hersfeld d"
replace nCH = 1 if gdename == "saint canut cdn"
replace nCH = 1 if gdename == "lima/peru"
replace nCH = 1 if gdename == "tel aviv/israel"
replace nCH = 1 if gdename == "opio f"
replace nCH = 1 if gdename == "great houghton, northampton gb"
replace nCH = 1 if gdename == "biberach (d)"
replace nCH = 1 if gdename == "albstadt ebingen d"
replace nCH = 1 if gdename == "siegburg (d)"
replace nCH = 1 if gdename == "brunn a"
replace nCH = 1 if gdename == "zwingenberg d"
replace nCH = 1 if gdename == "togliatti (rfsfr)"
replace nCH = 1 if gdename == "kensington usa"
replace nCH = 1 if gdename == "bozen i"
replace nCH = 1 if gdename == "dachau (d)"
replace nCH = 1 if gdename == "alleur b"
replace nCH = 1 if gdename == "seevetal d"
replace nCH = 1 if gdename == "cerfontaine b"
replace nCH = 1 if gdename == "bueren a.d.aare"
replace nCH = 1 if gdename == "muehlacker (d)"
replace nCH = 1 if gdename == "hojbjerg dk"
replace nCH = 1 if gdename == "vaerlose dk"
replace nCH = 1 if gdename == "johannesburg, suedafrika"
replace nCH = 1 if gdename == "vibeke d"
replace nCH = 1 if gdename == "dublino irl"
replace nCH = 1 if gdename == "julianstown, co. meath irl"
replace nCH = 1 if gdename == "gelnhausen d"
replace nCH = 1 if gdename == "weingarten (d)"
replace nCH = 1 if gdename == "swindon gb"
replace nCH = 1 if gdename == "thoune"
replace nCH = 1 if gdename == "nofels a"
replace nCH = 1 if gdename == "palermo i"
replace nCH = 1 if gdename == "firenze (i)"
replace nCH = 1 if gdename == "karlsruhe brd"
replace nCH = 1 if gdename == "la havana c"
replace nCH = 1 if gdename == "franklin park, illinois usa"
replace nCH = 1 if gdename == "lauchringen (d)"
replace nCH = 1 if gdename == "st. augustin d"
replace nCH = 1 if gdename == "weil schoenbuch d"
replace nCH = 1 if gdename == "grankulla (sf)"
replace nCH = 1 if gdename == "schuettorf d"
replace nCH = 1 if gdename == "ischgel a"
replace nCH = 1 if gdename == "deutschland (d)"
replace nCH = 1 if gdename == "umkirch d"
replace nCH = 1 if gdename == "el barraco e"
replace nCH = 1 if gdename == "bad herrenalb d"
replace nCH = 1 if gdename == "kinnelon usa"
replace nCH = 1 if gdename == "bisamberg a"
replace nCH = 1 if gdename == "brighton gb"
replace nCH = 1 if gdename == "kressbronn d"
replace nCH = 1 if gdename == "tranby n"
replace nCH = 1 if gdename == "leverkusen (d)"
replace nCH = 1 if gdename == "hoesbach d"
replace nCH = 1 if gdename == "vancouver (kanada)"
replace nCH = 1 if gdename == "voissins le bretonneux f"
replace nCH = 1 if gdename == "crondall gb"
replace nCH = 1 if gdename == "sanaa ymn"
replace nCH = 1 if gdename == "granvenhage nl"
replace nCH = 1 if gdename == "kufstein a"
replace nCH = 1 if gdename == "seewalchen a"
replace nCH = 1 if gdename == "ridderkerk slikkerveer n"
replace nCH = 1 if gdename == "chambery (f)"
replace nCH = 1 if gdename == "balingen (d)"
replace nCH = 1 if gdename == "javia e"
replace nCH = 1 if gdename == "djeddah as"
replace nCH = 1 if gdename == "osterfildern, ruit d"
replace nCH = 1 if gdename == "gand (b)"
replace nCH = 1 if gdename == "reiskirchen d"
replace nCH = 1 if gdename == "maassluis nl"
replace nCH = 1 if gdename == "mclean va usa"
replace nCH = 1 if gdename == "peronnas f"
replace nCH = 1 if gdename == "saint jean de tholome f"
replace nCH = 1 if gdename == "irvine ca usa"
replace nCH = 1 if gdename == "sesto s. giovanni i"
replace nCH = 1 if gdename == "pino torinese (i)"
replace nCH = 1 if gdename == "sark, channel island gb"
replace nCH = 1 if gdename == "dorking, surrey gb"
replace nCH = 1 if gdename == "brescia (i)"
replace nCH = 1 if gdename == "sint martenslennik b"
replace nCH = 1 if gdename == "isla de la palma, canarias e"
replace nCH = 1 if gdename == "meina, arona i"
replace nCH = 1 if gdename == "copponex f"
replace nCH = 1 if gdename == "brookvale n.s.w. aus"
replace nCH = 1 if gdename == "gerlingen bei stuttgart"
replace nCH = 1 if gdename == "bryanston rsa"
replace nCH = 1 if gdename == "santos (brasilien)"
replace nCH = 1 if gdename == "brunoy f"
replace nCH = 1 if gdename == "vaxholm s"
replace nCH = 1 if gdename == "boulogne billancourt (f)"
replace nCH = 1 if gdename == "st louis (f)"
replace nCH = 1 if gdename == "lindenberg im allgaeu d"
replace nCH = 1 if gdename == "ciquio i"
replace nCH = 1 if gdename == "encino usa"
replace nCH = 1 if gdename == "portoroz yu"
replace nCH = 1 if gdename == "toenisvorst d"
replace nCH = 1 if gdename == "willow grove, pennsylvania usa"
replace nCH = 1 if gdename == "tulsa usa"
replace nCH = 1 if gdename == "heverlee (b)"
replace nCH = 1 if gdename == "kerpen d"
replace nCH = 1 if gdename == "erfstadt d"
replace nCH = 1 if gdename == "porto cervo i"
replace nCH = 1 if gdename == "carrollton/georgia usa"
replace nCH = 1 if gdename == "portland usa"
replace nCH = 1 if gdename == "borgholzhausen d"
replace nCH = 1 if gdename == "gerdlstein d"
replace nCH = 1 if gdename == "ludwigshafen a. rh. d"
replace nCH = 1 if gdename == "selfkant d"
replace nCH = 1 if gdename == "mortsel b"
replace nCH = 1 if gdename == "hemel hempstead gb"
replace nCH = 1 if gdename == "meerhout b"
replace nCH = 1 if gdename == "hengersberg d"
replace nCH = 1 if gdename == "calw hirsau d"
replace nCH = 1 if gdename == "griesheim d"
replace nCH = 1 if gdename == "induno olona i"
replace nCH = 1 if gdename == "bobaeck sf"
replace nCH = 1 if gdename == "grenzach wyhlen (d)"
replace nCH = 1 if gdename == "leawood usa"
replace nCH = 1 if gdename == "cavallasca i"
replace nCH = 1 if gdename == "flemalle haute b"
replace nCH = 1 if gdename == "wermelskirchen (d)"
replace nCH = 1 if gdename == "blumberg d"
replace nCH = 1 if gdename == "st antoine f"
replace nCH = 1 if gdename == "kulmbach, bayern d"
replace nCH = 1 if gdename == "la celle saint cloud (f)"
replace nCH = 1 if gdename == "belleme f"
replace nCH = 1 if gdename == "weilstetten d"
replace nCH = 1 if gdename == "schopfheim brd"
replace nCH = 1 if gdename == "pittsburgh/usa"
replace nCH = 1 if gdename == "bad voelsau a"
replace nCH = 1 if gdename == "neuwied d"
replace nCH = 1 if gdename == "whitchurch on thames gb"
replace nCH = 1 if gdename == "rodenheim d"
replace nCH = 1 if gdename == "wupperthal d"
replace nCH = 1 if gdename == "focsani ro"
replace nCH = 1 if gdename == "esbo fin"
replace nCH = 1 if gdename == "santa monica (usa)"
replace nCH = 1 if gdename == "parigi"
replace nCH = 1 if gdename == "vougeot, cote d or f"
replace nCH = 1 if gdename == "ekeren b"
replace nCH = 1 if gdename == "hohentengen a.h. d"
replace nCH = 1 if gdename == "roquebrune cap martin f"
replace nCH = 1 if gdename == "oulu sf"
replace nCH = 1 if gdename == "fullerton, california usa"
replace nCH = 1 if gdename == "saltsjoe boo s"
replace nCH = 1 if gdename == "cernusco sul naviglio i"
replace nCH = 1 if gdename == "naarden"
replace nCH = 1 if gdename == "muenchberg/oberfranken d"
replace nCH = 1 if gdename == "toronto cnd"
replace nCH = 1 if gdename == "godalming gb"
replace nCH = 1 if gdename == "driebergen (nl)"
replace nCH = 1 if gdename == "evere b"
replace nCH = 1 if gdename == "la spezia i"
replace nCH = 1 if gdename == "uffheim/sierentz f"
replace nCH = 1 if gdename == "piacenza (i)"
replace nCH = 1 if gdename == "klippan s"
replace nCH = 1 if gdename == "condamine f"
replace nCH = 1 if gdename == "ebenhausenzell d"
replace nCH = 1 if gdename == "grimberger b"
replace nCH = 1 if gdename == "ferno i"
replace nCH = 1 if gdename == "barcelona es"
replace nCH = 1 if gdename == "dueren guerzenich d"
replace nCH = 1 if gdename == "segrate milano i"
replace nCH = 1 if gdename == "bibbona, livorno i"
replace nCH = 1 if gdename == "woerthsee d"
replace nCH = 1 if gdename == "loenen aan de vecht nl"
replace nCH = 1 if gdename == "st. clements, ontario cdn"
replace nCH = 1 if gdename == "st. louis, missouri usa"
replace nCH = 1 if gdename == "st. petersburg rus"
replace nCH = 1 if gdename == "wuppertal barmen d"
replace nCH = 1 if gdename == "hausen ag"
replace nCH = 1 if gdename == "san jose calif./usa"
replace nCH = 1 if gdename == "puettlingen d"
replace nCH = 1 if gdename == "hoelsingborg"
replace nCH = 1 if gdename == "wegscheid f"
replace nCH = 1 if gdename == "bour l"
replace nCH = 1 if gdename == "new york (ny. usa)"
replace nCH = 1 if gdename == "helvoirt nl"
replace nCH = 1 if gdename == "hemsbach d"
replace nCH = 1 if gdename == "ramat hasharon il"
replace nCH = 1 if gdename == "granger, indiana usa"
replace nCH = 1 if gdename == "marlow gb"
replace nCH = 1 if gdename == "tripoli (lt)"
replace nCH = 1 if gdename == "bagno a ripoli i"
replace nCH = 1 if gdename == "saint jean de gonville f"
replace nCH = 1 if gdename == "schliern/koeniz"
replace nCH = 1 if gdename == "gaggenau/ selbach d"
replace nCH = 1 if gdename == "sistrans a"
replace nCH = 1 if gdename == "armonk usa"
replace nCH = 1 if gdename == "yarze rl"
replace nCH = 1 if gdename == "eccles gb"
replace nCH = 1 if gdename == "friedrichshafen brd"
replace nCH = 1 if gdename == "lysaker n"
replace nCH = 1 if gdename == "lenzerheide/vaz/obervaz"
replace nCH = 1 if gdename == "liverpool gb"
replace nCH = 1 if gdename == "hemmingen d"
replace nCH = 1 if gdename == "taipei (taiwan)"
replace nCH = 1 if gdename == "dassendorf d"
replace nCH = 1 if gdename == "new canaan ct usa"
replace nCH = 1 if gdename == "tailfingen d"
replace nCH = 1 if gdename == "schoorl nl"
replace nCH = 1 if gdename == "senden ay d"
replace nCH = 1 if gdename == "manila (pi)"
replace nCH = 1 if gdename == "roethis a"
replace nCH = 1 if gdename == "renchen d"
replace nCH = 1 if gdename == "chaumont gistoux b"
replace nCH = 1 if gdename == "theheran ir"
replace nCH = 1 if gdename == "homborsund n"
replace nCH = 1 if gdename == "castenaso i"
replace nCH = 1 if gdename == "manchline/ayrshire gb"
replace nCH = 1 if gdename == "alzey d"
replace nCH = 1 if gdename == "route de la"
replace nCH = 1 if gdename == "troy usa"
replace nCH = 1 if gdename == "hillside usa"
replace nCH = 1 if gdename == "kokomo, indiana usa"
replace nCH = 1 if gdename == "siegsdorf d"
replace nCH = 1 if gdename == "stellenbosch za"
replace nCH = 1 if gdename == "albavilla i"
replace nCH = 1 if gdename == "todtnau d"
replace nCH = 1 if gdename == "kauniainen (sf)"
replace nCH = 1 if gdename == "drammen n"
replace nCH = 1 if gdename == "altadena usa"
replace nCH = 1 if gdename == "knowle, solihull gb"
replace nCH = 1 if gdename == "auckland (neuseeland)"
replace nCH = 1 if gdename == "heiligenhaus (d)"
replace nCH = 1 if gdename == "waldbronn brd"
replace nCH = 1 if gdename == "eddelsen/hittfeld d"
replace nCH = 1 if gdename == "weed nl"
replace nCH = 1 if gdename == "belfast irl"
replace nCH = 1 if gdename == "saint souplet f"
replace nCH = 1 if gdename == "mittelhausbergen f"
replace nCH = 1 if gdename == "tivoli i"
replace nCH = 1 if gdename == "marostica, vicenza i"
replace nCH = 1 if gdename == "isolella, corse f"
replace nCH = 1 if gdename == "st gingolph f"
replace nCH = 1 if gdename == "basking ridge, new jersey usa"
replace nCH = 1 if gdename == "bethesda usa"
replace nCH = 1 if gdename == "yokohama (j)"
replace nCH = 1 if gdename == "youghal irl"
replace nCH = 1 if gdename == "camberley (gb)"
replace nCH = 1 if gdename == "moraira e"
replace nCH = 1 if gdename == "spring, texas usa"
replace nCH = 1 if gdename == "goerwihl d"
replace nCH = 1 if gdename == "eppelhorn wiesbach d"
replace nCH = 1 if gdename == "southampton gb"
replace nCH = 1 if gdename == "frossasco i"
replace nCH = 1 if gdename == "vetraz monthoux (f)"
replace nCH = 1 if gdename == "windsor berkshire gb"
replace nCH = 1 if gdename == "cheadle, cheshire gb"
replace nCH = 1 if gdename == "st. helier, jersey gb"
replace nCH = 1 if gdename == "lonsle saunier f"
replace nCH = 1 if gdename == "eschwege d"
replace nCH = 1 if gdename == "freiberg/neckar d"
replace nCH = 1 if gdename == "great missenden gb"
replace nCH = 1 if gdename == "la louviere b"
replace nCH = 1 if gdename == "sao paulo (brasilien)"
replace nCH = 1 if gdename == "nelspruit za"
replace nCH = 1 if gdename == "elsmark dk"
replace nCH = 1 if gdename == "woodbury usa"
replace nCH = 1 if gdename == "sesto calende i"
replace nCH = 1 if gdename == "hartenfels (d)"
replace nCH = 1 if gdename == "hoofddorp nl"
replace nCH = 1 if gdename == "goefis a"
replace nCH = 1 if gdename == "mooloolah/australien"
replace nCH = 1 if gdename == "highland park usa"
replace nCH = 1 if gdename == "assling d"
replace nCH = 1 if gdename == "rodgau d"
replace nCH = 1 if gdename == "granges"
replace nCH = 1 if gdename == "great neck usa"
replace nCH = 1 if gdename == "hagenthal le haut f"
replace nCH = 1 if gdename == "leuven b"
replace nCH = 1 if gdename == "bergen (n)"
replace nCH = 1 if gdename == "belem br"
replace nCH = 1 if gdename == "kiryat ono il"
replace nCH = 1 if gdename == "sheffield gb"
replace nCH = 1 if gdename == "bickenbach d"
replace nCH = 1 if gdename == "guipavas f"
replace nCH = 1 if gdename == "cuckfield gb"
replace nCH = 1 if gdename == "elm grove, winsconsin usa"
replace nCH = 1 if gdename == "divonne les bains (f)"
replace nCH = 1 if gdename == "lampertheim f"
replace nCH = 1 if gdename == "nesbru n"
replace nCH = 1 if gdename == "kettwig/d"
replace nCH = 1 if gdename == "nordwick nl"
replace nCH = 1 if gdename == "essen brd"
replace nCH = 1 if gdename == "bramhall gb"
replace nCH = 1 if gdename == "cleveland, ohio usa"
replace nCH = 1 if gdename == "salerno i"
replace nCH = 1 if gdename == "glen ridge, new jersey usa"
replace nCH = 1 if gdename == "altena d"
replace nCH = 1 if gdename == "arbere f"
replace nCH = 1 if gdename == "carugate i"
replace nCH = 1 if gdename == "raincy f"
replace nCH = 1 if gdename == "abu dhabi sv"
replace nCH = 1 if gdename == "malmoe angdalavaegen s"
replace nCH = 1 if gdename == "samoens f"
replace nCH = 1 if gdename == "bourges b"
replace nCH = 1 if gdename == "kelsterbach d"
replace nCH = 1 if gdename == "osteras n"
replace nCH = 1 if gdename == "balibago rp"
replace nCH = 1 if gdename == "gronau d"
replace nCH = 1 if gdename == "montcherand vd"
replace nCH = 1 if gdename == "iver heath gb"
replace nCH = 1 if gdename == "alton hampshire gb"
replace nCH = 1 if gdename == "contamine sarzin f"
replace nCH = 1 if gdename == "nauders a"
replace nCH = 1 if gdename == "radlett hertfordshire gb"
replace nCH = 1 if gdename == "busto arsizio (i)"
replace nCH = 1 if gdename == "cosne sur loire, nievre f"
replace nCH = 1 if gdename == "saint cergues f"
replace nCH = 1 if gdename == "halifax, west yorkshire gb"
replace nCH = 1 if gdename == "montreal can"
replace nCH = 1 if gdename == "luxenburg"
replace nCH = 1 if gdename == "nutley, nj usa"
replace nCH = 1 if gdename == "aldridge gb"
replace nCH = 1 if gdename == "esslingen (d)"
replace nCH = 1 if gdename == "katlenburg lindau"
replace nCH = 1 if gdename == "ugine f"
replace nCH = 1 if gdename == "lingolsheim f"
replace nCH = 1 if gdename == "mon desir (sark, channel"
replace nCH = 1 if gdename == "lauingen an der donau d"
replace nCH = 1 if gdename == "hockenheim (d)"
replace nCH = 1 if gdename == "djedda sa"
replace nCH = 1 if gdename == "nagoya j"
replace nCH = 1 if gdename == "furtwangen d"
replace nCH = 1 if gdename == "evillers f"
replace nCH = 1 if gdename == "pluederhausen d"
replace nCH = 1 if gdename == "savyon il"
replace nCH = 1 if gdename == "rueckersdorf d"
replace nCH = 1 if gdename == "odense dk"
replace nCH = 1 if gdename == "sulz a"
replace nCH = 1 if gdename == "odenthal d"
replace nCH = 1 if gdename == "rochester (usa)"
replace nCH = 1 if gdename == "kingswood gb"
replace nCH = 1 if gdename == "englewood, new jersey usa"
replace nCH = 1 if gdename == "dar es salaam eaz"
replace nCH = 1 if gdename == "sennestadt (d)"
replace nCH = 1 if gdename == "macau"
replace nCH = 1 if gdename == "stadtbergen d"
replace nCH = 1 if gdename == "lohmar (d)"
replace nCH = 1 if gdename == "sainte adresse f"
replace nCH = 1 if gdename == "mohadessin giza et"
replace nCH = 1 if gdename == "la celle st cloud(f)"
replace nCH = 1 if gdename == "askim s"
replace nCH = 1 if gdename == "auckland nz"
replace nCH = 1 if gdename == "obersulmetigen d"
replace nCH = 1 if gdename == "charlottesville, virginia (usa)"
replace nCH = 1 if gdename == "inning d"
replace nCH = 1 if gdename == "bergen nl"
replace nCH = 1 if gdename == "vegsund n"
replace nCH = 1 if gdename == "santa ponsa/iles baleares e"
replace nCH = 1 if gdename == "geel b"
replace nCH = 1 if gdename == "bruessel/b"
replace nCH = 1 if gdename == "richmond (virginia) usa"
replace nCH = 1 if gdename == "chens le pont f"
replace nCH = 1 if gdename == "nicosia cy"
replace nCH = 1 if gdename == "serrieres de briord f"
replace nCH = 1 if gdename == "costamasnaga i"
replace nCH = 1 if gdename == "schoeneck d"
replace nCH = 1 if gdename == "villmar d"
replace nCH = 1 if gdename == "gilching d"
replace nCH = 1 if gdename == "aspen/colorado (usa)"
replace nCH = 1 if gdename == "sevenoaks (gb)"
replace nCH = 1 if gdename == "meinerzhagen d"
replace nCH = 1 if gdename == "neu isenburg"
replace nCH = 1 if gdename == "untersteinbach d"
replace nCH = 1 if gdename == "san mateo, california usa"
replace nCH = 1 if gdename == "gand b"
replace nCH = 1 if gdename == "weilheim ay d"
replace nCH = 1 if gdename == "rockville, maryland usa"
replace nCH = 1 if gdename == "singapour"
replace nCH = 1 if gdename == "darmstadt/d"
replace nCH = 1 if gdename == "schwarzach a"
replace nCH = 1 if gdename == "schwertberg a"
replace nCH = 1 if gdename == "toronto/ kanada"
replace nCH = 1 if gdename == "eppstein d"
replace nCH = 1 if gdename == "mbabane (swaziland)"
replace nCH = 1 if gdename == "n djamena, tschad"
replace nCH = 1 if gdename == "culford gb"
replace nCH = 1 if gdename == "ismaning (d)"
replace nCH = 1 if gdename == "wassenhaar nl"
replace nCH = 1 if gdename == "gundelsheim d"
replace nCH = 1 if gdename == "greenville (usa)"
replace nCH = 1 if gdename == "pesmes f"
replace nCH = 1 if gdename == "noisy le roy f"
replace nCH = 1 if gdename == "akron, ohio (usa)"
replace nCH = 1 if gdename == "sassenberg (d)"
replace nCH = 1 if gdename == "turino i"
replace nCH = 1 if gdename == "fie/bolzano i"
replace nCH = 1 if gdename == "ham sur heure b"
replace nCH = 1 if gdename == "wiehl (d)"
replace nCH = 1 if gdename == "north vancouver can"
replace nCH = 1 if gdename == "akron, ohio usa"
replace nCH = 1 if gdename == "sandy hook, conn usa"
replace nCH = 1 if gdename == "coburg (d)"
replace nCH = 1 if gdename == "laarne b"
replace nCH = 1 if gdename == "riihimaeki sf"
replace nCH = 1 if gdename == "kuenzelsau gaisbach d"
replace nCH = 1 if gdename == "rheinbach (d)"
replace nCH = 1 if gdename == "hallein (a)"
replace nCH = 1 if gdename == "telgte d"
replace nCH = 1 if gdename == "burbank usa"
replace nCH = 1 if gdename == "moskau ussr"
replace nCH = 1 if gdename == "seven hills, ohio usa"
replace nCH = 1 if gdename == "frankfurt/main (d)"
replace nCH = 1 if gdename == "rose/avry sur matran"
replace nCH = 1 if gdename == "bronnen d"
replace nCH = 1 if gdename == "villingen schwenningen"
replace nCH = 1 if gdename == "schortens d"
replace nCH = 1 if gdename == "noordwjk nl"
replace nCH = 1 if gdename == "campodarsego i"
replace nCH = 1 if gdename == "oirschot nl"
replace nCH = 1 if gdename == "westmount can"
replace nCH = 1 if gdename == "dallas (texas. usa)"
replace nCH = 1 if gdename == "dallas texas usa"
replace nCH = 1 if gdename == "tulsa (usa)"
replace nCH = 1 if gdename == "costa rica"
replace nCH = 1 if gdename == "garches (f)"
replace nCH = 1 if gdename == "hartfortshire gb"
replace nCH = 1 if gdename == "kalamazoo usa"
replace nCH = 1 if gdename == "luedenscheid (d)"
replace nCH = 1 if gdename == "damman (sa)"
replace nCH = 1 if gdename == "ville d avray, haute seine f"
replace nCH = 1 if gdename == "porto colom, mallorca e"
replace nCH = 1 if gdename == "dourdan f"
replace nCH = 1 if gdename == "bratislava csfr"
replace nCH = 1 if gdename == "vence (f)"
replace nCH = 1 if gdename == "baden a"
replace nCH = 1 if gdename == "otterp dk"
replace nCH = 1 if gdename == "morez f"
replace nCH = 1 if gdename == "maidstone gb"
replace nCH = 1 if gdename == "papendrecht nl"
replace nCH = 1 if gdename == "doha, quatar (rl)"
replace nCH = 1 if gdename == "rexdale"
replace nCH = 1 if gdename == "starnberg (d)"
replace nCH = 1 if gdename == "st andre nord f"
replace nCH = 1 if gdename == "thaur a"
replace nCH = 1 if gdename == "pocklington, yorkshire gb"
replace nCH = 1 if gdename == "bevel, nijlen b"
replace nCH = 1 if gdename == "villanova pa, usa"
replace nCH = 1 if gdename == "pougny f"
replace nCH = 1 if gdename == "verso f"
replace nCH = 1 if gdename == "peck f"
replace nCH = 1 if gdename == "joppe nl"
replace nCH = 1 if gdename == "chieti i"
replace nCH = 1 if gdename == "castino i"
replace nCH = 1 if gdename == "arhus dk"
replace nCH = 1 if gdename == "buxtehude d"
replace nCH = 1 if gdename == "maarn nl"
replace nCH = 1 if gdename == "graested dk"
replace nCH = 1 if gdename == "neuberg d"
replace nCH = 1 if gdename == "weilheim (d)"
replace nCH = 1 if gdename == "netmold d"
replace nCH = 1 if gdename == "hoquiam, wash. usa"
replace nCH = 1 if gdename == "loinsin f"
replace nCH = 1 if gdename == "garbagnate milanese i"
replace nCH = 1 if gdename == "winchester usa"
replace nCH = 1 if gdename == "lucinges f"
replace nCH = 1 if gdename == "carvalhos p"
replace nCH = 1 if gdename == "koblenz (d)"
replace nCH = 1 if gdename == "scarsdale, n.y. usa"
replace nCH = 1 if gdename == "pomezia i"
replace nCH = 1 if gdename == "brackel d"
replace nCH = 1 if gdename == "unter schleissheim d"
replace nCH = 1 if gdename == "halle d"
replace nCH = 1 if gdename == "macau (macao)"
replace nCH = 1 if gdename == "duelmen (d)"
replace nCH = 1 if gdename == "franklin lakes (new jersey. usa)"
replace nCH = 1 if gdename == "clairefontaine f"
replace nCH = 1 if gdename == "tampa usa"
replace nCH = 1 if gdename == "goetzis (a)"
replace nCH = 1 if gdename == "wertheim (d)"
replace nCH = 1 if gdename == "andresy (f)"
replace nCH = 1 if gdename == "mulfingen d"
replace nCH = 1 if gdename == "tasisn la dmei lune f"
replace nCH = 1 if gdename == "concord ma usa"
replace nCH = 1 if gdename == "wintzenheim f"
replace nCH = 1 if gdename == "elze d"
replace nCH = 1 if gdename == "erlangen (d)"
replace nCH = 1 if gdename == "stroud gb"
replace nCH = 1 if gdename == "aumuehle brd"
replace nCH = 1 if gdename == "sunndalsoera n"
replace nCH = 1 if gdename == "st. louis (usa)"
replace nCH = 1 if gdename == "boskoop nl"
replace nCH = 1 if gdename == "grossburgwedel d"
replace nCH = 1 if gdename == "laucherthal b.s. d"
replace nCH = 1 if gdename == "apeldorn nl"
replace nCH = 1 if gdename == "schmitten d"
replace nCH = 1 if gdename == "ontario (cdn)"
replace nCH = 1 if gdename == "st. martens latem b"
replace nCH = 1 if gdename == "neupre b"
replace nCH = 1 if gdename == "outremont cdn"
replace nCH = 1 if gdename == "la moraleja e"
replace nCH = 1 if gdename == "conakry (gb)"
replace nCH = 1 if gdename == "baton rouge usa"
replace nCH = 1 if gdename == "perugia i"
replace nCH = 1 if gdename == "douglas, isle of man gb"
replace nCH = 1 if gdename == "keltern ellmendingen d"
replace nCH = 1 if gdename == "wartenberg d"
replace nCH = 1 if gdename == "woerleschwang d"
replace nCH = 1 if gdename == "grigny (f)"
replace nCH = 1 if gdename == "villarbasse i"
replace nCH = 1 if gdename == "headley gb"
replace nCH = 1 if gdename == "genes (i)"
replace nCH = 1 if gdename == "bensheim (d)"
replace nCH = 1 if gdename == "guesnou f"
replace nCH = 1 if gdename == "kettwig d"
replace nCH = 1 if gdename == "kaarst (d)"
replace nCH = 1 if gdename == "guernesey gb"
replace nCH = 1 if gdename == "xanten marienbaum d"
replace nCH = 1 if gdename == "lindome s"
replace nCH = 1 if gdename == "grandcamp (f)"
replace nCH = 1 if gdename == "den haag/nl"
replace nCH = 1 if gdename == "birkeroed dk"
replace nCH = 1 if gdename == "livingston n.j. usa"
replace nCH = 1 if gdename == "schaebisch gmuend d"
replace nCH = 1 if gdename == "belmont, mass (usa)"
replace nCH = 1 if gdename == "balzheim d"
replace nCH = 1 if gdename == "tinqueux f"
replace nCH = 1 if gdename == "balingen d"
replace nCH = 1 if gdename == "bad duerkheim (d)"
replace nCH = 1 if gdename == "ingersheim d"
replace nCH = 1 if gdename == "shropshire gb"
replace nCH = 1 if gdename == "kronberg im taunus d"
replace nCH = 1 if gdename == "sark/channel islands gb"
replace nCH = 1 if gdename == "enebyberg s"
replace nCH = 1 if gdename == "rockwall, texas usa"
replace nCH = 1 if gdename == "sutton coldfield, west midlands gb"
replace nCH = 1 if gdename == "clinton, connecticut usa"
replace nCH = 1 if gdename == "wuelfrath (d)"
replace nCH = 1 if gdename == "zoersel b"
replace nCH = 1 if gdename == "scarsdale ny usa"
replace nCH = 1 if gdename == "larchmont, new york usa"
replace nCH = 1 if gdename == "rungsted kyst (daenemark)"
replace nCH = 1 if gdename == "idaho usa"
replace nCH = 1 if gdename == "brentford gb"
replace nCH = 1 if gdename == "tolmezzo/udine i"
replace nCH = 1 if gdename == "cavaion veronese i"
replace nCH = 1 if gdename == "quebec (cdn)"
replace nCH = 1 if gdename == "pittsburgh (pennsylvania. usa)"
replace nCH = 1 if gdename == "battenberg/eder d"
replace nCH = 1 if gdename == "cos cob usa"
replace nCH = 1 if gdename == "saint jorioz f"
replace nCH = 1 if gdename == "dallas, texas us"
replace nCH = 1 if gdename == "ljubliana slo"
replace nCH = 1 if gdename == "hamilton bermudes"
replace nCH = 1 if gdename == "dubai et"
replace nCH = 1 if gdename == "edenvale rsa"
replace nCH = 1 if gdename == "brande dk"
replace nCH = 1 if gdename == "buching d"
replace nCH = 1 if gdename == "spokane, washington usa"
replace nCH = 1 if gdename == "palatine usa"
replace nCH = 1 if gdename == "nairobi (eak)"
replace nCH = 1 if gdename == "woodside usa"
replace nCH = 1 if gdename == "champs sur marne f"
replace nCH = 1 if gdename == "bad laer d"
replace nCH = 1 if gdename == "kreuth d"
replace nCH = 1 if gdename == "peoria usa"
replace nCH = 1 if gdename == "sprimont b"
replace nCH = 1 if gdename == "stresa i"
replace nCH = 1 if gdename == "hof a, haag"
replace nCH = 1 if gdename == "kronberg/taunus d"
replace nCH = 1 if gdename == "punta gorda, fla usa"
replace nCH = 1 if gdename == "eaubonne f"
replace nCH = 1 if gdename == "madison usa"
replace nCH = 1 if gdename == "north salem, new york (usa)"
replace nCH = 1 if gdename == "neunburg vorm wald d"
replace nCH = 1 if gdename == "vaihingen enz d"
replace nCH = 1 if gdename == "shenfield gb"
replace nCH = 1 if gdename == "vetileon"
replace nCH = 1 if gdename == "hofheim (d)"
replace nCH = 1 if gdename == "needham/massachusetts usa"
replace nCH = 1 if gdename == "truro, cornwall gb"
replace nCH = 1 if gdename == "steinheim a. m"
replace nCH = 1 if gdename == "braunfels d"
replace nCH = 1 if gdename == "lorch waldhausen d"
replace nCH = 1 if gdename == "herzlyya il"
replace nCH = 1 if gdename == "geilenkirchen d"
replace nCH = 1 if gdename == "chicago illinois usa"
replace nCH = 1 if gdename == "todtnau (d)"
replace nCH = 1 if gdename == "new port richey"
replace nCH = 1 if gdename == "gueterloh d"
replace nCH = 1 if gdename == "st martin guernsey gb"
replace nCH = 1 if gdename == "palazzolo sul l oglio i"
replace nCH = 1 if gdename == "ostende b"
replace nCH = 1 if gdename == "marina die porto cervo i"
replace nCH = 1 if gdename == "agadir ma"
replace nCH = 1 if gdename == "pinneberg d"
replace nCH = 1 if gdename == "lincolnwood il usa"
replace nCH = 1 if gdename == "klettgau buehl d"
replace nCH = 1 if gdename == "leamington gb"
replace nCH = 1 if gdename == "gstaad/saanen"
replace nCH = 1 if gdename == "dnepropetrovsk uk"
replace nCH = 1 if gdename == "wildpoldsried d"
replace nCH = 1 if gdename == "niedernhausen d"
replace nCH = 1 if gdename == "dorking gb"
replace nCH = 1 if gdename == "union mills (isle of man)"
replace nCH = 1 if gdename == "surbiton (gb)"
replace nCH = 1 if gdename == "singapore, sgp"
replace nCH = 1 if gdename == "hambruecken d"
replace nCH = 1 if gdename == "monza, i"
replace nCH = 1 if gdename == "ratingen hoesel (d)"
replace nCH = 1 if gdename == "hoeselt b"
replace nCH = 1 if gdename == "woluwe saint pierre b"
replace nCH = 1 if gdename == "alzenau d"
replace nCH = 1 if gdename == "birkenau d"
replace nCH = 1 if gdename == "mochenwangen d"
replace nCH = 1 if gdename == "charenton le pont f"
replace nCH = 1 if gdename == "fontainebleau (f)"
replace nCH = 1 if gdename == "pembroke bermudas"
replace nCH = 1 if gdename == "elewijt b"
replace nCH = 1 if gdename == "velden nl"
replace nCH = 1 if gdename == "butzbach d"
replace nCH = 1 if gdename == "soedertaelje s"
replace nCH = 1 if gdename == "kalutara cl"
replace nCH = 1 if gdename == "rielasingen (d)"
replace nCH = 1 if gdename == "mons nimy b"
replace nCH = 1 if gdename == "golasecca i"
replace nCH = 1 if gdename == "vagen d"
replace nCH = 1 if gdename == "buttikon/schuebelbach"
replace nCH = 1 if gdename == "woluwe saint lambert b"
replace nCH = 1 if gdename == "lower gwynedd usa"
replace nCH = 1 if gdename == "strassbourg"
replace nCH = 1 if gdename == "pilsbach voecklabruck a"
replace nCH = 1 if gdename == "quebec can"
replace nCH = 1 if gdename == "rye usa"
replace nCH = 1 if gdename == "farringdon, oxon gb"
replace nCH = 1 if gdename == "blasimon f"
replace nCH = 1 if gdename == "nanterre f"
replace nCH = 1 if gdename == "guemligen/muri bei bern"
replace nCH = 1 if gdename == "atherton (usa)"
replace nCH = 1 if gdename == "gustavsberg s"
replace nCH = 1 if gdename == "paradise valley, arizona usa"
replace nCH = 1 if gdename == "laughton, essex gb"
replace nCH = 1 if gdename == "sulzburg d"
replace nCH = 1 if gdename == "peking rc"
replace nCH = 1 if gdename == "eschenbach d"
replace nCH = 1 if gdename == "cernusco lombardone i"
replace nCH = 1 if gdename == "brenzo sk"
replace nCH = 1 if gdename == "ettenhausen tg"
replace nCH = 1 if gdename == "valduggia i"
replace nCH = 1 if gdename == "st doulchard f"
replace nCH = 1 if gdename == "grand zapids usa"
replace nCH = 1 if gdename == "graefelfing (d)"
replace nCH = 1 if gdename == "safat (kwt)"
replace nCH = 1 if gdename == "kristiansand"
replace nCH = 1 if gdename == "nabek syr"
replace nCH = 1 if gdename == "uehlingen birkendorf d"
replace nCH = 1 if gdename == "dietzenbach/steinberg d"
replace nCH = 1 if gdename == "irvine ca, usa"
replace nCH = 1 if gdename == "doha, qatar"
replace nCH = 1 if gdename == "vaduz li"
replace nCH = 1 if gdename == "cesano maderno i"
replace nCH = 1 if gdename == "lahr (d)"
replace nCH = 1 if gdename == "leuven (b)"
replace nCH = 1 if gdename == "margate usa"
replace nCH = 1 if gdename == "dalfsen nl"
replace nCH = 1 if gdename == "aeussere"
replace nCH = 1 if gdename == "plochingen d"
replace nCH = 1 if gdename == "taino i"
replace nCH = 1 if gdename == "st jean d ardieres f"
replace nCH = 1 if gdename == "bettola/piacenza i"
replace nCH = 1 if gdename == "eching d"
replace nCH = 1 if gdename == "bischofsheim d"
replace nCH = 1 if gdename == "haarlem (nl)"
replace nCH = 1 if gdename == "tripolis lt"
replace nCH = 1 if gdename == "kerken eyll (d)"
replace nCH = 1 if gdename == "seattle (usa)"
replace nCH = 1 if gdename == "ramkweil a"
replace nCH = 1 if gdename == "dondelange l"
replace nCH = 1 if gdename == "le sappey f"
replace nCH = 1 if gdename == "etobicoke, ontario ca"
replace nCH = 1 if gdename == "gif sur yvette (f)"
replace nCH = 1 if gdename == "koweit city"
replace nCH = 1 if gdename == "park city utah usa"
replace nCH = 1 if gdename == "cluses (f)"
replace nCH = 1 if gdename == "londre gb"
replace nCH = 1 if gdename == "leusden (nl)"
replace nCH = 1 if gdename == "la chapelle reanville f"
replace nCH = 1 if gdename == "panama city beach, florida usa"
replace nCH = 1 if gdename == "kullavik s"
replace nCH = 1 if gdename == "sovere i"
replace nCH = 1 if gdename == "marietta ga usa"
replace nCH = 1 if gdename == "naples fl usa"
replace nCH = 1 if gdename == "wedel (d)"
replace nCH = 1 if gdename == "leicester gb"
replace nCH = 1 if gdename == "fort smith usa"
replace nCH = 1 if gdename == "east horsley gb"
replace nCH = 1 if gdename == "franklin lakes usa"
replace nCH = 1 if gdename == "bludenz"
replace nCH = 1 if gdename == "reutte a"
replace nCH = 1 if gdename == "palamos e"
replace nCH = 1 if gdename == "oostkamp b"
replace nCH = 1 if gdename == "murcia e"
replace nCH = 1 if gdename == "laab im walde a"
replace nCH = 1 if gdename == "klettgau griessen d"
replace nCH = 1 if gdename == "lyon (france)"
replace nCH = 1 if gdename == "sparta (usa)"
replace nCH = 1 if gdename == "plaisir (f)"
replace nCH = 1 if gdename == "alger (dz)"
replace nCH = 1 if gdename == "st. albans gb"
replace nCH = 1 if gdename == "kansas city usa"
replace nCH = 1 if gdename == "moncalieri i"
replace nCH = 1 if gdename == "haessleholm s"
replace nCH = 1 if gdename == "saint mande f"
replace nCH = 1 if gdename == "rio"
replace nCH = 1 if gdename == "lamberthurst gb"
replace nCH = 1 if gdename == "maisons laffitte (f)"
replace nCH = 1 if gdename == "wertheim waldenhausen (d)"
replace nCH = 1 if gdename == "sunnigdale gb"
replace nCH = 1 if gdename == "ahrensburg d"
replace nCH = 1 if gdename == "westport, connecticut usa"
replace nCH = 1 if gdename == "leun d"
replace nCH = 1 if gdename == "montreal/kanada"
replace nCH = 1 if gdename == "igny f"
replace nCH = 1 if gdename == "canada"
replace nCH = 1 if gdename == "aschaffenburg (d)"
replace nCH = 1 if gdename == "kiev uk"
replace nCH = 1 if gdename == "derby gb"
replace nCH = 1 if gdename == "guenwald d"
replace nCH = 1 if gdename == "caluire et cruire f"
replace nCH = 1 if gdename == "sark, chanel island"
replace nCH = 1 if gdename == "ljublijana yu"
replace nCH = 1 if gdename == "westport ct usa"
replace nCH = 1 if gdename == "olgiate comasco"
replace nCH = 1 if gdename == "worcester gb"
replace nCH = 1 if gdename == "forchtenberg d"
replace nCH = 1 if gdename == "kfar vradim il"
replace nCH = 1 if gdename == "hommersaak n"
replace nCH = 1 if gdename == "winnetka usa"
replace nCH = 1 if gdename == "gubkin rus"
replace nCH = 1 if gdename == "jh genderen nl"
replace nCH = 1 if gdename == "galati ro"
replace nCH = 1 if gdename == "quickborn d"
replace nCH = 1 if gdename == "stenungsund s"
replace nCH = 1 if gdename == "aalst b"
replace nCH = 1 if gdename == "frambouhans/charquemont f"
replace nCH = 1 if gdename == "villadossola i"
replace nCH = 1 if gdename == "san angelo, texas usa"
replace nCH = 1 if gdename == "hadley wood gb"
replace nCH = 1 if gdename == "waregem b"
replace nCH = 1 if gdename == "brentwood, essex gb"
replace nCH = 1 if gdename == "kisangani cgo"
replace nCH = 1 if gdename == "straelen d"
replace nCH = 1 if gdename == "ketsch am rhein d"
replace nCH = 1 if gdename == "wogeningen nl"
replace nCH = 1 if gdename == "south ascot (gb)"
replace nCH = 1 if gdename == "collonges f"
replace nCH = 1 if gdename == "ottawa (cdn)"
replace nCH = 1 if gdename == "phoenix/arizona usa"
replace nCH = 1 if gdename == "ramat aviv il"
replace nCH = 1 if gdename == "weald, sevenoaks gb"
replace nCH = 1 if gdename == "elancourt f"
replace nCH = 1 if gdename == "old greenwich usa"
replace nCH = 1 if gdename == "sandhurst gb"
replace nCH = 1 if gdename == "steyerberg d"
replace nCH = 1 if gdename == "edgware gb"
replace nCH = 1 if gdename == "siena i"
replace nCH = 1 if gdename == "bad mergentheim brd"
replace nCH = 1 if gdename == "poortugaal nl"
replace nCH = 1 if gdename == "taipeh rc"
replace nCH = 1 if gdename == "valbonne f"
replace nCH = 1 if gdename == "stratford upon avon gb"
replace nCH = 1 if gdename == "jona sg"
replace nCH = 1 if gdename == "singen"
replace nCH = 1 if gdename == "neu ulm d"
replace nCH = 1 if gdename == "somerset gb"
replace nCH = 1 if gdename == "pfullendorf"
replace nCH = 1 if gdename == "schoenwald d"
replace nCH = 1 if gdename == "bad driburg (d)"
replace nCH = 1 if gdename == "malbuisson f"
replace nCH = 1 if gdename == "nantes f"
replace nCH = 1 if gdename == "eppingen muehldach d"
replace nCH = 1 if gdename == "mission, kansas usa"
replace nCH = 1 if gdename == "norwood usa"
replace nCH = 1 if gdename == "kinshasa cgo"
replace nCH = 1 if gdename == "frankfurt dl"
replace nCH = 1 if gdename == "mansfield usa"
replace nCH = 1 if gdename == "solna (s)"
replace nCH = 1 if gdename == "st maur, val de marne f"
replace nCH = 1 if gdename == "voisins le bretonneux f"
replace nCH = 1 if gdename == "lutherville, maryland usa"
replace nCH = 1 if gdename == "dinslaken (d)"
replace nCH = 1 if gdename == "taif sa"
replace nCH = 1 if gdename == "la riviere drugeon f"
replace nCH = 1 if gdename == "senta, voijvodina jugosl"
replace nCH = 1 if gdename == "obernkirchen (d)"
replace nCH = 1 if gdename == "springe (d)"
replace nCH = 1 if gdename == "leighton buzzard gb"
replace nCH = 1 if gdename == "hampstead gb"
replace nCH = 1 if gdename == "gmund am tegernsee d"
replace nCH = 1 if gdename == "fourgueux (f)"
replace nCH = 1 if gdename == "lucerna"
replace nCH = 1 if gdename == "palm springs, kalifornien usa"
replace nCH = 1 if gdename == "hausen solscheid d"
replace nCH = 1 if gdename == "karlruhe d"
replace nCH = 1 if gdename == "new port richey usa"
replace nCH = 1 if gdename == "haelsingborg s"
replace nCH = 1 if gdename == "sonderborg dk"
replace nCH = 1 if gdename == "garner, north carolina (usa)"
replace nCH = 1 if gdename == "klein st. paul a"
replace nCH = 1 if gdename == "trecastagni i"
replace nCH = 1 if gdename == "alexandria et"
replace nCH = 1 if gdename == "mashad (ir)"
replace nCH = 1 if gdename == "manama brn"
replace nCH = 1 if gdename == "hinterkappelen/wohlen b. bern"
replace nCH = 1 if gdename == "neuilly / seine f"
replace nCH = 1 if gdename == "mainbernheim d"
replace nCH = 1 if gdename == "flintbek d"
replace nCH = 1 if gdename == "ouderkerk nl"
replace nCH = 1 if gdename == "dreielch d"
replace nCH = 1 if gdename == "freudenstadt d"
replace nCH = 1 if gdename == "atlanta, georgia (usa)"
replace nCH = 1 if gdename == "avon f"
replace nCH = 1 if gdename == "schiltigheim f"
replace nCH = 1 if gdename == "les fins f"
replace nCH = 1 if gdename == "bantry bay, cape town za"
replace nCH = 1 if gdename == "bickley kent gb"
replace nCH = 1 if gdename == "rechy chalais"
replace nCH = 1 if gdename == "villennes sur seine f"
replace nCH = 1 if gdename == "bethel usa"
replace nCH = 1 if gdename == "st martens latem b"
replace nCH = 1 if gdename == "pisa i"
replace nCH = 1 if gdename == "saint heand f"
replace nCH = 1 if gdename == "figino serenza i"
replace nCH = 1 if gdename == "boston, massachusetts (usa)"
replace nCH = 1 if gdename == "aubenas f"
replace nCH = 1 if gdename == "sinsheim d"
replace nCH = 1 if gdename == "messing, essex gb"
replace nCH = 1 if gdename == "edgefield south carolina usa"
replace nCH = 1 if gdename == "salem ind"
replace nCH = 1 if gdename == "troy (usa)"
replace nCH = 1 if gdename == "bonndorf d"
replace nCH = 1 if gdename == "campo s. martino i"
replace nCH = 1 if gdename == "alessandria i"
replace nCH = 1 if gdename == "nellingen d"
replace nCH = 1 if gdename == "zwjndrecht nl"
replace nCH = 1 if gdename == "notodden n"
replace nCH = 1 if gdename == "leiden (nl)"
replace nCH = 1 if gdename == "bensberg (d)"
replace nCH = 1 if gdename == "kaskil jounie rl"
replace nCH = 1 if gdename == "chappaqua usa"
replace nCH = 1 if gdename == "northbrook usa"
replace nCH = 1 if gdename == "messtetten (d)"
replace nCH = 1 if gdename == "mon desir gb"
replace nCH = 1 if gdename == "gaevle (s)"
replace nCH = 1 if gdename == "wertheim d"
replace nCH = 1 if gdename == "bad liebenzell d"
replace nCH = 1 if gdename == "montreuil f"
replace nCH = 1 if gdename == "perl d"
replace nCH = 1 if gdename == "sharjah uae"
replace nCH = 1 if gdename == "lengerich d"
replace nCH = 1 if gdename == "saulgau (d)"
replace nCH = 1 if gdename == "egg/rickenbach d"
replace nCH = 1 if gdename == "moelndal (s)"
replace nCH = 1 if gdename == "lesmo/milano i"
replace nCH = 1 if gdename == "great neck, new york usa"
replace nCH = 1 if gdename == "bognor regis gb"
replace nCH = 1 if gdename == "maisons lafitte f"
replace nCH = 1 if gdename == "lsernhagen d"
replace nCH = 1 if gdename == "gladbeck (d)"
replace nCH = 1 if gdename == "muenchen gruenwald d"
replace nCH = 1 if gdename == "verrieres le buisson f"
replace nCH = 1 if gdename == "chandolin val"
replace nCH = 1 if gdename == "grosse pointe park usa"
replace nCH = 1 if gdename == "trooz grivegnee b"
replace nCH = 1 if gdename == "cruz del eze ra"
replace nCH = 1 if gdename == "kingston (gb)"
replace nCH = 1 if gdename == "viry (f)"
replace nCH = 1 if gdename == "zeil am main d"
replace nCH = 1 if gdename == "zuckenriet (niederhelfenschwil)"
replace nCH = 1 if gdename == "frankfurt/d"
replace nCH = 1 if gdename == "haderslev dk"
replace nCH = 1 if gdename == "norrkoeping s"
replace nCH = 1 if gdename == "harsewinkel (d)"
replace nCH = 1 if gdename == "erlinsbach ag"
replace nCH = 1 if gdename == "krems a"
replace nCH = 1 if gdename == "hohentengen lienheim d"
replace nCH = 1 if gdename == "strasslach d"
replace nCH = 1 if gdename == "arnheim nl"
replace nCH = 1 if gdename == "evanston, illnois usa"
replace nCH = 1 if gdename == "manilva, malaga e"
replace nCH = 1 if gdename == "northeim d"
replace nCH = 1 if gdename == "maenttae (sf)"
replace nCH = 1 if gdename == "grand cayman (b.w.i.)"
replace nCH = 1 if gdename == "bad rappenau d"
replace nCH = 1 if gdename == "fuengirola e"
replace nCH = 1 if gdename == "bietigheim bissingen (d)"
replace nCH = 1 if gdename == "roseland, new jersey usa"
replace nCH = 1 if gdename == "berg a/d maas nl"
replace nCH = 1 if gdename == "hameln (d)"
replace nCH = 1 if gdename == "stains f"
replace nCH = 1 if gdename == "lavenna ponte tresa i"
replace nCH = 1 if gdename == "beaune f"
replace nCH = 1 if gdename == "les herbiers f"
replace nCH = 1 if gdename == "inverness usa"
replace nCH = 1 if gdename == "soderstorf d"
replace nCH = 1 if gdename == "beilstein d"
replace nCH = 1 if gdename == "chambery les vieux/savoi f"
replace nCH = 1 if gdename == "pullach (d)"
replace nCH = 1 if gdename == "sitges (e)"
replace nCH = 1 if gdename == "puntagorta, isola la palma e"
replace nCH = 1 if gdename == "st antonius zoersoel b"
replace nCH = 1 if gdename == "maderno i"
replace nCH = 1 if gdename == "new york, ny usa"
replace nCH = 1 if gdename == "alba le teil f"
replace nCH = 1 if gdename == "sorunda s"
replace nCH = 1 if gdename == "bodenheim d"
replace nCH = 1 if gdename == "cookham dean gb"
replace nCH = 1 if gdename == "ponti i"
replace nCH = 1 if gdename == "arnstorf d"
replace nCH = 1 if gdename == "bisuschio i"
replace nCH = 1 if gdename == "st albans gb"
replace nCH = 1 if gdename == "valsolda i"
replace nCH = 1 if gdename == "michelbach d"
replace nCH = 1 if gdename == "aalst (nl)"
replace nCH = 1 if gdename == "paisley gb"
replace nCH = 1 if gdename == "muehlhausen ehingen d"
replace nCH = 1 if gdename == "viersen (d)"
replace nCH = 1 if gdename == "boston, massachusetts usa"
replace nCH = 1 if gdename == "argenta, british columbia can"
replace nCH = 1 if gdename == "dietingen d"
replace nCH = 1 if gdename == "wassenaar/nl"
replace nCH = 1 if gdename == "sterrebeek b"
replace nCH = 1 if gdename == "farum dk"
replace nCH = 1 if gdename == "pittsburg (usa)"
replace nCH = 1 if gdename == "samarate i"
replace nCH = 1 if gdename == "hennef d"
replace nCH = 1 if gdename == "cernay la ville f"
replace nCH = 1 if gdename == "seib d"
replace nCH = 1 if gdename == "ingensdorf unterlindelbach d"
replace nCH = 1 if gdename == "lincoln, massachusetts"
replace nCH = 1 if gdename == "oberhaching d"
replace nCH = 1 if gdename == "st. veit a.d. glan a"
replace nCH = 1 if gdename == "weyhe d"
replace nCH = 1 if gdename == "kuweit"
replace nCH = 1 if gdename == "timonium (usa)"
replace nCH = 1 if gdename == "damaskus"
replace nCH = 1 if gdename == "stadtlohn d"
replace nCH = 1 if gdename == "tauberbischofsheim (d)"
replace nCH = 1 if gdename == "cologno al serio i"
replace nCH = 1 if gdename == "montlebon f"
replace nCH = 1 if gdename == "huenxe bruckhausen d"
replace nCH = 1 if gdename == "ris orangis"
replace nCH = 1 if gdename == "ryiad sa"
replace nCH = 1 if gdename == "michelstadt d"
replace nCH = 1 if gdename == "zwettl a"
replace nCH = 1 if gdename == "goslar d"
replace nCH = 1 if gdename == "datteln (d)"
replace nCH = 1 if gdename == "lymington gb"
replace nCH = 1 if gdename == "ainring d"
replace nCH = 1 if gdename == "louth, lincolnshire gb"
replace nCH = 1 if gdename == "cocquio trevisago i"
replace nCH = 1 if gdename == "quartier moran berre l etang f"
replace nCH = 1 if gdename == "sandweiler l"
replace nCH = 1 if gdename == "claygate, surrey gb"
replace nCH = 1 if gdename == "jacksonville, florida usa"
replace nCH = 1 if gdename == "arans and"
replace nCH = 1 if gdename == "wilmslow cheshire gb"
replace nCH = 1 if gdename == "kiew, ukraine"
replace nCH = 1 if gdename == "folgensbourg f"
replace nCH = 1 if gdename == "bath, avon gb"
replace nCH = 1 if gdename == "fellbach (d)"
replace nCH = 1 if gdename == "hausach d"
replace nCH = 1 if gdename == "santa margherita ligure i"
replace nCH = 1 if gdename == "grimbergen b"
replace nCH = 1 if gdename == "sliedrecht nl"
replace nCH = 1 if gdename == "wels thalheim a"
replace nCH = 1 if gdename == "uffheim f"
replace nCH = 1 if gdename == "nizza f"
replace nCH = 1 if gdename == "crowthorne gb"
replace nCH = 1 if gdename == "qartier moran berre l etang f"
replace nCH = 1 if gdename == "kuala lumpur (mal)"
replace nCH = 1 if gdename == "jakobsberg s"
replace nCH = 1 if gdename == "hurden sz"
replace nCH = 1 if gdename == "guidizzolo i"
replace nCH = 1 if gdename == "fleurance f"
replace nCH = 1 if gdename == "grubbenvorst (nl)"
replace nCH = 1 if gdename == "kronburg d"
replace nCH = 1 if gdename == "saratoga, kalifornien usa"
replace nCH = 1 if gdename == "siegburg seligenthal d"
replace nCH = 1 if gdename == "bedford usa"
replace nCH = 1 if gdename == "hemstede (nl)"
replace nCH = 1 if gdename == "holland usa"
replace nCH = 1 if gdename == "parainen sf"
replace nCH = 1 if gdename == "andorra and"
replace nCH = 1 if gdename == "hufzen nl"
replace nCH = 1 if gdename == "solingen aufderhoehe d"
replace nCH = 1 if gdename == "rungstedt dk"
replace nCH = 1 if gdename == "allendale, new jersey usa"
replace nCH = 1 if gdename == "huis ter heide nl"
replace nCH = 1 if gdename == "livonia, michigan usa"
replace nCH = 1 if gdename == "stutensee d"
replace nCH = 1 if gdename == "footscray aus"
replace nCH = 1 if gdename == "cincinnati, ohio (usa)"
replace nCH = 1 if gdename == "welwyn gb"
replace nCH = 1 if gdename == "snaroya n"
replace nCH = 1 if gdename == "milan/i"
replace nCH = 1 if gdename == "la boissiere ecole f"
replace nCH = 1 if gdename == "bastad s"
replace nCH = 1 if gdename == "reading, berkshire (gb)"
replace nCH = 1 if gdename == "oberndorf a"
replace nCH = 1 if gdename == "monte carlo cm"
replace nCH = 1 if gdename == "rottach egern d"
replace nCH = 1 if gdename == "jesteburg d"
replace nCH = 1 if gdename == "rm veen nl"
replace nCH = 1 if gdename == "koenigstein dl"
replace nCH = 1 if gdename == "saronni i"
replace nCH = 1 if gdename == "feldkirchen westerham (d)"
replace nCH = 1 if gdename == "bofors (s)"
replace nCH = 1 if gdename == "lanzarote e"
replace nCH = 1 if gdename == "boulogne sur seine (f)"
replace nCH = 1 if gdename == "malaga e"
replace nCH = 1 if gdename == "venray nl"
replace nCH = 1 if gdename == "hittisau a"
replace nCH = 1 if gdename == "manhattan beach usa"
replace nCH = 1 if gdename == "hungerford gb"
replace nCH = 1 if gdename == "peyrins f"
replace nCH = 1 if gdename == "le raincy f"
replace nCH = 1 if gdename == "watford"
replace nCH = 1 if gdename == "new orleans, louisiana (usa)"
replace nCH = 1 if gdename == "boynton beach, fl usa"
replace nCH = 1 if gdename == "bilthoven (nl)"
replace nCH = 1 if gdename == "rockford (usa)"
replace nCH = 1 if gdename == "dresano i"
replace nCH = 1 if gdename == "hammamet"
replace nCH = 1 if gdename == "lavagna i"
replace nCH = 1 if gdename == "dietzhoelztal ewersbach d"
replace nCH = 1 if gdename == "achim bierden d"
replace nCH = 1 if gdename == "ettlingen (d)"
replace nCH = 1 if gdename == "niort f"
replace nCH = 1 if gdename == "eiringen kirchen d"
replace nCH = 1 if gdename == "schilde (b)"
replace nCH = 1 if gdename == "cheam gb"
replace nCH = 1 if gdename == "saltsjoebaden"
replace nCH = 1 if gdename == "tecklenburg d"
replace nCH = 1 if gdename == "iphofen d"
replace nCH = 1 if gdename == "appenweier d"
replace nCH = 1 if gdename == "bois guillaume f"
replace nCH = 1 if gdename == "triest/i"
replace nCH = 1 if gdename == "sasbach obersasbach d"
replace nCH = 1 if gdename == "helfrantzkirch f"
replace nCH = 1 if gdename == "prevessins moens f"
replace nCH = 1 if gdename == "messel d"
replace nCH = 1 if gdename == "huntington ny usa"
replace nCH = 1 if gdename == "nelspruit (suedafrika)"
replace nCH = 1 if gdename == "saint didier au mont d or f"
replace nCH = 1 if gdename == "herbertingen d"
replace nCH = 1 if gdename == "kairo (aegypten)"
replace nCH = 1 if gdename == "ammerbuch d"
replace nCH = 1 if gdename == "panissieres(frankreich)"
replace nCH = 1 if gdename == "richmond va usa"
replace nCH = 1 if gdename == "ludwigshafen a. rh. (d)"
replace nCH = 1 if gdename == "sao paulo, brasilien"
replace nCH = 1 if gdename == "muelheim a.d. ruhr brd"
replace nCH = 1 if gdename == "london, gb"
replace nCH = 1 if gdename == "spanien"
replace nCH = 1 if gdename == "stans a"
replace nCH = 1 if gdename == "tucson arizona usa"
replace nCH = 1 if gdename == "capriasca i"
replace nCH = 1 if gdename == "mont sur marchienne b"
replace nCH = 1 if gdename == "thonon les bains (f)"
replace nCH = 1 if gdename == "vs schwenningen d"
replace nCH = 1 if gdename == "maison alfort f"
replace nCH = 1 if gdename == "riverside ct usa"
replace nCH = 1 if gdename == "st germain saviese"
replace nCH = 1 if gdename == "zwolle nl"
replace nCH = 1 if gdename == "limbourg b"
replace nCH = 1 if gdename == "chesterfield, missouri usa"
replace nCH = 1 if gdename == "junglinster l"
replace nCH = 1 if gdename == "bujumbura ru"
replace nCH = 1 if gdename == "aberdeen gb"
replace nCH = 1 if gdename == "franconville, val d oise f"
replace nCH = 1 if gdename == "mook nl"
replace nCH = 1 if gdename == "loerrach brd"
replace nCH = 1 if gdename == "teheran iran"
replace nCH = 1 if gdename == "well nl"
replace nCH = 1 if gdename == "darien connecticut usa"
replace nCH = 1 if gdename == "bettola i"
replace nCH = 1 if gdename == "viernheim d"
replace nCH = 1 if gdename == "taunton, mass. usa"
replace nCH = 1 if gdename == "county clare irl"
replace nCH = 1 if gdename == "wehr (d)"
replace nCH = 1 if gdename == "overyse b"
replace nCH = 1 if gdename == "carlsbad californie usa"
replace nCH = 1 if gdename == "ljubliana yu"
replace nCH = 1 if gdename == "plettenberg d"
replace nCH = 1 if gdename == "gex, dep. de l ain f"
replace nCH = 1 if gdename == "clarens/montreux"
replace nCH = 1 if gdename == "gentofte (dk)"
replace nCH = 1 if gdename == "bad aibling d"
replace nCH = 1 if gdename == "tjodalyng n"
replace nCH = 1 if gdename == "civitanova i"
replace nCH = 1 if gdename == "varel d"
replace nCH = 1 if gdename == "bischheim f"
replace nCH = 1 if gdename == "leibnitz (a)"
replace nCH = 1 if gdename == "marnes la coquette f"
replace nCH = 1 if gdename == "maison lafitte f"
replace nCH = 1 if gdename == "sulzbach laufen d"
replace nCH = 1 if gdename == "lochham d"
replace nCH = 1 if gdename == "nicosie cy"
replace nCH = 1 if gdename == "mailand(i)"
replace nCH = 1 if gdename == "phoenix, arizona usa"
replace nCH = 1 if gdename == "ste foy les lyon f"
replace nCH = 1 if gdename == "meerbuch d"
replace nCH = 1 if gdename == "pontarlier/doubs f"
replace nCH = 1 if gdename == "warendorf d"
replace nCH = 1 if gdename == "alpena usa"
replace nCH = 1 if gdename == "crowborough gb"
replace nCH = 1 if gdename == "carpi i"
replace nCH = 1 if gdename == "camberley, surrey gb"
replace nCH = 1 if gdename == "torino 1"
replace nCH = 1 if gdename == "djakarta (ri)"
replace nCH = 1 if gdename == "ljubljana slo"
replace nCH = 1 if gdename == "le lude f"
replace nCH = 1 if gdename == "east brunswick, new jersey usa"
replace nCH = 1 if gdename == "elixhausen a"
replace nCH = 1 if gdename == "neuilly s/ seine f"
replace nCH = 1 if gdename == "moscow gus"
replace nCH = 1 if gdename == "nonnenhorn d"
replace nCH = 1 if gdename == "hechingen (d)"
replace nCH = 1 if gdename == "juelich d"
replace nCH = 1 if gdename == "reute d"
replace nCH = 1 if gdename == "eupen b"
replace nCH = 1 if gdename == "st. ingbert d"
replace nCH = 1 if gdename == "hyogoken japan"
replace nCH = 1 if gdename == "algier"
replace nCH = 1 if gdename == "weubrugge nl"
replace nCH = 1 if gdename == "witten herbede d"
replace nCH = 1 if gdename == "ickenham (gb)"
replace nCH = 1 if gdename == "oosterbeek nl"
replace nCH = 1 if gdename == "dion valmont b"
replace nCH = 1 if gdename == "ikesa i"
replace nCH = 1 if gdename == "garmisch partenkirchen (d)"
replace nCH = 1 if gdename == "dubai eau"
replace nCH = 1 if gdename == "franklin lakes nj usa"
replace nCH = 1 if gdename == "porlezza (i)"
replace nCH = 1 if gdename == "bloemendaal (nl)"
replace nCH = 1 if gdename == "divonne f"
replace nCH = 1 if gdename == "almese i"
replace nCH = 1 if gdename == "new south wales aus"
replace nCH = 1 if gdename == "biedermannsdorf a"
replace nCH = 1 if gdename == "frohnleiten (a)"
replace nCH = 1 if gdename == "ibbenbueren d"
replace nCH = 1 if gdename == "riverside usa"
replace nCH = 1 if gdename == "cresskil, new jersey usa"
replace nCH = 1 if gdename == "argenteuil f"
replace nCH = 1 if gdename == "brussel"
replace nCH = 1 if gdename == "nokia sf"
replace nCH = 1 if gdename == "pekin cn"
replace nCH = 1 if gdename == "muehlheim ruhr d"
replace nCH = 1 if gdename == "schwaz (a)"
replace nCH = 1 if gdename == "salins les bains f"
replace nCH = 1 if gdename == "arcata, calif. usa"
replace nCH = 1 if gdename == "uetlibergstr. 286, hamburg d"
replace nCH = 1 if gdename == "hayward/kalif"
replace nCH = 1 if gdename == "sprendlingen (d)"
replace nCH = 1 if gdename == "ecully (f)"
replace nCH = 1 if gdename == "chabeuil f"
replace nCH = 1 if gdename == "trelleborg (s)"
replace nCH = 1 if gdename == "wembley gb"
replace nCH = 1 if gdename == "carobbio degli angeli i"
replace nCH = 1 if gdename == "wiltz l"
replace nCH = 1 if gdename == "zirl a"
replace nCH = 1 if gdename == "isle adam (f)"
replace nCH = 1 if gdename == "neufahrn d"
replace nCH = 1 if gdename == "bad harzburg (d)"
replace nCH = 1 if gdename == "deu delhi ind"
replace nCH = 1 if gdename == "yoko hama j"
replace nCH = 1 if gdename == "biberach/d"
replace nCH = 1 if gdename == "bombey"
replace nCH = 1 if gdename == "steenwyk nl"
replace nCH = 1 if gdename == "pl noorddwijkerhout (nl)"
replace nCH = 1 if gdename == "frankfurt m"
replace nCH = 1 if gdename == "bourges f"
replace nCH = 1 if gdename == "(springen d"
replace nCH = 1 if gdename == "lerici i"
replace nCH = 1 if gdename == "westhausen d"
replace nCH = 1 if gdename == "golden, colorado usa"
replace nCH = 1 if gdename == "cattolica (i)"
replace nCH = 1 if gdename == "wayzata, minnesota usa"
replace nCH = 1 if gdename == "rye/new york usa"
replace nCH = 1 if gdename == "fredensborg (dk)"
replace nCH = 1 if gdename == "bierden brd"
replace nCH = 1 if gdename == "vaestervik s"
replace nCH = 1 if gdename == "le rhoon nl"
replace nCH = 1 if gdename == "mateo usa"
replace nCH = 1 if gdename == "lohja sf"
replace nCH = 1 if gdename == "uedem d"
replace nCH = 1 if gdename == "prestbury gb"
replace nCH = 1 if gdename == "pinelands za"
replace nCH = 1 if gdename == "nacka s"
replace nCH = 1 if gdename == "bakou, aserbaidjan, udssr"
replace nCH = 1 if gdename == "brentwood essex gb"
replace nCH = 1 if gdename == "santa barbara, kalifornien usa"
replace nCH = 1 if gdename == "frankfurt a. m"
replace nCH = 1 if gdename == "loughton gb"
replace nCH = 1 if gdename == "albino i"
replace nCH = 1 if gdename == "uppsala (s)"
replace nCH = 1 if gdename == "bukarest ro"
replace nCH = 1 if gdename == "roswell usa"
replace nCH = 1 if gdename == "horgau d"
replace nCH = 1 if gdename == "kapelle op den bos b"
replace nCH = 1 if gdename == "burgbrohl (d)"
replace nCH = 1 if gdename == "berkeroed dk"
replace nCH = 1 if gdename == "esbjerg dk"
replace nCH = 1 if gdename == "baltimore usa"
replace nCH = 1 if gdename == "milwaukee (usa)"
replace nCH = 1 if gdename == "litzelstetten d"
replace nCH = 1 if gdename == "bad ischel a"
replace nCH = 1 if gdename == "bindlach d"
replace nCH = 1 if gdename == "wilmslow gb"
replace nCH = 1 if gdename == "cergy f"
replace nCH = 1 if gdename == "yokohoma j"
replace nCH = 1 if gdename == "tulsa, oklahoma usa"
replace nCH = 1 if gdename == "buggiano castello i"
replace nCH = 1 if gdename == "ridderkerk nl"
replace nCH = 1 if gdename == "av. des bains 20, mondorf les bains l"
replace nCH = 1 if gdename == "labro rieti i"
replace nCH = 1 if gdename == "fernwald d"
replace nCH = 1 if gdename == "tienebronn (d)"
replace nCH = 1 if gdename == "cuyahoga falls usa"
replace nCH = 1 if gdename == "thames ditton, surrey gb"
replace nCH = 1 if gdename == "court st etienne b"
replace nCH = 1 if gdename == "slough gb"
replace nCH = 1 if gdename == "frankenberg d"
replace nCH = 1 if gdename == "serravalle scrivia i"
replace nCH = 1 if gdename == "florence"
replace nCH = 1 if gdename == "deinze b"
replace nCH = 1 if gdename == "karben (d)"
replace nCH = 1 if gdename == "san isidro pe"
replace nCH = 1 if gdename == "bjeringbro dk"
replace nCH = 1 if gdename == "illhaeusern f"
replace nCH = 1 if gdename == "al khobar sa"
replace nCH = 1 if gdename == "kp veendam nl"
replace nCH = 1 if gdename == "umstadt (d)"
replace nCH = 1 if gdename == "tienen b"
replace nCH = 1 if gdename == "surry (gb)"
replace nCH = 1 if gdename == "maspalomas, gran canaria e"
replace nCH = 1 if gdename == "hillside, alderley edge, cheshire gb"
replace nCH = 1 if gdename == "huntsville, alabama usa"
replace nCH = 1 if gdename == "sesto fiorentino (i)"
replace nCH = 1 if gdename == "contamine sur arve f"
replace nCH = 1 if gdename == "scottsdale, arizona usa"
replace nCH = 1 if gdename == "kufstein (a)"
replace nCH = 1 if gdename == "garatshausen d"
replace nCH = 1 if gdename == "erding d"
replace nCH = 1 if gdename == "st louis usa"
replace nCH = 1 if gdename == "bludenz (a)"
replace nCH = 1 if gdename == "santa clara usa"
replace nCH = 1 if gdename == "atlanta (usa)"
replace nCH = 1 if gdename == "bretigny sur orges f"
replace nCH = 1 if gdename == "canmore, alberta cdn"
replace nCH = 1 if gdename == "treviso (i)"
replace nCH = 1 if gdename == "oakville cdn"
replace nCH = 1 if gdename == "rotterdam/nl"
replace nCH = 1 if gdename == "heidelberg dl"
replace nCH = 1 if gdename == "merzhausen d"
replace nCH = 1 if gdename == "haiger d"
replace nCH = 1 if gdename == "kingston jamaica"
replace nCH = 1 if gdename == "baech/freienbach"
replace nCH = 1 if gdename == "wolfsberg a"
replace nCH = 1 if gdename == "abu dhabi arabische emirate"
replace nCH = 1 if gdename == "ornex (f)"
replace nCH = 1 if gdename == "gennevilliers f"
replace nCH = 1 if gdename == "tassin la demi lune f"
replace nCH = 1 if gdename == "chembery f"
replace nCH = 1 if gdename == "dushanbe tj"
replace nCH = 1 if gdename == "amman (hkj)"
replace nCH = 1 if gdename == "frankendorf"
replace nCH = 1 if gdename == "wevelgem b"
replace nCH = 1 if gdename == "viterbo i"
replace nCH = 1 if gdename == "borgosesia i"
replace nCH = 1 if gdename == "couvin b"
replace nCH = 1 if gdename == "alpicat e"
replace nCH = 1 if gdename == "rohrdorf d"
replace nCH = 1 if gdename == "vezelay f"
replace nCH = 1 if gdename == "nottuln d"
replace nCH = 1 if gdename == "lanskroun cfsr"
replace nCH = 1 if gdename == "greeba, isle of man"
replace nCH = 1 if gdename == "noisy le roi (f)"
replace nCH = 1 if gdename == "weissach d"
replace nCH = 1 if gdename == "davesco sorangno"
replace nCH = 1 if gdename == "edinbourg"
replace nCH = 1 if gdename == "bagdad iro"
replace nCH = 1 if gdename == "puurs b"
replace nCH = 1 if gdename == "bagno a ripoli (i)"
replace nCH = 1 if gdename == "casale corte cerro i"
replace nCH = 1 if gdename == "lohne d"
replace nCH = 1 if gdename == "yanbu (sa)"
replace nCH = 1 if gdename == "langenargen d"
replace nCH = 1 if gdename == "neckarsulm d"
replace nCH = 1 if gdename == "cortina d ampezzo i"
replace nCH = 1 if gdename == "potenza i"
replace nCH = 1 if gdename == "porto ceresio i"
replace nCH = 1 if gdename == "monceaux le compte f"
replace nCH = 1 if gdename == "luxembourg kirchberg l"
replace nCH = 1 if gdename == "parigi f"
replace nCH = 1 if gdename == "guilderland, ny usa"
replace nCH = 1 if gdename == "boennigheim d"
replace nCH = 1 if gdename == "nurcia e"
replace nCH = 1 if gdename == "germignaga i"
replace nCH = 1 if gdename == "schwaz a"
replace nCH = 1 if gdename == "landskrona s"
replace nCH = 1 if gdename == "ancy le libre f"
replace nCH = 1 if gdename == "grankulla sf"
replace nCH = 1 if gdename == "kehl kork d"
replace nCH = 1 if gdename == "pietrasanta fraz marina i"
replace nCH = 1 if gdename == "solihull (gb)"
replace nCH = 1 if gdename == "san francisco, californie usa"
replace nCH = 1 if gdename == "bad berleburg d"
replace nCH = 1 if gdename == "accra (ghana)"
replace nCH = 1 if gdename == "kappel grafenhausen"
replace nCH = 1 if gdename == "backnang d"
replace nCH = 1 if gdename == "rosental/waengi"
replace nCH = 1 if gdename == "luven b"
replace nCH = 1 if gdename == "belphraon"
replace nCH = 1 if gdename == "torrance (california)"
replace nCH = 1 if gdename == "castlecrag aus"
replace nCH = 1 if gdename == "castelnuovo di baradegna i"
replace nCH = 1 if gdename == "dubai (uae)"
replace nCH = 1 if gdename == "schachen d"
replace nCH = 1 if gdename == "hohenau a"
replace nCH = 1 if gdename == "laguna niguel, kalifornien usa"
replace nCH = 1 if gdename == "medford usa"
replace nCH = 1 if gdename == "kolding dae"
replace nCH = 1 if gdename == "chatou (f)"
replace nCH = 1 if gdename == "collonge sous saleve (f)"
replace nCH = 1 if gdename == "marienheide d"
replace nCH = 1 if gdename == "steenbergen nb (nl)"
replace nCH = 1 if gdename == "edenbridge, kent gb"
replace nCH = 1 if gdename == "kist d"
replace nCH = 1 if gdename == "poue (f)"
replace nCH = 1 if gdename == "paratico i"
replace nCH = 1 if gdename == "reading, berkshire gb"
replace nCH = 1 if gdename == "sunnyvale /california, usa"
replace nCH = 1 if gdename == "chicago (illinois/usa)"
replace nCH = 1 if gdename == "albiolo i"
replace nCH = 1 if gdename == "hamilton/bermuda"
replace nCH = 1 if gdename == "signes f"
replace nCH = 1 if gdename == "valdagno i"
replace nCH = 1 if gdename == "epinal f"
replace nCH = 1 if gdename == "rescaldina i"
replace nCH = 1 if gdename == "guilford (gb)"
replace nCH = 1 if gdename == "standton za"
replace nCH = 1 if gdename == "southport usa"
replace nCH = 1 if gdename == "saltsjobaden (s)"
replace nCH = 1 if gdename == "ikeja wan"
replace nCH = 1 if gdename == "eppelheim d"
replace nCH = 1 if gdename == "hochstatt f"
replace nCH = 1 if gdename == "herk de stad b"
replace nCH = 1 if gdename == "schanghai china"
replace nCH = 1 if gdename == "tokyo jap"
replace nCH = 1 if gdename == "schwaig d"
replace nCH = 1 if gdename == "helecine b"
replace nCH = 1 if gdename == "stavanger n"
replace nCH = 1 if gdename == "hove gb"
replace nCH = 1 if gdename == "muelheim a.d.ruhr d"
replace nCH = 1 if gdename == "obertshausen (d)"
replace nCH = 1 if gdename == "illertissen d"
replace nCH = 1 if gdename == "sint genesius rode b"
replace nCH = 1 if gdename == "wohltorf d"
replace nCH = 1 if gdename == "wilmslow"
replace nCH = 1 if gdename == "lathi (sf)"
replace nCH = 1 if gdename == "tampere fin"
replace nCH = 1 if gdename == "kleineheubach d"
replace nCH = 1 if gdename == "missisauga, ontario can"
replace nCH = 1 if gdename == "croydon (surrey) gb"
replace nCH = 1 if gdename == "sechelt cdn"
replace nCH = 1 if gdename == "sasso fortino i"
replace nCH = 1 if gdename == "wesel d"
replace nCH = 1 if gdename == "istanbul tur"
replace nCH = 1 if gdename == "singapur (sgp)"
replace nCH = 1 if gdename == "seminole, florida usa"
replace nCH = 1 if gdename == "riad"
replace nCH = 1 if gdename == "malagan e"
replace nCH = 1 if gdename == "avala alabang rp"
replace nCH = 1 if gdename == "radda in chianti i"
replace nCH = 1 if gdename == "stone gb"
replace nCH = 1 if gdename == "anif a"
replace nCH = 1 if gdename == "sandton sa"
replace nCH = 1 if gdename == "obersiggenthal kirchdorf"
replace nCH = 1 if gdename == "muehlheim an der ruhr d"
replace nCH = 1 if gdename == "sevrier f"
replace nCH = 1 if gdename == "montluel f"
replace nCH = 1 if gdename == "yakutsk (repubblica sakha) rus"
replace nCH = 1 if gdename == "vorbasse dk"
replace nCH = 1 if gdename == "san mauro (i)"
replace nCH = 1 if gdename == "curepipe mus"
replace nCH = 1 if gdename == "rio claro brasilien"
replace nCH = 1 if gdename == "newport gb"
replace nCH = 1 if gdename == "migne auxances f"
replace nCH = 1 if gdename == "oberkirch bottenau d"
replace nCH = 1 if gdename == "calw stammheim (d)"
replace nCH = 1 if gdename == "bagdad (iro)"
replace nCH = 1 if gdename == "tarragona e"
replace nCH = 1 if gdename == "la tour de salvagny f"
replace nCH = 1 if gdename == "perrignier"
replace nCH = 1 if gdename == "sudbury, ma usa"
replace nCH = 1 if gdename == "moskau (su)"
replace nCH = 1 if gdename == "backnang (d)"
replace nCH = 1 if gdename == "gaildorf (d)"
replace nCH = 1 if gdename == "les bruyeres du haut, chiddes f"
replace nCH = 1 if gdename == "aerzen d"
replace nCH = 1 if gdename == "kingston usa"
replace nCH = 1 if gdename == "jersey (channel islands)"
replace nCH = 1 if gdename == "potomac usa"
replace nCH = 1 if gdename == "barcelona/e"
replace nCH = 1 if gdename == "racine wi usa"
replace nCH = 1 if gdename == "bindura rsr"
replace nCH = 1 if gdename == "wertheim am main d"
replace nCH = 1 if gdename == "euskirchen (d)"
replace nCH = 1 if gdename == "zagreb (kroatien)"
replace nCH = 1 if gdename == "bilzen b"
replace nCH = 1 if gdename == "njagan, russland"
replace nCH = 1 if gdename == "katy, texas usa"
replace nCH = 1 if gdename == "dewelier"
replace nCH = 1 if gdename == "brooks, alberta (kanada)"
replace nCH = 1 if gdename == "saltsjoebaden (schweden)"
replace nCH = 1 if gdename == "carlazzo i"
replace nCH = 1 if gdename == "st cloud, hauts de seine f"
replace nCH = 1 if gdename == "melrose park usa"
replace nCH = 1 if gdename == "espo sf"
replace nCH = 1 if gdename == "muehlbach am hkg a"
replace nCH = 1 if gdename == "mouy f"
replace nCH = 1 if gdename == "nimes f"
replace nCH = 1 if gdename == "puerto plata dom"
replace nCH = 1 if gdename == "schermbeck d"
replace nCH = 1 if gdename == "melun f"
replace nCH = 1 if gdename == "verdes usa"
replace nCH = 1 if gdename == "chicago, illinois (usa)"
replace nCH = 1 if gdename == "humington beach (ca) usa"
replace nCH = 1 if gdename == "ketsch d"
replace nCH = 1 if gdename == "monnetier mornex (f)"
replace nCH = 1 if gdename == "brentwood gb"
replace nCH = 1 if gdename == "croix (f)"
replace nCH = 1 if gdename == "koksjide b"
replace nCH = 1 if gdename == "pfaffenhofen/ilm d"
replace nCH = 1 if gdename == "illertissen (d)"
replace nCH = 1 if gdename == "booton usa"
replace nCH = 1 if gdename == "leimen d"
replace nCH = 1 if gdename == "deepdene aus"
replace nCH = 1 if gdename == "lubljana yu"
replace nCH = 1 if gdename == "manly aus"
replace nCH = 1 if gdename == "steinheim d"
replace nCH = 1 if gdename == "sundsvall s"
replace nCH = 1 if gdename == "illertissen betlinshausen d"
replace nCH = 1 if gdename == "muehlheim am main d"
replace nCH = 1 if gdename == "framersheim d"
replace nCH = 1 if gdename == "st maxime f"
replace nCH = 1 if gdename == "yaounde cm"
replace nCH = 1 if gdename == "effingham, surrey gb"
replace nCH = 1 if gdename == "trostberg d"
replace nCH = 1 if gdename == "brusimpiano i"
replace nCH = 1 if gdename == "new malden, surrey gb"
replace nCH = 1 if gdename == "spresiano i"
replace nCH = 1 if gdename == "kreuztal d"
replace nCH = 1 if gdename == "nassau (d)"
replace nCH = 1 if gdename == "comerio i"
replace nCH = 1 if gdename == "doha (katar)"
replace nCH = 1 if gdename == "dresher, pa usa"
replace nCH = 1 if gdename == "olot e"
replace nCH = 1 if gdename == "bredsten (dk)"
replace nCH = 1 if gdename == "fenis i"
replace nCH = 1 if gdename == "muscatine (usa)"
replace nCH = 1 if gdename == "djeddah"
replace nCH = 1 if gdename == "aichwald lobenrot d"
replace nCH = 1 if gdename == "ste martine, quebec (cdn)"
replace nCH = 1 if gdename == "epe nl"
replace nCH = 1 if gdename == "oxford mills cnd"
replace nCH = 1 if gdename == "weston, mass. usa"
replace nCH = 1 if gdename == "karibib/namibia"
replace nCH = 1 if gdename == "staines middx gb"
replace nCH = 1 if gdename == "megliadino san fidenzio i"
replace nCH = 1 if gdename == "kinshasa (cgo)"
replace nCH = 1 if gdename == "iwaki city j"
replace nCH = 1 if gdename == "south molton gb"
replace nCH = 1 if gdename == "bratislava cssr"
replace nCH = 1 if gdename == "huettenberg rechtenbach d"
replace nCH = 1 if gdename == "vedano al lambro i"
replace nCH = 1 if gdename == "ludesch a"
replace nCH = 1 if gdename == "providence, rhode island usa"
replace nCH = 1 if gdename == "santos br"
replace nCH = 1 if gdename == "memphis, tennessee usa"
replace nCH = 1 if gdename == "chislehurst, kent gb"
replace nCH = 1 if gdename == "arcen nl"
replace nCH = 1 if gdename == "westbrook, maine usa"
replace nCH = 1 if gdename == "dublin ir"
replace nCH = 1 if gdename == "hopkinton usa"
replace nCH = 1 if gdename == "wesel dl"
replace nCH = 1 if gdename == "zoagli i"
replace nCH = 1 if gdename == "gundelfingen d"
replace nCH = 1 if gdename == "portsmouth usa"
replace nCH = 1 if gdename == "spanien e"
replace nCH = 1 if gdename == "kuurne b"
replace nCH = 1 if gdename == "corning n.y. usa"
replace nCH = 1 if gdename == "curacao cu"
replace nCH = 1 if gdename == "bussum"
replace nCH = 1 if gdename == "asti i"
replace nCH = 1 if gdename == "aichwald aichschiess d"
replace nCH = 1 if gdename == "havre f"
replace nCH = 1 if gdename == "san francisco, calif. usa"
replace nCH = 1 if gdename == "kempen d"
replace nCH = 1 if gdename == "riyad/arab"
replace nCH = 1 if gdename == "gerona e"
replace nCH = 1 if gdename == "hof (d)"
replace nCH = 1 if gdename == "ocean front venice cal. usa"
replace nCH = 1 if gdename == "allaines f"
replace nCH = 1 if gdename == "marina di carrara i"
replace nCH = 1 if gdename == "broxbourne gb"
replace nCH = 1 if gdename == "djeddah hz"
replace nCH = 1 if gdename == "st ismier f"
replace nCH = 1 if gdename == "wimborne gb"
replace nCH = 1 if gdename == "heiligenhausen d"
replace nCH = 1 if gdename == "miami fl, usa"
replace nCH = 1 if gdename == "welwyn (gb)"
replace nCH = 1 if gdename == "odder dk"
replace nCH = 1 if gdename == "concord, mass. usa"
replace nCH = 1 if gdename == "laveno mombello i"
replace nCH = 1 if gdename == "st joseph, missouri us"
replace nCH = 1 if gdename == "la muraz f"
replace nCH = 1 if gdename == "beekbergen nl"
replace nCH = 1 if gdename == "les pavillons sous bois f"
replace nCH = 1 if gdename == "dardilly f"
replace nCH = 1 if gdename == "mon desir, sark gb"
replace nCH = 1 if gdename == "nacka (s)"
replace nCH = 1 if gdename == "toulon f"
replace nCH = 1 if gdename == "kerridge macclesfield gb"
replace nCH = 1 if gdename == "ivry sur seine f"
replace nCH = 1 if gdename == "neuenbuerg (d)"
replace nCH = 1 if gdename == "itzehoe d"
replace nCH = 1 if gdename == "st aubin des ormeaux f"
replace nCH = 1 if gdename == "bolzano i"
replace nCH = 1 if gdename == "scarsdale (usa)"
replace nCH = 1 if gdename == "westmount cdn"
replace nCH = 1 if gdename == "houston texas (usa)"
replace nCH = 1 if gdename == "wasselonne f"
replace nCH = 1 if gdename == "great falls usa"
replace nCH = 1 if gdename == "prevessin moens"
replace nCH = 1 if gdename == "waldegg a"
replace nCH = 1 if gdename == "st. georgen peterzell (d)"
replace nCH = 1 if gdename == "gaissau (a)"
replace nCH = 1 if gdename == "roodeport za"
replace nCH = 1 if gdename == "fleet, hampshire gb"
replace nCH = 1 if gdename == "heliopolis et"
replace nCH = 1 if gdename == "bethesda md usa"
replace nCH = 1 if gdename == "werne d"
replace nCH = 1 if gdename == "usingen kransberg d"
replace nCH = 1 if gdename == "sault ste marie (cdn)"
replace nCH = 1 if gdename == "vinkeveen nl"
replace nCH = 1 if gdename == "braeunlingen d"
replace nCH = 1 if gdename == "luebeck travemuende d"
replace nCH = 1 if gdename == "truchtersheim f"
replace nCH = 1 if gdename == "nassau bh"
replace nCH = 1 if gdename == "sesto san giovanni i"
replace nCH = 1 if gdename == "staldenwang d"
replace nCH = 1 if gdename == "fort myers usa"
replace nCH = 1 if gdename == "ruetlingen d"
replace nCH = 1 if gdename == "octe ville sur mer (f)"
replace nCH = 1 if gdename == "obergruenburg a"
replace nCH = 1 if gdename == "wuelfrath d"
replace nCH = 1 if gdename == "tripoli lybien"
replace nCH = 1 if gdename == "kuenzelsau (d)"
replace nCH = 1 if gdename == "acheres f"
replace nCH = 1 if gdename == "gonesse f"
replace nCH = 1 if gdename == "villeneuve de berg f"
replace nCH = 1 if gdename == "bangalore ind"
replace nCH = 1 if gdename == "uden (nl)"
replace nCH = 1 if gdename == "vacoas mus"
replace nCH = 1 if gdename == "thorens glieres f"
replace nCH = 1 if gdename == "gingelom b"
replace nCH = 1 if gdename == "weybridge (gb)"
replace nCH = 1 if gdename == "harare zw"
replace nCH = 1 if gdename == "santa clara, ca usa"
replace nCH = 1 if gdename == "hausenstamm (d)"
replace nCH = 1 if gdename == "essen werden d"
replace nCH = 1 if gdename == "birmingham (gb)"
replace nCH = 1 if gdename == "aelvsjoe s"
replace nCH = 1 if gdename == "horsholm (dk)"
replace nCH = 1 if gdename == "jamestown usa"
replace nCH = 1 if gdename == "zeeland usa"
replace nCH = 1 if gdename == "st. peter port, guernesey gb"
replace nCH = 1 if gdename == "burgdorf d"
replace nCH = 1 if gdename == "lecce i"
replace nCH = 1 if gdename == "chilly marazin f"
replace nCH = 1 if gdename == "boulogne billancourt, sceaux f"
replace nCH = 1 if gdename == "orebro s"
replace nCH = 1 if gdename == "arolsen (d)"
replace nCH = 1 if gdename == "bilbao (e)"
replace nCH = 1 if gdename == "chalfont st. giles gb"
replace nCH = 1 if gdename == "neviges (d)"
replace nCH = 1 if gdename == "piossasco i"
replace nCH = 1 if gdename == "ivrea (i)"
replace nCH = 1 if gdename == "luxenburg (l)"
replace nCH = 1 if gdename == "laveno ponte tresa i"
replace nCH = 1 if gdename == "seefeld (a)"
replace nCH = 1 if gdename == "buebingen d"
replace nCH = 1 if gdename == "niddatal d"
replace nCH = 1 if gdename == "budenheim d"
replace nCH = 1 if gdename == "giarre, catania i"
replace nCH = 1 if gdename == "crumstadt"
replace nCH = 1 if gdename == "monte pedreguer e"
replace nCH = 1 if gdename == "oostvoorne nl"
replace nCH = 1 if gdename == "bishops gb"
replace nCH = 1 if gdename == "berchen antwerpen b"
replace nCH = 1 if gdename == "lbach schwyz"
replace nCH = 1 if gdename == "wapenveld nl"
replace nCH = 1 if gdename == "wettswil a. a"
replace nCH = 1 if gdename == "ansbach d"
replace nCH = 1 if gdename == "gibraltar gbz"
replace nCH = 1 if gdename == "altdorf bei nuernberg d"
replace nCH = 1 if gdename == "westmount ca"
replace nCH = 1 if gdename == "cortina (i)"
replace nCH = 1 if gdename == "montale i"
replace nCH = 1 if gdename == "kayseri tr"
replace nCH = 1 if gdename == "graben neudorf d"
replace nCH = 1 if gdename == "arsizio"
replace nCH = 1 if gdename == "herouvillette f"
replace nCH = 1 if gdename == "columbus oh usa"
replace nCH = 1 if gdename == "liege (b)"
replace nCH = 1 if gdename == "dubai (vae)"
replace nCH = 1 if gdename == "sandhurst za"
replace nCH = 1 if gdename == "monkton, maryland usa"
replace nCH = 1 if gdename == "irwin, pennsyl. usa"
replace nCH = 1 if gdename == "newton, massachusetts usa"
replace nCH = 1 if gdename == "ho chi minh ville vn"
replace nCH = 1 if gdename == "mascate ma"
replace nCH = 1 if gdename == "cebu city rp"
replace nCH = 1 if gdename == "siegburg braschoss (d)"
replace nCH = 1 if gdename == "adelberg d"
replace nCH = 1 if gdename == "endingen d"
replace nCH = 1 if gdename == "shorewood mn usa"
replace nCH = 1 if gdename == "alleroed dk"
replace nCH = 1 if gdename == "kimberley dcn"
replace nCH = 1 if gdename == "dworp b"
replace nCH = 1 if gdename == "villaricca i"
replace nCH = 1 if gdename == "oud beijerland nl"
replace nCH = 1 if gdename == "instanbul (tr)"
replace nCH = 1 if gdename == "barlassina i"
replace nCH = 1 if gdename == "gilley f"
replace nCH = 1 if gdename == "bryanston ras"
replace nCH = 1 if gdename == "valmorea i"
replace nCH = 1 if gdename == "wangen i. allgaeu d"
replace nCH = 1 if gdename == "huntington gb"
replace nCH = 1 if gdename == "quay house (guernsey)"
replace nCH = 1 if gdename == "massa i"
replace nCH = 1 if gdename == "fuschl am see a"
replace nCH = 1 if gdename == "ferdinand"
replace nCH = 1 if gdename == "bidermannsdorf a"
replace nCH = 1 if gdename == "riolo terme i"
replace nCH = 1 if gdename == "brighton (gb)"
replace nCH = 1 if gdename == "leichlingen d"
replace nCH = 1 if gdename == "edingburgh"
replace nCH = 1 if gdename == "warrenhurst gb"
replace nCH = 1 if gdename == "westport, conn. (usa)"
replace nCH = 1 if gdename == "buschwiler f"
replace nCH = 1 if gdename == "duffel b"
replace nCH = 1 if gdename == "ranica i"
replace nCH = 1 if gdename == "cachen/seine f"
replace nCH = 1 if gdename == "fountain hills (usa)"
replace nCH = 1 if gdename == "felixtowe gb"
replace nCH = 1 if gdename == "gordes f"
replace nCH = 1 if gdename == "barnet gb"
replace nCH = 1 if gdename == "aalesund n"
replace nCH = 1 if gdename == "berlin dl"
replace nCH = 1 if gdename == "bethel (usa)"
replace nCH = 1 if gdename == "owingen d"
replace nCH = 1 if gdename == "ayze f"
replace nCH = 1 if gdename == "struer dk"
replace nCH = 1 if gdename == "riyadh"
replace nCH = 1 if gdename == "st malo f"
replace nCH = 1 if gdename == "randersacker d"
replace nCH = 1 if gdename == "muelheim a.d. ruhr (d)"
replace nCH = 1 if gdename == "bruegge b"
replace nCH = 1 if gdename == "los altos hills, ca, usa"
replace nCH = 1 if gdename == "chennai ind"
replace nCH = 1 if gdename == "eymines meylan f"
replace nCH = 1 if gdename == "st. nom la breteche f"
replace nCH = 1 if gdename == "vaerloese dk"
replace nCH = 1 if gdename == "east amhearst, new york usa"
replace nCH = 1 if gdename == "ganleigh gb"
replace nCH = 1 if gdename == "huskvarna s"
replace nCH = 1 if gdename == "lans f"
replace nCH = 1 if gdename == "collecchio i"
replace nCH = 1 if gdename == "lasne (b)"
replace nCH = 1 if gdename == "kremnica"
replace nCH = 1 if gdename == "perchtoldsdorf (a)"
replace nCH = 1 if gdename == "mexiko city mex"
replace nCH = 1 if gdename == "gruenwald b. muenchen (d)"
replace nCH = 1 if gdename == "northport usa"
replace nCH = 1 if gdename == "rampton gb"
replace nCH = 1 if gdename == "dietzenbach frankfurt d"
replace nCH = 1 if gdename == "tilburg (nl)"
replace nCH = 1 if gdename == "bad reichenhall d"
replace nCH = 1 if gdename == "oropesa del mar e"
replace nCH = 1 if gdename == "englewood, colorado usa"
replace nCH = 1 if gdename == "palma de mallorca (e)"
replace nCH = 1 if gdename == "puerto galea rp"
replace nCH = 1 if gdename == "southport, conn. (usa)"
replace nCH = 1 if gdename == "nykobing dk"
replace nCH = 1 if gdename == "langenhagen (d)"
replace nCH = 1 if gdename == "fenegro i"
replace nCH = 1 if gdename == "arezzo i"
replace nCH = 1 if gdename == "beersel b"
replace nCH = 1 if gdename == "villingen schwenningen (d)"
replace nCH = 1 if gdename == "willowdale, ontario cdn"
replace nCH = 1 if gdename == "haag d"
replace nCH = 1 if gdename == "giessen nl"
replace nCH = 1 if gdename == "kasara sud"
replace nCH = 1 if gdename == "harrow gb"
replace nCH = 1 if gdename == "krakow pl"
replace nCH = 1 if gdename == "bickliy kent gb"
replace nCH = 1 if gdename == "siebnen galgenen"
replace nCH = 1 if gdename == "weatherford ok usa"
replace nCH = 1 if gdename == "panama city pan"
replace nCH = 1 if gdename == "guilford, surrey gb"
replace nCH = 1 if gdename == "castelfranco veneto i"
replace nCH = 1 if gdename == "schwaig bei nuernberg d"
replace nCH = 1 if gdename == "chelles f"
replace nCH = 1 if gdename == "chichester, west sussex gb"
replace nCH = 1 if gdename == "crailsheim d"
replace nCH = 1 if gdename == "constantia za"
replace nCH = 1 if gdename == "schwaebisch hall"
replace nCH = 1 if gdename == "balboa usa"
replace nCH = 1 if gdename == "hall/tirol a"
replace nCH = 1 if gdename == "scottsdale (usa)"
replace nCH = 1 if gdename == "laveno i"
replace nCH = 1 if gdename == "beerzel b"
replace nCH = 1 if gdename == "san rafael (usa)"
replace nCH = 1 if gdename == "accra, ghana"
replace nCH = 1 if gdename == "umbertide i"
replace nCH = 1 if gdename == "wurmlingen d"
replace nCH = 1 if gdename == "bergen op zoom nl"
replace nCH = 1 if gdename == "isnu d"
replace nCH = 1 if gdename == "kirchheim unter teck d"
replace nCH = 1 if gdename == "saronno i"
replace nCH = 1 if gdename == "wolfsburg (d)"
replace nCH = 1 if gdename == "st julien en genevois"
replace nCH = 1 if gdename == "cambridgeshire gb"
replace nCH = 1 if gdename == "trollhaettan s"
replace nCH = 1 if gdename == "rungis f"
replace nCH = 1 if gdename == "ettenheim d"
replace nCH = 1 if gdename == "trezzano sul naviglio i"
replace nCH = 1 if gdename == "hattersheim (d)"
replace nCH = 1 if gdename == "volshsky (russland)"
replace nCH = 1 if gdename == "karlsbad d"
replace nCH = 1 if gdename == "wemmel/brabant b"
replace nCH = 1 if gdename == "erzingen d"
replace nCH = 1 if gdename == "weissensberg d"
replace nCH = 1 if gdename == "boxmeer nl"
replace nCH = 1 if gdename == "waldshut eschbach d"
replace nCH = 1 if gdename == "bad mergentheim d"
replace nCH = 1 if gdename == "hardenberg nl"
replace nCH = 1 if gdename == "middlesbrough cleveland gb"
replace nCH = 1 if gdename == "knutsford gb"
replace nCH = 1 if gdename == "old westbury usa"
replace nCH = 1 if gdename == "lagos (nigeria)"
replace nCH = 1 if gdename == "pforzsheim d"
replace nCH = 1 if gdename == "soedertaelje (s)"
replace nCH = 1 if gdename == "chatel f"
replace nCH = 1 if gdename == "pyrford, surrey (gb)"
replace nCH = 1 if gdename == "danderyd/s"
replace nCH = 1 if gdename == "naples usa"
replace nCH = 1 if gdename == "lidingoe (schweden)"
replace nCH = 1 if gdename == "hindhead gb"
replace nCH = 1 if gdename == "riemerling (d)"
replace nCH = 1 if gdename == "mailand it"
replace nCH = 1 if gdename == "sta agata li battiati i"
replace nCH = 1 if gdename == "ballaison f"
replace nCH = 1 if gdename == "treibach a"
replace nCH = 1 if gdename == "bonn bad godesberg d"
replace nCH = 1 if gdename == "bernin le village f"
replace nCH = 1 if gdename == "quartino/magadino"
replace nCH = 1 if gdename == "vienna, virginia usa"
replace nCH = 1 if gdename == "bromley gb"
replace nCH = 1 if gdename == "eching am ammersee (d)"
replace nCH = 1 if gdename == "olgiate olona i"
replace nCH = 1 if gdename == "port arthur tx usa"
replace nCH = 1 if gdename == "bailly (f)"
replace nCH = 1 if gdename == "le caire"
replace nCH = 1 if gdename == "menlo park, california usa"
replace nCH = 1 if gdename == "tunis t"
replace nCH = 1 if gdename == "klein st. paul (a)"
replace nCH = 1 if gdename == "mexiko"
replace nCH = 1 if gdename == "freetown, sierra leone"
replace nCH = 1 if gdename == "vivry/hte savoie f"
replace nCH = 1 if gdename == "shinjuku ku j"
replace nCH = 1 if gdename == "ventspils (lettland)"
replace nCH = 1 if gdename == "devon usa"
replace nCH = 1 if gdename == "ottobeuren (d)"
replace nCH = 1 if gdename == "ballycommon, tipperary irl"
replace nCH = 1 if gdename == "serroue"
replace nCH = 1 if gdename == "clencoe usa"
replace nCH = 1 if gdename == "aylesbury gb"
replace nCH = 1 if gdename == "albbruck d"
replace nCH = 1 if gdename == "embourg b"
replace nCH = 1 if gdename == "velin f"
replace nCH = 1 if gdename == "st. georges (bermuda)"
replace nCH = 1 if gdename == "breitenwang a"
replace nCH = 1 if gdename == "herkenbosch nl"
replace nCH = 1 if gdename == "lemwerder d"
replace nCH = 1 if gdename == "nischni novgorod rus"
replace nCH = 1 if gdename == "weddington usa"
replace nCH = 1 if gdename == "meer b"
replace nCH = 1 if gdename == "roethenbach allgaeu d"
replace nCH = 1 if gdename == "neustadt wied d"
replace nCH = 1 if gdename == "duissbergen sur mer b"
replace nCH = 1 if gdename == "bremen dl"
replace nCH = 1 if gdename == "geiselbach d"
replace nCH = 1 if gdename == "fuessen d"
replace nCH = 1 if gdename == "velden am woerthersee a"
replace nCH = 1 if gdename == "sannois f"
replace nCH = 1 if gdename == "hausen b. offenbach a. m (d)"
replace nCH = 1 if gdename == "toledo usa"
replace nCH = 1 if gdename == "bad saeckingen (d)"
replace nCH = 1 if gdename == "pelago i"
replace nCH = 1 if gdename == "l aigle f"
replace nCH = 1 if gdename == "sao paolo bra"
replace nCH = 1 if gdename == "dingy saint clair f"
replace nCH = 1 if gdename == "mont vermeil f"
replace nCH = 1 if gdename == "belluno i"
replace nCH = 1 if gdename == "new york city usa"
replace nCH = 1 if gdename == "bremen oberneuland d"
replace nCH = 1 if gdename == "marino rom i"
replace nCH = 1 if gdename == "lindlar schmitzhoehe d"
replace nCH = 1 if gdename == "mount hermon cal. usa"
replace nCH = 1 if gdename == "wilmette usa"
replace nCH = 1 if gdename == "las arenas e"
replace nCH = 1 if gdename == "villettes d anthon f"
replace nCH = 1 if gdename == "preussin f"
replace nCH = 1 if gdename == "frederiksberg kopenhagen dk"
replace nCH = 1 if gdename == "lottstetten (d)"
replace nCH = 1 if gdename == "bagno ripoli (i)"
replace nCH = 1 if gdename == "wormerveer nl"
replace nCH = 1 if gdename == "rotkreuz/rasch"
replace nCH = 1 if gdename == "bad urach d"
replace nCH = 1 if gdename == "etroubles i"
replace nCH = 1 if gdename == "pompano beach usa"
replace nCH = 1 if gdename == "radolfzell boehringen (brd)"
replace nCH = 1 if gdename == "floreal mus"
replace nCH = 1 if gdename == "tripolis lar"
replace nCH = 1 if gdename == "santa fe, new mexico usa"
replace nCH = 1 if gdename == "hammel dk"
replace nCH = 1 if gdename == "ewshot surrey gb"
replace nCH = 1 if gdename == "neuweiler d"
replace nCH = 1 if gdename == "esloev (s)"
replace nCH = 1 if gdename == "geldrop nl"
replace nCH = 1 if gdename == "engelkirchen d"
replace nCH = 1 if gdename == "lemineur, paris f"
replace nCH = 1 if gdename == "portland, oregon (usa)"
replace nCH = 1 if gdename == "suffolk gb"
replace nCH = 1 if gdename == "udine (i)"
replace nCH = 1 if gdename == "perrineville usa"
replace nCH = 1 if gdename == "sulzbach (d)"
replace nCH = 1 if gdename == "playa de aro e"
replace nCH = 1 if gdename == "middenbeemster nl"
replace nCH = 1 if gdename == "ohlsbach d"
replace nCH = 1 if gdename == "bailly/yvelines f"
replace nCH = 1 if gdename == "kamp lintfort a"
replace nCH = 1 if gdename == "lenox usa"
replace nCH = 1 if gdename == "brentwood (gb)"
replace nCH = 1 if gdename == "sleeuwijk, nl"
replace nCH = 1 if gdename == "hausen b. offenbach a.m (d)"
replace nCH = 1 if gdename == "san francisco, californien usa"
replace nCH = 1 if gdename == "st etienne f"
replace nCH = 1 if gdename == "upie f"
replace nCH = 1 if gdename == "zemst b"
replace nCH = 1 if gdename == "berkley usa"
replace nCH = 1 if gdename == "livigno i"

replace geo = -9 if nCH != 1
replace geo = -91 if nCH == 1
gen gdename_used = gdename
keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page gdename_used geo nCH
save "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo-9.dta", replace


* Append all sub-datasets
use "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo-99.dta", clear
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo1.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo21.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo22.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo3.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo41.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo42.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo43.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo91.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo92.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo93.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo-9.dta"

replace gdename = gdename_used if gdename_used != ""
rename geo geo_merge

bysort PID year ID_dupl: gen gde_options = _n  // indicator for observations with multiple options on gdename

keep PID year ID_dupl titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
geo_merge countPLZ4 gdenr_2018 gdename GdeNr_E_CNTR GdeNr_N_CNTR nCH gde_options options4a options4c options_manual

order PID year ID_dupl gde_options options4a options4c options_manual titlejob firstname lastname address address2 alternateaddress PLZ4 city companiesid functions page ///
geo_merge countPLZ4 gdename gdenr_2018 GdeNr_E_CNTR GdeNr_N_CNTR nCH 
 
label var nCH "foreign city" 
label var ID_dupl "Duplicates in original data (PID-year)"
label var options4a "potential gde-matches in 4a (after split)"
label var options4c "potential gde-matches in 4c (after split & gdename options w/o canton info)"
label var options_manual "potential gde-matches after manual coding"
label var gde_options "indicator for multiple options in gde"
label var geo_merge "basis info of merge"

/*
1  : PLZ match
21 : string match gdename & year with minimum of name standardization
22 : string match gdename & year with further standardization of name 
3  : extract PLZ from string and match on PLZ4
41 : no string match, split gdename, loop and merge on cell fragments (minimal standardization): manual checks
42 : no string match, split gdename, loop and merge on cell fragments (further standardization): manual checks
43 : no string match, split gdename, loop and merge on cell fragments (without cantonal information in gdenames): manual checks
91 : no matches on PLZ / gdename & year - manual correctoin of PLZ with subsequent merging on PLZ4
92 : no matches on PLZ / gdename & year - manual coding with subsequent merging on gdename
93 : no matches on PLZ / gdename & year - discard year: merge on gdename only
-9 : no merge after all procedures
-91: no merge because foreign city 
-99: no Geo-Info - no PLZ, no city info
*/

save "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person-Geo.dta", replace

tab year geo_merge

preserve
keep if geo_merge == -9 & nCH != 1
keep city gdename
gen count = -1
collapse (sum) count, by(city gdename)
order count gdename city
sort count city
save "$path\02_Processed_data\10_Directors_1934_2003\stillmissing.dta", replace
restore


erase "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_dupl0.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_dupl1.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_dupl2.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_dupl3.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_dupl4.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\unique_Gdenr2018-Geo.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\uPLZ_Geo.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo1_notmerged.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo2a_notmerged.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo2b_notmerged.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo3_notmerged.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo4_notmerged.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\reshaped_4rd-round.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo4abc.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo1.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo21.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo22.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo3.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo41.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo42.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo43.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo90.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo90a.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo90b.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo91.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo92.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo93.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo-99.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Persons-Geo-9.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\stillmissing.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\check_3.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\check_4a.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\check_4c.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\check_4c2.dta"



*erase "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp3.dta"
*erase "$path\02_Processed_data\10_Directors_1934_2003\comp_1934-2003_tmp3.dta"

