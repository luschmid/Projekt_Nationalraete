clear
cap log close
set more 1
version 17

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"


***** Add 297 missing observations from the east-coordinate error
* 1) Prepare GT2
* 2) Prepare subset of 297 obs
* 3) merge and code missing obs


*** Prepare "alternative Ground Truth" (GT2) for quality checks: use existing registery of mandates for elected National Councilors


** Read-in manually coded lists of mandates (Ladina Oester & Karin Gächter-Meile/Brunschweiler) --- > 1996 (no directory of national councilors)

foreach year in 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1997 1998 1998_1999 1999 2000 2001 2002 2003 {
	import excel "$path\01_Raw_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate Nationalraete.xlsx",  sheet(`year') firstrow clear
	gen year = "`year'"
	save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate Nationalraete_`year'.dta", replace
}	
use "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate Nationalraete_1983.dta", clear
foreach year in 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1997 1998 1998_1999 1999 2000 2001 2002 2003 { 
	append using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate Nationalraete_`year'.dta", force
}
foreach year in 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1997 1998 1998_1999 1999 2000 2001 2002 2003 {
 erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate Nationalraete_`year'.dta"
}
drop O P // empty cells
replace year = "1982" if year == "1983"
replace year = "1983" if year == "1984"
replace year = "1984" if year == "1985"
replace year = "1985" if year == "1986"
replace year = "1986" if year == "1987"
replace year = "1987" if year == "1988"
replace year = "1988" if year == "1989"
replace year = "1989" if year == "1990"
replace year = "1990" if year == "1991"
replace year = "1991" if year == "1992"
replace year = "1992" if year == "1993"
replace year = "1993" if year == "1994"
replace year = "1994" if year == "1995"
replace year = "1995" if year == "1996"  // no real changes made --> 1996 fehlt!
replace year = "1996" if year == "1997"
replace year = "1997" if year == "1998"
replace year = "1998" if year == "1998_1999"
destring year, replace
replace Party = "GP" if Lastname == "von Felten" & Firstname == "Margrith" & year == 1999  // von Felten left the SP in 1998 and joined a "links-grün-feministische Liste" (Wikipedia)
order year Lastname Firstname Title Party ZIPcode Placeofresidence Militaryrank Profession Politicaloffices Companyname Companylocation Signature Role Capital 

** Prepare data structure
gen company_norm = ""
replace company_norm = ustrtrim(Companyname)
replace company_norm = ustrlower(company_norm)
replace company_norm = usubinstr(company_norm, "ü", "ue",.)
replace company_norm = usubinstr(company_norm, "ä", "ae",.)
replace company_norm = usubinstr(company_norm, "ö", "oe",.)
replace company_norm = usubinstr(company_norm, "é", "e",.)
replace company_norm = usubinstr(company_norm, "è", "e",.)
replace company_norm = usubinstr(company_norm, "ê", "e",.)
replace company_norm = usubinstr(company_norm, "ç", "c",.)
replace company_norm = usubinstr(company_norm, `"""', "",.)
replace company_norm = usubinstr(company_norm, "«", "",.)
replace company_norm = usubinstr(company_norm, "»", "",.)
replace company_norm = usubinstr(company_norm, ".", "",.)
replace company_norm = usubinstr(company_norm, ",", "",.)
replace company_norm = usubinstr(company_norm, ";", "",.)
replace company_norm = usubinstr(company_norm, "-", "",.)
replace company_norm = usubinstr(company_norm, ":", "",.)
replace company_norm = usubinstr(company_norm, "'", "",.)
replace company_norm = usubinstr(company_norm, "  ", "",.)
gen polmand = .
replace polmand = 1 if (regexm(company_norm, "kommission")) | (regexm(company_norm, "commission"))
replace polmand = 1 if (regexm(company_norm, "delegation")) | (regexm(company_norm, "delegazion"))
gen stiftung = .
replace stiftung = 1 if (regexm(company_norm, "stiftung")) | (regexm(company_norm, "fondation")) | (regexm(company_norm, "fondazion"))
gen genossenschaft = .
replace genossenschaft = 1 if (regexm(company_norm, "genossenschaft")) | (regexm(company_norm, "cooperativ"))
gen verband = .
replace verband = 1 if (regexm(company_norm, "verband")) | (regexm(company_norm, "verbaende")) | (regexm(company_norm, "verein")) | (regexm(company_norm, "association")) | (regexm(company_norm, "associazione"))
foreach noAG in polmand stiftung genossenschaft verband {
	tab `noAG'
}
/*
// Check non-AG names
preserve
keep if polmand == 1 
gen n = 1
collapse (sum) n, by(Companyname)
order n Companyname
export excel using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Check_PolMandates.xlsx",  first(var) replace
restore
preserve
keep if stiftung == 1 
gen n = 1
collapse (sum) n, by(Companyname)
order n Companyname
export excel using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Check_Stiftung.xlsx",  first(var) replace
restore
preserve
keep if genossenschaft == 1 
gen n = 1
collapse (sum) n, by(Companyname)
order n Companyname
export excel using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Check_Genossenschaft.xlsx",  first(var) replace
restore
preserve
keep if verband == 1 
gen n = 1
collapse (sum) n, by(Companyname)
order n Companyname
export excel using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Check_Verband.xlsx",  first(var) replace
restore
*/

drop if stiftung == 1
drop if polmand == 1
gen genossenschaftkeep = .
replace genossenschaftkeep =  1 if Companyname == "BEG Bank Europäischer Genossenschaftsbanken"
replace genossenschaftkeep =  1 if Companyname == "BEG Bank Europäischer Genossenschaftsbanken."
replace genossenschaftkeep =  1 if Companyname == "Bank Europäischer Genossenschaftsbanken"
replace genossenschaftkeep =  1 if Companyname == "Genossenschaftliche Zentralbank"
replace genossenschaftkeep =  1 if Companyname == "Genossenschaftliche Zentralbank AG"
replace genossenschaftkeep =  1 if Companyname == "Genossenschaftliche Zentralbank AG."
replace genossenschaftkeep =  1 if Companyname == "Genossenschaftliche Zentralbank Aktiengesellschaft"
replace genossenschaftkeep =  1 if Companyname == "Internationale Genossenschaftsbank AG"
drop if genossenschaft == 1 & genossenschaftkeep != 1
gen verbandskeep = .
replace verbandskeep = 1 if Companyname == "Verbandsdruckerei AG Bern"
replace verbandskeep = 1 if Companyname == "Verbandsmolkerei Region Bern AG"
replace verbandskeep = 1 if Companyname == "Vereinigte Huttwil Bahnen (VHB)"
replace verbandskeep = 1 if Companyname == "Vereinigte Huttwil-Bahnen (VHB)."
replace verbandskeep = 1 if Companyname == "Vereinigte Schweizerische Rheinsalinen"
drop if verband == 1 & verbandskeep != 1

bysort year Lastname Firstname Party: gen nth_entry = _n // consecutively numbered entries
gen help1 = 1 if Companyname != "" & polmand != 1
bysort year Lastname Firstname Party: egen nbrmand = total(help1)
drop help1

// corrections to make unique entries on name firstname year
replace Party = "SVP" if Party == "" & year == 1986 & Lastname == "Fischer" & Firstname == "Theo"
replace Party = "SVP" if Party == "" & year == 1987 & Lastname == "Fischer" & Firstname == "Theo"
replace Lastname = "Fischer Sursee" if year >= 1983 & year <= 1994 & Lastname == "Fischer" & Firstname == "Theo" & Party == "CVP"  // "official" name in NC to distinguish from Fischer-Hägglingen (Wikipedia)
replace Lastname = "Fischer Hägglingen" if year >= 1982 & year <= 1996 & Lastname == "Fischer" & Firstname == "Theo" & Party == "SVP" // "official" name in NC to distinguish from Fischer-Sursee (Wikipedia)
replace Lastname = usubinstr(Lastname, "-", " ",.)
replace Firstname = usubinstr(Firstname, "-", " ",.)

// corrections to match names in NR data
replace Lastname = "Aeppli" if Lastname == "Aeppli Wartmann" & Firstname == "Regina"  // in NR data as "Regine Aeppli" (checked with Wikipedia)
replace Firstname = "Regine" if Lastname == "Aeppli" & Firstname == "Regina"
replace Lastname = "Ammann" if Lastname == "Ammann Schellenberg" & Firstname == "Ulrich" 
replace Firstname = "Manfred" if Lastname == "Aregger" & Firstname == "Manfred sen."
replace Lastname = "Aubry Moine" if Lastname == "Aubry" & Firstname == "Geneviève" & (year <= 1983 | year > 1987)
replace Lastname = "Aguet" if Lastname == "Auget" & Firstname == "Pierre" 
replace Lastname = "Baader" if Lastname == "Baader Buri" & Firstname == "Caspar" 
replace Firstname = "Alma" if Lastname == "Bacciarini" & Firstname == "Ama"
replace Firstname = "Käthi" if Lastname == "Bangerter" & Firstname == "Katharina"
replace Lastname = "Bangerter" if Lastname == "Bangerter Schober" & Firstname == "Katharina" 
replace Firstname = "Käthi" if Lastname == "Bangerter" & Firstname == "Katharina"
replace Firstname = "J. Alexander" if Lastname == "Baumann" & Firstname == "Joseph Alexander"
replace Lastname = "Baumann" if Lastname == "Baumann Bieri" & Firstname == "Ruedi" 

// keep information on political offices (do we need this?)
preserve
keep year Lastname Firstname Party Politicaloffices nth_entry
drop if Politicaloffices == "" & nth_entry > 1
bysort year Lastname Firstname Party: gen PolOffices = Politicaloffices[1] 
by year Lastname Firstname Party: replace PolOffices = PolOffices[_n-1] + "; " + Politicaloffices if _n > 1 
by year Lastname Firstname Party: replace PolOffices = PolOffices[_N]
bysort year Lastname Firstname Party: gen n = _n
keep if n == 1
drop Politicaloffices nth_entry n
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\PolOffices.dta", replace
restore

merge m:1 year Lastname Firstname Party using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\PolOffices.dta"
drop _merge
drop Politicaloffices

// keep only company information (and 1 entry for those w/o mandates) --> drop infos pertaining to political specializations (e.g. Staatspolitik, Rechtsfragen,...)
drop if Companyname == "" & nbrmand != 0
keep if nbrmand == 0 | Companylocation != "" | Signature != "" | Role != "" | Capital != .
drop if Companyname == "" & nth_entry>1
replace Companyname = "none" if Companyname == "" 
replace company_norm = "none" if Companyname == "none"
sort year Lastname Firstname Party

duplicates tag year Lastname Firstname Party Companyname Companylocation, gen(dupl)
*br year Lastname Firstname Party Companyname Companylocation Signature Role Capital if dupl>0
//make sure no info is lost (sometimes different info exist in duplicate entries)
bysort year Lastname Firstname Party Companyname: gen signature = Signature[1] 
by year Lastname Firstname Party Companyname: replace signature = signature[_n-1] + "; " + Signature if _n > 1 
by year Lastname Firstname Party Companyname: replace signature = signature[_N]
bysort year Lastname Firstname Party Companyname: gen role = Role[1] 
by year Lastname Firstname Party Companyname: replace role = role[_n-1] + "; " + Role if _n > 1 
by year Lastname Firstname Party Companyname: replace role = role[_N]
tostring Capital, replace
bysort year Lastname Firstname Party Companyname: gen capital = Capital[1] 
by year Lastname Firstname Party Companyname: replace capital = capital[_n-1] + "; " + Capital if _n > 1 
by year Lastname Firstname Party Companyname: replace capital = capital[_N]
bysort year Lastname Firstname Party Companyname: gen n = _n
keep if n == 1 
drop Signature Role Capital dupl n

// are firmnames unique?
duplicates report year Lastname Firstname Party Companyname
duplicates report year Lastname Firstname Party company_norm
duplicates report year Lastname Firstname Party company_norm Companylocation
duplicates tag year Lastname Firstname Party company_norm, gen(dupl1) // 10 firmnames not unique
duplicates tag year Lastname Firstname Party company_norm Companylocation, gen(dupl2) // 3 firmnames not unique

/*
br year Lastname Firstname Party company_norm Companyname Companylocation if dupl1>0
year	Lastname	Firstname	Party	company_norm	Companyname	Companylocation
1984	Röthlin Lieb	Walter	CVP	reinhard a. ag	Reinhard A. AG	Luzern
1984	Röthlin Lieb	Walter	CVP	reinhard a. ag	Reinhard, A., AG	Kerns
1986	Fischer Sursee	Theo	CVP	luzerner kantonalbank	Luzerner Kantonalbank	
1986	Fischer Sursee	Theo	CVP	luzerner kantonalbank	Luzerner Kantonalbank.	Luzern
1987	Fischer Hägglingen	Theo	SVP	realit treuhand ag	Realit Treuhand AG	Lenzburg
1987	Fischer Hägglingen	Theo	SVP	realit treuhand ag	Realit Treuhand AG.	Lenzburg
1988	Fischer Hägglingen	Theo	SVP	realit treuhand ag	Realit Treuhand AG	Lenzburg
1988	Fischer Hägglingen	Theo	SVP	realit treuhand ag	Realit Treuhand AG.	Lenzburg
1989	Fischer Hägglingen	Theo	SVP	realit treuhand ag	Realit Treuhand AG	Lenzburg
1989	Fischer Hägglingen	Theo	SVP	realit treuhand ag	Realit Treuhand AG.	Lenzburg
1993	Spoerry Toneatti	Verena	FDP	schweizerische kreditanstalt	Schweizerische Kreditanstalt	
1993	Spoerry Toneatti	Verena	FDP	schweizerische kreditanstalt	Schweizerische Kreditanstalt.	Zürich
1997	Engler	Rolf	CVP	swica gesundheitsorganisation	SWICA Gesundheitsorganisation	
1997	Engler	Rolf	CVP	swica gesundheitsorganisation	Swica Gesundheitsorganisation	Winterthur
1997	Lachat	François	CVP	rassemblement jurassien	Rassemblement jurassien	
1997	Lachat	François	CVP	rassemblement jurassien	Rassemblement jurassien.	Delémont
1998	Engler	Rolf	CVP	swica gesundheitsorganisation	SWICA Gesundheitsorganisation	
1998	Engler	Rolf	CVP	swica gesundheitsorganisation	Swica Gesundheitsorganisation	Winterthur
1998	Lachat	François	CVP	rassemblement jurassien	Rassemblement jurassien	
1998	Lachat	François	CVP	rassemblement jurassien	Rassemblement jurassien.	Delémont
*/
drop dupl1 dupl2 
duplicates drop year Lastname Firstname Party company_norm, force

drop if year == 2003 & Lastname == "" & Firstname == "" & Party == ""
rename Lastname name
rename Firstname firstname
sort name firstname year
order year name firstname Title Party ZIPcode Placeofresidence Militaryrank Profession PolOffices Companyname Companylocation signature role capital
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete.dta", replace
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\PolOffices.dta"


** Geo-coding Mandate dataset 

* Standardize municipality names in geo dataset **
use "$path\02_Processed_data\08_Municipalities\GdeNr2018_Geo_1931-2018.dta", clear
keep year ctn gdename gdenr gdenr_2018 E_CNTR N_CNTR
drop if year==. & gdename==""
duplicates report gdename year
duplicates drop gdename year, force // drop 1 observation
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

* Standardize municipality names in Mandate dataset **
use "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete.dta", clear
gen gdename = Placeofresidence
replace gdename = "missing" if gdename == ""

preserve
keep if gdename == "missing" & ZIPcode == .
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo0.dta", replace 
restore
drop if gdename == "missing"
replace gdename = usubinstr(gdename, "`", " ",.)
replace gdename = usubinstr(gdename, "-", " ",.)
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
sort gdename year

merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2

// matched sample from first round
preserve 
keep if _merge==3
gen geomerge = 1
drop _merge
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo1.dta", replace
restore

drop if _merge == 3
drop _merge gdenr gdenr_2018 ctn E_CNTR N_CNTR

preserve
use "$path\02_Processed_data\08_Municipalities\uniquePLZ_(Avg)GdeNr2018-Geo.dta", clear
duplicates report PLZ4 
sort PLZ4 
save "$path\02_Processed_data\10_Directors_1934_2003\uPLZ_Geo.dta", replace
restore

// merge on ZIP codes
rename ZIPcode PLZ4
sort PLZ4
merge m:1 PLZ4 using "$path\02_Processed_data\10_Directors_1934_2003\uPLZ_Geo.dta"  // merge on PLZ
drop if _merge == 2

// matched sample from second round (PLZ)
preserve 
keep if _merge==3
gen geomerge = 2
drop _merge
rename GdeNr_N_CNTR N_CNTR
rename GdeNr_E_CNTR E_CNTR
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo2.dta", replace
restore

drop if _merge == 3
drop _merge gdenr gdenr_2018 GdeNr_N_CNTR GdeNr_E_CNTR

// manual check remaining Place of residences
preserve
gen n = -1
collapse (sum) n, by(gdename)
sort n
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_nomatch_gdename-PLZ.dta", replace 
restore

// clear cases with only one option
replace gdename = "aesch bei birmensdorf" if gdename == "aesch b. birmensdorf"
replace gdename = "fluehli" if gdename == "soerenberg"
replace gdename = "ludiano" if gdename == "malvaglia chiesa"
replace gdename = "courrendlin" if gdename == "courredlin"
replace gdename = "wisen (so)" if gdename == "wissen"
replace gdename = "dompierre (vd)" if gdename == "dompierre"  // this is the case of Jean Pierre Berger
replace gdename = "lully (vd)" if gdename == "lully"  // this is the case of Charles Friderici

sort gdename year
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2

keep if _merge==3
gen geomerge = 3
drop _merge
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo3.dta", replace

use "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo0.dta", clear
forv i = 1(1)3 {
append using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo`i'.dta"
}

replace gdenr_2018 = 99999 if gdename == "missing"
replace E_CNTR = 9999999 if gdename == "missing"
replace N_CNTR = 9999999 if gdename == "missing"

foreach var in gdename gdenr ctn E_CNTR N_CNTR geomerge {  // gdenr_2018
	rename `var' `var'_mand
}

sort year name firstname Party Companyname
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo.dta", replace

forv i = 0(1)3 {
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo`i'.dta"
}
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_nomatch_gdename-PLZ.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\uPLZ_Geo.dta"


** Merge with National Council dataset

use "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta", clear
encode ID, gen(IDno)
xtset IDno year
keep if year>=1979 & year<=2003 // sample restriction larger than mandate data: NR data are election years!
rename gdenr_2018_w gdenr_2018
rename E_CNTR_w  E_CNTR
rename N_CNTR_w N_CNTR
replace name = usubinstr(name, "-", " ",.)
replace firstname = usubinstr(firstname, "-", " ",.)
replace name = "Fischer Sursee" if year >= 1979 & year <= 1994 & name == "Fischer" & firstname == "Theo" & listname_bfs == "CVP/PDC"
replace name = "Fischer Hägglingen" if year >= 1979 & year <= 1995 & name == "Fischer" & firstname == "Theo" & listname_bfs == "SVP/UDC"

// check for duplicates in merge variables (first/last names): none
preserve
duplicates report year name firstname gdenr_2018
duplicates tag year name firstname gdenr_2018, gen(dupl)
restore

gen eleyear = year
expand 4, gen(tmp_year)
bysort ID year: gen n=_n  // the first year of office is the one following the election year (elections in October): VR Verzeichnisse typically finished by eletion time
replace year = year + n
drop tmp_year n
keep if year>=1982 & year<=2003
sort name firstname year

save "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo_tmp.dta", replace

merge 1:m year name firstname gdenr_2018 using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo.dta"


// export matched NR-Mandate data
preserve
keep if _merge == 3
drop _merge
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NR_mandate_match.dta", replace
restore

// export non-matched from mandate dataset
preserve
keep if _merge == 2
drop _merge
keep year name firstname gdenr_2018
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\non-match_mandates.dta", replace
restore

// keep non-matched from both datasets (drop _merge == 3)
preserve
drop if _merge == 3
bysort name firstname year _merge: gen n = _n
keep if n == 1
rename _merge no_merges
sort name firstname year no_merges
keep ID year name firstname Party gdename* gdenr* gdenr_2018* ctn* E_CNTR N_CNTR E_CNTR_mand N_CNTR_mand no_merges
order ID year name firstname no_merges
/*
export excel using  "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\unique_name-year-NR-mandate.xlsx", first(var) replace
//  --- > inactive as manual coding is competeted
*/
restore

// manual corrections
use "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo.dta", clear
gen ID = ""
replace ID = "ZH-1991-0004" if ID == "" & year == 1996 & name == "Aeppli" & firstname == "Regine" & gdenr_2018 == 99999
replace ID = "ZH-1991-0004" if ID == "" & year == 1997 & name == "Aeppli" & firstname == "Regine" & gdenr_2018 == 99999
replace ID = "ZH-1991-0004" if ID == "" & year == 1998 & name == "Aeppli" & firstname == "Regine" & gdenr_2018 == 99999
replace ID = "SG-1995-0003" if ID == "" & year == 1996 & name == "Alder" & firstname == "Fredi" & gdenr_2018 == 99999
replace ID = "GR-1983-0001" if ID == "" & year == 1983 & name == "Aliesch" & firstname == "Peter" & gdenr_2018 == 3901
replace ID = "SG-1975-0002" if ID == "" & year == 1985 & name == "Ammann" & firstname == "Walter" & gdenr_2018 == 99999
replace ID = "VS-1995-0002" if ID == "" & year == 2002 & name == "Antille" & firstname == "Charles Albert" & gdenr_2018 == 99999
replace ID = "VS-1995-0002" if ID == "" & year == 2003 & name == "Antille" & firstname == "Charles Albert" & gdenr_2018 == 99999
replace ID = "BEJU-1975-0017" if ID == "" & year == 1996 & name == "Aubry Moine" & firstname == "Geneviève" & gdenr_2018 == 713
replace ID = "SO-1995-0001" if ID == "" & year == 2002 & name == "Bader" & firstname == "Elvira" & gdenr_2018 == 99999
replace ID = "SO-1995-0001" if ID == "" & year == 2003 & name == "Bader" & firstname == "Elvira" & gdenr_2018 == 99999
replace ID = "FR-1955-0002" if ID == "" & year == 1982 & name == "Barras" & firstname == "Louis" & gdenr_2018 == 2175
replace ID = "ZH-1975-0029" if ID == "" & year == 1991 & name == "Baumberger" & firstname == "Peter" & gdenr_2018 == 230
replace ID = "VD-1991-0005" if ID == "" & year == 2000 & name == "Beck" & firstname == "Serge" & gdenr_2018 == 99999
replace ID = "VD-1991-0005" if ID == "" & year == 2001 & name == "Beck" & firstname == "Serge" & gdenr_2018 == 99999
replace ID = "VD-1991-0005" if ID == "" & year == 2002 & name == "Beck" & firstname == "Serge" & gdenr_2018 == 99999
replace ID = "VD-1991-0005" if ID == "" & year == 2003 & name == "Beck" & firstname == "Serge" & gdenr_2018 == 99999
replace ID = "ZH-1967-0034" if ID == "" & year == 1982 & name == "Biel" & firstname == "Walter" & gdenr_2018 == 96
replace ID = "ZH-1967-0034" if ID == "" & year == 1983 & name == "Biel" & firstname == "Walter" & gdenr_2018 == 96
replace ID = "SG-1999-0011" if ID == "" & year == 2000 & name == "Biger" & firstname == "Elmar" & gdenr_2018 == 3297
replace ID = "SG-1999-0011" if ID == "" & year == 2002 & name == "Bigger" & firstname == "Elmar" & gdenr_2018 == 99999
replace ID = "SG-1999-0011" if ID == "" & year == 2003 & name == "Bigger" & firstname == "Elmar" & gdenr_2018 == 99999
replace ID = "ZH-1991-0051" if ID == "" & year == 1991 & name == "Binder Gäumann" & firstname == "Max" & gdenr_2018 == 296
replace ID = "ZH-1991-0051" if ID == "" & year == 1996 & name == "Binder Gäumann" & firstname == "Max" & gdenr_2018 == 296
replace ID = "ZH-1991-0051" if ID == "" & year == 1997 & name == "Binder Gäumann" & firstname == "Max" & gdenr_2018 == 296
replace ID = "ZH-1991-0051" if ID == "" & year == 1998 & name == "Binder Gäumann" & firstname == "Max" & gdenr_2018 == 296
replace ID = "ZH-1991-0055" if ID == "" & year == 1991 & name == "Bischof" & firstname == "Hardi" & gdenr_2018 == 261
replace ID = "ZH-1991-0055" if ID == "" & year == 1993 & name == "Bischof" & firstname == "Leonhard" & gdenr_2018 == 261
replace ID = "ZH-1991-0055" if ID == "" & year == 1994 & name == "Bischof" & firstname == "Leonhard" & gdenr_2018 == 261
replace ID = "OW-1987-0002" if ID == "" & year == 1987 & name == "Blatter" & firstname == "Hans Ueli" & gdenr_2018 == 1402
replace ID = "OW-1987-0002" if ID == "" & year == 1988 & name == "Blatter" & firstname == "Hans Ueli" & gdenr_2018 == 1402
replace ID = "ZH-1979-0051" if ID == "" & year == 1998 & name == "Blocher" & firstname == "Christoph" & gdenr_2018 == 152
replace ID = "ZH-1979-0051" if ID == "" & year == 1999 & name == "Blocher" & firstname == "Christoph" & gdenr_2018 == 152
replace ID = "SZ-1971-0001" if ID == "" & year == 1982 & name == "Blunschy" & firstname == "Elisabeth" & gdenr_2018 == 1372
replace ID = "SZ-1971-0001" if ID == "" & year == 1983 & name == "Blunschy" & firstname == "Elisabeth" & gdenr_2018 == 1372
replace ID = "SZ-1971-0001" if ID == "" & year == 1984 & name == "Blunschy" & firstname == "Elisabeth" & gdenr_2018 == 1372
replace ID = "SZ-1971-0001" if ID == "" & year == 1985 & name == "Blunschy" & firstname == "Elisabeth" & gdenr_2018 == 1372
replace ID = "SZ-1971-0001" if ID == "" & year == 1986 & name == "Blunschy" & firstname == "Elisabeth" & gdenr_2018 == 1372
replace ID = "BEJU-1983-0037" if ID == "" & year == 1983 & name == "Bonny" & firstname == "Jean Pierre" & gdenr_2018 == 353
replace ID = "VS-1983-0005" if ID == "" & year == 1987 & name == "Bonvin" & firstname == "Hubert" & gdenr_2018 == 6240
replace ID = "VS-1983-0005" if ID == "" & year == 1988 & name == "Bonvin" & firstname == "Hubert" & gdenr_2018 == 6240
replace ID = "NE-1979-0004" if ID == "" & year == 1982 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "NE-1979-0004" if ID == "" & year == 1983 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "NE-1979-0004" if ID == "" & year == 1988 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "NE-1979-0004" if ID == "" & year == 1989 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "NE-1979-0004" if ID == "" & year == 1990 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "NE-1979-0004" if ID == "" & year == 1991 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "NE-1979-0004" if ID == "" & year == 1992 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "NE-1979-0004" if ID == "" & year == 1993 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "NE-1979-0004" if ID == "" & year == 1994 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "NE-1979-0004" if ID == "" & year == 1996 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "NE-1979-0004" if ID == "" & year == 1997 & name == "Borel" & firstname == "François" & gdenr_2018 == 6407
replace ID = "SO-1991-0007" if ID == "" & year == 1991 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 1992 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 1993 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 1994 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 1996 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 1997 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 1998 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 1999 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 2000 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 2001 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 2002 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "SO-1991-0007" if ID == "" & year == 2003 & name == "Borer" & firstname == "Roland" & gdenr_2018 == 2403
replace ID = "TI-1991-9053" if ID == "" & year == 1991 & name == "Borradori" & firstname == "Marco" & gdenr_2018 == 5192
replace ID = "ZH-1987-0077" if ID == "" & year == 1991 & name == "Bortoluzzi" & firstname == "Anton" & gdenr_2018 == 2
replace ID = "ZH-1987-0077" if ID == "" & year == 1992 & name == "Bortoluzzi" & firstname == "Anton" & gdenr_2018 == 2
replace ID = "ZH-1987-0077" if ID == "" & year == 1993 & name == "Bortoluzzi" & firstname == "Anton" & gdenr_2018 == 2
replace ID = "ZH-1967-0046" if ID == "" & year == 1984 & name == "Braunschweig" & firstname == "Heinz" & gdenr_2018 == 191
replace ID = "AG-1999-0028" if ID == "" & year == 2002 & name == "Bruderer" & firstname == "Pascale" & gdenr_2018 == 99999
replace ID = "AG-1999-0028" if ID == "" & year == 2003 & name == "Bruderer" & firstname == "Pascale" & gdenr_2018 == 99999
replace ID = "SG-1995-0026" if ID == "" & year == 2001 & name == "Brunner" & firstname == "Toni" & gdenr_2018 == 99999
replace ID = "SG-1995-0026" if ID == "" & year == 2002 & name == "Brunner" & firstname == "Toni" & gdenr_2018 == 99999
replace ID = "SG-1995-0026" if ID == "" & year == 2003 & name == "Brunner" & firstname == "Toni" & gdenr_2018 == 99999
replace ID = "VD-1991-0015" if ID == "" & year == 2002 & name == "Bugnon Maillard" & firstname == "André" & gdenr_2018 == 5646
replace ID = "VD-1991-0015" if ID == "" & year == 2003 & name == "Bugnon Maillard" & firstname == "André" & gdenr_2018 == 5646
replace ID = "BS-1987-0013" if ID == "" & year == 1987 & name == "Burckhardt" & firstname == "Martin H." & gdenr_2018 == 2701
replace ID = "BEJU-1987-0065" if ID == "" & year == 1987 & name == "Bär Schwab" & firstname == "Rosmarie" & gdenr_2018 == 356
replace ID = "BEJU-1987-0071" if ID == "" & year == 1987 & name == "Bäumlin" & firstname == "Ursula" & gdenr_2018 == 351
replace ID = "LU-1991-0005" if ID == "" & year == 1991 & name == "Bühlmann" & firstname == "Cécile" & gdenr_2018 == 1061
replace ID = "SH-1991-0001" if ID == "" & year == 1991 & name == "Bührer" & firstname == "Gerold" & gdenr_2018 == 2920
replace ID = "SG-1971-0010" if ID == "" & year == 1982 & name == "Bürer Wildhaber" & firstname == "Kurt" & gdenr_2018 == 3298
replace ID = "SG-1971-0010" if ID == "" & year == 1983 & name == "Bürer Wildhaber" & firstname == "Kurt" & gdenr_2018 == 3298
replace ID = "SG-1971-0010" if ID == "" & year == 1984 & name == "Bürer Wildhaber" & firstname == "Kurt" & gdenr_2018 == 3298
replace ID = "SG-1971-0010" if ID == "" & year == 1985 & name == "Bürer Wildhaber" & firstname == "Kurt" & gdenr_2018 == 3298
replace ID = "SZ-1987-0004" if ID == "" & year == 1987 & name == "Bürgi" & firstname == "Jakob" & gdenr_2018 == 1321
replace ID = "TI-1971-0007" if ID == "" & year == 1987 & name == "Caccia" & firstname == "Fulvio" & gdenr_2018 == 5003
replace ID = "TG-1983-9032" if ID == "" & year == 1986 & name == "Camenzind Wüest" & firstname == "Margrit" & gdenr_2018 == 4566
replace ID = "TI-1975-0008" if ID == "" & year == 1991 & name == "Camponovo" & firstname == "Geo" & gdenr_2018 == 5250
replace ID = "SG-1991-0016" if ID == "" & year == 1991 & name == "Caspar Hutter" & firstname == "Elisabeth" & gdenr_2018 == 3204
replace ID = "SG-1991-0016" if ID == "" & year == 1992 & name == "Caspar Hutter" & firstname == "Elisabeth" & gdenr_2018 == 3204
replace ID = "SG-1991-0016" if ID == "" & year == 1993 & name == "Caspar Hutter" & firstname == "Elisabeth" & gdenr_2018 == 3204
replace ID = "SG-1991-0016" if ID == "" & year == 1994 & name == "Caspar Hutter" & firstname == "Elisabeth" & gdenr_2018 == 3204
replace ID = "NE-1975-0005" if ID == "" & year == 1982 & name == "Cavadini" & firstname == "Jean" & gdenr_2018 == 6454
replace ID = "NE-1975-0005" if ID == "" & year == 1983 & name == "Cavadini" & firstname == "Jean" & gdenr_2018 == 6454
replace ID = "TI-1971-0014" if ID == "" & year == 1999 & name == "Cavalli" & firstname == "Francesco" & gdenr_2018 == 5091
replace ID = "FR-1991-0010" if ID == "" & year == 2000 & name == "Chapuis" & firstname == "Liliane" & gdenr_2018 == 99999
replace ID = "FR-1991-0010" if ID == "" & year == 2001 & name == "Chapuis" & firstname == "Liliane" & gdenr_2018 == 99999
replace ID = "FR-1991-0010" if ID == "" & year == 2002 & name == "Chapuis" & firstname == "Liliane" & gdenr_2018 == 99999
replace ID = "FR-1991-0010" if ID == "" & year == 2003 & name == "Chapuis" & firstname == "Liliane" & gdenr_2018 == 99999
replace ID = "VD-1991-0024" if ID == "" & year == 1991 & name == "Chevallaz" & firstname == "Olivier" & gdenr_2018 == 5586
replace ID = "AG-1963-9073" if ID == "" & year == 1982 & name == "Chopard" & firstname == " Max" & gdenr_2018 == 4044
replace ID = "AG-1963-9073" if ID == "" & year == 1983 & name == "Chopard" & firstname == " Max" & gdenr_2018 == 4044
replace ID = "AG-1963-9073" if ID == "" & year == 1984 & name == "Chopard" & firstname == " Max" & gdenr_2018 == 4044
replace ID = "AG-1963-9073" if ID == "" & year == 1985 & name == "Chopard" & firstname == " Max" & gdenr_2018 == 4044
replace ID = "AG-1963-9073" if ID == "" & year == 1986 & name == "Chopard" & firstname == " Max" & gdenr_2018 == 4044
replace ID = "GE-1975-0008" if ID == "" & year == 1985 & name == "Christinat" & firstname == "Amélia" & gdenr_2018 == 99999
replace ID = "BEJU-1983-0069" if ID == "" & year == 1983 & name == "Clivaz" & firstname == "Jean" & gdenr_2018 == 616
replace ID = "VS-1971-0006" if ID == "" & year == 1998 & name == "Comby" & firstname == "Bernard J." & gdenr_2018 == 6141
replace ID = "VS-1971-0006" if ID == "" & year == 1999 & name == "Comby" & firstname == "Bernard J." & gdenr_2018 == 6141
replace ID = "VS-1971-0006" if ID == "" & year == 1991 & name == "Comby" & firstname == "Bernard Jules" & gdenr_2018 == 6141
replace ID = "VS-1971-0006" if ID == "" & year == 1992 & name == "Comby" & firstname == "Bernard Jules" & gdenr_2018 == 6141
replace ID = "VS-1971-0006" if ID == "" & year == 1993 & name == "Comby" & firstname == "Bernard Jules" & gdenr_2018 == 6141
replace ID = "VS-1971-0006" if ID == "" & year == 1994 & name == "Comby" & firstname == "Bernard Jules" & gdenr_2018 == 6141
replace ID = "VS-1971-0006" if ID == "" & year == 1996 & name == "Comby" & firstname == "Bernard Jules" & gdenr_2018 == 6141
replace ID = "VS-1971-0006" if ID == "" & year == 1997 & name == "Comby" & firstname == "Bernard Jules" & gdenr_2018 == 6141
replace ID = "FR-1963-0004" if ID == "" & year == 1983 & name == "Cottet" & firstname == "Joseph" & gdenr_2018 == 2323
replace ID = "TI-1983-9023" if ID == "" & year == 1983 & name == "Cotti" & firstname == "Flavio" & gdenr_2018 == 5113
replace ID = "NE-1987-0008" if ID == "" & year == 2000 & name == "Cuche" & firstname == "Fernand" & gdenr_2018 == 99999
replace ID = "NE-1987-0008" if ID == "" & year == 2001 & name == "Cuche" & firstname == "Fernand" & gdenr_2018 == 99999
replace ID = "NE-1987-0008" if ID == "" & year == 2002 & name == "Cuche" & firstname == "Fernand" & gdenr_2018 == 99999
replace ID = "NE-1987-0008" if ID == "" & year == 2003 & name == "Cuche" & firstname == "Fernand" & gdenr_2018 == 99999
replace ID = "VS-1975-9021" if ID == "" & year == 1982 & name == "Darbellay" & firstname == "Vital" & gdenr_2018 == 6136
replace ID = "VS-1975-9021" if ID == "" & year == 1983 & name == "Darbellay" & firstname == "Vital" & gdenr_2018 == 6136
replace ID = "NE-1971-0008" if ID == "" & year == 1985 & name == "Deneys" & firstname == "Heidi" & gdenr_2018 == 99999
replace ID = "SZ-1991-0005" if ID == "" & year == 1991 & name == "Dettling" & firstname == "Toni" & gdenr_2018 == 1372
replace ID = "ZH-1987-9190" if ID == "" & year == 1987 & name == "Diener" & firstname == "Verena" & gdenr_2018 == 24
replace ID = "ZH-1987-9190" if ID == "" & year == 1992 & name == "Diener Aeppli" & firstname == "Verena" & gdenr_2018 == 24
replace ID = "ZH-1987-9190" if ID == "" & year == 1993 & name == "Diener Aeppli" & firstname == "Verena" & gdenr_2018 == 24
replace ID = "ZH-1987-9190" if ID == "" & year == 1994 & name == "Diener Aeppli" & firstname == "Verena" & gdenr_2018 == 24
replace ID = "ZH-1987-9190" if ID == "" & year == 1996 & name == "Diener Aeppli" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1987-9190" if ID == "" & year == 1997 & name == "Diener Aeppli" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "BEJU-1987-0105" if ID == "" & year == 1987 & name == "Dietrich" & firstname == "Franz" & gdenr_2018 == 355
replace ID = "TI-1995-0027" if ID == "" & year == 1999 & name == "Donato" & firstname == "Franco" & gdenr_2018 == 5396
replace ID = "BEJU-1995-0109" if ID == "" & year == 2001 & name == "Donzé" & firstname == "Walter" & gdenr_2018 == 99999
replace ID = "BEJU-1995-0109" if ID == "" & year == 2002 & name == "Donzé" & firstname == "Walter" & gdenr_2018 == 99999
replace ID = "BEJU-1995-0109" if ID == "" & year == 2003 & name == "Donzé" & firstname == "Walter" & gdenr_2018 == 99999
replace ID = "LU-1987-0014" if ID == "" & year == 1987 & name == "Dormann" & firstname == "Rosmarie" & gdenr_2018 == 1040
replace ID = "ZH-1979-0115" if ID == "" & year == 1987 & name == "Dreher" & firstname == "Michael E." & gdenr_2018 == 154
replace ID = "VD-1983-0057" if ID == "" & year == 1983 & name == "Dubois" & firstname == "Marcel" & gdenr_2018 == 5607
replace ID = "GE-1975-0015" if ID == "" & year == 1987 & name == "Ducret" & firstname == "Dominique" & gdenr_2018 == 6621
replace ID = "BS-1999-0009" if ID == "" & year == 2000 & name == "Dunand" & firstname == "Jean Henri" & gdenr_2018 == 99999
replace ID = "BS-1999-0009" if ID == "" & year == 2001 & name == "Dunant" & firstname == "Jean Henri" & gdenr_2018 == 99999
replace ID = "OW-1995-0002" if ID == "" & year == 2000 & name == "Durrer" & firstname == "Adalbert" & gdenr_2018 == 1401
replace ID = "OW-1995-0002" if ID == "" & year == 2001 & name == "Durrer" & firstname == "Adalbert" & gdenr_2018 == 1401
replace ID = "VD-1979-0045" if ID == "" & year == 1991 & name == "Duvoisin" & firstname == "Pierre André" & gdenr_2018 == 5938
replace ID = "VD-1979-0045" if ID == "" & year == 1992 & name == "Duvoisin" & firstname == "Pierre André" & gdenr_2018 == 5938
replace ID = "VD-1979-0045" if ID == "" & year == 1993 & name == "Duvoisin" & firstname == "Pierre André" & gdenr_2018 == 5938
replace ID = "VD-1979-0045" if ID == "" & year == 1994 & name == "Duvoisin" & firstname == "Pierre André" & gdenr_2018 == 5938
replace ID = "SZ-1995-0006" if ID == "" & year == 1996 & name == "Eberhard" & firstname == "Toni" & gdenr_2018 == 1331
replace ID = "SZ-1995-0006" if ID == "" & year == 1997 & name == "Eberhard" & firstname == "Toni" & gdenr_2018 == 1331
replace ID = "SZ-1995-0006" if ID == "" & year == 1998 & name == "Eberhard" & firstname == "Toni" & gdenr_2018 == 1331
replace ID = "SZ-1995-0006" if ID == "" & year == 1999 & name == "Eberhard" & firstname == "Toni" & gdenr_2018 == 1331
replace ID = "BEJU-1979-0069" if ID == "" & year == 1992 & name == "Eggenberger" & firstname == "Georges" & gdenr_2018 == 363
replace ID = "BEJU-1979-0069" if ID == "" & year == 1993 & name == "Eggenberger" & firstname == "Georges" & gdenr_2018 == 363
replace ID = "BEJU-1979-0069" if ID == "" & year == 1994 & name == "Eggenberger" & firstname == "Georges" & gdenr_2018 == 363
replace ID = "GE-1975-9001" if ID == "" & year == 1985 & name == "Eggli" & firstname == "Jacques Simon" & gdenr_2018 == 99999
replace ID = "GE-1975-9001" if ID == "" & year == 1983 & name == "Eggly" & firstname == "Jacques Simon" & gdenr_2018 == 6621
replace ID = "ZH-1963-0079" if ID == "" & year == 1990 & name == "Eisenring Frick" & firstname == "Paul" & gdenr_2018 == 151
replace ID = "NW-1995-0002" if ID == "" & year == 2000 & name == "Engelberger" & firstname == "Eduard sen." & gdenr_2018 == 1509
replace ID = "NW-1995-0002" if ID == "" & year == 2001 & name == "Engelberger" & firstname == "Eduard sen." & gdenr_2018 == 1509
replace ID = "NW-1995-0002" if ID == "" & year == 2002 & name == "Engelberger" & firstname == "Eduard sen." & gdenr_2018 == 1509
replace ID = "NW-1995-0002" if ID == "" & year == 2003 & name == "Engelberger" & firstname == "Eduard sen." & gdenr_2018 == 1509
replace ID = "AI-1987-0002" if ID == "" & year == 1987 & name == "Engler" & firstname == "Rolf" & gdenr_2018 == 3101
replace ID = "VS-1991-0012" if ID == "" & year == 1991 & name == "Epiney" & firstname == "Simon" & gdenr_2018 == 6252
replace ID = "SG-1975-0018" if ID == "" & year == 1982 & name == "Eppenberger" & firstname == "Susi" & gdenr_2018 == 3360
replace ID = "SG-1975-0018" if ID == "" & year == 1983 & name == "Eppenberger" & firstname == "Susi" & gdenr_2018 == 3360
replace ID = "SG-1975-0018" if ID == "" & year == 1988 & name == "Eppenberger" & firstname == "Susi" & gdenr_2018 == 3360
replace ID = "SG-1975-0018" if ID == "" & year == 1989 & name == "Eppenberger" & firstname == "Susi" & gdenr_2018 == 3360
replace ID = "SG-1975-0018" if ID == "" & year == 1990 & name == "Eppenberger" & firstname == "Susi" & gdenr_2018 == 3360
replace ID = "LU-1999-0025" if ID == "" & year == 2000 & name == "Estermann" & firstname == "Heinrich" & gdenr_2018 == 99999
replace ID = "LU-1999-0025" if ID == "" & year == 2001 & name == "Estermann" & firstname == "Heinrich" & gdenr_2018 == 99999
replace ID = "BL-1979-0013" if ID == "" & year == 1997 & name == "Fankhauser" & firstname == "Angeline" & gdenr_2018 == 2771
replace ID = "BL-1979-0013" if ID == "" & year == 1998 & name == "Fankhauser" & firstname == "Angeline" & gdenr_2018 == 2771
replace ID = "BL-1979-0013" if ID == "" & year == 1999 & name == "Fankhauser" & firstname == "Angeline" & gdenr_2018 == 2771
replace ID = "FR-1991-0019" if ID == "" & year == 1991 & name == "Fasel" & firstname == "Hugo" & gdenr_2018 == 2304
replace ID = "VD-1999-0084" if ID == "" & year == 2000 & name == "Favre" & firstname == "Charles" & gdenr_2018 == 5586
replace ID = "VD-1999-0084" if ID == "" & year == 2001 & name == "Favre" & firstname == "Charles" & gdenr_2018 == 5586
replace ID = "VD-1999-0084" if ID == "" & year == 2002 & name == "Favre" & firstname == "Charles" & gdenr_2018 == 5586
replace ID = "VD-1999-0084" if ID == "" & year == 2003 & name == "Favre" & firstname == "Charles" & gdenr_2018 == 5586
replace ID = "BEJU-1971-0093" if ID == "" & year == 1983 & name == "Fehr" & firstname == "Hermann" & gdenr_2018 == 371
replace ID = "ZH-1991-0177" if ID == "" & year == 1991 & name == "Fehr" & firstname == "Lisbeth" & gdenr_2018 == 32
replace ID = "ZH-1995-0175" if ID == "" & year == 2000 & name == "Fehr" & firstname == "Mario" & gdenr_2018 == 99999
replace ID = "ZH-1995-0175" if ID == "" & year == 2001 & name == "Fehr" & firstname == "Mario" & gdenr_2018 == 99999
replace ID = "ZH-1991-0177" if ID == "" & year == 1996 & name == "Fehr Ehrensperger" & firstname == "Lisbeth" & gdenr_2018 == 32
replace ID = "ZH-1991-0177" if ID == "" & year == 1997 & name == "Fehr Ehrensperger" & firstname == "Lisbeth" & gdenr_2018 == 32
replace ID = "BL-1971-0014" if ID == "" & year == 1982 & name == "Feigenwinter" & firstname == "Hans Rud." & gdenr_2018 == 2773
replace ID = "BL-1971-0014" if ID == "" & year == 1984 & name == "Feigenwinter" & firstname == "Hans Rud." & gdenr_2018 == 2773
replace ID = "BL-1971-0014" if ID == "" & year == 1985 & name == "Feigenwinter" & firstname == "Hans Rud." & gdenr_2018 == 2773
replace ID = "BEJU-1951-0046" if ID == "" & year == 1982 & name == "Fischer" & firstname == " Otto" & gdenr_2018 == 351
replace ID = "AG-1971-0027" if ID == "" & year == 1997 & name == "Fischer" & firstname == "Theo" & gdenr_2018 == 4082
replace ID = "AG-1971-0027" if ID == "" & year == 1998 & name == "Fischer" & firstname == "Theo" & gdenr_2018 == 4082
replace ID = "AG-1971-0027" if ID == "" & year == 1999 & name == "Fischer" & firstname == "Theo" & gdenr_2018 == 4068
replace ID = "BL-1967-0006" if ID == "" & year == 1982 & name == "Flubacher Haas" & firstname == "Karl" & gdenr_2018 == 2852
replace ID = "BL-1967-0006" if ID == "" & year == 1983 & name == "Flubacher Haas" & firstname == "Karl" & gdenr_2018 == 2852
replace ID = "BL-1967-0006" if ID == "" & year == 1984 & name == "Flubacher Haas" & firstname == "Karl" & gdenr_2018 == 2852
replace ID = "BL-1967-0006" if ID == "" & year == 1985 & name == "Flubacher Haas" & firstname == "Karl" & gdenr_2018 == 2852
replace ID = "VD-1947-0032" if ID == "" & year == 1982 & name == "Forel" & firstname == "Armand Auguste" & gdenr_2018 == 5724
replace ID = "ZH-1983-0119" if ID == "" & year == 1987 & name == "Frey" & firstname == "Walter" & gdenr_2018 == 154
replace ID = "VD-1987-0057" if ID == "" & year == 1987 & name == "Friderici" & firstname == "Charles" & gdenr_2018 == 5639
replace ID = "VD-1987-0057" if ID == "" & year == 1994 & name == "Friderici" & firstname == "Charles E." & gdenr_2018 == 5639
replace ID = "VD-1987-0057" if ID == "" & year == 1996 & name == "Friderici" & firstname == "Charles E." & gdenr_2018 == 5639
replace ID = "VD-1987-0057" if ID == "" & year == 1997 & name == "Friderici" & firstname == "Charles E." & gdenr_2018 == 5639
replace ID = "VD-1987-0057" if ID == "" & year == 1998 & name == "Friderici" & firstname == "Charles E." & gdenr_2018 == 5639
replace ID = "VD-1987-0057" if ID == "" & year == 1999 & name == "Friderici" & firstname == "Charles E." & gdenr_2018 == 5639
replace ID = "BEJU-1975-0151" if ID == "" & year == 1983 & name == "Friedli" & firstname == "Valentine" & gdenr_2018 == 6711
replace ID = "BEJU-1975-0151" if ID == "" & year == 1984 & name == "Friedli" & firstname == "Valentine" & gdenr_2018 == 6711
replace ID = "BEJU-1975-0151" if ID == "" & year == 1985 & name == "Friedli" & firstname == "Valentine" & gdenr_2018 == 6711
replace ID = "BEJU-1975-0151" if ID == "" & year == 1986 & name == "Friedli" & firstname == "Valentine" & gdenr_2018 == 6711
replace ID = "AR-1975-0002" if ID == "" & year == 1982 & name == "Früh" & firstname == "Hans Rud." & gdenr_2018 == 3021
replace ID = "AR-1975-0002" if ID == "" & year == 1983 & name == "Früh" & firstname == "Hans Rudolf" & gdenr_2018 == 3021
replace ID = "AR-1975-0002" if ID == "" & year == 1988 & name == "Früh" & firstname == "Hans Rudolf" & gdenr_2018 == 3021
replace ID = "AR-1975-0002" if ID == "" & year == 1989 & name == "Früh" & firstname == "Hans Rudolf" & gdenr_2018 == 3021
replace ID = "AR-1975-0002" if ID == "" & year == 1990 & name == "Früh" & firstname == "Hans Rudolf" & gdenr_2018 == 3021
replace ID = "AR-1975-0002" if ID == "" & year == 1991 & name == "Früh" & firstname == "Hans Rudolf" & gdenr_2018 == 3021
replace ID = "LU-1987-0019" if ID == "" & year == 1987 & name == "Fäh" & firstname == "Paul" & gdenr_2018 == 1061
replace ID = "SG-1991-0027" if ID == "" & year == 1997 & name == "Fässler Osterwalder" & firstname == "Hildegard" & gdenr_2018 == 3273
replace ID = "SG-1991-0027" if ID == "" & year == 1998 & name == "Fässler Osterwalder" & firstname == "Hildegard" & gdenr_2018 == 3273
replace ID = "SG-1991-0027" if ID == "" & year == 1999 & name == "Fässler Osterwalder" & firstname == "Hildegard" & gdenr_2018 == 3273
replace ID = "SG-1991-0027" if ID == "" & year == 2000 & name == "Fässler Osterwalder" & firstname == "Hildegard" & gdenr_2018 == 99999
replace ID = "SG-1991-0027" if ID == "" & year == 2002 & name == "Fässler Osterwalder" & firstname == "Hildegard" & gdenr_2018 == 99999
replace ID = "SG-1991-0027" if ID == "" & year == 2003 & name == "Fässler Osterwalder" & firstname == "Hildegard" & gdenr_2018 == 99999
replace ID = "SZ-1995-0007" if ID == "" & year == 2000 & name == "Föhn" & firstname == "Peter" & gdenr_2018 == 99999
replace ID = "SZ-1995-0007" if ID == "" & year == 2001 & name == "Föhn" & firstname == "Peter" & gdenr_2018 == 99999
replace ID = "SO-1975-0006" if ID == "" & year == 1982 & name == "Füeg" & firstname == "Cornelia" & gdenr_2018 == 2502
replace ID = "GR-1987-0014" if ID == "" & year == 1996 & name == "Gadient" & firstname == "Brigitta M." & gdenr_2018 == 3901
replace ID = "GR-1987-0014" if ID == "" & year == 1997 & name == "Gadient" & firstname == "Brigitta M." & gdenr_2018 == 3901
replace ID = "GR-1987-0014" if ID == "" & year == 1998 & name == "Gadient" & firstname == "Brigitta M." & gdenr_2018 == 3901
replace ID = "GR-1987-0014" if ID == "" & year == 1999 & name == "Gadient" & firstname == "Brigitta M." & gdenr_2018 == 3901
replace ID = "GR-1987-0014" if ID == "" & year == 2000 & name == "Gadient" & firstname == "Brigitta M." & gdenr_2018 == 3901
replace ID = "GR-1987-0014" if ID == "" & year == 2001 & name == "Gadient" & firstname == "Brigitta M." & gdenr_2018 == 3901
replace ID = "GR-1987-0014" if ID == "" & year == 2002 & name == "Gadient" & firstname == "Brigitta M." & gdenr_2018 == 3901
replace ID = "GR-1987-0014" if ID == "" & year == 2003 & name == "Gadient" & firstname == "Brigitta M." & gdenr_2018 == 3901
replace ID = "BEJU-1987-0148" if ID == "" & year == 2000 & name == "Galli" & firstname == "Remo" & gdenr_2018 == 99999
replace ID = "BEJU-1987-0148" if ID == "" & year == 2001 & name == "Galli" & firstname == "Remo" & gdenr_2018 == 99999
replace ID = "NE-1999-0016" if ID == "" & year == 2000 & name == "Garbani" & firstname == "Valérie" & gdenr_2018 == 99999
replace ID = "NE-1999-0016" if ID == "" & year == 2001 & name == "Garbani" & firstname == "Valérie" & gdenr_2018 == 99999
replace ID = "NE-1999-0016" if ID == "" & year == 2002 & name == "Garbani" & firstname == "Valérie" & gdenr_2018 == 99999
replace ID = "NE-1999-0016" if ID == "" & year == 2003 & name == "Garbani" & firstname == "Valérie" & gdenr_2018 == 99999
replace ID = "VD-1987-0058" if ID == "" & year == 1992 & name == "Gardiol" & firstname == "Irène" & gdenr_2018 == 5590
replace ID = "VD-1987-0058" if ID == "" & year == 1993 & name == "Gardiol" & firstname == "Irène" & gdenr_2018 == 5590
replace ID = "BEJU-1975-0164" if ID == "" & year == 1982 & name == "Gehler" & firstname == " Jean Paul" & gdenr_2018 == 703
replace ID = "BEJU-1975-0164" if ID == "" & year == 1983 & name == "Gehler" & firstname == " Jean Paul" & gdenr_2018 == 703
replace ID = "BEJU-1975-0164" if ID == "" & year == 1984 & name == "Gehler" & firstname == " Jean Paul" & gdenr_2018 == 703
replace ID = "BEJU-1975-0164" if ID == "" & year == 1985 & name == "Gehler" & firstname == " Jean Paul" & gdenr_2018 == 703
replace ID = "BEJU-1975-0164" if ID == "" & year == 1986 & name == "Gehler" & firstname == " Jean Paul" & gdenr_2018 == 703
replace ID = "ZH-1987-0216" if ID == "" & year == 2000 & name == "Genner Ritzmann" & firstname == "Ruth" & gdenr_2018 == 261
replace ID = "ZH-1987-0216" if ID == "" & year == 2001 & name == "Genner Ritzmann" & firstname == "Ruth" & gdenr_2018 == 261
replace ID = "ZH-1987-0216" if ID == "" & year == 2002 & name == "Genner Ritzmann" & firstname == "Ruth" & gdenr_2018 == 261
replace ID = "ZH-1987-0216" if ID == "" & year == 2003 & name == "Genner Ritzmann" & firstname == "Ruth" & gdenr_2018 == 261
replace ID = "AG-1991-0042" if ID == "" & year == 1991 & name == "Giezendanner" & firstname == "Ulrich" & gdenr_2018 == 4282
replace ID = "SG-1967-0016" if ID == "" & year == 1983 & name == "Giger" & firstname == "Titus" & gdenr_2018 == 3295
replace ID = "VD-1971-0060" if ID == "" & year == 1982 & name == "Girard" & firstname == "Gertrude" & gdenr_2018 == 5889
replace ID = "FR-1983-0015" if ID == "" & year == 1991 & name == "Gobet" & firstname == "Alexis" & gdenr_2018 == 2113
replace ID = "BL-1987-0019" if ID == "" & year == 1997 & name == "Gonseth" & firstname == "Ruth" & gdenr_2018 == 2829
replace ID = "BL-1987-0019" if ID == "" & year == 1998 & name == "Gonseth" & firstname == "Ruth" & gdenr_2018 == 2829
replace ID = "BL-1987-0019" if ID == "" & year == 1999 & name == "Gonseth" & firstname == "Ruth" & gdenr_2018 == 2829
replace ID = "BL-1987-0019" if ID == "" & year == 2000 & name == "Gonseth" & firstname == "Ruth" & gdenr_2018 == 2829
replace ID = "BL-1987-0019" if ID == "" & year == 2001 & name == "Gonseth" & firstname == "Ruth" & gdenr_2018 == 2829
replace ID = "TI-1975-9026" if ID == "" & year == 1983 & name == "Grassi" & firstname == "Mario" & gdenr_2018 == 5196
replace ID = "TI-1975-9026" if ID == "" & year == 1987 & name == "Grassi" & firstname == "Mario P." & gdenr_2018 == 5196
replace ID = "TI-1975-9026" if ID == "" & year == 1988 & name == "Grassi" & firstname == "Mario P." & gdenr_2018 == 5196
replace ID = "TI-1975-9026" if ID == "" & year == 1989 & name == "Grassi" & firstname == "Mario P." & gdenr_2018 == 5196
replace ID = "TI-1975-9026" if ID == "" & year == 1990 & name == "Grassi" & firstname == "Mario P." & gdenr_2018 == 5196
replace ID = "ZH-1979-0179" if ID == "" & year == 1985 & name == "Grendelmeier" & firstname == "Verena" & gdenr_2018 == 99999
replace ID = "GE-1971-9097" if ID == "" & year == 1996 & name == "Grobet" & firstname == "Christian" & gdenr_2018 == 6621
replace ID = "GE-1971-9097" if ID == "" & year == 1997 & name == "Grobet" & firstname == "Christian" & gdenr_2018 == 6621
replace ID = "GE-1971-9097" if ID == "" & year == 1998 & name == "Grobet" & firstname == "Christian" & gdenr_2018 == 6621
replace ID = "GE-1971-9097" if ID == "" & year == 1999 & name == "Grobet" & firstname == "Christian" & gdenr_2018 == 6621
replace ID = "GE-1971-9097" if ID == "" & year == 2000 & name == "Grobet" & firstname == "Christian" & gdenr_2018 == 99999
replace ID = "GE-1971-9097" if ID == "" & year == 2001 & name == "Grobet" & firstname == "Christian" & gdenr_2018 == 99999
replace ID = "GE-1971-9097" if ID == "" & year == 2002 & name == "Grobet" & firstname == "Christian" & gdenr_2018 == 99999
replace ID = "GE-1971-9097" if ID == "" & year == 2003 & name == "Grobet" & firstname == "Christian" & gdenr_2018 == 99999
replace ID = "GE-1987-0021" if ID == "" & year == 1987 & name == "Gros" & firstname == "Jean Michel" & gdenr_2018 == 6638
replace ID = "TI-1983-0025" if ID == "" & year == 1983 & name == "Guidici" & firstname == "Luciano" & gdenr_2018 == 5113
replace ID = "NE-1987-0016" if ID == "" & year == 1987 & name == "Guinand" & firstname == "Jean" & gdenr_2018 == 6458
replace ID = "VD-1991-0064" if ID == "" & year == 1998 & name == "Guisan" & firstname == "Yves" & gdenr_2018 == 5886
replace ID = "VD-1991-0064" if ID == "" & year == 1999 & name == "Guisan" & firstname == "Yves" & gdenr_2018 == 5886
replace ID = "VD-1991-0064" if ID == "" & year == 2000 & name == "Guisan" & firstname == "Yves" & gdenr_2018 == 5886
replace ID = "VD-1991-0064" if ID == "" & year == 2001 & name == "Guisan" & firstname == "Yves" & gdenr_2018 == 5886
replace ID = "VD-1991-0064" if ID == "" & year == 2002 & name == "Guisan" & firstname == "Yves" & gdenr_2018 == 5886
replace ID = "VD-1991-0064" if ID == "" & year == 2003 & name == "Guisan" & firstname == "Yves" & gdenr_2018 == 5886
replace ID = "BL-1987-0021" if ID == "" & year == 1987 & name == "Gysin Buser" & firstname == "Hans Rudolf" & gdenr_2018 == 2831
replace ID = "BL-1987-0021" if ID == "" & year == 1988 & name == "Gysin Buser" & firstname == "Hans Rudolf" & gdenr_2018 == 2831
replace ID = "BL-1987-0021" if ID == "" & year == 1989 & name == "Gysin Buser" & firstname == "Hans Rudolf" & gdenr_2018 == 2831
replace ID = "BL-1987-0021" if ID == "" & year == 1990 & name == "Gysin Buser" & firstname == "Hans Rudolf" & gdenr_2018 == 2831
replace ID = "BL-1987-0021" if ID == "" & year == 1991 & name == "Gysin Buser" & firstname == "Hans Rudolf" & gdenr_2018 == 2831
replace ID = "BL-1987-0021" if ID == "" & year == 1992 & name == "Gysin Buser" & firstname == "Hans Rudolf" & gdenr_2018 == 2831
replace ID = "BL-1987-0021" if ID == "" & year == 1993 & name == "Gysin Buser" & firstname == "Hans Rudolf" & gdenr_2018 == 2831
replace ID = "BL-1987-0021" if ID == "" & year == 1994 & name == "Gysin Buser" & firstname == "Hans Rudolf" & gdenr_2018 == 2831
replace ID = "BL-1987-0021" if ID == "" & year == 1996 & name == "Gysin Buser" & firstname == "Hans Rudolf" & gdenr_2018 == 2831
replace ID = "BL-1987-0021" if ID == "" & year == 2003 & name == "Gysin Buser" & firstname == "Hans Rudolf" & gdenr_2018 == 2831
replace ID = "BEJU-1971-7031" if ID == "" & year == 1996 & name == "Günter" & firstname == "Paul" & gdenr_2018 == 593
replace ID = "BEJU-1971-7031" if ID == "" & year == 1997 & name == "Günter" & firstname == "Paul" & gdenr_2018 == 593
replace ID = "BEJU-1971-7031" if ID == "" & year == 1998 & name == "Günter" & firstname == "Paul" & gdenr_2018 == 593
replace ID = "BEJU-1971-7031" if ID == "" & year == 1999 & name == "Günter" & firstname == "Paul" & gdenr_2018 == 593
replace ID = "BEJU-1987-0189" if ID == "" & year == 1987 & name == "Hafner" & firstname == "Rudolf" & gdenr_2018 == 351
replace ID = "SH-1987-9014" if ID == "" & year == 1987 & name == "Hafner" & firstname == "Ursula" & gdenr_2018 == 2939
replace ID = "SH-1987-9014" if ID == "" & year == 1992 & name == "Hafner" & firstname == "Ursula" & gdenr_2018 == 2939
replace ID = "SH-1987-9014" if ID == "" & year == 1993 & name == "Hafner" & firstname == "Ursula" & gdenr_2018 == 2939
replace ID = "SH-1987-9014" if ID == "" & year == 1994 & name == "Hafner" & firstname == "Ursula" & gdenr_2018 == 2939
replace ID = "SH-1987-9014" if ID == "" & year == 1998 & name == "Hafner" & firstname == "Ursula" & gdenr_2018 == 2939
replace ID = "SH-1987-9014" if ID == "" & year == 1999 & name == "Hafner" & firstname == "Ursula" & gdenr_2018 == 2939
replace ID = "GR-1999-0018" if ID == "" & year == 2000 & name == "Hassler" & firstname == "Hansjörg" & gdenr_2018 == 99999
replace ID = "GR-1999-0018" if ID == "" & year == 2001 & name == "Hassler" & firstname == "Hansjörg" & gdenr_2018 == 99999
replace ID = "ZH-1987-0274" if ID == "" & year == 1991 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 1992 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 1993 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 1994 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 1996 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 1997 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 1998 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 1999 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 2000 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 2001 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 2002 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "ZH-1987-0274" if ID == "" & year == 2003 & name == "Heberlein" & firstname == "Beatrix" & gdenr_2018 == 160
replace ID = "SO-1995-0023" if ID == "" & year == 1998 & name == "Heim" & firstname == "Alex" & gdenr_2018 == 2404
replace ID = "SO-1995-0023" if ID == "" & year == 1999 & name == "Heim" & firstname == "Alex" & gdenr_2018 == 2404
replace ID = "SO-1995-0023" if ID == "" & year == 2000 & name == "Heim" & firstname == "Alex" & gdenr_2018 == 99999
replace ID = "SO-1995-0023" if ID == "" & year == 2001 & name == "Heim" & firstname == "Alex" & gdenr_2018 == 99999
replace ID = "SO-1995-0023" if ID == "" & year == 2002 & name == "Heim" & firstname == "Alex" & gdenr_2018 == 99999
replace ID = "SO-1995-0023" if ID == "" & year == 2003 & name == "Heim" & firstname == "Alex" & gdenr_2018 == 99999
replace ID = "TG-1987-9010" if ID == "" & year == 1987 & name == "Hess" & firstname == "Otto" & gdenr_2018 == 4431
replace ID = "VS-1987-0019" if ID == "" & year == 1987 & name == "Hildbrand" & firstname == "Franz" & gdenr_2018 == 6118
replace ID = "BEJU-1991-0218" if ID == "" & year == 2000 & name == "Hochreutener" & firstname == "Norbert" & gdenr_2018 == 355
replace ID = "BEJU-1991-0218" if ID == "" & year == 1996 & name == "Hochreutner" & firstname == "Norbert" & gdenr_2018 == 99999
replace ID = "BEJU-1991-0218" if ID == "" & year == 1997 & name == "Hochreutner" & firstname == "Norbert" & gdenr_2018 == 355
replace ID = "AG-1995-0081" if ID == "" & year == 2000 & name == "Hofmann" & firstname == "Urs" & gdenr_2018 == 99999
replace ID = "SG-1991-0039" if ID == "" & year == 1991 & name == "Hollenstein" & firstname == "Pia" & gdenr_2018 == 3203
replace ID = "BS-1959-0028" if ID == "" & year == 1997 & name == "Hubacher" & firstname == "Helmut" & gdenr_2018 == 6807
replace ID = "ZH-1991-0307" if ID == "" & year == 1996 & name == "Hubmann" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0307" if ID == "" & year == 1997 & name == "Hubmann" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0307" if ID == "" & year == 1998 & name == "Hubmann" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0307" if ID == "" & year == 1999 & name == "Hubmann" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0307" if ID == "" & year == 2000 & name == "Hubmann" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0307" if ID == "" & year == 2001 & name == "Hubmann" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0307" if ID == "" & year == 2002 & name == "Hubmann" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0307" if ID == "" & year == 2003 & name == "Hubmann" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "GR-1991-0018" if ID == "" & year == 1991 & name == "Hämmerle" & firstname == "Andrea" & gdenr_2018 == 3673
replace ID = "GR-1991-0018" if ID == "" & year == 2002 & name == "Hämmerle" & firstname == "Andrea" & gdenr_2018 == 99999
replace ID = "SO-1987-0023" if ID == "" & year == 1987 & name == "Hänggi" & firstname == "Peter" & gdenr_2018 == 2621
replace ID = "GL-1979-0003" if ID == "" & year == 1982 & name == "Hösli" & firstname == "Fritz" & gdenr_2018 == 1631
replace ID = "GL-1979-0003" if ID == "" & year == 1983 & name == "Hösli" & firstname == "Fritz" & gdenr_2018 == 1631
replace ID = "OW-2003-0002" if ID == "" & year == 2002 & name == "Imfeld" & firstname == "Adrian" & gdenr_2018 == 1407
replace ID = "OW-2003-0002" if ID == "" & year == 2003 & name == "Imfeld" & firstname == "Adrian" & gdenr_2018 == 1407
replace ID = "NW-1979-0003" if ID == "" & year == 1983 & name == "Iten" & firstname == "Josef" & gdenr_2018 == 1507
replace ID = "NW-1979-0003" if ID == "" & year == 1984 & name == "Iten" & firstname == "Josef" & gdenr_2018 == 1507
replace ID = "VD-1987-0084" if ID == "" & year == 1987 & name == "Jeanprêtre" & firstname == "Francine" & gdenr_2018 == 5642
replace ID = "BEJU-1987-0227" if ID == "" & year == 1991 & name == "Jenni" & firstname == "Peter" & gdenr_2018 == 888
replace ID = "VS-1995-9048" if ID == "" & year == 2000 & name == "Jossen Zinsstag" & firstname == "Peter" & gdenr_2018 == 99999
replace ID = "VS-1995-9048" if ID == "" & year == 2001 & name == "Jossen Zinsstag" & firstname == "Peter" & gdenr_2018 == 99999
replace ID = "SO-1991-0027" if ID == "" & year == 1991 & name == "Jäggi" & firstname == "Paul" & gdenr_2018 == 2524
replace ID = "ZH-1983-0240" if ID == "" & year == 2003 & name == "Keller" & firstname == "Robert" & gdenr_2018 == 177
replace ID = "ZH-1991-0354" if ID == "" & year == 1991 & name == "Kern" & firstname == "Armin" & gdenr_2018 == 230
replace ID = "SO-1991-0031" if ID == "" & year == 1996 & name == "Kofmel" & firstname == "Peter" & gdenr_2018 == 2516
replace ID = "SO-1991-0031" if ID == "" & year == 1997 & name == "Kofmel" & firstname == "Peter" & gdenr_2018 == 2516
replace ID = "SO-1991-0031" if ID == "" & year == 1998 & name == "Kofmel" & firstname == "Peter" & gdenr_2018 == 2516
replace ID = "SO-1991-0031" if ID == "" & year == 1999 & name == "Kofmel" & firstname == "Peter" & gdenr_2018 == 2516
replace ID = "AI-1971-9000" if ID == "" & year == 1986 & name == "Koller Brander" & firstname == "Arnold" & gdenr_2018 == 3101
replace ID = "AG-1991-9152" if ID == "" & year == 1999 & name == "Kuhn" & firstname == "Kattin" & gdenr_2018 == 4082
replace ID = "ZH-1963-0168" if ID == "" & year == 1986 & name == "Künzi" & firstname == "Hans Paul" & gdenr_2018 == 99999
replace ID = "ZH-1963-0168" if ID == "" & year == 1982 & name == "Künzi Girsberger" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "ZH-1963-0168" if ID == "" & year == 1983 & name == "Künzi Girsberger" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "ZH-1963-0168" if ID == "" & year == 1984 & name == "Künzi Girsberger" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "ZH-1963-0168" if ID == "" & year == 1985 & name == "Künzi Girsberger" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "SZ-1999-0014" if ID == "" & year == 2000 & name == "Lalive dÈpinay Buff" & firstname == "Maya" & gdenr_2018 == 1322
replace ID = "SZ-1999-0014" if ID == "" & year == 2001 & name == "Lalive dÈpinay Buff" & firstname == "Maya" & gdenr_2018 == 1322
replace ID = "SZ-1999-0014" if ID == "" & year == 2002 & name == "Lalive dÈpinay Buff" & firstname == "Maya" & gdenr_2018 == 1322
replace ID = "SZ-1999-0014" if ID == "" & year == 2003 & name == "Lalive dÈpinay Buff" & firstname == "Maya" & gdenr_2018 == 1322
replace ID = "VD-1987-0091" if ID == "" & year == 1996 & name == "Langenberger Jaeger" & firstname == "Christiane" & gdenr_2018 == 5645
replace ID = "VD-1987-0091" if ID == "" & year == 1997 & name == "Langenberger Jaeger" & firstname == "Christiane" & gdenr_2018 == 5645
replace ID = "VD-1987-0091" if ID == "" & year == 1998 & name == "Langenberger Jaeger" & firstname == "Christiane" & gdenr_2018 == 5645
replace ID = "VD-1987-0091" if ID == "" & year == 1999 & name == "Langenberger Jaeger" & firstname == "Christiane" & gdenr_2018 == 5645
replace ID = "LU-1983-0023" if ID == "" & year == 1983 & name == "Lanz" & firstname == "Fritz" & gdenr_2018 == 1061
replace ID = "ZH-1987-0421" if ID == "" & year == 1987 & name == "Ledergerber" & firstname == "Elmar" & gdenr_2018 == 261
replace ID = "LU-1991-0024" if ID == "" & year == 1996 & name == "Leu" & firstname == "Josef" & gdenr_2018 == 1032
replace ID = "LU-1991-0024" if ID == "" & year == 1997 & name == "Leu" & firstname == "Josef" & gdenr_2018 == 1032
replace ID = "LU-1991-0024" if ID == "" & year == 1998 & name == "Leu" & firstname == "Josef" & gdenr_2018 == 1032
replace ID = "LU-1991-0024" if ID == "" & year == 1999 & name == "Leu" & firstname == "Josef" & gdenr_2018 == 1032
replace ID = "LU-1991-0024" if ID == "" & year == 1991 & name == "Leu Morgenthaler" & firstname == "Josef" & gdenr_2018 == 1032
replace ID = "VD-1987-0095" if ID == "" & year == 1987 & name == "Leuba" & firstname == "Jean François" & gdenr_2018 == 5601
replace ID = "VD-1987-0095" if ID == "" & year == 1996 & name == "Leuba" & firstname == "Jean François" & gdenr_2018 == 5601
replace ID = "VD-1987-0095" if ID == "" & year == 1997 & name == "Leuba" & firstname == "Jean François" & gdenr_2018 == 5601
replace ID = "VD-1987-0095" if ID == "" & year == 1998 & name == "Leuba" & firstname == "Jean François" & gdenr_2018 == 5601
replace ID = "SO-1975-9021" if ID == "" & year == 1983 & name == "Leuenberger ( Ramseier)" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1984 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1985 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1986 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1987 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1988 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1989 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1990 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1991 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1992 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1993 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1994 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1996 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1997 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "SO-1975-9021" if ID == "" & year == 1998 & name == "Leuenberger Ramseier" & firstname == "Ernst" & gdenr_2018 == 2601
replace ID = "BL-1983-9011" if ID == "" & year == 1988 & name == "Leutenegger Oberholzer" & firstname == "Susanne" & gdenr_2018 == 2762
replace ID = "BL-1983-9011" if ID == "" & year == 1989 & name == "Leutenegger Oberholzer" & firstname == "Susanne" & gdenr_2018 == 2762
replace ID = "BL-1983-9011" if ID == "" & year == 1990 & name == "Leutenegger Oberholzer" & firstname == "Susanne" & gdenr_2018 == 2762
replace ID = "BL-1983-9011" if ID == "" & year == 2000 & name == "Leutenegger Oberholzer" & firstname == "Susanne" & gdenr_2018 == 99999
replace ID = "BL-1983-9011" if ID == "" & year == 2001 & name == "Leutenegger Oberholzer" & firstname == "Susanne" & gdenr_2018 == 99999
replace ID = "BL-1983-9011" if ID == "" & year == 2002 & name == "Leutenegger Oberholzer" & firstname == "Susanne" & gdenr_2018 == 99999
replace ID = "BL-1983-9011" if ID == "" & year == 2003 & name == "Leutenegger Oberholzer" & firstname == "Susanne" & gdenr_2018 == 99999
replace ID = "AG-1999-0125" if ID == "" & year == 2002 & name == "Leuthard Hausin" & firstname == "Doris" & gdenr_2018 == 4234
replace ID = "AG-1999-0125" if ID == "" & year == 2003 & name == "Leuthard Hausin" & firstname == "Doris" & gdenr_2018 == 4234
replace ID = "BEJU-1983-0233" if ID == "" & year == 1987 & name == "Loeb Fischer" & firstname == "François" & gdenr_2018 == 356
replace ID = "BEJU-1983-0233" if ID == "" & year == 1988 & name == "Loeb Fischer" & firstname == "François" & gdenr_2018 == 356
replace ID = "LU-1999-0054" if ID == "" & year == 2000 & name == "Lustenberger" & firstname == "Rudolf" & gdenr_2018 == 1007
replace ID = "LU-1999-0054" if ID == "" & year == 2001 & name == "Lustenberger" & firstname == "Rudolf" & gdenr_2018 == 1007
replace ID = "LU-1999-0054" if ID == "" & year == 2002 & name == "Lustenberger" & firstname == "Rudolf" & gdenr_2018 == 1007
replace ID = "LU-1999-0054" if ID == "" & year == 2003 & name == "Lustenberger" & firstname == "Rudolf" & gdenr_2018 == 1007
replace ID = "ZH-1959-0153" if ID == "" & year == 1982 & name == "Lüchinger" & firstname == "Hans G." & gdenr_2018 == 14
replace ID = "AR-1983-0004" if ID == "" & year == 1983 & name == "Maeder" & firstname == "Herbert" & gdenr_2018 == 3034
replace ID = "AR-1983-0004" if ID == "" & year == 1988 & name == "Maeder" & firstname == "Herbert" & gdenr_2018 == 3034
replace ID = "AR-1983-0004" if ID == "" & year == 1989 & name == "Maeder" & firstname == "Herbert" & gdenr_2018 == 3034
replace ID = "AR-1983-0004" if ID == "" & year == 1990 & name == "Maeder" & firstname == "Herbert" & gdenr_2018 == 3034
replace ID = "AR-1983-0004" if ID == "" & year == 1991 & name == "Maeder" & firstname == "Herbert" & gdenr_2018 == 3034
replace ID = "VD-1999-0137" if ID == "" & year == 2000 & name == "Maillard" & firstname == "Pierre Yves" & gdenr_2018 == 99999
replace ID = "VD-1999-0137" if ID == "" & year == 2001 & name == "Maillars" & firstname == "Pierre Yves" & gdenr_2018 == 99999
replace ID = "VD-1999-0137" if ID == "" & year == 2002 & name == "Maillars" & firstname == "Pierre Yves" & gdenr_2018 == 99999
replace ID = "VD-1999-0137" if ID == "" & year == 2003 & name == "Maillars" & firstname == "Pierre Yves" & gdenr_2018 == 99999
replace ID = "VD-1991-0091" if ID == "" & year == 1991 & name == "Mamie" & firstname == "Philippe" & gdenr_2018 == 5764
replace ID = "GL-1991-9001" if ID == "" & year == 1991 & name == "Marti" & firstname == "Werner" & gdenr_2018 == 1631
replace ID = "GL-1991-9001" if ID == "" & year == 1999 & name == "Marti" & firstname == "Werner" & gdenr_2018 == 99999
replace ID = "GL-1991-9001" if ID == "" & year == 2002 & name == "Marti" & firstname == "Werner" & gdenr_2018 == 1632
replace ID = "GL-1991-9001" if ID == "" & year == 2003 & name == "Marti" & firstname == "Werner" & gdenr_2018 == 1632
replace ID = "VD-1987-0102" if ID == "" & year == 1987 & name == "Martin" & firstname == "Paul René" & gdenr_2018 == 5586
replace ID = "ZH-1991-0431" if ID == "" & year == 2000 & name == "Marty Kälin" & firstname == "Barbara" & gdenr_2018 == 99999
replace ID = "ZH-1991-0431" if ID == "" & year == 2001 & name == "Marty Kälin" & firstname == "Barbara" & gdenr_2018 == 298
replace ID = "ZH-1991-0431" if ID == "" & year == 2002 & name == "Marty Kälin" & firstname == "Barbara" & gdenr_2018 == 298
replace ID = "ZH-1991-0431" if ID == "" & year == 2003 & name == "Marty Kälin" & firstname == "Barbara" & gdenr_2018 == 298
replace ID = "BS-1975-0042" if ID == "" & year == 1982 & name == "Mascarin" & firstname == "Ruth" & gdenr_2018 == 2701
replace ID = "BS-1975-0042" if ID == "" & year == 1983 & name == "Mascarin" & firstname == "Ruth" & gdenr_2018 == 2701
replace ID = "TI-1991-0042" if ID == "" & year == 1991 & name == "Maspoli" & firstname == "Flavio" & gdenr_2018 == 5118
replace ID = "VD-1967-0065" if ID == "" & year == 1985 & name == "Massy" & firstname == "Claude" & gdenr_2018 == 99999
replace ID = "AG-1991-0081" if ID == "" & year == 2000 & name == "Mathys" & firstname == "Hans Ulrich" & gdenr_2018 == 99999
replace ID = "AG-1991-0081" if ID == "" & year == 2001 & name == "Mathys" & firstname == "Hans Ulrich" & gdenr_2018 == 99999
replace ID = "NE-1987-0019" if ID == "" & year == 1987 & name == "Matthey" & firstname == "Francis" & gdenr_2018 == 6421
replace ID = "NE-1987-0019" if ID == "" & year == 1996 & name == "Matthey" & firstname == "Francis" & gdenr_2018 == 6421
replace ID = "AG-1975-0057" if ID == "" & year == 1982 & name == "Mauch" & firstname == "Ursula" & gdenr_2018 == 4073
replace ID = "AG-1975-0057" if ID == "" & year == 1983 & name == "Mauch" & firstname == "Ursula" & gdenr_2018 == 4073
replace ID = "AG-1975-0057" if ID == "" & year == 1984 & name == "Mauch" & firstname == "Ursula" & gdenr_2018 == 4073
replace ID = "AG-1975-0057" if ID == "" & year == 1985 & name == "Mauch" & firstname == "Ursula" & gdenr_2018 == 4073
replace ID = "AG-1975-0057" if ID == "" & year == 1986 & name == "Mauch" & firstname == "Ursula" & gdenr_2018 == 4073
replace ID = "AG-1975-0057" if ID == "" & year == 1987 & name == "Mauch" & firstname == "Ursula" & gdenr_2018 == 4073
replace ID = "ZH-1987-0461" if ID == "" & year == 1991 & name == "Maurer" & firstname == "Ulrich" & gdenr_2018 == 117
replace ID = "ZH-1987-0461" if ID == "" & year == 1992 & name == "Maurer" & firstname == "Ulrich" & gdenr_2018 == 117
replace ID = "ZH-1987-0461" if ID == "" & year == 1993 & name == "Maurer" & firstname == "Ulrich" & gdenr_2018 == 117
replace ID = "ZH-1987-0461" if ID == "" & year == 1994 & name == "Maurer" & firstname == "Ulrich" & gdenr_2018 == 117
replace ID = "ZH-1987-0461" if ID == "" & year == 1996 & name == "Maurer" & firstname == "Ulrich" & gdenr_2018 == 117
replace ID = "ZH-1987-0461" if ID == "" & year == 1997 & name == "Maurer" & firstname == "Ulrich" & gdenr_2018 == 117
replace ID = "ZH-1987-0461" if ID == "" & year == 1998 & name == "Maurer" & firstname == "Ulrich" & gdenr_2018 == 117
replace ID = "ZH-1987-0461" if ID == "" & year == 1999 & name == "Maurer" & firstname == "Ulrich" & gdenr_2018 == 117
replace ID = "ZH-1987-0461" if ID == "" & year == 2000 & name == "Maurer" & firstname == "Ulrich" & gdenr_2018 == 117
replace ID = "ZH-1987-0461" if ID == "" & year == 2001 & name == "Maurer" & firstname == "Ulrich" & gdenr_2018 == 117
replace ID = "GE-1983-9110" if ID == "" & year == 1983 & name == "Maître" & firstname == "Jean Philippe" & gdenr_2018 == 6616
replace ID = "GE-1983-9110" if ID == "" & year == 1985 & name == "Maître" & firstname == "Jean Philippe" & gdenr_2018 == 6621
replace ID = "GE-1983-9110" if ID == "" & year == 1996 & name == "Maître" & firstname == "Jean Philippe" & gdenr_2018 == 6621
replace ID = "GE-1983-9110" if ID == "" & year == 2000 & name == "Maître" & firstname == "Jean Philippe" & gdenr_2018 == 6616
replace ID = "GE-1983-9110" if ID == "" & year == 2001 & name == "Maître" & firstname == "Jean Philippe" & gdenr_2018 == 6616
replace ID = "GE-1983-9110" if ID == "" & year == 2002 & name == "Maître" & firstname == "Jean Philippe" & gdenr_2018 == 6616
replace ID = "GE-1983-9110" if ID == "" & year == 2003 & name == "Maître" & firstname == "Jean Philippe" & gdenr_2018 == 6621
replace ID = "LU-1971-0017" if ID == "" & year == 1982 & name == "Meier" & firstname == " Kaspar" & gdenr_2018 == 1061
replace ID = "ZH-1967-0219" if ID == "" & year == 1983 & name == "Meier" & firstname == "Fritz" & gdenr_2018 == 298
replace ID = "ZH-1967-0219" if ID == "" & year == 1984 & name == "Meier" & firstname == "Fritz" & gdenr_2018 == 298
replace ID = "ZH-1967-0219" if ID == "" & year == 1985 & name == "Meier" & firstname == "Fritz" & gdenr_2018 == 298
replace ID = "BEJU-1971-0253" if ID == "" & year == 1982 & name == "Meier" & firstname == "Werner" & gdenr_2018 == 630
replace ID = "SG-1995-0122" if ID == "" & year == 2000 & name == "Meier Schatz" & firstname == "Lucretia" & gdenr_2018 == 3378
replace ID = "SG-1995-0122" if ID == "" & year == 2001 & name == "Meier Schatz" & firstname == "Lucretia" & gdenr_2018 == 3378
replace ID = "SG-1995-0122" if ID == "" & year == 2002 & name == "Meier Schatz" & firstname == "Lucretia" & gdenr_2018 == 3378
replace ID = "SG-1995-0122" if ID == "" & year == 2003 & name == "Meier Schatz" & firstname == "Lucretia" & gdenr_2018 == 3378
replace ID = "AR-1975-0004" if ID == "" & year == 1982 & name == "Merz" & firstname == "Christian" & gdenr_2018 == 3032
replace ID = "TG-1991-0038" if ID == "" & year == 2000 & name == "Messmer" & firstname == "Werner" & gdenr_2018 == 4506
replace ID = "TG-1991-0038" if ID == "" & year == 2001 & name == "Messmer" & firstname == "Werner" & gdenr_2018 == 4506
replace ID = "TG-1991-0038" if ID == "" & year == 2002 & name == "Messmer" & firstname == "Werner" & gdenr_2018 == 4506
replace ID = "BL-1979-0032" if ID == "" & year == 1992 & name == "Meyer" & firstname == "Theodor" & gdenr_2018 == 2770
replace ID = "BL-1979-0032" if ID == "" & year == 1993 & name == "Meyer" & firstname == "Theodor" & gdenr_2018 == 2770
replace ID = "BL-1979-0032" if ID == "" & year == 1994 & name == "Meyer" & firstname == "Theodor" & gdenr_2018 == 2770
replace ID = "BL-1979-0032" if ID == "" & year == 1996 & name == "Meyer" & firstname == "Theodor" & gdenr_2018 == 2770
replace ID = "BL-1979-0032" if ID == "" & year == 1997 & name == "Meyer" & firstname == "Theodor" & gdenr_2018 == 2770
replace ID = "BL-1979-0032" if ID == "" & year == 1998 & name == "Meyer" & firstname == "Theodor" & gdenr_2018 == 2770
replace ID = "BL-1979-0032" if ID == "" & year == 1999 & name == "Meyer" & firstname == "Theodor" & gdenr_2018 == 2773
replace ID = "FR-1995-0051" if ID == "" & year == 2000 & name == "Meyer Kaelin" & firstname == "Thérèse" & gdenr_2018 == 99999
replace ID = "FR-1995-0051" if ID == "" & year == 2001 & name == "Meyer Kaelin" & firstname == "Thérèse" & gdenr_2018 == 99999
replace ID = "SO-1975-0014" if ID == "" & year == 1991 & name == "Misteli" & firstname == "Marguerite" & gdenr_2018 == 2601
replace ID = "AG-1991-0091" if ID == "" & year == 1991 & name == "Moser Stauber" & firstname == "René O." & gdenr_2018 == 4082
replace ID = "AG-1991-0091" if ID == "" & year == 1992 & name == "Moser Stauber" & firstname == "René O." & gdenr_2018 == 4082
replace ID = "AG-1991-0091" if ID == "" & year == 1993 & name == "Moser Stauber" & firstname == "René O." & gdenr_2018 == 4082
replace ID = "AG-1991-0091" if ID == "" & year == 1994 & name == "Moser Stauber" & firstname == "René O." & gdenr_2018 == 4082
replace ID = "AG-1991-0091" if ID == "" & year == 1996 & name == "Moser Stauber" & firstname == "René O." & gdenr_2018 == 4082
replace ID = "GE-1999-0060" if ID == "" & year == 2000 & name == "Mugny" & firstname == "Patrice" & gdenr_2018 == 99999
replace ID = "GE-1999-0060" if ID == "" & year == 2001 & name == "Mugny" & firstname == "Patrice" & gdenr_2018 == 99999
replace ID = "VD-1971-9031" if ID == "" & year == 2000 & name == "Ménétrey Savary" & firstname == "Anne Catherine" & gdenr_2018 == 99999
replace ID = "VD-1971-9031" if ID == "" & year == 2001 & name == "Ménétrey Savary" & firstname == "Anne Catherine" & gdenr_2018 == 99999
replace ID = "VD-1971-9031" if ID == "" & year == 2002 & name == "Ménétrey Savary" & firstname == "Anne Catherine" & gdenr_2018 == 5601
replace ID = "VD-1971-9031" if ID == "" & year == 2003 & name == "Ménétrey Savary" & firstname == "Anne Catherine" & gdenr_2018 == 5601
replace ID = "ZH-1991-0468" if ID == "" & year == 2000 & name == "Mörgeli" & firstname == "Christoph" & gdenr_2018 == 99999
replace ID = "ZH-1991-0468" if ID == "" & year == 2001 & name == "Mörgeli" & firstname == "Christoph" & gdenr_2018 == 99999
replace ID = "ZH-1991-0468" if ID == "" & year == 2002 & name == "Mörgeli" & firstname == "Christoph" & gdenr_2018 == 99999
replace ID = "ZH-1991-0468" if ID == "" & year == 2003 & name == "Mörgeli" & firstname == "Christoph" & gdenr_2018 == 99999
replace ID = "TG-1983-0027" if ID == "" & year == 1983 & name == "Mühlemann" & firstname == "Ernst" & gdenr_2018 == 4646
replace ID = "ZH-1983-0341" if ID == "" & year == 1983 & name == "Müller" & firstname == "Arnold" & gdenr_2018 == 81
replace ID = "ZH-1979-0402" if ID == "" & year == 1985 & name == "Müller" & firstname == "Kurt" & gdenr_2018 == 99999
replace ID = "ZH-1991-0481" if ID == "" & year == 1996 & name == "Müller Hemmi" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0481" if ID == "" & year == 1997 & name == "Müller Hemmi" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0481" if ID == "" & year == 1998 & name == "Müller Hemmi" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0481" if ID == "" & year == 1999 & name == "Müller Hemmi" & firstname == "Verena" & gdenr_2018 == 261
replace ID = "ZH-1991-0481" if ID == "" & year == 2000 & name == "Müller Hemmi" & firstname == "Verena" & gdenr_2018 == 99999
replace ID = "ZH-1991-0481" if ID == "" & year == 2001 & name == "Müller Hemmi" & firstname == "Verena" & gdenr_2018 == 99999
replace ID = "ZH-1991-0481" if ID == "" & year == 2002 & name == "Müller Hemmi" & firstname == "Verena" & gdenr_2018 == 99999
replace ID = "ZH-1991-0481" if ID == "" & year == 2003 & name == "Müller Hemmi" & firstname == "Verena" & gdenr_2018 == 99999
replace ID = "ZH-1975-0373" if ID == "" & year == 1987 & name == "Nabholz" & firstname == "Lili" & gdenr_2018 == 161
replace ID = "ZH-1975-0373" if ID == "" & year == 1996 & name == "Nabholz Haidegger" & firstname == "Lili" & gdenr_2018 == 161
replace ID = "ZH-1975-0373" if ID == "" & year == 1997 & name == "Nabholz Haidegger" & firstname == "Lili" & gdenr_2018 == 161
replace ID = "ZH-1975-0373" if ID == "" & year == 1998 & name == "Nabholz Haidegger" & firstname == "Lili" & gdenr_2018 == 161
replace ID = "ZH-1975-0373" if ID == "" & year == 1999 & name == "Nabholz Haidegger" & firstname == "Lili" & gdenr_2018 == 161
replace ID = "BL-1967-0016" if ID == "" & year == 1982 & name == "Nebiker Dubach" & firstname == "Hans Rudolf" & gdenr_2018 == 2884
replace ID = "BL-1967-0016" if ID == "" & year == 1983 & name == "Nebiker Dubach" & firstname == "Hans Rudolf" & gdenr_2018 == 2884
replace ID = "BL-1967-0016" if ID == "" & year == 1984 & name == "Nebiker Dubach" & firstname == "Hans Rudolf" & gdenr_2018 == 2884
replace ID = "BL-1967-0016" if ID == "" & year == 1985 & name == "Nebiker Dubach" & firstname == "Hans Rudolf" & gdenr_2018 == 2884
replace ID = "VD-1999-0163" if ID == "" & year == 2000 & name == "Neirynck" & firstname == "Jacques" & gdenr_2018 == 99999
replace ID = "VD-1999-0163" if ID == "" & year == 2001 & name == "Neirynck" & firstname == "Jacques" & gdenr_2018 == 99999
replace ID = "ZH-1967-0245" if ID == "" & year == 1984 & name == "Neuenschwander" & firstname == "Willi" & gdenr_2018 == 246
replace ID = "ZH-1967-0245" if ID == "" & year == 1985 & name == "Neuenschwander" & firstname == "Willi" & gdenr_2018 == 246
replace ID = "ZH-1967-0245" if ID == "" & year == 1986 & name == "Neuenschwander" & firstname == "Willi" & gdenr_2018 == 246
replace ID = "ZH-1967-0245" if ID == "" & year == 1987 & name == "Neuenschwander" & firstname == "Willy" & gdenr_2018 == 246
replace ID = "BL-1983-9011" if ID == "" & year == 1987 & name == "Oberholzer" & firstname == "Susanne" & gdenr_2018 == 2762
replace ID = "BEJU-1971-7032" if ID == "" & year == 1982 & name == "Oehen" & firstname == "Valentin" & gdenr_2018 == 5222
replace ID = "BEJU-1971-7032" if ID == "" & year == 1983 & name == "Oehen" & firstname == "Valentin" & gdenr_2018 == 5222
replace ID = "BEJU-1987-0322" if ID == "" & year == 1996 & name == "Oehrli" & firstname == "Fritz" & gdenr_2018 == 932
replace ID = "BEJU-1987-0322" if ID == "" & year == 1997 & name == "Oehrli" & firstname == "Fritz" & gdenr_2018 == 932
replace ID = "BEJU-1987-0322" if ID == "" & year == 1998 & name == "Oehrli" & firstname == "Fritz" & gdenr_2018 == 932
replace ID = "BEJU-1987-0322" if ID == "" & year == 1999 & name == "Oehrli" & firstname == "Fritz" & gdenr_2018 == 932
replace ID = "BEJU-1987-0322" if ID == "" & year == 2000 & name == "Oehrli" & firstname == "Fritz" & gdenr_2018 == 932
replace ID = "BEJU-1987-0322" if ID == "" & year == 2001 & name == "Oehrli" & firstname == "Fritz" & gdenr_2018 == 932
replace ID = "BEJU-1987-0322" if ID == "" & year == 2002 & name == "Oehrli" & firstname == "Fritz" & gdenr_2018 == 932
replace ID = "BEJU-1987-0322" if ID == "" & year == 2003 & name == "Oehrli" & firstname == "Fritz" & gdenr_2018 == 932
replace ID = "ZH-1963-0220" if ID == "" & year == 1982 & name == "Oester" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "ZH-1963-0220" if ID == "" & year == 1983 & name == "Oester" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "ZH-1963-0220" if ID == "" & year == 1984 & name == "Oester" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "ZH-1963-0220" if ID == "" & year == 1985 & name == "Oester" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "ZH-1963-0220" if ID == "" & year == 1986 & name == "Oester" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "ZH-1963-0220" if ID == "" & year == 1987 & name == "Oester" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "ZH-1963-0220" if ID == "" & year == 1988 & name == "Oester" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "ZH-1963-0220" if ID == "" & year == 1989 & name == "Oester" & firstname == "Hans" & gdenr_2018 == 261
replace ID = "BEJU-1979-0254" if ID == "" & year == 1984 & name == "Ogi" & firstname == "Adolf" & gdenr_2018 == 538
replace ID = "BEJU-1979-0254" if ID == "" & year == 1985 & name == "Ogi" & firstname == "Adolf" & gdenr_2018 == 538
replace ID = "BEJU-1979-0254" if ID == "" & year == 1986 & name == "Ogi" & firstname == "Adolf" & gdenr_2018 == 538
replace ID = "BEJU-1979-0254" if ID == "" & year == 1987 & name == "Ogi" & firstname == "Adolf" & gdenr_2018 == 538
replace ID = "BL-1963-0016" if ID == "" & year == 1982 & name == "Ott" & firstname == "Heinrich" & gdenr_2018 == 2769
replace ID = "BL-1963-0016" if ID == "" & year == 1983 & name == "Ott" & firstname == "Heinrich" & gdenr_2018 == 2769
replace ID = "BL-1963-0016" if ID == "" & year == 1985 & name == "Ott" & firstname == "Heinrich" & gdenr_2018 == 2767
replace ID = "BL-1963-0016" if ID == "" & year == 1986 & name == "Ott" & firstname == "Heinrich" & gdenr_2018 == 2767
replace ID = "BL-1963-0016" if ID == "" & year == 1987 & name == "Ott" & firstname == "Heinrich" & gdenr_2018 == 2767
replace ID = "VS-1987-0027" if ID == "" & year == 1987 & name == "Paccolat" & firstname == "Monique" & gdenr_2018 == 6211
replace ID = "TI-1995-0049" if ID == "" & year == 1996 & name == "Pelly" & firstname == "Fulvio" & gdenr_2018 == 5192
replace ID = "VD-1983-0125" if ID == "" & year == 1983 & name == "Perey" & firstname == "André" & gdenr_2018 == 5653
replace ID = "FR-1987-0029" if ID == "" & year == 1987 & name == "Philipona" & firstname == "Jean Nicolas" & gdenr_2018 == 2140
replace ID = "VD-1983-0128" if ID == "" & year == 1983 & name == "Pidoux" & firstname == "Philippe" & gdenr_2018 == 5586
replace ID = "TI-1975-9062" if ID == "" & year == 1993 & name == "Pini" & firstname == "Massimo" & gdenr_2018 == 5281
replace ID = "TI-1975-9062" if ID == "" & year == 1994 & name == "Pini" & firstname == "Massimo" & gdenr_2018 == 5281
replace ID = "TI-1975-9062" if ID == "" & year == 1998 & name == "Pini" & firstname == "Massimo" & gdenr_2018 == 5398
replace ID = "TI-1975-9062" if ID == "" & year == 1999 & name == "Pini" & firstname == "Massimo" & gdenr_2018 == 5398
replace ID = "GE-1991-0084" if ID == "" & year == 1991 & name == "Poncet" & firstname == "Charles" & gdenr_2018 == 6642
replace ID = "GR-1987-0026" if ID == "" & year == 1988 & name == "Portmann" & firstname == "Theodor C." & gdenr_2018 == 3901
replace ID = "GR-1987-0026" if ID == "" & year == 1989 & name == "Portmann" & firstname == "Theodor C." & gdenr_2018 == 3901
replace ID = "GR-1987-0026" if ID == "" & year == 1990 & name == "Portmann" & firstname == "Theodor C." & gdenr_2018 == 3901
replace ID = "GR-1987-0026" if ID == "" & year == 1987 & name == "Portmann" & firstname == "Theodor Christian" & gdenr_2018 == 3901
replace ID = "TG-1991-0050" if ID == "" & year == 1998 & name == "Raggenbass" & firstname == "Hans Ulrich A." & gdenr_2018 == 4426
replace ID = "TG-1991-0050" if ID == "" & year == 1999 & name == "Raggenbass" & firstname == "Hans Ulrich A." & gdenr_2018 == 4426
replace ID = "TG-1991-0050" if ID == "" & year == 2000 & name == "Raggenbass" & firstname == "Hans Ulrich A." & gdenr_2018 == 4461
replace ID = "TG-1991-0050" if ID == "" & year == 2001 & name == "Raggenbass" & firstname == "Hans Ulrich A." & gdenr_2018 == 4461
replace ID = "TG-1991-0050" if ID == "" & year == 2002 & name == "Raggenbass" & firstname == "Hans Ulrich A." & gdenr_2018 == 4461
replace ID = "TG-1991-0050" if ID == "" & year == 2003 & name == "Raggenbass" & firstname == "Hans Ulrich A." & gdenr_2018 == 4461
replace ID = "TG-1991-0050" if ID == "" & year == 1991 & name == "Raggenbass" & firstname == "Hansueli" & gdenr_2018 == 4426
replace ID = "BS-1995-0050" if ID == "" & year == 1996 & name == "Randegger" & firstname == "Johannes R." & gdenr_2018 == 2702
replace ID = "BS-1995-0050" if ID == "" & year == 1997 & name == "Randegger" & firstname == "Johannes R." & gdenr_2018 == 2702
replace ID = "BS-1995-0050" if ID == "" & year == 1998 & name == "Randegger" & firstname == "Johannes R." & gdenr_2018 == 2702
replace ID = "BS-1995-0050" if ID == "" & year == 1999 & name == "Randegger" & firstname == "Johannes R." & gdenr_2018 == 2702
replace ID = "BS-1995-0050" if ID == "" & year == 2000 & name == "Randegger" & firstname == "Johannes R." & gdenr_2018 == 2702
replace ID = "BS-1995-0050" if ID == "" & year == 2001 & name == "Randegger" & firstname == "Johannes R." & gdenr_2018 == 2702
replace ID = "BS-1995-0050" if ID == "" & year == 2002 & name == "Randegger" & firstname == "Johannes R." & gdenr_2018 == 2702
replace ID = "BS-1995-0050" if ID == "" & year == 2003 & name == "Randegger" & firstname == "Johannes R." & gdenr_2018 == 2702
replace ID = "GE-1983-0034" if ID == "" & year == 1983 & name == "Rebeaud" & firstname == "Laurent" & gdenr_2018 == 6643
replace ID = "GE-1983-0034" if ID == "" & year == 1988 & name == "Rebeaud" & firstname == "Laurent" & gdenr_2018 == 6643
replace ID = "SG-1983-0041" if ID == "" & year == 2000 & name == "Rechsteiner" & firstname == "Paul" & gdenr_2018 == 99999
replace ID = "SG-1983-0041" if ID == "" & year == 2001 & name == "Rechsteiner" & firstname == "Paul" & gdenr_2018 == 99999
replace ID = "BS-1991-0050" if ID == "" & year == 1996 & name == "Rechsteiner" & firstname == "Rudolf" & gdenr_2018 == 2701
replace ID = "BS-1991-0050" if ID == "" & year == 1997 & name == "Rechsteiner" & firstname == "Rudolf" & gdenr_2018 == 2701
replace ID = "BS-1991-0050" if ID == "" & year == 1998 & name == "Rechsteiner" & firstname == "Rudolf" & gdenr_2018 == 2701
replace ID = "BS-1991-0050" if ID == "" & year == 1999 & name == "Rechsteiner" & firstname == "Rudolf" & gdenr_2018 == 2701
replace ID = "BS-1991-0050" if ID == "" & year == 2000 & name == "Rechsteiner" & firstname == "Rudolf" & gdenr_2018 == 2701
replace ID = "BS-1991-0050" if ID == "" & year == 2001 & name == "Rechsteiner" & firstname == "Rudolf" & gdenr_2018 == 2701
replace ID = "BS-1991-0050" if ID == "" & year == 2002 & name == "Rechsteiner" & firstname == "Rudolf" & gdenr_2018 == 2701
replace ID = "BS-1991-0050" if ID == "" & year == 2003 & name == "Rechsteiner" & firstname == "Rudolf" & gdenr_2018 == 2701
replace ID = "SG-1983-0041" if ID == "" & year == 1986 & name == "Reichsteiner" & firstname == "Paul" & gdenr_2018 == 3203
replace ID = "SG-1983-0041" if ID == "" & year == 1987 & name == "Reichsteiner" & firstname == "Paul" & gdenr_2018 == 3203
replace ID = "SG-1983-0041" if ID == "" & year == 1988 & name == "Reichsteiner" & firstname == "Paul" & gdenr_2018 == 3203
replace ID = "SG-1983-0041" if ID == "" & year == 1989 & name == "Reichsteiner" & firstname == "Paul" & gdenr_2018 == 3203
replace ID = "SG-1983-0041" if ID == "" & year == 1990 & name == "Reichsteiner" & firstname == "Paul" & gdenr_2018 == 3203
replace ID = "SG-1983-0041" if ID == "" & year == 1991 & name == "Reichsteiner" & firstname == "Paul" & gdenr_2018 == 3203
replace ID = "SG-1983-0041" if ID == "" & year == 1992 & name == "Reichsteiner" & firstname == "Paul" & gdenr_2018 == 3203
replace ID = "SG-1983-0041" if ID == "" & year == 1993 & name == "Reichsteiner" & firstname == "Paul" & gdenr_2018 == 3203
replace ID = "SG-1983-0041" if ID == "" & year == 1994 & name == "Reichsteiner" & firstname == "Paul" & gdenr_2018 == 3203
replace ID = "AG-1987-0107" if ID == "" & year == 1987 & name == "Reimann" & firstname == "Maximilian" & gdenr_2018 == 4165
replace ID = "ZH-1971-0340" if ID == "" & year == 1982 & name == "Ribi" & firstname == "Martha" & gdenr_2018 == 261
replace ID = "ZH-1975-0416" if ID == "" & year == 2000 & name == "Riklin" & firstname == "Katharina" & gdenr_2018 == 99999
replace ID = "ZH-1975-0416" if ID == "" & year == 2001 & name == "Riklin" & firstname == "Katharina" & gdenr_2018 == 261
replace ID = "ZH-1975-0416" if ID == "" & year == 2002 & name == "Riklin" & firstname == "Katharina" & gdenr_2018 == 261
replace ID = "ZH-1975-0416" if ID == "" & year == 2003 & name == "Riklin" & firstname == "Katharina" & gdenr_2018 == 261
replace ID = "TI-1979-9061" if ID == "" & year == 1982 & name == "Robbiani" & firstname == "Dario" & gdenr_2018 == 5176
replace ID = "TI-1979-9061" if ID == "" & year == 1983 & name == "Robbiani" & firstname == "Dario" & gdenr_2018 == 5176
replace ID = "TI-1999-0053" if ID == "" & year == 2001 & name == "Robbiani" & firstname == "Meinrado" & gdenr_2018 == 5171
replace ID = "TI-1999-0053" if ID == "" & year == 2002 & name == "Robbiani" & firstname == "Meinrado" & gdenr_2018 == 5171
replace ID = "TI-1999-0053" if ID == "" & year == 2003 & name == "Robbiani" & firstname == "Meinrado" & gdenr_2018 == 5171
replace ID = "BEJU-1971-0313" if ID == "" & year == 1991 & name == "Robert Bächtold" & firstname == "Leni" & gdenr_2018 == 351
replace ID = "FR-1987-9000" if ID == "" & year == 1987 & name == "Rohrbasser" & firstname == "Bernard" & gdenr_2018 == 2325
replace ID = "SG-1971-0063" if ID == "" & year == 1985 & name == "Rohrer" & firstname == "Hans" & gdenr_2018 == 3273
replace ID = "VS-1995-0048" if ID == "" & year == 2000 & name == "Rossini" & firstname == "Stéphane" & gdenr_2018 == 99999
replace ID = "VS-1995-0048" if ID == "" & year == 2001 & name == "Rossini" & firstname == "Stéphane" & gdenr_2018 == 99999
replace ID = "VS-1995-0048" if ID == "" & year == 2002 & name == "Rossini" & firstname == "Stéphane" & gdenr_2018 == 99999
replace ID = "VS-1995-0048" if ID == "" & year == 2003 & name == "Rossini" & firstname == "Stéphane" & gdenr_2018 == 99999
replace ID = "AG-1963-0057" if ID == "" & year == 1982 & name == "Roth" & firstname == "Hans" & gdenr_2018 == 2503
replace ID = "SG-1983-0043" if ID == "" & year == 1983 & name == "Ruckstuhl" & firstname == "Hans" & gdenr_2018 == 3427
replace ID = "BEJU-1979-0292" if ID == "" & year == 1983 & name == "Ruf" & firstname == "Markus" & gdenr_2018 == 351
replace ID = "BEJU-1983-0323" if ID == "" & year == 1994 & name == "Rychen" & firstname == "Albert" & gdenr_2018 == 306
replace ID = "BEJU-1983-0323" if ID == "" & year == 1996 & name == "Rychen" & firstname == "Albert" & gdenr_2018 == 306
replace ID = "BEJU-1983-0323" if ID == "" & year == 1997 & name == "Rychen" & firstname == "Albert" & gdenr_2018 == 306
replace ID = "BEJU-1983-0323" if ID == "" & year == 1998 & name == "Rychen" & firstname == "Albert" & gdenr_2018 == 306
replace ID = "BEJU-1983-0323" if ID == "" & year == 1999 & name == "Rychen" & firstname == "Albert" & gdenr_2018 == 306
replace ID = "BEJU-1971-0328" if ID == "" & year == 1982 & name == "Räz Rutsch" & firstname == "Fritz" & gdenr_2018 == 310
replace ID = "OW-1971-0002" if ID == "" & year == 1982 & name == "Röthlin Lieb" & firstname == "Walter" & gdenr_2018 == 1404
replace ID = "OW-1971-0002" if ID == "" & year == 1983 & name == "Röthlin Lieb" & firstname == "Walter" & gdenr_2018 == 1404
replace ID = "OW-1971-0002" if ID == "" & year == 1984 & name == "Röthlin Lieb" & firstname == "Walter" & gdenr_2018 == 1404
replace ID = "OW-1971-0002" if ID == "" & year == 1985 & name == "Röthlin Lieb" & firstname == "Walter" & gdenr_2018 == 1404
replace ID = "OW-1971-0002" if ID == "" & year == 1986 & name == "Röthlin Lieb" & firstname == "Walter" & gdenr_2018 == 1404
replace ID = "TI-1983-0046" if ID == "" & year == 1983 & name == "Salvioni" & firstname == "Sergio" & gdenr_2018 == 5113
replace ID = "TI-1983-0046" if ID == "" & year == 1984 & name == "Salvioni" & firstname == "Sergio" & gdenr_2018 == 5113
replace ID = "VD-1995-0174" if ID == "" & year == 1996 & name == "Sandoz" & firstname == "Marcel" & gdenr_2018 == 5639
replace ID = "FR-1983-0026" if ID == "" & year == 1983 & name == "Savary" & firstname == "Jean" & gdenr_2018 == 2113
replace ID = "FR-1983-0026" if ID == "" & year == 1984 & name == "Savary" & firstname == "Jean" & gdenr_2018 == 2113
replace ID = "FR-1983-0026" if ID == "" & year == 1985 & name == "Savary" & firstname == "Jean" & gdenr_2018 == 2113
replace ID = "FR-1983-0026" if ID == "" & year == 1986 & name == "Savary" & firstname == "Jean" & gdenr_2018 == 2113
replace ID = "FR-1983-0026" if ID == "" & year == 1987 & name == "Savary" & firstname == "Jean" & gdenr_2018 == 2113
replace ID = "FR-1983-0026" if ID == "" & year == 1988 & name == "Savary" & firstname == "Jean" & gdenr_2018 == 2113
replace ID = "FR-1983-0026" if ID == "" & year == 1989 & name == "Savary" & firstname == "Jean" & gdenr_2018 == 2113
replace ID = "FR-1983-0026" if ID == "" & year == 1990 & name == "Savary" & firstname == "Jean" & gdenr_2018 == 2113
replace ID = "SO-1987-9005" if ID == "" & year == 1987 & name == "Scheidegger" & firstname == "Urs" & gdenr_2018 == 2601
replace ID = "BEJU-1987-0383" if ID == "" & year == 1987 & name == "Scherrer" & firstname == "Jürg" & gdenr_2018 == 303
replace ID = "NE-1991-0025" if ID == "" & year == 1991 & name == "Scheurer" & firstname == "Rémy" & gdenr_2018 == 6454
replace ID = "VS-1991-0030" if ID == "" & year == 1996 & name == "Schmid Zimmermann" & firstname == "Odilo" & gdenr_2018 == 6002
replace ID = "BEJU-1991-0437" if ID == "" & year == 1991 & name == "Schmied" & firstname == "Walter" & gdenr_2018 == 700
replace ID = "BEJU-1995-0453" if ID == "" & year == 2000 & name == "Schneider" & firstname == "Johann N." & gdenr_2018 == 329
replace ID = "VD-1987-0155" if ID == "" & year == 2000 & name == "Schwaab" & firstname == "Jean Jaques" & gdenr_2018 == 5613
replace ID = "VD-1987-0155" if ID == "" & year == 2001 & name == "Schwaab" & firstname == "Jean Jaques" & gdenr_2018 == 5613
replace ID = "VD-1987-0155" if ID == "" & year == 2002 & name == "Schwaab" & firstname == "Jean Jaques" & gdenr_2018 == 5613
replace ID = "VD-1987-0155" if ID == "" & year == 2003 & name == "Schwaab" & firstname == "Jean Jaques" & gdenr_2018 == 5613
replace ID = "BEJU-1987-0405" if ID == "" & year == 1987 & name == "Schwab" & firstname == "Heinz" & gdenr_2018 == 312
replace ID = "AG-1967-9100" if ID == "" & year == 1982 & name == "Schwarz Clavadetscher" & firstname == "Urs" & gdenr_2018 == 4289
replace ID = "AG-1967-9100" if ID == "" & year == 1983 & name == "Schwarz Clavadetscher" & firstname == "Urs" & gdenr_2018 == 4289
replace ID = "AG-1967-9100" if ID == "" & year == 1984 & name == "Schwarz Clavadetscher" & firstname == "Urs" & gdenr_2018 == 4289
replace ID = "AG-1967-9100" if ID == "" & year == 1985 & name == "Schwarz Clavadetscher" & firstname == "Urs" & gdenr_2018 == 4289
replace ID = "LU-1975-0040" if ID == "" & year == 1982 & name == "Schärfli" & firstname == "Hans" & gdenr_2018 == 1143
replace ID = "LU-1975-0040" if ID == "" & year == 1983 & name == "Schärli Ludin" & firstname == "Hans" & gdenr_2018 == 1143
replace ID = "LU-1975-0040" if ID == "" & year == 1984 & name == "Schärli Ludin" & firstname == "Hans" & gdenr_2018 == 1143
replace ID = "LU-1975-0040" if ID == "" & year == 1985 & name == "Schärli Ludin" & firstname == "Hans" & gdenr_2018 == 1143
replace ID = "SG-1979-0045" if ID == "" & year == 1982 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1983 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1984 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1985 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1986 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1987 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1988 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1989 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1990 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1991 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1992 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1993 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "SG-1979-0045" if ID == "" & year == 1994 & name == "Segmüller" & firstname == "Eva" & gdenr_2018 == 3203
replace ID = "GE-1971-9101" if ID == "" & year == 1987 & name == "Segond" & firstname == "Guy Olivier" & gdenr_2018 == 6621
replace ID = "BEJU-1987-0414" if ID == "" & year == 1987 & name == "Seiler" & firstname == "Hanspeter" & gdenr_2018 == 590
replace ID = "GR-1995-0042" if ID == "" & year == 1996 & name == "Semedani" & firstname == "Silva A." & gdenr_2018 == 3911
replace ID = "GR-1995-0042" if ID == "" & year == 1997 & name == "Semedani" & firstname == "Silva A." & gdenr_2018 == 3911
replace ID = "GR-1995-0042" if ID == "" & year == 1998 & name == "Semedani" & firstname == "Silva A." & gdenr_2018 == 3911
replace ID = "GR-1995-0042" if ID == "" & year == 1999 & name == "Semedani" & firstname == "Silva A." & gdenr_2018 == 3911
replace ID = "ZH-1991-0640" if ID == "" & year == 1991 & name == "Sieber" & firstname == "Ernst" & gdenr_2018 == 261
replace ID = "VD-1995-0181" if ID == "" & year == 1997 & name == "Simon" & firstname == "Jean Charles" & gdenr_2018 == 5606
replace ID = "VD-1995-0181" if ID == "" & year == 1998 & name == "Simon" & firstname == "Jean Charles" & gdenr_2018 == 5606
replace ID = "VD-1995-0181" if ID == "" & year == 1999 & name == "Simon" & firstname == "Jean Charles" & gdenr_2018 == 5606
replace ID = "GE-1967-9224" if ID == "" & year == 1983 & name == "Soldini" & firstname == "Mario" & gdenr_2018 == 6621
replace ID = "GE-1967-9224" if ID == "" & year == 1984 & name == "Soldini" & firstname == "Mario" & gdenr_2018 == 6621
replace ID = "BEJU-1999-0371" if ID == "" & year == 2000 & name == "Sommaruga" & firstname == "Simonetta" & gdenr_2018 == 99999
replace ID = "BEJU-1999-0371" if ID == "" & year == 2001 & name == "Sommaruga" & firstname == "Simonetta" & gdenr_2018 == 99999
replace ID = "ZH-1983-0471" if ID == "" & year == 1985 & name == "Spaelti" & firstname == "Peter" & gdenr_2018 == 221
replace ID = "ZH-1983-0467" if ID == "" & year == 1989 & name == "Spoerry Toneatti" & firstname == "Verena" & gdenr_2018 == 295
replace ID = "ZH-1983-0467" if ID == "" & year == 1990 & name == "Spoerry Toneatti" & firstname == "Verena" & gdenr_2018 == 295
replace ID = "ZH-1983-0467" if ID == "" & year == 1991 & name == "Spoerry Toneatti" & firstname == "Verena" & gdenr_2018 == 295
replace ID = "ZH-1983-0467" if ID == "" & year == 1992 & name == "Spoerry Toneatti" & firstname == "Verena" & gdenr_2018 == 295
replace ID = "ZH-1983-0467" if ID == "" & year == 1993 & name == "Spoerry Toneatti" & firstname == "Verena" & gdenr_2018 == 295
replace ID = "ZH-1983-0467" if ID == "" & year == 1994 & name == "Spoerry Toneatti" & firstname == "Verena" & gdenr_2018 == 295
replace ID = "ZH-1983-0467" if ID == "" & year == 1983 & name == "Spoerry Toneatti" & firstname == "Vreni" & gdenr_2018 == 295
replace ID = "TG-1999-0075" if ID == "" & year == 2000 & name == "Spuhler" & firstname == "Peter Chr." & gdenr_2018 == 4621
replace ID = "TG-1999-0075" if ID == "" & year == 2001 & name == "Spuhler" & firstname == "Peter Chr." & gdenr_2018 == 4621
replace ID = "TG-1999-0075" if ID == "" & year == 2002 & name == "Spuhler" & firstname == "Peter Chr." & gdenr_2018 == 4621
replace ID = "TG-1999-0075" if ID == "" & year == 2003 & name == "Spuhler" & firstname == "Peter Chr." & gdenr_2018 == 4621
replace ID = "ZH-1983-0471" if ID == "" & year == 1983 & name == "Spälti" & firstname == "Peter" & gdenr_2018 == 221
replace ID = "LU-1983-0036" if ID == "" & year == 1983 & name == "Stamm" & firstname == "Judith" & gdenr_2018 == 1061
replace ID = "AG-1991-0129" if ID == "" & year == 1991 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 4021
replace ID = "AG-1991-0129" if ID == "" & year == 1992 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 4021
replace ID = "AG-1991-0129" if ID == "" & year == 1993 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 4021
replace ID = "AG-1991-0129" if ID == "" & year == 1994 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 4021
replace ID = "AG-1991-0129" if ID == "" & year == 1996 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 4021
replace ID = "AG-1991-0129" if ID == "" & year == 1997 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 4021
replace ID = "AG-1991-0129" if ID == "" & year == 1998 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 4021
replace ID = "AG-1991-0129" if ID == "" & year == 1999 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 4021
replace ID = "AG-1991-0129" if ID == "" & year == 2000 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 99999
replace ID = "AG-1991-0129" if ID == "" & year == 2001 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 99999
replace ID = "AG-1991-0129" if ID == "" & year == 2002 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 4021
replace ID = "AG-1991-0129" if ID == "" & year == 2003 & name == "Stamm" & firstname == "Luzius" & gdenr_2018 == 4021
replace ID = "ZH-1971-0395" if ID == "" & year == 1986 & name == "Stappung" & firstname == "Josef" & gdenr_2018 == 247
replace ID = "ZH-1971-0398" if ID == "" & year == 1988 & name == "Steffen" & firstname == "Hans" & gdenr_2018 == 114
replace ID = "ZH-1971-0398" if ID == "" & year == 1989 & name == "Steffen" & firstname == "Hans" & gdenr_2018 == 114
replace ID = "ZH-1971-0398" if ID == "" & year == 1990 & name == "Steffen" & firstname == "Hans" & gdenr_2018 == 114
replace ID = "ZH-1971-0398" if ID == "" & year == 1991 & name == "Steffen" & firstname == "Hans" & gdenr_2018 == 114
replace ID = "ZH-1971-0398" if ID == "" & year == 1996 & name == "Steffen" & firstname == "Hans" & gdenr_2018 == 114
replace ID = "ZH-1971-0398" if ID == "" & year == 1997 & name == "Steffen" & firstname == "Hans" & gdenr_2018 == 114
replace ID = "ZH-1971-0398" if ID == "" & year == 1998 & name == "Steffen" & firstname == "Hans" & gdenr_2018 == 114
replace ID = "ZH-1971-0398" if ID == "" & year == 1999 & name == "Steffen" & firstname == "Hans" & gdenr_2018 == 114
replace ID = "UR-1983-0002" if ID == "" & year == 1982 & name == "Steinegger" & firstname == "Franz" & gdenr_2018 == 1207
replace ID = "UR-1983-0002" if ID == "" & year == 1983 & name == "Steinegger" & firstname == "Franz" & gdenr_2018 == 1207
replace ID = "SO-1991-0053" if ID == "" & year == 1994 & name == "Steiner" & firstname == "Rudolph" & gdenr_2018 == 2493
replace ID = "SO-1991-0053" if ID == "" & year == 1996 & name == "Steiner" & firstname == "Rudolph" & gdenr_2018 == 2493
replace ID = "SO-1991-0053" if ID == "" & year == 1997 & name == "Steiner" & firstname == "Rudolph" & gdenr_2018 == 2493
replace ID = "SO-1991-0053" if ID == "" & year == 1998 & name == "Steiner" & firstname == "Rudolph" & gdenr_2018 == 2493
replace ID = "SO-1991-0053" if ID == "" & year == 1999 & name == "Steiner" & firstname == "Rudolph" & gdenr_2018 == 2493
replace ID = "SO-1991-0053" if ID == "" & year == 2000 & name == "Steiner" & firstname == "Rudolph" & gdenr_2018 == 2493
replace ID = "SO-1991-0053" if ID == "" & year == 2001 & name == "Steiner" & firstname == "Rudolph" & gdenr_2018 == 2493
replace ID = "SO-1991-0053" if ID == "" & year == 2002 & name == "Steiner" & firstname == "Rudolph" & gdenr_2018 == 2493
replace ID = "SO-1991-0053" if ID == "" & year == 2003 & name == "Steiner" & firstname == "Rudolph" & gdenr_2018 == 2493
replace ID = "ZH-1987-9254" if ID == "" & year == 1987 & name == "Stocker" & firstname == "Monika" & gdenr_2018 == 261
replace ID = "ZH-1987-9254" if ID == "" & year == 1988 & name == "Stocker" & firstname == "Monika" & gdenr_2018 == 261
replace ID = "ZH-1987-9254" if ID == "" & year == 1989 & name == "Stocker" & firstname == "Monika" & gdenr_2018 == 261
replace ID = "ZH-1987-9254" if ID == "" & year == 1990 & name == "Stocker" & firstname == "Monika" & gdenr_2018 == 261
replace ID = "BEJU-1987-0450" if ID == "" & year == 1991 & name == "Strahm" & firstname == "Rudolf H." & gdenr_2018 == 354
replace ID = "AG-1971-0115" if ID == "" & year == 2000 & name == "Studer" & firstname == "Heinz" & gdenr_2018 == 4045
replace ID = "AG-1971-0115" if ID == "" & year == 2001 & name == "Studer" & firstname == "Heinz" & gdenr_2018 == 4045
replace ID = "AG-1991-0135" if ID == "" & year == 1996 & name == "Stumpf" & firstname == "Doris" & gdenr_2018 == 4045
replace ID = "AG-1991-0135" if ID == "" & year == 1997 & name == "Stumpf" & firstname == "Doris" & gdenr_2018 == 4045
replace ID = "AG-1991-0135" if ID == "" & year == 1998 & name == "Stumpf" & firstname == "Doris" & gdenr_2018 == 4045
replace ID = "AG-1991-0135" if ID == "" & year == 1999 & name == "Stumpf" & firstname == "Doris" & gdenr_2018 == 4045
replace ID = "BEJU-1987-0461" if ID == "" & year == 1991 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 1992 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 1993 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 1994 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 1996 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 1997 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 1998 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 1999 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 2000 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 2001 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 2002 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "BEJU-1987-0461" if ID == "" & year == 2003 & name == "Suter" & firstname == "Marc F." & gdenr_2018 == 371
replace ID = "AG-1987-0139" if ID == "" & year == 1987 & name == "Thür" & firstname == "Hanspeter" & gdenr_2018 == 4001
replace ID = "VD-1987-9011" if ID == "" & year == 2000 & name == "Tillmann" & firstname == "Pierre" & gdenr_2018 == 5586
replace ID = "VD-1987-9011" if ID == "" & year == 2001 & name == "Tillmann" & firstname == "Pierre" & gdenr_2018 == 5586
replace ID = "GE-1979-0045" if ID == "" & year == 1991 & name == "Tschopp" & firstname == "Peter" & gdenr_2018 == 6642
replace ID = "LU-1983-0038" if ID == "" & year == 1983 & name == "Tschuppert" & firstname == "Karl" & gdenr_2018 == 1128
replace ID = "BEJU-1991-0510" if ID == "" & year == 1991 & name == "Tschäppät" & firstname == "Alexander" & gdenr_2018 == 351
replace ID = "ZH-1971-0414" if ID == "" & year == 1984 & name == "Uchtenhagen" & firstname == "Lilian" & gdenr_2018 == 261
replace ID = "ZH-1971-0414" if ID == "" & year == 1985 & name == "Uchtenhagen" & firstname == "Lilian" & gdenr_2018 == 261
replace ID = "ZH-1971-0414" if ID == "" & year == 1986 & name == "Uchtenhagen" & firstname == "Lilian" & gdenr_2018 == 261
replace ID = "ZH-1971-0414" if ID == "" & year == 1987 & name == "Uchtenhagen" & firstname == "Lilian" & gdenr_2018 == 261
replace ID = "ZH-1971-0414" if ID == "" & year == 1988 & name == "Uchtenhagen" & firstname == "Lilian" & gdenr_2018 == 261
replace ID = "ZH-1971-0414" if ID == "" & year == 1989 & name == "Uchtenhagen" & firstname == "Lilian" & gdenr_2018 == 261
replace ID = "ZH-1971-0414" if ID == "" & year == 1990 & name == "Uchtenhagen" & firstname == "Lilian" & gdenr_2018 == 261
replace ID = "TG-1983-0043" if ID == "" & year == 1983 & name == "Uhlmann" & firstname == "Hans" & gdenr_2018 == 4951
replace ID = "SO-1983-0037" if ID == "" & year == 1987 & name == "Ulrich" & firstname == "Ursula" & gdenr_2018 == 2581
replace ID = "AR-1995-0009" if ID == "" & year == 1996 & name == "Vallender" & firstname == "Dorothea" & gdenr_2018 == 3025
replace ID = "AR-1995-0009" if ID == "" & year == 1997 & name == "Vallender" & firstname == "Dorothea" & gdenr_2018 == 3025
replace ID = "AR-1995-0009" if ID == "" & year == 1998 & name == "Vallender" & firstname == "Dorothea" & gdenr_2018 == 3025
replace ID = "AR-1995-0009" if ID == "" & year == 1999 & name == "Vallender" & firstname == "Dorothea" & gdenr_2018 == 3025
replace ID = "AR-1995-0009" if ID == "" & year == 2000 & name == "Vallender" & firstname == "Dorothea" & gdenr_2018 == 3025
replace ID = "AR-1995-0009" if ID == "" & year == 2001 & name == "Vallender" & firstname == "Dorothea" & gdenr_2018 == 3025
replace ID = "AR-1995-0009" if ID == "" & year == 2002 & name == "Vallender" & firstname == "Dorothea" & gdenr_2018 == 3025
replace ID = "AR-1995-0009" if ID == "" & year == 2003 & name == "Vallender Clausen" & firstname == "Dorothea" & gdenr_2018 == 3025
replace ID = "VS-1979-0042" if ID == "" & year == 1982 & name == "Vannay" & firstname == "Françoise" & gdenr_2018 == 6217
replace ID = "VS-1979-0042" if ID == "" & year == 1983 & name == "Vannay" & firstname == "Françoise" & gdenr_2018 == 6217
replace ID = "VS-1979-0042" if ID == "" & year == 1984 & name == "Vannay" & firstname == "Françoise" & gdenr_2018 == 6217
replace ID = "VS-1979-0042" if ID == "" & year == 1985 & name == "Vannay" & firstname == "Françoise" & gdenr_2018 == 99999
replace ID = "VS-1979-0042" if ID == "" & year == 1986 & name == "Vannay" & firstname == "Françoise" & gdenr_2018 == 6217
replace ID = "ZH-1991-0724" if ID == "" & year == 1991 & name == "Vetterli" & firstname == "Werner" & gdenr_2018 == 248
replace ID = "BEJU-1983-0398" if ID == "" & year == 1997 & name == "Waber" & firstname == "Christian" & gdenr_2018 == 99999
replace ID = "BEJU-1983-0398" if ID == "" & year == 1998 & name == "Waber" & firstname == "Christian" & gdenr_2018 == 99999
replace ID = "BEJU-1995-0541" if ID == "" & year == 2000 & name == "Wandfluh" & firstname == "Hansruedi" & gdenr_2018 == 99999
replace ID = "BEJU-1995-0541" if ID == "" & year == 2001 & name == "Wandfluh" & firstname == "Hansruedi" & gdenr_2018 == 99999
replace ID = "BEJU-1995-0541" if ID == "" & year == 2002 & name == "Wandfluh" & firstname == "Hansruedi" & gdenr_2018 == 99999
replace ID = "SO-1975-0029" if ID == "" & year == 1983 & name == "Wanner" & firstname == "Christian" & gdenr_2018 == 2457
replace ID = "AG-1991-0147" if ID == "" & year == 1996 & name == "Weber" & firstname == "Agnes" & gdenr_2018 == 99999
replace ID = "AG-1991-0147" if ID == "" & year == 1997 & name == "Weber" & firstname == "Agnes" & gdenr_2018 == 99999
replace ID = "AG-1991-0147" if ID == "" & year == 1998 & name == "Weber" & firstname == "Agnes" & gdenr_2018 == 99999
replace ID = "AG-1955-0065" if ID == "" & year == 1982 & name == "Weber Huber" & firstname == "Leo" & gdenr_2018 == 4236
replace ID = "AG-1955-0065" if ID == "" & year == 1983 & name == "Weber Huber" & firstname == "Leo" & gdenr_2018 == 4236
replace ID = "AG-1955-0065" if ID == "" & year == 1984 & name == "Weber Huber" & firstname == "Leo" & gdenr_2018 == 4236
replace ID = "AG-1955-0065" if ID == "" & year == 1985 & name == "Weber Huber" & firstname == "Leo" & gdenr_2018 == 4236
replace ID = "SZ-1975-9014" if ID == "" & year == 1982 & name == "Weber Wiget" & firstname == "Karl" & gdenr_2018 == 1372
replace ID = "SZ-1975-9014" if ID == "" & year == 1983 & name == "Weber Wiget" & firstname == "Karl" & gdenr_2018 == 1372
replace ID = "SZ-1975-9014" if ID == "" & year == 1984 & name == "Weber Wiget" & firstname == "Karl" & gdenr_2018 == 1372
replace ID = "SZ-1975-9014" if ID == "" & year == 1985 & name == "Weber Wiget" & firstname == "Karl" & gdenr_2018 == 1372
replace ID = "SG-1995-0177" if ID == "" & year == 1997 & name == "Weigelt" & firstname == "Peter" & gdenr_2018 == 3214
replace ID = "SG-1995-0177" if ID == "" & year == 1998 & name == "Weigelt" & firstname == "Peter" & gdenr_2018 == 3214
replace ID = "SG-1995-0177" if ID == "" & year == 1999 & name == "Weigelt" & firstname == "Peter" & gdenr_2018 == 3214
replace ID = "BS-1983-0069" if ID == "" & year == 1983 & name == "Wick" & firstname == "Hugo" & gdenr_2018 == 2701
replace ID = "BS-1983-0069" if ID == "" & year == 1984 & name == "Wick" & firstname == "Hugo" & gdenr_2018 == 2701
replace ID = "BS-1983-0069" if ID == "" & year == 1985 & name == "Wick" & firstname == "Hugo" & gdenr_2018 == 2701
replace ID = "BS-1983-0069" if ID == "" & year == 1986 & name == "Wick" & firstname == "Hugo" & gdenr_2018 == 2701
replace ID = "BS-1983-0069" if ID == "" & year == 1991 & name == "Wick" & firstname == "Hugo" & gdenr_2018 == 2701
replace ID = "SG-1987-0089" if ID == "" & year == 1987 & name == "Widrig" & firstname == "Hans Werner" & gdenr_2018 == 3291
replace ID = "ZH-1987-0773" if ID == "" & year == 1987 & name == "Wiederkehr" & firstname == "Roland" & gdenr_2018 == 241
replace ID = "BS-1987-0081" if ID == "" & year == 2002 & name == "Wirz" & firstname == "Christine" & gdenr_2018 == 2701
replace ID = "BS-1987-0081" if ID == "" & year == 2003 & name == "Wirz" & firstname == "Christine" & gdenr_2018 == 2701
replace ID = "SG-1987-0091" if ID == "" & year == 1991 & name == "Wittenwiler" & firstname == "Milli" & gdenr_2018 == 3379
replace ID = "SG-1987-0091" if ID == "" & year == 1992 & name == "Wittenwiler" & firstname == "Milli" & gdenr_2018 == 3379
replace ID = "SG-1987-0091" if ID == "" & year == 1993 & name == "Wittenwiler" & firstname == "Milli" & gdenr_2018 == 3379
replace ID = "SG-1987-0091" if ID == "" & year == 1994 & name == "Wittenwiler" & firstname == "Milli" & gdenr_2018 == 3379
replace ID = "BEJU-1991-0543" if ID == "" & year == 2000 & name == "Wyss" & firstname == "Ursula" & gdenr_2018 == 99999
replace ID = "BEJU-1991-0543" if ID == "" & year == 2001 & name == "Wyss" & firstname == "Ursula" & gdenr_2018 == 99999
replace ID = "BEJU-1991-0543" if ID == "" & year == 2002 & name == "Wyss" & firstname == "Ursula" & gdenr_2018 == 99999
replace ID = "BEJU-1991-0543" if ID == "" & year == 2003 & name == "Wyss" & firstname == "Ursula" & gdenr_2018 == 99999
replace ID = "BEJU-1975-0489" if ID == "" & year == 1999 & name == "Wyss" & firstname == "William A." & gdenr_2018 == 988
replace ID = "BEJU-1975-0489" if ID == "" & year == 1989 & name == "Wyss" & firstname == "William Andreas" & gdenr_2018 == 988
replace ID = "BEJU-1975-0489" if ID == "" & year == 1990 & name == "Wyss" & firstname == "William Andreas" & gdenr_2018 == 988
replace ID = "BEJU-1975-0489" if ID == "" & year == 1993 & name == "Wyss" & firstname == "William Andreas" & gdenr_2018 == 988
replace ID = "BEJU-1975-0489" if ID == "" & year == 1994 & name == "Wyss" & firstname == "William Andreas" & gdenr_2018 == 988
replace ID = "BEJU-1975-0489" if ID == "" & year == 1996 & name == "Wyss" & firstname == "William Andreas" & gdenr_2018 == 988
replace ID = "BEJU-1975-0489" if ID == "" & year == 1997 & name == "Wyss" & firstname == "William Andreas" & gdenr_2018 == 988
replace ID = "BEJU-1975-0489" if ID == "" & year == 1998 & name == "Wyss" & firstname == "William Andreas" & gdenr_2018 == 988
replace ID = "SO-1999-0084" if ID == "" & year == 2000 & name == "Zanetti" & firstname == "Roberto" & gdenr_2018 == 99999
replace ID = "SO-1999-0084" if ID == "" & year == 2001 & name == "Zanetti" & firstname == "Roberto" & gdenr_2018 == 99999
replace ID = "ZH-1987-0794" if ID == "" & year == 2000 & name == "Zapfl Helbling" & firstname == "Rosmarie" & gdenr_2018 == 191
replace ID = "ZH-1987-0794" if ID == "" & year == 2001 & name == "Zapfl Helbling" & firstname == "Rosmarie" & gdenr_2018 == 191
replace ID = "ZH-1987-0794" if ID == "" & year == 2002 & name == "Zapfl Helbling" & firstname == "Rosmarie" & gdenr_2018 == 191
replace ID = "ZH-1987-0794" if ID == "" & year == 2003 & name == "Zapfl Helbling" & firstname == "Rosmarie" & gdenr_2018 == 191
replace ID = "AG-1979-0122" if ID == "" & year == 1999 & name == "Zbinden" & firstname == "Hans" & gdenr_2018 == 4021
replace ID = "AG-1979-0122" if ID == "" & year == 2000 & name == "Zbinden" & firstname == "Hans" & gdenr_2018 == 99999
replace ID = "AG-1979-0122" if ID == "" & year == 2001 & name == "Zbinden" & firstname == "Hans" & gdenr_2018 == 99999
replace ID = "FR-1975-0027" if ID == "" & year == 1982 & name == "Zbinden" & firstname == "Paul" & gdenr_2018 == 2196
replace ID = "FR-1975-0027" if ID == "" & year == 1983 & name == "Zbinden" & firstname == "Paul" & gdenr_2018 == 2196
replace ID = "FR-1975-0027" if ID == "" & year == 1984 & name == "Zbinden" & firstname == "Paul" & gdenr_2018 == 2196
replace ID = "FR-1975-0027" if ID == "" & year == 1985 & name == "Zbinden" & firstname == "Paul" & gdenr_2018 == 2196
replace ID = "FR-1975-0027" if ID == "" & year == 1989 & name == "Zbinden" & firstname == "Paul" & gdenr_2018 == 2196
replace ID = "FR-1975-0027" if ID == "" & year == 1990 & name == "Zbinden" & firstname == "Paul" & gdenr_2018 == 2196
replace ID = "GE-1967-9094" if ID == "" & year == 1996 & name == "Ziegler" & firstname == "Jean" & gdenr_2018 == 6621
replace ID = "GE-1967-9094" if ID == "" & year == 1997 & name == "Ziegler" & firstname == "Jean" & gdenr_2018 == 6621
replace ID = "GE-1967-9094" if ID == "" & year == 1998 & name == "Ziegler" & firstname == "Jean" & gdenr_2018 == 6621
replace ID = "GE-1967-9094" if ID == "" & year == 1999 & name == "Ziegler" & firstname == "Jean" & gdenr_2018 == 6621
replace ID = "GE-1967-9094" if ID == "" & year == 1987 & name == "Ziegler" & firstname == "Josef" & gdenr_2018 == 6614
replace ID = "GE-1967-9094" if ID == "" & year == 1988 & name == "Ziegler" & firstname == "Josef" & gdenr_2018 == 6614
replace ID = "VD-1987-0167" if ID == "" & year == 1996 & name == "Zisyadis" & firstname == "Josef" & gdenr_2018 == 5586
replace ID = "ZH-1987-0805" if ID == "" & year == 2000 & name == "Zuppiger" & firstname == "Bruno" & gdenr_2018 == 117
replace ID = "ZH-1987-0805" if ID == "" & year == 2001 & name == "Zuppiger" & firstname == "Bruno" & gdenr_2018 == 117
replace ID = "ZH-1987-0805" if ID == "" & year == 2002 & name == "Zuppiger" & firstname == "Bruno" & gdenr_2018 == 117
replace ID = "ZH-1987-0805" if ID == "" & year == 2003 & name == "Zuppiger" & firstname == "Bruno" & gdenr_2018 == 117
replace ID = "BEJU-1979-0380" if ID == "" & year == 1991 & name == "Zwahlen" & firstname == "Jean Claude" & gdenr_2018 == 717
replace ID = "SG-1979-0054" if ID == "" & year == 1987 & name == "Zwingli" & firstname == "Walter" & gdenr_2018 == 3235
replace ID = "BEJU-1983-0431" if ID == "" & year == 1983 & name == "Zwygart" & firstname == "Otto" & gdenr_2018 == 352
replace ID = "BEJU-1983-0431" if ID == "" & year == 1996 & name == "Zwygart " & firstname == "Otto" & gdenr_2018 == 352
replace ID = "BEJU-1983-0431" if ID == "" & year == 1997 & name == "Zwygart " & firstname == "Otto" & gdenr_2018 == 352
replace ID = "BEJU-1983-0431" if ID == "" & year == 1998 & name == "Zwygart " & firstname == "Otto" & gdenr_2018 == 352
replace ID = "BEJU-1983-0431" if ID == "" & year == 1999 & name == "Zwygart " & firstname == "Otto" & gdenr_2018 == 352
replace ID = "BEJU-1983-0431" if ID == "" & year == 2000 & name == "Zwygart " & firstname == "Otto" & gdenr_2018 == 352
replace ID = "BEJU-1987-0529" if ID == "" & year == 1987 & name == "Zölch Balmer" & firstname == "Elisabeth" & gdenr_2018 == 351
replace ID = "SZ-1983-0009" if ID == "" & year == 1992 & name == "Züger" & firstname == "Arthur" & gdenr_2018 == 1349
replace ID = "SZ-1983-0009" if ID == "" & year == 1993 & name == "Züger" & firstname == "Arthur" & gdenr_2018 == 1349
replace ID = "SZ-1983-0009" if ID == "" & year == 1994 & name == "Züger" & firstname == "Arthur" & gdenr_2018 == 1349
replace ID = "GE-1991-0115" if ID == "" & year == 1991 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6621
replace ID = "GE-1991-0115" if ID == "" & year == 1992 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6621
replace ID = "GE-1991-0115" if ID == "" & year == 1993 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6621
replace ID = "GE-1991-0115" if ID == "" & year == 1994 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6621
replace ID = "GE-1991-0115" if ID == "" & year == 1996 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6621
replace ID = "GE-1991-0115" if ID == "" & year == 1997 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6621
replace ID = "GE-1991-0115" if ID == "" & year == 1998 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6612
replace ID = "GE-1991-0115" if ID == "" & year == 1999 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6612
replace ID = "GE-1991-0115" if ID == "" & year == 2000 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6612
replace ID = "GE-1991-0115" if ID == "" & year == 2001 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6612
replace ID = "GE-1991-0115" if ID == "" & year == 2002 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6612
replace ID = "GE-1991-0115" if ID == "" & year == 2003 & name == "de Darel" & firstname == "Jean Nils" & gdenr_2018 == 6612
replace ID = "BEJU-1995-0584" if ID == "" & year == 1996 & name == "von Allmen" & firstname == "Hans Ueli" & gdenr_2018 == 942
replace ID = "BEJU-1995-0584" if ID == "" & year == 1997 & name == "von Allmen" & firstname == "Hans Ueli" & gdenr_2018 == 947
replace ID = "BEJU-1995-0584" if ID == "" & year == 1998 & name == "von Allmen" & firstname == "Hans Ueli" & gdenr_2018 == 947
replace ID = "BEJU-1995-0584" if ID == "" & year == 1999 & name == "von Allmen" & firstname == "Hans Ueli" & gdenr_2018 == 947
replace ID = "BS-1991-0078" if ID == "" & year == 1991 & name == "von Felten" & firstname == "Margrith" & gdenr_2018 == 2701
drop if ID == ""
sort ID year
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\non-match_mandates_corrected.dta", replace

// merge with NR data
use "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo_tmp.dta", clear
sort ID year

merge 1:m ID year using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\non-match_mandates_corrected.dta"
drop if _merge == 1
drop _merge

append using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NR_mandate_match.dta"

keep ID year name firstname gdename gdenr_2018 canton Companyname Companylocation company_norm
duplicates tag ID year gdenr_2018 company_norm, gen(dupl)  // no duplicates
drop dupl 

rename company_norm firmname
rename ID NRid

sort NRid year firmname

foreach var in canton name firstname gdename gdenr_2018 Companyname Companylocation {
	rename `var' `var'_mand
}
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NR_ID_GT2_tmp.dta", replace

** Export unique NR-ID: who is in the sample?
duplicates drop NRid, force
keep NRid
sort NRid
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\uniqueNR_ID_GT2.dta", replace


*** Merge NR Mandates with RL-Sugarcube output after (FP/FN) post-processing

** Prepare 297 observations
/*
use "$path\02_Processed_data\10_Directors_1934_2003\RL_sugarcube_297.dta", clear
merge 1:m PID year using "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta"
keep if _merge == 3
drop _merge
duplicates tag PID NRid year CID, gen(dupl)
br if dupl>0
drop if PID == 15826 & year == 1965 & ctn != "sg"
drop if PID == 16629 & year == 1966 & ctn != "sg"
drop if PID == 103870 & year == 1990 & ctn != "sg"
duplicates report PID NRid year CID
keep PID NRid year CID 
gen east_miss = 1 
save "$path\02_Processed_data\10_Directors_1934_2003\RL_sugarcube_297_CID.dta", replace 
*/

** Sugarcube firm-links after post-processing

* merge on NRid year firmname

use "$path\02_Processed_data\10_Directors_1934_2003\RL_NR-Sugarcube_1934-2003.dta", clear 

merge m:1 PID NRid year CID using "$path\02_Processed_data\10_Directors_1934_2003\RL_sugarcube_297_CID.dta"
drop _merge

drop firmname // this is "old" normalization of firmnames in firmname-cleaning process.
gen firmname = ""
replace firmname = ustrtrim(cname_orig)
replace firmname = ustrlower(firmname)
replace firmname = usubinstr(firmname, "ü", "ue",.)
replace firmname = usubinstr(firmname, "ä", "ae",.)
replace firmname = usubinstr(firmname, "ö", "oe",.)
replace firmname = usubinstr(firmname, "é", "e",.)
replace firmname = usubinstr(firmname, "è", "e",.)
replace firmname = usubinstr(firmname, "ê", "e",.)
replace firmname = usubinstr(firmname, "ç", "c",.)
replace firmname = usubinstr(firmname, `"""', "",.)
replace firmname = usubinstr(firmname, "«", "",.)
replace firmname = usubinstr(firmname, "»", "",.)
replace firmname = usubinstr(firmname, ".", "",.)
replace firmname = usubinstr(firmname, ",", "",.)
replace firmname = usubinstr(firmname, ";", "",.)
replace firmname = usubinstr(firmname, "-", "",.)
replace firmname = usubinstr(firmname, ":", "",.)
replace firmname = usubinstr(firmname, "'", "",.)
replace firmname = usubinstr(firmname, "  ", "",.)

sort NRid year firmname
recast str1000 firmname

// restrict sample to NRs that are covered in NR-VR Verzeichnis
merge m:1 NRid using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\uniqueNR_ID_GT2.dta"
keep if _merge==3
drop _merge
drop if year < 1982
drop if year > 2003

replace firmname = "none" if firmname == ""

duplicates report NRid year PID CID  // no duplicates. BUT: we want to prepare the comparison of NRid's with firmnames from GT2 -- > look for duplicate firmnames (w/o PID/CID)
duplicates report NRid year firmname
duplicates tag NRid year firmname, gen(dupl)
duplicates drop NRid year firmname, force // looks as if these are genuine duplicates (also w.r.t. location, capital, residential municipality)
drop dupl 

keep if east_miss == 1  // just focus on the originally 297 missing observations

merge m:1 NRid year firmname using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NR_ID_GT2_tmp.dta"

// prepare for manual coding: GT2 subset
preserve
keep if _merge == 2
gen source = "GT2"
keep NRid canton year name_mand firstname_mand gdename_mand firmname Companyname_mand Companylocation_mand source
rename name_mand name
rename firstname_mand firstname 
rename gdename_mand gdename 
rename Companyname_mand firm_orig 
rename Companylocation_mand firmlocation
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_eastmiss.dta", replace
restore

// prepare for manual coding: RL subset
keep if _merge == 1
gen source = "RL"
keep NRid canton year name firstname gdename cname_orig firmname firmlocation firstname_bus name_bus gdename_bus source east_miss
rename name name_nr
rename firstname firstname_nr
rename gdename gdename_nr
rename name_bus name
rename firstname_bus firstname 
rename gdename_bus gdename
rename cname_orig firm_orig
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge1_eastmiss.dta", replace

// combine non-matched observations for manual evaluation & linkage
append using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_eastmiss.dta"

replace firm_orig = "none" if firmname == "none"

// Identify affected NRids
preserve
keep if east_miss == 1
duplicates drop NRid year, force
keep NRid year
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NRid_297.dta", replace
restore

merge m:1 NRid year using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NRid_297.dta"
keep if _merge == 3
drop _merge


// create alternating ID-Block variable to distinguish different blocks of NR-IDs (even/odd)
egen help1 = group(year NRid)
gen NROddBlocks = 0 
replace NROddBlocks = 1 if mod(help1, 2) == 1  // odd IDs
drop help1

// create variable to link the sources (create this variable for GT2)
bysort year NRid source (firmname): gen linkIDs = _n
replace linkIDs = . if source == "RL"

gen out = .

keep NROddBlocks source year NRid out linkIDs firm_orig firmlocation firstname name gdename
order NROddBlocks source year NRid out linkIDs firm_orig firmlocation firstname name gdename

sort year NRid firm_orig source

/*/ export information for coding by individuals  -- > deactivated as manual coding is completed
export excel using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_coding297.xlsx", first(var) replace
*/

* import manually coded RL-GT2 data from the FIRST ROUND!
import excel "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_coding_in.xlsx", firstrow clear

*identify observations affected by 297-error
merge m:1 NRid year using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NRid_297.dta"


erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NR_mandate_match.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\non-match_mandates.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\non-match_mandates_corrected.dta"
erase "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo_tmp.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\Mandate_Nationalraete_geo.dta"

erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge1_eastmiss.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_eastmiss.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\uniqueNR_ID_GT2.dta"

erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NR_ID_GT2_tmp.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NRid_297.dta"

