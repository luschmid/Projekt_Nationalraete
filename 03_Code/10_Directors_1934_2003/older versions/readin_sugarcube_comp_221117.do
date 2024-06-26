version 17
clear all
set more off
cap log close

*global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\"
global path "C:\Users\yanni\Dropbox\Projekt Nationalräte\"

******************************************
*** Import and Geo-Code Sugarcube Data ***
******************************************

** import from original csv files (saved as ANSI (1934-1963) else UTF-8 with BOM)
foreach val in 1934 1943 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_companies.csv", clear  // ANSI code, delimiter is "tab"
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\comp_`val'_tmp.dta", replace
}
foreach val in 1960 1962 1963 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_companies.csv", delimiters(";") clear // ANSI code, delimiter is ";"
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\comp_`val'_tmp.dta", replace
}
foreach val in 1964 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_companies.csv", delimiters(";") encoding("utf-8") clear
	replace year = year  
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\comp_`val'_tmp.dta", replace
}
foreach val in 1966 1967 1970 1973 1976 1979_1980 1980_1981 1998_1999 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_companies.csv", delimiters(";") encoding("utf-8") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\comp_`val'_tmp.dta", replace
}
forv val = 1982(1)1998 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_companies.csv", delimiters(";") encoding("utf-8") clear
	replace year = year - 1 
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\comp_`val'_tmp.dta", replace
}
forv val = 1999(1)2003 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_companies.csv", delimiters(";") encoding("utf-8") clear
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\comp_`val'_tmp.dta", replace
}
use "$path\02_Processed_data\10_Directors_1934_2003\comp_1934_tmp.dta", clear
foreach val in 1943 1960 1962 1963 1964 1966 1967 1970 1973 1976 1979_1980 1980_1981 1998_1999 {
	append using "$path\02_Processed_data\10_Directors_1934_2003\comp_`val'_tmp.dta", force
	erase "$path\02_Processed_data\10_Directors_1934_2003\comp_`val'_tmp.dta"
}
erase "$path\02_Processed_data\10_Directors_1934_2003\comp_1934_tmp.dta"
forv val = 1982(1)2003 {
	append using "$path\02_Processed_data\10_Directors_1934_2003\comp_`val'_tmp.dta", force
	erase "$path\02_Processed_data\10_Directors_1934_2003\comp_`val'_tmp.dta"
}
save "$path\02_Processed_data\10_Directors_1934_2003\comp_1934-2003_tmp1.dta", replace


** Cleaning of special characters, etc. **
use "$path\02_Processed_data\10_Directors_1934_2003\comp_1934-2003_tmp1.dta", clear
rename id CID
rename name cname
rename city gdename
rename address caddress
label var CID "Company ID per year"
destring CID, replace
destring year, replace

foreach str in cname caddress gdename {
	replace `str' = "" if `str' == "."
	gen `str'_orig = `str'
	replace `str' = ustrlower(`str')
	replace `str' = ustrtrim(`str')
}
save "$path\02_Processed_data\10_Directors_1934_2003\comp_1934-2003_tmp2.dta", replace


**** Geo-Coding 
*ssc install mdesc
mdesc cname caddress gdename owners capital function signature 
drop caddress // it is 100% missing

** Cleaning of special characters, etc. **
foreach str in gdename /*cname*/ {
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
	replace `str' = usubinstr(`str', "-", " ",.)
	replace `str' = usubinstr(`str', "~", " ",.)	
	replace `str' = usubinstr(`str', ":", " ",.)
	replace `str' = usubinstr(`str', "    ", " ",.)
	replace `str' = usubinstr(`str', "   ", " ",.)
	replace `str' = usubinstr(`str', "  ", " ",.)
	replace `str' = ustrtrim(`str')
}
recast str144 gdename  // instead of strL, which cannot be used in merger command
save "$path\02_Processed_data\10_Directors_1934_2003\comp_1934-2003_tmp3.dta", replace

*** String-merge with sugarcube data **

** Prepare Geo-Datasets 
* Standardize municipality names in geo dataset **
use "$path\02_Processed_data\08_Municipalities\GdeNr2018_Geo_1931-2018.dta", clear
keep year ctn gdename gdenr gdenr_2018 E_CNTR N_CNTR
drop if year==. & gdename==""
duplicates report gdename year
duplicates drop gdename year, force // drop 4 observations
replace gdename = usubinstr(gdename, "ä", "ae",.)
replace gdename = usubinstr(gdename, "ö", "oe",.)
replace gdename = usubinstr(gdename, "ü", "ue",.)
replace gdename = usubinstr(gdename, "È", "e",.)
replace gdename = usubinstr(gdename, "è", "e",.)
replace gdename = usubinstr(gdename, "é", "e",.)
replace gdename = usubinstr(gdename, "ê", "e",.)
replace gdename = usubinstr(gdename, "â", "a",.)
replace gdename = usubinstr(gdename, "ô", "o",.)
replace gdename = usubinstr(gdename, "î", "i",.)
replace gdename = usubinstr(gdename, "'", " ",.)
replace gdename = usubinstr(gdename, "`", " ",.)
replace gdename = usubinstr(gdename, "´", " ",.)
replace gdename = usubinstr(gdename, ",", " ",.)
replace gdename = usubinstr(gdename, "-", " ",.)
replace gdename = usubinstr(gdename, ":", " ",.)
replace gdename = ustrlower(gdename)
replace ctn = ustrlower(ctn)
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

*********************
** Geo-Merge steps **
*********************

/*
Description of merge strategy
1) Merge our company dataset with normalized company-gdename information
2) a) String merge on further standardized company-gdename info (step 1: Mark's code)
2) b) String merge on further standardized company-gdename info (step 2: 300 biggest cases)
2) c) String merge on further standardized company-gdename info (step 3)
3) Extraction of information from gdename var (succursalle, branch, etc.)
4) Extraction of information from capital var (to find where the company is located)
5) Extraction of information from cname var (cases where municipality ended up in company name)
*/

* 1) String merge on normalized company-gdename info
use "$path\02_Processed_data\10_Directors_1934_2003\comp_1934-2003_tmp3.dta", clear
sort gdename year
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2
gen geo_merge = 1 if _merge == 3

preserve
keep if _merge == 3
save "$path\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_1.dta", replace
restore

* 2) a) String merge on further standardized company-gdename info (step 1: Mark's code)

keep if _merge == 1
drop _merge

replace gdename = "geneve" if gdename == "genve"
replace gdename = "sierre" if gdename == "sierre/siders"
replace gdename = "zuerich" if gdename == "zrich"
replace gdename = "biel (be)" if gdename == "biel/bienne"
replace gdename = "ollon" if gdename == "ollon vd"
replace gdename = "murten" if gdename == "murten morat"
replace gdename = "duedingen" if gdename == "duedingen (guin)"
replace gdename = "brugg" if gdename == "brugg ag"
replace gdename = "st margrethen" if gdename == "margrethen sg"
replace gdename = "lamone" if gdename == "lamone cadempino"
replace gdename = "villars sur glane" if gdename == "villars sur glaene"
replace gdename = "gams" if gdename == "garns"
replace gdename = "wisen (so)" if gdename == "wiesen"
replace gdename = "biel (be)" if gdename == "bienne"
replace gdename = "le locle" if gdename == "le lode"
replace gdename = "muri bei bern" if gdename == "muri/bern"
replace gdename = "davos" if gdename == "davos platz"
replace gdename = "sainte croix" if gdename == "ste croix"
replace gdename = "saint imier" if gdename == "st imier"
replace gdename = "le chenit" if gdename == "le brassus"
replace gdename = "le chenit" if gdename == "le sentier"
replace gdename = "berneck" if gdename == "heerbrugg"
replace gdename = "vaz/obervaz" if gdename == "lenzerheide"
replace gdename = "pregny chambesy" if gdename == "pregny"
replace gdename = "kuessnacht (sz)" if gdename == "kuessnacht"
replace gdename = "kuesnacht (zh)" if gdename == "kiisnacht"
replace gdename = "evilard" if gdename == "leubringen"
replace gdename = "saint maurice" if gdename == "st maurice"
replace gdename = "lancy" if gdename == "petit lancy"
replace gdename = "rueschlikon" if gdename == "ruschlikon"
replace gdename = "saint imier" if gdename == "st. imier"
replace gdename = "kirchberg (sg)" if gdename == "bazenheid"
replace gdename = "muenchenstein" if gdename == "munchenstein"
replace gdename = "teufen" if gdename == "niederteufen"
replace gdename = "collonge bellerive" if gdename == "collonges bellerive"
replace gdename = "neuhausen am rheinfall" if gdename == "neuhausen"
replace gdename = "zollikon" if gdename == "zollikerberg"
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
replace gdename = "disentis/muster" if gdename == "disentis"
replace gdename = "muensingen" if gdename == "munsingen"
replace gdename = "littau" if gdename == "reussbuehl"
replace gdename = "chatel saint denis" if gdename == "chatel st denis"
replace gdename = "rueschlikon" if gdename == "riischlikon"
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
replace gdename = "granges pres marnand" if gdename == "granges marnand"
replace gdename = "murten" if gdename == "morat"
replace gdename = "ittigen" if gdename == "worblaufen"
replace gdename = "feldbrunnen st. niklaus" if gdename == "feldbrunnen"
replace gdename = "langnau im emmental" if gdename == "langnau i. e."
replace gdename = "langnau im emmental" if gdename == "langnau i. e"
replace gdename = "feusisberg" if gdename == "schindellegi"
replace gdename = "ruederswil" if gdename == "zollbrueck"
replace gdename = "delemont" if gdename == "delsberg"
replace gdename = "reichenbach im kandertal" if gdename == "reichenbach"
replace gdename = "altstaetten" if gdename == "alstaetten"
replace gdename = "muri bei bern" if gdename == "guemligen muri"
replace gdename = "marin epagnier" if gdename == "marin"
replace gdename = "vilters wangs" if gdename == "wangs"
replace gdename = "crans montana" if gdename == "crans/chermignon"
replace gdename = "diessbach bei bueren" if gdename == "diessbach"
replace gdename = "l abbaye" if gdename == "le pont"
replace gdename = "saint legier la chiesaz" if gdename == "st legier"
replace gdename = "bern" if gdename == "bern buempliz"
replace gdename = "bussigny" if gdename == "bussigny s/morges"
replace gdename = "chezard saint martin" if gdename == "chezard"
replace gdename = "nuglar st. pantaleon" if gdename == "nuglar"
replace gdename = "le chenit" if gdename == "orient"
replace gdename = "villars sur glane" if gdename == "villars/glane"
replace gdename = "ollon" if gdename == "villars/ollon"
replace gdename = "muri bei bern" if gdename == "gumligen"
replace gdename = "walenstadt" if gdename == "wallenstadt"
replace gdename = "crans montana" if gdename == "crans/sierre"
replace gdename = "kuesnacht (zh)" if gdename == "goldbach"
replace gdename = "montreux" if gdename == "montreux clarens"
replace gdename = "muri bei bern" if gdename == "muri bern"
replace gdename = "geneve" if gdename == "genf"
replace gdename = "schwyz" if gdename == "ibach"
replace gdename = "ittigen" if gdename == "papiermuehle"
replace gdename = "saint ursanne" if gdename == "st ursanne"
replace gdename = "belmont sur lausanne" if gdename == "belmont/lausanne"
replace gdename = "lugano" if gdename == "cassarate"
replace gdename = "davos" if gdename == "davos glaris"
replace gdename = "la tour de peilz" if gdename == "la tour de peilt"
replace gdename = "le mont sur lausanne" if gdename == "mont s/lausanne"
replace gdename = "crans montana" if gdename == "montana/vermala"
replace gdename = "muenchwilen (tg)" if gdename == "st. margarethen"
replace gdename = "saint legier la chiesaz" if gdename == "st legier la chiesaz"
replace gdename = "cossonay" if gdename == "cossonay gare"
replace gdename = "gruyeres" if gdename == "epagny"
replace gdename = "hospental" if gdename == "hospenthal"
replace gdename = "kuesnacht (zh)" if gdename == "kusnacht"
replace gdename = "koeniz" if gdename == "niederwangen"
replace gdename = "saint ursanne" if gdename == "st. ursanne"
replace gdename = "corsier sur vevey" if gdename == "corsier/vevey"
replace gdename = "kuesnacht (zh)" if gdename == "kosnacht"
replace gdename = "neukirch an der thur" if gdename == "neukirch"
replace gdename = "arth" if gdename == "oberarth"
replace gdename = "zuerich" if gdename == "ziirich"
replace gdename = "nendaz" if gdename == "aproz"
replace gdename = "chene bougeries" if gdename == "chenes bougeries"
replace gdename = "montreux" if gdename == "chernex"
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
replace gdename = "lancy" if gdename == "petit laney"
replace gdename = "rueegsau" if gdename == "rueegsauschachen"
replace gdename = "rueegsau" if gdename == "ruegsauschachen"
replace gdename = "sonceboz sombeval" if gdename == "sonceboz"
replace gdename = "boudry" if gdename == "areuse"
replace gdename = "boettstein" if gdename == "bottstein"
replace gdename = "duebenedorf" if gdename == "dubendorf"
replace gdename = "wigoltingen" if gdename == "hasli wigoltingen"
replace gdename = "kilchberg (zh)" if gdename == "kilchberg/zuerich"
replace gdename = "niederhelfenschwil" if gdename == "lenggenwil"
replace gdename = "montilliez" if gdename == "montilier"
replace gdename = "vaz/obervaz" if gdename == "vaz"
replace gdename = "koeniz" if gdename == "thoerishaus"
replace gdename = "pregny chambesy" if gdename == "chambesy pregny"
replace gdename = "montagny (fr)" if gdename == "cousset"
replace gdename = "crans montana" if gdename == "crans s/sierre"
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
replace gdename = "quarten" if gdename == "unterterzen"
replace gdename = "zollikon" if gdename == "zolliken"
replace gdename = "neuenegg" if gdename == "brueggelbach"
replace gdename = "bern" if gdename == "buempliz"
replace gdename = "bussigny" if gdename == "bussigny s/ morges"
replace gdename = "muri bei bern" if gdename == "giimligen"
replace gdename = "freienbach" if gdename == "hurden"
replace gdename = "kuessnacht (zh)" if gdename == "kilsnacht"
replace gdename = "kuessnacht (sz)" if gdename == "kussnacht"
replace gdename = "crans montana" if gdename == "montana station"
replace gdename = "ruemlang" if gdename == "rumlang"
replace gdename = "simplon" if gdename == "simplon dorf"
replace gdename = "sonceboz sombeval" if gdename == "sombeval"
replace gdename = "egnach" if gdename == "steinebrunn"
replace gdename = "saint imier" if gdename == "st lmier"
replace gdename = "taegerwilen" if gdename == "tagerwilen"
replace gdename = "glarus nord" if gdename == "ziegelbruecke"
replace gdename = "carouge (ge)" if gdename == "carouge geneve"
replace gdename = "ollon" if gdename == "chesieres"
replace gdename = "collombey muraz" if gdename == "collombey"
replace gdename = "crans pres celigny" if gdename == "crans/celigny"
replace gdename = "crans montana" if gdename == "crans chermignon"
replace gdename = "kuessnacht (sz)" if gdename == "immensee"
replace gdename = "lausanne" if gdename == "lau sanne"
replace gdename = "faido" if gdename == "lavorgo"
replace gdename = "collombey muraz" if gdename == "muraz"
replace gdename = "sumvitg" if gdename == "rabius"
replace gdename = "worb" if gdename == "ruefenacht worb"
replace gdename = "waltensburg/vuorz" if gdename == "waltensburg"
replace gdename = "adligenswil" if gdename == "adligenschwil"
replace gdename = "riedholz" if gdename == "attisholz"
replace gdename = "port valais" if gdename == "bouveret"
replace gdename = "buelach" if gdename == "bulach"
replace gdename = "bussigny" if gdename == "bussigny sur morges"
replace gdename = "schuebelbach" if gdename == "buttikon"
replace gdename = "ormont dessus" if gdename == "diablerets"
replace gdename = "sumiswald" if gdename == "gruenen/sumiswald"
replace gdename = "muri bei bern" if gdename == "gumlingen"
replace gdename = "koeniz" if gdename == "liebefeld/koeniz"
replace gdename = "luetzelflueh" if gdename == "lutzelflueh"
replace gdename = "neckertal" if gdename == "nassen"
replace gdename = "neuchatel" if gdename == "neuenburg"
replace gdename = "schoenenwerd" if gdename == "schonenwerd"
replace gdename = "la grande beroche" if gdename == "st aubin sauges"
replace gdename = "saint george" if gdename == "st george"
replace gdename = "thonex" if gdename == "thoenex"
replace gdename = "trun" if gdename == "truns"
replace gdename = "schwende" if gdename == "weissbad"
replace gdename = "wartau" if gdename == "weite"
replace gdename = "brig" if gdename == "brigue"
replace gdename = "chavannes pres renens" if gdename == "chavannes renens"
replace gdename = "chevenez" if gdename == "chevenaz"
replace gdename = "egg" if gdename == "esslingen"
replace gdename = "kuesnacht (zh)" if gdename == "kuesnacht"
replace gdename = "hermetschwil staffeln" if gdename == "bremgarten"
replace gdename = "wetzikon (zh)" if gdename == "kempten wetzikon"
replace gdename = "boettstein" if gdename == "kleindottingen"
replace gdename = "le chenit" if gdename == "l orient"
replace gdename = "evilard" if gdename == "macolin"
replace gdename = "montreux" if gdename == "montreux les planches"
replace gdename = "moerschwil" if gdename == "morschwil"
replace gdename = "lauterbrunnen" if gdename == "muerren"
replace gdename = "pratteln" if gdename == "prattein"
replace gdename = "schaffhausen" if gdename == "schauffhausen"
replace gdename = "montreux" if gdename == "territet montreux"
replace gdename = "zug" if gdename == "zug oberwil"
replace gdename = "berguen/bravuogn" if gdename == "berguen"
replace gdename = "bern" if gdename == "bern bumpliz"
replace gdename = "chezard saint martin" if gdename == "chezard st martin"
replace gdename = "corsier sur vevey" if gdename == "corsier s/vevey"
replace gdename = "davesco soragno" if gdename == "davesco"
replace gdename = "les geneveys sur coffrane" if gdename == "geneveys sur coffrane"
replace gdename = "emmen" if gdename == "gerliswil"
replace gdename = "opfikon" if gdename == "glattbrugg opfikon"
replace gdename = "reichenbach im kandertal" if gdename == "kandertal"
replace gdename = "feuerthalen" if gdename == "langwiesen"
replace gdename = "montreux chatelard" if gdename == "le chatelard montreux"
replace gdename = "ormont dessus" if gdename == "les diablerets"
replace gdename = "leuk" if gdename == "leuk susten"
replace gdename = "bregaglia" if gdename == "maloja"
replace gdename = "montagny pres yverdon" if gdename == "montagny/yverdon"
replace gdename = "nesslau" if gdename == "neu st johann"
replace gdename = "pfaeffikon" if gdename == "pfaffikon"
replace gdename = "bad ragaz" if gdename == "ragaz"
replace gdename = "zell (zh)" if gdename == "rikon"
replace gdename = "schuepfheim" if gdename == "schupfheim"
replace gdename = "sils im engadin/segl" if gdename == "sils maria"
replace gdename = "sion" if gdename == "sitten"
replace gdename = "taeuffelen" if gdename == "tauffelen"
replace gdename = "zofingen" if gdename == "zofmgen"
replace gdename = "bre aldesago" if gdename == "aldesago di bre"
replace gdename = "alpnach" if gdename == "alpnachstad"
replace gdename = "quinto" if gdename == "ambri di quinto"
replace gdename = "wartau" if gdename == "azmoos"
replace gdename = "beatenberg" if gdename == "beatenbucht"
replace gdename = "deisswil bei muenchenbuchsee" if gdename == "deisswil"
replace gdename = "marthalen" if gdename == "ellikon"
replace gdename = "emmen" if gdename == "emmenbruecke emmen"
replace gdename = "ennetbuergen" if gdename == "ennetburgen"
replace gdename = "flims" if gdename == "films"
replace gdename = "geneve" if gdename == "ganeve"
replace gdename = "reichenbach im kandertal" if gdename == "kiental"
replace gdename = "boettstein" if gdename == "kleindoettingen"
replace gdename = "le mont sur lausanne" if gdename == "le mont s/lausanne"
replace gdename = "lengnau (be)" if gdename == "lengnau bei biel"
replace gdename = "vaz/obervaz" if gdename == "lenzerheid"
replace gdename = "leuk" if gdename == "loeche"
replace gdename = "rochefort" if gdename == "montezillon rochefort"
replace gdename = "neuhausen am rheinfall" if gdename == "neuhausen a. r."
replace gdename = "otringen" if gdename == "of tringen"
replace gdename = "olten" if gdename == "often"
replace gdename = "buchrain" if gdename == "perlen"
replace gdename = "urtenen schoenbuehl" if gdename == "schoenbuehl"
replace gdename = "uetikon am see" if gdename == "uetikon a see"
replace gdename = "zumikon" if gdename == "umikon"
replace gdename = "les verrieres" if gdename == "verrieres"
replace gdename = "visp" if gdename == "vieques"
replace gdename = "ollon" if gdename == "villars s/ollon"
replace gdename = "daerstetten" if gdename == "weissenburg"
replace gdename = "wil (sg)" if gdename == "wil (st. g.)"
replace gdename = "niederhelfenschwil" if gdename == "zuckenriet"
replace gdename = "allschwil" if gdename == "alischwil"
replace gdename = "gossau (sg)" if gdename == "arnegg gossau"
replace gdename = "arth" if gdename == "arth goldau"
replace gdename = "wiesendangen" if gdename == "attikon"
replace gdename = "basel" if gdename == "bale"
replace gdename = "le chenit" if gdename == "brassus"
replace gdename = "ruete" if gdename == "bruelisau"
replace gdename = "bussigny pres lausanne" if gdename == "bussigny/lausanne"
replace gdename = "le cerneux pequignot" if gdename == "cerneux pequignot"
replace gdename = "chateau d oex" if gdename == "chaeteau d oex"
replace gdename = "rochefort" if gdename == "chambrelien rochefort"
replace gdename = "chateau d oex" if gdename == "chateau d tex"
replace gdename = "commugny" if gdename == "commungy"
replace gdename = "courrendlin" if gdename == "courrendin"
replace gdename = "egerkingen" if gdename == "egerkinden"
replace gdename = "einsiedeln" if gdename == "euthal"
replace gdename = "littau" if gdename == "fluhmuehle"
replace gdename = "kirchberg (sg)" if gdename == "gaehwil"
replace gdename = "grosshoechstetten" if gdename == "grosshochstetten"
replace gdename = "volketswil" if gdename == "hegnau"
replace gdename = "altstaetten" if gdename == "hinterforst"
replace gdename = "ittigen" if gdename == "ittigen bolligen"
replace gdename = "horw" if gdename == "kastanienbaum horw"
replace gdename = "kuettigen" if gdename == "kuttigen"
replace gdename = "saint martin (vs)" if gdename == "la luette"
replace gdename = "lausanne" if gdename == "lausannne"
replace gdename = "lenzburg" if gdename == "lenzbourg"
replace gdename = "luetzelflueh" if gdename == "luetzelffueh"
replace gdename = "teufen (ar)" if gdename == "teufen"
replace gdename = "muenchenstein" if gdename == "miinchenstein"
replace gdename = "troistorrents" if gdename == "morgins"
replace gdename = "chateau d oex" if gdename == "moulins pres chateau d oex"
replace gdename = "niederwil (ag)" if gdename == "niederwil"
replace gdename = "koeniz" if gdename == "niederwangen koeniz"
replace gdename = "zug" if gdename == "oberwil zug"
replace gdename = "ormont dessus" if gdename == "ormond dessus"
replace gdename = "val de ruz" if gdename == "paquier"
replace gdename = "perly certoux" if gdename == "perly"
replace gdename = "faellanden" if gdename == "pfaffhausen"
replace gdename = "quinto" if gdename == "piotta di quinto"
replace gdename = "sion" if gdename == "pont de la morge/sion"
replace gdename = "lancy" if gdename == "pt lancy"
replace gdename = "regensberg" if gdename == "regenberg"
replace gdename = "riehen" if gdename == "riehen basel"
replace gdename = "zell (zh)" if gdename == "rikon zell"
replace gdename = "romanshorn" if gdename == "romanshom"
replace gdename = "sennwald" if gdename == "salez"
replace gdename = "s chanf" if gdename == "scanfs"
replace gdename = "seengen" if gdename == "seegen"
replace gdename = "signy avenex" if gdename == "signy"
replace gdename = "stansstad" if gdename == "stanstad"
replace gdename = "sutz lattrigen" if gdename == "sutz"
replace gdename = "turbenthal" if gdename == "thurbenthal"
replace gdename = "tramelan" if gdename == "tramelan dessous"
replace gdename = "unteraegeri" if gdename == "unterageri"
replace gdename = "vechigen" if gdename == "utzigen"
replace gdename = "vira (gambarogno)" if gdename == "vira cambarogno"
replace gdename = "koeniz" if gdename == "wabern koeniz"
replace gdename = "wallisellen" if gdename == "walisellen"
replace gdename = "wallisellen" if gdename == "walliselen"
replace gdename = "wynau" if gdename == "winau"
replace gdename = "wolfhalden" if gdename == "wollhalden"
replace gdename = "basel" if gdename == "base]"
replace gdename = "belmont sur lausanne" if gdename == "belmont s lausanne"
replace gdename = "bueren an der aare" if gdename == "bueren a. a."
replace gdename = "burgdorf" if gdename == "buergdorf"
replace gdename = "duebendorf" if gdename == "diibendorf"
replace gdename = "altstaetten" if gdename == "altstetten"
replace gdename = "silenen" if gdename == "amsteg"
replace gdename = "avenches" if gdename == "avanches"
replace gdename = "bellevue" if gdename == "bellevue geneve"
replace gdename = "beromuenster" if gdename == "beromunster"
replace gdename = "affoltern am albis" if gdename == "affoltern a. a."
replace gdename = "bettlach" if gdename == "bettfach"
replace gdename = "adelboden" if gdename == "adelhoden"
replace gdename = "adliswil" if gdename == "adliswii"
replace gdename = "ruedtligen alchenflueh" if gdename == "alchenflueh"
replace gdename = "allschwil" if gdename == "allscbwil"
replace gdename = "allschwil" if gdename == "allschwill"
replace gdename = "anieres" if gdename == "anieres geneve"
replace gdename = "arlesheim" if gdename == "arlesheim bl"
replace gdename = "arth" if gdename == "artb"
replace gdename = "bad ragaz" if gdename == "bad ragan"
replace gdename = "freienbach" if gdename == "baech freienbach"
replace gdename = "basel" if gdename == "basal"
replace gdename = "basel" if gdename == "base)"
replace gdename = "baetterkinden" if gdename == "batterkinden"
replace gdename = "beinwil am see" if gdename == "beinwil/see"
replace gdename = "altstaetten" if gdename == "altstiitten"
replace gdename = "nyon" if gdename == "amex sur nyon"
replace gdename = "bremgarten (ag)" if gdename == "anglikon"
replace gdename = "boudry" if gdename == "areuse boudry"
replace gdename = "arlesheim" if gdename == "arlesheiin"
replace gdename = "balsthal" if gdename == "balstal"
replace gdename = "basel" if gdename == "basel basel"
replace gdename = "basel" if gdename == "basle"
replace gdename = "beckenried" if gdename == "bechenried"
replace gdename = "bellach" if gdename == "belach"
replace gdename = "bern" if gdename == "ber n"
replace gdename = "bern" if gdename == "bern bumplitz"
replace gdename = "l abbaye" if gdename == "abbaye"
replace gdename = "thal" if gdename == "altenrhein"
replace gdename = "sainte croix" if gdename == "auberson"
replace gdename = "bassersdorf" if gdename == "basserdorf"
replace gdename = "bern" if gdename == "berne"
replace gdename = "waldkirch" if gdename == "bernhardzell"
replace gdename = "crans montana" if gdename == "bluche/randogne"
replace gdename = "bourg saint pierre" if gdename == "bourg st pierre"
replace gdename = "bre aldesago" if gdename == "bre"
replace gdename = "ingenbohl" if gdename == "brunnen"
replace gdename = "brusino arsizio" if gdename == "brusino"
replace gdename = "thunstetten" if gdename == "buetzberg"
replace gdename = "bussigny" if gdename == "bussigny/morges"
replace gdename = "brusio" if gdename == "campascio"
replace gdename = "brusio" if gdename == "campocologno"
replace gdename = "celerina/schlarigna" if gdename == "celerina"
replace gdename = "vernier" if gdename == "chatelaine"
replace gdename = "la chaux de fonds" if gdename == "chaux de fonds"
replace gdename = "ollon" if gdename == "chesieres s/ollon"
replace gdename = "ollon" if gdename == "chesieres/ollon"
replace gdename = "chene bourg" if gdename == "chine bourg"
replace gdename = "montreux" if gdename == "clarens"
replace gdename = "montreux" if gdename == "clarens montreux"
replace gdename = "corcelles le jorat" if gdename == "corcelles/payerne"
replace gdename = "corzoneso" if gdename == "corzonesco"
replace gdename = "crans montana" if gdename == "crans"
replace gdename = "crans montana" if gdename == "crans lens"
replace gdename = "crans montana" if gdename == "crans sur sierre"
replace gdename = "crans montana" if gdename == "crans/lens"
replace gdename = "crans pres celigny" if gdename == "crans/nyon"
replace gdename = "davos" if gdename == "davos wolfgang"
replace gdename = "rapperswil (be)" if gdename == "dieterswil"
replace gdename = "illnau effretikon" if gdename == "effretikon"
replace gdename = "emmen" if gdename == "emmebruecke"
replace gdename = "emmen" if gdename == "emmenbruecke"
replace gdename = "feldis/veulden" if gdename == "feldis"
replace gdename = "meilen" if gdename == "feldmeilen"
replace gdename = "flims" if gdename == "flims dorf"
replace gdename = "opfikon" if gdename == "glattbrugg"
replace gdename = "arth" if gdename == "goldau"
replace gdename = "lancy" if gdename == "grand lancy"
replace gdename = "le grand saconnex" if gdename == "grand saconnex"
replace gdename = "saanen" if gdename == "gstaad"
replace gdename = "muri bei bern" if gdename == "guemligen"
replace gdename = "muehleberg" if gdename == "guemmenen"
replace gdename = "sigriswil" if gdename == "gunten"
replace gdename = "nendaz" if gdename == "haute nendaz"
replace gdename = "val de ruz" if gdename == "hauts geneveys"
replace gdename = "buchholterberg" if gdename == "heimenschwand"
replace gdename = "waengi" if gdename == "heiterschen"
replace gdename = "bremgarten (ag)" if gdename == "hermetschwil"
replace gdename = "weggis" if gdename == "hertenstein"
replace gdename = "hoelstein" if gdename == "holstein"
replace gdename = "warth weiningen" if gdename == "ittingen"
replace gdename = "kaiseraugst" if gdename == "kaisersaugst"
replace gdename = "mauensee" if gdename == "kaltbach"
replace gdename = "koelliken" if gdename == "kolliken"
replace gdename = "koeniz" if gdename == "koniz"
replace gdename = "kreuzlingen" if gdename == "kreuzungen"
replace gdename = "wittenbach" if gdename == "kronbuehl"
replace gdename = "kuessnacht (sz)" if gdename == "kuessnacht a. r."
replace gdename = "sainte croix" if gdename == "l auberson"
replace gdename = "langendorf" if gdename == "langerdorf"
replace gdename = "langnau im emmental" if gdename == "langnau i.e."
replace gdename = "l abbaye" if gdename == "les bioux"
replace gdename = "montreux" if gdename == "les planches"
replace gdename = "koeniz" if gdename == "liebefeld"
replace gdename = "luzern" if gdename == "liftau"
replace gdename = "altstaetten" if gdename == "luechingen"
replace gdename = "luetisburg" if gdename == "luetisberg"
replace gdename = "lugano" if gdename == "lugano cassarate"
replace gdename = "lugano" if gdename == "lugano paradiso"
replace gdename = "teufen (ar)" if gdename == "lustmuehle"
replace gdename = "malters" if gdename == "maltera"
replace gdename = "maennedorf" if gdename == "mannedorf"
replace gdename = "matten bei interlaken" if gdename == "matten"
replace gdename = "mellingen" if gdename == "melligen"
replace gdename = "sigriswil" if gdename == "merligen"
replace gdename = "mezzovico vira" if gdename == "mezzovico"
replace gdename = "moeriken wildegg" if gdename == "moeriken"
replace gdename = "mont sur rolle" if gdename == "mont s/rolle"
replace gdename = "le mont sur lausanne" if gdename == "mont sur lausanne"
replace gdename = "morges" if gdename == "morgen"
replace gdename = "moeriken wildegg" if gdename == "moriken wildegg"
replace gdename = "muenchenstein" if gdename == "muechenstein"
replace gdename = "quarten" if gdename == "murg"
replace gdename = "muri bei bern" if gdename == "muri b. bern"
replace gdename = "lauterbrunnen" if gdename == "murren"
replace gdename = "allschwil" if gdename == "neu allschwil"
replace gdename = "allschwil" if gdename == "neuallschwil"
replace gdename = "neuchatel" if gdename == "neuchaetel"
replace gdename = "uzwil" if gdename == "niederuzwil"
replace gdename = "pambio noranco" if gdename == "noranco"
replace gdename = "vaz/obervaz" if gdename == "obervaz"
replace gdename = "oftringen" if gdename == "otringen"
replace gdename = "pailly" if gdename == "pally"
replace gdename = "pully" if gdename == "puffy"
replace gdename = "pully" if gdename == "pullt"
replace gdename = "pully" if gdename == "putty"
replace gdename = "chalais" if gdename == "rechy"
replace gdename = "richenthal" if gdename == "richental"
replace gdename = "nesslau krummenau" if gdename == "rietbad"
replace gdename = "romanel sur lausanne" if gdename == "romanel s/lausanne"
replace gdename = "romanshorn" if gdename == "romanhorn"
replace gdename = "kuettigen" if gdename == "rombach"
replace gdename = "ronco sopra ascona" if gdename == "ronco"
replace gdename = "risch" if gdename == "rotkreuz"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolfstetten"
replace gdename = "worb" if gdename == "ruefenacht"
replace gdename = "worb" if gdename == "rufenacht"
replace gdename = "rubigen" if gdename == "ruhigen"
replace gdename = "savosa" if gdename == "savosa paese"
replace gdename = "saanen" if gdename == "schoenried"
replace gdename = "schoetz" if gdename == "schuetz"
replace gdename = "scuol" if gdename == "schuls"
replace gdename = "bischofszell" if gdename == "schweizerholz"
replace gdename = "sempach" if gdename == "sempach station"
replace gdename = "selzach" if gdename == "setzach"
replace gdename = "schuebelbach" if gdename == "siebnen"
replace gdename = "zweisimmen" if gdename == "simmenthal"
replace gdename = "koeniz" if gdename == "spiegel"
replace gdename = "koeniz" if gdename == "spiegel koeniz"
replace gdename = "saint blaise" if gdename == "st blaise"
replace gdename = "saint cergue" if gdename == "st cergue"
replace gdename = "saint george" if gdename == "st georges"
replace gdename = "saint legier la chiesaz" if gdename == "st legier/la chiesaz"
replace gdename = "saint leonard" if gdename == "st leonard"
replace gdename = "saint prex" if gdename == "st prex"
replace gdename = "staefa" if gdename == "staef a"
replace gdename = "staufen" if gdename == "stauten"
replace gdename = "flums" if gdename == "tannenbodenalp"
replace gdename = "torricella taverne" if gdename == "taverne torricella"
replace gdename = "tenero contra" if gdename == "tenero"
replace gdename = "montreux" if gdename == "territet"
replace gdename = "thonex" if gdename == "thbnex"
replace gdename = "thun" if gdename == "thun gwatt"
replace gdename = "thun" if gdename == "thunersee"
replace gdename = "wartau" if gdename == "truebbach"
replace gdename = "wartau" if gdename == "truebbach/wartau"
replace gdename = "uetikon am see" if gdename == "uetikon a. see"
replace gdename = "wohlen bei bern" if gdename == "uettligen"
replace gdename = "uitikon" if gdename == "uitikon a a"
replace gdename = "uitikon" if gdename == "uitikon a. a."
replace gdename = "uitikon" if gdename == "uitikon a. a"
replace gdename = "uitikon" if gdename == "uitikon am albis"
replace gdename = "uitikon" if gdename == "uitikon waldegg"
replace gdename = "sion" if gdename == "uvrier sion"
replace gdename = "pfaefers" if gdename == "vaettis"
replace gdename = "vaz/obervaz" if gdename == "valbella"
replace gdename = "vandoeuvres" if gdename == "vandaeuvres"
replace gdename = "bagnes" if gdename == "verbier"
replace gdename = "chalais" if gdename == "vercorin"
replace gdename = "villars sainte croix" if gdename == "villars ste croix"
replace gdename = "ollon" if gdename == "villars sur ollon"
replace gdename = "villaz saint pierre" if gdename == "villaz st pierre"
replace gdename = "vuarrens" if gdename == "vuarrengel"
replace gdename = "scuol" if gdename == "vulpera"
replace gdename = "koeniz" if gdename == "wabern"
replace gdename = "wangen an der aare" if gdename == "wangen a. a."
replace gdename = "vilters wangs" if gdename == "wangs vilters"
replace gdename = "sumiswald" if gdename == "wasen i. e."
replace gdename = "grabs" if gdename == "werdenberg"
replace gdename = "wetzikon (zh)" if gdename == "wetzikon"
replace gdename = "moeriken wildegg" if gdename == "wildegg"
replace gdename = "winterthur" if gdename == "winter thur"
replace gdename = "wolhusen" if gdename == "wohlhusen"
replace gdename = "zuerich" if gdename == "zilrich"
replace gdename = "zuerich" if gdename == "zuerieh"
replace gdename = "zuerich" if gdename == "zurich"
replace gdename = "renens (vd)" if gdename == "renens"
replace gdename = "neyruz (fr)" if gdename == "neyruz"
replace gdename = "waedenswil" if gdename == "schoenenberg"
replace gdename = "birmensdorf (zh)" if gdename == "birmensdorf"
replace gdename = "rekingen (ag)" if gdename == "rekingen"
replace gdename = "oberriet (sg)" if gdename == "oberriet"
replace gdename = "emmen" if gdename == "emmenbruecke/emmen"
replace gdename = "carouge (ge)" if gdename == "carouge"
replace gdename = "horw" if gdename == "kastanienbaum/horw"
replace gdename = "ringgenberg (be)" if gdename == "ringgenberg"
replace gdename = "ingenbohl" if gdename == "brunnen/ingenbohl"
replace gdename = "littau" if gdename == "reussbuehl/littau"
replace gdename = "lachen" if gdename == "lachen sz"
replace gdename = "luzern" if gdename == "lucerne"
replace gdename = "mesocco" if gdename == "san bernardino/mesocco"
replace gdename = "freienbach" if gdename == "pfaeffikon sz"
replace gdename = "teufenthal (ag)" if gdename == "teufenthal"
replace gdename = "reckingen (vs)" if gdename == "reckingen"
replace gdename = "ried bei moerel" if gdename == "ried bei morel"
replace gdename = "kuesnacht (zh)" if gdename == "kuessnacht (zh)"
replace gdename = "grub (ar)" if gdename == "grub"
replace gdename = "zug" if gdename == "oberwil/zug"
replace gdename = "veltheim (ag)" if gdename == "veltheim"
replace gdename = "gaiserwald" if gdename == "abtwil/gaiserwald"
replace gdename = "risch" if gdename == "rotkreuz/risch"
replace gdename = "lauterbrunnen" if gdename == "wengen/lauterbrunnen"
replace gdename = "motiers (ne)" if gdename == "motiers"
replace gdename = "treytorrens (payerne)" if gdename == "treytorrens"
replace gdename = "koeniz" if gdename == "schliern"
replace gdename = "ruethi (rheintal)" if gdename == "ruethi"
replace gdename = "langnau am albis" if gdename == "langnau a. a"
replace gdename = "wartau" if gdename == "weite/wartau"
replace gdename = "birmenstorf (ag)" if gdename == "birmenstorf"
replace gdename = "duebendorf" if gdename == "duebenedorf"
replace gdename = "pfaeffikon" if gdename == "pfaeffikon zh"
replace gdename = "saint sulpice (vd)" if gdename == "st sulpice vd"
replace gdename = "thun" if gdename == "gwatt"
replace gdename = "langenthal" if gdename == "langenthal (be)"
replace gdename = "faellanden" if gdename == "faellenden"
replace gdename = "vaz/obervaz" if gdename == "valbella/vaz/obervaz"
replace gdename = "flums" if gdename == "flumserberg"
replace gdename = "reute (ar)" if gdename == "reute"
replace gdename = "sierre" if gdename == "siders"
replace gdename = "ormont dessus" if gdename == "les diablerets/ormont dessus"
replace gdename = "werthenstein" if gdename == "schachen/werthenstein"
replace gdename = "corcelles cormondreche" if gdename == "corcelles ne"
replace gdename = "menzingen" if gdename == "edlibach/menzingen"
replace gdename = "gebenstorf" if gdename == "gebensdorf"
replace gdename = "kirchlindach" if gdename == "herrenschwanden/kirchlindach"
replace gdename = "oberstammheim" if gdename == "stammheim"
replace gdename = "riehen" if gdename == "riehen (basel)"
replace gdename = "uitikon" if gdename == "ultikon"
replace gdename = "morrens (vd)" if gdename == "morrens"
replace gdename = "rueschlikon" if gdename == "rueschlikon (zh)"
replace gdename = "hilterfingen" if gdename == "huenibach hilterfingen"
replace gdename = "kuessnacht (sz)" if gdename == "kuessnacht sz"
replace gdename = "freienbach" if gdename == "pfaeffikon/freienbach"
replace gdename = "koeniz" if gdename == "wabern/koeniz"
replace gdename = "hausen (ag)" if gdename == "hausen b. brugg"
replace gdename = "duedingen" if gdename == "mariahilf"
replace gdename = "horn" if gdename == "horn tg"
replace gdename = "schwyz" if gdename == "rickenbach/schwyz"
replace gdename = "illnau" if gdename == "ilinau effretikon"
replace gdename = "wuerenlos" if gdename == "wrenlos"
replace gdename = "saignelegier" if gdename == "saignelgier"
replace gdename = "muestair" if gdename == "mstair"
replace gdename = "hausen am albis" if gdename == "hausen a. a."
replace gdename = "unterbaech" if gdename == "unterbch"
replace gdename = "fluehli" if gdename == "soerenberg/fluehli"
replace gdename = "thal" if gdename == "staad/thal"
replace gdename = "dietikon" if gdename == "dietikon zh"
replace gdename = "interlaken" if gdename == "matten/interlaken"
replace gdename = "le locle" if gdename == "le locie"
replace gdename = "freienbach" if gdename == "freienbach sz"
replace gdename = "winterthur" if gdename == "winterthour"
replace gdename = "lohn ammannsegg" if gdename == "ammansegg"


sort gdename year
drop gdenr gdenr_2018 ctn E_CNTR N_CNTR geo_merge
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2
gen geo_merge = 2 if _merge == 3

preserve
keep if _merge == 3
save "$path\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_2.dta", replace
restore


* 2) b) String merge on further standardized company-gdename info (step 2: 300 biggest cases)
keep if _merge == 1
drop _merge

replace gdename = "yverdon les bains" if gdename == "yverdon"
replace gdename = "st gallen" if gdename == "gallen"
replace gdename = "montana" if gdename == "crans montana"
replace gdename = "chur" if gdename == "coira"
replace gdename = "kreuzlingen" if gdename == "kreuzlangen"
replace gdename = "wil sg" if gdename == "wi sg"
replace gdename = "neuchatel" if gdename == "neuchtel"
replace gdename = "charmey" if gdename == "charmey (gruyere)"
replace gdename = "oberglatt" if gdename == "oberglatt zh"
replace gdename = "corcelles cormondreche" if gdename == "corcelles cormondreche ne"
replace gdename = "eschlikon" if gdename == "eschlikon tg"
replace gdename = "bruegg" if gdename == "bruegg bei biel"
replace gdename = "baar" if gdename == "baer"
replace gdename = "zug" if gdename == "zug pr"
replace gdename = "geneve" if gdename == "geneve vp"
replace gdename = "wil sg" if gdename == "wii sg"
replace gdename = "biel be" if gdename == "biel bienne"
replace gdename = "geneve" if gdename == "geneve vr"
replace gdename = "chatel saint denis" if gdename == "chaetel st denis"
replace gdename = "solothurn" if gdename == "solothum"
replace gdename = "lugano" if gdename == "lugano vr"
replace gdename = "lugano" if gdename == "lugano pr"
replace gdename = "alterswil" if gdename == "alterswil fr"
replace gdename = "wuennewil flamatt" if gdename == "wuennewil/flamatt"
replace gdename = "sainte croix" if gdename == "st croix"
replace gdename = "saint cierges" if gdename == "st cierges"
replace gdename = "saint aubin fr" if gdename == "st aubin fr"
replace gdename = "sorengo" if gdename == "sorengo pr"
replace gdename = "sorengo" if gdename == "sorengo s"
replace gdename = "sorengo" if gdename == "sorengo vr"
replace gdename = "solothurn" if gdename == "solo thurn"
replace gdename = "steg" if gdename == "steg vs"
replace gdename = "luzern" if gdename == "luzem"
replace gdename = "cornaux" if gdename == "cornaux ne"
replace gdename = "augst" if gdename == "augst bl"
replace gdename = "studen" if gdename == "studen be"
replace gdename = "kuessnacht am rigi" if gdename == "kuessnacht (sz)"
replace gdename = "laax" if gdename == "laax gr"
replace gdename = "niederglatt" if gdename == "niederglatt zh"
replace gdename = "lenk" if gdename == "lenk im simmental"
replace gdename = "biel benken" if gdename == "biel benken bl"
replace gdename = "bruegg" if gdename == "bruegg be"
replace gdename = "baar" if gdename == "bear"
replace gdename = "laupen" if gdename == "laupen be"
replace gdename = "geneve" if gdename == "geneve pr"
replace gdename = "altstaetten" if gdename == "altstaetten sg"
replace gdename = "martigny" if gdename == "martigny bourg"
replace gdename = "martigny" if gdename == "martigny ville"
replace gdename = "fribourg" if gdename == "fri bourg"
replace gdename = "geneve" if gdename == "ge neve"
replace gdename = "biel be" if gdename == "biellbienne"
replace gdename = "stetten ag" if gdename == "steffen ag"
replace gdename = "tentlingen" if gdename == "tentlingen/tinterin fr"
replace gdename = "geneve" if gdename == "geneve s"
replace gdename = "charmey" if gdename == "charmey (gruyere)"
replace gdename = "egg" if gdename == "egg bei zuerich"
replace gdename = "castel san pietro" if gdename == "pietro"
replace gdename = "la sagne" if gdename == "la sagne ne"
replace gdename = "la tour de peilz" if gdename == "la tourde peilz"
replace gdename = "la tour de treme" if gdename == "la tourde treme"
replace gdename = "lugano" if gdename == "lugano vp"
replace gdename = "zuerich" if gdename == "zuench"
replace gdename = "lugano" if gdename == "lugano s"
replace gdename = "payerne" if gdename == "payeme"
replace gdename = "saint luc" if gdename == "st luc"
replace gdename = "la joux fr" if gdename == "lajoux"
replace gdename = "riva san vitale" if gdename == "vitale"
replace gdename = "lugano" if gdename == "castagnola"
replace gdename = "sant antonino" if gdename == "antonino"
replace gdename = "lens" if gdename == "de lens"
replace gdename = "hagenbuch" if gdename == "hagenbuch zh"
replace gdename = "balerna" if gdename == "balema"
replace gdename = "seewen" if gdename == "seewen so"
replace gdename = "le locle" if gdename == "locle"
replace gdename = "geneve" if gdename == "genese"
replace gdename = "biel be" if gdename == "biel!bienne"
replace gdename = "ebnat kappel" if gdename == "ebnat kappet"
replace gdename = "ebnat kappel" if gdename == "ebnet kappel"
replace gdename = "sirnach" if gdename == "simach"
replace gdename = "manno" if gdename == "manne"
replace gdename = "manno" if gdename == "manna"
replace gdename = "marly" if gdename == "manly"
replace gdename = "igis" if gdename == "/gis"
replace gdename = "reiden" if gdename == "beiden"
replace gdename = "icogne" if gdename == "(cogne"
replace gdename = "uster" if gdename == "oster"
replace gdename = "stalden vs" if gdename == "stalden"
replace gdename = "la sagne" if gdename == "la sagne ne"
replace gdename = "chiasso" if gdename == "chiasso vr"
replace gdename = "chiasso" if gdename == "chiasso pr"
replace gdename = "muenchwilen ag" if gdename == "muenchwilenag"
replace gdename = "muenchwilen tg" if gdename == "muenchwilentg"
replace gdename = "uster" if gdename == "lister"
replace gdename = "studen" if gdename == "studen bei bruegg"
replace gdename = "klosters serneus" if gdename == "klosters"
replace gdename = "klosters serneus" if gdename == "klosters dorf"
replace gdename = "klosters serneus" if gdename == "klosters platz"
replace gdename = "klosters serneus" if gdename == "klosters semeus"
replace gdename = "klosters serneus" if gdename == "klosterssemeus"
replace gdename = "klosters serneus" if gdename == "klostersserneus"
replace gdename = "illnau effretikon" if gdename == "illnau"
replace gdename = "wuennewil flamatt" if gdename == "wuennewil"
replace gdename = "st antoni" if gdename == "antoni"
replace gdename = "schoetz" if gdename == "schaetz"
replace gdename = "igis" if gdename == "gis"
replace gdename = "san vittore" if gdename == "vittore"
replace gdename = "hofstetten flueh" if gdename == "hofstetten flueh so"
replace gdename = "hofstetten so" if gdename == "hofstetten flueh"
replace gdename = "hofstetten so" if gdename == "hofsteffen flueh so"
replace gdename = "taeuffelen" if gdename == "taeuffelen gerolfingen"
replace gdename = "basel" if gdename == "base!"
replace gdename = "bussigny pres lausanne" if gdename == "bussigny"
replace gdename = "boudevilliers" if gdename == "de boudevilliers"
replace gdename = "boudry" if gdename == "de boudry"
replace gdename = "chalais" if gdename == "de chalais"
replace gdename = "chermignon" if gdename == "de chermignon"
replace gdename = "constantine" if gdename == "de constantine"
replace gdename = "corcelles cormondreche" if gdename == "de corcelles cormond"
replace gdename = "corcelles cormondreche" if gdename == "de corcelles cormondrec"
replace gdename = "corcelles cormondreche" if gdename == "de corcelles cormondrech"
replace gdename = "finhaut" if gdename == "de finhaut"
replace gdename = "la chaux de fonds" if gdename == "de fonds"
replace gdename = "l abbaye" if gdename == "de l abbaye"
replace gdename = "lutry" if gdename == "de lutry"
replace gdename = "marin epagnier" if gdename == "de marin epagnier"
replace gdename = "montana" if gdename == "de montana"
replace gdename = "montreux" if gdename == "de montreux"
replace gdename = "ollon" if gdename == "de ollon"
replace gdename = "randogne" if gdename == "de randogne"
replace gdename = "sion" if gdename == "de sion"
replace gdename = "saint aubin sauges" if gdename == "de st aubin sauges"
replace gdename = "saint imier" if gdename == "de st imier"
replace gdename = "saint martin" if gdename == "de st martin"
replace gdename = "sainte croix" if gdename == "de ste croix"
replace gdename = "duedingen" if gdename == "ddingen"
replace gdename = "duebendorf" if gdename == "dbendorf"
replace gdename = "dornach" if gdename == "darnach"
replace gdename = "daellikon" if gdename == "dallikon"
replace gdename = "daenikon" if gdename == "daenikon zh"
replace gdename = "daellikon" if gdename == "daeilikon"
replace gdename = "ormont dessus" if gdename == "d ormont dessus"
replace gdename = "ormont dessous" if gdename == "d ormont dessous"
replace gdename = "ollon" if gdename == "d ollon"
replace gdename = "icogne" if gdename == "d icogne"
replace gdename = "cumbel" if gdename == "cumbel/cumbels"
replace gdename = "freienstein teufen" if gdename == "freienstein"
replace gdename = "st gallen" if gdename == "sl gallen"
replace gdename = "sierre" if gdename == "slerre"
replace gdename = "sion" if gdename == "slon"
replace gdename = "fluehli" if gdename == "soerenberg"
replace gdename = "fluehli" if gdename == "soerenberg fluehli"
replace gdename = "solothurn" if gdename == "soiothurn"
replace gdename = "st peterzell" if gdename == "peterzell"
replace gdename = "st peterzell" if gdename == "peterszell"
replace gdename = "igis" if gdename == "landquart"
replace gdename = "thonex" if gdename == "thnex"
replace gdename = "thonex" if gdename == "thdnex"
replace gdename = "thalwil" if gdename == "thalwil pr"
replace gdename = "thalwil" if gdename == "thalwii"
replace gdename = "bilten" if gdename == "bitten"
replace gdename = "delemont" if gdename == "emont"
replace gdename = "evilard" if gdename == "evilard (leubringen)"
replace gdename = "langnau be" if gdename == "langhau im emmental"
replace gdename = "vaz/obervaz" if gdename == "lenzerheide vaz"
replace gdename = "vaz/obervaz" if gdename == "gde. vaz/obervaz"
replace gdename = "vaz/obervaz" if gdename == "lenzerheide vaz/obervaz"
replace gdename = "vaz/obervaz" if gdename == "lenzerheide/vaz/obervaz"
replace gdename = "bolligen" if gdename == "ostermundigen" //ostermundingen was part of bolligen (with Ittigen) until 1983, they became distinct mun after that year
replace gdename = "bolligen" if gdename == "ostermundigen bolligen" & year < 1984
replace gdename = "birmensdorf zh" if gdename == "birmensdorfzh"
replace gdename = "plaffeien" if gdename == "plaffeien schwarzsee"
replace gdename = "plaffeien" if gdename == "plaffeien/schwarzsee"
replace gdename = "plaffeien" if gdename == "plaffelen"
replace gdename = "plaffeien" if gdename == "plaffeyen"
replace gdename = "plaffeien" if gdename == "plaffeyen/schwarzsee"
replace gdename = "geneve" if gdename == "plainpalais"
replace gdename = "plan les ouates" if gdename == "plan les ouates pr"
replace gdename = "plan les ouates" if gdename == "plan les quates"
replace gdename = "plan les ouates" if gdename == "planes ouates"
replace gdename = "plan les ouates" if gdename == "plans les ouates"
replace gdename = "plan les ouates" if gdename == "plantes ouates"
replace gdename = "wahlern" if gdename == "schwarzenburg"
replace gdename = "wahlern" if gdename == "schwarzenburg wahlern"
replace gdename = "plaffeien" if gdename == "schwarzsee"
replace gdename = "plaffeien" if gdename == "schwarzsee plaffeien"
replace gdename = "plaffeien" if gdename == "schwarzsee/plaffeien"
replace gdename = "sils im engadin/segl" if gdename == "sils/segl im engadin"
replace gdename = "sils im engadin/segl" if gdename == "sils im engadin"
replace gdename = "sils im engadin/segl" if gdename == "sils/segi im engadin"
replace gdename = "st ursen" if gdename == "ursen"
replace gdename = "st gallen" if gdename == "galien"
replace gdename = "muenster vs" if gdename == "muenster"
replace gdename = "toffen" if gdename == "tollen"
replace gdename = "igis" if gdename == "(gis"
replace gdename = "igis" if gdename == "!gis"
replace gdename = "icogne" if gdename == "!cogne"
replace gdename = "geneve" if gdename == "ageneve"
replace gdename = "bern" if gdename == "agbern"
replace gdename = "basel" if gdename == "agbasel"
replace gdename = "zug" if gdename == "ag zug"
replace gdename = "zuerich" if gdename == "ag zuerich"
replace gdename = "luzern" if gdename == "ag luzern"
replace gdename = "bern" if gdename == "ag bern"
replace gdename = "basel" if gdename == "ag basel"
replace gdename = "fribourg" if gdename == "afribourg"
replace gdename = "affoltern am albis" if gdename == "affoltem am albis"
replace gdename = "affoltern am albis" if gdename == "affoitern am albis"
replace gdename = "affoltern am albis" if gdename == "affoftem am albis"
replace gdename = "affoltern am albis" if gdename == "affoftern am albis"
replace gdename = "affoltern am albis" if gdename == "affoltern a a"
replace gdename = "affoltern am albis" if gdename == "affoltern am aibis"
replace gdename = "affoltern am albis" if gdename == "affoltern am albfis"
replace gdename = "affoltern am albis" if gdename == "affoltern am albin"
replace gdename = "affoltern am albis" if gdename == "affoltern am albis vr"
replace gdename = "affoltern am albis" if gdename == "affottem am albis"
replace gdename = "affoltern am albis" if gdename == "affoltern/albis"
replace gdename = "affoltern am albis" if gdename == "affoltern am albs"
replace gdename = "buerglen tg" if gdename == "buergten tg"
replace gdename = "ruefenach" if gdename == "ruefenach ag"
replace gdename = "le chenit" if gdename == "le chenil"
replace gdename = "gsteig" if gdename == "gsteig bei gstaad"
replace gdename = "wettswil am albis" if gdename == "wettswil"
replace gdename = "st gallen" if gdename == "gallen pr"
replace gdename = "st gallenkappel" if gdename == "gallenkappel"
replace gdename = "arni ag" if gdename == "arni islisberg"
replace gdename = "arni be" if gdename == "arni bei biglen"
replace gdename = "arni be" if gdename == "arni"
replace gdename = "st stephan" if gdename == "stephan"
replace gdename = "kuesnacht zh" if gdename == "kuesnachtzh"
replace gdename = "morges" if gdename == "merges"
replace gdename = "winkel" if gdename == "winkel bei buelach"
replace gdename = "lausanne" if gdename == "lausanne pr"
replace gdename = "roggwil tg" if gdename == "roggwiltg"
replace gdename = "jaun" if gdename == "bellegarde"
replace gdename = "jaun" if gdename == "bellegarde/vilette"
replace gdename = "bellinzona" if gdename == "bel linzona"
replace gdename = "ottenbach" if gdename == "offenbach"
replace gdename = "ottenbach" if gdename == "offenbach zh"
replace gdename = "worb" if gdename == "warb"
replace gdename = "wassen" if gdename == "wassen ur"
replace gdename = "waedenswil" if gdename == "wdenswil"
replace gdename = "weinfelden" if gdename == "weinfeiden"
replace gdename = "weinfelden" if gdename == "weinfelder"
replace gdename = "weinfelden" if gdename == "weinfeldes"
replace gdename = "weinfelden" if gdename == "weinfellen"
replace gdename = "knutwil" if gdename == "bad knutwil"
replace gdename = "montreux" if gdename == "montreux chaetelard"
replace gdename = "montreux" if gdename == "montreux chatelard"
replace gdename = "montreux" if gdename == "montreux le chaetelard"
replace gdename = "montreux" if gdename == "montreux le chatelard"
replace gdename = "montreux" if gdename == "montreux planches"
replace gdename = "montreux" if gdename == "montreux territet"
replace gdename = "montreux" if gdename == "montreux territet (les planches)"
replace gdename = "montreux" if gdename == "montreux vd"
replace gdename = "montreux" if gdename == "montreux vr"
replace gdename = "moerel" if gdename == "morel"
replace gdename = "torricella taverne" if gdename == "taveme torricella"
replace gdename = "torricella taverne" if gdename == "taverne"
replace gdename = "torricella taverne" if gdename == "taverne toricella"
replace gdename = "torricella taverne" if gdename == "taverne torriceila"
replace gdename = "ruedtligen alchenflueh" if gdename == "ruedtligen"
replace gdename = "monteggio" if gdename == "molinazzo di monteggio"
replace gdename = "rickenbach zh" if gdename == "rickenbach attikon"
replace gdename = "ried bei brig" if gdename == "ried brig"
replace gdename = "neuhausen am rheinfall" if gdename == "rh"
replace gdename = "luthern" if gdename == "luthern dorf"
replace gdename = "sumvitg" if gdename == "sumvitg/somvix"
replace gdename = "leuk" if gdename == "susten"
replace gdename = "leuk" if gdename == "susten leuk"
replace gdename = "leuk" if gdename == "susten gmd leuk"
replace gdename = "leuk" if gdename == "susten, gde. leuk"
replace gdename = "leuk" if gdename == "susten. gde. leuk"
replace gdename = "diessenhofen" if gdename == "willisdorf"
replace gdename = "willisdorf" if gdename == "willisdorf diessenhofen"
replace gdename = "wil sg" if gdename == "wit sg"
replace gdename = "wil zh" if gdename == "wit zh"
replace gdename = "kappel am albis" if gdename == "kappet am albis"
replace gdename = "kappel so" if gdename == "kappet so"
replace gdename = "igis" if gdename == "lgis"
replace gdename = "le lieu" if gdename == "lieu"
replace gdename = "kaiserstuhl" if gdename == "kaiserstuhl ag"
replace gdename = "buettikon" if gdename == "buettikon ag"
replace gdename = "renan be" if gdename == "renan"
replace gdename = "ittigen" if gdename == "lttigen"
replace gdename = "morges" if gdename == "marges"
replace gdename = "hasle lu" if gdename == "hasle"
replace gdename = "haslen" if gdename == "haslen al"
replace gdename = "haslen" if gdename == "haslen gl"
replace gdename = "st niklaus" if gdename == "niklaus vs"
replace gdename = "saint martin vs" if gdename == "st martin vs"
replace gdename = "saint martin fr" if gdename == "st martin fr"
replace gdename = "st margrethen" if gdename == "st margrethen sg"
replace gdename = "saint livres" if gdename == "st livres"
replace gdename = "saint jean" if gdename == "st jean"
replace gdename = "saint jean" if gdename == "st jean vs"
replace gdename = "saint gingolph" if gdename == "st gingolph"
replace gdename = "saint prex" if gdename == "st prei"
replace gdename = "saint saphorin (lavaux)" if gdename == "st saphorin (lavaux)"
replace gdename = "rohrbach" if gdename == "rohrbach bei huttwil"
replace gdename = "brig glis" if gdename == "glis"
replace gdename = "st silvester" if gdename == "silvester"
replace gdename = "unterbaech" if gdename == "unterbaech vs"
replace gdename = "wila" if gdename == "wile"
replace gdename = "domat/ems" if gdename == "domat ems"
replace gdename = "brig glis" if gdename == "brig glas"
replace gdename = "sisseln" if gdename == "sisseln ag"
replace gdename = "san nazzaro" if gdename == "nazzaro"
replace gdename = "vella" if gdename == "vella gr"
replace gdename = "lossy formangueires" if gdename == "lossy"
replace gdename = "ipsach" if gdename == "lpsach"
replace gdename = "birmenstorf ag" if gdename == "birmenstorfag"
replace gdename = "mesocco" if gdename == "san bernardino"
replace gdename = "arni ag" if gdename == "arni lslisberg"
replace gdename = "boenigen" if gdename == "boenigen bei interlaken"
replace gdename = "schuepfen" if gdename == "schoepfen"
replace gdename = "schuepfheim" if gdename == "schoepfheim"
replace gdename = "sarnen" if gdename == "sannen"
replace gdename = "freienbach" if gdename == "pfaeffikon freienbach"
replace gdename = "saint sulpice ne" if gdename == "st sulpice ne"
replace gdename = "gipf oberfrick" if gdename == "gipf 0berfrick"
replace gdename = "gipf oberfrick" if gdename == "gipf oberfrack"
replace gdename = "gipf oberfrick" if gdename == "gipf obertrick"
replace gdename = "verscio" if gdename == "verscio pedemonte"
replace gdename = "liestal" if gdename == "vestal"
replace gdename = "vaz/obervaz" if gdename == "vaz obervaz"
replace gdename = "bolligen" if gdename == "ittigen" // Ittigen was part of Bolligen before 1983, year of their divisions
replace gdename = "lauperswil" if gdename == "gde. lauperswil"
replace gdename = "lauterbrunnen" if gdename == "gde. lauterbrunnen"
replace gdename = "wahlen" if gdename == "wahlen ag"
replace gdename = "meyrin" if gdename == "meyrin vr"
replace gdename = "meyrin" if gdename == "meyrin vp"
replace gdename = "meyrin" if gdename == "meyrin pr"
replace gdename = "moehlin" if gdename == "mhlin"
replace gdename = "muenchenstein" if gdename == "mnchenstein"
replace gdename = "muensingen" if gdename == "mnsingen"
replace gdename = "maennedorf" if gdename == "mnnedorf"
replace gdename = "kuessnacht am rigi" if gdename == "kuessnachtam rigi"
replace gdename = "kuessnacht am rigi" if gdename == "rigi"
replace gdename = "oberbueren" if gdename == "oberbaeren"
replace gdename = "lauenen" if gdename == "lauenen bei gstaad"
replace gdename = "schattdorf" if gdename == "schaftdorf"
replace gdename = "schattdorf" if gdename == "schaltdorf"
replace gdename = "zuerich" if gdename == "zuerich vp"
replace gdename = "thunstetten" if gdename == "gde. thunstetten"
replace gdename = "baeriswil" if gdename == "baeriswil be"
replace gdename = "forel (lavaux)" if gdename == "forel (lavauz)"
replace gdename = "prangins" if gdename == "frangins"
replace gdename = "arbon" if gdename == "frasnacht"
replace gdename = "frasnacht" if gdename == "frasnacht arbon"
replace gdename = "fribourg" if gdename == "freiburg"
replace gdename = "freienbach" if gdename == "frelenbach"
replace gdename = "oberhofen bei kreuzlingen" if gdename == "oberhofen bei kreuzlangen"
replace gdename = "saanen" if gdename == "gstaad saanen"
replace gdename = "niederurnen" if gdename == "niederumen"
replace gdename = "guettingen" if gdename == "goettingen"
replace gdename = "roche vd" if gdename == "roche"
replace gdename = "thal" if gdename == "thai"
replace gdename = "le chenit" if gdename == "du chenit"
replace gdename = "kreuzlingen" if gdename == "kreuzsingen"
replace gdename = "eppenberg woeschnau" if gdename == "eppenberg"

sort gdename year
drop gdenr gdenr_2018 ctn E_CNTR N_CNTR geo_merge
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2
gen geo_merge = 3 if _merge == 3

preserve
keep if _merge == 3
save "$path\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_3.dta", replace
restore

* 2) c) String merge on further standardized company-gdename info (step 3)
keep if _merge == 1
drop _merge

replace gdename = "sierre" if gdename == "granges vs"
replace gdename = "sierre" if gdename == "granges/sierre"
replace gdename = "granges pres marnand" if gdename == "granges pres mamand"
replace gdename = "granges pres marnand" if gdename == "granges/marnand"
replace gdename = "stadel" if gdename == "stadel bei niederglatt"
replace gdename = "stabio" if gdename == "stablo"
replace gdename = "stabio" if gdename == "stable"
replace gdename = "stabio" if gdename == "stabia"
replace gdename = "lausanne" if gdename == "lausanne vr"
replace gdename = "muttenz" if gdename == "mutlenz"
replace gdename = "peseux" if gdename == "peseux pr"
replace gdename = "st peter" if gdename == "peter"
replace gdename = "uzwil" if gdename == "henau uzwil"
replace gdename = "neftenbach" if gdename == "neffenbach"
replace gdename = "zug" if gdename == "a zug"
replace gdename = "lugano" if gdename == "a lugano"
replace gdename = "lausanne" if gdename == "a lausanne"
replace gdename = "geneve" if gdename == "a geneve"
replace gdename = "fribourg" if gdename == "a fribourg"
replace gdename = "la chaux de fonds" if gdename == "a chaux de fonds"
replace gdename = "coire" if gdename == "a chur"
replace gdename = "basel" if gdename == "a basel"
replace gdename = "saint maurice" if gdename == "maurice"
replace gdename = "toffen" if gdename == "totten"
replace gdename = "palezieux" if gdename == "palesieux"
replace gdename = "palezieux" if gdename == "palezieux gare"
replace gdename = "palezieux" if gdename == "palezieux vd"
replace gdename = "palezieux" if gdename == "palezieuxvd"
replace gdename = "pambio noranco" if gdename == "pambio noranco pr"
replace gdename = "pambio noranco" if gdename == "pamblo noranco"
replace gdename = "nuglar st pantaleon" if gdename == "pantaleon"
replace gdename = "rieden" if gdename == "rieden sg"
replace gdename = "rieden" if gdename == "riedem"
replace gdename = "ried bei moerel" if gdename == "ried moerel"
replace gdename = "stansstad" if gdename == "fuerigen"
replace gdename = "stansstad" if gdename == "fuerigen standsstad"
replace gdename = "stansstad" if gdename == "fuerigen stansstad"
replace gdename = "stansstad" if gdename == "fuerigen/stansstad"
replace gdename = "innertkirchen" if gdename == "lnnertkirchen"
replace gdename = "vilters wangs" if gdename == "vilters"
replace gdename = "vilters" if gdename == "vilters wangs" & year < 1996
replace gdename = "fribourg" if gdename == "fribourg vr"
replace gdename = "basel" if gdename == "basel vr"
replace gdename = "zuerich" if gdename == "zue rich"
replace gdename = "zug" if gdename == "zug vr"
replace gdename = "basel" if gdename == "basei"
replace gdename = "zuerich" if gdename == "zuerich s"
replace gdename = "basel" if gdename == "basel pr"
replace gdename = "chiasso" if gdename == "sachiasso"
replace gdename = "sachseln" if gdename == "sachsein"
replace gdename = "sachseln" if gdename == "sacheln"
replace gdename = "sachseln" if gdename == "sachslen"
replace gdename = "delemont" if gdename == "sadelemont"
replace gdename = "fribourg" if gdename == "safribourg"
replace gdename = "sainte croix" if gdename == "saint croix"
replace gdename = "saint blaise" if gdename == "sainte blaise"
replace gdename = "sales (gruyere)" if gdename == "sales" & year > 1976
replace gdename = "sales" if gdename == "sales (gruyere)" & year > 2000
replace gdename = "lugano" if gdename == "salugano"
replace gdename = "geneve" if gdename == "gereve"
replace gdename = "la chaux de fonds" if gdename == "la chaux de fond"
replace gdename = "lausanne" if gdename == "lausanne s"
replace gdename = "zug" if gdename == "zug[p]"
replace gdename = "lugano" if gdename == "a„ lugano"
replace gdename = "luzern" if gdename == "lu zern"
replace gdename = "geneve" if gdename == "geneva"
replace gdename = "geneve" if gdename == "geneve[p]"
replace gdename = "paradiso" if gdename == "paradiso vr"
replace gdename = "manno" if gdename == "manno pr"
replace gdename = "zug" if gdename == "zug vp"
replace gdename = "locarno" if gdename == "locarno s"
replace gdename = "zuerich" if gdename == "zu rich"
replace gdename = "locarno" if gdename == "locarno pr"
replace gdename = "massagno" if gdename == "massagno pr"
replace gdename = "luzern" if gdename == "luzern pr"
replace gdename = "bern" if gdename == "bern[p]"
replace gdename = "fribourg" if gdename == "fnbourg"
replace gdename = "lausanne" if gdename == "lausanne[p]"
replace gdename = "la chaux de fonds" if gdename == "la chaux de fonds pr"
replace gdename = "lausanne" if gdename == "lausanne vp"
replace gdename = "basel" if gdename == "basel vp"
replace gdename = "basel" if gdename == "basel s"
replace gdename = "la chaux de fonds" if gdename == "la chauxde fonds"
replace gdename = "kloten" if gdename == "kioten"
replace gdename = "zuerich" if gdename == "zuerich[p]"
replace gdename = "zuerich" if gdename == "zuerich vr"
replace gdename = "geneve" if gdename == "a„ geneve"
replace gdename = "bern" if gdename == "bern pr"
replace gdename = "zug" if gdename == "zug s"
replace gdename = "zuerich" if gdename == "zuericn"
replace gdename = "locarno" if gdename == "locamo"
replace gdename = "neuchatel" if gdename == "neu chatel"
replace gdename = "buerglen tg" if gdename == "buerglentg"
replace gdename = "sils im engadin/segl" if gdename == "/segl"
replace gdename = "igis" if gdename == "/ igis"
replace gdename = "biel be" if gdename == "/bienne"
replace gdename = "lens" if gdename == "/lens"
replace gdename = "zuerich" if gdename == "(vzm) zuerich"
replace gdename = "zuerich" if gdename == "(stahag) zuerich"
replace gdename = "geneve" if gdename == "(geneve"
replace gdename = "lausanne" if gdename == "(enusa) lausanne"
replace gdename = "vevey" if gdename == "(curateur) vevey"
replace gdename = "berneck" if gdename == "bemeck"
replace gdename = "glarus" if gdename == "glaris"
replace gdename = "glarus" if gdename == "glaru"
replace gdename = "illnau effretikon" if gdename == "lllnau effretikon"
replace gdename = "liestal" if gdename == "llestal"
replace gdename = "bole" if gdename == "boele"
replace gdename = "belp" if gdename == "beip"
replace gdename = "sierre" if gdename == "sierra"
replace gdename = "sierre" if gdename == "sierrelsiders"
replace gdename = "pfeffikon" if gdename == "pfeffikon lu"
replace gdename = "wahlern" if gdename == "wahlem"
replace gdename = "wolfhalden" if gdename == "wolf halden"
replace gdename = "chermignon" if gdename == "com. de chermignon"
replace gdename = "kuesnacht zh" if gdename == "ksnacht"
replace gdename = "scuol/schuls" if gdename == "scuol"
replace gdename = "chene bourg" if gdename == "chne bourg"
replace gdename = "neuhausen am rheinfall" if gdename == "rheinfall"
replace gdename = "fischingen" if gdename == "dussnang"
replace gdename = "delemont" if gdename == "delmont"
replace gdename = "chiasso" if gdename == "chlasso"
replace gdename = "chene bougeries" if gdename == "chne bougeries"
replace gdename = "kuessnacht am rigi" if gdename == "kussnacht am rigi"
replace gdename = "ruemikon" if gdename == "ruemikon ag"
replace gdename = "thal" if gdename == "staad"
replace gdename = "thal" if gdename == "staad thai"
replace gdename = "thal" if gdename == "staad thal"
replace gdename = "bertschikon" if gdename == "bertschikon bei attikon"
replace gdename = "muntelier" if gdename == "montilliez"
replace gdename = "mels" if gdename == "mets"
replace gdename = "saint sulpice vd" if gdename == "stsulpice vd"
replace gdename = "egnach" if gdename == "neukirch egnach"
replace gdename = "saint sulpice vd" if gdename == "st sulpicevd"
replace gdename = "eisch" if gdename == "rotkreuz risch"
replace gdename = "lancy" if gdename == "laney"
replace gdename = "le chenit" if gdename == "chenit"
replace gdename = "einsiedeln" if gdename == "einsiedein"
replace gdename = "einsiedeln" if gdename == "einsielden"
replace gdename = "elsau" if gdename == "eisau"
replace gdename = "bueren an der aare" if gdename == "baeren an der aare"
replace gdename = "celerina/schlarigna" if gdename == "celerinalschlarigna"
replace gdename = "biel be" if gdename == "bielibienne"
replace gdename = "delemont" if gdename == "del emont"
replace gdename = "les breuleux" if gdename == "breuleux"
replace gdename = "la chaux de fonds" if gdename == "la chaos de fonds"
replace gdename = "la chaux de fonds" if gdename == "la chaun de fonds"
replace gdename = "la chaux de fonds" if gdename == "la chauo de fonds"
replace gdename = "la chaux de fonds" if gdename == "la chaus de fonds"
replace gdename = "la chaux de fonds" if gdename == "la chaut de fonds"
replace gdename = "la chaux de fonds" if gdename == "la chaux de fond"
replace gdename = "la chaux de fonds" if gdename == "la chaux de fonds pr"
replace gdename = "la chaux de fonds" if gdename == "la chaux des fonds"
replace gdename = "la chaux de fonds" if gdename == "la chauxde fonds"
replace gdename = "le paquier ne" if gdename == "le paequier ne"
replace gdename = "le paquier fr" if gdename == "le paequier fr"
replace gdename = "le noirmont" if gdename == "le noir mont"
replace gdename = "geneve" if gdename == "s. a. geneve"
replace gdename = "lausanne" if gdename == "s. a. lausanne"
replace gdename = "geneve" if gdename == "sa geneve"
replace gdename = "lugano" if gdename == "salugano"
replace gdename = "tavetsch" if gdename == "salva tavetsch"
replace gdename = "salvan" if gdename == "salven"
replace gdename = "tujetsch" if gdename == "tavetsch"
replace gdename = "urtenen schoenbuehl" if gdename == "urtenen"
replace gdename = "urtenen" if gdename == "urtenen schoenbuehl" & year < 2003
replace gdename = "wohlen ag" if gdename == "wohien ag"
replace gdename = "abtwil" if gdename == "abtwil ag"
replace gdename = "le grand saconnex" if gdename == "le grandsaconnex"
replace gdename = "le grand saconnex" if gdename == "le grand saconnex pr"
replace gdename = "le grand saconnex" if gdename == "le grand saconnex s"
replace gdename = "le grand saconnex" if gdename == "le grand saconnex vp"
replace gdename = "le grand saconnex" if gdename == "le grand saxonnex"
replace gdename = "le grand saconnex" if gdename == "le grandsaconnex pr"
replace gdename = "vernier" if gdename == "vemier"
replace gdename = "le noirmont" if gdename == "noirmont"
replace gdename = "opfikon" if gdename == "opfikon pr"
replace gdename = "baar" if gdename == "saar"
replace gdename = "cham" if gdename == "chem"
replace gdename = "bad ragaz" if gdename == "bad ragez"
replace gdename = "marly" if gdename == "marty"
replace gdename = "vaz/obervaz" if gdename == "vazjobervaz"
replace gdename = "vaz/obervaz" if gdename == "vazlobervaz"
replace gdename = "vaz/obervaz" if gdename == "vazobervaz"
replace gdename = "vaz/obervaz" if gdename == "vaziobervaz"
replace gdename = "chermignon" if gdename == "crans s/chermignon"
replace gdename = "villars sur glane" if gdename == "villars sur giane"
replace gdename = "essert fr" if gdename == "essert"
replace gdename = "essertines sur rolle" if gdename == "essertines s/rolle"
replace gdename = "essertines sur rolle" if gdename == "essertines/rolle"
replace gdename = "saas" if gdename == "saas im praettigau"
replace gdename = "saas" if gdename == "saas im praettigau gr"
replace gdename = "saas" if gdename == "saas praettigau"
replace gdename = "saas balen" if gdename == "saas baten"
replace gdename = "schwarzenburg" if gdename == "schwarzenburg lu"
replace gdename = "oetwil an der limmat" if gdename == "limmat"
replace gdename = "santa maria im muenstertal" if gdename == "maria im muenstertal"
replace gdename = "koeniz" if gdename == "kniz"
replace gdename = "mesocco" if gdename == "di mesocco"
replace gdename = "ormont dessus" if gdename == "diablerets/ormont dessus"
replace gdename = "lens" if gdename == "crans s/lens"
replace gdename = "icogne" if gdename == "crans s/lcogne"
replace gdename = "icogne" if gdename == "crans s/icogne"
replace gdename = "lens" if gdename == "crans s lens"
replace gdename = "chermignon" if gdename == "crans s chermignon"
replace gdename = "sierre" if gdename == "crans s/ sierre"
replace gdename = "montana" if gdename == "crans s/montana"
replace gdename = "lens" if gdename == "crans sierre"
replace gdename = "lens" if gdename == "crans s/ sierre"
replace gdename = "chermignon" if gdename == "crans sur chermignon"
replace gdename = "basel" if gdename == "base"
replace gdename = "basel" if gdename == "base 1"
replace gdename = "basel" if gdename == "basef"
replace gdename = "wuennewil flamatt" if gdename == "flamatt wuennewil"
replace gdename = "wuennewil flamatt" if gdename == "flamatt/wuennewil"
replace gdename = "flawil" if gdename == "flawii"
replace gdename = "lonay" if gdename == "loney"
replace gdename = "st antoenien ascharina" if gdename == "antoenien ascharina"
replace gdename = "st antoenien castels" if gdename == "antoenien castels"
replace gdename = "sils im engadin/segl" if gdename == "siislsegl im engadin"
replace gdename = "signy" if gdename == "signy vr"
replace gdename = "monte carasso" if gdename == "montecarasso"
replace gdename = "estavayer le lac" if gdename == "esta vayer le lac"
replace gdename = "estavayer le lac" if gdename == "estavayer e lac"
replace gdename = "estavayer le lac" if gdename == "estavayer ie lac"
replace gdename = "estavayer le lac" if gdename == "estavayerae lac"
replace gdename = "estavayer le lac" if gdename == "estavayer4e lac"
replace gdename = "estavayer le lac" if gdename == "estavayerle lac"
replace gdename = "estavayer le lac" if gdename == "estavayer le lac vr"
replace gdename = "estavayer le lac" if gdename == "estavayer le lac s"
replace gdename = "romainmotier envy" if gdename == "romainmotier"
replace gdename = "romainmotier envy" if gdename == "romainmoetier"
replace gdename = "ebnat kappel" if gdename == "ebnat"
replace gdename = "ebnat kappel" if gdename == "ebnat koppel"
replace gdename = "ebnat kappel" if gdename == "ebnat kappes"
replace gdename = "ebnat" if gdename == "ebnat kappel" & year < 1965
replace gdename = "saint barthelemy vd" if gdename == "st barthelemy"
replace gdename = "saint barthelemy vd" if gdename == "st barthelemy vd"
replace gdename = "saint blaise" if gdename == "st biaise"
replace gdename = "saint brais" if gdename == "st brais"
replace gdename = "grono" if gdename == "grano"
replace gdename = "inkwil" if gdename == "lnkwil"
replace gdename = "altstaetten" if gdename == "altstaellen sg"
replace gdename = "altstaetten" if gdename == "altstaetten so"
replace gdename = "altstaetten" if gdename == "altstaetten/sg"
replace gdename = "altstaetten" if gdename == "altstatten"
replace gdename = "altstaetten" if gdename == "altstaettensg"
replace gdename = "altstaetten" if gdename == "altstatten sg"
replace gdename = "altstaetten" if gdename == "altsttten"
replace gdename = "val d illiez" if gdename == "vai d illiez"
replace gdename = "vallorbe" if gdename == "vailorbe"
replace gdename = "val d illiez" if gdename == "val d llliez"
replace gdename = "valeyres sous ursins" if gdename == "valeryres sous ursins"
replace gdename = "vallorbe" if gdename == "vallorbe vp"
replace gdename = "crans pres celigny" if gdename == "crans pres celigny pr"
replace gdename = "mendrisio" if gdename == "mendriso"
replace gdename = "rueti zh" if gdename == "ruti zh"
replace gdename = "lausanne" if gdename == "ausanne"
replace gdename = "chermignon" if gdename == "corn chermignon"
replace gdename = "biel be" if gdename == "blel/bienne"
replace gdename = "zuerich" if gdename == "zuerich[d]"
replace gdename = "chiasso" if gdename == "chiasso s"
replace gdename = "basel" if gdename == "basel[p]"
replace gdename = "chur" if gdename == "coire"
replace gdename = "baar" if gdename == "baar vr"
replace gdename = "chur" if gdename == "chur vr"
replace gdename = "biel be" if gdename == "bieli bienne"
replace gdename = "la tour de treme" if gdename == "la tour de trime"
replace gdename = "yverdon les bains" if gdename == "yverdon les bain"
replace gdename = "biel be" if gdename == "biel/blenne"
replace gdename = "lens" if gdename == "crans sur sierre/lens"
replace gdename = "fribourg" if gdename == "fribourg s"
replace gdename = "la chaux de fonds" if gdename == "lachaux de fonds"
replace gdename = "hofen" if gdename == "hofen sh"
replace gdename = "kilchberg zh" if gdename == "kirchberg zh"
replace gdename = "vernier" if gdename == "vernier pr"
replace gdename = "pfaeffikon" if gdename == "pfaffikon zh"
replace gdename = "le grand saconnex" if gdename == "gr saconnex"
replace gdename = "yverdon les bains" if gdename == "yverdon4es bains"
replace gdename = "baar" if gdename == "baar pr"
replace gdename = "la tour de treme" if gdename == "tour de treme"
replace gdename = "saint sulpice vd" if gdename == "stsulpicevd"
replace gdename = "lussy sur morges" if gdename == "lussy s/morges"
replace gdename = "grone" if gdename == "groene"
replace gdename = "gaiserwald" if gdename == "engelburg gaiserwald"
replace gdename = "neuchatel" if gdename == "neuchatei"
replace gdename = "lugano" if gdename == "lu gano"
replace gdename = "vandoeuvres" if gdename == "vandeeuvres"
replace gdename = "sierre" if gdename == "sierreisiders"
replace gdename = "nyon" if gdename == "nyon pr"
replace gdename = "st margrethen" if gdename == "margrethen"
replace gdename = "wil" if gdename == "wil/sg"
replace gdename = "paradiso" if gdename == "paradise"
replace gdename = "adliswil" if gdename == "adliswil pr"
replace gdename = "bellinzona" if gdename == "bellinzona vr"
replace gdename = "wohlen bei bern" if gdename == "wohlen bei bem"
replace gdename = "ziefen" if gdename == "zielen"
replace gdename = "safnern" if gdename == "safnem"
replace gdename = "les brenets" if gdename == "brenets"
replace gdename = "dornach" if gdename == "domach"
replace gdename = "lens" if gdename == "crans sur lens"
replace gdename = "saint legier la chiesaz" if gdename == "legier la chiesaz"
replace gdename = "chateau d oex" if gdename == "chteau d oex"
replace gdename = "flums" if gdename == "plums"
replace gdename = "littau" if gdename == "reussbhl"
replace gdename = "littau" if gdename == "reussbuehl / littau"
replace gdename = "littau" if gdename == "reussbuehl littau"
replace gdename = "littau" if gdename == "reussbuhl"
replace gdename = "st margrethen" if gdename == "margarethen"
replace gdename = "st margrethen" if gdename == "margarethen/sg"
replace gdename = "martigny" if gdename == "matigny"
replace gdename = "mattstetten" if gdename == "mattsteffen"
replace gdename = "chene bougeries" if gdename == "ch bougeries"
replace gdename = "arth" if gdename == "gde arth"
replace gdename = "bern" if gdename == "gde bern"
replace gdename = "bolligen" if gdename == "gde bolligen"
replace gdename = "freienbach" if gdename == "gde freienbach"
replace gdename = "uzwil" if gdename == "gde uzwil"
replace gdename = "geneve" if gdename == "gebeve"
replace gdename = "geneve" if gdename == "geheve"
replace gdename = "henau" if gdename == "uzwil"
replace gdename = "uzwil" if gdename == "uzswil"
replace gdename = "uzwil" if gdename == "henau" & year > 1962
replace gdename = "saint aubin sauges" if gdename == "saint aubin ne"
replace gdename = "geneve" if gdename == "genlve"
replace gdename = "geneve" if gdename == "geni ve"
replace gdename = "geneve" if gdename == "genevo"
replace gdename = "geneve" if gdename == "genevez"
replace gdename = "gansingen" if gdename == "gensingen"
replace gdename = "hasle lu" if gdename == "haste lu"
replace gdename = "lausanne" if gdename == "alausanne"
replace gdename = "haslen" if gdename == "hasten"
replace gdename = "tschierv" if gdename == "tschiery"
replace gdename = "ecublens vd" if gdename == "ecublensvd"
replace gdename = "schlatt haslen" if gdename == "schlatt bei appenzell"
replace gdename = "itingen" if gdename == "ltingen"
replace gdename = "biel be" if gdename == "biel: bienne"
replace gdename = "inwil" if gdename == "lnwil"
replace gdename = "interlaken" if gdename == "inter laken"
replace gdename = "mels" if gdename == "meis"
replace gdename = "malters" if gdename == "mailers"
replace gdename = "sion" if gdename == "bramois"
replace gdename = "zug" if gdename == "zua"
replace gdename = "ilanz" if gdename == "llanz"
replace gdename = "buchberg" if gdename == "buchberg sh"
replace gdename = "praroman" if gdename == "praroman le mouret"
replace gdename = "pratteln" if gdename == "prat teln"
replace gdename = "igis" if gdename == "landquart igis"
replace gdename = "igis" if gdename == "landquart lgis"
replace gdename = "fischbach" if gdename == "fischbach lu"
replace gdename = "au sg" if gdename == "heerbrugg au"
replace gdename = "balgach" if gdename == "heerbrugg balgach"
replace gdename = "berneck" if gdename == "heerbrugg berneck"
replace gdename = "au sg" if gdename == "heerbrugg/au"
replace gdename = "au sg" if gdename == "heerbrugg/au sg"
replace gdename = "littau" if gdename == "uttau"
replace gdename = "lueterkofen ichertswil" if gdename == "lueterkofen"
replace gdename = "lueterkofen ichertswil" if gdename == "lueterkofen lchertswil"
replace gdename = "lugano" if gdename == "lugaeno"
replace gdename = "lugano" if gdename == "lugaho"
replace gdename = "lugano" if gdename == "lugamo"
replace gdename = "lugano" if gdename == "lugana"
replace gdename = "lugano" if gdename == "luganc"
replace gdename = "lugano" if gdename == "lugand"
replace gdename = "lugano" if gdename == "lugani"
replace gdename = "lugano" if gdename == "lugano[p]"
replace gdename = "lugano" if gdename == "lugaro"
replace gdename = "randogne" if gdename == "vermala sur randogne"
replace gdename = "tujetsch" if gdename == "sedrun tujetsch"
replace gdename = "tavetsch" if gdename == "tujetsch" & year < 1976
replace gdename = "lugano" if gdename == "tugano"
replace gdename = "zuerich" if gdename == "tuerich"
replace gdename = "renens vd" if gdename == "renensvd"
replace gdename = "biel be" if gdename == "biei/bienne"
replace gdename = "farvagny" if gdename == "farvagny le grand"
replace gdename = "farvagny le grand" if gdename == "favargny le grand"
replace gdename = "zell zh" if gdename == "zeii zh"
replace gdename = "buelach" if gdename == "búlach"
replace gdename = "villa gr" if gdename == "villa"
replace gdename = "vella" if gdename == "villa gr" & year > 1986
replace gdename = "ollon" if gdename == "oiion vd"
replace gdename = "ollon" if gdename == "olion"
replace gdename = "biel be" if gdename == "biel/bienne vr"
replace gdename = "hauenstein ifenthal" if gdename == "hauenstein"
replace gdename = "zillis reischen" if gdename == "zillis"
replace gdename = "zofingen" if gdename == "zo fingen"
replace gdename = "zollikon" if gdename == "zoilikon"
replace gdename = "zollikon" if gdename == "zol likon"
replace gdename = "zollikon" if gdename == "zoliikon"
replace gdename = "zollikon" if gdename == "zolli kon"
replace gdename = "zollikon" if gdename == "zollikon pr"
replace gdename = "zollikon" if gdename == "zollikon zh"
replace gdename = "zollikofen" if gdename == "zollikoten"
replace gdename = "zollikon" if gdename == "zolllkon"
replace gdename = "zollikofen" if gdename == "zolllkofen"
replace gdename = "zollikon" if gdename == "zoltikon"
replace gdename = "zollikofen" if gdename == "zotlikofen"
replace gdename = "zollikon" if gdename == "zotlikon"
replace gdename = "zollikofen" if gdename == "zollikofen station"
replace gdename = "zollikofen" if gdename == "zollikof en"
replace gdename = "zollikofen" if gdename == "zolli kofen"
replace gdename = "zollikofen" if gdename == "zoliikofen"
replace gdename = "zollikofen" if gdename == "zoilikofen"
replace gdename = "zofingen" if gdename == "zotingen"
replace gdename = "erlen" if gdename == "riedt bei erlen"
replace gdename = "faellanden" if gdename == "fuellanden"
replace gdename = "urnaesch" if gdename == "umaesch"
replace gdename = "murten" if gdename == "murten moral"
replace gdename = "stetten sh" if gdename == "steffen sh"
replace gdename = "wil sg" if gdename == "wiisg"
replace gdename = "chavannes pres renens" if gdename == "chavannes par renens"
replace gdename = "chavannes pres renens" if gdename == "chavannes pres renons"
replace gdename = "chavannes pres renens" if gdename == "chavannes prs renens"
replace gdename = "chavannes pres renens" if gdename == "chavannes/renens"
replace gdename = "chavannes pres renens" if gdename == "chavannes pres repens"
replace gdename = "chavannes de bogis" if gdename == "chavennes de bogis"
replace gdename = "orsieres" if gdename == "ornieres"
replace gdename = "ayent" if gdename == "anzere ayent"
replace gdename = "oron le chatel" if gdename == "oron le chaetel"
replace gdename = "kirchberg sg" if gdename == "bazenheid kirchberg"
replace gdename = "kirchberg sg" if gdename == "bazenheid kirchberg sg"
replace gdename = "bellerive vd" if gdename == "bellerive"
replace gdename = "bellerive vd" if gdename == "bellerive (vaud)"
replace gdename = "bellerive vd" if gdename == "bellerive[d]"
replace gdename = "bellerive vd" if gdename == "bellerive[p]"
replace gdename = "bellerive vd" if gdename == "bellerivevd"
replace gdename = "bellinzona" if gdename == "bellin zona"
replace gdename = "bellinzona" if gdename == "bellinzona pr"
replace gdename = "bellinzona" if gdename == "bellinzona ti"
replace gdename = "bellinzona" if gdename == "bellinzonna"
replace gdename = "bellinzona" if gdename == "bellizona"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolfstetten friedlisbg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolfstetten frieslisberg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolfstten friedlisberg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolistetten friedlisberg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudoltstetten friedlisberg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudotfstetten friedlisberg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolfstetten friediisberg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudolf stetten friedlisberg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudoifstetten friedlisberg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudofstetten friedlisberg"
replace gdename = "rudolfstetten friedlisberg" if gdename == "rudoffstetten friedlisberg"
replace gdename = "zug" if gdename == "zue"
replace gdename = "zug" if gdename == "zig"
replace gdename = "luzern" if gdename == "zern"
replace gdename = "zermatt" if gdename == "zer matt"
replace gdename = "zollikon" if gdename == "zallikon"
replace gdename = "zuerich" if gdename == "z!rich"
replace gdename = "zuerich" if gdename == "z rich"
replace gdename = "yvorne" if gdename == "yvome"
replace gdename = "yvorne" if gdename == "yvonne"
replace gdename = "yvorne" if gdename == "yvorn e"
replace gdename = "yvonand" if gdename == "yvonand vr"
replace gdename = "yvonand" if gdename == "yvonand vd"
replace gdename = "yvonand" if gdename == "yvonnand"
replace gdename = "zell lu" if gdename == "zeii lu"
replace gdename = "wartau" if gdename == "truebbach / wartau"
replace gdename = "wartau" if gdename == "truebbach gmd. wartau"
replace gdename = "wartau" if gdename == "truebbach sg"
replace gdename = "wartau" if gdename == "truebbach wartau"
replace gdename = "feusisberg" if gdename == "schindeleggi feusisberg"
replace gdename = "feusisberg" if gdename == "schindellegi (feusisberg)"
replace gdename = "feusisberg" if gdename == "schindellegi feusisberg"
replace gdename = "feusisberg" if gdename == "schindeliegi"
replace gdename = "schlatt tg" if gdename == "schlaff tg"
replace gdename = "schlatt tg" if gdename == "schlaft tg"
replace gdename = "adlikon" if gdename == "adlikon bei andelfingen"
replace gdename = "biel be" if gdename == "biel/ bienne"
replace gdename = "oberdorf bl" if gdename == "bl oberwil"
replace gdename = "reinach bl" if gdename == "bl reinach"
replace gdename = "st blaise" if gdename == "blaise"

sort gdename year
drop gdenr gdenr_2018 ctn E_CNTR N_CNTR geo_merge
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2
gen geo_merge = 4 if _merge == 3

preserve
keep if _merge == 3
save "$path\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_4.dta", replace
restore

* 2) d) String merge on further standardized company-gdename info (step 4)

keep if _merge == 1
drop _merge

replace gdename = "roveredo gr" if gdename == "roveredo"
replace gdename = "zuerich" if gdename == "zuerich pr"
replace gdename = "hofstetten flueh" if gdename == "hofstetten so"
replace gdename = "hofstetten flueh" if gdename == "hofstetten fluh so"
replace gdename = "hausen bei brugg" if gdename == "hausen (ag)"
replace gdename = "schwarzenbach" if gdename == "schwarzenberg lu"
replace gdename = "cressier fr" if gdename == "cressier sur morat"
replace gdename = "lully vd" if gdename == "lully sur morges"
replace gdename = "villars sur glane" if gdename == "villars sur glne"
replace gdename = "malix" if gdename == "malfix"
replace gdename = "illnau" if gdename == "illnau effretikon"
replace gdename = "illnau effretikon" if gdename == "illnau effretikon vr"
replace gdename = "illnau" if gdename == "illnau effretikon vr" & year < 1974
replace gdename = "illnau effretikon" if gdename == "illnau effretlkon"
replace gdename = "illnau" if gdename == "illnau effretlkon" & year < 1974
replace gdename = "saint imier" if gdename == "imier"
replace gdename = "illnau effretikon" if gdename == "illnau eff"
replace gdename = "torricella taverne" if gdename == "torricella"
replace gdename = "torricella taverne" if gdename == "torricella taveme"
replace gdename = "ollon" if gdename == "oiion"
replace gdename = "celerina/schlarigna" if gdename == "celerina/ schlarigna"
replace gdename = "jona" if gdename == "jana"
replace gdename = "celerina/schlarigna" if gdename == "celerina schlarigna"
replace gdename = "sennwald" if gdename == "haag"
replace gdename = "sennwald" if gdename == "haag sennwald"
replace gdename = "sennwald" if gdename == "haag/sennwald"
replace gdename = "le landeron" if gdename == "landeron combes"
replace gdename = "langenthal" if gdename == "langen thal"
replace gdename = "haegendorf" if gdename == "hagendorf"
replace gdename = "cham" if gdename == "hagendorn"
replace gdename = "cham" if gdename == "hagendorn cham"
replace gdename = "egnach" if gdename == "steinebrunn (egnach)"
replace gdename = "egnach" if gdename == "steinebrunn egnach"
replace gdename = "wangen bei olten" if gdename == "wangen olten"
replace gdename = "wangen bei olten" if gdename == "wangen/olten"
replace gdename = "wangen bruettisellen" if gdename == "wangen bnittisellen"
replace gdename = "wangen bruettisellen" if gdename == "wangen bnuettisellen"
replace gdename = "wangen bruettisellen" if gdename == "wangen bruettiseilen"
replace gdename = "wangen bruettisellen" if gdename == "wangen bruettisel len"
replace gdename = "wangen bruettisellen" if gdename == "wangen bruettiselien"
replace gdename = "wangen bruettisellen" if gdename == "wangen bruettisellen vr"
replace gdename = "wangen bruettisellen" if gdename == "wangen bruttisellen"
replace gdename = "wangen bei olten" if gdename == "wangen b olten"
replace gdename = "wangen bei olten" if gdename == "wangen b/olten"
replace gdename = "roggwil tg" if gdename == "freidorf roggwil"
replace gdename = "oberentfelden" if gdename == "oberent felden"
replace gdename = "oberentfelden" if gdename == "oberentfeden"
replace gdename = "oberentfelden" if gdename == "oberentfeiden"
replace gdename = "oberentfelden" if gdename == "oberenttelden"
replace gdename = "obergesteln" if gdename == "obergestein"
replace gdename = "sumiswald" if gdename == "wasen"
replace gdename = "sumiswald" if gdename == "wasen im emmental"
replace gdename = "sumiswald" if gdename == "wasen sumiswald"
replace gdename = "sumiswald" if gdename == "wasen/sumiswald"
replace gdename = "bern" if gdename == "bein"
replace gdename = "ollon" if gdename == "villars sur 0llon"
replace gdename = "buch sh" if gdename == "buch"
replace gdename = "barbereche" if gdename == "pensier"
replace gdename = "bernex" if gdename == "bemex"
replace gdename = "chavornay" if gdename == "chavomay"
replace gdename = "uetikon am see" if gdename == "uetikon"
replace gdename = "pfaeffikon" if gdename == "pfffikon"
replace gdename = "paradiso" if gdename == "paradiso pr"
replace gdename = "liestal" if gdename == "uestal"
replace gdename = "langnau am albis" if gdename == "langnau a. albis"
replace gdename = "bernex" if gdename == "bernez"
replace gdename = "mannens grandsivaz" if gdename == "mannens"
replace gdename = "murten" if gdename == "marten morat"
replace gdename = "murten" if gdename == "marten moral"
replace gdename = "romont fr" if gdename == "ramant fr"
replace gdename = "rafz" if gdename == "rat"
replace gdename = "rafz" if gdename == "raft"
replace gdename = "oberriet sg" if gdename == "oberniet sg"
replace gdename = "morlon" if gdename == "morton"
replace gdename = "roveredo gr" if gdename == "roveredogr"
replace gdename = "bernex" if gdename == "bernei"
replace gdename = "ingenbohl" if gdename == "brunnen ingenbohl"
replace gdename = "sion" if gdename == "sian"
replace gdename = "mendrisio" if gdename == "mendrislo"
replace gdename = "thunstetten" if gdename == "buetzberg thunstetten"
replace gdename = "ollon" if gdename == "0llon"
replace gdename = "olten" if gdename == "0lten"
replace gdename = "ollon" if gdename == "011on"
replace gdename = "olten" if gdename == "01 ten"
replace gdename = "bern" if gdename == "bem"
replace gdename = "jona" if gdename == "jena"
replace gdename = "jona" if gdename == "jona sg"
replace gdename = "brig glis" if gdename == "brig" & year > 1971
replace gdename = "rickenbach zh" if gdename == "rickenbach bei winterthur"
replace gdename = "rickenbach bei wil" if gdename == "rickenbach b/wil"
replace gdename = "ollon" if gdename == "villars ollon"
replace gdename = "langnau im emmental" if gdename == "langnau i/e"
replace gdename = "langnau im emmental" if gdename == "langneu im emmental"
replace gdename = "langnau im emmental" if gdename == "langnau i e"
replace gdename = "holderbank ag" if gdename == "holderbankag"
replace gdename = "sierre" if gdename == "sierre siders"
replace gdename = "genolier" if gdename == "genoller"
replace gdename = "buerglen ur" if gdename == "buergten ur"

sort gdename year
drop gdenr gdenr_2018 ctn E_CNTR N_CNTR geo_merge
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2
gen geo_merge = 5 if _merge == 3

preserve
keep if _merge == 3
save "$path\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_5.dta", replace
restore

* 3) Extraction of information from gdename var (succursalle, branch, etc.)

keep if _merge == 1
drop _merge

gen gdename_extract = "branch" if (regexm(gdename), "succursale") // create a new var indicating if succursale in the gdename
replace gdename_extract = "branch" if (regexm(gdename), "branch") // create a new var indicating if succursale in the gdename
replace gdename_extract = "branch" if (regexm(gdename), "zweigniederlassung") // create a new var indicating if succursale in the gdename
replace gdename = ustrregexra(gdename,"succursale de","")  // delete "succursale de" when in gdename
replace gdename = ustrregexra(gdename,"succursale di","")  // delete "succursale di" when in gdename
replace gdename = ustrregexra(gdename,"succursale d","")  // delete "succursale di" when in gdename
replace gdename = ustrregexra(gdename,"branch of","")  // delete "branch of" when in gdename
replace gdename = ustrregexra(gdename," branch","")  // delete "branch" when in gdename
replace gdename = ustrregexra(gdename,"zweigniederlassung","")  // delete "zweigniederlassung" when in gdename
replace gdename = ustrlower(gdename)
replace gdename = ustrtrim(gdename)
replace gdename = "geneve" if gdename == "geneva"
replace gdename = "geneve" if gdename == "geneva fice"
replace gdename = "zuerich" if gdename == "zurich"
replace gdename = "zug" if gdename == "fice zug"

sort gdename year
drop gdenr gdenr_2018 ctn E_CNTR N_CNTR geo_merge
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2
gen geo_merge = 6 if _merge == 3

preserve
keep if _merge == 3
save "$path\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_6.dta", replace
restore

* 4) Extraction of information from capital var (to find where the company is located)
keep if _merge == 1
drop _merge

* 4.1) Gdename with zip code & real mun name in capital var

replace gdename_extract = regexs(1) if (regexm(gdename, "([0-9][0-9][0-9][0-9])")) // finds 4-digit numbers within string
replace gdename = ustrregexra(gdename,"[0-9][0-9][0-9][0-9]", capital)  
replace gdename = ustrlower(gdename)
replace gdename = ustrtrim(gdename)

* 4.2) "sa" in gdename & real mun name in capital var
replace gdename_extract = "sa" if gdename == "sa" // gdename is "sa", but it should be part of cname. The gdename ended up in capital var.
replace gdename = ustrregexra(gdename,"sa", capital) if gdename_extract == "sa" // Replace the "sa" part with the capital var so that we find the correct gdename 

* 4.3) "le" in gdename & real rest of mun name in capital var
gen capital_extract = ""
replace capital_extract = capital if gdename == "le"
replace gdename = "le " + capital_extract if gdename == "le"

foreach str in gdename /*cname*/ {
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
	replace `str' = usubinstr(`str', "-", " ",.)
	replace `str' = usubinstr(`str', "~", " ",.)	
	replace `str' = usubinstr(`str', ":", " ",.)
	replace `str' = usubinstr(`str', "    ", " ",.)
	replace `str' = usubinstr(`str', "   ", " ",.)
	replace `str' = usubinstr(`str', "  ", " ",.)
	replace `str' = ustrtrim(`str')
	replace `str' = ustrlower(`str')
	replace `str' = ustrtrim(`str')
}

sort gdename year
drop gdenr gdenr_2018 ctn E_CNTR N_CNTR geo_merge
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2
gen geo_merge = 7 if _merge == 3

preserve
keep if _merge == 3
save "$path\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_7.dta", replace
restore

* 5) Extraction of information from cname var (cases where municipality ended up in company name)

keep if _merge == 1
drop _merge

* 5.1) Extracting "Hergiswil am see" to distinguish if Hergiswil LU or NW
gen cname_extract = "hergiswil am see" if (regexm(cname, "hergiswil am see")) // create new var indicating if "Hergiswil am see" within cname
replace gdename = "hergiswil nw" if cname_extract == "hergiswil am see"

* 5.2) Extracting "buchs sg/lu/zh" to distinguish "buchs" if AG, LU, SG or ZH 
replace cname_extract = "buchs sg" if (regexm(cname, "buchs sg")) // create new var indicating if "buchs sg" within cname
replace gdename = "buchs sg" if cname_extract == "buchs sg"
replace cname_extract = "buchs lu" if (regexm(cname, "buchs lu")) // create new var indicating if "buchs lu" within cname
replace gdename = "buchs lu" if cname_extract == "buchs lu"
replace cname_extract = "buchs zh" if (regexm(cname, "buchs zh")) // create new var indicating if "buchs zh" within cname
replace gdename = "buchs zh" if cname_extract == "buchs zh"

* 5.3) Extracting "wil sg/zh" to distinguish if or 
replace cname_extract = "wil sg" if (regexm(cname, "wil sg")) 
replace gdename = "wil sg" if cname_extract == "wil sg"
replace cname_extract = "wil zh" if (regexm(cname, "wil zh")) 
replace gdename = "wil zh" if cname_extract == "wil zh"

* 5.4) Extracting "wohlen be" to distinguish if wohlen be or wohlen ag 
replace cname_extract = "wohlen be" if (regexm(cname, "wohlen be")) 
replace gdename = "wohlen be" if cname_extract == "wohlen be"

* 5.5) Extracting "gossau sg/zh" to distinguish if gossau sg or gossau zh
replace cname_extract = "gossau sg" if (regexm(cname, "gossau sg")) 
replace gdename = "gossau sg" if cname_extract == "gossau sg"
replace cname_extract = "gossau zh" if (regexm(cname, "gossau zh")) 
replace gdename = "gossau zh" if cname_extract == "gossau zh"

* 5.6) Extracting "rapperswil sg/be" to distinguish if be or sg
replace cname_extract = "rapperswil sg" if (regexm(cname, "rapperswil sg")) 
replace gdename = "rapperswil sg" if cname_extract == "rapperswil sg"
replace cname_extract = "rapperswil be" if (regexm(cname, "rapperswil be")) 
replace gdename = "rapperswil be" if cname_extract == "rapperswil be"

* 5.7) Extracting "kilchberg bl/zh" 
replace cname_extract = "kilchberg bl" if (regexm(cname, "kilchberg bl")) 
replace gdename = "kilchberg bl" if cname_extract == "kilchberg bl"
replace cname_extract = "kilchberg zh" if (regexm(cname, "kilchberg zh")) 
replace gdename = "kilchberg zh" if cname_extract == "kilchberg zh"

* 5.8 ) Extracting "altdorf sh/ur"
replace cname_extract = "altdorf sh" if (regexm(cname, "altdorf sh")) 
replace gdename = "altdorf sh" if cname_extract == "altdorf sh"
replace cname_extract = "altdorf ur" if (regexm(cname, "altdorf ur")) 
replace gdename = "altdorf ur" if cname_extract == "altdorf ur"

* 5.9) Extracting "romont be/fr" 
replace cname_extract = "romont be" if (regexm(cname, "romont be")) 
replace gdename = "romont be" if cname_extract == "romont be"
replace cname_extract = "romont fr" if (regexm(cname, "romont fr")) 
replace gdename = "romont fr" if cname_extract == "romont fr"

* 5.10) Extracting "muri ag / bei bern"  
replace cname_extract = "muri b" if (regexm(cname, "muri b")) 
replace gdename = "muri bei bern" if cname_extract == "muri b"

* 5.11) Extracting "kirchberg be/sg" 
replace cname_extract = "kirchberg be" if (regexm(cname, "kirchberg be")) 
replace gdename = "kirchberg be" if cname_extract == "kirchberg be"
replace cname_extract = "kirchberg sg" if (regexm(cname, "kirchberg sg")) 
replace gdename = "kirchberg sg" if cname_extract == "kirchberg sg"

* 5.12) Extracting "wald ar/zh"  
replace cname_extract = "wald ar" if (regexm(cname, "wald ar")) 
replace gdename = "wald ar" if cname_extract == "wald ar"
replace cname_extract = "wald zh" if (regexm(cname, "wald zh")) 
replace gdename = "wald zh" if cname_extract == "wald zh"

* 5.13) Extracting "aesch bl/lu/zh"  
replace cname_extract = "aesch bl" if (regexm(cname, "aesch bl")) 
replace gdename = "aesch bl" if cname_extract == "aesch bl"
replace cname_extract = "aesch lu" if (regexm(cname, "aesch lu")) 
replace gdename = "aesch lu" if cname_extract == "aesch lu"
replace cname_extract = "aesch zh" if (regexm(cname, "aesch zh")) 
replace gdename = "aesch zh" if cname_extract == "aesch zh"

* 5.14) Extracting "rueti gl/zh" 
replace cname_extract = "rueti gl" if (regexm(cname, "rüti gl")) 
replace gdename = "rueti gl" if cname_extract == "rueti gl"
replace cname_extract = "rueti zh" if (regexm(cname, "rüti zh")) 
replace gdename = "rueti zh" if cname_extract == "rueti zh"

* 5.15) Extracting "villeneuve fr/vd"  
replace cname_extract = "villeneuve fr" if (regexm(cname, "villeneuve fr")) 
replace gdename = "villeneuve fr" if cname_extract == "villeneuve fr"
replace cname_extract = "villeneuve vd" if (regexm(cname, "villeneuve vd")) 
replace gdename = "villeneuve vd" if cname_extract == "villeneuve vd"
replace cname_extract = "villeneuve vd" if (regexm(cname, "villeneuve vd")) 
replace gdename = "villeneuve vd" if cname_extract == "villeneuve vd"

* 5.16) Extracting "ecublens fr/vd"  
replace cname_extract = "ecublens fr" if (regexm(cname, "ecublens fr")) 
replace gdename = "ecublens fr" if cname_extract == "ecublens fr"
replace cname_extract = "ecublens vd" if (regexm(cname, "ecublens vd")) 
replace gdename = "ecublens vd" if cname_extract == "ecublens vd"

* 5.17) Extracting "oberwil ag/bei bueren/bl/tg/im simmental/lieli"  
replace cname_extract = "oberwil tg" if (regexm(cname, "oberwil tg")) 
replace gdename = "oberwil tg" if cname_extract == "oberwil tg"
replace cname_extract = "oberwil bl" if (regexm(cname, "oberwil bl")) 
replace gdename = "oberwil bl" if cname_extract == "oberwil bl"

* 5.18) Extracting "zell lu/zh"  
replace cname_extract = "zell lu" if (regexm(cname, "zell lu")) 
replace gdename = "zell lu" if cname_extract == "zell lu"
replace cname_extract = "zell lu" if (regexm(cname, "hüswil")) 
replace gdename = "zell lu" if cname_extract == "zell lu"
replace cname_extract = "zell zh" if (regexm(cname, "zell zh")) 
replace gdename = "zell zh" if cname_extract == "zell zh"

* 5.19) Extracting "wangen sz/zh/so/be", "wangen bruettisellen", "wangen bei olten", "wangen an der aare"  
replace cname_extract = "wangen sz" if (regexm(cname, "wangen sz")) 
replace gdename = "wangen sz" if cname_extract == "wangen sz"

replace cname_extract = "wangen sz" if (regexm(cname, "siebnen")) & gdename == "wangen"
replace gdename = "wangen sz" if cname_extract == "wangen sz"

replace cname_extract = "walliswil bei wangen" if (regexm(cname, "walliswil")) & gdename == "wangen"
replace gdename = "walliswil bei wangen" if cname_extract == "walliswil bei wangen"

replace cname_extract = "wangen zh" if (regexm(cname, "wangen zh")) 
replace gdename = "wangen zh" if cname_extract == "wangen zh"

replace cname_extract = "wangen zh" if (regexm(cname, "wangen brüttisellen")) 
replace gdename = "wangen zh" if cname_extract == "wangen zh"

replace cname_extract = "wangen bei olten" if (regexm(cname, "wangen bei olten")) 
replace gdename = "wangen bei olten" if cname_extract == "wangen bei olten"

replace cname_extract = "wangen bei olten" if (regexm(cname, "wangen so")) 
replace gdename = "wangen bei olten" if cname_extract == "wangen bei olten"

replace cname_extract = "wangen an der aare" if (regexm(cname, "wangen a")) & gdename == "a"
replace gdename = "wangen an der aare" if cname_extract == "wangen an der aare" & gdename == "a"
replace cname_extract = "wangen an der aare" if (regexm(cname, "wangen a")) & gdename == "aare"
replace gdename = "wangen an der aare" if cname_extract == "wangen an der aare" & gdename == "aare"
replace gdename = "bueren an der aare" if gdename == "aare" // previous step left only bueren an der aare as option when gdename is "aare" (double checked with cname)


* 5.20) Extracting from gdename = "a" or "albis"
replace cname_extract = "affoltern am albis" if (regexm(cname, "affoltern a")) & gdename == "a"
replace gdename = "affoltern am albis" if cname_extract == "affoltern am albis" & gdename == "a"
replace cname_extract = "affoltern am albis" if (regexm(cname, "affoltern a")) & gdename == "albis"
replace gdename = "affoltern am albis" if cname_extract == "affoltern am albis" & gdename == "albis"

replace cname_extract = "langnau am albis" if (regexm(cname, "langnau a")) & gdename == "a"
replace gdename = "langnau am albis" if cname_extract == "langnau am albis" & gdename == "a"
replace cname_extract = "langnau am albis" if (regexm(cname, "langnau a")) & gdename == "albis"
replace gdename = "langnau am albis" if cname_extract == "langnau am albis" & gdename == "albis"

replace cname_extract = "wettswil am albis" if (regexm(cname, "wettswil a")) & gdename == "a" 
replace gdename = "wettswil" if cname_extract == "wettswil am albis" & gdename == "a"
replace cname_extract = "wettswil am albis" if (regexm(cname, "wettswil a")) & gdename == "albis" 
replace gdename = "wettswil" if cname_extract == "wettswil am albis" & gdename == "albis"
replace gdename = "wettswil am albis" if gdename == "wettswil" & year > 1975

replace cname_extract = "aeugst am albis" if (regexm(cname, "aeugst a")) & gdename == "a" 
replace gdename = "aeugst" if cname_extract == "aeugst am albis" & gdename == "a"
replace cname_extract = "aeugst am albis" if (regexm(cname, "aeugst a")) & gdename == "albis" 
replace gdename = "aeugst" if cname_extract == "aeugst am albis" & gdename == "albis"
replace gdename = "aeugst am albis" if gdename == "aeugst" & year > 1975

replace cname_extract = "hausen am albis" if (regexm(cname, "hausen a")) & gdename == "a" 
replace gdename = "hausen am albis" if cname_extract == "hausen am albis" & gdename == "a"
replace cname_extract = "hausen am albis" if (regexm(cname, "hausen a")) & gdename == "albis" 
replace gdename = "hausen am albis" if cname_extract == "hausen am albis" & gdename == "albis"

replace cname_extract = "kappel am albis" if (regexm(cname, "kappel a")) & gdename == "a" 
replace gdename = "kappel am albis" if cname_extract == "kappel am albis" & gdename == "a"

replace cname_extract = "uitikon" if (regexm(cname, "uitikon a")) & gdename == "a" // also called uitikon am albis
replace gdename = "uitikon" if cname_extract == "uitikon" & gdename == "a"
replace cname_extract = "uitikon" if (regexm(cname, "uitikon a")) & gdename == "albis" // also called uitikon am albis
replace gdename = "uitikon" if cname_extract == "uitikon" & gdename == "albis"

* 5.21) Extracting from gdename = "au" and heerbrugg is specified in cname (we know its au sg)
replace cname_extract = "au sg" if (regexm(cname, "heerbrugg")) & gdename == "au"
replace gdename = "au sg" if cname_extract == "au sg" & gdename == "au"
replace cname_extract = "au sg" if (regexm(cname, "au sg")) & gdename == "au"
replace gdename = "au sg" if cname_extract == "au sg" & gdename == "au"

* 5.22) Extracting from gdname ="biel" and bienne or lengnau is specified in cname (it is biel be or lengnau be)
replace cname_extract = "biel be" if (regexm(cname, "bienne")) & gdename == "biel"
replace gdename = "biel be" if cname_extract == "biel be" & gdename == "biel"
replace cname_extract = "lengnau be" if (regexm(cname, "lengnau")) & gdename == "biel"
replace gdename = "lengnau be" if cname_extract == "lengnau be" & gdename == "biel"

* 5.23) Extracting from gdname = "e": it is (very) often "im emmental" cut in two (one part in cname and "e" in gdename)
replace cname_extract = "langnau im emmental" if (regexm(cname, "langnau i")) & gdename == "e"
replace gdename = "langnau im emmental" if cname_extract == "langnau im emmental" & gdename == "e"

replace cname_extract = "wasen im emmental" if (regexm(cname, "wasen i")) & gdename == "e"
replace gdename = "sumiswald" if cname_extract == "wasen im emmental" & gdename == "e" // wasen i.e. is a locality in sumiswald

replace cname_extract = "affoltern im emmental" if (regexm(cname, "affoltern i")) & gdename == "e"
replace gdename = "affoltern im emmental" if cname_extract == "affoltern im emmental" & gdename == "e"

replace cname_extract = "roethenbach im emmental" if (regexm(cname, "röthenbach i")) & gdename == "e"
replace gdename = "roethenbach im emmental" if cname_extract == "roethenbach im emmental" & gdename == "e"

replace cname_extract = "st stephan" if (regexm(cname, "stephan")) & gdename == "e"
replace gdename = "st stephan" if cname_extract == "st stephan" & gdename == "e"

* 5.24) Extracting from gdname = "e": other option is "im engadin" (one part in cname and "e" in gdename)
replace cname_extract = "sils im engadin/segl" if (regexm(cname, "sils")) & gdename == "e"
replace gdename = "sils im engadin/segl" if cname_extract == "sils im engadin/segl" & gdename == "e"

* 5.25) Extracting from gdename = "johann" : it is either "alt st johann" or "neu st johann" which is a locality in Nesslau
replace cname_extract = "alt st johann" if (regexm(cname, "alt st")) & gdename == "johann"
replace gdename = "alt st johann" if cname_extract == "alt st johann" & gdename == "johann"

* 5.26) Extracting from gdename = "busswil" : many possibilities, I focus on those we can 100% distinguish
replace cname_extract = "busswil bei bueren" if (regexm(cname, "busswil")) & gdename == "bueren"
replace gdename = "busswil bei bueren" if cname_extract == "busswil bei bueren" & gdename == "bueren"

replace cname_extract = "diessbach bei bueren" if (regexm(cname, "diessbach")) & gdename == "bueren"
replace gdename = "diessbach bei bueren" if cname_extract == "diessbach bei bueren" & gdename == "bueren"

replace cname_extract = "oberwil bei bueren" if (regexm(cname, "oberwil")) & gdename == "bueren"
replace gdename = "oberwil bei bueren" if cname_extract == "oberwil bei bueren" & gdename == "bueren"

replace cname_extract = "rueti bei bueren" if (regexm(cname, "rüti")) & gdename == "bueren"
replace gdename = "rueti bei bueren" if cname_extract == "rueti bei bueren" & gdename == "bueren"

replace cname_extract = "wengi" if (regexm(cname, "wengi")) & gdename == "bueren"
replace gdename = "wengi" if cname_extract == "wengi" & gdename == "bueren"
*
replace cname_extract = "busswil bei bueren" if (regexm(cname, "busswil")) & gdename == "buren"
replace gdename = "busswil bei bueren" if cname_extract == "busswil bei bueren" & gdename == "buren"

replace cname_extract = "diessbach bei bueren" if (regexm(cname, "diessbach")) & gdename == "buren"
replace gdename = "diessbach bei bueren" if cname_extract == "diessbach bei bueren" & gdename == "buren"

replace cname_extract = "rueti bei bueren" if (regexm(cname, "rüti")) & gdename == "buren"
replace gdename = "rueti bei bueren" if cname_extract == "rueti bei bueren" & gdename == "buren"

replace cname_extract = "wengi" if (regexm(cname, "wengi")) & gdename == "buren"
replace gdename = "wengi" if cname_extract == "wengi" & gdename == "buren"

* 5.27) Extracting "erlenbach zh" or "erlenbach im simmental"  
replace cname_extract = "erlenbach zh" if (regexm(cname, "erlenbach zh")) & gdename == "erlenbach"
replace gdename = "erlenbach zh" if cname_extract == "erlenbach zh" & gdename == "erlenbach"

replace cname_extract = "erlenbach im simmental" if (regexm(cname, "erlenbach i.s")) & gdename == "erlenbach"
replace gdename = "erlenbach im simmental" if cname_extract == "erlenbach im simmental" & gdename == "erlenbach"

* 5.28) Extracting from gdename = "lengnau": either be (also called lengnau bei biel) or ag   
replace cname_extract = "lengnau be" if (regexm(cname, "lengnau b.")) & gdename == "lengnau"
replace gdename = "lengnau be" if cname_extract == "lengnau be" & gdename == "lengnau" 

replace cname_extract = "lengnau be" if (regexm(cname, "lengnau bei")) & gdename == "lengnau"
replace gdename = "lengnau be" if cname_extract == "lengnau be" & gdename == "lengnau" 

* 5.29) Extracting from gdename = "colombier": either ne or vd
replace cname_extract = "colombier vd" if (regexm(cname, "colombier vd")) & gdename == "colombier"
replace gdename = "colombier vd" if cname_extract == "colombier vd" & gdename == "colombier" 

replace cname_extract = "colombier ne" if (regexm(cname, "colombier ne")) & gdename == "colombier"
replace gdename = "colombier ne" if cname_extract == "colombier ne" & gdename == "colombier"

* 5.30) Extracting from gdename = "eschenbach": sg or lu
replace cname_extract = "eschenbach sg" if (regexm(cname, "eschenbach sg")) & gdename == "eschenbach"
replace gdename = "eschenbach sg" if cname_extract == "eschenbach sg" & gdename == "eschenbach" 
replace cname_extract = "eschenbach sg" if (regexm(cname, "neuhaus")) & gdename == "eschenbach"
replace gdename = "eschenbach sg" if cname_extract == "eschenbach sg" & gdename == "eschenbach" 

replace cname_extract = "eschenbach lu" if (regexm(cname, "eschenbach lu")) & gdename == "eschenbach"
replace gdename = "eschenbach lu" if cname_extract == "eschenbach lu" & gdename == "eschenbach"

* 5.31) Extracting gdename = "buerglen" & "buergten": either ur or tg 
replace cname_extract = "buerglen tg" if (regexm(cname, "bürglen tg")) & gdename == "buerglen"
replace gdename = "buerglen tg" if cname_extract == "buerglen tg" & gdename == "buerglen" 

replace cname_extract = "buerglen ur" if (regexm(cname, "bürglen ur")) & gdename == "buerglen"
replace gdename = "buerglen ur" if cname_extract == "buerglen ur" & gdename == "buerglen"

replace cname_extract = "buerglen tg" if (regexm(cname, "bürglen tg")) & gdename == "buergten"
replace gdename = "buerglen tg" if cname_extract == "buerglen tg" & gdename == "buergten" 

replace cname_extract = "buerglen ur" if (regexm(cname, "bürglen ur")) & gdename == "buerglen"
replace gdename = "buerglen ur" if cname_extract == "buerglen ur" & gdename == "buerglen"
replace gdename = "buerglen ur" if gdename == "buergten ur"

* 5.32) Extracting gdename = "niklaus": either feldbrunnen st niklaus or st niklaus
replace cname_extract = "feldbrunnen st niklaus" if (regexm(cname, "feldbrunnen")) & gdename == "niklaus"
replace gdename = "st niklaus" if gdename == "niklaus" 

* 5.33) Extracting gdename = "roggwil": either be or tg  
replace cname_extract = "roggwil tg" if (regexm(cname, "roggwil tg")) & gdename == "roggwil"
replace gdename = "roggwil tg" if cname_extract == "roggwil tg" & gdename == "roggwil" 

replace cname_extract = "roggwil tg" if (regexm(cname, "freidorf")) & gdename == "roggwil" // freidorf is a locality
replace gdename = "roggwil tg" if cname_extract == "roggwil tg" & gdename == "roggwil" 

* 5.34) Extracting gdename = "oberdorf": bl, nw or so  
replace cname_extract = "oberdorf bl" if (regexm(cname, "oberdorf (bl)")) & gdename == "oberdorf"
replace cname_extract = "oberdorf bl" if (regexm(cname, "oberdorf bl")) & gdename == "oberdorf"
replace gdename = "oberdorf bl" if cname_extract == "oberdorf bl" & gdename == "oberdorf"

replace cname_extract = "oberdorf nw" if (regexm(cname, "oberdorf nidwalden")) & gdename == "oberdorf"
replace cname_extract = "oberdorf nw" if (regexm(cname, "oberdorf nw")) & gdename == "oberdorf"
replace cname_extract = "oberdorf nw" if (regexm(cname, "büren")) & gdename == "oberdorf" // büren is a locality
replace cname_extract = "oberdorf nw" if (regexm(cname, "stans")) & gdename == "oberdorf" // also called oberdorf bei stans
replace gdename = "oberdorf nw" if cname_extract == "oberdorf nw" & gdename == "oberdorf"

replace cname_extract = "oberdorf bl" if (regexm(cname, "oberdorf bl")) & gdename == "oberdorf"
replace gdename = "oberdorf bl" if cname_extract == "oberdorf bl" & gdename == "oberdorf"

* 5.35) Extracting gdename = "zuzwil": either be or sg  
replace cname_extract = "zuzwil sg" if (regexm(cname, "zuberwangen")) & gdename == "zuzwil"
replace cname_extract = "zuzwil sg" if (regexm(cname, "züberwangen")) & gdename == "zuzwil"
replace cname_extract = "zuzwil sg" if (regexm(cname, "zuzwil sg")) & gdename == "zuzwil"
replace gdename = "zuzwil sg" if cname_extract == "zuzwil sg" & gdename == "zuzwil" 

replace cname_extract = "zuzwil be" if (regexm(cname, "zuzwil be")) & gdename == "zuzwil"
replace gdename = "zuzwil be" if cname_extract == "zuzwil be" & gdename == "zuzwil"

*5.36) Extracting gdename = "rickenbach": bl, lu, so, tg (also bei wil), zh  
replace cname_extract = "rickenbach lu" if (regexm(cname, "rickenbach/lu")) & gdename == "rickenbach"
replace cname_extract = "rickenbach lu" if (regexm(cname, "rickenbach lu")) & gdename == "rickenbach"
replace gdename = "rickenbach lu" if cname_extract == "rickenbach lu" & gdename == "rickenbach" 

replace cname_extract = "rickenbach bei wil" if (regexm(cname, "rickenbach bei wil")) & gdename == "rickenbach"
replace gdename = "rickenbach bei wil" if cname_extract == "rickenbach bei wil" & gdename == "rickenbach"

*5.37) Extracting gdename = "rorbas freienstein": two distinct mun  
replace cname_extract = "rorbas" if (regexm(cname, "rorbas")) & gdename == "rorbas freienstein"
replace gdename = "rorbas" if cname_extract == "rorbas" & gdename == "rorbas freienstein" 

replace cname_extract = "freienstein" if (regexm(cname, "freienstein")) & gdename == "rorbas freienstein"
replace gdename = "freienstein" if cname_extract == "freienstein" & gdename == "rorbas freienstein"

*5.38) Extracting gdename = "villars"  
replace cname_extract = "ollon" if (regexm(cname, "bretaye")) & gdename == "villars" // Bretaye (or Col de Bretaye) is a high mountain pass of the Swiss Alps, located above Villars-sur-Ollon (mun = Ollon)
replace cname_extract = "ollon" if (regexm(cname, "011on")) & gdename == "villars"
replace cname_extract = "ollon" if (regexm(cname, "villars-sur-ollon")) & gdename == "villars"
replace cname_extract = "ollon" if (regexm(cname, "roc d'orsay")) & gdename == "villars"
replace gdename = "ollon" if cname_extract == "ollon" & gdename == "villars" 

*5.39) Extracting gdename = "schmitten": fr or gr
replace cname_extract = "schmitten fr" if (regexm(cname, "ried")) & gdename == "schmitten"
replace cname_extract = "schmitten fr" if (regexm(cname, "fillistorf")) & gdename == "schmitten"
replace gdename = "schmitten fr" if cname_extract == "schmitten fr" & gdename == "schmitten" 

*5.40) Extracting gdename = "langnau": langnau bei reiden, langnau im emmental, langnau am albis
replace cname_extract = "langnau im emmental" if (regexm(cname, "emmental")) & gdename == "langnau"
replace cname_extract = "langnau im emmental" if (regexm(cname, "emmenthal")) & gdename == "langnau"
replace cname_extract = "langnau im emmental" if (regexm(cname, "langnau i")) & gdename == "langnau"
replace gdename = "langnau im emmental" if cname_extract == "langnau im emmental" & gdename == "langnau" 

replace cname_extract = "langnau am albis" if (regexm(cname, "langnau a")) & gdename == "langnau"
replace gdename = "langnau am albis" if cname_extract == "langnau am albis" & gdename == "langnau"

* 5.41) Extracting from gdename = "s" or "see" 
replace cname_extract = "beinwil am see" if (regexm(cname, "beinwil a")) & gdename == "s"
replace gdename = "beinwil am see" if cname_extract == "beinwil am see" & gdename == "s" 

replace cname_extract = "erlenbach im simmental" if (regexm(cname, "erlenbach i")) & gdename == "s"
replace gdename = "erlenbach im simmental" if cname_extract == "erlenbach im simmental" & gdename == "s" 

replace cname_extract = "oberwil im simmental" if (regexm(cname, "oberwil i")) & gdename == "s"
replace gdename = "oberwil im simmental" if cname_extract == "oberwil im simmental" & gdename == "s" 

replace cname_extract = "lenk" if (regexm(cname, "lenk i")) & gdename == "s"
replace gdename = "lenk" if cname_extract == "lenk" & gdename == "s" 

replace cname_extract = "boltigen" if (regexm(cname, "boltigen i")) & gdename == "s"
replace gdename = "boltigen" if cname_extract == "boltigen" & gdename == "s" 

replace cname_extract = "oetwil am see" if (regexm(cname, "oetwil a")) & gdename == "s"
replace gdename = "oetwil am see" if cname_extract == "oetwil am see" & gdename == "s" 

replace cname_extract = "uetikon am see" if (regexm(cname, "uetikon a")) & gdename == "s"
replace gdename = "uetikon am see" if cname_extract == "uetikon am see" & gdename == "s"
replace gdename = "uetikon" if gdename == "uetikon am see" & year < 1977
**
replace cname_extract = "beinwil am see" if (regexm(cname, "beinwil a")) & gdename == "see"
replace gdename = "beinwil am see" if cname_extract == "beinwil am see" & gdename == "see" 

replace cname_extract = "oetwil am see" if (regexm(cname, "oetwil a")) & gdename == "see"
replace gdename = "oetwil am see" if cname_extract == "oetwil am see" & gdename == "see" 

replace cname_extract = "uetikon am see" if (regexm(cname, "uetikon a")) & gdename == "see"
replace gdename = "uetikon am see" if cname_extract == "uetikon am see" & gdename == "see"
replace gdename = "uetikon" if gdename == "uetikon am see" & year < 1977

* 5.42) Extracting gdename = "stein": ar, ag, sg (=toggenburg before 1994), sh (=am rhein)  
replace cname_extract = "stein ar" if (regexm(cname, "stein ar")) & gdename == "stein"
replace gdename = "stein ar" if cname_extract == "stein ar" & gdename == "stein" 

replace cname_extract = "stein toggenburg" if (regexm(cname, "toggenburg")) & gdename == "stein"
replace gdename = "stein toggenburg" if cname_extract == "stein toggenburg" & gdename == "stein"

* 5.43) Extracting gdename = "r" (rigi, rhein, etc)  
replace cname_extract = "kuessnacht am rigi" if (regexm(cname, "küsnacht")) & gdename == "r"
replace cname_extract = "kuessnacht am rigi" if (regexm(cname, "küssnacht")) & gdename == "r"
replace gdename = "kuessnacht am rigi" if cname_extract == "kuessnacht am rigi" & gdename == "r" 

replace cname_extract = "neuhausen am rheinfall" if (regexm(cname, "neuhausen")) & gdename == "r"
replace gdename = "neuhausen am rheinfall" if cname_extract == "neuhausen am rheinfall" & gdename == "r"

replace cname_extract = "stein am rhein" if (regexm(cname, "stein a")) & gdename == "r"
replace gdename = "stein am rhein" if cname_extract == "stein am rhein" & gdename == "r"

* 5.44) Extracting gdename = "holderbank"  
replace cname_extract = "holderbank ag" if (regexm(cname, "wildegg")) & gdename == "holderbank"
replace gdename = "holderbank ag" if cname_extract == "holderbank ag" & gdename == "holderbank" 

* 5.45) Extracting gdename = "willisau": willisau land, willisau stadt, hergiswil bei willisau  
replace cname_extract = "hergiswil bei willisau" if (regexm(cname, "hergiswil")) & gdename == "willisau"
replace gdename = "hergiswil bei willisau" if cname_extract == "hergiswil bei willisau" & gdename == "willisau" 

* 5.46) Extracting gdename = "corcelles"  
replace cname_extract = "corcelles cormondreche" if (regexm(cname, "peseux")) & gdename == "corcelles"
replace gdename = "corcelles cormondreche" if cname_extract == "corcelles cormondreche" & gdename == "corcelles" 

* 5.47) Extracting gdename = "hofstetten"  
replace cname_extract = "hofstetten so" if (regexm(cname, "flüh")) & gdename == "hofstetten"
replace gdename = "hofstetten so" if cname_extract == "hofstetten so" & gdename == "hofstetten" 
replace gdename = "hofstetten flueh" if gdename == "hofstetten so" & year > 1985

* 5.48) Extracting gdename = "forel"  
replace cname_extract = "forel lavaux" if (regexm(cname, "lavaux")) & gdename == "forel"
replace gdename = "forel lavaux" if cname_extract == "forel lavaux" & gdename == "forel" 

* 5.49) Extracting gdename = "oberhofen"  
replace cname_extract = "oberhofen bei kreuzlingen" if (regexm(cname, "lengwil")) & gdename == "oberhofen"
replace gdename = "oberhofen bei kreuzlingen" if cname_extract == "oberhofen bei kreuzlingen" & gdename == "oberhofen" 

/* Example
* 5.) Extracting gdename = ""  
replace cname_extract = "" if (regexm(cname, "")) & gdename == ""
replace gdename = "" if cname_extract == "" & gdename == "" 

replace cname_extract = "" if (regexm(cname, "")) & gdename == ""
replace gdename = "" if cname_extract == "" & gdename == ""
*/

foreach str in gdename /*cname*/ {
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
	replace `str' = usubinstr(`str', "-", " ",.)
	replace `str' = usubinstr(`str', "~", " ",.)	
	replace `str' = usubinstr(`str', ":", " ",.)
	replace `str' = usubinstr(`str', "    ", " ",.)
	replace `str' = usubinstr(`str', "   ", " ",.)
	replace `str' = usubinstr(`str', "  ", " ",.)
	replace `str' = ustrtrim(`str')
	replace `str' = ustrlower(`str')
	replace `str' = ustrtrim(`str')
}

sort gdename year
drop gdenr gdenr_2018 ctn E_CNTR N_CNTR geo_merge
merge m:1 gdename year using "$path\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"
drop if _merge == 2
gen geo_merge = 8 if _merge == 3

preserve
keep if _merge == 3
save "$path\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_8.dta", replace
restore


*Useful commands and example
keep if _merge == 1
drop _merge
save "$path\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_remaining.dta", replace



/*
gen obs = 1
collapse (sum) obs, by(gdename)
replace obs = obs * -1
sort obs
*/

/*
replace gdename = ustrregexra(gdename,"3oll","3011")  // replace "3oll" in string (bern)
gen PLZ4_2 = regexs(1) if (regexm(gdename, "([0-9][0-9][0-9][0-9])")) // finds 4-digit numbers within string
gen nomcap = regexs(0) if regexm(gdename, "\([0-9]+\.[0-9]+\)") // extracts information in parenthesis (the info on nominal capital is represented that way)
replace gdename = ustrregexra(gdename,"([0-9][0-9][0-9][0-9])","")  // deletes 4-digit numbers coded as PLZ4_2
replace days2 = 30 if(regexm(days, "[eE]very[ ]*[dD]ay[.]*")) // Replace "every day" with 30
*/


