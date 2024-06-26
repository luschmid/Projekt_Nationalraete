version 15
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\"
global dump "C:\Users\schelkem\Desktop\"

local date $S_DATE

*log using "$path\Sugar_log_Pers-Firm_`date'", replace

*** Merge Person and Company Info from Sugarcube

mkdir "$dump\SugarCompMerge_tmp\", public  // create directory for all the tmp files.

use "$path\Sugarcube_Person_CleanName-Gender-Geo.dta", clear

/*
generate random = runiform()
sort random
generate insample = _n <= 10000 
keep if insample == 1
save "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_randomsample.dta", replace
*/

// split function is very slow...
// alternative: read out relevant Obs with company ID as csv (coma delimited) and read back in (do separately for function)

export delim PID year ID_dupl gde_options companiesid using "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_CID_tmp.csv", delimiter(",") replace
export delim PID year ID_dupl gde_options functions using "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_functions_tmp.csv", delimiter(",") replace

* Reshape mandates
import delim using "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_CID_tmp.csv", delimiter(",") encoding("utf-8") bindquotes(nobind) stringcols(_all) clear
rename companiesid v5
for var *: rename X cX
rename cpid PID
rename cyear year
rename cid_dupl ID_dupl
rename cgde_options gde_options
describe _all 
local maxvC = r(k)
local noCID = `maxvC' - 4

/*/ check whether mandates are perfectly "successive" (i.e., one with an 8th entry must have a 7th): it has to be the case and IS THE CASE
gen checkmandateC = 1
forv i = 6(1)`maxvC'  {
	local j = `i' - 1
	replace checkmandateC = 0 if  cv`i' != "" & cv`j' == ""  // alternative approach (not via split)
}
tab checkmandateC
drop checkmandateC
*/

// Alternative to reshape: Create subsamples to reassemble entire dataset in long format
gen tmp_CID = cv5
preserve
drop cv*
save "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_mand_1.dta", replace
restore

forv i = 6(1)`maxvC'  {
	keep if cv`i' != ""
	preserve
	replace tmp_CID = cv`i' // alternative approach (not via split)
	drop cv*
	local h = `i' - 4
	save "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_mand_`h'.dta", replace
	restore
}
use "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_mand_`noCID'.dta", clear
local m = `noCID'  - 1
forv i = `m'(-1)1 { 
append using "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_mand_`i'.dta"
}

/*
// there are nonnumeric characters in the companies ID (all 1962)
*ssc install egenmore
egen CIDn= sieve(CID), omit(0123456789)
br * if CIDn != ""
*/

rename tmp_CID CID
replace CID = usubinstr(CID, `"""', "",.)
replace CID = "9999999999" if CID =="Philips AG"
replace CID = "9999999999" if CID =="Hewlett-Pack"
replace CID = "9999999999" if CID =="Herbert Ammann"
replace CID = "9999999999" if CID =="Skilift Müsella AG"
replace CID = "9999999999" if CID =="Schweiz. "
replace CID = "9999999999" if CID ==" Kumag AG Maschinenfabrik" 
replace CID = "9999999999" if CID ==" Ems-Gelsenberg"
replace CID = "9999999999" if CID =="Pavesa AG"
replace CID = "9999999999" if CID =="Buchdruckerei zur Froschau"
replace CID = "9999999999" if CID =="Myceta S.A."
replace CID = "9999999999" if CID =="Progress Lederwaren"
replace CID = "9999999999" if CID =="Rhenum Metal"
replace CID = "9999999999" if CID=="  Zurich (0.3)"
destring CID, replace
format CID %10.0f
destring PID, replace
format PID %10.0f
destring year, replace
format year %10.0f
destring ID_dupl, replace
format ID_dupl %10.0f
destring gde_options, replace
format gde_options %10.0f

bysort PID year ID_dupl gde_options: gen mandID = _n
label var mandID "Mandate #"
sort PID year ID_dupl gde_options mandID
save "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_mand_long.dta", replace

* Reshape functions
import delim using "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_functions_tmp.csv", delimiter(",") encoding("utf-8") bindquotes(nobind) stringcols(_all) clear
rename functions v5
for var *: rename X fX
rename fpid PID
rename fyear year
rename fid_dupl ID_dupl
rename fgde_options gde_options
describe _all 
local maxvF = r(k)
local noMandateF = `maxvF' - 4

/*/ check whether mandates are perfectly "successive" (i.e., one with an 8th entry must have a 7th): it has to be the case and IS THE CASE
gen checkmandateF = 1
forv i = 6(1)`maxvF'  {
	local j = `i' - 1
	replace checkmandateF = 0 if  fv`i' != "" & fv`j' == ""  // alternative approach (not via split)
}
tab checkmandateF
drop checkmandateF
*/

// Alternative to reshape: Create subsamples to reassemble entire dataset in long format
gen tmp_funct = fv5
preserve
drop fv*
save "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_fct_1.dta", replace
restore

forv i = 6(1)`maxvF'  {
	keep if fv`i' != ""
	preserve
	replace tmp_funct = fv`i' // alternative approach (not via split)
	drop fv*
	local h = `i' - 4
	save "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_fct_`h'.dta", replace
	restore
}
use "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_fct_`noMandateF'.dta", clear
local m = `noMandateF'  - 1
forv i = `m'(-1)1 {  
append using "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_fct_`i'.dta"
}
rename tmp_funct compFunct
replace compFunct = usubinstr(compFunct, `"""', "",.)
destring PID, replace
format PID %10.0f
destring year, replace
format year %10.0f
destring ID_dupl, replace
format ID_dupl %10.0f
destring gde_options, replace
format gde_options %10.0f

bysort PID year ID_dupl gde_options: gen mandID = _n
label var mandID "Mandate #"
sort PID year ID_dupl gde_options mandID
save "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_fct_long.dta", replace

* Merge mandates & functions, merge with original data
use "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_mand_long.dta", clear
merge 1:1 PID year ID_dupl gde_options mandID using "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_fct_long.dta"

preserve
keep if _merge == 1
gen n = 1
collapse (sum) n, by(PID year ID_dupl gde_options)
save "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_mandfct_nomatch.dta", replace
restore

gen flaghelp1 = 1 if _merge!=3  // some PID have more mandates than functions: identify those cases
bysort PID year ID_dupl gde_options: egen flaghelp2 = total(flaghelp1)
gen mandate_flag = 0
replace mandate_flag = 1 if flaghelp2 > 0
label var mandate_flag "IDs where # mandates != # functions"
drop _merge flaghelp1 flaghelp2

merge m:1 PID year ID_dupl gde_options using "$path\Sugarcube_Person_CleanName-Gender-Geo.dta"
drop _merge

order PID year ID_dupl gde_options options4a options4b options_manual titlejob firstname lastname male address address2 alternateaddress PLZ4 city ///
companiesid mandID functions page geo_merge countPLZ4 gdename gdenr_2018 GdeNr_E_CNTR GdeNr_N_CNTR nCH ctn language_w firstname_orig lastname_orig namecorrection ///
id_sug year_sug e_id_sug n_id_sug CID compFunct mandate_flag
label var CID "Company identifier for merge with company data"
label var compFunct "Function in company"
sort PID year ID_dupl gde_options CID
save "$path\Sugarcube_Person_CleanName-Gender-Geo_long.dta", replace

* Merge with company information
use "$path\Sugarcube_Person_CleanName-Gender-Geo_long.dta", clear

preserve
use "$path\comp_1934-2003_tmp1.dta", clear  // company data  -- > use latest version
rename id CID
rename city firmlocation
rename name firmname
drop if CID == "."
destring CID, replace
format CID %10.0f
destring year, replace
format year %4.0f
drop address // empty
sort CID year
duplicates report CID year
save "$dump\SugarCompMerge_tmp\comp_1934-2003_merge.dta", replace
restore

merge m:1 CID year using "$dump\SugarCompMerge_tmp\comp_1934-2003_merge.dta"
drop if _merge == 2
drop _merge
sort PID year ID_dupl gde_options CID

save "$path\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", replace

* Prepare dataset for Sandro RL: reshape into wide, single, coma-separated mandate variable
use "$path\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", clear
sort PID year ID_dupl gde_option mandID CID 
egen UPID = group(PID year ID_dupl gde_options)

preserve
keep UPID mandID firmname
replace firmname = usubinstr(firmname, ";", "",.)
bysort UPID (mandID): gen firmstr = firmname[1]
by UPID: replace firmstr = firmstr[_n-1] + ";" + firmname if _n > 1 
by UPID: replace firmstr = firmstr[_N]
bysort UPID: gen n = _n
keep if n == 1
drop n firmname
rename firmstr firmnames
save "$path\Sugarcube_Person_Company_Names_wide.dta", replace
restore

keep if mandID == 1
drop CID compFunct mandate_flag firmname firmlocation owners capital function signature // these are variables for the long-format
sort UPID
merge 1:1 UPID using "$path\Sugarcube_Person_Company_Names_wide.dta"
keep id_sug year_sug e_id_sug n_id_sug firstname lastname gdename male GdeNr_E_CNTR GdeNr_N_CNTR ctn language_w e_id_sug n_id_sug firmnames

duplicates report id_sug year_sug e_id_sug n_id_sug
order id_sug year_sug e_id_sug n_id_sug
sort id_sug year_sug e_id_sug n_id_sug 
save "$path\Sugarcube_RecordLinkage_Persons-Firmnames.dta", replace
export delimited using "$path\Sugarcube_RecordLinkage_Persons-Firmnames.csv", delimiter(",") replace

cap log close

forv i = 1(1)`noCID'  {
erase "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_mand_`i'.dta"
}
forv i = 1(1)`noMandateF'  {
erase "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_fct_`i'.dta"
}

erase "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_CID_tmp.csv"
erase "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_functions_tmp.csv"
erase "$dump\SugarCompMerge_tmp\comp_1934-2003_merge.dta"
erase "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_mandfct_nomatch.dta"
erase "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_fct_long.dta"
erase "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_mand_long.dta"
erase "$dump\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_Companies_RLpanel_IDs_tmp.dta"
rmdir "$dump\SugarCompMerge_tmp"
*erase "$path\Sugarcube_Person_CleanName-Gender-Geo_randomsample.dta"
