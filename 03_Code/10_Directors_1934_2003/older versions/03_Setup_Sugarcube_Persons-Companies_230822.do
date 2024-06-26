version 17
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\"
global temp "C:\Users\schelkem\Desktop\"

local date $S_DATE

log using "$path\Sugar_log_Pers-Firm_`date'", replace

*** Merge Person and Company Info from Sugarcube

mkdir "$temp\SugarCompMerge_tmp\", public  // create directory for all the tmp files.

use "$path\Sugarcube_Person_CleanName-Gender-Geo.dta", clear

** Steps
* a) Separate coma-separted mandate & function information and put in long format
* b) Merge mandate & function information in long format
* c) RL Post-Processing: Prepare company info for post-processing procedures  
* d) Create Person-Company Panel: merge with company information
* e) RL Panel ID: Prepare company info for procedures to establish Ground Truth


/*
generate random = runiform()
sort random
generate insample = _n <= 10000 
keep if insample == 1
save "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_randomsample.dta", replace
*/

* a) Separate coma-separted mandate and function information

// split function is very slow...
// alternative: read out relevant Obs with company ID as csv (coma delimited) and read back in (do separately for function)

export delim PID year ID_dupl gde_options companiesid using "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_CID_tmp.csv", delimiter(",") replace
export delim PID year ID_dupl gde_options functions using "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_functions_tmp.csv", delimiter(",") replace

* Reshape mandates
import delim using "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_CID_tmp.csv", delimiter(",") encoding("utf-8") bindquotes(nobind) stringcols(_all) clear
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
gen mandID = 1
preserve
drop cv*
save "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_mand_1.dta", replace
restore

forv i = 6(1)`maxvC'  {
	keep if cv`i' != ""
	preserve
	replace tmp_CID = cv`i' // alternative approach (not via split)
	replace mandID = `i' - 4
	drop cv*
	local h = `i' - 4
	save "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_mand_`h'.dta", replace
	restore
}
use "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_mand_`noCID'.dta", clear
local m = `noCID'  - 1
forv i = `m'(-1)1 { 
append using "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_mand_`i'.dta"
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

label var mandID "Mandate #"
sort PID year ID_dupl gde_options mandID
save "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_mand_long.dta", replace

* Reshape functions
import delim using "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_functions_tmp.csv", delimiter(",") encoding("utf-8") bindquotes(nobind) stringcols(_all) clear
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
gen mandID = 1
preserve
drop fv*
save "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_fct_1.dta", replace
restore

forv i = 6(1)`maxvF'  {
	keep if fv`i' != ""
	preserve
	replace tmp_funct = fv`i' // alternative approach (not via split)
	replace mandID = `i' - 4
	drop fv*
	local h = `i' - 4
	save "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_fct_`h'.dta", replace
	restore
}
use "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_fct_`noMandateF'.dta", clear
local m = `noMandateF'  - 1
forv i = `m'(-1)1 {  
append using "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_fct_`i'.dta"
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

label var mandID "Mandate #"
sort PID year ID_dupl gde_options mandID
save "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_fct_long.dta", replace


* b) Merge mandate & function information with original Person information in long format

use "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_mand_long.dta", clear
merge 1:1 PID year ID_dupl gde_options mandID using "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_fct_long.dta"

preserve
keep if _merge == 1
gen n = 1
collapse (sum) n, by(PID year ID_dupl gde_options)
save "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_mandfct_nomatch.dta", replace
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


/* c) RL Post-Processing: Prepare company info for post-processing procedures  

* Prepare dataset for RL post-processing and manual FP cleaning (Nationalrat-VR RL)  -- > prepared for Sandro Cilurzo May 2023

* Merge with company information
preserve
use "$path\comp_geo-merge_temp_20230502.dta", clear  // company data  -- > use latest version (Email Yannick, 02.05.2023)
rename gdename firmlocation
rename cname firmname
rename gdenr_2018 firm_gdenr_2018
rename E_CNTR firm_E_CNTR
rename N_CNTR firm_N_CNTR
rename ctn firm_ctn
rename geo_merge firm_geo_merge
drop if CID == .
format CID %10.0f
keep CID year firmname firmlocation firm_gdenr_2018 firm_E_CNTR firm_N_CNTR firm_ctn firm_geo_merge owners capital function signature cname_orig caddress_orig gdename_orig branch_dummy si_dummy ag_dummy
sort CID year
duplicates report CID year
save "$temp\SugarCompMerge_tmp\comp_1934-2003_merge.dta", replace
restore

merge m:1 CID year using "$temp\SugarCompMerge_tmp\comp_1934-2003_merge.dta"
drop if _merge == 2
drop _merge
sort PID year ID_dupl gde_options CID
save "$path\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", replace

use "$path\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", clear
sort PID year ID_dupl gde_option mandID CID 
egen UPID = group(PID year ID_dupl gde_options)

preserve
keep UPID mandID firmname
replace firmname = strproper(firmname)
replace firmname = usubinstr(firmname, "Ag", "AG",.)
replace firmname = usubinstr(firmname, "Sa", "SA",.)
replace firmname = usubinstr(firmname, "Si", "SI",.)
sort UPID firmname
replace firmname = usubinstr(firmname, ";", "",.)
bysort UPID (firmname): gen firmstr = firmname[1]  // alphabetical order of firm names
by UPID: replace firmstr = firmstr[_n-1] + "; " + firmname if _n > 1 
by UPID: replace firmstr = firmstr[_N]
bysort UPID: gen n = _n
keep if n == 1
drop n firmname
rename firmstr firmnames
save "$temp\SugarCompMerge_tmp\Sugarcube_Person_Company_Names_wide_alpha.dta", replace
restore

keep if mandID == 1
drop CID compFunct mandate_flag firmname firmlocation owners capital function signature // these are variables for the long-format
sort UPID
merge 1:1 UPID using "$temp\SugarCompMerge_tmp\Sugarcube_Person_Company_Names_wide_alpha.dta"
keep id_sug year_sug firstname lastname gdename male GdeNr_E_CNTR GdeNr_N_CNTR ctn language_w e_id_sug n_id_sug firmnames ID_dupl gde_options geo_merge

duplicates report id_sug year_sug e_id_sug n_id_sug
order id_sug year_sug e_id_sug n_id_sug
sort id_sug year_sug e_id_sug n_id_sug ID_dupl
*save "$path\Sugarcube_RLPostProcessing_Persons-Firmnames.dta", replace
export delimited using "$path\Sugarcube_RLPostProcessing_Persons-Firmnames.csv", delimiter(",") replace
*/


* d) Create Person-Company Panel: merge with company information

use "$path\Sugarcube_Person_CleanName-Gender-Geo_long.dta", clear
duplicates tag PID year ID_dupl gde_options CID, gen(CIDdupl)
tab CIDdupl
tab year if CIDdupl>0
*br if CIDdupl != 0
/*
------------+-----------------------------------
          0 |  8,290,836       99.87       99.87
          1 |     10,884        0.13      100.00
          2 |         48        0.00      100.00
          3 |         20        0.00      100.00
         17 |         18        0.00      100.00
         47 |         48        0.00      100.00
         87 |         88        0.00      100.00
------------+-----------------------------------
      Total |  8,301,942      100.00

The original Sucarcube-Persons data contain for 5640 cases duplicates in Firm-IDs per unique person (ID year ID_dupl gde_options)
*/
preserve
keep if CIDdupl >0
gen n = 1
collapse (sum) n, by(PID year ID_dupl gde_options firstname lastname city gdename gdenr_2018 page CID)
save "$temp\SugarCompMerge_tmp\Sugarcube_Person_duplicateMandates.dta", replace
restore

duplicates drop PID year ID_dupl gde_options CID, force  // (5,640 observations deleted)
drop CIDdupl

preserve
use "$path\comp_geo-merge_temp_20230810.dta", clear  // company data with geo-localization Stand 10.08.2023
rename gdename firmlocation
rename cname firmname
rename gdenr_2018 firm_gdenr_2018
rename E_CNTR firm_E_CNTR
rename N_CNTR firm_N_CNTR
rename ctn firm_ctn
*rename geo_merge firm_geo_merge
drop if CID == .
format CID %10.0f
keep CID year firmname firmlocation firm_gdenr_2018 firm_E_CNTR firm_N_CNTR firm_ctn owners capital function signature cname_orig caddress_orig gdename_orig branch_dummy
sort CID year
duplicates report CID year
save "$temp\SugarCompMerge_tmp\comp_1934-2003_merge.dta", replace
restore

merge m:1 CID year using "$temp\SugarCompMerge_tmp\comp_1934-2003_merge.dta"
drop if _merge == 2
drop _merge
sort PID year ID_dupl gde_options CID
save "$path\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", replace


* e) Panel ID: Prepare company info for procedures to establish Ground Truth

* Prepare Ground Truth for Person Panel ID (August 2023)
*use "$path\Sugarcube_Person_CleanName-Gender-Geo_Companies.dta", clear
sort PID year ID_dupl gde_option mandID CID 
egen UPID = group(PID year ID_dupl gde_options)

preserve
keep UPID mandID firmname firmlocation
replace firmname = strproper(firmname)
replace firmname = usubinstr(firmname, "&" , " & ",.)
replace firmname = " " + firmname + " "
replace firmname = usubinstr(firmname, " Ag ", " AG ",.)
replace firmname = usubinstr(firmname, " Sa ", " SA ",.)
replace firmname = usubinstr(firmname, " Si ", " SI ",.)
replace firmname = usubinstr(firmname, "«" , "",.)
replace firmname = usubinstr(firmname, "»" , "",.)
replace firmname = usubinstr(firmname, "  " , " ",.)
replace firmname = usubinstr(firmname, "  " , " ",.)
replace firmname = ustrtrim(firmname)
replace firmlocation = strproper(firmlocation)
foreach ctn in Ag Ai Ar Be Bl Bs Fr Ge Gl Gr Ju Lu Ne Nw Ow Sg Sh So Sz Tg Ti Ur Vd Vs Zg Zh {
	replace firmlocation = usubinstr(firmlocation, "(`ctn')" , ustrupper("(`ctn')"),.)
}
gen firmnamelocation = firmname + " " + "(" + firmlocation + ")"
replace firmnamelocation = usubinstr(firmnamelocation, "  " , " ",.)
drop firmname firmlocation
sort UPID firmnamelocation
replace firmnamelocation = usubinstr(firmnamelocation, ";", "",.)
bysort UPID (firmnamelocation): gen firmstr = firmnamelocation[1]  // alphabetical order of firm names
by UPID: replace firmstr = firmnamelocation[_n-1] + "; " + firmnamelocation if _n > 1 
by UPID: replace firmstr = firmstr[_N]
bysort UPID: gen n = _n
keep if n == 1
drop n firmnamelocation mandID
rename firmstr firmnamelocation
save "$temp\SugarCompMerge_tmp\Sugarcube_Person_CompanyNames-Location_wide_alpha.dta", replace
restore

keep if mandID == 1
drop CID compFunct mandate_flag firmname firmlocation owners capital function signature // these are variables for the long-format
sort UPID
merge 1:1 UPID using "$temp\SugarCompMerge_tmp\Sugarcube_Person_CompanyNames-Location_wide_alpha.dta"
keep PID year firstname lastname city gdename male GdeNr_E_CNTR GdeNr_N_CNTR ctn language_w firmnamelocation ID_dupl gde_options geo_merge

duplicates report PID year ID_dupl gde_options
order PID year ID_dupl gde_options
sort PID year ID_dupl gde_options
*save "$path\Sugarcube_PersPanelID_Firmnames_wide.dta", replace
export delimited using "$path\Sugarcube_PersPanelID_Firmnames_wide.csv", delimiter(",") replace


cap log close

forv i = 1(1)`noCID'  {
erase "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_mand_`i'.dta"
}
forv i = 1(1)`noMandateF'  {
erase "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_split_fct_`i'.dta"
}

erase "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_CID_tmp.csv"
erase "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_functions_tmp.csv"
erase "$temp\SugarCompMerge_tmp\comp_1934-2003_merge.dta"
erase "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_mandfct_nomatch.dta"
erase "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_fct_long.dta"
erase "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_mand_long.dta"
erase "$temp\SugarCompMerge_tmp\Sugarcube_Person_CleanName-Gender-Geo_Companies_RLpanel_IDs_tmp.dta"
erase "$temp\SugarCompMerge_tmp\Sugarcube_Person_CompanyNames-Location_wide_alpha.dta"
rmdir "$temp\SugarCompMerge_tmp"
*erase "$path\Sugarcube_Person_CleanName-Gender-Geo_randomsample.dta"
