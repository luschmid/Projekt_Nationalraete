version 17
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
global archive "E:\12. Cloud\Dropbox\Projekt Nationalräte\07_Archive\10_Directors_1934_2003\06_Old_person_files_(before_correction_in_nov_2020)\"
global temp "C:\Users\schelkem\Desktop\"

mkdir "$temp\SugarPersonsCheck_tmp\", public  // create directory for all the tmp files.

*** Compare the Companies ID's and functions from the original data delivery by Sugarcube with those after corrections @ UniLU

** Import original Sugarcube csv files (saved as UTF-8 with BOM) - before corrections @ UniLU
foreach val in 1934 1943 1960 1962 1963 1964 {
	import delim using "$archive\Kern_`val'_persons.csv", delimiters(";") encoding("utf-8") clear
	tostring *, replace
	save "$temp\SugarPersonsCheck_tmp\pers_`val'_orig.dta", replace
}
foreach val in 1966 1967 1970 1973 1976 1980_1981 {
	import delim using "$archive\Kern_`val'_persons.csv", delimiters(";") encoding("utf-8") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$temp\SugarPersonsCheck_tmp\pers_`val'_orig.dta", replace
}
foreach val in 1979_1980 {
	import delim using "$archive\Kern_`val'_persons.csv", encoding("utf-8") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$temp\SugarPersonsCheck_tmp\pers_`val'_orig.dta", replace
}
forv val = 1982(1)1985 {
	import delim using "$archive\Kern_`val'_persons.csv", delimiters(";") encoding("utf-8") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$temp\SugarPersonsCheck_tmp\pers_`val'_orig.dta", replace
}
use "$temp\SugarPersonsCheck_tmp\pers_1934_orig.dta", clear
foreach val in 1943 1960 1962 1963 1964 1966 1967 1970 1973 1976 1979_1980 1980_1981 1982 1983 1984 1985 {
	append using "$temp\SugarPersonsCheck_tmp\pers_`val'_orig.dta", force
	erase "$temp\SugarPersonsCheck_tmp\pers_`val'_orig.dta"
}
erase "$temp\SugarPersonsCheck_tmp\pers_1934_orig.dta"
keep id year companiesid functions firstname lastname
rename companiesid companiesid_o
rename functions functions_o
rename id PID
label var PID "Person ID per year"
destring PID, replace
destring year, replace
gen nbrcompIDs_o = length(companiesid_o) - length(subinstr(companiesid_o, ",", "", .)) + 1
gen nbrfunctions_o = length(functions_o) - length(subinstr(functions_o, ",", "", .)) + 1
gen d_compID_funct_o = nbrcompIDs_o - nbrfunctions_o
tab d_compID_funct_o 
tab nbrcompIDs_o if d_compID_funct_o != 0
sort PID year
duplicates report PID year
duplicates tag PID year, gen(dupl_o)
gen compmerge = companiesid
recast str2045 compmerge, force
replace compmerge = "" if dupl_o == 0
gen firstname_merge = ""
replace firstname_merge = firstname if dupl_o > 0
gen lastname_merge = ""
replace lastname_merge = lastname if dupl_o > 0
drop firstname lastname
save "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-1985_original.dta", replace

** Use corrected (@UniLU) raw data
foreach val in 1934 1943 1960 1962 1963 1964 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", delimiters(";") encoding("utf-8") clear
	tostring *, replace
	save "$temp\SugarPersonsCheck_tmp\pers_`val'_corr.dta", replace
}
foreach val in 1966 1967 1970 1973 1976 1980_1981 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", delimiters(";") encoding("utf-8") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$temp\SugarPersonsCheck_tmp\pers_`val'_corr.dta", replace
}
foreach val in 1979_1980 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", encoding("utf-8") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$temp\SugarPersonsCheck_tmp\pers_`val'_corr.dta", replace
}
forv val = 1982(1)1985 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", delimiters(";") encoding("utf-8") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$temp\SugarPersonsCheck_tmp\pers_`val'_corr.dta", replace
}
use "$temp\SugarPersonsCheck_tmp\pers_1934_corr.dta", clear
foreach val in 1943 1960 1962 1963 1964 1966 1967 1970 1973 1976 1979_1980 1980_1981 1982 1983 1984 1985 {
	append using "$temp\SugarPersonsCheck_tmp\pers_`val'_corr.dta", force
	erase "$temp\SugarPersonsCheck_tmp\pers_`val'_corr.dta"
}
erase "$temp\SugarPersonsCheck_tmp\pers_1934_corr.dta"
keep id year companiesid functions firstname lastname
rename id PID
label var PID "Person ID per year"
destring PID, replace
destring year, replace
gen nbrcompIDs = length(companiesid) - length(subinstr(companiesid , ",", "", .)) + 1
gen nbrfunctions = length(functions) - length(subinstr(functions, ",", "", .)) + 1
gen d_compID_funct = nbrcompIDs - nbrfunctions
tab d_compID_funct 
tab nbrcompIDs if d_compID_funct != 0
duplicates report PID year
duplicates tag PID year, gen(dupl_c)
gen compmerge = companiesid
recast str2045 compmerge, force
replace compmerge = "" if dupl_c == 0
gen firstname_merge = ""
replace firstname_merge = firstname if dupl_c > 0
gen lastname_merge = ""
replace lastname_merge = lastname if dupl_c > 0
drop firstname lastname
/*
--------------------------------------
   Copies | Observations       Surplus
----------+---------------------------
        1 |      1199278             0
        2 |          132            66
        3 |            3             2
--------------------------------------

# Obs: 1,199,413


ORIGINAL Sugarcube: Duplicates in terms of PID year

--------------------------------------
   Copies | Observations       Surplus
----------+---------------------------
        1 |      1199285             0
        2 |          126            63
        3 |            3             2
--------------------------------------

# Obs: 	1,199,414
*/
save "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-1985_corr_UniLU.dta", replace


** Identify difference between original and corrected functions:
* Merge original and corrected CID's, check functions

*use "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-1985_corr_UniLU.dta", clear
merge 1:1 PID year compmerge firstname_merge lastname_merge using "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-1985_original.dta" // we require additional ID vars due to ID duplicates in PID year

*br if _merge!=3
drop if _merge == 2

gen compid_diff = companiesid_o != companiesid // identify differences in CID's due to potential corrections @UniLU

tab d_compID_funct compid_diff
keep if d_compID_funct != 0 & compid_diff == 0
duplicates report PID year
keep PID year functions_o

save "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-1985_Functions_corrections.dta", replace


/*
replace functions = functions_o if d_compID_funct != 0 & compid_diff == 0  // replace functions with the original, if there were no corrections in the CID @UniLU -- > Excel-Errors?

drop nbrcompIDs* nbrfunctions* d_compID_funct*

gen nbrcompIDs = length(companiesid) - length(subinstr(companiesid, ",", "", .)) + 1
gen nbrfunctions = length(functions) - length(subinstr(functions, ",", "", .)) + 1
gen d_compID_funct = nbrcompIDs - nbrfunctions
tab d_compID_funct 
tab nbrcompIDs if d_compID_funct != 0
*/












