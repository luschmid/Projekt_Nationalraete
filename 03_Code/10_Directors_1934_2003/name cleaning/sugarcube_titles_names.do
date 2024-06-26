version 16
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"

*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"


************************************************************
*** Import, Geo-Code and Name Cleaning in Sugarcube Data ***
************************************************************

** import from original csv files (saved as UTF-8 with BOM)
foreach val in 1934 1943 1960 1962 1963 1964 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", delimiters(";") encoding("utf-8") clear
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta", replace
}
foreach val in 1966 1967 1970 1973 1976 1980_1981 1998_1999 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", delimiters(";") encoding("utf-8") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta", replace
}
foreach val in 1979_1980 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", encoding("utf-8") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta", replace
}
forv val = 1982(1)1998 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", delimiters(";") encoding("utf-8") clear
	replace year = year-1 // years of publication != years of survey
	tostring *, replace
	save "$path\02_Processed_data\10_Directors_1934_2003\pers_`val'_tmp.dta", replace
}
forv val = 1999(1)2003 {
	import delim using "$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_persons.csv", delimiters(";") encoding("utf-8") clear
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

*use "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp1.dta", clear
rename id PID
label var PID "Person ID per year"
destring PID, replace
destring year, replace
save "$path\02_Processed_data\10_Directors_1934_2003\pers_1934-2003_tmp2.dta", replace

* Collapse data for inspection
keep PID year titlejob firstname lastname

preserve
keep titlejob
gen count = -1
replace titlejob = ustrtrim(titlejob)
collapse (sum) count, by(titlejob)
sort titlejob
gen newtitlejob = ""
order count titlejob newtitlejob
export excel using "$path\02_Processed_data\10_Directors_1934_2003\01_Round1_Outfiles\titles_unique.xlsx", first(var) replace 
restore

preserve
keep firstname
gen count = -1
replace firstname = ustrtrim(firstname)
collapse (sum) count, by(firstname)
sort firstname
save "$path\02_Processed_data\10_Directors_1934_2003\01_Round1_Outfiles\firstname_unique_R0.dta", replace
restore

preserve
keep lastname
gen count = -1
replace lastname = ustrtrim(lastname)
collapse (sum) count, by(lastname)
sort lastname
save "$path\02_Processed_data\10_Directors_1934_2003\01_Round1_Outfiles\lastname_unique_R0.dta", replace
restore

* Prepare collapsed data for inspection
use "$path\02_Processed_data\10_Directors_1934_2003\01_Round1_Outfiles\firstname_unique_R0.dta", clear
sort firstname
gen first = usubstr(firstname, 1, 1)
replace first = ustrupper(first)
replace first = "A" if first == "Ä"
replace first = "A" if first == "À"
replace first = "A" if first == "Â"
replace first = "E" if first == "È"
replace first = "E" if first == "É"
replace first = "E" if first == "Ê"
replace first = "O" if first == "Ö"
replace first = "O" if first == "Ô"
replace first = "O" if first == "Ò"
replace first = "O" if first == "Ó"
replace first = "U" if first == "Ü"
replace first = "U" if first == "Û"
replace first = "U" if first == "Ù"
replace first = "I" if first == "Ï"
replace first = "I" if first == "Î"
replace first = "I" if first == "Ì"
replace first = "C" if first == "Ç"
replace first = "S" if first == "SS"
gen help = 0 
replace help = 1 if regexm(first, "[A-Z]")
replace first = "@" if help == 0
drop help

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z {
	preserve
	keep if first == "`i'"
	drop first
	gen newfirstname = ""
	order count firstname newfirstname
	export excel using "$path\02_Processed_data\10_Directors_1934_2003\01_Round1_Outfiles\firstname_unique`i'.xlsx", first(var) replace 
	restore
}

use "$path\02_Processed_data\10_Directors_1934_2003\01_Round1_Outfiles\lastname_unique_R0.dta", clear
sort lastname
gen first = usubstr(lastname, 1, 1)
replace first = ustrupper(first)
replace first = "A" if first == "Ä"
replace first = "A" if first == "À"
replace first = "A" if first == "Â"
replace first = "E" if first == "È"
replace first = "E" if first == "É"
replace first = "E" if first == "Ê"
replace first = "O" if first == "Ö"
replace first = "O" if first == "Ô"
replace first = "O" if first == "Ò"
replace first = "O" if first == "Ó"
replace first = "U" if first == "Ü"
replace first = "U" if first == "Û"
replace first = "U" if first == "Ù"
replace first = "I" if first == "Ï"
replace first = "I" if first == "Î"
replace first = "I" if first == "Ì"
replace first = "C" if first == "Ç"
replace first = "S" if first == "SS"
gen help = 0 
replace help = 1 if regexm(first, "[A-Z]")
replace first = "@" if help == 0
drop help

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z {
	preserve
	keep if first == "`i'"
	drop first
	gen newlastname = ""
	order count lastname newlastname
	export excel using "$path\02_Processed_data\10_Directors_1934_2003\01_Round1_Outfiles\lastname_unique`i'.xlsx", first(var) replace 
	restore
}

/*
* Prepare collapsed data for inspection
use "$path\02_Processed_data\10_Directors_1934_2003\firstname_unique.dta", clear
sort firstname
gen n = _n
gen N = _N
gen NoFiles = ceil(N/20000)
local NoFiles = NoFiles
gen file = .
gen upper = .
gen lower = .
forv i = 1(1)`NoFiles' {
	replace lower = (`i'-1) * 20000
	replace upper = `i' * 20000
	replace file = `i' if n >= lower & n < upper
	preserve
	keep if file == `i'
	drop upper lower NoFiles n N file
	gen newfirstname = ""
	order count firstname newfirstname
	export excel using "$path\02_Processed_data\10_Directors_1934_2003\firstname_unique`i'.xlsx", first(var) replace 
	restore
}

use "$path\02_Processed_data\10_Directors_1934_2003\lastname_unique.dta", clear
sort lastname
gen n = _n
gen N = _N
gen NoFiles = ceil(N/20000)
local NoFiles = NoFiles
gen file = .
gen upper = .
gen lower = .
forv i = 1(1)`NoFiles' {
	replace lower = (`i'-1) * 20000
	replace upper = `i' * 20000
	replace file = `i' if n >= lower & n < upper
	preserve
	keep if file == `i'
	drop upper lower NoFiles n N file
	gen newlastname = ""
	order count lastname newlastname
	export excel using "$path\02_Processed_data\10_Directors_1934_2003\lastname_unique`i'.xlsx", first(var) replace 
	restore
}
