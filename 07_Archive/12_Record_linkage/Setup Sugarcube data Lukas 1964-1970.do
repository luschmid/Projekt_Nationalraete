******************************************************
* (A) Set directories
******************************************************

clear
set more off

*global path "E:\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "C:\Dropbox\Projekt Nationalräte"
*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"
*global path "C:\Users\Lukas\Dropbox\Projekt Nationalräte"

******************************************************
* (A) MANUAL CORRECTIONS
******************************************************
/*local years "1934 1943 1960 1962 1963 1964 1966 1967 1970 1973 1976 1979_1980 1980_1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1998_1999 2000 2001 2002 2003"
foreach yr in `years'{
import delim ///
	"$path\01_Raw_data\10_Directors_1934_2003\Kern_`yr'_persons.csv", ///
	 clear delim(";")  encoding("utf-8") 

duplicates tag id, gen(dups)
keep if dups>0 | firstname=="" | lastname==""
drop dups 
capture tostring address
capture tostring address2 
capture tostring alternateaddress 
capture tostring postcode 
if "`yr'"=="1934" save "$path\02_Processed_data\10_Directors_1934_2003\temp.dta", replace
else {
append using "$path\02_Processed_data\10_Directors_1934_2003\temp.dta", force
save "$path\02_Processed_data\10_Directors_1934_2003\temp.dta", replace
}  
}*/



local vals "1964 1966 1967 1970"
foreach val in `vals'{

* (a) import company file 

import delim ///
	"$path\01_Raw_data\10_Directors_1934_2003\Kern_`val'_companies.csv", ///
	 clear delim(";")  encoding("utf-8") 

gen company=name+"("+city+")"
rename id  companiesid
sum year 
local yr=r(max)
rename year year_comp
keep companiesid company year_comp
save "$path\02_Processed_data\10_Directors_1934_2003\Lukas\temp_companies.dta", replace

* (b) import persons file 


use "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person-Geo.dta", clear
keep if year==(`yr'-1)

*drop PID
gen id_new=_n

order id year firstname lastname city titlejob address address2 alternateaddress
capture tostring address2, replace 
capture tostring alternateaddress, replace 
split companiesid, parse(",")
drop companiesid

preserve 
keep companiesid*
describe 
local companiesid_max = c(k)
di `companiesid_max'
restore

forvalues i = 1(1)`companiesid_max'{
preserve
local iplus1 =`i' +1
if `i'<`companiesid_max'{
keep if companiesid`i'!="" & companiesid`iplus1'=="" 
drop companiesid`iplus1'-companiesid`companiesid_max'
}
else{
keep if companiesid`i'!="" 
}
reshape long companiesid, i(id_new) j(number)
drop if companiesid==""
destring companiesid, replace
if `i'==1 save "$path\02_Processed_data\10_Directors_1934_2003\Lukas\temp.dta", replace 
else{
append using "$path\02_Processed_data\10_Directors_1934_2003\Lukas\temp.dta"
save "$path\02_Processed_data\10_Directors_1934_2003\Lukas\temp.dta", replace 
}
restore
}



* (c) merge persons and company file 

use "$path\02_Processed_data\10_Directors_1934_2003\Lukas\temp.dta", clear
merge m:1 companiesid using  "$path\02_Processed_data\10_Directors_1934_2003\Lukas\temp_companies.dta"

preserve
keep if _merge<3
if `yr'==1964 save "$path\02_Processed_data\10_Directors_1934_2003\Lukas\pers_comp_notmerged.dta", replace
else {
append using "$path\02_Processed_data\10_Directors_1934_2003\Lukas\pers_comp_notmerged.dta", force
save "$path\02_Processed_data\10_Directors_1934_2003\Lukas\pers_comp_notmerged.dta", replace
}
restore

drop if _merge<3

drop _merge
drop companiesid
drop number
duplicates drop id_new company, force
sort id_new company // sort companies alphabetically 
bysort id_new: gen number=_n
sum number 
local number_max=r(max)
reshape wide company, i(id_new) j(number)

local vars "company"
foreach var in `vars'{
local counter=1
capture drop `var'
forvalues i=1(1)`number_max'{
	//di "`var'"`i'
	capture confirm variable `var'`i'
		if !_rc {
			if `counter'==1{
			gen `var'=`var'`i'
			}
			else{
			qui replace `var'=`var'+", "+`var'`i' if `var'`i'!=""
			}
		}
		else {
		}
	local counter=`counter'+1
}
}

capture tostring address2 
capture tostring alternateaddress 

keep PID id_new year firstname lastname city titlejob address address2 alternateaddress company
order PID id_new year firstname lastname city titlejob address address2 alternateaddress company
rename company companies

if `yr'==1964 save "$path\02_Processed_data\10_Directors_1934_2003\Lukas\temp_allyears.dta", replace
else {
append using "$path\02_Processed_data\10_Directors_1934_2003\Lukas\temp_allyears.dta", force
save "$path\02_Processed_data\10_Directors_1934_2003\Lukas\temp_allyears.dta", replace
} 
}
** export
keep PID year companies
export delimited "$path\02_Processed_data\10_Directors_1934_2003\Lukas\Sugarcube_Companies_Lukas.csv", delimiter(";") replace


/* 

duplicates tag year firstname lastname city titlejob address address2 alternateaddress, gen(dups)

preserve
keep if dups==0
save "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_1.dta", replace
restore

preserve
keep if dups>0
save "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_2.dta", replace
restore

* (d) merge constructed file to initial file


use "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person-Geo.dta", clear
duplicates tag year firstname lastname city titlejob address address2 alternateaddress, gen(dups)

preserve
keep if dups>0
save "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person-Geo_temp.dta", replace
restore


keep if dups==0
merge 1:1 year firstname lastname city titlejob address address2 alternateaddress ///
	 using "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_1.dta", gen(merge_all1)

append using "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person-Geo_temp.dta"
append using "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_2.dta"

drop geo_merge countPLZ4 gdename gdenr_2018 GdeNr_E_CNTR GdeNr_N_CNTR nCH gde_options dups merge_all
 
save "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person-VRMandates.dta", replace


erase "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears.dta" 
erase "$path\02_Processed_data\10_Directors_1934_2003\temp_companies.dta" 
erase "$path\02_Processed_data\10_Directors_1934_2003\temp.dta" 
erase "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_1.dta" 
erase "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_2.dta" 
erase "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person-Geo_temp.dta" 


