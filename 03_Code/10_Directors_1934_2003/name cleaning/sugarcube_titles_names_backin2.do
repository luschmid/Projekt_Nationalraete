version 16
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"

*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"

*****************************************************************
*** Import corrected first and last name files (second round) ***
*****************************************************************

*** First names: R2 *** 

gen count = .
gen firstname = ""
gen newfirstname = ""
gen newlastname = ""

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z { 
	import excel using "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\firstname_unique`i'_c2_out.xlsx", clear firstrow
	save "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\firstname_unique`i'_R2.dta", replace
}
use "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\firstname_unique@_R2.dta", clear
foreach i in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z {
	append using "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\firstname_unique`i'_R2.dta", force
	save "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\firstname_unique_R2.dta", replace
}
rename firstname firstnameR2
rename newfirstname newfirstnameR2
rename newlastname newlastnameR2
rename errorcode indfirstR2
rename comment commentR2
rename serie serieR2
keep fid firstnameR2 newfirstnameR2 newlastnameR2 indfirstR2 serieR2 commentR2 
sort fid
save "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\firstname_unique_R2.dta", replace

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z { 
	erase "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\firstname_unique`i'_R2.dta"
}

** merge with first round **
use "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\firstname_unique_R1.dta"
sort fid
merge m:1 fid using "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\firstname_unique_R2.dta"
drop if firstname == ""

* checks
gen seriecheck = 1 if serie != serieR2 & _merge==3 // _merge == 3: just for subsample of second round checks
tab seriecheck
drop seriecheck serieR2
gen firstnamecheck = 1 if firstname != firstnameR2 & _merge==3 // _merge == 3: just for subsample of second round checks
tab firstnamecheck
// error is in firstnameR2 -- > coder renames it form "Dixon Quackenbuch" to "Dixon"
drop firstnamecheck firstnameR2

** indicators for newfirstname issues
gen indfirst1 = .
replace indfirst1 = 2 if newlastname != "" & newfirstname == ""  // lastname in firstname column
replace indfirst1 = 3 if newfirstname != ""  // corrections of only the first name (--> standard case)
replace indfirst1 = 4 if newfirstname != "" & newlastname != ""  // there is first and last name info in the firstname variable
replace indfirst1 = 1 if newfirstname == "0" | newfirstname == "1" | newfirstname == "?" | newfirstname == "Ortsfragment" | newfirstname == "ortsfragment" | newfirstname == "Namensfragment" | newfirstname == "namensfragment" | newfirstname == "Fragment" | newfirstname == "fragment" | newfirstname == "Firmeninformation" | newfirstname == "firmeninformation" | newfirstname == "Firmenbezeichnung" | newfirstname == "firmenbezeichnung" | newfirstname == "Firstname fragment" | newfirstname == "firstname fragment" | newfirstname == "Lastname fragment" | newfirstname == "lastname fragment"  | newlastname == "Lastname fragment" | newlastname == "lastname fragment" | newlastname == "lastnamefragment" // flagged cases by coders
tab indfirst indfirst1  // display changes from R2

replace indfirst = indfirst1  // replace the first version of the indicator (changed due to R2 corrections)
drop indfirst1

replace newfirstname = newfirstnameR2 if _merge==3
replace newlastname = newlastnameR2 if _merge==3
rename newfirstname newfirstname_f
rename newlastname newlastname_f


// check the errorcode in second round flag
*br * if indfirstR2 != .  // no action required
// check cases with commen in second round
*br * if commentR2 != ""  // no action required

// duplicates
duplicates report firstname
duplicates tag firstname, gen(dup)
*br * if dup == 1 // checked manually
drop if newfirstname_f == "" & dup == 1
drop dup
duplicates report firstname

// how many fragment cases? -- > indfirst == 1
egen helpcount = total(count) if indfirst==1
tab helpcount // 2528 cases of fragments to check in raw data (w/o adress fragments)

drop count newfirstnameR2 newlastnameR2 indfirstR2 commentR2 _merge fid serie helpcount

sort firstname
save "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\firstname_unique.dta", replace



*** Last names: R2 *** 

clear 
gen count = .
gen firstname = ""
gen newfirstname = ""
gen newlastname = ""

save "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique_R2.dta", replace

foreach i in @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z { 
	import excel using "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\lastname_unique`i'_c2_out.xlsx", clear firstrow
	save "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\lastname_unique`i'_R2.dta", replace
}
use "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\lastname_unique@_R2.dta", clear
foreach i in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z {  
	append using "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\lastname_unique`i'_R2.dta", force
	save "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\lastname_unique_R1.dta", replace
}
rename lastname lastnameR2
rename newfirstname newfirstnameR2
rename newlastname newlastnameR2
rename comment commentR2
rename serie serieR2
gen indlastR2 = 1 if commentR2 == "Lastname fragment" | commentR2 == "Fragment" | commentR2 == "Fragment "
*br * if errorcode != . // nothing notable

keep lid lastnameR2 newfirstnameR2 newlastnameR2 indlastR2 serieR2 
sort lid
save "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\lastname_unique_R2.dta", replace

foreach i in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z {
	erase "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\lastname_unique`i'_R2.dta"
}

** merge with first round **
use "$path\02_Processed_data\10_Directors_1934_2003\02_Round1_Infiles\lastname_unique_R1.dta"
sort lid
merge m:1 lid using "$path\02_Processed_data\10_Directors_1934_2003\04_Round2_Infiles\lastname_unique_R2.dta"
drop if lastname == ""

* checks
gen seriecheck = 1 if serie != serieR2 & _merge==3 // _merge == 3: just for subsample of second round checks
tab seriecheck
drop seriecheck serieR2
gen lastnamecheck = 1 if lastname != lastnameR2 & _merge==3 // _merge == 3: just for subsample of second round checks
tab lastnamecheck
drop lastnamecheck lastnameR2

** indicators for newfirstname issues
gen indlast1 = .
replace indlast1 = 2 if newlastname == "" & newfirstname != ""  // firstname in lastname column
replace indlast1 = 3 if newlastname != ""  // corrections of only the last name (--> standard case)
replace indlast1 = 4 if newlastname != "" & newfirstname != ""  // there is first and last name info in the lastname variable
replace indlast1 = 1 if newlastname == "0" | newlastname == "1" | newlastname == "?" | newlastname == "Ortsfragment" | newlastname == "ortsfragment" | newlastname == "newlastname" | newlastname == "firmeninformation" | newlastname == "Firmeninformation" | newlastname == "Firmeninfortmation" | newlastname == "Firmenbezeichnung" | newlastname == "firmenbezeichnung" | newfirstname == "Firstname fragment" | newfirstname == "firstname fragment"  // flagged cases by coders

tab indlast indlast1
replace indlast = indlast1
drop indlast1

replace newfirstname = newfirstnameR2 if _merge==3
replace newlastname = newlastnameR2 if _merge==3
rename newfirstname newfirstname_l
rename newlastname newlastname_l
replace indlast = indlastR2 if indlastR2 != .

// duplicates
duplicates report lastname
duplicates tag lastname, gen(dup)
bysort lastname newlastname_l: gen helpdup = _n
drop if helpdup > 1 & dup == 1
drop dup helpdup
duplicates tag lastname, gen(dup)
drop if newlastname_l == "" & dup == 1
drop dup 
duplicates report lastname

// how many fragment cases? -- > indlast == 1
egen helpcount = total(count) if indlast==1
tab helpcount // 2028 cases of fragments to check in raw data (w/o adress fragments)

drop count newfirstnameR2 newlastnameR2 indlastR2 _merge lid serie helpcount

sort lastname
save "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\lastname_unique.dta", replace



*** Merge back into the original raw data

use "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person-Geo.dta", clear

gen tempID = _n
gen firstname_orig = firstname
gen lastname_orig = lastname 

replace firstname = ustrtrim(firstname)
replace lastname = ustrtrim(lastname)

sort firstname

merge m:1 firstname using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\firstname_unique.dta",  gen(mergefirst)
*br * if mergefirst==2
replace newfirstname_f = "Paul Alois" if firstname == "Paul AloiB541s"
drop if mergefirst == 2

sort lastname
merge m:1 lastname using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\lastname_unique.dta",  gen(mergelast)
*br * if mergelast==2
drop if mergelast == 2

* checks

/*
Cases to check
indfirst: "Indicator firstname issue: 1(fragment), 2(lastname in first), 3(corr firstname only - std), 4(corr first&last)"
indlast: "Indicator lastname issue: 1(fragment), 2(firstname in last), 3(corr lastname only - std), 4(corr first&last)"

-- > Change cases w/o further checks: "3"
-- > Check cases: "1", "2", "4"

*/

tab indfirst
tab indlast 

* corrections

replace firstname = newfirstname_f if newfirstname_f != "" & newlastname_f == "" & newfirstname_l == ""  // no lastnames found in firstname files
replace lastname = newlastname_l if newlastname_l != "" & newfirstname_l == "" & newlastname_f == ""   // no lastnames found in firstname files

replace firstname = newfirstname_l if firstname == "" // replace empty original firstnames where we find first names in lastname files
replace lastname = newlastname_l if lastname == ""    // replace empty original lastnames where we find last names in firstname files

* sample to manually check
gen cases = 0 
replace cases = 1 if newfirstname_f != "" |  newlastname_f != "" | newfirstname_l != "" | newlastname_l != ""  // all cases with some correction
gen corrections = 0
replace corrections = 1 if newfirstname_f != "" & newlastname_f == "" & newfirstname_l == ""  // already corrected: no lastnames found in firstname files
replace corrections = 1 if newlastname_l != "" & newfirstname_l == "" & newlastname_f == ""  //  already corrected: no lastnames found in firstname files

gen manualcheck = 1 if cases == 1 & corrections == 0
tab manualcheck

preserve
keep if manualcheck == 1
keep firstname lastname tempID PID year ID_dupl
gen finalfirstname = ""
gen finallastname = ""
gen cases = _n
order cases tempID lastname firstname finallastname finalfirstname
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_1.xlsx" if cases<1500, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_2.xlsx" if cases>=1500 & cases<3000, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_3.xlsx" if cases>=3000 & cases<4500, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_4.xlsx" if cases>=4500 & cases<6000, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_5.xlsx" if cases>=6000 & cases<7500, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_6.xlsx" if cases>=7500 & cases<9000, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_7.xlsx" if cases>=9000 & cases<10500, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_8.xlsx" if cases>=10500 & cases<12000, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_9.xlsx" if cases>=12000 & cases<13500, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_10.xlsx" if cases>=13500 & cases<15000, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_11.xlsx" if cases>=15000 & cases<16500, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_12.xlsx" if cases>=16500 & cases<18000, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_13.xlsx" if cases>=18000 & cases<19500, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_14.xlsx" if cases>=19500 & cases<21000, first(var) replace
export excel using "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Original_check_15.xlsx" if cases>=21000 & cases<22600, first(var) replace
restore

/*
save "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\05_Round3_Outfiles\Sugarcube_Person_RevName-Geo.dta", replace
*/














