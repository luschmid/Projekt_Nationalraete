
******************************************************
* (0) Set directories
******************************************************

	clear all
	set more off


    capture cd "G:\004300-01\Verwaltungsräte"

	if _rc==0{ 
	global hauptpfad "G:\004300-01\Verwaltungsräte"
	global datapath "G:\004300-01\Verwaltungsräte"
	}


	
**********************************
* (A) Check Version 11, May 2018
**********************************

* (i) Read-in, rename, and label variables

import excel  "Check Sugarcube v11 Anita.xlsx"	, sheet("Tabelle1") clear cellrange(A3:U103)
save "Check Sugarcube v11.dta", replace
import excel  "Check Sugarcube v11 Gabriela.xlsx"	, sheet("Tabelle1") clear cellrange(A4:U100)
append using "Check Sugarcube v11.dta"
save "Check Sugarcube v11.dta", replace


rename A volume
rename B filename
rename C entry
rename D person_in_data
rename E person_in_personsfile
rename F lastname_correctfield
rename G firstname_correctfield
rename H address_correctfield
rename I lastname_correct
rename J firstname_correct
rename K address_correct
rename L no_directorships_book
rename M no_directorships_data
rename N selected_directorship
rename O directorship_in_companiesfile
rename P firmname_correctfield
rename Q firmaddress_correctfield
rename R firmcapital_correctfield
rename S firmname_correct
rename T firmaddress_correct
rename U firmcapital_correct

list if person_in_data==0

// recode all information as wrong if the information on higher level is wrong

local vars "lastname_correctfield firstname_correctfield address_correctfield lastname_correct firstname_correct address_correct no_directorships_data directorship_in_companiesfile firmname_correctfield firmaddress_correctfield firmcapital_correctfield firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if person_in_data==0 |  person_in_personsfile==0
}

local vars "lastname_correct firstname_correct address_correct"
foreach var of varlist `vars'{
replace `var'=0 if `var'field==0
}

local vars "firmname_correctfield firmaddress_correctfield firmcapital_correctfield firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if directorship_in_companiesfile==0 
}

local vars "firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if `var'field==0
}

gen share_recorded_directorships=no_directorships_data/no_directorships_book

label var volume 
label var filename
label var entry
label var person_in_data "Person recorded in persons or companies file"
label var person_in_personsfile "Person recorded in persons file"
label var lastname_correctfield "Last name in correct field"
label var firstname_correctfield "First name in correct field"
label var address_correctfield "Address in correct field"
label var lastname_correct "Last name correct"
label var firstname_correct "First name correct"
label var address_correct "Address name correct"
label var no_directorships_book "Number of directorships in book"
label var no_directorships_data "Number of directorships recorded"
label var selected_directorship 
label var directorship_in_companiesfile "Directorship in companies file"
label var firmname_correctfield "Firm name in correct field"
label var firmaddress_correctfield "Firm address in correct field"
label var firmcapital_correctfield "Firm capital in correct field"
label var firmname_correct "Firm name correct"
label var firmaddress_correct "Firm address correct"
label var firmcapital_correct "Firm capital correct"
label var share_recorded_directorships "Share of recorded directorships"


fsum person_in_data-address_correct  share_recorded_directorships directorship_in_companiesfile-firmcapital_correct, label format(%9.3f) 


***********************************************
* (B) Check Version 14: 1973-2003, August 2018
************************************************

* (i) Read-in, rename, and label variables

import excel  "Check Sugarcube v14 Niederberger_1973-2003.xlsx"	, sheet("Tabelle1") clear cellrange(A3:U154)
save "Check Sugarcube v14_1973-2003.dta", replace
import excel  "Check Sugarcube v14 Anita_1973-2003.xlsx"	, sheet("Tabelle1") clear cellrange(A4:U152)
append using "Check Sugarcube v14_1973-2003.dta"
save "Check Sugarcube v14_1973-2003.dta", replace


rename A volume
rename B filename
rename C entry
rename D person_in_data
rename E person_in_personsfile
rename F lastname_correctfield
rename G firstname_correctfield
rename H address_correctfield
rename I lastname_correct
rename J firstname_correct
rename K address_correct
rename L no_directorships_book
rename M no_directorships_data
rename N selected_directorship
rename O directorship_in_companiesfile
rename P firmname_correctfield
rename Q firmaddress_correctfield
rename R firmcapital_correctfield
rename S firmname_correct
rename T firmaddress_correct
rename U firmcapital_correct


// recode all information to missing if person is not in data b/c these are observations from the appendices (checked by SL and LS August 2018)

foreach var of varlist person_in_personsfile-firmcapital_correct person_in_data {
replace `var' = . if person_in_data == 0
}

// recode all information as wrong if the information on higher level is wrong

local vars "lastname_correctfield firstname_correctfield address_correctfield lastname_correct firstname_correct address_correct no_directorships_data directorship_in_companiesfile firmname_correctfield firmaddress_correctfield firmcapital_correctfield firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if person_in_data==0 |  person_in_personsfile==0
}

local vars "lastname_correct firstname_correct address_correct"
foreach var of varlist `vars'{
replace `var'=0 if `var'field==0
}

local vars "firmname_correctfield firmaddress_correctfield firmcapital_correctfield firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if directorship_in_companiesfile==0 
}

local vars "firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if `var'field==0
}

gen share_recorded_directorships=no_directorships_data/no_directorships_book


label var volume 
label var filename
label var entry
label var person_in_data "Person recorded in persons or companies file"
label var person_in_personsfile "Person recorded in persons file"
label var lastname_correctfield "Last name in correct field"
label var firstname_correctfield "First name in correct field"
label var address_correctfield "Address in correct field"
label var lastname_correct "Last name correct"
label var firstname_correct "First name correct"
label var address_correct "Address name correct"
label var no_directorships_book "Number of directorships in book"
label var no_directorships_data "Number of directorships recorded"
label var selected_directorship 
label var directorship_in_companiesfile "Directorship in companies file"
label var firmname_correctfield "Firm name in correct field"
label var firmaddress_correctfield "Firm address in correct field"
label var firmcapital_correctfield "Firm capital in correct field"
label var firmname_correct "Firm name correct"
label var firmaddress_correct "Firm address correct"
label var firmcapital_correct "Firm capital correct"
label var share_recorded_directorships "Share of recorded directorships"


fsum person_in_data-address_correct  share_recorded_directorships directorship_in_companiesfile-firmcapital_correct, label  format(%9.3f)




***********************************************
* (C) Check Version 14: 1934-1973, August 2018
************************************************

* (i) Read-in, rename, and label variables

import excel  "Check Sugarcube v14 Anita_1934-1970.xlsx"	, sheet("Tabelle1") clear cellrange(A3:U132)
save "Check Sugarcube v14_1934-1970.dta", replace
import excel  "Check Sugarcube v14 Gabriela_1934-1970.xlsx"	, sheet("Tabelle1") clear cellrange(A4:U143)
append using "Check Sugarcube v14_1934-1970.dta"
save "Check Sugarcube v14_1973-2003.dta", replace



rename A volume
rename B filename
rename C entry
rename D person_in_data
rename E person_in_personsfile
rename F lastname_correctfield
rename G firstname_correctfield
rename H address_correctfield
rename I lastname_correct
rename J firstname_correct
rename K address_correct
rename L no_directorships_book
rename M no_directorships_data
rename N selected_directorship
rename O directorship_in_companiesfile
rename P firmname_correctfield
rename Q firmaddress_correctfield
rename R firmcapital_correctfield
rename S firmname_correct
rename T firmaddress_correct
rename U firmcapital_correct

// recode all information to missing if person is not in data b/c these are observations from the appendices (checked by SL and LS August 2018)

foreach var of varlist person_in_personsfile-firmcapital_correct person_in_data {
replace `var' = . if person_in_data == 0
}

// recode all information as wrong if the information on higher level is wrong

local vars "lastname_correctfield firstname_correctfield address_correctfield lastname_correct firstname_correct address_correct no_directorships_data directorship_in_companiesfile firmname_correctfield firmaddress_correctfield firmcapital_correctfield firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if person_in_data==0 |  person_in_personsfile==0
}

local vars "lastname_correct firstname_correct address_correct"
foreach var of varlist `vars'{
replace `var'=0 if `var'field==0
}

local vars "firmname_correctfield firmaddress_correctfield firmcapital_correctfield firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if directorship_in_companiesfile==0 
}

local vars "firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if `var'field==0
}

gen share_recorded_directorships=no_directorships_data/no_directorships_book


label var volume 
label var filename
label var entry
label var person_in_data "Person recorded in persons or companies file"
label var person_in_personsfile "Person recorded in persons file"
label var lastname_correctfield "Last name in correct field"
label var firstname_correctfield "First name in correct field"
label var address_correctfield "Address in correct field"
label var lastname_correct "Last name correct"
label var firstname_correct "First name correct"
label var address_correct "Address name correct"
label var no_directorships_book "Number of directorships in book"
label var no_directorships_data "Number of directorships recorded"
label var selected_directorship 
label var directorship_in_companiesfile "Directorship in companies file"
label var firmname_correctfield "Firm name in correct field"
label var firmaddress_correctfield "Firm address in correct field"
label var firmcapital_correctfield "Firm capital in correct field"
label var firmname_correct "Firm name correct"
label var firmaddress_correct "Firm address correct"
label var firmcapital_correct "Firm capital correct"
label var share_recorded_directorships "Share of recorded directorships"


fsum person_in_data-address_correct  share_recorded_directorships directorship_in_companiesfile-firmcapital_correct, label  format(%9.3f)




***********************************************
* (D) Check Version 15: 1934-1973, December 2018
************************************************

* (i) Read-in, rename, and label variables

import excel  "Check Sugarcube v15 Anita_1934-1970.xlsx"	, sheet("Tabelle1") clear cellrange(A3:U132)
save "Check Sugarcube v15_1934-1970.dta", replace
import excel  "Check Sugarcube v15 Gabriela_1934-1970.xlsx"	, sheet("Tabelle1") clear cellrange(A3:U135)
append using "Check Sugarcube v14_1934-1970.dta"
save "Check Sugarcube v15_1973-2003.dta", replace



rename A volume
rename B filename
rename C entry
rename D person_in_data
rename E person_in_personsfile
rename F lastname_correctfield
rename G firstname_correctfield
rename H address_correctfield
rename I lastname_correct
rename J firstname_correct
rename K address_correct
rename L no_directorships_book
rename M no_directorships_data
rename N selected_directorship
rename O directorship_in_companiesfile
rename P firmname_correctfield
rename Q firmaddress_correctfield
rename R firmcapital_correctfield
rename S firmname_correct
rename T firmaddress_correct
rename U firmcapital_correct

// recode all information to missing if person is not in data b/c these are observations from the appendices (checked by SL and LS August 2018)

foreach var of varlist person_in_personsfile-firmcapital_correct person_in_data {
replace `var' = . if person_in_data == 0
}

// recode all information as wrong if the information on higher level is wrong

local vars "lastname_correctfield firstname_correctfield address_correctfield lastname_correct firstname_correct address_correct no_directorships_data directorship_in_companiesfile firmname_correctfield firmaddress_correctfield firmcapital_correctfield firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if person_in_data==0 |  person_in_personsfile==0
}

local vars "lastname_correct firstname_correct address_correct"
foreach var of varlist `vars'{
replace `var'=0 if `var'field==0
}

local vars "firmname_correctfield firmaddress_correctfield firmcapital_correctfield firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if directorship_in_companiesfile==0 
}

local vars "firmname_correct firmaddress_correct firmcapital_correct"
foreach var of varlist `vars'{
replace `var'=0 if `var'field==0
}

gen share_recorded_directorships=no_directorships_data/no_directorships_book


label var volume 
label var filename
label var entry
label var person_in_data "Person recorded in persons or companies file"
label var person_in_personsfile "Person recorded in persons file"
label var lastname_correctfield "Last name in correct field"
label var firstname_correctfield "First name in correct field"
label var address_correctfield "Address in correct field"
label var lastname_correct "Last name correct"
label var firstname_correct "First name correct"
label var address_correct "Address name correct"
label var no_directorships_book "Number of directorships in book"
label var no_directorships_data "Number of directorships recorded"
label var selected_directorship 
label var directorship_in_companiesfile "Directorship in companies file"
label var firmname_correctfield "Firm name in correct field"
label var firmaddress_correctfield "Firm address in correct field"
label var firmcapital_correctfield "Firm capital in correct field"
label var firmname_correct "Firm name correct"
label var firmaddress_correct "Firm address correct"
label var firmcapital_correct "Firm capital correct"
label var share_recorded_directorships "Share of recorded directorships"


fsum person_in_data-address_correct  share_recorded_directorships directorship_in_companiesfile-firmcapital_correct, label  format(%9.3f)


