version 17
clear all
set more off
cap log close


global data "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\11_Directors_1994_2018\"

use "${data}Bisnode_Person-Firmen_Geo.dta", clear
keep if gremium == "Verwaltungsrat"
gen eintritt=date(eintrittdatum,"YMD")
gen eintrittjahr = year(eintritt)
gen austritt=date(austrittdatum,"YMD")
gen austrittjahr = year(austritt)
gen dob=date(geburtstag,"YMD")
gen yob=geburtstag if dob==.
destring yob, replace
gen geburtsjahr = year(dob)
replace geburtsjahr = yob if geburtsjahr == . & yob > 0

gen altereintritt = eintrittjahr - geburtsjahr
gen alteraustritt = austrittjahr - geburtsjahr

hist altereintritt, freq 
graph export "${data}freq_Alter_VREintritt.pdf", as(pdf) replace
hist alteraustritt, freq
graph export "${data}freq_Alter_VRAustritt.pdf", as(pdf) replace


