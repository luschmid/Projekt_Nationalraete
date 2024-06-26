clear
cap log close
set more 1
version 17

*** Paths
global HZ "E:\12. Cloud\Dropbox\Projekt Nationalräte\05_Texts_and_presns\04_Knappheit_2023\"


****** Analysis
import excel "$HZ\closeness_2023.xlsx", first sheet("closeness_2023") clear

gen target = 0
replace target = 1 if firstname == "Zeno" & lastname == "Staub"
replace target = 1 if firstname == "Donato" & lastname == "Scognamiglio"
replace target = 1 if firstname == "Stefan" & lastname == "Brupbacher"
replace target = 1 if firstname == "Christof" & lastname == "Züger"
replace target = 1 if firstname == "Hans-Ulrich" & lastname == "Bigler"
replace target = 1 if firstname == "Michel" & lastname == "Matter"
replace target = 1 if firstname == "Martin" & lastname == "Rufer"
replace target = 1 if firstname == "Matthias Samuel" & lastname == "Jauslin"

bysort canton_num list_num: egen targetlist = max(target)
tab targetlist

bysort canton_num list_num: gen perslist = _N
bysort canton_num list_num (votes): gen lrank = (_n)*(-1)
bysort canton_num list_num: gen listperform = perslist + lrank + 1
drop perslist lrank
bysort canton_num list_num: egen seats = total(elected)
bysort canton_num list_num: gen ersatz = listperform - seats
replace ersatz = . if ersatz <= 0 
replace ersatz = . if seats == 0

keep if targetlist == 1

gen listnr = list_num
replace listnr = "3" if list_num == "3a"
replace listnr = "4" if list_num == "4a"
destring listnr, replace
gen mvotes = votes * (-1)
sort canton_num listnr mvotes
drop mvotes listnr
save "$HZ\closeness_2023_Valda.dta", replace