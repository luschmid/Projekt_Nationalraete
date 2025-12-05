global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\1 Data\1_data\14_jarckspuehlerwidmer\Cantonal Information"

import excel using "$path\Zürich.xlsx", sheet(Parteien) first clear
keep ID_PERSON_NEW PARTEIBEZEICHNUNG
rename PARTEIBEZEICHNUNG party
bysort ID_PERSON_NEW: gen tot = _N
drop if party == "" & tot > 1
bysort ID_PERSON_NEW: gen no = _n
reshape wide party, i(ID_PERSON_NEW) j(no)
gen party = party1
replace party = party + ", " + party2 if party2 != ""
replace party = party + ", " + party3 if party3 != ""
replace party = party + ", " + party4 if party4 != ""
replace party = party + ", " + party5 if party5 != ""
replace party = party + ", " + party6 if party6 != ""
replace party = party + ", " + party7 if party7 != ""
keep ID_PERSON_NEW party
save "$path\Zürich.dta", replace

import excel using "$path\Zürich.xlsx", sheet(Bürgerorte) first clear
keep ID_PERSON_NEW ORT
rename ORT origin
bysort ID_PERSON_NEW: gen tot = _N
drop if origin == "" & tot > 1
bysort ID_PERSON_NEW: gen no = _n
reshape wide origin, i(ID_PERSON_NEW) j(no)
gen origin = origin1
replace origin = origin + ", " + origin2 if origin2 != ""
replace origin = origin + ", " + origin3 if origin3 != ""
replace origin = origin + ", " + origin4 if origin4 != ""
merge 1:1 ID_PERSON_NEW using "$path\Zürich.dta"
keep ID_PERSON_NEW party origin
save "$path\Zürich.dta", replace

import excel using "$path\Zürich.xlsx", sheet(Wohnorte) first clear
keep ID_PERSON_NEW GEMEINDENAME
rename GEMEINDENAME municipality
bysort ID_PERSON_NEW: gen tot = _N
drop if municipality == "" & tot > 1
bysort ID_PERSON_NEW: gen no = _n
reshape wide municipality, i(ID_PERSON_NEW) j(no)
gen municipality = municipality1
replace municipality = municipality + ", " + municipality2 if municipality2 != ""
replace municipality = municipality + ", " + municipality3 if municipality3 != ""
replace municipality = municipality + ", " + municipality4 if municipality4 != ""
replace municipality = municipality + ", " + municipality5 if municipality5 != ""
replace municipality = municipality + ", " + municipality6 if municipality6 != ""
merge 1:1 ID_PERSON_NEW using "$path\Zürich.dta"
keep ID_PERSON_NEW party origin municipality
save "$path\Zürich.dta", replace

import excel using "$path\Zürich.xlsx", sheet(Personen) first clear
merge 1:1 ID_PERSON_NEW using "$path\Zürich.dta"
order ID_PERSON_NEW NACHNAME VORNAME GESCHLECHT BERUF AKAD_TITEL party ///
    municipality origin DATUM_GEBURT_TAG DATUM_GEBURT_MONAT DATUM_GEBURT_JAHR ///
    DATUM_TOD_TAG DATUM_TOD_MONAT DATUM_TOD_JAHR
keep ID_PERSON_NEW-DATUM_TOD_JAHR

rename ID_PERSON_NEW IDZH
rename NACHNAME name
rename VORNAME firstname
rename GESCHLECHT sex 
rename BERUF job
rename AKAD_TITEL title
rename DATUM_GEBURT_TAG bday
rename DATUM_GEBURT_MONAT bmonth
rename DATUM_GEBURT_JAHR byear
rename DATUM_TOD_TAG dday
rename DATUM_TOD_MONAT dmonth
rename DATUM_TOD_JAHR dyear
save "$path\Zürich.dta", replace

export excel using "$path\Zürich_processed.xlsx", firstrow(variables)



