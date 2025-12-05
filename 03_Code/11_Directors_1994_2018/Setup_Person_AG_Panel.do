clear
set more off
version 17

*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "C:\Current\Dropbox\Projekt Nationalräte"
global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

** Use Person-Firm data: drop non-AGs, build panel structure
use personenid anrede vorname nachname wohnort gdenr_2018_w nat_ch foreignWO duns handelsregisternummer rechtsform gremium funktion eintrittdatum austrittdatum kapitalnominalag kapitaleinbezahltag plzdomiziladresse firma using ///
"$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Person-Firmen_Geo.dta", clear

* keep only AG's and VR mandates (post-processing only considered those)
keep if gremium == "Verwaltungsrat" & rechtsform == "Aktiengesellschaft"
drop if inlist(funktion, "Aktuar/in (nicht Mitglied)", "Sekretär/in (nicht Mitglied)", "Beisitzer/in", "Ausseramtliche/r Konkursverwalter/in", "Generalsekretär/in (nicht Mitglied)", "Liquidator/in", "Kassier/in (nicht Mitglied)", "Protokollführer/in (nicht Mitglied)", "Verwalter/in (nicht Mitglied)")

duplicates report personenid duns handelsregisternummer rechtsform gremium funktion eintrittdatum austrittdatum kapitalnominalag kapitaleinbezahltag plzdomiziladresse firma
duplicates tag personenid duns handelsregisternummer rechtsform gremium funktion eintrittdatum austrittdatum kapitalnominalag kapitaleinbezahltag plzdomiziladresse firma, gen(dupl) // visual inspection  
*br * if dupl>0   // visual inspection suggests genuine duplicates
drop dupl
duplicates drop personenid duns handelsregisternummer rechtsform gremium funktion eintrittdatum austrittdatum kapitalnominalag kapitaleinbezahltag plzdomiziladresse firma, force

drop handelsregisternummer rechtsform gremium kapitalnominalag kapitaleinbezahltag plzdomiziladresse firma

bysort personenid duns funktion: gen mandid = _n
expand 25, gen(expy)
bysort personenid duns funktion mandid: gen year=_n + 1993
drop expy
sort personenid duns funktion mandid year
order year personenid duns funktion mandid

gen entrydate = date(eintrittdatum, "YMD")
format entrydate %td
gen exitdate = date(austrittdatum, "YMD")
format exitdate %td
gen flag = 1 if exitdate<entrydate & entrydate != . & exitdate != .
gen exitdate_tmp = exitdate
replace exitdate = entrydate if flag == 1 // assumption: confused entry/exit date
replace entrydate = exitdate_tmp if flag == 1
drop exitdate_tmp flag
gen entryear=year(entrydate)
gen exityear=year(exitdate)
tab entryear
tab exityear
replace entryear = 1993 if eintrittdatum == ""
replace exityear = 2019 if austrittdatum == ""
tab entryear
tab exityear
drop if year<entryear | year>exityear

duplicates report personenid duns funktion year 
duplicates tag personenid duns funktion year, gen(dupl)
*br * if dupl>0
duplicates drop personenid duns funktion year, force  // drop "true" duplicates ("funktion" included)
drop dupl

* Overlapping spells as President and Member --- > keep highest function (e.g., President)
bysort personenid duns year: gen no_duplMand = _n  
bysort personenid duns year: egen avg_duplMand = mean(no_duplMand)
gen duplMand = 0
bysort personenid duns year: replace duplMand = 1 if avg_duplMand>1
drop avg_duplMand

gen MandH = 9
replace MandH = 1 if funktion=="Präsident/in"
replace MandH = 2 if funktion=="Vizepräsident/in"
replace MandH = 3 if funktion=="Delegierte/r"
replace MandH = 4 if funktion=="Mitglied"
replace MandH = 5 if funktion=="Einziges Mitglied"
replace MandH = 6 if funktion=="Verwalter/in"
replace MandH = 7 if funktion=="Aktuar/in"
replace MandH = 8 if funktion=="Sekretär/in"
bysort personenid duns year: egen minMandH = min(MandH)
keep if MandH == minMandH  // keep the "highest" function
duplicates tag personenid duns year, gen(dupl)
*br * if dupl>0
duplicates drop personenid duns year, force
drop duplMand no_duplMand MandH minMandH dupl

keep if year>1993 & year<2018
sort personenid duns year

tab year

gen sex_ambig = 0
replace sex_ambig = 1 if anrede == ""
tab sex_ambig

preserve
drop if anrede != ""
bysort vorname: gen n = -1
collapse (sum) n, by(vorname)
sort n
save "$path\02_Processed_data\11_Directors_1994_2018\vornamen.dta", replace
restore

gen no_geo = 0
replace no_geo = 1 if gdenr_2018_w == .
replace no_geo = 3 if foreignWO == 1
tab no_geo


erase "$path\02_Processed_data\11_Directors_1994_2018\vornamen.dta"
