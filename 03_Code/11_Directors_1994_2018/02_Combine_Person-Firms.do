version 17
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"

**** Import Bisnode firm data
import delim using "$path\01_Raw_data\11_Directors_1994_2018\Uni Fribourg Firmendaten 20-07-2018.csv", delimiters(";") varnames(1) clear
duplicates report duns
rename gründungsjahr gruendungsjahr
rename aktivitätsstatus aktivitaetsstatus
tostring *, replace

save "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_raw.dta", replace

*** Correct errors in firm data
use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_raw.dta", clear
gen flag = 1 if v37!="" | v38!="" | v39!="" | v40!="" 
egen flagobs = total(flag)
di flagobs
drop flag flagobs

*** 1) errors to the left of nogabeschrieb (the variable that causes most problems)
tab nogacode
tab rechtsform
/*
CH |          1        0.00       99.98
Detailhandel mit Herrenbekleidung |          1        0.00       99.98
                          deutsch |        177        0.01      100.00
                       franzsisch |         19        0.00      100.00
                      italienisch |          2        0.00      100.00
*/
*br * if nogacode =="CH" | nogacode =="Detailhandel mit Herrenbekleidung" | nogacode =="deutsch" | nogacode =="franzsisch" | nogacode =="italienisch" 
keep if nogacode =="CH" | nogacode =="Detailhandel mit Herrenbekleidung" | nogacode =="deutsch" | nogacode =="französisch" | nogacode =="italienisch" | ///
rechtsform == "deutsch"


* a) Single error entry (subsample_1a)
preserve
keep if duns == "482467966"
gen codomiziladresse1 = ""
rename codomiziladresse domiziladresse1
rename domiziladresse zusatzdomiziladresse1
rename zusatzdomiziladresse plzdomiziladresse1
rename plzdomiziladresse ortdomiziladresse1
rename ortdomiziladresse kantondomiziladresse1
rename kantondomiziladresse landdomiziladresse1
rename landdomiziladresse copostadresse1
rename copostadresse postadresse1
rename postadresse zusatzpostadresse1
rename zusatzpostadresse plzpostadresse1
rename plzpostadresse ortpostadresse1
rename ortpostadresse kantonpostadresse1
rename kantonpostadresse landpostadresse1
rename landpostadresse postfachadresse1
rename postfachadresse plzpostfachadresse1
rename plzpostfachadresse ortpostfachadresse1
rename ortpostfachadresse kantonpostfachadresse1
rename kantonpostfachadresse landpostfachadresse1
rename landpostfachadresse sprache1
rename sprache nogacode1
rename nogacode nogabeschrieb1
rename nogabeschrieb rechtsform1
rename rechtsform filialeindikator1
rename filialeindikator datumhandelsregistereintrag1
rename datumhandelsregistereintrag gruendungsjahr1
rename gruendungsjahr kapitalnominalag1
rename kapitalnominalag kapitaleinbezahltag1
rename kapitaleinbezahltag kapitalnominalgmbh1
rename kapitalnominalgmbh kapitaleinbezahltgmbh1
rename kapitaleinbezahltgmbh aktivitaetsstatus1
drop aktivitaetsstatus 
local varlist codomiziladresse domiziladresse zusatzdomiziladresse plzdomiziladresse ortdomiziladresse kantondomiziladresse landdomiziladresse /// 
copostadresse postadresse zusatzpostadresse plzpostadresse ortpostadresse kantonpostadresse landpostadresse postfachadresse ///
plzpostfachadresse ortpostfachadresse kantonpostfachadresse landpostfachadresse sprache nogacode nogabeschrieb rechtsform filialeindikator ///
datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus
foreach var in `varlist' {
	rename `var'1 `var'
}
order duns handelsregisternummer uid firma handelsname `varlist'
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_1a.dta", replace
restore 
drop if duns == "482467966"

* b) domiziadresse is empty: errors due to one shift in variable
preserve
keep if domiziladresse == "" & zusatzdomiziladresse!= ""
drop domiziladresse
rename zusatzdomiziladresse domiziladresse
rename plzdomiziladresse zusatzdomiziladresse1
rename ortdomiziladresse plzdomiziladresse1
rename kantondomiziladresse ortdomiziladresse1
rename landdomiziladresse kantondomiziladresse1
rename copostadresse landdomiziladresse1
rename postadresse copostadresse1
rename zusatzpostadresse postadresse1
rename plzpostadresse zusatzpostadresse1
rename ortpostadresse plzpostadresse1
rename kantonpostadresse ortpostadresse1
rename landpostadresse kantonpostadresse1
rename postfachadresse landpostadresse1
rename plzpostfachadresse postfachadresse1
rename ortpostfachadresse plzpostfachadresse1
rename kantonpostfachadresse ortpostfachadresse1
rename landpostfachadresse kantonpostfachadresse1
rename sprache landpostfachadresse1
rename nogacode sprache1
rename nogabeschrieb nogacode1
rename rechtsform nogabeschrieb1
rename filialeindikator rechtsform1
rename datumhandelsregistereintrag filialeindikator1
rename gruendungsjahr datumhandelsregistereintrag1
rename kapitalnominalag gruendungsjahr1
rename kapitaleinbezahltag kapitalnominalag1
rename kapitalnominalgmbh kapitaleinbezahltag1
rename kapitaleinbezahltgmbh kapitalnominalgmbh1
rename aktivitaetsstatus kapitaleinbezahltgmbh1
rename v37 aktivitaetsstatus1
local varlist zusatzdomiziladresse plzdomiziladresse ortdomiziladresse kantondomiziladresse landdomiziladresse /// 
copostadresse postadresse zusatzpostadresse plzpostadresse ortpostadresse kantonpostadresse landpostadresse postfachadresse ///
plzpostfachadresse ortpostfachadresse kantonpostfachadresse landpostfachadresse sprache nogacode nogabeschrieb rechtsform filialeindikator ///
datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus
foreach var in `varlist' {
	rename `var'1 `var'
}
order duns handelsregisternummer uid firma handelsname domiziladresse `varlist'
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_1b.dta", replace
restore 

* c) domiziadresse is non-empty: still errors due to one shift in variable
keep if domiziladresse != ""
drop if duns == "482500287" | duns == "485308451"
replace handelsname = handelsname + ", " + codomiziladresse 
replace handelsname = "" if handelsname == ", "
replace codomiziladresse = ""
rename ortpostadresse plzpostadresse1
rename kantonpostadresse ortpostadresse1 
rename landpostadresse kantonpostadresse1
rename postfachadresse landpostadresse1 
rename plzpostfachadresse postfachadresse1 
rename ortpostfachadresse plzpostfachadresse1 
rename kantonpostfachadresse ortpostfachadresse1 
rename landpostfachadresse kantonpostfachadresse1 
rename sprache landpostfachadresse1 
rename nogacode sprache1 
rename nogabeschrieb nogacode1 
rename rechtsform nogabeschrieb1 
rename filialeindikator rechtsform1 
rename datumhandelsregistereintrag filialeindikator1 
rename gruendungsjahr datumhandelsregistereintrag1 
rename kapitalnominalag gruendungsjahr1 
rename kapitaleinbezahltag kapitalnominalag1 
rename kapitalnominalgmbh kapitaleinbezahltag1
rename kapitaleinbezahltgmbh kapitalnominalgmbh1 
rename aktivitaetsstatus kapitaleinbezahltgmbh1 
rename v37 aktivitaetsstatus1
drop plzpostadresse
gen v37 = ""
local varlist plzpostadresse ortpostadresse kantonpostadresse landpostadresse postfachadresse ///
plzpostfachadresse ortpostfachadresse kantonpostfachadresse landpostfachadresse sprache nogacode nogabeschrieb rechtsform filialeindikator ///
datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus 
foreach var in `varlist' {
	rename `var'1 `var'
}
order duns handelsregisternummer uid firma handelsname domiziladresse codomiziladresse domiziladresse zusatzdomiziladresse plzdomiziladresse ///
ortdomiziladresse kantondomiziladresse landdomiziladresse copostadresse postadresse zusatzpostadresse plzpostadresse `varlist' v37 v38 v39 v40
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_1c.dta", replace 

* d) specific case correction
use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_raw.dta", clear
keep if duns == "482500287"
replace handelsname = handelsname + ", " + codomiziladresse 
replace handelsname = "" if handelsname == ", "
replace codomiziladresse = ""
rename domiziladresse codomiziladresse1
rename plzdomiziladresse domiziladresse1
rename ortdomiziladresse zusatzdomiziladresse1
rename kantondomiziladresse plzdomiziladresse1
rename landdomiziladresse ortdomiziladresse1
rename copostadresse kantondomiziladresse1
rename postadresse landdomiziladresse1
rename zusatzpostadresse copostadresse1
rename plzpostadresse postadresse1
rename ortpostadresse zusatzpostadresse1
rename kantonpostadresse plzpostadresse1
rename landpostadresse ortpostadresse1
rename postfachadresse kantonpostadresse1
rename plzpostfachadresse landpostadresse1
rename ortpostfachadresse postfachadresse1
rename kantonpostfachadresse plzpostfachadresse1
rename landpostfachadresse ortpostfachadresse1
rename sprache kantonpostfachadresse1
rename nogacode landpostfachadresse1
rename nogabeschrieb sprache1
rename rechtsform nogacode1
rename filialeindikator nogabeschrieb1
rename datumhandelsregistereintrag rechtsform1
rename gruendungsjahr filialeindikator1
rename kapitalnominalag datumhandelsregistereintrag1
rename kapitaleinbezahltag gruendungsjahr1
rename kapitalnominalgmbh kapitalnominalag1
rename kapitaleinbezahltgmbh kapitaleinbezahltag1
rename aktivitaetsstatus kapitalnominalgmbh1
rename v37 kapitaleinbezahltgmbh1
rename v38 aktivitaetsstatus1
drop codomiziladresse zusatzdomiziladresse
gen v37 = ""
gen v38 = ""
local varlist codomiziladresse zusatzdomiziladresse domiziladresse plzdomiziladresse ortdomiziladresse kantondomiziladresse landdomiziladresse /// 
copostadresse postadresse zusatzpostadresse plzpostadresse ortpostadresse kantonpostadresse landpostadresse postfachadresse ///
plzpostfachadresse ortpostfachadresse kantonpostfachadresse landpostfachadresse sprache nogacode nogabeschrieb rechtsform filialeindikator ///
datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus 
foreach var in `varlist' {
	rename `var'1 `var'
}
order duns handelsregisternummer uid firma handelsname `varlist'
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_1d.dta", replace

* e) specific case correction
use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_raw.dta", clear
keep if duns == "485308451"
replace handelsname = handelsname + ", " + codomiziladresse 
replace handelsname = "" if handelsname == ", "
replace codomiziladresse = ""
rename zusatzdomiziladresse domiziladresse1
rename plzdomiziladresse zusatzdomiziladresse1
rename ortdomiziladresse plzdomiziladresse1
rename kantondomiziladresse ortdomiziladresse1
rename landdomiziladresse kantondomiziladresse1
rename copostadresse landdomiziladresse1
rename postadresse copostadresse1
rename zusatzpostadresse postadresse1
rename plzpostadresse zusatzpostadresse1
rename ortpostadresse plzpostadresse1
rename kantonpostadresse ortpostadresse1
rename landpostadresse kantonpostadresse1
rename postfachadresse landpostadresse1
rename plzpostfachadresse postfachadresse1
rename ortpostfachadresse plzpostfachadresse1
rename kantonpostfachadresse ortpostfachadresse1
rename landpostfachadresse kantonpostfachadresse1
rename sprache landpostfachadresse1
rename nogacode sprache1
rename nogabeschrieb nogacode1
rename rechtsform nogabeschrieb1
rename filialeindikator rechtsform1
rename datumhandelsregistereintrag filialeindikator1
rename gruendungsjahr datumhandelsregistereintrag1
rename kapitalnominalag gruendungsjahr1
rename kapitaleinbezahltag kapitalnominalag1
rename kapitalnominalgmbh kapitaleinbezahltag1
rename kapitaleinbezahltgmbh kapitalnominalgmbh1
rename aktivitaetsstatus kapitaleinbezahltgmbh1
rename v37 aktivitaetsstatus1
rename v38 v371
drop domiziladresse 
gen v38 = ""
local varlist zusatzdomiziladresse domiziladresse plzdomiziladresse ortdomiziladresse kantondomiziladresse landdomiziladresse /// 
copostadresse postadresse zusatzpostadresse plzpostadresse ortpostadresse kantonpostadresse landpostadresse postfachadresse ///
plzpostfachadresse ortpostfachadresse kantonpostfachadresse landpostfachadresse sprache nogacode nogabeschrieb rechtsform filialeindikator ///
datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus v37
foreach var in `varlist' {
	rename `var'1 `var'
}
order duns handelsregisternummer uid firma handelsname codomiziladresse `varlist'
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_1e.dta", replace

* f) specific case correction
use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_raw.dta", clear
keep if duns == "483675919"
rename ortdomiziladresse plzdomiziladresse1
rename kantondomiziladresse ortdomiziladresse1
rename landdomiziladresse kantondomiziladresse1
rename copostadresse landdomiziladresse1
rename postadresse copostadresse1
rename zusatzpostadresse postadresse1
rename plzpostadresse zusatzpostadresse1
rename ortpostadresse plzpostadresse1
rename kantonpostadresse ortpostadresse1
rename landpostadresse kantonpostadresse1
rename postfachadresse landpostadresse1
rename plzpostfachadresse postfachadresse1
rename ortpostfachadresse plzpostfachadresse1
rename kantonpostfachadresse ortpostfachadresse1
rename landpostfachadresse kantonpostfachadresse1
rename sprache landpostfachadresse1
rename nogacode sprache1
rename nogabeschrieb nogacode1
rename rechtsform nogabeschrieb1
rename filialeindikator rechtsform1
rename datumhandelsregistereintrag filialeindikator1
rename gruendungsjahr datumhandelsregistereintrag1
rename kapitalnominalag gruendungsjahr1
rename kapitaleinbezahltag kapitalnominalag1
rename kapitalnominalgmbh kapitaleinbezahltag1
rename kapitaleinbezahltgmbh kapitalnominalgmbh1
rename aktivitaetsstatus kapitaleinbezahltgmbh1
rename v37 aktivitaetsstatus1
drop plzdomiziladresse
gen v37 = ""
local varlist plzdomiziladresse ortdomiziladresse kantondomiziladresse landdomiziladresse /// 
copostadresse postadresse zusatzpostadresse plzpostadresse ortpostadresse kantonpostadresse landpostadresse postfachadresse ///
plzpostfachadresse ortpostfachadresse kantonpostfachadresse landpostfachadresse sprache nogacode nogabeschrieb rechtsform filialeindikator ///
datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus
foreach var in `varlist' {
	rename `var'1 `var'
}
order duns handelsregisternummer uid firma handelsname codomiziladresse domiziladresse zusatzdomiziladresse plzdomiziladresse `varlist' v37 v38 v39 v40
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_1f.dta", replace

* append subsamples a-f of problem case 1 (left of nogabeschrieb)
use "$path\02_Processed_data\11_Directors_1994_2018\subsample_1a.dta", clear
foreach i in b c d e f {
	append using "$path\02_Processed_data\11_Directors_1994_2018\subsample_1`i'.dta", force
}
foreach var of varlist _all {
	label var `var' ""
}
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_1.dta", replace
foreach i in a b c d e f {
	erase "$path\02_Processed_data\11_Directors_1994_2018\subsample_1`i'.dta"
}

** Reconstruct corrected Bisnode data: Sample 1 (left of nogabeschrieb)
* Drop all observations that will be corrected by batch of observations (subsamples)
* Append corrected subsamples

use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_raw.dta", clear
drop if nogacode =="CH" | nogacode =="Detailhandel mit Herrenbekleidung" | nogacode =="deutsch" | nogacode =="französisch" | nogacode =="italienisch" | rechtsform == "deutsch"
append using "$path\02_Processed_data\11_Directors_1994_2018\subsample_1.dta"
save "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_temp.dta", replace

*** 2) errors to the right of nogabeschrieb (the variable that causes most problems)
use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_temp.dta", clear
tab rechtsform

drop if rechtsform == "Aktiengesellschaft" | rechtsform == "Einfache Gesellschaft" | rechtsform == "Einzelfirma" | rechtsform == "Genossenschaft" | ///
rechtsform == "Gesellschaft mit beschränkter Haftung" | rechtsform == "Kollektivgesellschaft" | rechtsform == "Kommanditgesellschaft für kollektive Kapitalanlagen" | ///
rechtsform == "Kommanditgesellschaft" | rechtsform == "Stiftung" | rechtsform == "Verein" | rechtsform == "Zweigniederlassung" | ///
rechtsform == "Zweigniederlassung (Hauptsitz im Ausland)" | rechtsform == "Öffentlich-rechtliche Institution" | rechtsform == "Rechtsform unbekannt" | ///
rechtsform == "Anstalt" | rechtsform == "Duns Support Record" | rechtsform == "FL - Aktiengesellschaft" | rechtsform == "Haupt von Gemeinderschaften" | ///
rechtsform == "Investmentgesellschaft mit variablem Kapital (SICAV)" // last 5 entries are proper rechtsformen (not due to column shifts)

/*
relevant cases are 
a) if v37!="" & v38=="" & v39=="" & v40==""
b) if v37=="" & v38!="" & v39=="" & v40=="" 
c) if v37!="" & v38!="" & v39=="" & v40=="" 
d) if v37=="" & v38=="" & v39!="" & v40=="" 
e) one special case: if duns == "487483919"
f) if v37=="" & v38=="" & v39=="" & v40!="" 
g) two special cases: if duns == "480010742" | duns == "487049678"
*/

* a) one right shift (v37 non-empty, rest empty)
preserve
keep if v37!="" & v38=="" & v39=="" & v40=="" 
replace nogabeschrieb = nogabeschrieb + ", " + rechtsform
drop rechtsform
rename filialeindikator rechtsform1
rename datumhandelsregistereintrag filialeindikator1
rename gruendungsjahr datumhandelsregistereintrag1
rename kapitalnominalag gruendungsjahr1
rename kapitaleinbezahltag kapitalnominalag1
rename kapitalnominalgmbh kapitaleinbezahltag1
rename kapitaleinbezahltgmbh kapitalnominalgmbh1
rename aktivitaetsstatus kapitaleinbezahltgmbh1
rename v37 aktivitaetsstatus1
foreach var in rechtsform filialeindikator datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus {
	rename `var'1 `var'
}
tab rechtsform
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_2a.dta", replace
restore

* b) two right shifts (v38 non-empty, rest empty)
preserve
keep if v37=="" & v38!="" & v39=="" & v40==""
drop if duns == "480010742" | duns == "487049678"  // special case g)
replace nogabeschrieb = nogabeschrieb + ", " + rechtsform + ", " + filialeindikator
drop rechtsform filialeindikator
rename datumhandelsregistereintrag rechtsform1
rename gruendungsjahr filialeindikator1
rename kapitalnominalag datumhandelsregistereintrag1
rename kapitaleinbezahltag gruendungsjahr1
rename kapitalnominalgmbh kapitalnominalag1
rename kapitaleinbezahltgmbh kapitaleinbezahltag1
rename aktivitaetsstatus kapitalnominalgmbh1
rename v37 kapitaleinbezahltgmbh1
rename v38 aktivitaetsstatus1
foreach var in rechtsform filialeindikator datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus {
	rename `var'1 `var'
}
tab rechtsform
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_2b.dta", replace
restore

* c) two right shifts (v37 & v38 non-empty, rest empty)
preserve
keep if v37!="" & v38!="" & v39=="" & v40=="" 
replace nogabeschrieb = nogabeschrieb + ", " + rechtsform + ", " + filialeindikator
drop rechtsform filialeindikator
rename datumhandelsregistereintrag rechtsform1
rename gruendungsjahr filialeindikator1
rename kapitalnominalag datumhandelsregistereintrag1
rename kapitaleinbezahltag gruendungsjahr1
rename kapitalnominalgmbh kapitalnominalag1
rename kapitaleinbezahltgmbh kapitaleinbezahltag1
rename aktivitaetsstatus kapitalnominalgmbh1
rename v37 kapitaleinbezahltgmbh1
rename v38 aktivitaetsstatus1
foreach var in rechtsform filialeindikator datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus {
	rename `var'1 `var'
}
tab rechtsform
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_2c.dta", replace
restore

* d) three right shifts (v37 & v38 & v39 non-empty, v40 empty)
preserve
keep if v37=="" & v38=="" & v39!="" & v40==""
drop if duns == "487483919"
replace handelsname = handelsname + ", " + domiziladresse
rename zusatzdomiziladresse codomiziladresse1
rename ortdomiziladresse domiziladresse1
rename landdomiziladresse plzdomiziladresse1
rename plzdomiziladresse kantondomiziladresse1
rename zusatzpostadresse landdomiziladresse1
rename copostadresse ortdomiziladresse1
rename plzpostadresse postadresse1
rename kantonpostadresse plzpostadresse1
rename landpostadresse ortpostadresse1
rename postfachadresse kantonpostadresse1
rename plzpostfachadresse landpostadresse1
rename nogabeschrieb sprache1
rename rechtsform nogacode1
rename filialeindikator nogabeschrieb1
rename datumhandelsregistereintrag rechtsform1
rename gruendungsjahr filialeindikator1
rename kapitalnominalag datumhandelsregistereintrag1
rename kapitaleinbezahltag gruendungsjahr1
rename kapitalnominalgmbh kapitalnominalag1
rename kapitaleinbezahltgmbh kapitaleinbezahltag1
rename aktivitaetsstatus kapitalnominalgmbh1
rename v37 kapitaleinbezahltgmbh1
rename v39 aktivitaetsstatus1
drop codomiziladresse domiziladresse kantondomiziladresse postadresse ortpostadresse ortpostfachadresse kantonpostfachadresse landpostfachadresse sprache nogacode
foreach var in codomiziladresse domiziladresse plzdomiziladresse kantondomiziladresse landdomiziladresse ortdomiziladresse ///
postadresse plzpostadresse ortpostadresse kantonpostadresse landpostadresse sprache nogacode nogabeschrieb rechtsform filialeindikator datumhandelsregistereintrag /// 
gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus {
	rename `var'1 `var'
}
tab rechtsform
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_2d.dta", replace
restore

* e) special case in three right shifts (v37 & v38 & v39 non-empty, v40 empty)
preserve
keep if duns == "487483919"
replace nogabeschrieb = nogabeschrieb + ", " + rechtsform + ", " + filialeindikator
rename datumhandelsregistereintrag rechtsform1
rename gruendungsjahr filialeindikator1
rename kapitalnominalag datumhandelsregistereintrag1
rename kapitaleinbezahltag gruendungsjahr1
rename kapitalnominalgmbh kapitalnominalag1
rename kapitaleinbezahltgmbh kapitaleinbezahltag1
rename aktivitaetsstatus kapitalnominalgmbh1
rename v37 kapitaleinbezahltgmbh1
rename v39 aktivitaetsstatus1
drop rechtsform filialeindikator
foreach var in rechtsform filialeindikator datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus {
	rename `var'1 `var'
}
tab rechtsform
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_2e.dta", replace 
restore

* f) four right shifts (v40 non-empty, else empty)
preserve
keep if v37=="" & v38=="" & v39=="" & v40!="" 
replace codomiziladresse = codomiziladresse + ", " + domiziladresse + ", " + zusatzdomiziladresse + ", " + plzdomiziladresse
rename kantondomiziladresse domiziladresse1
rename landdomiziladresse zusatzdomiziladresse1
rename copostadresse plzdomiziladresse1
rename postadresse ortdomiziladresse1
rename zusatzpostadresse kantondomiziladresse1
rename plzpostadresse landdomiziladresse1
rename ortpostadresse copostadresse1
rename kantonpostadresse postadresse1
rename landpostadresse zusatzpostadresse1
rename postfachadresse plzpostadresse1
rename plzpostfachadresse ortpostadresse1
rename ortpostfachadresse kantonpostadresse1
rename kantonpostfachadresse landpostadresse1
rename landpostfachadresse postfachadresse1
rename sprache plzpostfachadresse1
rename nogacode ortpostfachadresse1
rename nogabeschrieb kantonpostfachadresse1
rename rechtsform landpostfachadresse1
rename filialeindikator sprache1
rename datumhandelsregistereintrag nogacode1
rename gruendungsjahr nogabeschrieb1
rename kapitalnominalag rechtsform1
rename kapitaleinbezahltag filialeindikator1
rename kapitalnominalgmbh datumhandelsregistereintrag1
rename kapitaleinbezahltgmbh gruendungsjahr1
rename aktivitaetsstatus kapitalnominalag1
rename v37 kapitaleinbezahltag1
rename v38 kapitalnominalgmbh1
rename v39 kapitaleinbezahltgmbh1
rename v40 aktivitaetsstatus1
drop codomiziladresse domiziladresse zusatzdomiziladresse plzdomiziladresse ortdomiziladresse
foreach var in domiziladresse zusatzdomiziladresse plzdomiziladresse ortdomiziladresse kantondomiziladresse landdomiziladresse copostadresse postadresse ///
zusatzpostadresse plzpostadresse ortpostadresse kantonpostadresse landpostadresse postfachadresse plzpostfachadresse ortpostfachadresse kantonpostfachadresse /// 
landpostfachadresse sprache nogacode nogabeschrieb rechtsform filialeindikator datumhandelsregistereintrag gruendungsjahr kapitalnominalag kapitaleinbezahltag ///
kapitalnominalgmbh kapitaleinbezahltgmbh aktivitaetsstatus  {
	rename `var'1 `var'
}
tab rechtsform
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_2f.dta", replace 
restore

* g) two special cases
preserve
keep if duns == "480010742" | duns == "487049678"
replace nogabeschrieb = nogabeschrieb + ", " + rechtsform
drop rechtsform
rename filialeindikator rechtsform1
rename datumhandelsregistereintrag filialeindikator1
rename gruendungsjahr datumhandelsregistereintrag1
rename kapitalnominalag gruendungsjahr1
rename v38 aktivitaetsstatus1
drop aktivitaetsstatus
foreach var in rechtsform filialeindikator datumhandelsregistereintrag gruendungsjahr aktivitaetsstatus {
	rename `var'1 `var'
}
tab rechtsform
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_2g.dta", replace 
restore

* append subsamples a-g of problem case 2 (right of nogabeschrieb)
use "$path\02_Processed_data\11_Directors_1994_2018\subsample_2a.dta", clear
foreach i in b c d e f g {
	append using "$path\02_Processed_data\11_Directors_1994_2018\subsample_2`i'.dta", force
}
foreach var of varlist _all {
	label var `var' ""
}
save "$path\02_Processed_data\11_Directors_1994_2018\subsample_2.dta", replace
foreach i in a b c d e f g {
	erase "$path\02_Processed_data\11_Directors_1994_2018\subsample_2`i'.dta"
}

** Reconstruct corrected Bisnode data: Sample 1 (left of nogabeschrieb)
* Drop all observations that will be corrected by batch of observations (subsamples)
* Append corrected subsamples

use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_temp.dta", clear
keep if rechtsform == "Aktiengesellschaft" | rechtsform == "Einfache Gesellschaft" | rechtsform == "Einzelfirma" | rechtsform == "Genossenschaft" | ///
rechtsform == "Gesellschaft mit beschränkter Haftung" | rechtsform == "Kollektivgesellschaft" | rechtsform == "Kommanditgesellschaft für kollektive Kapitalanlagen" | ///
rechtsform == "Kommanditgesellschaft" | rechtsform == "Stiftung" | rechtsform == "Verein" | rechtsform == "Zweigniederlassung" | ///
rechtsform == "Zweigniederlassung (Hauptsitz im Ausland)" | rechtsform == "Öffentlich-rechtliche Institution" | rechtsform == "Rechtsform unbekannt" | ///
rechtsform == "Anstalt" | rechtsform == "Duns Support Record" | rechtsform == "FL - Aktiengesellschaft" | rechtsform == "Haupt von Gemeinderschaften" | ///
rechtsform == "Investmentgesellschaft mit variablem Kapital (SICAV)" | /// // last 5 entries are proper rechtsformen (not due to column shifts) 
duns == "480382568" | duns == "480382600" | duns == "480395024" | duns == "480433239" | duns == "480568794" | duns == "480617278" | duns == "480665327" | /// 
duns == "480722862" | duns == "480957450" | duns == "481111370" | duns == "481197320" | duns == "481739449" | duns == "482193492" | duns == "482991861" | ///
duns == "483612797" | duns == "483840526" | duns == "483858205" | duns == "485719801" | duns == "486017577" | duns == "486320542" | duns == "486717341" | ///
duns == "487531501" | duns == "487802415" | duns == "487888885" // these are observations with missing information on rechtsform

append using "$path\02_Processed_data\11_Directors_1994_2018\subsample_2.dta"

* formatting
destring duns, replace
format duns %10.0f
sort duns
destring nogacode, replace
destring gruendungsjahr, replace
destring kapitalnominalag, replace
destring kapitaleinbezahltag, replace 
destring kapitalnominalgmbh, replace 
destring kapitaleinbezahltgmbh, replace

/*
destring datumhandelsregistereintrag, format(format) replace
destring plzdomiziladresse, replace 
destring plzpostadresse, replace 
destring plzpostfachadresse, replace
*/

drop v37 v38 v39 v40
sort duns
duplicates report duns
save "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_corr.dta", replace

erase "$path\02_Processed_data\11_Directors_1994_2018\subsample_1.dta"
erase "$path\02_Processed_data\11_Directors_1994_2018\subsample_2.dta"
erase "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_temp.dta"

use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_corr.dta", clear


*** Combine Bisnode Person data with firm data
use "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Person_Geo.dta", clear
sort duns

merge m:1 duns using "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Firmen_corr.dta"

*br personenid anrede titel vorname nachname nationalität funktion duns if _merge==1

/*
    Result                      Number of obs
    -----------------------------------------
    Not matched                            12
        from master                        12  (_merge==1)
        from using                          0  (_merge==2)

    Matched                         5,898,177  (_merge==3)
    -----------------------------------------

personenid	anrede	titel	vorname	nachname	nationalität	funktion	duns
3544532	männlich		Bernard	Terrier		Leiter der Zweigniederlassung	480223555
3544533	männlich		Ralph	Bobbiá		Leiter der Zweigniederlassung	480223558
3544534	weiblich		Iris	Klausmann		Leiter der Zweigniederlassung	480223559
3544535	weiblich		Nadine	Meyer		Leiter der Zweigniederlassung	480223560
3546801	weiblich		Erika	Peter		Leiter der Zweigniederlassung	480224704
3546808	männlich		Bruno	Kuster		Betriebsleiter/in	480224714
1914932	männlich		Peter	Muster	Schweiz	Mitglied	481295843
1914931	männlich		Werner	Muster	Schweiz	Präsident/in	481295843
1915058	männlich	Dipl. Ing.	Vaclav	Sobriecz	Polen	Direktor/in	481295843
1914933	männlich		Hans	Muster	Schweiz	Mitglied	481295843
2008142	weiblich		Sonja	Kvas	Deutschland	Inhaber/in	481996440
3038797	männlich		Federico	Nannini	Italien	Prokurist/in	488252610
*/

drop _merge
sort personenid duns
save "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Person-Firmen_Geo.dta", replace

preserve
keep personenid duns anrede vorname nachname geburtstag eintrittdatum austrittdatum E_CNTR_w N_CNTR_w rechtsform gremium funktion firma
save "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_Person-Firmen_Geo_reduced.dta", replace
restore


* Check information overlap: legal statute (Rechtsform) vs. board name (Gremium)
log using "$path\02_Processed_data\11_Directors_1994_2018\Check_Gremium-Rechtsform.smcl", replace
tab gremium

preserve
keep if gremium == "Verwaltungsrat"
tab rechtsform
restore

preserve
keep if gremium == "Gesellschafter"
tab rechtsform
restore

preserve
keep if gremium == "Inhaber"
tab rechtsform
restore

preserve
keep if gremium == "Vorstand"
tab rechtsform
restore

preserve
keep if gremium == "Stiftungsrat"
tab rechtsform
restore

preserve
keep if gremium == "Zeichnungsberechtigte"
tab rechtsform
restore

cap log close

