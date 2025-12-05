
******************************************************
* (0) Set directories
******************************************************

	clear all
	set more off


    capture cd "D:\SchmidLu\Dropbox\Projekt Nationalräte"

	if _rc==0{ 
	global hauptpfad "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
	global datapath "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
	}


cd "$hauptpfad/7_Verwaltungsratsmandate/Bisnode/Datenanreicherung20180410"

// Note: We sent Bisnode a file with all combinations of different names, firstnames, origin, municipality for a person in the long format
//       We received a file from Bisnode in the long format. Therefore, the directory information may be merged to multiple observations per
//	     person. In this do file we merge all possible mandates to a candidate file in the wide format and drop duplicate mandates. 


* (i) Prepare our Input file with panel ID

set excelxlsxlargefile on
import excel  "Uni Fribourg Output 10-04-2018.xlsx"	, sheet("Input Uni Fribourg") clear first cellrange(A1:I36841)
keep ID PanelID
save PanelID.dta, replace



* (ii) Read-in file from Bisnode from 10. April 2018 (only sheet with SHAB information)  
	
import excel  "Uni Fribourg Output 10-04-2018.xlsx"	, sheet("Personen - Mandate zusätzlich") clear first cellrange(A1:AI118327)
merge m:1 ID using PanelID.dta, gen(merge_panelid)
drop ID Nachname Vorname Geburtsjahr Heimatort No_Heimatort Wohnort No_Wohnort 
rename PanelID ID // note that in our bisnode_out file the panel ID ("FR-1987-XXXX") is called ID 

save Bisnode_Personen_temp.dta, replace

import excel  "Uni Fribourg Output 10-04-2018.xlsx"	, sheet("Firmendaten") clear first cellrange(A2:CS51979)
save Bisnode_Firmen_temp.dta, replace



* (iii) Merge candidate information in wide format


use "$hauptpfad/1_data/12 bisnode/bisnode_out_wide.dta", clear 
drop if ID==""

merge 1:m ID using Bisnode_Personen_temp.dta, gen(merge_bisnode_personen)
merge m:1 DUNS using Bisnode_Firmen_temp.dta, gen(merge_bisnode_firmen)

sort  DUNS merge_bisnode_firmen
br if merge_bisnode_firmen==1|  merge_bisnode_firmen==2


drop if merge_panelid==2  // drop individuals with no mandate
drop if merge_bisnode_firmen==2 // drop firms related to individuals of B2C database (here we use B2B database)


rename U Vorname_Bisnode
rename V Nachname_Bisnode
rename Bürgerort Bürgerort_Bisnode
rename AA Wohnort_Bisnode
rename Location_current Wohnort_aktuell_Bisnode
rename Street_current Strasse_aktuell_Bisnode
rename Anrede Geschlecht_Bisnode
rename EintrittDatum EintrittDatum_Bisnode	 
rename AustrittDatum AustrittDatum_Bisnode	
rename DUNS DUNS_Bisnode
rename Gremium Gremium_Bisnode	
rename Funktion	 Funktion_Bisnode
rename Unterschrift Unterschrift_Bisnode
rename Titel Titel_Bisnode
rename Nationalität Nationalität_Bisnode
rename AustrittIndikator AustrittIndikator_Bisnode
rename Kapitalanteil Kapitalanteil_Bisnode	
rename Firma Firma_Bisnode




gen Birthdate_matched_day=day(Birthdate_matched)
gen Birthdate_matched_month=month(Birthdate_matched)
gen Birthdate_matched_year=year(Birthdate_matched)

tostring Birthdate_matched_*, replace

gen Geburtsdatum_Bisnode=Birthdate_matched_day+"."+Birthdate_matched_month+"."+Birthdate_matched_year  if Birthdate_year_only==0
replace Geburtsdatum_Bisnode=Birthdate_matched_year  if Birthdate_year_only==1



* (iv) Drop duplicates

duplicates report ID Inland-Kapitalanteil 
duplicates drop ID Inland-Kapitalanteil, force

* (v) Create out file for research assistants

tostring birthyear*, replace

prog combine
gen `1'=`1'1
forvalues i=2(1)`2'{
replace `1'=`1'+"; "+`1'`i' if `1'`i'!="" & `1'`i'!="."
}
end

combine name 3
combine firstname 3
combine origin 8
combine municipality 5
combine birthyear 2

sort ID PersonenID

* (vi) Exclude too old and too young candidates


foreach var of varlist birthyear1 birthyear2 Birthdate_matched_year{
gen `var'_cs=real(substr(`var',1,1))+real(substr(`var',2,1))+real(substr(`var',3,1))+real(substr(`var',4,1))
}

destring birthyear1 birthyear2 Birthdate_matched_year, replace



drop if min(abs(Birthdate_matched_year-birthyear1), abs(Birthdate_matched_year-birthyear2))>10 & !missing(Birthdate_matched_year) & ///
		birthyear1_cs !=Birthdate_matched_year_cs & birthyear2_cs !=Birthdate_matched_year_cs


* (vii) Exclude too old and too young candidates

bysort ID PersonenID: egen EintrittDatum_Bisnode_min=min(year(EintrittDatum_Bisnode))
bysort ID PersonenID: egen EintrittDatum_Bisnode_max=max(year(EintrittDatum_Bisnode))
bysort ID PersonenID: egen AustrittDatum_Bisnode_max=max(year(AustrittDatum_Bisnode))

	* (a) Eintritt


	drop if EintrittDatum_Bisnode_min-min(birthyear1,birthyear2)<8 & !missing(EintrittDatum_Bisnode_min) // we drop all individuals who are below 18 (=mündig) and allow an error range of 10 years

	drop if EintrittDatum_Bisnode_max-max(birthyear1,birthyear2)>100 & !missing(EintrittDatum_Bisnode_max) // we drop all individuals who are above 110 and allow an error range of 10 years


	* (b) Austritt

	drop if AustrittDatum_Bisnode_max-max(birthyear1,birthyear2)>110 & !missing(AustrittDatum_Bisnode_max) // we drop all individuals who are above 120 and allow an error range of 10 years



* (viii) Exclude persons with no directorship


gen gremium_vr=0
replace gremium_vr=1 if Gremium=="Verwaltungsrat" 
bysort PersonenID: egen gremium_vr_max=max(gremium_vr)
keep if gremium_vr_max==1

br Gremium gremium_other_max if name1=="Blocher" 


* (ix) Mark individuals for which one Bisnode observation was added to multiple observations in the election data

preserve
keep ID PersonenID
duplicates drop ID PersonenID, force
bysort PersonenID: gen irgendwas=_n
reshape wide ID, i(PersonenID) j(irgendwas)
gen MultipleIDs="" if ID2==""
replace MultipleIDs=ID1+"; "+ID2 if ID2!=""
replace MultipleIDs=MultipleIDs+"; "+ID3 if ID3!=""
replace MultipleIDs=MultipleIDs+"; "+ID4 if ID4!=""
keep PersonenID MultipleIDs
save MultipleIDs.dta, replace
restore

merge m:1 PersonenID using MultipleIDs, gen(merge_multipleids)


* (x) Create out file for research assistants

gen Geburtsdatum_New=""
gen Todesdatum_New=""

local keepvars ID PersonenID PID_current name Vorname_Bisnode firstname Nachname_Bisnode origin Bürgerort_Bisnode municipality Wohnort_Bisnode /// 
      Wohnort_aktuell_Bisnode Strasse_aktuell_Bisnode birthyear Geburtsdatum_Bisnode Geburtsdatum_New Todesdatum_New EintrittDatum_Bisnode	AustrittDatum_Bisnode DUNS_Bisnode ///
      Firma_Bisnode Gremium_Bisnode	Funktion_Bisnode	Unterschrift_Bisnode AustrittIndikator_Bisnode Kapitalanteil_Bisnode sex Geschlecht_Bisnode job Titel_Bisnode /// 
      Nationalität_Bisnode MultipleIDs

order `keepvars'
keep `keepvars'


export excel  "Bisnode_Out_Checkfile.xlsx"	,replace    firstrow(variables)

unique ID PersonenID
unique ID

erase PanelID.dta 
erase Bisnode_Personen_temp.dta
erase Bisnode_Firmen_temp.dta
