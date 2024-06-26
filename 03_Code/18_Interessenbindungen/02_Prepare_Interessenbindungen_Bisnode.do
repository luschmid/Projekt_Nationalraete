clear
cap log close
set more 1
version 17

global dataNR "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\"
global dataRL "E:\12. Cloud\Dropbox\Record Linkage\02_Data\"


*** Prepare Interessenbindungen Dataset
use "$dataNR\18_Interessenbindungen\Interessenbindungen_1996-2014.dta", clear
keep if Council == "NR"
keep if LegalForm == "AG"
keep if inlist(Organ, "VR", "VR*", "V" "CA", "CdA")   // all mandates year > 2012 are lost
drop if inlist(Function, "A", "AM", "Adm." "Dir", "Gf." "Sekr") 
keep year canton mp_id mp_name mp_pname PartyS Profession Mandate Function Affiliation Place BirthDate partylas WohnAdr WohnPlz WohnOrt Sex Hreg
gen birthyear = year(BirthDate)

gen mp_lname = subinstr(mp_name, mp_pname,"", .) 
replace mp_lname = ustrtrim(mp_lname)
replace mp_lname = "Fischer" if mp_name == "Fischer Theo (Hägglingen)"
replace WohnOrt = "Hägglingen" if mp_name == "Fischer Theo (Hägglingen)"
replace mp_lname = "Schneider" if mp_name == "Schneider Johann Niklaus"
replace mp_pname = "Johann Niklaus" if mp_name == "Schneider Johann Niklaus"
replace mp_lname = "Schmid-Federer" if mp_name == "Schmid-Federer Barbar"
replace mp_lname = "Baader" if mp_name == "Baader-Buri Casper"
replace mp_lname = "Darbellay" if mp_name == "Darbellay Christoph"
replace mp_lname = "Loeb" if mp_name == "Loeb Francois"
replace mp_lname = "Lachat" if mp_name == "Lachat Francois"
replace mp_lname = "Nebiker" if mp_name == "Nebiker Hans Rudolf"
replace mp_lname = "Leuba" if mp_name == "Leuba Jean-Francois"
replace mp_lname = "Borer" if mp_name == "Borer Roland"
replace mp_pname = "Adrian" if mp_pname == "Adriano" & mp_lname == "Imfeld"  // official name is actually "Adriano" on Wiki & NR website -- > Name in NR dataset is "Adrian"
replace mp_pname = "Brigitta Maria" if mp_pname == "Brigitta M." & mp_lname == "Gadient"
replace mp_pname = "Bruno J." if mp_pname == "Bruno" & mp_lname == "Zuppiger"
replace mp_lname = "Markwalder" if mp_lname == "Markwalder Bär" & mp_pname == "Christa" & year == 2007
replace mp_pname = "Andy" if mp_pname == "Andreas" & mp_lname == "Tschümperlin"
replace mp_lname = "Haering Binder" if mp_lname == "Haering" & mp_pname == "Barbara" & (year == 2002 | year == 2003)
replace mp_lname = "Langenberger" if mp_lname == "Langenberger-Jaeger" & mp_pname == "Christiane" & year == 1999
replace mp_lname = "von Rotz-Spichtig" if mp_lname == "von Rotz" & mp_pname == "Christoph"
replace mp_pname = "Edi" if mp_pname == "Eduard" & mp_lname == "Engelberger" & year > 1999
replace mp_pname = "Guido A." if mp_pname == "Guido" & mp_lname == "Zäch"
replace mp_pname = "Hans Rudolf" if mp_pname == "Hans-Rudolf" & mp_lname == "Nebiker"
replace mp_lname = "Maître" if mp_lname == "Maitre" & mp_pname == "Jean-Philippe" & year >= 1996 & year <= 1999
replace mp_lname = "Schneider-Ammann" if mp_lname == "Schneider" & mp_pname == "Johann Niklaus" 
replace mp_pname = "Johann-Niklaus" if mp_lname == "Schneider-Ammann" & mp_pname == "Johann Niklaus" & year >= 2000 & year <= 2003
replace mp_pname = "Johann" if mp_lname == "Schneider-Ammann" & mp_pname == "Johann Niklaus" & year >= 2004 & year <= 2011
replace mp_pname = "Johannes Robert" if mp_pname == "Johannes" & mp_lname == "Randegger" 
replace mp_lname = "Leu-Morgenthaler" if mp_lname == "Leu" & mp_pname == "Josef" & year >= 1996 & year <= 1999
replace mp_lname = "Nabholz-Haidegger" if mp_lname == "Nabholz" & mp_pname == "Lili" & year == 2000
replace mp_lname = "Binder-Gäumann" if mp_lname == "Binder" & mp_pname == "Max" & year >= 2004 & year <= 2007
replace mp_pname = "Ruedi" if mp_lname == "Rechsteiner" & mp_pname == "Rudolf" & year >= 1996 & year <= 2007
replace mp_pname = "Rudolf" if mp_lname == "Noser" & mp_pname == "Ruedi" & year >= 2008 & year <= 2011
replace mp_pname = "Rudolf" if mp_lname == "Aeschbacher" & mp_pname == "Ruedi" & year >= 2000 & year <= 2003
replace mp_pname = "Anton" if mp_lname == "Eberhard" & mp_pname == "Toni" & year >= 1996 & year <= 1999
replace mp_lname = "Fischer" if mp_lname == "Fischer-Seengen" & mp_pname == "Ulrich" & year >= 2002 & year <= 2003
replace mp_lname = "Haller-Iseli" if mp_lname == "Haller" & mp_pname == "Ursula" & year >= 2000 & year <= 2007
replace mp_lname = "Kleiner-Schläpfer" if mp_lname == "Kleiner" & mp_pname == "Marianne" & year >= 2004 & year <= 2011

order year canton mp_pname mp_lname WohnOrt WohnAdr WohnPlz PartyS partylas Profession Sex birthyear Affiliation Place Hreg Function Mandate mp_id
rename mp_pname firstname
rename mp_lname name
sort firstname name year
duplicates report firstname name year  // there is a maximum of 30 mandates per person
tab year
save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\IG_AG_1996-2012.dta", replace

*** Prepare Nationalrat
use "$dataNR\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta", clear
keep ID canton year name firstname birthyear sex job elected gdename gdenr_2018_w list*
keep if year>=1991
gen eleyear = year
expand 4, gen(tmp_year)
bysort ID year: gen n=_n  // the first year of office is the one following the election year (elections in October): therefore, not _n-1
replace year = year + n
drop tmp_year n
sort firstname name year
duplicates report firstname name canton year
tab year
rename ID NRid
save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_BisnodeRL.dta", replace
keep if elected == 1
save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_BisnodeRL_elected.dta", replace

*** Merge IG & NR
merge 1:m firstname name canton year using "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\IG_AG_1996-2012.dta"
sort firstname name year
drop if year>2011
drop if year<1996
*br * if _merge == 1
*br * if _merge == 2
gen source = "NR" if _merge == 1
replace source = "IG" if _merge == 2
replace source = "NR&IG" if _merge == 3
drop _merge
gen issue = "ok" if source == "NR&IG"
replace issue = "Nachrücker" if firstname == "Adrian" & name == "Imfeld" & year == 2003
replace issue = "Nachrücker" if firstname == "Charles-Albert" & name == "Antille" & year == 1999
replace issue = "Nachrücker" if firstname == "Andreas" & name == "Brönnimann" & (year == 2010 | year == 2011)
replace issue = "Nachrücker" if firstname == "Christian" & name == "Waber" & year == 1999
replace issue = "Nachrücker" if firstname == "Andreas" & name == "Zeller" & year == 2007
replace issue = "Nachrücker" if firstname == "Caspar" & name == "Baader" & year == 1999
replace issue = "Nachrücker" if firstname == "Christine" & name == "Keller" & (year == 1998 | year == 1999)
replace issue = "Nachrücker" if firstname == "Christine" & name == "Wirz-von Planta" & (year == 2002 | year == 2003)
replace issue = "Nachrücker" if firstname == "Eric" & name == "Nussbaumer"
replace issue = "Nachrücker" if firstname == "Ernst" & name == "Hasler"
replace issue = "Nachrücker" if firstname == "Ernst" & name == "Schibli" & year == 2003
replace issue = "Nachrücker" if firstname == "Felix" & name == "Walker" & year >= 1999 & year <= 2003
replace issue = "Nachrücker" if firstname == "Gilbert" & name == "Debons" & year == 1999
replace issue = "Nachrücker" if firstname == "Hans" & name == "Rutschmann" & year >= 2004 & year <= 2007
replace issue = "Nachrücker" if firstname == "Hans" & name == "Stöckli" & year >= 2005 & year <= 2007
replace issue = "Nachrücker" if firstname == "Laurent" & name == "Favre" & year >= 2009 & year <= 2011
replace issue = "Nachrücker" if firstname == "Lieni" & name == "Füglistaller" & year >= 2006 & year <= 2007
replace issue = "Nachrücker" if firstname == "Markus" & name == "Hutter"
replace issue = "Nachrücker" if firstname == "Markus" & name == "Zemp" & year == 2007
replace issue = "Nachrücker" if firstname == "Norman" & name == "Gobbi" & year == 2011
replace issue = "Nachrücker" if firstname == "Paul-André" & name == "Roux" & year == 2011
replace issue = "Nachrücker" if firstname == "Peter" & name == "Flück" & year == 2011
replace issue = "Nachrücker" if firstname == "Pierre" & name == "Salvi" & year == 2007
replace issue = "Nachrücker" if firstname == "Pierre" & name == "Tillmanns" & year >= 2000 & year <= 2001
replace issue = "Nachrücker" if firstname == "René" & name == "Vaudroz" & year == 2007
replace issue = "Nachrücker" if firstname == "Sebastian" & name == "Frehner" & year == 2011
replace issue = "Nachrücker" if firstname == "Thomas" & name == "Müller" & year == 2007
replace issue = "Nachrücker" if firstname == "Ulrich" & name == "Schlüer" & year >= 2010 & year <= 2011
replace issue = "Nachrücker" if firstname == "Urs" & name == "Bernhardsgrütter" & year == 2007
replace issue = "Nachrücker" if firstname == "Urs" & name == "Hany" & year == 2007
replace issue = "Nachrücker" if firstname == "Urs" & name == "Schweizer" & year == 2007
replace issue = "Nachrücker" if firstname == "Walter" & name == "Bosshard" & year >= 1996 & year <= 1999
replace issue = "Nachrücker" if firstname == "Boris" & name == "Banga" & year >= 2000 & year <= 2003
replace issue = "Nachrücker" if firstname == "Corina" & name == "Eichenberger-Walther" & year >= 2008 & year <= 2011
replace issue = "Nachrücker" if firstname == "Fabio" & name == "Abate" & year == 2001
replace issue = "Nachrücker" if firstname == "Jacques-André" & name == "Maire" & year == 2011
replace issue = "Nachrücker" if firstname == "Margret" & name == "Kiener Nellen" & year == 2006
replace issue = "Nachrücker" if firstname == "Ruth" & name == "Kalbermatten" & year == 1999
replace issue = "Nachrücker" if firstname == "Viola" & name == "Amherd" & year >= 2006 & year <= 2007

* Include NR ID
replace NRid = "OW-2003-7055" if firstname == "Adrian" & name == "Imfeld" & year == 2003
replace NRid = "VS-1995-0002" if firstname == "Charles-Albert" & name == "Antille" & year == 1999
replace NRid = "BEJU-1999-0043" if firstname == "Andreas" & name == "Brönnimann" & (year == 2010 | year == 2011)
replace NRid = "BEJU-1983-0398" if firstname == "Christian" & name == "Waber" & year == 1999
replace NRid = "SG-2003-0157" if firstname == "Andreas" & name == "Zeller" & year == 2007
replace NRid = "BL-1995-0007" if firstname == "Caspar" & name == "Baader" & year == 1999
replace NRid = "BS-1995-0032" if firstname == "Christine" & name == "Keller" & (year == 1998 | year == 1999)
replace NRid = "BS-1987-0081" if firstname == "Christine" & name == "Wirz-von Planta" & (year == 2002 | year == 2003)
replace NRid = "BL-2003-0056" if firstname == "Eric" & name == "Nussbaumer"
replace NRid = "AG-1995-9071" if firstname == "Ernst" & name == "Hasler"
replace NRid = "ZH-1991-0583" if firstname == "Ernst" & name == "Schibli" & year == 2003
replace NRid = "SG-1999-0163" if firstname == "Felix" & name == "Walker" & year >= 1999 & year <= 2003
replace NRid = "VS-1995-0017" if firstname == "Gilbert" & name == "Debons" & year == 1999
replace NRid = "ZH-1987-0584" if firstname == "Hans" & name == "Rutschmann" & year >= 2004 & year <= 2007
replace NRid = "BEJU-1987-0459" if firstname == "Hans" & name == "Stöckli" & year >= 2005 & year <= 2007
replace NRid = "NE-2007-0021" if firstname == "Laurent" & name == "Favre" & year >= 2009 & year <= 2011
replace NRid = "AG-1995-0049" if firstname == "Lieni" & name == "Füglistaller" & year >= 2006 & year <= 2007
replace NRid = "ZH-1999-0374" if firstname == "Markus" & name == "Hutter"
replace NRid = "AG-1995-0212" if firstname == "Markus" & name == "Zemp" & year == 2007
replace NRid = "TI-1999-0032" if firstname == "Norman" & name == "Gobbi" & year == 2011
replace NRid = "VS-2007-0103" if firstname == "Paul-André" & name == "Roux" & year == 2011
replace NRid = "BEJU-2003-9058" if firstname == "Peter" & name == "Flück" & year == 2011
replace NRid = "VD-1999-0211" if firstname == "Pierre" & name == "Salvi" & year == 2007
replace NRid = "VD-1987-9011" if firstname == "Pierre" & name == "Tillmanns" & year >= 2000 & year <= 2001
replace NRid = "VD-1995-0191" if firstname == "René" & name == "Vaudroz" & year == 2007
replace NRid = "BS-2003-0014" if firstname == "Sebastian" & name == "Frehner" & year == 2011
replace NRid = "SG-1995-0126" if firstname == "Thomas" & name == "Müller" & year == 2007
replace NRid = "ZH-1971-9246" if firstname == "Ulrich" & name == "Schlüer" & year >= 2010 & year <= 2011
replace NRid = "SG-1991-0010" if firstname == "Urs" & name == "Bernhardsgrütter" & year == 2007
replace NRid = "ZH-1991-0259" if firstname == "Urs" & name == "Hany" & year == 2007
replace NRid = "BS-1999-0060" if firstname == "Urs" & name == "Schweizer" & year == 2007
replace NRid = "ZH-1987-0081" if firstname == "Walter" & name == "Bosshard" & year >= 1996 & year <= 1999
replace NRid = "SO-1995-0002" if firstname == "Boris" & name == "Banga" & year >= 2000 & year <= 2003
replace NRid = "AG-1995-0036" if firstname == "Corina" & name == "Eichenberger-Walther" & year >= 2008 & year <= 2011
replace NRid = "TI-1999-0001" if firstname == "Fabio" & name == "Abate" & year == 2001
replace NRid = "NE-2007-0036" if firstname == "Jacques-André" & name == "Maire" & year == 2011
replace NRid = "BEJU-2003-0227" if firstname == "Margret" & name == "Kiener Nellen" & year == 2006
replace NRid = "VS-1991-0018" if firstname == "Ruth" & name == "Kalbermatten" & year == 1999
replace NRid = "VS-2003-0002" if firstname == "Viola" & name == "Amherd" & year >= 2006 & year <= 2007

order source issue
sort NRid year Hreg
save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_IG_NRid_BisnodeRL.dta", replace

*** Merge Bisnode info to Generation 7 (Bisnode)
* Step 1) retain only RL links of NRid and VRids: use "clusterid" to construct data-key-dataset
* Step 2) merge NRids on Bisnode data -- > NR with firm connection (Bisnode)
* Step 3) merge IG data on NR_Bisnode sample via NRid

** Step 1): construct data key: NR <--> VR link
import delim using "$dataRL\07_Output_RL_Bisnode\RL-Results-02ca62d1-d407-43dd-9284-0e24a2a011db-generation_7.csv", delimiters(",") encoding("UTF-8") clear

preserve
keep clusterid sourcefile id firstname name e_cntr_w n_cntr_w birthyear sex
keep if sourcefile == 0
rename id NRid
rename firstname firstname_pol 
rename name name_pol
rename e_cntr_w e_cntr_w_pol
rename n_cntr_w n_cntr_w_pol
rename birthyear birthyear_pol
rename sex sex_pol
replace sex = "M" if sex == "1"
replace sex = "W" if sex == "0"
label var NRid "ID Nationalrat"
drop sourcefile
drop if clusterid == .
save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_NRid.dta", replace
restore

preserve
keep clusterid sourcefile id firstname name e_cntr_w n_cntr_w birthyear sex linkscore
keep if sourcefile == 1
rename firstname firstname_bis 
rename name name_bis
rename e_cntr_w e_cntr_w_bis
rename n_cntr_w n_cntr_w_bis
rename birthyear birthyear_bis
rename sex sex_bis
destring id, gen(VRid)
label var VRid "ID Person Bisnode"
drop sourcefile id
drop if clusterid == .
save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_VRid.dta", replace
restore

use "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_VRid.dta", clear
merge 1:1 clusterid using "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_NRid.dta"
drop if _merge == 2 // 5 observations only contain NRid without a link to a VRid
drop _merge
sort VRid
duplicates report VRid NRid
duplicates drop VRid NRid, force
bysort VRid: gen nmand = _n
bysort VRid: egen maxVRid = max(nmand)
sum maxVRid
local max = r(max)
order VRid NRid firstname_* name_* birthyear_* sex_* e_cntr_w_* n_cntr_w_*
save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_VRid_NRid.dta", replace

forv i = 1(1)`max' {
	preserve
	keep VRid NRid clusterid linkscore nmand maxVRid
	keep if nmand == `i'
	drop maxVRid nmand
	save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_VRid_NRid_`i'.dta", replace
	restore
}

** Step 2): merge NRids on Bisnode data -- > NRid with firm connection (Bisnode)
use "$dataNR\11_Directors_1994_2018\Bisnode_Person-Firmen_Geo.dta", clear
keep if gremium == "Verwaltungsrat"
keep if rechtsform == "Aktiengesellschaft"
drop if inlist(funktion, "Aktuar/in (nicht Mitglied)", "Sekretär/in (nicht Mitglied)", "Beisitzer/in", "Ausseramtliche/r Konkursverwalter/in", "Generalsekretär/in (nicht Mitglied)", "Liquidator/in", "Kassier/in (nicht Mitglied)", "Protokollführer/in (nicht Mitglied)", "Verwalter/in (nicht Mitglied)")
sort personenid
rename personenid VRid
keep VRid duns vorname nachname W_PLZ4 wohnort funktion eintrittdatum austrittdatum E_CNTR_w N_CNTR_w uid handelsregisternummer firma handelsname plzdomiziladresse ortdomiziladresse kantondomiziladresse nogacode unterschrift kapitalnominalag kapitaleinbezahltag
sort VRid

forv i = 1(1)`max' {  // `max'
	preserve
	merge m:1 VRid using "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_VRid_NRid_`i'.dta"
	drop if _merge == 1
	drop if _merge == 2 // these are links to non-AG's (with all observations, there are no _merge == 2)
	drop _merge
	save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_Bisnode-NRid_`i'.dta", replace
	restore	
}
use "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_Bisnode-NRid_1.dta", clear
forv i = 2(1)`max' {
	append using "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_Bisnode-NRid_`i'.dta"
}
order NRid VRid vorname nachname wohnort W_PLZ4 E_CNTR_w N_CNTR_w funktion unterschrift duns uid handelsregisternummer firma handelsname plzdomiziladresse ortdomiziladresse kantondomiziladresse nogacode kapitalnominalag kapitaleinbezahltag
save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_Bisnode-NRid.dta", replace

* Panel structure
*use "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_Bisnode-NRid.dta", clear
bysort VRid NRid duns funktion: gen mandid = _n
expand 25, gen(expy)
bysort VRid NRid duns funktion mandid: gen year=_n + 1993
drop expy
sort VRid NRid duns funktion mandid year
order year VRid NRid duns funktion mandid

gen entrydate = date(eintrittdatum, "YMD")
format entrydate %td
gen exitdate = date(austrittdatum, "YMD")
format exitdate %td
gen flag = 1 if exitdate<entrydate & entrydate != . & exitdate != .
gen exitdate_tmp = exitdate
replace exitdate = exitdate if flag == 1 // assumption: confused entry/exit date
replace entrydate = exitdate_tmp if flag == 1
drop exitdate_tmp 
gen entryear=year(entrydate)
gen exityear=year(exitdate)
gen entrymiss = 1 if eintrittdatum == ""
gen exitmiss = 1 if austrittdatum == ""
tab entryear
tab exityear
replace entryear = 1993 if entrymiss == 1
replace exityear = 2018 if exitmiss == 1
tab entryear
tab exityear

replace year = . if year<entryear | year>exityear
drop if year == .
drop eintrittdatum austrittdatum

duplicates report VRid NRid duns funktion year
duplicates drop VRid NRid duns funktion year, force
sort NRid year

rename vorname firstname_bis 
rename nachname name_bis
rename wohnort residence_bis

save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_Bisnode-NRid_1994-2018.dta", replace

* Merge with NR info (elected)
use "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_BisnodeRL.dta", clear
sort NRid year
merge 1:m NRid year using "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_Bisnode-NRid_1994-2018.dta"
drop if year>2015
drop if year<1994
tab _merge
drop if elected == 0
keep if _merge == 3
drop _merge
sort NRid year uid
gen source = "bisnode"
drop if linkscore < 0.24 // F1 score = 0.24
save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_NR-Bisnode_1994-2015.dta", replace

** Step 3): merge with IG information
* prepare NR-Bisnode data
*use "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_NR_Bisnode-1994-2015.dta", clear
rename handelsregisternummer HregNr // varname in IG dataset
rename funktion function
rename firma firm
gen name_NR = firstname + " " + name
gen fullname_VR = firstname_bis + " " + name_bis
rename ortdomiziladresse firm_place
drop firstname name plzdomiziladresse kantondomiziladresse nogacode kapitalnominalag kapitaleinbezahltag clusterid linkscore entrydate exitdate flag entryear exityear entrymiss exitmiss VRid mandid firstname_bis name_bis W_PLZ4 E_CNTR_w N_CNTR_w unterschrift sex job elected list listname_bfs gdenr_2018_w handelsname uid eleyear birthyear
gen dataset = "Bisnode"
order NRid year canton name_NR gdename dataset firm firm_place Hreg function duns

* prepare NR-IG data
preserve
use "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_IG_NRid_BisnodeRL.dta", clear
drop if source == "IG" // drop Durrer & Engelberger
gen name_NR = firstname + " " + name
gen HregNr = Hreg
replace HregNr = usubinstr(HregNr, "-", "",.) // normalize nummer to reflect Bisnode-format
replace HregNr = usubinstr(HregNr, ".", "",.) 
drop issue firstname name job WohnOrt WohnAdr WohnPlz PartyS partylas Sex Mandate mp_id BirthDate sex job elected list listname_bfs gdenr_2018_w birthyear eleyear Hreg job Profession
rename Affiliation firm
rename Place firm_place
rename mp_name fullname_VR
rename Function function
gen dataset = "IG"
order NRid year canton name_NR gdename dataset firm firm_place Hreg function 
save "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_IG_NRid_BisnodeRL_tmp.dta", replace
restore

append using "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_IG_NRid_BisnodeRL_tmp.dta"

sort NRid name_NR year firm
gen coding = ""
gen PY_mandID = ""
order NRid year canton name_NR gdename dataset coding PY_mandID firm firm_place HregNr function duns source fullname_VR residence_bis
drop if year < 1996 | year > 2011


export excel "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR-Bisnode-IG_1996-2011.xlsx",  firstrow(variables) replace


*** Erase files
forv i = 1(1)7 {  // `max'
	erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_VRid_NRid_`i'.dta"
	erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_Bisnode-NRid_`i'.dta"
}
erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\IG_AG_1996-2012.dta"
erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_Bisnode-NRid.dta"
erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_VRid_NRid.dta"
erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_VRid.dta"
erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_NRid.dta"
*erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_IG_NRid_BisnodeRL.dta"
erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_IG_NRid_BisnodeRL_tmp.dta"
*erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_NR-Bisnode_1994-2015.dta"
*erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\RL_Bisnode-NRid_1994-2018.dta"
erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_BisnodeRL.dta"
erase "$dataNR\18_Interessenbindungen\01_Check_Bisnode_Interessenbindungen\NR_BisnodeRL_elected.dta"


