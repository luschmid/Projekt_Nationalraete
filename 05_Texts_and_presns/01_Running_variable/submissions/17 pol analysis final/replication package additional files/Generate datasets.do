global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"

*---------------------------
* A) Dataset for Switzerland
*---------------------------

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear

drop ID_pers
egen ID_pers=group(ID)

* (i) Generate running variable, election in next election,
* 	   party variables, and cantonal dummy variable

xtset ID_pers ID_time
sort ID_pers ID_time


gen party_sp= listname_bfs=="SP/PS" if year>=1971 
gen party_cvp= listname_bfs=="CVP/PDC" if year>=1971 
gen party_fdp= listname_bfs=="FDP/PLR (PRD)" if year>=1971 
gen party_svp= listname_bfs=="SVP/UDC" if year>=1971 
replace party_svp=1 if listname_bfs=="BDP/PBD" & year>=1971 

tab canton, gen(canton_)

* (ii) Variable labeling

label var ID_pers "Person ID"
label var no_seats "Number of seats in canton"
label var votes "Candidate votes"
label var eligible_cant "Number of eligible voters"
label var firstname "Firstname"
label var name "Name"

label var canton_1 "Aargau"
label var canton_2 "Appenzell Innerrhoden"
label var canton_3 "Appenzell Ausserrhoden"
label var canton_4 "Bern"
label var canton_5 "Basel Landschaft"
label var canton_6 "Basel Stadt"
label var canton_7 "Fribourg"
label var canton_8 "Geneva"
label var canton_9 "Glarus"
label var canton_10 "Graubünden"
label var canton_11 "Jura"
label var canton_12 "Lucerne"
label var canton_13 "Neuchâtel"
label var canton_14 "Nidwalden"
label var canton_15 "Obwalden"
label var canton_16 "St. Gallen"
label var canton_17 "Schaffhausen"
label var canton_18 "Solothurn"
label var canton_19 "Schwyz"
label var canton_20 "Thurgau"
label var canton_21 "Ticino"
label var canton_22 "Uri"
label var canton_23 "Vaud"
label var canton_24 "Valais"
label var canton_25 "Zug"
label var canton_26 "Zurich"

label var party_sp "Social democrats (SP)"
label var party_cvp "Christian democrats (CVP)"
label var party_fdp "Liberals (FDP)"
label var party_svp "National convervatives (SVP)"

* (iii) Keep only relevant variables

keep ID_pers ID_time canton firstname name year alliance suballiance ///
 list suballiance votes pvotes elected eligible_cant eligible_cant canton canton_* ///
 party_* age sex year no_seats
 
order ID_pers ID_time canton year firstname name votes pvotes elected list alliance suballiance ///
  eligible_cant no_seats age sex canton_* ///
 party_*  

save "$path\05_Texts_and_presns\01_Running_variable\submissions\17 pol analysis final\replication package\data_switzerland", replace


*-------------------------
* B) Dataset for Honduras
*-------------------------

use "$path\02_Processed_data\15_Elections_Honduras\elections_hn_final.dta", clear

* (i) All vars to lower case

rename *, lower

* (ii) Generate numeric id for candidates and numeric votes variable

egen ID_pers=group(id_stata)
destring votes, replace
egen ID_time=group(year)

* (iii) Generate variable for party votes, party dummies, and
* 		department dummies

bysort department party year: egen pvotes=sum(votes)

gen party_pn= party=="PARTIDO NACIONAL DE HONDURAS" 
gen party_pl= party=="PARTIDO LIBERAL DE HONDURAS" 
gen party_libre= party=="PARTIDO LIBERTAD Y REFUNDACION" 

tab department, gen(department_)

gen sex_string=sex
drop sex

gen sex=.
replace sex=1 if sex_string=="H"
replace sex=0 if sex_string=="M"
drop sex_string

* (iv) Variable labeling

label var candidate "Candidate"
label var department "Department"
label var ID_pers "Person ID"
label var elected "Elected"
label var party "Party"
label var votes "Candidate votes"
label var pvotes "Party votes"
label var total_electoral "Number of eligible voters"
label var no_seats "Number of seats in department"
label var year "Year"
label var sex "Sex"

label var department_1 "ATLANTIDA"
label var department_2 "CHOLUTECA"
label var department_3 "COLON"
label var department_4 "COMAYAGUA"
label var department_5 "COPAN"
label var department_6 "CORTES"
label var department_7 "EL PARAISO"
label var department_8 "FRANCISCO MORAZAN"
label var department_9 "GRACIAS A DIOS"
label var department_10 "INTIBUCA"
label var department_11 "ISLAS DE LA BAHIA "
label var department_12 "LA PAZ"
label var department_13 "LEMPIRA"
label var department_14 "OCOTEPEQUE"
label var department_15 "OLANCHO"
label var department_16 "SANTA BARBARA"
label var department_17 "VALLE"
label var department_18 "YORO"

label var party_pn "National conservatives (PN)"
label var party_pl "Social liberals (PL)"
label var party_libre "Social democrats (Libre)"

* (v) Keep only relevant variables

keep ID_pers ID_time candidate  ///
department year sex party votes pvotes elected total_electoral no_seats department_* party_*

order ID_pers ID_time department year candidate votes pvotes elected party total_electoral ///
no_seats sex department_* party_*

save "$path\05_Texts_and_presns\01_Running_variable\submissions\17 pol analysis final\replication package\data_honduras.dta", replace

*-------------------------
* C) Dataset for Norway
*-------------------------

use "$path\01_Raw_data\13_Running_variable\NSD2405_en_F7.dta", clear

rename *, lower

* (i) Keep only relevant years

keep if year>=1953 & year<=1985

* Note: We keep 1985 b/c we use this data as outcome data. We drop it before
*		saving the data (see below in this file). 

* (ii) Small correction of data

replace pid=22893 if candidatename_orig=="Lars Fagerland" & year==1973
* Note: The original pid refers to the pid of another person. 
*		Fiva and Smith (2018) made a mistake by editing the original candidate 
*		name. It changed from Lars Fagerland to Styrk Lothe. We have to make 
*		this change, otherwise the candidates run in the same district and year 
*		in different parties



* (ii) Generate party dummies and district dummies

rename margin votemargin_fiva_smith
rename votes pvotes

tab party, gen(party_)
tab district, gen(district_)


* (iii) Generate unique pid
* Note: We do not know how and when Jon drops duplicate pids. Instead, 
* 		we use our way and assign a unique pid per district

gen double ID_pers=districtid*100000+pid
egen ID_time=group(year)
xtset ID_pers ID_time


* (iv) Create dummy for Fiva/Smith sample

gen marginal_candidate=""
replace marginal_candidate="Marginal candidate" if votemargin_fiva_smith!=.
replace marginal_candidate="Non-marginal candidate" if votemargin_fiva_smith==.

* (v) Keep only relevant variables 

label var ID_pers "Person ID"
label var districtid "District ID"
label var elected "Elected"
label var party "Party"
label var pvotes "Party votes"
label var electorate "Number of eligible voters"
label var firstname "Firstname"
label var lastname "Name"
label var rank "List position"
label var marginal_candidate "Marginal candidate"

keep ID_pers ID_time year districtid firstname lastname rank pvotes elected party electorate ///
 marginal_candidate
 
order ID_pers ID_time year districtid firstname lastname rank pvotes elected party electorate ///
 marginal_candidate 

 
* (vi) Keep only relevant years for estimation

save "$path\05_Texts_and_presns\01_Running_variable\submissions\17 pol analysis final\replication package\data_norway.dta", replace
