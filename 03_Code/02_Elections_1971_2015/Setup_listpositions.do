
* Import NRW candidate data including list positions 1975-2003 *
forv i=1975(4)2003 {
	import excel "$path\01_Raw_data\02_Elections_1971_2015\Kandidatstimmen 1975-2003.xlsx", sheet(`i') firstrow clear
	tostring *, replace
	gen year = `i'
	order year
	save "$path\02_Processed_data\02_Elections_1971_2015\Kandidatstimmen_`i'.dta", replace
}
use "$path\02_Processed_data\02_Elections_1971_2015\Kandidatstimmen_1975.dta", clear
forv i=1979(4)2003 {
	append using "$path\02_Processed_data\02_Elections_1971_2015\Kandidatstimmen_`i'.dta"
}
rename ID cand_no
rename KANTON cantonno
rename NAME_KANTON canton
rename CANDIDATE_NAME name
rename FIRST_NAME firstname
rename SEX sex
rename BIRTH birthyear
rename first_position listpos
rename CUMULATIVE  sdlistpos
rename ELECTED  elected
rename party_affiliation  party_no
rename Label_D_F  party_name
rename name_short  list_name
rename LIST_ID list_nr_off
rename VOTES_NOT_CHANGED  votes_unchanged
rename VOTES_CHANGED  votes_changed
rename VOTES_TOTAL votes
destring cantonno cand_no birthyear listpos sdlistpos party_no votes_unchanged votes_changed votes, replace
order year cantonno canton cand_no list_nr_off list_name name firstname sex birthyear listpos sdlistpos elected party_no party_name votes_unchanged votes_changed votes
sort canton year name firstname votes
save "$path\02_Processed_data\02_Elections_1971_2015\Listenpositionen_1975-2003.dta", replace

* Import NRW candidate data including list positions 2007-2015 *
import excel "$path\01_Raw_data\02_Elections_1971_2015\su-d-17.02.02.05.02.zc.2007.k.xlsx", cellrange(A5:S3138) clear
rename A year
rename B cantonno
rename C canton
rename D list_nr
rename E list_nr_off
rename F list_name
rename G party_no
rename H party_name
rename I cand_no
rename J name
rename K firstname
rename L sex
rename M birthyear
rename N listpos
rename O sdlistpos
rename P elected
rename Q votes_unchanged
rename R votes_changed
rename S votes
replace canton = "NW" if canton=="NW 5)"
label var year "Wahljahr"
label var cantonno "Kantons-Nr."
label var canton "Kanton"
label var list_nr "Listen-Nr. (numerisch)"
label var list_nr_off "Listen-Nr. (offiziell)"
label var list_name "Liste"
label var party_no "Partei-Nr."
label var party_name "Partei"
label var cand_no "Kandidaten-Nr."
label var name "Name"
label var firstname "Vorname"
label var sex "Geschlecht"
label var birthyear "Geburtsjahr"
label var listpos "Erste Position auf Wahlzettel"
label var sdlistpos "Zweite Position auf Wahlzettel (Kumulierung)"
label var elected "Gewählt"
label var votes_unchanged "Stimmen aus unveränderten Wahlzetteln"
label var votes_changed "Stimmen aus veränderten Wahlztteln"
label var votes "Total erhaltene Stimmen"
drop if canton=="" & year==. & firstname=="" & name=="" & votes==.
order year cantonno canton cand_no list_nr list_name name firstname sex birthyear listpos sdlistpos elected party_no party_name votes_unchanged votes_changed votes
sort year canton name firstname votes
save "$path\02_Processed_data\02_Elections_1971_2015\Kandidatstimmen_2007.dta", replace

import excel "$path\01_Raw_data\02_Elections_1971_2015\su-d-17.02.02.05.02.zb.2011.k.xlsx", cellrange(A5:S3510) clear
rename A year
rename B cantonno
rename C canton
rename D list_nr
rename E list_nr_off
rename F list_name
rename G party_no
rename H party_name
rename I cand_no
rename J name
rename K firstname
rename L sex
rename M birthyear
rename N listpos
rename O sdlistpos
rename P elected
rename Q votes_unchanged
rename R votes_changed
rename S votes
label var year "Wahljahr"
label var cantonno "Kantons-Nr."
label var canton "Kanton"
label var list_nr "Listen-Nr. (numerisch)"
label var list_nr_off "Listen-Nr. (offiziell)"
label var list_name "Liste"
label var party_no "Partei-Nr."
label var party_name "Partei"
label var cand_no "Kandidaten-Nr."
label var name "Name"
label var firstname "Vorname"
label var sex "Geschlecht"
label var birthyear "Geburtsjahr"
label var listpos "Erste Position auf Wahlzettel"
label var sdlistpos "Zweite Position auf Wahlzettel (Kumulierung)"
label var elected "Gewählt"
label var votes_unchanged "Stimmen aus unveränderten Wahlzetteln"
label var votes_changed "Stimmen aus veränderten Wahlztteln"
label var votes "Total erhaltene Stimmen"
drop if canton=="" & year==. & firstname=="" & name=="" & votes==.
destring listpos sdlistpos, replace 
order year cantonno canton cand_no list_nr list_name name firstname sex birthyear listpos sdlistpos elected party_no party_name votes_unchanged votes_changed votes
sort year canton name firstname votes
save "$path\02_Processed_data\02_Elections_1971_2015\Kandidatstimmen_2011.dta", replace

import excel "$path\01_Raw_data\02_Elections_1971_2015\su-d-17.02.02.05.02.za.2015.k.xlsx", cellrange(A5:S3810) clear
rename A year
rename B cantonno
rename C canton
rename D list_nr
rename E list_nr_off
rename F list_name
rename G party_no
rename H party_name
rename I cand_no
rename J name
rename K firstname
rename L sex
rename M birthyear
rename N listpos
rename O sdlistpos
rename P elected
rename Q votes_unchanged
rename R votes_changed
rename S votes
label var year "Wahljahr"
label var cantonno "Kantons-Nr."
label var canton "Kanton"
label var list_nr "Listen-Nr. (numerisch)"
label var list_nr_off "Listen-Nr. (offiziell)"
label var list_name "Liste"
label var party_no "Partei-Nr."
label var party_name "Partei"
label var cand_no "Kandidaten-Nr."
label var name "Name"
label var firstname "Vorname"
label var sex "Geschlecht"
label var birthyear "Geburtsjahr"
label var listpos "Erste Position auf Wahlzettel"
label var sdlistpos "Zweite Position auf Wahlzettel (Kumulierung)"
label var elected "Gewählt"
label var votes_unchanged "Stimmen aus unveränderten Wahlzetteln"
label var votes_changed "Stimmen aus veränderten Wahlztteln"
label var votes "Total erhaltene Stimmen"
drop if canton=="" & year==. & firstname=="" & name=="" & votes==.
destring listpos sdlistpos votes_changed, replace 
order year cantonno canton cand_no list_nr_off list_nr list_name name firstname sex birthyear listpos sdlistpos elected party_no party_name votes_unchanged votes_changed votes
sort year canton name firstname votes
save  "$path\02_Processed_data\02_Elections_1971_2015\Kandidatstimmen_2015.dta", replace

use "$path\02_Processed_data\02_Elections_1971_2015\Listenpositionen_1975-2003.dta", clear
forv i=2007(4)2015 {
	append using "$path\02_Processed_data\02_Elections_1971_2015\Kandidatstimmen_`i'.dta"
}
rename cand_no candidateno
rename list_nr_off listno
rename list_name list_listpos
rename name name_listpos
rename firstname firstname_listpos
rename party_name party_listpos
rename votes votes_listpos
label var list_nr "Numerical list number (NOT official)"
label var listpos "Position on list"
label var sdlistpos "Second position on list"
label var list_listpos "List name BfS"

// Adjust listno to match info in candidate data
replace listno = "1" if list_listpos == "Arthur Loepfe (CVP)" & canton == "AI" & year == 2007
replace listno = "1" if list_listpos == "Marianne Kleiner-S. (FDP)" & canton == "AR" & year == 2007
replace listno = "2" if list_listpos == "Edgar Bischof (Übrige)" & canton == "AR" & year == 2007
replace listno = "3" if list_listpos == "Ivo Müller (Übrige)" & canton == "AR" & year == 2007
replace listno = "4" if list_listpos == "M. Weisshaupt (Übrige)" & canton == "AR" & year == 2007
replace listno = "5" if list_listpos == "Jakob Freund (Übrige)" & canton == "AR" & year == 2007
replace listno = "1" if list_listpos == "Werner Marti (SP)" & canton == "GL" & year == 2007
replace listno = "2" if list_listpos == "Martin Dürst (Junge SVP)" & canton == "GL" & year == 2007
replace listno = "3" if list_listpos == "Julius Fäh (Übrige)" & canton == "GL" & year == 2007
replace listno = "1" if list_listpos == "Edi Engelberger (FDP)" & canton == "NW" & year == 2007
replace listno = "1" if list_listpos == "von Rotz-Spichtig (SVP)" & canton == "OW" & year == 2007
replace listno = "2" if list_listpos == "Beat von Wyl (SP)" & canton == "OW" & year == 2007
replace listno = "3" if list_listpos == "Lukas Gasser (Übrige)" & canton == "OW" & year == 2007
replace listno = "4" if list_listpos == "Patrick Imfeld (CVP)" & canton == "OW" & year == 2007
replace listno = "1" if list_listpos == "Gabi Huber (FDP)" & canton == "UR" & year == 2007

replace candidateno = 1 if list_listpos == "Arthur Loepfe (CVP)" & canton == "AI" & year == 2007
replace candidateno = 1 if list_listpos == "Marianne Kleiner-S. (FDP)" & canton == "AR" & year == 2007
replace candidateno = 1 if list_listpos == "Edgar Bischof (Übrige)" & canton == "AR" & year == 2007
replace candidateno = 1 if list_listpos == "Ivo Müller (Übrige)" & canton == "AR" & year == 2007
replace candidateno = 1 if list_listpos == "M. Weisshaupt (Übrige)" & canton == "AR" & year == 2007
replace candidateno = 1 if list_listpos == "Jakob Freund (Übrige)" & canton == "AR" & year == 2007
replace candidateno = 1 if list_listpos == "Werner Marti (SP)" & canton == "GL" & year == 2007
replace candidateno = 1 if list_listpos == "Martin Dürst (Junge SVP)" & canton == "GL" & year == 2007
replace candidateno = 1 if list_listpos == "Julius Fäh (Übrige)" & canton == "GL" & year == 2007
replace candidateno = 1 if list_listpos == "Edi Engelberger (FDP)" & canton == "NW" & year == 2007
replace candidateno = 1 if list_listpos == "von Rotz-Spichtig (SVP)" & canton == "OW" & year == 2007
replace candidateno = 1 if list_listpos == "Beat von Wyl (SP)" & canton == "OW" & year == 2007
replace candidateno = 1 if list_listpos == "Lukas Gasser (Übrige)" & canton == "OW" & year == 2007
replace candidateno = 1 if list_listpos == "Patrick Imfeld (CVP)" & canton == "OW" & year == 2007
replace candidateno = 1 if list_listpos == "Gabi Huber (FDP)" & canton == "UR" & year == 2007

replace listno = "2" if list_listpos == "Vereinzelte" & canton == "AI" & year == 2007
replace listno = "6" if list_listpos == "Vereinzelte" & canton == "AR" & year == 2007
replace listno = "4" if list_listpos == "Vereinzelte" & canton == "GL" & year == 2007
replace listno = "5" if list_listpos == "Vereinzelte" & canton == "OW" & year == 2007
replace listno = "missing" if list_listpos == "Vereinzelte" & canton == "UR" & year == 2007

replace listno = upper(listno)

keep year canton candidateno listno listpos list_nr name_listpos firstname_listpos
order year canton candidateno listno listpos list_nr
sort canton year listno candidateno
save "$path\02_Processed_data\02_Elections_1971_2015\Listenpositionen_1975-2015.dta", replace

* erase yearly data
erase "$path\02_Processed_data\02_Elections_1971_2015\Listenpositionen_1975-2003.dta"
forv i=1975(4)2015 {
	erase "$path\02_Processed_data\02_Elections_1971_2015\Kandidatstimmen_`i'.dta"
}


