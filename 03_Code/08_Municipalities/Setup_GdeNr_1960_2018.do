

****************************************************************
*** Create municipality panel: Take account of amalgamations ***
****************************************************************


*	Liste historique des communes
*	https://www.bfs.admin.ch/bfs/fr/home/bases-statistiques/repertoire-officiel-communes-suisse/liste-historisee-communes.assetdetail.4062823.html

import excel "$path/01_Raw_data/08_Municipalities/BfS_180829_Gde_1960_2018.xlsx", sheet("Liste") firstrow
* Add canton number
gen ctnnr = 1 if ctn=="ZH" 
replace ctnnr = 2 if ctn=="BE"
replace ctnnr = 3 if ctn=="LU"
replace ctnnr = 4 if ctn=="UR"
replace ctnnr = 5 if ctn=="SZ"
replace ctnnr = 6 if ctn=="OW"
replace ctnnr = 7 if ctn=="NW"
replace ctnnr = 8 if ctn=="GL"
replace ctnnr = 9 if ctn=="ZG"
replace ctnnr = 10 if ctn=="FR"
replace ctnnr = 11 if ctn=="SO"
replace ctnnr = 12 if ctn=="BS"
replace ctnnr = 13 if ctn=="BL"
replace ctnnr = 14 if ctn=="SH"
replace ctnnr = 15 if ctn=="AR"
replace ctnnr = 16 if ctn=="AI"
replace ctnnr = 17 if ctn=="SG"
replace ctnnr = 18 if ctn=="GR"
replace ctnnr = 19 if ctn=="AG"
replace ctnnr = 20 if ctn=="TG"
replace ctnnr = 21 if ctn=="TI"
replace ctnnr = 22 if ctn=="VD"
replace ctnnr = 23 if ctn=="VS"
replace ctnnr = 24 if ctn=="NE"
replace ctnnr = 25 if ctn=="GE"
replace ctnnr = 26 if ctn=="JU"

* Variable labels 
label var histnr "Gde historic id"
label var district_histnr "District historic id"
label var ctn "Canton"
label var ctnnr "Canton id"
label var gdenr "Gde BSF id"
label var gdename "Gde name"
label var gdename_sht "Gde name short"
label var type "Recording type 11 = commune politique, 12 = territoire non attribuÈ; 13 = partie cantonale de lac"
label var status "provis = 0; definitive = 1"
label var entry_mutnr "Entry id"
label var entry_type "Entry type (20-27)"
label var entry_yr "Entry year"
label var exit_mutnr "Exit id"
label var exit_type "Exit type (22-30)"
label var exit_yr "Exit year"
label var modif_yr "Modification year"

order histnr district_histnr ctn ctnnr gdenr gdename gdename_sht type status entry_mutnr ///
entry_type entry_yr exit_mutnr exit_type exit_yr modif_yr

* What type of municipality? political municipality => type == 11
keep if type == 11
drop type

* Only municipalities in final status (no "provisorisch") => status == 1
keep if status == 1
drop status

sort histnr

save "$path/02_Processed_data/08_Municipalities/municipalities_list_BfS.dta", replace


*	Construct the universe of existing municipalities on 31.12. of every year
********************************************************************************
* Produce yearly universe of municipalities

clear
quietly foreach i of num 1960/2018 {
clear
import excel "$path/01_Raw_data/08_Municipalities/`i'.xls", sheet("Data") firstrow

rename Numérodhistorisation histnr
rename Canton 				ctn
rename Numérodudistrict 	district_histnr
rename Nomdudistrict 		district
rename Numérodelacommune	gdenr 
rename Nomdelacommune 		gdename
rename Datedelinscription	entry_yr

gen year = `i'

save "$path/02_Processed_data/08_Municipalities/Gde_in_`i'.dta", replace
}

clear
use "$path/02_Processed_data/08_Municipalities/Gde_in_1960"
quietly foreach i of num 1961/2018 {
append using "$path/02_Processed_data/08_Municipalities/Gde_in_`i'.dta", force
}

sort gdenr year


gen ctnnr = 1 if ctn=="ZH" 
replace ctnnr = 2 if ctn=="BE"
replace ctnnr = 3 if ctn=="LU"
replace ctnnr = 4 if ctn=="UR"
replace ctnnr = 5 if ctn=="SZ"
replace ctnnr = 6 if ctn=="OW"
replace ctnnr = 7 if ctn=="NW"
replace ctnnr = 8 if ctn=="GL"
replace ctnnr = 9 if ctn=="ZG"
replace ctnnr = 10 if ctn=="FR"
replace ctnnr = 11 if ctn=="SO"
replace ctnnr = 12 if ctn=="BS"
replace ctnnr = 13 if ctn=="BL"
replace ctnnr = 14 if ctn=="SH"
replace ctnnr = 15 if ctn=="AR"
replace ctnnr = 16 if ctn=="AI"
replace ctnnr = 17 if ctn=="SG"
replace ctnnr = 18 if ctn=="GR"
replace ctnnr = 19 if ctn=="AG"
replace ctnnr = 20 if ctn=="TG"
replace ctnnr = 21 if ctn=="TI"
replace ctnnr = 22 if ctn=="VD"
replace ctnnr = 23 if ctn=="VS"
replace ctnnr = 24 if ctn=="NE"
replace ctnnr = 25 if ctn=="GE"
replace ctnnr = 26 if ctn=="JU"

order histnr gdenr gdename year ctnnr ctn

label var histnr "Gde historic id"
label var district_histnr "District historic id"
label var district "District"
label var ctn "Canton"
label var ctnnr "Canton id"
label var gdenr "Gde BSF id"
label var gdename "Gde name"
label var entry_yr "Entry year"
label var year "year"

save "$path/02_Processed_data/08_Municipalities/municipalities_1960_2018.dta", replace


* Merge data with data from BfS
clear
use "$path/02_Processed_data/08_Municipalities/municipalities_1960_2018.dta"

merge m:1 histnr using "$path/02_Processed_data/08_Municipalities/municipalities_list_BfS.dta"
sort histnr year

* some municipalities do not merge: 
/*
histnr	gdenr	gdename			year	ctnnr	ctn	district_histnr	district	entry_yr	gdename_sht			entry_mutnr	entry_type	exit_mutnr	exit_type	exit_yr			modif_yr	_merge
14538	5451	Avenches				22		VD	10052						01jul2006	Avenches			2337			26		2371		24			8/31/2006		8/31/2006	using only (2)
14539	5830	Villarzel				22		VD	10041						01jul2006	Villarzel			2338			26		2649		24			8/31/2006		8/31/2006	using only (2)
16088	6417	La Grande-BÈroche		24		NE	10320						01jan2018	La Grande-BÈroche	3590			21		3623		23			3/31/2018		3/31/2018	using only (2)

These municipalities changed inbetween 31.12.n and 31.12.n+1
*/
tab _merge
drop _merge

duplicates report

duplicates report histnr year

duplicates report gdenr year
duplicates list gdenr gdename year, sepby(gdenr) 
/*
  +---------------------------------------------------+
  | group:     obs:   gdenr            gdename   year |
  |---------------------------------------------------|
  |      1    64432    4471       Bischofszell   1995 |
  |      1   154674    4471       Bischofszell   1995 |
  |      1   154705    4471       Bischofszell   1995 |
  |---------------------------------------------------|
  |      2    76274    4486          Gottshaus   1995 |
  |      2   154673    4486          Gottshaus   1995 |
  |---------------------------------------------------|
  |      3    13399    4555          Wittenwil   1995 |
  |      3   154597    4555          Wittenwil   1995 |
  |---------------------------------------------------|
  |      4    72570    4566         Frauenfeld   1997 |
  |      4   155067    4566         Frauenfeld   1997 |
  |---------------------------------------------------|
  |      5    78561    4571           Gachnang   1997 |
  |      5   155068    4571           Gachnang   1997 |
  |---------------------------------------------------|
  |      6    98965    4915   Opfershofen (TG)   1994 |
  |      6   154337    4915   Opfershofen (TG)   1994 |
  |---------------------------------------------------|
  |      7    17416    4946         Weinfelden   1994 |
  |      7   154483    4946         Weinfelden   1994 |
  +---------------------------------------------------+

*/

save "$path/02_Processed_data/08_Municipalities/municipalities_UBpanel_1960_2018.dta", replace


clear
use "$path/02_Processed_data/08_Municipalities/municipalities_UBpanel_1960_2018.dta"

br if exit_mutnr==. // 2222 municipalities
sum if exit_mutnr==. & year==2018 // 2222 municipalities
gen histnr_2018 = histnr  if exit_mutnr==.
gen gdenr_2018 = gdenr  if exit_mutnr==.
gen ctnnr_2018 = ctnnr  if exit_mutnr==.

egen float id = group(histnr_2018) //  1 to 2222
label var id "id for municipalities in 2018"


* Process:
*	For every municipality existing in 2018
*		take historical number
*		take entry number
*		if entry number = exit number of other municipality, put historical number in "histnr_2018"

* Remarque: histnr_2018 correspond pour chaque commune au numÈro historisÈ en vigueur en 2018 (derniËre annÈe)
* Cette variable permettra de faire le collapse, i.e. de fusionner aritificiellement les communes qui disparaissent
* durant la pÈriode.

sum id //  1       2222

quietly forvalues i = 1/2222 { 	

sum histnr if id == `i' 
scalar hist2018 = r(mean)
sum gdenr if id == `i' 
scalar gdenr2018 = r(mean)
sum ctnnr if id == `i' 
scalar ctnnr2018 = r(mean)

sum entry_mutnr if id == `i'
scalar mut2018 = r(mean)

replace histnr_2018 = hist2018 	if exit_mutnr == mut2018 
replace gdenr_2018	= gdenr2018 if exit_mutnr == mut2018 
replace ctnnr_2018	= ctnnr2018 if exit_mutnr == mut2018 
}
scalar drop _all

* Generate an id_B for all municipality that were linked in first round
egen float id_B = group(histnr) if histnr_2018 != . & id==. //    1       1852
label var id_B "id for municipalities that are linked to 2018 after 1st iteration"
sum id_B

* 2nd iteration
quietly forvalues i = 1/1852 { 	
sum histnr_2018 if id_B == `i' 
scalar hist2018 = r(mean)
sum gdenr_2018 if id_B == `i' 
scalar gdenr2018 = r(mean)
sum ctnnr_2018 if id_B == `i' 
scalar ctnnr2018 = r(mean)

sum entry_mutnr if id_B == `i'
scalar mut2018 = r(mean)

replace histnr_2018 = hist2018 	if exit_mutnr == mut2018 
replace gdenr_2018	= gdenr2018 if exit_mutnr == mut2018
replace ctnnr_2018	= ctnnr2018 if exit_mutnr == mut2018
}
scalar drop _all

egen float id_C = group(histnr) if histnr_2018 != . & id==. & id_B==. //     1       1073
label var id_C "id for municipalities that are linked to 2018 after 2nd iteration"
sum id_C

* 3rd iteration
quietly forvalues i = 1/1073 { 	
sum histnr_2018 if id_C == `i' 
scalar hist2018 = r(mean)
sum gdenr_2018 if id_C == `i' 
scalar gdenr2018 = r(mean)
sum ctnnr_2018 if id_C == `i' 
scalar ctnnr2018 = r(mean)

sum entry_mutnr if id_C == `i'
scalar mut2018 = r(mean)

replace histnr_2018 = hist2018 	if exit_mutnr == mut2018 
replace gdenr_2018	= gdenr2018 if exit_mutnr == mut2018 
replace ctnnr_2018	= ctnnr2018 if exit_mutnr == mut2018 
}
scalar drop _all

egen float id_D = group(histnr) if histnr_2018 != .  & id==. & id_B==. & id_C==. //    1        475
label var id_D "id for municipalities that are linked to 2018 after 3rd iteration"
sum id_D

* 4th iteration
quietly forvalues i = 1/475 { 	
sum histnr_2018 if id_D == `i' 
scalar hist2018 = r(mean)
sum gdenr_2018 if id_D == `i' 
scalar gdenr2018 = r(mean)
sum ctnnr_2018 if id_D == `i' 
scalar ctnnr2018 = r(mean)

sum entry_mutnr if id_D == `i'
scalar mut2018 = r(mean)

replace histnr_2018 = hist2018 	if exit_mutnr == mut2018 
replace gdenr_2018	= gdenr2018 if exit_mutnr == mut2018 
replace ctnnr_2018	= ctnnr2018 if exit_mutnr == mut2018
}
scalar drop _all

egen float id_E = group(histnr) if histnr_2018 != .  & id==. & id_B==. & id_C==. & id_D==.  //	1         95
label var id_E "id for municipalities that are linked to 2018 after 4th iteration"
sum id_E

* 5th iteration
quietly forvalues i = 1/95 { 	
sum histnr_2018 if id_E == `i' 
scalar hist2018 = r(mean)
sum gdenr_2018 if id_E == `i' 
scalar gdenr2018 = r(mean)
sum ctnnr_2018 if id_E == `i' 
scalar ctnnr2018 = r(mean)

sum entry_mutnr if id_E == `i'
scalar mut2018 = r(mean)

replace histnr_2018 = hist2018 	if exit_mutnr == mut2018 
replace gdenr_2018	= gdenr2018 if exit_mutnr == mut2018 
replace ctnnr_2018	= ctnnr2018 if exit_mutnr == mut2018
}
scalar drop _all

egen float id_F = group(histnr) if histnr_2018 != .  & id==. & id_B==. & id_C==. & id_D==. & id_E==. //  
label var id_F "id for municipalities that are linked to 2018 after 5th iteration"
sum id_F

sum histnr histnr_2018
sum if histnr_2018==.

sort gdenr year

drop id id_B id_C id_D id_E id_F

label var histnr_2018 	"Historic id in 2018"
label var gdenr_2018 	"Gde BSF id in 2018"
label var ctnnr_2018	"Ctn id in 2018"

order histnr gdenr gdename gdename_sht

* Correction of BfS GemeindeNr due to exchanges of teritory between municipalities
drop if year == 1994 & gdename == "Berg (TG)" & gdenr == 4892 & gdenr_2018 == 4946 // 2 obs for Berg in 1994; this one is wrong.
replace gdenr_2018 = 4891 if gdenr_2018 == 4946 & gdename == "Berg (TG)" // Berg swapped territory with Weinfelden (4946); Berg is correct.
replace gdenr_2018 = 4891 if gdenr_2018 == 4946 & gdename == "Andhausen" // Andhausen merged with Berg (4891), Berg later swapped territory with Weinfelden (4946); Berg is correct.
replace gdenr_2018 = 4891 if gdenr_2018 == 4946 & gdename == "Graltshausen" // Ditto for Graltshausen
replace gdenr_2018 = 4891 if gdenr_2018 == 4946 & gdename == "Guntershausen bei Birwinken" // Ditto for Guntershausen bei Birwinken
replace gdenr_2018 = 4891 if gdenr_2018 == 4946 & gdename == "Mauren" // Ditto for Mauren
replace gdenr_2018 = 4781 if gdenr_2018 == 4551 & gdename == "Wängi" // Unclear how Wängi is connected to Aadrof (4551); connection seems to run from Wängi to Wittenwil to Aadorf; number for 2018 should be 4781 throughout.
replace gdenr_2018 = 4501 if gdenr_2018 == 4471 & gdename == "Buhwil" // Buhwil merged to Kradolf-Schönenberg (4501), Kradolf-Schönenberg later swapped territory with Bischofszell (4471); Kradolf-Schönenberg is correct.
replace gdenr_2018 = 4501 if gdenr_2018 == 4471 & gdename == "Kradolf" // Ditto for Kradolf
replace gdenr_2018 = 4501 if gdenr_2018 == 4471 & gdename == "Neukirch an der Thur" // Ditto for Neukirch an der Thur
replace gdenr_2018 = 4501 if gdenr_2018 == 4471 & gdename == "Schönenberg an der Thur" // Ditto for Schönenberg an der Thur
replace gdenr_2018 = 4911 if gdenr_2018 == 4506 & gdename == "Opfershofen (TG)" // Opfershofen (TG) is related to Sulgen (4506) through mergers or territorial swaps, but in the end it seems to have merged with Bürglen (TG) (4911)
replace gdenr_2018 = 4761 if gdenr_2018 == 4724 & gdename == "Sirnach" // Sirnach (4761) and Eschlikon (4724) swapped territory in 1997, but Sirnach should be coded as Sirnach.
replace gdenr_2018 = 4761 if gdenr_2018 == 4724 & gdename == "Busswil (TG)" // Busswil (TG) merged to Sirnach (4761); Sirnach and Eschlikon (4724) later swapped territory; Sirnach is correct.
replace gdenr_2018 = 4761 if gdenr_2018 == 4724 & gdename == "Horben" // Ditto for Horben
replace gdenr_2018 = 4761 if gdenr_2018 == 4724 & gdename == "Wiezikon" // Ditto for Wiezikon
replace gdenr_2018 = 4571 if gdenr_2018 == 4566 & gdename == "Oberwil (TG)" // Oberwil (TG), Niederwil (TG), Kefikon, Islikon, and Gachnang merged to Gachnang (4571), Gachnang and Frauenfeld later swapped territory; Gachnang is correct.
replace gdenr_2018 = 4571 if gdenr_2018 == 4566 & gdename == "Niederwil (TG)" // Ditto for Niederwil (TG)
replace gdenr_2018 = 4571 if gdenr_2018 == 4566 & gdename == "Kefikon" // Ditto for Kefikon
replace gdenr_2018 = 4571 if gdenr_2018 == 4566 & gdename == "Islikon" // Ditto for Islikon
replace gdenr_2018 = 4571 if gdenr_2018 == 4566 & gdename == "Gachnang" // Ditto for Gachnang
replace gdenr_2018 = 4486 if gdenr_2018 == 4471 & gdename == "Gottshaus" // Gottshaus swapped territory with Bischofszell (4471), but should be coded as Hauptwil-Gottshaus (4486)
replace gdenr_2018 = 4683 if gdenr_2018 == 4891 & gdename == "Lengwil" // Illighausen and Oberhofen bei Kreuzlingen merged to Lengwil (4683); number for 2018 should be 4683 throughout.
replace gdenr_2018 = 4683 if gdenr_2018 == 4891 & gdename == "Illighausen" // Ditto for Illighausen 
replace gdenr_2018 = 4683 if gdenr_2018 == 4891 & gdename == "Oberhofen bei Kreuzlingen" // Ditto for Oberhofen bei Kreuzlingen
replace gdenr_2018 = 4946 if gdenr_2018 == 4891 & gdename == "Weinfelden" // Weinfelden (4946) swapped territory with Berg (TG) (4891) in 1995; Weinfelden is correct.
replace gdenr_2018 = 4061 if gdenr_2018 == 4084 & gdename == "Arni-Islisberg" // Arni-Islisberg separated in 1983; Arni (AG) seems to be the more important secceeding municipality with around 3 times as many inhabitants as Islisberg; Arni's number for 2018 used.
replace gdenr_2018 = 1081 if gdenr_2018 == 1099 & gdename == "Beromünster" // Beromünster (1081) and Schenkon (1099) swapped territory in 2015, but Beromünster should be coded as Beromünster throughout.
replace gdenr_2018 = 1081 if gdenr_2018 == 1099 & gdename == "Gunzwil" // Gunzwil, Neudorf, and Schwarzenbach merged with Beromünster (1081), which later swapped territory with Schenkon (1099); Beromünster is correct.
replace gdenr_2018 = 1081 if gdenr_2018 == 1099 & gdename == "Neudorf" // Ditto for Neudorf
replace gdenr_2018 = 1081 if gdenr_2018 == 1099 & gdename == "Schwarzenbach" // Ditto for Neudorf
replace gdenr_2018 = 24 if gdenr_2018 == 223 & gdename == "Buch am Irchel" // Buch am Irchel (24) and Neftenbach (223) swapped territory in 2013, but Buch am Irchel should be coded as Buch am Irchel throughout.
replace gdenr_2018 = 5249 if gdenr_2018 == 5269 & gdename == "Casima" // Casima, Monte, and Castel San Pietro as well as one village of Caneggio merged to Castel San Pietro (5249); number in 2018 should therefore be 5249, not 5269 (Breggia) to which Caneggio later merged.
replace gdenr_2018 = 5249 if gdenr_2018 == 5269 & gdename == "Monte" // Ditto for Monte
replace gdenr_2018 = 5249 if gdenr_2018 == 5269 & gdename == "Castel San Pietro" // Ditto for Castel San Pietro
replace gdenr_2018 = 3951 if gdenr_2018 == 3953 & gdename == "Fläsch" // Common areas of Fläsch and Maienfeld are merged to the respective communities in 1977, but number for 2018 for Fläsch should be 3951 throughout.
replace gdenr_2018 = 623 if gdenr_2018 == 616 & gdename == "Rubigen" // Rubigen splits in 1993 into Trimstein (631), Allmendingen (630), and Rubigen (623), Trimstein later merged to Münsingen (616); Rubigen was by far the largest of the three successor municipalities with population figures around 5 times the ones of the other two municipalities.
replace gdenr_2018 = 608 if gdenr_2018 == 629 & gdename == "Schlosswil" // Oberhünigen (629) is separated from Schlosswil in 1980, but Schlosswil is the larger of the successor communities with arount twice the number of inhabitants; Schlosswil merged to Grosshöchstetten (608) in 2018.



sort year gdenr_2018

save "$path/02_Processed_data/08_Municipalities/Gde_1960_2018", replace



* Erase intermediate processed data
erase "$path/02_Processed_data/08_Municipalities/municipalities_1960_2018.dta"
erase "$path/02_Processed_data/08_Municipalities/municipalities_list_BfS.dta"
erase "$path/02_Processed_data/08_Municipalities/municipalities_UBpanel_1960_2018.dta"
forv v = 1960(1)2018 {
erase "$path/02_Processed_data/08_Municipalities/Gde_in_`v'.dta"
}
