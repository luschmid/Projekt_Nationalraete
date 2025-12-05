******************************************************
* (A) Set directories
******************************************************

clear all
set more off

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
global outpath "$path\02_Processed_data\15_Elections_Honduras"

******************************************************
* (B) READ-IN ELECTION RESULTS DATA (2009-2017)
******************************************************

* 1. Read in election results for 2009

global inpath "$path\01_Raw_data\15_Elections_Honduras\02_Elecciones_Generales_2009\"
program readin1
	import exc using "$inpath\04_Detalle de Votos\\`1'_`2'.xlsx", locale(es_HN) ///
	clear
	if A[1] != "Resultado detallado de votos por candidato" exit 688
	qui sum B
	if `r(sum)' != `3' {
		di "Vote total `r(sum)' instead of `3'"
		exit 688
	}
	* 2 tests to check if the files are corrupted
	g Year = 2009
	g Department = A[3]
	g Party = A[2]
	g Party_abbr = "`1'"
	drop if [_n] < 4
	rename A Candidate
	rename B Votes
	destring C, g(Elected)
	destring D, g(ID_TSE)
	rename E Sex
	if "`1'" == "PN" & "`2'" == "ATLANTIDA" save "$outpath\elections_hn.dta", ///
	replace
	else {
	append using "$outpath\elections_hn.dta"
	save "$outpath\elections_hn.dta", replace
	}
end

readin1	PN		ATLANTIDA			277121
readin1	PN		COLON				86608
readin1	PN		COMAYAGUA			253374
readin1	PN		COPAN				265176
readin1	PN		CORTES				2078565
readin1	PN		CHOLUTECA			395356
readin1	PN		EL_PARAISO			264504
readin1	PN		FRANCISCO_MORAZAN	3448399
readin1	PN		GRACIAS_A_DIOS		8022
readin1	PN		INTIBUCA			71860
readin1	PN		ISLAS_DE_LA_BAHIA	5315
readin1	PN		LA_PAZ				57931
readin1	PN		LEMPIRA				214400
readin1	PN		OCOTEPEQUE			28588
readin1	PN		OLANCHO				296367
readin1	PN		SANTA_BARBARA		409076
readin1	PN		VALLE				77177
readin1	PN		YORO				323738
readin1	PINU	ATLANTIDA			21497
readin1	PINU	COLON				4761
readin1	PINU	COMAYAGUA			21550
readin1	PINU	COPAN				14588
readin1	PINU	CORTES				234755
readin1	PINU	CHOLUTECA			37762
readin1	PINU	EL_PARAISO			15118
readin1	PINU	FRANCISCO_MORAZAN	594942
readin1	PINU	INTIBUCA			3452
readin1	PINU	ISLAS_DE_LA_BAHIA	109
readin1	PINU	LA_PAZ				3446
readin1	PINU	LEMPIRA				9629
readin1	PINU	OCOTEPEQUE			1551
readin1	PINU	OLANCHO				23814
readin1	PINU	SANTA_BARBARA		17711
readin1	PINU	VALLE				3332
readin1	PINU	YORO				23201
readin1	PDCH	ATLANTIDA			22260
readin1	PDCH	COLON				2800
readin1	PDCH	COMAYAGUA			11997
readin1	PDCH	COPAN				8099
readin1	PDCH	CORTES				199071
readin1	PDCH	CHOLUTECA			35643
readin1	PDCH	EL_PARAISO			11624
readin1	PDCH	FRANCISCO_MORAZAN	403736
readin1	PDCH	GRACIAS_A_DIOS		36
readin1	PDCH	INTIBUCA			960
readin1	PDCH	ISLAS_DE_LA_BAHIA	111
readin1	PDCH	LA_PAZ				2862
readin1	PDCH	LEMPIRA				3561
readin1	PDCH	OCOTEPEQUE			705
readin1	PDCH	OLANCHO				27352
readin1	PDCH	SANTA_BARBARA		20913
readin1	PDCH	VALLE				1860
readin1	PDCH	YORO				28961
readin1	PL		ATLANTIDA			187994
readin1	PL		COLON				67836
readin1	PL		COMAYAGUA			167950
readin1	PL		COPAN				155017
readin1	PL		CORTES				1267970
readin1	PL		CHOLUTECA			297090
readin1	PL		EL_PARAISO			173538
readin1	PL		FRANCISCO_MORAZAN	1646678
readin1	PL		GRACIAS_A_DIOS		5637
readin1	PL		INTIBUCA			31799
readin1	PL		ISLAS_DE_LA_BAHIA	3000
readin1	PL		LA_PAZ				36205
readin1	PL		LEMPIRA				111954
readin1	PL		OCOTEPEQUE			21931
readin1	PL		OLANCHO				167213
readin1	PL		SANTA_BARBARA		260305
readin1	PL		VALLE				64372
readin1	PL		YORO				271506
readin1	PUD		ATLANTIDA			18713
readin1	PUD		COLON				5347
readin1	PUD		COMAYAGUA			17646
readin1	PUD		COPAN				8277
readin1	PUD		CORTES				187583
readin1	PUD		EL_PARAISO			11017
readin1	PUD		FRANCISCO_MORAZAN	370202
readin1	PUD		GRACIAS_A_DIOS		116
readin1	PUD		INTIBUCA			1567
readin1	PUD		ISLAS_DE_LA_BAHIA	1393
readin1	PUD		LA_PAZ				3940
readin1	PUD		LEMPIRA				14854
readin1	PUD		OLANCHO				17647
readin1	PUD		SANTA_BARBARA		25962
readin1	PUD		VALLE				3099
readin1	PUD		YORO				36381
readin1	MIPP	ISLAS_DE_LA_BAHIA	3545

/* There are 5 candidates with invalid names. In 4 cases, it possible to identify
* the likely name of the candidate on the ballot. However, there are some dif-
* fer​ences between the lists in the results and the ballots. In one case, these
* differences are so large, that it is not possible to identify the most likely
* name with any confidence.
* Source: http://www.tse.hn/web/documentos/papeletas/Diputados.pdf; accessed
* February 7, 2020.
replace Candidate = "SILVIA BESSY AYALA FIGUEROA" if Party_abbr == "PUD" & ///
	Department == "CORTES" & regexm(Candidate, "NO SE ENCONTRO EN PADRON ELECTORAL") == 1
replace Candidate = "LESLIE ESAU CARVAJAL CONDE" if Party_abbr == "PL" & ///
	Department == "CORTES" & regexm(Candidate, "NO SE ENCONTRO EN PADRON ELECTORAL") == 1
*replace Candidate = "" if Party_abbr == "PINU" & ///
*	Department == "CORTES" & regexm(Candidate, "NO SE ENCONTRO EN PADRON ELECTORAL") == 1
replace Candidate = "FELIPE NERY MENDEZ SORIANO" if Party_abbr == "PL" & ///
	Department == "CHOLUTECA" & regexm(Candidate, "por definir") == 1
replace Candidate = "ARNOLDO FERMIN MAYORGA ROAA" if Party_abbr == "PL" & ///
	Department == "CHOLUTECA" & regexm(Candidate, "NO SE ENCONTRO EN PADRON ELECTORAL") == 1
*/

preserve
collapse (sum) Votes, by(Party)
list Party Votes
restore
* MOVIMIENTO INDEPENDIENTE POPULAR PROGRESISTA	3545
* PARTIDO DEMOCRATA CRISTIANO DE HONDURAS		782551
* PARTIDO INNOVACION Y UNIDAD SOCIAL DEMOCRATA	1031218
* PARTIDO LIBERAL DE HONDURAS					4937995
* PARTIDO NACIONAL DE HONDURAS					8561577
* PARTIDO UNIFICACION DEMOCRATICA				723744

duplicates tag Candidate, g(dups)
list Candidate dups if dups>0
desc

charlist Candidate

* 2. Read in election results for 2013
global inpath "$path\01_Raw_data\15_Elections_Honduras\03_Elecciones_Generales_2013\"
program readin2
	import exc using "$inpath\04_Results TSE\\`1'.xlsx", locale(es_HN) ///
	first clear
	qui sum Numero
	if `r(max)' != `2' {
		di "Number of candidates `r(max)' instead of `2'"
		exit 688
	}
	qui sum Marcas if Numero == `r(max)'
	if `r(mean)' != `3' {
		di "Votes of lowest placed candidate `r(mean)' instead of `3'"
		exit 688
	}
	qui sum Marcas
	if `r(sum)' != `4' {
		di "Vote total `r(sum)' instead of `4'"
		exit 688
	}
	* 2 tests to check if the files are corrupted
	g Year = 2013
	g Department = "`1'"
	rename Partido Party
	rename Candidato Candidate
	rename Identidad ID_TSE
	rename Marcas Votes
	rename Elegido Elected
	rename Sexo Sex
	if "`1'" == "ATLANTIDA" save "$outpath\elections_hn_2013.dta", ///
	replace
	else {
	append using "$outpath\elections_hn_2013.dta"
	save "$outpath\elections_hn_2013.dta", replace
	}
end

readin2	"ATLANTIDA"   	 	72	307		815131
readin2	"CHOLUTECA"   	 	82	169		1115497
readin2	"COLON"        		36	124		299874
readin2	"COMAYAGUA"    		63	650		847407
readin2	"COPAN"        		49	179		689777
readin2	"CORTES"       		180	1130	7195963
readin2	"EL PARAISO"   		54	176		866555
readin2	"FRANCISCO MORAZAN"	208	1163	10939086
readin2	"GRACIAS A DIOS"	8	87		19956
readin2	"INTIBUCA"     		24	120		185773
readin2	"ISLAS DE LA BAHIA"	9	20		15996
readin2	"LA PAZ"      	 	27	125		158066
readin2	"LEMPIRA"      		45	34		532744
readin2	"OCOTEPEQUE"   		18	40		73812
readin2	"OLANCHO"      		63	266		1030844
readin2	"SANTA BARBARA"		72	121		1245163
readin2	"VALLE"       		32	74		219870
readin2	"YORO"       		81	263		1279319

append using "$outpath\elections_hn.dta"
save "$outpath\elections_hn.dta", replace
erase "$outpath\elections_hn_2013.dta"

* 3. Read in election results for 2017
global inpath "$path\01_Raw_data\15_Elections_Honduras\04_Elecciones_Generales_2017\"
program readin3
	import exc using "$inpath\04_Results TSE\\`1'.xlsx", locale(es_HN) ///
	first clear
	qui sum Pos
	if `r(max)' != `2' {
		di "Number of candidates `r(max)' instead of `2'"
		exit 688
	}
	qui sum Marcas if Pos == `r(max)'
	if `r(mean)' != `3' {
		di "Votes of lowest placed candidate `r(mean)' instead of `3'"
		exit 688
	}
	* 2 tests to check if the files are corrupted
	g Year = 2017
	g Department = "`1'"
	rename Partido Party
	rename Candidato Candidate
	rename Identidad ID_TSE
	rename Marcas Votes
	rename Elegido Elected
	rename Sexo Sex
	if "`1'" == "ATLANTIDA" save "$outpath\elections_hn_2017.dta", ///
	replace
	else {
	append using "$outpath\elections_hn_2017.dta"
	save "$outpath\elections_hn_2017.dta", replace
	}
end

readin3	"ATLANTIDA"  	  	80	398
readin3	"CHOLUTECA"    		90	332
readin3	"COLON"        		40	171
readin3	"COMAYAGUA"    		70	285
readin3	"COPAN"        		70	159
readin3	"CORTES"       		200	1628
readin3	"EL PARAISO"   		60	197
readin3	"FRANCISCO MORAZAN"	230	2884
readin3	"GRACIAS A DIOS"	10	56
readin3	"INTIBUCA"     		27	101
readin3	"ISLAS DE LA BAHIA"	10	18
readin3	"LA PAZ"      	 	30	73
readin3	"LEMPIRA"      		50	52
readin3	"OCOTEPEQUE"   		20	51
readin3	"OLANCHO"      		70	159
readin3	"SANTA BARBARA"		90	164
readin3	"VALLE"        		36	57
readin3	"YORO"         		90	332

append using "$outpath\elections_hn.dta"
save "$outpath\elections_hn.dta", replace
erase "$outpath\elections_hn_2017.dta"

tab Year
* 2009:   629
* 2013: 1,123
* 2017: 1,273

replace Candidate = stritrim(Candidate)
replace Candidate = strtrim(Candidate)
replace Candidate = regexr(Candidate,"Ñ","N")
replace Candidate = regexr(Candidate,"Ñ","N")
replace Candidate = "LIDESKY MILENIA ARGENAL MARTINEZ" ///
	if Candidate == "LIDESKY MILENIA ARGE?AL MARTINEZ" & Year == 2013 & ///
	Department == "EL PARAISO"
format %13.0f ID_TSE
save "$outpath\elections_hn.dta", replace
drop if regexm(Candidate, "NO SE ENCONTRO EN PADRON ELECTORAL") == 1 | ///
	Candidate == "por definir"
use "$outpath\elections_hn.dta", clear

* 4. Harmonize party information

* a) Parties in 2009
* MOVIMIENTO INDEPENDIENTE POPULAR PROGRESISTA   -                              -
* PARTIDO DEMOCRATA CRISTIANO DE HONDURAS        -	                            Centre-right
* PARTIDO INNOVACION Y UNIDAD SOCIAL DEMOCRATA   Centroizquierda                Centre-left
* PARTIDO LIBERAL DE HONDURAS                    Centroizquierda                Centre to centre-left
* PARTIDO NACIONAL DE HONDURAS                   Derecha a extrema derecha      Right-wing to far-right
* PARTIDO UNIFICACION DEMOCRATICA                Izquierda                      Left-wing

* b) Additional parties in 2013
* FRENTE AMPLIO                                  -                              Centre-left
* INDEPENDIENTE SOCIALISTA                       -                              -
* PARTIDO LIBERTAD Y REFUNDACION                 Izquierda a extrema izquierda  Left-wing to far-left
* PARTIDO ANTI CORRUPCION                        Centroderecha                  Centre-right
* PARTIDO ALIANZA PATRIOTICA                     Centroderecha a derecha        Centre-right
* UD-FAPER                                       -                              -
* UNIDOS CHOLUTECA                               -                              -

replace Party = "PARTIDO DEMOCRATA CRISTIANO DE HONDURAS" if Party == ///
	"PARTIDO DEMOCRATA CRISTIANO" & Year == 2013
replace Party = "PARTIDO INNOVACION Y UNIDAD SOCIAL DEMOCRATA" if Party == ///
	"PARTIDO INNOVACION Y UNIDAD SD" & Year == 2013
replace Party = "PARTIDO FRENTE AMPLIO" if Party == ///
	"PARTIDO FRENTE AMPLIO POLITICO ELECTORAL" & Year == 2013

* c) Additional parties in 2017
* PARTIDO VA MOVIMIENTO SOLIDARIO                -                              -

replace Party = "PARTIDO ALIANZA PATRIOTICA" if Party == ///
	"PARTIDO ALIANZA PATRIOTICA HONDUREÑA" & Year == 2017
	
* d) Left-right coding
* Coding scheme:
* 1: Left (left-wing, far-left)
* 2: Centre-left (including "centre to centre-left")
* 3: Centre and unkown party position
* 4: Centre-right
* 5: Right (right-wing, far-right)

g Position = .
replace Position = 3 if Party == "CANDIDATURA INDEPENDIENTE SOCIALISTA_CIS"
* INDEPENDIENTE SOCIALISTA                       -                              -
* No information about position found. It is a single individual.
replace Position = 3 if Party == "MOVIMIENTO INDEPENDIENTE POPULAR PROGRESISTA"
* MOVIMIENTO INDEPENDIENTE POPULAR PROGRESISTA   -                              -
* No information about position found. It is a single individual.
replace Position = 4 if Party == "PARTIDO ALIANZA PATRIOTICA"
* PARTIDO ALIANZA PATRIOTICA                     Centroderecha a derecha        Centre-right
replace Position = 4 if Party == "PARTIDO ANTICORRUPCION"
* PARTIDO ANTI CORRUPCION                        Centroderecha                  Centre-right
replace Position = 4 if Party == "PARTIDO DEMOCRATA CRISTIANO DE HONDURAS"
* PARTIDO DEMOCRATA CRISTIANO DE HONDURAS        -	                            Centre-right
replace Position = 2 if Party == "PARTIDO FRENTE AMPLIO"
* FRENTE AMPLIO                                  -                              Centre-left
replace Position = 2 if Party == "PARTIDO INNOVACION Y UNIDAD SOCIAL DEMOCRATA"
* PARTIDO INNOVACION Y UNIDAD SOCIAL DEMOCRATA   Centroizquierda                Centre-left
replace Position = 2 if Party == "PARTIDO LIBERAL DE HONDURAS"
* PARTIDO LIBERAL DE HONDURAS                    Centroizquierda                Centre to centre-left
replace Position = 1 if Party == "PARTIDO LIBERTAD Y REFUNDACION"
* PARTIDO LIBERTAD Y REFUNDACION                 Izquierda a extrema izquierda  Left-wing to far-left
replace Position = 5 if Party == "PARTIDO NACIONAL DE HONDURAS"
* PARTIDO NACIONAL DE HONDURAS                   Derecha a extrema derecha      Right-wing to far-right
replace Position = 1 if Party == "PARTIDO UNIFICACION DEMOCRATICA"
* PARTIDO UNIFICACION DEMOCRATICA                Izquierda                      Left-wing
replace Position = 3 if Party == "PARTIDO VA MOVIMIENTO SOLIDARIO"
* PARTIDO VA MOVIMIENTO SOLIDARIO                -                              -
* Self-positioning: Christian social philosophy and against extremes (i.e., cen-
* ter); founded by former members of other parties such as PARTIDO DEMOCRATA 
* CRISTIANO DE HONDURAS (centre-right). Sources:
* https://www.laprensa.hn/honduras/elecciones2017/1115117-410/partido-vamos-quiere-captar-votos-ciudadanos-radicales-jovenes	
* https://www.elheraldo.hn/pais/928954-466/cruz-asensio-sale-de-la-dc-y-funda-el-partido-vamos
replace Position = 2 if Party == "ALIANZA UD-FAPER"
* ALIANZA UD-FAPER                               -                              -
* Coded as FAPER above
replace Position = 5 if Party == "UNIDOS POR CHOLUTECA"
* UNIDOS POR CHOLUTECA                           -                              -
* This is a single individual, JOSE CESAR ORTEGA SOUZA, who ran as an indepen-
* dent but who was a former member of the National Congress of Honduras and the
* Central American Parliament for the PARTIDO NACIONAL DE HONDURAS (right-wing 
* to far-right). He was also father-in-law of the representative of the PARTIDO 
* LIBERAL DE HONDURAS (centre-left). Because of his former affiliation, coded
* as right. Source:
* https://proceso.hn/mas-noticias/32-m%C3%A1s-noticias/fallece-exdiputado-al-parlacen-del-partido-nacional-cesar-ortega-souza.html

save "$outpath\elections_hn.dta", replace

* Pedro Ernesto Ferrera Sánchez: Transgender, born male, now female; coded as 
* female


* Output for record linkage: 

local var Year Department Party Position Candidate ID_TSE Votes Elected Sex
keep `var'
order `var'
save "$outpath\elections_hn.dta", replace

preserve
tab Sex
split Candidate, p("")
g Male = Sex == "H"
g Female = Sex == "M"
collapse (sum) Male Female, by(Candidate1)
table Candidate1 if Male > 0 & Female == 0, c(mean Male)
table Candidate1 if Male == 0 & Female > 0, c(mean Female)
table Candidate1 if Male > 0 & Female > 0, c(mean Male mean Female)
restore
split Candidate, p("")

* ANGEL: ANGEL LARISSA ESPINAL PINEL - correct
* CARMEN: CARMEN ALCIDES GUEVARA CORTES - correct
* FRANCIS: FRANCIS OMAR CABRERA MIRANDA - correct
* FRANCIS: FRANCIS RAFAEL CONTRERAS RIVERA - correct
* MERLIN: MERLIN CAROLINA MARTINEZ VIGIL - correct
* NORMAN: NORMAN PORTILLO VELASQUEZ - ???
* PEDRO: PEDRO ERNESTO FERRERA SANCHEZ - born male, transgender female; coded as
* female
* TRANCITO: TRANCITO MARIBEL CHINCHILLA SORTO - correct
* TRANCITO: TRANCITO ZUNIGA VARELA - correct
* YURI: YURI FERNANDO MELARA BERLIOZ - correct
* YURI: YURI NICOLL MANZANARES TORRES - correct

* Input from record linkage 
import delimited using ///
 "$outpath\MatchesHondurasAll.csv",  ///
 delimiters(";") clear
 
replace id_tse1="" if id_tse1=="NA"
destring id_tse1, replace  

gen equal_name_department=0
replace equal_name_department=1 if candidate1==candidate2 &  ///
department1==department2
 
save "$outpath\MatchesHondurasAll.dta", replace

reshape long candidate party year department id_tse ,i(pair) j(candidate_number)

bysort pair: egen id_tse_sd=sd(id_tse)

br if id_tse_sd==0 & prob<0.9

save "$outpath\MatchesHonduras_ChecksSimon.dta", replace

* Checks Simon

use "$outpath\elections_hn_fuzzy_merge.dta", clear

* False positives
egen id_dptcand = group(Department Candidate)
bysort id_Stata: egen sd_id_dptcand = sd(id_dptcand)
sum sd_id_dptcand //No variation in names within Stata id
bysort id_Stata: egen sd_ID_TSE = sd(ID_TSE)
sum sd_ID_TSE //No variation in in ID_TSE within Stata id

* False negatives
egen num_id_stata = group(id_Stata)
bysort id_dptcand: egen sd_id_Stata = sd(num_id_stata)
sum sd_id_Stata
list Department Candidate Year Party if sd_id_Stata>0 & sd_id_Stata<. //Only 3 candidates with unknown names in Cortes in 2009
drop sd_id_Stata
bysort ID_TSE: egen sd_id_Stata = sd(num_id_stata)
sum sd_id_Stata
list Department Candidate Year Party if sd_id_Stata>0 & sd_id_Stata<. & !missing(ID_TSE)
replace id_Stata = "EL PARAISO-2013-0009" if id_Stata == "EL PARAISO-2017-0029"
replace id_Stata = "FRANCISCO MORAZAN-2009-0090" if id_Stata == "FRANCISCO MORAZAN-2017-0181"
replace id_Stata = "OLANCHO-2013-0007" if id_Stata == "OLANCHO-2017-0007"
drop sd_id_Stata num_id_stata
egen num_id_stata = group(id_Stata)
bysort ID_TSE: egen sd_id_Stata = sd(num_id_stata)
sum sd_id_Stata
list id_Stata Department Candidate Year if sd_id_Stata>0 & sd_id_Stata<. & !missing(ID_TSE)

sort Department Candidate
g ident_id_Stata = 1 if id_Stata == id_Stata[_n-1]
order id_Stata ident_id_Stata
* Manually checked if similar names have identical Stata id.

save "$outpath\elections_hn_fuzzy_merge_final.dta", replace

erase "$outpath\MatchesHonduras_ChecksSimon.dta"


/* end

 
 keep if prob>0.5
 
 keep *1 pair
 rename candidate1 Candidate 
 rename party1 Party  
 rename year1 Year 
 rename department1 Department
 
save "$outpath\rl_pair1.dta", replace


import delimited using ///
 "\\Srw-hpc5\e\Projekt Nationalräte\02_Processed_data\15_Elections_Honduras\MatchesHondurasAll.csv",  ///
 delimiters(";") clear
 
 keep *2 pair
 rename candidate2 Candidate 
 rename party2 Party  
 rename year2 Year 
 rename department2 Department
 
save "$outpath\rl_pair2.dta", replace


* Combine dataset with information from record linkage: 
use "$outpath\elections_hn.dta", clear

merge 1:1 Candidate Party Year Department using "$outpath\rl_pair1.dta"


