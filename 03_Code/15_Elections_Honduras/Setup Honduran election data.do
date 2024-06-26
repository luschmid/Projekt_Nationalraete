******************************************************
* (A) Set directories
******************************************************
version 16
clear all
set more off

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
global outpath "$path\02_Processed_data\15_Elections_Honduras"

******************************************************
* (B) READ-IN ELECTION RESULTS DATA (2009-2017)
******************************************************

* Read-in results for 2009, 2013, and 2017
* Recode parties, add information on gender 

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

*save "$outpath\elections_hn.dta", replace

* ANGEL: ANGEL LARISSA ESPINAL PINEL - correct
* CARMEN: CARMEN ALCIDES GUEVARA CORTES - correct
* FRANCIS: FRANCIS OMAR CABRERA MIRANDA - correct
* FRANCIS: FRANCIS RAFAEL CONTRERAS RIVERA - correct
* MERLIN: MERLIN CAROLINA MARTINEZ VIGIL - correct
* PEDRO: PEDRO ERNESTO FERRERA SANCHEZ - born male, transgender female; coded as
* female
* TRANCITO: TRANCITO MARIBEL CHINCHILLA SORTO - correct
* TRANCITO: TRANCITO ZUNIGA VARELA - correct
* YURI: YURI FERNANDO MELARA BERLIOZ - correct
* YURI: YURI NICOLL MANZANARES TORRES - correct


* 5. Add missing information on sex of candidates

sort Year Department Party Candidate
list Candidate Department Year if Sex == "" 

preserve
tab Sex
split Candidate, p("")
g Male = Sex == "H"
g Female = Sex == "M"
collapse (sum) Male Female, by(Candidate1)
table Candidate1 if Male > 0 & Female == 0, c(mean Male)
table Candidate1 if Male == 0 & Female > 0, c(mean Female)
table Candidate1 if Male > 0 & Female > 0, c(mean Male mean Female)
* Result: Lukas checked these cases on 11.5.2020 and checked in what follows all unclear cases
restore

split Candidate, p("")
sort Year Department Party 
list Year Department Party Candidate if Candidate1=="NATIVIDAD" | Candidate1=="ONASIS" | Candidate1=="PAIM" | Candidate1=="YAUDET" 

list Year Department Party Candidate if Candidate1=="CHRIS" | Candidate1=="DARMY" | Candidate1=="DELSY" | Candidate1=="" | Candidate1=="DENICE" | Candidate1=="DUVIS" | Candidate1=="EBLIN" | Candidate1=="ESLEM" | Candidate1=="GLIRIAN" | Candidate1=="GRACIBEL" | Candidate1=="JEIMY" | Candidate1=="JISELA" | Candidate1=="KEILY" | Candidate1=="KEYLIN" | Candidate1=="KLENY" | Candidate1=="KLEYMER" | Candidate1=="LENNY" | Candidate1=="LIDESKY" | Candidate1=="NELIN" | Candidate1=="YILDIS" | Candidate1=="" | Candidate1=="YANOLY"

* Result: Lukas checked these cases with the pictures and found that all sexes are correctly coded

drop Candidate1 Candidate2 Candidate3 Candidate4 Candidate5 Candidate6

* Sex coded based on names for observations without pictures

replace Sex = "H" if Candidate == "NAHUM OMAR QUEZADA MURILLO" & Department == "ATLANTIDA" & Year == 2009
replace Sex = "H" if Candidate == "JAVIER ANTONIO RIVERA TABORA" & Department == "COPAN" & Year == 2009
replace Sex = "H" if Candidate == "JOSE ANGEL SAAVEDRA POSADAS" & Department == "COPAN" & Year == 2009
replace Sex = "H" if Candidate == "JOHNI FRANCISCO SALINAS RAMOS" & Department == "CORTES" & Year == 2009
replace Sex = "H" if Candidate == "JOSE ADALBERTO MALDONADO ACOSTA" & Department == "CORTES" & Year == 2009
replace Sex = "H" if Candidate == "OSCAR ALEJANDRO MARTINEZ PAVON" & Department == "CORTES" & Year == 2009
replace Sex = "H" if Candidate == "VICTORINO CARRANZA PONCE" & Department == "CORTES" & Year == 2009
replace Sex = "H" if Candidate == "JORGE ALFREDO AGUILAR MARTINEZ" & Department == "CORTES" & Year == 2009
replace Sex = "H" if Candidate == "JOSE NOEL PAZ FERNANDEZ" & Department == "CORTES" & Year == 2009
replace Sex = "H" if Candidate == "ANTONIO AGUILAR CHAVEZ" & Department == "CORTES" & Year == 2009
replace Sex = "H" if Candidate == "GERMAN GENARO GUTIERREZ RAPALO" & Department == "CORTES" & Year == 2009
replace Sex = "H" if Candidate == "HECTOR MANUEL CALDERON ORTIZ" & Department == "CORTES" & Year == 2009
replace Sex = "H" if Candidate == "LUCIO ERNESTO LOPEZ CHAVEZ" & Department == "CORTES" & Year == 2009
replace Sex = "H" if Candidate == "FREDY RAMON OLIVA CALDERON" & Department == "EL PARAISO" & Year == 2009
replace Sex = "H" if Candidate == "FRANKLIN ALEXIS MURILLO RODRIGUEZ" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex = "H" if Candidate == "JOSE GUILLERMO CANO VALLADARES" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex = "H" if Candidate == "MARIO ALFREDO MOYA LOVO" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex = "H" if Candidate == "MARIO AREVALO COLINDRES" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex = "H" if Candidate == "REYNALDO ANTONIO TILGUATH" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex = "H" if Candidate == "CELESTINO GONZALEZ VASQUEZ" & Department == "INTIBUCA" & Year == 2009
replace Sex = "H" if Candidate == "JOSUE NOE RUIZ PADILLA" & Department == "OLANCHO" & Year == 2009
replace Sex = "H" if Candidate == "JOSE LEONEL CACERES" & Department == "YORO" & Year == 2009
replace Sex = "H" if Candidate == "LUIS AUDON ARTEAGA MARADIAGA" & Department == "YORO" & Year == 2009
replace Sex = "H" if Candidate == "ADRIAN CRUZ REYES" & Department == "YORO" & Year == 2009
replace Sex = "H" if Candidate == "HECTOR DARIO MARTINEZ SOLIZ" & Department == "YORO" & Year == 2009
replace Sex = "H" if Candidate == "JULIO CESAR BUESO LUNA" & Department == "COLON" & Year == 2017
replace Sex = "H" if Candidate == "RICARDO ALBERTO SUAZO MEDINA" & Department == "CORTES" & Year == 2017
replace Sex = "H" if Candidate == "ALEX WILFREDO MEJIA HERNANDEZ" & Department == "FRANCISCO MORAZAN" & Year == 2017
replace Sex = "H" if Candidate == "CELESTINO GONZALEZ VASQUEZ" & Department == "INTIBUCA" & Year == 2017
replace Sex = "H" if Candidate == "OLVIN REYNIERY GONZALES DOMINGUEZ" & Department == "INTIBUCA" & Year == 2017
replace Sex = "H" if Candidate == "CARLOS ROBERTO PADILLA ANTUNEZ" & Department == "OLANCHO" & Year == 2017
replace Sex = "H" if Candidate == "JUSTO RUFINO DUARTE LOPEZ" & Department == "OLANCHO" & Year == 2017
replace Sex = "H" if Candidate == "LUIS ALONSO BU SABILLON" & Department == "SANTA BARBARA" & Year == 2017
replace Sex = "H" if Candidate == "DENIS ADALID PADILLA CRUZ" & Department == "YORO" & Year == 2017

replace Sex = "M" if Candidate == "MARIA ENMA MENDEZ BAJURTO" & Department == "COPAN" & Year == 2009
replace Sex = "M" if Candidate == "MARIA LEONOR GONZALES" & Department == "COPAN" & Year == 2009
replace Sex = "M" if Candidate == "XIOMARA LIZETH PAVON" & Department == "CORTES" & Year == 2009
replace Sex = "M" if Candidate == "WALESKA ANTONIA FERNANDEZ RAMIREZ" & Department == "CORTES" & Year == 2009
replace Sex = "M" if Candidate == "BLANCA AIDE RODAS CARDONA" & Department == "CORTES" & Year == 2009
replace Sex = "M" if Candidate == "YESENIA WALESKA ANDINO LOPEZ" & Department == "CORTES" & Year == 2009
replace Sex = "M" if Candidate == "ISABEL SANDOVAL" & Department == "EL PARAISO" & Year == 2009
replace Sex = "M" if Candidate == "PATRICIA ICELINA VALLECILLO ARDON" & Department == "EL PARAISO" & Year == 2009
replace Sex = "M" if Candidate == "CONSUELO DACOSTA MEMBRENO" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex = "M" if Candidate == "WENDY VANESA FLORES SAUCEDA" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex = "M" if Candidate == "ALEJANDRA MELISSA LAINEZ OBANDO" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex = "M" if Candidate == "DIANA PATRICIA VASQUEZ CARDENAS" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex = "M" if Candidate == "REINA ESTHER OCHOA ORTIZ" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex = "M" if Candidate == "MARTA GARCIA GONZALEZ" & Department == "INTIBUCA" & Year == 2009
replace Sex = "M" if Candidate == "REYNA AUXILIADORA MENDEZ GONZALEZ" & Department == "INTIBUCA" & Year == 2009
replace Sex = "M" if Candidate == "DORIS AMANDA MARTINEZ" & Department == "LEMPIRA" & Year == 2009
replace Sex = "M" if Candidate == "FRANCISCA GUADALUPE BREVE MEJIA" & Department == "OLANCHO" & Year == 2009
replace Sex = "M" if Candidate == "ROSA PATRICIA CRUZ MEJIA" & Department == "OLANCHO" & Year == 2009
replace Sex = "M" if Candidate == "SEIDY LILI SANTOS REYES" & Department == "OLANCHO" & Year == 2009
replace Sex = "M" if Candidate == "CLEOFAS ENAMORADO CABALLERO" & Department == "SANTA BARBARA" & Year == 2009
replace Sex = "M" if Candidate == "CLAUDIA SUYAPA MARTINEZ ENAMORADO" & Department == "YORO" & Year == 2009
replace Sex = "M" if Candidate == "MARIA SUYAPA GONZALES VELIZ" & Department == "YORO" & Year == 2009
replace Sex = "M" if Candidate == "KAREN LIZZETH MUNGUIA" & Department == "COLON" & Year == 2013
replace Sex = "M" if Candidate == "NINA ROSARIO RODRIGUEZ RAMOS" & Department == "COMAYAGUA" & Year == 2013
replace Sex = "M" if Candidate == "FLORA IDALMA LOPEZ VALLECILLO" & Department == "SANTA BARBARA" & Year == 2013
replace Sex = "M" if Candidate == "BEATRIZ ANGELICA CANALES ORDONEZ" & Department == "FRANCISCO MORAZAN" & Year == 2017
replace Sex = "M" if Candidate == "BERTHA ROSA AGUILAR AGUILAR" & Department == "INTIBUCA" & Year == 2017
replace Sex = "M" if Candidate == "MA SUSANA MARTINEZ" & Department == "LA PAZ" & Year == 2017
replace Sex = "M" if Candidate == "GUADALUPE ELIZABETH ZELAYA LOBO" & Department == "OLANCHO" & Year == 2017
replace Sex = "M" if Candidate == "YESSICA SULEMA SIERRA LAINEZ" & Department == "VALLE" & Year == 2017

* first names not on the list above

list Candidate Department Year if Sex == ""
replace Sex ="M" if Candidate == "SERGIA VASQUEZ SABILLON" & Department == "CORTES" & Year == 2009
replace Sex ="H" if Candidate == "LENING NAHUN URRUTIA ALVARADO" & Department == "CORTES" & Year == 2009 
* https://www.facebook.com/Lening-Nahun-Urrutia-Alvarado-1568680239857344/
replace Sex ="M" if Candidate == "ONIS WALDINA VELASQUEZ DOLMO" & Department == "CORTES" & Year == 2009
* https://es-es.facebook.com/oniswaldina.velasquezdolmo
replace Sex ="M" if Candidate == "GERALDINA BOGRAN RIVERA" & Department == "FRANCISCO MORAZAN" & Year == 2009
replace Sex ="M" if Candidate == "MARLENY FUNES VASQUEZ" & Department == "LA PAZ" & Year == 2009
replace Sex ="H" if Candidate == "KENEDY SANTOS GODFRY" & Department == "GRACIAS A DIOS" & Year == 2013
replace Sex ="M" if Candidate == "RUMUALDA YESENIA MOREL PAREDES" & Department == "SANTA BARBARA" & Year == 2017

* Output for record linkage: 

local var Year Department Party Position Candidate ID_TSE Votes Elected Sex
keep `var'
order `var'
save "$outpath\elections_hn_fuzzy_merge_input.dta", replace
erase "$outpath\elections_hn.dta"



******************************************************
* (C) READ-IN OUTPUT FROM RECORD LINKAGE
******************************************************

use "$outpath\elections_hn_fuzzy_merge_ouptut.dta", clear

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

keep id_Stata Candidate Year Party Votes Department

save "$outpath\elections_hn_fuzzy_merge_final.dta", replace


********************************************************************************
* (D) MERGE OUTPUT FROM RECORD LINKAGE AND RUNNING VARIABLE TO INITIAL DATASET
********************************************************************************

* (i) Save output from calculation of running variable in stata format	
* Note: The rv calculation is done in RV Honduras.R

preserve
import delimited using "$outpath\Margins_Honduras.csv", delim(";") clear
drop if candidate=="NA" // to do: these are 7 candidates with the same number of seats that cause one additional observation in the respective party
rename candidate Candidate
rename year Year 
rename votes_h Votes
rename district Department
rename party Party
rename votemargin Votemargin
keep Candidate Year Party Votes Department Votemargin
destring Votes Year, replace
save "$outpath\Margins_Honduras.dta", replace
restore

use "$outpath\elections_hn_fuzzy_merge_input.dta", clear
tostring   Candidate  Party, replace
merge 1:1 Candidate Year Party Votes Department using "$outpath\elections_hn_fuzzy_merge_final.dta", gen(merge_fuzzy_merge)
merge 1:1 Candidate Year Party Votes Department using "$outpath\Margins_Honduras.dta", gen(merge_rv)

order id_Stata Candidate Department Year Party Sex Elected Votes Votemargin

* (ii) Generate running variable (rescaled vote margin)

preserve
import delimited using "$path\01_Raw_data\15_Elections_Honduras\honduras_seats.csv", delim(";") clear
rename district Department
rename year Year
save "$outpath\Department_Aggregates.dta", replace
restore

merge m:1 Department Year using  "$outpath\Department_Aggregates.dta", gen(merge_votes_department)

gen votemargin_rel=Votemargin/total_electoral 

* (iii) Generate time variables

egen ID_time=group(Year)
egen ID_pers=group(id_Stata)

xtset ID_pers ID_time

bysort ID_pers: egen ID_time_min=min(ID_time) 

sum Year
local maxyear=r(max)
local minyear=r(min)

gen participation_L1=0
replace participation_L1=1 if L1.Elected!=. & Year!=`minyear'
replace participation_L1=. if Year==`minyear' // set to missing for first year of dataset

gen participation_F1=0
replace participation_F1=1 if F1.Elected!=. & Year!=`maxyear'
replace participation_F1=. if  Year==`maxyear' // set to missing for last year of dataset

sort ID_pers ID_time

gen Votemargin_L1=.
replace Votemargin_L1=L1.Votemargin if Year!=`minyear'
gen votaron_L1=.
replace votaron_L1=L1.votaron if Year!=`minyear'
gen total_electoral_L1=.
replace total_electoral_L1=L1.total_electoral if Year!=`minyear'
gen no_seats_L1=.
replace no_seats_L1=L1.no_seats if Year!=`minyear'

 
 

local vars "Elected" 
foreach var in `vars'{			
gen `var'_F1=F1.`var'
gen `var'_L1=L1.`var'
replace `var'_F1=0 if `var'_F1==. & Year!=`maxyear'
replace `var'_L1=0 if `var'_L1==. & Year!=`minyear'
}

/* (iv) Merge rv by Kotakorpi et al. (2017) and generate rescaled rv

preserve
import delimited using "$path\02_Processed_data\13_Running_variable\kotakorpietal2017_honduras.csv", delim(",") clear
duplicates drop id_stata year, force
rename id_stata id_Stata
rename year Year
save  "$outpath\kotakorpietal2017_honduras.dta", replace
restore

merge 1:1 id_Stata Year using "$outpath\kotakorpietal2017_honduras.dta", gen(merge_kotakorpi2017)
erase "$outpath\kotakorpietal2017_honduras.dta"
drop party_num


sort Department Year Party
egen party_num=group(Department Year Party) // numeric version of party

gsort party_num -Votes -Elected Candidate
qui bysort party_num: egen Votes_h_rank=rank(-Votes), unique // generate rank of candidate on party list
gsort party_num -p -Votes -Elected Candidate
qui bysort party_num: egen p_rank=rank(-p), unique // generate rank of candidate on party list

gen rank_diff=Votes_h_rank-p_rank
sum rank_diff
sum rank_diff if rank_diff!=0

qui bysort party_num: egen rank_diff_max=max(rank_diff)
sort party_num Votes_h
br ID canton year name firstname Elected Votes Votes_h_rank p p_rank if rank_diff_max!=0
br ID canton year name firstname Elected Votes Votes_h_rank p p_rank if rank_diff!=0

qui bysort party_num: egen p_min_help=min(p) if Elected==1
qui bysort party_num Elected: egen p_max_help=max(p) if Elected==0
qui bysort party_num: egen p_min=min(p_min_help) 
qui bysort party_num: egen p_max=min(p_max_help)
qui bysort party_num: gen number_Elected_help=_N if Elected==1
qui bysort party_num: gen number_notElected_help=_N if Elected==0
qui bysort party_num: egen number_Elected=min(number_Elected_help) 
qui bysort party_num: egen number_notElected=min(number_notElected_help)
drop *_help
replace number_Elected=0 if number_Elected==.
replace number_notElected=0 if number_notElected==.

gen p_marginal=(p_min+p_max)/2
replace p_marginal=1 if number_Elected==0
replace p_marginal=0 if number_notElected==0

gen p_transformed=p-p_marginal
hist p_transformed
bysort Elected: sum p_transformed

br ID party_num canton year name firstname Elected Votes p p_marginal number_* if p_transformed<0.2 & Elected==1
br ID party_num canton year name firstname Elected Votes p p_marginal number_* if p_transformed>0.2 & Elected==0
br ID party_num canton year name firstname Elected Votes p p_marginal number_* if party_num==1989

drop merge*
*/

save "$outpath\elections_hn_final.dta", replace

//erase "$outpath\elections_hn_fuzzy_merge_final.dta"
erase "$outpath\elections_hn_fuzzy_merge_input.dta"
erase "$outpath\Margins_Honduras.dta"
erase  "$outpath\Department_Aggregates.dta"
