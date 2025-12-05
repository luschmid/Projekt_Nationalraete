clear 	
cap log close
set more off	  

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"

*** Enrich NR Candidate data with Geo-Coordinates

* Note: This do-file was written by Mark. The steps (roman letters), the
* 		part "Create ID for geo coordinates", and "Check whether new dataset ..." 
*		was added by Lukas (July 2021).

* (i) Preparation of geo data: drop duplicate gdner year entries

use "$path\02_Processed_data\08_Municipalities\GdeNr2018_Geo_1931-2018.dta", clear
duplicates list gdenr year
duplicates drop gdenr year, force
sort gdenr year
save "$path\02_Processed_data\08_Municipalities\GdeNr2018_Geo_1931-2018_tmp.dta", replace

* (ii) Merge NR dataset and geo dataset on residence municipality

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
gen gdenr = municipalityno
*drop _merge
sort gdenr year
merge m:1 gdenr year using "$path\02_Processed_data\08_Municipalities\GdeNr2018_Geo_1931-2018_tmp.dta"
/*
    Result                           # of obs.
    -----------------------------------------
    not matched                       246,400
        from master                       279  (_merge==1)
        from using                    246,121  (_merge==2)

    matched                            41,346  (_merge==3)
    -----------------------------------------

	--> _merge==1: unmerged NR: only candidates living in foreign countries (BFS No = 9950)
*/	
drop if _merge==2
drop _merge
drop gdenr

drop ctnnr gdename_sht district_histnr district entry_yr entry_mutnr entry_type exit_mutnr exit_type ///
exit_yr modif_yr histnr_2018 ctnnr_2018 GMDNAME BZNR KTNR GRNR AREA_HA X_MIN X_MAX Y_MIN Y_MAX ///
X_CNTR Y_CNTR Z_MIN Z_MAX Z_AVG Z_MED Z_CNTR E_MIN E_MAX N_MIN N_MAX

rename E_CNTR E_CNTR_w
rename N_CNTR N_CNTR_w
rename gdenr_2018 gdenr_2018_w
label var E_CNTR_w "E-coordinate centroid Wohnort (BFS-GdeNr)"
label var N_CNTR_w "N-coordinate centroid Wohnort (BFS-GdeNr)"
label var gdenr_2018_w "Municipal number 2018 (Wohnort)"

* (iii) Merge NR dataset and geo dataset on origin municipality

* (a) corrections on BFS Numbers in Bürgerort 1 & 2

replace originno1 = 4941 if originno1 == 4942 & year == 1975 // Ottoberg 				4942 	1975	-- > Märstetten 		4941
replace originno1 = 5422 if originno1 == 5433 & year == 2011 // Pizy 					5433	2011	-- > Aubonne 			5422
replace originno1 = 5464 if originno1 == 5453 & year == 2011 // Chabrey 				5453	2011	-- > Vully-les-Lacs		5464
replace originno1 = 5451 if originno1 == 5461 & year == 2011 // Oleyres					5461	2011	-- > Avanches			5451
replace originno1 = 5464 if originno1 == 5463 & year == 2011 // Villars-le-Grand		5463	2011	-- > Vully-les-Lacs		5464
replace originno1 = 5540 if originno1 == 5517 & year == 2011 // Dommartin				5517	2011	-- > Moutilliez			5540
replace originno1 = 5540 if originno1 == 5528 & year == 2011 // Naz						5528	2011	-- > Moutilliez			5540
replace originno1 = 5804 if originno1 == 5538 & year == 2011 // Villars-Tiercelin		5538	2011	-- > Jorat-Menthue		5804
replace originno1 = 5613 if originno1 == 5603 & year == 2011 // Epesses					5603	2011	-- > Bourg-en-Lavaux	5613
replace originno1 = 5613 if originno1 == 5608 & year == 2011 // Riex					5608	2011	-- > Bourg-en-Lavaux	5613
replace originno1 = 5804 if originno1 == 5796 & year == 2011 // Peney-le-Jorat			5796	2011	-- > Jorat-Menthue		5804
replace originno1 = 5831 if originno1 == 5818 & year == 2011 // Granges-près-Marnand	5818	2011	-- > Valbroye			5813
replace originno1 = 6181 if originno1 == 6180 & year == 2003 // Ried-Mörel				6180	2003	-- > Riederalp			6181

replace originno2 = 5613 if originno2 == 5603 & year == 2011 // Epesses					5603	2011	-- > Bourg-en-Lavaux	5613

* (b) merge datasets on bürgerort
forv v =1(1)6 {
	gen gdenr = originno`v'
	sort gdenr year

	merge m:1 gdenr year using "$path\02_Processed_data\08_Municipalities\GdeNr2018_Geo_1931-2018_tmp.dta"
	tab gdenr originno`v' if gdenr!=. & _merge==1
	drop if _merge==2
	drop _merge 
	drop gdenr

	drop ctnnr gdename_sht district_histnr district entry_yr entry_mutnr entry_type exit_mutnr exit_type ///
	exit_yr modif_yr histnr_2018 ctnnr_2018 gdenr_2018 GMDNAME BZNR KTNR GRNR AREA_HA X_MIN X_MAX Y_MIN Y_MAX ///
	X_CNTR Y_CNTR Z_MIN Z_MAX Z_AVG Z_MED Z_CNTR E_MIN E_MAX N_MIN N_MAX

	rename E_CNTR E_CNTR_b`v'
	rename N_CNTR N_CNTR_b`v'
	label var E_CNTR_b`v' "E-coordinate centroid Bürgerort`v' (BFS-GdeNr)"
	label var N_CNTR_b`v' "N-coordinate centroid Bprgerort`v' (BFS-GdeNr)"
}


* (iv) add language information
preserve
use "$path\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018.dta", clear
duplicates drop gdenr_2018, force
sort gdenr_2018
save "$path\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018_temp.dta", replace
restore

rename gdenr_2018_w gdenr_2018
sort gdenr_2018 year
merge m:1 gdenr_2018 using "$path\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018_temp.dta"
drop if _merge == 2
/*
    Result                           # of obs.
    -----------------------------------------
    not matched                           515
        from master                       209  (_merge==1)
        from using                        306  (_merge==2)

    matched                            41,254  (_merge==3)
	
-- > master only (_merge == 1): all foreign municipalities 
*/
drop _merge gem_cd gem_name
rename gdenr_2018 gdenr_2018_w
rename language language_w
label var language_w "majority language in Wohnort"

* (v) Create ID for geo coordinates
* Note: This is necessary b/c both NR and Sugarcube data is not unique
* 		in terms of ID and to merge it with ground truth data. 

gen e_id_polit = int(E_CNTR_w) 
gen n_id_polit = int(N_CNTR_w)
label var e_id_polit "Geo-ID East for link w/ Ground Truth"
label var n_id_polit "Geo-ID Nort for link w/ Ground Truth" 

* (vi) Sort and save dataset

sort ID year

save "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta", replace	
erase "$path\02_Processed_data\08_Municipalities\GdeNr2018_Geo_1931-2018_tmp.dta"
erase "$path\02_Processed_data\08_Municipalities\BfS_Municipalities_Language_gdenr2018_temp.dta"

/* (vii) Check whether new dataset merges uniquely with ground truth
* Note: Lukas, July 2021.

* (a) Prepare NR-Geo data 

use "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo.dta", clear	

rename ID id_polit
duplicates drop id_polit e_id_polit n_id_polit, force

keep id_polit e_id_polit n_id_polit

save "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo_Nodups.dta", replace	


* (b) Prepare ground truth data and merge it to NR-Geo data data

use "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\Sugarcube_Ground_Truth_ID_Sandro.dta", clear
duplicates tag id_polit e_id_polit n_id_polit, gen(dups)
tab dups
duplicates tag id_polit e_id_polit n_id_polit id_sug year_sug e_id_sug n_id_sug, gen(dups2)
tab dups2
egen group_ids=group(id_polit e_id_polit n_id_polit id_sug year_sug e_id_sug n_id_sug)
br if dups2>0

merge m:1 id_polit e_id_polit n_id_polit using "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo_Nodups.dta"

* Result: All ground truth entries have an entry in politician data (no entry with _merge==1)

erase "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo_Nodups.dta"
