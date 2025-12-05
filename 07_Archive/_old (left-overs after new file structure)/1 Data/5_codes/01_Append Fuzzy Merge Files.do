
	/*** DEFINE MAIN PATH ***/
	capture cd "C:\Dropbox\Incumbency Advantage"
	if _rc==0{
		global hauptpfad "C:\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "C:\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "C:\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "C:\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "C:\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"

	}
	
    capture cd "D:\SchmidLu\Dropbox\Projekt Nationalräte"

	if _rc==0{ 
	global hauptpfad "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"
	}

	
	****************************************
	* (A) VASCO
	****************************************
		
	* (i) Start with AG
		
	set more off			
	cd "$datapath\5 fuzzy merge back in\vasco"
	import excel  "erf_vasco_AG.xlsx"	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	saveold "$datapath\5 fuzzy merge back in\data_federal_temp.dta", replace	version(12)


	* (ii) Merge all other files in folder endgültig fertig

	local dateien AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH
	foreach datei in `dateien'{
	cd "$datapath\5 fuzzy merge back in\vasco"
	import excel using erf_vasco_`datei'	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	cd "$datapath\5 fuzzy merge back in"
	append using data_federal_temp.dta
	saveold "$datapath\5 fuzzy merge back in\data_federal_temp.dta", replace	version(12)
	}
	
	gen Kanton=substr(id_Stata,1,2)
	
	local vars indicator candidateid_threshold_min candidateid_threshold_opt candidateid_threshold_max candidateid_     id_Stata  unsicher Kommentar   
	foreach var in `vars'{
	rename  `var' `var'_vasco
	}
	
	duplicates drop Vorname Nachname Geschlecht Geburtsjahr Listenbezeichnung Gemeinde Kanton Wahljahr, force // drop AG observations
	
	* (iii) Change birthyears and surnames for some observations (probably this information was correct by BfS after we gave the fuzzy merge data to Vasco Schelbert)
	
	replace Geburtsjahr=1942 if Vorname=="Brigitte" & Nachname=="Schmid" & Wahljahr==1983 & Kanton=="ZH" 
	replace Geburtsjahr=1942 if Vorname=="Brigitte" & Nachname=="Schmid" & Wahljahr==1987 & Kanton=="ZH" 
	replace Geburtsjahr=1951 if Vorname=="Christoph" & Nachname=="Keller" & Wahljahr==1979 & Kanton=="ZH" 
	replace Geburtsjahr=1909 if Vorname=="Christoph" & Nachname=="Wolfensberger" & Wahljahr==1975 & Kanton=="ZH" 
	replace Geburtsjahr=1960 if Vorname=="Daniel" & Nachname=="Holzreuter" & Wahljahr==1995 & Kanton=="ZH" 
	replace Geburtsjahr=1945 if Vorname=="Françoise" & Nachname=="Glatz-Studer" & Wahljahr==2015 & Kanton=="VD" 
	replace Geburtsjahr=1956 if Vorname=="Guy" & Nachname=="Morin" & Wahljahr==1991 & Kanton=="BS" 
	replace Geburtsjahr=1943 if Vorname=="Kurt" & Nachname=="Schreiber" & Wahljahr==1991 & Kanton=="ZH" 
	replace Geburtsjahr=1946 if Vorname=="Marianne" &  Nachname=="Jaccard" & Wahljahr==1987 & Kanton=="VD" 
	replace Geburtsjahr=1946 if Vorname=="Marianne" &  Nachname=="Jaccard" & Wahljahr==1991 & Kanton=="VD" 
	replace Geburtsjahr=1944 if Vorname=="Niklaus" &  Nachname=="Scherr" & Wahljahr==1979 & Kanton=="ZH" 
	replace Geburtsjahr=1941 if Vorname=="Rudolf" &  Nachname=="Bautz" & Wahljahr==1983 & Kanton=="ZH" 


	replace Vorname="Margrit" if Vorname=="Margrith" & Nachname=="Kolp" & Wahljahr==1999 & Kanton=="ZH" 
	replace Nachname="Wäfler"  if Vorname=="Markus" & Nachname=="Wäffler"  & Wahljahr==1983 & Kanton=="ZH" 
	replace Nachname="Sturny"  if Vorname=="Max" & Nachname=="Sturni"  & Wahljahr==1979 & Kanton=="ZH" 
	replace Vorname="Regine"  if Vorname=="Regina" & Nachname=="Aeppli"  & Wahljahr==1995 & Kanton=="ZH" 
	replace Vorname="Stefan"  if Vorname=="Stephan" & Nachname=="Dollenmeier"  & Wahljahr==2011 & Kanton=="ZH" 
	
	
	*(iv) Generate final IDs
	
	gen id_Final_vasco=id_Stata
	replace id_Final_vasco=id_recoded if id_recoded!="." & id_recoded!=""

	saveold erf_vasco_all.dta, version(12) replace
	
	
	****************************************
	* (B) JONAS
	****************************************
		
	* (i) Start with AG
		
	set more off			
	cd "$datapath\5 fuzzy merge back in\jonas"
	import excel  "erf_2015_AG.xlsx"	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	saveold "$datapath\5 fuzzy merge back in\data_federal_temp.dta", replace	version(12)


	* (ii) Merge all other files 

	local dateien AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH
	foreach datei in `dateien'{
	cd "$datapath\5 fuzzy merge back in\jonas"
	import excel using erf_2015_`datei'	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	cd "$datapath\5 fuzzy merge back in"
	append using data_federal_temp.dta
	saveold "$datapath\5 fuzzy merge back in\data_federal_temp.dta", replace	version(12)
	}
	
	gen Kanton=substr(id_Stata,1,2)

	rename id_recoded_2011 id_recoded 
	replace Kommentar=Kommentar_2011 if Kommentar==""
	drop *_2011
	
	local vars indicator candidateid_threshold_min candidateid_threshold_opt candidateid_threshold_max candidateid_     id_Stata unsicher Kommentar   
	foreach var in `vars'{
	rename  `var' `var'_jonas
	}
	
	
	
	duplicates drop Vorname Nachname Geschlecht Geburtsjahr Listenbezeichnung Gemeinde Kanton Wahljahr, force // drop AG observations
	
	
	* (iii) Change birthyears and surnames for some observations (probably this information was correct by BfS after we gave the fuzzy merge data to Jonas Röllin)
	
	replace Geburtsjahr=1942 if Vorname=="Brigitte" & Nachname=="Schmid" & Wahljahr==1983 & Kanton=="ZH" 
	replace Geburtsjahr=1942 if Vorname=="Brigitte" & Nachname=="Schmid" & Wahljahr==1987 & Kanton=="ZH" 
	replace Geburtsjahr=1951 if Vorname=="Christoph" & Nachname=="Keller" & Wahljahr==1979 & Kanton=="ZH" 
	replace Geburtsjahr=1909 if Vorname=="Christoph" & Nachname=="Wolfensberger" & Wahljahr==1975 & Kanton=="ZH" 
	replace Geburtsjahr=1960 if Vorname=="Daniel" & Nachname=="Holzreuter" & Wahljahr==1995 & Kanton=="ZH" 
	replace Geburtsjahr=1945 if Vorname=="Françoise" & Nachname=="Glatz-Studer" & Wahljahr==2015 & Kanton=="VD" 
	replace Geburtsjahr=1956 if Vorname=="Guy" & Nachname=="Morin" & Wahljahr==1991 & Kanton=="BS" 
	replace Geburtsjahr=1943 if Vorname=="Kurt" & Nachname=="Schreiber" & Wahljahr==1991 & Kanton=="ZH" 
	replace Geburtsjahr=1946 if Vorname=="Marianne" &  Nachname=="Jaccard" & Wahljahr==1987 & Kanton=="VD" 
	replace Geburtsjahr=1946 if Vorname=="Marianne" &  Nachname=="Jaccard" & Wahljahr==1991 & Kanton=="VD" 
	replace Geburtsjahr=1944 if Vorname=="Niklaus" &  Nachname=="Scherr" & Wahljahr==1979 & Kanton=="ZH" 
	replace Geburtsjahr=1941 if Vorname=="Rudolf" &  Nachname=="Bautz" & Wahljahr==1983 & Kanton=="ZH" 


	replace Vorname="Margrit" if Vorname=="Margrith" & Nachname=="Kolp" & Wahljahr==1999 & Kanton=="ZH" 
	replace Nachname="Wäfler"  if Vorname=="Markus" & Nachname=="Wäffler"  & Wahljahr==1983 & Kanton=="ZH" 
	replace Nachname="Sturny"  if Vorname=="Max" & Nachname=="Sturni"  & Wahljahr==1979 & Kanton=="ZH" 
	replace Vorname="Regine"  if Vorname=="Regina" & Nachname=="Aeppli"  & Wahljahr==1995 & Kanton=="ZH" 
	replace Vorname="Stefan"  if Vorname=="Stephan" & Nachname=="Dollenmeier"  & Wahljahr==2011 & Kanton=="ZH" 
	
	
	*(iv) Generate final IDs
	
	gen id_Final_jonas=id_Stata
	replace id_Final_jonas=id_recoded if id_recoded!="." & id_recoded!=""

	saveold erf_jonas_all.dta, version(12) replace
	
	
	
	
	
	*********************************************
	* (C) LUZERN (JANA JARCK UND ROXANE BRUENDLER)
	*********************************************
		
	* (i) Start with AG
		
	set more off			
	cd "$datapath\5 fuzzy merge back in\luzern"
	import excel  "erf_AG_unilu_RB.xlsx"	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	saveold "$datapath\5 fuzzy merge back in\data_federal_temp.dta", replace	version(12)


	* (ii) Merge all other files

	local dateien AI_unilu_RB AR_unilu_RB BEJU_unilu_RB BL_unilu_JJ BS_unilu_JJ FR_unilu_RB_JJ GE_unilu_JJ GL_unilu_RB GR_unilu_RB LU_unilu_RB NE_unilu_RB NW_unilu_RB OW_unilu_RB SG_unilu_RB  SH_unilu_JJ  SO_unilu_JJ SZ_unilu_JJ TG_unilu_JJ TI_unilu_JJ_RB UR_unilu_JJ VD_unilu_JJ VS_unilu_JJ ZG_JJ ZH_JJ
	foreach datei in `dateien'{
	cd "$datapath\5 fuzzy merge back in\luzern"
	display "`datei'"
	import excel using erf_`datei'	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	cd "$datapath\5 fuzzy merge back in"
	append using data_federal_temp.dta
	saveold "$datapath\5 fuzzy merge back in\data_federal_temp.dta", replace	version(13)
	}
	
	gen Kanton=substr(id_Stata,1,2)
	
	local vars indicator candidateid_threshold_min candidateid_threshold_opt candidateid_threshold_max candidateid_     id_Stata unsicher Kommentar   
	foreach var in `vars'{
	rename  `var' `var'_lu
	}
	
	duplicates drop Vorname Nachname Geschlecht Geburtsjahr Listenbezeichnung Gemeinde Kanton Wahljahr, force // drop AG observations

	recast str244 Beruf, force 
	
	*(iii) Generate final IDs
	
	gen id_Final_lu=id_Stata_lu
	replace id_Final_lu=id_recoded if  id_recoded!="." & id_recoded!=""
	
	saveold erf_luzern_all.dta,  replace version(12)

	

	**************************************************
	* (D) FRIBOURG (NICOLAS MAURI AND BENEDICTE DROZ)
	**************************************************
		
	* (i) Start with AG
		
	set more off			
	cd "$datapath\5 fuzzy merge back in\fribourg"
	import excel  "erf_AG.xlsx"	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	saveold "$datapath\5 fuzzy merge back in\data_federal_temp.dta", replace	version(12)


	* (ii) Merge all other files

	local dateien AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH
	foreach datei in `dateien'{
	cd "$datapath\5 fuzzy merge back in\fribourg"
	display "`datei'"
	import excel using erf_`datei'	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	cd "$datapath\5 fuzzy merge back in"
	append using data_federal_temp.dta
	saveold "$datapath\5 fuzzy merge back in\data_federal_temp.dta", replace	version(13)
	}
	
	gen Kanton=substr(id_Stata,1,2)
	
	local vars indicator candidateid_threshold_min candidateid_threshold_opt candidateid_threshold_max candidateid_     id_Stata  unsicher Kommentar   
	foreach var in `vars'{
	rename  `var' `var'_fr
	}
	
	duplicates drop Vorname Nachname Geschlecht Geburtsjahr Listenbezeichnung Gemeinde Kanton Wahljahr, force // drop AG observations

	recast str244 Beruf, force 
	
	// reverse all name changes made by coders in Fribourg
	
	replace Nachname="Stebler" if id_Stata=="BS-1955-0046" & Nachname=="Stehler"
	replace Nachname="Degen" if id_Stata=="BS-1967-0016" & Nachname=="Degea"
	replace Nachname="Corbat" if id_Stata=="GE-1963-0006" & Nachname=="Corbaz"
	replace Nachname="Hafter" if id_Stata=="TG-1967-0013" & Nachname=="Haffter"
	replace Nachname="Hugli" if id_Stata=="VD-1963-0039" & Nachname=="Huggli"
	replace Nachname="Günthard" if id_Stata=="ZH-1931-0053" & Nachname=="Günthart"
	replace Nachname="Bosshart" if id_Stata=="ZH-1935-0025" & Nachname=="Bosshardt"
	replace Nachname="Bänninger" if id_Stata=="ZH-1935-0034" & Nachname=="Bänniger"
	replace Nachname="König" if id_Stata=="ZH-1955-0152" & Nachname=="Kong"
	replace Nachname="Sturny" if id_Stata=="ZH-1971-0406" & Nachname=="Sturni"
	replace Nachname="Wüfler" if id_Stata=="ZH-1983-0553" & Nachname=="Wüffler"
	
	
	
	*(iii) Generate final IDs
	
	gen id_Final_fr=id_Stata_fr
	replace id_Final_fr=id_recoded if id_recoded!="." & id_recoded!=""
	
	
	rename Vorname  Vorname_fr
	rename Geburtsjahr Geburtsjahr_fr
	gen Nachname_fr=Nachname 
	
	
	
	saveold erf_fribourg_all.dta, version(12) replace
			
	
	
	****************************************
	* (E) MERGE ALL DATASETS
	****************************************
		
	* (i) Start with LU data
	
	cd "$datapath\5 fuzzy merge back in"
	use erf_luzern_all.dta, clear
	
	gen id_Stata_fr=id_Stata_lu
	

	* (ii) Merge FR data

	merge 1:1  id_Stata_fr Wahljahr Nachname using erf_fribourg_all.dta, gen(merge_fribourg)
	
	// check cases which were changed by Jana Jarck (this section was added later, 10.10.2017)
	
	gen changed=0
	replace changed=1 if Vorname_fr!=Vorname
	replace changed=1 if Nachname_fr!=Nachname
	replace changed=1 if Geburtsjahr_fr!=Geburtsjahr
	
	br Wahljahr Kanton  Vorname* Nachname* Geburt* if changed==1
	 
	
	
	* (iii) Rename first names (Patrick "Pat", Maxximiliano "Max" etc. as well as double entry Werner Müller)
	
	replace Vorname="Patrick" if  id_Stata_fr=="BS-2011-0068" &	Geschlecht=="M" &	Geburtsjahr==1983 & Wahljahr==2011
	replace Vorname="Massimiliano" if  id_Stata_fr=="TI-2007-0003" &	Geschlecht=="M" &	Geburtsjahr==1982 & Wahljahr==2011
	replace Vorname="Sergio" if  id_Stata_fr=="TI-2011-0001" &	Geschlecht=="M" &	Geburtsjahr==1965 & Wahljahr==2011
	replace Vorname="Werner 2" if  id_Stata_fr=="ZH-1975-0369" &	Geschlecht=="M" &	Geburtsjahr==1932 & Wahljahr==1975
	replace Vorname="Willi" if  id_Stata_fr=="ZH-1963-0298" &	Geschlecht=="M" &	Geburtsjahr==1922 & Wahljahr==1975

	* (iv) Merge data from Vasco and Jonas
	
	merge 1:1 Vorname Nachname Kanton Geschlecht Geburtsjahr Listenbezeichnung Wahljahr using erf_jonas_all.dta, gen(merge_jonas)
	merge 1:1 Vorname Nachname Kanton Geschlecht Geburtsjahr Listenbezeichnung Wahljahr using erf_vasco_all.dta, gen(merge_vasco)
	
	//br Vorname Nachname Kanton Geschlecht Geburtsjahr Listenbezeichnung Wahljahr  if merge_jonas==2
	//br  Vorname Nachname Kanton Geschlecht Geburtsjahr  Wahljahr if Nachname=="Dollenmeier" 

	/*
	erase data_federal_temp.dta 
	erase erf_luzern_all.dta 
	erase erf_fribourg_all.dta 
	erase erf_vasco_all.dta  
	erase erf_jonas_all.dta 
    */
	
	
	* (vii) Drop 53 Vereinzelte in cantons UR, AI, GL, OW
	
	drop if Vorname==""

	
	* (vii) Tag "unsicher" cases
	
	gen unsicher=""
	replace unsicher="lu" if unsicher_lu=="0" | unsicher_lu=="10"
	replace unsicher=unsicher+"fr" if unsicher_fr=="0" | unsicher_fr=="10" 
	replace unsicher=unsicher+"vasco" if unsicher_vasco=="0" | unsicher_vasco=="10"
	replace unsicher=unsicher+"jonas" if unsicher_jonas=="0" | unsicher_jonas=="10"
	
		
	local coders lu fr vasco jonas
	foreach coder in `coders'{
	gen helper`coder'=0
	replace helper`coder'=1 if  unsicher_`coder'!="." & unsicher_`coder'!=""
	bysort id_Final_`coder': egen unsicher_num1_`coder'=max(helper`coder')
	bysort id_Stata_`coder': egen unsicher_num2_`coder'=max(helper`coder')
	}
	
	drop helper*

	gen unsicher_all=unsicher_num1_lu+ unsicher_num1_fr +unsicher_num1_vasco +unsicher_num1_jonas+unsicher_num2_lu+ unsicher_num2_fr +unsicher_num2_vasco +unsicher_num2_jonas
	tab unsicher_all
	
	
	
	
	* (vii) Tag inconsistent cases
	
	local coders lu fr vasco jonas
	foreach coder in `coders'{
	egen id_Finalnum_`coder'=group(id_Final_`coder')
	}	
	
	
	local coders lu fr vasco jonas
	local coders2 lu fr vasco jonas
	foreach coder in `coders'{
	foreach coder2 in `coders2'{
	bysort id_Finalnum_`coder': egen id_Final_sd_`coder'_`coder2'=sd(id_Finalnum_`coder2')
	}
	}
	
	
	** replace all cases that involve vasco or jonas for the years before 1975
	
	bysort id_Final_lu: egen Wahljahr_min_lu=min(Wahljahr)
	bysort id_Final_fr: egen Wahljahr_min_fr=min(Wahljahr)
	
	local vars id_Final_sd_lu_vasco id_Final_sd_lu_jonas id_Final_sd_fr_vasco id_Final_sd_fr_jonas id_Final_sd_vasco_vasco id_Final_sd_vasco_jonas id_Final_sd_vasco_lu id_Final_sd_vasco_fr id_Final_sd_jonas_lu id_Final_sd_jonas_fr id_Final_sd_jonas_vasco id_Final_sd_jonas_jonas
	foreach var in `vars'{
	tab Wahljahr  if `var'>0 & !missing(`var')& Wahljahr_min_lu<1975 |Wahljahr_min_fr<1975 & id_Final_fr!="."& id_Final_lu!="."
	replace `var'=0 if `var'>0 & !missing(`var')& Wahljahr_min_lu<1975 |Wahljahr_min_fr<1975 & id_Final_fr!="."& id_Final_lu!="."
	}
	
	
	gen inconsistent=""
	local exclude 
	
	local coders lu fr vasco jonas
	local coders2 lu fr vasco jonas
	foreach coder in `coders'{
	display "`coders2'"
	foreach coder2 in `coders2'{
	bysort id_Final_`coder': egen helper=max(id_Final_sd_`coder'_`coder2')
	gen helper2="`coder'_`coder2'" if helper>0 & !missing(helper) & "`coder'"!="`coder2'"
	replace inconsistent=inconsistent+"/"+helper2 if helper>0 & !missing(helper) & "`coder'"!="`coder2'"
	drop helper*
	//local exclude `exclude' `coder' 
	//local coders2: list coders2 - exclude
	}
	}
	
	replace inconsistent=substr(inconsistent,2,.)
	replace inconsistent="" if inconsistent=="vasco_lu/vasco_fr/jonas_lu/jonas_fr" & Wahljahr<1975
	replace inconsistent="" if inconsistent=="fr_lu/vasco_lu/vasco_fr/jonas_lu/jonas_fr" & Wahljahr<1975
	replace inconsistent="" if inconsistent=="lu_fr/fr_lu/vasco_lu/vasco_fr/jonas_lu/jonas_fr" & Wahljahr<1975
	replace inconsistent="" if inconsistent=="fr_lu/vasco_lu/vasco_fr/jonas_lu/jonas_fr" & Wahljahr<1975
	


	
	split inconsistent, p("_" "/")	
	
	gen deviation_lu=0
	gen deviation_fr=0
	gen deviation_vasco=0
	gen deviation_jonas=0
	
	
	ds inconsistent* // ssc install ds
	foreach var in `r(varlist)'{
	replace deviation_lu=deviation_lu+1 if `var'=="lu"
	replace deviation_fr=deviation_fr+1 if `var'=="fr"
	replace deviation_vasco=deviation_vasco+1 if `var'=="vasco"
	replace deviation_jonas=deviation_jonas+1 if `var'=="jonas"
	}
	
	tab deviation_lu if Wahljahr>=1975
	tab deviation_fr if Wahljahr>=1975
	tab deviation_vasco if Wahljahr>=1975
	tab deviation_jonas if Wahljahr>=1975
	

	gen inconsistent_new=inconsistent
	split inconsistent_new, p("/")	

	gen deviation_lu_fr=0
	gen deviation_lu_vasco=0
	gen deviation_lu_jonas=0
	gen deviation_fr_vasco=0
	gen deviation_fr_jonas=0
	gen deviation_vasco_jonas=0
	
	ds inconsistent_new* // ssc install ds
	foreach var in `r(varlist)'{
	replace deviation_lu_fr=1 if `var'=="lu_fr" | `var'=="fr_lu"
	replace deviation_lu_vasco=1 if `var'=="lu_vasco" | `var'=="vasco_lu"
	replace deviation_lu_jonas=1 if `var'=="lu_jonas" | `var'=="jonas_lu"
	replace deviation_fr_vasco=1 if `var'=="fr_vasco" | `var'=="vasco_fr"
	replace deviation_fr_jonas=1 if `var'=="fr_jonas" | `var'=="jonas_fr"
	replace deviation_vasco_jonas=1 if `var'=="vasco_jonas" | `var'=="jonas_vasco"
	}
	
	tab deviation_lu_fr if Wahljahr>=1975
	tab deviation_lu_vasco  if Wahljahr>=1975
	tab deviation_lu_jonas  if Wahljahr>=1975
	tab deviation_fr_vasco  if Wahljahr>=1975
	tab deviation_fr_jonas  if Wahljahr>=1975
	tab deviation_vasco_jonas if Wahljahr>=1975
	
	
	gen deviation_all=0
	ds deviation_* // ssc install ds
	foreach var in `r(varlist)'{
	replace deviation_all=deviation_all+1 if `var'==1
	}
	
	
	order unsicher inconsistent deviation*  id_Final_lu id_Final_fr id_Final_jonas id_Final_vasco Vorname Nachname Wahljahr Geschlecht Geburtsjahr Listenbezeichnung Gemeinde  Kommentar* Beruf Kanton 

	gen check=0
	replace check=1 if inconsistent1!=""
	replace check=1 if unsicher!=""

	
	
	gen Vorname_3=substr(Vorname, 1,3)
	gen Nachname_3=substr(Nachname, 1,3)
	
	replace Listenbezeichnung = subinstr(Listenbezeichnung, "Liste ", "", .)
	replace Listenbezeichnung = subinstr(Listenbezeichnung, "Lista ", "", .)
	
	forvalues i = 1(1)32 {
	replace Listenbezeichnung = subinstr(Listenbezeichnung, "`i'", "", .)
	}
	
	replace Listenbezeichnung = subinstr(Listenbezeichnung, ".", "", .)
	
	//replace Listenbezeichnung=substr(Listenbezeichnung, 1,30)
	
	

	// Kanton finalIDs Partei Geburtsjahr Gemeinde klein
	
	rename id_Final_lu id1
	rename id_Final_fr id2
	rename id_Final_vasco id3
	rename id_Final_jonas  id4

	
	gsort  Kanton Vorname_3 Nachname_3 Wahljahr Geschlecht id1 
	gen unsicher_Final=""
	gen id_Final=""

	order check Kanton Vorname Nachname id_Final unsicher_Final id1 id2 id3 id4 Listenbezeichnung Geburtsjahr  Gemeinde Wahljahr unsicher inconsistent

	local i=1	
	local coders lu fr  vasco jonas
	foreach coder in `coders'{
	replace unsicher = subinstr(unsicher, "`coder'", "`i'", .)
	replace inconsistent = subinstr(inconsistent, "`coder'", "`i'", .)
	local i=`i'+1
	}

	
	// generate outfile for Jana Jarck and Roxane Bründler who worked on that from July until September 2017
	
	preserve 
	keep 	 check Kanton Vorname Nachname id_Final unsicher_Final id1 id2 id3 id4 Listenbezeichnung Geburtsjahr  Gemeinde Wahljahr unsicher inconsistent
	export excel "$datapath\6 fuzzy merge final check\Final_Check.xlsx", firstrow(variables) replace  
	restore 
	

	
	//br Vorname Nachname Kanton Geschlecht Geburtsjahr Listenbezeichnung Gemeinde Wahljahr if merge_jonas==2

			
	
	*******************************************************************
	* (F) AGGREGATE FINAL CHECK FILES BY JANA JARCK AND ROXANE BRÜNDER
	*******************************************************************

	
	* (i) Merge data by JJ and RB
	
	import excel "$datapath/9 fuzzy merge final check/Final_Check_2.0_JJ.xlsx", first clear
	
	gen id_Ultimate=id1
	replace id_Ultimate=id_Final if id_Final!=""
	rename  id_Ultimate id_Ultimate_JJ
	rename unsicher_Final unsicher_Ultimate_JJ
	rename Kommentar_Final Kommentar_Ultimate_JJ
	
	preserve 
	import excel "$datapath/9 fuzzy merge final check/Final_Check_RB.xlsx", first clear
	gen id_Ultimate=id1
	replace id_Ultimate=id_Final if id_Final!=""
	rename  id_Ultimate id_Ultimate_RB
	rename unsicher_Final unsicher_Ultimate_RB
	rename Kommentar_Final Kommentar_Ultimate_RB

	saveold "$datapath/9 fuzzy merge final check/Final_Check_RB.dta", replace version(13)
	restore 
	

	

	//drop Vorname Nachname Geburtsjahr // drop Vorname and Nachname by Jana Jarck b/c she changed those fields
	
	merge 1:1 id1 id2 Wahljahr using "$datapath/9 fuzzy merge final check/Final_Check_RB.dta"
	drop id1 id2 id3 id4 id_Final unsicher
	

	
	* (ii) Tag inconsistent cases 
	
	egen id_Ultimate_JJ_num=group(id_Ultimate_JJ)
	egen id_Ultimate_RB_num=group(id_Ultimate_RB)
	
	bysort id_Ultimate_JJ_num: egen help1=sd(id_Ultimate_RB_num)
	bysort id_Ultimate_RB_num: egen help2=sd(id_Ultimate_JJ_num)
	
	replace help1=1 if !missing(help1) & help1!=0
	replace help2=1 if !missing(help2) & help2!=0
	
	bysort id_Ultimate_JJ_num: egen check_JJ1=max(help1)
	bysort id_Ultimate_JJ_num: egen check_JJ2=max(help2)
	bysort id_Ultimate_RB_num: egen check_RB1=max(help1)
	bysort id_Ultimate_RB_num: egen check_RB2=max(help2)
	
	gen toCheck_inconsistent=0
	replace toCheck_inconsistent=1 if check_JJ1==1 |check_JJ2==1 |check_RB1==1 |check_RB2==1
		
	gen toCheck_unsicher=0
	replace toCheck_unsicher=1 if unsicher_Ultimate_JJ>0 & !missing(unsicher_Ultimate_JJ)
	replace toCheck_unsicher=1 if unsicher_Ultimate_RB>0 & !missing(unsicher_Ultimate_RB)

	order toCheck* id_Ul*  unsicher_Ultimate* Kommentar_Ultimate* 
	sort id_Ultimate_JJ	id_Ultimate_RB Wahljahr
	
	* (iii) Manual corrections
	
	br if toCheck_inconsistent==1 | toCheck_unsicher==1

	gen id_Ultimate=id_Ultimate_JJ  // set ultimate ID by Jana Jarck as benchmark
	
	// correct entries (RB is right)
	
	replace id_Ultimate="BEJU-1951-7005" if id_Ultimate_JJ=="BEJU-1963-7006" & Wahljahr==1963 // Adolf von Allmen
	replace id_Ultimate="BEJU-1959-7024" if id_Ultimate_JJ=="BEJU-1963-7025" & Wahljahr==1963 // Louis Zihlmann
	replace id_Ultimate="BEJU-2003-7011" if id_Ultimate_JJ=="BEJU-2007-7010" & Wahljahr==2007 // Andrea Kämpf/Wegmann
	replace id_Ultimate="BL-1931-9012" if id_Ultimate_JJ=="BL-1939-0019" & Wahljahr==1947 // (Theo) Gottlieb Rutschi
	replace id_Ultimate="BL-1931-9012" if id_Ultimate_JJ=="BL-1939-0019" & Wahljahr==1951 // (Theo) Gottlieb Rutschi
	replace id_Ultimate="BL-1971-7034" if id_Ultimate_JJ=="BL-1983-7035" & Wahljahr==1983 // Otto Amrein
	replace id_Ultimate="LU-1975-7053" if id_Ultimate_JJ=="LU-1979-7054" & Wahljahr==1979 // Maria Caminati
	replace id_Ultimate="VD-1991-7092" if id_Ultimate_JJ=="VD-2015-7093" & Wahljahr==2015 // Sylviane Klein
	replace id_Ultimate="ZH-1955-7149" if id_Ultimate_JJ=="ZH-1963-7148" & Wahljahr==1963 // Siegfried Pfister
	replace id_Ultimate="ZH-1955-7149" if id_Ultimate_JJ=="ZH-1963-7148" & Wahljahr==1967 // Siegfried Pfister
	replace id_Ultimate="ZH-1979-0651" if id_Ultimate_JJ=="ZH-2007-0793" & Wahljahr==2007 // Werner Zücker
	replace id_Ultimate="ZH-1979-9225" if id_Ultimate_JJ=="ZH-2007-0413" & Wahljahr==2007 // Hans Läubli
	replace id_Ultimate="ZH-1979-9225" if id_Ultimate_JJ=="ZH-2007-0413" & Wahljahr==2011 // Hans Läubli
	replace id_Ultimate="ZH-1979-9225" if id_Ultimate_JJ=="ZH-2007-0413" & Wahljahr==2015 // Hans Läubli
	replace id_Ultimate="ZH-1983-0484" if id_Ultimate_JJ=="ZH-2003-0824" & Wahljahr==2003 // Eugène Suter
	replace id_Ultimate="ZH-1987-7124" if id_Ultimate_JJ=="ZH-1991-7125" & Wahljahr==1991 // Gabrielle/Gabi Petri
	replace id_Ultimate="ZH-1987-7124" if id_Ultimate_JJ=="ZH-1991-7125" & Wahljahr==1995 // Gabrielle/Gabi Petri
	replace id_Ultimate="ZH-1987-7131" if id_Ultimate_JJ=="ZH-1995-7130" & Wahljahr==1995 // Hans Gregor Egloff
	replace id_Ultimate="ZH-1987-7131" if id_Ultimate_JJ=="ZH-1995-7130" & Wahljahr==1999 // Hans Gregor Egloff
	replace id_Ultimate="ZH-1991-0308" if id_Ultimate_JJ=="ZH-1995-0332" & Wahljahr==1995 // Markus Hug
	replace id_Ultimate="ZH-1983-0481" if id_Ultimate_JJ=="ZH-1991-9257" & Wahljahr==1991 // Rolf Strasser
	replace id_Ultimate="ZH-1983-0481" if id_Ultimate_JJ=="ZH-1991-9257" & Wahljahr==1995 // Rolf Strasser
	replace id_Ultimate="ZH-1995-7139" if id_Ultimate_JJ=="ZH-1999-7140" & Wahljahr==1999 // Linda Frauenfelder-Bornhauser
	replace id_Ultimate="ZH-1995-7139" if id_Ultimate_JJ=="ZH-1999-7140" & Wahljahr==2003 // Linda Frauenfelder-Bornhauser
	replace id_Ultimate="ZH-1991-9133" if id_Ultimate_JJ=="ZH-1995-9132" & Wahljahr==1995 // Sylvia Frei-Maurer
	replace id_Ultimate="ZH-1991-7105" if id_Ultimate_JJ=="ZH-1995-7106" & Wahljahr==1995 // Andreas Scheu
	replace id_Ultimate="ZH-1991-7105" if id_Ultimate_JJ=="ZH-1995-7106" & Wahljahr==2003 // Andreas Scheu
	replace id_Ultimate="ZH-1991-0347" if id_Ultimate_JJ=="ZH-2003-0435" & Wahljahr==2003 // Heidi Kaufmann
	replace id_Ultimate="ZH-1971-0211" if id_Ultimate_JJ=="ZH-2007-0355" & Wahljahr==2007 // Heinz Keller
	
	// change variable names for merge to election results
	
	rename Vorname firstname
	rename Nachname name
	rename Kanton canton
	rename Listenbezeichnung list
	rename Gemeinde municipality
	rename Geburtsjahr birthyear 
	rename Wahljahr year

	// correct wrong changes by Jana Jarck

	replace birthyear=1945 if canton=="VD" & year==1987 & name=="Jaccard" & firstname=="Marianne"
	replace birthyear=1945 if canton=="VD" & year==1991 & name=="Jaccard" & firstname=="Marianne"
	replace birthyear=1961 if canton=="ZH" & year==1995 & name=="Holzreuter" & firstname=="Daniel"
	replace birthyear=1873 if canton=="ZH" & year==1935 & name=="Ryffel-Schiess" & firstname=="Albert"
	replace firstname="Willy" if canton=="ZH" & year==1975 & name=="Walker" & firstname=="Willi"
		
	gen canton_merge=canton
	replace canton_merge="BEJU" if canton=="BE" | canton=="JU"

	gen name_merge=name // undo correct changes by Jana Jarck (in 1st check round as saved in  erf_luzern_all.dta) for merge to election results
	gen firstname_merge=firstname
	gen birthyear_merge=birthyear
		
	replace birthyear_merge=1955 if canton=="BS" & year==1991 & name=="Morin" & firstname=="Guy"
	replace birthyear_merge=1901 if canton=="GE" & year==1963 & name=="Revaclier" & firstname=="François"
	replace name_merge="Haffter" if canton=="TG" & year==1971 & name=="Hafter" & firstname=="Arthur"
	replace birthyear_merge=1890 if canton=="TG" & year==1939 & name=="Holliger" & firstname=="Hans"
	replace birthyear_merge=1874 if canton=="TI" & year==1943 & name=="Varesi" & firstname=="Giovanni"
	replace birthyear_merge=1946 if canton=="VD" & year==2015 & name=="Glatz-Studer" & firstname=="Françoise"
	replace firstname_merge="Regina" if canton=="ZH" & year==1995 & name=="Aeppli" & firstname=="Regine"
	replace birthyear_merge=1908 if canton=="ZH" & year==1951 & name=="Arnold" & firstname=="Max"
	replace birthyear_merge=1942 if canton=="ZH" & year==1983 & name=="Bautz" & firstname=="Rudolf"
	replace firstname_merge="Stephan" if canton=="ZH" & year==2011 & name=="Dollenmeier" & firstname=="Stefan"
	replace birthyear_merge=1898 if canton=="ZH" & year==1939 & name=="Hofmann" & firstname=="Paul"
	replace birthyear_merge=1949 if canton=="ZH" & year==1979 & name=="Keller" & firstname=="Christoph"
	replace firstname_merge="Margrith" if canton=="ZH" & year==1999 & name=="Kolp" & firstname=="Margrit"
	replace firstname_merge="Karl" if canton=="ZH" & year==1931 & name=="Läuchli" & firstname=="Carl"
	replace birthyear_merge=1942 if canton=="ZH" & year==1979 & name=="Scherr" & firstname=="Niklaus"
	replace birthyear_merge=1943 if canton=="ZH" & year==1983 & name=="Schmid" & firstname=="Brigitte"
	replace birthyear_merge=1943 if canton=="ZH" & year==1987 & name=="Schmid" & firstname=="Brigitte"
	replace birthyear_merge=1953 if canton=="ZH" & year==1991 & name=="Schreiber" & firstname=="Kurt"
	replace birthyear_merge=1908 if canton=="ZH" & year==1975 & name=="Wolfensberger" & firstname=="Christoph"
	replace name_merge="Degea" if canton=="BS" & year==1967 & name=="Degen" & firstname=="Benjamin"
	replace name_merge="Stehler" if canton=="BS" & year==1967 & name=="Stebler" & firstname=="Hans"
	replace name_merge="Corbaz" if canton=="GE" & year==1967 & name=="Corbat" & firstname=="Fernand"
	replace name_merge="Huggli" if canton=="VD" & year==1971 & name=="Hugli" & firstname=="Jean"
	replace name_merge="Bosshardt" if canton=="ZH" & year==1939 & name=="Bosshart" & firstname=="Eduard"
	replace name_merge="Bänniger" if canton=="ZH" & year==1935 & name=="Bänninger" & firstname=="Alfred"
	replace name_merge="Günthart" if canton=="ZH" & year==1939 & name=="Günthard" & firstname=="Alois"
	replace name_merge="Kong" if canton=="ZH" & year==1963 & name=="König" & firstname=="Adolf"
	replace name_merge="Sturni" if canton=="ZH" & year==1979 & name=="Sturny" & firstname=="Max"



	keep canton canton_merge firstname name list birthyear municipality year id_Ultimate ///
		birthyear_merge firstname_merge name_merge

	saveold "$datapath/FuzzyMerge_Ultimate.dta", replace version(13)
	
	
	

	
