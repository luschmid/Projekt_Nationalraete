
	/* TABLE OF CONTENTS 
	
	PART I: FUZZY MERGE INCLUDING 2015 DATA
	
	A: SET OVERALL PARAMETERS
	B: READ-IN ELECTION RESULTS DATA (UNTIL 2015)
	C: READ-IN AND MERGE LEFT-RIGHT-CODING
	D: DATA CLEANSING 2015 DATA
	E: GENERATE A DATASET FOR EACH CANTON
	F: SET PARAMETERS FOR EACH CANTON AND DO FUZZY MERGE
	G: ERASE INTERMEDIATE FILES
	
	PART II: MERGE 2015 DATA WITH 2011 DATA
	
	H: READ-IN RECODED DATA UNTIL 2011 (FUZZY MERGE)
	E: MERGE CANDIDATE FILE (UNTIL 2015) WITH FUZZY MERGE (UNTIL 2011) AND LEFT-RIGHT CODING  
	*/
	
	
	//////////////////////////////////////////////////////
	/// 	PART I: FUZZY MERGE INCLUDING 2015 DATA
	//////////////////////////////////////////////////////
	
	
	
	********************************
	* (A) SET OVERALL PARAMETERS
	********************************

	/*** DEFINE MAIN PATH ***/
	capture cd "C:\Dropbox\Incumbency Advantage"
	if _rc==0{
		global hauptpfad "C:\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "C:\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "C:\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "C:\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "C:\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"

	}
	
    capture cd "C:\SchmidLu\Dropbox\Projekt Nationalräte"

	if _rc==0{ 
	global hauptpfad "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"
	}


	* PATH OF FUZZY MERGE FILES
	local fuzzyPath .\fuzzy



	/*** SETTINGS: SET LOCALS ***/
	* the cantons name

	local cantons AG AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH


	**************************************
	* (B) READ-IN ELECTION RESULTS DATA (UNTIL 2015)
	*****************************************


	* (i) Read in information on votes from 1975 to 2015
	
	// Note: this dataset includes information on the types of votes a candidate received in a certain municipality
	//		 -> since we are not using this information, we collapse the dataset by candidate
	
		cd "$path_original_data"
		clear
		insheet using "NRW_KANDIDATENSTIMMEN including 2015.csv",  delimiter(";")
		
		gen year=substr(nationalratswahl,4	,4)
		destring year, replace
		collapse (sum) kandidst_unver_wz kandidst_ver_wz kandidst_tot (mean) year, by( nationalratswahl listennummer kandidatennummer canton_no)
		
		gen knr=canton_no
		merge m:1 knr using cantonal_translator,gen(merge_cantonaltranslator)
		gen kantonsnummer=knr
	
	
	* (ii) Read in information on a candidate's name, party, list, residence municipality, gender, previous position (incumbent yes/no)
		
		
		preserve
		clear
		import excel "NRW_KANDIDATEN including 2015.xlsx",  first // Note: I (Lukas) save the orginal csv file "NRW_KANDIDATEN including 2015.csv" as "NRW_KANDIDATEN including 2015.xlsx" b/c  the import excel command handles Umlaute correctly (unlike the command insheet)
		rename _all, lower // transform all variables to lower case variables

		gen year=substr(nationalratswahl,4	,4)
		destring year, replace

		saveold "NRW_KANDIDATEN including 2015.dta", replace version(12)
		restore
		
		
	* (iii)	Merge vote data with candidate information and check not merged cases
		
		merge 1:1 year kantonsnummer kandidatennummer listennummer using "NRW_KANDIDATEN including 2015.dta", gen(merge_candinfos)

		preserve
		keep if   merge_candinfos==1 |  merge_candinfos==2  // Note: these are probably "übrige" in the small cantons UR, OW, NW, GL, AI, UR
		saveold "NRW_KANDIDATEN including 2015 - NOT MERGED.dta", replace version(12)
		restore
		

		
	* (iv) Merge information on birthyear
		
		preserve
		cd "$path_original_data"
		clear
		import excel  "	including birthyear.xlsx",  first
		rename *, lower
		gen year=substr(nationalratswahl,4	,4)
		destring year, replace
		keep year kantonsnummer kandidatennummer listennummer geburtsjahr 
		saveold "NRW_KANDIDATEN including birthyear.dta", replace version(12)
		restore
		
		merge 1:1 year kantonsnummer kandidatennummer listennummer using "NRW_KANDIDATEN including birthyear.dta", gen(merge_birthyear)

	
		cd "$datapath"
		saveold "NRW_KANDIDATEN_ALL including 2015.dta", replace version(12)

	
	
	
	******************************************************
	* (C) READ-IN LEFT-RIGHT-CODING AND MERGE IT TO DATA
	*******************************************************
	
	
	use "NRW_KANDIDATEN_ALL including 2015.dta", clear


	* (i) Read-in left-right coding

	preserve
	cd "$datapath"
	clear
	import excel using "LeftRightCoding.xlsx", firstrow
	save LeftRightCoding.dta, replace
	restore
	
	* (ii) Merge votes to left-right coding
	
	gen partei_nr=parteinummer

	merge m:1 partei_nr using LeftRightCoding.dta, gen(merge_leftright)	
	
	// (iii) check merge
	 /*
	 bysort partei_nr: gen indi=_n
	 br  partei_nr if merge_leftright==1 &  indi==1
	 br  partei_nr parteiname listenname if merge_leftright==1 &  indi==1
	*/
	
	
	
		
	******************************************************
	* (D) DATA CLEANSING 2015 DATA
	*******************************************************
	
	
	gen gemeinde=gemeindename
	gen listen_nr_offiziell=listennummer
	gen wahljahr=year
	
	replace canton="BEJU" if canton=="BE" |canton=="JU"
	
	* (i) Correct entries which cannot be handled by STATA (nicknames in "") or are wrong entries (gender of candidate 
	
	replace vorname="Patrick" if canton=="BS" & wahljahr==2011 	& nachname=="Mächler" // this is the only Mächler in BS in 2011
	replace vorname="Sergio" if canton=="TI" & wahljahr==2011 & gemeinde=="Isone" & nachname=="Arigoni"
	replace vorname="Massimiliano" if canton=="TI" & wahljahr==2011 & gemeinde=="Bellinzona" & nachname=="Ay"
	replace geschlecht="F" if canton=="BS" & wahljahr==1987 &  nachname=="Frei" & vorname=="Saskia"
	
	replace vorname="Werner 2" if nachname=="Müller" & vorname =="Werner" & year==1975 & listen_nr_offiziell=="19" // Werner Müller born in 1932 and candidate for "Neue Demokrat. Bewegung" in ZH in the year 1975 

		
	* (ii) Drop stille Wahlen and Vereinzelte
	
	br if kandidst_unver_wz==. 		// Stille Wahlen
	drop if kandidst_unver_wz==.

	br canton wahljahr vorname nachname geschlecht gemeinde listenname kandidst_unver_wz kandidst_ver_wz if nachname=="" // Vereinzelte
	drop if nachname=="" 

	cd "$datapath"
	saveold "fuzzy_merge_2015_all.dta", replace version(12)
	
	
	
	
	
	*****************************************
	* (E) GENERATE A DATASET FOR EACH CANTON
	*****************************************
	
	
	cd "$datapath"
	use fuzzy_merge_2015_all.dta, clear

		
		// Note: We merge candidates in Berne and Jura since these two cantons were together until 1979

		local cantons AG AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH


		foreach can in `cantons'{

		preserve
		display "`can'"
		drop if nachname=="Vereinzelte"
		keep if canton=="`can'" & year>=1975
		cd "$path_fuzzy_data"
		rename nachname lastname
		rename vorname firstname 
		rename LinksRechts LeftRight
		rename geschlecht gender
		rename geburtsjahr birthyear
		//gen birthyear=.
		rename partei party	
		rename  gemeindenummer_bfs town_nr
		rename  gemeindename town
		rename listenname listname  
		rename listennummer list_nr  
		sort lastname firstname party listen_nr_offiziell kandidatennummer
		gen candidateid=_n // individual unique id
		save `can'_all.dta, replace
		restore

			}
			



		**************************************
		* (F) SET PARAMETERS FOR EACH CANTON
		*****************************************
		
	/*Guten Nachmittag Herr Schmid (Mail vom 09.07.2015 14:22)
			Hier eine kleine Liste wann es wo Stille Wahlen gab:
			NRW1979               AR          P            
			NRW1987               AR          P            
			NRW1999               OW          M                            
			NRW2007               NW          M            
			Freundliche Grüsse
			Corinne Straub*/
		

	

	* (ii) Do fuzzy merge canton by canton
	
	global canton "AG"
	do "$path_fuzzy/fuzzymergeSettings.do"
	cd "$path_fuzzy_data"
	
		//cd "$path_fuzzy"
	do "$path_fuzzy/fuzzymergeNEW_2015.do"
	
	* (iii) Do fuzzy merge in a loop 
	
	local cantons AG AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH
	foreach can in `cantons'{
	global canton `can'
	do "$path_fuzzy/fuzzymergeSettings.do"
	cd "$path_fuzzy_data"
	do "$path_fuzzy/fuzzymergeNEW_2015.do"
	}
	
	
	* (iv) Append all cantonal files
	
	cd "$path_fuzzy_data"
	import excel  "erf_AG.xlsx"	, sheet(Sheet1) clear first
	saveold temp.dta, version(12) replace
	
	local cantons AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH
	foreach can in `cantons'{
	cd "$path_fuzzy_data"
	import excel  "erf_`can'.xlsx"	, sheet(Sheet1) clear first
	append using temp.dta
	saveold temp.dta, version(12) replace
	}
	
	rename 	Vorname Nachname Wahljahr Gemeinde Listenbezeichnung  Geschlecht Geburtsjahr  , lower // transform  variables to lower case variables

	saveold fuzzy_merge_2015.dta, version(12) replace
	erase temp.dta
	


	
	

	**************************************
	* (G) ERASE INTERMEDIATE FILES
	*****************************************

	local cantons AG AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH
	foreach can in `cantons'{	
	foreach year in $years{
	display "`can'"
	cd "$path_fuzzy_data"
	erase `can'_all.dta
	erase `can'_`year'.dta
	erase data_`can'_`year'.dta
	erase `can'new_`year'.dta
	erase inter_`year'.dta
	erase inter_`year'_wide.dta
		}
		
	erase fuzzy_0.6_*
	erase fuzzy_0.7_*
	erase fuzzy_0.8_*
	erase fuzzy_0.9_*
	}
	}
	
	
		
	//////////////////////////////////////////////////////
	/// 	PART II: MERGE 2015 DATA WITH 2011 DATA
	//////////////////////////////////////////////////////
	

	
	
	****************************************
	* (H) READ-IN RECODED DATA UNTIL 2011
	****************************************
		
	* (i) Start with AG
		
	set more off			
	cd "$datapath\2 fuzzy merge back in\endgültig fertig"
	import excel  "erf_AG.xlsx"	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	saveold "$datapath\data_federal_temp.dta", replace	version(12)


	* (ii) Merge all other files in folder endgültig fertig

	cd "$datapath\2 fuzzy merge back in\endgültig fertig"
	quietly fs *		// ssc install fs
	display `r(files)'

	foreach datei in `r(files)'{ // all cantons except ZH
	//display  "`datei'"
	cd "$datapath\2 fuzzy merge back in\endgültig fertig"
	import excel  `datei'	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	cd "$datapath"
	append using data_federal_temp.dta
	saveold "$datapath\data_federal_temp.dta", replace	version(12)
	}


	cd "$datapath\2 fuzzy merge back in\unsicher einfügen"
	import excel  erf_ZH.xlsx	, sheet(Sheet1) clear first
	cd "$datapath"
	append using data_federal_temp.dta, force
	//saveold "$datapath\data_federal_temp.dta", replace	version(12)
	
	erase "$datapath\data_federal_temp.dta"


	* (iii) Data cleansing and variable generation

	 // drop canton AG (start canton in loop above)
	duplicates drop Nachname Vorname	Wahljahr Geschlecht	Geburtsjahr	Listenbezeichnung	Gemeinde, force 	


	* (iv) ID recoding

	//gen id_Stata_new=id_Stata
	//replace id_Stata_new=id_recoded if id_recoded!="" & id_recoded!="."
	
	* (v) lower case variables
	rename 	Vorname Nachname Wahljahr Gemeinde Listenbezeichnung  Geschlecht Geburtsjahr  , lower // transform  variables to lower case variables

	* (vi) generate canton variable
	
	gen canton=substr(id_Stata,1,2)
	
	* (vii) Merge old list names and associated list numbers (we need this for merging new vote data and fuzzy merge until 2011)

	gen liste=listenbezeichnung
	
	cd "$path_original_data"
	merge m:1 liste using ListNames_until2011.dta, gen(merge_listinfos)
	
	// check merge: 
	// br if merge_listinfos!=3
	// -> These are all observations from 1971
	
	drop if merge_listinfos!=3
	
	* (viii) Check for duplicates
	
	bysort canton wahljahr vorname nachname geschlecht gemeinde: gen indimax=_N
	br if indimax>1
	
	replace vorname="Werner 2" if id_Stata=="ZH-1975-0369" // Werner Müller born in 1932 and candidate for "Neue Demokrat. Bewegung" in ZH in the year 1975 

	
	* (ix) Replace entries from Bern and Jura
	
	replace canton="BEJU" if canton=="BE" |canton=="JU"
	
	* (x) Correct some names that have changed in the new version including 2015 data
	

	replace vorname="Šárka" if  vorname=="¿árka"
	replace nachname="Bläsi" if canton=="GE" & wahljahr==2011 & vorname=="Thomas"	
	replace vorname="Sergio" if canton=="TI" & wahljahr==2011 & gemeinde=="Isone" & nachname=="Arigoni"
	replace vorname="Massimiliano" if canton=="TI" & wahljahr==2011 & gemeinde=="Bellinzona" & nachname=="Ay"
	replace vorname="Patrick" if canton=="BS" & wahljahr==2011 	& nachname=="Mächler" // this is the only Mächler in BS in 2011
	
	
	local variables indicator	candidateid_threshold_opt	candidateid_threshold_min	candidateid_threshold_max	candidateid_	id_Stata	id_recoded	unsicher	Kommentar	listenbezeichnung	
	foreach var in `variables'{	
	rename `var' `var'_2011
	}
	
		
	
	cd "$datapath"
	saveold "ID_Recoding_until_2011.dta", replace version(12)

	
	
	
	*********************************************************************************************
	* (I) MERGE CANDIDATE FILE (UNTIL 2015) WITH FUZZY MERGE (UNTIL 2011) AND LEFT-RIGHT CODING  
	*********************************************************************************************
	
	
	* (i) Merge 2011 and 2015 data
	
	cd "$path_fuzzy_data"
	use fuzzy_merge_2015.dta, clear
	cd "$datapath"
	
	gen canton=substr(id_Stata,1,2)
	replace canton="BEJU" if canton=="BE"
	
	merge 1:1 canton wahljahr vorname nachname geschlecht gemeinde  using ID_Recoding_until_2011.dta, gen(merge_idrecoding)	
	
	
	* (ii) Figure out those cases that should be manually checked 
	
	gen tocheck=0
	
	* (a) Cases in 2015 for which candidate_id_min != candidate_id_max
	
	 replace tocheck=1 if indicator==1 & wahljahr==2015

	 * -> we also need to check all cases that are possibly related with these entries
  
	   gen help2=0
	   replace help2=1 if indicator==1 & wahljahr==2015
	   
	   local variables candidateid_threshold_min	candidateid_threshold_min_2011 candidateid_threshold_opt	candidateid_threshold_opt_2011	candidateid_threshold_max	candidateid_threshold_max_2011		 
	   foreach var in `variables'{	
	   bysort `var' canton: egen h_`var'=sum(help2)
	   replace tocheck=2 if h_`var'>0 & wahljahr!=2015
	   drop h_`var'
		}
		
	 
	 * (b) Cases before 2015 for which id_Stata changed because of the added 2015 data
	  
	  replace tocheck=3 if id_Stata!=id_Stata_2011 & wahljahr!=2015
	  
	  gen help1=0
	  replace help1=1 if id_Stata!=id_Stata_2011 & wahljahr!=2015
	 
	  
	 * -> we also need to check all cases that are possibly related with these entries
  
	
	   local variables candidateid_threshold_min	candidateid_threshold_min_2011 candidateid_threshold_opt	candidateid_threshold_opt_2011	candidateid_threshold_max	candidateid_threshold_max_2011		 
	   foreach var in `variables'{	
	   bysort `var' canton: egen h_`var'=sum(help1)
	   replace tocheck=4 if h_`var'>0 & tocheck!=3
	   drop h_`var'
		}
		
		drop help* candidateid_threshold_min_sorted
	  
	    tab tocheck
	  
	  
	    br id_Stata* id_recoded* canton wahljahr vorname nachname geschlecht gemeinde if tocheck==1   
		
	

	 * (iii) Generate Output for Jonas Röllin
	 
		 drop indicator
		 gen indicator=tocheck
	 
	 
		rename vorname Vorname
		rename nachname Nachname
		rename listenbezeichnung Listenbezeichnung
		rename geschlecht Geschlecht
		rename geburtsjahr Geburtsjahr
		rename wahljahr Wahljahr
		rename gemeinde Gemeinde
		rename listenbezeichnung_2011 Listenbezeichnung_2011
		
		drop id_recoded
		
		gen firstname_start_ = strlower(substr(Vorname,1,3))
		gen lastname_start_ = strlower(substr(Nachname,1,3))
	
		sort firstname_start_ lastname_start_  candidateid_threshold_opt
		egen candidateid_threshold_min_sorted=group(candidateid_threshold_min)
	

		sort canton  firstname_start_ lastname_start_  


		keep canton indicator candidateid_threshold_min candidateid_threshold_opt  candidateid_threshold_max candidateid_ Vorname Nachname id_recoded id_Stata unsicher Kommentar Geschlecht Geburtsjahr Listenbezeichnung Gemeinde Wahljahr canton indicator_2011 candidateid_threshold_opt_2011 candidateid_threshold_min_2011 candidateid_threshold_max_2011 candidateid__2011 id_Stata_2011 id_recoded_2011 unsicher_2011 Kommentar_2011 Listenbezeichnung_2011

		order canton indicator candidateid_threshold_min candidateid_threshold_opt candidateid_threshold_max candidateid_ candidateid_threshold_opt_2011	candidateid_threshold_min_2011	candidateid_threshold_max_2011	candidateid__2011 Vorname Nachname  id_Stata id_recoded_2011 unsicher Kommentar Geschlecht Geburtsjahr Listenbezeichnung Gemeinde Wahljahr canton indicator_2011 candidateid_threshold_opt_2011 candidateid_threshold_min_2011 candidateid_threshold_max_2011 candidateid__2011 id_Stata_2011  unsicher_2011 Kommentar_2011 Listenbezeichnung_2011
	

		*cd "$datapath"
	    local cantons AG AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH
		foreach can in `cantons'{	
		preserve
		keep if canton=="`can'"
		drop canton
		replace id_recoded_2011="" if id_recoded_2011=="."
		export excel using "./3 fuzzy merge to check 2015/erf_2015_`can'.xlsx", firstrow(variables) replace
		restore
		}


	
