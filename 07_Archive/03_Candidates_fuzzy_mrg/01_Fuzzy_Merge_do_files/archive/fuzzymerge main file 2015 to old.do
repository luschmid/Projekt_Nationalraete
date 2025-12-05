
	/* TABLE OF CONTENTS 
	
	A: SET OVERALL PARAMETERS
	B: READ-IN ELECTION RESULTS DATA (UNTIL 2015)
	C: READ-IN RECODED DATA UNTIL 2011 (FUZZY MERGE)
	D: READ-IN LEFT-RIGHT-CODING
	E: MERGE CANDIDATE FILE (UNTIL 2015) WITH FUZZY MERGE (UNTIL 2011) AND LEFT-RIGHT CODING  
	F: GENERATE A DATASET FOR EACH CANTON
	G: SET PARAMETERS FOR EACH CANTON
	H: ERASE INTERMEDIATE FILES
	*/
	
	
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
	
    capture cd "C:\SchmidLu\Dropbox\Incumbency Advantage"

	if _rc==0{ 
	global hauptpfad "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "C:\SchmidLuDropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"
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
		
		cd "$datapath"
		saveold "NRW_KANDIDATEN_ALL including 2015.dta", replace version(12)



	**************************************
	* (C) READ-IN RECODED DATA UNTIL 2011
	*****************************************
		
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
	saveold "$datapath\data_federal_temp.dta", replace	version(12)


	* (iii) Data cleansing and variable generation

	 // drop canton AG (start canton in loop above)
	duplicates drop Nachname Vorname	Wahljahr Geschlecht	Geburtsjahr	Listenbezeichnung	Gemeinde, force 	


	* (iv) ID recoding

	gen id_Stata_new=id_Stata
	replace id_Stata_new=id_recoded if id_recoded!="" & id_recoded!="."
	
	* (v) lower case variables
	rename 	Vorname Nachname Wahljahr Gemeinde Listenbezeichnung  Geschlecht, lower //Geburtsjahr , lower // transform  variables to lower case variables

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
	
	bysort canton wahljahr vorname nachname geschlecht gemeinde
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
	
	
	cd "$datapath"
	saveold "ID_Recoding_until_2011.dta", replace version(12)

	
	
	**************************************
	* (D) READ-IN LEFT-RIGHT-CODING
	*****************************************

	cd "$datapath"
	clear
	import excel using "LeftRightCoding.xlsx", firstrow
	save LeftRightCoding.dta, replace
	
	
	
	*********************************************************************************************
	* (E) MERGE CANDIDATE FILE (UNTIL 2015) WITH FUZZY MERGE (UNTIL 2011) AND LEFT-RIGHT CODING  
	*********************************************************************************************
	
	
	cd "$datapath"
	use "NRW_KANDIDATEN_ALL including 2015.dta", clear
	
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

	* (ii) Merge votes to previous fuzzy merge evaluation
	

	merge 1:1 canton wahljahr vorname nachname geschlecht gemeinde  using ID_Recoding_until_2011.dta, gen(merge_idrecoding)	
	
	
	/* Different rounds of Merge checks 
	
	// merge check 1: 
	br canton wahljahr vorname nachname geschlecht gemeinde listenbezeichnung if merge_idrecoding==2
 	br canton wahljahr vorname nachname geschlecht gemeinde listenbezeichnung if vorname=="Fritz" & nachname=="Bodenmann" & canton=="AG" & wahljahr==1975
 	br canton wahljahr vorname nachname geschlecht gemeinde listenbezeichnung if vorname=="Fritz" & nachname=="Hodler" & canton=="AG" & wahljahr==1975

	// -> Note: The variable listenbezeichnung was changed by BfS (with the 2015 data)
	// 			Solution: I merged the listennumber to the fuzzy merge (merge m:1 liste using ListNames_until2011.dta, gen(merge_listinfos))
	
	// merge check 2: 
	br canton wahljahr vorname nachname geschlecht gemeinde listenbezeichnung merge_idrecoding if merge_idrecoding==2
	br canton wahljahr vorname nachname geschlecht gemeinde merge_idrecoding if vorname=="Alphonse" &	nachname=="Froidevaux"
	br canton wahljahr vorname nachname geschlecht gemeinde merge_idrecoding if vorname=="Dominique" &	nachname=="Hubleur"
	
	// -> Note: These are mostly observations from Jura which are coded as BE in the fuzzy merge data but as JU in the new candidate data (until 2015)
	

	// merge check 3: 
	tab year if merge_idrecoding==1
	br canton wahljahr vorname nachname geschlecht gemeinde listenbezeichnung kandidst_unver_wz kandidst_ver_wz if merge_idrecoding==1& wahljahr<2015
	tab year if merge_idrecoding==2
	br canton wahljahr vorname nachname geschlecht gemeinde listenbezeichnung kandidst_unver_wz kandidst_ver_wz if merge_idrecoding==2

	
	// -> no observations, which means that all oservations from using data 
	
	br canton wahljahr vorname nachname geschlecht gemeinde listenbezeichnung merge_idrecoding if merge_idrecoding!=3 & year<2015
	
	*/
	
	
 	
	* (ii) Merge votes to left-right coding
	
	gen partei_nr=parteinummer

	merge m:1 partei_nr using LeftRightCoding.dta, gen(merge_leftright)	
	
	// check merge
	 /*
	 bysort partei_nr: gen indi=_n
	 br  partei_nr if merge_leftright==1 &  indi==1
	 br  partei_nr parteiname listenname if merge_leftright==1 &  indi==1
	*/
	
	
	
	* (iii) Drop stille Wahlen and Vereinzelte
	
	br if kandidst_unver_wz==. 		// Stille Wahlen
	drop if kandidst_unver_wz==.

	br canton wahljahr vorname nachname geschlecht gemeinde listenbezeichnung kandidst_unver_wz kandidst_ver_wz if nachname=="" // Vereinzelte
	drop if nachname=="" 

	
	* (iv) Save one dataset pre2011 and post2011
	
	/*preserve
	keep if wahljahr<=2011
	saveold fuzzy_merge_out_pre2011.dta, replace version(12)
	restore
	
	keep if wahljahr<=2011
	*/

	
	
	*****************************************
	* (F) GENERATE A DATASET FOR EACH CANTON
	*****************************************
	
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
	//rename geburtsjahr birthyear
	rename partei party	
	rename  gemeindenummer_bfs town_nr
	rename  gemeindename town
	rename liste listname  
	rename listennummer list_nr  
	sort lastname firstname party listen_nr_offiziell kandidatennummer
	gen candidateid=_n // individual unique id
	save `can'_all_until2015.dta, replace
	restore

		}



	**************************************
	* (G) SET PARAMETERS FOR EACH CANTON
	*****************************************
	
	/*Guten Nachmittag Herr Schmid (Mail vom 09.07.2015 14:22)
		Hier eine kleine Liste wann es wo Stille Wahlen gab:
		NRW1979               AR          P            
		NRW1987               AR          P            
		NRW1999               OW          M                            
		NRW2007               NW          M            
		NRW1971               ZG          P            
		Freundliche Grüsse
		Corinne Straub*/
	

	global canton VD


	if "$canton"=="AI"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	
	if "$canton"=="AR"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	
	if "$canton"=="AG"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	if "$canton"=="BEJU"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.8
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	
	if "$canton"=="BL"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.8
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	if "$canton"=="BS"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.8
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	if "$canton"=="FR"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}	
	
	if "$canton"=="GE"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.8
	global threshold_min 0.6
	global threshold_max 0.9
	}	
	
	if "$canton"=="GE"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}	
	
	if "$canton"=="GR"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.8
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	if "$canton"=="LU"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.8
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	if "$canton"=="NE"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
			
	if "$canton"=="NW"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	if "$canton"=="OW"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	if "$canton"=="SG"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.8
	global threshold_min 0.7
	global threshold_max 0.9
	}
	
	
		
	if "$canton"=="SH"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}



	
	if "$canton"=="SO"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}



	
	if "$canton"=="SZ"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	if "$canton"=="TG"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	
	if "$canton"=="TI"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	
	
	if "$canton"=="UR"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	if "$canton"=="VD"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.8
	global threshold_min 0.7
	global threshold_max 0.9
	}
	
	if "$canton"=="VS"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	
	
	if "$canton"=="ZG"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.7
	global threshold_min 0.6
	global threshold_max 0.9
	}
	
	
	if "$canton"=="ZH"{
	* cantonal years
	global years 2015
	* the matching variables
	global vars lastname firstname LeftRight gender town_nr birthyear year
	* the weights for the matching variables (for matches)
	global yesweights 15 15 12 15 13 15 1 
	* the weights for the matching variables (for non-matches)
	global noweights 8 10 11 15 8 12 1
	* define thresholds for merging at the end (colour cells to simplify task of recoders)
	global threshold_opt 0.8
	global threshold_min 0.6
	global threshold_max 0.9
	}

	cd "$path_fuzzy"
	do fuzzymergeNEW.do
	
	
	**************************************
	* (H) ERASE INTERMEDIATE FILES
	*****************************************


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
