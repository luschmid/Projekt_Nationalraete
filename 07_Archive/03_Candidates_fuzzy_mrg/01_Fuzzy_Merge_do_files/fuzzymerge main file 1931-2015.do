
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
	
    capture cd "D:\SchmidLu\Dropbox\Projekt Nationalräte"

	if _rc==0{ 
	global hauptpfad "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"
	}


	* PATH OF FUZZY MERGE FILES
	local fuzzyPath .\fuzzy



	/*** SETTINGS: SET LOCALS ***/
	* the cantons name

	local cantons AG AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH


	
	
	
	************************************************
	* (B) READ-IN ELECTION RESULTS DATA (1931-2015)
	************************************************

	cd "$path_fuzzy_data"
	use nationalraete_1931_2015_fuzzy.dta, clear
	
	replace canton="BEJU" if canton=="BE" | canton=="JU"
			
	
	*****************************************
	* (C) GENERATE A DATASET FOR EACH CANTON
	*****************************************
	
		
		// Note: We merge candidates in Berne and Jura since these two cantons were together until 1979

		local cantons AG AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH


		foreach can in `cantons'{

		preserve
		display "`can'"
		drop if name=="Vereinzelte"
		keep if canton=="`can'" & year>=1931
		rename name lastname
		rename leftright LeftRight
		rename sex gender
		rename  municipalityno town_nr
		rename  municipality town
		rename list listname  
		sort lastname firstname 
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
	
	/*global canton "AG"
	do "$path_fuzzy/fuzzymergeSettings_1931-2015.do"
	cd "$path_fuzzy_data"
	
		//cd "$path_fuzzy"
	do "$path_fuzzy/fuzzymergeNEW_1931-2015.do"
	*/
	
	* (iii) Do fuzzy merge in a loop 
	
	local cantons AG AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH
	foreach can in `cantons'{
	global canton `can'
	do "$path_fuzzy/fuzzymergeSettings_1931-2015.do"
	cd "$path_fuzzy_data"
	do "$path_fuzzy/fuzzymergeNEW_1931-2015.do"
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
	
	

	
