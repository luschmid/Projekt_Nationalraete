
	clear
	set more off

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




	******************************************************
	* (A) CHECK FIRST NAMES
	******************************************************


	
	
		// (a) by canton
	
	
		use "$hauptpfad\1_data\nr_data_1931_1975_jan2017.dta", clear 
		
		merge m:1 canton using "$path_original_data\cantonal_translator", gen(merge_cantonal_translator)
		
		drop if merge_cantonal_translator==2 // drop JU
		
		egen vorname_num=group(vorname)
		sum vorname_num
		
		bysort vorname_num: gen vorname_total=_N
		
		collapse  vorname_num vorname_total, by(vorname canton)
		
		gen vorname_unchanged=vorname
		
		reshape wide vorname  , i(vorname_num) j(canton) string
		
		gen vorname_cantons="AG" if vornameAG!=""
		
		foreach i in AI AR BE BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH {
        replace vorname_cantons= vorname_cantons +" " + "`i'" if vorname`i'!=""
        //display "`i'"
		}
			
		keep vorname_unchanged vorname_cantons vorname_total
		
		save "$hauptpfad\firstnames_canton_31_75.dta", replace 
		
		
		// (b) by year
		
		use "$hauptpfad\1_data\nr_data_1931_1975_jan2017.dta", clear 
		
		egen vorname_num=group(vorname)
		sum vorname_num
		
		bysort vorname_num: gen vorname_total=_N
		
		collapse  vorname_num vorname_total, by(vorname year)
		
		gen vorname_unchanged=vorname
		
		reshape wide vorname  , i(vorname_num) j(year) 
		
		gen vorname_year="1931" if vorname1931!=""
		
		forvalues i = 1935(4)1975{
        replace vorname_year= vorname_year +" " + "`i'" if vorname`i'!=""
		}
			
		keep vorname_unchanged vorname_year vorname_total
		
		
		// (c) merge canton and year 
		
		merge 1:1 vorname_unchanged using "$hauptpfad\firstnames_canton_31_75.dta", gen(merge_ktn_year) 

		
		rename vorname_unchanged Vorname 
		rename vorname_year Jahre
		rename vorname_cantons Kantone 
		rename vorname_total Anzahl_Beobachtungen
		
		order Vorname Anzahl_Beobachtungen Kantone Jahre
		keep  Vorname Anzahl_Beobachtungen Kantone Jahre
		sort Vorname Anzahl_Beobachtungen Kantone Jahre
		
		export excel "$hauptpfad\1_data\Vornamen_Check_Jan17.xlsx", firstrow(variables) replace
		erase "$hauptpfad\firstnames_canton_31_75.dta"



		******************************************************
		* (B) CHECK SURNAMES
		******************************************************


		// (a) by canton
	
	
		use "$hauptpfad\1_data\nr_data_1931_1975_jan2017.dta", clear 
		
		merge m:1 canton using "$path_original_data\cantonal_translator", gen(merge_cantonal_translator)
		
		drop if merge_cantonal_translator==2 // drop JU
		
		egen nachname_num=group(nachname)
		sum nachname_num
		
		bysort nachname_num: gen nachname_total=_N
		
		collapse  nachname_num nachname_total, by(nachname canton)
		
		gen nachname_unchanged=nachname
		
		reshape wide nachname  , i(nachname_num) j(canton) string
		
		gen nachname_cantons="AG" if nachnameAG!=""
		
		foreach i in AI AR BE BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH {
        replace nachname_cantons= nachname_cantons +" " + "`i'" if nachname`i'!=""
        //display "`i'"
		}
			
		keep nachname_unchanged nachname_cantons nachname_total
		
		save "$hauptpfad\surnames_canton_31_75.dta", replace 
		
		
		// (b) by year
		
		use "$hauptpfad\1_data\nr_data_1931_1975_jan2017.dta", clear 
		
		egen nachname_num=group(nachname)
		sum nachname_num
		
		bysort nachname_num: gen nachname_total=_N
		
		collapse  nachname_num nachname_total, by(nachname year)
		
		gen nachname_unchanged=nachname
		
		reshape wide nachname  , i(nachname_num) j(year) 
		
		gen nachname_year="1931" if nachname1931!=""
		
		forvalues i = 1935(4)1975{
        replace nachname_year= nachname_year +" " + "`i'" if nachname`i'!=""
		}
			
		keep nachname_unchanged nachname_year nachname_total
		
		
		// (c) merge canton and year 
		
		merge 1:1 nachname_unchanged using "$hauptpfad\surnames_canton_31_75.dta", gen(merge_ktn_year) 


		rename nachname_unchanged Nachname 
		rename nachname_year Jahre
		rename nachname_cantons Kantone 
		rename nachname_total Anzahl_Beobachtungen
		
		order Nachname Anzahl_Beobachtungen Kantone Jahre
		keep  Nachname Anzahl_Beobachtungen Kantone Jahre
		sort Nachname Anzahl_Beobachtungen Kantone Jahre
		
		
		export excel "$hauptpfad\1_data\Nachnamen_Check_Jan17.xlsx", firstrow(variables) replace
		erase "$hauptpfad\surnames_canton_31_75.dta"
