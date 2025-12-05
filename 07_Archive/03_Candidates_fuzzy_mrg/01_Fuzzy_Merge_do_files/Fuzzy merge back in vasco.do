
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

	
	****************************************
	* (A) VASCO
	****************************************
		
	* (i) Start with AG
		
	set more off			
	cd "$datapath\5 fuzzy merge back in vasco"
	import excel  "erf_vasco_AG.xlsx"	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	saveold "$datapath\data_federal_temp.dta", replace	version(12)


	* (ii) Merge all other files in folder endgültig fertig

	cd "$datapath\5 fuzzy merge back in vasco"
	quietly fs *		// ssc install fs
	display `r(files)'

	foreach datei in `r(files)'{ // all cantons except ZH
	//display  "`datei'"
	cd "$datapath\5 fuzzy merge back in vasco"
	import excel  `datei'	, sheet(Sheet1) clear first
	tostring id_recoded, replace
	tostring Kommentar, replace
	tostring unsicher, replace
	cd "$datapath"
	append using data_federal_temp.dta
	saveold "$datapath\data_federal_temp.dta", replace	version(12)
	}
	
	gen canton=substr(id_Stata,1,2)
	
	duplicates drop Vorname Nachname Geschlecht Geburtsjahr Listenbezeichnung Gemeinde canton Wahljahr, force // drop AG observations

	saveold erf_vasco_all.dta, version(12) replace
	
	
	****************************************
	* (B) MR JONES
	****************************************
		
	cd "$datapath\6 fuzzy merge back in 2015 jonas"
	
	import excel  erf_2015_BEJU.xlsx	, sheet(Sheet1) clear first
	
	local vars indicator candidateid_threshold_min candidateid_threshold_opt candidateid_threshold_max candidateid_ candidateid_threshold_opt_2011 candidateid_threshold_min_2011 candidateid_threshold_max_2011 candidateid__2011  id_Stata id_recoded_2011 unsicher Kommentar  indicator_2011 id_Stata_2011 unsicher_2011 Kommentar_2011 Listenbezeichnung_2011
	foreach var in `vars'{
	rename  `var' `var'_j
	}
	
	cd "$datapath"
	saveold erf_jonas_be.dta, version(12) replace
	
	****************************************
	* (C) MERGE BOTH DATASETS
	****************************************
		
	cd "$datapath"
	use erf_vasco_all.dta, clear
	
	merge 1:1 Vorname Nachname Geschlecht Geburtsjahr Listenbezeichnung Gemeinde Wahljahr using erf_jonas_be.dta, gen(merge_vasco_jonas)
	
	
	br Vorname Nachname Geschlecht Geburtsjahr Listenbezeichnung Gemeinde Wahljahr id_Stata id_recoded id_Stata_j id_recoded_2011_j  if merge_vasco_jonas==3

	
	
	br if id_Stata_j=="BEJU-2003-0445"
