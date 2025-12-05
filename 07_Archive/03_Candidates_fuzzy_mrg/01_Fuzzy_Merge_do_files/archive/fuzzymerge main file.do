
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
	else{ 
		global hauptpfad "C:\Users\08609901\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "C:\Users\08609901\Dropbox\Projekt Nationalräte\1 Data\1_data"\
		global path_original_data "C:\Users\08609901\Dropbox\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "C:\Users\08609901\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "C:\Users\08609901\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"

	}


	* PATH OF FUZZY MERGE FILES
	local fuzzyPath .\fuzzy



	/*** SETTINGS: SET LOCALS ***/
	* the cantons name

	local cantons AG AI AR BEJU BL BS FR GE GL GR LU NE NW OW SG SH SO SZ TG TI UR VD VS ZG ZH



	**************************************
	* (B) READ-IN DATA 
	*****************************************


	* (ii) Read in information on left-right position of parties and place of residence of candidates

	cd "$datapath"
	clear
	import excel using "LeftRightCoding.xlsx", firstrow
	save LeftRightCoding.dta, replace

	cd "$path_original_data"
	clear
	insheet using  "NRW_KANDIDATEN.csv", delimiter(";")
	cd "$datapath"
	save NRW_KANDIDATEN.dta, replace

	* (iii) Merge left-right position information 

	
	cd "$path_original_data"
	clear
	insheet using "KANDIDATENSTIMMEN (BfS).csv", delimiter(";")
	cd "$datapath"
	merge m:1 partei_nr using LeftRightCoding.dta, gen(merge_leftright)	


	* (iv) Merge place of residence

	gen nachname=name
	gen kandidatennummer=kandidaten_nr
	gen parteinummer=partei_nr
	gen listennummer=listen_nr_offiziell
	gen kantonsname=kanton
	
	drop if liste=="Vereinzelte"

	merge 1:1 nationalratswahl kantonsname listennummer 	kandidatennummer	vorname	nachname geschlecht using NRW_KANDIDATEN.dta, gen(merge_residence)
	
	* (v) Check non-merges 
	
	tab nationalratswahl if merge_residence==2 // most of non-merges are candidates who ran in 1971 for which we have no election results
	drop if nationalratswahl=="NRW1971"
	
	preserve // out for manual checking
	sort merge_residenc nationalratswahl kantonsname listennummer kandidatennummer vorname	nachname kandidst_tot
	order nationalratswahl listennummer kantonsname	kandidatennummer	vorname	nachname name geschlecht geburtsjahr partei liste kandidst_tot
	keep if inrange(merge_residence,1,2)
	export excel using residence_manualcheck.xlsx, firstrow(variables) replace
	restore
	
	/* CHECK THE FOLLOWING CANDIDATES (MAIL SENT TO BFS ON 8 JULY 2015 
	nationalratswahl	listennummer	kantonsname	kandidatennummer	vorname		nachname	geschlecht
	NRW1979				1				AR							1	Christian	Merz			M
	NRW1979				2				AR							1	Hans-Rudolf	Früh			M
	NRW1987				1				AR							1	Hans-Rudolf	Früh			M
	NRW1987				2				AR							1	Herbert		Maeder			M
	NRW1999				1				OW							1	Adalbert	Durrer			M
	NRW2007				1				NW							1	Edi			Engelberger		M
	*/


	sort kantons_nr name vorname nationalratswahl
	// br nationalratswahl kantons_nr kanton partei_nr name name vorname geburts* geschlecht gemeindenummer_bfs

	* (v) Generate canton-year substet of datasets 

	gen year=substr(nationalratswahl,4,4)
	destring year, replace

	g canton=kanton
	replace canton="BEJU" if kanton=="BE" | kanton=="JU" // take BE and JU together for merging procedure

	foreach can in `cantons'{

	preserve
	display "`can'"
	drop if name=="Vereinzelte"
	keep if canton=="`can'" & year>=1975
	cd "$path_fuzzy_data"
	rename name lastname
	rename vorname firstname 
	rename LinksRechts LeftRight
	rename geschlecht gender
	rename geburtsjahr birthyear
	rename partei party	
	rename  gemeindenummer_bfs town_nr
	rename  gemeindename town
	rename liste listname  
	rename listennummer list_nr  
	sort lastname firstname party listen_nr_offiziell kandidaten_nr
	gen candidateid=_n // individual unique id
	save `can'_all.dta, replace
	restore

		}




	**************************************
	* (C) SET PARAMETERS FOR EACH CANTON
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1983 1991 1995 1999 2003 2007 2011 // stille wahl: 1979, 1987
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2011
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
	global years 1975 1979 1983 1987 1991 1995 2003 2007 2011 // stille wahl 1999
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	global years 1975 1979 1983 1987 1991 1995 1999 2003 2007 2011
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
	* (D) ERASE INTERMEDIATE FILES
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
