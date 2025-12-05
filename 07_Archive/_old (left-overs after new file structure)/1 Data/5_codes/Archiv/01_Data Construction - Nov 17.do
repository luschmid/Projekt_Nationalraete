
******************************************************
* (0) Set directories
******************************************************

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
	
    capture cd "D:\SchmidLu\Dropbox\Projekt Nationalräte"

	if _rc==0{ 
	global hauptpfad "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"
	}


		


	******************************************************
	* (A) READ-IN ELECTION RESULTS DATA (1975 until 2015)
	******************************************************


	* (i) Read in individual election results votes from 1975 to 2015
	
	// Note: this dataset includes information on the types of votes a candidate received in a certain municipality
	//		 -> since we are not using this information, we collapse the dataset by candidate
	
		cd "$path_original_data"
		clear
		insheet using "NRW_KANDIDATENSTIMMEN including 2015.csv",  delimiter(";")
		
		gen year=substr(nationalratswahl,4	,4)
		destring year, replace
		collapse (sum) kandidst_unver_wz kandidst_ver_wz kandidst_tot (mean) year, by( nationalratswahl listennummer kandidatennummer canton_no)
		
	
	* (ii) Read in information on a candidate's name, party, list, residence municipality, gender, previous position (incumbent yes/no)
		
		
		preserve
		clear
		import excel "NRW_KANDIDATEN including birthyear.xlsx",  first // Note: I (Lukas) saved the orginal csv file "NRW_KANDIDATEN including 2015.csv" as "NRW_KANDIDATEN including 2015.xlsx" b/c  the import excel command handles Umlaute correctly (unlike the command insheet)
		rename _all, lower // transform all variables to lower case variables

		gen year=substr(nationalratswahl,4	,4)
		destring year, replace
		
		do "$hauptpfad/5_codes/01_Add Missing Municipalities.do" // add missing municipality information in bfs data (based on Florence's coding, Jan 17)

		saveold "NRW_KANDIDATEN including 2015.dta", replace version(12)
		restore
		
		
		
	* (iii)	Merge vote data with candidate information and check not merged cases
	
		
		gen kantonsnummer=canton_no
		merge 1:1 year kantonsnummer kandidatennummer listennummer using "NRW_KANDIDATEN including 2015.dta", gen(merge_candinfos)
		erase "NRW_KANDIDATEN including 2015.dta"
		
		
		/* merge_candinfos==1: 53 Stimmen für Vereinzelte in Kantonen mit einem NR-Sitz (Majorzkantone) 
		   merge_candinfos==2: 6 stille Wahlen
		*/

		gen knr=kantonsnummer
		merge m:1 knr using cantonal_translator,gen(merge_cantonaltranslator)
		
		drop canton_no knr // cantonal number is kantonsnummer from now on; we drop all other canton number variables

	
	* (iv) Merge left-right coding
	
		gen partei_nr=parteinummer
		
		merge m:1 partei_nr using "$hauptpfad\1_data\LeftRightCoding.dta", gen(merge_parteien)
		
	
	* (v) Corrections
	
		replace vorname="Willy" if nachname=="Walker" & vorname=="Willi" & year==1975 & canton=="ZH"
		replace gemeindenummer_bfs=4937 if  nachname=="Möckli" & year==1975 & vorname=="Gustav" & gemeindenummer_bfs==4631 & geburtsjahr==1926
		replace gemeindenummer_bfs=195 if  nachname=="Reich" & year==1975 & vorname=="Richard" & gemeindenummer_bfs==253 & geburtsjahr==1927
		replace gemeindenummer_bfs=195 if  nachname=="Reich" & year==1979 & vorname=="Richard" & gemeindenummer_bfs==253 & geburtsjahr==1927
		replace gemeindenummer_bfs=195 if  nachname=="Reich" & year==1983 & vorname=="Richard" & gemeindenummer_bfs==. & geburtsjahr==1927
		replace gemeindenummer_bfs=195 if  nachname=="Reich" & year==1987 & vorname=="Richard" & gemeindenummer_bfs==. & geburtsjahr==1927
		
		replace gemeindename="Maur" if  nachname=="Reich" & year==1975 & vorname=="Richard" & gemeindename=="Zürich"  & geburtsjahr==1927
		replace gemeindename="Maur" if  nachname=="Reich" & year==1979 & vorname=="Richard" & gemeindename=="Zürich" & geburtsjahr==1927


		
	* (vi) Rename variables in order to merge it to pre-1975 data
	
		
		gen elected=. 
		replace elected=1 if gewaehlt=="G"
		replace elected=0 if gewaehlt=="N"

	

		keep nachname vorname kandidatennummer kandidst_tot year kantonsnummer canton listenname gemeindenummer_bfs gemeindename geschlecht geburtsjahr elected status LinksRechts
		
		
		rename kandidatennummer candidateno
		rename kandidst_tot votes
		rename kantonsnummer cantonno
		rename listenname list
		rename gemeindenummer_bfs municipalityno
		rename gemeindename municipality
		rename LinksRechts leftright
		rename geschlecht sex
		rename geburtsjahr birthyear 
		rename vorname firstname
		rename nachname name
	
	
		saveold "NRW_KANDIDATEN_ALL including 2015.dta", replace version(12)


		

		******************************************************
		* (B) READ-IN ELECTION RESULTS DATA (1931 until 1975)
		******************************************************
		
		
		* (i) Import all excel files
		
		clear
		import excel "$hauptpfad/6_Bundesblatt/Bundesblatt_in_FR_checked/1931_Bundesblatt_backin.xlsx"	,  first // start with 1931
		gen year=1931
		replace canton=canton[_n-1] if canton==""
		saveold "$hauptpfad\1_data\temp.dta", replace
		
		forvalues y=1935(4)1975{
		clear
		import excel "$hauptpfad/6_Bundesblatt/Bundesblatt_in_FR_checked/`y'_Bundesblatt_backin.xlsx"	,  first
		display "`y'"
		gen year=`y'
		replace canton=canton[_n-1] if canton==""
		append using "$hauptpfad\1_data\temp.dta"
		save "$hauptpfad\1_data\temp.dta",  replace
		}
		
		erase "$hauptpfad\1_data\temp.dta"
		
		tab year
		
		gen order_initial=_n
		
		
		
		* (ii) Correct names of canton
		
		preserve
		clear
		import excel "$hauptpfad/1_data/Kantonsnamen_harmonisiert.xlsx"	,  first // Manuelle Vereinheitlichung der unterschiedlichen Kantonsnamen
		saveold "$hauptpfad\1_data\Kantonsnamen_harmonisiert.dta",  replace
		restore
		
		merge m:1 canton using "$hauptpfad\1_data\Kantonsnamen_harmonisiert.dta", gen(merge_kantonsname)
		erase "$hauptpfad\1_data\Kantonsnamen_harmonisiert.dta"
		
		drop canton 
		gen canton=kantonsname

		merge m:1 canton using "$path_original_data/cantonal_translator",gen(merge_cantonaltranslator)
		drop if canton=="JU"
		
		
		* (iii) Add municipality numbers and correct municipality names
		
				
		/*preserve						// Outfile for Simon Luechinger for correction of municipality names
		egen city_num=group(city)
		collapse  (first) city, by(city_num canton)	
		export excel "$hauptpfad\1_data\Municipalities All.xlsx", firstrow(variables) replace  
		restore
		*/ 
		
		gen gemeindename=city
		
		
		do "$hauptpfad/5_codes/01_Correct Municipalities"  // Correct municipality names  (written by Simon)
		rename city gemeindename_original
	
	
	
		preserve
		clear
		import excel "$hauptpfad/1_data/Gemeindebestand_1969.xls"	,  first sheet("exportList")  cellrange(A12:G3090) // Gemeindestand 1969
		rename Gemeindename gemeindename
		rename BFSGdenummer gemeindenummer_bfs
		saveold "$hauptpfad\1_data\Gemeindebestand_1969.dta",  replace version(12)
		restore
		
		merge m:1 gemeindename using "$hauptpfad\1_data\Gemeindebestand_1969.dta",gen(merge_gemeindebestand)
		drop if merge_gemeindebestand==2
		// 181 merge=2 all in 1975
		
		erase "$hauptpfad\1_data\Gemeindebestand_1969.dta"
		
		
		
		* (iv) Correct party lists
		

		/*preserve						// Outfile for Florence Stempfel
		egen party_num=group(party)
		collapse  (first) party, by(party_num canton)
		//export excel "$hauptpfad\1_data\Parteinamen Alle.xlsx", firstrow(variables) replace
		restore		
		*/
		
		egen party_num=group(party)
		
		preserve    // Read in party file harmonized by Florence and Lukas
		import excel "$hauptpfad\1_data\Parteinamen All Florence.xlsx", firstrow clear 
		duplicates drop canton Partei, force // drop two duplicate observations 
		rename Partei party
		saveold "$hauptpfad\1_data\Parteinamen All Florence.dta", version(12) replace
		import excel "$hauptpfad\1_data\LeftRightCoding Florence.xlsx", firstrow clear sheet("all parties recoded") 
		rename Partei party
		saveold "$hauptpfad\1_data\LeftRightCoding Florence.dta", version(12) replace
		restore
		
		merge m:1 canton party using "$hauptpfad\1_data\Parteinamen All Florence.dta", gen(merge_party_all)
		erase "$hauptpfad\1_data\Parteinamen All Florence.dta"

		// merge: 1,761 not merged from using data
		//		     32 from AI, NW, OW, UR 
		//    	  1,729 from 1975 
		
		
		split Code, p(",")
		destring Code1 Code2 Code3 Code4, replace
		rename Code Code_original
		gen varnew=_n
		reshape long Code , i(varnew) j(Code_Nr)
		drop Code_Nr
		rename party party_all
		merge m:1 Code using "$hauptpfad\1_data\LeftRightCoding Florence.dta", gen(merge_qualcheck)
		erase "$hauptpfad\1_data\LeftRightCoding Florence.dta"		
		drop if merge_qualcheck==2 // 61 merge=2 -> probably new parties (Grüne, ...)
		sort varnew
		collapse (first) no-Code_original (mean) LinksRechts, by(varnew)
		
		
		
		
		* (v) Change names to 1975-2015 format
		
		rename firstname vorname
		rename name nachname
		rename birth geburtsjahr
		rename votes kandidst_tot
		
		gen gewaehlt=.
		replace gewaehlt=1 if elected=="true" | elected=="1"
		replace gewaehlt=0  if elected=="false"| elected=="0"
		drop elected
		rename gewaehlt elected
		
		
		//save "$hauptpfad\1_data\nr_data_1931_1975_jan2017.dta", replace  // out file as a backup for the version in jan 2017

		
		* (vi) Corrections and additional information 
		
		replace elected=0 if year==1967 & canton=="VS" & party_all=="Liste 1. Liste socialiste" // four candidates misclassified as elected
		replace geburtsjahr=1906 if year==1947 & nachname=="Michlig" & vorname=="Meinard" & kandidst_tot==6468 // this information is according to Bundesblatt 1951

		
		
		
		* (vii) Quality checks
		
		
		// (a) Check number cantons per year
		
		preserve
		gen ones=1
		collapse (mean) ones , by(year canton)
		tab year
		restore
		
		// result: 25 cantons in all years
		
		tab canton if LinksRechts==. & year!=1975
		
		// result: LinksRechts is only missing for AI, NW, OW, UR
		
		tab year
		
		
		
		
		// (b) Check number candidates per year

		preserve
		drop if canton=="UR" | canton=="AI" | canton=="NW" |canton=="OW"
		drop if canton=="GL" & year>=1971
		bysort year: gen nobs=_N
		bysort year: gen indi=_n
		br nobs year if indi==1		
		restore
		
		// result: no. candidates equal to bfs data except for 1967 in which 1 candidate is missing in our data
		
		br nachname	vorname	geburtsjahr canton if  (canton=="SZ" | canton=="AR" |canton=="GL") & year==1967
		
		bysort year: egen elected_sum=sum(elected)
		bysort year: sum elected_sum
		
		/*
		gen ones=1
		keep if elected==1
		collapse (sum) ones , by(year canton)
		br year canton ones if year==1967
		*/
		
		
		br nachname vorname geburtsjahr canton kandidst_tot if year==1967 & canton=="VS" & elected==1
		sort party_all kandidst_tot 
		br nachname vorname geburtsjahr canton kandidst_tot elected party_all if year==1967 & canton=="VS" 

		// result: 1967 are 4 candidates more than on BfS list -> explanation: in one party (list 1) in VS all candidates were classified as elected
		
		preserve
		gen ones=1
		collapse (sum) ones (mean) knr , by(year canton party_all elected)
		sort year knr party_all elected
		br year knr party_all elected if year==1967
		
		gen party_all_num = regexs(2) if regexm(party_all, "^([^0-9]*)([0-9]+)([^0-9]*)$")
		destring party_all_num, replace
		sort year knr party_all_num elected 
		
		br   year knr party_all elected ones if year==1967 
		restore
		
		/* result: All differences to the BfS list (Anzahl Kandidierende bei den Nationalratswahlen 1928 – 2015)
				   can be explained by "stille Wahlen" and Majorzkantone except 1967. 
				   In 1967 we have one candidate less than on the BfS list. We checked the number of 
				   elected and not elected candidates on all lists and compared it to the original 
				   Bundesblatt file. We found no difference. 
		
		(first) |
		   year |      Freq.     Percent        Cum.
	------------+-----------------------------------
		   1931 |        777        5.59        5.59
		   1935 |        973        7.00       12.59
		   1939 |        783        5.63       18.22
		   1943 |      1,021        7.35       25.57
		   1947 |        998        7.18       32.75
		   1951 |      1,070        7.70       40.45
		   1955 |      1,089        7.83       48.28
		   1959 |      1,077        7.75       56.03
		   1963 |      1,200        8.63       64.66
		   1967 |      1,262        9.08       73.74
		   1971 |      1,696       12.20       85.94
		   1975 |      1,954       14.06      100.00
	------------+-----------------------------------
		  Total |     13,900      100.00
		  
		  */
		
		
			// (c) Check number of elected candidates per canton and year
	
			preserve
			gen ones=1
			keep if elected==1
			collapse (sum) ones (mean) knr , by(year canton )
			sort canton year
			
			merge 1:1 year canton using "$hauptpfad\1_data\quality checks\seats.dta", gen(merge_seats) 
			
			gen seats_diff=seats-ones
			tab seats_diff
			
			br if seats_diff!=0
			restore
			
			sort party_all
			br if canton=="FR" & elected==1 & year==1955
			
			// result: FR 1955: Pierre Pasquier misclassified as "elected" (corrected in original file)
			
			
			
			
			
			// (d) Check changes of lists
			
			preserve
			sort order_initial
			gen party_all_num = regexs(2) if regexm(party_all, "^([^0-9]*)([0-9]+)([^0-9]*)$")
			destring party_all_num, replace
			tsset order_initial
			
			gen party_all_num_diff=d1.party_all_num
			gen knr_diff=d1.knr
			
			tab party_all_num_diff
			gen indi_partydiff=0 
			replace indi_partydiff=1 if (party_all_num_diff<0 & knr_diff==0) | (party_all_num_diff>1 & knr_diff==0)
			
			gen indi_partydiff_d=d1.indi_partydiff
			
			br nachname	vorname	geburtsjahr canton year party_all if indi_partydiff==1
			restore

			// result: all checked, one party name changed in original lists and in "Parteinamen all Florence"
			
			
			
			
			// (e) Check decreasing vote numbers within a list
			
			preserve
			sort year canton party_all kandidst_tot 
			bysort year canton party_all: gen votes_order=_n
				
			tsset order_initial
			gen votes_order_diff=d1.votes_order
			
			gen party_all_num = regexs(2) if regexm(party_all, "^([^0-9]*)([0-9]+)([^0-9]*)$")
			destring party_all_num, replace
			gen party_all_num_diff=d1.party_all_num
			gen knr_diff=d1.knr
			

			tab votes_order_diff
	
			gen indi_votesdiff=0 
			replace indi_votesdiff=1 if votes_order_diff!=-1 & knr_diff==0 & party_all_num_diff==0
			br nachname	vorname	geburtsjahr canton year party_all votes_order kandidst_tot party_all_num_diff party_all_num knr_diff votes_order_diff if indi_votesdiff==1
			restore

			
			// (f) Check missings in ocd number
			
			 br nachname vorname geburtsjahr canton year party_all  kandidst_tot elected if ocd==.
			 
			 // result: all checked with Bundesblatt and data are correct
			 
			 
			 
			 
			 // (g) Check stille Wahlen
			 
			 br nachname vorname geburtsjahr canton year party_all  kandidst_tot elected if kandidst_tot==. & elected==0
			 
			 // result: no one was not elected when "stille Wahlen" took place
			 
			 
			 
			 // (h) Check whether not elected candidates have more votes than elected candidates on the same list
			 
			 bysort year canton party_all : egen help1=min(kandidst_tot) if elected==1
			 bysort year canton party_all: egen kandidst_tot_min=max(help1)
			 bysort year canton party_all : egen help2=max(kandidst_tot) if elected==0
			 bysort year canton party_all : egen kandidst_tot_max=max(help2)
			 drop help1 help2
			 
	 
			 sort canton year party_all elected kandidst_tot
	 		 

			 br nachname vorname geburtsjahr  canton year party_all elected kandidst_tot kandidst_tot_min kandidst_tot_max if kandidst_tot_min<=kandidst_tot_max 
			
			 // result: all checked and plausible

	

			// (h) Check whether party changes if canton changes

			sort year canton party_all kandidst_tot 
			bysort year canton party_all: gen votes_order=_n
				
			tsset order_initial
			gen votes_order_diff=d1.votes_order
			gen elected_diff=d1.elected
	
			gen party_all_num = regexs(2) if regexm(party_all, "^([^0-9]*)([0-9]+)([^0-9]*)$")
			destring party_all_num, replace
			gen party_all_num_diff=d1.party_all_num
			gen knr_diff=d1.knr
			
			br if party_all_num_diff==0 & knr_diff!=0
			
			// result: party always changes when canton changes
			
			
			
			
			// (i) Check whether there exist cases with no change in party list but change from "not elected" to "elected"

						
			br nachname vorname geburtsjahr  canton year party_all elected elected_diff knr_diff party_all_num_diff if party_all_num_diff==0 & elected_diff==1

			// result: no such cases 
			
			
			
			// (j) Check whether first list in a canton is list number 1 (and not 2)
			

			tab party_all_num if knr_diff!=0,  missing 
			
			// result: all first lists in a canton are list 1
			
			
		
		
		
			
			
			* (vii) Merge 1971 information on gender, municipality and municipality id
			
			
			
			preserve
			cd "$path_original_data"
			import excel using "Kovariate_1971.xlsx" , clear first
			foreach x of varlist _all {
			rename `x' `x'_1971
			} 
			rename Kantonsname_1971 canton
			rename Vorname_1971 vorname
			rename Nachname_1971 nachname 
			gen year=1971
			replace vorname="Elisabeth" if nachname=="Rüegsegger-Suter" & vorname=="Elisbeth"
			replace nachname="Memmishofer" if nachname=="Hemmishofer" & vorname=="Jean"
			replace vorname="François" if nachname=="Rossé" & vorname=="Francois"
			replace nachname="Kaiser" if nachname=="Kaiser-Kraemer" & vorname=="Gisela"
			replace nachname="Tuscher" if nachname=="Tüscher" & vorname=="Hans"
			replace nachname="Jenny" if nachname=="Jenni" & vorname=="Olga"
			replace vorname="Jürg" if nachname=="Niklaus" & vorname=="JÜrg" 
			replace nachname="Hernandez-Kartaschoff" if nachname=="Hernandez-Kartaschof" & vorname=="Rosemarie"
			replace vorname="Ingrid" if nachname=="Risch" & vorname=="Ingrid Henriette"
			replace nachname="Bäder" if nachname=="Baeder" & vorname=="Peter"
			replace nachname="Güntensperger-Gsell" if nachname=="Guentensperger-Gsell" & vorname=="Sibylla"
			replace vorname="André J.-R." if nachname=="Feignoux" & vorname=="André J.-R"
			replace nachname="Dietrich-Schellenberg" if nachname=="Dietrich-Schellenber" & vorname=="Erica"
			replace vorname="Hansjörg W." if nachname=="Furrer" & vorname=="Hansjörg W"
			replace nachname="De Lorenzo" if nachname=="de Lorenzo" & vorname=="Kurt"
			replace vorname="Martha" if nachname=="Müller-Ledergerber" & vorname=="Mahrta"
			replace nachname="à Porta" if nachname=="à.Porta" & vorname=="Reto"
			replace vorname="Rudolf-Christoph" if nachname=="Schellenberg" & vorname=="Rudolf Christoph"
			replace vorname="Willy" if nachname=="Walker" & vorname=="Willi"
			
			// change municipality ids
			
			
			replace Gemeindenummer_BFS_1971=4260 if nachname=="Schmid-Bruggisser" & vorname=="Elisabeth" & Gemeindenummer_BFS_1971==4036 
			replace Gemeindenummer_BFS_1971=352 if nachname=="Vetsch" & vorname=="Hans" & Gemeindenummer_BFS_1971==351 
			replace Gemeindenummer_BFS_1971=308 if nachname=="Weber" & vorname=="Hans" & Gemeindenummer_BFS_1971==588 
			replace Gemeindenummer_BFS_1971=511 if nachname=="Gerber" & vorname=="Isaac" & Gemeindenummer_BFS_1971==6800 
			replace Gemeindenummer_BFS_1971=445 if nachname=="Beuret" & vorname=="Jean-Pierre" & Gemeindenummer_BFS_1971==435 
			replace Gemeindenummer_BFS_1971=958 if nachname=="Blumenstein" & vorname=="Jürg" & Gemeindenummer_BFS_1971==955 
			replace Gemeindenummer_BFS_1971=352 if nachname=="Blatter" & vorname=="Rolf" & Gemeindenummer_BFS_1971==351 
			replace Gemeindenummer_BFS_1971=2701 if  nachname=="Weber" & vorname=="Rudolf" & Gemeindenummer_BFS_1971==1058 & canton=="BE"
			replace Gemeindenummer_BFS_1971=3236 if  nachname=="Gautschi" & vorname=="Ernst" & Gemeindenummer_BFS_1971==4746 
			replace Gemeindenummer_BFS_1971=3236 if  nachname=="Rosenberger" & vorname=="Hans-Peter" & Gemeindenummer_BFS_1971==4746 
			replace Gemeindenummer_BFS_1971=1346 if  nachname=="Diethelm-Dobler" & vorname=="Josef" & Gemeindenummer_BFS_1971==1349
			replace Gemeindenummer_BFS_1971=4937 if  nachname=="Saner"  & vorname=="Fritz" & Gemeindenummer_BFS_1971==4631
			replace Gemeindenummer_BFS_1971=4534 if  nachname=="Möckli"  & vorname=="Gustav" & Gemeindenummer_BFS_1971==4781 
			replace Gemeindenummer_BFS_1971=4762 if  nachname=="Fritschi"  & vorname=="Rudolf" & Gemeindenummer_BFS_1971==4724 
			replace Gemeindenummer_BFS_1971=4671 if  nachname=="Hangartner"  & vorname=="Ulrich" & Gemeindenummer_BFS_1971==4641 
			replace Gemeindenummer_BFS_1971=5172 if  nachname=="Monico"  & vorname=="Nice" & Gemeindenummer_BFS_1971==5192 
			replace Gemeindenummer_BFS_1971=1201 if  nachname=="Weber"  & vorname=="Alfred" & Gemeindenummer_BFS_1971==6460 
		    replace Gemeindenummer_BFS_1971=5502 if  nachname=="Huggli"  & vorname=="Jean" & Gemeindenummer_BFS_1971==5487 
		    replace Gemeindenummer_BFS_1971=6236 if  nachname=="Rey"  & vorname=="Alfred" & Gemeindenummer_BFS_1971==6248 
			replace Gemeindenummer_BFS_1971=6217 if  nachname=="Rouiller"  & vorname=="Claude" & Gemeindenummer_BFS_1971==6616 
			replace Gemeindenummer_BFS_1971=118 if  nachname=="Schenk"  & vorname=="Fritz" & Gemeindenummer_BFS_1971==142  
			replace Gemeindenummer_BFS_1971=118 if  nachname=="Zedi"  & vorname=="Gotthard" & Gemeindenummer_BFS_1971==142 
			replace Gemeindenummer_BFS_1971=118 if  nachname=="Rüegg"  & vorname=="Hans" & Gemeindenummer_BFS_1971==142 
			replace Gemeindenummer_BFS_1971=118 if  nachname=="Hungerbühler"  & vorname=="Hugo" & Gemeindenummer_BFS_1971==142 
			replace Gemeindenummer_BFS_1971=118 if  nachname=="Lienhard"  & vorname=="Konrad" & Gemeindenummer_BFS_1971==142 
			replace Gemeindenummer_BFS_1971=118 if  nachname=="Kunz"  & vorname=="Walter" & Gemeindenummer_BFS_1971==142 
			replace Gemeindenummer_BFS_1971=118 if  nachname=="Kyburz"  & vorname=="Walter" & Gemeindenummer_BFS_1971==142 
			replace Gemeindenummer_BFS_1971=63 if  nachname=="Furrer"  & vorname=="Hansjörg W." & Gemeindenummer_BFS_1971==62 
			replace Gemeindenummer_BFS_1971=192 if  nachname=="Basler"  & vorname=="Konrad" & Gemeindenummer_BFS_1971==55 
			replace Gemeindenummer_BFS_1971=192 if  nachname=="Meier"  & vorname=="Kurt" & Gemeindenummer_BFS_1971==55 
			replace Gemeindenummer_BFS_1971=195 if  nachname=="Reich"  & vorname=="Richard" & Gemeindenummer_BFS_1971==154 

			replace Gemeinde_1971="Bolligen" if nachname=="Vetsch" & vorname=="Hans" & Gemeinde_1971=="Bern" 
			replace Gemeinde_1971="Le Bémont (BE)" if nachname=="Gerber" & vorname=="Isaac" & Gemeinde_1971=="Aux Rouges-Terre" 
			replace Gemeinde_1971="Sonvilier" if nachname=="Beuret" & vorname=="Jean-Pierre" & Gemeinde_1971=="La Chaux d'Abel"
			replace Gemeinde_1971="Bolligen" if nachname=="Blatter" & vorname=="Rolf" & Gemeinde_1971=="Bern" 
			replace Gemeinde_1971="Basel" if  nachname=="Weber" & vorname=="Rudolf" & Gemeinde_1971=="St. Niklausen" & canton=="BE"
			replace Gemeinde_1971="Schübelbach" if  nachname=="Diethelm-Dobler" & vorname=="Josef" & Gemeinde_1971=="Wangen (SZ)" 
			replace Gemeinde_1971="Unterschlatt" if  nachname=="Möckli"  & vorname=="Gustav" & Gemeinde_1971=="Wängi" 
			replace Gemeinde_1971="Kreuzlingen" if  nachname=="Hangartner"  & vorname=="Ulrich" & Gemeinde_1971=="Altnau" 
			replace Gemeinde_1971="Castagnola" if  nachname=="Monico"  & vorname=="Nice" & Gemeinde_1971=="Lugano"
			replace Gemeinde_1971="Altdorf (UR)"  if  nachname=="Weber"  & vorname=="Alfred" & Gemeinde_1971=="Altdorf" 
		    replace Gemeinde_1971="Granges (VS)" if  nachname=="Rey"  & vorname=="Alfred" & Gemeinde_1971=="Sierre" 
			replace Gemeinde_1971="Lufingen"  if  nachname=="Furrer"  & vorname=="Hansjörg W." & Gemeinde_1971=="Kloten"  
			replace Gemeinde_1971="Maur" if  nachname=="Reich"  & vorname=="Richard" & Gemeinde_1971=="Küsnacht (ZH)" 

			save Kovariate_1971.dta, replace
			restore
			
	
	
			replace nachname="Döbeli" if nachname=="Dobeli" & vorname=="Ernst" & year==1971
			replace nachname="Schwarz-Knecht" if nachname=="Schwarz" & vorname=="Heidi" & year==1971
			replace nachname="Wieser-Nielsen" if nachname=="Wieser" & vorname=="Helga" & year==1971
			replace vorname="Marie-Luise" if nachname=="Fischer" & vorname=="Marie-Louise" & year==1971
			replace nachname="Schmidt-Brugger" if nachname=="Schmidt" & vorname=="Sonja" & year==1971
			replace vorname="Claire-Lise" if nachname=="Renggli-Bonsack" & vorname=="Claire Lise" & year==1971
			replace nachname="Hofmann" if nachname=="Hofman" & vorname=="Fritz" & year==1971
			replace nachname="Kaiser" if nachname=="Kaiser-Krämer" & vorname=="Gisela" & year==1971
			replace nachname="Isler" if nachname=="Gisler" & vorname=="Heinz" & year==1971
			replace vorname="Hélène" if nachname=="Hirschi-Jeanprêtre" & vorname=="Helene" & year==1971
			replace vorname="Isaac" if nachname=="Gerber" & vorname=="Isaak" & year==1971
			replace nachname="Düby" if nachname=="Duby" & vorname=="Hans" & year==1971
			replace vorname="Hansrudolf" if nachname=="Lutz" & vorname=="Hans Rudolf" & year==1971
			replace nachname="Schletti-Stössel" if nachname=="Schletti-Stossel" & vorname=="Lucie" & year==1971
			replace nachname="Häusermann-Tièche" if nachname=="Häusermann-Tïèche" & vorname=="Lucie" & year==1971
			replace nachname="Ammann" if nachname=="Amman" & vorname=="Ulrich" & year==1971
			replace nachname="Ammann" if nachname=="Amman" & vorname=="Ulrich" & year==1967
			replace nachname="Kunz-Aeschlimann" if nachname=="Kunz-Aeschlimarm" & vorname=="Vreni" & year==1971
			replace nachname="Thalmann" if nachname=="Thalmann-Leutenegger" & vorname=="Ruth" & year==1971
			replace nachname="Faust-Kübler" if nachname=="Faust-Kubler" & vorname=="Erika" & year==1971
			replace vorname="Georges" if nachname=="Degen" & vorname=="Georg" & year==1971
			replace nachname="Benz-Gürtler" if nachname=="Benz-Gurtler" & vorname=="Jean" & year==1971
			replace nachname="Zimmerli-Silbernagel" if nachname=="Silbernagel-Zimmerli" & vorname=="Marty" & year==1971
			replace nachname="Ruffieux-Overney" if nachname=="Ruffieux Overney" & vorname=="Monique" & year==1971
			replace vorname="Hanni" if nachname=="Schwab" & vorname=="Hanny" & year==1971
			replace vorname="Ingrid" if nachname=="Risch" & vorname=="I. Henriette" & year==1971
			replace nachname="Widmer-Schmid" if nachname=="Widmer" & vorname=="Ursula" & year==1971
			replace vorname="Ernst" if nachname=="Müller" & vorname=="Ernst, Jun." & year==1971
			replace nachname="Schlatter" if nachname=="Schlauer" & vorname=="Gaspard" & year==1971
			replace nachname="Klauser-Schwab" if nachname=="Klauser Schwab" & vorname=="Irma" & year==1971
			replace vorname="Jakob Jun." if nachname=="Bütler" & vorname=="Jakob, Jun." & year==1971
			replace nachname="Diethelm-Dobler" if nachname=="Diethelm" & vorname=="Josef" & year==1959
			replace nachname="Diethelm-Dobler" if nachname=="Diethelm" & vorname=="Josef" & year==1963
			replace nachname="Diethelm-Dobler" if nachname=="Diethelm" & vorname=="Josef" & year==1967
			replace nachname="Diethelm-Dobler" if nachname=="Diethelm" & vorname=="Josef" & year==1971
			replace nachname="Diethelm-Dobler" if nachname=="Diethelm" & vorname=="Josef" & year==1975
			replace nachname="Bommeli-Reutlinger" if nachname=="Bommeli" & vorname=="Elisabeth" & year==1971
			replace vorname="Franz-Norbert" if nachname=="Bommer" & vorname=="Franz Norbert" & year==1971
			replace vorname="Hanspeter" if nachname=="Fischer" & vorname=="Hans Peter" & year==1971
			replace nachname="Pinoja" if nachname=="Pinoia" & vorname=="Innocente" & year==1963
			replace nachname="Pinoja" if nachname=="Pinoia" & vorname=="Innocente" & year==1967
			replace nachname="Pinoja" if nachname=="Pinoia" & vorname=="Innocente" & year==1971
			replace vorname="Emile" if nachname=="Chevalier" & vorname=="Emile-Louis" & year==1971
			replace nachname="Ethenoz-Damond" if nachname=="Ethenoz Damond" & vorname=="Gabrielle" & year==1971
			replace nachname="Payot" if nachname=="Payot " & vorname=="Pierre" & year==1971
			replace vorname="Anton E." if nachname=="Schrafl" & vorname=="Anton" & year==1971
			replace vorname="Erhard junior" if nachname=="Spörri" & vorname=="Erhard, Jun." & year==1971
			replace vorname="Hans-Ulrich" if nachname=="Frei-Wohlgemuth" & vorname=="Hans Ulrich" & year==1971
			replace vorname="Heinrich" if nachname=="Müller" & vorname=="Heinrich C." & year==1971
			replace nachname="Graemiger" if nachname=="Graeminger" & vorname=="Peter" & year==1971
			replace nachname="Graemiger" if nachname=="Graeminger" & vorname=="Peter" & year==1975
			replace vorname="Rudolf-Christoph" if nachname=="Schellenberg" & vorname=="Rudolf C." & year==1971
			replace vorname="Walter Albert" if nachname=="Peter" & vorname=="Walter" & geburtsjahr==1938 & year==1971
			
			bysort year canton vorname nachname: gen indimax=_N
			br if indimax>1
			gen identifier_1971=.
			replace identifier_1971=geburtsjahr if indimax>1 
			
			replace vorname="Hans2" if nachname=="Pfister" & vorname=="Hans" & geburtsjahr==1903 & year==1963 & kandidst_tot==21526 // duplicate Hans Pfister in 1963
			merge 1:1 year canton vorname nachname identifier_1971 using Kovariate_1971.dta ,gen(merge_covariates)	
			erase Kovariate_1971.dta
			replace vorname="Hans" if nachname=="Pfister" & vorname=="Hans2" & geburtsjahr==1903 & year==1963 & kandidst_tot==21526 // duplicate Hans Pfister in 1963
			
			// Check cases with different municipality ids (Simon's coding vs. BfS coding)
			
			br nachname vorname geburtsjahr canton gemeindenummer_bfs gemeindename Gemeinde_1971 Gemeindenummer_BFS_1971 if gemeindenummer_bfs!=Gemeindenummer_BFS_1971 & year==1971
			
			//    Notes: Different municipality ids for Lajoux (BE later JU): we did not change anything  
			//											Mervelier (BE later JU): we did not change anything
			//											Hasle bei Burgdorf or Rüegsau: we did not change anything b/c two possibilities (Hasle-Rüegsau in original document)
			//											Wuppenau: we did not change anything
		
			
			
			
			
			gen str1 geschlecht="M" if year<1971
			replace geschlecht=Geschlecht_1971 if year==1971 


			replace gemeindename=Gemeinde_1971 if year==1971 
			replace gemeindenummer_bfs=Gemeindenummer_BFS_1971 if year==1971
			

			
			* (viii) Recode missing origin information
			
			
			replace origin=gemeindename if origin=="von und" | origin=="de et"	| origin=="i ed"

			
			
			preserve // create origin list
			bysort canton origin: gen indi=_n
			keep if indi==1
			keep origin canton
			export excel "$hauptpfad\1_data\8 origin data\OriginUnique_Out.xlsx", replace firstrow(varlabels)
			restore
	
	
			* (ix) Rename variables for merge with post-1971 data
			
			keep Kandidatennummer_1971 nachname vorname geburtsjahr canton knr kandidst_tot job party_all year gemeindename gemeindenummer_bfs LinksRechts	elected geschlecht
		
		
			rename Kandidatennummer_1971 candidateno
			rename kandidst_tot votes
			rename knr cantonno
			rename party_all list
			rename gemeindenummer_bfs municipalityno
			rename gemeindename municipality
			rename LinksRechts leftright
			rename geschlecht sex
			rename geburtsjahr birthyear 
			rename vorname firstname
			rename nachname name

			
		
			drop if year==1975 // drop 1975 b/c we use BFS data for this year

	
	
			

		
			
		******************************************************
		* (C) MERGE PRE- AND POST-1975-DATA
		******************************************************

		 * (i) Append all data and labeling

			append 	using "NRW_KANDIDATEN_ALL including 2015.dta"
			
			order canton* year name firstname birthyear sex votes elected candidateno municipality municipalityno leftright list job status

			label var canton "Canton"
			label var cantonno "Canton no."
			label var year "Year"
			label var name "Name"
			label var firstname "First name"
			label var birthyear "Year of birth"
			label var sex "Sex"
			label var votes "Votes"
			label var elected "Elected"
			label var candidateno "Candidate no. on list"
			label var municipality "Municipality"
			label var municipalityno "Municipality no."
			label var leftright "Left-right coding of list"						
			label var list "List"
			label var job "Job"
			label var status "Status"

		* (ii) Quality checks

		
			* (a) Is Municipality id constant per municipality? 
			
			bysort canton municipality: egen municipality_sd=sd(municipalityno) 			
			br if municipality_sd>0 & !missing(municipality_sd) // we saved this dataset to check (Quality Check Gemeinde für Florence.xlsx)		
			drop municipality_sd
			
			
			// result: 164 observations need correction, action: run do file with our corrections
			
			do "$hauptpfad/5_codes/01_Correct Municipalities 1975-2015 part I.do" //
		
		
		
			* (b) Check missing values
			
			// total observations in dataset: 41'606
			
			count if missing(cantonno) // result: 0 non-missings, action: nothing
			count if missing(year) // result: 0 non-missings, action: nothing
			count if missing(name) // result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
			tab canton if  missing(name)
			count if missing(firstname) // result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
			tab canton if  missing(firstname)
			count if missing(birthyear) // result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
			count if missing(sex) // result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
			count if missing(votes) // result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
			sort year canton name firstname
			count if missing(votes)
			br if missing(votes) & !missing(birthyear) // result: 92 observations with no votes, all b/c "Stille Wahlen" (57 individuals in 1939, 35 individuals in the period of 1945-2015)
			count if missing(elected)
			br if missing(elected) // result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
			count if missing(candidateno)
			tab year if missing(candidateno) // result: missing before 1971, action: nothing
			tab year if !missing(candidateno)
			count if municipality=="" // result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
			
			preserve
			count if missing(municipalityno)
			count if missing(municipalityno) & municipality!=""	
			//br if missing(municipalityno) // we collapsed this subsample (by canton and municipality) to construct the outfile for Florence (Municipalities All Florence.xlsx)
			keep if missing(municipalityno)
			collapse municipalityno, by(municipality canton)
			restore
			
			do "$hauptpfad/5_codes/01_Correct Municipalities 1975-2015 part II.do" 
			
			preserve    // Read in party file harmonized by Florence and Lukas
			import excel "$hauptpfad\1_data\Gemeindebestand_20151231.xls", firstrow clear cellrange(A12:F2336)
			rename BFSGdenummer municipalityno_2015
			rename Gemeindename municipality
			keep municipalityno_2015 municipality
			saveold "$hauptpfad\1_data\Gemeindebestand_20151231.dta", version(12) replace
			restore
			
			merge m:1 municipality using "$hauptpfad\1_data\Gemeindebestand_20151231.dta" // note: many mismatches from master and using ok
			drop if _merge==2
			erase "$hauptpfad\1_data\Gemeindebestand_20151231.dta"
						
			replace municipalityno=municipalityno_2015 if missing(municipalityno)
			drop _merge municipalityno_2015
			
			count if missing(municipalityno) // result: 53 Vereinzelte have missing municipalityno
			
			count if missing(leftright) 
			tab canton if missing(leftright) // result: all 85 observations with missing data are in AI, NW, OW, UR 
			br if  missing(leftright)        // action: we do not replace the missing leftright score b/c the files for these cantons are small (<39 observations)

			
			count if list=="" 
			br if list==""  & !missing(leftright) // result: same set of individuals that also lacks leftright coding
												  // action: no change
			
			count if job=="" 
			tab year if job==""  // result: no jobs for those>1975, only about 200 observations with missing job information for year<1975
								 // action: no change 
			
			tab year if status=="" // result: no status available before 1975
			tab canton if status=="" & year>=1975 // 53 "Vereinzelte"
			br if status=="" & year>=1975

			

		******************************************************
		* (D) SAVE DATASET FOR BISNODE
		******************************************************
		
		
		* (i) Test output for Bisnode 
		
			preserve
			
			gsort cantonno name firstname birthyear municipalityno year
			
			bysort name firstname birthyear municipalityno: gen indicator=_n
			bysort name firstname birthyear municipalityno: gen indicator_max=_N
			keep if indicator==indicator_max
			drop indicator*
			
			
			// duplicates drop cantonno name firstname birthyear municipalityno, force
			
			egen rank = rank(runiform()) if !missing(municipalityno)
			
			keep if rank<1001
			
			keep canton name firstname birthyear municipality municipalityno 
			sort canton name firstname birthyear municipality municipalityno 
			
			export excel "$hauptpfad\1_data\Bisnode_test.xlsx", replace firstrow(varlabels)

			restore
			
		* (ii) Main file for Bisnode 
		
			
			preserve
			
			gsort cantonno name firstname birthyear municipalityno year
			
			bysort name firstname birthyear municipalityno: gen indicator=_n
			bysort name firstname birthyear municipalityno: gen indicator_max=_N
			keep if indicator==indicator_max
			drop indicator*
						
			//duplicates drop cantonno name firstname birthyear municipality, force
			
			
			keep canton name firstname birthyear municipality municipalityno 
			sort canton name firstname birthyear municipality municipalityno 
			
			export excel "$hauptpfad\1_data\Bisnode_all.xlsx", replace firstrow(varlabels)
			
			
			restore
	

		******************************************************
		* (E) SAVE DATASET FOR FUZZY MERGE (1931 until 2015)
		******************************************************

		cd "$path_fuzzy_data"
		save nationalraete_1931_2015_fuzzy.dta, replace
			
		******************************************************
		* (F) MERGE IDS AND ELECTION RESULTS (1931 until 2015)
		******************************************************
		
		
		* (i) Read in electoral results
		
		cd "$path_fuzzy_data"
		use nationalraete_1931_2015_fuzzy.dta, clear
		
		gen canton_merge=canton
		replace canton_merge="BEJU" if canton=="BE" | canton=="JU"	

		
		* (ii) correct birthyears and names (from fuzzy merge file)
		
		replace firstname="Patrick" if   name=="Mächler"   &	birthyear==1983 & year==2011
		replace firstname="Massimiliano" if name=="Ay" 	 &	birthyear==1982 & year==2011
		replace firstname="Sergio" if  name=="Arigoni" &		birthyear==1965 & year==2011
		replace firstname="Werner 2" if name=="Müller" &	 	birthyear==1932 & year==1975 & firstname=="Werner"
		
		rename name name_merge // undo correct changes by Jana Jarck (in 1st check round as saved in  erf_luzern_all.dta) for merge to election results
		rename firstname firstname_merge
		rename birthyear birthyear_merge 
		
	
		* (ii) Merge electoral results and fuzzy merge
		
		merge 1:1 firstname_merge name_merge canton_merge  birthyear_merge municipality year using "$datapath/FuzzyMerge_Ultimate.dta", gen(merge_fuzzy_election_results)
		
		// result: 59 from master but not in using (6 candidates that were elected in "Stillen Wahlen" (see below) and 53 "Vereinzelte"
		
		
		*** replace IDs of those candidates that were elected in "stille Wahlen" after 1975
				
		replace id_Ultimate="AR-1975-0004" if canton=="AR" &name=="Merz" & firstname=="Christian" & year==1979 // Christian Merz (stille Wahlen)
		replace id_Ultimate="AR-1975-0002" if canton=="AR" &name=="Früh" & year==1979 // Hans-Rudolf Früh (stille Wahlen)
		replace id_Ultimate="AR-1975-0002" if canton=="AR" &name=="Früh" & year==1987 // Hans-Rudolf Früh (stille Wahlen)
		replace id_Ultimate="AR-1983-0004" if canton=="AR" &name=="Maeder" & year==1987 // Herbert Mäder (stille Wahlen)
		replace id_Ultimate="NW-1995-0002" if canton=="NW" &name=="Engelberger" & year==2007 // Edi Engelberger(stille Wahlen)
		replace id_Ultimate="OW-1995-0002" if canton=="OW" &name=="Durrer" & year==1999 // Adalbert Durrer(stille Wahlen)
	
		replace firstname="Werner" if name=="Müller" &	 	birthyear==1932 & year==1975 & firstname=="Werner 2"
	
	
	
		* (iii) Incorporate corrections and unclear cases from final round by JJ and RB
	
		// correct birthyears, names and firstnames (from unsicher coding of JJ and RB in the final round)

		replace birthyear=1961 if id_Ultimate=="AG-1991-0052" & year==1995
		replace birthyear=1878 if id_Ultimate=="BEJU-1935-0165" & year==1939
		replace birthyear=1878 if id_Ultimate=="BEJU-1935-0165" & year==1943
		replace birthyear=1915 if id_Ultimate=="BEJU-1947-0030" & year==1971
		replace firstname="Paul" if id_Ultimate=="BEJU-1951-9104" & year==1959 // wrong name; according to our cantonal dataset there was never an Emil Schorer in the Berner Grossrat
		replace birthyear=1950 if id_Ultimate=="GR-1987-0001" & year==2011 
		replace birthyear=1949 if id_Ultimate=="LU-1975-7053" & year==1975
		replace birthyear=1894 if id_Ultimate=="SG-1935-9026" & year==1935
		replace birthyear=1882 if id_Ultimate=="SG-1935-9026" & year==1935
		replace birthyear=1891 if id_Ultimate=="VD-1931-9026" & year==1931
		replace birthyear=1878 if id_Ultimate=="VS-1931-7097" & year==1943
		replace birthyear=1924 if id_Ultimate=="VS-1955-0027" & year==1967
		replace name="Pfister" if id_Ultimate=="ZH-1955-7149" & year==1955 // wrong name: Siegfried Pfister
		
		
		// unclear birthyears
	
		* Ernst Zubler: 1894 or 1884, years: 1935 and 1939,   id_Ultimate_JJ=="AG-1935-9027" 
		* Paul Meier: 1906 or 1907, years:  1955 and 1963,  id_Ultimate_JJ=="BEJU-1955-9070" 
		* Hans Beyeler: 1910 or 1911, years:  1971 and 1975,  id_Ultimate_JJ=="BEJU-1971-9027" 
		* Florian Gantenbein: 1896 or 1903, years: 1943 and 1947,  id_Ultimate_JJ=="-1943-9009"
		* Jean Curdy, 1922 or 1927, years: 1967 and 1971, id_Ultimate_JJ=="VD-1967-7089"
		* Hans Rohner, 1928 or 1932, years: 1983 and 1987, id_Ultimate_JJ=="ZH-1983-9116"
		* Markus Hug, 1964 or 1969, years: 1991 and 1995, id_Ultimate_JJ=="ZH-1991-0308"
		* Jürg (Giorgio) Bösiger, 1955 or 1963, years: 1999 and 2003, id_Ultimate_JJ=="ZH-1999-9184"
		
		// unclear coding	
		* Beat Balmer: 1958, years: 1991 and 2007, id_Ultimate_RB=="BEJU-1991-7004", different profession, different party
			
		
		* (iv) Harmonization of names, firstnames and birthyears

		// rules: 1. If more than two observations: follow majority spelling and/or correct spelling in Bundesblatt, cantonal data, or other external sources
		//        2. If only two observations: check spelling in Bundesblatt (or cantonal data and other external sources) and if they differ we do not change anything
		//        3. We do not harmonize Umlaute and double names
		//		  4. Change abbrevations to full names (Joh. Friedr. becomes Johann Friedrich) if possible
		//		  5. Omit "Jun./sen." in first name
		//		  6. Change first names in paranthesis to separate names (Andreas (Res) becomes Res)
		//        7. Change "v." to "von", "Von" to "von", "Van" to "van, "De" to "de"
		//		  8. No change in first name if first name is used in some years but not in others (e.g., Josef Bernhard and Bernhard)

			
		//br if toCheck==1 & substr(name,1,3)=="v. "
	
		replace name="von Waldkirch" if name=="v. Waldkirch"
		replace name="von Almen" if name=="v. Almen"
		replace name="von Grünigen" if name=="v. Grünigen"
		replace name="von Siebenthal" if name=="v. Siebenthal"
		replace name="von Greyerz" if name=="v. Greyerz"
		
		replace name="Baeschlin-Régis" if name=="Baeschlin" & year==1975 & canton=="AG"
		replace name="Baeschlin-Régis" if name=="Baeschlin" & year==1979 & canton=="AG"
		replace firstname="Jacqueline" if firstname=="Regis Jacqueline" & year==1975 & canton=="AG"
		replace firstname="Jacqueline" if firstname=="Regis Jacqueline" & year==1979 & canton=="AG"
		replace name="Rey" if name=="Key" & year==1955 & canton=="AG"
		replace name="Ehrensperger" if name=="Ehrensberger" & canton=="AG" & year==1939
		replace name="Kuhn" if name=="Kühn" & canton=="AG" & year==1947
		replace name="Brandenberger" if name=="Brandenburger" & canton=="AG" & year==1943
		replace name="Gujer" if name=="Guyer" & canton=="AG" & year==1955

		// include changes_simon here
		
		
		
		
		// include changes_lukas here
		
		
		replace firstname="Katrin" if name=="Spahr" & canton=="AG" & year==1995

		
		* (v) Keep only relevant variables
		
		drop merge_fuzzy_election_results canton_merge name_merge firstname_merge birthyear_merge
		
		rename id_Ultimate ID
		label var ID "Unique person id"

		order canton ID name firstname name birthyear
		
		saveold  "$hauptpfad\1_data\nationalraete_1931_2015.dta", version(13) replace
		
		
		* (v) Check if fixed characterstics are identical within ID
		
		egen canton_num=group(canton)	
		bysort ID: egen canton_num_sd=sd(canton_num)
		br if canton_num_sd>0 & !missing(canton_num_sd)
		tab canton if canton_num_sd>0 & !missing(canton_num_sd) & ID!=""
		// result: canton varies only in BE and JU
		
		bysort ID year: gen indi_ID_year=_n if ID!=""
		tab indi_ID_year
		// result: no ID is found twice in a year
		
		gen toCheck=0
		
		
		egen name_num=group(name)	
		bysort ID: egen name_num_sd=sd(name_num)
		br if name_num_sd>0 & !missing(name_num_sd)	
		replace toCheck=1 if name_num_sd>0 & !missing(name_num_sd)	
		// result: 2195 with non-identical names
		

		egen firstname_num=group(firstname)	
		bysort ID: egen firstname_num_sd=sd(firstname_num)
		br if firstname_num_sd>0 & !missing(firstname_num_sd)	
		replace toCheck=1 if firstname_num_sd>0 & !missing(firstname_num_sd)	
		// result: 1838 with non-identical first names
		
		// change wrong firstname
		
	    replace firstname="Walter" if firstname=="Stöckli" & canton=="UR" & year==1975 // note: change of lastname by Simon Lüchinger in his recodings
		
		
		// change other cases
	
	    replace firstname="Max E." if firstname=="Max" & canton=="AG" & year==1947 & name=="Pflüger"
	    replace firstname="Max E." if firstname=="Max" & canton=="AG" & year==1959 & name=="Pflüger"
	    replace firstname="Guido A." if firstname=="Guido" & canton=="AG" & year==1999 & name=="Zäch"
	    replace firstname="Hans R." if firstname=="Hans" & canton=="BE" & year==1967 & name=="Gaschen"
	    replace firstname="Aldo A." if firstname=="Aldo" & canton=="BE" & year==1975 & name=="Merazzi"
	    replace firstname="Aldo A." if firstname=="Aldo" & canton=="BE" & year==1991 & name=="Merazzi"
	    replace firstname="Christoph G." if firstname=="Christoph" & canton=="BE" & year==1975 & name=="Walliser"
	    replace firstname="Christoph G." if firstname=="Christoph" & canton=="BE" & year==1979 & name=="Walliser"
	    replace firstname="Christoph G." if firstname=="Christoph" & canton=="BE" & year==1983 & name=="Walliser"
	    replace firstname="Jacqueline M." if firstname=="Jacqueline" & canton=="BE" & year==1987 & name=="Gottschalk"
	    replace firstname="Johann Ulrich" if firstname=="Johann" & canton=="BE" & year==2011 & name=="Grädel"
	    replace firstname="Fritz Abraham" if firstname=="Fritz" & canton=="BE" & year==1991 & name=="Oehrli"
	    replace firstname="Rudolf H." if firstname=="Rudolf" & canton=="BE" & year==1991 & name=="Strahm"
	    replace firstname="Rudolf H." if firstname=="Rudolf" & canton=="BE" & year==1995 & name=="Strahm"
	    replace firstname="Rudolf H." if firstname=="Rudolf" & canton=="BE" & year==1999 & name=="Strahm"
	    replace firstname="Wilfried" if firstname=="Willfried" & canton=="BE" & year==1999 & name=="Gasser"
	    replace firstname="Ernst B." if firstname=="Ernst" & canton=="BE" & year==1995 & name=="Hügli"
	    replace firstname="Robert C." if firstname=="Robert" & canton=="BE" & year==1999 & name=="Meyer"
	    replace firstname="Robert C." if firstname=="Robert" & canton=="BE" & year==2007 & name=="Meyer"
	    replace firstname="Robert C." if firstname=="Robert" & canton=="BE" & year==2011 & name=="Meyer"
	    replace firstname="German E." if firstname=="German" & canton=="BE" & year==1995 & name=="Clénin"
	    replace firstname="Heinz Benjamin" if firstname=="Heinz" & canton=="BE" & year==1995 & name=="Zaugg"
	    replace firstname="Erich J." if firstname=="Erich" & canton=="BE" & year==1999 & name=="Hess"
		replace firstname="Erich J." if firstname=="Erich" & canton=="BE" & year==2003 & name=="Hess"
	    replace firstname="Erich J." if firstname=="Erich" & canton=="BE" & year==2011 & name=="Hess"
	    replace firstname="Erich J." if firstname=="Erich" & canton=="BE" & year==2015 & name=="Hess"
	    replace firstname="Bernard W." if firstname=="Bernard" & canton=="BE" & year==1999 & name=="Liebich"
	    replace firstname="Josef" if firstname=="Joseph" & canton=="BE" & year==2003 & name=="Rothenfluh"
		replace firstname="Nathalie Djamila" if firstname=="Nathalie" & canton=="BE" & year==2007 & name=="Conrad"
		replace firstname="Christian W." if firstname=="Christian" & canton=="BE" & year==2007 & name=="Hodel"
		replace firstname="Yves Sandro" if firstname=="Yves" & canton=="BE" & year==2007 & name=="Berger"
		replace firstname="Katja Ursina" if firstname=="Katja" & canton=="BE" & year==2011 & name=="Nilsen-Schenkel"
		replace firstname="Willi" if firstname=="Willy" & canton=="BL" & year==1975 & name=="Breitenstein"
		replace firstname="Peter René" if firstname=="Peter" & canton=="BL" & year==1983 & name=="Staub"
		replace firstname="Peter René" if firstname=="Peter" & canton=="BL" & year==2011 & name=="Staub"
		replace firstname="Lucas" if firstname=="Lukas" & canton=="BS" & year==1959 & name=="Bernoulli"
		replace firstname="Bruno August" if firstname=="Bruno" & canton=="BS" & year==1975 & name=="Weber"
		replace firstname="Frédéric Paul" if firstname=="Frédéric" & canton=="BS" & year==1983 & name=="Walthard"
		replace firstname="Christmuth Martin" if firstname=="Christmuth" & canton=="BS" & year==1983 
		replace firstname="Johannes Robert" if firstname=="Johannes" & canton=="BS" & year==1995 & name=="Randegger"
		replace firstname="Johannes Robert" if firstname=="Johannes R." & canton=="BS" & year==2003 &  name=="Randegger"
		replace firstname="Peter Andreas" if firstname=="Peter-Andreas" & canton=="BS" & year==1999   
		replace firstname="Peter Andreas" if firstname=="Peter A." & canton=="BS" & year==2003   
		replace firstname="Lukas R." if firstname=="Lukas" & canton=="BS" & year==2011   & name=="Michel"
		replace firstname="Lukas R." if firstname=="Lukas" & canton=="BS" & year==2015   & name=="Michel"
		replace firstname="Heinrich S." if firstname=="Heinrich" & canton=="BS" & year==2007   & name=="Ueberwasser"
		replace firstname="Andrea Elisabeth" if firstname=="Andrea" & canton=="BS" & year==2011   & name=="Knellwolf"
		replace firstname="Ralph Alexander" if firstname=="Ralph" & canton=="FR" & year==2015   & name=="Schmid"
		replace firstname="Jean Louis" if firstname=="Louis" & canton=="GE" & year==1935   & name=="Segessemann"
		replace firstname="Jean-Charles" if firstname=="Jean Charles" & canton=="GE" & year==1995   & name=="Rielle"



	    replace firstname="Yahya Hasan" if name=="Bajwa" & canton=="AG" & year==2015
	    replace firstname="Johann Friedrich" if name=="Keller" & firstname=="Joh. Friedr." & canton=="BE" & year==1943
	    replace firstname="Léopold" if name=="Christe" &  canton=="BE" & year==1951
	    replace firstname="Emanuel" if name=="Gehret" &  canton=="BE" & year==1963
	    replace firstname="Gerhart" if name=="Schürch" &  canton=="BE" & year==1959
	    replace firstname="Heinrich" if name=="Bäbler" &  canton=="BE" & year==1975
	    replace firstname="Gréty" if name=="Hoffmeyer" &  canton=="JU" & year==1987
	    replace firstname="Marie-José" if firstname=="Marie-Josée" &  canton=="BE" & year==1987
	    replace firstname="Richard Hermann" if firstname=="Richard.Herm" &  canton=="BE" & year==1991
	    replace firstname="Rudolf" if firstname=="Rudolf H." &  canton=="BE" & year==1987
	    replace firstname="Marc Frédéric" if firstname=="Marc F." &  canton=="BE" & year==1991
	    replace firstname="Marc Frédéric" if firstname=="Marc F." &  canton=="BE" & year==1995
	    replace firstname="Marc Frédéric" if firstname=="Marc F." &  canton=="BE" & year==1999
	    replace firstname="Marc Frédéric" if firstname=="Marc F." &  canton=="BE" & year==2003
	    replace firstname="Marc Frédéric" if firstname=="Marc F." &  canton=="BE" & year==2007
	    replace firstname="Susanne" if firstname=="Susanna" &  canton=="BE" & year==2015 & name=="Meierhans"
	    replace firstname="Johann-Niklaus" if firstname=="Johann Niklaus" &  canton=="BE" & year==1999 
	    replace firstname="Johann" if firstname=="Johann N." &  canton=="BE" & year==2003 
	    replace firstname="Johann" if firstname=="Johann N." &  canton=="BE" & year==2007 
	    replace firstname="Heidy" if firstname=="Heidi" &  canton=="BE" & year==2007 & name=="Wegmann"
	    replace firstname="Patrizia" if firstname=="Patrizia-Carola" &  canton=="BE" & year==2003 
		replace firstname="Maria Esther" if firstname=="Maria E." &  canton=="BE" & year==2011
		replace firstname="Georg Max" if firstname=="Georg Max (Jorgo)" &  canton=="BE" & year==2011
		replace firstname="Séraphine" if firstname=="Seraphine" &  canton=="BE" & year==2015
		replace firstname="Anna Magdalena" if firstname=="Anna-Magdalena" &  canton=="BE" & year==2015
		replace firstname="Res" if firstname=="Andreas (Res)" &  canton=="BE" & year==2011
		replace firstname="Johannes" if firstname=="Joh. " &  canton=="BL" & year==1931
		replace firstname="Josef" if firstname=="Joseph" &  canton=="BL" & year==1943 & name=="Tschopp"
		replace firstname="Eduard" if firstname=="E." &  canton=="BS" & year==1939 & name=="Strub"
		replace firstname="August" if firstname=="Aug." &  canton=="BS" & year==1939 & name=="Ursprung"
		replace firstname="Edwin" if firstname=="E." &  canton=="BS" & year==1939 & name=="Zweifel"
		replace firstname="Willi" if firstname=="Willy" &  canton=="BS" & year==1963 & name=="Grieder-Rychen"
		replace firstname="Carl" if firstname=="C." &  canton=="BS" & year==1939 & name=="Ludwig"
		replace firstname="Max" if firstname=="M." &  canton=="BS" & year==1939 & name=="Dannenberger"
		replace firstname="Erich" if firstname=="E." &  canton=="BS" & year==1939 & name=="Bolza"
		replace firstname="Alfred" if firstname=="A." &  canton=="BS" & year==1939 & name=="Würz"
		replace firstname="Peter" if firstname=="P." &  canton=="BS" & year==1939 & name=="Zschokke"
		replace firstname="Emmanuel" if firstname=="F. Emmanuel" &  canton=="BS" & year==1951 & name=="Iselin"
		replace firstname="Felix Emmanuel" if firstname=="F. Emmanuel" &  canton=="BS" & year==1955 & name=="Iselin"
		replace firstname="Christof" if firstname=="Christoph" &  canton=="BS" & year==1975 & name=="Dressler"
		replace firstname="Peter H." if firstname=="Peter.H" &  canton=="BS" & year==1991
		replace firstname="Timothée" if firstname=="Timothé" &  canton=="BS" & year==2007
		replace firstname="Marc Antoine" if firstname=="Marc-Antoine" &  canton=="FR" & year==2015
		replace firstname="Hermann Alfred" if firstname=="Hermann.Alfred" &  canton=="GE" & year==1991
		replace firstname="Otto" if firstname=="Johann (Otto)" &  canton=="GE" & year==2015
		replace firstname="Loli" if firstname=="Dolorès (Loli)" &  canton=="GE" & year==1995
		replace firstname="Loly" if firstname=="Dolores (Loly)" &  canton=="GE" & year==2011
		replace firstname="Marie" if firstname=="Marie (Maryelle)" &  canton=="GE" & year==2015
		replace firstname="Marie-Gabrielle" if firstname=="Maryelle Marie-Gabrielle" &  canton=="GE" & year==2003
		replace firstname="Saliha" if firstname=="Saliha (Salika)" &  canton=="GE" & year==2015
		replace firstname="Tobias" if firstname=="Tobia" &  canton=="GE" & year==2011
		replace firstname="Stefania" if firstname=="Stefania (Stéfanie)" &  canton=="GE" & year==2015
		replace firstname="Josef" if firstname=="Joseph" &  canton=="GR" & year==1943 & name=="Condrau"
		replace firstname="Robert C." if firstname=="Robert" &  canton=="GR" & year==1935 & name=="Ganzoni"
		replace firstname="Peider" if firstname=="Pieder" &  canton=="GR" & year==1975 & name=="Ganzoni"
		replace firstname="Brigitta Maria" if firstname=="Brigitta M." &  canton=="GR" & year==1991 & name=="Gadient"
		replace firstname="Brigitta Maria" if firstname=="Brigitta M." &  canton=="GR" & year==1995 & name=="Gadient"
		replace firstname="Brigitta Maria" if firstname=="Brigitta M." &  canton=="GR" & year==1999 & name=="Gadient"
		replace firstname="Brigitta Maria" if firstname=="Brigitta M." &  canton=="GR" & year==2003 & name=="Gadient"
		replace firstname="Brigitta Maria" if firstname=="Brigitta M." &  canton=="GR" & year==2007 & name=="Gadient"
		replace firstname="Veronika" if firstname=="Veronika (Nicky)" &  canton=="GR" & year==1991
		replace firstname="Nicky" if firstname=="Nicky (Veronika)" &  canton=="GR" & year==1995
		replace firstname="Lorenzo" if firstname=="Lorenzo (Lolo)" &  canton=="GR" & year==2003
		replace firstname="Christiana" if firstname=="Christina" &  canton=="GR" & year==2003 & name=="Flütsch"
		replace firstname="Renata" if firstname=="Renate" &  canton=="GR" & year==2011 & name=="Birrer"
		replace firstname="Jan" if firstname=="Jann" &  canton=="GR" & year==2011 & name=="Koch"
		replace firstname="Franz Josef" if firstname=="Franz Joseph" &  canton=="LU" & year==1963 & name=="Kurmann"
		replace firstname="Alfons" if firstname=="Alphons" &  canton=="LU" & year==1967 & name=="Müller"
		replace firstname="Josi J." if firstname=="Josi" &  canton=="LU" & year==1975 & name=="Meier"
		replace firstname="Josi J." if firstname=="Josi" &  canton=="LU" & year==1979 & name=="Meier"
		replace firstname="Anton Franz" if firstname=="Anton F." &  canton=="LU" & year==1983 & name=="Steffen"
		replace firstname="Helen" if firstname=="Helene" &  canton=="LU" & year==1987 & name=="Leumann-Würsch"
		replace firstname="Räto B." if firstname=="Räto" &  canton=="LU" & year==2007 & name=="Camenisch"
		replace firstname="Yasikaran" if firstname=="Yasi" &  canton=="LU" & year==2011 & name=="Manoharan"
		replace firstname="Fernand-Alfred" if firstname=="F.-Alfred" &  canton=="NE" & year==1955 & name=="Landry"
		replace firstname="Frédéric" if firstname=="Frederic" &  canton=="NE" & year==1967 & name=="Blaser"
		replace firstname="Maria Angela" if firstname=="Maria-Angela" &  canton=="NE" & year==2007 & name=="Guyot"
		replace firstname="Sonia Stella" if firstname=="Sonia" &  canton=="NE" & year==2011 & name=="Barbosa"
		replace firstname="Johann Konrad" if firstname=="Joh. Konrad" &  canton=="SG" & year==1931 & name=="Müller"
		replace firstname="Bernhard" if firstname=="Bernh." &  canton=="SG" & year==1947 & name=="Roth"
		replace firstname="Josef" if firstname=="J." &  canton=="SG" & year==1939 & name=="Scherrer"
		replace firstname="Johannes" if firstname=="Johann" &  canton=="SG" & year==1935 & name=="Duft"
		replace firstname="Josef Alfred" if firstname=="Josef" &  canton=="SG" & year==1935 & name=="Minikus"
		replace firstname="Josef Alfred" if firstname=="Josef" &  canton=="SG" & year==1943 & name=="Minikus"
		replace firstname="Josef Alfred" if firstname=="Josef" &  canton=="SG" & year==1947 & name=="Minikus"
		replace firstname="Peter" if firstname=="Peter Alfons" &  canton=="SG" & year==1939 & name=="Schwizer"
		replace firstname="Traugott" if firstname=="Traugott, Jun." &  canton=="SG" & year==1935 & name=="Walter"
		replace firstname="Thomas" if firstname=="Th." &  canton=="SG" & year==1939 & name=="Holenstein"
		replace firstname="Josef" if firstname=="J." &  canton=="SG" & year==1939 & name=="Riedener"
		replace firstname="Albert" if firstname=="Alb." &  canton=="SG" & year==1947 & name=="Spindler"
		replace firstname="Arnold" if firstname=="Arn." &  canton=="SG" & year==1939 & name=="Kappler"
		replace firstname="Gottlieb" if firstname=="G." &  canton=="SG" & year==1947 & name=="Graf"
		replace firstname="Gallus" if firstname=="G." &  canton=="SG" & year==1939 & name=="Eugster"
		replace firstname="Sebastian" if firstname=="Seb." &  canton=="SG" & year==1947 & name=="Engel"
		replace firstname="Josef" if firstname=="Jos." &  canton=="SG" & year==1947 & name=="Strässle"
		replace firstname="Hans" if firstname=="H." &  canton=="SG" & year==1947 & name=="Sturzenegger"
		replace firstname="Adolf" if firstname=="Ad." &  canton=="SG" & year==1947 & name=="Huber"
		replace firstname="Gottfried" if firstname=="Gottfr." &  canton=="SG" & year==1947 & name=="Münger"
		replace firstname="Johann" if firstname=="Joh." &  canton=="SG" & year==1947 & name=="Reich"
		replace firstname="Josef" if firstname=="J." &  canton=="SG" & year==1947 & name=="Schöbi"
		replace firstname="Wendelin" if firstname=="Wendel" &  canton=="SG" & year==2007 & name=="Rüttimann"
		replace firstname="Reto F." if firstname=="Reto" &  canton=="SG" & year==1999 & name=="Denoth"
		replace firstname="Daniel M." if firstname=="Daniel" &  canton=="SG" & year==1999 & name=="Häusermann"
		replace firstname="Carl Eugen" if firstname=="Carl E." &  canton=="SH" & year==1959 & name=="Scherrer"
		replace firstname="Willy" if firstname=="Willi" &  canton=="SO" & year==1979 & name=="Pfund"
		replace firstname="Peter J." if firstname=="Peter" &  canton=="SO" & year==1983 & name=="Aebi"
		replace firstname="Peter M." if firstname=="Peter" &  canton=="SO" & year==2015 & name=="Linz"
		replace firstname="Roland F." if firstname=="Roland" &  canton=="SO" & year==1991 & name=="Borer"
		replace firstname="Roland F." if firstname=="Roland" &  canton=="SO" & year==1995 & name=="Borer"
		replace firstname="Roland F." if firstname=="Roland" &  canton=="SO" & year==1999 & name=="Borer"
		replace firstname="Roland F." if firstname=="Roland" &  canton=="SO" & year==2003 & name=="Borer"
		replace firstname="Roland F." if firstname=="Roland" &  canton=="SO" & year==2007 & name=="Borer"
		replace firstname="Roland F." if firstname=="Roland" &  canton=="SO" & year==2015 & name=="Borer"
		replace firstname="Hansjörg" if firstname=="Hansjürg" &  canton=="SO" & year==2007 & name=="Stoll"
		replace firstname="Johannes" if firstname=="Johann" &  canton=="TG" & year==1931 & name=="Lymann"
		replace firstname="Robert" if firstname=="Rob." &  canton=="TG" & year==1939 & name=="Züllig"
		replace firstname="Alphons" if firstname=="Alphonse" &  canton=="TG" & year==1939 & name=="von Streng"
		replace firstname="Hansueli" if firstname=="Hansuli" &  canton=="TG" & year==1991 & name=="Raggenbass"
		replace firstname="Willy J." if firstname=="Willy" &  canton=="TG" & year==1995 & name=="Schmidhauser"
		replace firstname="Willy J." if firstname=="Willy" &  canton=="TG" & year==1999 & name=="Schmidhauser"
		replace firstname="Willy J." if firstname=="Willy" &  canton=="TG" & year==2003 & name=="Schmidhauser"
		replace firstname="Willy J." if firstname=="Willy" &  canton=="TG" & year==2007 & name=="Schmidhauser"
		replace firstname="Gabi" if firstname=="Gabriele (Gabi)" &  canton=="TG" & year==2007 & name=="Badertscher"
		replace firstname="Urs-Peter" if firstname=="Urs Peter" &  canton=="TG" & year==2007 & name=="Beerli"
		replace firstname="Ruggero" if firstname=="Ruggero  " &  canton=="TI" & year==1935 & name=="Dollfuss"
		replace firstname="Ruggero" if firstname=="Roggero" &  canton=="TI" & year==1931 & name=="Dollfuss"
		replace firstname="Giovanni" if firstname=="Giovanni  " &  canton=="TI" & year==1935 & name=="Polar"
		replace firstname="Riccardo" if firstname=="Ricardo" &  canton=="TI" & year==1939 & name=="Rossi"
		replace firstname="Giovanni Battista" if firstname=="Giovanni-Battista" &  canton=="TI" & year==1947 & name=="Rusca"
		replace firstname="Piero" if firstname=="Pietro  " &  canton=="TI" & year==1935 & name=="Pellegrini"
		replace firstname="Giancarlo" if firstname=="Gian Carlo" &  canton=="TI" & year==1979 & name=="Staffieri"
		replace firstname="Franco" if firstname=="Francesco (Franco)" &  canton=="TI" & year==1995 & name=="Cavalli"
		replace firstname="Francesco" if firstname=="Francesco (Franco)" &  canton=="TI" & year==1999 & name=="Cavalli"
		replace firstname="Francesco" if firstname=="Francesco (Franco)" &  canton=="TI" & year==2003 & name=="Cavalli"
		replace firstname="Gian Mario" if firstname=="Gianmario" &  canton=="TI" & year==1975 & name=="Pagani"
		replace firstname="Giovan Battista" if firstname=="Giovan-Battista" &  canton=="TI" & year==1979 & name=="Pedrazzini"
		replace firstname="Mario" if firstname=="Gianmario" &  canton=="TI" & year==1975 & name=="Grassi"
		replace firstname="Giovanni Maria" if firstname=="Giovan Maria" &  canton=="TI" & year==1991 & name=="Staffieri"
		replace firstname="Gianpiero" if firstname=="Gian Pietro" &  canton=="TI" & year==1995 & name=="Bernasconi"
		replace firstname="Pierluigi" if firstname=="Pier Luigi" &  canton=="TI" & year==1995 & name=="Zanchi"
		replace firstname="Fausto" if firstname=="Fausto (Gerri)" &  canton=="TI" & year==1999 & name=="Beretta-Piccoli"
		replace firstname="Gerri" if firstname=="Fausto (Gerry)" &  canton=="TI" & year==2007 & name=="Beretta-Piccoli"
		replace firstname="Norman" if firstname=="Norman (Vais)" &  canton=="TI" & year==2007 & name=="Gobbi"
		replace firstname="Elena" if firstname=="Elena (Bacche)" &  canton=="TI" & year==2007 & name=="Bacchetta"
		replace firstname="Henry" if firstname=="Henri" &  canton=="VD" & year==1935 & name=="Vallotton"
		replace firstname="Jean" if firstname=="Jean-Gabriel" &  canton=="VD" & year==1931 & name=="Vincent"
		replace firstname="Jules" if firstname=="Jules " &  canton=="VD" & year==1947 & name=="Berruex"
		replace firstname="Aimé" if firstname=="André" &  canton=="VD" & year==1947 & name=="Dormond"
		replace firstname="Aloïs" if firstname=="Alois" &  canton=="VD" & year==1959 & name=="Grob"
		replace firstname="Paul-Abram" if firstname=="Paul-A." &  canton=="VD" & year==1947 & name=="Meylan"
		replace firstname="Pierre-David" if firstname=="Pierre" &  canton=="VD" & year==1967 & name=="Candaux"
		replace firstname="Jeanlouis" if firstname=="Jean-Louis" &  canton=="VD" & year==1987 & name=="Cornuz"
		replace firstname="Charles Frédéric" if firstname=="Charles" &  canton=="VD" & year==1987 & name=="Imfeld"
		replace firstname="Alain-Valéry" if firstname=="Alain-Valery" &  canton=="VD" & year==1999 & name=="Poitry"
		replace firstname="Pierre Marcel" if firstname=="Pierre" &  canton=="VD" & year==1991 & name=="Würsch"
		replace firstname="Pierre Marcel" if firstname=="Pierre" &  canton=="VD" & year==1999 & name=="Würsch"
		replace firstname="Pierre Marcel" if firstname=="Pierre" &  canton=="VD" & year==2015 & name=="Wuersch"
		replace firstname="Robert" if firstname=="Robert (dit Ted Robert)" &  canton=="VD" & year==2003 & name=="Gurtner"
		replace firstname="Suzi" if firstname=="Susi" &  canton=="VD" & year==1999 & name=="Dulex"
		replace firstname="Marguerite" if firstname=="Margarida" &  canton=="VD" & year==2007 & name=="Vernay"
		replace firstname="André Francis" if firstname=="André" &  canton=="VD" & year==2003 & name=="Cattin"
		replace firstname="Jean-Baptiste" if firstname=="Jean-Baptise" &  canton=="VD" & year==2007 & name=="Blanc"
		replace firstname="Naima" if firstname=="Naime" &  canton=="VD" & year==2007 & name=="Topkiran"
		replace firstname="Pablo Gabriel" if firstname=="Pablo" &  canton=="VD" & year==2011 & name=="Gutierrez"
		replace firstname="Oskar" if firstname=="Oscar" &  canton=="VS" & year==1943 & name=="Schnyder"
		replace firstname="Meinrad" if firstname=="Meinard" &  canton=="VS" & year==1947 & name=="Michlig"
		replace firstname="Leo" if firstname=="Léo" &  canton=="VS" & year==1959 & name=="Stoffel"
		replace firstname="Alphons" if firstname=="Alfons" &  canton=="VS" & year==1971 & name=="Imhasly"
		replace firstname="Kevin" if firstname=="Kévin" &  canton=="VS" & year==2011 & name=="Follonier"
		replace firstname="Kevin" if firstname=="Kévin" &  canton=="VS" & year==2011 & name=="Pellouchoud"
		replace firstname="Andreas C." if firstname=="Andreas" &  canton=="ZG" & year==1971 & name=="Brunner"
		replace firstname="Andreas C." if firstname=="Andreas" &  canton=="ZG" & year==1975 & name=="Brunner"
		replace firstname="David Hirsch" if firstname=="David" &  canton=="ZH" & year==1935 & name=="Farbstein"
		replace firstname="Eduard Jakob" if firstname=="Eduard" &  canton=="ZH" & year==1935 & name=="Geilinger"
		replace firstname="Eduard Jakob" if firstname=="Eduard" &  canton=="ZH" & year==1939 & name=="Geilinger"
		replace firstname="Emil Johann" if firstname=="Emil J." &  canton=="ZH" & year==1935 & name=="Graf"
		replace firstname="Emil Johann" if firstname=="Emil J." &  canton=="ZH" & year==1939 & name=="Graf"
		replace firstname="Emil Johann" if firstname=="Emil J." &  canton=="ZH" & year==1943 & name=="Graf"
		replace firstname="Ernst Walter" if firstname=="Ernst" &  canton=="ZH" & year==1931 & name=="Högger"
		replace firstname="Ernst Walter" if firstname=="Ernst" &  canton=="ZH" & year==1947& name=="Högger"
		replace firstname="Ernst Walter" if firstname=="Ernst" &  canton=="ZH" & year==1951& name=="Högger"
		replace firstname="Ernst Walter" if firstname=="Ernst" &  canton=="ZH" & year==1955& name=="Högger"
		replace firstname="Ludwig Max" if firstname=="Ludwig" &  canton=="ZH" & year==1935& name=="Schneller"
		replace firstname="Emil J." if firstname=="Emil" &  canton=="ZH" & year==1931& name=="Walter"
		replace firstname="Jean Martin" if firstname=="Jean" &  canton=="ZH" & year==1935& name=="Zahner"
		replace firstname="Jean Martin" if firstname=="Jean" &  canton=="ZH" & year==1939& name=="Zahner"
		replace firstname="Jean Martin" if firstname=="Jean" &  canton=="ZH" & year==1943& name=="Zahner"
		replace firstname="Emil Theodor" if firstname=="Emil" &  canton=="ZH" & year==1943& name=="Albert"
		replace firstname="Jakob Fridolin" if firstname=="Jakob" &  canton=="ZH" & year==1939& name=="Büsser"
		replace firstname="Emil Otto" if firstname=="Emil" &  canton=="ZH" & year==1943& name=="Duft"
		replace firstname="Emil Otto" if firstname=="Emil" &  canton=="ZH" & year==1947& name=="Duft"
		replace firstname="Emil Otto" if firstname=="Emil" &  canton=="ZH" & year==1951& name=="Duft"
		replace firstname="Emil Otto" if firstname=="Emil" &  canton=="ZH" & year==1955& name=="Duft"
		replace firstname="Emil Otto" if firstname=="Emil" &  canton=="ZH" & year==1959& name=="Duft"
		replace firstname="Emil Otto" if firstname=="Emil" &  canton=="ZH" & year==1963& name=="Duft"
		replace firstname="Wilhelm" if firstname=="J. Wilhelm" &  canton=="ZH" & year==1939& name=="Dürsteler"
		replace firstname="Franz Mathäus" if firstname=="Franz" &  canton=="ZH" & year==1943& name=="Egger"
		replace firstname="Franz Mathäus" if firstname=="Franz" &  canton=="ZH" & year==1955& name=="Egger"
		replace firstname="Karl Johann" if firstname=="Karl" &  canton=="ZH" & year==1935& name=="Grob"
		replace firstname="Jules Frédéric" if firstname=="Jules" &  canton=="ZH" & year==1951& name=="Humbert-Droz"
		replace firstname="Jules Frédéric" if firstname=="Jules" &  canton=="ZH" & year==1955& name=="Humbert-Droz"
		replace firstname="Johann" if firstname=="Johannes" &  canton=="ZH" & year==1935& name=="Lienhard"
		replace firstname="Adolf" if firstname=="Adolf, Jun." &  canton=="ZH" & year==1943& name=="Morf"
		replace firstname="Otto" if firstname=="Otto, Jun." &  canton=="ZH" & year==1943& name=="Peter"
		replace firstname="Gustav Alois" if firstname=="Gustav" &  canton=="ZH" & year==1939& name=="Schwartz"
		replace firstname="Emil Walter" if firstname=="Emil" &  canton=="ZH" & year==1939& name=="Stocker"
		replace firstname="Hans Anton" if firstname=="Hans" &  canton=="ZH" & year==1935& name=="Pestalozzi"
		replace firstname="Xaver Alois" if firstname=="Xaver" &  canton=="ZH" & year==1939& name=="Arnet"
		replace firstname="Xaver Alois" if firstname=="Xaver" &  canton=="ZH" & year==1943& name=="Arnet"
		replace firstname="Xaver Alois" if firstname=="Xaver" &  canton=="ZH" & year==1947& name=="Arnet"
		replace firstname="Walter Albert" if firstname=="Walter" &  canton=="ZH" & year==1943& name=="Baechi"
		replace firstname="Walter Albert" if firstname=="Walter" &  canton=="ZH" & year==1947& name=="Baechi"
		replace firstname="Walter Albert" if firstname=="Walter" &  canton=="ZH" & year==1951& name=="Baechi"
		replace firstname="Walter Albert" if firstname=="Walter" &  canton=="ZH" & year==1963& name=="Baechi"
		replace firstname="Josef" if firstname=="Joseph" &  canton=="ZH" & year==1947& name=="Spichtig"
		replace firstname="Ernst Alfred" if firstname=="Ernst" &  canton=="ZH" & year==1943& name=="Stiefel"
		replace firstname="Walter" if firstname=="Walther" &  canton=="ZH" & year==1951& name=="Trüb"
		replace firstname="Walter" if firstname=="Walther" &  canton=="ZH" & year==1955& name=="Trüb"
		replace firstname="Rudolf" if firstname=="Rudolf A.M." &  canton=="ZH" & year==1959& name=="Huber"
		replace firstname="Robert Fridolin" if firstname=="Robert" &  canton=="ZH" & year==1947& name=="Schmidt"
		replace firstname="Robert Eduard" if firstname=="Robert" &  canton=="ZH" & year==1943& name=="Bühler"
		replace firstname="Robert Eduard" if firstname=="Robert" &  canton=="ZH" & year==1951& name=="Bühler"
		replace firstname="Robert Eduard" if firstname=="Robert" &  canton=="ZH" & year==1955& name=="Bühler"
		replace firstname="Robert Eduard" if firstname=="Robert" &  canton=="ZH" & year==1959& name=="Bühler"
		replace firstname="Robert Eduard" if firstname=="Robert" &  canton=="ZH" & year==1963& name=="Bühler"
		replace firstname="Robert Eduard" if firstname=="Robert" &  canton=="ZH" & year==1967& name=="Bühler"
		replace firstname="Hans Ulrich" if firstname=="Hans U." &  canton=="ZH" & year==1947& name=="Schläpfer"
		replace firstname="Hans Ulrich" if firstname=="Hans U." &  canton=="ZH" & year==1951& name=="Schlaepfer"
		replace firstname="Hans Ulrich" if firstname=="Hans U." &  canton=="ZH" & year==1955& name=="Schlaepfer"
		replace firstname="Peter Hans Jakob" if firstname=="Peter" &  canton=="ZH" & year==1955& name=="Schmidheiny"
		replace firstname="Peter Hans Jakob" if firstname=="Peter" &  canton=="ZH" & year==1959& name=="Schmidheiny"
		replace firstname="Josef Erhard" if firstname=="Josef" &  canton=="ZH" & year==1955& name=="Frey"
		replace firstname="Viktor" if firstname=="Victor" &  canton=="ZH" & year==1963& name=="Jent"
		replace firstname="Paul Anton" if firstname=="Paul" &  canton=="ZH" & year==1951& name=="Wilhelm"
		replace firstname="Adelrich Jakob" if firstname=="Adelrich" &  canton=="ZH" & year==1951& name=="Schuler"
		replace firstname="Adelrich Jakob" if firstname=="Alderich" &  canton=="ZH" & year==1955& name=="Schuler"
		replace firstname="Adelrich Jakob" if firstname=="Adelrich" &  canton=="ZH" & year==1963& name=="Schuler"
		replace firstname="Adelrich Jakob" if firstname=="Adelrich" &  canton=="ZH" & year==1971& name=="Schuler"
		replace firstname="Hans Jakob" if firstname=="Hansjakob" &  canton=="ZH" & year==1959& name=="Keller"
		replace firstname="Werner F." if firstname=="Werner" &  canton=="ZH" & year==1963& name=="Leutenegger"
		replace firstname="Emile Willy" if firstname=="Emile" &  canton=="ZH" & year==1951& name=="Naville"
		replace firstname="Emile Willy" if firstname=="Emile" &  canton=="ZH" & year==1959& name=="Naville"
		replace firstname="Emile Willy" if firstname=="Emil" &  canton=="ZH" & year==1963& name=="Naville"
		replace firstname="Werner A." if firstname=="Werner" &  canton=="ZH" & year==1955& name=="Stahel"
		replace firstname="Werner A." if firstname=="Werner" &  canton=="ZH" & year==1963& name=="Stahel"
		replace firstname="Rolf Arnold" if firstname=="Rolf" &  canton=="ZH" & year==1951& name=="Widmer"
		replace firstname="Rolf Arnold" if firstname=="Rolf A." &  canton=="ZH" & year==1959& name=="Widmer"
		replace firstname="Rolf Arnold" if firstname=="Rolf A." &  canton=="ZH" & year==1963& name=="Widmer"
		replace firstname="Rolf Arnold" if firstname=="Rolf A." &  canton=="ZH" & year==1967& name=="Widmer"
		replace firstname="Hans Edouard" if firstname=="Hans" &  canton=="ZH" & year==1955& name=="von Fischer"
		replace firstname="Hans Rudolf" if firstname=="Hans" &  canton=="ZH" & year==1951& name=="Schinz"
		replace firstname="Hans Rudolf" if firstname=="Hans R." &  canton=="ZH" & year==1955& name=="Schinz"
		replace firstname="Cyril" if firstname=="Cyrill" &  canton=="ZH" & year==1963& name=="Hegnauer"
		replace firstname="Willy" if firstname=="Willi" &  canton=="ZH" & year==1971& name=="Hochuli"
		replace firstname="Armin J." if firstname=="Armin Joh." &  canton=="ZH" & year==1955& name=="Huber"
		replace firstname="Emil Eduard" if firstname=="Emil" &  canton=="ZH" & year==1959& name=="Meier"
		replace firstname="Felix Josef" if firstname=="Felix" &  canton=="ZH" & year==1959& name=="Stoffel"
		replace firstname="Hans Ulrich" if firstname=="Hansulrich" &  canton=="ZH" & year==1967& name=="Fröhlich"
		replace firstname="Hans Ulrich" if firstname=="Hans U." &  canton=="ZH" & year==1971& name=="Fröhlich"
		replace firstname="Johannes Rudolf" if firstname=="Johannes" &  canton=="ZH" & year==1963& name=="Pfenninger"
		replace firstname="Franz" if firstname=="Francisco" &  canton=="ZH" & year==1963& name=="Raurisch"
		replace firstname="Jacques" if firstname=="Jaques" &  canton=="ZH" & year==1959& name=="Ruedin"
		replace firstname="Walter Paul" if firstname=="Walter" &  canton=="ZH" & year==1959& name=="Siegmann"
		replace firstname="Anatole" if firstname=="Anatol" &  canton=="ZH" & year==1971& name=="Toedtli"
		replace firstname="Rolf A." if firstname=="Rolf" &  canton=="ZH" & year==1967& name=="Balsiger"
		replace firstname="Rolf A." if firstname=="Rolf" &  canton=="ZH" & year==1971& name=="Balsiger"
		replace firstname="Rolf A." if firstname=="Rolf" &  canton=="ZH" & year==1975& name=="Balsiger"
		replace firstname="Rolf A." if firstname=="Rolf" &  canton=="ZH" & year==1979& name=="Balsiger"
		replace firstname="Willy" if firstname=="Willi" &  canton=="ZH" & year==1967& name=="Kaufmann"
		replace firstname="Hans Paul" if firstname=="Hans" &  canton=="ZH" & year==1971 & name=="Künzi"
		replace firstname="Hans Paul" if firstname=="Hans" &  canton=="ZH" & year==1975 & name=="Künzi"
		replace firstname="Hans Paul" if firstname=="Hans" &  canton=="ZH" & year==1979 & name=="Künzi"
		replace firstname="Hans Paul" if firstname=="Hans" &  canton=="ZH" & year==1983 & name=="Künzi"
		replace firstname="Urs Max" if firstname=="Urs" &  canton=="ZH" & year==1963 & name=="Lenzlinger"
		replace firstname="Walter P." if firstname=="Walter" &  canton=="ZH" & year==1963 & name=="Moser"
		replace firstname="Hans Christian" if firstname=="Hans" &  canton=="ZH" & year==1963 & name=="Oester"
		replace firstname="Hans Christian" if firstname=="Hans" &  canton=="ZH" & year==1971 & name=="Oester"
		replace firstname="Hans Christian" if firstname=="Hans" &  canton=="ZH" & year==1975 & name=="Oester"
		replace firstname="Hans Christian" if firstname=="Hans" &  canton=="ZH" & year==1979 & name=="Oester"
		replace firstname="Hans Christian" if firstname=="Hans" &  canton=="ZH" & year==1983 & name=="Oester"
		replace firstname="Hans Christian" if firstname=="Hans" &  canton=="ZH" & year==1987 & name=="Oester"
		replace firstname="Josef" if firstname=="Joseph" &  canton=="ZH" & year==1983 & name=="Landolt"
		replace firstname="Jakob Eugen" if firstname=="Jakob" &  canton=="ZH" & year==1971 & name=="Jaggi"
		replace firstname="Alois Josef" if firstname=="Alois" &  canton=="ZH" & year==1975 & name=="Kistler"
		replace firstname="Josef" if firstname=="Joseph" &  canton=="ZH" & year==1967 & name=="Schmid"
		replace firstname="James Eduard" if firstname=="James" &  canton=="ZH" & year==1971 & name=="Schwarzenbach"
		replace firstname="James Eduard" if firstname=="James" &  canton=="ZH" & year==1975 & name=="Schwarzenbach"
		replace firstname="Hans" if firstname=="Hans Rudolf" &  canton=="ZH" & year==1979 & name=="Weidmann"
		replace firstname="Hans Georg" if firstname=="Hans" &  canton=="ZH" & year==1975 & name=="Ramseier"
		replace firstname="Hans Georg" if firstname=="Hans" &  canton=="ZH" & year==1979 & name=="Ramseier"
		replace firstname="Eduard E." if firstname=="Eduard" &  canton=="ZH" & year==1975 & name=="Perret"
		replace firstname="Eduard E." if firstname=="Eduard" &  canton=="ZH" & year==1987 & name=="Perret"
		replace firstname="Josef H." if firstname=="Josef" &  canton=="ZH" & year==1975 & name=="Ammann"
		replace firstname="Bruno A." if firstname=="Bruno" &  canton=="ZH" & year==1971 & name=="Gloor"
		replace firstname="Bruno A." if firstname=="Bruno" &  canton=="ZH" & year==1987 & name=="Gloor"
		replace firstname="Marie-Theres" if firstname=="Maria-Theres" &  canton=="ZH" & year==1971 & name=="Larcher-Schelbert"
		replace firstname="Vera M.E." if firstname=="Vera" &  canton=="ZH" & year==1971 & name=="Obeid-Ruggli"
		replace firstname="Erhard" if firstname=="Erhard junior" &  canton=="ZH" & year==1971 & name=="Spörri"
		replace firstname="Georg" if firstname=="George" &  canton=="ZH" & year==1987 & name=="Ganz"
		replace firstname="Fredi" if firstname=="Fredy" &  canton=="ZH" & year==1987 & name=="Rüegg"
		replace firstname="Jürg" if firstname=="Georg (Jürg)" &  canton=="ZH" & year==2015 & name=="Schmid"
		replace firstname="Roy Alfred" if firstname=="Roy" &  canton=="ZH" & year==1987 & name=="Kunz"
		replace firstname="Roy Alfred" if firstname=="Roy M.A." &  canton=="ZH" & year==1991 & name=="Kunz"
		replace firstname="Wally Laura" if firstname=="Wally" &  canton=="ZH" & year==1975 & name=="Widmer"
		replace firstname="Werner" if firstname=="Ernst Werner" &  canton=="ZH" & year==1987 & name=="Külling"
		replace firstname="Johannes J." if firstname=="Johannes" &  canton=="ZH" & year==1983 & name=="Müller"
		replace firstname="Ursula Charlotte" if firstname=="Ursula" &  canton=="ZH" & year==1979 & name=="Gross"
		replace firstname="Ursula Charlotte" if firstname=="Ursula" &  canton=="ZH" & year==2003 & name=="Gross Leemann"
		replace firstname="Ursula Charlotte" if firstname=="Ursula" &  canton=="ZH" & year==2015 & name=="Gross Leemann"
		replace firstname="Gérald" if firstname=="Gerald" &  canton=="ZH" & year==1979 & name=="Werner"
		replace firstname="Nelly Flora" if firstname=="Nelly" &  canton=="ZH" & year==1983 & name=="Bucher"
		replace firstname="Nelly Flora" if firstname=="Nelly" &  canton=="ZH" & year==1991 & name=="Bucher"
		replace firstname="Bernhard Andreas" if firstname=="Bernhard" &  canton=="ZH" & year==1987 & name=="Gubler"
		replace firstname="Susanne" if firstname=="Suzanne" &  canton=="ZH" & year==1987 & name=="Lechleiter-Schreiber"
		replace firstname="Hans A." if firstname=="Hans" &  canton=="ZH" & year==1987 & name=="Muther"
		replace firstname="Anjuska" if firstname=="Anjuška" &  canton=="ZH" & year==2015 & name=="Weil"
		replace firstname="Hans-Ulrich" if firstname=="Hansulrich" &  canton=="ZH" & year==1983 & name=="Lehmann"
		replace firstname="Marie-Louise" if firstname=="Marie-Louise (Malu)" &  canton=="ZH" & year==2007 & name=="Dubach"
		replace firstname="Martha" if firstname=="Marta" &  canton=="ZH" & year==1995 & name=="Günthart"
		replace firstname="Max R." if firstname=="Max" &  canton=="ZH" & year==2007 & name=="Homberger"
		replace firstname="Fritz Andreas" if firstname=="Fritz" &  canton=="ZH" & year==1987 & name=="Jäckli"
		replace firstname="Beat W." if firstname=="Beat" &  canton=="ZH" & year==1987 & name=="Müller" & birthyear==1950
		replace firstname="Beat W." if firstname=="Beat W" &  canton=="ZH" & year==1991& name=="Müller" & birthyear==1950
		replace firstname="Christof" if firstname=="Christoph" &  canton=="ZH" & year==2007& name=="Wolfer" 
		replace firstname="Bruno J." if firstname=="Bruno" &  canton=="ZH" & year==1991& name=="Zuppiger" 
		replace firstname="Bruno J." if firstname=="Bruno" &  canton=="ZH" & year==1995& name=="Zuppiger" 
		replace firstname="Bruno J." if firstname=="Bruno" &  canton=="ZH" & year==1999& name=="Zuppiger" 
		replace firstname="Bruno J." if firstname=="Bruno" &  canton=="ZH" & year==2003& name=="Zuppiger" 
		replace firstname="Bruno J." if firstname=="Bruno" &  canton=="ZH" & year==2007& name=="Zuppiger" 
		replace firstname="Bruno J." if firstname=="Bruno" &  canton=="ZH" & year==2011& name=="Zuppiger" 
		replace firstname="Bruno S." if firstname=="Bruno" &  canton=="ZH" & year==1987& name=="Müller"  & birthyear==1952
		replace firstname="Hans-Peter" if firstname=="Hanspeter" &  canton=="ZH" & year==1987& name=="Portmann" 
		replace firstname="Jean E." if firstname=="Jean E" &  canton=="ZH" & year==1991& name=="Bollier" 
		replace firstname="Stephan Urs" if firstname=="Stefan" &  canton=="ZH" & year==1995& name=="Breu" 
		replace firstname="Andreas Walter" if firstname=="Andreas" &  canton=="ZH" & year==1991& name=="Geiser" 
		replace firstname="Ueli" if firstname=="Ulrich (Ueli)" &  canton=="ZH" & year==2015& name=="Isler" 
		replace firstname="Joe A." if firstname=="Joe" &  canton=="ZH" & year==1999& name=="Manser" 
		replace firstname="Joe A." if firstname=="Joe" &  canton=="ZH" & year==2003& name=="Manser" 
		replace firstname="Joe A." if firstname=="Joe" &  canton=="ZH" & year==2011& name=="Manser" 
		replace firstname="André" if firstname=="Andi" &  canton=="ZH" & year==1991& name=="Meier" 
		replace firstname="Blanca Irene" if firstname=="Blanca" &  canton=="ZH" & year==1991& name=="Ramer" 
		replace firstname="Blanca Irene" if firstname=="Blanca" &  canton=="ZH" & year==1999& name=="Ramer" 
		replace firstname="Blanca Irene" if firstname=="Blanca" &  canton=="ZH" & year==2003& name=="Ramer" 
		replace firstname="Blanca Irene" if firstname=="Blanca" &  canton=="ZH" & year==2015& name=="Ramer" 
		replace firstname="Andreas" if firstname=="Andi" &  canton=="ZH" & year==2003& name=="Scheu" 
		replace firstname="Paul Robert" if firstname=="Paul R." &  canton=="ZH" & year==1991& name=="Graf" 
		replace firstname="Paul Robert" if firstname=="Paul" &  canton=="ZH" & year==1995& name=="Graf" 
		replace firstname="Christine Juliana" if firstname=="Christine" &  canton=="ZH" & year==1991& name=="Renner" 
		replace firstname="Hans-Jacob" if firstname=="Hans-Jakob" &  canton=="ZH" & year==2011& name=="Heitz" 
		replace firstname="Olav" if firstname=="Olaf" &  canton=="ZH" & year==1999& name=="Brunner" 
		replace firstname="Christina Barbara" if firstname=="Christina B." &  canton=="ZH" & year==1999& name=="Furrer" 
		replace firstname="Christina Barbara" if firstname=="Christina" &  canton=="ZH" & year==2003& name=="Furrer" 
		replace firstname="Christina Barbara" if firstname=="Christina" &  canton=="ZH" & year==2007& name=="Furrer" 
		replace firstname="Christina Barbara" if firstname=="Christina" &  canton=="ZH" & year==2011& name=="Furrer" 
		replace firstname="Christina Barbara" if firstname=="Christina" &  canton=="ZH" & year==2015& name=="Furrer" 
		replace firstname="Victor" if firstname=="Viktor" &  canton=="ZH" & year==2003& name=="Furrer" 
		replace firstname="Chantal Juliane" if firstname=="Chantal" &  canton=="ZH" & year==1995& name=="Galladé" 
		replace firstname="Chantal Juliane" if firstname=="Chantal" &  canton=="ZH" & year==1999& name=="Galladé" 
		replace firstname="Chantal Juliane" if firstname=="Chantal" &  canton=="ZH" & year==2003& name=="Galladé" 
		replace firstname="Chantal Juliane" if firstname=="Chantal" &  canton=="ZH" & year==2007& name=="Galladé" 
		replace firstname="Chantal Juliane" if firstname=="Chantal" &  canton=="ZH" & year==2015& name=="Galladé" 
		replace firstname="Hans-Heinrich" if firstname=="Hansheinrich" &  canton=="ZH" & year==1999& name=="Heusser" 
		replace firstname="Michèl M." if firstname=="Michèl" &  canton=="ZH" & year==1995& name=="Hurt" 
		replace firstname="Ursula" if firstname=="Ursula (Ursi)" &  canton=="ZH" & year==2007& name=="Hänni-Hauser" 
		replace firstname="Anton E." if firstname=="Anton" &  canton=="ZH" & year==1995& name=="Melliger" 
		replace firstname="Beat Andreas" if firstname=="Beat" &  canton=="ZH" & year==1999& name=="Müller" 
		replace firstname="Eliane" if firstname=="Elian" &  canton=="ZH" & year==1999& name=="Oehler" 
		replace firstname="Daniel" if firstname=="Dani" &  canton=="ZH" & year==2003& name=="Schärer" 
		replace firstname="Peter" if firstname=="Peter (Jochi)" &  canton=="ZH" & year==2007& name=="Weil-Goldstein" 
		replace firstname="Peter" if firstname=="Peter (Jochi)" &  canton=="ZH" & year==2015& name=="Weil" 
		replace firstname="Erwin Kurt" if firstname=="Erwin-Kurt" &  canton=="ZH" & year==2003& name=="Widmer" 
		replace firstname="Susanna" if firstname=="Susanne" &  canton=="ZH" & year==1999& name=="Fassnacht" 
		replace firstname="Niklaus" if firstname=="Niklaus (Nik)" &  canton=="ZH" & year==2007& name=="Gugger" 
		replace firstname="Niklaus" if firstname=="Niklaus (Nik)" &  canton=="ZH" & year==2015& name=="Gugger" 
		replace firstname="Heinz Peter" if firstname=="Heinz" &  canton=="ZH" & year==1999& name=="Kyburz" 
		replace firstname="Heinz Peter" if firstname=="Heinz" &  canton=="ZH" & year==2003& name=="Kyburz" 
		replace firstname="Heinz Peter" if firstname=="Heinz" &  canton=="ZH" & year==2007& name=="Kyburz" 
		replace firstname="Heinz Peter" if firstname=="Heinz" &  canton=="ZH" & year==2015& name=="Kyburz" 
		replace firstname="Ernst Edwin" if firstname=="Ernst" &  canton=="ZH" & year==1999& name=="Rebsamen" 
		replace firstname="Ernst Edwin" if firstname=="Ernst" &  canton=="ZH" & year==2007& name=="Rebsamen" 
		replace firstname="Gregor A." if firstname=="Gregor" &  canton=="ZH" & year==2007& name=="Rutz" 
		replace firstname="Gregor A." if firstname=="Gregor" &  canton=="ZH" & year==2015& name=="Rutz" 
		replace firstname="Erika Verena" if firstname=="Erika V." &  canton=="ZH" & year==1999& name=="Rüedi-Meier" 
		replace firstname="Bruno A." if firstname=="Bruno" &  canton=="ZH" & year==1999& name=="Sauter" 
		replace firstname="Bruno A." if firstname=="Bruno" &  canton=="ZH" & year==2007& name=="Sauter" 
		replace firstname="René" if firstname=="Rene" &  canton=="ZH" & year==1999& name=="Schwengeler" 
		replace firstname="Katharina Eva" if firstname=="Katharina" &  canton=="ZH" & year==1999& name=="Wachter-Renfer" 
		replace firstname="Katharina Eva" if firstname=="Katharina" &  canton=="ZH" & year==2003& name=="Wachter" 
		replace firstname="Katharina Eva" if firstname=="Katharina" &  canton=="ZH" & year==2011& name=="Wachter" 
		replace firstname="Hansueli" if firstname=="Hansueli (Zulu)" &  canton=="ZH" & year==2007& name=="Züllig" 
		replace firstname="Gabriela" if firstname=="Gabriela (Gabi)" &  canton=="ZH" & year==2015& name=="Bienz-Meier" 
		replace firstname="Ulrich" if firstname=="Ulrich (Ueli)" &  canton=="ZH" & year==2015& name=="Brugger" 
		replace firstname="Oskar Ulrich" if firstname=="Oskar" &  canton=="ZH" & year==2003& name=="Denzler" 
		replace firstname="Claudia Fabiana" if firstname=="Claudia" &  canton=="ZH" & year==2007& name=="Gambacciani" 
		replace firstname="Claudia Fabiana" if firstname=="Claudia" &  canton=="ZH" & year==2011& name=="Gambacciani" 
		replace firstname="Christina Marion" if firstname=="Christina Marion (Chrigi)" &  canton=="ZH" & year==2015& name=="Hug" 
		replace firstname="Christina Marion" if firstname=="Christina" &  canton=="ZH" & year==2011& name=="Hug" 
		replace firstname="Christina Marion" if firstname=="Christina (Chrigi)" &  canton=="ZH" & year==2007& name=="Hug" 
		replace firstname="Christina Marion" if firstname=="Christina" &  canton=="ZH" & year==2003& name=="Hug" 
		replace firstname="Glenda Irene" if firstname=="Glenda" &  canton=="ZH" & year==2007& name=="Loebell-Ryan" 
		replace firstname="Johann Baptist" if firstname=="Johann B." &  canton=="ZH" & year==2003& name=="Lutz" 
		replace firstname="Raphael" if firstname=="Raphael J.-P." &  canton=="ZH" & year==2015& name=="Meyer" 
		replace firstname="Elisabeth" if firstname=="Elisabeth (Lisette)" &  canton=="ZH" & year==2015& name=="Müller-Jaag" 
		replace firstname="Rudolf" if firstname=="Rudolf (Ruedi)" &  canton=="ZH" & year==2007& name=="Noser" 
		replace firstname="Lena" if firstname=="Magdalena (Lena)" &  canton=="ZH" & year==2007& name=="Schneller" 
		replace firstname="Berta" if firstname=="Bertha" &  canton=="ZH" & year==2011& name=="Stocker" 
		replace firstname="Daniel Stefan" if firstname=="Daniel" &  canton=="ZH" & year==2003& name=="Suter" 
		replace firstname="Daniel Stefan" if firstname=="Daniel" &  canton=="ZH" & year==2007& name=="Suter" 
		replace firstname="Daniel Stefan" if firstname=="Daniel" &  canton=="ZH" & year==2015& name=="Suter" 
		replace firstname="Hans Ulrich" if firstname=="Hans Ulrich (Hanf Ueli)" &  canton=="ZH" & year==2007& name=="Flückiger" 
		replace firstname="Hans Ulrich" if firstname=="Hans Ulrich (Hanf Ueli)" &  canton=="ZH" & year==2015& name=="Flückiger" 
		replace firstname="Robert" if firstname=="Robert (Röbi)" &  canton=="ZH" & year==2015& name=="Brunner" 
		replace firstname="Paul Walter" if firstname=="Paul" &  canton=="ZH" & year==2011& name=="Eggimann" 
		replace firstname="Helen" if firstname=="Helena (Helen)" &  canton=="ZH" & year==2007& name=="Freiermuth" 
		replace firstname="Helena" if firstname=="Helena (Helen)" &  canton=="ZH" & year==2015& name=="Freiermuth" 
		replace firstname="Thomas Andreas" if firstname=="Thomas" &  canton=="ZH" & year==2011& name=="Frey" 
		replace firstname="Kathrin Anna" if firstname=="Kathrin Anna (Katharina)" &  canton=="ZH" & year==2015& name=="Gander" 
		replace firstname="Bastien" if firstname=="Bastien (Bas)" &  canton=="ZH" & year==2015& name=="Girod" 
		replace firstname="Matthias Francis" if firstname=="Matthias" &  canton=="ZH" & year==2011& name=="Herfeldt" 
		replace firstname="Gabriele" if firstname=="Gabriele (Gabi)" &  canton=="ZH" & year==2015& name=="Kisker" 
		replace firstname="Mattea Julia" if firstname=="Mattea" &  canton=="ZH" & year==2007& name=="Meyer" 
		replace firstname="Mattea Julia" if firstname=="Mattea" &  canton=="ZH" & year==2015& name=="Meyer" 
		replace firstname="Mischa" if firstname=="Mischa (Mis)" &  canton=="ZH" & year==2015& name=="Müller" 
		replace firstname="Martin Thomas" if firstname=="Martin" &  canton=="ZH" & year==2007& name=="Neukom" 
		replace firstname="Martin Thomas" if firstname=="Martin" &  canton=="ZH" & year==2011& name=="Neukom" 
		replace firstname="Maria Agnes" if firstname=="Maria" &  canton=="ZH" & year==2007& name=="Rohweder-Lischer" 
		replace firstname="Maria Agnes" if firstname=="Maria" &  canton=="ZH" & year==2011& name=="Rohweder-Lischer" 
		replace firstname="Clemens" if firstname=="Clemens (Clemi)" &  canton=="ZH" & year==2015& name=="Ruckstuhl" 
		replace firstname="Marc Roland" if firstname=="Marc" &  canton=="ZH" & year==2007& name=="Rudin" 
		replace firstname="Marc Roland" if firstname=="Marc" &  canton=="ZH" & year==2015& name=="Rudin" 
		replace firstname="Horst Roland" if firstname=="Horst" &  canton=="ZH" & year==2007& name=="Zbinden" 
		replace firstname="Horst Roland" if firstname=="Horst Roland (Roli)" &  canton=="ZH" & year==2015& name=="Zbinden" 
		replace firstname="Michael Gustav" if firstname=="Michael" &  canton=="ZH" & year==2011 &name=="Zeugin" 
		replace firstname="Michael Gustav" if firstname=="Michael (Michi)" &  canton=="ZH" & year==2015 &name=="Zeugin" 
		replace firstname="Urs Michael" if firstname=="Urs" &  canton=="ZH" & year==2007 &name=="Zollinger" 
		replace firstname="Kurt Edwin" if firstname=="Kurt" &  canton=="ZH" & year==2015 &name=="Artho" 
		replace firstname="Samuel" if firstname=="Samuel (Sami)" &  canton=="ZH" & year==2015 &name=="Dubno" 
		replace firstname="Evelyn Christina" if firstname=="Evelyn" &  canton=="ZH" & year==2011 &name=="Funkhouser" 
		replace firstname="Andreas Stephan" if firstname=="Andreas" &  canton=="ZH" & year==2011 &name=="Geering" 
		replace firstname="Martin A." if firstname=="Martin" &  canton=="ZH" & year==2015 &name=="Huber" 
		replace firstname="Peter" if firstname=="Peter (Peschä)" &  canton=="ZH" & year==2015 &name=="Häni" 
		replace firstname="Peter" if firstname=="Peter (Seegras)" &  canton=="ZH" & year==2015 &name=="Keel" 
		replace firstname="Stefan" if firstname=="Stefan (Stefano)" &  canton=="ZH" & year==2015 &name=="Kunz" 
		replace firstname="Thomas Paul" if firstname=="Thomas" &  canton=="ZH" & year==2011 &name=="Lamprecht" 
		replace firstname="Marianne" if firstname=="Marianna" &  canton=="ZH" & year==2011 &name=="Landolt-Rickenbach" 
		replace firstname="Benjamin Nepomuk" if firstname=="Benjamin Samuel Nepomuk (Nepi)" &  canton=="ZH" & year==2015 &name=="Lepri" 
		replace firstname="Maurice" if firstname=="Maurice (Mo)" &  canton=="ZH" & year==2015 &name=="Maggi" 
		replace firstname="Florian Markus" if firstname=="Florian" &  canton=="ZH" & year==2011 &name=="Maier" 
		replace firstname="Samuel Pablo" if firstname=="Samuel" &  canton=="ZH" & year==2011 &name=="Müller" 
		replace firstname="Samuel Pablo" if firstname=="Samuel Pablo (Sam)" &  canton=="ZH" & year==2015 &name=="Müller" 
		replace firstname="Matthias Y." if firstname=="Matthias" &  canton=="ZH" & year==2011 &name=="Reich" 
		replace firstname="Marionna Madeleine" if firstname=="Marionna" &  canton=="ZH" & year==2015 &name=="Schlatter-Schmid" 
		replace firstname="Eva Verena" if firstname=="Eva" &  canton=="ZH" & year==2015 &name=="Steiner" 
		replace firstname="Judith Anna" if firstname=="Judith" &  canton=="ZH" & year==2011 &name=="Stofer" 
		replace firstname="Hanna" if firstname=="Hanni" &  canton=="ZH" & year==2015 &name=="Stutz" 
		replace firstname="Marc Dominic" if firstname=="Marc" &  canton=="ZH" & year==2011 &name=="Thalmann" 
		replace firstname="Robert Hans" if firstname=="Robert" &  canton=="ZH" & year==2015 &name=="Wenger" 
		replace firstname="Mark Anthony" if firstname=="Mark" &  canton=="ZH" & year==2011 &name=="Wisskirchen" 
		replace firstname="Simone Angela" if firstname=="Simone" &  canton=="ZH" & year==2015 &name=="Jundt-Wisskirchen" 
		replace firstname="Tania Kristin" if firstname=="Tania" &  canton=="ZH" & year==2015 &name=="Woodhatch" 
		
		
		// rule: change obvious abbreviations, leave cases with additional names (e.g. Alfredo (Benjamin), Claudia (Dana) etc)
		
		split firstname, p("(")
		
		br 	name	firstname	birthyear	cantonno	year if firstname2!=""
			
		replace firstname="Abenda" if firstname=="Abena (Lina)"& year==2007
		replace firstname="Alexander" if firstname=="Alexander (Aleks)"& year==2015
		replace firstname="Alexander" if firstname=="Alexander (Sascha)"& year==1995
		replace firstname="Ana Belen" if firstname=="Ana Belen (Ana)"& year==2015
		replace firstname="Annette Yasuka" if firstname=="Annette (Annette Yasuka)"& year==2015
		replace firstname="Arthur" if firstname=="Arthur (Turi)"& year==2007
		replace firstname="Benjamin" if firstname=="Benjamin (Beni)"& year==2015
		replace firstname="Charles-Emmanuel" if firstname=="Charles-Emmanuel (Charles)"& year==2015
		replace firstname="Chrisoula" if firstname=="Chrisoula (Chris)"& year==2015
		replace firstname="Christian" if firstname=="Christian (Chris)"& year==2015
		replace firstname="Claude-Marie" if firstname=="Claude (Claude-Marie)"& year==2015
		replace firstname="Cornelia" if firstname=="Cornelia (Conny)"& year==2015
		replace firstname="Cristina" if firstname=="Cristina (Criss)"& year==2015
		replace firstname="Daniel" if firstname=="Daniel (Dani)"& year==2015
		replace firstname="Dorothea" if firstname=="Dorothea (Dorothe)"& year==2015
		replace firstname="Francisco" if firstname=="Francisco (Franco)"& year==2015
		replace firstname="Gabriela" if firstname=="Gabriela (Gabi)"& year==2007
		replace firstname="Gerhard" if firstname=="Gerhard (Gery)"& year==2015
		replace firstname="Hans-Rudolf" if firstname=="Hans-Rudolf (Hans-Ruedi)"& year==2015
		replace firstname="Josef" if firstname=="Josef (Sepp) junior"& year==1995
		replace firstname="Joseph" if firstname=="Joseph (Giuseppe)"& year==1991
		replace firstname="Josianne" if firstname=="Josianne (Josi)"& year==2015
		replace firstname="Jürg" if firstname=="Jürg (Giorgio)"& year==1999
		replace firstname="Jürg" if firstname=="Jürg (Giorgio)"& year==2003
		replace firstname="Katharina" if firstname=="Katharina (Kathrin)"& year==2015
		replace firstname="Kreshnik" if firstname=="Kreshnik (Kresh)"& year==2015
		replace firstname="Kurt" if firstname=="Kurt (Kuwi)"& year==2015
		replace firstname="Ladislaus" if firstname=="Ladislaus (Lazlo)"& year==2007
		replace firstname="Markus" if firstname=="Markus (Marc)"& year==2015
		replace firstname="Melanie" if firstname=="Melanie (Melä)"& year==2015
		replace firstname="Michael" if firstname=="Michael (Michi)"& year==2015
		replace firstname="Michael" if firstname=="Michael (Mike)"& year==2015
		replace firstname="Mohan" if firstname=="Mohan (Mo)"& year==2015
		replace firstname="Monika" if firstname=="Monika (Moni)"& year==2015
		replace firstname="Nenad" if firstname=="Nenad (Neno)"& year==2015
		replace firstname="Niklaus" if firstname=="Niklaus (Niggi)"& year==2015
		replace firstname="Noe" if firstname=="Noe (Noc)"& year==2015
		replace firstname="Robert" if firstname=="Robert (Röbi)"& year==2015
		replace firstname="Simon" if firstname=="Simon (Seimon)"& year==2015
		replace firstname="Simone" if firstname=="Simone (Mona)"& year==2015
		replace firstname="Stefan" if firstname=="Stefan (Jim Bob)"& year==2011
		replace firstname="Stefan" if firstname=="Stefan (Jim Bob)"& year==2015
		replace firstname="Thomas" if firstname=="Thomas (Tom)"& year==2015
		replace firstname="Vincenzo" if firstname=="Vincenzo (Enzo)"& year==2015
		replace firstname="Walter" if firstname=="Walter (Walti)"& year==2015
		replace firstname="Walter" if firstname=="Walter (Walti)"& year==2015
		

		// check cases with only one letter in first name (subfolder abbreviation in check names.xlsx) 
		
		br firstname name canton year birthyear list job if firstname=="A."|firstname=="B."|firstname=="C."|firstname=="D."|firstname=="E."|firstname=="F."|firstname=="G."|firstname=="H."|firstname=="I."|firstname=="J."|firstname=="K."|firstname=="L."|firstname=="M."|firstname=="N."|firstname=="O."|firstname=="P."|firstname=="Q."|firstname=="R."|firstname=="S."|firstname=="T."|firstname=="U."|firstname=="V."|firstname=="W."|firstname=="X."|firstname=="Y."|firstname=="Z."
		
		replace firstname="F. Emmanuel" if firstname=="F."& name=="Iselin"
		replace firstname="Alexander" if firstname=="A."&name=="Sarasin"
		replace firstname="Benedikt" if firstname=="B."&name=="Mani"
		replace firstname="Charles" if firstname=="C."&name=="Virchaux"
		
		br firstname name canton year birthyear list job if firstname=="A. "|firstname=="B. "|firstname=="C. "|firstname=="D. "|firstname=="E. "|firstname=="F. "|firstname=="G. "|firstname=="H. "|firstname=="I. "|firstname=="J. "|firstname=="K. "|firstname=="L. "|firstname=="M. "|firstname=="N. "|firstname=="O. "|firstname=="P. "|firstname=="Q. "|firstname=="R. "|firstname=="S. "|firstname=="T. "|firstname=="U. "|firstname=="V. "|firstname=="W. "|firstname=="X. "|firstname=="Y. "|firstname=="Z. "
	
		
		// check unique firstnames manually
		
		
		preserve
		bysort firstname: gen indi=_n
		keep if indi==1
		br firstname
		restore
		
		// -> correct and add information on special cases
		

		replace firstname="Elmar Theodor" if firstname=="Elmar Th."
		replace firstname="Gebhard" if firstname=="Gebhard Jun."
		replace firstname="Hans-Jürg" if firstname=="H. J."
		replace firstname="Hans R." if firstname=="Hans R"
		replace firstname="Hans" if firstname=="Hans, Jun."
		replace firstname="Hans U." if firstname=="Hans-U."
		replace firstname="Hans-Jürg" if firstname=="Hans.Jürg"
		replace firstname="Hans Ulrich" if firstname=="Hans.Ulrich"
		replace firstname="Heinz P." if firstname=="Heinz.P"
		replace firstname="Johann Jakob" if firstname=="J. J."
		replace firstname="Jean-Jacques" if firstname=="J.-Jacques"
		replace firstname="Jakob" if firstname=="Jakob Jun."
		replace firstname="Jeanine Renée" if firstname=="Jeanine.Renée"
		replace firstname="Johann Christ." if firstname=="Joh. Christ."
		replace firstname="Johann Otto" if firstname=="Joh. Otto"
		replace firstname="Josef A." if firstname=="Jos. A."
		replace firstname="Josef" if firstname=="Josef (Sepp) junior"
		replace firstname="Josef" if firstname=="Josef, Sen."
		replace firstname="Maurice Nicolas" if firstname=="Maurice Nic."
		replace firstname="Max Johann" if firstname=="Max Joh."
		replace firstname="Otto" if firstname=="Otto, Jun."
		replace firstname="Roland W." if firstname=="Roland.W"

		

		bysort ID: egen birthyear_num_sd=sd(birthyear)
		br if birthyear_num_sd>0 & !missing(birthyear_num_sd)	
		replace toCheck=1 if birthyear_num_sd>0 & !missing(birthyear_num_sd)
		
		br  canton name firstname ID birthyear birthyear_num_sd
		// result:  with non-identical birthyears


