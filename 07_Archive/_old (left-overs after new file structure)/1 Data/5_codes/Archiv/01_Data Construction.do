
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
	
    capture cd "C:\SchmidLu\Dropbox\Projekt Nationalräte"

	if _rc==0{ 
	global hauptpfad "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
		global datapath "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
		global path_original_data "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\0_original_data"
		global path_fuzzy_data "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\fuzzy merge datasets"
		global path_fuzzy "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\3_fuzzy merge do files"
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
		
		
		
		/* merge_candinfos==1: 53 Stimmen für Vereinzelte in Kantonen mit einem NR-Sitz (Majorzkantone) 
		   merge_candinfos==2: 6 stille Wahlen
		*/

		gen knr=kantonsnummer
		merge m:1 knr using cantonal_translator,gen(merge_cantonaltranslator)
	

		
	* (iv) Merge information on birthyear
		
		preserve
		cd "$path_original_data"
		clear
		import excel  "NRW_KANDIDATEN including birthyear.xlsx",  first
		rename *, lower
		gen year=substr(nationalratswahl,4	,4)
		destring year, replace
		keep year kantonsnummer kandidatennummer listennummer geburtsjahr 
		saveold "NRW_KANDIDATEN including birthyear.dta", replace version(12)
		restore
		
		merge 1:1 year kantonsnummer kandidatennummer listennummer using "NRW_KANDIDATEN including birthyear.dta", gen(merge_birthyear)

	
		saveold "NRW_KANDIDATEN_ALL including 2015.dta", replace version(12)

		
		
		
	 * (iv) Append information on lists and merge ist to individual dataset
		 
			cd "$hauptpfad"

			preserve
			forvalues i=1975(4)2015{
			clear
			import excel using ".\0_original_data\NRW 1971_2015ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet("`i'") cellrange(A5:L445)
			tostring ListenNroffiziell, replace
			save "$hauptpfad\1_data\partyvotes`i'.dta", replace
			}



			use "$hauptpfad\1_data\partyvotes1971.dta", clear
			forvalues i=1975(4)2015{
			append using "$hauptpfad\1_data\partyvotes`i'.dta", force
			erase "$hauptpfad\1_data\partyvotes`i'.dta"
			}


			gen nationalratswahl=Wahljahr
			gen canton=Kanton
			gen listennummer=ListenNroffiziell
			
			drop if KantonsNr==. // drop empty cells
			save "$hauptpfad\1_data\partyvotes_all.dta", replace
			restore
			
			merge m:1 nationalratswahl listennummer  canton using "$hauptpfad\1_data\partyvotes_all.dta", gen(merge_lists)
			
			/* 
			
			Result                           # of obs.
			-----------------------------------------
			not matched                           167
				from master                         6  (merge_lists==1) -> THESE ARE ALL FROM MAJORITY VOTING CANTONS: AR (2 x 1979, 2 x 1987) AND OW (1999) AND NW (2007) 
				from using                        161  (merge_lists==2) -> THESE ARE ALL OBSERVATIONS FROM 1971

			matched                            29,654  (merge_lists==3)
			-----------------------------------------
			*/
			
			
		   	sort nationalratswahl canton_no listennummer kandidatennummer 
			
			drop if nationalratswahl=="NRW1971"
			
			save "$hauptpfad\nr_data_1975_2015.dta", replace


		******************************************************
		* (B) READ-IN ELECTION RESULTS DATA (1931 until 1971)
		******************************************************
		
		
		* (i) Import all excel files
		
		clear
		import excel "$hauptpfad/6_Bundesblatt/Bundesblatt in/1931_Bundesblatt_backin.xlsx"	,  first // start with 1931
		gen Wahljahr=1931
		replace canton=canton[_n-1] if canton==""
		saveold "$hauptpfad\1_data\temp.dta", version(12) replace
		
		forvalues y=1935(4)1971{
		clear
		// display "`y'"
		import excel "$hauptpfad/6_Bundesblatt/Bundesblatt in/`y'_Bundesblatt_backin.xlsx"	,  first
		replace canton=canton[_n-1] if canton==""
		append using "$hauptpfad\1_data\temp.dta"
		replace Wahljahr=`y' if Wahljahr==.
		saveold "$hauptpfad\1_data\temp.dta", version(13) replace
		}
		
		tab Wahljahr
		
		
		* (ii) Correct names of canton
		
		preserve
		clear
		import excel "$hauptpfad/1_data/Kantonsnamen_harmonisiert.xlsx"	,  first
		saveold "$hauptpfad\1_data\Kantonsnamen_harmonisiert.dta", version(13) replace
		restore
		
		merge m:1 canton using "$hauptpfad\1_data\Kantonsnamen_harmonisiert.dta", gen(merge_kantonsname)
		
		drop canton 
		gen canton=kantonsname

		merge m:1 canton using "$path_original_data/cantonal_translator",gen(merge_cantonaltranslator)
		drop if canton=="JU"

		
		
		* (iii) Correct party lists
		

		preserve
		egen party_num=group(party)
		collapse  (first) party, by(party_num canton)
		//export excel "$hauptpfad\1_data\Parteinamen Alle.xlsx", firstrow(variables) replace
		restore		
		
		
		* (iv) Add municipality numbers
		
		
		preserve
		egen city_num=group(city)
		collapse  (first) city, by(city_num canton)	
		export excel "$hauptpfad\1_data\Municipalities All.xlsx", firstrow(variables) replace  
		restore
		
		
		
		* (iv) Adapt names to 
		
		rename firstname vorname
		rename name nachname
		rename birth gebursjahr
		rename votes kandidst_tot
		
		gen gewaehlt=""
		replace gewaehlt="G" if elected="true"
		replace gewaehlt="N" if elected="false"
	
		
		

		******************************************************
		* (C) EXECUTE FUZZY MERGE (1931 until 2015)
		******************************************************
		
		
		
		******************************************************
		* (D) DATA PREPARATION
		******************************************************
		
		* (v)  Generate running variable
		
		
		* (a) Define votes and elected variables
		
		egen list_nr=group(nationalratswahl canton_no listennummer)
		egen canton_id=group(knr)
		gen individual_votes=kandidst_tot		
		gen list_votes=ErhalteneStimmen
		bysort list_nr: egen individual_votes_sum=total(individual_votes), missing
		
		gen individual_votes_perc= individual_votes/individual_votes_sum*100 // generate individual vote share in terms of total list votes
		gen elected=.
		replace elected=1 if gewaehlt=="G"
		replace elected=0 if gewaehlt=="N"
		
		corr list_votes individual_votes_sum
		
		sum list_votes if  list_votes==individual_votes_sum
 		sum list_votes if  list_votes!=individual_votes_sum

	 	br vorname nachname nationalratswahl canton listennummer list_votes individual_votes_sum 
		
		gen list_votes_diff=list_votes - individual_votes_sum 
	
		
		* (b) Sort candidates on list and define marginal candidate		
	
		gsort year canton_id  list_nr -individual_votes
		bysort year canton_id  list_nr: gen candidate_rank=_n // rank of candidate on list

		bysort year canton_id  list_nr: egen elected_perlist= sum(elected) // indicator if anybody is elected on list
		bysort year canton_id : egen elected_ct= sum(elected) 
				
		sort  year canton_id  list_nr candidate_rank
		gen elect_rank=elected*candidate_rank		// interaction rank and elected
		bysort year canton_id  list_nr: egen elect_rank_max=max(elect_rank)
		gen candidate_dist=candidate_rank-elect_rank_max // distance to marginal guy in positions

		gen marginal=0 // marginal guy from below (for all non-elected)
		replace marginal=1 if candidate_dist==0

		gen marginal_above=0 // marginal guy from above (for all elected)
		replace marginal_above=1 if candidate_dist==1

		
		* (c) Generate vote variable of marginal candidate per list
		
		
		gen individual_votes_marg=. 
		replace individual_votes_marg=individual_votes if marginal==1
		gen individual_votes_marg_above=. 
		replace individual_votes_marg_above=individual_votes if marginal_above==1
	
		gen individual_votes_perc_marg=. 
		replace individual_votes_perc_marg=individual_votes_perc if marginal==1
		gen individual_votes_perc_marg_above=. 
		replace individual_votes_perc_marg_above=individual_votes_perc if marginal_above==1
		
		bysort year canton_id  list_nr:egen  individual_votes_marg_all=max(individual_votes_marg) 
		bysort year canton_id  list_nr:egen  individual_votes_mar_ab_all=min(individual_votes_marg_above) 
		bysort year canton_id  list_nr:egen  individual_votes_perc_marg_all=max(individual_votes_perc_marg) 
		bysort year canton_id  list_nr:egen  individual_votes_perc_mar_ab_all=min(individual_votes_perc_marg_above) 
			
		
		* (d) Take difference in votes to marginal candidate 


		gen individual_diff_marginal=individual_votes-individual_votes_marg_all if elected==0
		replace individual_diff_marginal=individual_votes-individual_votes_mar_ab_all if elected==1

		gen individual_diff_marginal_perc=individual_votes_perc-individual_votes_perc_marg_all if elected==0
		replace individual_diff_marginal_perc=individual_votes_perc-individual_votes_perc_mar_ab_all if elected==1
		
		 br  if inrange(individual_diff_marginal,-10,10)
		 
		 sum if inrange(individual_diff_marginal,-100,100)


		
		// check 
		// br year canton_id  list_nr individual_diff_marginal individual_votes individual_votes_marg_all  elected
		
		sum elected_perlist if individual_diff_marginal==.
		
		br  year canton_id  canton list_nr individual_diff_marginal individual_votes individual_votes_marg_all  elected elected_perlist if individual_diff_marginal==. &  elected_perlist>0  

		
		// USE WITH ID
		
		sort id_Stata_new_num time_variable
		gen id_Stata_new_num=group(id_Stata_new)
		bysort id_Stata_new_num: gen nr_participation=_n
	
		gen candidate_votes_perc_lag=L1.candidate_votes_perc

		sort id_Stata_new_num time_variable
		gen candidate_diff_marginal_lag=L1.candidate_diff_marginal

		gen candidate_diff_marginal_lag_elec=candidate_diff_marginal_lag*elected
		
	
