

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
	
	local canton ZH  

	
	cd "$path_fuzzy_data"
	use `canton'_all.dta, clear // infile
	gen yearNew=year
	
	local vars lastname firstname gender birthyear town
	foreach var in `vars'{
	rename `var' `var'_
	}
	
	merge 1:1 yearNew  lastname_ firstname_ gender_ birthyear_ using fuzzy_out_`canton'.dta, gen(merge_original)
	// note: replace using dataset fuzzy_out_`canton'.dta with the one processed by typist
	
	
	
	br if candidateid_==550
	br if candidateid_==551
	br if candidateid_==552
	br if lastname_=="Schelbert"
