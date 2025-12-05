


	global path "C:/Schmidlu/Dropbox/Record Linkage"
	global path_projekt_nationalraete "C:/Schmidlu/Dropbox/Projekt Nationalräte"
	
	** Aggregate Unsicher Check_Negatives done by Nemo Krüger and Laura Decet
	
	
		
	**********************
	* (A) Out files
	**********************	

	/*
	* Nemo Krüger
	local first 1
	
	local files : dir "$path/Code/Metatesting/Nemo Krueger/Bisnode" files "*.xlsx"
		foreach file in `files' {
		display	"`file'"
		import excel using "$path/Code/Metatesting/Nemo Krueger/Bisnode/`file'", first clear
		gen filename="`file'"
		capture rename K Sicher
		capture rename L Unsicher
		capture rename M Kommentare		
		capture tostring Sicher,replace
		capture tostring Unsicher,replace
		capture tostring Kommentare,replace
		
		if `first'==0{
					append using "$path/Code/Metatesting/Nemo Krueger/Bisnode/help.dta", force
					}
		save "$path/Code/Metatesting/Nemo Krueger/Bisnode/help.dta", replace
		local first 0
	}
	
	replace Sicher = subinstr(Sicher, " ", "", .)
	replace Unsicher = subinstr(Unsicher, " ", "", .)
	split Sicher,p("/")
	drop Sicher
	keep ID Sicher*
	bysort ID: gen indi=_n
	egen ID_num=group(ID)
	xtset ID_num indi
	reshape long Sicher*, i(ID_num) j(Mandate_No)
	
	rename Sicher Sicher_NemoKrueger
	rename Unsicher Unsicher_NemoKrueger
	
	keep ID Sicher_NemoKrueger Unsicher_NemoKrueger filename
	
	bysort filename: tab Sicher_NemoKrueger
	bysort filename: tab Unsicher_NemoKrueger
	
	save "$path/Code/Metatesting/Nemo Krueger/Bisnode/NemoKruegerAll.dta", replace
	
	
	* Laura Decet
	
	local files : dir "$path/Code/Metatesting/Laura Decet/Bisnode" files "*.xlsx"
		foreach file in `files' {
		display	"`file'"
		import excel using "$path/Code/Metatesting/Laura Decet/Bisnode/`file'", first clear
		gen filename="`file'"
		capture rename K Sicher
		capture rename L Unsicher
		capture rename M Kommentare
		if `first'==0{
					append using "$path/Code/Metatesting/Laura Decet/Bisnode/help.dta", force
					}
		save "$path/Code/Metatesting/Laura Decet/Bisnode/help.dta", replace
		local first 0
	}
	
	rename Sicher Sicher_LauraDecet
	rename Unsicher Unsicher_LauraDecet
	
	merge 1:1 ID using "$path/Code/Metatesting/Nemo Krueger/Bisnode/NemoKruegerAll.dta", gen(merge_krueger_dect)
	
	order ID */
	
	
	**********************
	* (B) Read in results
	**********************	
	
	** Nemo Krüger
	
		local first 1

		local files : dir "$path/Code/Metatesting/Nemo Krueger/Bisnode" files "*.xlsx"
		foreach file in `files' {
		display	"`file'"
		import excel using "$path/Code/Metatesting/Nemo Krueger/Bisnode/`file'", first clear
		capture tostring Unsicher, replace
		gen filename="`file'"
		gen double batch= real(regexs(1)) if regexm(filename,"([0-9]+)")	
		capture rename K Sicher
		capture rename L Unsicher
		capture rename M Kommentare		
		capture tostring Sicher,replace
		capture tostring Unsicher,replace
		capture tostring Kommentare,replace
		
		if `first'==0{
					append using "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_nemokrueger.dta", force
					}
		save "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_nemokrueger.dta", replace
		local first 0
	}
	
	preserve // This is for Arbeitsteilung 3. März 2020
	drop Sicher Unsicher filename batch Kommentare N O
	collapse (first) Kanton Geschlecht Name Vorname Geburtsjahr Wohngemeinde Heimatort Beruf Parteiname, by(ID)
	save "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_nemokrueger_all.dta", replace
	restore
	
	keep ID Sicher Unsicher batch filename Kommentare N O
	
	foreach x of varlist _all {
	rename `x' `x'_Krueger
	} 
	rename ID_Krueger ID
	rename Sicher_Krueger Sicher
	
	bysort ID: gen nobs=_N
	tostring batch, gen(batchstring)
	replace ID=ID+"-"+batchstring if nobs>1
	drop nobs

	replace Sicher = subinstr(Sicher, "\", "/",.) 
	replace Sicher = subinstr(Sicher, " ", "",.) 
	replace Unsicher_Krueger = subinstr(Unsicher_Krueger, "\", "/",.) 
	replace Unsicher_Krueger = subinstr(Unsicher_Krueger, " ", "",.) 
	split Sicher, p("/")
	drop Sicher
	
	preserve	
	keep ID filename_Krueger Sicher* Unsicher_Krueger batch_Krueger
	reshape long Sicher, i(ID) j(mandate_number)
	drop if Sicher==""
	duplicates drop ID Sicher, force
	gen double Sicher_num= real(regexs(1)) if regexm(Sicher,"([0-9]+)")	
	drop Sicher
	rename Sicher_num Sicher
	save "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_nemokrueger.dta", replace
	restore
	
	
	
	** Laura Decet
	
		clear
		local first 1

		local files : dir "$path/Code/Metatesting/Laura Decet/Bisnode" files "*.xlsx"
		foreach file in `files' {
		display	"`file'"
		import excel using "$path/Code/Metatesting/Laura Decet/Bisnode/`file'", first clear
		capture tostring Unsicher, replace
		gen filename="`file'"
		gen double batch= real(regexs(1)) if regexm(filename,"([0-9]+)")	
		capture rename K Sicher
		capture rename L Unsicher
		capture rename M Kommentare		
		capture tostring Sicher,replace
		capture tostring Unsicher,replace
		capture tostring Kommentare,replace
		
		if `first'==0{
					append using "$path/Code/Metatesting/Laura Decet/Bisnode/help_lauradecet.dta", force
					}
		save "$path/Code/Metatesting/Laura Decet/Bisnode/help_lauradecet.dta", replace
		local first 0
	}
	
	keep ID Sicher Unsicher batch filename Kommentare N O
	
	foreach x of varlist _all {
	rename `x' `x'_Decet
	} 
	rename ID_Decet ID
	rename Sicher_Decet Sicher
	
	bysort ID: gen nobs=_N
	tostring batch, gen(batchstring)	
	replace ID=ID+"-"+batchstring if nobs>1
	drop nobs
	
	replace Sicher = subinstr(Sicher, "\", "/",.) 
	replace Sicher = subinstr(Sicher, " ", "",.) 
	replace Unsicher_Decet = subinstr(Unsicher_Decet, "\", "/",.) 
	replace Unsicher_Decet = subinstr(Unsicher_Decet, " ", "",.) 
	split Sicher, p("/")
	drop Sicher
	
	preserve	
	keep ID filename_Decet Sicher* Unsicher_Decet batch_Decet
	reshape long Sicher, i(ID) j(mandate_number)
	drop if Sicher==""
	duplicates drop ID Sicher, force
	gen double Sicher_num= real(regexs(1)) if regexm(Sicher,"([0-9]+)")	
	drop Sicher
	rename Sicher_num Sicher
	save "$path/Code/Metatesting/Laura Decet/Bisnode/help_lauradecet.dta", replace
	restore
		

	
	* Arbeitsteilung 3. März 2020
	
	use "$path/Code/Metatesting/Laura Decet/Bisnode/help_lauradecet.dta", clear
	merge 1:1 ID Sicher using "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_nemokrueger.dta", gen(merge_krueger_decet)
	
	gen batch=batch_Krueger
	replace batch=batch_Decet if batch==.
	
	order ID Sicher merge_krueger_decet Unsicher_Decet Unsicher_Krueger	
	bysort ID: egen merge_krueger_decet_min=min(merge_krueger_decet)
	keep if merge_krueger_decet_min<3

	gen Source_Connection=""
	replace Source_Connection="Decet" if merge_krueger_decet==1
	replace Source_Connection="Krueger" if merge_krueger_decet==2
	replace Source_Connection="Beide" if merge_krueger_decet==3
	
	drop if batch<=10
	keep ID Source_Connection Sicher Unsicher_Decet Unsicher_Krueger batch	
	
	rename ID ID_new
	gen ID=substr(ID_new,1,12)
	replace ID=substr(ID_new,1,14) if substr(ID_new,1,4)=="BEJU"

	merge m:1 ID using "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_nemokrueger_all.dta", gen(merge_info_all)
	
	keep if merge_info_all==3

	gsort batch ID -Source_Connection Sicher
	order ID ID_new Unsicher_Decet Unsicher_Krueger
	
	export excel "$path_projekt_nationalraete/02_Processed_data/12_Record_linkage/nationalraete_tocheck_batch_11_15.xlsx",  firstrow(variables) replace
	
	
	
	
	** Evaluate results
	
	use "$path/Code/Metatesting/Laura Decet/Bisnode/help_lauradecet.dta", clear
	merge 1:1 ID Sicher using "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_nemokrueger.dta", gen(merge_krueger_decet)
	
	gen batch=batch_Krueger
	replace batch=batch_Decet if batch==.
	
	order ID Sicher merge_krueger_decet Unsicher_Decet Unsicher_Krueger	
	bysort ID: egen merge_krueger_decet_min=min(merge_krueger_decet)
	gsort batch ID -merge_krueger_decet Sicher
	br if merge_krueger_decet_min<3
	
	// Simon und Lukas coding
	
	drop if ID=="BEJU-1995-0573" & Sicher==1515414
	drop if ID=="BL-2007-0028" & Sicher==1111362
	drop if ID=="GR-1999-0020" & Sicher==283019
	drop if ID=="ZH-1995-0489" & Sicher==3268564
	drop if ID=="GE-2011-9152" & Sicher==2501893
	drop if ID=="BEJU-2007-0251" & Sicher==1221557
	drop if ID=="GR-1999-0011-5" & Sicher==79946
	drop if ID=="LU-2003-0074" & Sicher==1640691
	drop if ID=="SG-1999-0090" & Sicher==2543929
	drop if ID=="ZH-1999-0768" & Sicher==1041924
	drop if ID=="AG-2003-0079" & Sicher==3091985
	drop if ID=="TG-1995-0068" & Sicher==1224184
	drop if ID=="TG-1995-0068" & Sicher==1224186
	drop if ID=="TG-1995-0068" & Sicher==1269231
	drop if ID=="ZH-1999-0053" & Sicher==748614
	drop if ID=="ZH-2007-0178" & Sicher==1838491
	drop if ID=="ZH-2007-0178" & Sicher==2277664
	drop if ID=="ZH-2007-0662-7" & Sicher==541055
	drop if ID=="BEJU-2011-0160" & Sicher==1250458
	drop if ID=="FR-2007-0057" & Sicher==1409473
	drop if ID=="ZH-1991-0416-8" & Sicher==73565
	drop if ID=="AG-2015-0005" & Sicher==2744073
	drop if ID=="VS-2003-0047" & Sicher==3131028
	drop if ID=="BEJU-1999-0306" & Sicher==3148398
	drop if ID=="SG-1991-0034" & Sicher==918820
	drop if ID=="ZH-1991-0308" & Sicher==3209695
	drop if ID=="AG-1995-0068" & Sicher==1092632
	drop if ID=="BEJU-1991-0177" & Sicher==1172878
	drop if ID=="BEJU-1999-0305" & Sicher==1236988
	drop if ID=="BEJU-1999-0305" & Sicher==42742
	drop if ID=="LU-2003-0021" & Sicher==232900
	drop if ID=="SG-1991-0010" & Sicher==1025655
	drop if ID=="SG-2003-0116" & Sicher==485396
	drop if ID=="SG-2003-0116" & Sicher==3292133
	drop if ID=="SG-2011-0174" & Sicher==1305972
	drop if ID=="SZ-1995-0012" & Sicher==1255659
	drop if ID=="TG-2003-0052" & Sicher==3582764

	// Batch Lukas (after Arbeitsteilung 3. März 2020)
	drop if ID=="ZH-1991-0071" & Sicher==1533786
	drop if ID=="ZH-2003-0952" & Sicher==2115637
	drop if ID=="AG-1983-0074" & Sicher==233982
	drop if ID=="AG-1983-0074" & Sicher==803284
	drop if ID=="AG-1983-0074" & Sicher==810274
	drop if ID=="AG-1983-0074" & Sicher==1324707
	drop if ID=="AG-1999-0218" & Sicher==2487401
	drop if ID=="BEJU-1991-0425" & Sicher==1230496
	drop if ID=="GR-1999-0017" & Sicher==260505
	drop if ID=="SG-2003-0133" & Sicher==249938
	drop if ID=="ZG-2007-0014" & Sicher==1845155
	drop if ID=="ZH-1995-0487" & Sicher==1767449
	drop if ID=="ZH-1995-0632" & Sicher==3520451
	drop if ID=="AG-2011-0120" & Sicher==271526
	drop if ID=="AR-1995-0004" & Sicher==268825
	drop if ID=="FR-2015-0041" & Sicher==1406563

	// Batch Simon (after Arbeitsteilung 3. März 2020)

	drop if ID=="NE-1979-0012" & Sicher==392202
	drop if ID=="VD-1999-0142" & Sicher==3265982	
	drop if ID=="ZH-1991-0416-13" & Sicher==73565
	drop if ID=="ZH-2003-0286" & Sicher==3213298
	drop if ID=="ZH-2003-0571" & Sicher==605477
	drop if ID=="ZH-2003-0571" & Sicher==3561086
	drop if ID=="ZH-2003-0937" & Sicher==2125677
	drop if ID=="BEJU-1987-0288" & Sicher==1561615
	drop if ID=="BEJU-1987-0288" & Sicher==3289192
	drop if ID=="ZH-1991-0042" & Sicher==1692319
	drop if ID=="ZH-1995-0290" & Sicher==166375
	drop if ID=="ZH-1995-0290" & Sicher==252337
	drop if ID=="ZH-1995-0290" & Sicher==952575
	drop if ID=="ZH-1995-0290" & Sicher==964996
	drop if ID=="ZH-2003-0543" & Sicher==3362897
	drop if ID=="ZH-2011-0214" & Sicher==1138887
	drop if ID=="ZH-2011-0688" & Sicher==2198159
	drop if ID=="LU-1995-0006" & Sicher==308761
	drop if ID=="LU-1995-0006" & Sicher==665506
	drop if ID=="LU-1995-0006" & Sicher==928142
	
	* check match between Decet and Krüger
	tab merge_krueger_decet 
	* Result:  Initially we have 1,375 matched and 145 not matched connections. 
	*	    (102 connections only by Decet, 43 only by Krüger)
	* 	    Of these 145 non-matched connections, we deleted 70 connections 
	* 	    and kept 75 (48 by Decet, 26 by Krüger)

	
	
	
	
	* Check doubletten
	
	/*
	
	split ID,p("-")
	
	rename Sicher personenid
	rename ID ID_new
	gen ID=substr(ID_new,1,12)
	replace ID=substr(ID_new,1,14) if substr(ID_new,1,4)=="BEJU"
	
	keep if ID4!=""
 	keep ID ID_new  personenid
	sort ID personenid
	bysort ID: gen numbering=_n 
	reshape wide personenid ID_new, i(ID) j(numbering)
	
	
	order ID personenid*
	
	br if personenid1!=personenid2
	
	* Zwei Abweichungen der Duplikate: 
	* VD-1991-0168 695376: Gleicher Wohn- und Bürgerort, obwohl sehr weit 
	* auseinder; scheint gleiche Person zu sein: 			https://www.letemps.ch/suisse/vaud-polemique-autour-president-club-economique-socialiste

	* VS-2007-0073	895896
	* Gleiche Emailadresse: 
	* http://gemeinde.buerchen.ch/freizeitsport/freizeit/vereine.php?e_19_22_event=grosses-wildbuffet-in-buerchen
	* http://oberwallis.cvp.ch/menschen/thomas-lehner/
	br ID ID_new personenid 
	*/
	
	rename ID ID_new
	
	gen ID=substr(ID_new,1,12)
	replace ID=substr(ID_new,1,14) if substr(ID_new,1,4)=="BEJU"
	rename Sicher personenid

	keep ID personenid
	duplicates list ID personenid
	duplicates report ID personenid
	duplicates drop ID personenid, force

	//br if ID=="VS-2007-0073" & personenid==895896
	
	
	save "$path_projekt_nationalraete/02_Processed_data/12_Record_linkage/02_Check_False_Negatives_in/connections.dta", replace
	// Note: these are all individuals (IDs) from the NR dataset with (possibly)
	// 		 more than one connection (personenid) in the Bisnodedate 
	
	
	* Merge connenctions and non-connections
	
	use "$path/Code/Metatesting/Nemo Krueger/Bisnode/NemoKruegerAll.dta", clear
	replace ID="ZH-1999-0701" if regex(ID,"ZH-1999-0701")
	
	keep ID
	duplicates drop ID, force
	
	merge 1:m ID using 	 "$path_projekt_nationalraete/02_Processed_data/12_Record_linkage/02_Check_False_Negatives_in/connections.dta", gen(merge_bisnode)
	keep ID personenid
	
	outsheet ID personenid using "$path_projekt_nationalraete/02_Processed_data/12_Record_linkage/bisnode_falsenegatives.csv",  delimiter(";") replace


	
	*** erase all temp files
	erase "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_nemokrueger_all.dta"
	erase "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_nemokrueger.dta"
	erase "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_lauradecet.dta"

	
	/*
	use "$path/Code/Metatesting/Nemo Krueger/Bisnode/help_nemokrueger_all.dta", replace
	
	preserve
	import delimited "$path_projekt_nationalraete/02_Processed_data/12_Record_linkage/help_check_batch_11_15.csv", clear
	keep id batch sicher
	gen ID=substr(id,1,12)
	replace ID=substr(id,1,14) if substr(id,1,4)=="BEJU"
	collapse sicher batch, by(ID)
	save "$path_projekt_nationalraete/02_Processed_data/12_Record_linkage/help_check_batch_11_15.dta", replace
	restore
	
	bysort ID: gen indi=_n
	drop if indi>1
	
	merge 1:1 ID using "$path_projekt_nationalraete/02_Processed_data/12_Record_linkage/help_check_batch_11_15.dta"
	
	keep if _merge==3
	drop indi _merge
	sort batch ID
	 
	
	erase "$path_projekt_nationalraete/02_Processed_data/12_Record_linkage/help_check_batch_11_15.dta"
	


	
	

	


	


	 
	
	


	

	