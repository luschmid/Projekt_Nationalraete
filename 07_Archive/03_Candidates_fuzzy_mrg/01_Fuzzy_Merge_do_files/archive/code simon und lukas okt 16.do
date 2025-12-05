
	
	//  CODE VON SIMON UND LUKAS, 20.10.2016
	
	
	* (i) 2011 data
	
	cd "$datapath"
	use data_federal_temp.dta, clear
	
	egen Vorname_num=group(Vorname)
	
	bysort id_Stata: egen sd_year=sd(Geburtsjahr)
	bysort id_Stata: egen sd_Vorname_num=sd(Vorname_num)
	bysort id_Stata: egen sum_indicator=sum(indicator)

	sort id_Stata
	br candidateid_threshold_* indicator id_Stata  Vorname Nachname Geschlecht Geburtsjahr Gemeinde Wahljahr sd_year if sd_year!=0 & !missing(sd_Vorname_num) & sum_indicator==0 & !missing(sd_indicator) 
	br candidateid_threshold_* indicator id_Stata  Vorname Nachname Geschlecht Geburtsjahr Gemeinde Wahljahr sd_year if sd_year!=0 & !missing(sd_year) & sum_indicator==0 & !missing(sd_indicator) 

	
	* (ii) Old 2011 checked files

	preserve
	cd "$datapath\2 fuzzy merge back in\endgültig fertig"
	import excel  "erf_AG.xlsx"	, sheet(Sheet1) clear first
	
	local variables indicator	candidateid_threshold_opt	candidateid_threshold_min	candidateid_threshold_max	candidateid_	id_Stata	id_recoded	unsicher	Kommentar	Listenbezeichnung	
	foreach var in `variables'{	
	rename `var' `var'_2011
	saveold erf_AG_upto2011.dta, version(12) replace
	}
	restore
	
	cd "$datapath\2 fuzzy merge back in\endgültig fertig"
	merge 1:1 Vorname	Nachname Geschlecht	Geburtsjahr Gemeinde Wahljahr using erf_AG_upto2011.dta, gen(merge_fuzzy)
	
	
	list id_Stata id_Stata_2011 Vorname Nachname Geschlecht Geburtsjahr Gemeinde Wahljahr if id_Stata!=id_Stata_2011 & Wahljahr<2015
	list id_Stata id_Stata_2011 Vorname Nachname Geschlecht Geburtsjahr Gemeinde Wahljahr if id_Stata!=id_Stata_2011 & Wahljahr<2015
	
	bysort id_Stata: egen sd_year=sd(Geburtsjahr)
	bysort id_Stata: egen sd_indicator=sd(indicator)

	sort id_Stata
	br candidateid_threshold_* indicator id_Stata  Vorname Nachname Geschlecht Geburtsjahr Gemeinde Wahljahr sd_year if sd_year!=0 & !missing(sd_year) & sd_indicator==0 & !missing(sd_indicator) 
	

	
	
	
	gen indicator_different_id=0
	replace indicator_different_id=1 if id_Stata != id_Stata_2011 
	
	bysort candidateid_threshold_opt: egen indicator_different_id_mean1=mean(indicator_different_id)
	bysort candidateid_threshold_opt_2011: egen indicator_different_id_mean2=mean(indicator_different_id)
	
	gen indicator_different_id_mean=indicator_different_id_mean1+indicator_different_id_mean1
	
	br id_Stata id_Stata_2011 Vorname Nachname Geschlecht Geburtsjahr Gemeinde Wahljahr if indicator_different_id_mean!=0
	
	gen tocheck=0
	
	ENDE CODE VON SIMON UND LUKAS, 20.10.2016
	
	*/
	

	* (ii) Merge votes to previous fuzzy merge evaluation

	
	
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
	
	
