		* (vi) Out for Heimatorte (Florence)
		
		preserve
		keep ID canton cantonno year name_merge firstname_merge birthyear_merge firstname name birthyear  job list municipalityno municipality sex
		drop if ID==""
		order ID canton cantonno year name_merge firstname_merge birthyear_merge firstname name birthyear  job list municipalityno municipality sex
		sort ID year 
		compress
		saveold  "$hauptpfad\1_data\8 origin data\nationalraete_origin.dta", replace version(13)
		restore
