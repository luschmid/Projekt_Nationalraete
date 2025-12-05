
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
