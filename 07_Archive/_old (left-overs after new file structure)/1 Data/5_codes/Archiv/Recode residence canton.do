		* (iii) Recode residence canton
		
		
		// 1931-1971
		
		gen gemeindenummer_bfs=municipalityno
		merge m:1 gemeindenummer_bfs using "$hauptpfad\1_data\Gemeindebestand_1969.dta", gen(gemeindemerge)
		tab year if gemeindemerge==1 
		br if gemeindemerge==1 & year==1971 // result: this is a person from Moutier in the canton of Bern (Moutier vs. Moutier BE)
		gen residencecanton=Kanton 
		replace residencecanton="BE" if year==1971 & name=="Saucy" & 	firstname=="Bernard" &	birthyear==1938
		drop  HistNummer Kanton Bezirksnummer Bezirksname gemeindename DatumderAufnahme gemeindemerge
		
		// 1975-2015
		
		gen municipalityno_2015=municipalityno	
		merge m:1 municipalityno_2015 using "$hauptpfad\1_data\Gemeindebestand_20151231.dta", gen(gemeindemerge)
		tab year if gemeindemerge==1 
		br if gemeindemerge==1 
