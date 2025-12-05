

	cd "C:/Schmidlu/Dropbox/Projekt Nationalräte"

	import excel using "02_Processed_data/12_Record_linkage/02_Check_False_Negatives_in/Duplicates-Check.xlsx", first sheet("Tabelle1") clear

	drop if ID=="LU-2015-0151"
	drop if ID=="ZH-1999-0279"
	drop if ID=="ZH-2007-0331"
	drop if ID=="ZH-2011-0115" & Source_File=="Check_False_Negatives15.xlsx"
	drop if ID=="BEJU-2011-0035" & Source_File=="Check_False_Negatives15.xlsx"


	duplicates tag ID Worker, gen(dups)
	tab dups
	//br if dups>0
	drop if dups>0

	gen ID_Bisnode=Sicher+"/"+Unsicher
	split ID_Bisnode, p("/" "\\" "\" "//")
	rename  ID_Bisnode ID_Bis_old
	gen obs=_n
	reshape long ID_Bisnode,i(obs) j(num)  string
	bysort ID Worker: gen first=_n
	drop if ID_Bisnode=="" & first>1
	drop first
	gsort ID Worker -ID_Bisnode
	bysort ID Worker: gen num2=_n
	tostring num2, replace
	drop num
	sort ID_Bisnode
	reshape wide ID_Bisnode,i(obs) j(num2)  string
	gen ID_Bisnode_harm = ID_Bisnode1+"/"+ID_Bisnode2+"/"+ID_Bisnode3+"/"+ID_Bisnode4+"/"+ID_Bisnode5+"/"+ID_Bisnode6
	unique ID_Bisnode_harm,by(ID) gen(ID_Bisnode_unique)
	bysort ID: egen ID_Bisnode_unique_min=min(ID_Bisnode_unique)

	bysort Worker: tab ID_Bisnode_unique_min  

	sort ID Worker
	order ID ID_Bisnode1-ID_Bisnode6
	br ID ID_Bisnode1-ID_Bisnode6 Beruf Geburtsjahr Geschlecht Heimatort ID Kanton Kommentare Vorname Name Parteiname Wohngemeinde if ID_Bisnode_unique_min>1

	tab ID

	preserve
	collapse (first) ID_Bisnode_unique_min, by(ID)
	tab ID_Bisnode_unique_min
	restore

	keep ID ID_Bisnode_unique_min Worker
	bysort ID: gen num3=_n
	reshape wide ID_Bisnode_unique_min Worker, i(ID) j(num3)
	gen Pair=Worker1+" " +Worker2+" " +Worker3
	sort Pair ID
	gen correct=.
	replace correct=0 if ID_Bisnode_unique_min1>1
	replace correct=1 if ID_Bisnode_unique_min1==1
	keep ID Pair correct
	

	bysort Pair: tab correct
