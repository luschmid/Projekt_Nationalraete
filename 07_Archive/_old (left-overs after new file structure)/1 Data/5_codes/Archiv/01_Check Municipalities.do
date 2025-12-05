
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
	* (A) CHECK MUNICIPALITIES
	******************************************************


	
	
	
	
		use "$hauptpfad\1_data\nr_data_1931_1975_jan2017.dta", clear 
		
		
		list nachname	vorname year if city=="Arlesheim" & canton=="BS"
		list nachname	vorname year if city=="Basel" & canton=="BE"
		list nachname	vorname year if city=="Bellinzona" & canton=="VD"
		
		list nachname	vorname year if city=="Bern" & canton=="ZH"
		list nachname	vorname year if city=="Berna" & canton=="TI"
		list nachname	vorname year if city=="Berne" & canton=="FR"
		list nachname	vorname year if city=="Berne" & canton=="NE"
		list nachname	vorname year if city=="Berne" & canton=="VD"
		list nachname	vorname year if city=="Berne" & canton=="VS"

	
		list nachname	vorname year if city=="Bex" & canton=="VS"
		list nachname	vorname year if city=="Biberstein" & canton=="BE"
		list nachname	vorname year if city=="Binningen" & canton=="BE"
		list nachname	vorname year if city=="Boudevilliers" & canton=="GE"
		list nachname	vorname year if city=="Boudevilliers" & canton=="VD"
		list nachname	vorname year if city=="Busswil" & canton=="BE"

		list nachname	vorname year if city=="Cheseaux" & canton=="VD"
		
		list nachname	vorname year if city=="Davos Platz" & canton=="ZH"
		list nachname	vorname year if city=="Estavayer" & canton=="FR"
		list nachname	vorname year if city=="Galmiz" & canton=="FR"
		list nachname	vorname year if city=="Genf" & canton=="ZH"
		list nachname	vorname year if city=="Genève" & canton=="VD"
		list nachname	vorname year if city=="Guntershausen" & canton=="TG"

		list nachname	vorname year if city=="Köniz" & canton=="VD"
		list nachname	vorname year if city=="Köniz (Berne)" & canton=="VD"

		
		list nachname	vorname year if city=="Lausanne" & canton=="BE"
		list nachname	vorname year if city=="Lausanne" & canton=="NE"
		list nachname	vorname year if city=="Lausanne" & canton=="VS"
		list nachname	vorname year if city=="Luzern" & canton=="BE"

	
		list nachname	vorname year if city=="Mies (Vaud)" & canton=="GE"
		list nachname	vorname year if city=="Montier" & canton=="BE"
		list nachname	vorname year if city=="Muri BE" & canton=="GE"
		list nachname	vorname year if city=="Nods et Moutier" & canton=="BE"
		list nachname	vorname year if city=="Oberhünigen/Niederhünigen" & canton=="BE"

		list nachname	vorname year if city=="Romanel" & canton=="VD"
		list nachname	vorname year if city=="Röthenbach" & canton=="BE"
		list nachname	vorname year if city=="Rüti/Bülach" & canton=="ZH"
		
		
		list nachname	vorname year if city=="Siebnen" & canton=="SZ"
		list nachname	vorname city year if nachname=="Stähli" & vorname=="Fritz"
		list nachname	vorname city year if nachname=="Diethelm" & vorname=="Josef"
		list nachname	vorname city year gebursjahr canton party if nachname=="Kürzi" & vorname=="Josef"

		list nachname	vorname year if city=="Steffisberg" & canton=="BE"
		list nachname	vorname year if city=="Thun-Steffisburg" & canton=="BE"
		list nachname	vorname year if city=="Tona SG" & canton=="SG"
		list nachname	vorname year if city=="Visperterminen/Zermatt" & canton=="VS"

		list nachname	vorname year if city=="Willerzwil" & canton=="SZ"
		list nachname	vorname year if city=="Willisau" & canton=="LU"
	
		list nachname	vorname year gebursjahr if city=="Zollbruck" & canton=="BE"
		list nachname	vorname year gebursjahr if city=="Zollbrück" & canton=="BE"
	
		list nachname	vorname year if city=="Zollikerberg" & canton=="ZH"

			
		list nachname	vorname year if city=="Zurich" & canton=="VD"
		list nachname	vorname year if city=="Zurich" & canton=="VS"
		list nachname	vorname year if city=="Zurigo" & canton=="TI"
		list nachname	vorname year if city=="Zürich" & canton=="NE"
		list nachname	vorname year if city=="evey" & canton=="VD"
		list nachname	vorname year if city=="in Ascona" & canton=="BS"

		
		list nachname	vorname year if city=="in Basel" & canton=="AG"
		list nachname	vorname year if city=="in Basel" & canton=="BE"
		list nachname	vorname year if city=="in Basel" & canton=="LU"
		list nachname	vorname year if city=="in Basel" & canton=="SG"
		list nachname	vorname year if city=="in Basel" & canton=="SH"
		list nachname	vorname year if city=="in Basel" & canton=="ZH"	
		list nachname	vorname year if city=="in Basel-Riehen" & canton=="ZH"	
	
		list nachname	vorname year if city=="in Bellinzona" & canton=="GR"
		list nachname	vorname year if city=="in Bellinzona Canton de Vaud" & canton=="TI"	
		list nachname	vorname year if city=="in Bern-Wyler" & canton=="BE"	
		list nachname	vorname year if city=="in Biel-Benken" & canton=="BL"	

		list nachname	vorname year city if nachname=="Wyss" & vorname=="Paul"	

		list nachname	vorname year if city=="in Bienen" & canton=="BS"	
		list nachname	vorname year city if nachname=="Senn" & vorname=="Karl"	

		list nachname	vorname year if city=="in Brissago" & canton=="BS"	
	
		list nachname	vorname year if city=="in Bümpliz" & canton=="TI"	
		list nachname	vorname year if city=="in Büren" & canton=="BE"	
		list nachname	vorname year if city=="in Deisswil" & canton=="BE"	
		list nachname	vorname year if city=="in Genf" & canton=="BE"	
		list nachname	vorname year if city=="in Genf" & canton=="ZH"	
		list nachname	vorname year if city=="in Heiden" & canton=="BL"	
		list nachname	vorname year if city=="in Holzhof" & canton=="TG"	

		list nachname	vorname year if city=="in Lausanne" & canton=="BE"	
		list nachname	vorname year if city=="in Losanna" & canton=="TI"
		
		list nachname	vorname year if city=="in Matten" & canton=="BE"	
		list nachname	vorname year if city=="in Reussbühl-Littau" & canton=="LU"	
		list nachname	vorname year if city=="in Riehen" & canton=="BL"	

		list nachname	vorname year if city=="in Romanshorn" & canton=="AG"
		
		list nachname	vorname year if city=="in Rüschlikon" & canton=="BE"	
		list nachname	vorname year if city=="in Rüschlikon" & canton=="SG"
		list nachname	vorname year if city=="in Rüschlikon (Zürich)" & canton=="BE"	
		tab canton if nachname=="Duttweiler" & vorname=="Gottlieb"
		
		
		list nachname	vorname year if city=="in Schinznach" & canton=="AG"	
		list nachname	vorname year if city=="in St Gallen" & canton=="BE"	
		list nachname	vorname year if city=="in St. Gallen" & canton=="AG"	
		list nachname	vorname year if city=="in St. Gallen" & canton=="BE"	

		list nachname	vorname year if city=="in Tann-Rüti" & canton=="ZH"	

		list nachname	vorname year if city=="in Trümmelbach-Kleine Scheidegg" & canton=="BE"	
		list nachname	vorname year if city=="in Trümmelbach-Scheidegg" & canton=="BE"	
		list nachname	vorname year if city=="in Trümmelbach/Scheidegg" & canton=="BE"	
		
		list nachname	vorname year if city=="in Unterengstringen" & canton=="AG"	
		list nachname	vorname year if city=="in Veltheim" & canton=="SG"	
		
		list nachname	vorname year if city=="in Waldenburg" & canton=="BS"		

		list nachname	vorname year if city=="in Wiggiswil-Münchenbuchsee" & canton=="BE"		
	 
		list nachname	vorname year if city=="in Zollbrück" & canton=="BE"		
	 
		list nachname	vorname year if city=="in Zollbrück" & canton=="BE"		
	 


		
		* (ii) To check
		
		list nachname	vorname year if city=="Bellinzona" & canton=="VD"
		list nachname	vorname year if city=="Luzern" & canton=="BE"
		list nachname	vorname year if city=="Biel-Benken BL" & canton=="BL"
		list nachname	vorname year if city=="in Romanshorn" & canton=="AG"
		