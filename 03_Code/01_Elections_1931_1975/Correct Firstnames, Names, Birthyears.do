		
		
		******************************************************
		* (I) FIRSTNAMES
		******************************************************
		
		// (a) change wrong firstname
		
	    replace firstname="Walter" if firstname=="Stöckli" & canton=="UR" & year==1975 // note: change of lastname by Simon Lüchinger in his recodings
		

		// (b) change other cases for which firstname differs within the same ID
	
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
		replace firstname="Katrin" if name=="Spahr" & canton=="AG" & year==1995



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
		replace firstname="Felix Emmanuel" if firstname=="Felix Emanuel" &  canton=="BS" & year==1947 & name=="Iselin"
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

		
		
		// (c) change cases with parenthesis in firstname	
		// rule: change obvious abbreviations, leave cases with additional names (e.g. Alfredo (Benjamin), Claudia (Dana) etc)
		

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
		replace firstname="Walter" if firstname=="Walter (Walti)"& year==2015 // 2 changes, both are correct
		
		
		// (d) correct and add information on special cases (after browsing all unique firstname)
		
		

		replace firstname="Elmar Theodor" if firstname=="Elmar Th."
		replace firstname="Gebhard" if firstname=="Gebhard Jun."
		replace firstname="Hans-Jürg" if firstname=="H. J."
		replace firstname="Hans R." if firstname=="Hans R"
		replace firstname="Hans" if firstname=="Hans, Jun."
		replace firstname="Hans U." if firstname=="Hans-U."
		replace firstname="Hans-Jürg" if firstname=="Hans.Jürg"
		replace firstname="Hans Ulrich" if firstname=="Hans.Ulrich"
		replace firstname="Heinz P." if firstname=="Heinz.P"
		replace firstname="Johann Jakob" if firstname=="J. J." // three changes, all are correct
		replace firstname="Jean-Jacques" if firstname=="J.-Jacques"
		replace firstname="Jakob" if firstname=="Jakob Jun."
		replace firstname="Jeanine Renée" if firstname=="Jeanine.Renée"
		replace firstname="Johann Christ." if firstname=="Joh. Christ."
		replace firstname="Johann Otto" if firstname=="Joh. Otto"
		replace firstname="Josef A." if firstname=="Jos. A."
		replace firstname="Josef" if firstname=="Josef, Sen."
		replace firstname="Maurice Nicolas" if firstname=="Maurice Nic."
		replace firstname="Max Johann" if firstname=="Max Joh."
		replace firstname="Otto" if firstname=="Otto, Jun."
		replace firstname="Roland W." if firstname=="Roland.W"

		// (e) check cases with only one letter in first name (subfolder abbreviation in check names.xlsx) 
		
		replace firstname="F. Emmanuel" if firstname=="F."& name=="Iselin"
		replace firstname="Alexander" if firstname=="A."&name=="Sarasin"
		replace firstname="Benedikt" if firstname=="B."&name=="Mani"
		replace firstname="Charles" if firstname=="C."&name=="Virchaux"		
		
		
		
		******************************************************
		* (II) NAMES 
		******************************************************
		
		
		replace name="von Waldkirch" if name=="v. Waldkirch"
		replace name="von Almen" if name=="v. Almen"
		replace name="von Grünigen" if name=="v. Grünigen"
		replace name="von Siebenthal" if name=="v. Siebenthal" // 2 changes, both are correct
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
		
		replace name="Baenziger" if name=="Bänzinger" & canton=="AG" & year==1959
		replace name="Heuberger-Kellenberger" if name=="Heuberger-Keller" & canton=="AG" & year==1999
		replace name="Bouchat" if name=="Bouchât" & canton=="BE" & year==1935
		replace name="Fell" if name=="Pell" & canton=="BE" & year==1931
		replace name="Hueber" if name=="Huber" & firstname=="Alfred" & canton=="BE" & year==1959
		replace name="Mertenat" if name=="Mertenet" & canton=="BE" & year==1955
		replace name="Schlapbach" if name=="Schläpbach" & canton=="BE" & year==1951
		replace name="Düby" if name=="Dübi" & canton=="BE" & year==1967
		replace name="Rutishauser" if name=="Rütishauser" & canton=="BE" & year==1955
		replace name="Borel" if name=="Borei" & canton=="BE" & year==1963
		replace name="Rosselet" if name=="Rosseiet" & canton=="BE" & year==1963
		replace name="Sauthier" if name=="Sautier" & canton=="BE" & year==1967
		replace name="Wermeille" if name=="Vermeille" & canton=="BE" & year==1967
		replace name="Brotbeck" if name=="Brodbeck" & canton=="BE" & year==1967
		replace name="Kästli" if name=="Kästii" & canton=="BE" & year==1967
		replace name="Casetti" if name=="Cassetti" & canton=="BE" & year==1967
		replace name="Marchand" if name=="Marchant" & canton=="BE" & year==1967
		replace name="Biétry" if name=="Bietry" & canton=="JU" & year==1979
		replace name="Guillaume-Gentil" if name=="Guillaume-Gentil (dit Gentil)" & canton=="JU" & year==1991
		replace name="Tanner-Gerber" if name=="Tanner Gerber" & canton=="BE" & year==2003
		replace name="Haefely" if name=="Haefeli" & canton=="BE" & year==2007
		replace name="Perrella" if name=="Perella" & canton=="BE" & year==2011
		replace name="von Kaenel" if name=="Von Kaenel" & canton=="BE" & year==2015
		replace name="von Blarer" if name=="Von Blarer" & canton=="BL" & year==1931
		replace name="Rutschi" if name=="Butschi" & canton=="BL" & year==1947
		replace name="van Baerle" if name=="Van Baerle" & canton=="BL" & year==1943
		replace name="Jenni" if name=="Jenny" & canton=="BL" & year==1947
		replace name="Jaques" if name=="Jacques" & canton=="BL" & year==1975
		replace name="Hanhart" if name=="Hankart" & canton=="BS" & year==1935
		replace name="Dürrenmatt" if name=="Dürremantt" & canton=="BS" & year==1955
		replace name="Haefeli" if name=="Haefely" & canton=="BS" & year==1959
		replace name="Jaquier" if name=="Jacquier" & canton=="FR" & year==1959
		replace name="Steinmann" if name=="Steimann" & canton=="FR" & year==1951
		replace name="Steinmann" if name=="Steimann" & canton=="FR" & year==1959
		replace name="Lanthemann" if name=="Lanthmann" & canton=="FR" & year==1967
		replace name="Friedly" if name=="Friedli" & canton=="FR" & year==1967
		replace name="Herrli" if name=="Heerli" & canton=="FR" & year==1971
		replace name="Dusseiller" if name=="Dussciller" & canton=="GE" & year==1947
		replace name="Déthiollaz" if name=="Dethiollaz" & canton=="GE" & year==1935
		replace name="Trüb" if name=="Trab" & canton=="GE" & year==1951
		replace name="Trüb" if name=="Trub" & canton=="GE" & year==1959
		replace name="Trüb" if name=="Trub" & canton=="GE" & year==1963
		replace name="Thévenaz" if name=="Thévénaz" & canton=="GE" & year==1975
		replace name="Thévenaz" if name=="Thevenaz" & canton=="GE" & year==1983
		replace name="Félix" if name=="Felix" & canton=="GE" & year==1979
		replace name="Maulini-Dreyfus" if name=="Maulini (Dreyfus)" & canton=="GE" & year==2015
		replace name="de Tassigny" if name=="De Tassigny" & canton=="GE" & year==2003
		replace name="Pétroz" if name=="Petroz" & canton=="GE" & year==2007
		replace name="Kasteler-Budde" if name=="Kasteler (Kasteler-Budde)" & canton=="GE" & year==2015
		replace name="Stebler Guenat" if name=="Stebler Guenat (Stebler)" & canton=="GE" & year==2015
		replace name="Mettier" if name=="Mettler" & canton=="GR" & year==1943
		replace name="Borgula" if name=="Borgola" & canton=="LU" & year==1983
		replace name="Pekdemir-Demirci" if name=="Pekdemir-Demirici" & canton=="LU" & year==2011
		replace name="Krügel" if name=="Krugel" & canton=="NE" & year==1931
		replace name="Krügel" if name=="Krugel" & canton=="NE" & year==1935
		replace name="Beyer" if name=="Beyer " & canton=="SG" & year==1935
		replace name="Kuhn" if name=="Kühn" & canton=="SG" & year==1939
		replace name="Wälter" if name=="Walter" & canton=="SG" & year==1935
		replace name="Eichenberger" if name=="Eicbenberger" & canton=="SG" & year==1935
		replace name="Ritter" if name=="Bitter" & canton=="SG" & year==1939
		replace name="Klingler" if name=="Klinger" & canton=="SG" & year==1947
		replace name="Hoby" if name=="Hobi" & canton=="SG" & year==1947
		replace name="Künzler" if name=="Kümzler" & canton=="SG" & year==1959
		replace name="Spreiter" if name=="Streiter" & canton=="SG" & year==1959
		replace name="Siegrist" if name=="Sigrist" & canton=="SO" & year==1931
		replace name="Stebler" if name=="Stehler" & canton=="SO" & year==1959
		replace name="Kräuchi" if name=="Krauchi" & canton=="SO" & year==1963
		replace name="Schürmann" if name=="Schurmann" & canton=="SO" & year==1963
		replace name="Bühler" if name=="Buhler" & canton=="TG" & year==1951
		replace name="Schümperli" if name=="Schumperli" & canton=="TG" & year==1951
		replace name="Rodel" if name=="Robel" & canton=="TG" & year==1959
		replace name="Labhart" if name=="Labhard" & canton=="TG" & year==1963
		replace name="Thür" if name=="Thur" & canton=="TG" & year==1951
		replace name="Thür" if name=="Thur" & canton=="TG" & year==1963
		replace name="Dollfus" if name=="Dollfuss" & canton=="TI" & year==1931
		replace name="Dollfus" if name=="Dollfuss" & canton=="TI" & year==1935
		replace name="Krähenbühl" if name=="Krähenbuhl" & canton=="TI" & year==1975
		replace name="von Wyttenbach" if name=="Von Wyttenbach" & canton=="TI" & year==2003
		replace name="Stöckli" if name=="Walter" & canton=="UR" & year==1975
		replace name="Viquerat" if name=="Viquerat " & canton=="VD" & year==1935
		replace name="von der Aa" if name=="Von der Aaa" & canton=="VD" & year==1931
		replace name="von der Aa" if name=="Von der Aa" & canton=="VD" & year==1943
		replace name="Fauquex" if name=="Fauquez" & canton=="VD" & year==1935
		replace name="Fauquex" if name=="Fauquez" & canton=="VD" & year==1939
		replace name="von Arx" if name=="Von Arx" & canton=="VD" & year==1951
		replace name="Debétaz" if name=="Debetaz" & canton=="VD" & year==1963
		replace name="Décosterd" if name=="Decosterd" & canton=="VD" & year==1959
		replace name="Décosterd" if name=="Decosterd" & canton=="VD" & year==1967
		replace name="Décosterd" if name=="Decosterd" & canton=="VD" & year==1971
		replace name="Baechtold" if name=="Beachtold" & canton=="VD" & year==1959
		replace name="Hediger" if name=="Hédiger" & canton=="VD" & year==1963
		replace name="Miéville" if name=="Mieville" & canton=="VD" & year==1979 //2 changes, both changes are correct
		replace name="Menétrey" if name=="Menetrey" & canton=="VD" & year==1971
		replace name="Menétrey" if name=="Menetrey" & canton=="VD" & year==1979
		replace name="Menétrey" if name=="Menetrey" & canton=="VD" & year==1983
		replace name="Menétrey Savary" if name=="Ménétrey Savary" & canton=="VD" & year==1999
		replace name="Clément" if name=="Clement" & canton=="VD" & year==1979
		replace name="Ganière" if name=="Ganiere" & canton=="VD" & year==1975
		replace name="Ganière" if name=="Ganiere" & canton=="VD" & year==1983
		replace name="Brélaz" if name=="Brelaz" & canton=="VD" & year==1975
		replace name="Brélaz" if name=="Brelaz" & canton=="VD" & year==1979
		replace name="van Singer" if name=="Van Singer" & canton=="VD" & year==2007
		replace name="van Singer" if name=="Van Singer" & canton=="VD" & year==2011
		replace name="van Singer" if name=="Van Singer" & canton=="VD" & year==2015
		replace name="Décosterd" if name=="Decosterd" & canton=="VD" & year==1995
		replace name="Cherix" if name=="Chérix" & canton=="VD" & year==2011
		replace name="Cretin-Meylan" if name=="Crétin-Meylan" & canton=="VD" & year==2011
		replace name="Métraux" if name=="Metraux" & canton=="VD" & year==1979
		replace name="Montangero" if name=="Montangéro" & canton=="VD" & year==1999
		replace name="de Preux" if name=="De Preux" & canton=="VD" & year==2007
		replace name="de Preux" if name=="De Preux" & canton=="VD" & year==2011
		replace name="von Siebenthal" if name=="Von Siebenthal" & canton=="VD" & year==2007
		replace name="Bréchet" if name=="Brechet" & canton=="VD" & year==2007
		replace name="Misiego" if name=="Misiégo" & canton=="VD" & year==2011
		replace name="Walter" if name=="Walther" & canton=="VS" & year==1931
		replace name="Anthamatten" if name=="Anthammatten" & canton=="VS" & year==1955
		replace name="de Riedmatten" if name=="De Riedmatten" & canton=="VS" & year==2003
		replace name="Hoppeler" if name=="Hoppeier" & canton=="ZH" & year==1935
		replace name="Helbig" if name=="Heibig" & canton=="ZH" & year==1935
		replace name="Kappeler" if name=="Kappeier" & canton=="ZH" & year==1935
		replace name="Lechleiter" if name=="Lechleitner" & canton=="ZH" & year==1935
		replace name="Kunz" if name=="Künz" & canton=="ZH" & year==1951
		replace name="Munz" if name=="Münz" & canton=="ZH" & year==1955
		replace name="Günthart" if name=="Günthardt" & canton=="ZH" & year==1943
		replace name="Günthart" if name=="Günthardt" & firstname == "Gottfried" & canton=="ZH" & year==1967
		replace name="Elber" if name=="Eiber" & canton=="ZH" & year==1955
		replace name="Frei" if name=="Frey" & firstname == "Josef Erhard" & canton=="ZH" & year==1955
		replace name="Bosshard" if name=="Bossard" & canton=="ZH" & year==1963
		replace name="Flueler" if name=="Flüeler" & canton=="ZH" & year==1955
		replace name="König" if name=="Köng" & canton=="ZH" & year==1967
		replace name="Kloter" if name=="Klotter" & canton=="ZH" & year==1967
		replace name="Raurich" if name=="Raurisch" & canton=="ZH" & year==1963
		replace name="Rüegg" if name=="Ruegg" & canton=="ZH" & year==1963
		replace name="Vonrufs" if name=="von Rufs" & canton=="ZH" & year==1967
		replace name="Budliger" if name=="Budlinger" & canton=="ZH" & year==1983
		replace name="Stöckli" if name=="Stöcklin" & canton=="ZH" & year==1979
		replace name="Berginz " if name=="Bergienz" & canton=="ZH" & year==1983
		replace name="Günthardt" if name=="Günthart" & canton=="ZH" & year==1983
		replace name="Wäfler" if name=="Wäffler" & canton=="ZH" & year==1983
		replace name="Gutzwiller" if name=="Gutzwiler" & canton=="ZH" & year==1991
		replace name="Nüssli" if name=="Nüssli jun." & canton=="ZH" & year==1995
		replace name="Ilija-Wildi" if name=="Jlija-Wildi" & canton=="ZH" & year==2003
		replace name="von Felten" if name=="Von Felten" & canton=="ZH" & year==2003
		replace name="Sangines" if name=="Sanginés" & canton=="ZH" & year==2011
		replace name="Umbricht-Décosterd" if name=="Unbricht-Décosterd" & canton=="AG" & year==1979
	

		******************************************************
		* (III) BIRTHYEARS 
		******************************************************
		
				
		replace birthyear=1884 if birthyear==1894 & name=="Zubler" & firstname=="Ernst" & year==1935
		replace birthyear=1886 if birthyear==1887 & name=="Glaser" & firstname=="Alfred" & year==1939
		replace birthyear=1886 if birthyear==1866 & name=="Held" & firstname=="Alfred" & year==1931
		replace birthyear=1897 if birthyear==1895 & name=="Weber" & firstname=="Max" & year==1947
		replace birthyear=1897 if birthyear==1895 & name=="Weber" & firstname=="Max" & year==1951
		replace birthyear=1898 if birthyear==1889 & name=="Brawand" & firstname=="Samuel" & year==1943
		replace birthyear=1901 if birthyear==1900 & name=="Burren" & firstname=="Ernst" & year==1939
		replace birthyear=1892 if birthyear==1890 & name=="Spieler-Gerster" & firstname=="Joseph" & year==1943
		replace birthyear=1869 if birthyear==1870 & name=="Grimaître" & firstname=="Alcide" & year==1939
		replace birthyear=1886 if birthyear==1896 & name=="Kasser" & firstname=="Walter" & year==1947
		replace birthyear=1899 if birthyear==1900 & name=="Schmid" & firstname=="Ernst" & year==1947
		replace birthyear=1903 if birthyear==1908 & name=="Tschumi" & firstname=="Ernst" & year==1939
		replace birthyear=1921 if birthyear==1920 & name=="Macquat" & firstname=="Roger" & year==1963
		replace birthyear=1900 if birthyear==1901 & name=="Rupp" & firstname=="August" & year==1947
		replace birthyear=1902 if birthyear==1903 & name=="Bohrer" & firstname=="Franz" & year==1955
		replace birthyear=1922 if birthyear==1908 & name=="Zihlmann" & firstname=="Louis" & year==1959
		replace birthyear=1918 if birthyear==1910 & name=="Veya" & firstname=="Raymond" & year==1967
		replace birthyear=1958 if birthyear==1957 & name=="Hofer-Scherer" & firstname=="Marc" & year==1995
		replace birthyear=1949 if birthyear==1948 & name=="Seydoux" & firstname=="André" & year==1991
		replace birthyear=1958 if birthyear==1959 & name=="Hirt" & firstname=="Ursula" & year==1991
		replace birthyear=1949 if birthyear==1947 & name=="Walz" & firstname=="Brigitte" & year==1991
		replace birthyear=1899 if birthyear==1898 & name=="Rutschi" & firstname=="Gottlieb" & year==1935
		replace birthyear=1903 if birthyear==1904 & name=="Ryser" & firstname=="Albert" & year==1951
		replace birthyear=1919 if birthyear==1924 & name=="Frey" & firstname=="Hans" & year==1959
		replace birthyear=1926 if birthyear==1925 & name=="Breitenstein" & firstname=="Willi" & year==1975
		replace birthyear=1939 if birthyear==1921 & name=="Gasser" & firstname=="Thomas" & year==1983
		replace birthyear=1904 if birthyear==1914 & name=="Würz" & firstname=="Alfred" & year==1951
		replace birthyear=1898 if birthyear==1899 & name=="Zschokke" & firstname=="Peter" & year==1943
		replace birthyear=1954 if birthyear==1935 & name=="Inglin-Buomberger" & firstname=="Beatrice" & year==1991
		replace birthyear=1954 if birthyear==1957 & name=="Oppliger-Schenker" & firstname=="Silvia" & year==1995
		replace birthyear=1883 if birthyear==1888 & name=="Pochon" & firstname=="Marc" & year==1931
		replace birthyear=1892 if birthyear==1891 & name=="Pasquier" & firstname=="Albert" & year==1947
		replace birthyear=1941 if birthyear==1951 & name=="Krauskopf" & firstname=="Eveline" & year==1999
		replace birthyear=1898 if birthyear==1895 & name=="Pugin" & firstname=="Antoine" & year==1939
		replace birthyear=1921 if birthyear==1927 & name=="Jacquiard-Renevier" & firstname=="Jacqueline" & year==1987
		replace birthyear=1881 if birthyear==1887 & name=="Branger" & firstname=="Erhard" & year==1931
		replace birthyear=1916 if birthyear==1914 & name=="Mayer" & firstname=="Oskar" & year==1951
		replace birthyear=1939 if birthyear==1938 & name=="Obrist" & firstname=="Robert" & year==1999
		replace birthyear=1961 if birthyear==1963 & name=="Flütsch" & firstname=="Christiana" & year==2003
		replace birthyear=1863 if birthyear==1868 & name=="Ackermann" & firstname=="Isidor" & year==1931
		replace birthyear=1873 if birthyear==1878 & name=="Grünenfelder" & firstname=="Emil" & year==1931
		replace birthyear=1883 if birthyear==1888 & name=="Duft" & firstname=="Johannes" & year==1931
		replace birthyear=1894 if birthyear==1882 & name=="Eichenberger" & firstname=="Rudolf" & year==1935
		replace birthyear=1889 if birthyear==1888 & name=="Kappler" & firstname=="Arnold" & year==1939
		replace birthyear=1903 if birthyear==1896 & name=="Gantenbein" & firstname=="Florian" & year==1943
		replace birthyear=1913 if birthyear==1931 & name=="Weber" & firstname=="Joachim" & year==1971
		replace birthyear=1883 if birthyear==1893 & name=="Borella" & firstname=="Francesco" & year==1943
		replace birthyear=1888 if birthyear==1883 & name=="Gasparini" & firstname=="Amilcare" & year==1935
		replace birthyear=1897 if birthyear==1896 & name=="Bordoni" & firstname=="Rodolfo" & year==1947
		replace birthyear=1918 if birthyear==1919 & name=="Stefani" & firstname=="Alberto" & year==1963
		replace birthyear=1887 if birthyear==1987 & name=="Muheim" & firstname=="Karl" & year==1935
		replace birthyear=1897 if birthyear==1887 & name=="Arnold" & firstname=="Franz" & year==1947
		replace birthyear=1877 if birthyear==1879 & name=="Masson" & firstname=="Eugène" & year==1931
		replace birthyear=1883 if birthyear==1888 & name=="Perrin-Cherbuin" & firstname=="Ernest" & year==1931
		replace birthyear=1890 if birthyear==1891 & name=="Miéville" & firstname=="Adrien" & year==1947
		replace birthyear=1927 if birthyear==1922 & name=="Curdy" & firstname=="Jean" & year==1967
		replace birthyear=1945 if birthyear==1946 & name=="Jaccard" & firstname=="Marianne" & year==1979
		replace birthyear=1961 if birthyear==1966 & name=="Buffat" & firstname=="Marc-Olivier" & year==2015
		replace birthyear=1909 if birthyear==1919 & name=="Germanier" & firstname=="Francis" & year==1959
		replace birthyear=1916 if birthyear==1918 & name=="Escher" & firstname=="Alfred" & year==1963
		replace birthyear=1873 if birthyear==1878 & name=="Ryffel" & firstname=="Albert" & year==1931
		replace birthyear=1893 if birthyear==1898 & name=="Stiefel" & firstname=="Heinrich" & year==1931
		replace birthyear=1905 if birthyear==1903 & name=="Pfister" & firstname=="Paul" & year==1963
		replace birthyear=1951 if birthyear==1961 & name=="Goldberger" & firstname=="Liliane" & year==1987
		replace birthyear=1961 if birthyear==1960 & name=="Holzreuter" & firstname=="Daniel" & year==1991
		replace birthyear=1964 if birthyear==1969 & name=="Hug" & firstname=="Markus" & year==1991

		