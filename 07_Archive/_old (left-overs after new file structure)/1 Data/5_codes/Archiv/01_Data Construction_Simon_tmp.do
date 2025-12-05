use "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\1 Data\1_data\nationalraete_1931_2015.dta", clear

		* (iv) Harmonization of names, firstnames and birthyears
		
		// rules: 1. If more than two observations: follow majority spelling and/or correct spelling in Bundesblatt, cantonal data, or other external sources
		//		  2. If only two observations: check spelling in Bundesblatt (or cantonal data and other external sources) and if they differ we do not change anything
		//		  3. We do not harmonize Umlaute and double names
		//		  4. Change "v." to "von", "Von" to "von", "Van" to "van", "De" to "de"

		//br if toCheck==1 & substr(name,1,3)=="v. "
		
		replace name="von Waldkirch" if name=="v. Waldkirch"
		replace name="von Almen" if name=="v. Almen"
		replace name="von Grünigen" if name=="v. Grünigen"
		replace name="von Siebenthal" if name=="v. Siebenthal"
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
		replace name="Frei" if name=="Frey" & firstname == "Josef" & canton=="ZH" & year==1955
		replace name="Bosshard" if name=="Bossard" & canton=="ZH" & year==1963
		replace name="Flueler" if name=="Flüeler" & canton=="ZH" & year==1955
		replace name="König" if name=="Köng" & canton=="ZH" & year==1967
		replace name="Kloter" if name=="Klotter" & canton=="ZH" & year==1967
		replace name="Raurich" if name=="Raurisch" & canton=="ZH" & year==1963
		replace name="Rüegg" if name=="Ruegg" & canton=="ZH" & year==1963
		replace name="von Rufs" if name=="von Rufs" & canton=="ZH" & year==1967
		replace name="Budliger" if name=="Budlinger" & canton=="ZH" & year==1983
		replace name="Stöckli" if name=="Stöckli" & canton=="ZH" & year==1979
		replace name="Berginz " if name=="Bergienz" & canton=="ZH" & year==1983
		replace name="Günthardt" if name=="Günthart" & canton=="ZH" & year==1983
		replace name="Wäfler" if name=="Wäffler" & canton=="ZH" & year==1983
		replace name="Gutzwiller" if name=="Gutzwiler" & canton=="ZH" & year==1991
		replace name="Nüssli" if name=="Nüssli jun." & canton=="ZH" & year==1995
		replace name="Ilija-Wildi" if name=="Jlija-Wildi" & canton=="ZH" & year==2003
		replace name="von Felten" if name=="Von Felten" & canton=="ZH" & year==2003
		replace name="Sangines" if name=="Sanginés" & canton=="ZH" & year==2011

		gen toCheck=0
		
		egen name_num=group(name)	
		bysort ID: egen name_num_sd=sd(name_num)
		sort ID year
		br if name_num_sd>0 & !missing(name_num_sd)
		replace toCheck=1 if name_num_sd>0 & !missing(name_num_sd)
		// result: 1779 (2178) with non-identical names
