// This file was written by Simon Luechinger (see Municipalities All Simon in Projekt Nationalräte\1 Data\1_data
// The file corrects old municipality names to the official municipality structure in 1975

replace gemeindename = "Aarwangen" if canton == "BE" & city == "Aarwangen-Hard" 
replace gemeindename = "Seegräben" if canton == "ZH" & city == "Aathal" // Aathal ist Teil von Seegräben.
replace gemeindename = "Seegräben" if canton == "ZH" & city == "Aathal-Seegräben" // Aathal ist Teil von Seegräben.
replace gemeindename = "Gaiserwald" if canton == "SG" & city == "Abtwil-Gaiserwald" // Abtwil ist Teil von Gaiserwald.
replace gemeindename = "Aesch (BL)" if canton == "BL" & city == "Aesch" 
replace gemeindename = "Aeschi bei Spiez" if canton == "BE" & city == "Aeschi" 
replace gemeindename = "Affoltern am Albis" if canton == "ZH" & city == "Affoltern a. Albis" 
replace gemeindename = "Affoltern am Albis" if canton == "ZH" & city == "Affoltern a.A." 
replace gemeindename = "Affoltern am Albis" if canton == "ZH" & city == "Affolterna.A." 
replace gemeindename = "Rubigen" if canton == "BE" & city == "Allmendingen" // Allmendingen war Teil von Rubigen.
replace gemeindename = "Altdorf (UR)" if canton == "UR" & city == "Altdorf" 
replace gemeindename = "Altstätten" if canton == "SG" & city == "Altstatten" 
replace gemeindename = "Grossandelfingen" if canton == "ZH" & city == "Andelfingen" // Früherer Name.
replace gemeindename = "Appenzell" if canton == "AI" & city == "Appenzell-Unterrhein" // Ortschaft bzw. Ortsteil nur in Verzeichnis aus 1818 gefunden. Eintrag in Bundesblatt bezieht sich auf Arnold Koller in 1975; 1979 ist dieser als in Appenzell wohnhaft vermerkt.
replace gemeindename = "Boudry" if canton == "NE" & city == "Areuse/ Boudry" 
replace gemeindename = "Boudry" if canton == "NE" & city == "Areuse/Boudry" 
replace gemeindename = "Arlesheim" if canton == "BL" & city == "Ariesheim" 
replace gemeindename = "Au (SG)" if canton == "SG" & city == "Au" 
replace gemeindename = "Wädenswil" if canton == "ZH" & city == "Au-Wädenswil" 
replace gemeindename = "Au (SG)" if canton == "SG" & city == "AuSG" 
replace gemeindename = "Lufingen" if canton == "ZH" & city == "Augwil/ Kloten" // Augwil ist Teil von Klotens Nachbargemeinde Lufingen.
replace gemeindename = "Montreux" if canton == "VD" & city == "Avants-sur-Montreux" 
replace gemeindename = "Basel" if canton == "BS" & city == "Base]" 
replace gemeindename = "Basel" if canton == "BS" & city == "Basel '" 
replace gemeindename = "Bassecourt" if canton == "BE" & city == "Berlincourt" // Berlincourt ist Teil von Bassecourt.
replace gemeindename = "Bern" if canton == "BE" & city == "Bern-Bethlehem" // Bethlehem ist Teil von Bern.
replace gemeindename = "Bern" if canton == "BE" & city == "Bern-Bümpliz" // Bümpliz ist Teil von Bern.
replace gemeindename = "Bern" if canton == "BE" & city == "Bern-Eymatt" // Eymatt ist Teil von Bern.
replace gemeindename = "Wohlen bei Bern" if canton == "BE" & city == "Bern-Hinterkappelen" // Hinterkappelen ist Teil von Wohlen bei Bern.
replace gemeindename = "Bern" if canton == "TI" & city == "Berna" // Kontrolle notwendig!
replace gemeindename = "Bern" if canton == "BE" & city == "Berne" 
replace gemeindename = "Bern" if canton == "FR" & city == "Berne" // Kontrolle notwendig!
replace gemeindename = "Bern" if canton == "NE" & city == "Berne" // Kontrolle notwendig!
replace gemeindename = "Bern" if canton == "VD" & city == "Berne" // Kontrolle notwendig!
replace gemeindename = "Bern" if canton == "VS" & city == "Berne" // Kontrolle notwendig!
replace gemeindename = "Biel (BE)" if canton == "BE" & city == "Biel" 
replace gemeindename = "Biel (BL)" if canton == "BL" & city == "Biel BL" 
replace gemeindename = "Biel (BE)" if canton == "BE" & city == "Bienne" 
replace gemeindename = "Walkringen" if canton == "BE" & city == "Bigenthal" // Bigenthal ist Teil von Walkringen.
replace gemeindename = "Nürensdorf" if canton == "ZH" & city == "Birchwil/Nürendorf" // Birchwil ist Teil von Nürensdorf.
replace gemeindename = "Birmensdorf (ZH)" if canton == "ZH" & city == "Birmensdorf" 
replace gemeindename = "Birmenstorf (AG)" if canton == "AG" & city == "Birmenstorf" 
replace gemeindename = "Les Bois" if canton == "BE" & city == "Bois" 
replace gemeindename = "Bolligen" if canton == "BE" & city == "Bolligen-Dorf" 
replace gemeindename = "Satigny" if canton == "GE" & city == "Bourdigny" // Bourdigny-Dessous ist Teil von Satigny.
replace gemeindename = "Satigny" if canton == "GE" & city == "Bourdigny-Satigny" // Bourdigny-Dessous ist Teil von Satigny.
replace gemeindename = "Sion" if canton == "VS" & city == "Bramois" // Fusionierte am 01.01.1968 mit Sion.
replace gemeindename = "Le Chenit" if canton == "VD" & city == "Brassus" // Le Brassus ist Teil von Le Chenit.
replace gemeindename = "Wynigen" if canton == "BE" & city == "Breitenegg bei Wynigen" // Breitenegg ist Teil von Wynigen.
replace gemeindename = "Bremgarten (AG)" if canton == "AG" & city == "Bremgarten" 
replace gemeindename = "Les Brenets" if canton == "NE" & city == "Brenets" 
replace gemeindename = "Les Breuleux" if canton == "BE" & city == "Breuleux" 
replace gemeindename = "Brienz (BE)" if canton == "BE" & city == "Brienz" 
replace gemeindename = "Brig" if canton == "VS" & city == "Brigue" 
replace gemeindename = "Raron" if canton == "VS" & city == "Brigue/Rarogne" // Rarogne ist frz. Name von Raron in der Nähe von Brig.
replace gemeindename = "Brügg" if canton == "BE" & city == "Brügg bei Biel" 
replace gemeindename = "Wangen (ZH)" if canton == "ZH" & city == "Brüttisellen" // Spätere Namensänderung zu Wangen-Brüttisellen.
replace gemeindename = "Luzein" if canton == "GR" & city == "Buchen i.P." // Buchen im Prättigau ist Teil von Luzein.
replace gemeindename = "Buchs (AG)" if canton == "AG" & city == "Buchs" 
replace gemeindename = "Buchs (SG)" if canton == "SG" & city == "Buchs" 
replace gemeindename = "Buchs (ZH)" if canton == "ZH" & city == "Buchs" 
replace gemeindename = "Buchs (SG)" if canton == "SG" & city == "Buchs SG" 
replace gemeindename = "Buchs (ZH)" if canton == "ZH" & city == "Buchs ZH" 
replace gemeindename = "Neukirch an der Thur" if canton == "TG" & city == "Buhl-Neukirch an der Thur" 
replace gemeindename = "Bünzen" if canton == "AG" & city == "Bunzen" // Umlaut in Bundesblatt teilweise nicht ersichtlich. Fehler jedoch ersichtlich bei Kandiderenden wie Herwig Abt aus "Bünzen" in einigen Jahren und aus "Bunzen" in anderen.
replace gemeindename = "Burg im Leimental" if canton == "BE" & city == "Burg i.L." // Burg im Leimental wechselte 1994 von BL zu BE.
replace gemeindename = "Bussigny-près-Lausanne" if canton == "VD" & city == "Bussigny s. Morges" // Früherer Name war Bussigny-sur-Morges.
replace gemeindename = "Busswil bei Büren" if canton == "BE" & city == "Busswil" // Kontrolle notwendig! 2 Busswil in BE (Busswil bei Büren und Busswil bei Melchnau). Vermutlich nur relevant bei Kandidat Thomas Rosa, Fensterfabrikant in Busswil bei Büren (vgl. Jubiläumsprospekt von ROSA Fenster: http://www.rosafenster.ch/userfiles/file/90-Jahre-Rosa-Fenster-Prospekt.pdf; 25.01.2017).
replace gemeindename = "Diemtigen" if canton == "BE" & city == "Bächlen/Oey-Diemtigen" 
replace gemeindename = "Le Bémont (BE)" if canton == "BE" & city == "Bémont" 
replace gemeindename = "Biel (BE)" if canton == "BE" & city == "Bözingen" // Bözingen ist Teil von Biel (BE).
replace gemeindename = "Bern" if canton == "BE" & city == "Bümpliz" // Bümpliz ist Teil von Bern.
replace gemeindename = "Büren an der Aare" if canton == "BE" & city == "Büren a. A." 
replace gemeindename = "Büren an der Aare" if canton == "BE" & city == "Büren a.A." 
replace gemeindename = "Büren an der Aare" if canton == "BE" & city == "Büren a.d.A." 
replace gemeindename = "Thunstetten" if canton == "BE" & city == "Bützberg" // Bützberg ist Teil von Thunstetten.
replace gemeindename = "Camorino" if canton == "TI" & city == "Camerino" 
replace gemeindename = "Carouge (GE)" if canton == "GE" & city == "Carouge" 
replace gemeindename = "Carouge (GE)" if canton == "GE" & city == "Carouge GE" 
replace gemeindename = "Carrouge (VD)" if canton == "VD" & city == "Carrouge VD" 
replace gemeindename = "Lausanne" if canton == "VD" & city == "Chailly sur Lausanne" // Chailly ist Teil von Lausanne.
replace gemeindename = "Gempenach" if canton == "FR" & city == "Champagny (Gempenach)" // Champagny ist frz. Name.
replace gemeindename = "Chapelle (Glâne)" if canton == "FR" & city == "Chapelle-sur-Gillarens" // Früherer Name.
replace gemeindename = "Chapelle (Glâne)" if canton == "FR" & city == "Chapelle-sur-Glâne" 
replace gemeindename = "Chavannes-sur-Moudon" if canton == "VD" & city == "Chavannes s/Moudon" 
replace gemeindename = "Chavannes-près-Renens" if canton == "VD" & city == "Chavannes-Renens" 
replace gemeindename = "Chavannes-près-Renens" if canton == "VD" & city == "Chavannes/Renens" 
replace gemeindename = "Chavannes-le-Veyron" if canton == "VD" & city == "Chavannesle-Veyron" 
replace gemeindename = "Le Chenit" if canton == "VD" & city == "Chenit" 
replace gemeindename = "Montreux" if canton == "VD" & city == "Chernex" // Chernex ist Teil von Montreux.
replace gemeindename = "Montreux" if canton == "VD" & city == "Chernex-sur-Montreux" // Chernex ist Teil von Montreux.
replace gemeindename = "Chesalles-sur-Moudon" if canton == "VD" & city == "Chesalles s. Moudon" 
replace gemeindename = "Cheseaux-Noréaz" if canton == "VD" & city == "Cheseaux" // Kontrolle notwendig! 2 Cheseaux in VD (Cheseaux-Noréaz und Cheseaux-sur-Lausanne). Cheseaux ist Weiler in Cheseaux-Noréaz, daher diese Gemeinde verwendet. Vermutlich nur relevant für Henri Cottier, Bauer.
replace gemeindename = "Ollon" if canton == "VD" & city == "Chesières" // Chesières ist Teil von Ollon.
replace gemeindename = "Anières" if canton == "GE" & city == "Chevrens s Anières" // Chevrens ist Teil von Anières.
replace gemeindename = "Anières" if canton == "GE" & city == "Chevrens-Anières" // Chevrens ist Teil von Anières.
replace gemeindename = "Giffers" if canton == "FR" & city == "Chevrilles (Giffers)" // Chevrilles ist frz. Name.
replace gemeindename = "Gorgier" if canton == "NE" & city == "Chez-le-Bart (commune de Gorgier)" 
replace gemeindename = "Kerzers" if canton == "FR" & city == "Chiètres" // Chiètres ist frz. Name.
replace gemeindename = "Vandoeuvres" if canton == "GE" & city == "Chougny-Vandœuvres" // Chougny ist Teil von Vandoeuvres.
replace gemeindename = "Château-d'Oex" if canton == "VD" & city == "Château d'Œx" 
replace gemeindename = "Château-d'Oex" if canton == "VD" & city == "Château-d'Œx" 
replace gemeindename = "Château-d'Oex" if canton == "VD" & city == "Châteaux-d'Oex" 
replace gemeindename = "Châtel-Saint-Denis" if canton == "FR" & city == "Châtel-St-Denis" 
replace gemeindename = "Montreux" if canton == "VD" & city == "Châtelard" // Kontrolle notwendig! Châtelard ohne Zusatz gibt es in Bundesblattt gar nicht; immer mit als Châtelard-Montreux oder Montreux-Châtelard. Eigenständige Gemeinde nur bis 1961.
replace gemeindename = "Montreux" if canton == "VD" & city == "Châtelard-Montreux" 
replace gemeindename = "Ollon" if canton == "VD" & city == "Chésières" // Chesières ist Teil von Ollon.
replace gemeindename = "Montreux" if canton == "VD" & city == "Clarens" // Clarens ist Teil von Montreux.
replace gemeindename = "Montreux" if canton == "VD" & city == "Clarens-Montreux" // Clarens ist Teil von Montreux.
replace gemeindename = "Montreux" if canton == "VD" & city == "Clarens/ Montreux" // Clarens ist Teil von Montreux.
replace gemeindename = "Montreux" if canton == "VD" & city == "Clarens/Montreux" // Clarens ist Teil von Montreux.
replace gemeindename = "Collonge-Bellerive" if canton == "GE" & city == "Collonge-Bellerive GE" 
replace gemeindename = "Combremont-le-Grand" if canton == "VD" & city == "Combremontle-Grand" 
replace gemeindename = "Genève" if canton == "GE" & city == "Conches (Genève)" // Conches ist Teil von Genève.
replace gemeindename = "Corcelles-Cormondrèche" if canton == "NE" & city == "Corcelles" 
replace gemeindename = "Corcelles-près-Payerne" if canton == "VD" & city == "Corcelles p Payerne" 
replace gemeindename = "Corcelles-près-Payerne" if canton == "VD" & city == "Corcelles p. Payerne" 
replace gemeindename = "Corcelles-près-Payerne" if canton == "VD" & city == "Corcelles sur Payerne" 
replace gemeindename = "Corcelles-Cormondrèche" if canton == "NE" & city == "Cormondrèche" 
replace gemeindename = "Corsier (GE)" if canton == "GE" & city == "Corsier" 
replace gemeindename = "Lutry" if canton == "VD" & city == "Corsy s Lutry" // Corsy ist Teil von Lutry.
replace gemeindename = "Lutry" if canton == "VD" & city == "Corsy s. Lutry" // Corsy ist Teil von Lutry.
replace gemeindename = "Courroux" if canton == "BE" & city == "Courcelon" // Courcelon ist Teil von Courroux.
replace gemeindename = "Courtételle" if canton == "BE" & city == "Courtetelle" 
replace gemeindename = "Montagny-les-Monts" if canton == "FR" & city == "Cousset" // Cousset war Teil von Montagny-les-Monts (heute Montagny); siehe: http://www.hls-dhs-dss.ch/textes/d/D826.php; 25.01.2017.
replace gemeindename = "Crans-près-Céligny" if canton == "VD" & city == "Crans" // Crans-près-Céligny scheint teilweise auch als Crans-sur-Nyon bezeichnet zu werden.
replace gemeindename = "Crans-près-Céligny" if canton == "VD" & city == "Crans s. Nyon" // Crans-près-Céligny scheint teilweise auch als Crans-sur-Nyon bezeichnet zu werden.
replace gemeindename = "Cressier (NE)" if canton == "NE" & city == "Cressier" 
replace gemeindename = "Crissier" if canton == "VD" & city == "Crisier" 
replace gemeindename = "Bardonnex" if canton == "GE" & city == "Croix-de-Rozon (Bardonnex)" // Croix-de-Rozon ist Teil von Bardonnex.
replace gemeindename = "Le Crêt" if canton == "FR" & city == "Crêt" 
replace gemeindename = "La Chaux-de-Fonds" if canton == "NE" & city == "Crêt-du-Locle" // Crêt-du-Locle ist Teil von La Chaux-de-Fonds.
replace gemeindename = "Genève" if canton == "GE" & city == "Crêts de Champel" // Crêts-de-Champel ist Teil von Genève.
replace gemeindename = "Cugy (VD)" if canton == "VD" & city == "Cugy" 
replace gemeindename = "Mettmenstetten" if canton == "ZH" & city == "Dachelsen/ Mettmenstetten" 
replace gemeindename = "Davos" if canton == "ZH" & city == "Davos Platz" // Kontrolle notwendig!
replace gemeindename = "Davos" if canton == "GR" & city == "Davos-Platz" 
replace gemeindename = "Stettlen" if canton == "BE" & city == "Deisswil-Stettlen" // Deisswil ist Teil von Stettlen.
replace gemeindename = "Unterschlatt" if canton == "TG" & city == "Dickihof" // Dickihof war Teil von Unterschlatt (dann Schlatt bei Diessenhofen, heute Schlatt).
replace gemeindename = "Bauma" if canton == "ZH" & city == "Dillhaus-Bauma" // Dillhaus ist Teil von Bauma.
replace gemeindename = "Dinhard" if canton == "ZH" & city == "Dinhard-Welsikon" // Welsikon ist Teil von Dinhard.
replace gemeindename = "Dompierre (FR)" if canton == "FR" & city == "Dompierre" 
replace gemeindename = "Dompierre (VD)" if canton == "VD" & city == "Dompierre" 
replace gemeindename = "Dübendorf" if canton == "ZH" & city == "Dubendorf" 
replace gemeindename = "Dürnten" if canton == "ZH" & city == "Durnten" 
replace gemeindename = "Stallikon" if canton == "ZH" & city == "Dägerst, Stallikon" // Tägerst ist Teil von Stallikon.
replace gemeindename = "Maur" if canton == "ZH" & city == "Ebmatingen" // Ebmatingen ist Teil von Maur.
replace gemeindename = "Ebnat-Kappel" if canton == "SG" & city == "Ebnat" // Bis 1965 eigenständige Gemeinde.
replace gemeindename = "Ecublens (VD)" if canton == "VD" & city == "Ecublens" 
replace gemeindename = "Ecublens (VD)" if canton == "VD" & city == "Ecublens (Vaud)" 
replace gemeindename = "Ecublens (VD)" if canton == "VD" & city == "Ecublens VD" 
replace gemeindename = "Ecublens (VD)" if canton == "VD" & city == "Ecublèns" 
replace gemeindename = "Illnau" if canton == "ZH" & city == "Effretikon" // Früherer Name für Illnau-Effretikon.
replace gemeindename = "Egg" if canton == "ZH" & city == "Egg ZH" 
replace gemeindename = "Embrach" if canton == "ZH" & city == "Einbrach" 
replace gemeindename = "Ellikon an der Thur" if canton == "ZH" & city == "Ellikon a. d. Thur" 
replace gemeindename = "Ellikon an der Thur" if canton == "ZH" & city == "Ellikona.d. Thur" 
replace gemeindename = "Emmen" if canton == "LU" & city == "Emmenbrücke" // Emmenbrücke ist Teil von Emmen.
replace gemeindename = "Lauperswil" if canton == "BE" & city == "Emmenmatt" // Emmenmatt ist Teil von Lauperswil.
replace gemeindename = "Lauperswil" if canton == "BE" & city == "Emmenmatt im Emmental" // Emmenmatt ist Teil von Lauperswil.
replace gemeindename = "Worb" if canton == "BE" & city == "Enggistein" // Enggistein ist Teil von Worb.
replace gemeindename = "Untersiggenthal" if canton == "AG" & city == "Ennetturgi Gde. Untersiggenthal" 
replace gemeindename = "Epesses" if canton == "VD" & city == "Epésses" 
replace gemeindename = "Erlenbach (ZH)" if canton == "ZH" & city == "Erlenbach" 
replace gemeindename = "Essert-sous-Champvent" if canton == "VD" & city == "Esseri/Champvent" 
replace gemeindename = "Essert (FR)" if canton == "FR" & city == "Essert" 
replace gemeindename = "Essert-sous-Champvent" if canton == "VD" & city == "Essert s. Champvent" 
replace gemeindename = "Essert-sous-Champvent" if canton == "VD" & city == "Essert/ Champvent" 
replace gemeindename = "Essert-sous-Champvent" if canton == "VD" & city == "Essert/Champvent" 
replace gemeindename = "Essertines-sur-Yverdon" if canton == "VD" & city == "Essertines sur Yverdon" 
replace gemeindename = "Essert-sous-Champvent" if canton == "VD" & city == "Essertsous-Champvent" 
replace gemeindename = "Egg" if canton == "ZH" & city == "Esslingen ZH" // Esslingen ist Teil von Egg.
replace gemeindename = "Egg" if canton == "ZH" & city == "Esslingen-Egg" // Esslingen ist Teil von Egg.
replace gemeindename = "Estavayer-le-Lac" if canton == "FR" & city == "Estavayer" // Kontrolle notwendig! 2 Estavayer in FR (Estavayer-le-Gibloux und Estavayer-le-Lac). Betrifft vernutlich nur Armand Droz, in anderen Jahren in Estavayer-le-Lac wohnhaft.
replace gemeindename = "Estavayer-le-Lac" if canton == "FR" & city == "Estavayer le-Lac" 
replace gemeindename = "Estavayer-le-Lac" if canton == "FR" & city == "Estavayerle-Lac" 
replace gemeindename = "Wetzikon (ZH)" if canton == "ZH" & city == "Ettenhausen/ Wetzikon" // Ettenhausen ist Teil von Wetzikon.
replace gemeindename = "Farvagny-le-Grand" if canton == "FR" & city == "Farvagny-leGrand" 
replace gemeindename = "Affoltern am Albis" if canton == "ZH" & city == "Fehrenbach-Affoltern a.A." // Fehrenbach ist Teil von Affoltern am Albis.
replace gemeindename = "Meilen" if canton == "ZH" & city == "Feldmeilen" // Feldmeilen ist Teil von Meilen.
replace gemeindename = "Bolligen" if canton == "BE" & city == "Ferenberg" // Ferenberg ist Teil von Bolligen.
replace gemeindename = "Ferlens (VD)" if canton == "VD" & city == "Ferlens" 
replace gemeindename = "Gsteig" if canton == "BE" & city == "Feutersoey bei Gstaad" // Feutersoey ist Teil von Gsteig.
replace gemeindename = "Wünnewil" if canton == "FR" & city == "Flamatt" // Flamatt ist Teil von Wünnewil-Flamatt (bis 1974 nur Wünnewil).
replace gemeindename = "Fontaines-sur-Grandson" if canton == "VD" & city == "Fontaines" 
replace gemeindename = "Maur" if canton == "ZH" & city == "Forch/Maur" 
replace gemeindename = "Seedorf (BE)" if canton == "BE" & city == "Frienisberg" // Frienisberg ist Teil von Seedorf.
replace gemeindename = "Frutigen" if canton == "BE" & city == "Frutìgen" 
replace gemeindename = "Galmiz" if canton == "FR" & city == "Galmiz/Guin" // Kontrolle notwendig! Galmiz und Guin (Düdingen) sind eigentlich 2 versch. Gemeinden.
replace gemeindename = "Gams" if canton == "SG" & city == "Garns" 
replace gemeindename = "Les Genevez (BE)" if canton == "BE" & city == "Genevez" 
replace gemeindename = "Genève" if canton == "ZH" & city == "Genf" // Kontrolle notwendig!
replace gemeindename = "Genève" if canton == "GE" & city == "Genèver" 
replace gemeindename = "Täuffelen" if canton == "BE" & city == "Gerolfingen" // Gerolfingen ist Teil von Täuffelen.
replace gemeindename = "Fischenthal" if canton == "ZH" & city == "Gibswil" // Gibswil ist Teil von Fischenthal.
replace gemeindename = "Opfikon" if canton == "ZH" & city == "Glattbrugg" // Glattbrugg ist Teil von Opfikon.
replace gemeindename = "Opfikon" if canton == "ZH" & city == "Glattbrugg/Opfikon" // Glattbrugg ist Teil von Opfikon.
replace gemeindename = "Montreux" if canton == "VD" & city == "Glion/Montreux" // Glion ist Teil von Montreux.
replace gemeindename = "Arth" if canton == "SZ" & city == "Goldau" // Goldau ist Teil von Arth.
replace gemeindename = "Thun" if canton == "BE" & city == "Goldiwil" // Goldiwil ist Teil von Thun.
replace gemeindename = "Gossau (SG)" if canton == "SG" & city == "Gossau" 
replace gemeindename = "Gossau (ZH)" if canton == "ZH" & city == "Gossau" 
replace gemeindename = "Gossau (SG)" if canton == "SG" & city == "Gossau (St. Gallen)" 
replace gemeindename = "Lancy" if canton == "GE" & city == "Grand-Lancy" // Grand-Lancy ist Teil von Lancy.
replace gemeindename = "Le Grand-Saconnex" if canton == "GE" & city == "Grand-Saconnex" 
replace gemeindename = "Granges (VS)" if canton == "VS" & city == "Granges" 
replace gemeindename = "Seeberg" if canton == "BE" & city == "Grasswil" // Grasswil ist Teil von Seeberg.
replace gemeindename = "Grellingen" if canton == "BE" & city == "Grellingue" // Grellingue ist frz. von Grellingen.
replace gemeindename = "Saanen" if canton == "BE" & city == "Grund bei Gstaad" // Grund bei Gstaad ist Teil von Saanen.
replace gemeindename = "Sumiswald" if canton == "BE" & city == "Grünen-Sumiswald" // Grünen ist Teil von Sumiswald.
replace gemeindename = "Lützelflüh" if canton == "BE" & city == "Grünenmatt" // Grünenmatt ist Teil von Lützelflüh.
replace gemeindename = "Gossau (ZH)" if canton == "ZH" & city == "Grüt b. Wetzikon" // Grüt ist Teil von Gossau.
replace gemeindename = "Gossau (ZH)" if canton == "ZH" & city == "Grüt/Gossau" // Grüt ist Teil von Gossau.
replace gemeindename = "Saanen" if canton == "BE" & city == "Gstaad" // Gstaad ist Teil von Saanen.
replace gemeindename = "Gsteig" if canton == "BE" & city == "Gsteig bei Gstaad" 
replace gemeindename = "Düdingen" if canton == "FR" & city == "Guin" // Guin ist frz. Name von Düdingen.
replace gemeindename = "Sigriswil" if canton == "BE" & city == "Gunten" // Gunten ist Teil von Sigriswil.
replace gemeindename = "Guntershausen bei Aadorf" if canton == "TG" & city == "Guntershausen" // Kontrolle notwendig! 2 Guntershausen in TG (Guntershausen bei Aadorf und Guntershausen bei Birwinken); betrifft nur Benedikt Beer aus vermutlich Guntershausen bei Aadord; vgl. http://www.aadorf.ch/xml_1/internet/de/application/d1/d190/f193.cfm; 25.01.2017.
replace gemeindename = "Volketswil" if canton == "ZH" & city == "Gutenswil" // Gutenwil ist Teil von Volketswil.
replace gemeindename = "Muri bei Bern" if canton == "BE" & city == "Gümligen" // Gümligen ist Teil von Muri bei Bern.
replace gemeindename = "Felben" if canton == "TG" & city == "Haldenhof-Felben" 
replace gemeindename = "Erlinsbach" if canton == "AG" & city == "Hard-Erlinsbach" // Hard ist Teil von Erlinsbach.
replace gemeindename = "Hasle bei Burgdorf" if canton == "BE" & city == "Hasle" 
replace gemeindename = "Hasle (LU)" if canton == "LU" & city == "Hasle" 
replace gemeindename = "Hasle bei Burgdorf" if canton == "BE" & city == "Hasle b/Burgdorf" 
replace gemeindename = "Hasle bei Burgdorf" if canton == "BE" & city == "Hasle-Rüegsau" // Hasle bei Burgdorf und Rüegsau sind eigentlich 2 versch. Gemeinden.
replace gemeindename = "Hausen am Albis" if canton == "ZH" & city == "Hausen a.A." 
replace gemeindename = "Volketswil" if canton == "ZH" & city == "Hegnau" // Hegnau ist Teil von Volketswil.
replace gemeindename = "Volketswil" if canton == "ZH" & city == "Hegnau-Volketswil" // Hegnau ist Teil von Volketswil.
replace gemeindename = "Volketswil" if canton == "ZH" & city == "Hegnau/ Volketswil" // Hegnau ist Teil von Volketswil.
replace gemeindename = "Neuenegg" if canton == "BE" & city == "Heitern-Neuenegg" // Heitern ist Teil von Neuenegg.
replace gemeindename = "Hütten" if canton == "ZH" & city == "Hengerten-Hutten" // Hengerten ist Teil von Hütten.
replace gemeindename = "Rüeggisberg" if canton == "BE" & city == "Hinterfultigen" // Hinterfultigen ist Teil von Rüeggisberg.
replace gemeindename = "Wohlen bei Bern" if canton == "BE" & city == "Hinterkappelen" // Hinterkappelen ist Teil von Wohlen bei Bern.
replace gemeindename = "Holderbank (AG)" if canton == "AG" & city == "Holderbank" 
replace gemeindename = "Münchringen" if canton == "BE" & city == "Holzmühle, Münchringen" // Holzmühle ist Teil von Münchringen.
replace gemeindename = "Spiez" if canton == "BE" & city == "Hondrich" // Hondrich ist Teil von Spiez.
replace gemeindename = "Roggwil (TG)" if canton == "TG" & city == "Häuslen-Roggwil" // Häuslen ist Teil von Roggwil.
replace gemeindename = "Hilterfingen" if canton == "BE" & city == "Hünibach" // Hünibach ist Teil von Hilterfingen.
replace gemeindename = "Hilterfingen" if canton == "BE" & city == "Hünibach bei Thun" // Hünibach ist Teil von Hilterfingen.
replace gemeindename = "Hilterfingen" if canton == "BE" & city == "Hünibach/Thun" // Hünibach ist Teil von Hilterfingen.
replace gemeindename = "Schwyz" if canton == "SZ" & city == "Ibach" // Ibach ist Teil von Schwyz.
replace gemeindename = "Wohlen bei Bern" if canton == "BE" & city == "Innerberg bei Säriswil" // Innerberg und Säriswil sind Teil  von Wohlen bei Bern.
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "Irgenhausen" // Irgenhausen ist Teil von Pfäffikon.
replace gemeindename = "Bolligen" if canton == "BE" & city == "Ittigen" // Ittigen gehörte zu Bolligen.
replace gemeindename = "Jona" if canton == "SG" & city == "Jona SG" 
replace gemeindename = "Ebnat-Kappel" if canton == "SG" & city == "Kappel" // Bis 1965 eigenständige Gemeinde.
replace gemeindename = "Lindau" if canton == "ZH" & city == "Kempttal" // Kemptthal ist Teil von Lindau.
replace gemeindename = "Kilchberg (ZH)" if canton == "ZH" & city == "Kilchberg" 
replace gemeindename = "Kilchberg (ZH)" if canton == "ZH" & city == "Kilchberg ZH" 
replace gemeindename = "Volketswil" if canton == "ZH" & city == "Kindhausen-Volketswil" // Kindhausen ist Teil von Volketswil.
replace gemeindename = "Kirchberg (BE)" if canton == "BE" & city == "Kirchberg" 
replace gemeindename = "Kirchberg (BE)" if canton == "BE" & city == "Kirchberg BE" 
replace gemeindename = "Kirchdorf (BE)" if canton == "BE" & city == "Kirchdorf" 
replace gemeindename = "Böttstein" if canton == "AG" & city == "Kleindöttingen-Böttstein" // Kleindöttingen ist Teil von Böttstein.
replace gemeindename = "Zell (ZH)" if canton == "ZH" & city == "Kollbrunn-Zell" // Kollbrunn ist Teil von Zell.
replace gemeindename = "Konolfingen" if canton == "BE" & city == "Konolflngen" 
replace gemeindename = "Köniz" if canton == "VD" & city == "Köniz (Berne)" // Kontrolle notwendig!
replace gemeindename = "Köniz" if canton == "BE" & city == "Köniz-Liebefeld" // Liebefeld ist Teil von Köniz.
replace gemeindename = "Küsnacht (ZH)" if canton == "ZH" & city == "Küsnacht" 
replace gemeindename = "Küsnacht (ZH)" if canton == "ZH" & city == "Küsnacht/Itschnach" // Itschnach ist Teil von Küsnacht.
replace gemeindename = "Küssnacht am Rigi" if canton == "SZ" & city == "Küssnacht" 
replace gemeindename = "Küttigen" if canton == "AG" & city == "Küttigen-Rombach" // Rombach ist Teil von Küttigen.
replace gemeindename = "Château-d'Oex" if canton == "VD" & city == "L'Etivaz" // L'Etivaz ist Teil von Château-d'Oex.
replace gemeindename = "Le Chenit" if canton == "VD" & city == "L'Orient" // L'Orient ist Teil von Le Chenit.
replace gemeindename = "Martigny" if canton == "VS" & city == "La Bâtiaz" // La Bâtiaz ist Teil von Martigny.
replace gemeindename = "La Chaux (Cossonay)" if canton == "VD" & city == "La Chaux" 
replace gemeindename = "La Chaux (Cossonay)" if canton == "VD" & city == "La Chaux près Cossonay" 
replace gemeindename = "La Chaux (Cossonay)" if canton == "VD" & city == "La Chaux s Cossonay" 
replace gemeindename = "Sonvilier" if canton == "BE" & city == "La Chaux-d'Abel" // La Chaux-d'Abel ist Teil von Sonvilier.
replace gemeindename = "La Chaux-de-Fonds" if canton == "NE" & city == "La Chaux-deFonds" 
replace gemeindename = "La Chaux-de-Fonds" if canton == "NE" & city == "La Chaux-déFonds" 
replace gemeindename = "La Chaux-de-Fonds" if canton == "NE" & city == "La Chauxde-Fonds" 
replace gemeindename = "Saint-Légier-La Chiésaz" if canton == "VD" & city == "La Chiésaz" // La Chiésaz ist Teil von Saint-Légier-La Chiésaz.
replace gemeindename = "Bardonnex" if canton == "GE" & city == "La Croix-de-Rozon" // La Croix-de-Rozon ist Teil von Bardonnex.
replace gemeindename = "Les Clées" if canton == "VD" & city == "La Russille" // La Russille ist Teil von Les Clées.
replace gemeindename = "La Tour-de-Peilz" if canton == "VD" & city == "La Tour-dePeilz" 
replace gemeindename = "Lajoux (BE)" if canton == "BE" & city == "Lajoux" 
replace gemeindename = "Langnau im Emmental" if canton == "BE" & city == "Langnau" 
replace gemeindename = "Langnau am Albis" if canton == "ZH" & city == "Langnau a.A." 
replace gemeindename = "Langnau im Emmental" if canton == "BE" & city == "Langnau i. E." 
replace gemeindename = "Langnau im Emmental" if canton == "BE" & city == "Langnau i.E." 
replace gemeindename = "Langnau im Emmental" if canton == "BE" & city == "Langnaui.E." 
replace gemeindename = "Läufelfingen" if canton == "BL" & city == "Laufelfingen" 
replace gemeindename = "Laufen" if canton == "BE" & city == "Laufon" // Laufon ist frz. Name von Laufen (ab 1994 in BL).
replace gemeindename = "Lavey-Morcles" if canton == "VD" & city == "Lavey-Village" // Lavey-Village ist Teil von Lavey-Morcles.
replace gemeindename = "Le Mont-sur-Lausanne" if canton == "VD" & city == "Le Mont s/Lausanne" 
replace gemeindename = "L'Abbaye" if canton == "VD" & city == "Le Pont" // Le Pont ist Teil von L'Abbaye.
replace gemeindename = "Lengnau (AG)" if canton == "AG" & city == "Lengnau" 
replace gemeindename = "Lengnau (BE)" if canton == "BE" & city == "Lengnau" 
replace gemeindename = "Vaz/Obervaz" if canton == "GR" & city == "Lenzerheide/ Lai" // Lai ist rm. Name von Lenzerheide; Lenzerheide ist Teil von Vaz/Obervaz.
replace gemeindename = "Vaz/Obervaz" if canton == "GR" & city == "Lenzerheide/Lai" // Lai ist rm. Name von Lenzerheide; Lenzerheide ist Teil von Vaz/Obervaz.
replace gemeindename = "Ormont-Dessus" if canton == "VD" & city == "Les Diablerets" // Les Diablerets ist Teil von Ormont-Dessus.
replace gemeindename = "Château-d'Oex" if canton == "VD" & city == "Les Moulins" // Les Moulins ist Teil von Château-d'Oex.
replace gemeindename = "Evilard" if canton == "BE" & city == "Leubringen" // Leubringen ist dt. Name von Evilard.
replace gemeindename = "Leuk" if canton == "VS" & city == "Leuk-Stadt" 
replace gemeindename = "Köniz" if canton == "BE" & city == "Liebefeld" // Liebefeld ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "Liebefeld-Köniz" // Liebefeld ist Teil von Köniz.
replace gemeindename = "Locarno" if canton == "TI" & city == "Locamo" 
replace gemeindename = "Le Locle" if canton == "NE" & city == "Locle" 
replace gemeindename = "Lurtigen" if canton == "FR" & city == "Lourtens" // Lourtens ist frz. Name für Lurtigen.
replace gemeindename = "Bagnes" if canton == "VS" & city == "Lourtier-Bagnes" // Lourtier ist Teil von Bagnes.
replace gemeindename = "Leuk" if canton == "VS" & city == "Louèche-Ville" // Loèche ist frz. Name von Leuk.
replace gemeindename = "Leuk" if canton == "VS" & city == "Loèche" 
replace gemeindename = "Leuk" if canton == "VS" & city == "Loèche Ville" 
replace gemeindename = "Leuk" if canton == "VS" & city == "Loèche-Souste" // Souste ist Teil von Leuk.
replace gemeindename = "Eischoll" if canton == "VS" & city == "Loèche/Eischoll" // Eischoll ist eigenständige Gemeinde im Bezirk Leuk.
replace gemeindename = "Lully (FR)" if canton == "FR" & city == "Lully" 
replace gemeindename = "Lutry" if canton == "VD" & city == "Lutry " 
replace gemeindename = "Häggenschwil" if canton == "SG" & city == "Lömmenschwil-Häggenschwil" // Lömmenschwil ist Teil von Häggenschwil.
replace gemeindename = "Evilard" if canton == "BE" & city == "Magglingen" // Magglingen ist Teil von Evilard.
replace gemeindename = "Marbach (SG)" if canton == "SG" & city == "Marbach" 
replace gemeindename = "Martigny" if canton == "VS" & city == "Martigny-Bourg" // Martigny-Bourg war eigenständige Gemeinde bis 1993.
replace gemeindename = "Martigny" if canton == "VS" & city == "Martigny-Ville" // Martigny-Ville war eigenständige Gemeinde bis 1993.
replace gemeindename = "Matten bei Interlaken" if canton == "BE" & city == "Matten b. Interlaken" 
replace gemeindename = "Mies" if canton == "GE" & city == "Mies (Vaud)" // Kontrolle notwendig!
replace gemeindename = "Köniz" if canton == "BE" & city == "Mittelhäusern" // Mittelhäusern ist Teil von Köniz.
replace gemeindename = "Saint-Ursanne" if canton == "BE" & city == "Monnat, St-Ursanne" // Monnat ist vermutlich Teil von Saint-Ursanne.
replace gemeindename = "Saint-Imier" if canton == "BE" & city == "Mont-Soleil" // Mont-Soleil ist vermutlich Teil von Saint-Imier.
replace gemeindename = "Montagny-près-Yverdon" if canton == "VD" & city == "Montagny" 
replace gemeindename = "Montagny-près-Yverdon" if canton == "VD" & city == "Montagny près Yverdon" 
replace gemeindename = "Moutier" if canton == "BE" & city == "Montier" // Kontrolle notwendig! Vermutlich Schreibfehler.
replace gemeindename = "Montreux" if canton == "VD" & city == "Montreux-Châtelard" // Eigenständige Gemeinde nur bis 1961.
replace gemeindename = "Murten" if canton == "FR" & city == "Morat" // Morat ist frz. Name für Murten.
replace gemeindename = "Quarten" if canton == "SG" & city == "Murg-Quarten" // Murg ist Teil von Quarten.
replace gemeindename = "Muri (AG)" if canton == "AG" & city == "Muri" 
replace gemeindename = "Muri bei Bern" if canton == "BE" & city == "Muri" 
replace gemeindename = "Muri bei Bern" if canton == "BE" & city == "Muri BE" 
replace gemeindename = "Muri bei Bern" if canton == "GE" & city == "Muri BE" // Kontrolle notwendig!
replace gemeindename = "Aristau" if canton == "AG" & city == "Murimoos (Gemeinde Aristau)" 
replace gemeindename = "Biel (BE)" if canton == "BE" & city == "Mâche-Bienne" // Mâche ist Teil von Bienne.
replace gemeindename = "Mézières (VD)" if canton == "VD" & city == "Mézières" 
replace gemeindename = "Mézières (VD)" if canton == "VD" & city == "Mézières (Vaud)" 
replace gemeindename = "Môtiers (NE)" if canton == "NE" & city == "Môtiers" 
replace gemeindename = "Mühledorf (SO)" if canton == "SO" & city == "Mühledorf" 
replace gemeindename = "Weisslingen" if canton == "ZH" & city == "Neschwil/Weisslingen" // Neschwil ist Teil von Weisslingen.
replace gemeindename = "Allschwil" if canton == "BL" & city == "Neu-Allschwil" // Neuallschwil ist vermutlich Teil von Allschwil.
replace gemeindename = "Neuchâtel" if canton == "NE" & city == "Neuchatel" 
replace gemeindename = "Rüegsau" if canton == "BE" & city == "Neuegg/Sumiswald" // Neuegg ist vermutlich Teil von Rüegsau, Nachbargemeinde von Sumiswald.
replace gemeindename = "Neuhausen am Rheinfall" if canton == "SH" & city == "Neuhausen a.Rhf." 
replace gemeindename = "Egnach" if canton == "TG" & city == "Neukirch-Egnach" // Neukirch ist Teil von Egnach.
replace gemeindename = "Neukirch an der Thur" if canton == "TG" & city == "Neukirch/ Thur" 
replace gemeindename = "La Neuveville" if canton == "BE" & city == "Neuveville" 
replace gemeindename = "Niederdorf" if canton == "BL" & city == "Niederdorf BL" 
replace gemeindename = "Niederried bei Kallnach" if canton == "BE" & city == "Niederried b. Kallnach" 
replace gemeindename = "Köniz" if canton == "BE" & city == "Niederscherli" // Niederscherli ist Teil von Köniz.
replace gemeindename = "Uzwil" if canton == "SG" & city == "Niederstetten-Henau" // Niederstetten und Henau sind Teil von Uzwil.
replace gemeindename = "Uzwil" if canton == "SG" & city == "Niederstetten-Uzwil" // Niederstetten ist Teil von Uzwil.
replace gemeindename = "Uzwil" if canton == "SG" & city == "Niederuzwil-Henau" // Niederuzwi und Henau sind Teil von Uzwil.
replace gemeindename = "Uzwil" if canton == "SG" & city == "Niederuzwil-Uzwil" // Niederuzwi ist Teil von Uzwil.
replace gemeindename = "Yvonand" if canton == "VD" & city == "Niédens" // Niédens-Dessous und Niédens-Dessus sind Teil von Yvonand.
replace gemeindename = "Yvonand" if canton == "VD" & city == "Niédens/Yvonand" // Niédens-Dessous und Niédens-Dessus sind Teil von Yvonand.
replace gemeindename = "Moutier" if canton == "BE" & city == "Nods et Moutier" // Kontrolle notwendig! Nods und Moutier sind 2 eigenständige Gemeinden. Betrifft vermutlich nur Francis Erard. Moutier zugeordnet.
replace gemeindename = "Le Noirmont" if canton == "BE" & city == "Noirmont" 
replace gemeindename = "Uster" if canton == "ZH" & city == "Nänikon-Uster" // Nänikon ist Teil von Uster.
replace gemeindename = "Wetzikon (ZH)" if canton == "ZH" & city == "Ober-Wetzikon" // Oberwetzikon ist Teil von Wetzikon.
replace gemeindename = "Wohlen bei Bern" if canton == "BE" & city == "Oberdettigen" // Oberdettigen ist Teil von Wohlen bei Bern.
replace gemeindename = "Wohlen bei Bern" if canton == "BE" & city == "Oberdettigen b/Wohlen" // Oberdettigen ist Teil von Wohlen bei Bern.
replace gemeindename = "Oberdorf (BL)" if canton == "BL" & city == "Oberdorf" 
replace gemeindename = "Oberhofen am Thunersee" if canton == "BE" & city == "Oberhofen" 
replace gemeindename = "Schlosswil" if canton == "BE" & city == "Oberhünigen" // Oberhünigen geörte bis 1979 zu Schlosswil.
replace gemeindename = "Schlosswil" if canton == "BE" & city == "Oberhünigen/Niederhünigen" // Kontrolle notwendig! Oberhünigen waren nie in selber Gemeinde. Jedoch Oberhünigen zugeordnet, da Johann Geissbühler in anderem Jahr in Oberhünigen. Oberhünigen geörte bis 1979 zu Schlosswil.
replace gemeindename = "Meilen" if canton == "ZH" & city == "Obermeilen" // Obermeilen ist Teil von Meilen.
replace gemeindename = "Oberriet (SG)" if canton == "SG" & city == "Oberriet" 
replace gemeindename = "Köniz" if canton == "BE" & city == "Oberscherli" // Oberscherli ist Teil von Köniz.
replace gemeindename = "Steinmaur" if canton == "ZH" & city == "Obersteinmaur" // Obersteinmaur ist Teil von Steinmaur.
replace gemeindename = "Köniz" if canton == "BE" & city == "Oberwangen" // Oberwangen ist Teil von Köniz .
replace gemeindename = "Köniz" if canton == "BE" & city == "Oberwangen BE" // Oberwangen ist Teil von Köniz .
replace gemeindename = "Oberwil (BL)" if canton == "BL" & city == "Oberwil" 
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "Oberwil-Pfäffikon" // Oberwil ist Teil von Pfäffikon.
replace gemeindename = "Oetwil am See" if canton == "ZH" & city == "Oetwil a. See" 
replace gemeindename = "Oetwil am See" if canton == "ZH" & city == "Oetwil a.S." 
replace gemeindename = "Oetwil an der Limmat" if canton == "ZH" & city == "Oetwil a.d.L." 
replace gemeindename = "Diemtigen" if canton == "BE" & city == "Oey-Diemtigen" // Oey ist Teil von Diemtigen.
replace gemeindename = "Onnens (VD)" if canton == "VD" & city == "Onnens" 
replace gemeindename = "Ormont-Dessus" if canton == "VD" & city == "Ormont-dessus" 
replace gemeindename = "Oberburg" if canton == "BE" & city == "Oschwand" // Oschwand ist Teil von Oberburg.
replace gemeindename = "Bolligen" if canton == "BE" & city == "Ostermundigen" // Ostermundigen gehörte zu Bolligen.
replace gemeindename = "Oulens-sur-Lucens" if canton == "VD" & city == "Oulens/Lucens" 
replace gemeindename = "Bolligen" if canton == "BE" & city == "Papiermühle BE" // Papiermühle ist Teil von Ittigen, Ittigen gehörte zu Bolligen.
replace gemeindename = "Buchrain" if canton == "LU" & city == "Perlen" // Perlen ist Teil von Buchrain.
replace gemeindename = "Lancy" if canton == "GE" & city == "Petit-Lancy" // Petit-Lancy ist Teil von Lancy.
replace gemeindename = "Genève" if canton == "GE" & city == "Petit-Saconnex" // Petit-Saconnex ist Teil von Genève.
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "Pfäffikon ZH" 
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "PfäffikonZH" 
replace gemeindename = "Port" if canton == "BE" & city == "Port bei Nidau" // Port ist Nachbargemeinde von Nidau.
replace gemeindename = "Port" if canton == "BE" & city == "Port/Nidau" // Port ist Nachbargemeinde von Nidau.
replace gemeindename = "Vully-le-Bas" if canton == "FR" & city == "Praz (Vully)" // Praz war Teil von Vully-le-Bas (heute Mont-Vully).
replace gemeindename = "Porrentruy" if canton == "BE" & city == "Pruntrut" // Pruntrut ist dt. Name für Porrentruy.
replace gemeindename = "Le Pâquier (FR)" if canton == "FR" & city == "Pâquier" 
replace gemeindename = "Le Pâquier (NE)" if canton == "NE" & city == "Pâquier" 
replace gemeindename = "Langnau im Emmental" if canton == "BE" & city == "Raingutberg/Langnau im Emmental" // Raingutberg ist vermutlich Teil von Langnau im Emmental.
replace gemeindename = "Rapperswil (SG)" if canton == "SG" & city == "Rapperswil" 
replace gemeindename = "Rapperswil (BE)" if canton == "BE" & city == "Rapperswil BE" 
replace gemeindename = "Rapperswil (SG)" if canton == "SG" & city == "Rapperswil SG" 
replace gemeindename = "Raron" if canton == "VS" & city == "Rarogne" // Rarogne ist frz. Name von Raron.
replace gemeindename = "Reckingen (VS)" if canton == "VS" & city == "Reckingen" 
replace gemeindename = "Meiringen" if canton == "BE" & city == "Reichenbach bei Meiringen" // Reichenbach ist vermutlich Teil von Meiringen.
replace gemeindename = "Reinach (AG)" if canton == "AG" & city == "Reinach" 
replace gemeindename = "Reinach (BL)" if canton == "BL" & city == "Reinach" 
replace gemeindename = "Reinach (BL)" if canton == "BL" & city == "Reinach BL" 
replace gemeindename = "Renan (BE)" if canton == "BE" & city == "Renan" 
replace gemeindename = "Renens (VD)" if canton == "VD" & city == "Renens" 
replace gemeindename = "Renens (VD)" if canton == "VD" & city == "Renens " 
replace gemeindename = "Renens (VD)" if canton == "VD" & city == "Renens VD" 
replace gemeindename = "Luzern" if canton == "LU" & city == "Reussbühl" // Reussbühl ist Teil von Luzern.
replace gemeindename = "Ursenbach" if canton == "BE" & city == "Richisberg-Ursenbach" // Richisberg ist Teil von Ursenbach.
replace gemeindename = "Rickenbach (SO)" if canton == "SO" & city == "Rickenbach" 
replace gemeindename = "Schwyz" if canton == "SZ" & city == "Rickenbach" // Rickenbach ist Teil von Schwyz.
replace gemeindename = "Rickenbach bei Wil" if canton == "TG" & city == "Rickenbach" 
replace gemeindename = "Ried bei Kerzers" if canton == "FR" & city == "Ried près Chiètres" // Chiètres ist frz. Name für Kerzers.
replace gemeindename = "Riedholz" if canton == "SO" & city == "Riedholz-Attisholz" // Attisholz ist Teil von Riedholz.
replace gemeindename = "Seeberg" if canton == "BE" & city == "Riedtwil" // Riedtwil  ist Teil von Seeberg.
replace gemeindename = "Roche-d'Or" if canton == "BE" & city == "Roche d'Or" 
replace gemeindename = "Roggwil (BE)" if canton == "BE" & city == "Roggwil" 
replace gemeindename = "Roggwil (BE)" if canton == "BE" & city == "Roggwil BE" 
replace gemeindename = "Rohr (AG)" if canton == "AG" & city == "Rohr" 
replace gemeindename = "Romanel-sur-Lausanne" if canton == "VD" & city == "Romanel" // Kontrolle notwendig! 2 Romanel in VD (Romanel-sur-Lausanne und Romanel-sur-Morges), jedoch: Bei Wohnort nur Romanel-sur-Lausanne in Bundesblatt gefunden.
replace gemeindename = "Romanel-sur-Lausanne" if canton == "VD" & city == "Romanel près Lausanne" 
replace gemeindename = "Romont (FR)" if canton == "FR" & city == "Romont" 
replace gemeindename = "Mühleberg" if canton == "BE" & city == "Rosshäusern" // Rosshäusern ist Teil von Mühleberg.
replace gemeindename = "Roveredo (GR)" if canton == "GR" & city == "Roveredo" 
replace gemeindename = "Rudolfstetten-Friedlisberg" if canton == "AG" & city == "Rudolfstetten" // Rudolfstetten ist Teil von Rudolfstetten-Friedlisberg.
replace gemeindename = "Worb" if canton == "BE" & city == "Rufenacht-Worb" // Rüfenacht ist Teil von Worb.
replace gemeindename = "Rüschlikon" if canton == "ZH" & city == "Ruschlikon" 
replace gemeindename = "Röthenbach im Emmental" if canton == "BE" & city == "Röthenbach" // Kontrolle notwendig! 2 Röthenbach in BE (Röthenbach im Emmental und Röthenbach bei Herzogenbuchsee). Vermutlich nur relevant bei Kandidat Rudolf Rüegsegger, zu anderem Zeitpunkt in Röthenbach im Emmental.
replace gemeindename = "Röthenbach im Emmental" if canton == "BE" & city == "Röthenbach i.E." 
replace gemeindename = "Wynigen" if canton == "BE" & city == "Rüedisbach" // Rüedisbach gehört zu Wynigen.
replace gemeindename = "Worb" if canton == "BE" & city == "Rüfenacht" // Rüfenacht ist Teil von Worb.
replace gemeindename = "Worb" if canton == "BE" & city == "Rüfenacht/Worb" // Rüfenacht ist Teil von Worb.
replace gemeindename = "Rüschegg" if canton == "BE" & city == "RüscheggGraben" // Rüschegg-Graben ist Teil von Rüschegg.
replace gemeindename = "Rüschegg" if canton == "BE" & city == "Rüschegggraben" // Rüschegg-Graben ist Teil von Rüschegg.
replace gemeindename = "Rüti (ZH)" if canton == "ZH" & city == "Rüti" 
replace gemeindename = "Rüti (ZH)" if canton == "ZH" & city == "Rüti ZH" 
replace gemeindename = "Rüti (ZH)" if canton == "ZH" & city == "Rüti/Bülach" // Kontrolle notwendig! Rüti und Bülach sind 2 eigenständige und weit auseinadnerliegende Gemeinden.
replace gemeindename = "Saas Fee" if canton == "VS" & city == "Saas-Fee" 
replace gemeindename = "Saignelégier" if canton == "BE" & city == "Saigne-légier" 
replace gemeindename = "Chézard-Saint-Martin" if canton == "NE" & city == "Saint-Martin" // Saint-Martin ist Teil von Chézard-Saint-Martin.
replace gemeindename = "Saint-Saphorin (Lavaux)" if canton == "VD" & city == "Saint-Saphorin/Lavaux" 
replace gemeindename = "Saint-Sulpice (VD)" if canton == "VD" & city == "Saint-Sulpice" 
replace gemeindename = "St. Silvester" if canton == "FR" & city == "Saint-Sylvestre" 
replace gemeindename = "Bauma" if canton == "ZH" & city == "Saland" // Saland ist Teil von Bauma.
replace gemeindename = "Constantine" if canton == "VD" & city == "Salavaux" // Salavaux ist Teil von Constantine.
replace gemeindename = "Sennwald" if canton == "SG" & city == "Salez" // Salez ist Teil von Sennwald.
replace gemeindename = "Sennwald" if canton == "SG" & city == "Salez-Sennwald" // Salez ist Teil von Sennwald.
replace gemeindename = "Salgesch" if canton == "VS" & city == "Salquenen" // Salquenen ist frz. Name für Salgesch.
replace gemeindename = "Salvenach" if canton == "FR" & city == "Salvagny" // Salvagny ist frz. Name für Salvenach.
replace gemeindename = "Horw" if canton == "LU" & city == "Sankt Niklausen" // St. Niklausen ist Teil von Horw.
replace gemeindename = "Savigny" if canton == "VD" & city == "Savigny " 
replace gemeindename = "Sennwald" if canton == "SG" & city == "Sax-Sennwald" // Sax ist Teil von Sennwald.
replace gemeindename = "Schlatt" if canton == "ZH" & city == "Schlatt b/Winterthur" // Eigentlich bis 1999 Schlatt bei Winterthur.
replace gemeindename = "Schmitten (FR)" if canton == "FR" & city == "Schmitten" 
replace gemeindename = "Wahlern" if canton == "BE" & city == "Schwarzenburg" // Schwarzenburg war Teil von Wahlern (heute Schwarzenburg).
replace gemeindename = "Unterlangenegg" if canton == "BE" & city == "Schwarzenegg" // Schwarzenegg ist Teil von Unterlangenegg.
replace gemeindename = "Sempach" if canton == "LU" & city == "Sempach-Station" // Sempach-Station ist Teil von Sempach.
replace gemeindename = "Sierre" if canton == "VS" & city == "Siders" // Siders ist dt. Name für Sierre.
replace gemeindename = "Schübelbach" if canton == "SZ" & city == "Siebnen" // Kontrolle notwendig! Siebnen ist auf 3 Gemeinden verteilt. Bei 2 Kandidaten (Fritz Stähli und Johann Wattenhofer) werden Siebnen und Siebnen-Schübelbach erwähnt; bei 2 anderen (Josef Diethelm und Josef Kürzi) nur Siebnen. Josef Diethelm ist gemäss HLS aus Schübelbach (Geburts- und Todesort)
replace gemeindename = "Sierre" if canton == "VS" & city == "Sierre/Bagnes" // Sierre und Bagnes sind 2 eigenständige Gemeinden.
replace gemeindename = "Signy-Avenex" if canton == "VD" & city == "Signy" 
replace gemeindename = "Sonceboz-Sombeval" if canton == "BE" & city == "Sonceboz" // Sonceboz ist Teil von Sonceboz-Sombeval.
replace gemeindename = "Köniz" if canton == "BE" & city == "Spiegel bei Bern" // Spiegel ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "Spiegel-Köniz" // Spiegel ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "Spiegel/Berne" // Spiegel ist Teil von Köniz.
replace gemeindename = "St. Margrethen" if canton == "SG" & city == "St Margrethen" 
replace gemeindename = "Saint-Aubin (FR)" if canton == "FR" & city == "St-Aubin" 
replace gemeindename = "Saint-Blaise" if canton == "NE" & city == "St-Blaise" 
replace gemeindename = "Saint-Brais" if canton == "BE" & city == "St-Brais" 
replace gemeindename = "Saint-Imier" if canton == "BE" & city == "St-Imier" 
replace gemeindename = "Saint-Maurice" if canton == "VS" & city == "St-Maurice" 
replace gemeindename = "Saint-Saphorin (Lavaux)" if canton == "VD" & city == "St-Saphorin (Lavaux)" 
replace gemeindename = "Saint-Saphorin (Lavaux)" if canton == "VD" & city == "St-Saphorin (Lx)" 
replace gemeindename = "Saint-Ursanne" if canton == "BE" & city == "St-Ursanne" 
replace gemeindename = "Knutwil" if canton == "LU" & city == "St. Erhard" // St. Erhard ist Teil von Knutwil.
replace gemeindename = "St. Margrethen" if canton == "SG" & city == "St. Margrethen SG" 
replace gemeindename = "Horw" if canton == "LU" & city == "St. Niklausen" // St. Niklausen ist Teil von Horw.
replace gemeindename = "Pfaffnau" if canton == "LU" & city == "St. Urban" // St. Urban ist Teil von Pfaffnau.
replace gemeindename = "Saint-Imier" if canton == "BE" & city == "St.-Imier" 
replace gemeindename = "St. Gallen" if canton == "SG" & city == "St.Gallen" 
replace gemeindename = "Horw" if canton == "LU" & city == "St.Niklausen" // St. Niklausen ist Teil von Horw.
replace gemeindename = "Horw" if canton == "LU" & city == "St.Niklausen/ Horw" // St. Niklausen ist Teil von Horw.
replace gemeindename = "Thal" if canton == "SG" & city == "Staad-Thal" // Staad ist Teil von Thal.
replace gemeindename = "Stadel" if canton == "ZH" & city == "Stadel bei Niederglatt" 
replace gemeindename = "Sainte-Croix" if canton == "VD" & city == "Ste-Croix" 
replace gemeindename = "Steffisburg" if canton == "BE" & city == "Steffisberg" // Kontrolle notwendig! Vermnutlich Schreibfehler bei Alfred Kunz sonst Steffisburg.
replace gemeindename = "Steffisburg" if canton == "BE" & city == "Steffisburg-Dorf" 
replace gemeindename = "Steffisburg" if canton == "BE" & city == "Steffisburg-Thun" 
replace gemeindename = "Fischenthal" if canton == "ZH" & city == "Steg/Fischenthal" // Steg ist ein Teil von Fischenthal.
replace gemeindename = "Stein (AG)" if canton == "AG" & city == "Stein" 
replace gemeindename = "Egnach" if canton == "TG" & city == "Steinebrunn" // Steinebrunn ist Teil von Egnach.
replace gemeindename = "Gadmen" if canton == "BE" & city == "Steingletscher am Sustenpass" // Steingletscher am Sustenpass ist Teil von Gadmen.
replace gemeindename = "Sâles (Gruyère)" if canton == "FR" & city == "Sâles" 
replace gemeindename = "Dürnten" if canton == "ZH" & city == "Tann" // Tann ist Teil von Dürnten.
replace gemeindename = "Dürnten" if canton == "ZH" & city == "Tann-Dürnten" // Tann ist Teil von Dürnten.
replace gemeindename = "Tafers" if canton == "FR" & city == "Tavel" // Tavel ist frz. Name für Tafers.
replace gemeindename = "Torricella-Taverne" if canton == "TI" & city == "Taverne" 
replace gemeindename = "Thônex" if canton == "GE" & city == "Thonex" 
replace gemeindename = "Thun" if canton == "BE" & city == "Thun-Dürrenast" // Dürrenast ist Teil von Thun.
replace gemeindename = "Thun" if canton == "BE" & city == "Thun-Goldiwil" // Goldiwil ist Teil von Thun.
replace gemeindename = "Thun" if canton == "BE" & city == "Thun-Steffisburg" // Kontrolle notwendig! Thun und Steffisburg sind 2 eigenständige Gemeinden. Betrifft vermutlich nur Otto Wirth, sonst nur Thun.
replace gemeindename = "Hilterfingen" if canton == "BE" & city == "Thun/Hünibach" // Hünibach ist Teil von Hilterfingen im Bezerik Thun.
replace gemeindename = "Köniz" if canton == "BE" & city == "Thörishaus" // Thörishaus ist Teil von Köniz.
replace gemeindename = "Jona" if canton == "SG" & city == "Tona SG" // Kontrolle notwendig! Vermutlich Schreibfehler.
replace gemeindename = "Torricella-Taverne" if canton == "TI" & city == "Torricella" 
replace gemeindename = "La Tour-de-Trême" if canton == "FR" & city == "Tour-de-Trême" 
replace gemeindename = "Turtmann" if canton == "VS" & city == "Tourtemagne" // Tourtemagne ist frz. Name für Turtmann.
replace gemeindename = "Tramelan" if canton == "BE" & city == "Tramelan-Dessus" // Tramelan-Dessus ist Teil von Tramelan.
replace gemeindename = "Tramelan" if canton == "BE" & city == "Tramelan-dessus" // Tramelan-Dessus ist Teil von Tramelan.
replace gemeindename = "Tramelan" if canton == "BE" & city == "Tramelandessus" // Tramelan-Dessus ist Teil von Tramelan.
replace gemeindename = "Rubigen" if canton == "BE" & city == "Trimstein" // Treimstein gehörte zu Rubigen.
replace gemeindename = "Trub" if canton == "BE" & city == "Trüb" 
replace gemeindename = "Wartau" if canton == "SG" & city == "Trübbach" // Trübbach ist Teil von Wartau.
replace gemeindename = "Uetikon" if canton == "ZH" & city == "Uetikon am See" 
replace gemeindename = "Uitikon" if canton == "ZH" & city == "Uitikon-Waldegg" // Waldegg ist Teil von Uitikon.
replace gemeindename = "Alt St. Johann" if canton == "SG" & city == "Unterwasser Alt St. Johann" // Unterwasser war Teil von Alt St. Johann (heute Wildhaus-Alt St. Johann).
replace gemeindename = "Urtenen" if canton == "BE" & city == "Urtenen-Schönbühl" // Früher nur Urtenen (heute Urtenen-Schönbühl).
replace gemeindename = "Bösingen" if canton == "FR" & city == "Uttewil" // Uttewil ist Teil von Bösingen.
replace gemeindename = "Vandoeuvres" if canton == "GE" & city == "Vandœuvres" 
replace gemeindename = "Vaux-sur-Morges" if canton == "VD" & city == "Vaux s. Morges" 
replace gemeindename = "Lausanne" if canton == "VD" & city == "Vernand-Dessus" // Vernand-Dessous ist Teil von Lausanne.
replace gemeindename = "Les Verrières" if canton == "NE" & city == "Verrières" 
replace gemeindename = "Ormont-Dessus" if canton == "VD" & city == "Vers l'Eglise" // Vers-l'Eglise ist Teil von Ormont-Dessus.
replace gemeindename = "Ormont-Dessus" if canton == "VD" & city == "Vers l'Eglise :" // Vers-l'Eglise ist Teil von Ormont-Dessus.
replace gemeindename = "Ormont-Dessus" if canton == "VD" & city == "Vers-l'Eglise" // Vers-l'Eglise ist Teil von Ormont-Dessus.
replace gemeindename = "Veyrier" if canton == "GE" & city == "Vessy-Veyrier" // Vessy ist Teil von Veyrier.
replace gemeindename = "Vicques" if canton == "BE" & city == "Vieques" // Schreibfehler.
replace gemeindename = "Villars-sur-Glâne" if canton == "FR" & city == "Villars s. Glâne" 
replace gemeindename = "Villars-sur-Glâne" if canton == "FR" & city == "Villars s/Glane" 
replace gemeindename = "Villars-sur-Glâne" if canton == "FR" & city == "Villars s/Glâne" 
replace gemeindename = "Ollon" if canton == "VD" & city == "Villars-sur-Ollon" // Villars-sur-Ollon ist Teil von Ollon.
replace gemeindename = "Villarsel-sur-Marly" if canton == "FR" & city == "Villarsel s/Marly" 
replace gemeindename = "Villars-sur-Glâne" if canton == "FR" & city == "Villarssur-Glâne" 
replace gemeindename = "Villaz-Saint-Pierre" if canton == "FR" & city == "Villaz-St-Pierre" 
replace gemeindename = "Villaz-Saint-Pierre" if canton == "FR" & city == "VillazSt-Pierre" 
replace gemeindename = "Villeneuve (VD)" if canton == "VD" & city == "Villeneuve" 
replace gemeindename = "Zermatt" if canton == "VS" & city == "Visperterminen/Zermatt" // Kontrolle notwendig! Visperterminen und Zermatt sind 2 eigenständige Gemeinden. Betrifft vermutlich nur Walter Zimmermann, sonst in Zermatt.
replace gemeindename = "Visp" if canton == "VS" & city == "Viège" // Viège ist frz. Name für Visp.
replace gemeindename = "Vugelles-La Mothe" if canton == "VD" & city == "Vugelles-la-Mothe" 
replace gemeindename = "Vuiteboeuf" if canton == "VD" & city == "Vuitebœuf" 
replace gemeindename = "Collonge-Bellerive" if canton == "GE" & city == "Vésenaz" // Vésenaz ist Teil von Collonge-Bellerive.
replace gemeindename = "Köniz" if canton == "BE" & city == "Wabern" // Wabern ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "Wabern-Köniz" // Wabern ist Teil von Köniz.
replace gemeindename = "Wallisellen" if canton == "ZH" & city == "WalUsellen" 
replace gemeindename = "Wald (ZH)" if canton == "ZH" & city == "Wald" 
replace gemeindename = "St. Peterzell" if canton == "SG" & city == "Wald-Schönengrund" // Wald-Schönengrund war ein Teil von St. Peterzell (heute Neckertal).
replace gemeindename = "Walenstadt" if canton == "SG" & city == "Walenstadtberg" // Walenstadtberg ist Teil von Walenstadt.
replace gemeindename = "Wangen an der Aare" if canton == "BE" & city == "Wangen" 
replace gemeindename = "Wangen (ZH)" if canton == "ZH" & city == "Wangen" 
replace gemeindename = "Wangen an der Aare" if canton == "BE" & city == "Wangen a.d. Aare" 
replace gemeindename = "Regensdorf" if canton == "ZH" & city == "Watt-Regensdorf" // Watt ist Tail von Regensdorf.
replace gemeindename = "Weiningen (ZH)" if canton == "ZH" & city == "Weinigen" // Vermutlich Schreibfehler.
replace gemeindename = "Weiningen (ZH)" if canton == "ZH" & city == "Weiningen" 
replace gemeindename = "Wartau" if canton == "SG" & city == "Weite-Wartau" // Weite ist Teil von Wartau.
replace gemeindename = "Dinhard" if canton == "ZH" & city == "Welsikon-Dinhard" // Welsikon ist Teil von Dinhard.
replace gemeindename = "Lauterbrunnen" if canton == "BE" & city == "Wengen" // Wengen ist Teil  von Lauterbrunnen.
replace gemeindename = "Grabs" if canton == "SG" & city == "Werdenberg Grabs" 
replace gemeindename = "Uster" if canton == "ZH" & city == "Wermatswil/Uster" // Wermatswil ist Teil von Uster.
replace gemeindename = "Wettswil" if canton == "ZH" & city == "Wettswil a.A." 
replace gemeindename = "Wetzikon (ZH)" if canton == "ZH" & city == "Wetzikon" 
replace gemeindename = "Wil (SG)" if canton == "SG" & city == "Wil" 
replace gemeindename = "Wil (SG)" if canton == "SG" & city == "Wil SG" 
replace gemeindename = "Wil (SG)" if canton == "SG" & city == "WilSG" 
replace gemeindename = "Möriken-Wildegg" if canton == "AG" & city == "Wildegg" // Wildegg ist Teil von Möriken-Wildegg.
replace gemeindename = "Einsiedeln" if canton == "SZ" & city == "Willerzwil" // Kontrolle notwendig! Ort gibt es nicht. Vermmutlich Schreibfehler: Willerzell. Betrifft nur Hans Fuch, sonst Willerzell. Willerzell ist Teil von Einsiedeln.
replace gemeindename = "Willisau Land" if canton == "LU" & city == "Willisau" // Kontrolle notwendig! 2 Willisau (Willisau Stadt und Willisau Land). Betrifft 3 Kandidaten: Franz Jos. Kurmann, sonst Willisau-Land, sowie Erwin Muff und Fritz Grüter.
replace gemeindename = "Willisau Land" if canton == "LU" & city == "Willisau-Land" 
replace gemeindename = "Winterthur" if canton == "ZH" & city == "Winterthnr" 
replace gemeindename = "Wohlen (AG)" if canton == "AG" & city == "Wohlen" 
replace gemeindename = "Wölflinswil" if canton == "AG" & city == "Wolflinswil" 
replace gemeindename = "Bolligen" if canton == "BE" & city == "Worblaufen" // Worblaufen war Teil von Bolligen (heute Ittigen).
replace gemeindename = "Winterthur" if canton == "ZH" & city == "Wïnterthur" 
replace gemeindename = "Zell (LU)" if canton == "LU" & city == "Zeli" 
replace gemeindename = "Zollikofen" if canton == "BE" & city == "Zoffikofen" // Vermutlich Schreibfehler.
replace gemeindename = "Lauperswil" if canton == "BE" & city == "Zollbruck" & (name == "Geissbühler" | name == "Badertscher") // Kontrolle notwendig! Zollbrück ist Teil von Rüderswil und Lauperswil. Betrifft vermutlich drei Kandidaten: Franz Badertscher ohne weitere Angabe, Fritz Geissbühler, sonst Lauperswil, Ernst Hirsbrunner, Vater in Rüderswil.
replace gemeindename = "Lauperswil" if canton == "BE" & city == "Zollbrück" & (name == "Geissbühler" | name == "Badertscher") // Kontrolle notwendig! Zollbrück ist Teil von Rüderswil und Lauperswil. Betrifft vermutlich drei Kandidaten: Franz Badertscher ohne weitere Angabe, Fritz Geissbühler, sonst Lauperswil, Ernst Hirsbrunner, Vater in Rüderswil.
replace gemeindename = "Rüderswil" if canton == "BE" & city == "Zollbrück" & name == "Hirsbrunner" // Kontrolle notwendig! Zollbrück ist Teil von Rüderswil und Lauperswil. Betrifft vermutlich drei Kandidaten: Franz Badertscher ohne weitere Angabe, Fritz Geissbühler, sonst Lauperswil, Ernst Hirsbrunner, Vater in Rüderswil.
replace gemeindename = "Zollikon" if canton == "ZH" & city == "Zollikerberg" // Zollikerberg ist Teil von Zollikon.
replace gemeindename = "Zollikon" if canton == "ZH" & city == "Zolükon" 
replace gemeindename = "Zürich" if canton == "VD" & city == "Zurich" // Kontrolle notwendig!
replace gemeindename = "Zürich" if canton == "VS" & city == "Zurich" // Kontrolle notwendig!
replace gemeindename = "Zürich" if canton == "TI" & city == "Zurigo" // Kontrolle notwendig!
replace gemeindename = "Glattfelden" if canton == "ZH" & city == "Zweidlen" // Zweidlen ist Teil von Glattfelden.
replace gemeindename = "Glattfelden" if canton == "ZH" & city == "Zweidlen-Glattfelden" // Zweidlen ist Teil von Glattfelden.
replace gemeindename = "Affoltern am Albis" if canton == "ZH" & city == "Zwillikon" // Zwillikon ist Teil von Affoltern am Albis.
replace gemeindename = "Zuzwil (SG)" if canton == "SG" & city == "Züberwangen-Zuzwil" // Züberwangen ist Teil von Zuzwil.
replace gemeindename = "Le Bémont (BE)" if canton == "BE" & city == "aux Rouges-Terres" // Rouges-Terres ist Teil von Le Bémont (BE).
replace gemeindename = "Vevey" if canton == "VD" & city == "evey" // Kontrolle notwendig! Vermutlich Schreibfehler.
replace gemeindename = "Erlenbach (ZH)" if canton == "ZH" & city == "im Jungholz, Erlenbach" 
replace gemeindename = "Aadorf" if canton == "TG" & city == "in Aadorf" 
replace gemeindename = "Aarau" if canton == "AG" & city == "in Aarau" 
replace gemeindename = "Aarberg" if canton == "BE" & city == "in Aarberg" 
replace gemeindename = "Aarburg" if canton == "AG" & city == "in Aarburg" 
replace gemeindename = "Aarwangen" if canton == "BE" & city == "in Aarwangen" 
replace gemeindename = "Seegräben" if canton == "ZH" & city == "in Aathal" // Aathal ist Teil von Seegräben.
replace gemeindename = "Seegräben" if canton == "ZH" & city == "in Aathal-Seegräben" // Aathal ist Teil von Seegräben.
replace gemeindename = "Seegräben" if canton == "ZH" & city == "in Aathal/Seegräben" // Aathal ist Teil von Seegräben.
replace gemeindename = "Adelboden" if canton == "BE" & city == "in Adelboden" 
replace gemeindename = "Adliswil" if canton == "ZH" & city == "in Adliswil" 
replace gemeindename = "Aesch (BL)" if canton == "BL" & city == "in Aesch" 
replace gemeindename = "Aesch (BL)" if canton == "BL" & city == "in Aesch (Basel-Landschaft)" 
replace gemeindename = "Aesch bei Birmensdorf" if canton == "ZH" & city == "in Aesch bei Birmensdorf" 
replace gemeindename = "Aeschi bei Spiez" if canton == "BE" & city == "in Aeschi" 
replace gemeindename = "Affoltern am Albis" if canton == "ZH" & city == "in Aeugsterthal" // Aeugstertal ist Teil von Aeugst am Albis.
replace gemeindename = "Affoltern am Albis" if canton == "ZH" & city == "in Affoltern a. A." 
replace gemeindename = "Affoltern am Albis" if canton == "ZH" & city == "in Affoltern a. A. :" 
replace gemeindename = "Affoltern am Albis" if canton == "ZH" & city == "in Affoltern am Albis" 
replace gemeindename = "Affoltern am Albis" if canton == "ZH" & city == "in Affoltern amAlbis" 
replace gemeindename = "Alle" if canton == "BE" & city == "in Alle" 
replace gemeindename = "Allschwil" if canton == "BL" & city == "in Allschwil" 
replace gemeindename = "Allschwil" if canton == "BL" & city == "in Allsehwil" 
replace gemeindename = "Alpnach" if canton == "OW" & city == "in Alpnach" 
replace gemeindename = "Alpnach" if canton == "OW" & city == "in Alpnachstad" // Alpnachstad ist Teil von Alpnach.
replace gemeindename = "Alt St. Johann" if canton == "SG" & city == "in Alt St. Johann" 
replace gemeindename = "Altdorf (UR)" if canton == "UR" & city == "in Altdorf" 
replace gemeindename = "Altdorf (UR)" if canton == "UR" & city == "in Altdorf " 
replace gemeindename = "Altnau" if canton == "TG" & city == "in Altnau" 
replace gemeindename = "Altstätten" if canton == "SG" & city == "in Altstätten" 
replace gemeindename = "Amden" if canton == "SG" & city == "in Amden" 
replace gemeindename = "Ammerswil" if canton == "AG" & city == "in Ammerswil" 
replace gemeindename = "Amriswil" if canton == "TG" & city == "in Amriswil" 
replace gemeindename = "Grossandelfingen" if canton == "ZH" & city == "in Andelfingen" // Früherer Name.
replace gemeindename = "Appenzell" if canton == "AI" & city == "in Appenzell" 
replace gemeindename = "Appenzell" if canton == "AI" & city == "in Appenzell mit Stimmen" 
replace gemeindename = "Arbon" if canton == "TG" & city == "in Arbon" 
replace gemeindename = "Arlesheim" if canton == "BL" & city == "in Ariesheim" 
replace gemeindename = "Aristau" if canton == "AG" & city == "in Aristau (Murimoos)" 
replace gemeindename = "Arlesheim" if canton == "BL" & city == "in Arlesheim" 
replace gemeindename = "Arlesheim" if canton == "BS" & city == "in Arlesheim" 
replace gemeindename = "Arni" if canton == "BE" & city == "in Arnisäge" // Arnisäge ist Teil von Arni.
replace gemeindename = "Arosa" if canton == "GR" & city == "in Arosa" 
replace gemeindename = "Arth" if canton == "SZ" & city == "in Arth" 
replace gemeindename = "Aesch bei Birmensdorf" if canton == "ZH" & city == "in Asch" 
replace gemeindename = "Ascona" if canton == "BS" & city == "in Ascona" // Kontrolle notwendig!
replace gemeindename = "Au (SG)" if canton == "SG" & city == "in Au" 
replace gemeindename = "Au (SG)" if canton == "SG" & city == "in Au (St. G.)" 
replace gemeindename = "Au (SG)" if canton == "SG" & city == "in Au (St. Gallen)" 
replace gemeindename = "Wädenswil" if canton == "ZH" & city == "in Au-Wädenswil" 
replace gemeindename = "Auw" if canton == "AG" & city == "in Auw" 
replace gemeindename = "Bad Ragaz" if canton == "SG" & city == "in Bad Ragaz" 
replace gemeindename = "Baden" if canton == "AG" & city == "in Baden" 
replace gemeindename = "Baldingen" if canton == "AG" & city == "in Baldingen" 
replace gemeindename = "Balerna" if canton == "TI" & city == "in Balerna" 
replace gemeindename = "Balgach" if canton == "SG" & city == "in Balgach" 
replace gemeindename = "Balsthal" if canton == "SO" & city == "in Balsthal" 
replace gemeindename = "Bassersdorf" if canton == "ZH" & city == "in BaltenswilBassersdorf" // Baltenswil ist Teil von Bassersdorf.
replace gemeindename = "Murgenthal" if canton == "AG" & city == "in Balzenwil-Murgenthal" // Balzenwil ist Teil von Murgenthal.
replace gemeindename = "Bangerten" if canton == "BE" & city == "in Bangerten" 
replace gemeindename = "Rapperswil (SG)" if canton == "SG" & city == "in Bapperswil" 
replace gemeindename = "Basadingen" if canton == "TG" & city == "in Basadingen" 
replace gemeindename = "Basel" if canton == "BS" & city == "in Base]" 
replace gemeindename = "Basel" if canton == "AG" & city == "in Basel" // Kontrolle notwendig!
replace gemeindename = "Basel" if canton == "BE" & city == "in Basel" // Kontrolle notwendig!
replace gemeindename = "Basel" if canton == "BS" & city == "in Basel" 
replace gemeindename = "Basel" if canton == "LU" & city == "in Basel" // Kontrolle notwendig!
replace gemeindename = "Basel" if canton == "SG" & city == "in Basel" // Kontrolle notwendig!
replace gemeindename = "Basel" if canton == "SH" & city == "in Basel" // Kontrolle notwendig!
replace gemeindename = "Basel" if canton == "ZH" & city == "in Basel" // Kontrolle notwendig!
replace gemeindename = "Basel" if canton == "BS" & city == "in Basel :" 
replace gemeindename = "Riehen" if canton == "ZH" & city == "in Basel-Riehen" // Kontrolle notwendig!
replace gemeindename = "Basel" if canton == "TI" & city == "in Basilea" // Kontrolle notwendig!
replace gemeindename = "Bassersdorf" if canton == "ZH" & city == "in Bassersdorf" 
replace gemeindename = "Bassersdorf" if canton == "ZH" & city == "in Bassersdorf-Baltenswil" // Baltenswil ist Teil von Bassersdorf.
replace gemeindename = "Bauma" if canton == "ZH" & city == "in Bauma" 
replace gemeindename = "Kirchberg (SG)" if canton == "SG" & city == "in Bazenheid" // Bazenheid ist Teil von Kirchberg.
replace gemeindename = "Beckenried" if canton == "NW" & city == "in Beckenried" 
replace gemeindename = "Beinwil am See" if canton == "AG" & city == "in Beinwil a. S." 
replace gemeindename = "Beinwil am See" if canton == "AG" & city == "in Beinwil a. See" 
replace gemeindename = "Beinwil am See" if canton == "AG" & city == "in Beinwil am See" 
replace gemeindename = "Bellinzona" if canton == "GR" & city == "in Bellinzona" // Kontrolle notwendig!
replace gemeindename = "Bellinzona" if canton == "TI" & city == "in Bellinzona" 
replace gemeindename = "Bellinzona" if canton == "TI" & city == "in Bellinzona (Daro)" // Daro ist Teil von Bellinzona.
replace gemeindename = "Bellinzona" if canton == "TI" & city == "in Bellinzona Canton de Vaud" // Kontrolle notwendig! Eintrag scheint fehlerhaft.
replace gemeindename = "Belp" if canton == "BE" & city == "in Belp" 
replace gemeindename = "Belprahon" if canton == "BE" & city == "in Belprahon" 
replace gemeindename = "Benken (SG)" if canton == "SG" & city == "in Benken" 
replace gemeindename = "Bennwil" if canton == "BL" & city == "in Bennwil" 
replace gemeindename = "Bergün/Bravuogn" if canton == "GR" & city == "in Bergün" 
replace gemeindename = "Oberhelfenschwil" if canton == "SG" & city == "in Berlig-Oberhelfenschwil" // Berlig ist vermutlich Teil von Oberhelfenschwil.
replace gemeindename = "Bern" if canton == "BE" & city == "in Bern" 
replace gemeindename = "Bern" if canton == "BE" & city == "in Bern-Bümpliz" 
replace gemeindename = "Bern" if canton == "BE" & city == "in Bern-Eymatt" 
replace gemeindename = "Köniz" if canton == "BE" & city == "in Bern-Liebefeld" // Liebefeld ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "in Bern-Spiegel" // Spiegel ist Teil von Köniz.
replace gemeindename = "Bern" if canton == "BE" & city == "in Bern-Wyler" // Kontrolle notwendig! Zusatz Wyler nicht gefunden. Betrifft vermutlich nur Hans Kästli, sonst Bern.
replace gemeindename = "Berneck" if canton == "SG" & city == "in Berneck" 
replace gemeindename = "Gossau (ZH)" if canton == "ZH" & city == "in Bertschikon-Gossau" // Bertschikon ist Teil von Gossau.
replace gemeindename = "Bettlach" if canton == "SO" & city == "in Bettlach" 
replace gemeindename = "Bettwil" if canton == "AG" & city == "in Bettwil" 
replace gemeindename = "Biasca" if canton == "TI" & city == "in Biasca" 
replace gemeindename = "Biberist" if canton == "SO" & city == "in Biberist" 
replace gemeindename = "Bichelsee" if canton == "TG" & city == "in Bichelsee" 
replace gemeindename = "Biel (BE)" if canton == "BE" & city == "in Biel" 
replace gemeindename = "Biel (BL)" if canton == "BL" & city == "in Biel (Baselland)" 
replace gemeindename = "Biel (BL)" if canton == "BL" & city == "in Biel-Benken" // Kontrolle notwendig! Bis 1971 2 eigenständige Gemeinden. Betrifft vermutlich nur 2 Kandidaten: Chrstoph Brodbeck, sonst Biel, und Paul Wyss.
replace gemeindename = "Biel (BE)" if canton == "BE" & city == "in Biel-Mett" // Mett ist Teil von Biel.
replace gemeindename = "Riehen" if canton == "BS" & city == "in Bienen" // Kontrolle notwendig! Vermutlich Schreibfehler.
replace gemeindename = "Biel (BE)" if canton == "BE" & city == "in Bienne" 
replace gemeindename = "Biglen" if canton == "BE" & city == "in Biglen" 
replace gemeindename = "Binningen" if canton == "BL" & city == "in Binningen" 
replace gemeindename = "Bioggio" if canton == "TI" & city == "in Bioggio" 
replace gemeindename = "Birmensdorf (ZH)" if canton == "ZH" & city == "in Birmensdorf" 
replace gemeindename = "Birsfelden" if canton == "BL" & city == "in Birsfelden" 
replace gemeindename = "Bischofszell" if canton == "TG" & city == "in Bischofszell" 
replace gemeindename = "Zweisimmen" if canton == "BE" & city == "in Blankenburg" // Blankenburg ist Teil von Zweisimmen.
replace gemeindename = "Lichtensteig" if canton == "SG" & city == "in Blatten-Lichtensteig" // Blatten ist Teil von Lichtensteig.
replace gemeindename = "Blauen" if canton == "BE" & city == "in Blauen" 
replace gemeindename = "Ebnat-Kappel" if canton == "SG" & city == "in Blomberg-Kappel" // Blomber ist vermutlich Teil von Kappel bzw später Ebnat-Kappel.
replace gemeindename = "Bodio" if canton == "TI" & city == "in Bodio" 
replace gemeindename = "Bolligen" if canton == "BE" & city == "in Bolligen" 
replace gemeindename = "Boltigen" if canton == "BE" & city == "in Boltigen" 
replace gemeindename = "Bonaduz" if canton == "GR" & city == "in Bonaduz" 
replace gemeindename = "Bonstetten" if canton == "ZH" & city == "in Bonstetten" 
replace gemeindename = "Boswil" if canton == "AG" & city == "in Boswil" 
replace gemeindename = "Flawil" if canton == "SG" & city == "in Botsberg-Flawil" // Botsberg ist Teil von Flawil.
replace gemeindename = "Bottmingen" if canton == "BL" & city == "in Bottmingen" 
replace gemeindename = "Bowil" if canton == "BE" & city == "in Bowil" 
replace gemeindename = "Zernez" if canton == "GR" & city == "in Brail-Zernez" // Brail ist Teil von Zernez.
replace gemeindename = "Breganzona" if canton == "TI" & city == "in Breganzona" 
replace gemeindename = "Breitenbach" if canton == "SO" & city == "in Breitenbach" 
replace gemeindename = "Bremgarten (AG)" if canton == "AG" & city == "in Bremgarten" 
replace gemeindename = "Bremgarten bei Bern" if canton == "BE" & city == "in Bremgarten" 
replace gemeindename = "Bremgarten bei Bern" if canton == "BE" & city == "in Bremgarten bei Bern" 
replace gemeindename = "Brienz (BE)" if canton == "BE" & city == "in Brienz" 
replace gemeindename = "Brislach" if canton == "BE" & city == "in Brislach" 
replace gemeindename = "Brissago" if canton == "BS" & city == "in Brissago" // Kontrolle notwendig!
replace gemeindename = "Brittnau" if canton == "AG" & city == "in Brittnau" 
replace gemeindename = "Brugg" if canton == "AG" & city == "in Brugg" 
replace gemeindename = "Brunnadern" if canton == "SG" & city == "in Brunnadern" 
replace gemeindename = "Brusio" if canton == "GR" & city == "in Brusio" 
replace gemeindename = "Brügg" if canton == "BE" & city == "in Brügg" 
replace gemeindename = "Brügg" if canton == "BE" & city == "in Brügg bei Biel" 
replace gemeindename = "Brüttelen" if canton == "BE" & city == "in Brüttelen" 
replace gemeindename = "Wangen (ZH)" if canton == "ZH" & city == "in Brüttisellen" // Spätere Namensänderung zu Wangen-Brüttisellen.
replace gemeindename = "Bubikon" if canton == "ZH" & city == "in Bubikon" 
replace gemeindename = "Uesslingen" if canton == "TG" & city == "in Buch-Uesslingen" // Heutiger Name ist Uesslingen-Buch.
replace gemeindename = "Uesslingen" if canton == "TG" & city == "in Buch-Üsslingen" // Heutiger Name ist Uesslingen-Buch.
replace gemeindename = "Buchs (AG)" if canton == "AG" & city == "in Buchs" 
replace gemeindename = "Buchs (SG)" if canton == "SG" & city == "in Buchs" 
replace gemeindename = "Buchs (ZH)" if canton == "ZH" & city == "in Buchs" 
replace gemeindename = "Buchs (AG)" if canton == "AG" & city == "in Buchs (Aargau)" 
replace gemeindename = "Buchs (ZH)" if canton == "ZH" & city == "in Buchs (Zch.)" 
replace gemeindename = "Buchs (ZH)" if canton == "ZH" & city == "in Buchs (Zürich)" 
replace gemeindename = "Schaffhausen" if canton == "SH" & city == "in Buchthalen" // Buchthalen ist Teil von Schaffhausen.
replace gemeindename = "Buochs" if canton == "NW" & city == "in Buochs" 
replace gemeindename = "Büren an der Aare" if canton == "BE" & city == "in Buren a A" 
replace gemeindename = "Büren an der Aare" if canton == "BE" & city == "in Buren a. A." 
replace gemeindename = "Burg (AG)" if canton == "AG" & city == "in Burg" 
replace gemeindename = "Burg im Leimental" if canton == "BE" & city == "in Burg" 
replace gemeindename = "Burg (AG)" if canton == "AG" & city == "in Burg (Aargau)" 
replace gemeindename = "Burgdorf" if canton == "BE" & city == "in Burgdorf" 
replace gemeindename = "Burgistein" if canton == "BE" & city == "in Burgistein" 
replace gemeindename = "Bürglen (TG)" if canton == "TG" & city == "in Burglen" 
replace gemeindename = "Schübelbach" if canton == "SZ" & city == "in Buttikon" // Buttikon ist Teil von Schübelbach.
replace gemeindename = "Buus" if canton == "BL" & city == "in Buus" 
replace gemeindename = "Diemtigen" if canton == "BE" & city == "in Bächlen (Simmental)" // Bächlen ist Teil von Diemtigen.
replace gemeindename = "Diemtigen" if canton == "BE" & city == "in Bächlen i. S." // Bächlen ist Teil von Diemtigen.
replace gemeindename = "Bäretswil" if canton == "ZH" & city == "in Bäretswil" 
replace gemeindename = "Bätterkinden" if canton == "BE" & city == "in Bätterkinden" 
replace gemeindename = "Bönigen" if canton == "BE" & city == "in Bönigen" 
replace gemeindename = "Böttstein" if canton == "AG" & city == "in BöttT stein" 
replace gemeindename = "Rüegsau" if canton == "BE" & city == "in Büegsau" 
replace gemeindename = "Büetigen" if canton == "BE" & city == "in Büetigen" 
replace gemeindename = "Bülach" if canton == "ZH" & city == "in Bülach" 
replace gemeindename = "Bern" if canton == "BE" & city == "in Bümpliz" // Bümpliz ist Teil von Bern.
replace gemeindename = "Bern" if canton == "TI" & city == "in Bümpliz" // Kontrolle notwendig! Bümpliz ist Teil von Bern.
replace gemeindename = "Bünzen" if canton == "AG" & city == "in Bünzen" 
replace gemeindename = "Büren an der Aare" if canton == "BE" & city == "in Büren" // Kontrolle notwendig! 2 Büren in BE (Büren an der Aare und Büren zum Hof). Betrifft vermutlich nur Werner Mülchi aus Büren an der Aara: http://www.e-periodica.ch/cntmng?pid=geo-004:1965:63::711; 26.01.2017.
replace gemeindename = "Bürglen (TG)" if canton == "TG" & city == "in Bürglen" 
replace gemeindename = "Büsserach" if canton == "SO" & city == "in Büsserach" 
replace gemeindename = "Schübelbach" if canton == "SZ" & city == "in Büttikon" // Buttikon ist Teil von Schübelbach.
replace gemeindename = "Camorino" if canton == "TI" & city == "in Camerino" 
replace gemeindename = "Camorino" if canton == "TI" & city == "in Camorino" 
replace gemeindename = "Capolago" if canton == "TI" & city == "in Capolago" 
replace gemeindename = "Lugano" if canton == "TI" & city == "in Castagnola" // Castagnola-Cassarate ist Teil von Lugano.
replace gemeindename = "Celerina/Schlarigna" if canton == "GR" & city == "in Celerina" 
replace gemeindename = "Chiasso" if canton == "TI" & city == "in Chiasso" 
replace gemeindename = "Chur" if canton == "GR" & city == "in Chur" 
replace gemeindename = "Chur" if canton == "GR" & city == "in Chur " 
replace gemeindename = "Coglio" if canton == "TI" & city == "in Coglio" 
replace gemeindename = "Cortébert" if canton == "BE" & city == "in Cortébert" 
replace gemeindename = "Courroux" if canton == "BE" & city == "in Courcelon, Gemeinde Courroux" 
replace gemeindename = "Courrendlin" if canton == "BE" & city == "in Courrendlin" 
replace gemeindename = "Dagmersellen" if canton == "LU" & city == "in Dagmersellen" 
replace gemeindename = "Davos" if canton == "GR" & city == "in Davos" 
replace gemeindename = "Davos" if canton == "GR" & city == "in Davos-Clavadel" 
replace gemeindename = "Davos" if canton == "GR" & city == "in Davos-Dorf" 
replace gemeindename = "Davos" if canton == "GR" & city == "in Davos-Monstein" 
replace gemeindename = "Davos" if canton == "GR" & city == "in Davos-Platz" 
replace gemeindename = "Degersheim" if canton == "SG" & city == "in Degersheim" 
replace gemeindename = "Stettlen" if canton == "BE" & city == "in Deisswil" // Kontrolle notwendig! Ohne Zusatz nur bei Hans Winzenried, sonst Deisswil bei Stettlen. Deisswil ist Teil von Stettlen.
replace gemeindename = "Stettlen" if canton == "BE" & city == "in Deisswil b. Stettlen" // Deisswil ist Teil von Stettlen.
replace gemeindename = "Stettlen" if canton == "BE" & city == "in Deisswil bei Stettlen" // Deisswil ist Teil von Stettlen.
replace gemeindename = "Stettlen" if canton == "BE" & city == "in Deisswil/Stettlen" // Deisswil ist Teil von Stettlen.
replace gemeindename = "Delémont" if canton == "BE" & city == "in Delsberg" // Delsberg ist dt. Name für Delémont.
replace gemeindename = "Delémont" if canton == "BE" & city == "in Delsberg'" // Delsberg ist dt. Name für Delémont.
replace gemeindename = "Derendingen" if canton == "SO" & city == "in Derendingen" 
replace gemeindename = "Dielsdorf" if canton == "ZH" & city == "in Dielsdorf" 
replace gemeindename = "Rapperswil (BE)" if canton == "BE" & city == "in Dieterswil" // Dieterswil ist Teil von Rapperswil.
replace gemeindename = "Dietikon" if canton == "ZH" & city == "in Dietikon" 
replace gemeindename = "Dietlikon" if canton == "ZH" & city == "in Dietlikon" 
replace gemeindename = "Dietwil" if canton == "AG" & city == "in Dietwil" 
replace gemeindename = "Bauma" if canton == "ZH" & city == "in Dillhaus-Bauma" // Dillhaus ist Teil von Bauma.
replace gemeindename = "Disentis/Mustér" if canton == "GR" & city == "in Disentis" 
replace gemeindename = "Disentis/Mustér" if canton == "GR" & city == "in Disentis/ Muster" 
replace gemeindename = "Disentis/Mustér" if canton == "GR" & city == "in Disentis/ Mustér" 
replace gemeindename = "Disentis/Mustér" if canton == "GR" & city == "in Disentis/Muster" 
replace gemeindename = "Disentis/Mustér" if canton == "GR" & city == "in Disentis/Mustèr" 
replace gemeindename = "Dittingen" if canton == "BE" & city == "in Dittingen" 
replace gemeindename = "Dorf" if canton == "ZH" & city == "in Dorf bei Andelfingen" 
replace gemeindename = "Dornach" if canton == "SO" & city == "in Dornach" 
replace gemeindename = "Dottikon" if canton == "AG" & city == "in Dottikon" 
replace gemeindename = "Lütisburg" if canton == "SG" & city == "in Dufertswil-Lütisburg" // Tufertschwil ist Tei von Lütisburg.
replace gemeindename = "Duggingen" if canton == "BE" & city == "in Duggingen" 
replace gemeindename = "Wädenswil" if canton == "ZH" & city == "in Dächenwies-Wädenswil" // Dächenwies ist vermutlich Teil von Wädenswil.
replace gemeindename = "Dägerlen" if canton == "ZH" & city == "in Dägerlen" 
replace gemeindename = "Dällikon" if canton == "ZH" & city == "in Dällikon" 
replace gemeindename = "Dänikon" if canton == "ZH" & city == "in Dänikon" 
replace gemeindename = "Därstetten" if canton == "BE" & city == "in Därstetten" 
replace gemeindename = "Dättlikon" if canton == "ZH" & city == "in Dättlikon" 
replace gemeindename = "Dübendorf" if canton == "ZH" & city == "in Dübendorf" 
replace gemeindename = "Dürnten" if canton == "ZH" & city == "in Dürnten" 
replace gemeindename = "Dürrenäsch" if canton == "AG" & city == "in Dürrenäsch" 
replace gemeindename = "Hittnau" if canton == "ZH" & city == "in Dürstelen-Hittnau" // Dürstelen ist Teil von Hittnau.
replace gemeindename = "Elgg" if canton == "ZH" & city == "in EIgg" 
replace gemeindename = "Hausen am Albis" if canton == "ZH" & city == "in Ebertswil/Hausen a. A." // Ebertswil ist Teil von Hausen am Albis.
replace gemeindename = "Ebnat-Kappel" if canton == "SG" & city == "in Ebnat" // Bis 1965 eigenständige Gemeinde.
replace gemeindename = "Saanen" if canton == "BE" & city == "in Ebnit/ Saanen" // Ebnit ist Teil von Saanen.
replace gemeindename = "Saanen" if canton == "BE" & city == "in Ebnit/Saanen" // Ebnit ist Teil von Saanen.
replace gemeindename = "Reinach (BL)" if canton == "BL" & city == "in Eeinach" 
replace gemeindename = "Illnau" if canton == "ZH" & city == "in Effretikon" // Früherer Name für Illnau-Effretikon.
replace gemeindename = "Eggiwil" if canton == "BE" & city == "in Eggiwil" 
replace gemeindename = "Eglisau" if canton == "ZH" & city == "in Eglisau" 
replace gemeindename = "Uetendorf" if canton == "BE" & city == "in Eichberg bei Uetendorf" // Eichberg ist Teil von Uetendorf.
replace gemeindename = "Uetendorf" if canton == "BE" & city == "in Eichberg/ Uetendorf" // Eichberg ist Teil von Uetendorf.
replace gemeindename = "Richterswil" if canton == "ZH" & city == "in Eichterswil" 
replace gemeindename = "Riehen" if canton == "BS" & city == "in Eiehen" 
replace gemeindename = "Eiken" if canton == "AG" & city == "in Eiken" 
replace gemeindename = "Embrach" if canton == "ZH" & city == "in Einbrach" 
replace gemeindename = "Einsiedeln" if canton == "SZ" & city == "in Einsiedeln" 
replace gemeindename = "Kilchberg (ZH)" if canton == "ZH" & city == "in Elchberg (Zch.)" 
replace gemeindename = "Elfingen" if canton == "AG" & city == "in Elfingen" 
replace gemeindename = "Elgg" if canton == "ZH" & city == "in Elgg" 
replace gemeindename = "Ellikon an der Thur" if canton == "ZH" & city == "in Ellikon a. d. Th." 
replace gemeindename = "Ellikon an der Thur" if canton == "ZH" & city == "in Ellikon a. d. Thur" 
replace gemeindename = "Embrach" if canton == "ZH" & city == "in Embrach" 
replace gemeindename = "Embrach" if canton == "ZH" & city == "in Embrach (Zürich)" 
replace gemeindename = "Emmen" if canton == "LU" & city == "in Emmen" 
replace gemeindename = "Emmen" if canton == "LU" & city == "in Emmenbrücke" // Emmenbrücke ist Teil von Emmen.
replace gemeindename = "Kreuzlingen" if canton == "TG" & city == "in Emmishofen" // Emmishofen ist Teil von Kreuzlingen.
replace gemeindename = "Endingen" if canton == "AG" & city == "in Endingen" 
replace gemeindename = "Engelberg" if canton == "OW" & city == "in Engelberg" 
replace gemeindename = "Worb" if canton == "BE" & city == "in Enggistein" // Enggistein ist Teil von Worb.
replace gemeindename = "Worb" if canton == "BE" & city == "in Enggistein b. Worb" // Enggistein ist Teil von Worb.
replace gemeindename = "Worb" if canton == "BE" & city == "in Enggistein bei Worb" // Enggistein ist Teil von Worb.
replace gemeindename = "Worb" if canton == "BE" & city == "in Enggistein/Worb" // Enggistein ist Teil von Worb.
replace gemeindename = "Engwang" if canton == "TG" & city == "in Engwang" 
replace gemeindename = "Ennenda" if canton == "GL" & city == "in Ennenda" 
replace gemeindename = "Ennetbaden" if canton == "AG" & city == "in Ennetbaden" 
replace gemeindename = "Ennetbürgen" if canton == "NW" & city == "in Ennetbürgen" 
replace gemeindename = "Ennetbürgen" if canton == "NW" & city == "in Ennetbürgen " 
replace gemeindename = "Entlebuch" if canton == "LU" & city == "in Entlebuch" 
replace gemeindename = "Rorschach" if canton == "SG" & city == "in Eorschach" 
replace gemeindename = "Eriz" if canton == "BE" & city == "in Eriz-Rufenen" // Rufenen ist Teil von Eriz.
replace gemeindename = "Erlach" if canton == "BE" & city == "in Erlach" 
replace gemeindename = "Erlenbach (ZH)" if canton == "ZH" & city == "in Erlenbach" 
replace gemeindename = "Erlenbach (ZH)" if canton == "ZH" & city == "in Erlenbach, Kt. Zürich" 
replace gemeindename = "Erschwil" if canton == "SO" & city == "in Erschwil" 
replace gemeindename = "Ersigen" if canton == "BE" & city == "in Ersigen" 
replace gemeindename = "Eschenz" if canton == "TG" & city == "in Eschenz" 
replace gemeindename = "Eschlikon" if canton == "TG" & city == "in Eschlikon" 
replace gemeindename = "Escholzmatt" if canton == "LU" & city == "in Escholzmatt" 
replace gemeindename = "Egg" if canton == "ZH" & city == "in Esslingen" // Esslingen ist Teil von Egg.
replace gemeindename = "Ettenhausen" if canton == "TG" & city == "in Ettenhausen" 
replace gemeindename = "Etzelkofen" if canton == "BE" & city == "in Etzelkofen" 
replace gemeindename = "Fahrwangen" if canton == "AG" & city == "in Fahrwangen" 
replace gemeindename = "Unterengstringen" if canton == "ZH" & city == "in Fahrweid-Unterengstringen" // Fahrweid ist vermutlich Teil von Unterengstringen (Hinweis: Kloster Fahr).
replace gemeindename = "Fehraltorf" if canton == "ZH" & city == "in Fehraltorf" 
replace gemeindename = "Hombrechtikon" if canton == "ZH" & city == "in Feldbach-Oberschirmensee" // Feldbach und Oberschirmensee sind Teile von Hombrechtikon.
replace gemeindename = "Meilen" if canton == "ZH" & city == "in Feldmeilen" // Feldmeilen ist Teil von Meilen.
replace gemeindename = "Felsberg" if canton == "GR" & city == "in Felsberg" 
replace gemeindename = "Feuerthalen" if canton == "ZH" & city == "in Feuerthalen" 
replace gemeindename = "Gsteig" if canton == "BE" & city == "in Feutersoey" // Feutersoey ist Teil von Gsteig.
replace gemeindename = "Fischenthal" if canton == "ZH" & city == "in Fischenthal" 
replace gemeindename = "Fisibach" if canton == "AG" & city == "in Fisibach" 
replace gemeindename = "Fisibach" if canton == "AG" & city == "in Fisibach-Waldhausen" // Waldhausen ist Teil von Fisibach.
replace gemeindename = "Fisibach" if canton == "AG" & city == "in Fisisbach" 
replace gemeindename = "Flawil" if canton == "SG" & city == "in Flawil" 
replace gemeindename = "Flerden" if canton == "GR" & city == "in Flerden" 
replace gemeindename = "Flims" if canton == "GR" & city == "in Flims" 
replace gemeindename = "Flums" if canton == "SG" & city == "in Flums" 
replace gemeindename = "Flurlingen" if canton == "ZH" & city == "in Flurlingen" 
replace gemeindename = "Flüelen" if canton == "UR" & city == "in Flüelen" 
replace gemeindename = "Forst" if canton == "BE" & city == "in Forst" 
replace gemeindename = "Frasnacht" if canton == "TG" & city == "in Frasnacht" 
replace gemeindename = "Frauenfeld" if canton == "TG" & city == "in Frauenfeld" 
replace gemeindename = "Frenkendorf" if canton == "BL" & city == "in Frenkendorf" 
replace gemeindename = "Frick" if canton == "AG" & city == "in Frick" 
replace gemeindename = "Seedorf (BE)" if canton == "BE" & city == "in Frienisberg-Seedorf" // Frienisberg ist Teil von Seedorf.
replace gemeindename = "Sumiswald" if canton == "BE" & city == "in Fritzenhaus (Wasen i. E.)" // Fritzenhaus und Wasen sind Teil von Sumiswald.
replace gemeindename = "Frutigen" if canton == "BE" & city == "in Frutigen" 
replace gemeindename = "Sennwald" if canton == "SG" & city == "in Frümsen" // Frümsen ist Teil von Sennwald.
replace gemeindename = "Sennwald" if canton == "SG" & city == "in Frümsen-Sennwald" // Frümsen ist Teil von Sennwald.
replace gemeindename = "Fulenbach" if canton == "SO" & city == "in Fulenbach" 
replace gemeindename = "Rüti (ZH)" if canton == "ZH" & city == "in Fägswil-Rüti" // Fägswil ist Teil von Rüti.
replace gemeindename = "Fällanden" if canton == "ZH" & city == "in Fällanden" 
replace gemeindename = "Füllinsdorf" if canton == "BL" & city == "in Füllinsdorf" 
replace gemeindename = "Gais" if canton == "AR" & city == "in Gais" 
replace gemeindename = "Gams" if canton == "SG" & city == "in Gams" 
replace gemeindename = "Thalwil" if canton == "ZH" & city == "in Gattikon-Thalwil" // Gattikon ist Teil von Thalwil.
replace gemeindename = "Gebenstorf" if canton == "AG" & city == "in Gebenstorf" 
replace gemeindename = "Gelterkinden" if canton == "BL" & city == "in Gelterkinden" 
replace gemeindename = "Genève" if canton == "BE" & city == "in Genf" // Kontrolle notwendig!
replace gemeindename = "Genève" if canton == "ZH" & city == "in Genf" // Kontrolle notwendig!
replace gemeindename = "Gerlafingen" if canton == "SO" & city == "in Gerlafingen" 
replace gemeindename = "Gerra (Gambarogno)" if canton == "TI" & city == "in Gerra Gambarogno" // Gerra ist Teil von Gambarogno.
replace gemeindename = "Gerzensee" if canton == "BE" & city == "in Gerzensee" 
replace gemeindename = "Giebenach" if canton == "BL" & city == "in Giebenach" 
replace gemeindename = "Giornico" if canton == "TI" & city == "in Giornico" 
replace gemeindename = "Hinwil" if canton == "ZH" & city == "in Girenbad, Hinwil" // Girenbad ist Teil von Hinwil.
replace gemeindename = "Giswil" if canton == "OW" & city == "in Giswil" 
replace gemeindename = "Giubiasco" if canton == "TI" & city == "in Giubiasco" 
replace gemeindename = "Glarus" if canton == "GL" & city == "in Glarus" 
replace gemeindename = "Flawil" if canton == "SG" & city == "in Glattal bei Flawil '" // Glattal ist vermutlich Teil von Flawil.
replace gemeindename = "Opfikon" if canton == "ZH" & city == "in Glattbrugg" // Glattbrugg ist Teil von Opfikon.
replace gemeindename = "Glattfelden" if canton == "ZH" & city == "in Glattfelden" 
replace gemeindename = "Goldach" if canton == "SG" & city == "in Goldach" 
replace gemeindename = "Arth" if canton == "SZ" & city == "in Goldau" // Goldau ist Teil von Arth.
replace gemeindename = "Küsnacht (ZH)" if canton == "ZH" & city == "in Goldbach-Küsnacht" // Goldbach ist Teil von Küsnacht.
replace gemeindename = "Hasliberg" if canton == "BE" & city == "in Goldern/Hasliberg" // Goldern ist Teil von Hasliberg.
replace gemeindename = "Gommiswald" if canton == "SG" & city == "in Gommiswald" 
replace gemeindename = "Gondiswil" if canton == "BE" & city == "in Gondiswil" 
replace gemeindename = "Gontenschwil" if canton == "AG" & city == "in Gontenschwil" 
replace gemeindename = "Gossau (SG)" if canton == "SG" & city == "in Gossau" 
replace gemeindename = "Gossau (ZH)" if canton == "ZH" & city == "in Gossau" 
replace gemeindename = "Gossau (SG)" if canton == "SG" & city == "in Gossau (St Gallen)" 
replace gemeindename = "Gossau (SG)" if canton == "SG" & city == "in Gossau (St. Gallen)" 
replace gemeindename = "Gossau (ZH)" if canton == "ZH" & city == "in Gossau (Zch.)" 
replace gemeindename = "Graben" if canton == "BE" & city == "in Graben b. H." 
replace gemeindename = "Grabs" if canton == "SG" & city == "in Grabs" 
replace gemeindename = "Grabs" if canton == "SG" & city == "in Grabsberg" // Grabsberg ist vermutlich Teil von Grabs.
replace gemeindename = "Grafenried" if canton == "BE" & city == "in Grafenried" 
replace gemeindename = "Seeberg" if canton == "BE" & city == "in Grasswil" // Grasswil ist Teil von Seeberg.
replace gemeindename = "Seeberg" if canton == "BE" & city == "in Grasswil-Seeberg" // Grasswil ist Teil von Seeberg.
replace gemeindename = "Gravesano" if canton == "TI" & city == "in Gravesano" 
replace gemeindename = "Grellingen" if canton == "BE" & city == "in Grellingen" 
replace gemeindename = "Grenchen" if canton == "SO" & city == "in Grenchen" 
replace gemeindename = "Greppen" if canton == "LU" & city == "in Greppen" 
replace gemeindename = "Grindelwald" if canton == "BE" & city == "in Grindelwald" 
replace gemeindename = "Grod" if canton == "SO" & city == "in Grod b. Däniken" 
replace gemeindename = "Grono" if canton == "GR" & city == "in Grono" 
replace gemeindename = "Grossaffoltern" if canton == "BE" & city == "in Grossaffoltern" 
replace gemeindename = "Grossandelfingen" if canton == "ZH" & city == "in Grossandelfingen" 
replace gemeindename = "Grosshöchstetten" if canton == "BE" & city == "in Grosshöchstetten" 
replace gemeindename = "Grosswangen" if canton == "LU" & city == "in Grosswangen" 
replace gemeindename = "Gränichen" if canton == "AG" & city == "in Gränichen" 
replace gemeindename = "Saanen" if canton == "BE" & city == "in Gstaad" // Gstaad ist Teil von Saanen.
replace gemeindename = "Saanen" if canton == "BE" & city == "in Gsteig bei Saanen" // Gstaad ist Teil von Saanen.
replace gemeindename = "Guggisberg" if canton == "BE" & city == "in Guggisberg" 
replace gemeindename = "Sigriswil" if canton == "BE" & city == "in Gunten" // Gunten ist Teil von Sigriswil.
replace gemeindename = "Guntershausen bei Birwinken" if canton == "TG" & city == "in Guntershausen-Leimbach" // Leimbach ist Teil von Guntershausen bei Birwinken (heute Guntershausen bei Berg).
replace gemeindename = "Gurbrü" if canton == "BE" & city == "in Gurbrü" 
replace gemeindename = "Gurzelen" if canton == "BE" & city == "in Gurzelen" 
replace gemeindename = "Thun" if canton == "BE" & city == "in Gwatt bei Thun" // Gwatt ist Teil von Thun.
replace gemeindename = "Hinwil" if canton == "ZH" & city == "in Gyrenbad/Hinwil" // Girenbad ist Teil von Hinwil.
replace gemeindename = "Muri bei Bern" if canton == "BE" & city == "in Gümligen" // Gümligen ist Teil von Muri bei Bern.
replace gemeindename = "Hagenbuch" if canton == "ZH" & city == "in Hagenbuch" 
replace gemeindename = "Frauenfeld" if canton == "TG" & city == "in Haldenhof-Herten" // Herten ist Teil von Frauenfeld.
replace gemeindename = "Oberriet (SG)" if canton == "SG" & city == "in Hard-Oberriet" // Hard ist Teil von Oberriet.
replace gemeindename = "Hasle bei Burgdorf" if canton == "BE" & city == "in Hasle" 
replace gemeindename = "Hasle bei Burgdorf" if canton == "BE" & city == "in Hasle b. B." 
replace gemeindename = "Hasle bei Burgdorf" if canton == "BE" & city == "in Hasle-Burgdorf" 
replace gemeindename = "Haslen" if canton == "GL" & city == "in Haslen" 
replace gemeindename = "Hasliberg" if canton == "BE" & city == "in Hasliberg" 
replace gemeindename = "Hausen am Albis" if canton == "ZH" & city == "in Hausen a. A" 
replace gemeindename = "Hausen am Albis" if canton == "ZH" & city == "in Hausen a. A." 
replace gemeindename = "Hedingen" if canton == "ZH" & city == "in Hedingen" 
replace gemeindename = "Au (SG)" if canton == "SG" & city == "in Heerbrugg" // Heerbrugg ist Teil von Au.
replace gemeindename = "Egnach" if canton == "TG" & city == "in Hegi-Winden" // Hegi und Winden sind Teil von Egnach.
replace gemeindename = "Heiden" if canton == "AR" & city == "in Heiden" 
replace gemeindename = "Heiden" if canton == "BL" & city == "in Heiden" // Kontrolle notwendig!
replace gemeindename = "Heiligenschwendi" if canton == "BE" & city == "in Heiligenschwendi" 
replace gemeindename = "Heimberg" if canton == "BE" & city == "in Heimberg" 
replace gemeindename = "Neuenegg" if canton == "BE" & city == "in Heitern/Neuenegg" // Heitern ist Teil von Neuenegg.
replace gemeindename = "Uzwil" if canton == "SG" & city == "in Henau" // Henau ist Teil von Uzwil.
replace gemeindename = "Herdern" if canton == "TG" & city == "in Herdern" 
replace gemeindename = "Hergiswil bei Willisau" if canton == "LU" & city == "in Hergiswil" 
replace gemeindename = "Hergiswil (NW)" if canton == "NW" & city == "in Hergiswil" 
replace gemeindename = "Hergiswil bei Willisau" if canton == "LU" & city == "in Hergiswil b. W." 
replace gemeindename = "Herisau" if canton == "AR" & city == "in Herisau" 
replace gemeindename = "Herrliberg" if canton == "ZH" & city == "in Herrliberg" 
replace gemeindename = "Herrliberg" if canton == "ZH" & city == "in Herrliberg-Wetzwil" // Wetzwil ist Teil von Herrliberg.
replace gemeindename = "Herzogenbuchsee" if canton == "BE" & city == "in Herzogenbuchsee" 
replace gemeindename = "Aeschi bei Spiez" if canton == "BE" & city == "in Heustrich-Emdtal" // Heustrich in Emdthal ist Teil von Aeschi bei Spiez.
replace gemeindename = "Hilterfingen" if canton == "BE" & city == "in Hilterfingen" 
replace gemeindename = "Hindelbank" if canton == "BE" & city == "in Hindelbank" 
replace gemeindename = "Eichberg" if canton == "SG" & city == "in Hinterforst-Eichberg" // Hinterforst ist Teil von Eichberg.
replace gemeindename = "Hinwil" if canton == "ZH" & city == "in Hinwil" 
replace gemeindename = "Hinwil" if canton == "ZH" & city == "in Hinwil-Hadlikon" // Hadlikon ist Teil von Hinwil.
replace gemeindename = "Hirzel" if canton == "ZH" & city == "in Hirzel" 
replace gemeindename = "Hittnau" if canton == "ZH" & city == "in Hittnau" 
replace gemeindename = "Hitzkirch" if canton == "LU" & city == "in Hitzkirch" 
replace gemeindename = "Hochdorf" if canton == "LU" & city == "in Hochdorf" 
replace gemeindename = "Hofstetten (SO)" if canton == "SO" & city == "in Hofstetten" 
replace gemeindename = "Hofstetten bei Brienz" if canton == "BE" & city == "in Hofstetten b. B." 
replace gemeindename = "Hohenrain" if canton == "LU" & city == "in Hohenrain" 
replace gemeindename = "Griesenberg" if canton == "TG" & city == "in Holzhof" // Kontrolle notwendig! Holzhof ohne Zusatz Fimmelsberg betrifft vermutlich nur Otto Wartmann. Fimmelsberg war Teil von Griesenberg (heute Amlikon-Bissegg).
replace gemeindename = "Griesenberg" if canton == "TG" & city == "in Holzhof-Fimmelsberg" // Fimmelsberg war Teil von Griesenberg (heute Amlikon-Bissegg).
replace gemeindename = "Holziken" if canton == "AG" & city == "in Holziken" 
replace gemeindename = "Homberg" if canton == "BE" & city == "in Homberg" 
replace gemeindename = "Hombrechtikon" if canton == "ZH" & city == "in Hombrechtikon" 
replace gemeindename = "Horgen" if canton == "ZH" & city == "in Horgen" 
replace gemeindename = "Horn" if canton == "TG" & city == "in Horn" 
replace gemeindename = "Hornussen" if canton == "AG" & city == "in Hornussen" 
replace gemeindename = "Horw" if canton == "LU" & city == "in Horw" 
replace gemeindename = "Gams" if canton == "SG" & city == "in Hub-Garns" // Hub ist vermutlich Teil von Gams.
replace gemeindename = "Humlikon" if canton == "ZH" & city == "in Humlikon" 
replace gemeindename = "Hunzenschwil" if canton == "AG" & city == "in Hunzenschwil" 
replace gemeindename = "Huttwil" if canton == "BE" & city == "in Huttwil" 
replace gemeindename = "Häggenschwil" if canton == "SG" & city == "in Häggenschwil" 
replace gemeindename = "Hägglingen" if canton == "AG" & city == "in Hägglingen" 
replace gemeindename = "Trub" if canton == "BE" & city == "in Häusern b. Trub" 
replace gemeindename = "Trub" if canton == "BE" & city == "in Häusern bei Trub" 
replace gemeindename = "St. Stephan" if canton == "BE" & city == "in Häusern/St. Stephan" // Häusern ist Tei von St. Stephan.
replace gemeindename = "Roggwil (TG)" if canton == "TG" & city == "in Häuslen-Boggwil" // Häuslen ist Teil von Roggwil.
replace gemeindename = "Roggwil (TG)" if canton == "TG" & city == "in Häuslen-Roggwil" // Häuslen ist Teil von Roggwil.
replace gemeindename = "Zürich" if canton == "ZH" & city == "in Höngg" // Höngg ist Teil von Zürich.
replace gemeindename = "Horn" if canton == "TG" & city == "in Hörn" 
replace gemeindename = "Hüntwangen" if canton == "ZH" & city == "in Hüntwangen" 
replace gemeindename = "Ilanz" if canton == "GR" & city == "in Ilanz" 
replace gemeindename = "Illnau" if canton == "ZH" & city == "in Illnau-Effretikon" // Früherer Name für Illnau-Effretikon.
replace gemeindename = "Innertkirchen" if canton == "BE" & city == "in Innertkirchen" 
replace gemeindename = "Ins" if canton == "BE" & city == "in Ins" 
replace gemeindename = "Interlaken" if canton == "BE" & city == "in Interlaken" 
replace gemeindename = "Inwil" if canton == "LU" & city == "in Inwil" 
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "in Irgenhausen-Pfäffikon" // Irgenhausen ist Teil von Pfäffikon.
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "in Irgenhausen-Pfäffikon (Zürich)" // Irgenhausen ist Teil von Pfäffikon.
replace gemeindename = "Islikon" if canton == "TG" & city == "in Islikon" 
replace gemeindename = "Itingen" if canton == "BL" & city == "in Itingen" 
replace gemeindename = "Küsnacht (ZH)" if canton == "ZH" & city == "in Itschnach/Küsnacht" // Itschnach ist Teil von Küsnacht.
replace gemeindename = "Bolligen" if canton == "BE" & city == "in Ittigen" // Ittigen gehörte zu Bolligen.
replace gemeindename = "Bolligen" if canton == "BE" & city == "in Ittigen bei Bern" // Ittigen gehörte zu Bolligen.
replace gemeindename = "Jegenstorf" if canton == "BE" & city == "in Jegenstorf" 
replace gemeindename = "Jona" if canton == "SG" & city == "in Jona" 
replace gemeindename = "Jonschwil" if canton == "SG" & city == "in Jonschwil" 
replace gemeindename = "Kallern" if canton == "AG" & city == "in Kallern" 
replace gemeindename = "Kaltbrunn" if canton == "SG" & city == "in Kaltbrunn" 
replace gemeindename = "Frutigen" if canton == "BE" & city == "in Kanderbrück" // Kanderbrück ist Teil von Frutigen.
replace gemeindename = "Kandersteg" if canton == "BE" & city == "in Kandersteg" 
replace gemeindename = "Kappel (SO)" if canton == "SO" & city == "in Kappel (Solothurn)" 
replace gemeindename = "Kappel am Albis" if canton == "ZH" & city == "in Kappel am Albis" 
replace gemeindename = "Wetzikon (ZH)" if canton == "ZH" & city == "in Kempten" // Kempten ist Teil von Wetzikon.
replace gemeindename = "Wetzikon (ZH)" if canton == "ZH" & city == "in Kempten (Zürich)" // Kempten ist Teil von Wetzikon.
replace gemeindename = "Wetzikon (ZH)" if canton == "ZH" & city == "in Kempten-Wetzikon" // Kempten ist Teil von Wetzikon.
replace gemeindename = "Wetzikon (ZH)" if canton == "ZH" & city == "in Kempten/Wetzikon" // Kempten ist Teil von Wetzikon.
replace gemeindename = "Wetzikon (ZH)" if canton == "ZH" & city == "in KemptenWetzikon" // Kempten ist Teil von Wetzikon.
replace gemeindename = "Lindau" if canton == "ZH" & city == "in Kempttal" // Kemptthal ist Teil von Lindau.
replace gemeindename = "Lindau" if canton == "ZH" & city == "in Kempttal-Lindau" // Kemptthal ist Teil von Lindau.
replace gemeindename = "Lindau" if canton == "ZH" & city == "in Kemptthal" // Kemptthal ist Teil von Lindau.
replace gemeindename = "Kestenholz" if canton == "SO" & city == "in Kestenholz" 
replace gemeindename = "Kilchberg (ZH)" if canton == "ZH" & city == "in Kilchberg" 
replace gemeindename = "Kilchberg (ZH)" if canton == "ZH" & city == "in Kilchberg (Zürich)" 
replace gemeindename = "Kilchberg (ZH)" if canton == "ZH" & city == "in Kilchberg b. Zch" 
replace gemeindename = "Kilchberg (ZH)" if canton == "ZH" & city == "in Kilchberg b.Zch" 
replace gemeindename = "Kilchberg (ZH)" if canton == "ZH" & city == "in Kilchberg bei Zürich" 
replace gemeindename = "Kirchberg (BE)" if canton == "BE" & city == "in Kirchberg" 
replace gemeindename = "Kirchberg (SG)" if canton == "SG" & city == "in Kirchberg" 
replace gemeindename = "Kirchberg (SG)" if canton == "SG" & city == "in Kirchberg (St. G.)" 
replace gemeindename = "Kirchberg (SG)" if canton == "SG" & city == "in Kirchberg (St. Gallen)" 
replace gemeindename = "Kirchleerau" if canton == "AG" & city == "in Kirchleerau" 
replace gemeindename = "Kirchlindach" if canton == "BE" & city == "in Kirchlindach" 
replace gemeindename = "Kleinandelfingen" if canton == "ZH" & city == "in Klein-Andelfingen" 
replace gemeindename = "Kleinandelfingen" if canton == "ZH" & city == "in Kleinandelfingen" 
replace gemeindename = "Klosters" if canton == "GR" & city == "in Klosters" 
replace gemeindename = "Kloten" if canton == "ZH" & city == "in Kloten" 
replace gemeindename = "Knutwil" if canton == "LU" & city == "in Knutwil" 
replace gemeindename = "Zell (ZH)" if canton == "ZH" & city == "in Kollbrunn" // Kollbrunn ist Teil von Zell.
replace gemeindename = "Zell (ZH)" if canton == "ZH" & city == "in Kollbrunn-Zell" // Kollbrunn ist Teil von Zell.
replace gemeindename = "Konolfingen" if canton == "BE" & city == "in Konolfingen" 
replace gemeindename = "Konolfingen" if canton == "BE" & city == "in Konolfingen-Stalden" // Stalden ist Teil von Konolfingen.
replace gemeindename = "Kradolf" if canton == "TG" & city == "in Kradolf" 
replace gemeindename = "Kreuzlingen" if canton == "TG" & city == "in Kreuzlingen" 
replace gemeindename = "Kreuzlingen" if canton == "TG" & city == "in Kreuzungen" 
replace gemeindename = "Kriens" if canton == "LU" & city == "in Kriens" 
replace gemeindename = "Wängi" if canton == "TG" & city == "in Krillberg-Wängi" // Krillberg ist Teil von Wängi.
replace gemeindename = "Bätterkinden" if canton == "BE" & city == "in Kräyligen" // Kräyligen ist Teil von Bätterkinden.
replace gemeindename = "Känerkinden" if canton == "BL" & city == "in Känerkinden" 
replace gemeindename = "Kölliken" if canton == "AG" & city == "in Kölliken" 
replace gemeindename = "Köniz" if canton == "BE" & city == "in Köniz" 
replace gemeindename = "Köniz" if canton == "BE" & city == "in Köniz-Liebefeld" // Liebefeld ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "in Köniz-Wabern" // Wabern ist Teil von Köniz.
replace gemeindename = "Künten" if canton == "AG" & city == "in Künten" 
replace gemeindename = "Küsnacht (ZH)" if canton == "ZH" & city == "in Küsnacht" 
replace gemeindename = "Küsnacht (ZH)" if canton == "ZH" & city == "in Küsnacht (Zch.)" 
replace gemeindename = "Küsnacht (ZH)" if canton == "ZH" & city == "in Küsnacht (Zürich)" 
replace gemeindename = "Küsnacht (ZH)" if canton == "ZH" & city == "in Küsnacht b. Zch" 
replace gemeindename = "Küsnacht (ZH)" if canton == "ZH" & city == "in Küsnacht bei Zürich" 
replace gemeindename = "Küsnacht (ZH)" if canton == "ZH" & city == "in Küsnacht, Kt. Zürich" 
replace gemeindename = "Küssnacht am Rigi" if canton == "SZ" & city == "in Küssnacht" 
replace gemeindename = "Küssnacht am Rigi" if canton == "SZ" & city == "in Küssnacht (Schwyz)" 
replace gemeindename = "Küssnacht am Rigi" if canton == "SZ" & city == "in Küssnacht a E" 
replace gemeindename = "Küssnacht am Rigi" if canton == "SZ" & city == "in Küssnacht a. R." 
replace gemeindename = "Küssnacht am Rigi" if canton == "SZ" & city == "in Küssnacht am Bigi" 
replace gemeindename = "Küttigkofen" if canton == "SO" & city == "in Küttigkofen" 
replace gemeindename = "La Chaux-de-Fonds" if canton == "BE" & city == "in La Chaux-de-Fonds" // Kontrolle notwendig!
replace gemeindename = "Laax" if canton == "GR" & city == "in Laax" 
replace gemeindename = "Lachen" if canton == "SZ" & city == "in Lachen" 
replace gemeindename = "Landiswil" if canton == "BE" & city == "in Landiswil" 
replace gemeindename = "Landiswil" if canton == "BE" & city == "in Landiswil-Obergoldbach" // Obergoldbach ist Teil von Landiswil.
replace gemeindename = "Igis" if canton == "GR" & city == "in Landquart" // Landquart war Teil von Igis (heute Landquart).
replace gemeindename = "Landschlacht" if canton == "TG" & city == "in Landschlacht" 
replace gemeindename = "Langenbruck" if canton == "BL" & city == "in Langenbruck" 
replace gemeindename = "Langenthal" if canton == "BE" & city == "in Langenthal" 
replace gemeindename = "Langnau im Emmental" if canton == "BE" & city == "in Langnau" 
replace gemeindename = "Langnau am Albis" if canton == "ZH" & city == "in Langnau a. A" 
replace gemeindename = "Langnau am Albis" if canton == "ZH" & city == "in Langnau a. A." 
replace gemeindename = "Langnau am Albis" if canton == "ZH" & city == "in Langnau amAlbis" 
replace gemeindename = "Langnau im Emmental" if canton == "BE" & city == "in Langnau i. E." 
replace gemeindename = "Feuerthalen" if canton == "ZH" & city == "in Langwiesen" // Langwiesen ist Teil von Feuerthalen.
replace gemeindename = "Läufelfingen" if canton == "BL" & city == "in Laufeltingen" 
replace gemeindename = "Laufen" if canton == "BE" & city == "in Laufen" 
replace gemeindename = "Laufen-Uhwiesen" if canton == "ZH" & city == "in Laufen-Uhwiesen" 
replace gemeindename = "Laufenburg" if canton == "AG" & city == "in Laufenburg" 
replace gemeindename = "Laupen" if canton == "BE" & city == "in Laupen" 
replace gemeindename = "Wald (ZH)" if canton == "ZH" & city == "in Laupen-Wald" // Laupen ist Teil von Wald.
replace gemeindename = "Lauperswil" if canton == "BE" & city == "in Lauperswil" 
replace gemeindename = "Lausanne" if canton == "BE" & city == "in Lausanne" // Kontrolle notwendig!
replace gemeindename = "Lausen" if canton == "BL" & city == "in Lausen" 
replace gemeindename = "Lauterbrunnen" if canton == "BE" & city == "in Lauterbrunnen" 
replace gemeindename = "Mühleberg" if canton == "BE" & city == "in Ledi bei Rosshäusern" // Ledi in Rosshäusern ist Teil von Mühleberg.
replace gemeindename = "Affoltern im Emmental" if canton == "BE" & city == "in Lehn-Weier i. E." // Weier im Emmental ist Teil von Affoltern im Emmental.
replace gemeindename = "Lengnau (BE)" if canton == "BE" & city == "in Lengnau" 
replace gemeindename = "Lenk" if canton == "BE" & city == "in Lenk (Bern)" 
replace gemeindename = "Lenk" if canton == "BE" & city == "in Lenk i. S." 
replace gemeindename = "Lenzburg" if canton == "AG" & city == "in Lenzburg" 
replace gemeindename = "Evilard" if canton == "BE" & city == "in Leubringen" // Leubringen ist dt. Name von Evilard.
replace gemeindename = "Leuggern" if canton == "AG" & city == "in Leuggern" 
replace gemeindename = "Leuzigen" if canton == "BE" & city == "in Leuzigen" 
replace gemeindename = "Lichtensteig" if canton == "SG" & city == "in Lichtensteig" 
replace gemeindename = "Köniz" if canton == "BE" & city == "in Liebefeld" // Liebefeld ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "in Liebefeld-Bern" // Liebefeld ist Teil von Köniz.
replace gemeindename = "Liesberg" if canton == "BE" & city == "in Liesberg" 
replace gemeindename = "Liestal" if canton == "BL" & city == "in Liestal" 
replace gemeindename = "Lindau" if canton == "ZH" & city == "in Lindau" 
replace gemeindename = "Linden" if canton == "BE" & city == "in Linden" 
replace gemeindename = "Littau" if canton == "LU" & city == "in Littau" 
replace gemeindename = "Locarno" if canton == "TI" & city == "in Locamo" 
replace gemeindename = "Locarno" if canton == "TI" & city == "in Locano" 
replace gemeindename = "Locarno" if canton == "TI" & city == "in Locarno" 
replace gemeindename = "Waldkirch" if canton == "SG" & city == "in Locherhof, Waldkirch" // Locherhof ist Teil von Waldkirch.
replace gemeindename = "Waldkirch" if canton == "SG" & city == "in Locherhof-Waldkirch" // Locherhof ist Teil von Waldkirch.
replace gemeindename = "Lausanne" if canton == "TI" & city == "in Losanna" // Kontrolle notwendig! Losanna ist it. Name für Lausanne.
replace gemeindename = "Lostorf" if canton == "SO" & city == "in Lostorf" 
replace gemeindename = "Ludiano" if canton == "TI" & city == "in Ludiano" 
replace gemeindename = "Lugano" if canton == "TI" & city == "in Lugano" 
replace gemeindename = "Lumbrein" if canton == "GR" & city == "in Lumbrein" 
replace gemeindename = "Luterbach" if canton == "SO" & city == "in Luterbach" 
replace gemeindename = "Luzein" if canton == "GR" & city == "in Luzein" 
replace gemeindename = "Luzern" if canton == "LU" & city == "in Luzern" 
replace gemeindename = "Lyss" if canton == "BE" & city == "in Lyss" 
replace gemeindename = "Lyssach" if canton == "BE" & city == "in Lyssach" 
replace gemeindename = "Lüsslingen" if canton == "SO" & city == "in Lüsslingen" 
replace gemeindename = "Lütisburg" if canton == "SG" & city == "in Lütisburg" 
replace gemeindename = "Lützelflüh" if canton == "BE" & city == "in Lützelflüh" 
replace gemeindename = "Madiswil" if canton == "BE" & city == "in Madiswil" 
replace gemeindename = "Magden" if canton == "AG" & city == "in Magden" 
replace gemeindename = "Maienfeld" if canton == "GR" & city == "in Maienfeld" 
replace gemeindename = "Malans" if canton == "GR" & city == "in Malans" 
replace gemeindename = "Malix" if canton == "GR" & city == "in Malix" 
replace gemeindename = "Malters" if canton == "LU" & city == "in Malters" 
replace gemeindename = "Mammern" if canton == "TG" & city == "in Mammern" 
replace gemeindename = "Mannenbach" if canton == "TG" & city == "in Mannenbach" 
replace gemeindename = "Marthalen" if canton == "ZH" & city == "in Marthalen" 
replace gemeindename = "Maschwanden" if canton == "ZH" & city == "in Maschwanden" 
replace gemeindename = "Massagno" if canton == "TI" & city == "in Massagno" 
replace gemeindename = "Matten bei Interlaken" if canton == "BE" & city == "in Matten" // Kontrolle notwendig! 2 Matten in BE (Matten im Simmental/Matten bei St. Stephan und Matten bei Interlagen). Matten ohne Zusatz betrifft vermutlich nur 2 Kandidaten: Walter Dürig, sonst Interlaken, Eduard Sterchi.
replace gemeindename = "St. Stephan" if canton == "BE" & city == "in Matten bei St Stephan" // Matte ist Teil von St. Stephan.
replace gemeindename = "St. Stephan" if canton == "BE" & city == "in Matten i. S." // Matte ist Teil von St. Stephan.
replace gemeindename = "Mattstetten" if canton == "BE" & city == "in Mattstetten" 
replace gemeindename = "Maur" if canton == "ZH" & city == "in Maur" 
replace gemeindename = "Meilen" if canton == "ZH" & city == "in Meilen" 
replace gemeindename = "Mellingen" if canton == "AG" & city == "in Meilingen" 
replace gemeindename = "Meiringen" if canton == "BE" & city == "in Meiringen" 
replace gemeindename = "Meisterschwanden" if canton == "AG" & city == "in Meisterschwanden" 
replace gemeindename = "Mellingen" if canton == "AG" & city == "in Mellingen" 
replace gemeindename = "Mels" if canton == "SG" & city == "in Mels" 
replace gemeindename = "Mendrisio" if canton == "TI" & city == "in Mendrisio" 
replace gemeindename = "Meiringen" if canton == "BE" & city == "in Meningen" // Kontrolle notwendig! Vermutlich Schreibfehler.
replace gemeindename = "Menziken" if canton == "AG" & city == "in Menziken" 
replace gemeindename = "Merenschwand" if canton == "AG" & city == "in Merenschwand" 
replace gemeindename = "Sigriswil" if canton == "BE" & city == "in Merligen" // Merligen ist Teil von Sigriswil.
replace gemeindename = "Merzligen" if canton == "BE" & city == "in Merzligen" 
replace gemeindename = "Mesocco" if canton == "GR" & city == "in Mesocco" 
replace gemeindename = "Biel (BE)" if canton == "BE" & city == "in Mett-Biel" // Mett ist Teil von Biel.
replace gemeindename = "Gossau (SG)" if canton == "SG" & city == "in Mettendorf-Gossau" // Mettendorf ist Teil von Gossau.
replace gemeindename = "Mettmenstetten" if canton == "ZH" & city == "in Mettmenstetten" 
replace gemeindename = "Minusio" if canton == "TI" & city == "in Minusio" 
replace gemeindename = "Mogelsberg" if canton == "SG" & city == "in Mogeisberg" 
replace gemeindename = "Mogelsberg" if canton == "SG" & city == "in Mogelsberg" 
replace gemeindename = "Oberriet (SG)" if canton == "SG" & city == "in Montlingen" // Montlingen ist Teil von Oberriet.
replace gemeindename = "Oberriet (SG)" if canton == "SG" & city == "in Montlingen-Oberriet" // Montlingen ist Teil von Oberriet.
replace gemeindename = "Moosleerau" if canton == "AG" & city == "in Moosleerau" 
replace gemeindename = "Münchenbuchsee" if canton == "BE" & city == "in Moospinte, Münchenbuchsee" 
replace gemeindename = "Muhen" if canton == "AG" & city == "in Muhen" 
replace gemeindename = "Mumpf" if canton == "AG" & city == "in Mumpf" 
replace gemeindename = "Muralto" if canton == "TI" & city == "in Muralto" 
replace gemeindename = "Quarten" if canton == "SG" & city == "in Murg" // Murg ist Teil von Quarten.
replace gemeindename = "Murgenthal" if canton == "AG" & city == "in Murgenthal" 
replace gemeindename = "Muri (AG)" if canton == "AG" & city == "in Muri" 
replace gemeindename = "Muri bei Bern" if canton == "BE" & city == "in Muri" 
replace gemeindename = "Muri bei Bern" if canton == "BE" & city == "in Muri b. B." 
replace gemeindename = "Muri bei Bern" if canton == "BE" & city == "in Muri b. Bern" 
replace gemeindename = "Muri bei Bern" if canton == "BE" & city == "in Muri bei Bern" 
replace gemeindename = "Muttenz" if canton == "BL" & city == "in Muttenz" 
replace gemeindename = "Muttenz" if canton == "BL" & city == "in Muttenz-Freidorf" // Freidorf ist Teil von Muttenz.
replace gemeindename = "Muttenz" if canton == "BL" & city == "in Muttenz-Rüttihardt" // Rüttihardt ist Teil von Muttenz.
replace gemeindename = "Männedorf" if canton == "ZH" & city == "in Männedorf" 
replace gemeindename = "Möhlin" if canton == "AG" & city == "in Möhlin" 
replace gemeindename = "Mönchaltorf" if canton == "ZH" & city == "in Mönchaltorf" 
replace gemeindename = "Möriken-Wildegg" if canton == "AG" & city == "in Möriken-Wildegg" 
replace gemeindename = "Wohlen bei Bern" if canton == "BE" & city == "in Möriswil" // Möriswil ist Teil von Wohlen bei Bern.
replace gemeindename = "Mörschwil" if canton == "SG" & city == "in Mörschwil" 
replace gemeindename = "Mühleberg" if canton == "BE" & city == "in Mühleberg" 
replace gemeindename = "Mühledorf (SO)" if canton == "SO" & city == "in Mühledorf" 
replace gemeindename = "Reichenbach im Kandertal" if canton == "BE" & city == "in Mülenen" // Mülenen ist Teil von Reichenbach im Kandertal.
replace gemeindename = "Münchenstein" if canton == "BL" & city == "in Münchenstein" 
replace gemeindename = "Münchenstein" if canton == "BL" & city == "in Münchenstein Post Neuewelt" 
replace gemeindename = "Kirchberg (SG)" if canton == "SG" & city == "in Münchwilen" // Münchwilen ist Teil von Kirchberg.
replace gemeindename = "Münchwilen (TG)" if canton == "TG" & city == "in Münchwilen" 
replace gemeindename = "Münsingen" if canton == "BE" & city == "in Müneingen" 
replace gemeindename = "Münsingen" if canton == "BE" & city == "in Münsingen" 
replace gemeindename = "Tägerwilen" if canton == "TG" & city == "in Nagelshausen" // Nagelshausen ist Teil von Tägerwilen.
replace gemeindename = "Mogelsberg" if canton == "SG" & city == "in Nassen" // Nassen war Teil von Mogelsberg (heute Neckertal).
replace gemeindename = "Niederhasli" if canton == "ZH" & city == "in Nassenwil-Niederhasli" // Nassenwil ist Teil von Niederhasli.
replace gemeindename = "Mogelsberg" if canton == "SG" & city == "in Necker" // Necker war Teil von Mogelsberg (heute Neckertal).
replace gemeindename = "Neftenbach" if canton == "ZH" & city == "in Neftenbach" 
replace gemeindename = "Nenzlingen" if canton == "BE" & city == "in Nenzlingen" 
replace gemeindename = "Nesslau" if canton == "SG" & city == "in Nesslau" 
replace gemeindename = "Netstal" if canton == "GL" & city == "in Netstal" 
replace gemeindename = "Nesslau" if canton == "SG" & city == "in Neu St. Johann" // Neu St. Johann ist Teil von Nesslau.
replace gemeindename = "Rüegsau" if canton == "BE" & city == "in Neuegg" // Neuegg ist vermutlich Teil von Rüegsau, Nachbargemeinde von Sumiswald.
replace gemeindename = "Rüegsau" if canton == "BE" & city == "in Neuegg-Rüegsau" // Neuegg ist vermutlich Teil von Rüegsau, Nachbargemeinde von Sumiswald.
replace gemeindename = "Rüegsau" if canton == "BE" & city == "in Neuegg-Sumiswald" // Neuegg ist vermutlich Teil von Rüegsau, Nachbargemeinde von Sumiswald.
replace gemeindename = "Neuenegg" if canton == "BE" & city == "in Neuenegg" 
replace gemeindename = "Neuenegg" if canton == "BE" & city == "in Neuenegg/Heitern" // Heitern ist Teil von Neuenegg.
replace gemeindename = "Neuenkirch" if canton == "LU" & city == "in Neuenkirch" 
replace gemeindename = "Neuhausen am Rheinfall" if canton == "SH" & city == "in Neuhausen" 
replace gemeindename = "Neuhausen am Rheinfall" if canton == "SH" & city == "in Neuhausen am Rheinfall" 
replace gemeindename = "Egnach" if canton == "TG" & city == "in Neukirch-Egnach" // Neukirch ist Teil von Egnach.
replace gemeindename = "Bäretswil" if canton == "ZH" & city == "in Neuthal-Bäretswil" 
replace gemeindename = "Wald (ZH)" if canton == "ZH" & city == "in Neuthal-Wald (Zürich)" 
replace gemeindename = "Nidau" if canton == "BE" & city == "in Nidau" 
replace gemeindename = "Sommeri" if canton == "TG" & city == "in Nieder-Sommeri" // Niedersommeri ist Teil von Sommeri.
replace gemeindename = "Uster" if canton == "ZH" & city == "in Nieder-Uster" // Niederuster ist Teil von Uster.
replace gemeindename = "Niederbipp" if canton == "BE" & city == "in Niederbipp" 
replace gemeindename = "Niederbüren" if canton == "SG" & city == "in Niederbüren" 
replace gemeindename = "Niederdorf" if canton == "BL" & city == "in Niederdorf" 
replace gemeindename = "Niedererlinsbach" if canton == "SO" & city == "in Niedererlinsbach (Solothurn)" 
replace gemeindename = "Gerlafingen" if canton == "SO" & city == "in Niedergerlafingen" // Früherer Name.
replace gemeindename = "Niederglatt" if canton == "ZH" & city == "in Niederglatt" 
replace gemeindename = "Niederlenz" if canton == "AG" & city == "in Niederlenz" 
replace gemeindename = "Niederrohrdorf" if canton == "AG" & city == "in Niederrohrdorf" 
replace gemeindename = "Köniz" if canton == "BE" & city == "in Niederscherli" // Niederscherli ist Teil von Köniz.
replace gemeindename = "Uzwil" if canton == "SG" & city == "in Niederuzwil" // Niederuzwil ist Teil von Uzwil.
replace gemeindename = "Niederweningen" if canton == "ZH" & city == "in Niederweningen" 
replace gemeindename = "Niederönz" if canton == "BE" & city == "in Niederönz" 
replace gemeindename = "Nunningen" if canton == "SO" & city == "in Nunningen" 
replace gemeindename = "Kirchberg (SG)" if canton == "SG" & city == "in Nutenwil-Bazenheid" // Nutenwil ist vermutlich Teil von Bazenheid, Bazenheid ist Teil von Kirchberg.
replace gemeindename = "Näfels" if canton == "GL" & city == "in Näfels" 
replace gemeindename = "Uster" if canton == "ZH" & city == "in Nänikon" // Nänikon ist Teil von Uster.
replace gemeindename = "Uster" if canton == "ZH" & city == "in Nänikon-Uster" // Nänikon ist Teil von Uster.
replace gemeindename = "Nürensdorf" if canton == "ZH" & city == "in Nürensdorf" 
replace gemeindename = "Oberembrach" if canton == "ZH" & city == "in Ober-Embrach" 
replace gemeindename = "Seuzach" if canton == "ZH" & city == "in Ober-Ohringen bei Seuzach" // Ober-Ohringen ist Teil von Seuzach.
replace gemeindename = "Uster" if canton == "ZH" & city == "in Ober-Uster" // Oberuster-Nossikon ist Teil von Uster.
replace gemeindename = "Winterthur" if canton == "ZH" & city == "in Ober-Winterthur" // Oberwinterthur ist Teil von Winterthur.
replace gemeindename = "Arth" if canton == "SZ" & city == "in Oberarth" // Oberarth ist Teil von Arth.
replace gemeindename = "Oberbalm" if canton == "BE" & city == "in Oberbalm" 
replace gemeindename = "Oberbipp" if canton == "BE" & city == "in Oberbipp" 
replace gemeindename = "Oberbuchsiten" if canton == "SO" & city == "in Oberbuchsiten" 
replace gemeindename = "Oberburg" if canton == "BE" & city == "in Oberburg" 
replace gemeindename = "Wohlen bei Bern" if canton == "BE" & city == "in Oberdettigen" // Oberdettigen ist Teil von Wohlen bei Bern.
replace gemeindename = "Oberdiessbach" if canton == "BE" & city == "in Oberdiessbach" 
replace gemeindename = "Bäretswil" if canton == "ZH" & city == "in OberdorfBäretswil" // Oberdorf ist vermutlich Teil von Bäretswil.
replace gemeindename = "Oberembrach" if canton == "ZH" & city == "in Oberembrach" 
replace gemeindename = "Endingen" if canton == "AG" & city == "in Oberendingen" // Früherer Name.
replace gemeindename = "Oberengstringen" if canton == "ZH" & city == "in Oberengstringen" 
replace gemeindename = "Oberentfelden" if canton == "AG" & city == "in Oberentfelden" 
replace gemeindename = "Oberglatt" if canton == "ZH" & city == "in Oberglatt" 
replace gemeindename = "Obergösgen" if canton == "SO" & city == "in Obergösgen" 
replace gemeindename = "Niederhasli" if canton == "ZH" & city == "in Oberhasli" // Oberhasli ist Teil von Niederhasli.
replace gemeindename = "Oberhelfenschwil" if canton == "SG" & city == "in Oberhelfenschwil" 
replace gemeindename = "Oberhof" if canton == "AG" & city == "in Oberhof" 
replace gemeindename = "Oberhofen am Thunersee" if canton == "BE" & city == "in Oberhofen" 
replace gemeindename = "Münchwilen (TG)" if canton == "TG" & city == "in Oberhofen-Münchwilen" // Oberhofen ist Teil von Münchwilen.
replace gemeindename = "Oberkulm" if canton == "AG" & city == "in Oberkulm" 
replace gemeindename = "Oberlunkhofen" if canton == "AG" & city == "in Oberlunkhofen" 
replace gemeindename = "Oberried am Brienzersee" if canton == "BE" & city == "in Oberried a. Brienzersee" 
replace gemeindename = "Oberrieden" if canton == "ZH" & city == "in Oberrieden" 
replace gemeindename = "Oberriet (SG)" if canton == "SG" & city == "in Oberriet" 
replace gemeindename = "Oberriet (SG)" if canton == "SG" & city == "in Oberriet (St. Gallen)" 
replace gemeindename = "Hombrechtikon" if canton == "ZH" & city == "in Oberschirmensee-Feldbach" // Feldbach und Oberschirmensee sind Teile von Hombrechtikon.
replace gemeindename = "Hombrechtikon" if canton == "ZH" & city == "in Oberschirmensee-Hombrechtikon" // Feldbach und Oberschirmensee sind Teile von Hombrechtikon.
replace gemeindename = "Oberstammheim" if canton == "ZH" & city == "in Oberstammheim" 
replace gemeindename = "Obersteckholz" if canton == "BE" & city == "in Obersteckholz" 
replace gemeindename = "Uzwil" if canton == "SG" & city == "in Oberstetten-Hena" // Oberstetten mit Henau ist Teil von Uzwil.
replace gemeindename = "Uster" if canton == "ZH" & city == "in Oberuster" // Oberuster-Nossikon ist Teil von Uster.
replace gemeindename = "Oberuzwil" if canton == "SG" & city == "in Oberuzwil" 
replace gemeindename = "Köniz" if canton == "BE" & city == "in Oberwangen" // Oberwangen ist Teil von Köniz .
replace gemeindename = "Oberweningen" if canton == "ZH" & city == "in Oberweningen" 
replace gemeindename = "Oberwichtrach" if canton == "BE" & city == "in Oberwichtrach" 
replace gemeindename = "Oberwil (BL)" if canton == "BL" & city == "in Oberwil" 
replace gemeindename = "Oberwil (BL)" if canton == "BL" & city == "in Oberwil (Baselland)" 
replace gemeindename = "Oberwil im Simmental" if canton == "BE" & city == "in Oberwil i. S." 
replace gemeindename = "Dägerlen" if canton == "ZH" & city == "in Oberwil-Dägerlen" // Oberwil ist Teil von Dägerlen.
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "in Oberwil-Pfäffikon" // Oberwil ist Teil von Pfäffikon.
replace gemeindename = "Winterthur" if canton == "ZH" & city == "in Oberwinterthur" // Oberwinterthur ist Teil von Winterthur.
replace gemeindename = "Obfelden" if canton == "ZH" & city == "in Obfelden" 
replace gemeindename = "Oensingen" if canton == "SO" & city == "in Oensingen" 
replace gemeindename = "Zürich" if canton == "ZH" & city == "in Oerlikon" // Oerlikon ist Teil von Zürich.
replace gemeindename = "Koppigen" if canton == "BE" & city == "in Oeschberg-Koppigen" // Oeschberg ist Teil von Koppigen.
replace gemeindename = "Diemtigen" if canton == "BE" & city == "in Oey" // Oey ist Teil von Diemtigen.
replace gemeindename = "Diemtigen" if canton == "BE" & city == "in Oey-Diemtigen" // Oey ist Teil von Diemtigen.
replace gemeindename = "Oftringen" if canton == "AG" & city == "in Oftringen" 
replace gemeindename = "Seuzach" if canton == "ZH" & city == "in Ohringen-Seuzach" // Ober- und Unterohringen sind Teil von Seuzach.
replace gemeindename = "Olten" if canton == "SO" & city == "in Olten" 
replace gemeindename = "Opfikon" if canton == "ZH" & city == "in Opfikon" 
replace gemeindename = "Oppligen" if canton == "BE" & city == "in Oppligen" 
replace gemeindename = "Steffisburg" if canton == "BE" & city == "in Ortbühl bei Steffisburg" // Ortbühl ist Teil von Steffisburg.
replace gemeindename = "Ossingen" if canton == "ZH" & city == "in Ossingen" 
replace gemeindename = "Bolligen" if canton == "BE" & city == "in Ostermundigen" // Ostermundigen gehörte zu Bolligen.
replace gemeindename = "Bolligen" if canton == "BE" & city == "in Ostermundigen-Bern" // Ostermundigen gehörte zu Bolligen.
replace gemeindename = "Otelfingen" if canton == "ZH" & city == "in Otelfingen" 
replace gemeindename = "Ottenbach" if canton == "ZH" & city == "in Ottenbach" 
replace gemeindename = "Gossau (ZH)" if canton == "ZH" & city == "in Ottikon-Gossau" // Ottikon ist Teil von Gossau.
replace gemeindename = "Unterschlatt" if canton == "TG" & city == "in Paradies-Schlatt" // Paradies war vermutlich Teil von Unterschlatt (heute Schlatt); http://www.hls-dhs-dss.ch/textes/d/D1890.php; 26.01.2017.
replace gemeindename = "Pfyn" if canton == "TG" & city == "in Pfyn" 
replace gemeindename = "Pfäfers" if canton == "SG" & city == "in Pfäfers" 
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "in Pfäffikon" 
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "in Pfäffikon (Zch.)" 
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "in Pfäffikon (Zürich)" 
replace gemeindename = "Pfäffikon" if canton == "ZH" & city == "in Pfäffikon-Zch." 
replace gemeindename = "Pieterlen" if canton == "BE" & city == "in Pieterlen" 
replace gemeindename = "Pohlern" if canton == "BE" & city == "in Pohlern" 
replace gemeindename = "Porrentruy" if canton == "BE" & city == "in Porrentruy" 
replace gemeindename = "Port" if canton == "BE" & city == "in Port" 
replace gemeindename = "Port" if canton == "BE" & city == "in Port bei Nidau" // Port ist Nachbargemeinde von Nidau.
replace gemeindename = "Poschiavo" if canton == "GR" & city == "in Poschiavo" 
replace gemeindename = "Pratteln" if canton == "BL" & city == "in Pratteln" 
replace gemeindename = "Porrentruy" if canton == "BE" & city == "in Pruntrut" // Pruntrut ist dt. Name für Porrentruy.
replace gemeindename = "Rafz" if canton == "ZH" & city == "in Rafz" 
replace gemeindename = "Rapperswil (SG)" if canton == "SG" & city == "in Rapperswil" 
replace gemeindename = "Rapperswil (SG)" if canton == "SG" & city == "in Rapperswil " 
replace gemeindename = "Rapperswil (SG)" if canton == "SG" & city == "in Rapperswil (St. Gallen)" 
replace gemeindename = "Rapperswil (BE)" if canton == "BE" & city == "in RapperswilWierezwil" // Wierezwil ist Teil von Rapperswil.
replace gemeindename = "Weinfelden" if canton == "TG" & city == "in Ratwies-Weinfelden" // Ratwies ist Teil von Weinfelden.
replace gemeindename = "Rebstein" if canton == "SG" & city == "in Rebstein" 
replace gemeindename = "Reconvilier" if canton == "BE" & city == "in Reconvilier" 
replace gemeindename = "Regensberg" if canton == "ZH" & city == "in Regensberg" 
replace gemeindename = "Regensdorf" if canton == "ZH" & city == "in Regensdorf" 
replace gemeindename = "Rehetobel" if canton == "AR" & city == "in Rehetobel" 
replace gemeindename = "Meiringen" if canton == "BE" & city == "in Reichenbach b. M." // Reichenbach ist vermutlich Teil von Meiringen.
replace gemeindename = "Meiringen" if canton == "BE" & city == "in Reichenbach bei Meiringen" // Reichenbach ist vermutlich Teil von Meiringen.
replace gemeindename = "Reichenbach im Kandertal" if canton == "BE" & city == "in Reichenbach im Kandertal" 
replace gemeindename = "Reichenburg" if canton == "SZ" & city == "in Reichenburg" 
replace gemeindename = "Reiden" if canton == "LU" & city == "in Reiden" 
replace gemeindename = "Reigoldswil" if canton == "BL" & city == "in Reigoldswil" 
replace gemeindename = "Reinach (AG)" if canton == "AG" & city == "in Reinach" 
replace gemeindename = "Reinach (BL)" if canton == "BL" & city == "in Reinach" 
replace gemeindename = "Reitnau" if canton == "AG" & city == "in Reitnau" 
replace gemeindename = "Rekingen (AG)" if canton == "AG" & city == "in Rekingen" 
replace gemeindename = "Reinach (BL)" if canton == "BL" & city == "in Remach" 
replace gemeindename = "Ramosch" if canton == "GR" & city == "in Remus" // Remüs ist dt. Name von Ramosch.
replace gemeindename = "Ramosch" if canton == "GR" & city == "in Remüs" // Remüs ist dt. Name von Ramosch.
replace gemeindename = "Luzern" if canton == "LU" & city == "in Reussbühl-Littau" // Kontrolle notwendig! Reussbühl ist Teil von Luzern; Luzern und Littau waren 2 eigenständige Gemeinden. Betrifft vermutlich jedoch nur Peter Brünisholz.
replace gemeindename = "Reute (AR)" if canton == "AR" & city == "in Reute" 
replace gemeindename = "Rheinau" if canton == "ZH" & city == "in Rheinau" 
replace gemeindename = "Rheineck" if canton == "SG" & city == "in Rheineck" 
replace gemeindename = "Rheinfelden" if canton == "AG" & city == "in Rheinfelden" 
replace gemeindename = "Ursenbach" if canton == "BE" & city == "in Richisberg-Ursenbach" // Richisberg ist Teil von Ursenbach.
replace gemeindename = "Richterswil" if canton == "ZH" & city == "in Richterswil" 
replace gemeindename = "Wattwil" if canton == "SG" & city == "in Ricken-Wattwil" // Ricken ist z.T. Teil von Wattwil.
replace gemeindename = "Rickenbach (SO)" if canton == "SO" & city == "in Rickenbach" 
replace gemeindename = "Rickenbach bei Wil" if canton == "TG" & city == "in Rickenbach" 
replace gemeindename = "Rickenbach (ZH)" if canton == "ZH" & city == "in Rickenbach" 
replace gemeindename = "Winterthur" if canton == "ZH" & city == "in Ricketwil-Winterthur" // Ricketwil ist Teil von Winterthur.
replace gemeindename = "Wallisellen" if canton == "ZH" & city == "in Rieden-Wallisellen" // Rieden ist Teil von Wallisellen.
replace gemeindename = "Riehen" if canton == "BL" & city == "in Riehen" // Kontrolle notwendig!
replace gemeindename = "Riehen" if canton == "BS" & city == "in Riehen" 
replace gemeindename = "Riehen" if canton == "BS" & city == "in Riehen (Basel-Stadt)" 
replace gemeindename = "Riehen" if canton == "BS" & city == "in Riehen-Basel" 
replace gemeindename = "Riehen" if canton == "BS" & city == "in Riehen/Basel" 
replace gemeindename = "Pratval" if canton == "GR" & city == "in Rietberg" // Rietberg ist vermutlich Teil von Pratval.
replace gemeindename = "Neftenbach" if canton == "ZH" & city == "in Riethof-Neftenbach" // Riethof ist Teil von Neftenbach.
replace gemeindename = "Guggisberg" if canton == "BE" & city == "in Riffenmatt" // Riffenmatt ist Teil von Guggisberg.
replace gemeindename = "Guggisberg" if canton == "BE" & city == "in Riffenmatt/Schwendi" // Riffenmatt und Schwendi sind Teil von Guggisberg.
replace gemeindename = "Riggisberg" if canton == "BE" & city == "in Riggisberg" 
replace gemeindename = "Zell (ZH)" if canton == "ZH" & city == "in Rikon" // Rikon ist Teil von Zell.
replace gemeindename = "Zell (ZH)" if canton == "ZH" & city == "in Rikon-Zell" // Rikon ist Teil von Zell.
replace gemeindename = "Zell (ZH)" if canton == "ZH" & city == "in Rikon/Zell" // Rikon ist Teil von Zell.
replace gemeindename = "Ringgenberg (BE)" if canton == "BE" & city == "in Ringgenberg" 
replace gemeindename = "Riniken" if canton == "AG" & city == "in Riniken" 
replace gemeindename = "Roches (BE)" if canton == "BE" & city == "in Roches" 
replace gemeindename = "Roggwil (BE)" if canton == "BE" & city == "in Roggwil" 
replace gemeindename = "Roggwil (TG)" if canton == "TG" & city == "in Roggwil" 
replace gemeindename = "Rohr (AG)" if canton == "AG" & city == "in Rohr" 
replace gemeindename = "Rohrbach" if canton == "BE" & city == "in Rohrbach" 
replace gemeindename = "Schönholzerswilen" if canton == "TG" & city == "in Rohren-Toos" // Toos ist Teil von Schönholzerswilen.
replace gemeindename = "Schönholzerswilen" if canton == "TG" & city == "in Rohren-Toss" // Toos ist Teil von Schönholzerswilen.
replace gemeindename = "Rorschach" if canton == "SG" & city == "in Rohrschach" 
replace gemeindename = "Romanshorn" if canton == "TG" & city == "in Romanshorn" 
replace gemeindename = "Romoos" if canton == "LU" & city == "in Romoos" 
replace gemeindename = "Root" if canton == "LU" & city == "in Root" 
replace gemeindename = "Rorbas" if canton == "ZH" & city == "in Rorbas" 
replace gemeindename = "Rorschach" if canton == "SG" & city == "in Rorschach" 
replace gemeindename = "Diemtigen" if canton == "BE" & city == "in Rothbad-Diemtigen" // Rothbad ist Teil von Diemtigen.
replace gemeindename = "Rothenbrunnen" if canton == "GR" & city == "in Rothenbrunnen" 
replace gemeindename = "Rothenburg" if canton == "LU" & city == "in Rothenburg" 
replace gemeindename = "Rothenthurm" if canton == "SZ" & city == "in Rothenthurm" 
replace gemeindename = "Rothrist" if canton == "AG" & city == "in Rothrist" 
replace gemeindename = "Roveredo (GR)" if canton == "GR" & city == "in Roveredo" 
replace gemeindename = "Trüllikon" if canton == "ZH" & city == "in Rudolfingen" // Rudolfingen ist Teil von Trüllikon.
replace gemeindename = "Russikon" if canton == "ZH" & city == "in Russikon" 
replace gemeindename = "Buchs (SG)" if canton == "SG" & city == "in Räfis-Buchs" // Räfis ist Teil von Buchs.
replace gemeindename = "Zell (ZH)" if canton == "ZH" & city == "in Rämismühle-Zell" // Rämismühle ist Teil von Zell.
replace gemeindename = "Elsau" if canton == "ZH" & city == "in Räterschen (Zürich)" // Räterschen ist Teil von Elsau.
replace gemeindename = "Elsau" if canton == "ZH" & city == "in Räterschen-Elsau" // Räterschen ist Teil von Elsau.
replace gemeindename = "Röschenz" if canton == "BE" & city == "in Röschenz" 
replace gemeindename = "Rüderswil" if canton == "BE" & city == "in Rüderswil" 
replace gemeindename = "Rüegsau" if canton == "BE" & city == "in Rüegsau" 
replace gemeindename = "Rüegsau" if canton == "BE" & city == "in Rüegsau-Schachen" // Rüegsauschachen ist Teil von Rüegsau.
replace gemeindename = "Rüegsau" if canton == "BE" & city == "in Rüegsauschachen" // Rüegsauschachen ist Teil von Rüegsau.
replace gemeindename = "Rüegsau" if canton == "BE" & city == "in Rüegsbach" // Rüegsbach ist Teil von Rüegsau.
replace gemeindename = "Worb" if canton == "BE" & city == "in Rüfenacht" // Rüfenacht ist Teil von Worb.
replace gemeindename = "Rüschlikon" if canton == "ZH" & city == "in Rüschliko" 
replace gemeindename = "Rüschlikon" if canton == "BE" & city == "in Rüschlikon" // Kontrolle notwendig!
replace gemeindename = "Rüschlikon" if canton == "SG" & city == "in Rüschlikon" // Kontrolle notwendig!
replace gemeindename = "Rüschlikon" if canton == "ZH" & city == "in Rüschlikon" 
replace gemeindename = "Rüschlikon" if canton == "BE" & city == "in Rüschlikon (Zürich)" // Kontrolle notwendig!
replace gemeindename = "Rüti (ZH)" if canton == "ZH" & city == "in Rüt" 
replace gemeindename = "Rüthi (Rheintal)" if canton == "SG" & city == "in Rüthi (St Gallen)" 
replace gemeindename = "Rüthi (Rheintal)" if canton == "SG" & city == "in Rüthi (St. Gallen)" 
replace gemeindename = "Rüti (GL)" if canton == "GL" & city == "in Rüti" 
replace gemeindename = "Rüti (ZH)" if canton == "ZH" & city == "in Rüti" 
replace gemeindename = "Rüti (ZH)" if canton == "ZH" & city == "in Rüti (Zh)" 
replace gemeindename = "Rüti (ZH)" if canton == "ZH" & city == "in Rüti (Zürich)" 
replace gemeindename = "Rüti bei Büren" if canton == "BE" & city == "in Rüti b. B." 
replace gemeindename = "Rüti bei Büren" if canton == "BE" & city == "in Rüti bei Bern" // Kontrolle notwendig! Rüti bei Bern in Bundesblatt nicht gefunden.
replace gemeindename = "Rüti bei Büren" if canton == "BE" & city == "in Rüti bei Büren" 
replace gemeindename = "Rüti (ZH)" if canton == "ZH" & city == "in Rüti, Kt. Zürich" 
replace gemeindename = "Winkel" if canton == "ZH" & city == "in Rüti-Winkel" // Rüti ist Teil von Winkel.
replace gemeindename = "Saanen" if canton == "BE" & city == "in Saanen" 
replace gemeindename = "Sachseln" if canton == "OW" & city == "in Sachseln" 
replace gemeindename = "Safenwil" if canton == "AG" & city == "in Safenwil" 
replace gemeindename = "Bauma" if canton == "ZH" & city == "in Saland-Bauma" // Saland ist Teil von Bauma.
replace gemeindename = "Salenstein" if canton == "TG" & city == "in Salenstein" 
replace gemeindename = "Sennwald" if canton == "SG" & city == "in Salez" // Salez ist Teil von Sennwald.
replace gemeindename = "Sennwald" if canton == "SG" & city == "in Salez-Sennwald" // Salez ist Teil von Sennwald.
replace gemeindename = "Salouf" if canton == "GR" & city == "in Salouf" 
replace gemeindename = "Samedan" if canton == "GR" & city == "in Samaden" 
replace gemeindename = "Sargans" if canton == "SG" & city == "in Sargans" 
replace gemeindename = "Sarn" if canton == "GR" & city == "in Sarn" 
replace gemeindename = "Sarnen" if canton == "OW" & city == "in Sarnen" 
replace gemeindename = "Schaffhausen" if canton == "SH" & city == "in Schaff hausen" 
replace gemeindename = "Schaffhausen" if canton == "SH" & city == "in Schaffhausen" 
replace gemeindename = "Schafisheim" if canton == "AG" & city == "in Schafisheim" 
replace gemeindename = "Schänis" if canton == "SG" & city == "in Schanis" 
replace gemeindename = "Reichenbach im Kandertal" if canton == "BE" & city == "in Scharnachthal" // Scharnachtal ist Teil von Reichenbach im Kandertal.
replace gemeindename = "Schiers" if canton == "GR" & city == "in Schiers" 
replace gemeindename = "Schinznach Dorf" if canton == "AG" & city == "in Schinznach" // Kontrolle notwendig! 2 Schinznach in AG (Schinznach Bad und Schinznach Dorf). Betrifft vermmutlich nur Otto Rinicker, sonst Schinznach Dorf.
replace gemeindename = "Schinznach Dorf" if canton == "AG" & city == "in Schinznach-Dorf" 
replace gemeindename = "Hombrechtikon" if canton == "ZH" & city == "in Schirmensee-Feldbach" // Feldbach und Oberschirmensee sind Teile von Hombrechtikon.
replace gemeindename = "Schlatt" if canton == "ZH" & city == "in Schlatt" 
replace gemeindename = "Schlieren" if canton == "ZH" & city == "in Schlieren" 
replace gemeindename = "Schlossrued" if canton == "AG" & city == "in Schlossrued" 
replace gemeindename = "Schmerikon" if canton == "SG" & city == "in Schmerikon" 
replace gemeindename = "Schmiedrued" if canton == "AG" & city == "in Schmidrued" 
replace gemeindename = "Schmiedrued" if canton == "AG" & city == "in Schmiedrued" 
replace gemeindename = "Schnottwil" if canton == "SO" & city == "in Schnottwil" 
replace gemeindename = "Schwadernau" if canton == "BE" & city == "in Schwadernau" 
replace gemeindename = "Schwanden (GL)" if canton == "GL" & city == "in Schwanden" 
replace gemeindename = "Wahlern" if canton == "BE" & city == "in Schwarzenburg" // Schwarzenburg war Teil von Wahlern (heute Schwarzenburg).
replace gemeindename = "Unterlangenegg" if canton == "BE" & city == "in Schwarzenegg" // Schwarzenegg ist Teil von Unterlangenegg.
replace gemeindename = "Schwerzenbach" if canton == "ZH" & city == "in Schwerzenbach" 
replace gemeindename = "Schwyz" if canton == "SZ" & city == "in Schwyz" 
replace gemeindename = "Schöftland" if canton == "AG" & city == "in Schöftland" 
replace gemeindename = "Schönenberg (ZH)" if canton == "ZH" & city == "in Schönenberg" 
replace gemeindename = "Schönenwerd" if canton == "SO" & city == "in Schönenwerd" 
replace gemeindename = "Schötz" if canton == "LU" & city == "in Schötz" 
replace gemeindename = "Schüpfen" if canton == "BE" & city == "in Schüpfen" 
replace gemeindename = "Scuol/Schuls" if canton == "GR" & city == "in Scuol/Schuls" 
replace gemeindename = "Zürich" if canton == "ZH" & city == "in Seebach" // Seebach ist Teil von Zürich.
replace gemeindename = "Winterthur" if canton == "ZH" & city == "in Seen-Winterthur" // Seen ist ein Teil von Winterthur.
replace gemeindename = "Schwyz" if canton == "SZ" & city == "in Seewen (Schwyz)" // Seewen ist Teil von Schwyz,
replace gemeindename = "Seftigen" if canton == "BE" & city == "in Seftigen" 
replace gemeindename = "Selzach" if canton == "SO" & city == "in Selzach" 
replace gemeindename = "Semione" if canton == "TI" & city == "in Semione" 
replace gemeindename = "Sennwald" if canton == "SG" & city == "in Sennwald" 
replace gemeindename = "Seon" if canton == "AG" & city == "in Seon" 
replace gemeindename = "Klosters" if canton == "GR" & city == "in Serneus" // Serneus ist Teil von Klosters (heute Klosters-Serneus).
replace gemeindename = "Seuzach" if canton == "ZH" & city == "in Seuzach" 
replace gemeindename = "Sevelen" if canton == "SG" & city == "in Sevelen" 
replace gemeindename = "Schübelbach" if canton == "SZ" & city == "in Siebnen" // Kontrolle notwendig! Siebnen ist auf 3 Gemeinden verteilt. Bei 2 Kandidaten (Fritz Stähli und Johann Wattenhofer) werden Siebnen und Siebnen-Schübelbach erwähnt; bei 2 anderen (Josef Diethelm und Josef Kürzi) nur Siebnen.
replace gemeindename = "Schübelbach" if canton == "SZ" & city == "in Siebnen '" // Kontrolle notwendig! Siebnen ist auf 3 Gemeinden verteilt. Bei 2 Kandidaten (Fritz Stähli und Johann Wattenhofer) werden Siebnen und Siebnen-Schübelbach erwähnt; bei 2 anderen (Josef Diethelm und Josef Kürzi) nur Siebnen.
replace gemeindename = "Schübelbach" if canton == "SZ" & city == "in Siebnen-Schübelbach" 
replace gemeindename = "Sils im Domleschg" if canton == "GR" & city == "in Sils i. D." 
replace gemeindename = "Sils im Engadin/Segl" if canton == "GR" & city == "in Sils i. E. / Segl-Fex" 
replace gemeindename = "Sils im Engadin/Segl" if canton == "GR" & city == "in Sils i. E. /Segl-Fex" 
replace gemeindename = "Sils im Engadin/Segl" if canton == "GR" & city == "in Sils-Fex" 
replace gemeindename = "Sirnach" if canton == "TG" & city == "in Sirnach" 
replace gemeindename = "Sissach" if canton == "BL" & city == "in Sissach" 
replace gemeindename = "Sisseln" if canton == "AG" & city == "in Sisseln" 
replace gemeindename = "Sitterdorf" if canton == "TG" & city == "in Sitterdorf" 
replace gemeindename = "Solothurn" if canton == "SO" & city == "in Solothurn" 
replace gemeindename = "Sommeri" if canton == "TG" & city == "in Sommeri" 
replace gemeindename = "Speicher" if canton == "AR" & city == "in Speicher" 
replace gemeindename = "Köniz" if canton == "BE" & city == "in Spiegel (Köniz)" // Spiegel ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "in Spiegel, Köniz" // Spiegel ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "in Spiegel/Bern" // Spiegel ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "in Spiegel/Köniz" // Spiegel ist Teil von Köniz.
replace gemeindename = "Spiez" if canton == "BE" & city == "in Spiez" 
replace gemeindename = "Aarberg" if canton == "BE" & city == "in Spins" // Spins ist Teil von Aarberg.
replace gemeindename = "Aarberg" if canton == "BE" & city == "in Spins b. Aarberg" // Spins ist Teil von Aarberg.
replace gemeindename = "Aarberg" if canton == "BE" & city == "in Spins bei Aarberg" // Spins ist Teil von Aarberg.
replace gemeindename = "Aarberg" if canton == "BE" & city == "in Spins-Aarberg" // Spins ist Teil von Aarberg.
replace gemeindename = "St. Gallen" if canton == "BE" & city == "in St Gallen" // Kontrolle notwendig!
replace gemeindename = "St. Gallen" if canton == "SG" & city == "in St Gallen" 
replace gemeindename = "St. Gallen" if canton == "SG" & city == "in St Gallen-Bruggen" // Bruggen ist Teil von St. Gallen.
replace gemeindename = "St. Gallen" if canton == "SG" & city == "in St Galleu" 
replace gemeindename = "St. Moritz" if canton == "GR" & city == "in St Moritz" 
replace gemeindename = "St. Gallen" if canton == "AG" & city == "in St. Gallen" // Kontrolle notwendig!
replace gemeindename = "St. Gallen" if canton == "BE" & city == "in St. Gallen" // Kontrolle notwendig!
replace gemeindename = "St. Gallen" if canton == "SG" & city == "in St. Gallen" 
replace gemeindename = "St. Gallen" if canton == "SG" & city == "in St. Gallen-Winkeln" // Winkeln ist Teil von St. Gallen.
replace gemeindename = "Saint-Imier" if canton == "BE" & city == "in St. Immer" // Sankt Immer ist dt. Name für Saint-Imier.
replace gemeindename = "St. Margrethen" if canton == "SG" & city == "in St. Margrethen" 
replace gemeindename = "St. Margrethen" if canton == "SG" & city == "in St. Margrethen (St. Gallen)" 
replace gemeindename = "St. Moritz" if canton == "GR" & city == "in St. Moritz" 
replace gemeindename = "St. Stephan" if canton == "BE" & city == "in St. Stephan" 
replace gemeindename = "Thal" if canton == "SG" & city == "in Staad" // Staad ist Teil von Thal.
replace gemeindename = "Staffelbach" if canton == "AG" & city == "in Staffelbach" 
replace gemeindename = "Stans" if canton == "NW" & city == "in Stans" 
replace gemeindename = "Starrkirch-Wil" if canton == "SO" & city == "in Starrkirch-Wil" 
replace gemeindename = "Staufen" if canton == "AG" & city == "in Staufen" 
replace gemeindename = "Steckborn" if canton == "TG" & city == "in Steckborn" 
replace gemeindename = "Steffisburg" if canton == "BE" & city == "in Steffisburg" 
replace gemeindename = "Steffisburg" if canton == "BE" & city == "in Steffisburg, Schwäbis" // Schwäbis ist Teil von Steffisburg.
replace gemeindename = "Stein (AG)" if canton == "AG" & city == "in Stein" 
replace gemeindename = "Stein am Rhein" if canton == "SH" & city == "in Stein am Rhein" 
replace gemeindename = "Steinach" if canton == "SG" & city == "in Steinach" 
replace gemeindename = "Steinmaur" if canton == "ZH" & city == "in Steinmaur" 
replace gemeindename = "Stetten (AG)" if canton == "AG" & city == "in Stetten" 
replace gemeindename = "Stettlen" if canton == "BE" & city == "in Stettlen" 
replace gemeindename = "Strengelbach" if canton == "AG" & city == "in Strengelbach" 
replace gemeindename = "Stäfa" if canton == "ZH" & city == "in Stäfa" 
replace gemeindename = "Suhr" if canton == "AG" & city == "in Suhr" 
replace gemeindename = "Sulz (AG)" if canton == "AG" & city == "in Sulz" 
replace gemeindename = "Sulz (LU)" if canton == "LU" & city == "in Sulz" 
replace gemeindename = "Rickenbach (ZH)" if canton == "ZH" & city == "in Sulz-Rickenbach" // Sulz ist Teil von Rickenbach.
replace gemeindename = "Sumiswald" if canton == "BE" & city == "in Sumiswald" 
replace gemeindename = "Surava" if canton == "GR" & city == "in Surava" 
replace gemeindename = "Sursee" if canton == "LU" & city == "in Sursee" 
replace gemeindename = "Lindau" if canton == "ZH" & city == "in Tagelswangen" // Tagelswangen ist Teil von Lindau.
replace gemeindename = "Dürnten" if canton == "ZH" & city == "in Tann-Dürnten" // Tann ist Teil von Dürnten.
replace gemeindename = "Dürnten" if canton == "ZH" & city == "in Tann-Rüti" // Kontrolle notwendig! Dürnten und Rüti sind 2 eigenständige Gemeinden. Betrifft vermutlich nur Max Hertig, sonst Dürnten, und Kurt Wick.
replace gemeindename = "Tegerfelden" if canton == "AG" & city == "in Tegerfelden" 
replace gemeindename = "Tenniken" if canton == "BL" & city == "in Tenniken" 
replace gemeindename = "Teufen (AR)" if canton == "AR" & city == "in Teufen" 
replace gemeindename = "Freienstein-Teufen" if canton == "ZH" & city == "in Teufen" // Teufen ist Teil von Freienstein-Teufen.
replace gemeindename = "Freienstein-Teufen" if canton == "ZH" & city == "in Teufen-Freienstein" 
replace gemeindename = "Freienstein-Teufen" if canton == "ZH" & city == "in Teufen/Freienstein" 
replace gemeindename = "Teufenthal (AG)" if canton == "AG" & city == "in Teufenthal" 
replace gemeindename = "Thalwil" if canton == "ZH" & city == "in Thalwil" 
replace gemeindename = "Thalwil" if canton == "ZH" & city == "in Thalwil-Gattikon" // Gattikon ist Teil von Thalwil.
replace gemeindename = "Therwil" if canton == "BL" & city == "in Therwil" 
replace gemeindename = "Thierachern" if canton == "BE" & city == "in Thierachern" 
replace gemeindename = "Krauchthal" if canton == "BE" & city == "in Thorberg" // Thorberg ist Teil von Krauchthal.
replace gemeindename = "Thun" if canton == "BE" & city == "in Thun" 
replace gemeindename = "Thusis" if canton == "GR" & city == "in Thusis" 
replace gemeindename = "Trachselwald" if canton == "BE" & city == "in Trachselwald" 
replace gemeindename = "Tramelan" if canton == "BE" & city == "in Tramelan-dessus" // Tramelan-Dessus ist Teil von Tramelan.
replace gemeindename = "Treiten" if canton == "BE" & city == "in Treiten" 
replace gemeindename = "Trimbach" if canton == "SO" & city == "in Trimbach" 
replace gemeindename = "Trimmis" if canton == "GR" & city == "in Trimmis" 
replace gemeindename = "Wartau" if canton == "SG" & city == "in Trubbach" // Trübbach ist Teil von Wartau.
replace gemeindename = "Trub" if canton == "BE" & city == "in Trüb" 
replace gemeindename = "Trüllikon" if canton == "ZH" & city == "in Trüllikon" 
replace gemeindename = "Lauterbrunnen" if canton == "BE" & city == "in Trümmelbach-Kleine Scheidegg" // Kontrolle notwendig! Trümmelbach und Kleine Scheidegg sind vermutlich Teil von Lauterbrunnen.
replace gemeindename = "Lauterbrunnen" if canton == "BE" & city == "in Trümmelbach-Scheidegg" // Kontrolle notwendig! Trümmelbach und Kleine Scheidegg sind vermutlich Teil von Lauterbrunnen.
replace gemeindename = "Lauterbrunnen" if canton == "BE" & city == "in Trümmelbach/Scheidegg" // Kontrolle notwendig! Trümmelbach und Kleine Scheidegg sind vermutlich Teil von Lauterbrunnen.
replace gemeindename = "Turbenthal" if canton == "ZH" & city == "in Turbenthal" 
replace gemeindename = "Turgi" if canton == "AG" & city == "in Turg" 
replace gemeindename = "Turgi" if canton == "AG" & city == "in Turgi" 
replace gemeindename = "Twann" if canton == "BE" & city == "in Twann" 
replace gemeindename = "Tägerwilen" if canton == "TG" & city == "in Tägerwilen" 
replace gemeindename = "Tübach" if canton == "SG" & city == "in Tübach" 
replace gemeindename = "Tüscherz-Alfermée" if canton == "BE" & city == "in Tüscherz" 
replace gemeindename = "Uebeschi" if canton == "BE" & city == "in Uebeschi" 
replace gemeindename = "Stäfa" if canton == "ZH" & city == "in Uerikon-Stäfa" // Uerikon ist Teil von Stäfa.
replace gemeindename = "Stäfa" if canton == "ZH" & city == "in Uerikon/Stäfa" // Uerikon ist Teil von Stäfa.
replace gemeindename = "Uetendorf" if canton == "BE" & city == "in Uetendorf" 
replace gemeindename = "Uetikon" if canton == "ZH" & city == "in Uetikon" 
replace gemeindename = "Uetikon" if canton == "ZH" & city == "in Uetikon a. S." 
replace gemeindename = "Uetikon" if canton == "ZH" & city == "in Uetikon a. See" 
replace gemeindename = "Wohlen bei Bern" if canton == "BE" & city == "in Uettligen" // Uettligen ist Teil von Wohlen bei Bern.
replace gemeindename = "Uffikon" if canton == "LU" & city == "in Uffikon" 
replace gemeindename = "Uitikon" if canton == "ZH" & city == "in Uitikon am Albis" 
replace gemeindename = "Umiken" if canton == "AG" & city == "in Umiken" 
replace gemeindename = "Undervelier" if canton == "BE" & city == "in Undervelier" 
replace gemeindename = "Hittnau" if canton == "ZH" & city == "in Unter-Hittnau" // Under-Hitnau ist Teil von Hitnau.
replace gemeindename = "Unterstammheim" if canton == "ZH" & city == "in Unter-Stammheim" 
replace gemeindename = "Meiringen" if canton == "BE" & city == "in Unterbach/Meiringen" // Unterbachist Teil von Meiringen.
replace gemeindename = "Unterengstringen" if canton == "AG" & city == "in Unterengstringen" // Kontrolle notwendig! Möglich, da Kloster Fahr = Enklave von AG in Unterengstringen; siehe http://www.unterengstringen.ch/de/links/; 27.01.2017.
replace gemeindename = "Unterengstringen" if canton == "ZH" & city == "in Unterengstringen" 
replace gemeindename = "Unterentfelden" if canton == "AG" & city == "in Unterentfelden" 
replace gemeindename = "Hallau" if canton == "SH" & city == "in Unterhallau" // Früherer Name.
replace gemeindename = "Unterkulm" if canton == "AG" & city == "in Unterkulm" 
replace gemeindename = "Unterseen" if canton == "BE" & city == "in Unterseen" 
replace gemeindename = "Untersiggenthal" if canton == "AG" & city == "in Untersiggenthal" 
replace gemeindename = "Untervaz" if canton == "GR" & city == "in Untervaz" 
replace gemeindename = "Urdorf" if canton == "ZH" & city == "in Urdorf" 
replace gemeindename = "Urtenen" if canton == "BE" & city == "in Urtenen" 
replace gemeindename = "Urtenen" if canton == "BE" & city == "in Urtenen-Schönbühl" // Schönbühl ist Teil von Urtenen (heute Urtenen-Schönbühl).
replace gemeindename = "Uster" if canton == "ZH" & city == "in Uster" 
replace gemeindename = "Uttigen" if canton == "BE" & city == "in Uttigen" 
replace gemeindename = "Utzenstorf" if canton == "BE" & city == "in Utzenstorf" 
replace gemeindename = "Vechigen" if canton == "BE" & city == "in Utzigen-Boll" // Utzigen und Boll sind Teil von Vechigen.
replace gemeindename = "Uznach" if canton == "SG" & city == "in Uznach" 
replace gemeindename = "Uzwil" if canton == "SG" & city == "in Uzwil" 
replace gemeindename = "Vacallo" if canton == "TI" & city == "in Vacallo" 
replace gemeindename = "Winterthur" if canton == "SG" & city == "in Veltheim" // Kontrolle notwendig! Kein Veltheim in St. Gallen gefunden. Ist vermutlich Stadtkreis Veltheim in Winterthur. Betrifft vermutlich nur Hans Mosimann, sonst Winterthur.
replace gemeindename = "Worb" if canton == "BE" & city == "in Vielbringen (Worb)" // Vielbringen ist Teil von Worb.
replace gemeindename = "Worb" if canton == "BE" & city == "in Vielbringen b. Worb" // Vielbringen ist Teil von Worb.
replace gemeindename = "Viganello" if canton == "TI" & city == "in Viganello" 
replace gemeindename = "Vigens" if canton == "GR" & city == "in Vigens" 
replace gemeindename = "Villmergen" if canton == "AG" & city == "in Villmergen" 
replace gemeindename = "Vinelz" if canton == "BE" & city == "in Vinelz" 
replace gemeindename = "Vitznau" if canton == "LU" & city == "in Vitznau" 
replace gemeindename = "Vordemwald" if canton == "AG" & city == "in Vordemwald" 
replace gemeindename = "Pfäfers" if canton == "SG" & city == "in Vättis" // Vättis ist Teil von Pfäfers.
replace gemeindename = "Köniz" if canton == "BE" & city == "in Wabern" // Wabern ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "in Wabern b.B." // Wabern ist Teil von Köniz.
replace gemeindename = "Köniz" if canton == "BE" & city == "in Wabern-Bern" // Wabern ist Teil von Köniz.
replace gemeindename = "Wagenhausen" if canton == "TG" & city == "in Wagenhausen" 
replace gemeindename = "Wahlen" if canton == "BE" & city == "in Wahlen" 
replace gemeindename = "Wald (ZH)" if canton == "ZH" & city == "in Wald" 
replace gemeindename = "Wald (ZH)" if canton == "ZH" & city == "in Wald (Zch.)" 
replace gemeindename = "Wald (ZH)" if canton == "ZH" & city == "in Wald (Zürich )" 
replace gemeindename = "Wald (ZH)" if canton == "ZH" & city == "in Wald (Zürich)" 
replace gemeindename = "Waldenburg" if canton == "BL" & city == "in Waldenburg" 
replace gemeindename = "Waldenburg" if canton == "BS" & city == "in Waldenburg" // Kontrolle notwendig!
replace gemeindename = "Waldkirch" if canton == "SG" & city == "in Waldkirch" 
replace gemeindename = "Walenstadt" if canton == "SG" & city == "in Walenstadt" 
replace gemeindename = "Walenstadt" if canton == "SG" & city == "in Walenstadtberg" // Walenstadtberg gehört zu Walenstadt.
replace gemeindename = "Walenstadt" if canton == "SG" & city == "in Wallenstadt" 
replace gemeindename = "Wallisellen" if canton == "ZH" & city == "in Wallisellen" 
replace gemeindename = "Walliswil bei Wangen" if canton == "BE" & city == "in Walliswil-Wangen" 
replace gemeindename = "Waltenschwil" if canton == "AG" & city == "in Waltenschwil" 
replace gemeindename = "Walzenhausen" if canton == "AR" & city == "in Walzenhausen" 
replace gemeindename = "Wangen (ZH)" if canton == "ZH" & city == "in Wangen" 
replace gemeindename = "Wangen (ZH)" if canton == "ZH" & city == "in Wangen (Zürich)" 
replace gemeindename = "Wangen an der Aare" if canton == "BE" & city == "in Wangen a. A." 
replace gemeindename = "Wangen an der Aare" if canton == "BE" & city == "in Wangen a. d. A." 
replace gemeindename = "Wangen (ZH)" if canton == "ZH" & city == "in Wangen b. Dübendorf" 
replace gemeindename = "Wangen (ZH)" if canton == "ZH" & city == "in Wangen bei Dübendorf" 
replace gemeindename = "Wangenried" if canton == "BE" & city == "in Wangenried" 
replace gemeindename = "Warth" if canton == "TG" & city == "in Warth bei Märwil" 
replace gemeindename = "Hasliberg" if canton == "BE" & city == "in Wasserwendi bei Hasleberg" // Wasserwendi ist Teil von Hasliberg.
replace gemeindename = "Hasliberg" if canton == "BE" & city == "in Wasserwendi-Hasleberg" // Wasserwendi ist Teil von Hasliberg.
replace gemeindename = "Hasliberg" if canton == "BE" & city == "in Wasserwendi/Hasleberg" // Wasserwendi ist Teil von Hasliberg.
replace gemeindename = "Regensdorf" if canton == "ZH" & city == "in Watt-Regensdorf" // Watt ist Tail von Regensdorf.
replace gemeindename = "Wattenwil" if canton == "BE" & city == "in Wattenwil" 
replace gemeindename = "Wattwil" if canton == "SG" & city == "in Wattwil" 
replace gemeindename = "Wattwil" if canton == "SG" & city == "in Wattwil " 
replace gemeindename = "Weesen" if canton == "SG" & city == "in Weesen" 
replace gemeindename = "Wegenstetten" if canton == "AG" & city == "in Wegenstetten" 
replace gemeindename = "Weggis" if canton == "LU" & city == "in Weggis" 
replace gemeindename = "Weiach" if canton == "ZH" & city == "in Weiach" 
replace gemeindename = "Affoltern im Emmental" if canton == "BE" & city == "in Weier i. E." // Weier im Emmental ist Teil von Affoltern im Emmental.
replace gemeindename = "Affoltern im Emmental" if canton == "BE" & city == "in Weier/Affoltern" // Weier im Emmental ist Teil von Affoltern im Emmental.
replace gemeindename = "Weinfelden" if canton == "TG" & city == "in Weinfelden" 
replace gemeindename = "Weiningen (ZH)" if canton == "ZH" & city == "in Weiningen" 
replace gemeindename = "Weiningen (ZH)" if canton == "ZH" & city == "in Weiningen (Zürich)" 
replace gemeindename = "Boltigen" if canton == "BE" & city == "in Weissenbach i. S." // Weissenbach ist Teil von Boltigen.
replace gemeindename = "Boltigen" if canton == "BE" & city == "in Weissenbach/Boltigen" // Weissenbach ist Teil von Boltigen.
replace gemeindename = "Weisslingen" if canton == "ZH" & city == "in Weisslingen" 
replace gemeindename = "Wartau" if canton == "SG" & city == "in Weite-Wartau" // Weite ist Teil von Wartau.
replace gemeindename = "Lauterbrunnen" if canton == "BE" & city == "in Wengen" // Wengen ist Teil  von Lauterbrunnen.
replace gemeindename = "Wenslingen" if canton == "BL" & city == "in Wenslingen" 
replace gemeindename = "Grabs" if canton == "SG" & city == "in Werdenberg" // Werdenberg ist Teil von Grabs.
replace gemeindename = "Wettingen" if canton == "AG" & city == "in Wettingen" 
replace gemeindename = "Wetzikon (ZH)" if canton == "ZH" & city == "in Wetzikon" 
replace gemeindename = "Widnau" if canton == "SG" & city == "in Widnau" 
replace gemeindename = "Wiedlisbach" if canton == "BE" & city == "in Wiedlisbach" 
replace gemeindename = "Uzwil" if canton == "SG" & city == "in Wiesbühl-Uzwil" // Wiesbühl ist vermutlich ein Bauernhof in Uzwil.
replace gemeindename = "Wiggiswil" if canton == "BE" & city == "in Wiggiswil" 
replace gemeindename = "Wiggiswil" if canton == "BE" & city == "in Wiggiswil-Münchenbuchsee" // Kontrolle notwendig! Wiggiswil und Münchenbuchsee sind 2 eigenständige Gemeinden. Betrifft vermutlich nur Otto Häberli, sonste einmal nur Wiggiswil und einmal nur Münchenbuchsee. --> Wiggiswil zugeordnet.
replace gemeindename = "Wikon" if canton == "LU" & city == "in Wikon" 
replace gemeindename = "Wil (SG)" if canton == "SG" & city == "in Wil" 
replace gemeindename = "Wil (SG)" if canton == "SG" & city == "in Wil (St. G.)" 
replace gemeindename = "Wil (SG)" if canton == "SG" & city == "in Wil (St. Gallen)" 
replace gemeindename = "Wil (SG)" if canton == "SG" & city == "in Wil-St Gallen" 
replace gemeindename = "Wila" if canton == "ZH" & city == "in Wila" 
replace gemeindename = "Wilderswil" if canton == "BE" & city == "in Wilderswil" 
replace gemeindename = "Wildhaus" if canton == "SG" & city == "in Wildhaus" 
replace gemeindename = "Sarnen" if canton == "OW" & city == "in Wilen-Sarnen" // Wilen ist Teil von Sarnen.
replace gemeindename = "Sarnen" if canton == "OW" & city == "in Wilen-Sarner" 
replace gemeindename = "Seedorf (BE)" if canton == "BE" & city == "in Wiler-Seedorf" // Wiler ist Teil von Seedorf.
replace gemeindename = "Einsiedeln" if canton == "SZ" & city == "in Willerzell" // Willerzell ist Teil von Einsiedeln.
replace gemeindename = "Willisau Land" if canton == "LU" & city == "in Willisau" // Kontrolle notwendig! 2 Willisau (Willisau Stadt und Willisau Land). Betrifft 3 Kandidaten: Franz Jos. Kurmann, sonst Willisau-Land, sowie Erwin Muff und Fritz Grüter.
replace gemeindename = "Windisch" if canton == "AG" & city == "in Windisch" 
replace gemeindename = "Winkel" if canton == "ZH" & city == "in Winkel-Rüti" // Rüti ist Teil von Winkel.
replace gemeindename = "St. Gallen" if canton == "SG" & city == "in Winkeln-St Gallen" // Winkeln ist Teil von St. Gallen.
replace gemeindename = "Lindau" if canton == "ZH" & city == "in Winterberg/Lindau" // Winterberg ist Teil von Lindau.
replace gemeindename = "Winterthur" if canton == "ZH" & city == "in Winterthur" 
replace gemeindename = "Winterthur" if canton == "ZH" & city == "in Winterthur-Seen" // Seen ist ein Teil von Winterthur.
replace gemeindename = "Winterthur" if canton == "ZH" & city == "in Winterthur-Veitheim" // Veltheim ist Teil von Winterthur.
replace gemeindename = "Winterthur" if canton == "ZH" & city == "in Winterthur-Wülflingen" // Wülflingen ist Teil von Winterthur.
replace gemeindename = "Zürich" if canton == "ZH" & city == "in Witikon" // Witikon ist Teil von Zürich (seit 1934).
replace gemeindename = "Wittenbach" if canton == "SG" & city == "in Wittenbach" 
replace gemeindename = "Wittnau" if canton == "AG" & city == "in Wittnau" 
replace gemeindename = "Wohlen (AG)" if canton == "AG" & city == "in Wohlen" 
replace gemeindename = "Wolfhalden" if canton == "AR" & city == "in Wolfhalden" 
replace gemeindename = "Wolfisberg" if canton == "BE" & city == "in Wolfisberg" 
replace gemeindename = "Wölflinswil" if canton == "AG" & city == "in Wolflinswi" 
replace gemeindename = "Wollerau" if canton == "SZ" & city == "in Wollerau" 
replace gemeindename = "Worb" if canton == "BE" & city == "in Worb" 
replace gemeindename = "Worb" if canton == "BE" & city == "in Worb-Rüfenacht" // Rüfenacht ist Teil von Worb.
replace gemeindename = "Bolligen" if canton == "BE" & city == "in Worblaufen" // Worblaufen war Teil von Bolligen (heute Ittigen).
replace gemeindename = "Wynau" if canton == "BE" & city == "in Wynau" 
replace gemeindename = "Wynigen" if canton == "BE" & city == "in Wynigen" 
replace gemeindename = "Wädenswil" if canton == "ZH" & city == "in Wädenswil" 
replace gemeindename = "Wängi" if canton == "TG" & city == "in Wängi" 
replace gemeindename = "Wölflinswil" if canton == "AG" & city == "in Wölflinswil" 
replace gemeindename = "Wölflinswil" if canton == "AG" & city == "in Wölflinwil" 
replace gemeindename = "Winterthur" if canton == "ZH" & city == "in Wülflingen-Winterthur" // Wülflingen ist Teil von Winterthur.
replace gemeindename = "Würenlingen" if canton == "AG" & city == "in Würenlingen" 
replace gemeindename = "Würenlos" if canton == "AG" & city == "in Würenlos" 
replace gemeindename = "Zell (LU)" if canton == "LU" & city == "in Zell" 
replace gemeindename = "Zell (ZH)" if canton == "ZH" & city == "in Zell" 
replace gemeindename = "Zernez" if canton == "GR" & city == "in Zernez" 
replace gemeindename = "Zetzwil" if canton == "AG" & city == "in Zetzwil" 
replace gemeindename = "Zimmerwald" if canton == "BE" & city == "in Zimmerwald" 
replace gemeindename = "Zizers" if canton == "GR" & city == "in Zizers" 
replace gemeindename = "Zofingen" if canton == "AG" & city == "in Zofingen" 
replace gemeindename = "Zollikon" if canton == "ZH" & city == "in Zolikerberg" // Zollikerberg ist Teil von Zollikon.
replace gemeindename = "Lauperswil" if canton == "BE" & city == "in Zollbrück" & name == "Geissbühler" // Kontrolle notwendig! Zollbrück ist Teil von Rüderswil und Lauperswil. Betrifft vermutlich drei Kandidaten: Franz Badertscher ohne weitere Angabe (jedoch Heimatort = Lauperswil), Fritz Geissbühler, sonst Lauperswil, Ernst Hirsbrunner, Vater in Rüderswil.
replace gemeindename = "Zollikon" if canton == "ZH" & city == "in Zollikerberg" // Zollikerberg ist Teil von Zollikon.
replace gemeindename = "Zollikofen" if canton == "BE" & city == "in Zollikofen" 
replace gemeindename = "Zollikofen" if canton == "BE" & city == "in Zollikofen *" 
replace gemeindename = "Zollikon" if canton == "ZH" & city == "in Zollikon" 
replace gemeindename = "Zuchwil" if canton == "SO" & city == "in Zuchwil" 
replace gemeindename = "Zug" if canton == "ZG" & city == "in Zug" 
replace gemeindename = "Zumikon" if canton == "ZH" & city == "in Zumikon" 
replace gemeindename = "Zunzgen" if canton == "BL" & city == "in Zunzgen" 
replace gemeindename = "Zuoz" if canton == "GR" & city == "in Zuoz" 
replace gemeindename = "Zürich" if canton == "TI" & city == "in Zurigo" // Kontrolle notwendig! Zurigo ist it. Name für Zürich.
replace gemeindename = "Zurzach" if canton == "AG" & city == "in Zurzach" 
replace gemeindename = "Zuzwil (SG)" if canton == "SG" & city == "in Zuzwil" 
replace gemeindename = "Glattfelden" if canton == "ZH" & city == "in Zweidlen-Glattfelden" // Zweidlen ist Teil von Glattfelden.
replace gemeindename = "Zweisimmen" if canton == "BE" & city == "in Zweisimmen" 
replace gemeindename = "Zwingen" if canton == "BE" & city == "in Zwingen" 
replace gemeindename = "Zürich" if canton == "ZH" & city == "in Zürch" 
replace gemeindename = "Zürich" if canton == "AG" & city == "in Zürich" // Kontrolle notwendig!
replace gemeindename = "Zürich" if canton == "BE" & city == "in Zürich" // Kontrolle notwendig!
replace gemeindename = "Zürich" if canton == "LU" & city == "in Zürich" // Kontrolle notwendig!
replace gemeindename = "Zürich" if canton == "SG" & city == "in Zürich" // Kontrolle notwendig!
replace gemeindename = "Zürich" if canton == "ZH" & city == "in Zürich" 
replace gemeindename = "Zürich" if canton == "ZH" & city == "in Zürich-Oerlikon" // Oerlikon ist Teil von Zürich.
replace gemeindename = "Zürich" if canton == "ZH" & city == "in Zürich-Seebach" // Seebach ist Teil von Zürich.
replace gemeindename = "Brugg" if canton == "AG" & city == "in brugg" 
replace gemeindename = "Klingnau" if canton == "AG" & city == "in Klingnau" 
replace gemeindename = "Langenthal" if canton == "BE" & city == "inLangenthal" 
replace gemeindename = "Unterengstringen" if canton == "ZH" & city == "inUnterengstringen" 
replace gemeindename = "L'Abbaye" if canton == "VD" & city == "l'Abbaye" 
replace gemeindename = "Bourrignon" if canton == "BE" & city == "la Bürgisberg" // La Bürgisberg ist dt. Name für Bourrignon.
replace gemeindename = "Diesse" if canton == "BE" & city == "la Montagne de Diesse" // Kontrolle notwendig! La Montagne de Diesse (dt. Tessenberg) umfasst 4 (evtl. 5) Gemeinden: Diesse, Lamboing, Prêles (3 ehemalige Gemeinden sind heute Plateau de Diesse) und Nods (sowie evtl. Lignières). Betrifft vermutlich nur Georges Luterbacher. Arbiträre Zuordnung zu Diesse.
replace gemeindename = "La Rippe" if canton == "VD" & city == "la Rippe" 
replace gemeindename = "Montmelon" if canton == "BE" & city == "la Saigne-du-Milieu (Montmelon)" 
