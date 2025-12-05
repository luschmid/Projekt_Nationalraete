clear
cap log close
set more off
version 14

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\1 Data\1_data"
global path "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"


*** Prepare information on party votes


use "${path}\16_alliances\Listenverbindungen_1931-2015.dta", clear
gen liste_norm = Listenname
duplicates list canton year liste_norm
/*
  +-------------------------------------------------+
  | group:   obs:   canton   year        Listenname |
  |-------------------------------------------------|
  |      1   2167       NW   1943   kein Listenname |
  |      1   2168       NW   1943   kein Listenname |
  |      2   2193       OW   1943   kein Listenname |
  |      2   2194       OW   1943   kein Listenname |
  +-------------------------------------------------+
  
*/
replace liste_norm = "kein Listenname 1" if year == 1943 & canton == "NW" & Listennummer_num == 1 // the one party receiving a seat is now "kein Listenname 1"
replace liste_norm = "kein Listenname 2" if year == 1943 & canton == "NW" & Listennummer_num == 2 
replace liste_norm = "kein Listenname 1" if year == 1943 & canton == "OW" & Listennummer_num == 1 // the one party receiving a seat is now "kein Listenname 1"
replace liste_norm = "kein Listenname 2" if year == 1943 & canton == "OW" & Listennummer_num == 2
replace liste_norm = "kein Listenname" if Listenname == ""

sort canton year liste_norm
duplicates list canton year liste_norm

replace liste_norm = upper(liste_norm)
replace liste_norm = usubinstr(liste_norm, ".", " ", .)
replace liste_norm = usubinstr(liste_norm, ",", " ", .)
replace liste_norm = usubinstr(liste_norm, ";", " ", .)
replace liste_norm = usubinstr(liste_norm, ":", " ", .)
replace liste_norm = usubinstr(liste_norm, "-", " ", .)
replace liste_norm = usubinstr(liste_norm, "Ä", "A", .)
replace liste_norm = usubinstr(liste_norm, "Ä", "A", .)
replace liste_norm = usubinstr(liste_norm, "Ö", "O", .)
replace liste_norm = usubinstr(liste_norm, "Ö", "O", .)
replace liste_norm = usubinstr(liste_norm, "Ü", "U", .)
replace liste_norm = usubinstr(liste_norm, "Ü", "U", .)
replace liste_norm = usubinstr(liste_norm, "É", "E", .)
replace liste_norm = usubinstr(liste_norm, "É", "E", .)
replace liste_norm = usubinstr(liste_norm, "È", "E", .)
replace liste_norm = usubinstr(liste_norm, "Ê", "E", .)
replace liste_norm = usubinstr(liste_norm, "Â", "A", .)
replace liste_norm = usubinstr(liste_norm, "À", "A", .)
replace liste_norm = usubinstr(liste_norm, "ä", "A", .)
replace liste_norm = usubinstr(liste_norm, "ö", "O", .)
replace liste_norm = usubinstr(liste_norm, "ü", "U", .)
replace liste_norm = usubinstr(liste_norm, "é", "E", .)
replace liste_norm = usubinstr(liste_norm, "è", "E", .)
replace liste_norm = usubinstr(liste_norm, "ê", "E", .)
replace liste_norm = usubinstr(liste_norm, "à", "A", .)
replace liste_norm = usubinstr(liste_norm, "â", "A", .)
replace liste_norm = usubinstr(liste_norm, "  ", " ", .)
replace liste_norm = usubinstr(liste_norm, "«", "", .)
replace liste_norm = usubinstr(liste_norm, "»", "", .)
replace liste_norm = strtrim(liste_norm)

/*
merge 1:m canton year liste_norm using "${path2}nationalraete_1931_2015_normalizedListNames.dta"
br canton year liste_norm if _merge==1
br canton year liste_norm if _merge==2
*--- > see Excel "CompareUnmatchedParties.xlsx" with mismaches and code
*/
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "AG" & year == 1935 & liste_norm =="FREISINNIG DEMOKRATISCHE VOLKSPARTEI"
replace liste_norm = "SCHWEIZ BAUERNHEIMATBEWEGUNG JUNGBAUERN" if canton == "AG" & year == 1935 & liste_norm =="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG JUNGBAUERN"
replace liste_norm = "SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG (JUNGBAUERN)" if canton == "AG" & year == 1943 & liste_norm =="SCHWEIZ BAUERNHEIMATBEWEGUNG JUNGBAUERN"
replace liste_norm = "FREIE STIMMBURGERFUR DIE AUFHEBUNG DES STIMMZWANGS" if canton == "AG" & year == 1959 & liste_norm =="FREIE STIMMBURGER FUR DIE AUFHEBUNG DES STIMMZWANGES"
replace liste_norm = "FREISINNIG DEMOKRATISCHE VOLKSPARTEI UND JUNGLIBERALE BEWEGUNG" if canton == "AG" & year == 1963 & liste_norm =="FREISINNIG DEMOKRATISCHE VOLKSPARTEI UND JUNGLIHERALE BEWEGUNG"
replace liste_norm = "BAUERN  GEWERBE UND BURGERPARTEI (BGB MITTELSTANDSIISTE)" if canton == "AG" & year == 1967 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI (BGB MITTELSTANDSLISTE)"
replace liste_norm = "TEAM 67 LIBERALE AARGAUER FUR EINE MODERNE SCHWEIZ" if canton == "AG" & year == 1967 & liste_norm =="TEAM 1967 LIBERALE AARGAUER FUR EINE MODERNE SCHWEIZ"
replace liste_norm = "FEDERATION LIBERALE POPULAIRE JURASSIENNE" if canton == "BE" & year == 1931 & liste_norm =="JURASSISCHE LIBERALE PARTEI"
replace liste_norm = "PARTI DEMOCRATIQUE CATHOLIQUE DU CANTON DE BERNE" if canton == "BE" & year == 1931 & liste_norm =="KATHOLISCHE VOLKSPARTEI"
replace liste_norm = "BAUERN  GEWERBE UND BURGERPARTEI" if canton == "BE" & year == 1935 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI DES KANTONS BERN"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "BE" & year == 1935 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN UND JUNG LIBERALE BEWEGUNG"
replace liste_norm = "SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG JUNGBAUERN" if canton == "BE" & year == 1935 & liste_norm =="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "BE" & year == 1935 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS BERN"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "BE" & year == 1939 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN"
replace liste_norm = "KOMMUNISTISCHE PARTEI" if canton == "BE" & year == 1939 & liste_norm =="KOMMUNISTISCHE PARTEI DES KANTONS BERN"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "BE" & year == 1939 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS BERN"
replace liste_norm = "BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI EMMENTAL MITTELLAND OBERAARGAU SEELAND" if canton == "BE" & year == 1943 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI EMMENTAL MITTELLAND OBERAARGAU  SEELAND"
replace liste_norm = "FREISINNIG DEMOKRATISCHE VOLKSPARTEI OBERAARGAU/EMMENTAL" if canton == "BE" & year == 1943 & liste_norm =="FREISINNIG DEMOKRATISCHE VOLKSLISTE OBERAARGAU EMMENTAL"
replace liste_norm = "LANDESTEILVERBAND OBERLAND DER BERN BAUERN  GEWERBE UND BURGERPARTEI" if canton == "BE" & year == 1943 & liste_norm =="LANDESTEILVERBAND OBERLAND DER BAUERN  GEWERBE UND BURGERPARTEI"
replace liste_norm = "PARTI DEMOCRATIQUE CATHOLIQUE DU CANTON DE BERNE" if canton == "BE" & year == 1943 & liste_norm =="PARTI DEMOCRATIQUE CATHOLIQUE"
replace liste_norm = "SCHWEIZERISCHE BAUERN HEIMATBEWEGUNG" if canton == "BE" & year == 1943 & liste_norm =="SCHWEIZ BAUERNHEIMATBEWEGUNG"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI DES KANTONS BERN" if canton == "BE" & year == 1943 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm = "BERN BAUERN  GEWERBE UND BURGERPARTEI EMMENTAL MITTELLAND OBERAARGAU SEELAND" if canton == "BE" & year == 1947 & liste_norm =="LANDESTEILVERBAND EMMENTAL MITTELLAND OBERAARGAU SEELAND DER BERN BAUERN  GEWERBE UND BURGERPARTEI"
replace liste_norm = "BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI EMMENTAL JURA MITTELLAND OBERAARGAU SEELAND" if canton == "BE" & year == 1951 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI EMMENTAL JURA MITTELLAND OBERAARGAU SEELAND"
replace liste_norm = "BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI LANDESTEIL OBERLAND" if canton == "BE" & year == 1951 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI LANDESTEIL OBERLAND"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEILVERBAND OBERLAND" if canton == "BE" & year == 1951 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI LANDESTEIL VERBAND OBERLAND"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEILE EMMENTAL MITTELLAND OBERAARGAN SEELAND" if canton == "BE" & year == 1951 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI LANDESTEILE EMMENTAL MITTELLAND OBERAARGAU SEELAND"
replace liste_norm = "LIBERAL SOZIALISTISCHE PARTEI DES KANTONS BERN" if canton == "BE" & year == 1951 & liste_norm =="LIBERAL SOZIALISTISCHE PARTEI"
replace liste_norm = "PARTI SOCIALISTE JURASSIE" if canton == "BE" & year == 1951 & liste_norm =="PARTI SOCIALISTE JURASSIEN"
replace liste_norm = "BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI LANDESTEILVERBAND OBERLAND FREIE DEMOKRATISCHE MITTELSTANDSPARTEI" if canton == "BE" & year == 1955 & liste_norm =="BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI LANDESTEIL OBERLAND FREIE DEMOKRATISCHE MITTELSTANDSPARTEI"
replace liste_norm = "PARTI LIBERAL RADICAL JURASSIEN PARTI NATIONAL ROMAND DE BIENNE GROUPE RADICAL ROMAND DE BERNE" if canton == "BE" & year == 1959 & liste_norm =="PARTI LIBERAL RADICAL JURASSIEN"
replace liste_norm = "LANDESRING DER UNABHANGIGEN ALLIANCE DES INDEPENDANTS" if canton == "BE" & year == 1963 & liste_norm =="LANDESRING DER UNABHANGIGEN"
replace liste_norm = "PARTI LIBERAL RADICAL JURASSIEN FREISINNIGE PARTEI DES JURA" if canton == "BE" & year == 1963 & liste_norm =="PARTI LIBERAL RADICAL JURASSIEN (FREISINNIGE PARTEI DES JURA)"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI DES KANTONS BERN PARTI SOCIALISTE DU CANTON DE BERNE" if canton == "BE" & year == 1963 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS BERN"
replace liste_norm = "BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI FREIE DEMOKRATISCHE MITTELSTANDSPARTEI LANDESTEIL OBERLAND" if canton == "BE" & year == 1967 & liste_norm =="BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI FREIE DEMOKRATISCHE MITTELSTANDSPARTEI LANDESTEILVERBAND OBERLAND"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEIL BERN MITTELLAND" if canton == "BE" & year == 1967 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEIL MITTELLAND"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEILE SEELAND  LAUFENTAL" if canton == "BE" & year == 1967 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEILE SEELAND LAUFENTAL"
replace liste_norm = "KONSERVATIV CHRISTLICHSOZIALE VOLKSPARTEI DES KANTONS BERN" if canton == "BE" & year == 1967 & liste_norm =="KONSERVATIV CHRISTLICHSOZIALE VOLKSPARTEI"
replace liste_norm = "PARTI LIBERAL RADICAL JURASSIEN (FREISINNIGE PARTEI DES JURA) PARTI NATIONAL ROMAND DE BIENNE" if canton == "BE" & year == 1967 & liste_norm =="PARTI LIBERAL RADICAL JURASSIEN (FREISINNIGE PARTEI DES JURA)"
replace liste_norm = "BERNISCHE BGB OBERAARGAU (BAUERN  GEWERBE UND BARGERPARTEI)" if canton == "BE" & year == 1971 & liste_norm =="BERNISCHE BGB OBERAARGAU (BAUERN  GEWERBE UND BURGERPARTEI)"
replace liste_norm = "PARTI JURASSIEN DES PAYSANS ARTISANS ET BOURGEOIS(PAB)" if canton == "BE" & year == 1971 & liste_norm =="PAB JURA PARTI JURASSIEN DES PAYSANS ARTISANS ET BOURGEOIS"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "BL" & year == 1931 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI BASELLAND"
replace liste_norm = "KATHOLISCHE VOLKSPARTEI UND CHRISTLICHSOZIALE VEREINIGUNG BASELLAND" if canton == "BL" & year == 1931 & liste_norm =="KATHOLISCHE VOLKSPARTEI UND CHRISTLICH SOZIALE VEREINIGUNG"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "BL" & year == 1931 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI BASELLAND"
replace liste_norm = "KATHOLISCHE VOLKSPARTEI UND CHRISTLICH SOZIALE VEREINIGUNG BASELLAND" if canton == "BL" & year == 1935 & liste_norm =="KATHOLISCHE VOLKSPARTEI UND CHRISTLICHSOZIALE VEREINIGUNG BASELLAND"
replace liste_norm = "SCHWEIZ BAUERNHEIMATBEWEGUNG JUNGBAUERN" if canton == "BL" & year == 1935 & liste_norm =="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG JUNGBAUERN"
replace liste_norm = "KATHOLISCHE VOLKSPARTEI UND CHRISTLICHSOZIALE VEREINIGUNG BASELLAND" if canton == "BL" & year == 1939 & liste_norm =="KATHOLISCHE VOLKSPARTEI UND CHRISTLICH SOZIALE VEREINIGUNG BASELLAND"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1939 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI BASEL LANDSCHAFT"
replace liste_norm = "DEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1951 & liste_norm =="DEMOKRATISCHE PARTEI"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1951 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1951 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm = "DEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1955 & liste_norm =="DEMOKRATISCHE PARTEI"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1955 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1955 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI BASEL LAND" if canton == "BL" & year == 1959 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI BASEL LAND" if canton == "BL" & year == 1959 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm = "WIEDERVEREINIGUNGSFREUNDLICHELISTE AKTION KANTON BASEL" if canton == "BL" & year == 1959 & liste_norm =="WIEDERVEREINIGUNGSFREUNDLICHE LISTE AKTION KANTON BASEL"
replace liste_norm = "BAUERN  GEWERBE UND BURGERPARTEI BASELLAND" if canton == "BL" & year == 1967 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI"
replace liste_norm = "CHRISTLICHSOZIALE VOLKSPARTEI BASELLAND" if canton == "BL" & year == 1967 & liste_norm =="CHRISTLICHSOZIALE VOLKSPARTEI"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1967 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1967 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm = "SOZIALDEMOKRATCN UND GEWERKSCHAFTER" if canton == "BL" & year == 1971 & liste_norm =="SOZIALDEMOKRATEN UND GEWERKSCHAFTER"
replace liste_norm = "RADIKALDEMOKRATISCHE PARTEI BASEL" if canton == "BS" & year == 1935 & liste_norm =="RADIKAL DEMOKRATISCHE PARTEI BASEL"
replace liste_norm = "KATHOLISCHE VOLKSPARTEI" if canton == "BS" & year == 1939 & liste_norm =="KATHOLISCHE VOLKSPARTEI BASEL"
replace liste_norm = "KOMMUNISTISCHE PARTEI" if canton == "BS" & year == 1939 & liste_norm =="KOMMUNISTISCHE PARTEI BASEL"
replace liste_norm = "LIBERALE PARTEI" if canton == "BS" & year == 1939 & liste_norm =="LIBERALE PARTEI BASEL"
replace liste_norm = "NATIONALE VOLKSPARTEI" if canton == "BS" & year == 1939 & liste_norm =="NATIONALE VOLKSPARTEI BASEL"
replace liste_norm = "RADIKAL DEMOKRATISCHE PARTEI" if canton == "BS" & year == 1939 & liste_norm =="RADIKAL DEMOKRATISCHE PARTEI BASEL"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "BS" & year == 1939 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI BASEL"
replace liste_norm = "RADIKAL DEMOKRATISCHE PARTEI" if canton == "BS" & year == 1943 & liste_norm =="RADIKAL DEMOKRATISCHE PARTEI BASEL"
replace liste_norm = "BURGER UND GEWERBE PARTEI" if canton == "BS" & year == 1947 & liste_norm =="BURGER UND GEWERBEPARTEI"
replace liste_norm = "RADIKAL DEMOKRATISCHE PARTEI" if canton == "BS" & year == 1951 & liste_norm =="RADIKALDEMOKRATISCHE PARTEI"
replace liste_norm = "RADIKAL DEMOKRATISCHE PARTEI" if canton == "BS" & year == 1955 & liste_norm =="RADIKALDEMOKRATISCHE PARTEI"
replace liste_norm = "LISTE7 KATHOLISCHE VOLKSPARTEI" if canton == "BS" & year == 1959 & liste_norm =="KATHOLISCHE VOLKSPARTEI"
replace liste_norm = "KATHOLISCHE UND CHRISTLICHSOZIALE VOLKSPARTEI" if canton == "BS" & year == 1963 & liste_norm =="KATHOLISCHE UND CHRISTLICHSOZIALE VOLKSPARTEI KANTONALE ORGANISATION DER KONSERVATIV CHRISTLICHSOZIALEN VOLKSPARTEI DER SCHWEIZ"
replace liste_norm = "LIBERAL DEMOKRATISCHE BURGERPARTEI" if canton == "BS" & year == 1963 & liste_norm =="LIBERAL DEMOKRATISCHE BURGERPARTEI BASEL STADT"
replace liste_norm = "RADIKAL DEMOKRATISCHE PARTEI" if canton == "BS" & year == 1963 & liste_norm =="RADIKAL DEMOKRATISCHE PARTEI BASEL STADT KANTONALE ORGANISATION DER FREISINNIG DEMOKRATISCHEN PARTEI DER SCHWEIZ"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "BS" & year == 1963 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI BASEL STADT"
replace liste_norm = "CHRISTLICHDEMOKRATISCHE VOLKSPARTEI BASEL STADT (CVF)" if canton == "BS" & year == 1971 & liste_norm =="CHRISTLICHDEMOKRATISCHE VOLKSPARTEI BASEL STADT (CVP)"
replace liste_norm = "NATIONALE AKTION GEGEN DIE UBERFREMDUNG VONVOLK UND HEIMAT" if canton == "BS" & year == 1971 & liste_norm =="NATIONALE AKTION GEGEN DIE UBERFREMDUNG VON VOLK UND HEIMAT"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI BASEL STADT" if canton == "BS" & year == 1971 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI BS"
replace liste_norm = "LISTE CONSERVATRICE CHRETIENNE SOCIALE" if canton == "FR" & year == 1963 & liste_norm =="LISTE CONSERVATRICE CHRETIENNE SOCIALE  KONSERVATIV CHRISTLICHSOZIALE LISTE"
replace liste_norm = "LISTE DES PAYSANS ARTISANS ET INDEPENDANTS" if canton == "FR" & year == 1963 & liste_norm =="LISTE DES PAYSANS ARTISANS ET INDEPENDANTS  BAUERN GEWERBE UND BURGERPARTEI"
replace liste_norm = "LISTE RADICALE DEMOCRATIQUE" if canton == "FR" & year == 1963 & liste_norm =="LISTE RADICALE DEMOCRATIQUE  FREISINNIG DEMOKRATISCHE LISTE"
replace liste_norm = "LISTE SOCIALISTE" if canton == "FR" & year == 1963 & liste_norm =="LISTE SOCIALISTE  SOZIALISTISCHE LISTE"
replace liste_norm = "RADICALE DEMOCRATIQUE" if canton == "FR" & year == 1967 & liste_norm =="LISTE RADICALE DEMOCRATIQUE  RADIKAL DEMOKRATISCHE LISTE"
replace liste_norm = "SOCIALISTE" if canton == "FR" & year == 1967 & liste_norm =="LISTE SOCIALISTE  SOZIALISTISCHE LISTE"
replace liste_norm = "PARTI CONSERVATEUR CHRETIEN SOCIAL ET PARTI INDEPENDANT CHRETIEN SOCIAL" if canton == "FR" & year == 1967 & liste_norm =="PARTI CONSERVATEUR CHRETIEN SOCIAL ET DU PARTI INDEPENDANT CHRETIEN SOCIAL  KONSERVATIVCHRISTLICHSOZIALE VOLKSPARTEI UND DER UNABHANGIG CHRISTLICHSOZIALEN PARTEI"
replace liste_norm = "PARTI FRIBOURGEOIS DES PAYSANS ARTISANS ET DES INDEPENDANTS" if canton == "FR" & year == 1967 & liste_norm =="PARTI FRIBOURGEOIS DES PAYSANS ARTISANS ET DES INDEPENDANTS  BAUERN  GEWERBE UND BURGERPARTEI DES KANTONS FREIBURG"
replace liste_norm = "PARTI SOCIALISTE" if canton == "GE" & year == 1931 & liste_norm =="PARTI SOCIALISTE GENEVOIS"
replace liste_norm = "PARTI INDEPENDANT ET CHRETIEN SOCIAL GENEVOIS" if canton == "GE" & year == 1939 & liste_norm =="PARTI INDEPENDANT CHRETIEN SOCIAL"
replace liste_norm = "PART RADICAL" if canton == "GE" & year == 1939 & liste_norm =="PARTI RADICAL"
replace liste_norm = "PARTI SOCIALISTE DE GENEVE" if canton == "GE" & year == 1939 & liste_norm =="PARTI SOCIALISTE DE GENÈVE"
replace liste_norm = "PARTI SOCIALISTE GENEVOIS" if canton == "GE" & year == 1955 & liste_norm =="PARTI SOCIALISTE"
replace liste_norm = "PARTI SOCIALISTE GENEVOIS" if canton == "GE" & year == 1959 & liste_norm =="PARTI SOCIALISTE"
replace liste_norm = "KONSERVATIV DEMOKRATISCHEPARTEI" if canton == "GR" & year == 1935 & liste_norm =="KONSERVATIV DEMOKRATISCHE PARTEI"
replace liste_norm = "FREISINNIG DEMOKRATISCH PARTEI" if canton == "GR" & year == 1951 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm = "KONSERVATIVE UND CHRISTLICHSOZIALE VOLKSPARTEI" if canton == "GR" & year == 1963 & liste_norm =="KONSERVATIV CHRISTLICHSOZIALE VOLKSPARTEI"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI (ARBEITERUNION)" if canton == "LU" & year == 1931 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm = "KONSERVATIVE UND CHRISTLICHSOZIALE VOLKSPARTEI DES KANTONS LUZERN" if canton == "LU" & year == 1943 & liste_norm =="KONSERVATIVE UND CHRISTLICHSOZIALE VOLKSPARTEI"
replace liste_norm = "LIBERALE PARTEI DES KANTONS LUZERN" if canton == "LU" & year == 1943 & liste_norm =="LIBERALE PARTEI"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI DES KANTONS LUZERN" if canton == "LU" & year == 1943 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI UND GEWERKSCHAFTSKARTELL" if canton == "LU" & year == 1951 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI UND GEWERKSCHAFTSKARTEL"
replace liste_norm = "LIBERATE PARTEI" if canton == "LU" & year == 1963 & liste_norm =="LIBERALE PARTEI"
replace liste_norm = "LIBERALE PARTEI" if canton == "LU" & year == 1967 & liste_norm =="LIBERATE PARTEI"
replace liste_norm = "PARTI DEMOCRATE POPULAIRE NEUCHATELOIS" if canton == "NE" & year == 1931 & liste_norm =="PARTI DEMOCRATE POPULAIRE NEUCHÂTELOIS"
replace liste_norm = "LIBERALE" if canton == "NE" & year == 1935 & liste_norm =="LISTE LIBERALE"
replace liste_norm = "PROGRESSISTE NATIONALE" if canton == "NE" & year == 1935 & liste_norm =="LISTE PROGRESSISTE NATIONALE"
replace liste_norm = "RADICALE" if canton == "NE" & year == 1935 & liste_norm =="LISTE RADICALE"
replace liste_norm = "SOCIALISTE" if canton == "NE" & year == 1935 & liste_norm =="LISTE SOCIALISTE"
replace liste_norm = "PARTI SOCIALISTE NEUCHATELOIS" if canton == "NE" & year == 1943 & liste_norm =="PARTI SOCIALISTE NEUCHÂTELOIS"
replace liste_norm = "LISTE DU PARTI LIBERAL" if canton == "NE" & year == 1947 & liste_norm =="PARTI LIBERAL"
replace liste_norm = "LISTE DU PARTI OUVRIER ET POPULAIRE" if canton == "NE" & year == 1947 & liste_norm =="PARTI OUVRIER ET POPULAIRE"
replace liste_norm = "LISTE DU PARTI RADICAL" if canton == "NE" & year == 1947 & liste_norm =="PARTI RADICAL"
replace liste_norm = "LISTE DU PARTI SOCIALISTE" if canton == "NE" & year == 1947 & liste_norm =="PARTI SOCIALISTE"
replace liste_norm = "LISTE LIBERALE" if canton == "NE" & year == 1951 & liste_norm =="PARTI LIBERAL"
replace liste_norm = "LISTE DU PARTI OUVRIER ET POPULAIRE" if canton == "NE" & year == 1951 & liste_norm =="PARTI OUVRIER ET POPULAIRE NEUCHATELOIS"
replace liste_norm = "LISTE RADICALE" if canton == "NE" & year == 1951 & liste_norm =="PARTI RADICAL"
replace liste_norm = "LISTE SOCIALISTE" if canton == "NE" & year == 1951 & liste_norm =="PARTI SOCIALISTE"
replace liste_norm = "PARTI OUVRIER ET POPULAIRE" if canton == "NE" & year == 1963 & liste_norm =="PARTI OUVRIER ET POPULAIRE NEUCHATELOIS"
replace liste_norm = "SCHWEIZ BAUERNHEIMATBEWEGUNG JUNGBAUERN" if canton == "SG" & year == 1935 & liste_norm =="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG JUNGBAUERN"
replace liste_norm = "SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG (JUNGBAUERN)" if canton == "SG" & year == 1939 & liste_norm =="SCHWEIZ BAUERNHEIMATBEWEGUNG (JUNGBAUERN)"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG O" if canton == "SG" & year == 1959 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG" if canton == "SG" & year == 1963 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI UND DER JUNGLIBERALEN BEWEGUNG"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI UND GEWERKSCHAFTEN" if canton == "SG" & year == 1963 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI UND DER GEWERKSCHAFTEN"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI UND GEWERKSCHAFTEN" if canton == "SG" & year == 1967 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI UNDJUNGLIBERALE BEWEGUNG LISTE NORD" if canton == "SG" & year == 1971 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG LISTE NORD"
replace liste_norm = "PROGRESSIVE ORGANISATIONEN ST GALLEN (POSG)" if canton == "SG" & year == 1971 & liste_norm =="PROGRESSIVE ORGANISATION ST GALLEN (POSG)"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "SH" & year == 1931 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS SCHAFFHAUSEN"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "SH" & year == 1931 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS SCHAFFHAUSEN"
replace liste_norm = "SCHAFFHAUSER BAUERNPARTEI" if canton == "SH" & year == 1943 & liste_norm =="BAUERNPARTEI"
replace liste_norm = "LISTE DER BAUERNPARTEI" if canton == "SH" & year == 1947 & liste_norm =="BAUERNPARTEI"
replace liste_norm = "LISTE DER FREISINNIG DEMOKRATISCHEN PARTEI" if canton == "SH" & year == 1947 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm = "LISTE DER KATHOLISCHEN VOLKSPARTEI" if canton == "SH" & year == 1947 & liste_norm =="KATHOLISCHE VOLKSPARTEI"
replace liste_norm = "LISTE DER SOZIALISTISCHEN ARBEITERPARTEI" if canton == "SH" & year == 1947 & liste_norm =="SOZIALISTISCHE ARBEITERPARTEI"
replace liste_norm = "SOZIALISTISCHEARBEITERPARTEI" if canton == "SH" & year == 1959 & liste_norm =="SOZIALISTISCHE PARTEI"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "SO" & year == 1931 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS SOLOTHURN"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "SO" & year == 1931 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS SOLOTHURN"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "SO" & year == 1935 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS SOLOTHURN"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "SO" & year == 1935 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS SOLOTHURN"
replace liste_norm = "SOLOTHURNISCHE BAUERN  GEWERBE UND BURGERPARTEI" if canton == "SO" & year == 1947 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI"
replace liste_norm = "SOLOTHURNISCHE VOLKSPARTEI" if canton == "SO" & year == 1947 & liste_norm =="VOLKSPARTEI"
replace liste_norm = "SOLOTHURNISCNE VOLKSPARTEI UND CHRISTLICHSOZIALE" if canton == "SO" & year == 1951 & liste_norm =="SOLOTHURNISCHE VOLKSPARTEI UND CHRISTLICHSOZIALE"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG" if canton == "SO" & year == 1955 & liste_norm =="PREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PAXTEI UND JUNGLIBERALE BEWEGUNG" if canton == "SO" & year == 1959 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG"
replace liste_norm = "ARBEITERPARTEI" if canton == "SZ" & year == 1931 & liste_norm =="ARBEITERPARTEI DES KANTONS SCHWYZ"
replace liste_norm = "KONSERVATIVE VOLKSPARTEI" if canton == "SZ" & year == 1931 & liste_norm =="KONSERVATIVE VOLKSPARTEI DES KANTONS SCHWYZ"
replace liste_norm = "LIBERALE VOLKSPARTEI" if canton == "SZ" & year == 1931 & liste_norm =="LIBERALE VOLKSPARTEI DES KANTONS SCHWYZ"
replace liste_norm = "ARBEITERUNION" if canton == "SZ" & year == 1935 & liste_norm =="ARBEITERUNION DES KANTONS SCHWYZ"
replace liste_norm = "KONSERVATIVE VOLKSPARTEI" if canton == "SZ" & year == 1935 & liste_norm =="KONSERVATIVE VOLKSPARTEI DES KANTONS SCHWYZ"
replace liste_norm = "LIBERALE VOLKSPARTEI" if canton == "SZ" & year == 1935 & liste_norm =="LIBERALE VOLKSPARTEI DES KANTONS SCHWYZ"
replace liste_norm = "BAUERNVEREINIGUNG" if canton == "SZ" & year == 1947 & liste_norm =="BAUERN VEREINIGUNG"
replace liste_norm = "ARBEITERUNION" if canton == "SZ" & year == 1959 & liste_norm =="LISTE DER ARBEITERUNION"
replace liste_norm = "CHRISTLICHSOZIALE PARTEI" if canton == "SZ" & year == 1959 & liste_norm =="LISTE DER CHRISTLICHSOZIALEN PARTEI"
replace liste_norm = "KONSERVATIVE VOLKSPARTEI" if canton == "SZ" & year == 1959 & liste_norm =="LISTE DER KONSERVATIVEN VOLKSPARTEI"
replace liste_norm = "LIBERALE VOLKSPARTEI UND JUNGLIBERALE BEWEGUNG" if canton == "SZ" & year == 1959 & liste_norm =="LISTE DER LIBERALEN VOLKSPARTEI UND JUNGLIBERALEN BEWEGUNG"
replace liste_norm = "ARBEITER AND ANGESTELLTENUNION" if canton == "SZ" & year == 1963 & liste_norm =="ARBEITER UND ANGESTELLTENUNION"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "TG" & year == 1935 & liste_norm =="LISTE DER FREISINNIG DEMOKRATISCHEN PARTEI"
replace liste_norm = "KATHOLISCHE VOLKSPARTEI" if canton == "TG" & year == 1935 & liste_norm =="LISTE DER KATHOLISCHEN VOLKSPARTEI"
replace liste_norm = "BAUERNHEIMATBEWEGUNG JUNGBAUERN" if canton == "TG" & year == 1935 & liste_norm =="LISTE DER SCHWEIZERISCHEN BAUERNHEIMATBEWEGUNG (JUNGBAUERN LISTE)"
replace liste_norm = "JUNG THURGAU" if canton == "TG" & year == 1935 & liste_norm =="LISTE JUNG THURGAU"
replace liste_norm = "SOZIALDEMOKRATISCHE LISTE" if canton == "TG" & year == 1939 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm = "BAUERNPARTEI" if canton == "TG" & year == 1947 & liste_norm =="BAUERN PARTEI"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "TG" & year == 1947 & liste_norm =="SOZIALDEMOKRATISCH GEWERKSCHAFTLICHE PARTEI"
replace liste_norm = "FREISINNIG DEMOKRATISCHE LISTE" if canton == "TG" & year == 1951 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm = "SOZIALDEMOKRATISCHE GEWERKSCHAFTLICHELISTE" if canton == "TG" & year == 1951 & liste_norm =="SOZIALDEMOKRATISCHE UND GEWERKSCHAFTLICHE LISTE"
replace liste_norm = "BAUERNPARTEI" if canton == "TG" & year == 1955 & liste_norm =="BAUERNLISTE"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "TG" & year == 1955 & liste_norm =="FREISINNIG DEMOKRATISCHE LISTE"
replace liste_norm = "CHRISTLICHSOZIALE PARTEI" if canton == "TG" & year == 1955 & liste_norm =="LISTE DER CHRISTLICHSOZIALEN"
replace liste_norm = "SOZIALDEMOKRATISCHE UND GEWERKSCHAFTLICHE PARTEI" if canton == "TG" & year == 1955 & liste_norm =="LISTE DER SOZIALDEMOKRATISCHEN PARTEI UND GEWERKSCHAFTEN"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI UND GEWERKSCHAFTEN" if canton == "TG" & year == 1959 & liste_norm =="LISTE DER SOZIALDEMOKRATISCHEN PARTEI UND GEWERKSCHAFTEN"
replace liste_norm = "CHRISTLICH DEMOKRATISCHE VOLKSPARTEI" if canton == "TG" & year == 1971 & liste_norm =="CHRISTLICH DEMOKRATISCHE VOKSPARTEI"
replace liste_norm = "PARTITO COMUNISTA" if canton == "TI" & year == 1931 & liste_norm =="PARTITO COMMUNISTA"
replace liste_norm = "LISTA DEL GRUPPO AGRARIO POPOLARE TICINESE" if canton == "TI" & year == 1947 & liste_norm =="GRUPPO AGRARIO POPOLARE TICINESE"
replace liste_norm = "LISTA DEL GRUPPO CONSERVATORE DEMOCRATICO" if canton == "TI" & year == 1947 & liste_norm =="GRUPPO CONSERVATORE DEMOCRATICO"
replace liste_norm = "LISTA DEL PARTITO OPERAIO E CONTADINO TICINESE" if canton == "TI" & year == 1947 & liste_norm =="GRUPPO OPERAIO E CONTADINO TICINESE"
replace liste_norm = "LISTA DEL GRUPPO LIBERALE RADICALE" if canton == "TI" & year == 1947 & liste_norm =="PARTITO LIBERALE RADICALE"
replace liste_norm = "LISTA DEL PARTITO SOCIALISTE" if canton == "TI" & year == 1947 & liste_norm =="PARTITO SOCIALISTA"
replace liste_norm = "LISTA DEL GRUPPO CONSERVATORE DEMOCRATICO" if canton == "TI" & year == 1951 & liste_norm =="GRUPPO CONSERVATORE DEMOCRATICO"
replace liste_norm = "LISTA LIBERALE RADICALE" if canton == "TI" & year == 1951 & liste_norm =="PARTITO LIBERALE RADICALE"
replace liste_norm = "LISTA DEL PARTITO SOCIALISTA TICINESE" if canton == "TI" & year == 1951 & liste_norm =="PARTITO SOCIALISTA TICINESE"
replace liste_norm = "CONSERVATORE DEMOCRATICO" if canton == "TI" & year == 1955 & liste_norm =="GRUPPO CONSERVATORE DEMOCRATICO"
replace liste_norm = "LIBERALE RADICALE" if canton == "TI" & year == 1955 & liste_norm =="PARTITO LIBERALE RADICALE"
replace liste_norm = "CONSERVATORE DEMOCRATICO" if canton == "TI" & year == 1959 & liste_norm =="GRUPPO CONSERVATORE DEMOCRATICO"
replace liste_norm = "LIBERALE RADICALE" if canton == "TI" & year == 1959 & liste_norm =="PARTITO LIBERALE RADICALE"
replace liste_norm = "PARTITO AGRARI ARTIGIANI ET PATRIZI" if canton == "TI" & year == 1963 & liste_norm =="PARTITO AGRARI ARTIGIANI E PATRIZI"
replace liste_norm = "LIBERALE RADICALE" if canton == "TI" & year == 1967 & liste_norm =="PARTITO LIBERALE RADICALE"
replace liste_norm = "PARTITE DEL LAVORO" if canton == "TI" & year == 1971 & liste_norm =="PARTITO DEL LAVORO"
replace liste_norm = "PARTITO LIBERATE RADICALE TICINESE" if canton == "TI" & year == 1971 & liste_norm =="PARTITO LIBERALE RADICALE TICINESE"
replace liste_norm = "PARTITE POPOLARE DEMOCRATICO PPD" if canton == "TI" & year == 1971 & liste_norm =="PARTITO POPOLARE DEMOCRATICO (PPD)"
replace liste_norm = "LIBERALE DEMOCRATIQUE" if canton == "VD" & year == 1935 & liste_norm =="LISTE LIBERALE DEMOCRATIQUE"
replace liste_norm = "LISTE DE L'ALLIANCE DES INDEPENDANTS" if canton == "VD" & year == 1943 & liste_norm =="LISTE DE L'ALLIANCE DES INDEPENDANTS VAUDOIS"
replace liste_norm = "LISTE NATIONALE PAYSANNE (LISTE DU PARTI NATIONAL DES PAYSANS ARTISANS ET BOURGEOIS)" if canton == "VD" & year == 1943 & liste_norm =="LISTE NATIONALE PAYSANNE"
replace liste_norm = "LISTE DA PARTI OUVRIER ET POPULAIRE VAUDOIS" if canton == "VD" & year == 1951 & liste_norm =="LISTE DU PARTI OUVRIER ET POPULAIRE VAUDOIS"
replace liste_norm = "PARTI SOCIALISTE VAUDOIS" if canton == "VD" & year == 1955 & liste_norm =="PARTI SOCIALISTE"
replace liste_norm = "PARTI RADICAL DEMOCRATIQUE" if canton == "VD" & year == 1963 & liste_norm =="PARTI RADICAL DEMOCRATIQUE VAUDOIS"
replace liste_norm = "PARTI VAUDOIS DES PAYSANS ARTISANS ET INDEPENDANTS" if canton == "VD" & year == 1963 & liste_norm =="PARTI VAUDOIS DES PAYSANS"
replace liste_norm = "PARTI LIBERALE DEMOCRATIQUE" if canton == "VD" & year == 1967 & liste_norm =="PARTI LIBERAL DEMOCRATIQUE"
replace liste_norm = "PARTI RADICALE DEMOCRATIQUE" if canton == "VD" & year == 1967 & liste_norm =="PARTI RADICAL DEMOCRATIQUE"
replace liste_norm = "PARTI VAUDOIS DES PAYSANS ARTISANS ET INDEPENDANTS (PAI)" if canton == "VD" & year == 1967 & liste_norm =="PARTI VAUDOIS DES PAYSANS ARTISANS ET INDEPENDANTS"
replace liste_norm = "MOUVEMEAT NATIONAL D'ACTION REPUBLICAINE ET SOCIALE" if canton == "VD" & year == 1971 & liste_norm =="MOUVEMENT NATIONAL D'ACTION REPUBLICAINE ET SOCIALE"
replace liste_norm = "LISTE OUVRIERE ET PAYSANNE" if canton == "VS" & year == 1935 & liste_norm =="LISTE OUVRIÈRE ET PAYSANNE"
replace liste_norm = "LISTE OUVRIERE ET PAYSANNE" if canton == "VS" & year == 1939 & liste_norm =="LISTE OUVRIÈRE ET PAYSANNE"
replace liste_norm = "CHRISTLICH SOZIALE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1955 & liste_norm =="LISTE DER CHRISTLICH SOZIALEN VOLKSPARTEI DES OBERWALLIS"
replace liste_norm = "KONSERVATIVE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1955 & liste_norm =="LISTE DER KONSERVATIVEN VOLKSPARTEI DES OBERWALLIS"
replace liste_norm = "MOUVEMENT SOCIAL PAYSAN INDEPENDANT" if canton == "VS" & year == 1955 & liste_norm =="LISTE DU MOUVEMENT SOCIAL PAYSAN INDEPENDANT"
replace liste_norm = "PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND" if canton == "VS" & year == 1955 & liste_norm =="LISTE DU PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND"
replace liste_norm = "PARTI RADICAL DEMOCRATIQUE" if canton == "VS" & year == 1955 & liste_norm =="LISTE DU PARTI RADICAL DEMOCRATIQUE"
replace liste_norm = "CHRISTLICHSOZIALE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1959 & liste_norm =="LISTE DER CHRISTLICH SOZIALEN VOLKSPARTEI OBERWALLIS"
replace liste_norm = "KONSERVATIVE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1959 & liste_norm =="LISTE DER KONSERVATIVEN VOLKSPARTEI OBERWALLIS"
replace liste_norm = "MOUVEMENT SOCIAL DES PAYSANS OUVRIERS ET INDEPENDANTS" if canton == "VS" & year == 1959 & liste_norm =="LISTE DU MOUVEMENT SOCIAL DES PAYSANS OUVRIERS ET INDEPENDANTS"
replace liste_norm = "PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND" if canton == "VS" & year == 1959 & liste_norm =="LISTE DU PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND"
replace liste_norm = "PARTI RADICAL DEMOCRATIQUE" if canton == "VS" & year == 1959 & liste_norm =="LISTE DU PARTI RADICAL DEMOCRATIQUE"
replace liste_norm = "PARTI SOCIALISTE" if canton == "VS" & year == 1959 & liste_norm =="LISTE DU PARTI SOCIALISTE"
replace liste_norm = "CHRISTLICHSOZIALE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1963 & liste_norm =="LISTE DER CHRISTLICH SOZIALEN VOLKSPARTEI OBERWALLIS"
replace liste_norm = "KONSERVATIVE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1963 & liste_norm =="LISTE DER KONSERVATIVEN VOLKSPARTEI OBERWALLIS"
replace liste_norm = "PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND" if canton == "VS" & year == 1963 & liste_norm =="LISTE DU PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND"
replace liste_norm = "PARTI RADICAL DEMOCRATIQUE" if canton == "VS" & year == 1963 & liste_norm =="LISTE DU PARTI RADICAL DEMOCRATIQUE"
replace liste_norm = "PARTI SOCIALISTE" if canton == "VS" & year == 1963 & liste_norm =="LISTE SOCIALISTE  SOZIALDEMOKRATISCHE LISTE"
replace liste_norm = "ALLIANCE DES INDEPENDANTS" if canton == "VS" & year == 1967 & liste_norm =="ALLIANCE DES INDEPENDANTS  LANDESRING DER UNABHANGIGEN"
replace liste_norm = "LISTE SOCIALISTE" if canton == "VS" & year == 1967 & liste_norm =="LISTE SOCIALISTE  SOZIALDEMOKRATISCHE LISTE"
replace liste_norm = "SOCIALISTE POPULAIRE" if canton == "VS" & year == 1967 & liste_norm =="LISTE SOCIALISTE POPULAIRE  SOZIALISTISCHE VOLKSPARTEI"
replace liste_norm = "MOUVEMENT SOCIAL INDEPENDENT" if canton == "VS" & year == 1967 & liste_norm =="MOUVEMENT SOCIAL INDEPENDANT"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "ZG" & year == 1935 & liste_norm =="FREISINNIG DEMOKRATISCHE LISTE"
replace liste_norm = "KONSERVATIVE VOLKS UND ARBEITERPARTEI" if canton == "ZG" & year == 1935 & liste_norm =="KONSERVATIVE VOLKS UND ARBEITERLISTE"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "ZG" & year == 1935 & liste_norm =="SOZIALDEMOKRATISCHE LISTE"
replace liste_norm = "KONSERVATIVE VOLKS ARBEITERPARTEI" if canton == "ZG" & year == 1943 & liste_norm =="KONSERVATIVE VOLKS UND ARBEITERPARTEI"
replace liste_norm = "KONSERVATIV CHRISTLICHSOZIALE PARTEI" if canton == "ZG" & year == 1955 & liste_norm =="KONSERVATIV CHRISTLICHSOZIALE LISTE"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "ZG" & year == 1955 & liste_norm =="SOZIALDEMOKRATISCHE LISTE"
replace liste_norm = "CONSERVATIV CHRISTLICHSOZIALE PARTEI" if canton == "ZG" & year == 1959 & liste_norm =="CONSERVATIV CHRISTLICHSOZIALE LISTE"
replace liste_norm = "FREISINNIG DEMOKRATISCHE PARTEI" if canton == "ZG" & year == 1959 & liste_norm =="FREISINNIG DEMOKRATISCHE LISTE"
replace liste_norm = "SOZIALDEMOKRATISCHE PARTEI" if canton == "ZG" & year == 1959 & liste_norm =="SOZIALDEMOKRATISCHE LISTE"
replace liste_norm = "KONSERVATIV CHRISTLICHSOZIALE PARTEI" if canton == "ZG" & year == 1967 & liste_norm =="CONSERVATIV CHRISTLICHSOZIALE PARTEI"
replace liste_norm = "FREISINNIGE PARTEI" if canton == "ZH" & year == 1931 & liste_norm =="FREISINNIGE LISTE"
replace liste_norm = "CHRISTLICHSOZIALE LISTE" if canton == "ZH" & year == 1935 & liste_norm =="CHRISTLICH SOZIALE LISTE"
replace liste_norm = "LISTE DER SCHWEIZERISCHEN BAUERNHEIMATBEWEGUNG (JUNGBAUERN)" if canton == "ZH" & year == 1935 & liste_norm =="LISTE DER SCHWEIZERISCHEN BAUERNHEIMATBEWEGUNG JUNGBAUERN"
replace liste_norm = "LISTE DER KANTONALEN BAUERNPARTEI BAUERLICH GEWERBLICHBURGERLICHE LISTE" if canton == "ZH" & year == 1939 & liste_norm =="LISTE DER KANTONALEN BAUERNPARTEI BAUERLICH GEWERBLICH BURGERLICHE LISTE"
replace liste_norm = "SOZIALDEMOKRATISCHELISTE" if canton == "ZH" & year == 1939 & liste_norm =="SOZIALDEMOKRATISCHE LISTE"
replace liste_norm = "CHRISTLICHSOZIALE LISTE" if canton == "ZH" & year == 1943 & liste_norm =="LISTE DER CHRISTLICH SOZIALEN PARTEI"
replace liste_norm = "LISTE DER SCHWEIZ BAUERN HEIMATBEWEGUNG" if canton == "ZH" & year == 1943 & liste_norm =="LISTE DER SCHWEIZERISCHEN BAUERNHEIMATBEWEGUNG"
replace liste_norm = "EVANGELISCHEN VOLKSPARTEI" if canton == "ZH" & year == 1951 & liste_norm =="EVANGELISCHE VOLKSPARTEI"
replace liste_norm = "FREISINNIGE LISTE ZURICH STADT" if canton == "ZH" & year == 1951 & liste_norm =="PREISINNIGE LISTE ZURICH STADT"
replace liste_norm = "LISTE DER BAUERN  GEWERBE UND BURGERPARTEI" if canton == "ZH" & year == 1955 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI"
replace liste_norm = "LISTE DER EVANGELISCHEN VOLKSPARTEI" if canton == "ZH" & year == 1955 & liste_norm =="EVANGELISCHE VOLKSPARTEI"
replace liste_norm = "LISTE DER PARTEI DER ARBEIT" if canton == "ZH" & year == 1955 & liste_norm =="PARTEI DER ARBEIT"
replace liste_norm = "FREISINNIGE LISTE ZURICH LAND" if canton == "ZH" & year == 1955 & liste_norm =="PREISINNIGE LISTE ZURICH LAND"
replace liste_norm = "LISTE EVANGELISCHE VOLKSPARTEI" if canton == "ZH" & year == 1959 & liste_norm =="EVANGELISCHE VOLKSPARTEI"
replace liste_norm = "FREISINNIGE LISTE STADT ZURICH" if canton == "ZH" & year == 1963 & liste_norm =="FREISINNIGE LISTE ZURICH STADT"
replace liste_norm = "LISTE DER BAUERN  GEWERBE UND BURGERPARTEI (MITTELSTANDSLISTE)" if canton == "ZH" & year == 1963 & liste_norm =="LISTE DER BAUERN  GEWERBE UND BURGERPARTEI (MITTELSTANDS LISTE)"
replace liste_norm = "LISTE 13 LISTE DER UBERPARTEILICHEN UNION" if canton == "ZH" & year == 1963 & liste_norm =="LISTE DER UBERPARTEILICHEN UNION"
replace liste_norm = "BGB MITTELSTANDSLISTE ZURICH LAND" if canton == "ZH" & year == 1967 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI ZURICH LAND (MITTELSTANDSLISTE)"
replace liste_norm = "LISTE DER BAUERN  GEWERBE UND BURGERPARTEI ZURICH STADT (MITTELSTANDSLISTE)" if canton == "ZH" & year == 1967 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI ZURICH STADT (MITTELSTANDSLISTE)"
replace liste_norm = "LISTE DER SOZIALDEMOKRATEN GEWERKSCHAFTER UND ANGESTELLTEN" if canton == "ZH" & year == 1967 & liste_norm =="SOZIALDEMOKRATEN GEWERKSCHAFTER UND ANGESTELLTE"
replace liste_norm = "ERWA BUND" if canton == "ZH" & year == 1971 & liste_norm =="ERWA BUND (KAMPF FUR RECHT UND UMWELTSCHUTZ)"
replace liste_norm = "LANDSEKTIONEN DER NATIONALEN AKTION GEGEN DIE UBERFREMDUNG VON VOLK UND HEIMAT" if canton == "ZH" & year == 1971 & liste_norm =="LANDSEKTION DER NATIONALEN AKTION GEGEN DIE UBERFREMDUNG VON VOLK UND HEIMAT"
replace liste_norm = "BGB MITTELSTANDSPARTEI" if canton == "ZH" & year == 1971 & liste_norm =="LISTE DER BGB MITTELSTANDSPARTEI"
replace liste_norm = "EVANGELISCHE VOLKSPARTEI (EVP)" if canton == "ZH" & year == 1971 & liste_norm =="LISTE EVANGELISCHE VOLKSPARTEI (EVP)"
replace liste_norm = "JUNGE MITTE" if canton == "ZH" & year == 1971 & liste_norm =="LISTE JUNGE MITTE"
sort canton year liste_norm
save "${path}\16_alliances\Listenverbindungen_1931-2015_Norm.dta", replace
