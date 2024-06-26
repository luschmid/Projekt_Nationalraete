version 17
clear all
set more off
cap log close

*global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"
global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\" /* Office */
global dump "C:\Users\schmutzy\Desktop\Projekt Nationalräte\" /* Office Dumpster */
global ED "C:\Users\schmutzy\Dropbox\Travaux_EmilieDousse\" /* Office ED */

********************************************************************************
*Standardization of Firm names in Bisnode*
use "$dump\02_Processed_data\10_Directors_1934_2003\18_Bisnode\Bisnode_Firmen_corr.dta", clear
drop if rechtsform == "Einzelfirma" 
drop if rechtsform == "Kollektivgesellschaft"
drop if rechtsform == "Anstalt" 
drop if rechtsform == "Duns Support Record" 
drop if rechtsform == "Einfache Gesellschaft"
drop if rechtsform == "Genossenschaft"
drop if rechtsform == "Gesellschaft mit beschränkter Haftung"
drop if rechtsform == "Haupt von Gemeinderschaften"
drop if rechtsform == "Kommanditgesellschaft"
drop if rechtsform == "Kommanditgesellschaft für kollektive Kapitalanlagen"
drop if rechtsform == "Rechtsform unbekannt"
drop if rechtsform == "Stiftung"
drop if rechtsform == "Verein"
drop if rechtsform == "Öffentlich-rechtliche Institution"
drop if missing(rechtsform)

* Simple & obvious character corrections
gen firma_orig = firma
*Same as in sugarcube
	replace firma = ustrlower(firma)
	replace firma = usubinstr(firma, "ä", "ae",.)
	replace firma = usubinstr(firma, "ö", "oe",.)
	replace firma = usubinstr(firma, "ü", "ue",.)
	replace firma = usubinstr(firma, "è", "e",.)
	replace firma = usubinstr(firma, "é", "e",.)
	replace firma = usubinstr(firma, "ê", "e",.)
	replace firma = usubinstr(firma, "ë", "e",.)
	replace firma = usubinstr(firma, "à", "a",.)
	replace firma = usubinstr(firma, "â", "a",.)
	replace firma = usubinstr(firma, "û", "u",.)
	replace firma = usubinstr(firma, "ù", "u",.)
	replace firma = usubinstr(firma, "ô", "o",.)
	replace firma = usubinstr(firma, "ò", "o",.)
	replace firma = usubinstr(firma, "ó", "o",.)	
	replace firma = usubinstr(firma, "î", "i",.)
	replace firma = usubinstr(firma, "ì", "i",.)
	replace firma = usubinstr(firma, "ï", "i",.)
	replace firma = usubinstr(firma, "ç", "c",.)
	replace firma = usubinstr(firma, " und ", " & ",.)
	replace firma = usubinstr(firma, " et ", " & ",.)
	replace firma = usubinstr(firma, "& co.", "& co",.)
	replace firma = usubinstr(firma, "& cie.", "& co",.)
	replace firma = usubinstr(firma, "& cie", "& co",.)
	replace firma = usubinstr(firma, "et cie", "& co",.)
	replace firma = usubinstr(firma, "&cie", "& co",.)
	replace firma = usubinstr(firma, "vormals", "vorm",.)
	replace firma = usubinstr(firma, ",", " ",.)
	replace firma = usubinstr(firma, "/", " ",.)
	replace firma = usubinstr(firma, "+", " ",.)
	replace firma = usubinstr(firma, "%", " ",.)
	replace firma = usubinstr(firma, ":", " ",.)
	replace firma = usubinstr(firma, "'", " ",.)
	replace firma = usubinstr(firma, "succursale de", "",.)
	replace firma = usubinstr(firma, "succursale", "",.)
	replace firma = usubinstr(firma, "branch of", "",.)
	replace firma = usubinstr(firma, "branch", "",.)
	replace firma = usubinstr(firma, "zweigniederlassung", "",.)
	replace firma = usubinstr(firma, "ag fuer", "ag",.)
	replace firma = usubinstr(firma, "ag fur", "ag",.)
	replace firma = usubinstr(firma, "ag.", "ag",.)
	replace firma = usubinstr(firma, "aktiengesellschaft", "ag",.)
	replace firma = usubinstr(firma, "(ag)", "ag",.)
	replace firma = usubinstr(firma, "a.-g.", "ag",.)
	replace firma = usubinstr(firma, "a.g.", "ag",.)
	replace firma = usubinstr(firma, " sa ", " ",.)
	replace firma = usubinstr(firma, "s.a.", "sa",.)
	replace firma = usubinstr(firma, "s. a.", "sa",.)
	replace firma = usubinstr(firma, "s.a", "sa",.)
	replace firma = usubinstr(firma, "sa.", "sa",.)
	replace firma = usubinstr(firma, "schweiz.", "schweizerische",.)
	replace firma = usubinstr(firma, "holding", "",.)
	replace firma = usubinstr(firma, "vorm.", "vorm",.)
	replace firma = usubinstr(firma, "compagnie", "",.)
	replace firma = usubinstr(firma, "cie.", "& co",.)
	replace firma = usubinstr(firma, "co.", "& co",.)
	replace firma = usubinstr(firma, "-", " ",.)
	replace firma = usubinstr(firma, "    ", " ",.)
	replace firma = usubinstr(firma, "   ", " ",.)
	replace firma = usubinstr(firma, "  ", " ",.)
	replace firma = ustrtrim(firma)
	
*specific to Bisnode
	replace firma = usubinstr(firma, "geloescht gemaess liquidation", "",.)
	replace firma = usubinstr(firma, "en liquidati", "",.)
	replace firma = usubinstr(firma, "en liquidation", "",.)
	replace firma = usubinstr(firma, "i. liquidation", "",.)
	replace firma = usubinstr(firma, "en. liquidation", "",.)
	replace firma = usubinstr(firma, "e liquidation", "",.)
	replace firma = usubinstr(firma, "en liquidatation", "",.)
	replace firma = usubinstr(firma, "en liquidaiton", "",.)
	replace firma = usubinstr(firma, "enliquidation", "",.)
	replace firma = usubinstr(firma, "en liquidaton", "",.)
	replace firma = usubinstr(firma, "en liquidateur", "",.)
	replace firma = usubinstr(firma, "inliquidation", "",.)
	replace firma = usubinstr(firma, "in liquidation", "",.)
	replace firma = usubinstr(firma, "im liquidation", "",.)
	replace firma = usubinstr(firma, "in liquidaton", "",.)
	replace firma = usubinstr(firma, "in liquidiation", "",.)
	replace firma = usubinstr(firma, "in liquidat.", "",.)
	replace firma = usubinstr(firma, "in liq.", "",.)
	replace firma = usubinstr(firma, "in liquid.", "",.)
	replace firma = usubinstr(firma, "in liquidaiton", "",.)
	replace firma = usubinstr(firma, "in liquidtion", "",.)
	replace firma = usubinstr(firma, "in liquidazione", "",.)
	replace firma = usubinstr(firma, "en liquidazione", "",.)
	replace firma = usubinstr(firma, "in liquidaz.", "",.)
	replace firma = usubinstr(firma, "in liquidazi", "",.)
	replace firma = usubinstr(firma, "e liquidazione", "",.)
	replace firma = usubinstr(firma, "en faillite", "",.)
	replace firma = usubinstr(firma, "in fallimento", "",.)
	replace firma = usubinstr(firma, "in nachlassliquidation", "",.)
	replace firma = usubinstr(firma, "in nachlass liquidation", "",.)
	replace firma = usubinstr(firma, "in nachlassliquidaton", "",.)
	replace firma = usubinstr(firma, "in konkursliquidation", "",.)
	replace firma = usubinstr(firma, "in konursliquidation", "",.)
	replace firma = usubinstr(firma, "in konkursliquidati on", "",.)
	replace firma = usubinstr(firma, "in konkursliquid.", "",.)
	replace firma = usubinstr(firma, "in konkusliquidation", "",.)
	replace firma = usubinstr(firma, "in konkurssliquidation", "",.)
	replace firma = usubinstr(firma, "in konkursliquidaton", "",.)
	replace firma = usubinstr(firma, "in konkursliquidatin", "",.)
	replace firma = usubinstr(firma, "dissoute sans liquidation", "",.)
	replace firma = usubinstr(firma, "radiee liquidee faillite cloturee", "",.)
	replace firma = usubinstr(firma, "en liquidqtion par suite de faillite", "",.)
	replace firma = usubinstr(firma, "en liquiditation par suite de faillite", "",.)
	replace firma = usubinstr(firma, "dissoute sans liquida tion le 14.12.1995", "",.)
	replace firma = usubinstr(firma, "dissolution sans liquidation dissoute", "",.)
	replace firma = usubinstr(firma, "(liquidiert)", "",.)
	replace firma = usubinstr(firma, "in liquidatio", "",.)
	replace firma = usubinstr(firma, "(liquid.)", "",.)
	replace firma = usubinstr(firma, "liquidiert", "",.)
	replace firma = usubinstr(firma, "liquidateur", "",.)
	replace firma = usubinstr(firma, "konkursliquidation", "",.)
	replace firma = usubinstr(firma, "par suite de faillite", "",.)
	replace firma = usubinstr(firma, "liquidation faillite", "",.)
	replace firma = usubinstr(firma, "en liquidat", "",.)
	replace firma = usubinstr(firma, "in konkurslliquidation", "",.)
	replace firma = usubinstr(firma, "il liquidazione", "",.)
	replace firma = usubinstr(firma, "en liquidtion", "",.)
	replace firma = usubinstr(firma, "liquidee par faillite", "",.)
	replace firma = usubinstr(firma, "nachlassliquidation", "",.)
	replace firma = usubinstr(firma, "en liquida.", "",.)
	replace firma = usubinstr(firma, "in konkutsliquidation", "",.)
	replace firma = usubinstr(firma, "in liquida.", "",.)
	replace firma = usubinstr(firma, "in liquiditation", "",.)
	replace firma = usubinstr(firma, "in liquidazone e in falliment", "",.)
	replace firma = usubinstr(firma, "in liquidastion", "",.)
	replace firma = usubinstr(firma, "en liquid", "",.)
	replace firma = usubinstr(firma, " liquidazione ", "",.)
	replace firma = usubinstr(firma, " liquidation", "",.)
	replace firma = usubinstr(firma, " in liq", "",.)
	replace firma = usubinstr(firma, "radiee faillite cloturee", "",.)
	replace firma = usubinstr(firma, "faillite cloturee", "",.)
	replace firma = usubinstr(firma, "suite a la faillite", "",.)
	replace firma = usubinstr(firma, "par suite defaillite", "",.)
	replace firma = usubinstr(firma, "(faillitee)", "",.)
	replace firma = usubinstr(firma, "ein faillite", "",.)
	replace firma = usubinstr(firma, "faillite radiee", "",.)
	replace firma = usubinstr(firma, "dissoute par suite faillite", "",.)
	replace firma = usubinstr(firma, "dissoute faillite", "",.)
	replace firma = usubinstr(firma, "par suite da faillite", "",.)
	replace firma = usubinstr(firma, "pa suite de faillite", "",.)
	replace firma = usubinstr(firma, "radiee faillite suspendue defaut d actif", "",.)
	replace firma = usubinstr(firma, "radiee suite cloture judiciaire de faillite", "",.)
	replace firma = usubinstr(firma, "par suite faillite", "",.)
	replace firma = usubinstr(firma, "par faillite", "",.)
	replace firma = usubinstr(firma, "en liq", "",.)
	replace firma = usubinstr(firma, "per suite de faillite", "",.)
	replace firma = usubinstr(firma, "radiee faillite suspendue faute d actif", "",.)
	replace firma = usubinstr(firma, " faillite", "",.)
	replace firma = usubinstr(firma, "(radiee)", "",.)
	replace firma = usubinstr(firma, " radiee", "",.)
	replace firma = usubinstr(firma, "(radiee suite fusion)", "",.)
	replace firma = usubinstr(firma, "s. i.", "si",.)
	replace firma = usubinstr(firma, ".", " ",.)
	replace firma = usubinstr(firma, "()", "",.)
	replace firma = usubinstr(firma, "    ", " ",.)
	replace firma = usubinstr(firma, "   ", " ",.)
	replace firma = usubinstr(firma, "  ", " ",.)
	replace firma = ustrtrim(firma)

save "$dump\02_Processed_data\10_Directors_1934_2003\18_Bisnode\Bisnode_Firmen_corr_clean.dta", replace

*Geolocalisation in Bisnode
use "$dump\02_Processed_data\10_Directors_1934_2003\18_Bisnode\Bisnode_Firmen_corr_clean.dta", clear
** Cleaning of special characters, etc. **

foreach str in ortdomiziladresse {
	replace `str' = usubinstr(`str', "ä", "ae",.)
	replace `str' = usubinstr(`str', "ö", "oe",.)
	replace `str' = usubinstr(`str', "ü", "ue",.)
	replace `str' = usubinstr(`str', "è", "e",.)
	replace `str' = usubinstr(`str', "é", "e",.)
	replace `str' = usubinstr(`str', "ê", "e",.)
	replace `str' = usubinstr(`str', "ë", "e",.)
	replace `str' = usubinstr(`str', "œ", "oe",.)
	replace `str' = usubinstr(`str', "à", "a",.)
	replace `str' = usubinstr(`str', "â", "a",.)
	replace `str' = usubinstr(`str', "û", "u",.)
	replace `str' = usubinstr(`str', "ù", "u",.)
	replace `str' = usubinstr(`str', "ô", "o",.)
	replace `str' = usubinstr(`str', "ò", "o",.)
	replace `str' = usubinstr(`str', "ó", "o",.)	
	replace `str' = usubinstr(`str', "î", "i",.)
	replace `str' = usubinstr(`str', "ì", "i",.)
	replace `str' = usubinstr(`str', "ï", "i",.)
	replace `str' = usubinstr(`str', "ç", "c",.)
	replace `str' = usubinstr(`str', "ß", "ss",.)
	replace `str' = usubinstr(`str', "$", "s",.)
	replace `str' = usubinstr(`str', "£", " ",.)	
	replace `str' = usubinstr(`str', "^", " ",.)	
	replace `str' = usubinstr(`str', "'", " ",.)
	replace `str' = usubinstr(`str', "`", " ",.)
	replace `str' = usubinstr(`str', "´", " ",.)
	replace `str' = usubinstr(`str', "-", " ",.)
	replace `str' = usubinstr(`str', "~", " ",.)	
	replace `str' = usubinstr(`str', ":", " ",.)
	replace `str' = usubinstr(`str', ";", " ",.)
	replace `str' = usubinstr(`str', "    ", " ",.)
	replace `str' = usubinstr(`str', "   ", " ",.)
	replace `str' = usubinstr(`str', "  ", " ",.)
	replace `str' = ustrtrim(`str')
	replace `str' = ustrlower(`str')
}

rename ortdomiziladresse gdename
*rename gruendungsjahr year
recast str40 gdename
*replace year = 1931 if year < 1931
*merge m:1 gdename year using "$dump\02_Processed_data\10_Directors_1934_2003\Gde_GdeNr2018_Centroid_1931-2018_tmp.dta"

sort firma plzdomiziladresse
bysort firma plzdomiziladresse:  gen dup = cond(_N==1,0,_n)
order duns handelsregisternummer uid firma dup
drop if dup > 1

save "$dump\02_Processed_data\10_Directors_1934_2003\18_Bisnode\Bisnode_Firmen_corr_gde_clean.dta", replace

*Tests 
use "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\complete_crosssections_1934_2003.dta", clear
gen matchit_id = _n
save "$dump\02_Processed_data\10_Directors_1934_2003\09_Deduplicates\complete_crosssections_1934_2003_fuzzymatch.dta", replace

ssc install matchit
ssc install freqindex

matchit matchit_id cname using "$dump\02_Processed_data\10_Directors_1934_2003\18_Bisnode\Bisnode_Firmen_corr_gde_clean.dta", idu(duns) txtu(firma)
 replace cname = firma if year ==.
 replace kantondomiziladresse = ustrlower(kantondomiziladresse)
 replace ctn = kantondomiziladresse if year ==.
 replace year = gruendungsjahr if missing(year)
 replace year = 1987 if missing(year)
  drop if year > 1988
  replace year = 1987
  

 
bysort gdename: strgroup cname, gen(ID_sug_bis_1) thresh(0.00000000000001) force



