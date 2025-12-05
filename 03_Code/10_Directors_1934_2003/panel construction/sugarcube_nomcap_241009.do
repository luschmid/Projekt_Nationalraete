version 17
clear all
set more off
cap log close

*global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"
global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\" /*office*/
global dump "C:\Users\schmutzy\Desktop\Projekt Nationalräte\" /*office dumpster*/
*global path "C:\Users\yanni\Dropbox\Projekt Nationalräte\" /*home*/
*global dump "C:\Users\yanni\Desktop\Projekt Nationalräte\" /*home dumpster*/

* Step 1

use "$dump\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_complete.dta", clear

rename capital capital_orig
gen capital = capital_orig
replace capital = ustrlower(capital)
replace capital = ustrregexra(capital,"[a-z]","") 
replace capital = usubinstr(capital, "44.1803)", "",.)
replace capital = usubinstr(capital, " 15 6300     55 6300", "",.)
replace capital = usubinstr(capital, " 10  8704", "",.)
replace capital = usubinstr(capital, "41.1422)", "",.)
replace capital = usubinstr(capital, "38.1807)", "",.)
replace capital = usubinstr(capital, "10.1242)", "",.)   
replace capital = usubinstr(capital, "18.1723)", "",.)  
replace capital = usubinstr(capital, "32.1801)", "",.)
replace capital = usubinstr(capital, "26.1814)", "",.)     
replace capital = usubinstr(capital, "26.1806)", "",.)   
replace capital = usubinstr(capital, "21.1220)", "",.)   
replace capital = usubinstr(capital, "38.1093)", "",.)
replace capital = usubinstr(capital, "9.1400)", "",.)            
replace capital = usubinstr(capital, "6.1752)", "",.)
replace capital = usubinstr(capital, "12.1752)", "",.)
replace capital = usubinstr(capital, "27.1234)", "",.)
replace capital = usubinstr(capital, "18.1214)", "",.)
replace capital = usubinstr(capital, "8.1253)", "",.)
replace capital = usubinstr(capital, ")     34  1223", "",.)
replace capital = usubinstr(capital, ")      12", "",.)
replace capital = usubinstr(capital, ")   6234", "",.)
replace capital = usubinstr(capital, ")     1415", "",.) 
replace capital = usubinstr(capital, "88 8618", "",.)       
replace capital = usubinstr(capital, "&", "",.)
replace capital = usubinstr(capital, "ß", "",.)
replace capital = usubinstr(capital, ") °", "",.)
replace capital = usubinstr(capital, "°", "",.)
replace capital = usubinstr(capital, ")2000â", "",.)
replace capital = usubinstr(capital, ")„", "",.)
replace capital = usubinstr(capital, ")~", "",.)
replace capital = usubinstr(capital, ")236340", "",.)
replace capital = usubinstr(capital, ")10.05", "",.)
replace capital = usubinstr(capital, ")„", "",.)
replace capital = usubinstr(capital, ")236340", "",.)
replace capital = usubinstr(capital, ")01", "",.)
replace capital = usubinstr(capital, ")02", "",.)
replace capital = usubinstr(capital, ")03", "",.)
replace capital = usubinstr(capital, ")04", "",.)
replace capital = usubinstr(capital, ")05", "",.)
replace capital = usubinstr(capital, ")06", "",.)
replace capital = usubinstr(capital, ")07", "",.)
replace capital = usubinstr(capital, ")08", "",.)
replace capital = usubinstr(capital, ")09", "",.)
replace capital = usubinstr(capital, ")001", "",.)
replace capital = usubinstr(capital, ")002", "",.)
replace capital = usubinstr(capital, ")003", "",.)
replace capital = usubinstr(capital, ")004", "",.)
replace capital = usubinstr(capital, ")005", "",.)
replace capital = usubinstr(capital, ")006", "",.)
replace capital = usubinstr(capital, ")007", "",.)
replace capital = usubinstr(capital, ")008", "",.)
replace capital = usubinstr(capital, ")009", "",.)
replace capital = usubinstr(capital, ") 1", "",.)
replace capital = usubinstr(capital, ") 2", "",.)
replace capital = usubinstr(capital, ") 3", "",.)
replace capital = usubinstr(capital, ") 4", "",.)
replace capital = usubinstr(capital, ") 5", "",.)
replace capital = usubinstr(capital, ") 6", "",.)
replace capital = usubinstr(capital, ") 7", "",.)
replace capital = usubinstr(capital, ") 8", "",.)
replace capital = usubinstr(capital, ") 9", "",.) 
replace capital = usubinstr(capital, ") ...", "",.)
replace capital = usubinstr(capital, ") ..", "",.)
replace capital = usubinstr(capital, ") .", "",.)
replace capital = usubinstr(capital, ")  ...", "",.)
replace capital = usubinstr(capital, ")  ..", "",.)
replace capital = usubinstr(capital, ")  .", "",.) 
replace capital = usubinstr(capital, ")4", "",.)
replace capital = usubinstr(capital, ")1", "",.)
replace capital = usubinstr(capital, ")2", "",.)
replace capital = usubinstr(capital, ")3", "",.)
replace capital = usubinstr(capital, ")5", "",.)
replace capital = usubinstr(capital, ")6", "",.)
replace capital = usubinstr(capital, ")7", "",.)
replace capital = usubinstr(capital, ")8", "",.)
replace capital = usubinstr(capital, ")9", "",.)
replace capital = usubinstr(capital, ")*", "",.)
replace capital = usubinstr(capital, ".)", "",.)
replace capital = usubinstr(capital, "ç", "",.)
replace capital = usubinstr(capital, "í", "",.)
replace capital = usubinstr(capital, "ì", "",.)
replace capital = usubinstr(capital, "'", "",.)
replace capital = usubinstr(capital, "-", "",.)
replace capital = usubinstr(capital, ",", ".",.)
replace capital = usubinstr(capital, "!", "",.)
replace capital = usubinstr(capital, "`", "",.)
replace capital = usubinstr(capital, "�", "",.)
replace capital = usubinstr(capital, "/", "",.)
replace capital = usubinstr(capital, "?", "",.)
replace capital = usubinstr(capital, ":", "",.)
replace capital = usubinstr(capital, "è", "",.)
replace capital = usubinstr(capital, "é", "",.)
replace capital = usubinstr(capital, "ê", "",.)
replace capital = usubinstr(capital, "ó", "",.)
replace capital = usubinstr(capital, "ö", "",.)
replace capital = usubinstr(capital, "ô", "",.)
replace capital = usubinstr(capital, "ü", "",.)
replace capital = usubinstr(capital, "ú", "",.)
replace capital = usubinstr(capital, "ù", "",.)
replace capital = usubinstr(capital, "►", "",.) 
replace capital = usubinstr(capital, "á", "",.)
replace capital = usubinstr(capital, "â", "",.)
replace capital = usubinstr(capital, "à", "",.)
replace capital = usubinstr(capital, "ä", "",.)
replace capital = usubinstr(capital, ">", "",.)
replace capital = usubinstr(capital, "<", "",.)
replace capital = usubinstr(capital, "+", "",.)
replace capital = usubinstr(capital, "*", "",.)
replace capital = usubinstr(capital, ").  .", "",.)
replace capital = usubinstr(capital, ") .", "",.)
replace capital = usubinstr(capital, ".. ).   . ", "",.)
replace capital = usubinstr(capital, ").   . ", "",.)
replace capital = usubinstr(capital, "23.1860) ", "",.)
replace capital = usubinstr(capital, ").    . ", "",.)
replace capital = usubinstr(capital, ") „", "",.)
replace capital = usubinstr(capital, "„", "",.)
replace capital = usubinstr(capital, ")  .", "",.)
replace capital = usubinstr(capital, ")   ..", "",.)
replace capital = usubinstr(capital, "704.1295", "",.)
replace capital = usubinstr(capital, ").   .. ", "",.)
replace capital = usubinstr(capital, ").     . ", "",.)
replace capital = usubinstr(capital, ") 0 2", "",.)
replace capital = usubinstr(capital, "28.1293) ", "",.)
replace capital = usubinstr(capital, " )   .    . .    . ", "",.)
replace capital = usubinstr(capital, ") ~", "",.)
replace capital = usubinstr(capital, "~", "",.)
replace capital = usubinstr(capital, "). .    . ", "",.) 
replace capital = usubinstr(capital, "). .. .   . ", "",.)
replace capital = usubinstr(capital, "). .  .  ", "",.)
replace capital = usubinstr(capital, ")   . ", "",.)
replace capital = usubinstr(capital, " .   . ", "",.)
replace capital = usubinstr(capital, " .    .     .   . ", "",.)
replace capital = usubinstr(capital, "36.1227)  ", "",.)
replace capital = usubinstr(capital, "97.1618) ", "",.)
replace capital = usubinstr(capital, ").    .. ", "",.)
replace capital = usubinstr(capital, "80.1224) ", "",.)
replace capital = usubinstr(capital, "8.7000) ", "",.)
replace capital = usubinstr(capital, "14.1254) ", "",.)
replace capital = usubinstr(capital, "5.1291) ", "",.)
replace capital = usubinstr(capital, "18.1291) ", "",.)
replace capital = usubinstr(capital, "11.1232) ", "",.)
replace capital = usubinstr(capital, "29.1023) ", "",.)
replace capital = usubinstr(capital, ")      . ", "",.)

replace capital = usubinstr(capital, "", "",.)
replace capital = usubinstr(capital, "", "",.)
replace capital = usubinstr(capital, "", "",.)
replace capital = usubinstr(capital, "", "",.)
replace capital = usubinstr(capital, "", "",.)
replace capital = usubinstr(capital, "", "",.)
replace capital = usubinstr(capital, "", "",.)
replace capital = usubinstr(capital, "", "",.)


replace capital = usubinstr(capital, ")  .", "",.)
replace capital = usubinstr(capital, ").", "",.)
replace capital = usubinstr(capital, ")", "",.)

replace capital = usubinstr(capital, "   ", " ",.)
replace capital = usubinstr(capital, "  ", " ",.)
replace capital = usubinstr(capital, " ", "",.)
replace capital = ustrtrim(capital)
destring capital, replace force 
replace capital = 0.5 if capital_orig == "0.S"
replace capital = 0.5 if capital_orig == "0.5)."
replace capital = 0.05 if capital_orig == "0:05"
replace capital = 0.3 if capital_orig == "0:3"
replace capital = 0.3 if capital_orig == "0.3) 5"
replace capital = 0.1 if capital_orig == "0-1"
replace capital = 0.1 if capital_orig == "0.'"
replace capital = 0.051 if capital_orig == "0.05l"
replace capital = 0.5 if capital_orig == "32.6602).5"
replace capital = 0.1 if capital_orig == "/0,1"
replace capital = 145.0 if capital_orig == "145:0"
replace capital = 0.1 if capital_orig == "0.1) 5"
replace capital = 0.98 if capital_orig == "O,98"
replace capital = 0.06 if capital_orig == "0.C6"
replace capital = . if capital_orig == "Liesta140,1"
replace capital = . if capital_orig == ",4,0"
replace capital = . if capital_orig == "41,"
replace capital = 0.1 if capital_orig == "0.i"
replace capital = 0.05 if capital_orig == "0..05"
replace capital = 0.05 if capital_orig == "0.0.5"
replace capital = 0.05 if capital_orig == "0.05) Pr.-"
replace capital = 0.1 if capital_orig == "0.1) Vr."
replace capital = 0.1 if capital_orig == "0.1."
replace capital = 0.15 if capital_orig == "0.15."
replace capital = 0.1 if capital_orig == "0.25)' t..."
replace capital = 0.05 if capital_orig == "0.05) S-"
replace capital = 0.05 if capital_orig == "0.05) VrAbi"
replace capital = 0.1 if capital_orig == "0.1) Vrk"
replace capital = 0.05 if capital_orig == "005) VrVr Techinco AG, Vaduz, Succursale di Lugano. Lugano"
replace capital = 0.1 if capital_orig == "0.1)Grienbachstrasse 15, 6300 Zug Abicht Ursula, Industriestrasse 55, 6300"
replace capital = 0.1 if capital_orig =="0.1) 31, 8640 PrRapperswil SG"
replace capital = 4.25 if capital_orig =="4.25)e"
replace capital = 0.15 if capital_orig =="0.15)k"
replace capital = 0.15 if capital_orig =="0.15) Colombière Del131, 2900 Porrentruyk"
replace capital = 0.1 if capital_orig =="0.1)k"
replace capital = 0.5 if capital_orig =="0.5)in Willigen, Amt Oberhasli,"

save "$dump\02_Processed_data\10_Directors_1934_2003\CID_nomcap_temp.dta", replace

* Step 2 : add corrections

*Put back into cross-sections
foreach y in 1934 1943 1960 1962 1963 1964 1965 1966 1969 1972 1975 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 {
merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\21_Capital\comp_`y'_missing_capital_mun_corrected.dta", keep (1 3) nogen
replace capital = capital_corrected if missing(capital)
drop capital_corrected
}

* Step 3 : keeping only relevant vars

keep CID year capital_orig capital

drop if CID == .

save "$dump\02_Processed_data\10_Directors_1934_2003\CID_nomcap.dta", replace