******************************************************
* (A) Set directories
******************************************************

clear
set more off
version 17

global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "C:\Dropbox\Projekt Nationalräte"
*global path "C:\Current\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"
*global path "C:\Users\Lukas\Dropbox\Projekt Nationalräte"

******************************************************
* (B) READ-IN ELECTION RESULTS DATA (1975 until 2015)
******************************************************

* (i) Read in individual election results from 1975 to 2015
* Note: This dataset includes information on the types of votes a candidate re-
* ceived in a certain municipality. Since we are not using this information, we
* collapse the dataset by candidate.

insheet using ///
	"$path\01_Raw_data\02_Elections_1971_2015\NRW_KANDIDATENSTIMMEN.csv", ///
	delim(";") clear

gen year=substr(nationalratswahl,4,4)
destring year, replace
collapse (sum) kandidst_unver_wz kandidst_ver_wz kandidst_tot (mean) year, ///
	by( nationalratswahl listennummer kandidatennummer canton_no)

* (ii) Read in information on a candidate's name, party, list, residence munici-
* pality, gender, and status

preserve
import exc "$path\01_Raw_data\02_Elections_1971_2015\NRW_KANDIDATEN.xlsx", ///
	first clear
* Note: We saved the orginal csv file "NRW_KANDIDATEN including 2015.csv" as 
* "NRW_KANDIDATEN.xlsx" b/c the import excel command handles Umlaute correctly

rename _all, lower
* Note: Transforms all variables to lower case variables.
gen year=substr(nationalratswahl,4,4)
destring year, replace

do "$path\03_Code\08_Municipalities\Add_missing_municipalities.do"
* Note: Add missing municipality information in BfS data; based on Florence's
* coding, January 2017

save "$path\02_Processed_data\02_Elections_1971_2015\NRW_KANDIDATEN.dta", ///
	replace
restore

* (iii)	Merge vote data with candidate information and check not merged cases

gen kantonsnummer=canton_no
merge 1:1 year kantonsnummer kandidatennummer listennummer using ///
	"$path\02_Processed_data\02_Elections_1971_2015\NRW_KANDIDATEN.dta", ///
	gen(merge_candinfos)
* merge_candinfos==1: 53 Stimmen für Vereinzelte in Kantonen mit einem NR-Sitz
* (Majorzkantone) 
* merge_candinfos==2: 6 stille Wahlen

replace canton_no = kantonsnummer if canton_no == .
merge m:1 canton_no using ///
	"$path\01_Raw_data\07_Cantons\cantonal_translator.dta", nogen

* (iv) Merge left-right coding

preserve
import exc "$path\01_Raw_data\06_Parties\LeftRightCoding.xlsx", first ///
	sheet("bfs") clear
rename partei_nr parteinummer
save "$path\02_Processed_data\06_Parties\LeftRightCoding.dta", replace
restore

merge m:1 parteinummer using ///
	"$path\02_Processed_data\06_Parties\LeftRightCoding.dta", ///
	gen(merge_parteien)

* (v) Corrections

replace vorname="Willy" if nachname=="Walker" & vorname=="Willi" & year==1975 ///
	& canton=="ZH"
replace gemeindenummer_bfs=4937 if nachname=="Möckli" & year==1975 ///
	& vorname=="Gustav" & gemeindenummer_bfs==4631 & geburtsjahr==1926
replace gemeindenummer_bfs=195 if nachname=="Reich" & year==1975 ///
	& vorname=="Richard" & gemeindenummer_bfs==253 & geburtsjahr==1927
replace gemeindenummer_bfs=195 if nachname=="Reich" & year==1979 ///
	& vorname=="Richard" & gemeindenummer_bfs==253 & geburtsjahr==1927
replace gemeindenummer_bfs=195 if nachname=="Reich" & year==1983 ///
	& vorname=="Richard" & gemeindenummer_bfs==. & geburtsjahr==1927
replace gemeindenummer_bfs=195 if nachname=="Reich" & year==1987 ///
	& vorname=="Richard" & gemeindenummer_bfs==. & geburtsjahr==1927

replace gemeindename="Maur" if nachname=="Reich" & year==1975 ///
	& vorname=="Richard" & gemeindename=="Zürich" & geburtsjahr==1927
replace gemeindename="Maur" if nachname=="Reich" & year==1979 ///
	& vorname=="Richard" & gemeindename=="Zürich" & geburtsjahr==1927

* (vi) Rename variables in order to merge it to pre-1975 data

gen elected=. 
replace elected=1 if gewaehlt=="G"
replace elected=0 if gewaehlt=="N"

keep nachname vorname kandidatennummer kandidst_tot year kantonsnummer canton ///
	parteiname listenname listennummer gemeindenummer_bfs gemeindename ///
	geschlecht geburtsjahr elected status LinksRechts

rename kandidatennummer candidateno
rename kandidst_tot votes
rename kantonsnummer cantonno
rename listenname list
rename listennummer listno
rename parteiname partyname
rename gemeindenummer_bfs municipalityno
rename gemeindename municipality
rename LinksRechts leftright
rename geschlecht sex
rename geburtsjahr birthyear 
rename vorname firstname
rename nachname name

save "$path\02_Processed_data\02_Elections_1971_2015\NRW_KANDIDATEN 1975-2015.dta", ///
	replace
erase "$path\02_Processed_data\06_Parties\LeftRightCoding.dta"
erase "$path\02_Processed_data\02_Elections_1971_2015\NRW_KANDIDATEN.dta"


******************************************************
* (C) READ-IN ELECTION RESULTS DATA (1931 until 1975)
******************************************************

* (i) Import all excel files

import exc ///
	"$path\01_Raw_data\01_Elections_1931_1975\1931_Bundesblatt_backin.xlsx", ///
	first clear
gen year=1931
replace canton=canton[_n-1] if canton==""
save "$path\02_Processed_data\01_Elections_1931_1975\temp.dta", replace

forvalues y=1935(4)1975 {
import exc ///
	"$path\01_Raw_data\01_Elections_1931_1975\\`y'_Bundesblatt_backin.xlsx", ///
	first clear
display "`y'"
gen year=`y'
replace canton=canton[_n-1] if canton==""
append using "$path\02_Processed_data\01_Elections_1931_1975\temp.dta"
save "$path\02_Processed_data\01_Elections_1931_1975\temp.dta", replace
}

tab year
gen order_initial=_n

* (ii) Correct names of canton

preserve
import exc ///
	"$path\01_Raw_data\01_Elections_1931_1975\Kantonsnamen_harmonisiert.xlsx", ///
	first clear
* Manual harmonization of canton names
save "$path\02_Processed_data\01_Elections_1931_1975\Kantonsnamen_harmonisiert.dta", ///
	replace 
restore

merge m:1 canton using ///
	"$path\02_Processed_data\01_Elections_1931_1975\Kantonsnamen_harmonisiert.dta", ///
	gen(merge_kantonsname)

drop canton 
gen canton=kantonsname

merge m:1 canton using "$path\01_Raw_data\07_Cantons\cantonal_translator.dta", ///
	gen(merge_cantonaltranslator)
* merge_cantonaltranslator==1: one obs for canton of JU
drop if canton=="JU"

* (iii) Add municipality numbers and correct municipality names

/* Outfile for Simon Luechinger for correction of municipality names
preserve
egen city_num=group(city)
collapse (first) city, by(city_num canton)	
export exc "$path\07_Archive\08_Municipalities\Municipalities All.xlsx", ///
	first(variables) replace 
restore
*/ 
gen gemeindename=city

do "$path\03_Code\08_Municipalities\Correct Municipalities.do"
* Correct municipality names; written by Simon.
rename city gemeindename_original

preserve
import exc "$path\01_Raw_data\08_Municipalities\Gemeindebestand_1969.xls", ///
	first sheet("exportList")  cellra(A12:G3090) clear
rename Gemeindename gemeindename
rename BFSGdenummer gemeindenummer_bfs
save "$path\02_Processed_data\08_Municipalities\Gemeindebestand_1969.dta", ///
	replace
restore

merge m:1 gemeindename using ///
	"$path\02_Processed_data\08_Municipalities\Gemeindebestand_1969.dta", ///
	gen(merge_gemeindebestand)
* 181 merge=2 all in 1975
drop if merge_gemeindebestand==2

* (iv) Correct party lists

/* Outfile for Florence
preserve
egen party_num=group(party)
collapse (first) party, by(party_num canton)
export exc "$path\07_Archive\06_Parties\Parteinamen Alle.xlsx", ///
	first(variables) replace
restore		
*/

egen party_num=group(party)

* Read in party file harmonized by Florence and Lukas
preserve
import exc "$path\01_Raw_data\06_Parties\Parteinamen All Florence.xlsx", first ///
	clear 
duplicates drop canton Partei, force
* Drop 12 duplicate observations 
rename Partei party
save "$path\02_Processed_data\06_Parties\Parteinamen All Florence.dta", replace
import exc "$path\01_Raw_data\06_Parties\LeftRightCoding.xlsx", first ///
	sheet("all parties recoded") clear
rename Parteiname party
save "$path\02_Processed_data\06_Parties\LeftRightCoding.dta", replace
restore

merge m:1 canton party using ///
	"$path\02_Processed_data\06_Parties\Parteinamen All Florence.dta", ///
	gen(merge_party_all)

* merge: 1,761 not merged from using data
*           32 from AI, NW, OW, UR
*        1,729 from 1975

split Code, p(",")
destring Code1 Code2 Code3 Code4, replace
rename Code Code_original
reshape long Code, i(order_initial) j(Code_Nr)
drop Code_Nr
rename party party_all
merge m:1 Code using "$path\02_Processed_data\06_Parties\LeftRightCoding.dta", ///
	gen(merge_qualcheck)
drop if merge_qualcheck==2
* 60 merge=2, mostly new parties (Grüne, ...)
collapse (first) no-Code_original (firstnm) party (mean) LinksRechts, ///
	by(order_initial)

* (v) Change names to 1975-2015 format

rename firstname vorname
rename name nachname
rename birth geburtsjahr
rename votes kandidst_tot

gen gewaehlt=.
replace gewaehlt=1 if elected=="true" | elected=="1"
replace gewaehlt=0 if elected=="false"| elected=="0"
drop elected
rename gewaehlt elected

* (vi) Corrections and additional information 

replace elected=0 if year==1967 & canton=="VS" & ///
	party_all=="Liste 1. Liste socialiste"
* 4 candidates misclassified as elected
replace geburtsjahr=1906 if year==1947 & nachname=="Michlig" & ///
	vorname=="Meinard" & kandidst_tot==6468
* This information is according to Bundesblatt 1951

* (vii) Quality checks

* (a) Check number cantons per year

preserve
gen ones=1
collapse (mean) ones , by(year canton)
tab year
restore
* Result: 25 cantons in all years

tab canton if LinksRechts==. & year!=1975
* Result: LinksRechts is only missing for AI, NW, OW, UR

tab year

* (b) Check number candidates per year

preserve
drop if canton=="UR" | canton=="AI" | canton=="NW" |canton=="OW"
drop if canton=="GL" & year>=1971
bysort year: gen nobs=_N
bysort year: gen indi=_n
list nobs year if indi==1
restore
* Result: No. candidates equal to BfS data except for 1967 in which 1 candidate
* is missing in our data; see 07_Archive\01_Elections_1931_1975\Quality_checks
* anzahl kandidierende bfs.xlsx

list nachname vorname geburtsjahr canton if ///
	(canton=="SZ" | canton=="AR" |canton=="GL") & year==1967

bysort year: egen elected_sum=sum(elected)
bysort year: sum elected_sum

preserve
gen ones=1
keep if elected==1
collapse (sum) ones, by(year canton)
table canton if year == 1967, statistic(mean ones)
restore

*br nachname vorname geburtsjahr canton kandidst_tot if year==1967 & ///
*	canton=="VS" & elected==1
*br nachname vorname geburtsjahr canton kandidst_tot elected party_all if ///
*	year==1967 & canton=="VS" 
* Result: 1967 are 4 candidates more than on BfS list; explanation: in one party
* (list 1) in VS all candidates were classified as elected (correction above)

preserve
gen ones=1
collapse (sum) ones (mean) canton_no, by(year canton party_all elected)
sort year canton_no party_all elected
*br year canton_no party_all elected if year==1967
gen party_all_num = regexs(2) ///
	if regexm(party_all, "^([^0-9]*)([0-9]+)([^0-9]*)$")
destring party_all_num, replace
sort year canton_no party_all_num elected 
*br year canton_no party_all elected ones if year==1967 
restore

* Result: All differences to the BfS list (Anzahl Kandidierende bei den Natio-
* nalratswahlen 1928 – 2015) can be explained by "stille Wahlen" and Majorzkan-
* tone except 1967. In 1967 we have one candidate less than on the BfS list. We 
* checked the number of elected and not elected candidates on all lists and com-
* pared it to the original Bundesblatt file. We found no difference. 

*(first) |
*   year |      Freq.     Percent        Cum.
*------------+-----------------------------------
*   1931 |        777        5.59        5.59
*   1935 |        973        7.00       12.59
*   1939 |        783        5.63       18.22
*   1943 |      1,021        7.35       25.57
*   1947 |        998        7.18       32.75
*   1951 |      1,070        7.70       40.45
*   1955 |      1,089        7.83       48.28
*   1959 |      1,077        7.75       56.03
*   1963 |      1,200        8.63       64.66
*   1967 |      1,262        9.08       73.74
*   1971 |      1,696       12.20       85.94
*   1975 |      1,954       14.06      100.00
*------------+-----------------------------------
*  Total |     13,900      100.00

* (c) Check number of elected candidates per canton and year

preserve
gen ones=1
keep if elected==1
collapse (sum) ones (mean) canton_no, by(year canton )
sort canton year
merge 1:1 year canton using ///
	"$path\01_Raw_data\01_Elections_1931_1975\seats.dta", gen(merge_seats) 
gen seats_diff=seats-ones
tab seats_diff
list if seats_diff!=0
restore

sort party_all
* br if canton=="FR" & elected==1 & year==1955
* Result: FR 1955: Pierre Pasquier misclassified as "elected" (corrected in ori-
* ginal file)

* (d) Check changes of lists

preserve
sort order_initial
gen party_all_num = regexs(2) ///
	if regexm(party_all, "^([^0-9]*)([0-9]+)([^0-9]*)$")
destring party_all_num, replace
tsset order_initial
gen party_all_num_diff=d1.party_all_num
gen canton_no_diff=d1.canton_no
tab party_all_num_diff
gen indi_partydiff=0 
replace indi_partydiff=1 ///
	if (party_all_num_diff<0 & canton_no_diff==0) | ///
	(party_all_num_diff>1 & canton_no_diff==0)
gen indi_partydiff_d=d1.indi_partydiff
list nachname vorname geburtsjahr canton year party_all if indi_partydiff==1
restore
* Result: All checked, one party name changed in original lists and in "Partei-
* namen all Florence"

* (e) Check decreasing vote numbers within a list

preserve
sort year canton party_all kandidst_tot order_initial
bysort year canton party_all: gen votes_order=_n
tsset order_initial
gen votes_order_diff=d1.votes_order
gen party_all_num = regexs(2) ///
	if regexm(party_all, "^([^0-9]*)([0-9]+)([^0-9]*)$")
destring party_all_num, replace
gen party_all_num_diff=d1.party_all_num
gen canton_no_diff=d1.canton_no
tab votes_order_diff
gen indi_votesdiff=0 
replace indi_votesdiff=1 ///
	if votes_order_diff!=-1 & canton_no_diff==0 & party_all_num_diff==0
*br nachname vorname geburtsjahr canton year party_all votes_order kandidst_tot ///
*	party_all_num_diff party_all_num canton_no_diff votes_order_diff ///
*	if indi_votesdiff==1
* Result: 339 cases; all correct. Most cases refer to party lists on which 2 or
* more candidates have the same number of votes, some to "stille Wahlen" (1939)
* and others are genuine errors in ordering in the Bundesblatt.
restore

* (f) Check missings in ocd number

*br nachname vorname geburtsjahr canton year party_all kandidst_tot elected ///
*	if ocd==.
* Result: All checked with Bundesblatt and data are correct.

* (g) Check stille Wahlen

*br nachname vorname geburtsjahr canton year party_all kandidst_tot elected ///
*	if kandidst_tot==. & elected==0
* Result: No one was not elected when "stille Wahlen" took place

* (h) Check whether not elected candidates have more votes than elected candi-
* dates on the same list

bysort year canton party_all: g listsize = _N
bysort year canton party_all: egen elected_votes_rank=min(kandidst_tot) if elected==1
bysort year canton party_all: egen not_elected_votes_rank=max(kandidst_tot) if elected==0
*bysort year canton party_all: egen kandidst_tot_max=max(help2)
*drop help1 help2
*sort canton year party_all elected kandidst_tot
*br nachname vorname geburtsjahr canton year party_all elected kandidst_tot ///
*	kandidst_tot_min kandidst_tot_max if kandidst_tot_min<=kandidst_tot_max ///
*	& listsize > 1
* Result: 74 cases; all checked and plausible; "stille Wahlen" or party lists
* with all candidates elected.

* (i) Check whether party changes if canton changes

sort year canton party_all kandidst_tot 
bysort year canton party_all: gen votes_order=_n

tsset order_initial
gen votes_order_diff=d1.votes_order
gen elected_diff=d1.elected

gen party_all_num = regexs(2) ///
	if regexm(party_all, "^([^0-9]*)([0-9]+)([^0-9]*)$")
destring party_all_num, replace
gen party_all_num_diff=d1.party_all_num
gen canton_no_diff=d1.canton_no

*br if party_all_num_diff==0 & canton_no_diff!=0
* Result: Party always changes when canton changes.

* (j) Check whether there exist cases with no change in party list but change 
* from "not elected" to "elected"
	
*br nachname vorname geburtsjahr canton year party_all elected elected_diff ///
*	canton_no_diff party_all_num_diff if party_all_num_diff==0 & elected_diff==1
* Result: No such cases.

* (k) Check whether first list in a canton is list number 1 (and not 2)

tab party_all_num if canton_no_diff!=0, missing 
* Result: All first lists in a canton are list 1

* (vii) Merge 1971 and 1975 information on gender, municipality and municipality
* id

* 1971
* Note that 1971 Kovariate is information from BfS which is merged to Bundes-
* blatt data

preserve
import exc using "$path\01_Raw_data\01_Elections_1931_1975\Kovariate_1971.xlsx", ///
	first clear
foreach x of varlist _all {
rename `x' `x'_1971
} 
rename Kantonsname_1971 canton
rename Vorname_1971 vorname
rename Nachname_1971 nachname 
gen year=1971
replace vorname="Elisabeth" if nachname=="Rüegsegger-Suter" & vorname=="Elisbeth"
replace nachname="Memmishofer" if nachname=="Hemmishofer" & vorname=="Jean"
replace vorname="François" if nachname=="Rossé" & vorname=="Francois"
replace nachname="Kaiser" if nachname=="Kaiser-Kraemer" & vorname=="Gisela"
replace nachname="Tuscher" if nachname=="Tüscher" & vorname=="Hans"
replace nachname="Jenny" if nachname=="Jenni" & vorname=="Olga"
replace vorname="Jürg" if nachname=="Niklaus" & vorname=="JÜrg" 
replace nachname="Hernandez-Kartaschoff" if nachname=="Hernandez-Kartaschof" ///
	& vorname=="Rosemarie"
replace vorname="Ingrid" if nachname=="Risch" & vorname=="Ingrid Henriette"
replace nachname="Bäder" if nachname=="Baeder" & vorname=="Peter"
replace nachname="Güntensperger-Gsell" if nachname=="Guentensperger-Gsell" ///
	& vorname=="Sibylla"
replace vorname="André J.-R." if nachname=="Feignoux" & vorname=="André J.-R"
replace nachname="Dietrich-Schellenberg" if nachname=="Dietrich-Schellenber" ///
	& vorname=="Erica"
replace vorname="Hansjörg W." if nachname=="Furrer" & vorname=="Hansjörg W"
replace nachname="De Lorenzo" if nachname=="de Lorenzo" & vorname=="Kurt"
replace vorname="Martha" if nachname=="Müller-Ledergerber" & vorname=="Mahrta"
replace nachname="à Porta" if nachname=="à.Porta" & vorname=="Reto"
replace vorname="Rudolf-Christoph" if nachname=="Schellenberg" ///
	& vorname=="Rudolf Christoph"
replace vorname="Willy" if nachname=="Walker" & vorname=="Willi"

* Change municipality ids
replace Gemeindenummer_BFS_1971=4260 if nachname=="Schmid-Bruggisser" ///
	& vorname=="Elisabeth" & Gemeindenummer_BFS_1971==4036 
replace Gemeindenummer_BFS_1971=352 if nachname=="Vetsch" & vorname=="Hans" ///
	& Gemeindenummer_BFS_1971==351 
replace Gemeindenummer_BFS_1971=308 if nachname=="Weber" & vorname=="Hans" ///
	& Gemeindenummer_BFS_1971==588 
replace Gemeindenummer_BFS_1971=511 if nachname=="Gerber" & vorname=="Isaac" ///
	& Gemeindenummer_BFS_1971==6800 
replace Gemeindenummer_BFS_1971=445 if nachname=="Beuret" ///
	& vorname=="Jean-Pierre" & Gemeindenummer_BFS_1971==435 
replace Gemeindenummer_BFS_1971=958 if nachname=="Blumenstein" & vorname=="Jürg" ///
	& Gemeindenummer_BFS_1971==955 
replace Gemeindenummer_BFS_1971=352 if nachname=="Blatter" & vorname=="Rolf" ///
	& Gemeindenummer_BFS_1971==351 
replace Gemeindenummer_BFS_1971=2701 if nachname=="Weber" & vorname=="Rudolf" ///
	& Gemeindenummer_BFS_1971==1058 & canton=="BE"
replace Gemeindenummer_BFS_1971=3236 if nachname=="Gautschi" & vorname=="Ernst" ///
	& Gemeindenummer_BFS_1971==4746 
replace Gemeindenummer_BFS_1971=3236 if nachname=="Rosenberger" ///
	& vorname=="Hans-Peter" & Gemeindenummer_BFS_1971==4746 
replace Gemeindenummer_BFS_1971=1346 if nachname=="Diethelm-Dobler" ///
	& vorname=="Josef" & Gemeindenummer_BFS_1971==1349
replace Gemeindenummer_BFS_1971=4937 if nachname=="Saner" & vorname=="Fritz" ///
	& Gemeindenummer_BFS_1971==4631
replace Gemeindenummer_BFS_1971=4534 if nachname=="Möckli" & vorname=="Gustav" ///
	& Gemeindenummer_BFS_1971==4781 
replace Gemeindenummer_BFS_1971=4762 if nachname=="Fritschi" & vorname=="Rudolf" ///
	& Gemeindenummer_BFS_1971==4724 
replace Gemeindenummer_BFS_1971=4671 if nachname=="Hangartner" ///
	& vorname=="Ulrich" & Gemeindenummer_BFS_1971==4641 
replace Gemeindenummer_BFS_1971=5172 if nachname=="Monico" & vorname=="Nice" ///
	& Gemeindenummer_BFS_1971==5192 
replace Gemeindenummer_BFS_1971=1201 if nachname=="Weber" & vorname=="Alfred" ///
	& Gemeindenummer_BFS_1971==6460 
replace Gemeindenummer_BFS_1971=5502 if nachname=="Huggli" & vorname=="Jean" ///
	& Gemeindenummer_BFS_1971==5487 
replace Gemeindenummer_BFS_1971=6236 if nachname=="Rey" & vorname=="Alfred" ///
	& Gemeindenummer_BFS_1971==6248 
replace Gemeindenummer_BFS_1971=6217 if nachname=="Rouiller" ///
	& vorname=="Claude" & Gemeindenummer_BFS_1971==6616 
replace Gemeindenummer_BFS_1971=118 if nachname=="Schenk" ///
	& vorname=="Fritz" & Gemeindenummer_BFS_1971==142 
replace Gemeindenummer_BFS_1971=118 if nachname=="Zedi" & vorname=="Gotthard" ///
	& Gemeindenummer_BFS_1971==142 
replace Gemeindenummer_BFS_1971=118 if nachname=="Rüegg" & vorname=="Hans" ///
	& Gemeindenummer_BFS_1971==142 
replace Gemeindenummer_BFS_1971=118 if nachname=="Hungerbühler" ///
	& vorname=="Hugo" & Gemeindenummer_BFS_1971==142 
replace Gemeindenummer_BFS_1971=118 if nachname=="Lienhard" & vorname=="Konrad" ///
	& Gemeindenummer_BFS_1971==142 
replace Gemeindenummer_BFS_1971=118 if nachname=="Kunz" & vorname=="Walter" ///
	& Gemeindenummer_BFS_1971==142 
replace Gemeindenummer_BFS_1971=118 if nachname=="Kyburz" & vorname=="Walter" ///
	& Gemeindenummer_BFS_1971==142 
replace Gemeindenummer_BFS_1971=63 if nachname=="Furrer" ///
	& vorname=="Hansjörg W." & Gemeindenummer_BFS_1971==62 
replace Gemeindenummer_BFS_1971=192 if nachname=="Basler" & vorname=="Konrad" ///
	& Gemeindenummer_BFS_1971==55 
replace Gemeindenummer_BFS_1971=192 if nachname=="Meier" & vorname=="Kurt" ///
	& Gemeindenummer_BFS_1971==55 
replace Gemeindenummer_BFS_1971=195 if nachname=="Reich" & vorname=="Richard" ///
	& Gemeindenummer_BFS_1971==154 

replace Gemeinde_1971="Bolligen" if nachname=="Vetsch" & vorname=="Hans" ///
	& Gemeinde_1971=="Bern" 
replace Gemeinde_1971="Le Bémont (BE)" if nachname=="Gerber" & vorname=="Isaac" ///
	& Gemeinde_1971=="Aux Rouges-Terre" 
replace Gemeinde_1971="Sonvilier" if nachname=="Beuret" & vorname=="Jean-Pierre" ///
	& Gemeinde_1971=="La Chaux d'Abel"
replace Gemeinde_1971="Bolligen" if nachname=="Blatter" & vorname=="Rolf" ///
	& Gemeinde_1971=="Bern" 
replace Gemeinde_1971="Basel" if nachname=="Weber" & vorname=="Rudolf" ///
	& Gemeinde_1971=="St. Niklausen" & canton=="BE"
replace Gemeinde_1971="Schübelbach" if nachname=="Diethelm-Dobler" ///
	& vorname=="Josef" & Gemeinde_1971=="Wangen (SZ)" 
replace Gemeinde_1971="Unterschlatt" if nachname=="Möckli" & vorname=="Gustav" ///
	& Gemeinde_1971=="Wängi" 
replace Gemeinde_1971="Kreuzlingen" if nachname=="Hangartner" ///
	& vorname=="Ulrich" & Gemeinde_1971=="Altnau" 
replace Gemeinde_1971="Castagnola" if nachname=="Monico" & vorname=="Nice" ///
	& Gemeinde_1971=="Lugano"
replace Gemeinde_1971="Altdorf (UR)" if nachname=="Weber" & vorname=="Alfred" ///
	& Gemeinde_1971=="Altdorf" 
replace Gemeinde_1971="Granges (VS)" if nachname=="Rey" & vorname=="Alfred" ///
	& Gemeinde_1971=="Sierre" 
replace Gemeinde_1971="Lufingen" if nachname=="Furrer" & vorname=="Hansjörg W." ///
	& Gemeinde_1971=="Kloten" 
replace Gemeinde_1971="Maur" if nachname=="Reich" & vorname=="Richard" ///
	& Gemeinde_1971=="Küsnacht (ZH)" 

save "$path\02_Processed_data\01_Elections_1931_1975\Kovariate_1971.dta", replace
restore

replace nachname="Döbeli" if nachname=="Dobeli" & vorname=="Ernst" & year==1971
replace nachname="Schwarz-Knecht" if nachname=="Schwarz" & vorname=="Heidi" ///
	& year==1971
replace nachname="Wieser-Nielsen" if nachname=="Wieser" & vorname=="Helga" ///
	& year==1971
replace vorname="Marie-Luise" if nachname=="Fischer" & vorname=="Marie-Louise" ///
	& year==1971
replace nachname="Schmidt-Brugger" if nachname=="Schmidt" & vorname=="Sonja" ///
	& year==1971
replace vorname="Claire-Lise" if nachname=="Renggli-Bonsack" ///
	& vorname=="Claire Lise" & year==1971
replace nachname="Hofmann" if nachname=="Hofman" & vorname=="Fritz" & year==1971
replace nachname="Kaiser" if nachname=="Kaiser-Krämer" & vorname=="Gisela" ///
	& year==1971
replace nachname="Isler" if nachname=="Gisler" & vorname=="Heinz" & year==1971
replace vorname="Hélène" if nachname=="Hirschi-Jeanprêtre" & vorname=="Helene" ///
	& year==1971
replace vorname="Isaac" if nachname=="Gerber" & vorname=="Isaak" & year==1971
replace nachname="Düby" if nachname=="Duby" & vorname=="Hans" & year==1971
replace vorname="Hansrudolf" if nachname=="Lutz" & vorname=="Hans Rudolf" ///
	& year==1971
replace nachname="Schletti-Stössel" if nachname=="Schletti-Stossel" ///
	& vorname=="Lucie" & year==1971
replace nachname="Häusermann-Tièche" if nachname=="Häusermann-Tïèche" ///
	& vorname=="Lucie" & year==1971
replace nachname="Ammann" if nachname=="Amman" & vorname=="Ulrich" & year==1971
replace nachname="Ammann" if nachname=="Amman" & vorname=="Ulrich" & year==1967
replace nachname="Kunz-Aeschlimann" if nachname=="Kunz-Aeschlimarm" ///
	& vorname=="Vreni" & year==1971
replace nachname="Thalmann" if nachname=="Thalmann-Leutenegger" ///
	& vorname=="Ruth" & year==1971
replace nachname="Faust-Kübler" if nachname=="Faust-Kubler" & vorname=="Erika" ///
	& year==1971
replace vorname="Georges" if nachname=="Degen" & vorname=="Georg" & year==1971
replace nachname="Benz-Gürtler" if nachname=="Benz-Gurtler" & vorname=="Jean" ///
	& year==1971
replace nachname="Zimmerli-Silbernagel" if nachname=="Silbernagel-Zimmerli" ///
	& vorname=="Marty" & year==1971
replace nachname="Ruffieux-Overney" if nachname=="Ruffieux Overney" ///
	& vorname=="Monique" & year==1971
replace vorname="Hanni" if nachname=="Schwab" & vorname=="Hanny" & year==1971
replace vorname="Ingrid" if nachname=="Risch" & vorname=="I. Henriette" ///
	& year==1971
replace nachname="Widmer-Schmid" if nachname=="Widmer" & vorname=="Ursula" ///
	& year==1971
replace vorname="Ernst" if nachname=="Müller" & vorname=="Ernst, Jun." ///
	& year==1971
replace nachname="Schlatter" if nachname=="Schlauer" & vorname=="Gaspard" ///
	& year==1971
replace nachname="Klauser-Schwab" if nachname=="Klauser Schwab" ///
	& vorname=="Irma" & year==1971
replace vorname="Jakob Jun." if nachname=="Bütler" & vorname=="Jakob, Jun." ///
	& year==1971
replace nachname="Diethelm-Dobler" if nachname=="Diethelm" & vorname=="Josef" ///
	& year==1959
replace nachname="Diethelm-Dobler" if nachname=="Diethelm" & vorname=="Josef" ///
	& year==1963
replace nachname="Diethelm-Dobler" if nachname=="Diethelm" & vorname=="Josef" ///
	& year==1967
replace nachname="Diethelm-Dobler" if nachname=="Diethelm" & vorname=="Josef" ///
	& year==1971
replace nachname="Diethelm-Dobler" if nachname=="Diethelm" & vorname=="Josef" ///
	& year==1975
replace nachname="Bommeli-Reutlinger" if nachname=="Bommeli" ///
	& vorname=="Elisabeth" & year==1971
replace vorname="Franz-Norbert" if nachname=="Bommer" & vorname=="Franz Norbert" ///
	& year==1971
replace vorname="Hanspeter" if nachname=="Fischer" & vorname=="Hans Peter" ///
	& year==1971
replace nachname="Pinoja" if nachname=="Pinoia" & vorname=="Innocente" ///
	& year==1963
replace nachname="Pinoja" if nachname=="Pinoia" & vorname=="Innocente" ///
	& year==1967
replace nachname="Pinoja" if nachname=="Pinoia" & vorname=="Innocente" ///
	& year==1971
replace vorname="Emile" if nachname=="Chevalier" & vorname=="Emile-Louis" ///
	& year==1971
replace nachname="Ethenoz-Damond" if nachname=="Ethenoz Damond" ///
	& vorname=="Gabrielle" & year==1971
replace nachname="Payot" if nachname=="Payot " & vorname=="Pierre" & year==1971
replace vorname="Anton E." if nachname=="Schrafl" & vorname=="Anton" ///
	& year==1971
replace vorname="Erhard junior" if nachname=="Spörri" & vorname=="Erhard, Jun." ///
	& year==1971
replace vorname="Hans-Ulrich" if nachname=="Frei-Wohlgemuth" ///
	& vorname=="Hans Ulrich" & year==1971
replace vorname="Heinrich" if nachname=="Müller" & vorname=="Heinrich C." ///
	& year==1971
replace nachname="Graemiger" if nachname=="Graeminger" & vorname=="Peter" ///
	& year==1971
replace nachname="Graemiger" if nachname=="Graeminger" & vorname=="Peter" ///
	& year==1975
replace vorname="Rudolf-Christoph" if nachname=="Schellenberg" ///
	& vorname=="Rudolf C." & year==1971
replace vorname="Walter Albert" if nachname=="Peter" & vorname=="Walter" ///
	& geburtsjahr==1938 & year==1971

bysort year canton vorname nachname: gen indimax=_N
*br if indimax>1
gen identifier_1971=.
replace identifier_1971=geburtsjahr if indimax>1 

replace vorname="Hans2" if nachname=="Pfister" & vorname=="Hans" ///
	& geburtsjahr==1903 & year==1963 & kandidst_tot==21526
* Duplicate Hans Pfister in 1963
merge 1:1 year canton vorname nachname identifier_1971 ///
	using "$path\02_Processed_data\01_Elections_1931_1975\Kovariate_1971.dta", ///
	gen(merge_covariates_1971)
tab year if merge_covariates_1971==3
replace vorname="Hans" if nachname=="Pfister" & vorname=="Hans2" ///
	& geburtsjahr==1903 & year==1963 & kandidst_tot==21526
* Duplicate Hans Pfister in 1963
* Check cases with different municipality ids (Simon's coding vs. BfS coding).

*br nachname vorname geburtsjahr canton gemeindenummer_bfs gemeindename ///
*	Gemeinde_1971 Gemeindenummer_BFS_1971 ///
* if gemeindenummer_bfs!=Gemeindenummer_BFS_1971 & year==1971
* Notes: Different municipality ids for
* Lajoux (BE later JU): we did not change anything 
* Mervelier (BE later JU): we did not change anything
* Hasle bei Burgdorf or Rüegsau: we did not change anything b/c two possibili-
* ties (Hasle-Rüegsau in original document)
* Wuppenau: we did not change anything

gen str1 geschlecht="M" if year<1971
replace geschlecht=Geschlecht_1971 if year==1971 

replace gemeindename=Gemeinde_1971 if year==1971 
replace gemeindenummer_bfs=Gemeindenummer_BFS_1971 if year==1971

* (viii) Recode missing origin information

replace origin=gemeindename if origin=="von und" | origin=="de et"	| ///
	origin=="i ed"
	
* (ix) Replace party names by harmonized BfS partynames

replace party=Parteiname_1971 if year==1971

* (x) Rename variables for merge with post-1971 data

keep Kandidatennummer_1971 nachname vorname geburtsjahr canton origin canton_no ///
	kandidst_tot job party party_all year gemeindename gemeindenummer_bfs ///
	LinksRechts elected geschlecht 

rename Kandidatennummer_1971 candidateno
rename kandidst_tot votes
rename canton_no cantonno
rename party_all list
rename gemeindenummer_bfs municipalityno
rename gemeindename municipality
rename LinksRechts leftright
rename geschlecht sex
rename geburtsjahr birthyear 
rename vorname firstname
rename nachname name
rename party partyname		
	
* (xi) 

* 1975
* Note that we save 1975 covariate information from the Bundesblatt which we la-
* ter merge to BfS data
				
replace firstname = "Adolf Kurt" if firstname =="Adolf" & name =="Leemann" ///
	& birthyear ==1927 & canton =="SZ" & votes==482
replace birthyear = 1931 if firstname =="Adolf Kurt" & name =="Leemann" ///
	& birthyear ==1927 & canton =="SZ" & votes==482
replace name = "Müller" if firstname =="Adrian" & name =="Muller" ///
	& birthyear ==1945 & canton =="BL" & votes==3215
replace name = "Ehrhard" if firstname =="Alexander" & name =="Erhard" ///
	& birthyear ==1911 & canton =="ZH" & votes==20213
replace name = "Steiner-Derrer" if firstname =="Alice" & name =="Steiner" ///
	& birthyear ==1916 & canton =="BE" & votes==20059
replace name = "Müller" if firstname =="Andreas" & name =="Mülller" ///
	& birthyear ==1934 & canton =="AG" & votes==9950
replace firstname = "Andre" if firstname =="André" & name =="Friche" ///
	& birthyear ==1940 & canton =="BE" & votes==3097
replace name = "Castelli-Rérat" if firstname =="Andrée" & name =="Castelli" ///
	& birthyear ==1945 & canton =="BE" & votes==1032
replace name = "Keller-Della" if firstname =="Anina" & name =="Keller" ///
	& birthyear ==1939 & canton =="BE" & votes==5009
replace firstname = "Casa Anina" if firstname =="Anina" & name =="Keller-Della" ///
	& birthyear ==1939 & canton =="BE" & votes==5009
replace name = "Deforel-Andreoli" if firstname =="Anna" & name =="Deforel" ///
	& birthyear ==1927 & canton =="BE" & votes==930
replace name = "Albrecht" if firstname =="August" & name =="Albrecht " ///
	& birthyear ==1907 & canton =="NW" & votes==5240
replace name = "Philippe-Waldener" if firstname =="Carla" & name =="Philippe" ///
	& birthyear ==1937 & canton =="BE" & votes==2150
replace name = "Simon-Vermot" if firstname =="Claude" & name =="Simon-Vermod" ///
	& birthyear ==1914 & canton =="NE" & votes==2484
replace name = "Baumgartner" if firstname =="David" & name =="Baumgartner " ///
	& birthyear ==1908 & canton =="GL" & votes==4554
replace name = "Deluc-Schildknecht" if firstname =="Dora" & name =="Deluc" ///
	& birthyear ==1928 & canton =="BE" & votes==6600
replace name = "Wyss-Nüssli" if firstname =="Dora" & name =="Wyss" ///
	& birthyear ==1937 & canton =="BE" & votes==8938
replace name = "Binz-Gehring" if firstname =="Doris" & name =="Binz" ///
	& birthyear ==1941 & canton =="BE" & votes==16647
replace name = "Kläusli" if firstname =="Egon" & name =="Klausli" ///
	& birthyear ==1945 & canton =="ZH" & votes==16341
replace name = "Oehrli-Weber" if firstname =="Elisabeth" & name =="Oehrli" ///
	& birthyear ==1926 & canton =="BE" & votes==16126
replace firstname = "Ernest" if firstname =="Ernst" & name =="Herren" ///
	& birthyear ==1914 & canton =="FR" & votes==1937
replace name = "Ganiere" if firstname =="François" & name =="Ganière" ///
	& birthyear ==1941 & canton =="VD" & votes==5318
replace name = "Marquis-Choulat" if firstname =="Françoise" & name =="Marquis" ///
	& birthyear ==1949 & canton =="BE" & votes==2314
replace birthyear = 1943 if firstname =="Gabrielle" & name =="Nanchen" ///
	& birthyear ==1934 & canton =="VS" & votes==23450
replace firstname = "Gebhard Jun." if firstname =="Gebhard, Jun." ///
	& name =="Forrer" & birthyear ==1949 & canton =="ZH" & votes==4303
replace name = "Merz-Steffen" if firstname =="Gertrud" & name =="Merz" ///
	& birthyear ==1922 & canton =="BE" & votes==11586
replace name = "Kaiser-Kraemer" if firstname =="Gisela" & name =="Kaiser" ///
	& birthyear ==1924 & canton =="BE" & votes==13080
replace firstname = "Gisela" if firstname =="Gisela " & name =="Traub" ///
	& birthyear ==1944 & canton =="BS" & votes==14124
replace name = "Hoffmeyer-Schori" if firstname =="Gréty" & name =="Hoffmeyer" ///
	& birthyear ==1934 & canton =="BE" & votes==2039
replace name = "Lindt-Loosli" if firstname =="Hanna" & name =="Lindt" ///
	& birthyear ==1926 & canton =="BE" & votes==23516
replace firstname = "Hans-Rudolf" if firstname =="Hans Rudolf" & name =="Haegi" ///
	& birthyear ==1938 & canton =="ZH" & votes==32723
replace firstname = "Hans-Ulrich" if firstname =="Hans Ulrich" & name =="Büschi" ///
	& birthyear ==1940 & canton =="BE" & votes==26365
replace firstname = "Hans-Peter" if firstname =="Hanspeter" & name =="Buob" ///
	& birthyear ==1932 & canton =="SG" & votes==7541
replace firstname = "Heinrich Jun." if firstname =="Heinrich, Jun." ///
	& name =="Bäbler" & birthyear ==1951 & canton =="BE" & votes==2061
replace name = "Moser-Wälti" if firstname =="Helene" & name =="Moser" ///
	& birthyear ==1927 & canton =="BE" & votes==21264
replace name = "de Goumoens" if firstname =="Henri" & name =="de Goumoëns" ///
	& birthyear ==1913 & canton =="VD" & votes==13814
replace birthyear = 1941 if firstname =="Herbert" & name =="Dirren" ///
	& birthyear ==1942 & canton =="VS" & votes==12370
replace firstname = "Hermann Ferdinand" if firstname =="Hermann" ///
	& name =="Haller" & birthyear ==1904 & canton =="ZH" & votes==796
replace name = "Mäder-Lüthi" if firstname =="Hertha" & name =="Mäder" ///
	& birthyear ==1921 & canton =="BE" & votes==23654
replace name = "Plomb" if firstname =="Hugues" & name =="Plomp" ///
	& birthyear ==1944 & canton =="BE" & votes==2238
replace name = "Aegerter-Merk" if firstname =="Irene" & name =="Aegerter" ///
	& birthyear ==1940 & canton =="BE" & votes==18490
replace firstname = "Regis Jacqueline" if firstname =="Jacqueline" ///
	& name =="Baeschlin-Régis" & birthyear ==1937 & canton =="AG" & votes==16090
replace name = "Baeschlin" if firstname =="Regis Jacqueline" ///
	& name =="Baeschlin-Régis" & birthyear ==1937 & canton =="AG" & votes==16090
replace name = "Jacquier" if firstname =="Jean-Jacques" & name =="Jaquier" ///
	& birthyear ==1941 & canton =="BE" & votes==9943
replace firstname = "Maxime Jules Xavier" if firstname =="Jules Xavier Maxime" ///
	& name =="Joly" & birthyear ==1912 & canton =="BE" & votes==2842
replace name = "Robert-Bächtold" if firstname =="Leni" & name =="Robert" ///
	& birthyear ==1936 & canton =="BE" & votes==22906
replace birthyear = 1917 if firstname =="Louis" & name =="Barras" ///
	& birthyear ==1914 & canton =="FR" & votes==25099
replace name = "Jacques" if firstname =="Madeleine" & name =="Jaques" ///
	& birthyear ==1949 & canton =="BL" & votes==2761
replace name = "Cattin-Aubry" if firstname =="Mady" & name =="Cattin" ///
	& birthyear ==1930 & canton =="BE" & votes==2106
replace name = "Krähenbuhl" if firstname =="Marco" & name =="Krähenbühl" ///
	& birthyear ==1941 & canton =="TI" & votes==6226
replace name = "Zurlinden-Wymann" if firstname =="Marianne" & name =="Zurlinden" ///
	& birthyear ==1922 & canton =="BE" & votes==12315
replace name = "Hitz-Droz" if firstname =="Marie-Louise" & name =="Hitz" ///
	& birthyear ==1944 & canton =="BE" & votes==7977
replace name = "Oeuvray" if firstname =="Martin" & name =="Œuvray" ///
	& birthyear ==1934 & canton =="BE" & votes==17263
replace firstname = "Matthias" if firstname =="Mathias" & name =="Zimmermann" ///
	& birthyear ==1942 & canton =="ZH" & votes==1916
replace name = "Willfratt" if firstname =="Max" & name =="Willfrat" ///
	& birthyear ==1952 & canton =="ZH" & votes==4055
replace name = "Bettex-Galland" if firstname =="Micheline" & name =="Bettex" ///
	& birthyear ==1928 & canton =="BE" & votes==17313
replace name = "de Lorenzi-Rossini" if firstname =="Milena" ///
	& name =="De Lorenzi-Rossini" & birthyear ==1948 & canton =="TI" ///
	& votes==5532
replace name = "Schaffner-Farine" if firstname =="Monique" & name =="Schaffner" ///
	& birthyear ==1940 & canton =="BE" & votes==3182
replace name = "Bretscher-Bickel" if firstname =="Odette" & name =="Bretscher" ///
	& birthyear ==1921 & canton =="BE" & votes==14373
replace name = "Küng-Gloor" if firstname =="Olga" & name =="Küng" ///
	& birthyear ==1929 & canton =="BE" & votes==9670
replace firstname = "Reymond" if firstname =="Raymond" & name =="Maître" ///
	& birthyear ==1927 & canton =="BE" & votes==2320
replace name = "Gaume-Gigandet" if firstname =="Raymonde" & name =="Gaume" ///
	& birthyear ==1951 & canton =="BE" & votes==3078
replace firstname = "René" if firstname =="Rene" & name =="Agassis" ///
	& birthyear ==1921 & canton =="VD" & votes==12849
replace name = "Clément" if firstname =="René" & name =="Clement" ///
	& birthyear ==1922 & canton =="VD" & votes==27033
replace name = "Jürg" if firstname =="Robert" & name =="Jörg" ///
	& birthyear ==1934 & canton =="ZH" & votes==2306
replace name = "Düberdorfer" if firstname =="Rudolf" & name =="Dübendorfer" ///
	& birthyear ==1908 & canton =="ZH" & votes==760
replace name = "Naegeli-Baur" if firstname =="Ruth" & name =="Naegeli" ///
	& birthyear ==1923 & canton =="BE" & votes==11154
replace name = "Schläppy" if firstname =="Rémy" & name =="Schlappy" ///
	& birthyear ==1917 & canton =="NE" & votes==18128
replace firstname = "Suzanne" if firstname =="Susanne" & name =="Tschäppät" ///
	& birthyear ==1927 & canton =="ZH" & votes==4357
replace name = "Jacquet" if firstname =="Suzanne" & name =="Jaquet" ///
	& birthyear ==1929 & canton =="BE" & votes==9482
replace name = "Lutz-Courant" if firstname =="Suzanne" & name =="Lutz" ///
	& birthyear ==1915 & canton =="BE" & votes==8205
replace name = "Friedli-Grässli" if firstname =="Valentine" & name =="Friedli" ///
	& birthyear ==1929 & canton =="BE" & votes==3741
replace firstname = "Vera M.E." if firstname =="Vera M. E." & name =="Obeid" ///
	& birthyear ==1932 & canton =="ZH" & votes==20114
replace name = "Weber" if firstname =="Vreni" & name =="Weber-Brändli" ///
	& birthyear ==1946 & canton =="SG" & votes==1098
replace firstname = "Wally" if firstname =="Wally " & name =="Widmer" ///
	& birthyear ==1925 & canton =="ZH" & votes==12409
replace name = "Huber-Ravazzi" if firstname =="Walter" & name =="Huber Ravazzi" ///
	& birthyear ==1938 & canton =="ZH" & votes==1252
replace name = "Jäger-Stamm" if firstname =="Walter" & name =="Jaeger-Stamm" ///
	& birthyear ==1911 & canton =="GR" & votes==2617
replace firstname = "Stöckli" if firstname =="Walter" & name =="Stöckli" ///
	& birthyear ==1942 & canton =="UR" & votes==1596
replace name = "Walter" if firstname =="Stöckli" & name =="Stöckli" ///
	& birthyear ==1942 & canton =="UR" & votes==1596
replace firstname = "Willy" if firstname =="Willi" & name =="Walker" ///
	& birthyear ==1922 & canton =="ZH" & votes==44829

preserve
keep if year==1975
keep year canton name firstname birthyear origin job votes
rename votes votes_bundesblatt1975
merge 1:1 year canton name firstname birthyear using ///
	"$path\02_Processed_data\02_Elections_1971_2015\NRW_KANDIDATEN 1975-2015.dta", ///
	gen(merge_1975)
save "$path\02_Processed_data\02_Elections_1971_2015\NRW_KANDIDATEN 1975-2015.dta", ///
	replace
tab year if merge_1975==2
tab canton if merge_1975==2 & year==1975
restore

drop if year==1975
* Drop 1975 b/c we use BFS data for this year.

erase "$path\02_Processed_data\01_Elections_1931_1975\temp.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\Kantonsnamen_harmonisiert.dta"
erase "$path\02_Processed_data\06_Parties\LeftRightCoding.dta"
erase "$path\02_Processed_data\06_Parties\Parteinamen All Florence.dta"
erase "$path\02_Processed_data\08_Municipalities\Gemeindebestand_1969.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\Kovariate_1971.dta"


******************************************************
* (D) MERGE PRE- AND POST-1975-DATA
******************************************************

* (i) Append all data and labeling

append 	using ///
	"$path\02_Processed_data\02_Elections_1971_2015\NRW_KANDIDATEN 1975-2015.dta"
* Append BfS data including the previously added information on origin and job 
* from Bundesblatt 1975

order canton* year name firstname birthyear sex votes elected candidateno ///
	municipality municipalityno leftright list job status

label var canton "Canton"
label var cantonno "Canton no."
label var year "Year"
label var name "Name"
label var firstname "First name"
label var birthyear "Year of birth"
label var sex "Sex"
label var votes "Votes"
label var elected "Elected"
label var candidateno "Candidate no. on list"
label var municipality "Municipality"
label var municipalityno "Municipality no."
label var leftright "Left-right coding of list"						
label var list "List"
label var listno "List number"
label var job "Job"
label var status "Status"
label var origin "Origin"

* (ii) Quality checks

* (a) Is Municipality id constant per municipality? 

bysort canton municipality: egen municipality_sd=sd(municipalityno) 			
*br if municipality_sd>0 & !missing(municipality_sd)
* Note: We saved this dataset to check (Quality Check Gemeinde für Florence.xlsx)		
drop municipality_sd
* Result: 164 observations need correction, action: run do file with our correc-
* tions

do "$path\03_Code\08_Municipalities\Correct Municipalities 1975-2015 part I.do"
* Note: Corrects municipality numbers.

* (b) Check missing values
* Total observations in dataset: 41'606

count if missing(cantonno)
* Result: 0 non-missings, action: nothing
count if missing(year)
* Result: 0 non-missings, action: nothing
count if missing(name)
* Result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
tab canton if missing(name)
count if missing(firstname)
* Result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
tab canton if missing(firstname)
count if missing(birthyear)
* Result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
count if missing(sex)
* Result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
count if missing(votes)
sort year canton name firstname
count if missing(votes)
*br if missing(votes) & !missing(birthyear)
* Result: 92 observations with no votes, all b/c "Stille Wahlen" (57 individuals)
* in 1939, 35 individuals in the period of 1945-2015)
count if missing(elected)
*br if missing(elected)
* Result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"
count if missing(candidateno)
tab year if missing(candidateno)
* Result: missing before 1971, action: nothing
tab year if !missing(candidateno)
count if municipality==""
* Result: 53 non-missings, action: nothing b/c these are "Vereinzelte Kandidaten"

preserve
count if missing(municipalityno)
count if missing(municipalityno) & municipality!=""	
*br if missing(municipalityno)
* We collapsed this subsample (by canton and municipality) to construct the out-
* file for Florence (Municipalities All Florence.xlsx)
keep if missing(municipalityno)
collapse municipalityno, by(municipality canton)
restore

do "$path\03_Code\08_Municipalities\Correct Municipalities 1975-2015 part II.do"
* Note: Corrects municipality names.

preserve
* Read in Gemeindestand by Simon Lüchinger
import exc "$path\01_Raw_data\08_Municipalities\Gemeindebestand_20151231.xls", ///
	first clear  cellra(A12:F2336)
rename BFSGdenummer municipalityno_2015
rename Gemeindename municipality
keep municipalityno_2015 municipality
save "$path\02_Processed_data\08_Municipalities\Gemeindebestand_20151231.dta", ///
	replace
restore

merge m:1 municipality using ///
	"$path\02_Processed_data\08_Municipalities\Gemeindebestand_20151231.dta"
* Note: Many mismatches from master and using ok.
drop if _merge==2

replace municipalityno=municipalityno_2015 if missing(municipalityno)
drop _merge municipalityno_2015

count if missing(municipalityno)
* Result: 53 Vereinzelte have missing municipalityno.
count if missing(leftright) 
tab canton if missing(leftright)
* Result: all 85 observations with missing data are in AI, NW, OW, UR. 
*br if missing(leftright)
* Action: We do not replace the missing leftright score b/c the files for these
* cantons are small (<39 observations).
count if list=="" 
*br if list=="" & !missing(leftright)
* Result: Same set of individuals that also lacks leftright coding; action: no
* change.
count if job=="" 
tab year if job==""
* Result: No jobs for those>1975, only about 200 observations with missing job
* information for year<1975; action: no change .
tab year if status==""
* Result: No status available before 1975.
tab canton if status=="" & year>=1975
* 53 "Vereinzelte"
*br if status=="" & year>=1975

* (c) Compare votes in 1975 across sources (BfS vs. Bundesblatt)

gen votes_diff=votes-votes_bundesblatt1975
tab votes_diff
* Result: 82% with identical result, two large differences of 5096 and 1616
* votes, all other observations <110 votes difference.
	
erase "$path\02_Processed_data\08_Municipalities\Gemeindebestand_20151231.dta"

******************************************************
* (E) SAVE DATASET FOR FUZZY MERGE (1931 until 2015)
******************************************************

*save "07_Archive\03_Candidates_fuzzy_mrg\nationalraete_1931_2015_fuzzy.dta", ///
*	replace


******************************************************
* (F) MERGE IDS AND ELECTION RESULTS (1931 until 2015)
******************************************************

* (i) Merge BE and JU to one canton

gen canton_merge=canton
replace canton_merge="BEJU" if canton=="BE" | canton=="JU"	

* (ii) Correct birthyears and names (from fuzzy merge file)

replace firstname="Patrick" if name=="Mächler" & birthyear==1983 & year==2011
replace firstname="Massimiliano" if name=="Ay" & birthyear==1982 & year==2011
replace firstname="Sergio" if name=="Arigoni" & birthyear==1965 & year==2011
replace firstname="Werner 2" if name=="Müller" & birthyear==1932 & year==1975 ///
	& firstname=="Werner"

rename name name_merge
* Note: Undo correct changes by Jana Jarck (in 1st check round as saved in 
* erf_luzern_all.dta) for merge to election results.
rename firstname firstname_merge
rename birthyear birthyear_merge 

* (ii) Merge electoral results and fuzzy merge

merge 1:1 firstname_merge name_merge canton_merge birthyear_merge municipality ///
	year using ///
	"$path\01_Raw_data\03_Candidates_fuzzy_mrg\FuzzyMerge_Ultimate.dta", ///
	gen(merge_fuzzy_election_results)
* Result: 59 from master but not in using (6 candidates that were elected in 
* "Stillen Wahlen" (see below) and 53 "Vereinzelte".

* Replace IDs of those candidates that were elected in "stille Wahlen" after
* 1975

replace id_Ultimate="AR-1975-0004" if canton=="AR" ///
	& name_merge=="Merz" & firstname_merge=="Christian" & year==1979
* Christian Merz (stille Wahlen)
replace id_Ultimate="AR-1975-0002" if canton=="AR" & name_merge=="Früh" ///
	& year==1979
* Hans-Rudolf Früh (stille Wahlen)
replace id_Ultimate="AR-1975-0002" if canton=="AR" & name_merge=="Früh" ///
	& year==1987
* Hans-Rudolf Früh (stille Wahlen)
replace id_Ultimate="AR-1983-0004" if canton=="AR" & name_merge=="Maeder" ///
	& year==1987
* Herbert Mäder (stille Wahlen)
replace id_Ultimate="NW-1995-0002" if canton=="NW" & name_merge=="Engelberger" ///
	& year==2007
* Edi Engelberger(stille Wahlen)
replace id_Ultimate="OW-1995-0002" if canton=="OW" & name_merge=="Durrer" ///
	& year==1999
* Adalbert Durrer(stille Wahlen)
replace firstname="Werner" if name_merge=="Müller" & birthyear==1932 ///
	& year==1975 & firstname_merge=="Werner 2"

* (iii) Incorporate corrections and unclear cases from final round by JJ and RB

* Correct birthyears, names and firstnames (from unsicher coding of JJ and RB in
* the final round)

replace birthyear=1961 if id_Ultimate=="AG-1991-0052" & year==1995
replace birthyear=1878 if id_Ultimate=="BEJU-1935-0165" & year==1939
replace birthyear=1878 if id_Ultimate=="BEJU-1935-0165" & year==1943
replace birthyear=1915 if id_Ultimate=="BEJU-1947-0030" & year==1971
replace firstname="Paul" if id_Ultimate=="BEJU-1951-9104" & year==1959
* Wrong name; according to our cantonal dataset there was never an Emil Schorer
* in the Berner Grossrat.
replace birthyear=1950 if id_Ultimate=="GR-1987-0001" & year==2011 
replace birthyear=1949 if id_Ultimate=="LU-1975-7053" & year==1975
replace birthyear=1894 if id_Ultimate=="SG-1935-9026" & year==1935
replace birthyear=1882 if id_Ultimate=="SG-1935-9026" & year==1935
replace birthyear=1891 if id_Ultimate=="VD-1931-9026" & year==1931
replace birthyear=1878 if id_Ultimate=="VS-1931-7097" & year==1943
replace birthyear=1924 if id_Ultimate=="VS-1955-0027" & year==1967
replace name="Pfister" if id_Ultimate=="ZH-1955-7149" & year==1955
* Wrong name: Siegfried Pfister

* Unclear birthyears:
* Ernst Zubler: 1894 or 1884, years: 1935 and 1939,
* id_Ultimate_JJ=="AG-1935-9027" 
* Paul Meier: 1906 or 1907, years: 1955 and 1963,
* id_Ultimate_JJ=="BEJU-1955-9070" 
* Hans Beyeler: 1910 or 1911, years: 1971 and 1975,
* id_Ultimate_JJ=="BEJU-1971-9027" 
* Florian Gantenbein: 1896 or 1903, years: 1943 and 1947,
* id_Ultimate_JJ=="-1943-9009"
* Jean Curdy, 1922 or 1927, years: 1967 and 1971, id_Ultimate_JJ=="VD-1967-7089"
* Hans Rohner, 1928 or 1932, years: 1983 and 1987,
* id_Ultimate_JJ=="ZH-1983-9116"
* Markus Hug, 1964 or 1969, years: 1991 and 1995,
* id_Ultimate_JJ=="ZH-1991-0308"
* Jürg (Giorgio) Bösiger, 1955 or 1963, years: 1999 and 2003,
* id_Ultimate_JJ=="ZH-1999-9184"

* Unclear coding	
* Beat Balmer: 1958, years: 1991 and 2007, id_Ultimate_RB=="BEJU-1991-7004",
* different profession, different party

* (iv) Harmonization of firstnames, names, and birthyears

* Rules: 1. If more than two observations: follow majority spelling and/or cor-
*        rect spelling in Bundesblatt, cantonal data, or other external sources
*        2. If only two observations: check spelling in Bundesblatt (or cantonal
*        data and other external sources) and if they differ we do not change
*        anything
*        3. We do not harmonize Umlaute and double names
*		 4. Change abbrevations to full names (Joh. Friedr. becomes Johann
*        Friedrich) if possible
*		 5. Omit "Jun./sen." in first name
*		 6. Change first names in paranthesis to separate names (Andreas (Res)
*        becomes Res)
*        7. Change "v." to "von", "Von" to "von", "Van" to "van, "De" to "de"
*		 8. No change in first name if first name is used in some years but not
*        in others (e.g., Josef Bernhard and Bernhard)

* (I) FIRST NAMES 

rename id_Ultimate ID

* (0) Several pre-tests

* Check whether canton varies within ID

egen canton_num=group(canton)	
bysort ID: egen canton_num_sd=sd(canton_num)
*br if canton_num_sd>0 & !missing(canton_num_sd)
tab canton if canton_num_sd>0 & !missing(canton_num_sd) & ID!=""
* Result: Canton varies only in BE and JU.

* Check whether we find the same ID twice in a year

bysort ID year: gen indi_ID_year=_n if ID!=""
tab indi_ID_year
* Result: No ID is found twice in a year.

* Check whether we find the same ID twice in a year

* (a) Check wrong entries
* Result: only Walter Stöckli (UR, 1975)

* (b) Check cases for which firstname differs within the same ID

egen firstname_num=group(firstname)	
bysort ID: egen firstname_num_sd=sd(firstname_num)
*br if firstname_num_sd>0 & !missing(firstname_num_sd)
* Result: 1842 with non-identical first names

* (c) Change cases with parenthesis in firstname	
* Rule: change obvious abbreviations, leave cases with additional names (e.g.
* Alfredo (Benjamin), Claudia (Dana) etc)

split firstname, p("(")
*br name firstname birthyear cantonno year if firstname2!=""

* (d) Check unique firstnames manually

preserve
bysort firstname: gen indi=_n
keep if indi==1
*br firstname
restore

* (e) Check cases with only one letter in first name (sheet abbreviation in
* check names.xlsx) 

*br firstname name canton year birthyear list job if firstname=="A." | ///
*	firstname=="B." | firstname=="C." | firstname=="D." | firstname=="E." | ///
*	firstname=="F." | firstname=="G." | firstname=="H." | firstname=="I." | ///
*	firstname=="J." | firstname=="K." | firstname=="L." | firstname=="M." | ///
*	firstname=="N." | firstname=="O." | firstname=="P." | firstname=="Q." | ///
*	firstname=="R." | firstname=="S." | firstname=="T." | firstname=="U." | ///
*	firstname=="V." | firstname=="W." | firstname=="X." | firstname=="Y." | ///
*	firstname=="Z."

* (II) NAMES

egen name_num=group(name)
bysort ID: egen name_num_sd=sd(name_num)
*br if name_num_sd>0 & !missing(name_num_sd)
* Result: 2212 cases with non-identical names.

* (III) BIRTHYEARS 
* Note: changes are stored in birthyears.xlsx

bysort ID: egen birthyear_num_sd=sd(birthyear)
*br if birthyear_num_sd>0 & !missing(birthyear_num_sd)	
* Result: 233 with non-identical birthyears

do "$path\03_Code\01_Elections_1931_1975\Correct Firstnames, Names, Birthyears.do" 

* (v) Add information about six cases elected in stille Wahlen (not included in
* BfS data and thus not used for fuzzy merge)

replace firstname="Hans-Rudolf" if name_merge=="Früh" & canton=="AR" ///
	& year==1979
replace firstname="Hans-Rudolf" if name_merge=="Früh" & canton=="AR" ///
	& year==1987
replace name="Früh" if name_merge=="Früh" & canton=="AR" & year==1979
replace name="Früh" if name_merge=="Früh" & canton=="AR" & year==1987
replace birthyear=1936 if name_merge=="Früh" & canton=="AR" & year==1979
replace birthyear=1936 if name_merge=="Früh" & canton=="AR" & year==1987

replace firstname="Christian" if name_merge=="Merz" & canton=="AR" & year==1979
replace name="Merz" if name_merge=="Merz" & canton=="AR" & year==1979
replace birthyear=1943 if name_merge=="Merz" & canton=="AR" & year==1979

replace firstname="Herbst" if name_merge=="Maeder" & canton=="AR" & year==1987
replace name="Maeder" if name_merge=="Maeder" & canton=="AR" & year==1987
replace birthyear=1930 if name_merge=="Maeder" & canton=="AR" & year==1987

replace firstname="Adalbert" if name_merge=="Durrer" & canton=="OW" & year==1999
replace name="Durrer" if name_merge=="Durrer" & canton=="OW" & year==1999
replace birthyear=1950 if name_merge=="Durrer" & canton=="OW" & year==1999
					
replace firstname="Edi" if name_merge=="Engelberger" & canton=="NW" ///
	& year==2007
replace name="Engelberger" if name_merge=="Engelberger" & canton=="NW" ///
	& year==2007
replace birthyear=1940 if name_merge=="Engelberger" & canton=="NW" & year==2007
		
*br if firstname==""
*Note: Browse those with missing observations in firstname, name, and birthyear.
*br if name==""
*br if birthyear==.

* (vi) Save final dataset

keep ID canton cantonno year name firstname birthyear sex votes elected ///
	candidateno municipality municipalityno partyname leftright list listno job ///
	origin status votes_bundesblatt1975
order ID canton cantonno year name firstname birthyear sex votes elected ///
	candidateno municipality municipalityno partyname leftright list listno job ///
	origin status votes_bundesblatt1975
sort canton year ID

* (vii) Merge Heimatorte and Jobs from Bundesblatt (see email from Florence
* 8.12.2017)

merge 1:1 ID year canton using ///
	"$path\01_Raw_data\01_Elections_1931_1975\origins_and_jobs_tomerge.dta"
* Note: According to email from Florence from 8.12.2017 there are two missing 
* candidates in the BfS data for which we add Bundesblatt information.

replace cantonno=4 if ID=="UR-1995-9001"
replace name="Merminod" if ID=="UR-1995-9001"
replace firstname="Yves" if ID=="UR-1995-9001"
replace sex="M" if ID=="UR-1995-9001"
replace birthyear=1947 if ID=="UR-1995-9001"
replace votes=25 if ID=="UR-1995-9001"
replace elected=0 if ID=="UR-1995-9001"
replace municipality="Neuchâtel" if ID=="UR-1995-9001"
replace municipalityno=6458 if ID=="UR-1995-9001"
replace job="Jurist" if ID=="UR-1995-9001"

replace cantonno=6 if ID=="OW-1995-9001"
replace name="Durrer" if ID=="OW-1995-9001"
replace firstname="Rudolf O." if ID=="OW-1995-9001"
replace birthyear=1951 if ID=="OW-1995-9001"
replace sex="M" if ID=="OW-1995-9001"
replace votes=136 if ID=="OW-1995-9001"
replace elected=0 if ID=="OW-1995-9001"
replace municipality="Sarnen" if ID=="OW-1995-9001"
replace municipalityno=1407 if ID=="OW-1995-9001"
replace job="Student" if ID=="OW-1995-9001"

* Note: Add Rudolf Günthardt, ZH, 1987
replace job="Einkaufschef" if ID=="ZH-1987-0261" & year==1987
replace votes_BU=1442 if ID=="ZH-1987-0261" & year==1987
replace Gemeindename1="Küsnacht (ZH)" if ID=="ZH-1987-0261" & year==1987
replace BFSGdenummer1=154 if ID=="ZH-1987-0261" & year==1987

forv i=1(1)6{
rename Gemeindename`i' origin`i' 
rename BFSGdenummer`i' originno`i' 
}

replace votes_BU = votes_bundesblatt1975 if year==1975 
replace votes_BU = votes if year==1975 & votes_BU==.
* Note: Add Bundeblatt information for Vereinzelte in AR,UR,OW,NW,GL in 1975.

replace job=job_BU if job_BU!="" & job==""

drop birthyear_BU
* Note: We already corrected birthyear (see above).
drop job_BU votes_bundesblatt1975 _merge
drop origin
* Note: Drop pre-1975 origin information from Bundesblatt b/c we have added ori-
* gin information for all years (by Florence)

sort ID year

* (viii) Municipality corrections by Mark (Email 4.10.2018)

replace municipalityno = 3404 if year == 1931 & ID == "SG-1931-0024"
replace municipalityno = 3404 if year == 1931 & ID == "SG-1931-0052"
replace municipalityno = 3404 if year == 1935 & ID == "SG-1935-0035"
replace municipalityno = 3404 if year == 1935 & ID == "SG-1931-0024"
replace municipalityno = 3404 if year == 1935 & ID == "SG-1935-0055"
replace municipalityno = 3404 if year == 1935 & ID == "SG-1935-0004"
replace municipalityno = 3404 if year == 1939 & ID == "SG-1931-0024"
replace municipalityno = 3404 if year == 1939 & ID == "SG-1939-0009"
replace municipalityno = 3404 if year == 1939 & ID == "SG-1939-0044"
replace municipalityno = 3404 if year == 1939 & ID == "SG-1935-0055"
replace municipalityno = 3404 if year == 1943 & ID == "SG-1931-0024"
replace municipalityno = 3404 if year == 1943 & ID == "SG-1935-0004"
replace municipalityno = 3404 if year == 1943 & ID == "SG-1935-0055"
replace municipalityno = 3404 if year == 1943 & ID == "SG-1939-0009"
replace municipalityno = 3404 if year == 1947 & ID == "SG-1947-0017"
replace municipalityno = 3404 if year == 1947 & ID == "SG-1947-0019"
replace municipalityno = 3404 if year == 1947 & ID == "SG-1935-0004"
replace municipalityno = 3404 if year == 1947 & ID == "SG-1947-0012"
replace municipalityno = 3404 if year == 1947 & ID == "SG-1935-0055"
replace municipalityno = 3404 if year == 1947 & ID == "SG-1947-0049"
replace municipalityno = 3404 if year == 1951 & ID == "SG-1947-0019"
replace municipalityno = 3404 if year == 1951 & ID == "SG-1947-0017"
replace municipalityno = 3404 if year == 1951 & ID == "SG-1935-0055"
replace municipalityno = 3404 if year == 1955 & ID == "SG-1935-0055"
replace municipalityno = 3404 if year == 1955 & ID == "SG-1947-0017"
replace municipalityno = 3404 if year == 1955 & ID == "SG-1947-0019"
replace municipalityno = 3404 if year == 1959 & ID == "SG-1935-0055"
replace municipalityno = 3404 if year == 1963 & ID == "SG-1963-0003"
replace municipalityno = 3404 if year == 1963 & ID == "SG-1963-0033"
replace municipalityno = 695 if year == 1971 & ID == "BEJU-1971-0333"
replace municipalityno = 2805 if year == 1975 & ID == "BS-1975-0069"
replace municipalityno = 2811 if year == 1975 & ID == "ZH-1975-0578"
replace municipalityno = 2812 if year == 1975 & ID == "BS-1975-0043"
replace municipalityno = 4507 if year == 1975 & ID == "TG-1975-0032"
replace municipalityno = 4518 if year == 1975 & ID == "TG-1975-0034"
replace municipalityno = 511 if year == 1975 & ID == "BEJU-1971-0127"
replace municipalityno = 515 if year == 1975 & ID == "BEJU-1975-0047"
replace municipalityno = 2203 if year == 1979 & ID == "FR-1979-0002"
replace municipalityno = 4507 if year == 1979 & ID == "TG-1975-0032"
replace municipalityno = 4507 if year == 1979 & ID == "TG-1979-0041"
replace municipalityno = 4534 if year == 1979 & ID == "TG-1971-0022"
replace municipalityno = 4762 if year == 1979 & ID == "TG-1979-0022"
replace municipalityno = 4762 if year == 1979 & ID == "TG-1979-0023"
replace municipalityno = 252 if year == 1983 & ID == "ZH-1983-0300"
replace municipalityno = 252 if year == 1983 & ID == "ZH-1983-0052"
replace municipalityno = 252 if year == 1983 & ID == "ZH-1963-2224"
replace municipalityno = 3425 if year == 1983 & ID == "SG-1983-0043"
replace municipalityno = 3776 if year == 1983 & ID == "GR-1983-0010"
replace municipalityno = 4912 if year == 1983 & ID == "TG-1983-0003"
replace municipalityno = 4507 if year == 1983 & ID == "TG-1979-0041"
replace municipalityno = 4507 if year == 1983 & ID == "TG-1975-0032"
replace municipalityno = 4534 if year == 1983 & ID == "TG-1971-0022"
replace municipalityno = 623 if year == 1987 & ID == "BEJU-1987-0337"
replace municipalityno = 854 if year == 1987 & ID == "BEJU-1987-0510"
replace municipalityno = 854 if year == 1987 & ID == "BEJU-1987-0161"
replace municipalityno = 2809 if year == 1987 & ID == "SO-1987-0053"
replace municipalityno = 3377 if year == 1987 & ID == "SG-1987-0030"
replace municipalityno = 4912 if year == 1987 & ID == "TG-1983-0003"
replace municipalityno = 4507 if year == 1987 & ID == "TG-1975-0032"
replace municipalityno = 4637 if year == 1987 & ID == "TG-1987-0026"
replace municipalityno = 4752 if year == 1987 & ID == "TG-1987-0031"
replace municipalityno = 4752 if year == 1987 & ID == "TG-1987-0022"
replace municipalityno = 5116 if year == 1987 & ID == "TI-1987-0038"
replace municipalityno = 6473 if year == 1987 & ID == "NE-1987-0018"
replace municipalityno = 171 if year == 1991 & ID == "ZH-1991-0188"
replace municipalityno = 4541 if year == 1991 & ID == "TG-1991-0042"
replace municipalityno = 4637 if year == 1991 & ID == "TG-1987-0026"
replace municipalityno = 4762 if year == 1991 & ID == "TG-1991-0041"
replace municipalityno = 4762 if year == 1991 & ID == "TG-1991-0028"
replace municipalityno = 4762 if year == 1991 & ID == "TG-1991-0027"
replace municipalityno = 3943 if year == 2007 & ID == "GR-2007-0053"
replace municipalityno = 306 if year == 2011 & ID == "BEJU-2003-0371"
replace municipalityno = 3943 if year == 2011 & ID == "GR-2011-0009"
replace municipalityno = 3943 if year == 2011 & ID == "GR-2011-0015"
replace municipalityno = 3943 if year == 2011 & ID == "GR-2011-0049"
replace municipalityno = 3943 if year == 2011 & ID == "GR-2007-0053"
replace municipalityno = 875 if year == 2015 & ID == "FR-2015-0049"
replace municipalityno = 1630 if year == 2015 & ID == "GL-2011-0003"
replace municipalityno = 1631 if year == 2015 & ID == "GL-2015-9000"

replace municipality = " Götighofen " if year == 1975 & ID == "TG-1975-0032"
replace municipality = "Zihlschacht" if year == 1975 & ID == "TG-1975-0034"
replace municipality = "Lossy-Formangueires" if year == 1979 ///
	& ID == "FR-1979-0002"
replace municipality = "Götighofen" if year == 1979 & ID == "TG-1975-0032"
replace municipality = "Götighofen" if year == 1979 & ID == "TG-1979-0041"
replace municipality = "Unterschlatt" if year == 1979 & ID == "TG-1971-0022"
replace municipality = "Vicosprano" if year == 1983 & ID == "GR-1983-0010"
replace municipality = "Donzhausen" if year == 1983 & ID == "TG-1983-0003"
replace municipality = "Götighofen" if year == 1983 & ID == "TG-1979-0041"
replace municipality = "Götighofen" if year == 1983 & ID == "TG-1975-0032"
replace municipality = "Unterschlatt" if year == 1983 & ID == "TG-1971-0022"
replace municipality = "Donzhausen" if year == 1987 & ID == "TG-1983-0003"
replace municipality = "Götighofen" if year == 1987 & ID == "TG-1975-0032"
replace municipality = "Siegershausen" if year == 1987 & ID == "TG-1987-0026"
replace municipality = "Wilen bei Will" if year == 1987 & ID == "TG-1987-0031"
replace municipality = "Wilen bei Will" if year == 1987 & ID == "TG-1987-0022"
replace municipality = "Magadino" if year == 1987 & ID == "TI-1987-0038"
replace municipality = "Chézard-Saint-Martin" if year == 1987 ///
	& ID == "NE-1987-0018"
replace municipality = "Diessenhofen" if year == 1991 & ID == "TG-1991-0042"
replace municipality = "Siegershausen" if year == 1991 & ID == "TG-1987-0026"

* (ix) Correct IDs from fuzzy merge failures (based on Simon's birthdate and
* death date check, Sep-Nov 2018 and status recoding by Simon/Lukas May 2019)

replace ID="AG-1955-0007" if ID=="AG-1971-9103"
replace ID="BL-1931-9003" if ID=="BL-1935-9004"
replace ID="GE-1931-0001" if ID=="GE-1935-9024"
replace ID="GE-1931-9012" if ID=="GE-1939-9045"
replace ID="GE-1931-9012" if ID=="GE-1947-0021"
replace ID="GE-1963-0026" if ID=="GE-1967-9021"
replace ID="SG-1935-0018" if ID=="SG-1939-0022"
replace ID="ZH-1983-9113" if ID=="ZH-1991-9112"
replace ID="AG-1959-0050" if ID=="AG-1963-0056" 
replace ID="GE-1951-9034" if ID=="GE-1963-9077"
replace ID="GE-1939-9040" if ID=="GE-1943-0006"
replace ID="GE-1947-0004" if ID=="GE-1959-9036"
replace ID="GE-1947-9055" if ID=="GE-1963-9216"
replace ID="SO-1931-9001" if ID=="SO-1943-0033"
replace ID="ZH-1963-9058" if ID=="ZH-1967-0223"

replace municipality="Murten" if ID=="FR-1971-0004"
replace job="pubblicista" if ID=="TI-1935-0015"
replace birthyear=1973 if ID=="VD-2011-0003"

erase "$path\02_Processed_data\02_Elections_1971_2015\NRW_KANDIDATEN 1975-2015.dta"

****************************************************************
* (G) MERGE INFORMATION ON DEATH DATES
****************************************************************

preserve
import exc ///
	"$path\01_Raw_data\09_DoB_and_DoD\GeburtsTodesdaten_Check_Round_1_and_2_20181203.xlsx", ///
	first  cellra(A1:X3386) clear

drop if ID=="AG-1971-9103"
drop if ID=="BL-1935-9004"
drop if ID=="GE-1935-9024"
drop if ID=="GE-1939-9045"
drop if ID=="GE-1947-0021"
drop if ID=="GE-1967-9021"
drop if ID=="SG-1939-0022"
drop if ID=="ZH-1991-9112"
drop if ID=="GE-1943-0006"
drop if ID=="GE-1959-9036"
drop if ID=="GE-1963-9216"
drop if ID=="SO-1943-0033"
* Note: Drop all observations that we have recoded just above, we checked all
* birth and death information and found that these are equal (in (F)).

keep ID bday bmonth byear QuelleGeburtsdatumNeu lebend dday dmonth dyear ///
	QuelleTodesdatumNeu
save "$path\02_Processed_data\09_DoB_and_DoD\GeburtsTodesdaten_Check_Round_1_and_2.dta", ///
	replace
restore	

merge m:1 ID using ///
	"$path\02_Processed_data\09_DoB_and_DoD\GeburtsTodesdaten_Check_Round_1_and_2.dta", ///
	gen(Merge_Round_1_2)

gen death_info =0
replace death_info =1 if lebend==1 | dyear!=.

erase "$path\02_Processed_data\09_DoB_and_DoD\GeburtsTodesdaten_Check_Round_1_and_2.dta"


*******************************************************************
* (H) MERGE INFORMATION ON ALLIANCES, SUBALLIANCES, AND PARTY VOTES
*******************************************************************
		
* Normalize list names

sort list
gen liste_norm = list
forv i = 1(1)40 {
	replace liste_norm = usubinstr(liste_norm, "Liste `i'. ", "", .)
	replace liste_norm = usubinstr(liste_norm, "Lista `i'. ", "", .)
	replace liste_norm = usubinstr(liste_norm, "Liste `i'.", "", .) 
	replace liste_norm = usubinstr(liste_norm, "Lista `i'.", "", .) 
	}
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
replace liste_norm = "LISTE EIDGENOSSISCHE FRONT" if ID =="ZH-1931-0035" ///
	& canton == "ZH" & year == 1931
replace liste_norm = "LISTE EIDGENOSSISCHE FRONT" if ID =="ZH-1931-0024" ///
	& canton == "ZH" & year == 1931
replace liste_norm = "LISTE EIDGENOSSISCHE FRONT" if ID =="ZH-1931-0064" ///
	& canton == "ZH" & year == 1931
replace liste_norm = "LISTE EIDGENOSSISCHE FRONT" if ID =="ZH-1931-0134" ///
	& canton == "ZH" & year == 1931
replace liste_norm = "LISTE EIDGENOSSISCHE FRONT" if ID =="ZH-1931-0168" ///
	& canton == "ZH" & year == 1931
replace liste_norm = "KEIN LISTENNAME" if list == ""

replace liste_norm=name if year==1943 & (canton=="NW" | canton=="OW")
* Note: Replace missing list information with candidate names for four cases.

replace liste_norm="VEREINZELTE" if ID=="" & year==2011 & canton=="AI" & ///
	liste_norm=="KEIN LISTENNAME"
drop if ID=="" & year==2011 & canton=="AI" & liste_norm=="KEIN LISTENNAME" & ///
	votes==.
replace liste_norm="VEREINZELTE" if ID=="" & year==2011 & canton=="AR" & ///
	liste_norm=="KEIN LISTENNAME"
drop if ID=="" & year==2011 & canton=="AR" & liste_norm=="KEIN LISTENNAME" & ///
	votes==.
replace liste_norm="VEREINZELTE" if ID=="" & year==2011 & canton=="GL" & ///
	liste_norm=="KEIN LISTENNAME"
drop if ID=="" & year==2011 & canton=="GL" & liste_norm=="KEIN LISTENNAME" & ///
	votes==.
replace liste_norm="VEREINZELTE" if ID=="" & year==2011 & canton=="UR" & ///
	liste_norm=="KEIN LISTENNAME"
drop if ID=="" & year==2011 & canton=="UR" & liste_norm=="KEIN LISTENNAME" & ///
	votes==.

preserve
do "$path\03_Code\05_Alliances\NR_Listenstimmen.do"
restore

merge m:1 canton year liste_norm using ///
	"$path\02_Processed_data\05_Alliances\Listenverbindungen_1931-2015.dta", ///
	gen(merge_alliance)
drop if merge_alliance==2

* 118 obseravations in master data only (merge==1): candidate observations without 
* cantonal information.
* Reason: Candidate data provide individual list names, while cantonal election 
* information (BfS) contains no list name -- > "KEIN LISTENNAME"
* 12 observation in using data only (merge==2): years without a list called 
* "KEIN LISTENNAME" in both candidate and cantonal datasets.
* Cantons affected: AI, AR, GL, NW, OW, UR, ZG
* Statement of problem: cantonal information on no of eligibles, no of voters, 
* valid votes, empty votes, etc. are not merged on all observations in candidate data
* Solution: identify missing canton-year observations and merge cantonal voting 
* information separately

preserve 
// identify affected canton-years 
keep if merge_alliance != 3
keep canton year
duplicates drop canton year, force
save "$path\02_Processed_data\05_Alliances\nomatch_temp.dta", replace
// extract cantonal voting information
use "$path\02_Processed_data\05_Alliances\Listenverbindungen_1931-2015.dta", clear
merge m:1 canton year using "$path\02_Processed_data\05_Alliances\nomatch_temp.dta", ///
	gen(merge_nomatch)
keep if merge_nomatch == 3
drop merge_nomatch
keep canton year seats_cant eligible_cant voters_cant invalid_cant empty_cant valid_cant
rename seats_cant seats_cant1
rename eligible_cant eligible_cant1 
rename voters_cant voters_cant1
rename invalid_cant invalid_cant1
rename empty_cant empty_cant1 
rename valid_cant valid_cant1
duplicates drop canton year, force
sort canton year
save "$path\02_Processed_data\05_Alliances\Listenverbindungen_nomatch_temp.dta", replace
restore

// insert missing cantonal voting data to candidate data
sort canton year  
merge m:1 canton year using "$path\02_Processed_data\05_Alliances\Listenverbindungen_nomatch_temp.dta", ///
	gen(merge_cantinfo)

foreach var in seats_cant eligible_cant voters_cant invalid_cant empty_cant valid_cant {
	replace `var' = `var'1 if `var' == .
	drop `var'1
}	
drop merge_cantinfo
erase "$path\02_Processed_data\05_Alliances\nomatch_temp.dta"
erase "$path\02_Processed_data\05_Alliances\Listenverbindungen_nomatch_temp.dta"

* Manual corrections for OW & NW error       
replace partyname="Arnold Durrer (Freisinnig-Demokratische Partei)" if ID=="OW-1971-0001" & year==1971           
replace partyname="Vereinzelte" if ID=="NW-1955-0001" & year==1967
replace list="ARNOLD DURRER (FREISINNIG-DEMOKRATISCHE PARTEI)" if ID=="OW-1971-0001" & year==1971
replace list="VEREINZELTE" if ID=="NW-1955-0001" & year==1967

replace pvotes=votes if pvotes==. // We replace party votes for all candidates in one-seat cantons (see Email by Mark 30.1.2019)


****************************************************************
* (I) READ-IN STATUS CHECK BY ANITA AND GABRIELA
****************************************************************

preserve
import exc "$path\01_Raw_data\04_Candidates_status\StatusCheck_AnitaIn.xlsx", ///
	first clear
drop if ID==""
keep ID cand_1928
save "$path\02_Processed_data\04_Candidates_status\StatusCheck_AnitaIn.dta", ///
	replace 

import exc "$path\01_Raw_data\04_Candidates_status\StatusCheck_GabrielaIn.xlsx", ///
	first clear
keep ID cand_1925
drop if ID==""
merge 1:1 ID using ///
	"$path\02_Processed_data\04_Candidates_status\StatusCheck_AnitaIn.dta", ///
	gen(merge_sc)

gen cand_before1931=0
replace cand_before1931=1 if cand_1925==1 | cand_1928==1

keep ID cand_before1931		
save "$path\02_Processed_data\04_Candidates_status\StatusCheck.dta", replace 
restore

merge m:1 ID ///
	using "$path\02_Processed_data\04_Candidates_status\StatusCheck.dta", ///
	gen(merge_status)
drop if merge_status==2 
* Note: These are observations for which we merged IDs after the status check
*		(point (F)): AG-1963-0056, GE-1943-0006, GE-1959-9036, GE-1963-9216, 
*		and SO-1943-0033

replace cand_before1931=0 if cand_before1931==.
drop merge_status

erase "$path\02_Processed_data\04_Candidates_status\StatusCheck_AnitaIn.dta"
erase "$path\02_Processed_data\04_Candidates_status\StatusCheck.dta"


****************************************************************
* (J) OUTFILE FÜR BISNODE, ANITA, TODESDATEN, STATUS CHECK 
****************************************************************
	
* (i) Create file with all variants of name, firstname, birthyear, origin, and
* municipality



local vars "name firstname birthyear job partyname"
* Note: Create a file with all different name/firstname spellings and birthyears
* which we then append to original dataset.

foreach var in `vars'{	
preserve
bysort ID `var': gen indi_`var'=_n		
keep if indi_`var'==1
* Note: Keep only one observation per ID-var group.

bysort ID: gen observation_number =_n
keep ID `var' observation_number
reshape wide `var', i(ID) j(observation_number) 
save "$path\02_Processed_data\01_Elections_1931_1975\\`var'", replace
restore
}

preserve
reshape long origin originno,i(ID canton year) j(number)

drop if origin==""		
bysort ID origin originno: gen indi=_n
keep if indi==1
bysort ID: gen observation_number =_n
keep ID origin* originno* observation_number
reshape wide origin originno,i(ID) j(observation_number)		

save "$path\02_Processed_data\01_Elections_1931_1975\origin.dta", replace
restore
		
preserve

bysort ID municipality municipalityno: gen indi=_n
keep if indi==1
bysort ID: gen observation_number =_n
keep ID municipality* municipalityno* observation_number
reshape wide municipality municipalityno,i(ID) j(observation_number)		

save "$path\02_Processed_data\01_Elections_1931_1975\municipality.dta", replace
restore

preserve

drop name firstname birthyear municipality municipalityno origin* originno*
* Note: Drop original timevariant variables.

merge m:1 ID using "$path\02_Processed_data\01_Elections_1931_1975\name.dta", ///
	gen(merge_firstname)
merge m:1 ID using "$path\02_Processed_data\01_Elections_1931_1975\firstname.dta", ///
	gen(merge_name)
merge m:1 ID using "$path\02_Processed_data\01_Elections_1931_1975\birthyear.dta", ///
	gen(merge_birthyear)
merge m:1 ID using "$path\02_Processed_data\01_Elections_1931_1975\job.dta", ///
	gen(merge_job)
merge m:1 ID using "$path\02_Processed_data\01_Elections_1931_1975\partyname.dta", ///
	gen(merge_partyname)
merge m:1 ID using "$path\02_Processed_data\01_Elections_1931_1975\origin.dta", ///
	gen(merge_origin)
merge m:1 ID ///
	using "$path\02_Processed_data\01_Elections_1931_1975\municipality.dta", ///
	gen(merge_municipality)

erase "$path\02_Processed_data\01_Elections_1931_1975\name.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\firstname.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\birthyear.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\job.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\partyname.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\origin.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\municipality.dta"

bysort ID: gen indi=_n
* Note: Keep one observation per ID.
keep if indi==1

drop merge* indi*

*export exc "$path\07_Archive\11_Directors_1994_2018\Bisnode_out_wide.xlsx", ///
*	replace first(varlabels)

drop originno* municipalityno*
keep ID canton sex name1-name3 firstname1-firstname3 birthyear1 birthyear2 ///
	job1-job10 partyname1-partyname5 origin1-origin8 municipality1-municipality5 
gen GeburtsdatumNeu=.
gen TodesdatumNeu=.
gen WohnorteNeu=.
gen Bemerkungen=.
*export exc ///
*	"$path\07_Archive\09_DoB_and_DoD\JarckSpuehlerWidmer_out_wide.xlsx", ///
*	replace first(varlabels)

keep if birthyear1<=1912 | birthyear2<=1912

sort canton name1 firstname1

drop name3 job9 job10 origin5-origin8 municipality4 municipality5 ///
	GeburtsdatumNeu TodesdatumNeu WohnorteNeu Bemerkungen

gen cand_1925=""
gen cand_1928=""
		
gen name=name1
replace name=name+", "+name2 if name2!=""

gen firstname=firstname1
replace firstname=firstname+", "+firstname2 if firstname2!=""

tostring birthyear1 birthyear2, replace

gen birthyear=birthyear1
replace birthyear=birthyear+", "+birthyear2 if birthyear2!="."

gen job=job1
replace job=job+", "+job2 if job2!=""
replace job=job+", "+job3 if job3!=""
replace job=job+", "+job4 if job4!=""
replace job=job+", "+job5 if job5!=""
replace job=job+", "+job6 if job6!=""
replace job=job+", "+job7 if job7!=""
replace job=job+", "+job8 if job8!=""

gen partyname=partyname1
replace partyname=partyname+", "+partyname2 if partyname2!=""
replace partyname=partyname+", "+partyname3 if partyname3!=""
replace partyname=partyname+", "+partyname4 if partyname4!=""
replace partyname=partyname+", "+partyname5 if partyname5!=""

gen origin=origin1
replace origin=origin+", "+origin2 if origin2!=""
replace origin=origin+", "+origin3 if origin3!=""
replace origin=origin+", "+origin4 if origin4!=""

gen municipality=municipality1								
replace municipality=municipality+", "+municipality2 if municipality2!=""
replace municipality=municipality+", "+municipality3 if municipality3!=""

merge m:1 canton using "$path\01_Raw_data\07_Cantons\cantonal_translator.dta", ///
	gen(merge1)
	
tab canton if merge1==2

sort canton_no name firstname 
*drop name1-municipality3 cantonno merge1

*export exc "$path\07_Archive\09_DoB_and_DoDStatus\Status_Check_out_wide.xlsx", ///
*	first(varlabels) replace

restore	

* (ii) Reshape everything into long format

preserve

use "$path\02_Processed_data\11_Directors_1994_2018\bisnode_out_wide.dta", clear

local vars_done ""

local vars "name firstname birthyear"
foreach var in `vars'{	
reshape long `var', i(ID "`vars_done'") j(number_`var')
capture drop if `var'==.
capture drop if `var'==""
local vars_done `vars_done' "`var'"
}

local vars_done ""

local vars "origin municipality"
foreach var in `vars'{	
reshape long `var' `var'no, i(ID name firstname birthyear "`vars_done'") ///
	j(number_`var')
capture drop if `var'==.
capture drop if `var'==""
local vars_done `vars_done' "`var' `var'no"
}

save "$path\02_Processed_data\11_Directors_1994_2018\bisnode_out_long.dta", ///
	replace
*export exc "$path\02_Processed_data\11_Directors_1994_2018\Bisnode_out_long.xlsx", ///
*	first(varlabels) replace

restore

****************************************************************
* (K) READ-IN NACHRÜCKER
****************************************************************

* Mail von Madeleine Schneider, 15.10.2018 
* Guten Tag Herr Schmid
* Die Variable «Status» zeigt, ob eine Kandidatin/ein Kandidat in den National-
* rat gewählt wurde bzw. bereits früher Einsitz im Rat hatte.
* N=neu gewählt
* B=bisher
* E=ehemalig
* A=andere, d.h. noch nie gewählt

preserve
import exc "$path\01_Raw_data\04_Candidates_status\Ratsmitglieder_1848_DE.xlsx", ///
	first clear
keep if CouncilName=="Nationalrat" | (FirstName=="Edouard" & LastName=="Debétaz" & DateOfBirth=="16.09.1917")
replace DateLeaving="30.05.2019" if DateLeaving==""
* Note: Replace current date as leaving date for active politicians.

* Corrections after comparison with BfS status variable
replace DateJoining="01.12.1975" if FirstName=="Mario" & LastName=="Soldini" ///
	& DateOfBirth=="02.02.1913"
replace DateJoining="29.11.1971" if FirstName=="Franz" & LastName=="Jaeger" ///
	& DateOfBirth=="04.12.1941"
replace DateJoining="07.12.1959" if FirstName=="Roger" & LastName=="Dafflon" ///
	& DateOfBirth=="02.12.1914"	
expand 2 if FirstName=="Rolf" & LastName=="Seiler" & DateOfBirth=="12.03.1932", ///
	gen(obs_dup)
replace DateLeaving="25.11.1979" if FirstName=="Rolf" & LastName=="Seiler" & ///
	DateOfBirth=="12.03.1932"& DateOfBirth=="12.03.1932"	& obs_dup==0
replace DateJoining="28.11.1983" if FirstName=="Rolf" & LastName=="Seiler" & ///
	DateOfBirth=="12.03.1932"& DateOfBirth=="12.03.1932"	& obs_dup==1
drop obs_dup

expand 2 if FirstName=="Felix" & LastName=="Moeschlin" & DateOfBirth=="31.07.1882" , ///
	gen(obs_dup)
replace CantonAbbreviation="BS" if FirstName=="Felix" & LastName=="Moeschlin" & ///
	DateOfBirth=="31.07.1882"	& obs_dup==1
drop obs_dup
expand 2 if FirstName=="Edouard" & LastName=="Debétaz" & DateOfBirth=="16.09.1917" , ///
	gen(obs_dup)
replace DateJoining="02.12.1956" if  FirstName=="Edouard" & LastName=="Debétaz" &  ///
	DateOfBirth=="16.09.1917" & obs_dup==1
replace DateLeaving="30.11.1975" if FirstName=="Edouard" & LastName=="Debétaz" & ///
	DateOfBirth=="16.09.1917" & obs_dup==1
drop if  FirstName=="Edouard" & LastName=="Debétaz" & ///
	DateOfBirth=="16.09.1917" & obs_dup==0 // drop Ständerat spell of this person
replace CouncilName="Nationalrat" if FirstName=="Edouard" & LastName=="Debétaz" & ///
	DateOfBirth=="16.09.1917"
drop obs_dup	
replace DateJoining="02.12.1935" if FirstName=="Heinrich" & LastName=="Schnyder" & ///
	DateOfBirth=="08.08.1897" 

replace DateJoining="04.12.1939" if FirstName=="Felix" & LastName=="Moeschlin" & ///
	DateOfBirth=="31.07.1882"	
replace DateJoining="12.09.1984" if FirstName=="Armand" & LastName=="Magnin" & ///
	DateOfBirth=="07.02.1920" 	
	

* (i) check deviations in firstname and lastname 

egen FirstName_num=group(FirstName)
egen LastName_num=group(LastName)
		
bysort LastName GenderAsString DateOfBirth : egen FirstName_sd=sd(FirstName_num)
bysort FirstName GenderAsString DateOfBirth : egen LastName_sd=sd(LastName_num)

list FirstName LastName DateOfBirth if FirstName_sd>0 & FirstName_sd!=.
list FirstName LastName DateOfBirth if LastName_sd>0 & LastName_sd!=.

replace FirstName="Eduard" if FirstName=="Edi" & LastName=="Engelberger" ///
	& DateOfBirth=="26.01.1940"
replace FirstName="Gregor A." if FirstName=="Gregor" & LastName=="Rutz" ///
	& DateOfBirth=="12.10.1972"
replace FirstName="Victor Emil" if FirstName=="Victor" & LastName=="Scherer" ///
	& DateOfBirth=="03.10.1881"
replace FirstName="Silva" if FirstName=="Silva Anita" & LastName=="Semadeni" ///
	& DateOfBirth=="08.02.1952"
replace LastName="Haering" if FirstName=="Barbara" & LastName=="Haering Binder" ///
	& DateOfBirth=="20.09.1953"
replace LastName="Markwalder" if FirstName=="Christa" ///
	& LastName=="Markwalder Bär" & DateOfBirth=="27.07.1975"
replace LastName="Leuthard" if FirstName=="Doris" & LastName=="Leuthard Hausin" ///
	& DateOfBirth=="10.04.1963"
replace LastName="Hutter" if FirstName=="Jasmin" & LastName=="Hutter-Hutter" ///
	& DateOfBirth=="11.06.1978"
replace LastName="Schneider-Ammann" if FirstName=="Johann N." ///
	& LastName=="Schneider" & DateOfBirth=="18.02.1952"
replace LastName="Dormond" if FirstName=="Marlyse" ///
	& LastName=="Dormond Béguelin" & DateOfBirth=="31.10.1949"
replace LastName="Haller" if FirstName=="Ursula" & LastName=="Haller Vannini" ///
	& DateOfBirth=="04.11.1948"

* (ii) Check possible cases in which both first and last name differ per indi-
* vidual

drop FirstName_num LastName_num
egen FirstName_num=group(FirstName)
egen LastName_num=group(LastName)
		
bysort CantonAbbreviation DateOfBirth : egen FirstName_sd2=sd(FirstName_num)
bysort CantonAbbreviation DateOfBirth: egen LastName_sd2=sd(LastName_num)
sort CantonAbbreviation DateOfBirth
list CantonAbbreviation FirstName LastName DateOfBirth if FirstName_sd2>0 ///
	& FirstName_sd2!=.
list CantonAbbreviation FirstName LastName DateOfBirth if LastName_sd2>0 ///
	& LastName_sd2!=.
* Result: All these cases are different people who are running in the same can-
* ton and have the same birthday	
	
drop FirstName_sd LastName_sd FirstName_num LastName_num FirstName_sd2 ///
	LastName_sd2

* (iii) Eliminate duplicates

duplicates tag, g(surplus)
list FirstName LastName CantonAbbreviation PartyAbbreviation DateJoining ///
	DateLeaving surplus if surplus>0
duplicates drop Active - DateOfDeath, force

* (iv) Calculate date of joining per individual

gen EDateJoining = date(DateJoining, "DMY")
gen EDateLeaving = date(DateLeaving, "DMY")

replace CantonAbbreviation="BEJU" if CantonAbbreviation=="JU" | ///
	CantonAbbreviation=="BE"

egen id_nachruecker=group(CantonAbbreviation FirstName LastName DateOfBirth)
* Note: Install dataset as a panel.
sort id_nachruecker EDateJoining 
bysort id_nachruecker: gen id_time=_n
xtset id_nachruecker id_time
	
local vars "Joining Leaving" 
foreach var in `vars'{			
gen EYear`var' = year(EDate`var')
gen EMonth`var' = month(EDate`var')		
bysort id_nachruecker: egen EDate`var'_Min=min(EDate`var')
bysort id_nachruecker: egen EDate`var'_Max=max(EDate`var')
bysort id_nachruecker: egen EYear`var'_Min=min(EYear`var')
bysort id_nachruecker: egen EYear`var'_Max=max(EYear`var')
gen EDate`var'_F1=F1.EDate`var'
gen EDate`var'_L1=L1.EDate`var'
}

drop if EYearLeaving_Max<=1931
* Note: Keep only relevant information

gen EDate_Diff=EDateJoining-EDateLeaving_L1
bysort id_nachruecker: egen EDate_Diff_Max=max(EDate_Diff)


* a) individuals with one single spell

gen EDateJoining1=EDateJoining_Min if EDate_Diff_Max==1 | EDate_Diff_Max==.
gen EDateLeaving1=EDateLeaving_Max if EDate_Diff_Max==1 | EDate_Diff_Max==.
	

* b) individuals with more than one spell

sort id_nachruecker id_time
*br FirstName LastName DateOfBirth id* DateJoining DateLeaving EDate_Diff ///
*	EDateJoining1 EDateLeaving1

forvalues i = 2(1)10 {
generate EDateJoining`i' = .
generate EDateLeaving`i' = .
}

bysort id_nachruecker: egen id_time_max=max(id_time)
levelsof id_nachruecker, local(id_levels), if id_time_max>1
* Note: Loop over all individuals with more than one spell.
foreach l of local id_levels {
local x=1
levelsof id_time, local(time_levels), if id_nachruecker == `l'
foreach t of local time_levels {
	sum EDate_Diff if id_nachruecker == `l' & id_time==`t'
	local Spell_Diff=r(max)
	if `t'==1{ // first observation
			replace EDateJoining1=EDateJoining if id_nachruecker == `l' ///
			& id_time==1
			replace EDateLeaving1=EDateLeaving if id_nachruecker == `l' ///
			& id_time==1
			}
	else if `Spell_Diff'==1 { // no breaks between spells
			sum EDateLeaving if id_nachruecker == `l' & id_time==`t'
			replace EDateLeaving`x'=r(max) if id_nachruecker == `l' & id_time==1
			}
	else if `Spell_Diff'>1 & `Spell_Diff'!=.{ // breaks between spells 
			local x=`x'+1
			sum EDateLeaving if id_nachruecker == `l' & id_time==`t'
			replace EDateLeaving`x'=r(max) if id_nachruecker == `l' & id_time==1				
			sum EDateJoining if id_nachruecker == `l' & id_time==`t'
			replace EDateJoining`x'=r(max) if id_nachruecker == `l' & id_time==1
			}
}
}

* (v) Check outfile
sort id_nachruecker id_time
*br FirstName LastName DateOfBirth id* DateJoining DateLeaving EDate_Diff ///
*	EDateJoining EDateLeaving EDateJoining1 EDateLeaving1 EDateJoining2 ///
*	EDateLeaving2 EDateJoining3 EDateLeaving3 EDateJoining4 EDateLeaving4 

* (vi) save outfile

gen birthyear=substr(DateOfBirth,-4,.)
destring birthyear, replace

rename LastName name 
rename FirstName firstname
rename GenderAsString sex
rename CantonAbbreviation canton_tomerge
replace sex="F" if sex=="f"
replace sex="M" if sex=="m"

replace name="ab Yberg" if name=="Ab Yberg" & firstname=="Alois" & canton_tomerge=="SZ" ///
	& sex=="M" & birthyear==1878
replace firstname="Roman" if name=="Abt" & firstname=="H. Roman" & canton_tomerge=="AG" ///
	& sex=="M" & birthyear==1883
replace name="Aeppli" if name=="Aeppli Wartmann" & firstname=="Regine" ///
	& canton_tomerge=="ZH" & sex=="F" & birthyear==1952
replace firstname="Rudolf" if name=="Aeschbacher" & firstname=="Ruedi" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1941
replace name="Amacker" if name=="Amacker-Amann" & firstname=="Kathrin" ///
	& canton_tomerge=="BL" & sex=="F" & birthyear==1962
replace name="Baumann-Bieri" if name=="Baumann" & firstname=="Stephanie" ///
	& canton_tomerge=="BEJU" & sex=="F" & birthyear==1951
replace firstname="Hans-Ulrich" if name=="Baumberger" & firstname=="Hans Ulrich" ///
	& canton_tomerge=="AR" & sex=="M" & birthyear==1932
replace firstname="Pierre" if name=="Benninger" & firstname=="Peter" ///
	& canton_tomerge=="FR" & sex=="M" & birthyear==1879
replace firstname="Ueli" if name=="Blatter" & firstname=="Ulrich" ///
	& canton_tomerge=="OW" & sex=="M" & birthyear==1940
replace name="Börlin" if name=="Boerlin" & firstname=="Ernst" & canton_tomerge=="BL" ///
	& sex=="M" & birthyear==1905
replace firstname="Francesco" if name=="Borella" & firstname=="Francesco-Nino" ///
	& canton_tomerge=="TI" & sex=="M" & birthyear==1883
replace name="Böschung" if name=="Boschung" & firstname=="Franz" & canton_tomerge=="FR" ///
	& sex=="M" & birthyear==1868
replace name="Bremi" if name=="Bremi-Forrer" & firstname=="Ulrich" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1929
replace firstname="Walter" if name=="Bringolf" & firstname=="Walther" ///
	& canton_tomerge=="SH" & sex=="M" & birthyear==1895
replace firstname="Albert" if name=="Broger" ///
	& firstname=="Johann Baptist Albert" & canton_tomerge=="AI" & sex=="M" ///
		& birthyear==1897
replace name="Bruegger" if name=="Brügger" & firstname=="Cyrill" ///
	& canton_tomerge=="FR" & sex=="M" & birthyear==1938
replace firstname="Martin H." if name=="Burckhardt" & firstname=="Martin" ///
	& canton_tomerge=="BS" & sex=="M" & birthyear==1921
replace firstname="Henri" if name=="Burrus" & firstname=="Henry" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1882
replace firstname="Alfred J." if name=="Büchi" & firstname=="Alfred" ///
	& canton_tomerge=="ZH" & sex=="M"& birthyear==1879
replace firstname="Robert Eduard" if name=="Bühler" & firstname=="Robert" ///
	& canton_tomerge=="ZH" & sex=="M"& birthyear==1902
replace firstname="Rolf" if name=="Bühler" & firstname=="Rolf Theodor" ///
	& canton_tomerge=="SG" & sex=="M"& birthyear==1903
replace firstname="Alexandre" if name=="Cailler" ///
	& firstname=="Alexandre-François-Louis" & canton_tomerge=="FR" & sex=="M" ///
	& birthyear==1866
replace firstname="Germain" if name=="Carnat" & firstname=="Germain-Joseph" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1894
replace firstname="Andrea Claudio" if name=="Caroni" & firstname=="Andrea" ///
	& canton_tomerge=="AR" & sex=="M" & birthyear==1980
replace firstname="Francesco" if name=="Cavalli" & firstname=="Franco" ///
	& canton_tomerge=="TI" & sex=="M" & birthyear==1942
replace firstname="Charles" if name=="Dellberg" & firstname=="Karl" ///
	& canton_tomerge=="VS" & sex=="M" & birthyear==1886
replace name="Diethelm-Dobler" if name=="Diethelm" & firstname=="Josef" ///
	& canton_tomerge=="SZ" & sex=="M" & birthyear==1914
replace firstname="Emil Otto" if name=="Duft" & firstname=="Emil" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1895
replace firstname="Kurt Rudolf" if name=="Düby" & firstname=="Kurt RU" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1900
replace firstname="Anton" if name=="Eberhard" & firstname=="Toni" ///
	& canton_tomerge=="SZ" & sex=="M" & birthyear==1949
replace firstname="Mathias" if name=="Eggenberger" & firstname=="Matthias" ///
	& canton_tomerge=="SG" & sex=="M" & birthyear==1905
replace firstname="Franz Mathäus" if name=="Egger" & firstname=="Franz" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1899
replace name="Egger" if name=="Egger-Wyss" & firstname=="Esther" & canton_tomerge=="AG" ///
	& sex=="F" & birthyear==1952
replace name="Eichenberger" if name=="Eichenberger-Walther" ///
	& firstname=="Corina" & canton_tomerge=="AG" & sex=="F" & birthyear==1954
replace firstname="Edi" if name=="Engelberger" & firstname=="Eduard" ///
	& canton_tomerge=="NW" & sex=="M" & birthyear==1940
replace firstname="Joseph" if name=="Escher" & firstname=="Josef" ///
	& canton_tomerge=="VS" & sex=="M" & birthyear==1885
replace name="Eugster-Züst" if name=="Eugster" & firstname=="Howard" ///
	& canton_tomerge=="AR" & sex=="M" & birthyear==1861
replace firstname="Fritz" if name=="Eymann" & firstname=="Fritz Henry" ///
	& canton_tomerge=="NE" & sex=="M" & birthyear==1880
replace firstname="David Hirsch" if name=="Farbstein" & firstname=="David" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1868
replace firstname="Charles" if name=="Favrod-Coune" ///
	& firstname=="Charles-Albert" & canton_tomerge=="VD" & sex=="M" & birthyear==1877
replace name="Flückiger" if name=="Flückiger-Bäni" & firstname=="Sylvia" ///
	& canton_tomerge=="AG" & sex=="F" & birthyear==1952
replace firstname="Armand" if name=="Forel" & firstname=="Armand-Auguste" ///
	& canton_tomerge=="VD" & sex=="M" & birthyear==1920
replace birthyear=1878 if name=="Frank" & firstname=="Ferdinand" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1873
replace name="Fässler" if name=="Fässler-Osterwalder" & firstname=="Hildegard" ///
	& canton_tomerge=="SG" & sex=="F" & birthyear==1951
replace name="Füeg-Hitz" if name=="Füeg" & firstname=="Cornelia" & canton_tomerge=="SO" ///
	& sex=="F" & birthyear==1941
replace firstname="Brigitta Maria" if name=="Gadient" & firstname=="Brigitta M." ///
	& canton_tomerge=="GR" & sex=="F" & birthyear==1960
replace firstname="Chantal Juliane" if name=="Galladé" & firstname=="Chantal" ///
	& canton_tomerge=="ZH" & sex=="F" & birthyear==1972
replace firstname="Remo" if name=="Galli" & firstname=="Remo Giosué" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1943
replace firstname="Raymund" if name=="Gamma" & firstname=="Raymond" ///
	& canton_tomerge=="UR" & sex=="M" & birthyear==1919
replace firstname="August" if name=="Gattiker" & firstname=="Oskar August" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1873
replace firstname="Andrea" if name=="Geissbühler" & firstname=="Andrea Martina" ///
	& canton_tomerge=="BEJU" & sex=="F" & birthyear==1976
replace firstname="Rudolf" if name=="Gelpke" & firstname=="Rudolf Arnold" ///
	& canton_tomerge=="BS" & sex=="M" & birthyear==1873
replace name="Geser-Rhoner" if name=="Geser" & firstname=="Albert" ///
	& canton_tomerge=="SG" & sex=="M" & birthyear==1868
replace name="Girard-Montet" if name=="Girard" & firstname=="Gertrude" ///
	& canton_tomerge=="VD" & sex=="F" & birthyear==1913
replace name="Glauser" if name=="Glauser-Zufferey" & firstname=="Alice" ///
	& canton_tomerge=="VD" & sex=="F" & birthyear==1954
replace name="Gmür" if name=="Gmür-Schönenberger" & firstname=="Andrea" ///
	& canton_tomerge=="LU" & sex=="F" & birthyear==1964
replace firstname="Emil Johann" if name=="Graf" & firstname=="Emil" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1890
replace firstname="Otto" if name=="Graf" & firstname=="Ernst Otto" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1877
replace firstname="Hans Ulrich" if name=="Graf" & firstname=="Hans-Ulrich" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1922
replace firstname="Achille" if name=="Grospierre" & firstname=="Achille-Tell" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1872
replace firstname="Niklaus" if name=="Gugger" & firstname=="Niklaus-Samuel" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1970
replace firstname="Wilfried" if name=="Gusset" & firstname=="Wilfried Ernest" ///
	& canton_tomerge=="TG" & sex=="M" & birthyear==1946
replace name="Gyr" if name=="Gyr-Steiner" & firstname=="Josy" & canton_tomerge=="SZ" ///
	& sex=="F" & birthyear==1949
replace firstname="Hans" if name=="Gysin" & firstname=="Hans Rudolf" ///
	& canton_tomerge=="BL" & sex=="M" & birthyear==1940
replace firstname="Ueli" if name=="Götsch" & firstname=="Ulrich" & canton_tomerge=="ZH" ///
	& sex=="M" & birthyear==1925
replace birthyear=1874 if name=="Hartmann" & firstname=="Georg" & canton_tomerge=="GR" ///
	& sex=="M" & birthyear==1873
replace firstname="Fritz" if name=="Hauser" & firstname=="Fritz E." ///
	& canton_tomerge=="BS" & sex=="M" & birthyear==1884
replace firstname="Erich J." if name=="Hess" & firstname=="Erich" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1981
replace firstname="Otto" if name=="Hess" & firstname=="Otto (senior)" ///
	& canton_tomerge=="TG" & sex=="M" & birthyear==1897
replace firstname="Franz" if name=="Hildbrand" & firstname=="Franz-Joseph" ///
	& canton_tomerge=="VS" & sex=="M" & birthyear==1942
replace firstname="Walter" if name=="Hofer" & firstname=="Walther" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1920
replace firstname="René" if name=="Houriet" & firstname=="René-Albert" ///
	& canton_tomerge=="VD" & sex=="M" & birthyear==1920
replace name="Humbel Näf" if name=="Humbel" & firstname=="Ruth" & canton_tomerge=="AG" ///
	& sex=="F" & birthyear==1957
replace firstname="Jules Frédéric" if name=="Humbert-Droz" & firstname=="Jules" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1891
replace name="Häberli" if name=="Häberli-Koller" & firstname=="Brigitte" ///
	& canton_tomerge=="TG" & sex=="F" & birthyear==1958
replace name="Hoesli" if name=="Hösli" & firstname=="Fritz" & canton_tomerge=="GL" ///
	& sex=="M" & birthyear==1922
replace name="Jacottet" if name=="Jaccottet" & firstname=="Henri-Louis" ///
	& canton_tomerge=="VD" & sex=="M" & birthyear==1881
replace firstname="Henri" if name=="Jacottet" & firstname=="Henri-Louis" ///
	& canton_tomerge=="VD" & sex=="M" & birthyear==1881
replace birthyear=1906 if name=="Jacquod" & firstname=="René" & canton_tomerge=="VS" ///
	& sex=="M" & birthyear==1905
replace name="Jäger" if name=="Jaeger" & firstname=="Walter" & canton_tomerge=="BS" ///
	& sex=="M" & birthyear==1911
replace name="Jaquet" if name=="Jaquet-Berger" & firstname=="Christiane" ///
	& canton_tomerge=="VD" & sex=="F" & birthyear==1937
replace firstname="Matthias" if name=="Jauslin" & firstname=="Matthias Samuel" ///
	& canton_tomerge=="AG" & sex=="M" & birthyear==1962
replace name="Jenny-Schuler" if name=="Jenny" & firstname=="Heinrich" ///
	& canton_tomerge=="GL" & sex=="M" & birthyear==1861
replace name="Jossen" if name=="Jossen-Zinsstag" & firstname=="Peter" ///
	& canton_tomerge=="VS" & sex=="M" & birthyear==1955
replace name="Junod-Leder" if name=="Junod" & firstname=="Georges" ///
	& canton_tomerge=="VD" & sex=="M" & birthyear==1880
replace name="Kleiner-Schläpfer" if name=="Kleiner" & firstname=="Marianne" ///
	& canton_tomerge=="AR" & sex=="F" & birthyear==1947
replace firstname="Marcel" if name=="Krügel" & firstname=="Marcel-René" ///
	& canton_tomerge=="NE" & sex=="M" & birthyear==1893
replace firstname="Franz Josef" if name=="Kurmann" & firstname=="Franz-Josef" ///
	& canton_tomerge=="LU" & sex=="M" & birthyear==1917
replace firstname="Jakob" if name=="Kägi" & firstname=="E. Jakob" & canton_tomerge=="ZH" ///
	& sex=="M" & birthyear==1886
replace name="Kaempfen" if name=="Kämpfen" & firstname=="Moritz" & canton_tomerge=="VS" ///
	& sex=="M" & birthyear==1907
replace firstname="Hans Paul" if name=="Künzi" & firstname=="Hans" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1924
replace name="Landolt-Rast" if name=="Landolt" & firstname=="Franz" ///
	& canton_tomerge=="GL" & sex=="M" & birthyear==1901
replace firstname="Erwin A." if name=="Lang" & firstname=="Erwin" & canton_tomerge=="ZH" ///
	& sex=="M" & birthyear==1908
replace name="Lang" if name=="Lang-Gehri" & firstname=="Hedi" & canton_tomerge=="ZH" ///
	& sex=="F" & birthyear==1931
replace firstname="Hajo" if name=="Leutenegger" & firstname=="Hansjakob" ///
	& canton_tomerge=="ZG" & sex=="M" & birthyear==1944
replace name="Leutenegger" if name=="Leutenegger Oberholzer" ///
	& firstname=="Susanne" & canton_tomerge=="BL" & sex=="F" & birthyear==1948
replace firstname="Hans Georg" if name=="Lüchinger" & firstname=="Hans-Georg" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1927
replace name="Maître" if name=="Maitre" & firstname=="Yves" & canton_tomerge=="GE" ///
	& sex=="M" & birthyear==1917
replace firstname="Fritz" if name=="Marbach" & firstname=="Friedrich" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1892
replace name="Marty Kaelin" if name=="Marty Kälin" & firstname=="Barbara" ///
	& canton_tomerge=="ZH" & sex=="F" & birthyear==1954
replace name="Mascarin" if name=="Mascarin-Bircher" & firstname=="Ruth" ///
	& canton_tomerge=="BS" & sex=="F" & birthyear==1945
replace firstname="Paul Ulrich" if name=="Meierhans" & firstname=="Paul" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1895
replace name="Menétrey" if name=="Menétrey-Savary" ///
	& firstname=="Anne-Catherine" & canton_tomerge=="VD" & sex=="F" & birthyear==1938
replace firstname="John" if name=="Mermod" & firstname=="John-Edouard" ///
	& canton_tomerge=="VD" & sex=="M" & birthyear==1882
replace firstname="Rodolphe" if name=="Metry" & firstname=="Rudolf" ///
	& canton_tomerge=="VS" & sex=="M" & birthyear==1885
replace firstname="Ludwig Friedrich" if name=="Meyer" ///
	& firstname=="Ludwig-Friedrich" & canton_tomerge=="LU" & sex=="M" & birthyear==1872
replace firstname="Mattea Julia" if name=="Meyer" ///
	& firstname=="Mattea" & canton_tomerge=="ZH" & sex=="F" & birthyear==1987
replace firstname="Carl" if name=="Miville" & firstname=="Carl jun." ///
	& canton_tomerge=="BS" & sex=="M" & birthyear==1921
replace firstname="Gion Rudolf" if name=="Mohr" & firstname=="Gian Rudolf" ///
	& canton_tomerge=="GR" & sex=="M" & birthyear==1885
replace name="von Moos" if name=="Moos" & firstname=="Herbert" & canton_tomerge=="BEJU" ///
	& sex=="M" & birthyear==1893
replace firstname="Hans" if name=="Müller" & firstname=="Hans Gottfried" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1893
replace firstname="Heinrich" if name=="Müller" & firstname=="Heinrich C." ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1911
replace firstname="Johann" if name=="Müller" & firstname=="Johannes" ///
	& canton_tomerge=="SG" & sex=="M" & birthyear==1883
replace name="Müller" if name=="Müller-Altermatt" & firstname=="Stefan" ///
	& canton_tomerge=="SO" & sex=="M" & birthyear==1976
replace name="Muri" if name=="Müri" & firstname=="Hermann" & canton_tomerge=="AG" ///
	& sex=="M" & birthyear==1874
replace name="Nägeli" if name=="Naegeli" & firstname=="Wilfried" & canton_tomerge=="TG" ///
	& sex=="M" & birthyear==1932
replace firstname="Hans Rudolf" if name=="Nebiker" & firstname=="Hans-Rudolf" ///
	& canton_tomerge=="BL" & sex=="M" & birthyear==1929
replace firstname="Rudolf" if name=="Noser" & firstname=="Ruedi" & canton_tomerge=="ZH" ///
	& sex=="M" & birthyear==1961
replace birthyear=1892 if name=="Nüesch" & firstname=="Jakob" & canton_tomerge=="SG" ///
	& sex=="M" & birthyear==1893
replace firstname="Josef" if name=="Odermatt" & firstname=="Joseph" ///
	& canton_tomerge=="NW" & sex=="M" & birthyear==1892
replace firstname="Hans Christian" if name=="Oester" & firstname=="Hans" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1931
replace sex="F" if name=="Paccolat" & firstname=="Monique" & canton_tomerge=="VS" ///
	& sex=="M" & birthyear==1954
replace birthyear=1878 if name=="Pfister" & firstname=="Eduard" & canton_tomerge=="TG" ///
	& sex=="M" & birthyear==1873
replace birthyear=1901 if name=="Philippe" & firstname=="Etienne" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1902
replace firstname="Albert" if name=="Picot" & firstname=="Albert-Edouard" ///
	& canton_tomerge=="GE" & sex=="M" & birthyear==1882
replace name="Piller" if name=="Piller Carrard" & firstname=="Valérie" ///
	& canton_tomerge=="FR" & sex=="F" & birthyear==1978
replace firstname="Johannes Robert" if name=="Randegger" & firstname=="Johannes" ///
	& canton_tomerge=="BS" & sex=="M" & birthyear==1941
replace firstname="Paul" if name=="Randon" & firstname=="Georges-Paul" ///
	& canton_tomerge=="GE" & sex=="M" & birthyear==1901
replace firstname="Rudolf" if name=="Reichling" & firstname=="Rudolf jun." ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1924
replace firstname="Katharina" if name=="Riklin" & firstname=="Kathy" ///
	& canton_tomerge=="ZH" & sex=="F" & birthyear==1952
replace firstname="John" if name=="Rochaix" & firstname=="John-Marc" ///
	& canton_tomerge=="GE" & sex=="M" & birthyear==1879
replace firstname="Max" if name=="Rohr" & firstname=="Max Albert" & canton_tomerge=="AG" ///
	& sex=="M" & birthyear==1890
replace firstname="François" if name=="Rossiaud" & firstname=="François Joseph" ///
	& canton_tomerge=="GE" & sex=="M" & birthyear==1880
replace name="Roth Bernasconi" if name=="Roth-Bernasconi" & firstname=="Maria" ///
	& canton_tomerge=="GE" & sex=="F" & birthyear==1955
replace firstname="Rebecca" if name=="Ruiz" & firstname=="Rebecca Ana" ///
	& canton_tomerge=="VD" & sex=="F" & birthyear==1982
replace firstname="Giovan Battista" if name=="Rusca" ///
	& firstname=="Giovan-Battista" & canton_tomerge=="TI" & sex=="M" & birthyear==1881
replace firstname="Henry" if name=="Sandoz" & firstname=="Henri Auguste" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1878
replace name="Oppliger-Schenker" if name=="Schenker" & firstname=="Silvia" ///
	& canton_tomerge=="BS" & sex=="F" & birthyear==1954
replace firstname="Carl Eugen" if name=="Scherrer" & firstname=="Carl" ///
	& canton_tomerge=="SH" & sex=="M" & birthyear==1909
replace firstname="Josef" if name=="Scherrer-Brisig" & firstname=="Joseph Anton" ///
	& canton_tomerge=="SG" & sex=="M" & birthyear==1891
replace name="Scherrer" if name=="Scherrer-Brisig" & firstname=="Josef" ///
	& canton_tomerge=="SG" & sex=="M" & birthyear==1891
replace firstname="August" if name=="Schirmer" & firstname=="L. August" ///
	& canton_tomerge=="SG" & sex=="M" & birthyear==1881
replace firstname="Ueli" if name=="Schlüer" & firstname=="Ulrich" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1944
replace name="Schmid" if name=="Schmid-Käser" & firstname=="Rudolf" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1888
replace firstname="Johann" if name=="Schneider-Ammann" & firstname=="Johann N." ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1952
replace name="Schneider" if name=="Schneider-Schneiter" & firstname=="Elisabeth" ///
	& canton_tomerge=="BL" & sex=="F" & birthyear==1964
replace firstname="Ludwig Max" if name=="Schneller" & firstname=="Ludwig" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1879
replace firstname="Adelrich Jakob" if name=="Schuler" ///
	& firstname=="Adelrich Jacob" & canton_tomerge=="ZH" & sex=="M" & birthyear==1922
replace birthyear=1898 if name=="Schuler" & firstname=="Hans" & canton_tomerge=="GL" ///
	& sex=="M" & birthyear==1899
replace firstname="James Eduard" if name=="Schwarzenbach" & firstname=="James" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1911
replace firstname="Jakob" if name=="Schwendener" & firstname=="Johann Jakob" ///
	& canton_tomerge=="SG" & sex=="M" & birthyear==1888
replace name="Schaerli" if name=="Schärli" & firstname=="Hans" & canton_tomerge=="LU" ///
	& sex=="M" & birthyear==1925
replace firstname="Fritz" if name=="Segessenmann" & firstname=="Friedrich" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1897
replace name="Segmüller-Weber" if name=="Segmüller" & firstname=="Eva" ///
	& canton_tomerge=="SG" & sex=="F" & birthyear==1932
replace firstname="Adolf" if name=="Seiler" & firstname=="Gustav Adolf" ///
	& canton_tomerge=="BL" & sex=="M" & birthyear==1875
replace name="Seiler-Graf" if name=="Seiler Graf" & firstname=="Priska" ///
	& canton_tomerge=="ZH" & sex=="F" & birthyear==1968
replace firstname="Walter Paul" if name=="Siegmann" & firstname=="Walter" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1910
replace name="Sigerist-Schalch" if name=="Sigerist" & firstname=="Heinrich" ///
	& canton_tomerge=="SH" & sex=="M" & birthyear==1881
replace firstname="Hans Konrad" if name=="Sonderegger" ///
	& firstname=="Hans-Konrad" & canton_tomerge=="BL" & sex=="M" & birthyear==1891
replace firstname="Hannes" if name=="Steffen" & firstname=="Hans" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1931
replace firstname="Edouard" if name=="Steinmetz" & firstname=="Edouard Ch." ///
	& canton_tomerge=="GE" & sex=="M" & birthyear==1865
replace firstname="Ernst Alfred" if name=="Stiefel" & firstname=="Ernst" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1887
replace name="Stocker" if name=="Stocker-Meier" & firstname=="Monika" ///
	& canton_tomerge=="ZH" & sex=="F" & birthyear==1948
replace name="Streiff" if name=="Streiff-Feller" & firstname=="Marianne" ///
	& canton_tomerge=="BEJU" & sex=="F" & birthyear==1957
replace firstname="Carl" if name=="Sulzer" & firstname=="Carl J." ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1865
replace name="Sulzer-Schmid" if name=="Sulzer" & firstname=="Carl" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1865
replace firstname="Johannes" if name=="Surbeck" & firstname=="Johann" ///
	& canton_tomerge=="BL" & sex=="M" & birthyear==1892
replace name="Thorens" if name=="Thorens Goumaz" & firstname=="Adèle" ///
	& canton_tomerge=="VD" & sex=="F" & birthyear==1971
replace name="Travelletti" if name=="Traveletti" & firstname=="Adolphe" ///
	& canton_tomerge=="VS" & sex=="M" & birthyear==1914
replace firstname="Walter" if name=="Trüb" & firstname=="Walther" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1883
replace firstname="Klemens" if name=="Ulrich" & firstname=="Clemenz" ///
	& canton_tomerge=="SZ" & sex=="M" & birthyear==1875
replace firstname="Henry" if name=="Vallotton" ///
	& firstname=="Henri-François-Jules" & canton_tomerge=="VD" & sex=="M" & birthyear==1891
replace firstname="Georges" if name=="Wander" & firstname=="Georg" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1898
replace firstname="Rudolf" if name=="Weber" & firstname=="Jakob Rudolf" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1887
replace name="Welti-Preiswerk" if name=="Welti" & firstname=="Franz" ///
	& canton_tomerge=="BS" & sex=="M" & birthyear==1879
replace firstname="Otto" if name=="Wenger" & firstname=="Otto M." ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1910
replace firstname="Max" if name=="Wey" & firstname=="Max Sigmund" ///
	& canton_tomerge=="LU" & sex=="M" & birthyear==1892
replace firstname="Hans" if name=="Widmer" & firstname=="Hans-Gottlieb" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1889
replace name="Widmer" if name=="Widmer-Kunz" & firstname=="Walter" ///
	& canton_tomerge=="AG" & sex=="M" & birthyear==1897
replace firstname="Karl" if name=="Wunderli" & firstname=="Karl J." ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1881
replace firstname="Fritz" if name=="Wüthrich" & firstname=="Fritz-Charles" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1893
replace birthyear=1943 if name=="Zbinden" & firstname=="Hans" & canton_tomerge=="AG" ///
	& sex=="M" & birthyear==1945
replace firstname="Bruno J." if name=="Zuppiger" & firstname=="Bruno" ///
	& canton_tomerge=="ZH" & sex=="M" & birthyear==1952
replace firstname="Otto" if name=="Zwygart" & firstname=="Otto (senior)" ///
	& canton_tomerge=="BEJU" & sex=="M" & birthyear==1911
replace firstname="Guido A." if name=="Zäch" & firstname=="Guido" ///
	& canton_tomerge=="AG" & sex=="M" & birthyear==1935
replace name="Zölch-Balmer" if name=="Zölch" & firstname=="Elisabeth" ///
	& canton_tomerge=="BEJU" & sex=="F" & birthyear==1951
replace firstname="Hans" if name=="von Matt" & firstname=="Hans jun." ///
	& canton_tomerge=="NW" & sex=="M" & birthyear==1869
replace name="von Rotz-Spichtig" if name=="von Rotz" & firstname=="Christoph" ///
	& canton_tomerge=="OW" & sex=="M" & birthyear==1966
replace firstname="Alphons" if name=="von Streng" & firstname=="Alfons" ///
	& canton_tomerge=="TG" & sex=="M" & birthyear==1852
replace firstname="Carl" if name=="von Weber" & firstname=="Karl" ///
	& canton_tomerge=="SZ" & sex=="M" & birthyear==1879

keep if id_time==1
keep name firstname canton_tomerge sex birthyear EDateLeaving1-EDateLeaving6 ///
	EDateJoining1-EDateJoining6 PartyName DateOfBirth ///
	DateOfDeath
order name firstname canton_tomerge sex birthyear EDateJoining1 EDateLeaving1 ///
	EDateJoining2 EDateLeaving2 EDateJoining3 EDateLeaving3 ///
	EDateJoining4 EDateLeaving4 EDateJoining5 EDateLeaving5 ///
	EDateJoining6 EDateLeaving6 PartyName DateOfBirth ///
	DateOfDeath

save "$path\02_Processed_data\04_Candidates_status\nachruecker_spells.dta", replace

restore

* Check duplicate in terms of canton, firstname, name, and birthyear
* duplicates drop ID, force
* duplicates tag name firstname canton sex birthyear, gen(dups)
* br name firstname canton sex birthyear job status partyname municipality elected if dups>0

* drop if dups>0
* Note: duplicates that are not the same persons in election data in terms of 
* canton, firstname, name, and birthyear
*ID				year	canton	firstname	name		birthyear	job							partyname
*BEJU-1991-0018	1991	BE		Beat		Balmer		1958		Konstrukteur				FPS/PSL
*BEJU-2007-0026	2007	BE		Beat		Balmer		1958		Bankangestellter			Übrige/Autres
*BEJU-2011-0213	2011	BE		Rudolf		Hofer		1948		Redaktor Amt. Bulletin		CVP/PDC
*BEJU-1971-0170	1999	BE		Rudolf		Hofer		1948		Fraktionssekretär LdU/EVP	LdU/AdI
*BEJU-1947-0032	1947	BE		Werner		Bärtschi	1911		Lehrer und Sektionschef		Freisinnig-Demokratische Partei / Parti radical-démocratique
*BEJU-1975-0084	1979	BE		Werner		Bärtschi	1911		PD Dr. med., Dr. biol.		LdU/AdI
*ZH-1955-0209	1963	ZH		Hans		Pfister		1903		Maschineningenieur			Radikale / Freisinnige
*ZH-1963-0229	1963	ZH		Hans		Pfister		1903		Molkereiverwalter			Bauern-,Gewerbe- und Bürgerpartei
*ZH-1991-0598	1991	ZH		Peter		Schmid		1962		Steuerberater				SVP/UDC
*ZH-1987-0615	1987	ZH		Peter		Schmid		1962		FEAM						POCH
*
* Lukas checked all these not unique names in file "Ratsmitglieder_1848_DE.xlsx"
* Result: Only Peter Schmid appears as a former member of the National Council
* but in different canton (TG), with different birthyear (1940) and different 
* party (Grüne)


* (v) Merge and aggregate over all observations per ID

gen canton_tomerge=canton
replace canton_tomerge="BEJU" if canton_tomerge=="BE" | canton_tomerge=="JU"

merge m:1 name firstname canton_tomerge sex birthyear using ///
	"$path\02_Processed_data\04_Candidates_status\nachruecker_spells.dta", gen(merge_nachruecker)
	
* Four individuals not merged from using data. These are politicians who never
*   stood for election and entered as replacement candidates. 
* Despland	Gabriel: Von 1941 bis 1943 im Nationalrat, danach Ständerat. 
* Gabriel Theodor: Rückte nach von 1932 bis 1935. 
* Rappard William: Von 1941 bis 1943 im Nationalrat. 
* Ruh Jakob: Von 1932 bis 1935 im Nationalrat. 

forvalues i = 1(1)6 {
bysort ID: egen helpvar1= max(EDateJoining`i') 
bysort ID: egen helpvar2= max(EDateLeaving`i') 
list ID year if EDateJoining`i'!=helpvar1 & EDateJoining`i'!=. & merge_nachruecker!=2
list ID year if EDateLeaving`i'!=helpvar2 & EDateLeaving`i'!=. & merge_nachruecker!=2
replace EDateJoining`i'=helpvar1 if EDateJoining`i'==. & ID!=""
replace EDateLeaving`i'=helpvar2 if EDateLeaving`i'==. & ID!=""
drop helpvar1 helpvar2
}

preserve
import exc "$path\01_Raw_data\01_Elections_1931_1975\Wahltermine.xlsx", ///
	first clear
save  "$path\02_Processed_data\01_Elections_1931_1975\Wahltermine.dta", ///
	replace
restore

merge m:1 year using "$path\02_Processed_data\01_Elections_1931_1975\Wahltermine.dta", ///
	gen(merge_electiondate)
	

gen incumbent=0
gen incumbent_L1=0
gen incumbent_F1=0

gen tenure=0
	
forvalues i = 1(1)6 {
replace incumbent=1 if  election_date>=EDateJoining`i' & election_date<=EDateLeaving`i'
replace incumbent_L1=1 if  election_date_L1>=EDateJoining`i' & election_date_L1<=EDateLeaving`i'
replace incumbent_F1=1 if  election_date_F1>=EDateJoining`i' & election_date_F1<=EDateLeaving`i'
replace tenure=tenure+election_date-EDateJoining`i' if election_date>EDateJoining`i' & election_date<EDateLeaving`i'
}

format EDate* %dM/d/CY

* Quality check: compare our incumbent coding with BfS coding
* tab incumbent status
*           |                   Status
* incumbent |         A          B          E          N |     Total
*-----------+--------------------------------------------+----------
*         0 |    27,106         12        101        638 |    27,857 
*         1 |        13      1,736          1          0 |     1,750 
*-----------+--------------------------------------------+----------
*     Total |    27,119      1,748        102        638 |    29,607 

* Case 1: 
* br if incumbent==1 & status=="A" 
* canton	year	ID				firstname	name		incumbent	status
*AG			2015	AG-2011-0061	Beat		Flach			1		A
*AG			2015	AG-2007-0112	Hansjörg	Knecht			1		A
*AG			2015	AG-1991-0042	Ulrich		Giezendanner	1		A
*AG			2015	AG-1995-0036	Corina		Eichenberger	1		A
*AG			2015	AG-1995-0086	Ruth		Humbel Näf		1		A
*AG			2015	AG-2007-0226	Cédric		Wermuth			1		A
*AG			2015	AG-2007-9167	Yvonne		Feri			1		A
*AG			2015	AG-1999-9180	Philipp		Müller			1		A
*AG			2015	AG-1995-9072	Max			Chopard-Acklin	1		A
*AG			2015	AG-1991-0129	Luzi		Stamm			1		A
*AG			2015	AG-2011-0083	Bernhard	Guhl			1		A
*AG			2015	AG-1987-0107	Maximilian	Reimann			1		A
*AG			2015	AG-1999-0060	Sylvia		Flückiger-Bäni	1		A
*
* Result: All observations are wrongly coded in BfS data

replace status="B" if status=="A" & incumbent==1 & canton=="AG" & year==2015	

* Case 2: 
* br canton year ID firstname name incumbent status if incumbent==0 & status=="B"
*
* Case 2a: Wrong in our parlament data (corrected in  (F) READ-IN NACHRÜCKER)
*GE		1979	GE-1967-9224	Mario		Soldini		0	B
*GE		1983	GE-1967-9224	Mario		Soldini		0	B
*SG		1975	SG-1971-0039	Franz		Jaeger		0	B
*SG		1979	SG-1971-0039	Franz		Jaeger		0	B
*SG		1987	SG-1971-0039	Franz		Jaeger		0	B
*SG		1983	SG-1971-0039	Franz		Jaeger		0	B
*SG		1991	SG-1971-0039	Franz		Jaeger		0	B

* Case 2b: The following individuals were previously elected in a differnet canton
*GR		1975	GR-1975-0015	Walter		Jäger-Stamm	0	B 
* Note: Elected in 1971 in BS
*ZH		1999	ZH-1999-0875	Jean		Ziegler		0	B
* Note: Elected in GE in 1995
*ZH		1975	ZH-1975-0374	Wilfried	Naegeli		0	B
* Not: Elected in TG in 1971

* Case 2c: Other definition of incumbency as BfS
*GE		1983	GE-1975-0047	Robert		Tochon		0	B
* Note: Leaving on 01/09/1983
*BE		1979	BEJU-1955-0173	Emil		Schaffer	0	B
* Note: Leaving on 25/01/1979

* Case 3: Incumbent in our data but "Ehemalig" in BfS data 
* br canton year ID 		firstname 	name 	incumbent status if incumbent==1 & status=="E"
*ZH	1983	ZH-1975-0483	Rolf		Seiler	1			vE


* Quality check:

gen spell_check1=1 if elected==1
gen spell_check2=0 if elected==0

forvalues i = 1(1)6 {
replace spell_check1=0 if inrange(election_date+60,EDateJoining`i', EDateLeaving`i') ///
	& elected==1 & EDateJoining`i'!=. & EDateLeaving`i'!=.
replace spell_check2=1 if inrange(election_date+60,EDateJoining`i', EDateLeaving`i') ///
	& elected==0 & EDateJoining`i'!=. & EDateLeaving`i'!=.
}

*br ID year name firstname election_date EDateJoining* EDateLeaving* if spell_check1==1
*88 cases
*AG-1931-0016	1943 	Killer	Karl: Ständerat (Replacement: Adolf Aeschbach)
*AG-1931-0033	1935	Siegrist Rudolf: Wahlannahmeverzicht (Replacement: Adolf Welti)
*AG-1931-0033	1939	Siegrist Rudolf: Wahlannahmeverzicht (Replacement: Adolf Gloor)
*AG-1947-0036	1955	Richner	Adolf: Zurückgetreten (Replacement: Werner Allemann)
*AG-1947-0036	1959	Richner	Adolf: Zurückgetreten (Replacement: Arthur Schmid)
*AG-1947-0054	1955	Stöckli	Xaver: Wahlannahmeverzicht (Replacement: Robert Reimann)
*AG-1951-0040	1963	Reimann	Robert: Wahlannahmeverzicht (Replacement: Julius Binder)
*AG-1987-0107	1995	Reimann	Maximilian: Wahlannahmeverzicht (Replacement: Ernst Hasler)
*AG-1995-0033	2007	Egerszegi-Obrist Christine: Ständerat (Replacement: Corina Eichengerger-Walther)
*AG-1999-0028	2011	Bruderer Wyss	Pascale: Ständerat (Replacement: Yvonne Feri)
*AG-1999-9180	2015	Müller	Philipp: Ständerat (Replacement: Matthias Samuel Jauslin)
**BEJU-1931-01371951	Weber Max: Bundesrat (Replacement later in January 1952 by Karl Geissbühler)	
*BEJU-1931-0138	1935	Weber Rudolf: Wahlannahmeverzicht (Replacement: Alfred Held)
**BEJU-1935-00071935	Anliker	Ernst -> Check
*BEJU-1935-0023	1947	Brawand	Samuel: s.w. Wahlannahmeverzicht weil Regierungrat (Replacement: Albert Fawer)
**BEJU-1935-00501951	Feldmann Markus: Bundesrat (Replacement: Hans Stähli end of January 1952)
*BEJU-1947-0014	1955	Bauder	Robert: Wahlannahmeverzicht (Replacement: Georg Rutishauser)
**BEJU-1979-02541987	Ogi	Adolf: Bundesrat (Replacement later in February 1988 replaced by Susanna Daepp-Heiniger)
*BEJU-1987-0459	2011	Stöckli	Hans: Ständerat (Replacement Alexander Tschäppät)
*BEJU-1999-0371	2003	Sommaruga Simonetta: Ständerat (Replacement Margret Kiener Nellen)
**BL-1931-0003	1951	Blunschi Jules: Verzicht auf das Nationalratsmandat (later January 1952 replaced by Josef Tschopp)
*BL-1939-0010	1939	Hilfiker Walter: Wahlannahmeverzicht (Replacement: Leo Mann)
*BL-1987-0029	2007	Janiak Claude: Ständerat (Replacement Eric Nussbaumer)
*BS-1943-0024	1943	Moeschlin Felix: previously elected in ZH (1939), then in BS (1943) 
*BS-1983-0022	2003	Fetz Anita: Ständerat (Replacement Silvia Schenker)
*FR-1931-0007	1931	Chassot	Charles: Wahlannahmeverzicht (Replacement: Joseph Delatena)
*FR-1935-0013	1939	Musy Jean-Marie: Wahlannahmeverzicht (Replacement: Robert Colliard)
*FR-1943-0004	1943	Bondallaz Paul: Unbekanntes Ausscheiden (Replacement 1944: Jakob Meyer)
*FR-1955-0013	1955	Ducotterd Georges: Wahlannahmeverzicht (Replacement: Robert Colliard, mandate incompability)
*GE-1951-9068	1951	Treina Jean: Wahlannahmeverzicht weil Beamte (Replacement: Georges Borel)
*GE-1959-7045	1983	Duboule	Gilbert: Gestorben, 6.11.1983 (Replacement: Revaclier	Jean)
*GE-1987-0006	1995	Brunner	Christiane: Ständerat (Replacement: Maria Roth Bernasconi)
*LU-1943-0007	1955	Clavadetscher Christian: Ständerat (Replacement: Niklaus Honauer)
*NE-2003-0009	2007	Burkhalter	Didier: Ständerat (Replacement: Laurent Favre)
*SG-1935-0008	1935	Duttweiler	Gottlieb: Wahlannahmeverzicht (Replacement: Ulrich Eggenberger)
*SG-1947-0019	1971	Eggenberger	Mathias: Ständerat (Replacement: Sahlfeld-Singer Hanna)
**SG-1951-0014	1971	Furgler	Kurt: Bundesrat (Replacement later in February 1972 replaced by Remigius Kaufmann)
*SG-1983-0012	1999	David Eugen: Ständerat (Replacement: Felix Walker)
*SG-1983-0041	2011	Rechsteiner	Paul: Ständerat (Replacement: Barbara Gysi)
*SO-1955-0017	1963	Ritschard Willi: Wahlannahmeverzicht (Replacement: Otto Stich)
*SO-1975-9021	1999	Leuenberger	Ernst: Ständerat (Replacement: Boris Banga)
*SO-1987-0005	2011	Bischof	Pirmin: Ständerat (Replacement: Urs Schläfli-Kocher)
*TG-1931-0018	1943	Roth August:  Annahmeverzicht(Replacement: Schümperli)
*TG-1939-0011	1951	Holliger Hans: Wahlannahmeverzicht(Replacement: Walter Tuchschmid)
*TG-1999-0043	2011	Häberli-Koller Brigitte: Ständerat (Replacement: Christian Lohr)
*TI-1931-9002	1935	Olgiati	Camillo: Unbekannte Gründe (Replacement: Francesco Rusca)
*TI-1955-9010	1963	Stefani	Alberto: Ständerat (Replacement: Ugo Gianella)
*TI-1971-0029	1975	Martinelli Pietro: Wahlannameverzeicht (Replacement: Werner Carobbio)
*VD-1931-0043	1931	Perret Paul: Rücktritt 22. Dez 1931 wegen Staatsrat (Replacement: Pierre Rochat)
*VD-1947-0014	1963	Bussey Alfred: Wahlannameverzeicht (Replacement: Marcel Bravand)
*VD-1947-0055	1947	Muret André: Unbekannte Gründe (ev. wegen Stadtrat Lausanne, Replacement: Armand Forel)
*VD-1947-0058	1947	Peitrequin Jean: Unbekannte Gründe (ev. wegen Stadtrat Lausanne, Replacement: Ulysse Péclard)
**VD-1975-7090	1983	Delamuraz Jean-Pascal: Bundesrat (Replacement later in March 1984 by Pierre Savary)
*VD-1979-0067	1987	Jaggi Yvette: Ständerat (Replacement: Pierre Aguet)
*VD-1983-0033	1999	Béguelin Michel: Ständerat (Replacement: Pierre Tillmanns)
*VD-1983-0033	2003	Béguelin Michel: Ständerat (Replacement: Pierre Salvi)
*VD-1983-0083	2007	Huguenin Marianne: Wahlannameverzeicht (Replacement: Josef Zisyadis)
*VD-1983-0135	2007	Recordon Luc: Ständerat (Replacement: Christian van Singer)
*VD-1983-0135	2011	Recordon Luc: Ständerat (Replacement: Christian van Singer)
*VD-1983-0135	2015	Recordon Luc: Wahlannameverzeicht (Replacement:  Adèle Thorens Goumaz)
*VD-1987-0091	1999	Langenberger Christiane: Ständerat (Replacement: René Vaudroz)
*VD-1987-0091	2003	Langenberger Christiane: Ständerat (Replacement: René Vaudroz)
*VD-1999-0216	2007	Savary Géraldine: Ständerat (Replacement: Ada Mara)
*VD-1999-0216	2011	Savary Géraldine: Ständerat (Replacement: Jean Christoph Schwaab)
*VD-1999-0216	2015	Savary Géraldine: Ständerat (Replacement: Jean Christoph Schwaab)
*VD-2007-0071	2015	Français Olivier: Ständerat (Replacement: Fathi Derber)
**VS-1931-0004	1931	Escher Joseph: Unbekannte Gründe (ev. wegen Staatsrat, Replacement: s.w. Rudolf Metry weil keine Ersatzkandidaten auf Liste, die nicht gewählt wurden)
*VS-1931-0008	1935	Kuntschen Joseph: Unbekannte Gründe (ev. wegen Gemeindepräsident Sion, Replacement: Rodolphe Metry)
**VS-1951-0011	1951	Guntern	Leo: Wahlannameverzeicht (Replacement: Leo Stoffel end of January 1952)
*ZH-1931-0110	1943	Nobs Ernst: Bundesrat (Replacement: later in March 1944 by Valentin Gitermann)
*ZH-1935-0045	1935	Duttweiler Gottlieb: Wahlannahmeverzicht (Replacement: Heinrich Schnyder)
**ZH-1935-0045	1939	Duttweiler Gottlieb: Wahlannahmeverzicht (Replacement: Felix Moeschlin)
*ZH-1935-0145	1951	Meier Rudolf: Wahlannahmeverzicht (Replacement: Heinrich Brändli)
*ZH-1947-0129	1955	König Walter: Wahlannahmeverzicht (Replacement: Rudolf Schmid)
*ZH-1975-0554	1987	Weber Monika: Wahlannahmeverzicht (Replacement: Roland Wiederkehr)
**ZH-1979-0051	2003	Blocher	Christoph: Bundesrat (Replacement, later in March 2003 replaced by Hans Rutschmann)
*ZH-1987-0274	2003	Heberlein Trix: Ständerat (Replacement: Markus Hutter)
*ZH-1987-9190	2007	Diener	Verena: Ständerat (Replacement: Thomas Weibel)
*ZH-1991-0250	2007	Gutzwiller	Felix: Ständerat (Replacement: Markus Hutter)
*ZH-2003-0420	2015	Jositsch Daniel: Ständerat (Replacement: Angelo Barrile)
*ZH-2003-0602	2015	Noser Ruedi: Ständerat (Replacement: Hans-Ulrich Bigler)


*br ID year name firstname election_date EDateJoining* EDateLeaving*  if spell_check2==1

* 67 cases
*AG-1931-0041	1935	Welti	Adolf  (Replacement of Rudolf Siegrist)
*AG-1935-0019	1939	Gloor	Adolf (Replacement of Rudolf Siegrist)
*AG-1939-0003	1943	Aeschbach	Adolf (Replacement of Karl Killer)
*AG-1939-0004	1955	Allemann	Werner (Replacement of Adolf Richner)
*AG-1951-0040	1955	Reimann	Robert (Replacement of Stöckli Xaver)
*AG-1955-0007	1963	Binder	Julius (Replacement of Reimann	Robert)
*AG-1963-9003	1959	Schmid	Arthur (Replacement of Adolf Richne)
*AG-1995-0036	2007	Eichenberger-Walther Corina (Replacement of Egerszegi-Obrist Christine)
*AG-1995-9071	1995	Hasler	Ernst (Replacement of Reimann Maximilian)
*AG-2007-9167	2011	Feri	Yvonne (Replacement of Pascale Bruderer)
*AG-2011-0119	2015	Jauslin	Matthias Samuel (Replacement of Müller Philipp)
*BEJU-1931-0056	1935	Held	Alfred (Replacement of Weber Rudolf)
*BEJU-1943-0053	1947	Fawer	Albert (Replacement of Samuel Brawand)
*BEJU-1955-0172	1955	Rutishauser	Georg (Replacement of Bauder Robert)
*BEJU-1991-0510	2011	Tschäppät	Alexander (Replacement of StöckliHans)
*BEJU-2003-0227	2003	Kiener Nellen Margret (Replacement of Simonetta Sommaruga)
*BL-1931-0009	1939	Mann	Leo (Replacement of Hilfiker Walter)
*BL-2003-0056	2007	Nussbaumer	Eric (Replacement of Claude Janiak)
*BS-1995-9024	2003	Schenker	Silvia (Replacement of Anita Fetz)
*FR-1931-0009	1931	Delatena	Joseph (Replacement of Chassot	Charles)
*FR-1931-9003	1939	Colliard	Robert (Replacement of Musy Jean-Marie)
*FR-1931-9003	1955	Colliard	Robert (Replacement of Ducotterd Georges)
*GE-1939-9007	1951	Borel	Georges (Replacement of Treina Jean)
**GE-1955-0016	1983	Magnin	Armand (Replacement of Roger Dafflon, but in 1984 -> corrected above)
*GE-1975-9062	1983	Revaclier	Jean (Replacement of Duboule Gilbert)
*GE-1991-0091	1995	Roth Bernasconi	Maria (Replacement of Brunner	Christiane)
*LU-1951-0009	1955	Honauer	Niklaus (Replacement of Clavadetscher Christian)
*NE-2007-0021	2007	Favre	Laurent (Replacement of Burkhalter	Didier)
*SG-1931-0011	1935	Eggenberger	Ulrich (Replacement of Duttweiler Gottlieb)
*SG-1971-0066	1971	Sahlfeld-Singer	Hanna (Replacement of Eggenberger Mathias)
*SG-1995-0071	2011	Gysi	Barbara (Replacement of Rechsteiner	Paul)
*SG-1999-0163	1999	Walker	Felix (Replacement of Eugen David)
*SO-1959-0016	1963	Stich	Otto (Replacement of Willi Ritschard)
*SO-1995-0002	1999	Banga	Boris (Replacement of Leuenberger	Ernst)
*SO-2011-7069	2011	Schläfli-Kocher	Urs (Replacement of Bischof	Pirmin)
*TG-1935-7079	1951	Tuchschmid	Walter (Replacement of Holliger Hans)
*TG-1943-0024	1943	Schümperli	Rudolf (Replacement of Roth August)
*TG-2003-0039	2011	Lohr	Christian (Replacement of Brigitte Häberli-Koller)
*TI-1931-0021	1935	Rusca	Giovanni Battista (Replacement of OlgiatiCamillo)
*TI-1963-0006	1975	Carobbio	Werner (Replacement of Martinelli Pietro)
*TI-1963-0014	1963	Gianella	Ugo (Replacement of Stefani	Alberto)
*VD-1931-0050	1931	Rochat	Pierre (Replacement of Paul Perret)
*VD-1935-0010	1963	Brawand	Marcel (Replacement of Bussey Alfred)
*VD-1943-0044	1947	Péclard	Ulysse (Replacement of Peitrequin Jean)
*VD-1947-0032	1947	Forel	Armand (Replacement of André Muret)
*VD-1983-7086	2007	van Singer	Christian (Replacement of Recordon Luc)
*VD-1983-7086	2011	van Singer	Christian (Replacement of Recordon Luc)
*VD-1987-0002	1987	Aguet	Pierre (Replacement of Jaggi Yvette)
*VD-1987-0167	2007	Zisyadis	Josef (Replacement of Huguenin Marianne)
*VD-1987-9011	1999	Tillmanns	 (Replacement of Béguelin Michel)
*VD-1995-0191	1999	Vaudroz	René (Replacement of Langenberger Christiane)
*VD-1995-0191	2003	Vaudroz	René (Replacement of Langenberger Christiane)
*VD-1999-0139	2007	Marra	Ada (Replacement of Savary Géraldine)
*VD-1999-0211	2003	Salvi	Pierre (Replacement of Béguelin Michel)
*VD-2003-0172	2011	Schwaab	Jean Christophe (Replacement of Savary Géraldine)
*VD-2003-0172	2015	Schwaab	Jean Christophe (Replacement of Savary Géraldine)
*VD-2003-0181	2015	Thorens Goumaz	Adèle (Replacement of Recordon Luc)
*VD-2011-0093	2015	Derder	Fathi (Replacement of Français Olivier)
*VS-1935-0011	1935	Metry	Rodolphe (Replacement of Kuntschen Joseph)
*ZH-1939-0029	1951	Brändli	Heinrich (Replacement of Meier Rudolf)
*ZH-1943-0193	1955	Schmid	Rudolf (Replacement of König Walter)
*ZH-1987-0773	1987	Wiederkehr	Roland (Replacement of Weber Monika)
*ZH-1999-0374	2003	Hutter	Markus (Replacement of Trix Heberlein)
*ZH-1999-0374	2007	Hutter	Markus (Replacement of Felix Gutzwiler)
*ZH-1999-9179	2015	Bigler	Hans-Ulrich (Replacement of Ruedi Noser)
*ZH-2003-0894	2007	Weibel	Thomas (Replacement of Diener Verena)
*ZH-2015-0038	2015	Barrile	Angelo (Replacement of Daniel Jositsch)

* (xi) Outfile for ground truth 1 (June 2024)

preserve
clear
set obs 73
gen year=1930+_n
save "$path\02_Processed_data\02_Elections_1971_2015\temp5.dta", replace
restore

preserve
duplicates drop ID, force
keep ID EDate*
cross using "$path\02_Processed_data\02_Elections_1971_2015\temp5.dta"
gen office_year=0
gen year_start=mdy(1,1,year)
gen year_end=mdy(12,31,year)

forvalues i = 1(1)6 {
replace office_year=1 if  year_start>=EDateJoining`i' & year_end<=EDateLeaving`i'
}
keep if office_year==1
keep ID year
rename ID id_0
rename year year_1
save "$path\02_Processed_data\02_Elections_1971_2015\OfficeYears_GT1.dta", replace
restore

erase "$path\02_Processed_data\02_Elections_1971_2015\temp5.dta"


* (xii) Outfile for ground truth 2 (June 2024)

preserve
clear
set obs 22
gen year=1981+_n
save "$path\02_Processed_data\02_Elections_1971_2015\temp5.dta"
restore

preserve
duplicates drop ID, force
keep ID EDate*
cross using "$path\02_Processed_data\02_Elections_1971_2015\temp5.dta"
gen office_year=0
gen year_start=mdy(1,1,year)
gen year_end=mdy(12,31,year)

forvalues i = 1(1)6 {
replace office_year=1 if  year_start>=EDateJoining`i' & year_end<=EDateLeaving`i'
}
keep if office_year==1
keep ID year
rename year Year
rename ID NRid
save "$path\02_Processed_data\02_Elections_1971_2015\OfficeYears_GT2.dta", replace
restore

erase "$path\02_Processed_data\02_Elections_1971_2015\temp5.dta"



******************************************************
* (L) ADD LIST OFFICIAL POSITION
******************************************************

egen listIDhelp = group(canton year list)  // generate artificial unique list ID for all years withouth this info (year<1971)
replace listIDhelp = listIDhelp + 900000
tostring listIDhelp, replace
gen listno_orig = listno  // store original values of variables
gen candidateno_orig = candidateno 
replace listno = listIDhelp if listno=="" // temporary addition of artificial list no. to avoid missing values (year<1971) for subsequent merging
bysort listIDhelp: gen listhelp=_n
replace candidateno = listhelp if candidateno==.
duplicates report canton year listno candidateno

preserve
do "$path\03_Code\02_Elections_1971_2015\Setup_listpositions.do"
restore

merge 1:1 canton year listno candidateno using ///
"$path\02_Processed_data\02_Elections_1971_2015\Listenpositionen_1975-2015.dta", ///
	gen(merge_listpos)
drop if merge_listpos==2  // only "vereinzelte" of UR 2007
/*
Master only: merge_listpos == 1: 11,959 obs
-- > except for 13 obs (listed below), all from period before 1975
-- > checked in list position data: there is no information on list positions
canton	year	firstname	name
AR		1979	Christian	Merz
AR		1979	Hans-Rudolf	Früh
AR		1987	Hans-Rudolf	Früh
AR		1987	Herbst		Maeder
OW		1995	Rudolf O.	Durrer
OW		1999	Adalbert	Durrer
UR		1995	Yves		Merminod
*/

drop listIDhelp listhelp
replace listno = listno_orig  // restore the original variable (including missings)
replace candidateno = candidateno_orig
drop listno_orig candidateno_orig
drop list_nr name_listpos firstname_listpos merge_listpos

****************************************************************
* (M) GENERATE TIME VARIABLES AND ADDITIONAL COVARIATES 
****************************************************************

* (i) Generate time variables

egen ID_time=group(year)
egen ID_pers=group(ID)

drop if year==2019 
* Note: one observation from 2019 (from date of election) is dropped 
drop if ID_pers==. 
* Note: Drop 53 Vereinzelte and four replacement candidates. 
drop if votes==.
* Note: drop 92 candidates that were elected in "stille Wahlen"


xtset ID_pers ID_time

bysort ID_pers: egen year_min=min(year) 
gen first_participation = year == year_min
replace first_participation = 0 if cand_before1931==1

sum year
local maxyear=r(max)
local minyear=r(min)

gen participation_L1=0
replace participation_L1=1 if L1.elected!=. & year!=`minyear'
replace participation_L1=. if year==`minyear' // set to missing for first year of dataset

gen participation_F1=0
replace participation_F1=1 if F1.elected!=. & year!=`maxyear'
replace participation_F1=. if  year==`maxyear' // set to missing for last year of dataset



local vars "elected" 
foreach var in `vars'{			
gen `var'_F1=F1.`var'
gen `var'_L1=L1.`var'
replace `var'_F1=0 if `var'_F1==. & year!=`maxyear'
replace `var'_L1=0 if `var'_L1==. & year!=`minyear'
}


* check timing variabls

gen group_timing=.
replace group_timing=1 if elected==0 & incumbent_F1==0 & participation_F1==0 & elected_F1==0
replace group_timing=2 if elected==0 & incumbent_F1==0 & participation_F1==0 & elected_F1==1 // this case should not be possible to be elected with no participation
replace group_timing=3 if elected==0 & incumbent_F1==0 & participation_F1==1 & elected_F1==0
replace group_timing=4 if elected==0 & incumbent_F1==0 & participation_F1==1 & elected_F1==1 
replace group_timing=5 if elected==0 & incumbent_F1==1 & participation_F1==0 & elected_F1==0
replace group_timing=6 if elected==0 & incumbent_F1==1 & participation_F1==0 & elected_F1==1 // this case should not be possible to be elected with no participation
replace group_timing=7 if elected==0 & incumbent_F1==1 & participation_F1==1 & elected_F1==0
replace group_timing=8 if elected==0 & incumbent_F1==1 & participation_F1==1 & elected_F1==1
replace group_timing=9 if elected==1 & incumbent_F1==0 & participation_F1==0 & elected_F1==0
replace group_timing=10 if elected==1 & incumbent_F1==0 & participation_F1==0 & elected_F1==1 // this case should not be possible to be elected with no participation
replace group_timing=11 if elected==1 & incumbent_F1==0 & participation_F1==1 & elected_F1==0
replace group_timing=12 if elected==1 & incumbent_F1==0 & participation_F1==1 & elected_F1==1
replace group_timing=13 if elected==1 & incumbent_F1==1 & participation_F1==0 & elected_F1==0
replace group_timing=14 if elected==1 & incumbent_F1==1 & participation_F1==0 & elected_F1==1 // this case should not be possible to be elected with no participation
replace group_timing=15 if elected==1 & incumbent_F1==1 & participation_F1==1 & elected_F1==0
replace group_timing=16 if elected==1 & incumbent_F1==1 & participation_F1==1 & elected_F1==1

tab group_timing, missing // note: 2,6,10, and 14 should not be possible

label define group_timings  1 "el=0 & in_F1=0 & pa_F1=0 & el_F1=0" ///
							2 "el=0 & in_F1=0 & pa_F1=0 & el_F1=1" ///
							3 "el=0 & in_F1=0 & pa_F1=1 & el_F1=0" ///
							4 "el=0 & in_F1=0 & pa_F1=1 & el_F1=1" /// 
							5 "el=0 & in_F1=1 & pa_F1=0 & el_F1=0" ///
							6 "el=0 & in_F1=1 & pa_F1=0 & el_F1=1" ///
							7 "el=0 & in_F1=1 & pa_F1=1 & el_F1=0" ///
							8 "el=0 & in_F1=1 & pa_F1=1 & el_F1=1" ///
							9 "el=1 & in_F1=0 & pa_F1=0 & el_F1=0" ///
							10 "el=1 & in_F1=0 & pa_F1=0 & el_F1=1" /// 
							11 "el=1 & in_F1=0 & pa_F1=1 & el_F1=0" /// 
							12 "el=1 & in_F1=0 & pa_F1=1 & el_F1=1" ///
							13 "el=1 & in_F1=1 & pa_F1=0 & el_F1=0" ///
							14 "el=1 & in_F1=1 & pa_F1=0 & el_F1=1" /// 
							15 "el=1 & in_F1=1 & pa_F1=1 & el_F1=0" ///
							16 "el=1 & in_F1=1 & pa_F1=1 & el_F1=1"

label value group_timing group_timings

bysort year: tab elected_F1

xtset ID_pers ID_time

forvalues i=1931(4)2011{
display "Year: `i'"
tab elected_F1 if year==`i'
sum elected if year==(`i'+4) & L.sex==. & elected==1 
}

* Result: Number of elected_F1=1 corresponds with total elected candidates in 
* next election minus those who did not run in current election

* (ii) Generate additional covariates

sort ID year
bysort ID: gen candidacyno=_n

preserve
import exc "$path\01_Raw_data\08_Municipalities\Municipalities_Population_1981_2015.xlsx", ///
	first clear
destring municipalityno year,replace
rename population populationmun
keep populationmun year municipalityno
save "$path\02_Processed_data\08_Municipalities\Municipalities_Population_1981_2015.dta", replace
restore

merge m:1 municipalityno year using ///
	"$path\02_Processed_data\08_Municipalities\Municipalities_Population_1981_2015.dta", gen(mergemun)

tab year if mergemun==1
drop if mergemun==2

replace sex="1" if sex=="M"
replace sex="0" if sex=="F"
destring sex, replace force
gen age=year-birthyear

* (iii) Add information on running variable 

save "$path\02_Processed_data\nationalraete_1931_2015_temp.dta", replace

*do "$path\03_Code\13_Running_variable\01_Vote_margin_analytical.do"
*do "$path\03_Code\13_Running_variable\02_Vote_margin_simulated.do"

import delimited  "$path\02_Processed_data\13_Running_variable\rv_analytical_r.csv", ///
	delim(";") clear
rename id ID 	
keep ID votemargin year
save "$path\02_Processed_data\13_Running_variable\rv_analytical_r.dta", replace



use "$path\02_Processed_data\nationalraete_1931_2015_temp.dta", clear
merge 1:1 ID year using ///
	"$path\02_Processed_data\13_Running_variable\rv_analytical_r.dta", gen(merge_running) // note: this is the version of june 2020 (corrected for ZH 1967)
* br ID first name elected votemargin_rel  if (elected==1 &  votemargin_rel<0) | (elected==0 &  votemargin_rel>0) & !missing(votemargin_rel)
gen votemargin_rel=votemargin/eligible_cant

bysort canton year: egen votes_sum=sum(votes)
bysort canton year: egen no_seats=sum(elected)
gen votemargin_rel_alt=votemargin/valid_cant 
bysort canton year list: gen indi=_n 
bysort canton year: egen help1=sum(pvotes) if indi==1
bysort canton year: egen pvotes_sum=max(help1)  
gen votemargin_rel_alt_2=votemargin/(votes_sum/no_seats)
drop help1 votes_sum pvotes_sum seats_cant indi

corr votemargin_rel votemargin_rel_alt votemargin_rel_alt_2

sum year
local maxyear=r(max)
local minyear=r(min)

sort ID_pers ID_time

gen votemargin_L1=.
replace votemargin_L1= L1.votemargin if  year!=`minyear'
gen voters_cant_L1=.
replace voters_cant_L1= L1.voters_cant if  year!=`minyear'
gen eligible_cant_L1=.
replace eligible_cant_L1= L1.eligible_cant if  year!=`minyear'
gen no_seats_L1=.
replace no_seats_L1= L1.no_seats if  year!=`minyear'

	
* (iv) Labeling and save dataset

label var ID "Candidate ID"
label var cand_before1931 "Candidate participated in 1925 and/or 1928 election"
label var populationmun "Population in municipality of candidate"
label var elected_F1 "Elected in next election"
label var elected_L1 "Elected in previous election"
label var participation_F1 "Participation in next election"
label var participation_L1 "Participation in previous election"
label var ID_time "Numeric ID variable for time dimension"
label var year_min "First observed participation of a candidate (since 1931)"
label var first_participation "First participation of a candidate (since 1925)"
label var ID_pers "Numeric ID variable for candidate dimension"
label var tenure "Number of days in office at election date"
label var incumbent "Incumbent dummy at election date"
*label var partymarg "Party margin"
*label var votesdiff "Candidate margin"
label var votemargin "Vote margin"
label var votemargin_rel "Running variable (relative vote margin)"
label var age "Age"
label var no_seats "Total no. of seats in canton"
label var votemargin_L1 "Vote margin in previous election"
label var voters_cant_L1 "No of actual voters in previous election"
label var eligible_cant_L1 "No of eligble voters in previous election"
label var alliance "Alliance (Listenverbindung)"
label var suballiance "Suballiance (Unterlistenverbindung)"
label var pvotes "Party votes"
label var firstname "First name of candidate"
label var name "Surname of candidate"
label var list "Party list"
label var originno1 "Municipality no. of first origin municipality"
label var originno2 "Municipality no. of second origin municipality"
label var originno3 "Municipality no. of third origin municipality"
label var originno4 "Municipality no. of fourth origin municipality"
label var originno5 "Municipality no. of fifth origin municipality"
label var originno6 "Municipality no. of sixth origin municipality"
label var EDateJoining1 "Date of Joining National Council, First Spell"
label var EDateLeaving1 "Date of Leaving National Council, First Spell"
label var EDateJoining2 "Date of Joining National Council, Second Spell"
label var EDateLeaving2 "Date of Leaving National Council, Second Spell"
label var EDateJoining3 "Date of Joining National Council, Third Spell"
label var EDateLeaving3 "Date of Leaving National Council, Third Spell"
label var EDateJoining4 "Date of Joining National Council, Fourth Spell"
label var EDateLeaving4 "Date of Leaving National Council, Fourth Spell"
label var EDateJoining5 "Date of Joining National Council, Fifth Spell"
label var EDateLeaving5 "Date of Leaving National Council, Fifth Spell"
label var EDateJoining6 "Date of Joining National Council, Sixth Spell"
label var EDateLeaving6 "Date of Leaving National Council, Sixth Spell"



* (v) Final corrections

list name firstname year if elected==1 & votemargin<0
list name firstname year if elected==1 & votemargin<0

replace elected=0 if year==1939 & name=="Musy" & firstname=="Jean-Marie" // Initially wrong results delivered by the cantonal government of FR, later corrected (see file "Spezialfälle.pdf")
replace elected=1 if year==1939 & name=="Colliard" & firstname=="Robert"

* Note: as a consequence, we also need to re-create the elected_F1 variables
*       (April 5, 2023)

drop elected_F1 elected_L1

local vars "elected" 
foreach var in `vars'{			
gen `var'_F1=F1.`var'
gen `var'_L1=L1.`var'
replace `var'_F1=0 if `var'_F1==. & year!=`maxyear'
replace `var'_L1=0 if `var'_L1==. & year!=`minyear'
}


replace partyname="" if year<1971
rename partyname listname_bfs
label var listname_bfs "Party list name according to BfS (>=1971)"

* (vi) Keep only selected variables

rename cantonno ID_canton
local vars ID year canton /// // e_id_sug n_id_sug 
	  name firstname birthyear age sex job /// 
	  votes elected pvotes votemargin votemargin_rel incumbent tenure ///
	  list listname_bfs alliance suballiance ///
	  eligible_cant voters_cant no_seats populationmun ///
	  cand_before1931 ///
	  ID_time ID_pers ID_canton year_min first_participation ///
	  participation_L1 participation_F1 elected_F1 elected_L1 ///
	  votemargin_L1 voters_cant_L1 eligible_cant_L1 ///
	  municipalityno originno1 originno2 originno3 originno4 originno5 listpos ///
	  originno6 EDateJoining1-EDateJoining6 EDateLeaving1-EDateLeaving6
	  
order `vars'
keep `vars'

save "$path\02_Processed_data\nationalraete_1931_2015.dta", replace

erase "$path\02_Processed_data\04_Candidates_status\nachruecker_spells.dta"
erase "$path\02_Processed_data\01_Elections_1931_1975\Wahltermine.dta"
erase "$path\02_Processed_data\08_Municipalities\Municipalities_Population_1981_2015.dta"
erase "$path\02_Processed_data\nationalraete_1931_2015_temp.dta"
erase "$path\02_Processed_data\13_Running_variable\rv_analytical_r.dta"

/*

******************************************************
* (M) RECORD LINKAGE CHECKS
******************************************************

*  (i) Checks false negatives

use "$path\02_Processed_data\11_Directors_1994_2018\bisnode_out_wide.dta", clear

tostring birthyear1 birthyear2, replace
replace birthyear2="" if birthyear2=="."

local vars "name firstname birthyear job partyname municipality origin"
foreach var in `vars'{
local counter=1
capture drop `var'
forvalues i=1(1)10{
	di "`var'"`i'
	capture confirm variable `var'`i'
		if !_rc {
			if `counter'==1{
			gen `var'=`var'`i'
			}
			else{
			replace `var'=`var'+", "+`var'`i' if `var'`i'!=""
			}
		}
		else {
		}
	local counter=`counter'+1
}
}


set seed 1234
drop if year<1991
sort ID year
generate random = runiform() 
sort random
egen random_rank=rank(random)	

order ID canton sex name firstname birthyear municipality origin job partyname random_rank 
keep ID canton sex name firstname birthyear municipality origin job partyname random_rank

label var canton "Kanton"
label var sex "Geschlecht"
label var name "Name"
label var firstname "Vorname"
label var birthyear "Geburtsjahr"
label var municipality "Wohngemeinde"
label var origin "Heimatort"
label var job "Beruf"
label var partyname "Parteiname"


* Ladina Oester and Julian Koller
/*
forvalues i=1(1)10{
preserve
keep if inrange(random_rank,(`i'-1)*100+1,(`i')*100)
drop random_rank
sort ID
export excel  using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives`i'.xlsx",  replace  firstrow(varlabels)
restore
} 



* Peter Wallmüller

forvalues i=1(1)5{
preserve
keep if inrange(random_rank,(`i'-1)*90+1000,(`i')*90+1001) | random_rank==(`i'-1)*100+11 | random_rank==(`i'-1)*100+22 | random_rank==(`i'-1)*100+55 | random_rank==(`i'-1)*100+66 | random_rank==(`i'-1)*100+77 | random_rank==(`i'-1)*100+511 | random_rank==(`i'-1)*100+522 | random_rank==(`i'-1)*100+555 | random_rank==(`i'-1)*100+566 | random_rank==(`i'-1)*100+577
drop random_rank
sort ID
export excel  using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives1`i'.xlsx",  replace  firstrow(varlabels)
restore
preserve
keep if random_rank==(`i'-1)*100+11 | random_rank==(`i'-1)*100+22 | random_rank==(`i'-1)*100+55 | random_rank==(`i'-1)*100+66 | random_rank==(`i'-1)*100+77 | random_rank==(`i'-1)*100+511 | random_rank==(`i'-1)*100+522 | random_rank==(`i'-1)*100+555 | random_rank==(`i'-1)*100+566 | random_rank==(`i'-1)*100+577
capture append using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\overlap_wallmueller_and_oesterkoller.dta"
save "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Overlap_Wallmueller_and_OesterKoller.dta", replace
restore 
}  

use "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\overlap_wallmueller_and_oesterkoller.dta", clear
sort ID
*export excel  using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Overlap_Wallmueller_and_OesterKoller.xlsx",  replace  firstrow(varlabels)
erase "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\overlap_wallmueller_and_oesterkoller.dta"
*/ 

preserve
import excel "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Overlap_Wallmueller_and_OesterKoller.xlsx", clear first
keep ID
save "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Overlap_Wallmueller_and_OesterKoller.dta", replace
restore

merge 1:1 ID using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Overlap_Wallmueller_and_OesterKoller.dta"
erase "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Overlap_Wallmueller_and_OesterKoller.dta"

drop if _merge==3 // drop those that were already double checked in the first round
drop _merge

sum random_rank if inrange(random_rank,1001,1020) // look for 20 observations in each set
sum random_rank if inrange(random_rank,478,497)
sum random_rank if inrange(random_rank,950,971)


preserve // Julian Koller: Check_False_Negatives16
keep if inrange(random_rank,1001,1020) | inrange(random_rank,478,499) | inrange(random_rank,1601,1658) 
drop random_rank
sort ID
export excel  using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives16.xlsx",  replace  firstrow(varlabels)
*save "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives16.dta", replace
restore

preserve // Peter Wallmüller: Check_False_Negatives17
keep if inrange(random_rank,478,497)| inrange(random_rank,950,975)  | inrange(random_rank,1661,1716) 
drop random_rank
sort ID
export excel  using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives17.xlsx",  replace  firstrow(varlabels)
*save "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives17.dta", replace
restore

 // check whether doubletten has worked
/*
preserve
import excel "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives11.xlsx", clear first
merge 1:1 ID using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives16.dta",gen(merge1)
merge 1:1 ID using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives17.dta",gen(merge2)
import excel "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives5.xlsx", clear first
merge 1:1 ID using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives16.dta",gen(merge1)
merge 1:1 ID using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives17.dta",gen(merge2)
restore

erase "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives16.dta"
erase "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Check_False_Negatives17.dta"
*/

preserve
keep if inrange(random_rank,478,497)| inrange(random_rank,950,975) | inrange(random_rank,1001,1020)
gen Overlap="Oester-Koller-Wallmüller" if inrange(random_rank,478,497)
replace Overlap="Koller-Wallmüller" if inrange(random_rank,950,975) | inrange(random_rank,1001,1020)
export excel  using "$path\02_Processed_data\12_Record_linkage\01_Check_False_Negatives_out\Overlap_Wallmueller_and_OesterKoller_Round2.xlsx",  replace  firstrow(varlabels)
restore


*  (iii) Checks false sugarcube data (Nemo Krüger, Laura Decet und Peter Wallmüller)

use "$path\02_Processed_data\11_Directors_1994_2018\bisnode_out_wide.dta", clear

drop if (birthyear1!=. & birthyear1>=1986) | (birthyear2!=. & birthyear2>=1986)
* drop those that were 18 years or below in 2003 (last year of sugarcube data)
* and are thus unlikely to have a directorship

tostring birthyear1 birthyear2, replace
replace birthyear2="" if birthyear2=="."

local vars "name firstname birthyear job partyname municipality origin"
foreach var in `vars'{
local counter=1
capture drop `var'
forvalues i=1(1)10{
	di "`var'"`i'
	capture confirm variable `var'`i'
		if !_rc {
			if `counter'==1{
			gen `var'=`var'`i'
			}
			else{
			replace `var'=`var'+", "+`var'`i' if `var'`i'!=""
			}
		}
		else {
		}
	local counter=`counter'+1
}
}

preserve
use "$path\02_Processed_data\02_Elections_1971_2015\nationalraete_1931_2015_Geo", clear
tostring E_CNTR_w N_CNTR_w, replace
gen GdeNr_CNTR=E_CNTR_w+"/"+N_CNTR_w

local vars "GdeNr_CNTR"
* Note: Create a E_CNTR_w N_CNTR_w
* which we then append to original dataset.

foreach var in `vars'{	
bysort ID `var': gen indi_`var'=_n		
keep if indi_`var'==1
* Note: Keep only one observation per ID-var group.

bysort ID: gen observation_number =_n
keep ID `var' observation_number
reshape wide `var', i(ID) j(observation_number) 
save "$path\02_Processed_data\02_Elections_1971_2015\\`var'", replace
}
restore

 
merge m:1 ID using "$path\02_Processed_data\02_Elections_1971_2015\GdeNr_CNTR.dta", ///
	gen(merge_gdenr_cntr)
erase "$path\02_Processed_data\02_Elections_1971_2015\GdeNr_CNTR.dta"
gen GdeNr_CNTR=GdeNr_CNTR1
replace GdeNr_CNTR=GdeNr_CNTR+", "+GdeNr_CNTR2 if GdeNr_CNTR2!=""
replace GdeNr_CNTR=GdeNr_CNTR+", "+GdeNr_CNTR3 if GdeNr_CNTR3!=""
replace GdeNr_CNTR=GdeNr_CNTR+", "+GdeNr_CNTR4 if GdeNr_CNTR4!=""

drop GdeNr_CNTR1-GdeNr_CNTR4


save "$path\02_Processed_data\11_Directors_1994_2018\sugarcube_rl_all_input.dta", replace // input for search for sugarcubedata



set seed 12345678

sort ID year
generate random = runiform() 
sort random
egen random_rank=rank(random)	

order ID canton sex name firstname birthyear municipality GdeNr_CNTR origin job partyname random_rank 
keep ID canton sex name firstname birthyear municipality GdeNr_CNTR origin job partyname random_rank

label var canton "Kanton"
label var sex "Geschlecht"
label var name "Name"
label var firstname "Vorname"
label var birthyear "Geburtsjahr"
label var municipality "Wohngemeinde"
label var origin "Heimatort"
label var job "Beruf"
label var partyname "Parteiname"
labe var GdeNr_CNTR "GdeNr_CNTR"

* Krueger and Decet

local no_packages=10 // number of packages
local no_coders=2 // number of coders
local no_packages_percoder=`no_packages'/`no_coders' // number of packages per coder
local overlap=10 // overlap per package
local package_size=100 // size of package
local rest_package=`package_size'-`overlap' // size of package
local start_rest=`overlap'*`no_packages'/`no_coders' // size of package

di "`rest_package'"
di "`no_packages_percoder'"


local coder_no=0
forvalues i=1(1)`no_packages'{
if mod(`i',`no_packages_percoder')==1 & `i'>`no_packages_percoder'{
local coder_no=`coder_nr'+1
}
local overlap_no=`i'-`no_packages_percoder'*`coder_no'-1
di "`coder_no'"
di "`overlap_no'"

preserve
keep if inrange(random_rank,`start_rest'+(`i'-1)*`rest_package'+1,`start_rest'+(`i')*`rest_package') ///
	| inrange(random_rank,(`overlap_no')*`overlap'+1,(`overlap_no'+1)*`overlap')
drop random_rank
sort ID
export excel  using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\01_Nationalrat_Random_Sample\Check_False_Negatives_Sugarcube`i'.xlsx",  replace  firstrow(varlabels)
restore
} 

br if random_rank<=5
br if inrange(random_rank,41,50)

* Peter Wallmüller

local no_packages=10 // number of packages
local no_coders=2 // number of coders
local no_packages_percoder=`no_packages'/`no_coders' // number of packages per coder
local overlap=10 // overlap per package
local package_size=100 // size of package
local rest_package=`package_size'-`overlap' // size of package
local start_rest=`overlap'*`no_packages'/`no_coders' // size of package

di "`rest_package'"
di "`no_packages_percoder'"


local coder_no=0
forvalues i=1(1)`no_packages'{
if mod(`i',`no_packages_percoder')==1 & `i'>`no_packages_percoder'{
local coder_no=`coder_nr'+1
}
local overlap_no=`i'-`no_packages_percoder'*`coder_no'-1
di "`coder_no'"
di "`overlap_no'"

preserve
keep if inrange(random_rank,1000+`start_rest'+(`i'-1)*`rest_package'+1,1000+`start_rest'+(`i')*`rest_package') ///
	| inrange(random_rank,(`overlap_no')*`overlap'+1,(`overlap_no'+1)*`overlap')
drop random_rank
gen Sicher=""
gen Unsicher=""
gen Kommentare=""
sort ID
export excel  using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\01_Nationalrat_Random_Sample\Check_False_Negatives_Sugarcube1`i'.xlsx",  replace  firstrow(varlabels)
restore
} 




******************************************************
* (M) RD ESTIMATION
******************************************************

use "$path\02_Processed_data\nationalraete_1931_2015.dta" , clear 

* (i) Sharp RD 

rdrobust  elected_F1 votemargin_rel, p(2)
rdrobust  elected_F1 votemargin_rel, p(3)
rdrobust  elected_F1 votemargin_rel, p(2) h(0.05)
rdrobust  elected_F1 votemargin_rel, p(3) h(0.05)
rdrobust  elected_F1 votemargin_rel, p(2) h(0.01)
rdrobust  elected_F1 votemargin_rel, p(3) h(0.01)


* (ii) Fuzzy RD 

bysort elected: tab  incumbent_F1 elected_F1

gen votemargin_rel_sq=votemargin_rel^2
reg incumbent_F1 elected votemargin_rel  votemargin_rel_sq if elected_F1!=. & inrange(votemargin_rel,-0.01,0.01)

rdrobust  elected_F1 votemargin_rel, p(2) fuzzy(incumbent_F1)
rdrobust  elected_F1 votemargin_rel, p(3) fuzzy(incumbent_F1)
rdrobust  elected_F1 votemargin_rel, p(2) fuzzy(incumbent_F1)  h(0.05)
rdrobust  elected_F1 votemargin_rel, p(3) fuzzy(incumbent_F1)  h(0.05)
rdrobust  elected_F1 votemargin_rel, p(2) fuzzy(incumbent_F1)  h(0.01)
rdrobust  elected_F1 votemargin_rel, p(3) fuzzy(incumbent_F1)  h(0.01)

* (iii) Information for newspaper article


rdplot elected_F1 votemargin_rel, p(2) fuzzy(incumbent_F1), if inrange(votemargin_rel,-0.01,0.01)

rdplot elected_F1 votemargin_rel, p(2) , if inrange(votemargin_rel,-0.001,0.001)


rdrobust  elected_F1 votemargin_rel, p(2) fuzzy(incumbent_F1)  h(0.01) 



rdrobust  elected_F1 votemargin_rel, p(2) h(0.01)
di e(tau_cl_l)
di e(tau_cl_r)

tab canton year if inrange(votemargin_rel,-0.01,0.01) & elected_F1!=.


ttest age, by(incumbent)
ttest incumbent_F1, by(incumbent)
ttest elected_F1, by(incumbent)



******************************************************
* (O) TOMAS TURNER-ZWINKELS OUT
******************************************************
/*
use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear

keep if ID!=""

sort canton ID year name firstname birthyear sex votes
egen ID_numeric=group(ID)

keep canton year ID_numeric name firstname birthyear sex elected
order canton year ID_numeric name firstname birthyear sex elected

keep if elected==1
keep if inrange(year,1951,2015) 

sort canton ID_numeric year

export exc "$path\OutTomasTurnerZwinkels.xlsx", replace first(variables)

*/

*---------------------------
* A) Dataset for Switzerland
*---------------------------

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear

drop ID_pers
egen ID_pers=group(ID)

* (i) Generate running variable, election in next election,
* 	   party variables, and cantonal dummy variable

xtset ID_pers ID_time
sort ID_pers ID_time

* (iii) Keep only relevant variables

keep ID_pers year job
 
order ID_pers year job
 
export delimited "$path\Luechinger_Schelker_Schmid_2023_Occupations.csv", replace delim(";")





