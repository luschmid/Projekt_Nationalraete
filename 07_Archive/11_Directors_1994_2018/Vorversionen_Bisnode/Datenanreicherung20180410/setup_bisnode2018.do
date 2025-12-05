clear
cap log close
set more 1
version 14

global bisnode "E:\12. Cloud\Dropbox\Projekt Nationalräte\2 Data Collection\Verwaltungsratsmandate\Bisnode\"

log using "${bisnode}setup_bisnode2018.smcl", replace

************************
**** Bisnode data  *****
************************


/*** Data import 18.03.2018 ***
set excelxlsxlargefile on
import excel "${bisnode}Datenanreicherung20180329\Uni Fribourg Output 29-03-2018.xlsx", sheet(Input Uni Fribourg) cellrange(A1:I36841) firstrow clear
save "${bisnode}Datenanreicherung20180329\bisnode_InputUniFR2018.dta", replace
clear
import excel "${bisnode}Datenanreicherung20180329\Uni Fribourg Output 29-03-2018.xlsx", sheet(Matching Daten) cellrange(A1:AG36841) firstrow clear
destring PLZ Zip_matched PID_matched MPID_matched PID_current Zip_current, replace
rename PLZ PLZ_Wohnort
save "${bisnode}Datenanreicherung20180329\bisnode_MatchingDaten2018.dta", replace
clear
import excel "${bisnode}Datenanreicherung20180329\Uni Fribourg Output 29-03-2018.xlsx", sheet(Personen - Mandate) cellrange(A1:I20583) firstrow clear
save "${bisnode}Datenanreicherung20180329\bisnode_Personen-Mandate2018.dta", replace
clear
import excel "${bisnode}Datenanreicherung20180329\Uni Fribourg Output 29-03-2018.xlsx", sheet(Firmendaten) cellrange(A2:CS51979) firstrow clear
destring DUNS PLZDomiziladresse PLZPostfachadresse NOGACode, replace
save "${bisnode}Datenanreicherung20180329\bisnode_Firmendaten2018.dta", replace
clear
import excel "${bisnode}Datenanreicherung20180329\Uni Fribourg Output 29-03-2018.xlsx", sheet(Personen - Mandate zusätzlich) cellrange(A1:AB118333) firstrow clear
destring DUNS, replace
rename N firstname
rename O name
rename T residence
save "${bisnode}Datenanreicherung20180329\bisnode_Personen-Mandate-zus2018.dta", replace
*/

/*** Data import 10.04.2018 ***
set excelxlsxlargefile on
import excel "${bisnode}Datenanreicherung20180410\Uni Fribourg Output 10-04-2018.xlsx", sheet(Input Uni Fribourg) cellrange(A1:I36841) firstrow clear
save "${bisnode}Datenanreicherung20180410\bisnode_InputUniFR2018.dta", replace
clear
import excel "${bisnode}Datenanreicherung20180410\Uni Fribourg Output 10-04-2018.xlsx", sheet(Matching Daten) cellrange(A1:AG36841) firstrow clear
destring PLZ Zip_matched PID_matched MPID_matched PID_current Zip_current, replace
rename PLZ PLZ_Wohnort
save "${bisnode}Datenanreicherung20180410\bisnode_MatchingDaten2018.dta", replace
clear
import excel "${bisnode}Datenanreicherung20180410\Uni Fribourg Output 10-04-2018.xlsx", sheet(Personen - Mandate) cellrange(A1:I20583) firstrow clear
save "${bisnode}Datenanreicherung20180410\bisnode_Personen-Mandate2018.dta", replace
clear
import excel "${bisnode}Datenanreicherung20180410\Uni Fribourg Output 10-04-2018.xlsx", sheet(Firmendaten) cellrange(A2:CS51979) firstrow clear
destring DUNS PLZDomiziladresse PLZPostfachadresse NOGACode Gründungsjahr, replace
save "${bisnode}Datenanreicherung20180410\bisnode_Firmendaten2018.dta", replace
clear
import excel "${bisnode}Datenanreicherung20180410\Uni Fribourg Output 10-04-2018.xlsx", sheet(Personen - Mandate zusätzlich) cellrange(A1:AI118327) firstrow clear
destring DUNS Zip_current, replace
rename U firstname
rename V name
rename AA residence
save "${bisnode}Datenanreicherung20180410\bisnode_Personen-Mandate-zus2018.dta", replace
*/

global bisnode "E:\12. Cloud\Dropbox\Projekt Nationalräte\2 Data Collection\Verwaltungsratsmandate\Bisnode\Datenanreicherung20180410\"

*** Data comparisons ***
* Matching Daten
use "${bisnode}bisnode_MatchingDaten2018.dta", clear
order ID PID_current PID_matched MPID_matched First_name_matched Vorname Last_name_matched Nachname
sort PID_current
codebook PID_current 
* 27'010 distinct observations; 36'840 observations in "Matching Daten"  

* Personen - Mandate 
use "${bisnode}bisnode_Personen-Mandate2018.dta", clear
sort PID_current
codebook PID_current DUNS
* PID_current 3'382 distinct observations; 20'582 observations in "Personen - Mandate"  
* DUNS 12'490 distinct observations  

* Personen - Mandate zusätzlich 
use "${bisnode}bisnode_Personen-Mandate-zus2018.dta", clear
sort PersonenID
codebook PID_current PersonenID DUNS
* PID_current 10'615 distinct observations; 118'326 observations; 49,605 missing observations in "Personen - Mandate zusätzlich"  
* PersonenID 27'410 distinct observations; 118'326 observations in "Personen - Mandate zusätzlich"  
* DUNS 49'594 distinct observations

* Firmendaten 
use "${bisnode}bisnode_Firmendaten2018.dta", clear
sort DUNS
codebook DUNS
* DUNS 51'977 distinct observations


* Compare "Matching Daten" with "Personen - Mandate zusätzlich"
//PID: "Matching Daten" contains 27'010 unique PID's, "Personen - Mandate zusätzlich" contains 27'410 unique PID's
use "${bisnode}bisnode_MatchingDaten2018.dta", clear
duplicates drop PID_current, force
keep PID_current
sort PID_current
save "${bisnode}bisnode_MatchingDaten2018_temp.dta", replace
use "${bisnode}bisnode_Personen-Mandate-zus2018.dta", clear
duplicates drop PID_current, force
keep PID_current
sort PID_current
save "${bisnode}bisnode_Personen-Mandate-zus2018_temp.dta", replace
merge 1:1 PID_current using "${bisnode}bisnode_MatchingDaten2018_temp.dta"
tab _merge

/*
_merge			Freq.	Percent	Cum.	
master only (1)	6,744	19.98	19.98		--> only "Personen - Mandate zusätzlich"
using only (2)	23,139	68.55	88.53		--> only "Matching Daten"
matched (3)		3,872	11.47	100.00			
Total			33,755	100.00
*/

erase "${bisnode}bisnode_MatchingDaten2018_temp.dta"
erase "${bisnode}bisnode_Personen-Mandate-zus2018_temp.dta"

//DUNS: "bisnode_Personen - Mandate" contains 12'490 unique DUNS's, "Personen - Mandate zusätzlich" contains 49'594 unique DUNS's
use "${bisnode}bisnode_Personen-Mandate2018.dta", clear
duplicates drop DUNS, force
keep DUNS
sort DUNS
save "${bisnode}bisnode_Personen-Mandate2018_temp.dta", replace
use "${bisnode}bisnode_Personen-Mandate-zus2018.dta", clear
duplicates drop DUNS, force
keep DUNS 
sort DUNS
save "${bisnode}bisnode_Personen-Mandate-zus2018_temp.dta", replace
use "${bisnode}bisnode_Firmendaten2018.dta", clear
duplicates drop DUNS, force
keep DUNS 
sort DUNS
save "${bisnode}bisnode_Firmendaten2018_temp.dta", replace

use "${bisnode}bisnode_Personen-Mandate-zus2018_temp.dta", clear
merge 1:1 DUNS using "${bisnode}bisnode_Personen-Mandate2018_temp.dta"
tab _merge

/*
_merge			Freq.	Percent	Cum.		
master only (1)	39,487	75.97	75.97		--> only "Personen - Mandate zusätzlich"
using only (2)	2,383	4.58	80.55		--> only "Matching Daten"
matched (3)		10,107	19.45	100.00		
Total			51,977	100.00
*/

rename _merge merge1
sort DUNS
merge 1:1 DUNS using "${bisnode}bisnode_Firmendaten2018_temp.dta"
tab _merge

/*
_merge		Freq.	Percent	Cum.		
matched (3)	51,977	100.00	100.00		
Total		51,977	100.00
*/

erase "${bisnode}bisnode_Personen-Mandate2018_temp.dta"
erase "${bisnode}bisnode_Personen-Mandate-zus2018_temp.dta"
erase "${bisnode}bisnode_Firmendaten2018_temp.dta"

cap log close

