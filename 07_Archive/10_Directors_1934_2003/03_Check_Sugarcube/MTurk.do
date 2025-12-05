
******************************************************
* (0) Set directories
******************************************************

cd "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\8 MTurk"


* (i) Start with AG
	

import excel  "Gesamtdatei_Ausreisser_zusammengetragen_180309.xlsx"	, sheet(Tabelle1) clear first


tab RICHTIGFALSCH

gen persons_all=persons
replace persons_all=persons_corrected if  RICHTIGFALSCH=="FALSCH " |  RICHTIGFALSCH=="0"
replace persons_all=persons_corrected if  persons_corrected!=.

gen ones=1
bysort ones: egen persons_total=sum(persons_all)
sum persons_total
