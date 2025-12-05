
*** Read-in Bisnode data ***

global datapath "E:\12. Cloud\Dropbox\Projekt Nationalräte\2 Data Collection\Verwaltungsratsmandate\Bisnode"

set more off
cd "$datapath\"
import excel  "Uni Fribourg Output 14-08-2017.xlsx", sheet(Matching Daten) clear first
destring PID_current, replace
sort PID_current
duplicates report PID_current
save "$datapath\Bisnode_matchingdata_temp.dta", replace
import excel  "Uni Fribourg Output 14-08-2017.xlsx", sheet(Personen - Mandate) clear first
destring PID_current, replace
sort PID_current
duplicates report PID_current
save "$datapath\Bisnode_person-mandate_temp.dta", replace
import excel  "Uni Fribourg Output 14-08-2017.xlsx", sheet(Firmendaten) cellrange(A2:CS25539) clear first
destring DUNS, replace
sort DUNS
duplicates report DUNS
save "$datapath\Bisnode_firms_temp.dta", replace

use "$datapath\Bisnode_matchingdata_temp.dta", clear
merge m:m PID_current using "$datapath\Bisnode_person-mandate_temp.dta"
drop _merge
sort DUNS
merge m:1 DUNS using "$datapath\Bisnode_firms_temp.dta"
drop _merge
save "$datapath\Bisnode.dta", replace
erase "$datapath\Bisnode_matchingdata_temp.dta"
erase "$datapath\Bisnode_person-mandate_temp.dta"
erase "$datapath\Bisnode_firms_temp.dta"

