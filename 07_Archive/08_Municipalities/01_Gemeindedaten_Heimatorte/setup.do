global path "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\Gemeindedaten Heimatorte\"

* Read in and append the municipality data for 1960 (first year available) and 
* all subsequent election years

import exc using "${path}\Gemeindestand 1960.xls", sh("Daten") first clear
g Jahr = 1931
g str10 Gemeindestand = "01.01.1960"
save "${path}\Gemeindestand.dta", replace

forv yr = 1935(4)2015 {
if `yr' < 1960 {
import exc using "${path}\Gemeindestand 1960.xls", sh("Daten") first clear
g Jahr = `yr'
g str10 Gemeindestand = "01.01.1960"
}
else {
import exc using "${path}\Gemeindestand `yr'.xls", sh("Daten") first clear
g Jahr = `yr'
g str10 Gemeindestand = "01.01.`yr'"
}
append using "${path}\Gemeindestand.dta"
sort Jahr BFS BFSGdenummer
save "${path}\Gemeindestand.dta", replace
}
