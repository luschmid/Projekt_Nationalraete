******************************************************
* (A) Set directories
******************************************************

clear
set more off

global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "C:\Dropbox\Projekt Nationalräte"
*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"
*global path "C:\Users\Lukas\Dropbox\Projekt Nationalräte"
global outpath "$path\02_Processed_data\15_Elections_Honduras"

******************************************************
* (B) SWITZERLAND
******************************************************

*(i) Number of lists per canton and year

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
gen ones=1
collapse ones, by(canton year liste_norm)
collapse (sum) ones, by(canton year)

sum ones

******************************************************
* (C) HONDURAS
******************************************************

*(i) Number of lists per department and year

use "$outpath\elections_hn_final.dta", replace

gen ones=1
collapse (first) ones (sum) Elected, by(Department Year Party)
collapse (sum) ones Elected, by(Department Year )

sum ones Elected	
