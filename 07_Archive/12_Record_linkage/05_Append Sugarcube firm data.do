******************************************************
* (A) Set directories
******************************************************

clear
set more off

*global path "E:\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "C:\Dropbox\Projekt Nationalräte"
*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"
*global path "C:\Users\Lukas\Dropbox\Projekt Nationalräte"

use "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_34_63.dta", clear

append using "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_64_70.dta"

append using "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_71_98.dta"

append using "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_99_00.dta"

append using "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_01.dta"

append using "$path\02_Processed_data\10_Directors_1934_2003\temp_allyears_02_03.dta"

keep PID year companies
export delimited using ///
"$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Companies.csv", delimiter(";") replace