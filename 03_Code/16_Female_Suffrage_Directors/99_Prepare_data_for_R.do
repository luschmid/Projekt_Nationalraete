
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte\"

use "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person_CleanName-Gender-Geo.dta" 

egen ID_new=group(PID year firstname lastname)
replace ID_new= 5000000+_n if ID_new==. // 427 obs
sum ID_new

keep PID year firstname lastname male companiesid language_w /// 
GdeNr_E_CNTR GdeNr_N_CNTR PLZ4 city language_w ctn functions

split companiesid, parse(",")


outsheet using Sugarcube_Person_CleanName-Gender-Geo.csv , delimiter(";") nolabel