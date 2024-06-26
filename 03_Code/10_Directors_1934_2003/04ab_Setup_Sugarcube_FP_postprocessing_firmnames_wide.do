
version 17
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\10_Directors_1934_2003\"

use "$path\Sugarcube_RLPostProcessing_Persons-Firmnames_230503.dta", clear

* last manual geo-codes
replace e_id_sug = 2667300 if id_sug == 106906 & year_sug == 1992
replace e_id_sug = 2613300 if id_sug == 112413 & year_sug == 1991
replace e_id_sug = 2630000 if id_sug == 158857 & year_sug == 1989
replace e_id_sug = 2638300 if id_sug == 2436 & year_sug == 1995
replace e_id_sug = 2654100 if id_sug == 3948 & year_sug == 1986
replace e_id_sug = 2685800 if id_sug == 45255 & year_sug == 1969
replace e_id_sug = 2675700 if id_sug == 56200 & year_sug == 1992
replace e_id_sug = 2683100 if id_sug == 77810 & year_sug == 1995
replace e_id_sug = 2604900 if id_sug == 87936 & year_sug == 1988
replace e_id_sug = 2604900 if id_sug == 93820 & year_sug == 1989
replace e_id_sug = 2611300 if id_sug == 9682 & year_sug == 1969

replace n_id_sug = 1249100 if id_sug == 106906 & year_sug == 1992
replace n_id_sug = 1258800 if id_sug == 112413 & year_sug == 1991
replace n_id_sug = 1236100 if id_sug == 158857 & year_sug == 1989
replace n_id_sug = 1237700 if id_sug == 2436 & year_sug == 1995
replace n_id_sug = 1259400 if id_sug == 3948 & year_sug == 1986
replace n_id_sug = 1243800 if id_sug == 45255 & year_sug == 1969
replace n_id_sug = 1209500 if id_sug == 56200 & year_sug == 1992
replace n_id_sug = 1247100 if id_sug == 77810 & year_sug == 1995
replace n_id_sug = 1231000 if id_sug == 87936 & year_sug == 1988
replace n_id_sug = 1231000 if id_sug == 93820 & year_sug == 1989
replace n_id_sug = 1267600 if id_sug == 9682 & year_sug == 1969

replace gdename = "kuenten" if id_sug == 106906 & year_sug == 1992
replace gdename = "dornach" if id_sug == 112413 & year_sug == 1991
replace gdename = "fulenbach" if id_sug == 158857 & year_sug == 1989
replace gdename = "zofingen" if id_sug == 2436 & year_sug == 1995
replace gdename = "unterboezberg" if id_sug == 3948 & year_sug == 1986
replace gdename = "zollikon" if id_sug == 45255 & year_sug == 1969
replace gdename = "weggis" if id_sug == 56200 & year_sug == 1992
replace gdename = "zuerich" if id_sug == 77810 & year_sug == 1995
replace gdename = "oberdorf" if id_sug == 87936 & year_sug == 1988
replace gdename = "oberdorf" if id_sug == 93820 & year_sug == 1989
replace gdename = "basel" if id_sug == 9682 & year_sug == 1969

save "$path\Sugarcube_RLPostProcessing_Persons-Firmnames_orig_geo.dta"