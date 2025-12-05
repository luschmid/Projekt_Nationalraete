
clear
set more off

global path_nr "D:\SchmidLu\Dropbox\Projekt Nationalräte"
global path "D:\SchmidLu\Dropbox\Record Linkage"
*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"


* (i) Read in data

import delimited "$path\Data\experiments\experiment3\found_matches_experiment_3_full_data.csv", ///
encoding("utf-8") delimiters(",") clear
drop if id_polit==""
* Note: drop observations with no id_polit

* (ii) Drop all observations with missing geo data

drop if e_cntr_w_bis==. & n_cntr_w_bis==. & e_cntr_b_bis==. & n_cntr_b_bis ==. 


* (iii) Check whether there is variation in id_bis per id_polit
*count if (e_cntr_w_bis==. & n_cntr_w_bis==.) |  (e_cntr_b_bis==. & n_cntr_b_bis ==.) 

bysort id_polit: egen id_bis_sd=sd(id_bis)





* (iii) Tag observations with the same birthyear (plausibility check 1)

gen equalbirthyear=0
replace equalbirthyear=1 if birthyear_polit==birthyear_bis

* (iii) Tag observations with the same sex (plausibility check 2)

replace sex_bis="F" if sex_bis=="W"

gen equalsex=0
replace equalsex=1 if sex_bis==sex_polit

* (iv)



collapse (first) name_polit firstname_polit birthyear_polit sex_polit e_cntr_w_polit ///
	n_cntr_w_polit e_cntr_b_polit n_cntr_b_polit id_bis name_bis firstname_bis ///
	birthyear_bis sex_bis e_cntr_w_bis n_cntr_w_bis e_cntr_b_bis n_cntr_b_bis ///
, by(id_polit)


export delimited "$path\Data\experiments\experiment3\found_matches_experiment_3_full_data_unique_samebirthyearandsex.csv", replace


*********************
* Recursive matching
*********************

use "$path_nr\02_Processed_data\nationalraete_1931_2015.dta", clear
rename ID id_polit
merge m:1 id_polit using "$path\Data\experiments\experiment3\found_matches_experiment_3_full_data_unique.dta", gen(merge_nr)
