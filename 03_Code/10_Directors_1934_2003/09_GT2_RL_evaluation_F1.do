

*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "C:\Dropbox\Projekt Nationalräte"
*global path "C:\Current\Dropbox\Projekt Nationalräte"
global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"
*global path "C:\Users\Lukas\Dropbox\Projekt Nationalräte"


use "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2.dta", clear
tab Source

preserve
collapse (sum) TP FP FN
gen precision = (TP / (TP + FP))
gen recall = (TP / (TP + FN))
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1
restore

* Metrics per Year
preserve
collapse (sum) TP FP FN, by(Year)
gen precision = (TP / (TP + FP))
gen recall = (TP / (TP + FN))
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
*export excel using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\F1_perYear.xlsx", first(var) replace
restore

* Metrics per Gender
preserve 
use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
duplicates drop ID, force
keep ID sex
rename ID NRid
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\ID_sex.dta", replace
restore

preserve
merge m:1 NRid using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\ID_sex.dta"
keep if _merge == 3
collapse (sum) TP FP FN, by(sex)
gen precision = (TP / (TP + FP))
gen recall = (TP / (TP + FN))
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
*export excel using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\F1_perSex.xlsx", first(var) replace
restore



* Exclude some common names
preserve 
use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
duplicates drop ID, force
keep ID name
rename name nameNR
rename ID NRid
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\ID_Names.dta", replace
restore

preserve
merge m:1 NRid using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\ID_Names.dta"
keep if _merge == 3
drop _merge

drop if inlist(nameNR, "Müller", "Schmid", "Meier", "Keller", "Huber", "Weber", "Meyer", "Fischer")
drop if inlist(nameNR, "Graf", "Frey", "Schneider", "Brunner", "Widmer", "Kunz", "Baumann")
collapse (sum) TP FP FN
gen precision = (TP / (TP + FP))
gen recall = (TP / (TP + FN))
gen F1 = (2 * (TP / (TP + FP)) * (TP / (TP + FN))) / ((TP / (TP + FP)) + (TP / (TP + FN)))
sum precision recall F1
restore