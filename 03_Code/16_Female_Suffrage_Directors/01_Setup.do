
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte\"
global gb_path "C:\Schmidlu\Dropbox\Gender_boards_politics\"


use "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person_CleanName-Gender-Geo.dta", clear

* (i) Drop observations with missing Geo code and/or missing gender information

drop if missing(gdenr_2018)
drop if missing(male)

// temp: assign Allmendingen and Ittingen to the canton of Bern

replace ctn="be" if gdenr_2018==362 | gdenr_2018==630


* (ii) Get number of companies per entry

gen nbrcompIDs = length(companiesid) - length(subinstr(companiesid , ",", "", .)) + 1

* (iii) Collapse dataset

gen ones=1

collapse (sum) ones nbrcompIDs , by(ctn year male)

* (iv) Reshape to wide format

reshape wide ones nbrcompIDs, i(ctn year) j(male)

rename ones0 fem_dir
rename ones1 male_dir
rename nbrcompIDs0 fem_man 
rename nbrcompIDs1 male_man

gen fem_share_dir=fem_dir/(fem_dir+male_dir)
gen fem_share_man=fem_man/(fem_man+male_man)

* (v) Merge control variables (from regulation project by Simon and Mark)

merge 1:1 ctn year using "$path\02_Processed_data\19_Cantonal_Covariates\ctn_cov_regdata.dta"

/*
 Result                      Number of obs
    -----------------------------------------
    Not matched                         2,039
        from master                         0  (_merge==1)
        from using                      2,039  (_merge==2)

    Matched                               925  (_merge==3)
    -----------------------------------------
*/
* Note: JU before 1979 missing, some gap years missing

*drop if _merge==2
drop _merge

* (vi) Label var and save dataset 

label var fem_share_dir "Directors: share females"
label var fem_share_man "Mandates: share females"
label var fem_dir "Directors: female"
label var male_dir "Directors: male"
label var fem_man "Mandates: female"
label var male_man "Mandates: male"

save "$gb_path\data\CH\cantonal\genderboards_ctn_ch.dta", replace
