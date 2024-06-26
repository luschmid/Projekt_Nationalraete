version 16
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"

*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"

******************************************************************************************
*** Import corrected first and last name files from raw-data corrections (third round) ***
******************************************************************************************


forv i = 1(1)15 { 
	import excel using "$path\02_Processed_data\10_Directors_1934_2003\06_Round3_Infiles\Original_check_`i'_out.xlsx", clear firstrow
	save "$path\02_Processed_data\10_Directors_1934_2003\06_Round3_Infiles\Original_check_`i'_out.dta", replace
}
use "$path\02_Processed_data\10_Directors_1934_2003\06_Round3_Infiles\Original_check_1_out.dta", clear
forv i = 2(1)15 {
	append using "$path\02_Processed_data\10_Directors_1934_2003\06_Round3_Infiles\Original_check_`i'_out.dta", force
	save "$path\02_Processed_data\10_Directors_1934_2003\06_Round3_Infiles\Original_check_corrections.dta", replace
}

rename firstname firstname_ch3
rename lastname lastname_ch3

sort tempID
merge 1:m tempID using "$path\02_Processed_data\10_Directors_1934_2003\06_Round3_Infiles\tempID_PID_year_ID_dupl.dta", gen(mergeIDs)  // import original ID's from original data
drop if mergeIDs == 2
drop mergeIDs

save "$path\02_Processed_data\10_Directors_1934_2003\06_Round3_Infiles\Original_check_corrections.dta", replace


/*
-- > all of that in setup do-file 


use "$path\02_Processed_data\10_Directors_1934_2003\05_Round3_Outfiles\Sugarcube_Person_RevName-Geo.dta", clear
sort tempID

merge m:1 tempID using "$path\02_Processed_data\10_Directors_1934_2003\06_Round3_Infiles\Original_check_corrections.dta", gen(mergenamefinal)

/*
   Result                           # of obs.
    -----------------------------------------
    not matched                     4,585,478
        from master                 4,585,478  (_merge==1)
        from using                          0  (_merge==2)

    matched                            22,583  (_merge==3)
    -----------------------------------------
all of the observations given for manual inspection (22,583) are merged back in again
*/

drop G // two cells containing only a space

replace firstname = finalfirstname if mergenamefinal == 3 & finalfirstname != ""
replace lastname = finallastname if mergenamefinal == 3 & finallastname != ""

// indicator for sort of name correction

gen namecorrection = 0
replace namecorrection = 1 if firstname_orig != firstname & lastname_orig == lastname // firstname corrected
replace namecorrection = 2 if lastname_orig != lastname & firstname_orig == firstname // lastname corrected
replace namecorrection = 3 if firstname_orig != firstname & lastname_orig != lastname // firstname and lastname corrected
tab namecorrection


drop tempID newfirstname_f newlastname_f indfirst mergefirst newlastname_l newfirstname_l indlast mergelast cases corrections manualcheck lastname_ch3 firstname_ch3 finallastname finalfirstname mergenamefinal

label var firstname_orig "First name as coded originally by sugarcube"
label var lastname_orig "Last name as coded originally by sugarcube"
label var namecorrection "Name corrections: 0: no correction; 1: first only; 2: lastname only; 3: first & last"

save "$path\02_Processed_data\10_Directors_1934_2003\Sugarcube_Person_CleanName-Geo.dta", replace

forv i = 1(1)15 {
erase "$path\02_Processed_data\10_Directors_1934_2003\06_Round3_Infiles\Original_check_`i'_out.dta"
}
