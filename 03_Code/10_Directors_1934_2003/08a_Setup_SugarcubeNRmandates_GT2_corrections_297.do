

use "$path\02_Processed_data\10_Directors_1934_2003\RL_NR-Sugarcube_1934-2003.dta", clear 

merge m:1 PID NRid year CID using "$path\02_Processed_data\10_Directors_1934_2003\RL_sugarcube_297_CID.dta" // these are 297 NRid-year-PID observations (including respective companies) that were missing in first attempt to compare with GT2
drop _merge

drop firmname // this is "old" normalization of firmnames in firmname-cleaning process.
gen firmname = ""
replace firmname = ustrtrim(cname_orig)
replace firmname = ustrlower(firmname)
replace firmname = usubinstr(firmname, "ü", "ue",.)
replace firmname = usubinstr(firmname, "ä", "ae",.)
replace firmname = usubinstr(firmname, "ö", "oe",.)
replace firmname = usubinstr(firmname, "é", "e",.)
replace firmname = usubinstr(firmname, "è", "e",.)
replace firmname = usubinstr(firmname, "ê", "e",.)
replace firmname = usubinstr(firmname, "ç", "c",.)
replace firmname = usubinstr(firmname, `"""', "",.)
replace firmname = usubinstr(firmname, "«", "",.)
replace firmname = usubinstr(firmname, "»", "",.)
replace firmname = usubinstr(firmname, ".", "",.)
replace firmname = usubinstr(firmname, ",", "",.)
replace firmname = usubinstr(firmname, ";", "",.)
replace firmname = usubinstr(firmname, "-", "",.)
replace firmname = usubinstr(firmname, ":", "",.)
replace firmname = usubinstr(firmname, "'", "",.)
replace firmname = usubinstr(firmname, "  ", "",.)

sort NRid year firmname
recast str1000 firmname

// restrict sample to NRs that are covered in NR-VR Verzeichnis
merge m:1 NRid using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\uniqueNR_ID_GT2.dta"
keep if _merge==3
drop _merge
drop if year < 1982
drop if year > 2003

replace firmname = "none" if firmname == ""

duplicates report NRid year PID CID  // no duplicates. BUT: we want to prepare the comparison of NRid's with firmnames from GT2 -- > look for duplicate firmnames (w/o PID/CID)
duplicates report NRid year firmname
duplicates tag NRid year firmname, gen(dupl)
duplicates drop NRid year firmname, force // looks as if these are genuine duplicates (also w.r.t. location, capital, residential municipality)
drop dupl 

keep if east_miss == 1  // just focus on the originally 297 missing observations

merge m:1 NRid year firmname using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NR_ID_GT2_tmp.dta"

// prepare for manual coding: GT2 subset
preserve
keep if _merge == 2
gen source = "GT2"
keep NRid canton year name_mand firstname_mand gdename_mand firmname Companyname_mand Companylocation_mand source
rename name_mand name
rename firstname_mand firstname 
rename gdename_mand gdename 
rename Companyname_mand firm_orig 
rename Companylocation_mand firmlocation
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_eastmiss.dta", replace
restore

// prepare for manual coding: RL subset
keep if _merge == 1
gen source = "RL"
keep NRid canton year name firstname gdename cname_orig firmname firmlocation firstname_bus name_bus gdename_bus source east_miss
rename name name_nr
rename firstname firstname_nr
rename gdename gdename_nr
rename name_bus name
rename firstname_bus firstname 
rename gdename_bus gdename
rename cname_orig firm_orig
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge1_eastmiss.dta", replace

// combine non-matched observations for manual evaluation & linkage
append using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_eastmiss.dta"

replace firm_orig = "none" if firmname == "none"

// Identify affected NRids
preserve
keep if east_miss == 1
duplicates drop NRid year, force
keep NRid year
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NRid_297.dta", replace
restore

merge m:1 NRid year using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\NRid_297.dta"
keep if _merge == 3
drop _merge


// create alternating ID-Block variable to distinguish different blocks of NR-IDs (even/odd)
egen help1 = group(year NRid)
gen NROddBlocks = 0 
replace NROddBlocks = 1 if mod(help1, 2) == 1  // odd IDs
drop help1

// create variable to link the sources (create this variable for GT2)
bysort year NRid source (firmname): gen linkIDs = _n
replace linkIDs = . if source == "RL"

gen out = .

keep NROddBlocks source year NRid out linkIDs firm_orig firmlocation firstname name gdename
order NROddBlocks source year NRid out linkIDs firm_orig firmlocation firstname name gdename

sort year NRid firm_orig source

/*/ export information for coding by manual coders  -- > deactivated as manual coding is completed
export excel using "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_coding297.xlsx", first(var) replace
*/
import excel "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_coding297_in.xlsx", firstrow clear
drop if out == "1" // "Butterzentrale Gossau"
save "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_coding297_in.dta", replace


erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge1_eastmiss.dta"
erase "$path\02_Processed_data\10_Directors_1934_2003\11_Elected_Mandates_GT2\RL_GT2_merge2_eastmiss.dta"



