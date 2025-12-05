
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
global serverpath "\\Srw-hpc5\e\Projekt Nationalräte"

import delimited "$path\02_Processed_data\13_Running_variable\kotakorpietal2017_honduras.csv", clear
duplicates drop id_stata year, force
rename id_stata id_Stata
rename year Year
save "$path\02_Processed_data\13_Running_variable\kotakorpietal2017_honduras.dta", replace

use "$path\02_Processed_data\15_Elections_Honduras\elections_hn_final.dta", clear
merge 1:1 id_Stata Year using "$path\02_Processed_data\13_Running_variable\kotakorpietal2017_honduras.dta", gen(kotakorpietal2017_honduras)

bysort Elected: sum p

bysort Department Year party_num: egen help1=min(p) if Elected==1
bysort Department Year party_num: egen help2=max(p) if Elected==0
bysort Department Year party_num: egen p_min=max(help1) 
bysort Department Year party_num: egen p_max=max(help2) 
gen p_diff=p_min-p_max
tab p_diff
unique Department Year party_num if p_diff<0


* Result: For 14 parties and 184 cases, the smallest p of the elected candidates
* on a party list is smaller than the highest p of the non-elected candidates on the same list. 
