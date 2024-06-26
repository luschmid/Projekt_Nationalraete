*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path_rl "C:\Schmidlu\Dropbox\Record Linkage"

global path "C:\Current\Dropbox\Projekt Nationalräte"
global path_rl "C:\Current\Dropbox\Record Linkage"

***********************************
* Search for gaps in Sugarcube data
***********************************

* (i) Create panel 

* (a) Get unique politicians

use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear

duplicates drop ID, force
keep ID
rename ID id_0

save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\unique_pol_ids.dta", replace


clear 
set obs 36
gen time_period=_n

* (b) Cross unique politician data 

cross using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\unique_pol_ids.dta" 

gen year_1=.
replace year_1= 1933 if time_period==1 
replace year_1= 1942 if time_period==2 
replace year_1= 1959 if time_period==3 
replace year_1= 1961 if time_period==4 
replace year_1= 1962 if time_period==5 
replace year_1= 1963 if time_period==6 
replace year_1= 1965 if time_period==7 
replace year_1= 1966 if time_period==8 
replace year_1= 1969 if time_period==9 
replace year_1= 1972 if time_period==10
replace year_1= 1975 if time_period==11
replace year_1= 1979 if time_period==12
replace year_1= 1980 if time_period==13
replace year_1= 1981 if time_period==14
replace year_1= 1982 if time_period==15
replace year_1= 1983 if time_period==16
replace year_1= 1984 if time_period==17
replace year_1= 1985 if time_period==18
replace year_1= 1986 if time_period==19
replace year_1= 1987 if time_period==20
replace year_1= 1988 if time_period==21
replace year_1= 1989 if time_period==22
replace year_1= 1990 if time_period==23
replace year_1= 1991 if time_period==24
replace year_1= 1992 if time_period==25
replace year_1= 1993 if time_period==26
replace year_1= 1994 if time_period==27
replace year_1= 1995 if time_period==28
replace year_1= 1996 if time_period==29
replace year_1= 1997 if time_period==30
replace year_1= 1998 if time_period==31
replace year_1= 1999 if time_period==32
replace year_1= 2000 if time_period==33
replace year_1= 2001 if time_period==34
replace year_1= 2002 if time_period==35
replace year_1= 2003 if time_period==36

* (ii) Merge correspondence table

merge 1:m id_0 year_1 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\11_RL_Output_Round1_In\RL_Output_Round1b_fp_maj_IDs.dta"

* (iii) Merge Sugarcube data 

rename e_cntr_w_1 e_id_sug
rename n_cntr_w_1 n_id_sug
rename year_1 year_sug
rename id_1 id_sug
rename e_cntr_w_0 e_id_polit
rename n_cntr_w_0 n_id_polit

destring id_sug, replace

gen ID_dupl=1

merge m:1 id_sug year_sug e_id_sug n_id_sug ID_dupl using "$path/02_Processed_data/10_Directors_1934_2003/Sugarcube_RLPostProcessing_Persons-Firmnames.dta", gen(merge_all)

drop if merge_all==2


rename lastname lastname_sug
rename firstname firstname_sug

* (iv) Merge politician data 

preserve
use "$path/02_Processed_data/02_Elections_1971_2015/nationalraete_1931_2015_Geo.dta",clear
*rename E_CNTR_w e_id_polit 
*rename N_CNTR_w  n_id_polit
keep ID name firstname e_id_polit n_id_polit
duplicates drop ID e_id_polit n_id_polit, force
save "$path/02_Processed_data/02_Elections_1971_2015/nationalraete_1931_2015_Geo_temp.dta", replace
restore

gen ID=id_0

merge m:1 ID e_id_polit n_id_polit using "$path/02_Processed_data/02_Elections_1971_2015/nationalraete_1931_2015_Geo_temp.dta", gen(merge2)
erase "$path/02_Processed_data/02_Elections_1971_2015/nationalraete_1931_2015_Geo_temp.dta"

drop if merge2==2

* (v) Drop duplicates 
* Note: remove two sources of municipality duplicates: 
* 1. same politician who lives in different municipalities and gets 
*    the same mandate in the same year (rapperswil and uznach case). 
* 2. same politician but not clear municipality in sugarcube data
*    (rüti and buchs)


duplicates drop id_0 id_sug year_sug, force

* (vi) Keep only relevant variables and sort

keep id_0 id_sug year_sug firstname_sug firstname lastname_sug name firmnames time_period

sort id_0 year_sug id_sug
order id_0 id_sug year_sug firstname* lastname* name firmnames


* (vii) Sort out those with only 0 or 1 entries 

gen entry=0
replace entry=1 if id_sug!=.

bysort id_0: egen no_entries=sum(entry)
drop if no_entries==1 | no_entries==0

* (viii) Sort out those with two subsequent entries 

gen time_period_temp=time_period
replace time_period_temp=.  if entry==0

bysort id_0: egen time_period_min=min(time_period_temp)
bysort id_0: egen time_period_max=max(time_period_temp)

drop if no_entries==2 & time_period_min + 1 == time_period_max

* (ix) Keep only relevant variables

egen id_0_num=group(id_0)

keep id_0_num year_sug firstname name id_sug id_0  time_period firstname_sug lastname_sug  firmnames 
order id_0_num year_sug firstname name id_sug id_0  time_period firstname_sug lastname_sug  firmnames  


* (x) Create random id and outfiles in batches

* (a) Create random id and merge it to initial dataset

preserve
duplicates drop id_0, force
keep id_0
set seed 1234
generate random = runiform()
sort random 
generate id_new = _n 
save "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\12_RL_Output_Round_2_Out\RL_Output_Gaps_randomsample.dta",replace
restore

merge m:1 id_0 using "$path\02_Processed_data\12_Record_linkage\02_Sugarcube\12_RL_Output_Round_2_Out\RL_Output_Gaps_randomsample.dta", gen(merge_rs)


* (b) Outfile in batches of 250 ids

sort id_0 time_period

local batch_size 250
forvalues i = 1(1)27 {
preserve
keep if inrange(id_0_num,(`i'-1)*`batch_size'+1,`i'*`batch_size')
export excel using 	"$path\02_Processed_data\12_Record_linkage\02_Sugarcube\12_RL_Output_Round_2_Out\RL_Output_Gaps_`i'.xlsx", first(var) replace
restore
}



