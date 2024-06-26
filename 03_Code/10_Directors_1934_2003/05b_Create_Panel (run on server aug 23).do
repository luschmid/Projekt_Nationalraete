global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
global path_rl "C:\Schmidlu\Dropbox\Record Linkage"

use "$path/02_Processed_data/10_Directors_1934_2003/Sugarcube_PersPanelID_Firmnames_wide.dta", clear

split lastname, p("-") gen(lastname_split)

replace lastname_split2=lastname_split2+"-"+lastname_split3 if lastname_split3!=""
replace lastname_split2=lastname_split2+"-"+lastname_split4 if lastname_split4!=""
replace lastname_split2=lastname_split2+"-"+lastname_split5 if lastname_split5!=""

keep id_sug year_sug e_id_sug n_id_sug ID_dupl gde_options firstname lastname male gdename firmnames lastname_split1 lastname_split2

gen firstname_stub=substr(firstname,1,15)
gen lastname_stub=substr(lastname,1,15)
gen lastname_split2_stub=substr(lastname_split2,1,15)

export delim "$path/02_Processed_data/10_Directors_1934_2003/12_Panel_Round1/Sugarcube_Panel_start.csv", delim(",") replace
