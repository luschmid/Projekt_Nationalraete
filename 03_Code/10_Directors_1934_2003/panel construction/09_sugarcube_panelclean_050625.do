clear all

*global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"
*global path "C:\Users\schelkem\Dropbox\Projekt Nationalräte\"
global path "C:\Users\schmutzy\Dropbox\Projekt Nationalräte\" /* Office */
global dump "C:\Users\schmutzy\Desktop\Projekt Nationalräte\" /* Office Dumpster */
*global path "C:\Users\yanni\Dropbox\Projekt Nationalräte\" /* Home */
*global dump "C:\Users\yanni\Desktop\Projekt Nationalräte\" /* Home Dumpster */

********************************************************************************
*Loading final panel data at the end of "sugarcube_gaps_281124.do"
use "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\firm_panel_1934_2003_finale_3.dta", clear
* Finale panel / raw (rematch with initial dataset to have clean names, etc.)

keep CID UCID_1934_2003 year cname_orig gdenr gdenr_2012 gdenr_2018 gdename gdename_orig owners periode si_dummy_2 CID_str capital function signature page caddress_orig geo_merge ctn E_CNTR N_CNTR nb_owners
rename cname_orig cname_orig_temp

merge 1:1 CID year using "$dump\02_Processed_data\10_Directors_1934_2003\comp_geo-merge_complete.dta", keep(3)
replace cname_orig_temp = cname_orig if missing(cname_orig_temp)

drop cname cname_orig branch_dummy gdename_extract capital_extract cname_extract _merge
rename cname_orig_temp cname_orig

order CID CID_str UCID_1934_2003 year periode cname_orig gdenr gdenr_2012 gdenr_2018 gdename gdename_orig ctn caddress_orig E_CNTR N_CNTR owners nb_owners capital function signature page

rename capital capital_temp

merge 1:1 CID year using "$path\02_Processed_data\10_Directors_1934_2003\CID_nomcap.dta", keep(3)
drop capital_temp

order CID CID_str UCID_1934_2003 year periode cname_orig gdenr gdenr_2012 gdenr_2018 gdename gdename_orig ctn caddress_orig E_CNTR N_CNTR owners nb_owners capital capital_orig function signature page

drop _merge

drop gdenr_2012 gdenr_2018

replace gdenr = 261 if gdenr == 253 & year > 1989

merge m:1 gdenr year using "C:\Users\schmutzy\Desktop\Projekt Nationalräte\02_Processed_data\08_Municipalities\Gde_1931_2003_2012_2018"

sort gdenr year

drop if _merge == 2

********************************************************************************
* Some mistakes show up

replace gdenr = . if _merge == 1
drop _merge

********************************************************************************
save "$dump\02_Processed_data\10_Directors_1934_2003\23_Panel\Firm_Panel_1934_2003.dta", replace
save "$path\02_Processed_data\10_Directors_1934_2003\Firm_Panel_1934_2003.dta", replace