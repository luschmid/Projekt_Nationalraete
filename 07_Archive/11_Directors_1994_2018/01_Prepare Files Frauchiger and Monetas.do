
******************************************************
* (0) Set directories
******************************************************

	clear
	set more off

	capture cd "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
		if _rc==0{
		global datapath "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data"
		}
		

		
********************************************************
* (1) Prepare file for Monetas
********************************************************		
		
* (i) Check missings in terms of birthyear and municipality
	
cd "$datapath"
use "NRW_KANDIDATEN_ALL including 2015.dta", clear

tab year if gemeindename==""
tab year if geburtsjahr==.
	
br  vorname  nachname geburtsjahr gemeindename if gemeindename==""
br  vorname  nachname geburtsjahr gemeindename if geburtsjahr==.


	
* (ii) Drop duplicates in terms of br  vorname  nachname geburtsjahr gemeindename if gemeindename==""

duplicates tag vorname  nachname geburtsjahr gemeindename, gen(dups)

tab dups

sort vorname  nachname geburtsjahr gemeindename year

br  year vorname  nachname geburtsjahr gemeindename dups if dups==1
br year vorname  nachname geburtsjahr gemeindename if vorname=="Adalbert" &	nachname=="Durrer"

drop if dups>0

* (iii) Save out-files

keep vorname  nachname geburtsjahr gemeindename
export excel NRW_out_Monetas.xlsx, replace
		
		
********************************************************
* (2) Prepare file for Monetas
********************************************************		

cd "$datapath"
use "NRW_KANDIDATEN_ALL including 2015.dta", clear


* (ii) Prepare data on vote lists

preserve
forvalues i=1975(4)2015{
clear
import excel using "$data\0_original_data\NRW 1971_2015ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet("`i'") cellrange(A5:L166)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes`i'.dta", replace
}
restore

save "$data\1_data\partyvotes`i'.dta", replace


use "$data\1_data\partyvotes1971.dta", clear
forvalues i=1975(4)2011{
append using "$data\1_data\partyvotes`i'.dta"
erase "$data\1_data\partyvotes`i'.dta"
}


gen Nationalratswahl=Wahljahr
save "$data\1_data\partyvotes_all.dta", replace
restore
 


