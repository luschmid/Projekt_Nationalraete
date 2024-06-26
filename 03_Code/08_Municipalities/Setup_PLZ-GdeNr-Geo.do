
*******************************************
*** Merge GdeNr - PLZ - Geo-Coordinates ***
*******************************************

* Join GdeNr/PLZ with geo coordinates

* Read-in Geo/PLZ-Data (2018 cross-section: GdeNr & PLZ & GeoCoordinates 2018)
clear
import excel "$path/01_Raw_data/08_Municipalities/PerimeterGDE.xlsx", sheet("PLZO_CSV_LV03") firstrow clear  // data from the Amtliches Ortschaftenverzeichnis (Catastre.ch)

// include the missing Oberlangenegg BfS 0935, PLZ 3616, E 623006, N 183248
set obs `=_N+1'
count
local Nobs = r(N)
replace Ortschaftsname = "Oberlangenegg" in `Nobs'
replace PLZ = 3616 in `Nobs'
replace Zusatzziffer = 0 in `Nobs'
replace Gemeindename = "Oberlangenegg" in `Nobs'
replace BFSNr = 0935 in `Nobs'
replace Kantonskürzel = "BE" in `Nobs'
replace E = 623006 in `Nobs'
replace N = 183248 in `Nobs'
gen gdenr_2018 = BFSNr
sort gdenr_2018 PLZ Zusatzziffer
save "$path/02_Processed_data/08_Municipalities/PerimeterGDE-PLZ.dta", replace

* Read-in Geo/GdeNr BfS (2018 cross-section: GdeNr & GeoCoordinates)
clear 
import excel "$path/01_Raw_data/08_Municipalities/be-b-00.03-gg18.xls", sheet("g1g18") firstrow clear
rename GMDNR gdenr_2018
duplicates list, gdenr_2018
sort gdenr_2018
save "$path/02_Processed_data/08_Municipalities/PerimeterGDE-GdeNR.dta", replace

* Standardize Municipalities datset
use "$path/02_Processed_data/08_Municipalities/Gde_1960_2018", clear
order year gdename gdenr gdenr_2018 ctnnr ctn
sort year gdenr 
duplicates list year gdenr
duplicates drop year gdenr, force
save "$path/02_Processed_data/08_Municipalities/muni_temp1.dta", replace

* Generate GdeNr for years prior to 1960 based on 1960 and append to 1960-2018
drop if year!=1960
sort gdenr
drop year
save "$path/02_Processed_data/08_Municipalities/muni_temp2.dta", replace
use "$path/01_Raw_data/08_Municipalities/AddYears.dta", clear
cross using "$path/02_Processed_data/08_Municipalities/muni_temp2.dta"
sort year gdenr
save "$path/02_Processed_data/08_Municipalities/Gde_1931_1959.dta", replace
append using "$path/02_Processed_data/08_Municipalities/muni_temp1.dta"
save "$path/02_Processed_data/08_Municipalities/Gde_1931_2018.dta", replace
erase "$path/02_Processed_data/08_Municipalities/muni_temp1.dta"
erase "$path/02_Processed_data/08_Municipalities/muni_temp2.dta"

* Merge with Geo-Data at municipality level (based on GdeNr BfS)
use "$path/02_Processed_data/08_Municipalities/Gde_1931_2018", clear
sort gdenr_2018 gdenr year
merge m:1 gdenr_2018 using "$path/02_Processed_data/08_Municipalities/PerimeterGDE-GdeNR.dta"
// 3 obs only in using: Staatswald Galm, C'za Cadenazzo/Monteceneri, C'za Capriasca/Lugano (keine Gemeinden)
drop _merge
save "$path/02_Processed_data/08_Municipalities/GdeNr2018_Geo_1931-2018.dta", replace


* Merge with Geo-Data at PLZ level

// Municip.Data: same GdeNr_2018 for various GdeNr and various years (one GdeNr per year, but many GdeNr_2018 per year)
// Geo-PLZ Data: same GdeNr for various PLZ (PLZ only for 2018)

// Loop through years, keep only one municipality per year and merge PLZ on it
// Save yearly datasets
// Append all yearly data

use "$path/02_Processed_data/08_Municipalities/GdeNr2018_Geo_1931-2018.dta", clear
sort year gdenr_2018 gdenr
quietly by year gdenr_2018:  gen dup = cond(_N==1,0,_n)
sort year gdenr_2018 gdenr
save "$path/02_Processed_data/08_Municipalities/muni_1931_2018_temp", replace

forv y=1931(1)2018  {
	use "$path/02_Processed_data/08_Municipalities/muni_1931_2018_temp", clear
	keep if year == `y'
	sum dup,
	local max = r(max)
	save "$path/02_Processed_data/08_Municipalities/muni_`y'", replace
	forv i = 0(1)`max' {
		use "$path/02_Processed_data/08_Municipalities/muni_`y'", clear
		keep if dup == `i'
		sort gdenr_2018
		merge 1:m gdenr_2018 using "$path/02_Processed_data/08_Municipalities/PerimeterGDE-PLZ.dta"
		keep if _merge==3
		save "$path/02_Processed_data/08_Municipalities/muni_`y'`i'", replace	
	}
	use "$path/02_Processed_data/08_Municipalities/muni_`y'0", clear
	forv j = 1(1)`max' {
		append using "$path/02_Processed_data/08_Municipalities/muni_`y'`j'"
	}
	sort gdenr_2018
	save "$path/02_Processed_data/08_Municipalities/Muni_Geo-PLZ_`y'.dta", replace
	forv k = 0(1)`max' {
		erase "$path/02_Processed_data/08_Municipalities/muni_`y'`k'.dta"
	}
}
use "$path/02_Processed_data/08_Municipalities/Muni_Geo-PLZ_1931.dta", clear
forv v=1932(1)2018  {
	append using "$path/02_Processed_data/08_Municipalities/Muni_Geo-PLZ_`v'.dta"
}
rename E_CNTR GdeNr_E_CNTR  // LV95
rename N_CNTR GdeNr_N_CNTR
rename E PLZ_E_LV03  // LV03
rename N PLZ_N_LV03
drop _merge
save "$path/02_Processed_data/08_Municipalities/GdeNr-Geo_GdeNr-PLZ_1931-2018.dta", replace

* erase
forv y=1931(1)2018 {
	erase "$path/02_Processed_data/08_Municipalities/muni_`y'.dta"
	erase "$path/02_Processed_data/08_Municipalities/Muni_Geo-PLZ_`y'.dta"
}



**** Generate Panel: PLZ - Geo-Coordnates ****

*** Add contemporary PLZ for years prior to 2018: Use PLZ information of Swiss Post and add geo-coordinates accoring to GdeNr-2018PLZ pairings 

** Prepare Swiss Post data: merge it with GdeNr and respective geo-coordinates
import excel "$path/01_Raw_data/08_Municipalities/PLZ-Liste Stand_28.09.2018.xlsx", sheet("Basisdaten") firstrow clear

gen PLZ4 = substr(PLZ6,1,4)
gen PLZ2 = substr(PLZ6,5,2)
destring PLZ4, replace
destring PLZ2, replace
destring PLZ6, replace
order PLZ6 PLZ4 PLZ2

replace PLZABSDatum = date("30-09-2018","DMY") if PLZABSDatum==.  
replace GMDEABSDatumGemeinde = date("30-09-2018","DMY") if GMDEABSDatumGemeinde==.

* Generate year of IBS & ABS (if date in 3 quarter -- > + 1 year)
* PLZ
gen yPLZin = year(PLZIBSDatum) 
gen yPLZout = year(PLZABSDatum) 
sum yPLZin
local yplzinmin = r(min) 
local yplzinmax = r(max)
sum yPLZout
local yplzoutmin = r(min) 
local yplzoutmax = r(max)
gen yPLZin_t =.
gen yPLZout_t =.
foreach z in IBS ABS {  // define year of entry/exit according to rule: if 3rd quarter --> next year
	gen qdate`z'PLZ = qofd(PLZ`z'Datum)  
	format qdate`z'PLZ %tq 
	generate qdate`z'PLZstr = string(qdate`z'PLZ, "%tq")
	replace qdate`z'PLZstr = substr(qdate`z'PLZstr,6,1)
	destring qdate`z'PLZstr, replace
	drop qdate`z'PLZ
	rename qdate`z'PLZstr qdate`z'PLZ
}
forv y = `yplzinmin'(1)`yplzinmax' {
replace yPLZin_t = yPLZin + 1 if qdateIBSPLZ>2 & yPLZin== `y'
}
replace yPLZin = yPLZin_t if yPLZin_t!=.

forv y = `yplzoutmin'(1)`yplzoutmax' {
replace yPLZout_t = yPLZout + 1 if qdateABSPLZ>2 & yPLZout== `y'
}
replace yPLZout = yPLZout_t if yPLZout_t!=.
drop yPLZin_t yPLZout_t qdateIBSPLZ qdateABSPLZ

* GdeNr
gen yGdeNrin = year(GMDEIBSDatumGemeinde) 
gen yGdeNrout = year(GMDEABSDatumGemeinde) 
sum yGdeNrin
local ygdenrinmin = r(min) 
local ygdenrinmax = r(max)
sum yGdeNrout
local ygdenroutmin = r(min) 
local ygdenroutmax = r(max)
gen yGdeNrin_t =.
gen yGdeNrout_t =.
foreach z in IBS ABS {  // define year of entry/exit according to rule: if 3rd quarter --> next year
	gen qdate`z'GdeNr = qofd(GMDE`z'DatumGemeinde)
	format qdate`z'GdeNr %tq
	generate qdate`z'GdeNrstr = string(qdate`z'GdeNr, "%tq")
	replace qdate`z'GdeNrstr = substr(qdate`z'GdeNrstr,6,1)
	destring qdate`z'GdeNrstr, replace
	drop qdate`z'GdeNr
	rename qdate`z'GdeNrstr qdate`z'GdeNr
}
forv y = `ygdenrinmin'(1)`ygdenrinmax' {
replace yGdeNrin_t = yGdeNrin + 1 if qdateIBSGdeNr>2 & yGdeNrin== `y'
}
replace yGdeNrin = yGdeNrin_t if yGdeNrin_t!=.

forv y = `ygdenroutmin'(1)`ygdenroutmax' {
replace yGdeNrout_t = yGdeNrout + 1 if qdateABSGdeNr>2 & yGdeNrout== `y'
}
replace yGdeNrout = yGdeNrout_t if yGdeNrout_t!=.
drop yGdeNrin_t yGdeNrout_t qdateIBSGdeNr qdateABSGdeNr

*br * if yPLZin > yGdeNrout  // no overlap of PLZ and GdeNR -- > no observations
*br * if yPLZout < yGdeNrin  // no overlap of PLZ and GdeNR -- > 572 Obs. -- > inconsistent
*br * if yPLZin < yGdeNrin & yGdeNrout < yPLZout  //  partial overlap: GdeNr changes within PLZ window: 80 Obs.
*br * if yGdeNrin < yPLZin & yPLZout < yGdeNrout  //  partial overlap: PLZ changes within GdeNr window: 8019 Obs.
drop if yPLZout < yGdeNrin

gen gdenr = GMDEBFSNr
sort gdenr
save "$path/02_Processed_data/08_Municipalities/PLZ-Post-Panel_temp.dta", replace

** Use GdeNr panel and merge on Post data
use "$path/02_Processed_data/08_Municipalities/GdeNr-Geo_GdeNr-PLZ_1931-2018.dta", clear 
bysort gdenr: egen ymin = min(year)  // define how many years one specific contemporary GdeNr was in use
bysort gdenr: egen ymax = max(year)
duplicates list gdenr year
duplicates drop gdenr, force
drop year
sort gdenr
save "$path/02_Processed_data/08_Municipalities/GdeNr_Geo-GdeNr_temp.dta", replace  // cross-section: every contemporary GdeNr once, with "in-use" time period

* merge
use "$path/02_Processed_data/08_Municipalities/PLZ-Post-Panel_temp.dta", clear
merge m:1 gdenr using "$path/02_Processed_data/08_Municipalities/GdeNr_Geo-GdeNr_temp.dta"
// _merge == 1 -- > Lichtenstein & some muni in Germany
// _merge == 2 -- > 373 unique observations. many gdenr exited before 1980, some afterwards.
drop if _merge!=3
drop _merge
save "$path/02_Processed_data/08_Municipalities/SwissPost-PLZ_GdeNr2018-Geo.dta", replace  // every PLZ from Post-data receives its GdeNr_2018-Coordinates, some PLZ appear more than once (in different Gde)

* Produce panel structure
sort PLZ6
gen long id = _n
order id
gen yRangePLZ = yPLZout - yPLZin + 1
gen yRangeGdeNr = yGdeNrout - yGdeNrin + 1
expand yRangePLZ
sort id
bysort id: gen idy = _n
bysort id: gen year = yPLZin + idy - 1
order id idy year

duplicates report PLZ6 gdenr year

* which PLZ - GdeNr combinations to keep
/*bysort PLZ6 year: egen sd_gdenr = sd(gdenr)
duplicates report PLZ6 year if sd_gdenr == 0
duplicates drop PLZ6 year if sd_gdenr == 0, force  // drop all observations that are the same within PLZ6 and year and always feature same specific GdeNr
drop sd_gdenr
duplicates report PLZ6 year 
bysort PLZ4 year: egen sd_gdenr = sd(gdenr)
duplicates report PLZ4 year if sd_gdenr == 0
duplicates drop PLZ4 year if sd_gdenr == 0, force  // drop all duplicate observations that are the same within PLZ4 and year and always feature same specific GdeNr
duplicates report PLZ4 year
*/

duplicates report PLZ4 gdenr year
duplicates drop PLZ4 gdenr year, force  // drop all duplicate observations that are the same within PLZ4 and year and always feature same specific GdeNr
duplicates report PLZ4 year
duplicates report PLZ4 gdenr year
duplicates report PLZ4 gdenr_2018 year


* Create contemporary PLZ-year panel with Geo-coordinates (centroids) according to GdeNr2018
// WHAT TO DO ABOUT MULTIPLE GdeNr PER PLZ4?  WE CALCULATE "AVERAGE COORDINATES" per PLZ & year.
preserve 
sort PLZ4 year
gen count = 1
bysort PLZ4 year: egen countPLZ4 = total(count)  // create indicator: 1 if only one GdeNR -- > Coordinates are not averaged
collapse (mean) countPLZ4 GdeNr_E_CNTR GdeNr_N_CNTR (first) gdenr_2018, by (PLZ4 year)
sort PLZ4 year
save "$path/02_Processed_data/08_Municipalities/uniquePLZ-Year-(Avg)GdeNr2018-Geo.dta", replace	// dataset with some averaged coordinates
collapse (mean) countPLZ4 GdeNr_E_CNTR GdeNr_N_CNTR (first) gdenr_2018, by (PLZ4)		// collapse on PLZ without year -- > Bisnode data might require average coordinates per PLZ
sort PLZ4 																				// --> changes of PLZ coordinates != durations of mandates in Bisnode data
save "$path/02_Processed_data/08_Municipalities/uniquePLZ_(Avg)GdeNr2018-Geo.dta", replace 	
restore

* merge geo-coordinates back into the expanded Post dataset
sort PLZ4 year
keep year PLZ6 PLZ4 PLZ2 gdenr gdenr_2018 ctnnr ctn histnr 
merge m:1 PLZ4 year using "$path/02_Processed_data/08_Municipalities/uniquePLZ-Year-(Avg)GdeNr2018-Geo.dta"
drop _merge
save "$path/02_Processed_data/08_Municipalities/PLZ-Year_(Avg)GdeNr2018-Geo.dta", replace

* Create dataset with unique municipality name and its corresponding geo-coordinates
use "$path/02_Processed_data\08_Municipalities/GdeNr2018_Geo_1931-2018.dta", clear
sort gdename
replace gdename = usubinstr(gdename, "'", " ",.)
replace gdename = usubinstr(gdename, "-", " ",.)
replace gdename = usubinstr(gdename, "Ä", "ae",.)
replace gdename = usubinstr(gdename, "Ö", "oe",.)
replace gdename = usubinstr(gdename, "Ü", "ue",.)
replace gdename = usubinstr(gdename, "ä", "ae",.)
replace gdename = usubinstr(gdename, "ö", "oe",.)
replace gdename = usubinstr(gdename, "ü", "ue",.)
replace gdename = usubinstr(gdename, "È", "e",.)
replace gdename = usubinstr(gdename, "É", "e",.)
replace gdename = usubinstr(gdename, "Ê", "e",.)
replace gdename = usubinstr(gdename, "è", "e",.)
replace gdename = usubinstr(gdename, "é", "e",.)
replace gdename = usubinstr(gdename, "ê", "e",.)
replace gdename = usubinstr(gdename, "À", "a",.)
replace gdename = usubinstr(gdename, "Â", "a",.)
replace gdename = usubinstr(gdename, "à", "a",.)
replace gdename = usubinstr(gdename, "â", "a",.)
replace gdename = usubinstr(gdename, "Û", "u",.)
replace gdename = usubinstr(gdename, "Û", "u",.)
replace gdename = usubinstr(gdename, "û", "u",.)
replace gdename = usubinstr(gdename, "ù", "u",.)
replace gdename = usubinstr(gdename, "Ô", "o",.)
replace gdename = usubinstr(gdename, "Ô", "o",.)
replace gdename = usubinstr(gdename, "ô", "o",.)
replace gdename = usubinstr(gdename, "ò", "o",.)
replace gdename = usubinstr(gdename, "Î", "i",.)
replace gdename = usubinstr(gdename, "Î", "i",.)
replace gdename = usubinstr(gdename, "î", "i",.)
replace gdename = usubinstr(gdename, "ì", "i",.)
replace gdename = usubinstr(gdename, "ç", "c",.)
replace gdename = ustrlower(gdename)
replace gdename = ustrtrim(gdename)
duplicates report gdename gdenr_2018
duplicates drop gdename gdenr_2018, force
duplicates list gdename
duplicates drop gdename, force
drop if gdename == ""
keep gdenr_2018 gdename E_CNTR N_CNTR
save "$path/02_Processed_data\08_Municipalities/unique_Gdename_Gdenr2018-Geo.dta", replace


erase "$path/02_Processed_data/08_Municipalities/Gde_1931_1959.dta"
erase "$path/02_Processed_data/08_Municipalities/Gde_1931_2018.dta"
erase "$path/02_Processed_data/08_Municipalities/GdeNr_Geo-GdeNr_temp.dta"
erase "$path/02_Processed_data/08_Municipalities/muni_1931_2018_temp.dta"
erase "$path/02_Processed_data/08_Municipalities/PLZ-Post-Panel_temp.dta"
erase "$path/02_Processed_data/08_Municipalities/PerimeterGDE-GdeNR.dta"
erase "$path/02_Processed_data/08_Municipalities/PerimeterGDE-PLZ.dta"
erase "$path/02_Processed_data/08_Municipalities/PLZ-Year_(Avg)GdeNr2018-Geo.dta"
erase "$path/02_Processed_data/08_Municipalities/GdeNr-Geo_GdeNr-PLZ_1931-2018.dta"

















