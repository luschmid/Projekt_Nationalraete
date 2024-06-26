

forv i=1931(4)1967  {
	foreach j in ZH BE LU UR SZ OW NW GL ZG FR SO BS BL SH AR AI SG GR AG TG TI VD VS NE GE {
	import excel "$path\01_Raw_data\05_Alliances\Bundesblatt_`i'_Listenverbindungen.xlsx", sheet(`j') cellrange(A14:H34) firstrow clear
	drop if Listenname==""
	gen Wahljahr = `i'
	gen Kantonsname = "`j'"
	tostring Listenname Listenverbindung Unterlistenverbindung, replace
	order Kantonsname Wahljahr
	save "$path\02_Processed_data\05_Alliances\Liste_`i'_`j'.dta", replace
	clear
	}
	use "$path\02_Processed_data\05_Alliances\Liste_`i'_ZH.dta", clear
	foreach j in BE LU UR SZ OW NW GL ZG FR SO BS BL SH AR AI SG GR AG TG TI VD VS NE GE {
	append using "$path\02_Processed_data\05_Alliances\Liste_`i'_`j'.dta"
	}
	save "$path\02_Processed_data\05_Alliances\Liste_`i'.dta", replace
	foreach k in ZH BE LU UR SZ OW NW GL ZG FR SO BS BL SH AR AI SG GR AG TG TI VD VS NE GE {
	erase "$path\02_Processed_data\05_Alliances\Liste_`i'_`k'.dta"
	}
}
use "$path\02_Processed_data\05_Alliances\Liste_1931.dta", clear
forv l = 1935(4)1967 {
	append using "$path\02_Processed_data\05_Alliances\Liste_`l'.dta"
}
sort Kantonsname Wahljahr Listennummer_num
tostring Listennummer Parteiname, replace
save "$path\02_Processed_data\05_Alliances\BdBlatt_Listenverbindungen_1931-1967.dta", replace
forv m = 1931(4)1967 {
	erase "$path\02_Processed_data\05_Alliances\Liste_`m'.dta"
}

*** Merge Information from the Bundesblätter (1931-1967) whith BFS data (1971-2015) ***
import excel "$path\01_Raw_data\05_Alliances\BfS_Listenverbindungen_1971-2015.xlsx", sheet(Listenverbindungen_1971-2015) cellrange(A8:K3073) firstrow clear
sort Kantonsname Wahljahr Listennummer_num
save "$path\02_Processed_data\05_Alliances\BfS_Listenverbindungen_1971-2015.dta", replace
append using "$path\02_Processed_data\05_Alliances\BdBlatt_Listenverbindungen_1931-1967.dta"
drop Kantonsnummer
sort Kantonsname Wahljahr Listennummer_num
save "$path\02_Processed_data\05_Alliances\Listenverbindungen_1931-2015_temp1.dta", replace


*** Import additional information from Bundesblätter 1931-1967 ***
forv i=1931(4)1967  {
	foreach j in ZH BE LU UR SZ OW NW GL ZG FR SO BS BL SH AR AI SG GR AG TG TI VD VS NE GE {
	import excel "$path\01_Raw_data\05_Alliances\Bundesblatt_`i'_Listenverbindungen.xlsx", sheet(`j') cellrange(A5:B10) clear
	gen numb = 1
	replace A = "A" if A == "Zahl der Sitze"
	replace A = "B" if A == "Stimmberechtigte"
	replace A = "C" if A == "Stimmende"
	replace A = "D" if A == "Ungültige Wahlzettel"
	replace A = "E" if A == "Leere Wahlzettel"
	replace A = "F" if A == "Gültige Wahlzettel" 
	encode A, gen(AA)
	drop A
	reshape wide B, i(numb) j(AA)
	drop numb
	rename B1 seats_cant
	rename B2 eligible_cant
	rename B3 voters_cant
	rename B4 invalid_cant
	rename B5 empty_cant
	rename B6 valid_cant
	gen Wahljahr = `i'
	gen Kantonsname = "`j'"

	save "$path\02_Processed_data\05_Alliances\BdBlatt_`i'_`j'.dta", replace
	clear
	}
	use "$path\02_Processed_data\05_Alliances\BdBlatt_`i'_ZH.dta", clear
	foreach j in BE LU UR SZ OW NW GL ZG FR SO BS BL SH AR AI SG GR AG TG TI VD VS NE GE {
	append using "$path\02_Processed_data\05_Alliances\BdBlatt_`i'_`j'.dta"
	}
	save "$path\02_Processed_data\05_Alliances\BdBlatt_`i'.dta", replace
	foreach k in ZH BE LU UR SZ OW NW GL ZG FR SO BS BL SH AR AI SG GR AG TG TI VD VS NE GE {
	erase "$path\02_Processed_data\05_Alliances\BdBlatt_`i'_`k'.dta"
	}
}
use "$path\02_Processed_data\05_Alliances\BdBlatt_1931.dta", clear
forv l = 1935(4)1967 {
	append using "$path\02_Processed_data\05_Alliances\BdBlatt_`l'.dta"
}

label var seats_cant "No of seats in National Council"
label var eligible_cant "No eligible voters"
label var voters_cant "No of actual voters in election"
label var invalid_cant "No of invalid votes"
label var empty_cant "No. of empty votes"
label var valid_cant "No. of valid votes"

order Kantonsname Wahljahr
 
*** check if typed information on valid votes is equivalent to calculated number **

*br Kantonsname Wahljahr Sitze Stimmberechtigte Stimmende ungueltige leere gueltige gueltig_calc diffgueltige if diffgueltige!=0 & diffgueltige!=.
/*
Kantonsname	Wahljahr	Sitze	Stimmberechtigte	Stimmende	ungueltige	leere	gueltige	guelt	dguelt	
SO	1935	7	43523	36213	409	480	 35354	35324	30		correct typing of all numbers, some error in official data from Bundesblatt
SG	1939	13	77171	64614	557	2867 61192	61190	2		correct typing of all numbers, some error in official data from Bundesblatt
BS	1943	8	53333	35358	87	94	 34508	35177	-669	correct typing of all numbers, some error in official data from Bundesblatt
BS	1955	8	65204	40286	27	127	 40002	40132	-130	correct typing of all numbers, some error in official data from Bundesblatt
BS	1959	8	67346	40765	37	154	 40442	40574	-132	correct typing of all numbers, some error in official data from Bundesblatt
LU	1959	9	69773	59301	185	766	 58351	58350	1		correct typing of all numbers, some error in official data from Bundesblatt
TG	1959	6	43424	32179	114	675	 31408	31390	18		correct typing of all numbers, some error in official data from Bundesblatt
BS	1967	8	66719	33441	22	96	 33040	33323	-283	correct typing of all numbers, some error in official data from Bundesblatt
*/

order Kantonsname Wahljahr seats_cant
sort Kantonsname Wahljahr
save "$path\02_Processed_data\05_Alliances\NRWahl_ZusatzInfo_BdBlatt_1931-1967.dta", replace

forv m = 1931(4)1967 {
	erase "$path\02_Processed_data\05_Alliances\BdBlatt_`m'.dta"
}


*** Import additional election information 1971-2015 ***
forv v=1971(4)2015 {
import excel "$path\01_Raw_data\05_Alliances\su-d-17.02.02.04.02.xlsx", sheet(`v') cellrange(A10:H36) firstrow clear
	foreach char in A B C D E F G H {
	tostring *, replace
		replace `char' = "" if `char'=="…"
		replace `char' = "" if `char'=="*"
	}
	rename A cantonname
	rename B voters_cant1
	rename C empty_cant1
	rename D invalid_cant1
	rename E valid_cant1
	rename F unchanged_cant1
	rename G changed_cant1
	rename H noparty_cant1
	destring *, replace
	gen canton = ""
	gen year = `v'
	replace canton = "ZH" if cantonname == "Zürich"
	replace canton = "BE" if cantonname == "Bern"
	replace canton = "LU" if cantonname == "Luzern"
	replace canton = "UR" if cantonname == "Uri 1)"
	replace canton = "SZ" if cantonname == "Schwyz"
	replace canton = "OW" if cantonname == "Obwalden 1)"
	replace canton = "NW" if cantonname == "Nidwalden 1)"
	replace canton = "GL" if cantonname == "Glarus 1)"
	replace canton = "ZG" if cantonname == "Zug"
	replace canton = "FR" if cantonname == "Freiburg"
	replace canton = "SO" if cantonname == "Solothurn"
	replace canton = "BS" if cantonname == "Basel-Stadt"
	replace canton = "BL" if cantonname == "Basel-Landschaft"
	replace canton = "SH" if cantonname == "Schaffhausen"
	replace canton = "AR" if cantonname == "Appenzell A. Rh. 1)"
	replace canton = "AI" if cantonname == "Appenzell I. Rh. 1)"
	replace canton = "SG" if cantonname == "St. Gallen"
	replace canton = "GR" if cantonname == "Graubünden"
	replace canton = "AG" if cantonname == "Aargau"
	replace canton = "TG" if cantonname == "Thurgau"
	replace canton = "TI" if cantonname == "Tessin"
	replace canton = "VD" if cantonname == "Waadt"
	replace canton = "VS" if cantonname == "Wallis"
	replace canton = "NE" if cantonname == "Neuenburg"
	replace canton = "GE" if cantonname == "Genf"
	replace canton = "JU" if cantonname == "Jura"
	replace canton = "AR" if cantonname == "Appenzell A. Rh."
	replace canton = "OW" if cantonname == "Obwalden 2)"
	replace canton = "AI" if cantonname == "Appenzell I. Rh. 2)"
	replace canton = "GL" if cantonname == "Glarus 2)"
	replace canton = "UR" if cantonname == "Uri 2)"
	replace canton = "NW" if cantonname == "Nidwalden 2)"
	replace canton = "GL" if cantonname == "Glarus 2)"
	replace canton = "UR" if cantonname == "Uri 2)"
	replace canton = "OW" if cantonname == "Obwalden 1), 2)"
	replace canton = "NW" if cantonname == "Nidwalden 1) 2)"
	replace canton = "ZG" if cantonname == "Zug 2)"
	drop if cantonname == ""
	
save "$path\02_Processed_data\05_Alliances\eletemp1_`v'.dta", replace
}
use "$path\02_Processed_data\05_Alliances\eletemp1_1971.dta"
forv v=1975(4)2015 {
append using "$path\02_Processed_data\05_Alliances\eletemp1_`v'.dta"
}
keep canton year voters_cant1 invalid_cant1 empty_cant1 valid_cant1
sort canton year
order canton year
rename canton Kantonsname
rename year Wahljahr
save "$path\02_Processed_data\05_Alliances\votes_1971-2015.dta", replace
forv x=1971(4)2015 {
	erase "$path\02_Processed_data\05_Alliances\eletemp1_`x'.dta"
}


*** Import voter information for elections 1919-2015 ***
forv v=1971(4)2015 {
import excel "$path\01_Raw_data\05_Alliances\je-d-17.02.02.04.01.xlsx", sheet(`v') cellrange(A8:C34) firstrow clear
	foreach char in A B C {
	tostring *, replace
		replace `char' = "" if `char'=="…"
		replace `char' = "" if `char'=="*"
	}
	rename A cantonname
	rename B eligible_cant2
	rename C voters_cant2
	destring *, replace
	gen canton = ""
	gen year = `v'
	replace canton = "ZH" if cantonname == "Zürich"
	replace canton = "BE" if cantonname == "Bern"
	replace canton = "LU" if cantonname == "Luzern"
	replace canton = "UR" if cantonname == "Uri 1)"
	replace canton = "SZ" if cantonname == "Schwyz"
	replace canton = "OW" if cantonname == "Obwalden 1)"
	replace canton = "NW" if cantonname == "Nidwalden 1)"
	replace canton = "GL" if cantonname == "Glarus 1)"
	replace canton = "ZG" if cantonname == "Zug"
	replace canton = "FR" if cantonname == "Freiburg"
	replace canton = "SO" if cantonname == "Solothurn"
	replace canton = "BS" if cantonname == "Basel-Stadt"
	replace canton = "BL" if cantonname == "Basel-Landschaft"
	replace canton = "SH" if cantonname == "Schaffhausen"
	replace canton = "AR" if cantonname == "Appenzell A. Rh. 1)"
	replace canton = "AI" if cantonname == "Appenzell I. Rh. 1)"
	replace canton = "SG" if cantonname == "St. Gallen"
	replace canton = "GR" if cantonname == "Graubünden"
	replace canton = "AG" if cantonname == "Aargau"
	replace canton = "TG" if cantonname == "Thurgau"
	replace canton = "TI" if cantonname == "Tessin"
	replace canton = "VD" if cantonname == "Waadt"
	replace canton = "VS" if cantonname == "Wallis"
	replace canton = "NE" if cantonname == "Neuenburg"
	replace canton = "GE" if cantonname == "Genf"
	replace canton = "JU" if cantonname == "Jura"
	replace canton = "NW" if cantonname == "Nidwalden 1) 3)"
	replace canton = "OW" if cantonname == "Obwalden 1), 2)"	
	replace canton = "UR" if cantonname == "Uri 2)"	
	replace canton = "OW" if cantonname == "Obwalden 2)"
	replace canton = "NW" if cantonname == "Nidwalden 2)"
	replace canton = "GL" if cantonname == "Glarus 2)"
	replace canton = "AR" if cantonname == "Appenzell A. Rh."
	replace canton = "AR" if cantonname == "Appenzell A. Rh. 2)"
	replace canton = "AI" if cantonname == "Appenzell I. Rh. 2)"
	replace canton = "ZG" if cantonname == "Zug 2)"
	drop if cantonname == ""
	
save "$path\02_Processed_data\05_Alliances\eletemp2_`v'.dta", replace
}
use "$path\02_Processed_data\05_Alliances\eletemp2_1971.dta"
forv v=1975(4)2015 {
append using "$path\02_Processed_data\05_Alliances\eletemp2_`v'.dta"
}
keep canton year eligible_cant2 voters_cant2
sort canton year
order canton year
rename canton Kantonsname
rename year Wahljahr
save "$path\02_Processed_data\05_Alliances\eligible_1971-2015.dta", replace
forv x=1971(4)2015 {
	erase "$path\02_Processed_data\05_Alliances\eletemp2_`x'.dta"
}

forv v=1919(3)1931 {
import excel "$path\01_Raw_data\05_Alliances\je-d-17.02.02.04.01.xlsx", sheet(`v') cellrange(B8:D33) firstrow clear
	foreach char in B C D {
	tostring *, replace
		replace `char' = "" if `char'=="…"
		replace `char' = "" if `char'=="*"
	}
	rename B cantonname
	rename C eligible_cant2
	rename D voters_cant2
	destring *, replace
	gen canton = ""
	gen year = `v'
	replace canton = "ZH" if cantonname == "Zürich"
	replace canton = "BE" if cantonname == "Bern"
	replace canton = "LU" if cantonname == "Luzern"
	replace canton = "UR" if cantonname == "Uri"
	replace canton = "SZ" if cantonname == "Schwyz"
	replace canton = "OW" if cantonname == "Obwalden"
	replace canton = "NW" if cantonname == "Nidwalden"
	replace canton = "GL" if cantonname == "Glarus"
	replace canton = "ZG" if cantonname == "Zug"
	replace canton = "FR" if cantonname == "Freiburg"
	replace canton = "SO" if cantonname == "Solothurn"
	replace canton = "BS" if cantonname == "Basel-Stadt"
	replace canton = "BL" if cantonname == "Basel-Landschaft"
	replace canton = "SH" if cantonname == "Schaffhausen"
	replace canton = "AR" if cantonname == "Appenzell A.Rh."
	replace canton = "AI" if cantonname == "Appenzell I.Rh."
	replace canton = "SG" if cantonname == "St. Gallen"
	replace canton = "GR" if cantonname == "Graubünden"
	replace canton = "AG" if cantonname == "Aargau"
	replace canton = "TG" if cantonname == "Thurgau"
	replace canton = "TI" if cantonname == "Tessin"
	replace canton = "VD" if cantonname == "Waadt"
	replace canton = "VS" if cantonname == "Wallis"
	replace canton = "NE" if cantonname == "Neuenburg"
	replace canton = "GE" if cantonname == "Genf"
	drop if cantonname == ""
	
save "$path\02_Processed_data\05_Alliances\eletemp3_`v'.dta", replace
}
use "$path\02_Processed_data\05_Alliances\eletemp3_1919.dta"
forv v=1922(3)1931 {
append using "$path\02_Processed_data\05_Alliances\eletemp3_`v'.dta"
}
keep canton year eligible_cant2 voters_cant2
sort canton year
order canton year
rename canton Kantonsname
rename year Wahljahr
save "$path\02_Processed_data\05_Alliances\eligible_1919-1931.dta", replace
forv x=1919(3)1931 {
	erase "$path\02_Processed_data\05_Alliances\eletemp3_`x'.dta"
}

forv v=1935(4)1967 {
import excel "$path\01_Raw_data\05_Alliances\je-d-17.02.02.04.01.xlsx", sheet(`v') cellrange(B8:D33) firstrow clear
	foreach char in B C D {
	tostring *, replace
		replace `char' = "" if `char'=="…"
		replace `char' = "" if `char'=="*"
	}
	rename B cantonname
	rename C eligible_cant2
	rename D voters_cant2
	destring *, replace
	gen canton = ""
	gen year = `v'
	replace canton = "ZH" if cantonname == "Zürich"
	replace canton = "BE" if cantonname == "Bern"
	replace canton = "LU" if cantonname == "Luzern"
	replace canton = "UR" if cantonname == "Uri"
	replace canton = "SZ" if cantonname == "Schwyz"
	replace canton = "OW" if cantonname == "Obwalden"
	replace canton = "NW" if cantonname == "Nidwalden"
	replace canton = "GL" if cantonname == "Glarus"
	replace canton = "ZG" if cantonname == "Zug"
	replace canton = "FR" if cantonname == "Freiburg"
	replace canton = "SO" if cantonname == "Solothurn"
	replace canton = "BS" if cantonname == "Basel-Stadt"
	replace canton = "BL" if cantonname == "Basel-Landschaft"
	replace canton = "SH" if cantonname == "Schaffhausen"
	replace canton = "AR" if cantonname == "Appenzell A.Rh."
	replace canton = "AI" if cantonname == "Appenzell I.Rh."
	replace canton = "SG" if cantonname == "St. Gallen"
	replace canton = "GR" if cantonname == "Graubünden"
	replace canton = "AG" if cantonname == "Aargau"
	replace canton = "TG" if cantonname == "Thurgau"
	replace canton = "TI" if cantonname == "Tessin"
	replace canton = "VD" if cantonname == "Waadt"
	replace canton = "VS" if cantonname == "Wallis"
	replace canton = "NE" if cantonname == "Neuenburg"
	replace canton = "GE" if cantonname == "Genf"
	drop if cantonname == ""
	
save "$path\02_Processed_data\05_Alliances\eletemp3_`v'.dta", replace
}
use "$path\02_Processed_data\05_Alliances\eletemp3_1935.dta"
forv v=1939(4)1967 {
append using "$path\02_Processed_data\05_Alliances\eletemp3_`v'.dta"
}
keep canton year eligible_cant2 voters_cant2
sort canton year
order canton year
rename canton Kantonsname
rename year Wahljahr
save "$path\02_Processed_data\05_Alliances\eligible_1935-1967.dta", replace
forv x=1935(4)1967 {
	erase "$path\02_Processed_data\05_Alliances\eletemp3_`x'.dta"
}

use "$path\02_Processed_data\05_Alliances\eligible_1919-1931.dta", clear
append using "$path\02_Processed_data\05_Alliances\eligible_1935-1967.dta"
append using "$path\02_Processed_data\05_Alliances\eligible_1971-2015.dta"
drop if Wahljahr<1931
sort Kantonsname Wahljahr
save "$path\02_Processed_data\05_Alliances\eligible_1931-2015.dta", replace

erase "$path\02_Processed_data\05_Alliances\eligible_1919-1931.dta"
erase "$path\02_Processed_data\05_Alliances\eligible_1935-1967.dta"
erase "$path\02_Processed_data\05_Alliances\eligible_1971-2015.dta"


*** Combine Datasets: "Listenverbindungen_1931-2015_temp" and NRWahl_ZusatzInfo_BdBlatt_1931-1967 & BfS Date with eligible voters, invalid votes, etc. ***
use "$path\02_Processed_data\05_Alliances\Listenverbindungen_1931-2015_temp1.dta", clear
merge m:1 Kantonsname Wahljahr using "$path\02_Processed_data\05_Alliances\NRWahl_ZusatzInfo_BdBlatt_1931-1967.dta"
tab _merge
drop _merge 

merge m:1 Kantonsname Wahljahr using "$path\02_Processed_data\05_Alliances\votes_1971-2015.dta"
replace voters_cant = voters_cant1 if Wahljahr>=1971
replace invalid_cant = invalid_cant1 if Wahljahr>=1971
replace empty_cant = empty_cant1 if Wahljahr>=1971
replace valid_cant = valid_cant1 if Wahljahr>=1971
drop voters_cant1 invalid_cant1 empty_cant1 valid_cant1
tab _merge
drop _merge
*save "$path\02_Processed_data\05_Alliances\Listenverbindungen_1931-2015_temp3.dta", replace

merge m:1 Kantonsname Wahljahr using "$path\02_Processed_data\05_Alliances\eligible_1931-2015.dta"
gen diffelig = eligible_cant - eligible_cant2
bysort Wahljahr Kantonsname: g indi_Wahljahr_Kanton=_n 
*br Wahljahr Kantonsname eligible_cant eligible_cant2 diffelig if abs(diffelig)>0 & indi_Wahljahr_Kanton==1 & diffelig!=. // there are no difference between BfS data and Bundesblatt data -- > take BfS
replace eligible_cant = eligible_cant2 if eligible_cant2!=. 
gen diffvoters = voters_cant - voters_cant2
*br Wahljahr Kantonsname voters_cant voters_cant2 diffvoters if abs(diffvoters)!=0 & indi_Wahljahr_Kanton==1 & diffvoters!=. // there are some differences between BfS data and Bundesblatt data: -- > take BfS
replace voters_cant = voters_cant2 if voters_cant2!=.
drop voters_cant2 eligible_cant2
tab _merge
drop _merge
save "$path\02_Processed_data\05_Alliances\Listenverbindungen_1931-2015_temp2.dta", replace

*** Add missing seat numbers ***
use "$path\01_Raw_data\05_Alliances\NRSeatsPerCanton_1931-2015.dta", clear
rename canton_abr Kantonsname
rename year Wahljahr
sort Kantonsname Wahljahr
merge 1:m Kantonsname Wahljahr using "$path\02_Processed_data\05_Alliances\Listenverbindungen_1931-2015_temp2.dta"
tab _merge

/* BFS does not include observations for several cantons with one single seat (and majority rule):
canton	Wahljahr
Appenzell I.Rh.	1971 - 2007, 2015
Appenzell A.Rh.	1979, 1987, 2003, 2007, 2015
Glarus			1971 - 2007, 2015
Nidwalden 		1971 - 2007, 2015
Obwalden		1971 - 2007, 2015
Uri				1971 - 2007, 2015
Zug				1971 (two NR seats!)
*/

drop _merge
replace seats_cant = seats if seats_cant==.
drop seats
rename canton cantonname // adjust variable names to match main datasets
rename Kantonsname canton
rename cantonname Kantonsname
rename Wahljahr year


rename Listenverbindung alliance
rename Listenstimmen pvotes
rename Unterlistenverbindung suballiance

gen liste_norm=Listenname
replace liste_norm = upper(liste_norm)
replace liste_norm = usubinstr(liste_norm, ".", " ", .)
replace liste_norm = usubinstr(liste_norm, ",", " ", .)
replace liste_norm = usubinstr(liste_norm, ";", " ", .)
replace liste_norm = usubinstr(liste_norm, ":", " ", .)
replace liste_norm = usubinstr(liste_norm, "-", " ", .)
replace liste_norm = usubinstr(liste_norm, "Ä", "A", .)
replace liste_norm = usubinstr(liste_norm, "Ä", "A", .)
replace liste_norm = usubinstr(liste_norm, "Ö", "O", .)
replace liste_norm = usubinstr(liste_norm, "Ö", "O", .)
replace liste_norm = usubinstr(liste_norm, "Ü", "U", .)
replace liste_norm = usubinstr(liste_norm, "Ü", "U", .)
replace liste_norm = usubinstr(liste_norm, "É", "E", .)
replace liste_norm = usubinstr(liste_norm, "É", "E", .)
replace liste_norm = usubinstr(liste_norm, "È", "E", .)
replace liste_norm = usubinstr(liste_norm, "Ê", "E", .)
replace liste_norm = usubinstr(liste_norm, "Â", "A", .)
replace liste_norm = usubinstr(liste_norm, "À", "A", .)
replace liste_norm = usubinstr(liste_norm, "ä", "A", .)
replace liste_norm = usubinstr(liste_norm, "ö", "O", .)
replace liste_norm = usubinstr(liste_norm, "ü", "U", .)
replace liste_norm = usubinstr(liste_norm, "é", "E", .)
replace liste_norm = usubinstr(liste_norm, "è", "E", .)
replace liste_norm = usubinstr(liste_norm, "ê", "E", .)
replace liste_norm = usubinstr(liste_norm, "à", "A", .)
replace liste_norm = usubinstr(liste_norm, "â", "A", .)
replace liste_norm = usubinstr(liste_norm, "  ", " ", .)
replace liste_norm = usubinstr(liste_norm, "«", "", .)
replace liste_norm = usubinstr(liste_norm, "»", "", .)
replace liste_norm = strtrim(liste_norm)
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "AG" & year == 1935 & liste_norm =="FREISINNIG DEMOKRATISCHE VOLKSPARTEI"
replace liste_norm ="SCHWEIZ BAUERNHEIMATBEWEGUNG JUNGBAUERN" if canton == "AG" & year == 1935 & liste_norm =="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG JUNGBAUERN"
replace liste_norm ="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG (JUNGBAUERN)" if canton == "AG" & year == 1943 & liste_norm =="SCHWEIZ BAUERNHEIMATBEWEGUNG JUNGBAUERN"
replace liste_norm ="FREIE STIMMBURGERFUR DIE AUFHEBUNG DES STIMMZWANGS" if canton == "AG" & year == 1959 & liste_norm =="FREIE STIMMBURGER FUR DIE AUFHEBUNG DES STIMMZWANGES"
replace liste_norm ="FREISINNIG DEMOKRATISCHE VOLKSPARTEI UND JUNGLIBERALE BEWEGUNG" if canton == "AG" & year == 1963 & liste_norm =="FREISINNIG DEMOKRATISCHE VOLKSPARTEI UND JUNGLIHERALE BEWEGUNG"
replace liste_norm ="BAUERN  GEWERBE UND BURGERPARTEI (BGB MITTELSTANDSIISTE)" if canton == "AG" & year == 1967 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI (BGB MITTELSTANDSLISTE)"
replace liste_norm ="TEAM 67 LIBERALE AARGAUER FUR EINE MODERNE SCHWEIZ" if canton == "AG" & year == 1967 & liste_norm =="TEAM 1967 LIBERALE AARGAUER FUR EINE MODERNE SCHWEIZ"
replace liste_norm ="FEDERATION LIBERALE POPULAIRE JURASSIENNE" if canton == "BE" & year == 1931 & liste_norm =="JURASSISCHE LIBERALE PARTEI"
replace liste_norm ="PARTI DEMOCRATIQUE CATHOLIQUE DU CANTON DE BERNE" if canton == "BE" & year == 1931 & liste_norm =="KATHOLISCHE VOLKSPARTEI"
replace liste_norm ="BAUERN  GEWERBE UND BURGERPARTEI" if canton == "BE" & year == 1935 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI DES KANTONS BERN"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "BE" & year == 1935 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN UND JUNG LIBERALE BEWEGUNG"
replace liste_norm ="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG JUNGBAUERN" if canton == "BE" & year == 1935 & liste_norm =="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "BE" & year == 1935 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS BERN"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "BE" & year == 1939 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN"
replace liste_norm ="KOMMUNISTISCHE PARTEI" if canton == "BE" & year == 1939 & liste_norm =="KOMMUNISTISCHE PARTEI DES KANTONS BERN"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "BE" & year == 1939 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS BERN"
replace liste_norm ="BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI EMMENTAL MITTELLAND OBERAARGAU SEELAND" if canton == "BE" & year == 1943 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI EMMENTAL MITTELLAND OBERAARGAU  SEELAND"
replace liste_norm ="FREISINNIG DEMOKRATISCHE VOLKSPARTEI OBERAARGAU/EMMENTAL" if canton == "BE" & year == 1943 & liste_norm =="FREISINNIG DEMOKRATISCHE VOLKSLISTE OBERAARGAU EMMENTAL"
replace liste_norm ="LANDESTEILVERBAND OBERLAND DER BERN BAUERN  GEWERBE UND BURGERPARTEI" if canton == "BE" & year == 1943 & liste_norm =="LANDESTEILVERBAND OBERLAND DER BAUERN  GEWERBE UND BURGERPARTEI"
replace liste_norm ="PARTI DEMOCRATIQUE CATHOLIQUE DU CANTON DE BERNE" if canton == "BE" & year == 1943 & liste_norm =="PARTI DEMOCRATIQUE CATHOLIQUE"
replace liste_norm ="SCHWEIZERISCHE BAUERN HEIMATBEWEGUNG" if canton == "BE" & year == 1943 & liste_norm =="SCHWEIZ BAUERNHEIMATBEWEGUNG"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI DES KANTONS BERN" if canton == "BE" & year == 1943 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm ="BERN BAUERN  GEWERBE UND BURGERPARTEI EMMENTAL MITTELLAND OBERAARGAU SEELAND" if canton == "BE" & year == 1947 & liste_norm =="LANDESTEILVERBAND EMMENTAL MITTELLAND OBERAARGAU SEELAND DER BERN BAUERN  GEWERBE UND BURGERPARTEI"
replace liste_norm ="BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI EMMENTAL JURA MITTELLAND OBERAARGAU SEELAND" if canton == "BE" & year == 1951 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI EMMENTAL JURA MITTELLAND OBERAARGAU SEELAND"
replace liste_norm ="BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI LANDESTEIL OBERLAND" if canton == "BE" & year == 1951 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI LANDESTEIL OBERLAND"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEILVERBAND OBERLAND" if canton == "BE" & year == 1951 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI LANDESTEIL VERBAND OBERLAND"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEILE EMMENTAL MITTELLAND OBERAARGAN SEELAND" if canton == "BE" & year == 1951 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI LANDESTEILE EMMENTAL MITTELLAND OBERAARGAU SEELAND"
replace liste_norm ="LIBERAL SOZIALISTISCHE PARTEI DES KANTONS BERN" if canton == "BE" & year == 1951 & liste_norm =="LIBERAL SOZIALISTISCHE PARTEI"
replace liste_norm ="PARTI SOCIALISTE JURASSIE" if canton == "BE" & year == 1951 & liste_norm =="PARTI SOCIALISTE JURASSIEN"
replace liste_norm ="BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI LANDESTEILVERBAND OBERLAND FREIE DEMOKRATISCHE MITTELSTANDSPARTEI" if canton == "BE" & year == 1955 & liste_norm =="BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI LANDESTEIL OBERLAND FREIE DEMOKRATISCHE MITTELSTANDSPARTEI"
replace liste_norm ="PARTI LIBERAL RADICAL JURASSIEN PARTI NATIONAL ROMAND DE BIENNE GROUPE RADICAL ROMAND DE BERNE" if canton == "BE" & year == 1959 & liste_norm =="PARTI LIBERAL RADICAL JURASSIEN"
replace liste_norm ="LANDESRING DER UNABHANGIGEN ALLIANCE DES INDEPENDANTS" if canton == "BE" & year == 1963 & liste_norm =="LANDESRING DER UNABHANGIGEN"
replace liste_norm ="PARTI LIBERAL RADICAL JURASSIEN FREISINNIGE PARTEI DES JURA" if canton == "BE" & year == 1963 & liste_norm =="PARTI LIBERAL RADICAL JURASSIEN (FREISINNIGE PARTEI DES JURA)"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI DES KANTONS BERN PARTI SOCIALISTE DU CANTON DE BERNE" if canton == "BE" & year == 1963 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS BERN"
replace liste_norm ="BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI FREIE DEMOKRATISCHE MITTELSTANDSPARTEI LANDESTEIL OBERLAND" if canton == "BE" & year == 1967 & liste_norm =="BERNISCHE BAUERN  GEWERBE UND BURGERPARTEI FREIE DEMOKRATISCHE MITTELSTANDSPARTEI LANDESTEILVERBAND OBERLAND"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEIL BERN MITTELLAND" if canton == "BE" & year == 1967 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEIL MITTELLAND"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEILE SEELAND  LAUFENTAL" if canton == "BE" & year == 1967 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS BERN LANDESTEILE SEELAND LAUFENTAL"
replace liste_norm ="KONSERVATIV CHRISTLICHSOZIALE VOLKSPARTEI DES KANTONS BERN" if canton == "BE" & year == 1967 & liste_norm =="KONSERVATIV CHRISTLICHSOZIALE VOLKSPARTEI"
replace liste_norm ="PARTI LIBERAL RADICAL JURASSIEN (FREISINNIGE PARTEI DES JURA) PARTI NATIONAL ROMAND DE BIENNE" if canton == "BE" & year == 1967 & liste_norm =="PARTI LIBERAL RADICAL JURASSIEN (FREISINNIGE PARTEI DES JURA)"
replace liste_norm ="BERNISCHE BGB OBERAARGAU (BAUERN  GEWERBE UND BARGERPARTEI)" if canton == "BE" & year == 1971 & liste_norm =="BERNISCHE BGB OBERAARGAU (BAUERN  GEWERBE UND BURGERPARTEI)"
replace liste_norm ="PARTI JURASSIEN DES PAYSANS ARTISANS ET BOURGEOIS(PAB)" if canton == "BE" & year == 1971 & liste_norm =="PAB JURA PARTI JURASSIEN DES PAYSANS ARTISANS ET BOURGEOIS"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "BL" & year == 1931 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI BASELLAND"
replace liste_norm ="KATHOLISCHE VOLKSPARTEI UND CHRISTLICHSOZIALE VEREINIGUNG BASELLAND" if canton == "BL" & year == 1931 & liste_norm =="KATHOLISCHE VOLKSPARTEI UND CHRISTLICH SOZIALE VEREINIGUNG"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "BL" & year == 1931 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI BASELLAND"
replace liste_norm ="KATHOLISCHE VOLKSPARTEI UND CHRISTLICH SOZIALE VEREINIGUNG BASELLAND" if canton == "BL" & year == 1935 & liste_norm =="KATHOLISCHE VOLKSPARTEI UND CHRISTLICHSOZIALE VEREINIGUNG BASELLAND"
replace liste_norm ="SCHWEIZ BAUERNHEIMATBEWEGUNG JUNGBAUERN" if canton == "BL" & year == 1935 & liste_norm =="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG JUNGBAUERN"
replace liste_norm ="KATHOLISCHE VOLKSPARTEI UND CHRISTLICHSOZIALE VEREINIGUNG BASELLAND" if canton == "BL" & year == 1939 & liste_norm =="KATHOLISCHE VOLKSPARTEI UND CHRISTLICH SOZIALE VEREINIGUNG BASELLAND"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1939 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI BASEL LANDSCHAFT"
replace liste_norm ="DEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1951 & liste_norm =="DEMOKRATISCHE PARTEI"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1951 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1951 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm ="DEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1955 & liste_norm =="DEMOKRATISCHE PARTEI"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1955 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1955 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI BASEL LAND" if canton == "BL" & year == 1959 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI BASEL LAND" if canton == "BL" & year == 1959 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm ="WIEDERVEREINIGUNGSFREUNDLICHELISTE AKTION KANTON BASEL" if canton == "BL" & year == 1959 & liste_norm =="WIEDERVEREINIGUNGSFREUNDLICHE LISTE AKTION KANTON BASEL"
replace liste_norm ="BAUERN  GEWERBE UND BURGERPARTEI BASELLAND" if canton == "BL" & year == 1967 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI"
replace liste_norm ="CHRISTLICHSOZIALE VOLKSPARTEI BASELLAND" if canton == "BL" & year == 1967 & liste_norm =="CHRISTLICHSOZIALE VOLKSPARTEI"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1967 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI BASELLAND" if canton == "BL" & year == 1967 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm ="SOZIALDEMOKRATCN UND GEWERKSCHAFTER" if canton == "BL" & year == 1971 & liste_norm =="SOZIALDEMOKRATEN UND GEWERKSCHAFTER"
replace liste_norm ="RADIKALDEMOKRATISCHE PARTEI BASEL" if canton == "BS" & year == 1935 & liste_norm =="RADIKAL DEMOKRATISCHE PARTEI BASEL"
replace liste_norm ="KATHOLISCHE VOLKSPARTEI" if canton == "BS" & year == 1939 & liste_norm =="KATHOLISCHE VOLKSPARTEI BASEL"
replace liste_norm ="KOMMUNISTISCHE PARTEI" if canton == "BS" & year == 1939 & liste_norm =="KOMMUNISTISCHE PARTEI BASEL"
replace liste_norm ="LIBERALE PARTEI" if canton == "BS" & year == 1939 & liste_norm =="LIBERALE PARTEI BASEL"
replace liste_norm ="NATIONALE VOLKSPARTEI" if canton == "BS" & year == 1939 & liste_norm =="NATIONALE VOLKSPARTEI BASEL"
replace liste_norm ="RADIKAL DEMOKRATISCHE PARTEI" if canton == "BS" & year == 1939 & liste_norm =="RADIKAL DEMOKRATISCHE PARTEI BASEL"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "BS" & year == 1939 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI BASEL"
replace liste_norm ="RADIKAL DEMOKRATISCHE PARTEI" if canton == "BS" & year == 1943 & liste_norm =="RADIKAL DEMOKRATISCHE PARTEI BASEL"
replace liste_norm ="BURGER UND GEWERBE PARTEI" if canton == "BS" & year == 1947 & liste_norm =="BURGER UND GEWERBEPARTEI"
replace liste_norm ="RADIKAL DEMOKRATISCHE PARTEI" if canton == "BS" & year == 1951 & liste_norm =="RADIKALDEMOKRATISCHE PARTEI"
replace liste_norm ="RADIKAL DEMOKRATISCHE PARTEI" if canton == "BS" & year == 1955 & liste_norm =="RADIKALDEMOKRATISCHE PARTEI"
replace liste_norm ="LISTE7 KATHOLISCHE VOLKSPARTEI" if canton == "BS" & year == 1959 & liste_norm =="KATHOLISCHE VOLKSPARTEI"
replace liste_norm ="KATHOLISCHE UND CHRISTLICHSOZIALE VOLKSPARTEI" if canton == "BS" & year == 1963 & liste_norm =="KATHOLISCHE UND CHRISTLICHSOZIALE VOLKSPARTEI KANTONALE ORGANISATION DER KONSERVATIV CHRISTLICHSOZIALEN VOLKSPARTEI DER SCHWEIZ"
replace liste_norm ="LIBERAL DEMOKRATISCHE BURGERPARTEI" if canton == "BS" & year == 1963 & liste_norm =="LIBERAL DEMOKRATISCHE BURGERPARTEI BASEL STADT"
replace liste_norm ="RADIKAL DEMOKRATISCHE PARTEI" if canton == "BS" & year == 1963 & liste_norm =="RADIKAL DEMOKRATISCHE PARTEI BASEL STADT KANTONALE ORGANISATION DER FREISINNIG DEMOKRATISCHEN PARTEI DER SCHWEIZ"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "BS" & year == 1963 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI BASEL STADT"
replace liste_norm ="CHRISTLICHDEMOKRATISCHE VOLKSPARTEI BASEL STADT (CVF)" if canton == "BS" & year == 1971 & liste_norm =="CHRISTLICHDEMOKRATISCHE VOLKSPARTEI BASEL STADT (CVP)"
replace liste_norm ="NATIONALE AKTION GEGEN DIE UBERFREMDUNG VONVOLK UND HEIMAT" if canton == "BS" & year == 1971 & liste_norm =="NATIONALE AKTION GEGEN DIE UBERFREMDUNG VON VOLK UND HEIMAT"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI BASEL STADT" if canton == "BS" & year == 1971 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI BS"
replace liste_norm ="LISTE CONSERVATRICE CHRETIENNE SOCIALE" if canton == "FR" & year == 1963 & liste_norm =="LISTE CONSERVATRICE CHRETIENNE SOCIALE  KONSERVATIV CHRISTLICHSOZIALE LISTE"
replace liste_norm ="LISTE DES PAYSANS ARTISANS ET INDEPENDANTS" if canton == "FR" & year == 1963 & liste_norm =="LISTE DES PAYSANS ARTISANS ET INDEPENDANTS  BAUERN GEWERBE UND BURGERPARTEI"
replace liste_norm ="LISTE RADICALE DEMOCRATIQUE" if canton == "FR" & year == 1963 & liste_norm =="LISTE RADICALE DEMOCRATIQUE  FREISINNIG DEMOKRATISCHE LISTE"
replace liste_norm ="LISTE SOCIALISTE" if canton == "FR" & year == 1963 & liste_norm =="LISTE SOCIALISTE  SOZIALISTISCHE LISTE"
replace liste_norm ="RADICALE DEMOCRATIQUE" if canton == "FR" & year == 1967 & liste_norm =="LISTE RADICALE DEMOCRATIQUE  RADIKAL DEMOKRATISCHE LISTE"
replace liste_norm ="SOCIALISTE" if canton == "FR" & year == 1967 & liste_norm =="LISTE SOCIALISTE  SOZIALISTISCHE LISTE"
replace liste_norm ="PARTI CONSERVATEUR CHRETIEN SOCIAL ET PARTI INDEPENDANT CHRETIEN SOCIAL" if canton == "FR" & year == 1967 & liste_norm =="PARTI CONSERVATEUR CHRETIEN SOCIAL ET DU PARTI INDEPENDANT CHRETIEN SOCIAL  KONSERVATIVCHRISTLICHSOZIALE VOLKSPARTEI UND DER UNABHANGIG CHRISTLICHSOZIALEN PARTEI"
replace liste_norm ="PARTI FRIBOURGEOIS DES PAYSANS ARTISANS ET DES INDEPENDANTS" if canton == "FR" & year == 1967 & liste_norm =="PARTI FRIBOURGEOIS DES PAYSANS ARTISANS ET DES INDEPENDANTS  BAUERN  GEWERBE UND BURGERPARTEI DES KANTONS FREIBURG"
replace liste_norm ="PARTI SOCIALISTE" if canton == "GE" & year == 1931 & liste_norm =="PARTI SOCIALISTE GENEVOIS"
replace liste_norm ="PARTI INDEPENDANT ET CHRETIEN SOCIAL GENEVOIS" if canton == "GE" & year == 1939 & liste_norm =="PARTI INDEPENDANT CHRETIEN SOCIAL"
replace liste_norm ="PART RADICAL" if canton == "GE" & year == 1939 & liste_norm =="PARTI RADICAL"
replace liste_norm ="PARTI SOCIALISTE DE GENEVE" if canton == "GE" & year == 1939 & liste_norm =="PARTI SOCIALISTE DE GENÈVE"
replace liste_norm ="PARTI SOCIALISTE GENEVOIS" if canton == "GE" & year == 1955 & liste_norm =="PARTI SOCIALISTE"
replace liste_norm ="PARTI SOCIALISTE GENEVOIS" if canton == "GE" & year == 1959 & liste_norm =="PARTI SOCIALISTE"
replace liste_norm ="KONSERVATIV DEMOKRATISCHEPARTEI" if canton == "GR" & year == 1935 & liste_norm =="KONSERVATIV DEMOKRATISCHE PARTEI"
replace liste_norm ="FREISINNIG DEMOKRATISCH PARTEI" if canton == "GR" & year == 1951 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm ="KONSERVATIVE UND CHRISTLICHSOZIALE VOLKSPARTEI" if canton == "GR" & year == 1963 & liste_norm =="KONSERVATIV CHRISTLICHSOZIALE VOLKSPARTEI"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI (ARBEITERUNION)" if canton == "LU" & year == 1931 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm ="KONSERVATIVE UND CHRISTLICHSOZIALE VOLKSPARTEI DES KANTONS LUZERN" if canton == "LU" & year == 1943 & liste_norm =="KONSERVATIVE UND CHRISTLICHSOZIALE VOLKSPARTEI"
replace liste_norm ="LIBERALE PARTEI DES KANTONS LUZERN" if canton == "LU" & year == 1943 & liste_norm =="LIBERALE PARTEI"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI DES KANTONS LUZERN" if canton == "LU" & year == 1943 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI UND GEWERKSCHAFTSKARTELL" if canton == "LU" & year == 1951 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI UND GEWERKSCHAFTSKARTEL"
replace liste_norm ="LIBERATE PARTEI" if canton == "LU" & year == 1963 & liste_norm =="LIBERALE PARTEI"
replace liste_norm ="LIBERALE PARTEI" if canton == "LU" & year == 1967 & liste_norm =="LIBERATE PARTEI"
replace liste_norm ="PARTI DEMOCRATE POPULAIRE NEUCHATELOIS" if canton == "NE" & year == 1931 & liste_norm =="PARTI DEMOCRATE POPULAIRE NEUCHÂTELOIS"
replace liste_norm ="LIBERALE" if canton == "NE" & year == 1935 & liste_norm =="LISTE LIBERALE"
replace liste_norm ="PROGRESSISTE NATIONALE" if canton == "NE" & year == 1935 & liste_norm =="LISTE PROGRESSISTE NATIONALE"
replace liste_norm ="RADICALE" if canton == "NE" & year == 1935 & liste_norm =="LISTE RADICALE"
replace liste_norm ="SOCIALISTE" if canton == "NE" & year == 1935 & liste_norm =="LISTE SOCIALISTE"
replace liste_norm ="PARTI SOCIALISTE NEUCHATELOIS" if canton == "NE" & year == 1943 & liste_norm =="PARTI SOCIALISTE NEUCHÂTELOIS"
replace liste_norm ="LISTE DU PARTI LIBERAL" if canton == "NE" & year == 1947 & liste_norm =="PARTI LIBERAL"
replace liste_norm ="LISTE DU PARTI OUVRIER ET POPULAIRE" if canton == "NE" & year == 1947 & liste_norm =="PARTI OUVRIER ET POPULAIRE"
replace liste_norm ="LISTE DU PARTI RADICAL" if canton == "NE" & year == 1947 & liste_norm =="PARTI RADICAL"
replace liste_norm ="LISTE DU PARTI SOCIALISTE" if canton == "NE" & year == 1947 & liste_norm =="PARTI SOCIALISTE"
replace liste_norm ="LISTE LIBERALE" if canton == "NE" & year == 1951 & liste_norm =="PARTI LIBERAL"
replace liste_norm ="LISTE DU PARTI OUVRIER ET POPULAIRE" if canton == "NE" & year == 1951 & liste_norm =="PARTI OUVRIER ET POPULAIRE NEUCHATELOIS"
replace liste_norm ="LISTE RADICALE" if canton == "NE" & year == 1951 & liste_norm =="PARTI RADICAL"
replace liste_norm ="LISTE SOCIALISTE" if canton == "NE" & year == 1951 & liste_norm =="PARTI SOCIALISTE"
replace liste_norm ="PARTI OUVRIER ET POPULAIRE" if canton == "NE" & year == 1963 & liste_norm =="PARTI OUVRIER ET POPULAIRE NEUCHATELOIS"
replace liste_norm ="SCHWEIZ BAUERNHEIMATBEWEGUNG JUNGBAUERN" if canton == "SG" & year == 1935 & liste_norm =="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG JUNGBAUERN"
replace liste_norm ="SCHWEIZERISCHE BAUERNHEIMATBEWEGUNG (JUNGBAUERN)" if canton == "SG" & year == 1939 & liste_norm =="SCHWEIZ BAUERNHEIMATBEWEGUNG (JUNGBAUERN)"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG O" if canton == "SG" & year == 1959 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG" if canton == "SG" & year == 1963 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI UND DER JUNGLIBERALEN BEWEGUNG"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI UND GEWERKSCHAFTEN" if canton == "SG" & year == 1963 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI UND DER GEWERKSCHAFTEN"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI UND GEWERKSCHAFTEN" if canton == "SG" & year == 1967 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI UNDJUNGLIBERALE BEWEGUNG LISTE NORD" if canton == "SG" & year == 1971 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG LISTE NORD"
replace liste_norm ="PROGRESSIVE ORGANISATIONEN ST GALLEN (POSG)" if canton == "SG" & year == 1971 & liste_norm =="PROGRESSIVE ORGANISATION ST GALLEN (POSG)"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "SH" & year == 1931 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS SCHAFFHAUSEN"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "SH" & year == 1931 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS SCHAFFHAUSEN"
replace liste_norm ="SCHAFFHAUSER BAUERNPARTEI" if canton == "SH" & year == 1943 & liste_norm =="BAUERNPARTEI"
replace liste_norm ="LISTE DER BAUERNPARTEI" if canton == "SH" & year == 1947 & liste_norm =="BAUERNPARTEI"
replace liste_norm ="LISTE DER FREISINNIG DEMOKRATISCHEN PARTEI" if canton == "SH" & year == 1947 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm ="LISTE DER KATHOLISCHEN VOLKSPARTEI" if canton == "SH" & year == 1947 & liste_norm =="KATHOLISCHE VOLKSPARTEI"
replace liste_norm ="LISTE DER SOZIALISTISCHEN ARBEITERPARTEI" if canton == "SH" & year == 1947 & liste_norm =="SOZIALISTISCHE ARBEITERPARTEI"
replace liste_norm ="SOZIALISTISCHEARBEITERPARTEI" if canton == "SH" & year == 1959 & liste_norm =="SOZIALISTISCHE PARTEI"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "SO" & year == 1931 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS SOLOTHURN"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "SO" & year == 1931 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS SOLOTHURN"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "SO" & year == 1935 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI DES KANTONS SOLOTHURN"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "SO" & year == 1935 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI DES KANTONS SOLOTHURN"
replace liste_norm ="SOLOTHURNISCHE BAUERN  GEWERBE UND BURGERPARTEI" if canton == "SO" & year == 1947 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI"
replace liste_norm ="SOLOTHURNISCHE VOLKSPARTEI" if canton == "SO" & year == 1947 & liste_norm =="VOLKSPARTEI"
replace liste_norm ="SOLOTHURNISCNE VOLKSPARTEI UND CHRISTLICHSOZIALE" if canton == "SO" & year == 1951 & liste_norm =="SOLOTHURNISCHE VOLKSPARTEI UND CHRISTLICHSOZIALE"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG" if canton == "SO" & year == 1955 & liste_norm =="PREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PAXTEI UND JUNGLIBERALE BEWEGUNG" if canton == "SO" & year == 1959 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI UND JUNGLIBERALE BEWEGUNG"
replace liste_norm ="ARBEITERPARTEI" if canton == "SZ" & year == 1931 & liste_norm =="ARBEITERPARTEI DES KANTONS SCHWYZ"
replace liste_norm ="KONSERVATIVE VOLKSPARTEI" if canton == "SZ" & year == 1931 & liste_norm =="KONSERVATIVE VOLKSPARTEI DES KANTONS SCHWYZ"
replace liste_norm ="LIBERALE VOLKSPARTEI" if canton == "SZ" & year == 1931 & liste_norm =="LIBERALE VOLKSPARTEI DES KANTONS SCHWYZ"
replace liste_norm ="ARBEITERUNION" if canton == "SZ" & year == 1935 & liste_norm =="ARBEITERUNION DES KANTONS SCHWYZ"
replace liste_norm ="KONSERVATIVE VOLKSPARTEI" if canton == "SZ" & year == 1935 & liste_norm =="KONSERVATIVE VOLKSPARTEI DES KANTONS SCHWYZ"
replace liste_norm ="LIBERALE VOLKSPARTEI" if canton == "SZ" & year == 1935 & liste_norm =="LIBERALE VOLKSPARTEI DES KANTONS SCHWYZ"
replace liste_norm ="BAUERNVEREINIGUNG" if canton == "SZ" & year == 1947 & liste_norm =="BAUERN VEREINIGUNG"
replace liste_norm ="ARBEITERUNION" if canton == "SZ" & year == 1959 & liste_norm =="LISTE DER ARBEITERUNION"
replace liste_norm ="CHRISTLICHSOZIALE PARTEI" if canton == "SZ" & year == 1959 & liste_norm =="LISTE DER CHRISTLICHSOZIALEN PARTEI"
replace liste_norm ="KONSERVATIVE VOLKSPARTEI" if canton == "SZ" & year == 1959 & liste_norm =="LISTE DER KONSERVATIVEN VOLKSPARTEI"
replace liste_norm ="LIBERALE VOLKSPARTEI UND JUNGLIBERALE BEWEGUNG" if canton == "SZ" & year == 1959 & liste_norm =="LISTE DER LIBERALEN VOLKSPARTEI UND JUNGLIBERALEN BEWEGUNG"
replace liste_norm ="ARBEITER AND ANGESTELLTENUNION" if canton == "SZ" & year == 1963 & liste_norm =="ARBEITER UND ANGESTELLTENUNION"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "TG" & year == 1935 & liste_norm =="LISTE DER FREISINNIG DEMOKRATISCHEN PARTEI"
replace liste_norm ="KATHOLISCHE VOLKSPARTEI" if canton == "TG" & year == 1935 & liste_norm =="LISTE DER KATHOLISCHEN VOLKSPARTEI"
replace liste_norm ="BAUERNHEIMATBEWEGUNG JUNGBAUERN" if canton == "TG" & year == 1935 & liste_norm =="LISTE DER SCHWEIZERISCHEN BAUERNHEIMATBEWEGUNG (JUNGBAUERN LISTE)"
replace liste_norm ="JUNG THURGAU" if canton == "TG" & year == 1935 & liste_norm =="LISTE JUNG THURGAU"
replace liste_norm ="SOZIALDEMOKRATISCHE LISTE" if canton == "TG" & year == 1939 & liste_norm =="SOZIALDEMOKRATISCHE PARTEI"
replace liste_norm ="BAUERNPARTEI" if canton == "TG" & year == 1947 & liste_norm =="BAUERN PARTEI"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "TG" & year == 1947 & liste_norm =="SOZIALDEMOKRATISCH GEWERKSCHAFTLICHE PARTEI"
replace liste_norm ="FREISINNIG DEMOKRATISCHE LISTE" if canton == "TG" & year == 1951 & liste_norm =="FREISINNIG DEMOKRATISCHE PARTEI"
replace liste_norm ="SOZIALDEMOKRATISCHE GEWERKSCHAFTLICHELISTE" if canton == "TG" & year == 1951 & liste_norm =="SOZIALDEMOKRATISCHE UND GEWERKSCHAFTLICHE LISTE"
replace liste_norm ="BAUERNPARTEI" if canton == "TG" & year == 1955 & liste_norm =="BAUERNLISTE"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "TG" & year == 1955 & liste_norm =="FREISINNIG DEMOKRATISCHE LISTE"
replace liste_norm ="CHRISTLICHSOZIALE PARTEI" if canton == "TG" & year == 1955 & liste_norm =="LISTE DER CHRISTLICHSOZIALEN"
replace liste_norm ="SOZIALDEMOKRATISCHE UND GEWERKSCHAFTLICHE PARTEI" if canton == "TG" & year == 1955 & liste_norm =="LISTE DER SOZIALDEMOKRATISCHEN PARTEI UND GEWERKSCHAFTEN"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI UND GEWERKSCHAFTEN" if canton == "TG" & year == 1959 & liste_norm =="LISTE DER SOZIALDEMOKRATISCHEN PARTEI UND GEWERKSCHAFTEN"
replace liste_norm ="CHRISTLICH DEMOKRATISCHE VOLKSPARTEI" if canton == "TG" & year == 1971 & liste_norm =="CHRISTLICH DEMOKRATISCHE VOKSPARTEI"
replace liste_norm ="PARTITO COMUNISTA" if canton == "TI" & year == 1931 & liste_norm =="PARTITO COMMUNISTA"
replace liste_norm ="LISTA DEL GRUPPO AGRARIO POPOLARE TICINESE" if canton == "TI" & year == 1947 & liste_norm =="GRUPPO AGRARIO POPOLARE TICINESE"
replace liste_norm ="LISTA DEL GRUPPO CONSERVATORE DEMOCRATICO" if canton == "TI" & year == 1947 & liste_norm =="GRUPPO CONSERVATORE DEMOCRATICO"
replace liste_norm ="LISTA DEL PARTITO OPERAIO E CONTADINO TICINESE" if canton == "TI" & year == 1947 & liste_norm =="GRUPPO OPERAIO E CONTADINO TICINESE"
replace liste_norm ="LISTA DEL GRUPPO LIBERALE RADICALE" if canton == "TI" & year == 1947 & liste_norm =="PARTITO LIBERALE RADICALE"
replace liste_norm ="LISTA DEL PARTITO SOCIALISTE" if canton == "TI" & year == 1947 & liste_norm =="PARTITO SOCIALISTA"
replace liste_norm ="LISTA DEL GRUPPO CONSERVATORE DEMOCRATICO" if canton == "TI" & year == 1951 & liste_norm =="GRUPPO CONSERVATORE DEMOCRATICO"
replace liste_norm ="LISTA LIBERALE RADICALE" if canton == "TI" & year == 1951 & liste_norm =="PARTITO LIBERALE RADICALE"
replace liste_norm ="LISTA DEL PARTITO SOCIALISTA TICINESE" if canton == "TI" & year == 1951 & liste_norm =="PARTITO SOCIALISTA TICINESE"
replace liste_norm ="CONSERVATORE DEMOCRATICO" if canton == "TI" & year == 1955 & liste_norm =="GRUPPO CONSERVATORE DEMOCRATICO"
replace liste_norm ="LIBERALE RADICALE" if canton == "TI" & year == 1955 & liste_norm =="PARTITO LIBERALE RADICALE"
replace liste_norm ="CONSERVATORE DEMOCRATICO" if canton == "TI" & year == 1959 & liste_norm =="GRUPPO CONSERVATORE DEMOCRATICO"
replace liste_norm ="LIBERALE RADICALE" if canton == "TI" & year == 1959 & liste_norm =="PARTITO LIBERALE RADICALE"
replace liste_norm ="PARTITO AGRARI ARTIGIANI ET PATRIZI" if canton == "TI" & year == 1963 & liste_norm =="PARTITO AGRARI ARTIGIANI E PATRIZI"
replace liste_norm ="LIBERALE RADICALE" if canton == "TI" & year == 1967 & liste_norm =="PARTITO LIBERALE RADICALE"
replace liste_norm ="PARTITE DEL LAVORO" if canton == "TI" & year == 1971 & liste_norm =="PARTITO DEL LAVORO"
replace liste_norm ="PARTITO LIBERATE RADICALE TICINESE" if canton == "TI" & year == 1971 & liste_norm =="PARTITO LIBERALE RADICALE TICINESE"
replace liste_norm ="PARTITE POPOLARE DEMOCRATICO PPD" if canton == "TI" & year == 1971 & liste_norm =="PARTITO POPOLARE DEMOCRATICO (PPD)"
replace liste_norm ="LIBERALE DEMOCRATIQUE" if canton == "VD" & year == 1935 & liste_norm =="LISTE LIBERALE DEMOCRATIQUE"
replace liste_norm ="LISTE DE L'ALLIANCE DES INDEPENDANTS" if canton == "VD" & year == 1943 & liste_norm =="LISTE DE L'ALLIANCE DES INDEPENDANTS VAUDOIS"
replace liste_norm ="LISTE NATIONALE PAYSANNE (LISTE DU PARTI NATIONAL DES PAYSANS ARTISANS ET BOURGEOIS)" if canton == "VD" & year == 1943 & liste_norm =="LISTE NATIONALE PAYSANNE"
replace liste_norm ="LISTE DA PARTI OUVRIER ET POPULAIRE VAUDOIS" if canton == "VD" & year == 1951 & liste_norm =="LISTE DU PARTI OUVRIER ET POPULAIRE VAUDOIS"
replace liste_norm ="PARTI SOCIALISTE VAUDOIS" if canton == "VD" & year == 1955 & liste_norm =="PARTI SOCIALISTE"
replace liste_norm ="PARTI RADICAL DEMOCRATIQUE" if canton == "VD" & year == 1963 & liste_norm =="PARTI RADICAL DEMOCRATIQUE VAUDOIS"
replace liste_norm ="PARTI VAUDOIS DES PAYSANS ARTISANS ET INDEPENDANTS" if canton == "VD" & year == 1963 & liste_norm =="PARTI VAUDOIS DES PAYSANS"
replace liste_norm ="PARTI LIBERALE DEMOCRATIQUE" if canton == "VD" & year == 1967 & liste_norm =="PARTI LIBERAL DEMOCRATIQUE"
replace liste_norm ="PARTI RADICALE DEMOCRATIQUE" if canton == "VD" & year == 1967 & liste_norm =="PARTI RADICAL DEMOCRATIQUE"
replace liste_norm ="PARTI VAUDOIS DES PAYSANS ARTISANS ET INDEPENDANTS (PAI)" if canton == "VD" & year == 1967 & liste_norm =="PARTI VAUDOIS DES PAYSANS ARTISANS ET INDEPENDANTS"
replace liste_norm ="MOUVEMEAT NATIONAL D'ACTION REPUBLICAINE ET SOCIALE" if canton == "VD" & year == 1971 & liste_norm =="MOUVEMENT NATIONAL D'ACTION REPUBLICAINE ET SOCIALE"
replace liste_norm ="LISTE OUVRIERE ET PAYSANNE" if canton == "VS" & year == 1935 & liste_norm =="LISTE OUVRIÈRE ET PAYSANNE"
replace liste_norm ="LISTE OUVRIERE ET PAYSANNE" if canton == "VS" & year == 1939 & liste_norm =="LISTE OUVRIÈRE ET PAYSANNE"
replace liste_norm ="CHRISTLICH SOZIALE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1955 & liste_norm =="LISTE DER CHRISTLICH SOZIALEN VOLKSPARTEI DES OBERWALLIS"
replace liste_norm ="KONSERVATIVE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1955 & liste_norm =="LISTE DER KONSERVATIVEN VOLKSPARTEI DES OBERWALLIS"
replace liste_norm ="MOUVEMENT SOCIAL PAYSAN INDEPENDANT" if canton == "VS" & year == 1955 & liste_norm =="LISTE DU MOUVEMENT SOCIAL PAYSAN INDEPENDANT"
replace liste_norm ="PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND" if canton == "VS" & year == 1955 & liste_norm =="LISTE DU PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND"
replace liste_norm ="PARTI RADICAL DEMOCRATIQUE" if canton == "VS" & year == 1955 & liste_norm =="LISTE DU PARTI RADICAL DEMOCRATIQUE"
replace liste_norm ="CHRISTLICHSOZIALE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1959 & liste_norm =="LISTE DER CHRISTLICH SOZIALEN VOLKSPARTEI OBERWALLIS"
replace liste_norm ="KONSERVATIVE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1959 & liste_norm =="LISTE DER KONSERVATIVEN VOLKSPARTEI OBERWALLIS"
replace liste_norm ="MOUVEMENT SOCIAL DES PAYSANS OUVRIERS ET INDEPENDANTS" if canton == "VS" & year == 1959 & liste_norm =="LISTE DU MOUVEMENT SOCIAL DES PAYSANS OUVRIERS ET INDEPENDANTS"
replace liste_norm ="PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND" if canton == "VS" & year == 1959 & liste_norm =="LISTE DU PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND"
replace liste_norm ="PARTI RADICAL DEMOCRATIQUE" if canton == "VS" & year == 1959 & liste_norm =="LISTE DU PARTI RADICAL DEMOCRATIQUE"
replace liste_norm ="PARTI SOCIALISTE" if canton == "VS" & year == 1959 & liste_norm =="LISTE DU PARTI SOCIALISTE"
replace liste_norm ="CHRISTLICHSOZIALE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1963 & liste_norm =="LISTE DER CHRISTLICH SOZIALEN VOLKSPARTEI OBERWALLIS"
replace liste_norm ="KONSERVATIVE VOLKSPARTEI OBERWALLIS" if canton == "VS" & year == 1963 & liste_norm =="LISTE DER KONSERVATIVEN VOLKSPARTEI OBERWALLIS"
replace liste_norm ="PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND" if canton == "VS" & year == 1963 & liste_norm =="LISTE DU PARTI CONSERVATEUR CHRETIEN SOCIAL DU VALAIS ROMAND"
replace liste_norm ="PARTI RADICAL DEMOCRATIQUE" if canton == "VS" & year == 1963 & liste_norm =="LISTE DU PARTI RADICAL DEMOCRATIQUE"
replace liste_norm ="PARTI SOCIALISTE" if canton == "VS" & year == 1963 & liste_norm =="LISTE SOCIALISTE  SOZIALDEMOKRATISCHE LISTE"
replace liste_norm ="ALLIANCE DES INDEPENDANTS" if canton == "VS" & year == 1967 & liste_norm =="ALLIANCE DES INDEPENDANTS  LANDESRING DER UNABHANGIGEN"
replace liste_norm ="LISTE SOCIALISTE" if canton == "VS" & year == 1967 & liste_norm =="LISTE SOCIALISTE  SOZIALDEMOKRATISCHE LISTE"
replace liste_norm ="SOCIALISTE POPULAIRE" if canton == "VS" & year == 1967 & liste_norm =="LISTE SOCIALISTE POPULAIRE  SOZIALISTISCHE VOLKSPARTEI"
replace liste_norm ="MOUVEMENT SOCIAL INDEPENDENT" if canton == "VS" & year == 1967 & liste_norm =="MOUVEMENT SOCIAL INDEPENDANT"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "ZG" & year == 1935 & liste_norm =="FREISINNIG DEMOKRATISCHE LISTE"
replace liste_norm ="KONSERVATIVE VOLKS UND ARBEITERPARTEI" if canton == "ZG" & year == 1935 & liste_norm =="KONSERVATIVE VOLKS UND ARBEITERLISTE"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "ZG" & year == 1935 & liste_norm =="SOZIALDEMOKRATISCHE LISTE"
replace liste_norm ="KONSERVATIVE VOLKS ARBEITERPARTEI" if canton == "ZG" & year == 1943 & liste_norm =="KONSERVATIVE VOLKS UND ARBEITERPARTEI"
replace liste_norm ="KONSERVATIV CHRISTLICHSOZIALE PARTEI" if canton == "ZG" & year == 1955 & liste_norm =="KONSERVATIV CHRISTLICHSOZIALE LISTE"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "ZG" & year == 1955 & liste_norm =="SOZIALDEMOKRATISCHE LISTE"
replace liste_norm ="CONSERVATIV CHRISTLICHSOZIALE PARTEI" if canton == "ZG" & year == 1959 & liste_norm =="CONSERVATIV CHRISTLICHSOZIALE LISTE"
replace liste_norm ="FREISINNIG DEMOKRATISCHE PARTEI" if canton == "ZG" & year == 1959 & liste_norm =="FREISINNIG DEMOKRATISCHE LISTE"
replace liste_norm ="SOZIALDEMOKRATISCHE PARTEI" if canton == "ZG" & year == 1959 & liste_norm =="SOZIALDEMOKRATISCHE LISTE"
replace liste_norm ="KONSERVATIV CHRISTLICHSOZIALE PARTEI" if canton == "ZG" & year == 1967 & liste_norm =="CONSERVATIV CHRISTLICHSOZIALE PARTEI"
replace liste_norm ="FREISINNIGE PARTEI" if canton == "ZH" & year == 1931 & liste_norm =="FREISINNIGE LISTE"
replace liste_norm ="LISTE EIDGENOSSISCHE FRONT" if canton == "ZH" & year == 1931 & liste_norm =="LISTE EIDGENÖSSISCHE FRONT"
replace liste_norm ="CHRISTLICHSOZIALE LISTE" if canton == "ZH" & year == 1935 & liste_norm =="CHRISTLICH SOZIALE LISTE"
replace liste_norm ="LISTE DER SCHWEIZERISCHEN BAUERNHEIMATBEWEGUNG (JUNGBAUERN)" if canton == "ZH" & year == 1935 & liste_norm =="LISTE DER SCHWEIZERISCHEN BAUERNHEIMATBEWEGUNG JUNGBAUERN"
replace liste_norm ="LISTE DER KANTONALEN BAUERNPARTEI BAUERLICH GEWERBLICHBURGERLICHE LISTE" if canton == "ZH" & year == 1939 & liste_norm =="LISTE DER KANTONALEN BAUERNPARTEI BAUERLICH GEWERBLICH BURGERLICHE LISTE"
replace liste_norm ="SOZIALDEMOKRATISCHELISTE" if canton == "ZH" & year == 1939 & liste_norm =="SOZIALDEMOKRATISCHE LISTE"
replace liste_norm ="CHRISTLICHSOZIALE LISTE" if canton == "ZH" & year == 1943 & liste_norm =="LISTE DER CHRISTLICH SOZIALEN PARTEI"
replace liste_norm ="LISTE DER SCHWEIZ BAUERN HEIMATBEWEGUNG" if canton == "ZH" & year == 1943 & liste_norm =="LISTE DER SCHWEIZERISCHEN BAUERNHEIMATBEWEGUNG"
replace liste_norm ="EVANGELISCHEN VOLKSPARTEI" if canton == "ZH" & year == 1951 & liste_norm =="EVANGELISCHE VOLKSPARTEI"
replace liste_norm ="FREISINNIGE LISTE ZURICH STADT" if canton == "ZH" & year == 1951 & liste_norm =="PREISINNIGE LISTE ZURICH STADT"
replace liste_norm ="LISTE DER BAUERN  GEWERBE UND BURGERPARTEI" if canton == "ZH" & year == 1955 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI"
replace liste_norm ="LISTE DER EVANGELISCHEN VOLKSPARTEI" if canton == "ZH" & year == 1955 & liste_norm =="EVANGELISCHE VOLKSPARTEI"
replace liste_norm ="LISTE DER PARTEI DER ARBEIT" if canton == "ZH" & year == 1955 & liste_norm =="PARTEI DER ARBEIT"
replace liste_norm ="FREISINNIGE LISTE ZURICH LAND" if canton == "ZH" & year == 1955 & liste_norm =="PREISINNIGE LISTE ZURICH LAND"
replace liste_norm ="LISTE EVANGELISCHE VOLKSPARTEI" if canton == "ZH" & year == 1959 & liste_norm =="EVANGELISCHE VOLKSPARTEI"
replace liste_norm ="FREISINNIGE LISTE STADT ZURICH" if canton == "ZH" & year == 1963 & liste_norm =="FREISINNIGE LISTE ZURICH STADT"
replace liste_norm ="LISTE DER BAUERN  GEWERBE UND BURGERPARTEI (MITTELSTANDSLISTE)" if canton == "ZH" & year == 1963 & liste_norm =="LISTE DER BAUERN  GEWERBE UND BURGERPARTEI (MITTELSTANDS LISTE)"
replace liste_norm ="LISTE 13 LISTE DER UBERPARTEILICHEN UNION" if canton == "ZH" & year == 1963 & liste_norm =="LISTE DER UBERPARTEILICHEN UNION"
replace liste_norm ="BGB MITTELSTANDSLISTE ZURICH LAND" if canton == "ZH" & year == 1967 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI ZURICH LAND (MITTELSTANDSLISTE)"
replace liste_norm ="LISTE DER BAUERN  GEWERBE UND BURGERPARTEI ZURICH STADT (MITTELSTANDSLISTE)" if canton == "ZH" & year == 1967 & liste_norm =="BAUERN  GEWERBE UND BURGERPARTEI ZURICH STADT (MITTELSTANDSLISTE)"
replace liste_norm ="LISTE DER SOZIALDEMOKRATEN GEWERKSCHAFTER UND ANGESTELLTEN" if canton == "ZH" & year == 1967 & liste_norm =="SOZIALDEMOKRATEN GEWERKSCHAFTER UND ANGESTELLTE"
replace liste_norm ="ERWA BUND" if canton == "ZH" & year == 1971 & liste_norm =="ERWA BUND (KAMPF FUR RECHT UND UMWELTSCHUTZ)"
replace liste_norm ="LANDSEKTIONEN DER NATIONALEN AKTION GEGEN DIE UBERFREMDUNG VON VOLK UND HEIMAT" if canton == "ZH" & year == 1971 & liste_norm =="LANDSEKTION DER NATIONALEN AKTION GEGEN DIE UBERFREMDUNG VON VOLK UND HEIMAT"
replace liste_norm ="BGB MITTELSTANDSPARTEI" if canton == "ZH" & year == 1971 & liste_norm =="LISTE DER BGB MITTELSTANDSPARTEI"
replace liste_norm ="EVANGELISCHE VOLKSPARTEI (EVP)" if canton == "ZH" & year == 1971 & liste_norm =="LISTE EVANGELISCHE VOLKSPARTEI (EVP)"
replace liste_norm ="JUNGE MITTE" if canton == "ZH" & year == 1971 & liste_norm =="LISTE JUNGE MITTE"
replace liste_norm = "KEIN LISTENNAME" if Listenname == ""
replace liste_norm = "KEIN LISTENNAME" if liste_norm == "KEINE LISTENNAMEN"

replace liste_norm="Wagner" if pvotes==2256 & year==1943 ///
	& canton=="NW"
* Notes: Replace missing list information with candidate names for four cases.
replace liste_norm="Deschwanden" if pvotes==1506 & year==1943 ///
	& canton=="NW"
replace liste_norm="Odermatt" if pvotes==2421 & year==1943 ///
	& canton=="OW"
replace liste_norm="Infanger" if pvotes==1642 & year==1943 ///
	& canton=="OW"

	
keep canton year liste_norm alliance suballiance pvotes seats_cant eligible_cant voters_cant invalid_cant empty_cant valid_cant
order year canton liste_norm alliance suballiance pvotes seats_cant eligible_cant voters_cant invalid_cant empty_cant valid_cant
sort canton year
save "$path\02_Processed_data\05_Alliances\Listenverbindungen_1931-2015.dta", replace


erase "$path\02_Processed_data\05_Alliances\Listenverbindungen_1931-2015_temp1.dta"
erase "$path\02_Processed_data\05_Alliances\Listenverbindungen_1931-2015_temp2.dta"
erase "$path\02_Processed_data\05_Alliances\BdBlatt_Listenverbindungen_1931-1967.dta"
erase "$path\02_Processed_data\05_Alliances\BfS_Listenverbindungen_1971-2015.dta"
erase "$path\02_Processed_data\05_Alliances\NRWahl_ZusatzInfo_BdBlatt_1931-1967.dta"
erase "$path\02_Processed_data\05_Alliances\eligible_1931-2015.dta"
erase "$path\02_Processed_data\05_Alliances\votes_1971-2015.dta"

