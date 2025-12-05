
******************************************************
* (0) Set directories
******************************************************

clear
clear matrix
set mem 500m
set more off

capture cd "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
	if _rc==0{
		global data "C:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data"
	}
	

* (i) Start with AG
	
set more off			
cd "$data\1_data\2 fuzzy merge back in\endgültig fertig"
import excel  "erf_AG.xlsx"	, sheet(Sheet1) clear first
tostring id_recoded, replace
tostring Kommentar, replace
tostring unsicher, replace
save "$data\1_data\data_federal_temp.dta", replace	


* (ii) Merge all other files in folder endgültig fertig

cd "$data\1_data\2 fuzzy merge back in\endgültig fertig"
quietly fs *		// ssc install fs
display `r(files)'

foreach datei in `r(files)'{ // all cantons except ZH
//display  "`datei'"
cd "$data\1_data\2 fuzzy merge back in\endgültig fertig"
import excel  `datei'	, sheet(Sheet1) clear first
tostring id_recoded, replace
tostring Kommentar, replace
tostring unsicher, replace
cd "$data\1_data"
append using data_federal_temp.dta
save "data_federal_temp.dta", replace	
}


cd "$data\1_data\2 fuzzy merge back in\unsicher einfügen"
import excel  erf_ZH.xlsx	, sheet(Sheet1) clear first
cd "$data\1_data"
append using data_federal_temp.dta, force
save "$data\1_data\data_federal_temp.dta", replace	


* (iii) Data cleansing and variable generation

 // drop canton AG (start canton in loop above)
duplicates drop Nachname Vorname	Wahljahr Geschlecht	Geburtsjahr	Listenbezeichnung	Gemeinde, force 	


* (iv) ID recoding

gen id_Stata_new=id_Stata
replace id_Stata_new=id_recoded if id_recoded!="" & id_recoded!="."

* (v) Read out "unsicher" cases

preserve
keep if unsicher==.
cd "$data\1_data"
export excel "$data\1_data\unsicher_cases_check.xlsx", replace	
restore




* (vi) Merge additional information: Total list votes, individual votes and municipality

** Procedure: Prepare list data (in (a)) and municipality information (in (b))
* 			  Merge them with individual vote results in (c) and with ID-recoding in (d)

* (a) Dataset on votes per list

preserve
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(1971) cellrange(A5:L166)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes1971.dta", replace
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(1975) cellrange(A5:L187)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes1975.dta", replace
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(1979) cellrange(A5:L182)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes1979.dta", replace
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(1983) cellrange(A5:L202)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes1983.dta", replace
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(1987) cellrange(A5:L241)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes1987.dta", replace
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(1991) cellrange(A5:L266)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes1991.dta", replace
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(1995) cellrange(A5:L298)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes1995.dta", replace
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(1999) cellrange(A5:L285)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes1999.dta", replace
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(2003) cellrange(A5:L289)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes2003.dta", replace
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(2007) cellrange(A5:L335)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes2007.dta", replace
clear
import excel using "$data\0_original_data\NRW 1971_2011ErhalteneStimmen_StärkeWahllisten_nach Kantonen (BfS).xls", first sheet(2011) cellrange(A5:L393)
tostring ListenNroffiziell, replace
save "$data\1_data\partyvotes2011.dta", replace


use "$data\1_data\partyvotes1971.dta", clear
forvalues i=1975(4)2011{
append using "$data\1_data\partyvotes`i'.dta"
erase "$data\1_data\partyvotes`i'.dta"
}


gen Nationalratswahl=Wahljahr
save "$data\1_data\partyvotes_all.dta", replace
restore
 
 
 
* (b) Municipality dataset

preserve
clear
import excel using "$data\0_original_data\NRW_KANDIDATEN.xlsx", first
//insheet using "$data\0_original_data\NRW_KANDIDATEN.csv", delimiter(";")
rename Gemeindename Gemeinde
// rename candidates with equal characteristics in terms of "vorname name kantonsname listennummer   nationalratswahl gewaehlt"
replace Nachname="Bärtschi 2" if Parteinummer==4 & Gemeindenummer_BFS=="538" & Kandidatennummer== 9  & Kantonsnummer==2 &  Nationalratswahl=="NRW1979" // Bärtschi Jakob
replace Listennummer="1011" if Parteinummer==31 & Kandidatennummer== 15 & Kantonsnummer==1 &  Nationalratswahl=="NRW2007" // Stalder Martin
save "$data\1_data\candidates_municipalities.dta", replace
restore




* (c) Dataset on individual vote counts (with merged municipality  and list information)

preserve
clear
import excel using "$data\0_original_data\KANDIDATENSTIMMEN (BfS).xlsx", first
gen kantonsname=Kanton
gen listennummer=Listen_Nr_offiziell
gen kandidatennummer=Kandidaten_nr
gen parteinummer=Partei_nr
gen kantonsnummer=Kantons_nr
// rename candidates with equal characteristics in terms of "vorname name kantonsname listennummer   nationalratswahl gewaehlt"
replace Name="Bärtschi 2" if parteinummer==4 & kandidatennummer== 9  & KandidSt_ver_WZ==21092 & kantonsnummer==2 &  Nationalratswahl=="NRW1979" // Bärtschi Jakob
replace listennummer="1011" if parteinummer==31 & kandidatennummer== 15 & kantonsnummer==1 &  Nationalratswahl=="NRW2007" // Stalder Martin

rename Name Nachname 
gen Wahljahr=substr(Nationalratswahl,4,4)
destring Wahljahr, replace
gen Kantonsname=Kanton
rename Listen_Nr_offiziell Listennummer

// temporary: drop two candidates with the exact same name on the ballot (Martin Stalder, GLP ZH in 2007)

bysort Vorname Nachname Kantonsname Listennummer   Nationalratswahl: gen indimax=_N
drop if indimax>1

merge 1:1 Vorname Nachname Kantonsname Listennummer   Nationalratswahl  using "$data\1_data\candidates_municipalities.dta", gen(merge2)
keep if merge2==3
/* NOTE: 1749 OBSERVATIONS DO NOT MERGE (1689 OF THEM ARE OBSERVATIONS FROM 1971) -> WE SHOULD CHECK THIS IN A LATER STAGE OF THE PROJECT
    Result                           # of obs.
    -----------------------------------------
    not matched                         1,749
        from master                        54  (merge2==1)
        from using                      1,695  (merge2==2)

    matched                            25,806  (merge2==3)
    -----------------------------------------
*/

gen ListenNroffiziell=Listennummer

merge m:1 ListenNroffiziell  Nationalratswahl Kanton using "$data\1_data\partyvotes_all.dta", gen(merge3) force
keep if merge3==3

/* NOTE: SIMILAR PROBLEM AS ABOVE WITH 1971 ELECTION DATA

    Result                           # of obs.
    -----------------------------------------
    not matched                         1,844
        from master                     1,695  (merge3==1)
        from using                        149  (merge3==2)

    matched                            25,860  (merge3==3)
    -----------------------------------------
*/

drop if Vorname=="" | Vorname=="Vereinzelte"| Vorname=="vereinzelte"

save candidates_votes_and_municipalities.dta, replace
restore

* (d) Merge municipality, list votes and individuals votes to ID data

gen Liste=Listenbezeichnung

merge 1:1 Vorname Nachname  Geburtsjahr Wahljahr Gemeinde Liste using "$data\1_data\candidates_votes_and_municipalities.dta", gen(merge4) force


* (e) Generate running variable and other variables


gen canton_id=Kanton
gen canton_id_num=Kanton

sort canton_id listennummer
gen list_nr=ListenNrnumerisch
gen year=Wahljahr

gen elected=.
replace elected=1 if Gewaehlt=="G"
replace elected=0 if Gewaehlt=="N"

gen individual_votes=KandidSt_tot
gen list_votes=ErhalteneStimmen

egen time_variable=group(Wahljahr)
egen id_Stata_new_num=group(id_Stata_new)

// temporary: drop duplicates
bysort id_Stata_new_num time_variable: gen indimax2=_N
drop if indimax2>1

xtset id_Stata_new_num time_variable
gen elected_lag=L1.elected


gen birthyear=Geburtsjahr

gen age=year-birthyear
gen age_sq=age^2
gen sex=.
replace sex=1 if Geschlecht=="M"
replace sex=0 if Geschlecht=="F"

// Generate running variable

bysort id_Stata_new_num: gen nr_participation=_n

gsort year canton_id  list_nr -individual_votes
bysort year canton_id  list_nr: gen candidate_rank=_n // rank of candidate on list


bysort year canton_id  list_nr: egen elected_perlist= sum(elected) // indicator if anybody is elected on list
bysort year canton_id : egen elected_ct= sum(elected) 
		
sort  year canton_id  list_nr candidate_rank
gen elect_rank=elected*candidate_rank		// interaction rank and elected
bysort year canton_id  list_nr: egen elect_rank_max=max(elect_rank)
gen candidate_dist=candidate_rank-elect_rank_max // distance to marginal guy in positions

gen marginal=0 // marginal guy from below (for all non-elected)
replace marginal=1 if candidate_dist==0

gen marginal_above=0 // marginal guy from above (for all elected)
replace marginal_above=1 if candidate_dist==1


gen candidate_votes_perc= individual_votes/list_votes*100 // generate individual vote share in terms of total list votes
sort id_Stata_new_num time_variable
gen candidate_votes_perc_lag=L1.candidate_votes_perc



gen candidate_votes_perc_marg=. 
replace candidate_votes_perc_marg=candidate_votes_perc if marginal==1

gen candidate_votes_perc_marg_above=. 
replace candidate_votes_perc_marg_above=candidate_votes_perc if marginal_above==1


bysort year canton_id  list_nr:egen  candidate_votes_perc_marg_all=max(candidate_votes_perc_marg) 
bysort year canton_id  list_nr:egen  candidate_votes_perc_mar_ab_all=min(candidate_votes_perc_marg_above) 

gen candidate_diff_marginal=candidate_votes_perc-candidate_votes_perc_marg_all if elected==0
replace candidate_diff_marginal=candidate_votes_perc-candidate_votes_perc_mar_ab_all if elected==1

sort id_Stata_new_num time_variable
gen candidate_diff_marginal_lag=L1.candidate_diff_marginal

gen candidate_diff_marginal_lag_elec=candidate_diff_marginal_lag*elected

saveold nr_analysis.dta, replace version(12)




///////// TEMPORARY: FIRST REGRESSION ANALYSIS
/////////         (for cantonal tax administrations and data protection officers)

cd"$data\1_data"
use nr_analysis.dta, clear

preserve
keep if canton_id=="BE"
keep if inrange(year,1999,2011)


reg candidate_votes_perc  elected_lag candidate_diff_marginal_lag candidate_diff_marginal_lag_elec i.year, cl(id_Stata_new_num)
	outreg2 using overall1.xls , excel  stats(coef se  ) noaster cttop("") dec(2) replace	
											
reg candidate_votes_perc  elected_lag  age sex  candidate_diff_marginal_lag candidate_diff_marginal_lag_elec i.year, cl(id_Stata_new_num)
	outreg2 using overall1.xls , excel  stats(coef se  ) noaster cttop("") dec(2) append
												
xtreg candidate_votes_perc  elected_lag age sex candidate_diff_marginal_lag candidate_diff_marginal_lag_elec i.year, cl(id_Stata_new_num)
	outreg2 using overall1.xls , excel  stats(coef se  ) noaster cttop("") dec(2) append
	
xtreg candidate_votes_perc  elected_lag age sex candidate_diff_marginal_lag candidate_diff_marginal_lag_elec i.year i. ParteiNr, cl(id_Stata_new_num)
	outreg2 using overall1.xls , excel  stats(coef se  ) noaster cttop("") dec(2) append	

sum candidate_votes_perc if e(sample)==1 & elected==0
	
	
xtreg candidate_votes_perc  elected_lag age sex candidate_diff_marginal_lag candidate_diff_marginal_lag_elec i.year i. ParteiNr, cl(id_Stata_new_num), ///
		if inrange(candidate_diff_marginal_lag,-2.5,2.5)
	outreg2 using overall1.xls , excel  stats(coef se  ) noaster cttop("") dec(2) append	

sum candidate_votes_perc if e(sample)==1 & elected==0 & inrange(candidate_diff_marginal_lag,-2.5,2.5)


restore
