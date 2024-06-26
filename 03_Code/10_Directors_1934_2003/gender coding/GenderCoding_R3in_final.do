version 17
clear all
set more off
cap log close

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\"


*****************
* Gender Coding *
*****************

* Read in first (R1) and third (R3) round of gender coding (second round contained in R3)
* Merge R1 and R3 (R2 included in R3)
* Checks with BfS data

* Read in R3 
import excel using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R3_in.xlsx", firstrow cellrange(A1:AI1630) clear
save "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R3.dta", replace

gen Comment_R2 = Comment_R2_1 + "; " + Comment_R2_2
replace Comment_R2 = "" if Comment_R2 == "; "
drop Comment_R2_1 Comment_R2_2
gen R = 3
keep cleanname male_R3 undecided_R3 dontknow_R3 Comment_R2 Comment_R3
sort cleanname

* Read in R1
preserve
use "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R1.dta", clear
gen male_R1 = male_R1_1 if Diff_R1 != 1		
gen undecided_R1 = undecided_R1_1 if Diff_R1 != 1	
gen dontknow_R1 = dontknow_R1_1 if Diff_R1 != 1	
gen Comment_R1 = Comment_1 + "; " + Comment_2
replace Comment_R1 = "" if Comment_R1 == "; "
keep nobs cleanname male_R1 undecided_R1 dontknow_R1 Comment_R1 Diff_R1
gen R = 1
replace R = . if nobs > -4
sort cleanname
save "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R1_ok.dta", replace
restore 

merge 1:1 cleanname using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecode_R1_ok.dta"

gen gender = 1 if male_R1 == 1 | male_R3 == 1
replace gender = 2 if male_R1 == 2 | male_R3 == 2
label var gender "male = 1; female = 2"

gen undecided = 1 if undecided_R1 == 1 | undecided_R3 == 1
gen dontknow = 1 if dontknow_R1 == 1 | dontknow_R3 == 1

label var undecided "male and female first names"
label var dontknow "unclear gender from first names"
label var R "1: R1 manual coding, 2: R2+R3 manual coding, 4: BfS"

gen Comments = Comment_R1 + "; " + Comment_R2 + "; " + Comment_R3
replace Comments = "" if Comments == "; " | Comments == "; ; "

keep nobs cleanname gender undecided dontknow Comments R


* Compare to BfS names
sort cleanname 
merge 1:1 cleanname using "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\BfS_Vornamen2015-19.dta"
drop if _merge == 2
drop _merge

egen totBfS = rowtotal(maleBfS femaleBfS)
replace totBfS = . if totBfS == 0
gen mbalBfS = maleBfS / totBfS
gen genderBfS = 1 if mbalBfS >= 0.9 & mbalBfS != .  // 90% cut-off
replace genderBfS = 2 if mbalBfS <= 0.1 & mbalBfS != .
gen undecidedBfS = 1 if mbalBfS != . & (mbalBfS > 0.1 & mbalBfS < 0.9) 

gen Flag_BfS = 1 if gender != genderBfS	& genderBfS != . 	// flag difference in coding and BfS
replace Flag_BfS = 1 if  undecided != undecidedBfS & undecidedBfS != . 
*br * if Flag_BfS == 1

* Check coding quality wrt to BfS first name count 2015 - 2019
preserve  // 90% cut-off for BfS gender coding
drop if R == .
gen corr = 0
replace corr = 1 if gender == genderBfS & genderBfS != .
replace corr = 1 if undecided == undecidedBfS & undecidedBfS != . 
replace corr = 1 if dontknow == undecidedBfS & undecidedBfS != .
replace corr = . if genderBfS == . & undecidedBfS == .
drop if corr == .
egen pctcorr = mean(corr)
di pctcorr " correctly coded if 1) name contained in BfS and 2) 90% cut-off for gender coding is applied"
restore

preserve // gender coded according to majority in BfS?
drop if R == .
replace genderBfS = 1 if maleBfS > femaleBfS
replace genderBfS = 2 if femaleBfS > maleBfS
gen corr = 0
replace corr = 1 if gender == genderBfS & genderBfS != .
replace corr = 1 if undecided == undecidedBfS & undecidedBfS != .
replace corr = 1 if dontknow == undecidedBfS & undecidedBfS != .
replace corr = . if genderBfS == . & undecidedBfS == .
drop if corr == .
egen pctcorr = mean(corr)
di pctcorr " correctly coded if 1) name contained in BfS and 2) gender coding is based on majority cut-off"
restore

* Corrections based on BfS counts
replace gender = genderBfS if Flag_BfS == 1
replace undecided = undecidedBfS if Flag_BfS == 1
replace R = 4 if Flag_BfS == 1 

preserve // % coded genders
gen coded = -1 * nobs if R != .
gen mancoded = -1 * nobs if nobs > -4
egen tmancoded = total(mancoded)
gen err = -1 * nobs if Flag_BfS == 1 & nobs > -4
egen pcterr = total(err)
replace pcterr = pcterr / tmancoded * 100
gen gen = -1 * nobs if gender == 1 | gender == 2
gen undec = -1 * nobs if undecided == 1
gen nBfS = -1 * nobs if R == 4
egen tot = total(nobs)
replace tot = tot * -1
egen pctcoded = total(coded) 
replace pctcoded = (pctcoded) / tot * 100
egen pctgen = total(gen)
replace pctgen = pctgen / tot * 100
egen pctundec = total(undec)
replace pctundec = pctundec / tot* 100
egen pctBfS = total(nBfS)
replace pctBfS = pctBfS / tot * 100
di "% coded names:            " pctcoded
di "% known genders:          " pctgen
di "% undecided genders:      " pctundec
di "% codings based on BfS:   " pctBfS
di "% manual coding errors    " pcterr
restore

* Coding Harmonisierung

replace gender = 0 if gender == 2
label var gender "0: female, 1: male"
drop if cleanname == ""
keep cleanname gender

sort cleanname
save "$path\02_Processed_data\10_Directors_1934_2003\07_GenderCoding\Gender_ToRecoded_final.dta", replace
