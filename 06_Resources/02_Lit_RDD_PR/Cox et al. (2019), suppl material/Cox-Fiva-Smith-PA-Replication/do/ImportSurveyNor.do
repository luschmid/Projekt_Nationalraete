******************************************************************************************************************************************************
********************* NORWAY SURVEY ******************************************************************************************************************
******************************************************************************************************************************************************
use dta/1965.dta, clear 
rename v002 knr

gen visitNKP=0
replace visitNKP=1 if v148==1
replace visitNKP=1 if v148==2

gen visitSF=0
replace visitSF=1 if v149==1
replace visitSF=1 if v149==2

gen visitDNA=0
replace visitDNA=1 if v150==1
replace visitDNA=1 if v150==2

gen visitV=0
replace visitV=1 if v151==1
replace visitV=1 if v151==2

gen visitKRF=0
replace visitKRF=1 if v152==1
replace visitKRF=1 if v152==2

gen visitSP=0
replace visitSP=1 if v153==1
replace visitSP=1 if v153==2

gen visitH=0
replace visitH=1 if v154==1
replace visitH=1 if v154==2

egen visitTOT=rowtotal(visitNKP-visitH)
gen visitANY=visitTOT>0
gen cnr=floor(knr/100)

sort cnr
merge cnr using dta/1965_magnitude.dta
drop if _merge==2

rename v266 age
tab v267, gen(d_married)
gen married=0
replace married=1 if v267==1
gen female=0
replace female=1 if v265==1
tab v288, gen(d_educ)
keep magnitude visitANY cnr knr d_* married age female
gen year=1965
save dta/1965_temp.dta, replace

******************'
use dta/1969.dta, clear

gen female=0
replace female=1 if v379==1
tab v384, gen(d_educ)
gen temp=v380
replace temp=1 if v380==2 /* code married as 1, as in 1965 */
replace temp=2 if v380==1 /* code unmarried as 2, as in 1965 */
tab temp, gen(d_married)
gen married=0
replace married=1 if v380==1  /* notice discrepancy in coding 1965, 1969 */

gen age=1969-omk4

gen visitANY=0
replace visitANY=1 if v222==1
replace visitANY=1 if v222==2
replace visitANY=. if v222==9
rename v002 knr

gen cnr=floor(knr/100)

sort cnr
merge cnr using dta/1965_magnitude.dta
drop if _merge==2

gen year=1969
keep year magnitude visitANY cnr knr d_* age female
append using dta/1965_temp.dta
tab cnr
tab year, gen(d_year)
gen age2=age^2

save dta/ImportSurveyNor.dta, replace
