clear

import excel "..\dta\Retriever.xlsx", firstrow
drop if Year==2016
drop if Id==.

/* monthly */
/* 
foreach var in Eiendomsskatt Kommunaleavgifter Barnehage Barnehagedekning Skole Grunnskole Eldreomsorg Sosialhjelp Barnevern Helsevesen KulturANDIdrett Samferdsel Kommunaladministrasjon {
scatter `var' Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(9 57 105, lpattern(dash)) xscale(range(25 110)) ytitle("") xtitle("") xlabel(none) subtitle(`var') connect(l)
graph save ..\figures\gph\_`var'.gph, replace
}
*/


/* quarterly */
collapse (sum) Eiendomsskatt Kommunaleavgifter Barnehage Barnehagedekning Skole Grunnskole Eldreomsorg Sosialhjelp Barnevern Helsevesen KulturANDIdrett Samferdsel Kommunaladministrasjon, by(Id)


/* REGRESSION ANALYSIS */
gen LocalElection=0
replace LocalElection=1 if Month==8 & Year==2007
replace LocalElection=1 if Month==8 & Year==2011
replace LocalElection=1 if Month==8 & Year==2015
replace LocalElection=1 if Month==9 & Year==2007
replace LocalElection=1 if Month==9 & Year==2011
replace LocalElection=1 if Month==9 & Year==2015

gen NationalElection=0
replace NationalElection=1 if Month==8 & Year==2005
replace NationalElection=1 if Month==8 & Year==2009
replace NationalElection=1 if Month==8 & Year==2013

tab Month, gen(d_Month)

log using ..\log\RetrieverRegressionOutput.txt, text replace
foreach var in Eiendomsskatt Kommunaleavgifter Barnehagedekning Grunnskole Eldreomsorg Sosialhjelp KulturANDIdrett Samferdsel Kommunaladministrasjon {
reg `var' LocalElection NationalElection d_Month* Id
}
log close
