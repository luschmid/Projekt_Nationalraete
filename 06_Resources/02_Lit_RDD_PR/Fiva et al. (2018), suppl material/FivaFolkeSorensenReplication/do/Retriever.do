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

*foreach var in Eiendomsskatt Kommunaleavgifter Barnehage Barnehagedekning Skole Grunnskole Eldreomsorg Sosialhjelp Barnevern Helsevesen KulturANDIdrett Samferdsel Kommunaladministrasjon {
*scatter `var' Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(3 19 35, lpattern(dash)) xlabel(3 "Q3-07" 19 "Q3-11" 35 "Q3-15")  xscale(range(1 36)) ytitle("") xtitle("") subtitle(`var') connect(l) scale(0.85)
*graph save ..\figures\gph\_`var'.gph, replace
*}

scatter Eiendomsskatt Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(3 19 35, lpattern(dash)) xscale(off) xmtick(,nolabels)  xscale(range(1 36)) ytitle("") xtitle("") subtitle(Property tax) connect(l) scale(0.85) graphregion(fcolor(white)) scheme(s2mono)
graph save ..\figures\gph\_Eiendomsskatt.gph, replace

scatter Kommunaleavgifter Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(3 19 35, lpattern(dash)) xscale(off) xmtick(,nolabels)  xscale(range(1 36)) ytitle("") xtitle("") subtitle(User charges) connect(l) scale(0.85) graphregion(fcolor(white)) scheme(s2mono)
graph save ..\figures\gph\_Kommunaleavgifter.gph, replace

scatter Barnehagedekning Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(3 19 35, lpattern(dash)) xscale(off) xmtick(,nolabels)  xscale(range(1 36)) ytitle("") xtitle("") subtitle(Child care) connect(l) scale(0.85) graphregion(fcolor(white)) scheme(s2mono)
graph save ..\figures\gph\_Barnehagedekning.gph, replace

scatter Grunnskole Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(3 19 35, lpattern(dash)) xscale(off) xmtick(,nolabels)  xscale(range(1 36)) ytitle("") xtitle("") subtitle(Education) connect(l) scale(0.85) graphregion(fcolor(white)) scheme(s2mono)
graph save ..\figures\gph\_Grunnskole.gph, replace

scatter Eldreomsorg Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(3 19 35, lpattern(dash)) xscale(off) xmtick(,nolabels)  xscale(range(1 36)) ytitle("") xtitle("") subtitle(Elderly care) connect(l) scale(0.85) graphregion(fcolor(white)) scheme(s2mono)
graph save ..\figures\gph\_Eldreomsorg.gph, replace

scatter Sosialhjelp Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(3 19 35, lpattern(dash)) xscale(off) xmtick(,nolabels)  xscale(range(1 36)) ytitle("") xtitle("") subtitle(Health care) connect(l) scale(0.85) graphregion(fcolor(white)) scheme(s2mono)
graph save ..\figures\gph\_Sosialhjelp.gph, replace

scatter KulturANDIdrett Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(3 19 35, lpattern(dash)) xlabel(3 "Q3-07" 19 "Q3-11" 35 "Q3-15")  xscale(range(1 36)) ytitle("") xtitle("") subtitle(Culture) connect(l) scale(0.85) graphregion(fcolor(white)) scheme(s2mono)
graph save ..\figures\gph\_KulturANDIdrett.gph, replace

scatter Samferdsel Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(3 19 35, lpattern(dash)) xlabel(3 "Q3-07" 19 "Q3-11" 35 "Q3-15")  xscale(range(1 36)) ytitle("") xtitle("") subtitle(Transport) connect(l) scale(0.85) graphregion(fcolor(white)) scheme(s2mono)
graph save ..\figures\gph\_Samferdsel.gph, replace

scatter Kommunaladministrasjon Id, mstyle(Oh) msymbol(circle) mcolor(white) mlcolor(black) msize(medlarge) xline(3 19 35, lpattern(dash)) xlabel(3 "Q3-07" 19 "Q3-11" 35 "Q3-15")  xscale(range(1 36)) ytitle("") xtitle("") subtitle(Administration) connect(l) scale(0.85) graphregion(fcolor(white)) scheme(s2mono)
graph save ..\figures\gph\_Kommunaladministrasjon.gph, replace

cd ..\figures\gph
graph combine _Eiendomsskatt.gph _Kommunaleavgifter.gph _Barnehagedekning.gph _Grunnskole.gph _Eldreomsorg.gph _Sosialhjelp.gph _KulturANDIdrett.gph _Samferdsel.gph _Kommunaladministrasjon.gph, rows(3) cols(3) graphregion(fcolor(white)) scheme(s2mono)
graph save Figure4.gph, replace
graph export ..\Figure4.tif, replace
graph export ..\Figure4.eps, replace
cd ..\..\do\
