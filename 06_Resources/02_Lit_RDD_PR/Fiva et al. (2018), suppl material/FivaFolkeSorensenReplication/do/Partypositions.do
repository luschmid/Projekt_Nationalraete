use ..\dta\PolicyIndexLongAccurate, clear

cd ..\figures\gph

rename DPTAX PTAXpref
rename Fees FEEpref
rename centraladm Administration
rename education Education
rename childcare ChildCare
rename elderlycare ElderlyCare
rename healthsocial HealthCare
rename transport Transportation
rename culture Culture
rename other Other

foreach depvar in PTAXpref {
twoway (scatter   `depvar' LeftRight if Party=="DNA", mcolor(black) msize(15) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(0,1)) ylabel(0(0.2)1)) ///
(scatter   `depvar' LeftRight if Party=="SV", mcolor(black) msize(7) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(0,1)) ylabel(0(0.2)1)) ///
(scatter   `depvar' LeftRight if Party=="RV", mcolor(black) msize(2) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(0,1)) ylabel(0(0.2)1)) ///
(scatter   `depvar' LeftRight if Party=="PP", mcolor(gray*0.55) msize(2) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(0,1)) ylabel(0(0.2)1)) ///
(scatter   `depvar' LeftRight if Party=="V", mcolor(gray*0.55) msize(5.5) mlabel(Party) mlabposition(6) mlabcolor(black) msymbol(oh) yscale(r(0,1)) ylabel(0(0.2)1)) ///
(scatter   `depvar' LeftRight if Party=="SP", mcolor(gray*0.55) msize(8) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(0,1)) ylabel(0(0.2)1)) ///
(scatter   `depvar' LeftRight if Party=="KRF", mcolor(gray*0.55) msize(7.5) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(0,1)) ylabel(0(0.2)1)) ///
(scatter   `depvar' LeftRight if Party=="H", mcolor(gray*0.55) msize(10.5) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(0,1)) ylabel(0(0.2)1)) ///
(scatter   `depvar' LeftRight if Party=="FRP", mcolor(gray*0.55) msize(8.5) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(0,1)) ylabel(0(0.2)1)), legend(off) graphregion(fcolor(white)) scheme(s2mono) saving(`depvar'LR, replace) ytitle("") xtitle("") subtitle(`depvar') xscale(off) xmtick(,nolabels)
}

foreach depvar in FEEpref ElderlyCare Education ChildCare HealthCare {
twoway (scatter   `depvar' LeftRight if Party=="DNA", mcolor(black) msize(15) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="SV", mcolor(black) msize(7) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="RV", mcolor(black) msize(2) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="PP", mcolor(gray*0.55) msize(2) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="V", mcolor(gray*0.55) msize(5.5) mlabel(Party) mlabposition(6) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="SP", mcolor(gray*0.55) msize(8) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="KRF", mcolor(gray*0.55) msize(7.5) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="H", mcolor(gray*0.55) msize(10.5) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="FRP", mcolor(gray*0.55) msize(8.5) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)), legend(off) graphregion(fcolor(white)) scheme(s2mono) saving(`depvar'LR, replace) ytitle("") xtitle("") subtitle(`depvar') xscale(off) xmtick(,nolabels)
}

foreach depvar in Administration Transportation Culture Other {
twoway (scatter   `depvar' LeftRight if Party=="DNA", mcolor(black) msize(15) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="SV", mcolor(black) msize(7) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="RV", mcolor(black) msize(2) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="PP", mcolor(gray*0.55) msize(2) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="V", mcolor(gray*0.55) msize(5.5) mlabel(Party) mlabposition(6) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="SP", mcolor(gray*0.55) msize(8) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="KRF", mcolor(gray*0.55) msize(7.5) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="H", mcolor(gray*0.55) msize(10.5) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)) ///
(scatter   `depvar' LeftRight if Party=="FRP", mcolor(gray*0.55) msize(8.5) mlabel(Party) mlabposition(12) mlabcolor(black) msymbol(oh) yscale(r(-0.6,0.6)) ylabel(-0.6(0.2)0.6)), legend(off) graphregion(fcolor(white)) scheme(s2mono) saving(`depvar'LR, replace) ytitle("") xtitle("") subtitle(`depvar')
}


graph combine PTAXprefLR.gph FEEprefLR.gph ChildCareLR.gph EducationLR.gph ElderlyCareLR.gph HealthCareLR.gph CultureLR.gph TransportationLR.gph AdministrationLR.gph, rows(3) cols(3) graphregion(fcolor(white)) scheme(s2mono)
graph save Figure2.gph, replace
*graph play SurveyEdit
graph export ..\Figure2.png, replace
graph export ..\Figure2.eps, replace
*graph export ..\..\..\Parties\figuresFINAL\Figure2.eps, replace

cd ..\..\do\
