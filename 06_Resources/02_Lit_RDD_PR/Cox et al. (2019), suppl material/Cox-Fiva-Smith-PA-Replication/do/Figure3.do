use dta/ImportSurveyEur.dta, clear
binscatter canvassing M, title("") ytitle("Hours canvassing") xtitle("M") msymbol(Oh) mcolor(black*0.5) lcolor(black) xlabel(1(19)20 20(20)120) yscale(range(0 0.2)) ylabel(0(2)10) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xsize(7) ysize(7) 
graph play do/PlayThick.grec
graph export figures/Figure3.pdf, replace
