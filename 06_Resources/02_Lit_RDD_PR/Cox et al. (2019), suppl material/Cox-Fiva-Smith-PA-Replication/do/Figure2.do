use dta/ImportSurveyNor.dta, clear
binscatter visitANY magnitude, controls(d_year* d_edu* female age d_married*) title("Norway") ytitle("Visit by campaign worker") xtitle("M") msymbol(Oh) mcolor(black*0.5) lcolor(black) yscale(range(0 0.1)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xsize(7) ysize(7) 
graph play do/PlayThick.grec
graph export figures/Figure2a.pdf, replace

use dta/ImportSurveySwi.dta, clear
binscatter ci14_use S, title("Switzerland") controls(sex age d_year* d_educ* d_mari*) ytitle("Conversation with candidate") xtitle("M") msymbol(Oh) mcolor(black*0.5) lcolor(black) yscale(range(0 0.2)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) 	xlabel(1(4)5 10(5)35) xsize(7) ysize(7) 
graph play do/PlayThick.grec
graph export figures/Figure2b.pdf, replace
