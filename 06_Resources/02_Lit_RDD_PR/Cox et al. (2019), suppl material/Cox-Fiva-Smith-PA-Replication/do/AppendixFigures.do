*****************************************************************************************
************************* Appendix Figure B1 ********************************************
*****************************************************************************************
use dta/ImportSurveyNor.dta, clear
binscatter visitANY magnitude, linetype(qfit) controls(d_year* d_edu* female age d_married*) title("Norway") ytitle("Visit by campaign worker") xtitle("M") msymbol(Oh) mcolor(black*0.5) lcolor(black) yscale(range(0 0.1)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xsize(7) ysize(7) 
graph play do/PlayThick.grec
graph export figures/FigureB1a.pdf, replace

use dta/ImportSurveySwi.dta, clear
binscatter ci14_use S, linetype(qfit) title("Switzerland") controls(sex age d_year* d_educ* d_mari*) ytitle("Conversation with candidate") xtitle("M") msymbol(Oh) mcolor(black*0.5) lcolor(black) yscale(range(0 0.2)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xlabel(1(4)5 10(5)35) xsize(7) ysize(7) 
graph play do/PlayThick.grec
graph export figures/FigureB1b.pdf, replace


*****************************************************************************************
************************* Appendix Figure B2 ********************************************
*****************************************************************************************
use dta/ImportSurveyEur.dta, clear
binscatter scalable M, title("") ytitle("Hours") xtitle("M") msymbol(Oh) mcolor(black*0.5) lcolor(black) yscale(range(0 16)) xlabel(1(19)20 20(20)120)  ylabel(0(4)16) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xsize(7) ysize(7) 
graph play do/PlayThick.grec
graph export figures/FigureB2.pdf, replace


*****************************************************************************************
************************* Appendix Figure B3 ********************************************
*****************************************************************************************
use dta/AlternativeMargin2, clear
sort year SMD_district
save dta/AlternativeMargin2_sort, replace

use dta/Turnout_District_Level, clear 
sort year SMD_district
drop _merge
merge year SMD_district using dta/AlternativeMargin2_sort

drop if _merge==2 /* PR-districts not in sample */
drop if year>1927
replace margin=min_distance1 if year>1920

gen marginBL=margin*magnitude  /* Blais-Lago normalize with district magnitude */
set obs 645
replace magnitude = 2 in 645 /* to add M=2 to x-axis */

*************************************************************************************
*************************************************************************************
*************************************************************************************
graph box C2, over(magnitude) title("") ytitle("Grofman-Selb alternative measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white))  scheme(s2mono) legend(off) yscale(range(0 1)) ylabel(0(0.2)1)
graph save figures/gph/NOR_GS_alt_M.gph, replace

cd figures/gph/
graph combine NOR_GS_alt_M.gph, plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) cols(1) ycommon title(Norway)
graph play ../../do/PlayM.grec
graph save NOR_Magnitude.gph, replace
cd ../../

*************************************************************************************
*************************************************************************************
*************************************************************************************
use dta/Panel.dta, clear
/* PUT ALL M>10 AS M=10 */
replace m=10 if m>10
keep if country=="SWI"

graph box C2, over(m) title("") ytitle("Grofman-Selb alternative measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) legend(off) yscale(range(0 1)) ylabel(0(0.2)1)
graph save figures/gph/SWI_GS_alt_M.gph, replace

cd figures/gph/
graph combine SWI_GS_alt_M.gph, plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) cols(3) ycommon title(Switzerland)
graph play ../../do/Play10.grec
graph play ../../do/PlayM.grec
graph save SWI_Magnitude.gph, replace

graph combine NOR_Magnitude.gph SWI_Magnitude.gph, rows(1) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) ycommon
graph export ../FigureB3.pdf, replace
cd ../../

*****************************************************************************************
************************* Appendix Figure B4 ********************************************
*****************************************************************************************

use dta/Panel.dta, clear
/* PUT ALL M>10 AS M=10 */
replace m=10 if m>10
keep if country=="SPA"

set obs 468
replace m = 2 in 468  /* to make m=2 show in figure */

graph box margin, over(m) title("") ytitle("Traditional measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) legend(off) ysize(7) xsize(7) yscale(range(0 1)) ylabel(0(0.2)1)
graph save figures/gph/SPA_CFS_M.gph, replace
graph play do/PlaySpaED.grec
graph export figures/FigureB4a.pdf, replace
graph box marginBL, over(m) title("") ytitle("Blais-Lago measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) legend(off) ysize(7) xsize(7) yscale(range(0 1)) ylabel(0(0.2)1)
graph save figures/gph/SPA_BL_M.gph, replace
graph play do/PlaySpaED.grec
graph export  figures/FigureB4b.pdf, replace
graph box C1, over(m) title("") ytitle("Grofman-Selb measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) legend(off) ysize(7) xsize(7) yscale(range(0 1)) ylabel(0(0.2)1)
graph save figures/gph/SPA_GS_M.gph, replace
graph play do/PlaySpaED.grec
graph export figures/FigureB4c.pdf, replace


*****************************************************************************************
************************* Appendix Figure B5 ********************************************
*****************************************************************************************
use dta/Panel.dta, clear

keep if country=="SPA"

twoway (scatter turnout margin if m==1, msymbol(oh)) (lfit turnout margin if m==1, lwidth(thick)), ytitle("Turnout") xtitle("Traditional measure") title(SMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1)) ylabel(0(0.2)1) legend(off)
graph save figures/gph/SPA_SMD_CFS.gph, replace
twoway (scatter turnout margin if m>1, msymbol(oh)) (lfit turnout margin if m>1, lwidth(thick)), ytitle("Turnout") xtitle("Traditional measure") title(MMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1)) ylabel(0(0.2)1) legend(off)
graph save figures/gph/SPA_MMD_CFS.gph, replace

twoway (scatter turnout marginBL if m==1, msymbol(oh)) (lfit turnout marginBL if m==1, lwidth(thick)), ytitle("Turnout") xtitle("Blais-Lago measure") title(SMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1)) ylabel(0(0.2)1) legend(off)
graph save figures/gph/SPA_SMD_BL.gph, replace
twoway (scatter turnout marginBL if m>1, msymbol(oh)) (lfit turnout marginBL if m>1, lwidth(thick)), ytitle("Turnout") xtitle("Blais-Lago measure") title(MMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1)) ylabel(0(0.2)1) legend(off)
graph save figures/gph/SPA_MMD_BL.gph, replace

twoway (scatter turnout C1 if m==1, msymbol(oh)) (lfit turnout C1 if m==1, lwidth(thick)), ytitle("Turnout") xtitle("Grofman-Selb measure") title(SMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1)) ylabel(0(0.2)1) legend(off)
graph save figures/gph/SPA_SMD_GS.gph, replace
twoway (scatter turnout C1 if m>1, msymbol(oh)) (lfit turnout C1 if m>1, lwidth(thick)), ytitle("Turnout") xtitle("Grofman-Selb measure") title(MMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1)) ylabel(0(0.2)1) legend(off)
graph save figures/gph/SPA_MMD_GS.gph, replace

cd figures/gph/
graph combine SPA_SMD_CFS.gph SPA_SMD_BL.gph SPA_SMD_GS.gph SPA_MMD_CFS.gph SPA_MMD_BL.gph SPA_MMD_GS.gph, ycommon xcommon cols(3) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) title("")
graph export ../FigureB5.pdf, replace 
cd ../../

*****************************************************************************************
************************* Appendix Figure B6 ********************************************
*****************************************************************************************
use dta/AlternativeMargin2, clear
sort year SMD_district
save dta/AlternativeMargin2_sort, replace

use dta/Turnout_District_Level, clear 
sort year SMD_district
drop _merge
merge year SMD_district using dta/AlternativeMargin2_sort

drop if _merge==2 /* PR-districts not in sample */
drop if year>1927
replace margin=min_distance1 if year>1920

gen marginBL=margin*magnitude  /* Blais-Lago normalize with district magnitude */

*************************************************************************************
*************************************************************************************
*************************************************************************************
twoway (scatter turnout C2 if year<1920, msymbol(oh)) (lfit turnout C2 if year<1920, lwidth(thick)), ytitle("Turnout") xtitle("Grofman-Selb alternative measure") title(Pre-reform (SMDs)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1)) ylabel(0(0.2)1) legend(off)
graph save figures/gph/NOR_SMD_GS_alt.gph, replace
twoway (scatter turnout C2 if year>1920, msymbol(oh)) (lfit turnout C2 if year>1920, lwidth(thick)), ytitle("Turnout") xtitle("Grofman-Selb alternative measure") title(Post-reform (MMDs)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1)) ylabel(0(0.2)1) legend(off)
graph save figures/gph/NOR_MMD_GS_alt.gph, replace

cd figures/gph/
graph combine NOR_SMD_GS_alt.gph NOR_MMD_GS_alt.gph, ycommon xcommon cols(2) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) title(Norway)
graph save NOR_Turnout.gph, replace 
cd ../../

*************************************************************************************
*************************************************************************************
*************************************************************************************
use dta/Panel.dta, clear
keep if country=="SWI"

twoway (scatter turnout C2 if m==1, msymbol(oh)) (lfit turnout C2 if m==1, lwidth(thick)), ytitle("Turnout") xtitle("Grofman-Selb alternative measure") title(SMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1)) ylabel(0(0.2)1) legend(off)
graph save figures/gph/SWI_SMD_GS_alt.gph, replace
twoway (scatter turnout C2 if m>1, msymbol(oh)) (lfit turnout C2 if m>1, lwidth(thick)), ytitle("Turnout") xtitle("Grofman-Selb alternative measure") title(MMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1)) ylabel(0(0.2)1) legend(off)
graph save figures/gph/SWI_MMD_GS_alt.gph, replace

cd figures/gph/
graph combine SWI_SMD_GS_alt.gph SWI_MMD_GS_alt.gph, ycommon xcommon cols(2) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) title(Switzerland)
graph save SWI_Turnout.gph, replace 
cd ../../

cd figures/gph/
graph combine NOR_Turnout.gph SWI_Turnout.gph, ycommon xcommon cols(1) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) 
graph export ../FigureB6.pdf, replace 
cd ../../
