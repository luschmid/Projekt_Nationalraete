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
graph box margin, over(magnitude) title("") ytitle("Traditional measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white))  scheme(s2mono) legend(off)
graph save figures/gph/NOR_CFS_M.gph, replace
graph box marginBL, over(magnitude) title("") ytitle("Blais-Lago measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white))  scheme(s2mono) legend(off)
graph save figures/gph/NOR_BL_M.gph, replace
graph box C1, over(magnitude) title("") ytitle("Grofman-Selb measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white))  scheme(s2mono) legend(off)
graph save figures/gph/NOR_GS_M.gph, replace

cd figures/gph/
graph combine NOR_CFS_M.gph NOR_BL_M.gph NOR_GS_M.gph, plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) cols(3) ycommon title(Norway)
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

graph box margin, over(m) title("") ytitle("Traditional measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) legend(off)
graph save figures/gph/SWI_CFS_M.gph, replace
graph box marginBL, over(m) title("") ytitle("Blais-Lago measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) legend(off)
graph save figures/gph/SWI_BL_M.gph, replace
graph box C1, over(m) title("") ytitle("Grofman-Selb measure") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) legend(off)
graph save figures/gph/SWI_GS_M.gph, replace

cd figures/gph/
graph combine SWI_CFS_M.gph SWI_BL_M.gph SWI_GS_M.gph, plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) cols(3) ycommon title(Switzerland)
graph play ../../do/Play10.grec
graph play ../../do/PlayM.grec
graph save SWI_Magnitude.gph, replace

graph combine NOR_Magnitude.gph SWI_Magnitude.gph, rows(2) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono)
graph export ../Figure5.pdf, replace
cd ../../
