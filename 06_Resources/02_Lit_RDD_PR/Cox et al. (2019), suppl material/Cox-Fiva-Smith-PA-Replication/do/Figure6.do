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
twoway (scatter turnout margin if year<1920, msymbol(oh)) (lfit turnout margin if year<1920, lwidth(thick)), ytitle("Turnout") xtitle("Traditional measure") title(Pre-reform (SMDs)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/NOR_SMD_CFS.gph, replace
twoway (scatter turnout margin if year>1920, msymbol(oh)) (lfit turnout margin if year>1920, lwidth(thick)), ytitle("Turnout") xtitle("Traditional measure") title(Post-reform (MMDs)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/NOR_MMD_CFS.gph, replace

twoway (scatter turnout marginBL if year<1920, msymbol(oh)) (lfit turnout marginBL if year<1920, lwidth(thick)), ytitle("Turnout") xtitle("Blais-Lago measure") title(Pre-reform (SMDs)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/NOR_SMD_BL.gph, replace
twoway (scatter turnout marginBL if year>1920, msymbol(oh)) (lfit turnout marginBL if year>1920, lwidth(thick)), ytitle("Turnout") xtitle("Blais-Lago measure") title(Post-reform (MMDs)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/NOR_MMD_BL.gph, replace

twoway (scatter turnout C1 if year<1920, msymbol(oh)) (lfit turnout C1 if year<1920, lwidth(thick)), ytitle("Turnout") xtitle("Grofman-Selb measure") title(Pre-reform (SMDs)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/NOR_SMD_GS.gph, replace
twoway (scatter turnout C1 if year>1920, msymbol(oh)) (lfit turnout C1 if year>1920, lwidth(thick)), ytitle("Turnout") xtitle("Grofman-Selb measure") title(Post-reform (MMDs)) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/NOR_MMD_GS.gph, replace

cd figures/gph/
graph combine NOR_SMD_CFS.gph NOR_SMD_BL.gph NOR_SMD_GS.gph NOR_MMD_CFS.gph NOR_MMD_BL.gph NOR_MMD_GS.gph, ycommon xcommon cols(3) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) title(Norway)
graph save NOR_Turnout.gph, replace 
graph export ../Figure6a.pdf, replace 
cd ../../

*************************************************************************************
*************************************************************************************
*************************************************************************************
use dta/Panel.dta, clear
keep if country=="SWI"

twoway (scatter turnout margin if m==1, msymbol(oh)) (lfit turnout margin if m==1, lwidth(thick)), ytitle("Turnout") xtitle("Traditional measure") title(SMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/SWI_SMD_CFS.gph, replace
twoway (scatter turnout margin if m>1, msymbol(oh)) (lfit turnout margin if m>1, lwidth(thick)), ytitle("Turnout") xtitle("Traditional measure") title(MMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/SWI_MMD_CFS.gph, replace

twoway (scatter turnout marginBL if m==1, msymbol(oh)) (lfit turnout marginBL if m==1, lwidth(thick)), ytitle("Turnout") xtitle("Blais-Lago measure") title(SMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/SWI_SMD_BL.gph, replace
twoway (scatter turnout marginBL if m>1, msymbol(oh)) (lfit turnout marginBL if m>1, lwidth(thick)), ytitle("Turnout") xtitle("Blais-Lago measure") title(MMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/SWI_MMD_BL.gph, replace

twoway (scatter turnout C1 if m==1, msymbol(oh)) (lfit turnout C1 if m==1, lwidth(thick)), ytitle("Turnout") xtitle("Grofman-Selb measure") title(SMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/SWI_SMD_GS.gph, replace
twoway (scatter turnout C1 if m>1, msymbol(oh)) (lfit turnout C1 if m>1, lwidth(thick)), ytitle("Turnout") xtitle("Grofman-Selb measure") title(MMDs) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) xscale(range(0 0.9)) yscale(range(0 1.0)) legend(off)
graph save figures/gph/SWI_MMD_GS.gph, replace

cd figures/gph/
graph combine SWI_SMD_CFS.gph SWI_SMD_BL.gph SWI_SMD_GS.gph SWI_MMD_CFS.gph SWI_MMD_BL.gph SWI_MMD_GS.gph, ycommon xcommon cols(3) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) title(Switzerland)
graph save SWI_Turnout.gph, replace 
graph export ../Figure6b.pdf, replace 
cd ../../
