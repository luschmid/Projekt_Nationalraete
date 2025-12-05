use dta/Turnout_District_Level, clear 

sort year PR_district
drop _merge
merge year PR_district using dta/AlternativeMargin

drop if _merge==2 /* PR-districts not in sample */

keep if year>1920 /* ONLY MMDs */
***********************************************
***********************************************
***********************************************
**** collapse to PR level
collapse min_distance min_distance1 min_distance_iXv_i min_distance1_iXv_i, by(year PR_district)

***********************************************
***********************************************
***********************************************
**** weight or use minima?

*(I) Multi-party (minima-based) on single-party (minima-based)
twoway (scatter min_distance min_distance1, msymbol(Oh)) (function y=x, range(0 0.15)), xsize(7) ysize(7) ///
yscale(range(0 0.14)) ylabel(0(0.02)0.14) xlabel(0(0.02)0.14) xscale(range(0 0.14)) legend(off) ///
ytitle(Multi-party measure (minima)) xtitle(Single-party measure (minima)) ///
plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) 
graph play do/PlayThick.grec
graph save figures/gph/I.gph, replace

*(II) Multi-party (weighted average) on single-party (minima-based)

twoway (scatter min_distance_iXv_i min_distance1, msymbol(Oh)) (function y=x, range(0 0.15)), xsize(7) ysize(7) ///
yscale(range(0 0.14)) ylabel(0(0.02)0.14) xlabel(0(0.02)0.14) xscale(range(0 0.14)) legend(off) ///
ytitle(Multi-party measure (weighted)) xtitle(Single-party measure (minima)) ///
plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) 
graph play do/PlayThick.grec
graph save figures/gph/II.gph, replace

*(III) Single-party (weighted average) on single-party (minima-based)
twoway (scatter min_distance1_iXv_i min_distance1, msymbol(Oh)) (function y=x, range(0 0.15)), xsize(7) ysize(7) ///
yscale(range(0 0.14)) ylabel(0(0.02)0.14) xlabel(0(0.02)0.14) xscale(range(0 0.14)) legend(off) ///
ytitle(Single-party measure (weighted)) xtitle(Single-party measure (minima)) ///
plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) 
graph play do/PlayThick.grec
graph save figures/gph/III.gph, replace

cd figures/gph/
graph combine I.gph II.gph III.gph, ycommon xcommon plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) cols(3) title(Norway)
graph play ../../do/PlayAspectSize.grec
graph save Figure4_NOR.gph, replace
cd ../../


use dta/Switzerland_CFS.dta, clear
keep if m>1 /* only MMDs */

**** collapse to PR level
collapse min_distance min_distance1 min_distance_iXv_i min_distance1_iXv_i, by(year id)

***********************************************
***********************************************
***********************************************
**** weight or use minima?

*(I) Multi-party (minima-based) on single-party (minima-based)
twoway (scatter min_distance min_distance1, msymbol(Oh)) (function y=x, range(0 0.24)), xsize(7) ysize(7) ///
yscale(range(0 0.24)) ylabel(0(0.04)0.24) xlabel(0(0.04)0.240) xscale(range(0 0.24)) legend(off) ///
ytitle(Multi-party measure (minima)) xtitle(Single-party measure (minima)) ///
plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) 
graph play do/PlayThick.grec
graph save figures/gph/I.gph, replace

*(II) Multi-party (weighted average) on single-party (minima-based)

twoway (scatter min_distance_iXv_i min_distance1, msymbol(Oh)) (function y=x, range(0 0.24)), xsize(7) ysize(7) ///
yscale(range(0 0.24)) ylabel(0(0.04)0.24) xlabel(0(0.04)0.240) xscale(range(0 0.24)) legend(off) ///
ytitle(Multi-party measure (weighted)) xtitle(Single-party measure (minima)) ///
plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) 
graph play do/PlayThick.grec
graph save figures/gph/II.gph, replace

*(III) Single-party (weighted average) on single-party (minima-based)
twoway (scatter min_distance1_iXv_i min_distance1, msymbol(Oh)) (function y=x, range(0 0.24)), xsize(7) ysize(7) ///
yscale(range(0 0.24)) ylabel(0(0.04)0.24) xlabel(0(0.04)0.240) xscale(range(0 0.24)) legend(off) ///
ytitle(Single-party measure (weighted)) xtitle(Single-party measure (minima)) ///
plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) 
graph play do/PlayThick.grec
graph save figures/gph/III.gph, replace

cd figures/gph/
graph combine I.gph II.gph III.gph, ycommon xcommon plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) cols(3) title(Switzerland)
graph play ../../do/PlayAspectSize.grec
graph save Figure4_SWI.gph, replace
cd ../../

cd figures/gph/
graph combine Figure4_NOR.gph Figure4_SWI.gph, rows(2) plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono)
graph export ../Figure4.pdf, replace
cd ../../
