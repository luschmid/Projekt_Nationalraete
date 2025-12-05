use dta/illustration.dta, clear

label var mag "M"
keep if mag<16

scatter voteshareneeded magnitude, subtitle((i) Share of votes B needs to gain a seat)  ytitle("Vote share") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) msymbol(Oh) yscale(range(0 0.5)) ylabel(0(0.1)0.5) xscale(range(1(10)15)) xlabel(1(1)15) xsize(7) ysize(7)  
graph save figures/gph/vs.gph, replace
graph export figures/Figure1a.pdf, replace

scatter rawneeded magnitude, subtitle((ii) Raw votes B needs to gain a seat) ytitle("Votes") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) yscale(range(0 1000)) ylabel(0(250)1000) ymtick(0(250)1000)  msymbol(Oh) xscale(range(1(10)15)) xlabel(1(1)15) xsize(7) ysize(7)  
graph save figures/gph/raw.gph, replace
graph export figures/Figure1b.pdf, replace

gen z=0.1
gen adsneeded=voteshareneeded/z

gen y=0.5
gen doorsneeded=rawneed/y

scatter adsneeded magnitude, subtitle((iii) Ads B needs to run in order to gain a seat) ytitle("Ads") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) yscale(range(0 5) titlegap(*3)) ylabel(0(1)5) msymbol(Oh) xscale(range(1(10)15)) xlabel(1(1)15) xsize(7) ysize(7)  
graph save figures/gph/ads.gph, replace
graph export figures/Figure1c.pdf, replace

scatter doorsneeded magnitude, subtitle((iv) Contacts B needs to make to gain a seat) ytitle("Contacts") plotregion(lcolor(white) ilcolor(white)) graphregion(color(white)) scheme(s2mono) yscale(range(0 2000) titlegap(*5)) ylabel(0(500)2000) ymtick(0(500)2000)  msymbol(Oh) xscale(range(1(10)15)) xlabel(1(1)15) xsize(7) ysize(7)  
graph save figures/gph/doors.gph, replace
graph export figures/Figure1d.pdf, replace
