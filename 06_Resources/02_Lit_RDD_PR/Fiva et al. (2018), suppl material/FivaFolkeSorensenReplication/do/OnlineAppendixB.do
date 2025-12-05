use ..\dta\MergedData, clear


****************************
****************************

sort knr ele
merge knr ele using "../dta/val med simuleringar DH4.dta"
drop _merge
bysort ele: sum dist_coal_margin
sort knr ele
merge knr ele using "../dta/val med simuleringar4.dta", update
bysort ele: sum dist_coal_margin

gen dist_coal_marginxMajLeft=dist_coal_margin*MajLeft
keep if abs(dist_coal_margin)<.5  /* limiting sample to only 1122 obs. where we also have simulations for distance to threshold */

gen SeatShareRIGHT=1-SeatShareLEFT
gen LVoteShareRIGHT=1-LVoteShareLEFT
gen MajRight=1-MajLeft

gen ss_coal= SeatShareRIGHT*100
gen vs_coal= LVoteShareRIGHT*100
gen ss_maj= SeatShareRIGHT*100
gen vs_maj= LVoteShareRIGHT*100
gen sim_dist=dist_coal_margin *100

****************************
****************************
append using "..\dta\online zeroes.dta"
vsss2dis MajRight ss_coal  2 "ss_coal!=. &ss_coal>40 & ss_coal<60" "MajRight!=0 & ss_coal>50 " "MajRight!=1  & ss_coal<50" 1000 "Observations" "Right-Wing Seat Share""Right Wing Majority" "Right Wing Minority" 
*graph play WhiteBackground
graph save ../figures/FIGURE_B1.gph, replace
graph export ../figures/FIGURE_B1.eps, replace
graph export ../figures/FIGURE_B1.tif, replace
*graph export ../../Parties/figuresFINAL/FIGURE_B1.eps, replace

vsss2dis ss_maj vs_maj 2 "vs_maj>40 & vs_maj<60" "MajRight!=0" "MajRight!=1" 1000 "Observations" "Right Wing Vote Share""Right Wing Majority" "Right Wing Minority" 
*graph play WhiteBackground
graph save ../figures/FIGURE_B2.gph, replace
graph export ../figures/FIGURE_B2.eps, replace
graph export ../figures/FIGURE_B2.tif, replace
*graph export ../../Parties/figuresFINAL/FIGURE_B2.eps, replace

twoway (scatter ss_maj vs_maj  if ss_coal!=. &vs_maj>40 & vs_maj<60, msize(small)),  ytitle(Right Wing Seat Share) xtitle(Right Wing Vote Share)legend(off) xline(50, lwidth(medthick) lcolor(gray))
*graph play WhiteBackground
graph save ../figures/FIGURE_B3.gph, replace
graph export ../figures/FIGURE_B3.eps, replace
graph export ../figures/FIGURE_B3.tif, replace
*graph export ../../Parties/figuresFINAL/FIGURE_B3.eps, replace

vsss2sim LVoteShareLEFT sim_dist  2 "SeatShareLEFT!=. " "MajRight==1" "MajRight==0" 10 "Observations" "Distance to Majority Change from Left to Right""Right Wing Majority" "Left Wing Majority"
*graph play WhiteBackground
graph save ../figures/FIGURE_B5.gph, replace
graph export ../figures/FIGURE_B5.eps, replace
graph export ../figures/FIGURE_B5.tif, replace
*graph export ../../Parties/figuresFINAL/FIGURE_B5.eps, replace

*** NOTE: Figure B4 is not included in the replication file
