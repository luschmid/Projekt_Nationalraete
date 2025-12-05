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

** OLD: vsss2dis MajRight ss_coal  2 "ss_coal!=. &ss_maj>40 & ss_maj<60" "MajRight==1" "MajRight==0" 1000 "Observations" "Coalition Seat Share""Seat Minority" "Seat Majority" 
vsss2dis MajRight ss_coal  2 "ss_coal!=. &ss_maj>40 & ss_maj<60" "MajRight==1" "MajRight==0" 1000 "Observations" "Right-Wing Seat Share""Right Wing Majority" "Right Wing Minority" 
graph save ../figures/ss_obs.gph, replace
graph export ../figures/ss_obs.eps, replace
graph export ../figures/ss_obs.tif, replace

** OLD: vsss2dis ss_coal vs_maj 2 "ss_coal!=. &vs_maj>40 & vs_maj<60" "MajRight==1" "MajRight==0" 1000 "Observations" "Coalition Vote Share""Seat Minority" "Seat Majority" 
vsss2dis ss_coal vs_maj 2 "ss_coal!=. &vs_maj>40 & vs_maj<60" "MajRight==1" "MajRight==0" 1000 "Observations" "Right Wing Vote Share""Right Wing Majority" "Right Wing Minority" 
graph save ../figures/vs_obs.gph, replace
graph export ../figures/vs_obs.eps, replace
graph export ../figures/vs_obs.tif, replace

** JHF: vsssalt doesn't work for me.. Update graph programs?
*vsssalt MajRight vs_maj 2 "ss_coal!=. &vs_maj>45 & vs_maj<55"  1000 "Seat Majority" "Coalition Vote Share"
*graph save ../figures/vs_maj.gph, replace
*graph export ../figures/vs_maj.eps, replace
*graph export ../figures/vs_maj.tif, replace

** OLD: twoway (scatter ss_maj vs_maj  if ss_coal!=. &vs_maj>40 & vs_maj<60, msize(small)),  ytitle(Left Wing Seat Share) xtitle(Left Wing Vote Share)legend(off) xline(50, lwidth(medthick) lcolor(gray) lpattern(dash))
twoway (scatter ss_maj vs_maj  if ss_coal!=. &vs_maj>40 & vs_maj<60, msize(small)),  ytitle(Right Wing Seat Share) xtitle(Right Wing Vote Share)legend(off) xline(50, lwidth(medthick) lcolor(gray))
graph save figa3.gph, replace
graph export figa3.tif, replace
graph export figa3.eps, replace

** OLD vsss2sim LVoteShareLEFT sim_dist  2 "SeatShareLEFT!=. " "MajRight==1" "MajRight==0" 10 "Observations" "Distance to Left Wing Seat Majority""Left Wing Minority" "Left Wing Majority"
vsss2sim LVoteShareLEFT sim_dist  2 "SeatShareLEFT!=. " "MajRight==1" "MajRight==0" 10 "Observations" "Distance to Left Wing Seat Majority""Right Wing Majority" "Left Wing Majority"
graph save figa5.gph, replace
graph export figa5.tif, replace
graph export figa5.eps, replace
