
gen ss_coal= 100-SeatShareLEFT*100
gen vs_coal= 100-LVoteShareLEFT*100
gen ss_maj= 100-SeatShareLEFT*100
gen vs_maj= 100-LVoteShareLEFT*100
gen sim_dist=dist_coal_margin *100
gen MajRight= 1-MajLeft

vsss2dis MajRight ss_coal  2 "ss_coal!=. &ss_maj>40 & ss_maj<60" "MajRight==0" "MajRight==1" 1000 "Observations" "Right Wing Seat Share""Right Wing Seat Minority" "Right Wing Seat Majority" 
graph save ss_obs, replace

vsss2dis  MajRight vs_maj 2 "ss_coal!=. &vs_maj>40 & vs_maj<60" "MajRight==0" "MajRight==1" 1000 "Observations" "Right Wing Vote Share""Right Wing Seat Minority" "Right Wing Seat Majority" 
graph save vs_obs, replace


twoway (scatter ss_maj vs_maj  if ss_coal!=. &vs_maj>40 & vs_maj<60, msize(small)),  ytitle(Right Wing Seat Share) xtitle(Right Wing Vote Share)legend(off) xline(50, lwidth(medthick) lcolor(black))

graph save figa3, replace
graph export figa3.eps, replace

vsss2sim   MajRight sim_dist  2 "SeatShareLEFT!=. " "MajRight==0" "MajRight==1" 10 "Observations" "Distance to Right Wing Seat Majority""Right Wing Seat Minority" "Right Wing Seat Majority"

graph save figa5, replace
graph export figa5.eps, replace


