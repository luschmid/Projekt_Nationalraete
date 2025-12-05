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

keep if abs(dist_coal_margin)<.5  /* limiting sample to only 1122 obs. where we also have simulations for distance to threshold */

keep knr electionperiod dist_coal_margin SeatShareRIGHT LVoteShareRIGHT

binscatter SeatShareRIGHT dist_coal_margin

save ..\dta\___Export.dta, replace
