use ..\dta\MergedData, clear

sort knr ele
merge knr ele using "../dta/val med simuleringar DH4.dta"
drop _merge
bysort ele: sum dist_coal_margin
sort knr ele
merge knr ele using "../dta/val med simuleringar4.dta", update
bysort ele: sum dist_coal_margin

keep if abs(dist_coal_margin)<.5  /* limiting sample to only 1122 obs. where we also have simulations for distance to threshold */

renpfix seatshare_ SeatShare

/* Table A.1: Descriptive statistics: Seat shares in the local council */

sutex SeatShare*, minmax title(Descriptive Statistics: Seatshares \label{DescrStatsParties}) ///
file(../tables/DescrStatsParties.tex) replace
