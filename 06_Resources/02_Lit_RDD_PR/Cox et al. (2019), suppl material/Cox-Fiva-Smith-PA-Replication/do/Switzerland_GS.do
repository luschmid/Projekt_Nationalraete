use dta/CH71-03parties.dta, clear

rename S m
rename kt id
drop s ss v_j s_j nop pv

/* INDEX OF COMPETITION, Swiss National Council elections, 1971-2007 */

/* NOTATION: id district, m district magnitude, v_i party i's vote share, s_i party i's no of seats */

/* Maximum vote share required for party i to gain an additional seat */

gen T_E = 1/(m+1)    /* T_E Threshold of exclusion */

gen x_G =(s_i+1)/(m+1)-v_i
replace x_G =. if s_i==m|x_G>T_E


/* Identify runner up j's vote share v_j and number of seats s_j */

gen s=0    /* counts seats allocated */
gen ss=0    /* counts sum of seats allocated */
gen pv=v_i    /* electoral quotient */
gen v_j=.
gen s_j=.
bys id year: egen nop=count(pv)    /*nop number of parties */
set more off
forval i=1/35 {        /* maximum number of rounds to allocate all the seats is 35 (=m) */
bys id year: egen r=rank(pv), u    /* this ranks the electoral quotients */
replace s=s+1 if r==nop&ss<m
gsort id year -r
replace v_j=v_i[_n+1] if r==nop&s==s_i&s_i>0
replace s_j=s_i[_n+1] if r==nop&s==s_i&s_i>0
replace pv=v_i/(s+1) if r==nop&ss<m
replace ss=ss+1 if ss<m
drop r
}

drop s-pv


/* Minimum vote share that could cost party i a seat */

gen x_L = (-s_i*v_j+s_j*v_i+v_i)/(s_i+s_j+1)


/* Indices of competition (C1, C2) */

gen c1_i = max(T_E-x_G, T_E-x_L) / T_E
sum c1_i

gen c2_i = 1-2*min(x_G, x_L)
sum c2_i

************ Aggregating to the district level ***************
gen C1=c1_i*v_i
gen C2=c2_i*v_i

collapse (mean) turnout (mean) m (sum) C1 (sum) C2, by(id year)
save dta/Switzerland_GS.dta, replace
