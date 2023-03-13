
program rdd_estimation

* (a) Summary statistics and bandwidth definition

sum $covs

rdrobust elected_F1 `1', p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1) 
global bw_opt_p1 =e(h_l)
global bbw_opt_p1 =e(b_l)
di $bw_opt_p1
di $bbw_opt_p1

rdrobust elected_F1 `1', p(2) bwselect(mserd) kernel($kernel)  ///
	vce($secl ID_pers) all rho(1)
global bw_opt_p2 =e(h_l)
global bbw_opt_p2 =e(b_l)
di $bw_opt_p2
di $bbw_opt_p2

* (b) RD plot

foreach z of varlist $outcomes {

local label: var label `z' 

cap rdplot `z' `1'  ///
	if -$bw_opt_p1 <= `1' & `1' <= $bw_opt_p1 & elected_F1 < ., ///
	p(1) kernel($kernel) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`2'_`1'_`z'_RD_Plot.pdf", ///
	as(pdf) replace
}

* (c) Bandwidth figures

foreach z of varlist $outcomes {
qui {
mat pointestimate = J(1,15,.)
mat CI = J(2,15,.)
local col=1
forv i = 0.1(0.1)1.6 {
local bw_opt_p1_tmp=$bw_opt_p1*`i'
local bbw_opt_p1_tmp=$bbw_opt_p1*`i'
rdrobust  `z' `1' if elected_F1 < ., p(1) h(`bw_opt_p1_tmp')  ///
	kernel($kernel) vce($secl ID_pers) all rho(1)
mat pointestimate[1,`col']=e(tau_cl)
mat CI[1,`col'] = e(ci_l_cl)
mat CI[2,`col'] = e(ci_r_cl)
local ++col
}
mat colnames pointestimate = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
mat colnames CI = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
coefplot matrix(pointestimate), ci(CI) vertical nolabel yline(0) ///
	mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`2'_`1'_`z'_Bw_Plot.pdf", ///
	as(pdf) replace
clear matrix

}
}


* (d) Balance table

local cp1_5 = 0
local cp1_10 = 0
local cp1_div2_5 = 0
local cp1_div2_10 = 0
local cp2_5 = 0
local cp2_10 = 0
local cp2_div2_5 = 0
local cp2_div2_10 = 0

foreach z of varlist $vars_all {
qui {
local label: var label `z' 
rdrobust  `z' `1' if elected_F1 < ., p(1) h($bw_opt_p1)  rho(1) ///
	kernel($kernel) vce($secl ID_pers) all
est store p1
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1" local ++cp1_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_10
local bw_opt_p1_div2=$bw_opt_p1/2
local bbw_opt_p1_div2=$bbw_opt_p1/2
rdrobust  `z' `1' if elected_F1 < ., p(1) h(`bw_opt_p1_div2')  rho(1) ///
	kernel($kernel) vce($secl ID_pers) all 
est store p1_div2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_div2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp1_div2_10
rdrobust  `z' `1' if elected_F1 < ., p(2) h($bw_opt_p2)  rho(1) ///
	kernel($kernel) vce($secl ID_pers) all
est store p2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_10
local bw_opt_p2_div2=$bw_opt_p2/2
local bbw_opt_p2_div2=$bbw_opt_p2/2
rdrobust  `z' `1' if elected_F1 < ., p(2) h(`bw_opt_p2_div2')  rho(1) ///
	kernel($kernel) vce($secl ID_pers) all
est store p2_div2
if `e(pv_rb)' < .05 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_div2_5
if `e(pv_rb)' < .1 & "`z'" != "elected_F1" & "`z'" != "participation_F1"  local ++cp2_div2_10

if "`z'" == "elected_F1" {
	esttab p1 p1_div2 p2 p2_div2 using "$path\05_Texts_and_presns\01_Running_variable\tables\\`2'_`1'_balance.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`label'") varwidth(21) star(* 0.10 ** 0.05 *** 0.01) fragment ///
	sca("h_l Bandwidth" "N_h_l Left obs" "N_h_r Right obs")  ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes replace

	}

else {
	esttab p1 p1_div2 p2 p2_div2 using "$path\05_Texts_and_presns\01_Running_variable\tables\\`2'_`1'_balance.tex", ///
	keep(Robust) b(%9.3f) se(%9.3f) coef(Robust "`label'") varwidth(21) fragment star(* 0.10 ** 0.05 *** 0.01) ///
	sfmt(%9.3f %9.0f %9.0f) noobs nogap noli nomti nonum nonotes append
	}
}
}

mat rrates = J(2,4,.)
mat rrates[1,1] =`cp1_5'
mat rrates[2,1] =`cp1_10'
mat rrates[1,2] =`cp1_div2_5'
mat rrates[2,2] =`cp1_div2_10'
mat rrates[1,3] =`cp2_5'
mat rrates[2,3] =`cp2_10'
mat rrates[1,4] =`cp2_div2_5'
mat rrates[2,4] =`cp2_div2_10'
mat rownames rrates = "5\%" "10\%"
mat colnames rrates = "\hspace{0.1cm}" "\hspace{0.2cm}" "\hspace{0.3cm}" "\hspace{0.4cm}"
*outtable using "$path\05_Texts_and_presns\01_Running_variable\tables\\`1'_balance", ///
*	mat(rrates) append nobox
esttab matrix(rrates) using "$path\05_Texts_and_presns\01_Running_variable\tables\\`2'_`1'_balance.tex",  ///
	nomti nonum append fragment


end