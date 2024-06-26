clear
cap log close
set more 1
version 17

global path "E:\12. Cloud\Dropbox\Projekt Nationalräte\02_Processed_data\"


*** Data Analysis: Bisnode ***

** Steps
* Collapse on the NR candidate level and create different firm aggregates (e.g. all AG's, large AG's, ...)
* Create temporal structure of analysis (lags & leads)
* Run regressions

use "$path\11_Directors_1994_2018\RL_G7_NR-Bisnode_1994-2018_allMandates.dta", clear

** Collapse at the NR candidate level


gen party_sp= listname_bfs=="SP/PS" if year>=1971 
gen party_cvp= listname_bfs=="CVP/PDC" if year>=1971 
gen party_fdp= listname_bfs=="FDP/PLR (PRD)" if year>=1971 
gen party_svp= listname_bfs=="SVP/UDC" if year>=1971 
replace party_svp=1 if listname_bfs=="BDP/PBD" & year>=1971 
gen age = year - birthyear

local first canton name firstname birthyear age sex job votes candidate elected pvotes votemargin votemargin_rel ///
incumbent tenure list listname_bfs alliance suballiance eligible_cant voters_cant no_seats populationmun ///
gdename gdenr_2018_w ctn language_w eleyear VRid duns funktion mandid firstname_bis name_bis W_PLZ4 residence_bis ///
wohnland nationalität geburtstag 

gen sd_kapitalnominalag = kapitalnominalag
gen sd_kapitaleinbezahltag = kapitaleinbezahltag

collapse (first) `first' (sum) mandate VRAG_1plus VRAG_NC150k HighFct OpFct (mean) kapitalnominalag kapitaleinbezahltag (sd) sd_kapitalnominalag sd_kapitaleinbezahltag, by(NRid year)
rename kapitalnominalag mean_capnom
rename kapitaleinbezahltag mean_cappayed
rename sd_kapitalnominalag sd_capnom
rename sd_kapitaleinbezahltag sd_cappayed
replace mandate = . if year<1994
replace VRAG_1plus = . if year<1994
replace VRAG_NC150k = . if year<1994
replace HighFct = . if year<1994
replace OpFct = . if year<1994

/* Collapse legislative term

gen eleperiod = .
local i = 1
forv y = 1980(4)2018 {
	local z = `y' + 3
	replace eleperiod = `i' if inrange(year, `y', `z')
	local i = `i' + 1
}

collapse (first) `first' year (mean) mandate VRAG_1plus VRAG_NC150k HighFct OpFct mean_capnom mean_cappayed sd_capnom sd_cappayed, by(NRid eleperiod)
rename year legtermy1
label var legtermy1 "First year of legislative term"
*/

*save "$path\11_Directors_1994_2018\RL_G7_NR-Bisnode_1994-2018_allMandates_regs.dta", replace

** Create temporal structure
*use "$path\11_Directors_1994_2018\RL_G7_NR-Bisnode_1994-2018_allMandates_regs.dta", clear

** Analysis
encode canton, gen(ctnr)
encode NRid, gen(NRidNo)
xtset NRidNo year

/** simple OLS
foreach dep in mandate VRAG_1plus VRAG_NC150k  { // mean_capnom mean_cappayed
	xtreg `dep' L1.elected i.year, fe cluster(NRidNo)
	xtreg `dep' F4.elected F3.elected F2.elected F1.elected elected L1.elected i.year, fe cluster(NRidNo)
	*xtreg `dep' F4.elected F3.elected F2.elected F1.elected elected L1.elected L2.elected L3.elected L4.elected L5.elected L6.elected L7.elected L8.elected i.year, fe cluster(NRidNo)
}
*/
* Covariates
sum year
local min = r(min)
local max = r(max)
forv y = `min'(1)`max' {
	gen dy_`y' = 0
	replace dy_`y' = 1 if `y' == year
}
sum ctnr
local min = r(min)
local max = r(max)
forv i = `min'(1)`max' {
	gen ctn_`i' = 0
	replace ctn_`i' = 1 if `i' == ctnr
}

* Lags and leads
forv l = 1(1)8 {
	bysort NRidNo: gen votemargin_rel_l`l' = votemargin_rel[_n-`l']
	bysort NRidNo: gen elected_l`l' = elected[_n-`l']
}
forv f = 4(-1)1 {
	bysort NRidNo: gen votemargin_rel_f`f' = votemargin_rel[_n+`f']
	bysort NRidNo: gen elected_f`f' = elected[_n-`f']
}

* Year of Term
gen yterm = year-eleyear

* RDD
global covs sex age no_seats // sex age party_sp party_cvp party_fdp party_svp no_seats
	
global ydummy dy_1996 dy_1997 dy_1998 dy_1999 dy_2000 dy_2001 dy_2002 dy_2003 dy_2004 dy_2005 ///
	dy_2006 dy_2007 dy_2008 dy_2009 dy_2010 dy_2011 dy_2012 dy_2013 dy_2014 dy_2015 dy_2016 dy_2017 dy_2018 

global ctndummy ctn_1 ctn_2 ctn_3 ctn_4 ctn_5 ctn_6 ctn_7 ctn_8 ctn_9 ctn_10  ctn_12 ctn_13 ///
	ctn_14 ctn_15 ctn_16 ctn_17 ctn_18 ctn_19 ctn_20 ctn_21 ctn_22 ctn_23 ctn_24 ctn_25 ctn_26 /// ctn_11: drop JU

foreach dep in mandate VRAG_1plus VRAG_NC150k HighFct OpFct { // mean_capnom mean_cappayed
	forv i = 1(1)4 {
		cap drop `dep'_weight
		rdrobust `dep' votemargin_rel if yterm==`i' & year>1994, p(1) kernel(triangular) bwselect(mserd) covs($covs $ctndummy) all  // $ydummy
		local bw_opt_`dep' = e(h_l) 
		gen `dep'_weight = 1 - abs(votemargin_rel/`bw_opt_`dep'') if abs(votemargin_rel)<=`bw_opt_`dep'' & year>1994
		sum `dep' if abs(votemargin_rel)<=`bw_opt_`dep'' & year>1994
		reg `dep' i.elected##c.votemargin_rel $covs i.year i.ctnr [aweight=`dep'_weight] if abs(votemargin_rel)<=`bw_opt_`dep'' & year>1994, robust
	}
	/* Graphs with optimal BW
		rdplot `dep' votemargin_rel_l1 if abs(votemargin_rel_l1)<=`bw_opt_`dep'' & year>1994, ///
		c(0) p(1) kernel(triangular) h(`bw_opt_`dep'' `bw_opt_`dep'') nbins(25) binselect = "esmv" ///
		graph_options(legend(off) ///
		xtitle("RV", size(large)) ///
		ytitle("`dep'", size(large)) ///
		yscale(range(0 1)) ylabel(0(0.2)1, labsize(medlarge)) ///
		graphregion(color(white) lwidth(large))) */
}

* RDD for first election only





/*
* Dynamic: RDD
local ydummy dy_1996 dy_1997 dy_1998 dy_1999 dy_2000 dy_2001 dy_2002 dy_2003 dy_2004 dy_2005 ///
 dy_2006 dy_2007 dy_2008 dy_2009 dy_2010 dy_2011 dy_2012 dy_2013 dy_2014 dy_2015 dy_2016 dy_2017 dy_2018 
foreach dep in mandate VRAG_1plus VRAG_NC150k { // mean_capnom mean_cappayed
	rdrobust `dep' votemargin_rel if year>1994, p(1) kernel(triangular) bwselect(mserd) covs(`ydummy') all
	forv i=1(1)4 {
		rdrobust `dep' votemargin_rel_l`i' if year>1994, p(1) kernel(triangular) bwselect(mserd) covs(`ydummy') all
	}
}

* Dynamic: RDD with Blinder-Oaxaca decomposition
local ydummy dy_1996 dy_1997 dy_1998 dy_1999 dy_2000 dy_2001 dy_2002 dy_2003 dy_2004 dy_2005 ///
 dy_2006 dy_2007 dy_2008 dy_2009 dy_2010 dy_2011 dy_2012 dy_2013 dy_2014 dy_2015 dy_2016 dy_2017 dy_2018 
foreach dep in mandate VRAG_1plus VRAG_NC150k { // mean_capnom mean_cappayed
	cap drop `dep'_weight
	rdrobust `dep' votemargin_rel if year>1994, p(1) kernel(triangular) bwselect(mserd) covs(`ydummy') all 
	display e(N_h_l)+e(N_h_r)
	display (e(N_h_l)+e(N_h_r))/e(N)
	local bw_opt=e(h_l) 				 
	local bbw_opt=e(b_l)  
	local bw_opt_`dep' = `bw_opt'
	gen `dep'_weight=0
	replace `dep'_weight=1-abs(votemargin_rel_l1/`bw_opt_`dep'') if (votemargin_rel_l1 > -`bw_opt_`dep'' & votemargin_rel_l1 < `bw_opt_`dep'')
	reg `dep' i.elected##c.votemargin_rel i.year [iweight=`dep'_weight] if (votemargin_rel_l1 > -`bw_opt_`dep'' & votemargin_rel_l1 < `bw_opt_`dep'') & year>1994, robust
	est store est_`dep'_0
	forv i = 1(1)4 {
		predict resid_`dep'_`i'
		reg resid_`dep'_`i' i.elected_l`i'##c.votemargin_rel_l`i' i.year [iweight=`dep'_weight] if (votemargin_rel_l1 > -`bw_opt_`dep'' & votemargin_rel_l1 < `bw_opt_`dep'') & year>1994, robust
	est store est_`dep'_`i'
	}	
	esttab est_`dep'_0 est_`dep'_1 est_`dep'_2 est_`dep'_3 est_`dep'_4, b(%9.3f) se(%9.3f) ///
    star(* 0.10 ** 0.05 *** 0.01) sca(rkf widstat) ///
    title(Dependent variable: `dep') label nomtitles ///
    keep(1.elected 1.elected_l*) 
}
eststo clear
drop _est_* resid_* *_weight


erase "$path\11_Directors_1994_2018\RL_G7_NR-Bisnode_1994-2018_allMandates_regs.dta"




