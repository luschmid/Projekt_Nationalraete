******************************************************
* (A) Set directories
******************************************************

clear
set more off

global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
*global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

******************************************************
* (B) Read in data on margins
******************************************************

import delimited using ///
	"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981.csv", ///
	encoding("utf-8") clear 
save "$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981.dta", ///
	replace

********************************************************************
* (C) Read in data from Fiva and Smith (2018) and generate margins
********************************************************************

use "$path\01_Raw_data\13_Running_variable\FivaSmithJune2019.dta", clear

* (i) Correction

replace pid=22893 if candidatename_orig=="Lars Fagerland" & year==1973
* Note: The original pid refers to the pid of another person. 
*		Fiva and Smith (2018) made a mistake by editing the original candidate 
*		name. It changed from Lars Fagerland to Styrk Lothe. We have to make 
*		this change, otherwise the candidates run in the same district and year 
*		in different parties

* (ii) Focus on years 1953-1985

keep if inrange(year,1953,1985)
*drop if candidatename_orig=="Knudsen, A. G." & year==1912
* note: drop duplicate

* (iii) Merge vote margins for Norway and generate running variable

merge 1:1 year partyname districtid pid using ///
	"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981.dta", ///
	gen(merge1)
drop if merge1 == 2
* Drop one votemargin w/out candidate
drop merge1

gen votemargin_vote=votemargin/approvedvotesoverall
gen votemargin_elig=votemargin/electorate
rename margin votemargin_fiva

***************************************
* (D) Destats and comparison of margins
***************************************

* (i) Get number of observations

tab pid if votemargin==. & votemargin_fiva!=. &  year>=1953 & year<=1981
* result: no observation with Fiva-Smith margin but not with our margin
tab pid if votemargin!=. & votemargin_fiva==. &  year>=1953 & year<=1981
* result: 8,274 observations with our margin but with no Fiva-Smith margin
tab pid if votemargin!=. & votemargin_fiva!=. &  year>=1953 & year<=1981
* result: 1,871 observations with both margins

* (ii) Inspect differences 

cor votemargin_vote votemargin_fiva 

*br margin votemargin votemargin_vote if year>=1953 & year<=1981

gen votemargin_vote_diff=votemargin_fiva-votemargin_vote 

hist votemargin_vote_diff

* result: it appears that the rescalation with voters is more similar to ours
* 		  in absolute terms

sum votemargin_vote_diff, d
tab votemargin_vote_diff

* (iii) Calculate absolut vote margin for Fiva-Smith variables

gen margin_absolute=votemargin_fiva*approvedvotesoverall
* result: both numbers are not integers

* (iv) where are the largest differences

scatter votemargin_vote_diff rank
tab rank if votemargin_vote_diff < .

******************
* (E) Application 
******************

* (i) Gen unique pid
* Note: We do not know how and when Jon drops duplicate pids. Instead, 
* 		we use our way and assign a unique pid per district

gen double pid_distid=districtid*100000+pid
egen time_id=group(year)
xtset pid_distid time_id

sort pid_distid time_id
gen elected_F1=F1.elected
replace elected_F1=0 if elected_F1==.

save "$path\02_Processed_data\13_Running_variable\FivaSmithJune2019_processed.dta", replace
* Note: For R file 16_Norway_Histogram

/* (ii) Sample restrictions by Fiva and Smith (2018)
* Note: See DataPrep.do line 18 to 35.

* (a) Focus on districts in which candidates achieved highest rank 

bysort pid year: egen min_rank=min(rank)
keep if rank==min_rank  /* keeping only the highest ranked entry for each candidate */

drop if candidatename_ed=="Anders Lange" & district=="vest-agder" /* Anders Lange was ranked first in both Oslo and Vest-Agder, but performed best in Oslo */

* (b) Focus on main parties

keep if party=="nkp"|party=="sv"|party=="dna"|party=="sp"|party=="v"|party=="krf"|party=="h"|party=="frp"
*/

* (iii) Define margin to use, kernel, covariates, clsuter

keep if inrange(year,1953,1981)
global margin votemargin_elig
global kernel triangular // Options: triangular (default), epanechnikov, uniform
global secl cluster // Options: cluster nncluster
global p 1

* (iv) Estimations

rdrobust elected_F1 votemargin_fiva, p($p) bwselect(mserd) kernel($kernel)  ///
	vce($secl pid_distid) all

rdrobust elected_F1 $margin, p($p) bwselect(mserd) kernel($kernel)  ///
	vce($secl pid_distid) all, if votemargin_fiva!=.
global bw_sample_fiva =e(h_l)
global bbw_sample_fiva =e(b_l)
	
rdrobust elected_F1 $margin, p($p) bwselect(mserd) kernel($kernel)  ///
	vce($secl pid_distid) all
global bw_full_sample =e(h_l)
global bbw_full_sample =e(b_l)

rdplot elected_F1 $margin ///
	if -$bw_sample_fiva <= $margin & $margin<= $bw_sample_fiva ///
	& votemargin_fiva < . & elected_F1 < ., ///
	p($p) kernel($kernel) h($bw_sample_fiva $bw_sample_fiva) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\norway_fiva_sample.pdf", ///
	as(pdf) replace	
	
rdplot elected_F1 $margin ///
	if -$bw_full_sample <= $margin & $margin<= $bw_full_sample & elected_F1 < ., ///
	p($p) kernel($kernel) h($bw_full_sample $bw_full_sample) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\norway_full_sample.pdf", ///
	as(pdf) replace
	
* (v) Bandwidth figures


*qui {
mat pointestimate = J(1,15,.)
mat CI = J(2,15,.)
local col=1
forv i = 0.1(0.1)1.6 {
local bw_opt_p1_tmp=$bw_full_sample*`i'
local bbw_opt_p1_tmp=$bbw_full_sample*`i'
rdrobust  elected_F1 $margin if elected_F1 < ., p(1) h(`bw_opt_p1_tmp') b(`bbw_opt_p1_tmp') ///
	kernel($kernel) vce($secl pid_distid) all
mat pointestimate[1,`col']=e(tau_bc)
mat CI[1,`col'] = e(ci_l_rb)
mat CI[2,`col'] = e(ci_r_rb)
local ++col
}
mat colnames pointestimate = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
mat colnames CI = ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" 
coefplot matrix(pointestimate), ci(CI) vertical nolabel yline(0) ///
	mcolor(black) ciopts(lc(black)) ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\norway_bw.pdf", ///
	as(pdf) replace
clear matrix
*}

