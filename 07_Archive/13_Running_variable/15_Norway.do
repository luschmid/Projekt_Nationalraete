******************************************************
* (A) Set directories
******************************************************

clear
set more off

*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte"
global path "C:\Schmidlu\Dropbox\Projekt Nationalräte"
*global path "E:/12. Cloud/Dropbox/Projekt Nationalräte"

******************************************************
* (B) Read in data on margins
******************************************************

import delimited using ///
	"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981.csv", ///
	encoding("utf-8") clear 
save "$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981.dta", ///
	replace
	
import delimited using ///
"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981_redis.csv", ///
encoding("utf-8") clear 
save ///
"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981_redis.dta", ///
	replace

import delimited using ///
	"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981_Easiest.csv", ///
	encoding("utf-8") clear 
keep votemargin year partyname districtid pid
rename votemargin votemargin_easiest
save "$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981_Easiest.dta", ///
	replace
	
********************************************************************
* (C) Read in data from Fiva and Smith (2018) and generate margins
********************************************************************

use "$path\01_Raw_data\13_Running_variable\FivaSmithJune2019.dta", clear

keep if inrange(year,1953,1985)
*drop if candidatename_orig=="Knudsen, A. G." & year==1912
* note: drop duplicate

merge 1:1 year partyname districtid pid using ///
	"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981.dta", ///
	gen(merge1)
merge 1:1 year partyname districtid rank using ///
"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981_redis.dta", ///
	gen(merge2)
merge 1:1 year partyname districtid pid using ///
"$path\02_Processed_data\13_Running_variable\RV_Norway_1953_1981_Easiest.dta", ///
	gen(merge3)
	
gen votemargin_vote=votemargin/approvedvotesoverall
gen votemargin_redis_vote=votemargin_redis/approvedvotesoverall
gen votemargin_vote_easiest=votemargin_easiest/approvedvotesoverall

tab pid if votemargin==. & margin!=. &  year>=1953 & year<=1981
* result: no observation with Fiva-Smith margin but not with our margin
tab pid if votemargin!=. & margin==. &  year>=1953 & year<=1981
* result: 8,274 observations with our margin but with no Fiva-Smith margin
tab pid if votemargin!=. & margin!=. &  year>=1953 & year<=1981
* result: 1,871 observations with both margins


******************************************************
* (D) Compare margins
******************************************************

* (i) Inspect differences 

cor votemargin_vote votemargin_vote_easiest votemargin_redis_vote margin 

*br margin votemargin votemargin_vote if year>=1953 & year<=1981

gen votemargin_vote_diff=margin-votemargin_vote 
gen votemargin_redis_vote_diff=margin-votemargin_redis_vote
gen votemargin_vote_easiest_diff=margin-votemargin_vote_easiest

hist votemargin_vote_diff
hist votemargin_redis_vote_diff

* result: it appears that the rescalation with voters is more similar to ours
* 		  in absolute terms

sum votemargin_vote_diff, d
tab votemargin_vote_diff

sum votemargin_redis_vote_diff, d
tab votemargin_redis_vote_diff

sum votemargin_vote_easiest_diff, d
tab votemargin_redis_vote_diff


* (ii) Calculate absolut vote margin for Fiva-Smith variables

gen margin_absolute=margin*approvedvotesoverall
* result: both numbers are not integers

* (iii) where are the largest differences?

scatter votemargin_vote_diff rank
tab rank if votemargin_vote_diff < .

scatter votemargin_redis_vote_diff rank
tab rank if votemargin_redis_vote_diff < .

*******************************************
* (E) Replication for Fiva and Smith (2018)
*******************************************

* (i) Correction

replace pid=22893 if candidatename_orig=="Lars Fagerland" & year==1973
* Note: The original pid refers to the pid of another person. 
*		Fiva and Smith (2018) made a mistake by editing the original candidate 
*		name. It changed from Lars Fagerland to Styrk Lothe. We have to make 
*		this change, otherwise the candidates run in the same district and year 
*		in different parties

* (ii) Gen unique pid
* Note: We do not know how and when Jon drops duplicate pids. Instead, 
* 		we use our way and assign a unique pid per district

gen double pid_unique=districtid*100000+pid
sort pid year
egen time_id=group(year)
xtset pid_unique time_id

sort pid_unique time_id
gen elected_F1=F1.elected
replace elected_F1=0 if elected_F1==.

save "$path\02_Processed_data\13_Running_variable\FivaSmithJune2019_processed.dta", replace

* (iii) Sample restriction by Fiva and Smith (2018)
* Note: See DataPrep.do line 18 to 35.

bysort pid year: egen min_rank=min(rank)
keep if rank==min_rank  /* keeping only the highest ranked entry for each candidate */
drop if party=="v" & year==2013 & district=="aust-agder" /* Venstre ran with identical lists in Vest-Agder and Aust-Agder in 2013. Venstre were 8.5pp away from winning a regular seat in Vest-Agder, and 13.1pp away from winning a regular seat in Aust-Agder, so we keep Vest-Agder */

drop if candidatename_ed=="Anders Lange" & district=="vest-agder" /* Anders Lange was ranked first in both Oslo and Vest-Agder, but performed best in Oslo */
drop if candidatename_ed=="Svend Haakon Jacobsen" & district=="rogaland" /* runs in both districts as hopeless candidate but same rank, arbitrary exclude one */
drop if candidatename_ed=="Kristin Dalehamn" & district=="rogaland" /* runs in both districts as hopeless candidate but same rank, arbitrary exclude one */
drop if candidatename_ed=="Steinar Bastesen" & district=="nord-trøndelag" /* runs in two districts in 2001, and wins in Nordland */

gen main=0
replace main=1 if party=="nkp"|party=="sv"|party=="dna"|party=="sp"|party=="v"|party=="krf"|party=="h"|party=="frp"
keep if main==1

* (iv) Define margin to use, kernel, covariates, clsuter

rename margin votemargin_fiva
global margin votemargin_vote
global kernel triangular // Options: triangular (default), epanechnikov, uniform
global secl cluster // Options: cluster nncluster

* (iv) Estimations

rdrobust elected_F1 votemargin_fiva, p(1) bwselect(mserd) kernel($kernel)  ///
	all vce($secl pid_unique), if inrange(year,1953,1981)
global bw_opt_p1 =e(h_l)
global bbw_opt_p1 =e(b_l)

rdrobust elected_F1 $margin, p(1) bwselect(mserd) kernel($kernel)  ///
	all vce($secl pid_unique), if inrange(year,1953,1981) &  votemargin_fiva!=. 
global bw_sample_fiva =e(h_l)
global bbw_sample_fiva =e(b_l)
	
rdrobust elected_F1 $margin, p(1) bwselect(mserd) kernel($kernel)  ///
	all vce($secl pid_unique), if inrange(year,1953,1981) 
global bw_full_sample =e(h_l)
global bbw_full_sample =e(b_l)

local bw_full_sample_div2=$bw_full_sample/2
local bbw_full_sample_div2=$bbw_full_sample/2
rdrobust  elected_F1 $margin if elected_F1 < ., p(1) h(`bw_full_sample_div2') ///
	b(`bbw_full_sample_div2') kernel($kernel) vce($secl pid_unique) all 

rdplot elected_F1 votemargin_fiva ///
	if -$bw_opt_p1 <= votemargin_fiva & votemargin_fiva<= $bw_opt_p1 & elected_F1 < ., ///
	p(1) kernel($kernel) h($bw_opt_p1 $bw_opt_p1) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig2.pdf", ///
	as(pdf) replace	

rdplot elected_F1 $margin ///
	if -$bw_sample_fiva <= $margin & $margin<= $bw_sample_fiva ///
	& votemargin_fiva!=. & elected_F1 < ., ///
	p(1) kernel($kernel) h($bw_sample_fiva $bw_sample_fiva) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig2.pdf", ///
	as(pdf) replace	
	
rdplot elected_F1 $margin ///
	if -$bw_full_sample <= $margin & $margin<= $bw_full_sample & elected_F1 < ., ///
	p(1) kernel($kernel) h($bw_full_sample $bw_full_sample) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off)) 
graph export "$path\05_Texts_and_presns\01_Running_variable\figures\\`margin'_`z'_Fig2.pdf", ///
	as(pdf) replace	
	
sum elected if -$bw_sample_fiva <= $margin & $margin<= $bw_sample_fiva & votemargin_fiva!=.
sum elected if -$bw_sample_fiva <= $margin & $margin<= $bw_sample_fiva 

gen votemargin_vote_abs=abs(votemargin_vote)
sum votemargin_vote_abs if votemargin_fiva!=. & votemargin_vote>=0
sum votemargin_vote_abs if votemargin_fiva==. & votemargin_vote>=0, d
um votemargin_vote_abs if votemargin_fiva!=. & votemargin_vote<0
sum votemargin_vote_abs if votemargin_fiva==. & votemargin_vote<0

/* Trickkiste

gen sample_fiva_smith=0
replace sample_fiva_smith=1 if $margin>=$bw_opt_p1 & $margin<=$bw_opt_p1

gen votemargin_vote_elected=elected*votemargin_vote
gen margin_elected=elected*margin

gen weight=1-abs(votemargin_vote/$bw_opt_p1)
gen weight2=1-abs(margin/$bw_opt_p1)

reg elected_F1 elected margin margin_elected ///
[pweight=weight2] if sample_fiva_smith==1

reg elected_F1 elected votemargin_vote votemargin_vote_elected ///
[pweight=weight] if sample_fiva_smith==1

rdrobust elected_F1 votemargin_vote, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl pid_unique) all, if   inrange(year,1953,1981)

rdrobust elected_F1 $margin, p(1) bwselect(mserd) kernel($kernel)  ///
	vce($secl pid_unique) all, if inrange(year,1953,1981)
	
/*
rdplot `z' $margin if elected_F1 < ., ///
	kernel($kernel) binselect(qs) ///
	graph_options(ylabel(, angle(horizontal) gsty(dot)) ///
	graphregion(fcolor(white) lcolor(white)) legend(off) ///
	ti(`label', position(11) orientation(horizontal) si(huge))) ///
	plotregion(lcolor(black))

