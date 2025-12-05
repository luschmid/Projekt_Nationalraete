*******************************************
** REPLICATION FILES FOR
** Curto-Grau, Marta, Solé-Ollé, Albert, and Sorribas-Navarro, Pilar. 2018. "Does electoral competition curb party favoritism?" AEJ: Applied Economics
** Tables in the Appendix
*******************************************

cd "...\"

clear all
set more off
set matsize 800


*************
*TABLE A1 and A2. Summary statistics
*************

use "db_main.dta", clear

su tk if ab==0

su ab  dist1 dab debt  tipo vcp  pob pob_mes6591 pob_5_1491 extr_noeu91 unempl income density localp  regre regde meanden educ2 presscirc regterm ecs1 ecs2 ecs3


use "db_full.dta", clear
su abcd bloc



************
*TABLE A10. Balance test
************
clear all

use "db_main.dta", clear

global regFE "i.codccaa"
global controls1 "dist1_ecs1 dist1 ecs1"

estimates clear 
foreach y in debt tipo vcp pob  density pob_mes6591 pob_5_1491 extr_noeu91 unempl income  {

rdrobust `y' dist1  , fuzzy(ab) vce(cluster codiine) kernel(uni)  //Select optimal bandwithd
eststo r`y':  reg `y' dab dist1 vda $regFE  if abs(dist1)<e(h_r), vce(cluster codiine)
}
	
esttab r*, keep(dab) se

estimates clear	
foreach y in presscirc {
rdrobust `y' dist1  , fuzzy(ab) vce(cluster cprov) kernel(uni) //Select optimal bandwithd
eststo r`y': reg `y' dab dist1 vda $regFE  if abs(dist1)<e(h_r), vce(cluster cprov)
 }	
esttab r*, keep(dab) se

 estimates clear
 foreach y in regre regde meanden educ2  regterm ecs1{

rdrobust `y' dist1  , fuzzy(ab) vce(cluster codccaa) kernel(uni) 
eststo r`y': reg `y' dab dist1 vda $regFE  if abs(dist1)<e(h_r), vce(cluster codccaa)

 }	
  
esttab r*, keep(dab) se


************
*TABLE A11. LATE with controls
************

use "db_main.dta", clear

global controls "lpob density debt vcp tipo"
global regFE "i.codccaa"

bysort codiine: center tk ab dt1 dt2  lpob density debt vcp tipo //For DiD
global controlsc "c_lpob c_density c_debt c_vcp c_tipo"


//Panel A (2nd stage) 

estimates clear

eststo r1: qui ivregress 2sls tk (ab vsa vsa2 = dab vda vda2) dist1 dist2 $regFE $controls, vce(cluster codiine) 
qui ivregress 2sls tk (ab vsa vsa2 = dab vda vda2) dist1 dist2  i.codccaa  $controls, vce(cluster codccaa) 
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r2: qui ivregress 2sls tk (ab vsa = dab vda) dist1       $regFE   $controls   if  abs(dist1)<0.386, vce(cluster codiine) 
qui ivregress 2sls tk (ab vsa = dab vda) dist1    i.codccaa    $controls     if  abs(dist1)<0.386, vce(cluster codccaa) 
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r3: qui ivregress 2sls tk (ab vsa = dab vda) dist1       $regFE  $controls    if  abs(dist1)<0.193, vce(cluster codiine) 
qui ivregress 2sls tk (ab vsa = dab vda) dist1    i.codccaa  $controls       if  abs(dist1)<0.193, vce(cluster codccaa) 
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r4: qui ivregress 2sls tk (ab vsa = dab vda) dist1      $regFE   $controls  if  abs(dist1)<0.0965, vce(cluster codiine) 
qui ivregress 2sls tk (ab vsa = dab vda) dist1     i.codccaa   $controls     if  abs(dist1)<0.0965, vce(cluster codccaa) 
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r5: qui ivregress 2sls tk (ab vsa = dab vda) dist1      $regFE   $controls if  abs(dist1)<0.048, vce(cluster codiine) 
qui ivregress 2sls tk (ab vsa = dab vda) dist1      i.codccaa   $controls    if  abs(dist1)<0.048, vce(cluster codccaa) 
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r6: qui reg c_tk c_ab   c_dt1 c_dt2 $controlsc, cluster(codiine) 
qui reg c_tk c_ab c_dt1 c_dt2  $controlsc, cluster(codccaa) 
boottest c_ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

esttab r1 r2 r3 r4 r5 r6, ar2 keep(ab c_ab ) se


//Panel B (1st stage) - with controls

estimates clear

eststo r1: qui reg ab dab dist1 dist2 vda vda2 $regFE $controls , vce(cluster codiine)
reg ab dab dist1 dist2 vda vda2 $regFE $controls , vce(cluster codccaa)
boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r2: qui reg ab dab dist1 vda $regFE $controls if abs(dist1)<0.386, vce(cluster codiine)
qui reg ab dab dist1 vda $regFE $controls  if abs(dist1)<0.386, vce(cluster codccaa)
boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r3: qui reg ab dab dist1 vda $regFE  $controls  if abs(dist1)<0.193, vce(cluster codiine)
qui reg ab dab dist1 vda $regFE $controls  if abs(dist1)<0.193, vce(cluster codccaa)
boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull


eststo r4: qui reg ab dab dist1  vda $regFE $controls  if  abs(dist1)<0.0965, vce(cluster codiine)
qui reg ab dab dist1  vda $regFE $controls  if  abs(dist1)<0.0965, vce(cluster codccaa)
boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull


eststo r5: qui reg ab dab dist1  vda $regFE $controls  if abs(dist1)<0.048,  vce(cluster codiine)
qui reg ab dab dist1  vda $regFE $controls  if abs(dist1)<0.048,  vce(cluster codccaa)
boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

esttab r1 r2 r3 r4 r5 , ar2 keep(dab) se



**********
**TABLE A12 - Polynomial Order
**********

ge dist3=dist1^3
ge vda3=dist3*dab
ge vsa3=dist3*ab

//Panel A (2nd stage) 
estimates clear

eststo r1: qui ivregress 2sls tk (ab vsa = dab vda) dist1  $regFE , vce(cluster codiine)
qui ivregress 2sls tk (ab vsa = dab vda) dist1  $regFE , vce(cluster codccaa)
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r2: qui ivregress 2sls tk (ab vsa vsa2  = dab vda vda2) dist1 dist2  $regFE, vce(cluster codiine)
qui ivregress 2sls tk (ab vsa vsa2  = dab vda vda2) dist1 dist2  $regFE, vce(cluster codccaa)
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull
  
eststo r3: qui ivregress 2sls tk (ab vsa vsa2 vsa3 = dab vda vda2 vda3) dist1 dist2 dist3 $regFE, vce(cluster codiine)
qui ivregress 2sls tk (ab vsa vsa2 vsa3 = dab vda vda2 vda3) dist1 dist2 dist3 $regFE, vce(cluster codccaa)
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

esttab r1 r2 r3  ,  keep(ab) se 


//Panel B (1st stage) 
eststo r1: qui reg ab dab dist1 vda  $regFE, vce(cluster codiine)
qui reg ab dab dist1 vda  $regFE, vce(cluster codccaa)
boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull
 
eststo r2:  qui reg ab dab dist1 dist2 vda vda2  $regFE, vce(cluster codiine)
qui reg ab dab dist1 dist2 vda vda2  $regFE, vce(cluster codccaa)
boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull
 
eststo r3: qui  reg ab dab dist1 dist2 dist3 vda vda2 vda3 $regFE, vce(cluster codiine)
qui  reg ab dab dist1 dist2 dist3 vda vda2 vda3 $regFE, vce(cluster codccaa)
boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull
 
esttab r1 r2 r3  , r2 keep(dab) se 
 

 
************
*TABLE A13. Placebo
************
use "db_main.dta", clear

preserve
keep if t<4

bysort codiine: center tk ab  dt1 lead_ab //FOR DID

eststo r1: qui reg c_tk c_ab c_lead_ab  c_dt1 if t<4, cluster(codiine)
boottest c_ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull
boottest c_lead_ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

esttab r1 , r2 keep(c_ab c_lead_ab) se 

restore


preserve
keep if t>2

bysort codiine: center tk ab  dt2  lag_ab //FOR DID

eststo r2: qui reg c_tk c_ab c_lag_ab  c_dt2 , cluster(codiine)
boottest c_ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull
boottest c_lag_ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

esttab r2 , r2 keep(c_ab c_lag_ab) se 
restore
 



************
*TABLE A14. Robustness
************

clear all

use "db_full.dta", clear

keep if tk!=.

global controls1 "dist1_ecs1  dist1 ecs1"

global regFE "i.codccaa"

//Col. 1 - Alternative forcing var
preserve

keep if ab!=. & dist_alt!=.

su compreg1
gen ecs1 = (compreg1 - r(mean))

ge esas1= ecs1*ab
ge edas1_bis= ecs1*dab_bis
ge vsalt_ecs1=vsa_alt*ecs1
ge vdalt_ecs1=vda_alt*ecs1
ge distalt_ecs1=dist_alt*ecs1


rdrobust tk dist_alt  , fuzzy(ab) vce(cluster codiine) kernel(uni) // Get optimal bandwidth
local b= e(h_r)

eststo dist_L: qui ivregress 2sls tk (ab vsa_alt  = dab_bis vda_alt ) dist_alt $regFE  if  abs(dist_alt)<`b', vce(cluster codiine)
qui ivregress 2sls tk (ab vsa_alt  = dab_bis vda_alt ) dist_alt $regFE  if  abs(dist_alt)<`b', vce(cluster codccaa)
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonul

eststo dist_H: qui ivregress 2sls tk (dca_ab* esas1 vsalt_ecs1  dca_vsa*  = dca_dbis* edas1_bis vdalt_ecs1  dca_vd_alt* ) dist_alt  distalt_ecs1 ecs1 $regFE   if  abs(dist_alt)<`b', vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul

esttab dist_L dist_H, r2 keep(ab esas1 ecs1 ) se 
restore


//Col. 2 & 7 -  No regional parties
preserve

keep if ab!=. & dist1!=. & reg==1

su compreg1
gen ecs1 = (compreg1 - r(mean))

ge esas1= ecs1*ab
ge edas1= ecs1*dab
ge vsa_ecs1=vsa*ecs1
ge vda_ecs1=vda*ecs1
ge dist1_ecs1=dist1*ecs1


bysort codiine: center tk ab  dt1 dt2 esas1 ecs1 dca* //For DID
rdrobust tk dist1  , fuzzy(ab) vce(cluster codiine) kernel(uni) // Get optimal bandwidth
local b= e(h_r)

eststo reg_L: qui ivregress 2sls tk (ab vsa  = dab vda ) dist1  $regFE  if  abs(dist1)<`b', vce(cluster codiine)
qui ivregress 2sls tk (ab vsa  = dab vda ) dist1  $regFE  if  abs(dist1)<`b', vce(cluster codccaa)
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonul

eststo reg_H: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1  dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* ) $controls1   $regFE if  abs(dist1)<`b' & dca5==0, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul

eststo reg_LD: qui reg c_tk c_ab   c_dt1 c_dt2 , cluster(codiine) 
qui reg c_tk c_ab   c_dt1 c_dt2 , cluster(codccaa) 
boottest c_ab ,  seed(987654321) boottype(score) reps(10000) nograph nonul

eststo reg_HD: qui regress c_tk c_dca_ab* c_esas1 c_ecs1 c_dt1 c_dt2  , cluster(codccaa)
boottest c_esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul

esttab reg_L reg_H reg_LD reg_HD, r2 keep(ab esas1 c_ab c_esas1 ecs1 c_ecs1) se 

restore	


//Col. 3 & 8 -  No local parties
preserve

keep if ab!=. & dist1!=. & indep==1

su compreg1
gen ecs1 = (compreg1 - r(mean))
ge esas1= ecs1*ab
ge edas1= ecs1*dab
ge vsa_ecs1=vsa*ecs1
ge vda_ecs1=vda*ecs1
ge dist1_ecs1=dist1*ecs1

bysort codiine: center tk ab  dt1 dt2 esas1 ecs1 dca* //FOR DID

rdrobust tk dist1  , fuzzy(ab) vce(cluster codiine) kernel(uni) // Get optimal bandwidth
local b= e(h_r)

eststo indep_L: qui ivregress 2sls tk (ab vsa = dab vda ) dist1 $regFE  if  abs(dist1)<`b', vce(cluster codiine)
qui ivregress 2sls tk (ab vsa = dab vda ) dist1 $regFE  if  abs(dist1)<`b', vce(cluster codccaa)
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo indep_H: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1  dca_vsa* = dca_dab* edas1 vda_ecs1  dca_vda* ) $controls1   $regFE if abs(dist1)<`b' & dca5==0 & dca6==0, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonull
boottest ecs1,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo indep_LD: qui reg c_tk c_ab   c_dt1 c_dt2 , cluster(codiine) 
qui reg c_tk c_ab   c_dt1 c_dt2 , cluster(codccaa) 
boottest c_ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo indep_HD: qui regress c_tk c_dca_ab* c_esas1 c_ecs1 c_dt1 c_dt2  , cluster(codccaa)
boottest c_esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonull
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonull

esttab indep_L indep_H indep_LD indep_HD, r2 keep(ab esas1 c_ab c_esas1 ecs1 c_ecs1) se 

restore	


//Col. 4 & 9 -  Partner alignment
preserve

keep if abcd!=. & dist1!=.

su compreg1
gen ecs1 = (compreg1 - r(mean))
ge esas1= ecs1*ab
ge edas1= ecs1*dab

ge vsa_ecs1=vsa*ecs1
ge vda_ecs1=vda*ecs1
ge dist1_ecs1=dist1*ecs1
ge esas1_bis= ecs1*abcd

ge vsbcd_ecs1=vsbcd*ecs1

bysort codiine: center tk abcd  dt1 dt2 esas1_bis ecs1 dca* //FOR DID

rdrobust tk dist1  , fuzzy(abcd) vce(cluster codiine) kernel(uni) // Get optimal bandwidth
local b = e(h_r)

eststo abcd_L: qui ivregress 2sls tk (abcd vsbcd = dab vda ) dist1  $regFE if abs(dist1)<`b', vce(cluster codiine)
qui ivregress 2sls tk (abcd vsbcd = dab vda ) dist1  $regFE if abs(dist1)<`b', vce(cluster codccaa)
boottest abcd ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo abcd_H:  ivregress 2sls tk (dca_bcd* esas1_bis vsbcd_ecs1 dca_vsbcd* = dca_dab* edas1 vda_ecs1  dca_vda* ) $controls1 $regFE  if abs(dist1)<`b', vce(cluster codccaa)
boottest esas1_bis ,  seed(987654321) boottype(score) reps(10000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo abcd_LD: qui reg c_tk c_abcd   c_dt1 c_dt2 , cluster(codiine) 
qui reg c_tk c_abcd   c_dt1 c_dt2 , cluster(codccaa) 
boottest c_abcd ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo abcd_HD: qui regress c_tk c_dca_bcd* c_esas1_bis c_ecs1 c_dt1 c_dt2, cluster(codccaa)
boottest c_esas1_bis ,  seed(987654321) boottype(score) reps(10000) nograph nonull
boottest  c_ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonull

esttab abcd_L abcd_H abcd_LD abcd_HD, r2 keep(abcd esas1_bis c_abcd c_esas1_bis ecs1 c_ecs1) se 

restore	


//Col. 5 & 10 - Bloc alignment
preserve

keep if bloc!=. & dist1!=.

su compreg1
gen ecs1 = (compreg1 - r(mean))

ge esas1= ecs1*ab
ge edas1= ecs1*dab
ge vsa_ecs1=vsa*ecs1
ge vda_ecs1=vda*ecs1
ge dist1_ecs1=dist1*ecs1
ge esas1_bisbis=ecs1*bloc
ge edas1_bis= ecs1*dab_bis
ge vsbloc_ecs1=vsbloc*ecs1

bysort codiine: center tk bloc  dt1 dt2 esas1_bisbis ecs1 dca* //FOR DID

rdrobust tk dist1  , fuzzy(bloc) vce(cluster codiine) kernel(uni) // Get optimal bandwidth
local b=e(h_r)

eststo bloc_L: qui ivregress 2sls tk (bloc vsbloc  = dab vda ) dist1 $regFE   if abs(dist1)<`b', vce(cluster codiine)
qui ivregress 2sls tk (bloc vsbloc  = dab vda ) dist1 $regFE   if abs(dist1)<`b', vce(cluster codccaa)
boottest bloc ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo bloc_H: qui ivregress 2sls tk (dca_bloc* esas1_bisbis vsbloc_ecs1  dca_vsbloc*  = dca_dab* edas1 vda_ecs1 dca_vda*) $controls1 $regFE   if abs(dist1)<`b', vce(cluster codccaa)
boottest esas1_bisbis ,  seed(987654321) boottype(score) reps(10000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo bloc_LD: qui reg c_tk c_bloc   c_dt1 c_dt2 , cluster(codiine) 
qui reg c_tk c_bloc   c_dt1 c_dt2 , cluster(codccaa) 
boottest c_bloc ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo bloc_HD: qui regress c_tk c_dca_bloc* c_esas1_bisbis c_ecs1 c_dt1 c_dt2, cluster(codccaa)
boottest c_esas1_bisbis ,  seed(987654321) boottype(score) reps(10000) nograph nonull
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonull

esttab bloc_L bloc_H bloc_LD bloc_HD, r2 keep(bloc esas1_bisbis c_bloc c_esas1_bisbis ecs1 c_ecs1) se 
restore	


//Col. 6 & 11 - Concurrent

preserve

keep if ab!=. & dist1!=. & concurrent==1

su compreg1
gen ecs1 = (compreg1 - r(mean))
ge esas1= ecs1*ab
ge edas1= ecs1*dab
ge vsa_ecs1=vsa*ecs1
ge vda_ecs1=vda*ecs1
ge dist1_ecs1=dist1*ecs1

bysort codiine: center tk ab  dt1 dt2 esas1 ecs1 dca* //FOR DID
rdrobust tk dist1  if concurrent==1, fuzzy(ab) vce(cluster codiine) kernel(uni) // Get optimal bandwidth
local b=e(h_r)

eststo conc_L: qui ivregress 2sls tk (ab vsa  = dab vda ) dist1  $regFE  if  abs(dist1)<`b', vce(cluster codiine)
qui ivregress 2sls tk (ab vsa  = dab vda ) dist1  $regFE  if  abs(dist1)<`b', vce(cluster codccaa)
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonul

eststo conc_H: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1  dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* ) $controls1   $regFE if  abs(dist1)<`b' & dca6==0, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul


eststo conc_LD: qui reg c_tk c_ab   c_dt1 c_dt2 , cluster(codiine) 
qui reg c_tk c_ab   c_dt1 c_dt2 , cluster(codccaa) 
boottest c_ab ,  seed(987654321) boottype(score) reps(10000) nograph nonul

eststo conc_HD: qui regress c_tk c_dca_ab* c_esas1 c_ecs1 c_dt1 c_dt2  , cluster(codccaa)
boottest c_esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul

esttab conc_L conc_H conc_LD conc_HD, r2 keep(ab esas1 c_ab c_esas1 ecs1 c_ecs1) se 
restore	



************
*TABLE A15. Alternative measures of the Regional seat margin 
************
use "db_full.dta", clear

reg tk ab dist1
keep if e(sample)

global controls1 "dist1_ecs1 dist1  ecs1"
global controls2 "dist1_ecs1 dist2_ecs1 dist1 dist2  ecs1"

global regFE "i.codccaa"

//Panel A. Regional competition = President’s v. Opposition’s blocs
preserve

su compreg2
gen ecs1 = (compreg2 - r(mean))
su ecs1
ge esas1= ecs1*ab
ge edas1= ecs1*dab

ge vsa_ecs1=vsa*ecs1
ge vsa2_ecs1=vsa2*ecs1

ge vda_ecs1=vda*ecs1
ge vda2_ecs1=vda2*ecs1

ge dist1_ecs1=dist1*ecs1
ge dist2_ecs1=dist2*ecs1


bysort codiine: center tk ab  dt1 dt2 esas1 ecs1 dca* //FOR DID

eststo H: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1 vsa2_ecs1 dca_vsa* dca_2vsa* = dca_dab* edas1 vda_ecs1 vda2_ecs1 dca_vda* dca_2vda*) $controls2   $regFE, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul

eststo H2: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1  dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* )  $controls1 $regFE if abs(dist1)<0.193 , vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul

eststo HD: qui regress c_tk c_dca_ab* c_esas1 c_ecs1 c_dt1 c_dt2  , cluster(codccaa) 
boottest c_esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul

esttab H H2 HD, r2 keep(esas1 ecs1 c_esas1 c_ecs1) se 

restore	


//Panel B. Regional competition = Main 2 parties considered
preserve

su compreg3
gen ecs1 = (compreg3 - r(mean))

ge esas1= ecs1*ab
ge edas1= ecs1*dab

ge vsa_ecs1=vsa*ecs1
ge vsa2_ecs1=vsa2*ecs1

ge vda_ecs1=vda*ecs1
ge vda2_ecs1=vda2*ecs1

ge dist1_ecs1=dist1*ecs1
ge dist2_ecs1=dist2*ecs1

bysort codiine: center tk ab  dt1 dt2 esas1 ecs1 dca* //FOR DID

eststo H: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1 vsa2_ecs1 dca_vsa* dca_2vsa* = dca_dab* edas1 vda_ecs1 vda2_ecs1 dca_vda* dca_2vda*) $controls2   $regFE, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul

eststo H2: ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1  dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* )  $controls1 $regFE if abs(dist1)<0.193 , vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul

eststo HD: qui regress c_tk c_dca_ab* c_esas1 c_ecs1 c_dt1 c_dt2  , cluster(codccaa) 
boottest c_esas1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(10000) nograph nonul

esttab H H2 HD, r2 keep(esas1 ecs1 c_esas1 c_ecs1) se 

restore	


************
*TABLE A16. Controlling for time varying covariates. Global RD estimates.
************
use "db_main.dta", clear

global controls2 "dist1_ecs1 dist2_ecs1 dist1 dist2  ecs1"
global regFE "i.codccaa"	

estimates clear
bysort codiine: center tk ab esas1 ecs1 dt1 dt2  resa regre desa regde denssa meandens termssa regterm educ2sa educ2 presssa presscirc dca_ab* //For DiD

// RD ESTIMATES, 2nd order pol
*Regional revenues
eststo r1: qui ivregress 2sls tk (dca_ab* esas1 resa vsa_ecs1 vsa2_ecs1 vsa_re vsa2_re dca_vsa* dca_2vsa* = dca_dab* edas1 vda_ecs1 vda2_ecs1 dca_vda* dca_2vda* reda  vda_re vda2_re) $controls2 $regFE  dist1_regre dist2_regre regre, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest resa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Regional debt
eststo r2: qui ivregress 2sls tk (dca_ab* esas1 desa vsa_ecs1 vsa2_ecs1 vsa_de vsa2_de dca_vsa* dca_2vsa*= dca_dab* edas1 vda_ecs1 vda2_ecs1 dca_vda* dca_2vda* deda  vda_de vda2_de) $controls2 $regFE  dist1_regde dist2_regde regde, vce(cluster codccaa) 
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest desa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Regional population density
eststo r3: qui ivregress 2sls tk (dca_ab* esas1 denssa vsa_ecs1 vsa2_ecs1 vsa_dens vsa2_dens dca_vsa* dca_2vsa* = dca_dab* edas1  densda vda_ecs1 vda2_ecs1 vda_dens vda2_dens dca_vda* dca_2vda*) $controls2 $regFE  dist1_meanden dist2_meanden meanden, vce(cluster codccaa) 
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest denssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Tenure in office
eststo r4: qui ivregress 2sls tk (dca_ab* esas1 termssa vsa_ecs1 vsa2_ecs1 vsa_te vsa2_te dca_vsa* dca_2vsa* = dca_dab* edas1 vda_ecs1 vda2_ecs1 dca_vda* dca_2vda* termsda vda_te vda2_te) $controls2 $regFE  dist1_regterm dist2_regterm regterm, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest termssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Press
eststo r5:qui  ivregress 2sls tk (dca_ab* esas1 presssa vsa_ecs1 vsa2_ecs1 vsa_pr vsa2_pr dca_vsa* dca_2vsa* = dca_dab* edas1 vda_ecs1 vda2_ecs1 dca_vda* dca_2vda* pressda vda_pr vda2_pr) $controls2 $regFE  dist1_presscirc dist2_presscirc presscirc, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest presssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*% Primary and secondary educ
eststo r6: qui ivregress 2sls tk (dca_ab* esas1 educ2sa vsa_ecs1 vsa2_ecs1 vsa_edu vsa2_edu dca_vsa* dca_2vsa* = dca_dab* edas1 vda_ecs1 vda2_ecs1 dca_vda* dca_2vda* educ2da vda_edu vda2_edu) $controls2 $regFE  dist1_educ dist2_educ educ2, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest educ2sa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

* All
eststo r7: qui ivregress 2sls tk (dca_ab* esas1 resa desa denssa termssa presssa educ2sa  vsa_ecs1 vsa2_ecs1 vsa_re vsa2_re vsa_dens vsa2_dens vsa_te vsa2_te vsa_pr vsa2_pr vsa_edu vsa2_edu dca_vsa* dca_2vsa* = dca_dab* edas1 reda deda densda vda_ecs1 vda2_ecs1 vda_re vda2_re vda_dens vda2_dens dca_vda* dca_2vda*  termsda pressda educ2da vda_te vda2_te vda_pr vda2_pr  vda_edu vda2_edu   ) $controls2 $regFE  dist1_regre dist2_regre regre regde dist1_regde dist2_regde meanden dist1_meanden dist2_meanden dist1_regterm dist2_regterm regterm dist1_presscirc dist2_presscirc presscirc  dist1_educ dist2_educ educ2  , vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest resa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest desa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest denssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest termssa  ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest presssa,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest educ2sa,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

esttab r1 r2 r3 r4 r5 r6  r7,  keep(esas1 resa desa denssa termssa presssa educ2sa ecs1) se 


// DID ESTIMATES

*Regional revenues
eststo r1: qui regress c_tk c_dca_ab* c_esas1 c_resa c_ecs1 c_dt1 c_dt2 c_regre , cluster(codccaa) 
boottest c_esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_resa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Regional debt
eststo r2: qui regress c_tk c_dca_ab* c_esas1 c_desa c_ecs1 c_dt1 c_dt2 c_regde , cluster(codccaa) 
boottest c_esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_desa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Regional population density
eststo r3: qui regress c_tk c_dca_ab* c_esas1 c_denssa c_ecs1 c_dt1 c_dt2 c_meanden , cluster(codccaa) 
boottest c_esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_denssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Tenure in office
eststo r4: qui regress c_tk c_dca_ab* c_esas1 c_termssa c_ecs1 c_dt1 c_dt2 c_regterm , cluster(codccaa) 
boottest c_esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_termssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Press
eststo r5: qui regress c_tk c_dca_ab* c_esas1 c_presssa c_ecs1 c_dt1 c_dt2 c_presscirc , cluster(codccaa) 
boottest c_esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_presssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Education
eststo r6:  qui regress c_tk c_dca_ab* c_esas1 c_educ2sa c_ecs1 c_dt1 c_dt2 c_educ2, cluster(codccaa) 
boottest c_esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_educ2sa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*All
eststo r7:  qui regress c_tk c_dca_ab* c_esas1 c_resa  c_desa  c_denssa  c_termssa  c_presssa  c_educ2sa  c_ecs1 c_dt1 c_dt2 c_regre c_regde c_meanden c_regterm c_presscirc c_educ2 , cluster(codccaa) 
boottest c_esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_resa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_desa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_denssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_termssa  ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_presssa,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_educ2sa,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest c_ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

esttab r1 r2 r3 r4 r5 r6  r7,  keep(c_esas1 c_resa c_desa c_denssa c_termssa c_presssa c_educ2sa c_ecs1) se 
