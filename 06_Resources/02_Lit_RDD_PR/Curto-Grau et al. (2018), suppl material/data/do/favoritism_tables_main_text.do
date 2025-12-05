*******************************************
** REPLICATION FILES FOR
** Curto-Grau, Marta, Solé-Ollé, Albert, and Sorribas-Navarro, Pilar. 2018. "Does electoral competition curb party favoritism?" AEJ: Applied Economics
** Tables in the Main Text
*******************************************

cd "..."

clear all
set more off
set matsize 800


use "db_main.dta", clear

global controls2 "dist1_ecs1 dist2_ecs1 dist1 dist2  ecs1"
global controls1 "dist1_ecs1  dist1 ecs1"

global regFE "i.codccaa"


bysort codiine: center tk ab esas1 ecs1 dt1 dt2  resa regre desa regde denssa meandens termssa regterm educ2sa educ2 presssa presscirc dca_ab* //For DiD

rdrobust tk dist1  , fuzzy(ab) vce(cluster codiine) kernel(uni) // Get optimal bandwidth


**********************
*TABLE 1. Panel A (2nd stage) - LATE
**********************

estimates clear

eststo r1: qui ivregress 2sls tk (ab vsa vsa2 = dab vda vda2) dist1 dist2 $regFE, vce(cluster codiine) 
qui ivregress 2sls tk (ab vsa vsa2 = dab vda vda2) dist1 dist2  i.codccaa, vce(cluster codccaa) 
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r2: qui ivregress 2sls tk (ab vsa = dab vda) dist1       $regFE     if  abs(dist1)<0.386, vce(cluster codiine) 
qui ivregress 2sls tk (ab vsa = dab vda) dist1    i.codccaa        if  abs(dist1)<0.386, vce(cluster codccaa) 
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r3: qui ivregress 2sls tk (ab vsa = dab vda) dist1       $regFE     if  abs(dist1)<0.193, vce(cluster codiine) 
qui ivregress 2sls tk (ab vsa = dab vda) dist1    i.codccaa        if  abs(dist1)<0.193, vce(cluster codccaa) 
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r4: qui ivregress 2sls tk (ab vsa = dab vda) dist1      $regFE    if  abs(dist1)<0.0965, vce(cluster codiine) 
qui ivregress 2sls tk (ab vsa = dab vda) dist1     i.codccaa       if  abs(dist1)<0.0965, vce(cluster codccaa) 
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r5: qui ivregress 2sls tk (ab vsa = dab vda) dist1      $regFE   if  abs(dist1)<0.048, vce(cluster codiine) 
qui ivregress 2sls tk (ab vsa = dab vda) dist1      i.codccaa      if  abs(dist1)<0.048, vce(cluster codccaa) 
boottest ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

eststo r6: qui reg c_tk c_ab   c_dt1 c_dt2, cluster(codiine) 
qui reg c_tk c_ab c_dt1 c_dt2 , cluster(codccaa) 
boottest c_ab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

esttab r1 r2 r3 r4 r5 r6 ,  keep(ab c_ab) se



**********************
*Table 1. Panel b (1st stage)
**********************
estimates clear

	eststo r1: qui reg ab dab dist1 dist2 vda vda2 $regFE, vce(cluster codiine)
	qui reg ab dab dist1 dist2 vda vda2 $regFE, vce(cluster codccaa)
	boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

	eststo r2: qui reg ab dab dist1 vda $regFE if abs(dist1)<0.386, vce(cluster codiine)
	qui reg ab dab dist1 vda $regFE if abs(dist1)<0.386, vce(cluster codccaa)
	boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

	eststo r3: qui reg ab dab dist1 vda $regFE if abs(dist1)<0.193, vce(cluster codiine)
	qui reg ab dab dist1 vda $regFE if abs(dist1)<0.193, vce(cluster codccaa)
	boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull


	eststo r4: qui reg ab dab dist1  vda $regFE if  abs(dist1)<0.0965, vce(cluster codiine)
	qui reg ab dab dist1  vda $regFE if  abs(dist1)<0.0965, vce(cluster codccaa)
	boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull


	eststo r5: qui reg ab dab dist1  vda $regFE if abs(dist1)<0.048,  vce(cluster codiine)
	qui reg ab dab dist1  vda $regFE if abs(dist1)<0.048,  vce(cluster codccaa)
	boottest dab ,  seed(987654321) boottype(score) reps(10000) nograph nonull

esttab r1 r2 r3 r4 r5 , ar2 keep(dab ) se



**********************
*Table 2. HLATE
**********************

estimates clear
	
	eststo r1: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1 vsa2_ecs1 dca_vsa* dca_2vsa* = dca_dab* edas1 vda_ecs1 vda2_ecs1 dca_vda* dca_2vda*) $controls2   $regFE, vce(cluster codccaa)
	boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
	boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
	
	eststo r2: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1  dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* )  $controls1 $regFE if abs(dist1)<0.386 , vce(cluster codccaa) 
	boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
	boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
		
	eststo r3: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1  dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* )  $controls1 $regFE if abs(dist1)<0.193 , vce(cluster codccaa)
	boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
	boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
	
	eststo r4: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1  dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* )  $controls1  $regFE if abs(dist1)<0.0965 & dca5==0 & dca6==0, vce(cluster codccaa) 
	boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
	boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
	
	eststo r5: qui ivregress 2sls tk  (dca_ab* esas1 vsa_ecs1  dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* )  $controls1   $regFE if abs(dist1)<0.048 & dca5==0 & dca6==0 & dca14==0, vce(cluster codccaa)
	boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
	boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
		
	eststo r6: qui  regress c_tk c_dca_ab* c_esas1 c_ecs1 c_dt1 c_dt2  , cluster(codccaa) 
	boottest c_esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
	boottest c_ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
	
esttab r1 r2 r3 r4 r5 r6, keep(esas1 ecs1 c_esas1  c_ecs1) se
	
	
********************
** Table 3: Controlling for potential confounders
********************

*Regional revenues
eststo r1: qui ivregress 2sls tk (dca_ab* esas1 resa vsa_ecs1 vsa_re dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* reda  vda_re ) $controls1 $regFE  dist1_regre regre if abs(dist1)<0.193, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest resa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Regional debt
eststo r2: qui ivregress 2sls tk (dca_ab* esas1 desa vsa_ecs1 vsa_de dca_vsa* = dca_dab* edas1 vda_ecs1 dca_vda* deda  vda_de ) $controls1 $regFE  dist1_regde regde if abs(dist1)<0.193, vce(cluster codccaa) 
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest desa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Regional population density
eststo r3: qui ivregress 2sls tk (dca_ab* esas1 denssa vsa_ecs1 vsa_dens  dca_vsa*  = dca_dab* edas1  densda vda_ecs1 vda_dens dca_vda* ) $controls1 $regFE  dist1_meanden meanden if abs(dist1)<0.193, vce(cluster codccaa) 
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest denssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Tenure in office
eststo r4: qui ivregress 2sls tk (dca_ab* esas1 termssa vsa_ecs1 vsa_te dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* termsda vda_te ) $controls1 $regFE  dist1_regterm regterm  if abs(dist1)<0.193, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest termssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*Press
eststo r5: qui ivregress 2sls tk (dca_ab* esas1 presssa vsa_ecs1 vsa_pr dca_vsa* = dca_dab* edas1 vda_ecs1 dca_vda* pressda vda_pr ) $controls1 $regFE  dist1_presscirc presscirc if abs(dist1)<0.193, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest presssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*College educated
eststo r6: qui ivregress 2sls tk (dca_ab* esas1 educ2sa vsa_ecs1 vsa_edu dca_vsa*  = dca_dab* edas1 vda_ecs1 dca_vda* educ2da vda_edu ) $controls1 $regFE  dist1_educ educ2 if abs(dist1)<0.193, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest educ2sa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull

*All
eststo r7: qui ivregress 2sls tk (dca_ab* esas1 resa desa denssa termssa presssa educ2sa vsa_ecs1 vsa_re vsa_dens vsa_te vsa_pr vsa_edu dca_vsa* = dca_dab* edas1 reda deda densda vda_ecs1 vda_re vda_dens dca_vda* termsda pressda educ2da  vda_te vda_pr vda_edu ) $controls1 $regFE  dist1_regre regre regde dist1_regde meanden dist1_meanden dist1_regterm regterm dist1_presscirc presscirc  dist1_educ educ2  if abs(dist1)<0.193, vce(cluster codccaa)
boottest esas1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest resa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest desa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest denssa ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest termssa  ,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest presssa,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest educ2sa,  seed(987654321) boottype(score) reps(1000) nograph nonull
boottest ecs1 ,  seed(987654321) boottype(score) reps(1000) nograph nonull


esttab r1 r2 r3 r4 r5 r6 r7, keep(esas1 resa desa denssa termssa presssa educ2sa ecs1) se



////////END
