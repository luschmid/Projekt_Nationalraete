
/**************************************Program definitions************************************************/
/**************************************Program definitions************************************************/

/****************************************************/
/* define program to define minimal distance to seat change          */
/****************************************************/

capture program drop calmindiff
program calmindiff
  args namelist namelist2
  
  foreach parti in `namelist'{
    gen spr`parti' =`parti'proc
	gen m`parti' = `parti'mandat
  }

   foreach parti in `namelist'{
   	gen p`parti'=m`parti'>0 & m`parti'!=.
		gen p1`parti'=m`parti'==1

		}

**define comparison numbers*******
   foreach parti in `namelist'{
		gen j`parti'=spr`parti'*(1-p`parti')/1.4 +spr`parti'*p`parti'/(1+2*m`parti')
		}

   foreach parti in `namelist'{
		gen j2`parti'=spr`parti'*(p1`parti')/1.4 +spr`parti'*(1-p1`parti')/(1+2*(m`parti'-1)) if m`parti'>0
		}

	

***define minimal and maximal comparison numbers*******
	gen maxj1=0 
	 foreach parti in `namelist'{
		replace maxj1= j`parti' if j`parti'>maxj1 & j`parti'!=.
		}

	gen minj2=100000000000000
	foreach parti in `namelist'{
		replace minj2= j2`parti' if j2`parti'<minj2 & j2`parti'!=0.
		}

		*****************ALTERNATIVE 1***************************************************
***define distance to seat gain in own votes*****
	foreach parti in `namelist'{
		gen `parti'diffjminj2= (minj2 -j`parti') * ( (1-p`parti') * 1.4 + p`parti'*(1+2*m`parti') )
		replace `parti'diffjminj2=.  if `parti'diffjminj2 <0
	}

***define distance to seat loss in own votes*****	
	foreach parti in `namelist'{
		gen `parti'diffj2maxj= (j2`parti'-maxj1)* (p1`parti' *1.4 +(1-p1`parti')*(1+2*(m`parti'-1))) if m`parti'>0 & m`parti'!=.
		replace `parti'diffj2maxj=. if `parti'diffj2maxj <0
	}

	foreach parti in `namelist'{
		gen mindiff`parti'p1= `parti'diffjminj2 
		gen mindiff`parti'n1=`parti'diffj2maxj   if m`parti'>0 & m`parti'!=.
	}

*****************ALTERNATIVE 2***************************************************
****define smallest distance in own votes for any party gaining a seat*********
	gen mindiffjminj2=100000000000000
	foreach parti in `namelist'{
		replace mindiffjminj2= `parti'diffjminj2 if `parti'diffjminj2<mindiffjminj2 & `parti'diffjminj2!=0 & `parti'diffjminj2!=.
		}	

***define smallest distance in own votes for any party loosing a seat*********
	gen mindiffj2maxj=100000000000000
	foreach parti in `namelist'{
		replace mindiffj2maxj= `parti'diffj2maxj if `parti'diffj2maxj<mindiffj2maxj & `parti'diffj2maxj!=0 & `parti'diffj2maxj!=.
	}	


*** Define distance for a party gaining a seat through gaining votes such that j=maxj1 and another party loosing votes****
	
	foreach parti in `namelist'{
		gen mindiff`parti'p2= mindiffj2maxj  + (maxj1 -j`parti')*((1-p`parti')*1.4 +p`parti'*(1+2*m`parti'))
		replace mindiff`parti'p2=. if maxj1 < j`parti'
	}

*** Define distance for a party loosing a seat through loosing votes such that j2=minj2 and another party gaining votes****
	foreach parti in `namelist'{
		gen mindiff`parti'n2=mindiffjminj2+ (j2`parti'-minj2 ) * (p1`parti'*1.4 + (1-p1`parti') * (1+2*(m`parti'-1) ) )  if m`parti'>0 & m`parti'!=.
		replace mindiff`parti'n2=. if j2`parti'< minj2
	}

*****************ALTERNATIVE 3***************************************************
		
*** Define how many votes  party2 would need to loose to such that party2j<party1j****
	foreach parti in `namelist'{
		
		foreach parti2 in `namelist2'{ 
			gen `parti'`parti2'jdiff=(j`parti2' -j`parti') *((1-p`parti2')*1.4 +p`parti2'*(1+2*m`parti2')) 
			replace `parti'`parti2'jdiff =. if `parti'`parti2'jdiff < 0
		}
	}
*** Define total vote losses for other parties such that  party1j=j1max*********
 
	foreach parti in `namelist'{
		gen	`parti'sumjdiff=0 
		foreach parti2 in `namelist2'{
			replace `parti'sumjdiff= `parti'sumjdiff +`parti'`parti2'jdiff if `parti'`parti2'jdiff !=.
		}
	}

*** Define how many votes  party2 would need to loose such that party2j2<party1j1****

	foreach parti in `namelist'{
		foreach parti2 in `namelist2'{ 
			gen `parti'`parti2'diffj2j= (j2`parti2' -j`parti') * (p1`parti2'*1.4 +(1-p1`parti2')*(1+2*(m`parti2'-1))) if m`parti2'>0 &  m`parti2'!=.
			replace `parti'`parti2'diffj2j =. if j`parti'==j`parti2' |`parti'`parti2'diffj2j<0
		}
	}

*** Define vote minimal loss for any party2  such that party2j2<party1j1****
	foreach parti in `namelist'{	
	gen `parti'mindiffj2j=100000000000000
		foreach parti2 in `namelist2'{ 
			replace `parti'mindiffj2j= `parti'`parti2'diffj2j if `parti'`parti2'diffj2j<`parti'mindiffj2j & `parti'`parti2'diffj2j!=.
		}	
	}

**** Define combination of minimal loss for any party2  such that party2j2<party1j1 and ****
	**** total vote losses for other parties such that  party1j<j1max to define disance to seat gain*****
	
foreach parti in `namelist'{
	gen mindiff`parti'p3=  `parti'sumjdiff+`parti'mindiffj2j
	}
 


*** Define how many votes  party2 would need to gain to such that party2j2>party1j2****
	foreach parti in `namelist'{	
		foreach parti2 in `namelist2'{ 
			gen `parti'`parti2'j2diff=(j2`parti'-j2`parti2') *(p1`parti2'*1.4 +(1-p1`parti2')*(1+2*(m`parti2'-1))) if m`parti2'>0 & m`parti2'!=. 
			replace `parti'`parti2'j2diff=. if `parti'`parti2'j2diff<0
		}
	}

*** Define total vote gains for other parties such that  party1j2=j2min*********
 
	foreach parti in `namelist'{
		gen	`parti'sumj2diff=0 
		foreach parti2 in `namelist2'{
			replace `parti'sumj2diff= `parti'sumj2diff +`parti'`parti2'j2diff if `parti'`parti2'j2diff !=.
		}
	}

*** Define how many votes  party2 would need to gain such that party2j>party1j2****

	foreach parti in `namelist'{
		foreach parti2 in `namelist2'{ 
			gen `parti'`parti2'diffjj2=(j2`parti' -j`parti2') * ((1-p`parti2')*1.4 + p`parti2'*(1+2*m`parti2')) 
			replace `parti'`parti2'diffjj2 =. if j`parti'==j`parti2' |`parti'`parti2'diffjj2<0
		}
	}
*** Define vote minimal gain for any party2  such that party2j>party1j2****
	foreach parti in `namelist'{	
	gen `parti'mindiffjj2=100000000000000
		foreach parti2 in `namelist2'{ 
			replace `parti'mindiffjj2= `parti'`parti2'diffjj2 if `parti'`parti2'diffjj2<`parti'mindiffjj2 & `parti'`parti2'diffjj2!=.
		}	
	}
	
**** Define combination of minimal gain for any party2  such that party2j>party1j2 and ****
	**** total vote gains for other parties such that  party1j2=j2min to define disance to seat gain*****

	foreach parti in `namelist'{
	gen mindiff`parti'n3= (`parti'sumj2diff+ `parti'mindiffjj2) if m`parti'>0 & m`parti'!=.
	}


*****Define minimal distance to seat change******
	foreach parti in `namelist'{
	egen mindiff`parti'p= rowmin (mindiff`parti'p1 mindiff`parti'p2 mindiff`parti'p3 ) 
	replace mindiff`parti'p=. if mindiff`parti'p<0
	egen mindiff`parti'n= rowmin (mindiff`parti'n1 mindiff`parti'n2 mindiff`parti'n3)  if m`parti'>0 & m`parti'!=. 
	replace mindiff`parti'n=. if mindiff`parti'n<0
	}
end




/****************************************************/
/****************************************************/
/* define program to dummies from minimum distance to change  calculated measured in vote share         */
/****************************************************/

capture program drop dum_mindiffpartycalvs
program dum_mindiffpartycalvs
args namelist namelist2


foreach lim in `namelist2' {  
gen sum_`lim'=0 
   foreach parti in `namelist'{
			gen `parti'diffpcalvs`lim'=  (mindiff`parti'p)<0.0`lim'

			gen `parti'diffncalvs`lim'=  (mindiff`parti'n)<0.0`lim'



			gen `parti'diffcalvs`lim' =(`parti'diffncalvs`lim'-`parti'diffpcalvs`lim')*.5
			replace sum_`lim'= sum_`lim' +`parti'diffcalvs`lim'
			drop `parti'diffncalvs`lim'  `parti'diffpcalvs`lim' 


		}
} 

foreach lim in `namelist2' {   
   foreach parti in `namelist'{




			replace `parti'diffcalvs`lim' =0 if sum_`lim'==0


		}
} 
end








/****************************************************/
/* define program to generated weighted vote share           */
/****************************************************/

capture program drop genws
program genws
args namelist 


foreach parti in `namelist'{
	gen `parti'p=`parti'mandat>0 
}
foreach parti in `namelist'{

	gen `parti'procpos=`parti'proc*`parti'p
}
egen sumprocpos = rowtotal (mprocpos-lok4procpos)

foreach parti  in `namelist'{
	gen `parti'procposk=`parti'proc/sumprocpos
	gen `parti'procposksq=`parti'procposk^2
	gen `parti'procposkcub=`parti'procposk^3
	gen `parti'procposk4p=`parti'procposk^4
}

foreach parti  in  `namelist'{
	gen `parti'procposkam=`parti'procposk*antalmandat
}
foreach parti  in  `namelist'{
	gen `parti'procposksqam=`parti'procposksq*antalmandat
}
foreach parti  in  `namelist'{
	gen `parti'procposcubkam=`parti'procposkcub*antalmandat
}
foreach parti  in  `namelist'{
	gen `parti'procposk4pam=`parti'procposk4p*antalmandat
	drop `parti'p `parti'procpos  `parti'procposk-`parti'procposk4p
}

drop sumprocpos

end
/****************************************************/
/****************************************************/
/* define program to generate alertantiv  weighted vote share           */
/****************************************************/

capture program drop altgenws
program altgenws
args namelist 

foreach parti in `namelist'{
gen `parti'p=`parti'mandat>0 
}
foreach parti in `namelist'{
gen `parti'procpos=`parti'proc*`parti'p
}

egen sumprocpos = rowtotal (mprocpos-lok4procpos)

foreach parti in `namelist'{
gen `parti'procposk=`parti'proc/sumprocpos
gen `parti'procposksq=`parti'procposk^2
gen `parti'procposkcub=`parti'procposk^3
gen `parti'procposk4p=`parti'procposk^4

}
foreach parti in  `namelist'{
	gen `parti'procposkamalt=`parti'procposk*antalmandat
}
foreach parti in  `namelist'{
	gen `parti'procposksqamalt=`parti'procposksq*antalmandat
}
foreach parti in  `namelist'{
	gen `parti'procposcubkamalt=`parti'procposkcub*antalmandat
}
foreach parti in  `namelist'{
	gen `parti'procposk4pamalt=`parti'procposk4p*antalmandat
	drop `parti'p `parti'procpos  `parti'procposk-`parti'procposk4p
}
end
/****************************************************/

/****************************************************/
/****************************************************/
/* define program to define dummies for random interval  weighted vote share*/
/****************************************************/
capture program drop dum_diffpartyclosewvs
program dum_diffpartyclosewvs
args namelist namelist2
foreach lim in `namelist2'{   
   foreach parti in `namelist'{
 	gen `parti'closewvs`lim'= `parti'diffcal`lim'!=0
	}
}
end
/****************************************************/

/****************************************************/
/* define program to define dummies for random interval  weighted vote share*/
/****************************************************/
capture program drop dum_diffpartyclosevs
program dum_diffpartyclosevs
args namelist namelist2
foreach lim in `namelist2' {   
   foreach parti in `namelist'{
 	gen `parti'closevs`lim'= abs(`parti'diffcalvs`lim')
	}
}

end
/****************************************************/


**************************************Main Program************************************************/

clear matrix
*cd "G:\simulering\"
set more off
clear
set mem 400M
use "electoral data basic.dta"

********program to mesure distance to threshold***********
calmindiff "m c fp kd mp s v sd nd lok1 lok2 lok3 lok4""m c fp kd mp s v sd nd lok1 lok2 lok3 lok4" "

drop sprm- mindifflok4n3

foreach parti in m c fp kd mp s v nd sd lok1 lok2 lok3 lok4{
replace mindiff`parti'n=. if  mindiff`parti'n <0
replace mindiff`parti'p=. if  mindiff`parti'p <0
}
foreach var of varlist mindiffmp- mindifflok4n{
gen `var'proc=`var'
}
*generate seatshare variable***************** mandat= seat proc=share

foreach var of varlist lok2mandat lok1mandat sdmandat lok3mandat lok4mandat smandat- vmandat ndmandat{
gen `var'proc=`var'/antalmandat
}


*******generate tertment and control variables******
dum_mindiffpartycalvs "m c fp kd mp s v sd nd lok1 lok2 lok3 lok4"  "05  025 01""


****generate vote share controls**********

genws "m c fp kd mp s v sd nd lok1 lok2 lok3 lok4"  ""

*collapse data to municipal level***********
collapse  (max)ar valkrets (rawsum)rostam mrostam- lok4rostam antalroster lok1rost lok2rost lok3rost lok4rost srost mrost crost fprost kdrost mprost vrost sdrost ndrost lok2mandat sdmandat lok1mandat lok3mandat lok4mandat smandat- vmandat ndmandat mprocam- lok4procam  mprocposkam- lok4procposk4pam    mdiffcalvs05- lok4diffcalvs01 antalmandat (min) mindiffmp- mindifflok4n , by(kommun kommunkod valar)

****dum_diffpartyclosewvs "m c fp kd mp s v sd nd lok1 lok2 lok3 lok4"  "3 25 2 15 1 05""
dum_diffpartyclosevs "m c fp kd mp s v sd nd lok1 lok2 lok3 lok4"  "05  025 01""

**caluclulate vote share, seat share,weighted vote shares etc..........
foreach var of varlist lok2mandat lok1mandat sdmandat lok3mandat lok4mandat smandat- vmandat ndmandat  mprocam-  lok4procposk4pam    mdiffcalvs05- lok4diffcalvs01 mclosevs05-lok4closevs01   {
gen `var'proc =`var'/antalmandat
}

foreach parti in m c fp kd mp s v nd sd lok1 lok2 lok3 lok4{
gen `parti'proc=`parti'rost/antalroster
}

foreach parti in m c fp kd mp s v nd sd lok1 lok2 lok3 lok4{
gen `parti'amproc=`parti'rostam/rostam
}

merge kommunkod valar using "electoral data block.dta"  "municipal data.dta", sort unique
gen placinvcap1= (placinv+1)/befolkning

foreach var of varlist Kommunalskatt befolkning  placinvcap1 { 
gen ln`var'=ln(`var')
}





foreach lim in 05  025 01 {

gen lokdiffcalvs`lim'proc=lok1diffcalvs`lim'proc+lok2diffcalvs`lim'proc+lok3diffcalvs`lim'proc+lok4diffcalvs`lim'proc
gen lokclosevs`lim'=(lok1closevs`lim'+lok2closevs`lim'+lok3closevs`lim'+lok4closevs`lim')
gen lokclosevs`lim'proc=(lok1closevs`lim'proc+lok2closevs`lim'proc+lok3closevs`lim'proc+lok4closevs`lim'proc)
}
gen lokprocposkamproc=lok1procposkamproc+lok2procposkamproc+lok3procposkamproc+lok4procposkamproc
gen lokmandatproc=lok1mandatproc+lok2mandatproc+lok3mandatproc+lok4mandatproc
gen lokprocposksqamproc= lokprocposkamproc^2
gen lokprocposcubkamproc=lokprocposkamproc^3
gen lokprocposk4pamproc=lokprocposkamproc^4
gen lokproc=lok1proc+lok2proc+lok3proc+lok4proc

xtset kommunkod

sort  kommunkod valar



keep if valar<2005



gen hoger_ss=mmandatproc+cmandatproc+fpmandatproc+kdmandatproc
gen hoger_maj=hoger_ss>.5
gen vanster_ss=smandatproc+vmandatproc
gen vanster_maj=vanster_ss>.5

gen mp_vag=hoger_maj==0 & vanster_maj==0 &((hoger_ss+mpmandatproc)>.5|(vanster_ss+mpmandatproc)>.5)
gen nd_vag=hoger_maj==0 & vanster_maj==0 &((hoger_ss+ndmandatproc)>.5|(vanster_ss+ndmandatproc)>.5)
gen odef_maj=hoger_maj==0 & vanster_maj==0 & nd_vag==0 & mp_vag==0

*****generate outcome variables with correct names*********
gen tax = Kommunalskatt 
gen env= miljorankingproc
gen imm = lnplacinvcap1 

	foreach lim in 05 025 01{
	foreach parti in m c fp kd {
	
	 gen `parti'_`lim'vag= (hoger_maj==1 & (hoger_ss-	2*`parti'diffcalvs`lim'proc<50))|(hoger_maj==0 & (hoger_ss-	2*`parti'diffcalvs`lim'proc>50))
	
	}
	}
		foreach lim in 05 025 01{
	foreach parti in v s {
	
	gen `parti'_`lim'vag= (vanster_maj==1 & (vanster_ss-	2*`parti'diffcalvs`lim'proc<50))|(vanster_maj==0 & (vanster_ss-	2*`parti'diffcalvs`lim'proc>50))
	
	}
	}
		
    foreach lim in 05 025 01{
	foreach parti in mp nd lok lok1 lok2 lok3 lok4 {
	
 gen	`parti'_`lim'vag= (vanster_maj==1 & (vanster_ss+2*`parti'diffcalvs`lim'proc<50))|(vanster_maj==0 & (vanster_ss+	2*`parti'diffcalvs`lim'proc>50)) |(hoger_maj==1 & (hoger_ss+2*`parti'diffcalvs`lim'proc<50))|(hoger_maj==0 & (hoger_ss+	2*`parti'diffcalvs`lim'proc>50))
	
	}
	}

gen befolkning3=befolkning/1000
gen pop14andelpe=pop0_14andel*100
******generated lagged outcome variables***********
sort kommunkod valar
foreach var of varlist lnplacinvcap1 Kommunalskatt  miljorankingproc  {
	gen `var'_lag = `var'[_n-1] if kommunkod==kommunkod[_n-1]
	
 }
 
******generate index variables**********
gen m_pos_tax=8.8					
gen m_imp_tax=0.45
gen m_pos_env=3.7
gen m_imp_env=0.01
gen m_pos_imm=4.6
gen m_imp_imm=0.00


gen c_pos_tax=5.8				
gen c_imp_tax=0.10
gen c_pos_env=6.9	
gen c_imp_env=0.42
gen c_pos_imm=5.4
gen c_imp_imm=0.01

					
gen fp_pos_tax=5.8
gen fp_imp_tax=0.15
gen fp_pos_env=5.0
gen fp_imp_env=0.01
gen fp_pos_imm=7.1
gen fp_imp_imm=0.05

gen kd_pos_tax=6.5				
gen kd_imp_tax=0.02
gen kd_pos_env=4.8	
gen kd_imp_env=0.01
gen kd_pos_imm=6.5
gen kd_imp_imm=0.01

gen mp_pos_tax=3.9				
gen mp_imp_tax=0.03
gen mp_pos_env=8.9
gen mp_imp_env=0.79
gen mp_pos_imm=5.8	
gen mp_imp_imm=0.00

gen s_pos_tax=3.2					
gen s_imp_tax=0.10
gen s_pos_env=5.1
gen s_imp_env=0.01
gen s_pos_imm=6.6
gen s_imp_imm=0.02

			
gen v_pos_tax=1.3
gen v_imp_tax=0.33	
gen v_pos_env=5.5	
gen v_imp_env=0.07
gen v_pos_imm=6.1
gen v_imp_imm=0.01

gen nd_pos_tax=7.5				
gen nd_imp_tax=0.19
gen nd_pos_env=3.0
gen nd_imp_env=0.01
gen nd_pos_imm=1.3	
gen nd_imp_imm=0.46

foreach policy in tax env imm{
egen mean_pos_`policy'= rowmean( m_pos_`policy' c_pos_`policy'  fp_pos_`policy'  kd_pos_`policy'  mp_pos_`policy'  s_pos_`policy'  v_pos_`policy'  nd_pos_`policy')
}

foreach policy in tax env imm{
	foreach lim in 05 025 01{
		gen `policy'_t_pos_`lim'=0
		gen `policy'_c_pos_`lim'=0
		gen `policy'_t_pos_imp_`lim'=0
		gen `policy'_c_pos_imp_`lim'=0
	
		foreach parti in m c fp kd mp s v nd{

			replace  `policy'_t_pos_`lim'= `policy'_t_pos_`lim'+`parti'diffcalvs`lim'proc* (`parti'_pos_`policy'-mean_pos_`policy')
			replace  `policy'_c_pos_`lim'= `policy'_c_pos_`lim'+`parti'closevs`lim'proc* (`parti'_pos_`policy'-mean_pos_`policy')
			replace  `policy'_t_pos_imp_`lim'= `policy'_t_pos_imp_`lim'+`parti'diffcalvs`lim'proc* (`parti'_pos_`policy'-mean_pos_`policy')* `parti'_imp_`policy'
			replace  `policy'_c_pos_imp_`lim'= `policy'_c_pos_imp_`lim'+`parti'closevs`lim'proc* (`parti'_pos_`policy'-mean_pos_`policy')* `parti'_imp_`policy'
		}
	}
}

foreach policy in tax env imm{
	gen `policy'_m_pos=0
	gen `policy'_m_pos_imp=0
	gen `policy'_m_pos_imp2=0
	
	foreach parti in m c fp kd mp s v nd{

		replace  `policy'_m_pos= `policy'_m_pos+`parti'mandatproc* (`parti'_pos_`policy'-mean_pos_`policy')
		replace  `policy'_m_pos_imp= `policy'_m_pos_imp+`parti'mandatproc*  (`parti'_pos_`policy'-mean_pos_`policy')* `parti'_imp_`policy'
		replace  `policy'_m_pos_imp2= `policy'_m_pos_imp2+`parti'mandatproc*  (`parti'_pos_`policy'-mean_pos_`policy')* (`parti'_imp_`policy'+.1)
	
	}
	egen `policy'_zm_pos=std(`policy'_m_pos)
	egen `policy'_zm_pos_imp=std(`policy'_m_pos_imp)
	egen `policy'_zm_pos_imp2=std(`policy'_m_pos_imp2)
}


foreach policy in tax env imm{
	gen `policy'_vs_pos=0
	gen `policy'_vs_pos_imp=0
	gen `policy'_vs_pos_imp2=0
	
	foreach parti in m c fp kd mp s v nd{

		replace  `policy'_vs_pos= `policy'_m_pos+`parti'procposkamproc* (`parti'_pos_`policy'-mean_pos_`policy')
		replace  `policy'_vs_pos_imp= `policy'_m_pos_imp+`parti'procposkamproc*  (`parti'_pos_`policy'-mean_pos_`policy')* `parti'_imp_`policy'
		replace  `policy'_vs_pos_imp2= `policy'_m_pos_imp2+`parti'procposkamproc*  (`parti'_pos_`policy'-mean_pos_`policy')* (`parti'_imp_`policy'+.1)
	
	}
	
	egen `policy'_zvs_pos=std(`policy'_vs_pos)
	egen `policy'_zvs_pos_imp=std(`policy'_vs_pos_imp)

}


*****Table 3**********

reg lnplacinvcap1  antalmandat, cluster (kommunkod) robust   
	outreg2   antalmandat using ind_A5_bal, ctitle ("Model start") excel nocons dec(4) se replace
foreach var of varlist lnplacinvcap1 Kommunalskatt  miljorankingproc { 

set more off
	

foreach lim in 025{
	
	xtreg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc  using ind_A5_bal, ctitle ("Model `var'`lim'fec") excel nocons dec(2) se append


}
}

reg realinkomst  antalmandat, cluster (kommunkod) robust   
	outreg2  antalmandat using index_A5_bal, ctitle ("Model start")  nocons dec(4)excel se replace
foreach var of varlist tax env imm  { 


	


foreach lim in 025 {

	xtivreg2 `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   ffirst 
	outreg2 `var'_zm_pos_imp   using index_A5_bal, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append

	xtivreg2 `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   ffirst
	outreg2 `var'_zm_pos using index_A5_bal, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append

}


}

