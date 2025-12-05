
/**************************************Program definitions************************************************/
/**************************************Program definitions************************************************/

/****************************************************/
/* define program to define minimal distance to seat change          */
/****************************************************/

/* Notes Simon, Lukas 5. January 2023
Alternative 1: nur zusätzliche Stimmen analog zu unserer Lösung

Alternative 2: minimale Zahl der Stimmen für einen Sitzverlust über alle Parteien hinweg plus eine Nebenbedingung: zusätzliche Stimmen für Fokuspartei, um auf die anfänglich knappste Partei für einen Sitzgewinn aufzuholen (Partei 3 in Teilfall d)

Problem: Die von der verlierenden Partei abgezogenen Stimmen werden nicht der Fokuspartei zugerechnet. Somit ist die Nebenbedingung schwerer erfüllbar. 

Folgerung: keine Umverteilung, nur eine Nebenbedingung beachtet

Alternative 3: Berechnet die Summe (i) der Stimmen, die die Partei mit dem knappsten Sitz verlieren müsste, sodass ihre Vergleichszahl mit der Vergleichszahl der Fokuspartei übereinstimmt und (ii) das Total der Stimmen, die diejenigen Parteien verlieren müssten, die näher an einem Sitzgewinn sind als die Fokuspartei, sodass deren Vergleichszahl gleich gross ist wie jene der Fokuspartei. In der Grafik A.2 entspricht (i) der Bewegung von C_2 zu C_1 und (ii) der Bewegung von C_3 zu C_1. 
F
olgerung: Hier scheint Folke alle Nebenbedingung zu beachten, jedoch findet nach wie vor keine Umverteilung statt
*/

*use "C:/Schmidlu/Dropbox/Projekt Nationalräte/06_Resources/02_Lit_RDD_PR/Folke (2014), suppl material/sigtuna1991.dta", clear
insheet using "C:/Schmidlu/Dropbox/Projekt Nationalräte/06_Resources/02_Lit_RDD_PR/Folke (2014), suppl material/Example Table A1 extended.csv", clear delimiter(";")

*capture program drop calmindiff // SLLS 5 Jan 2023
*program calmindiff // SLLS 5 Jan 2023
* args namelist namelist2 // SLLS 5 Jan 2023

* keep variables from appendix
gen spr=voteshares
keep parti spr mandat
rename mandat m

gen temp=1
reshape wide spr m, i(temp) j(parti) string
  


 /* SLLS: not necessary in our example
  foreach parti in `namelist'{
    gen spr`parti' =`parti'proc
	gen m`parti' = `parti'mandat
  }
  */
  
  local namelist "m c fp kd mp s v nd sd" // SLLS 5 Jan 2023 

   foreach parti in `namelist'{
   	gen p`parti'=m`parti'>0 & m`parti'!=. // Dummy variable for parties with at least one seat
		gen p1`parti'=m`parti'==1 // Dummy variable for parties with exactly one seat

		}

**define comparison numbers*******
   foreach parti in `namelist'{
		gen j`parti'=spr`parti'*(1-p`parti')/1.4 +spr`parti'*p`parti'/(1+2*m`parti')
		}
// Sainte-Laguë ratio for next seat
		
		
   foreach parti in `namelist'{
		gen j2`parti'=spr`parti'*(p1`parti')/1.4 +spr`parti'*(1-p1`parti')/(1+2*(m`parti'-1)) if m`parti'>0
		}
// Sainte-Laguë ratio for previous seat


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

// Comment do-file Fiva et al. (2018): "( a og c i olle-appendix)"
// Comment LSSL: In this alternative, Olle only adds vote to the focus party 

***define distance to seat gain in own votes*****
	foreach parti in `namelist'{
		gen `parti'diffjminj2= (minj2 -j`parti') * ( (1-p`parti') * 1.4 + p`parti'*(1+2*m`parti') )
		replace `parti'diffjminj2=.  if `parti'diffjminj2 <0
	}
	
* Result: diffjminj2 is almost equivalent to the row "diff gain ind %" in Table 
* 		  of Folke's Online Appendix.

***define distance to seat loss in own votes*****	
	foreach parti in `namelist'{
		gen `parti'diffj2maxj= (j2`parti'-maxj1)* (p1`parti' *1.4 +(1-p1`parti')*(1+2*(m`parti'-1))) if m`parti'>0 & m`parti'!=.
		replace `parti'diffj2maxj=. if `parti'diffj2maxj <0
	}
	
* Sainte-Laguë number of focus party's currenct seat minus highest Saint-Laguë numbers  
* among those Saint-Laguë numbers without seat. This is then retransformed into actual
* votes. 

	foreach parti in `namelist'{
		gen mindiff`parti'p1= `parti'diffjminj2 // Mininum distance to gain a seat
		gen mindiff`parti'n1=`parti'diffj2maxj   if m`parti'>0 & m`parti'!=.  // Mininum distance to lose a seat
	}

*****************ALTERNATIVE 2***************************************************

* Comment do-file Fiva et al. (2018): "( tilfelle b og d)"
* Comment LSSL: In this alternative, Olle seems to only subtracts votes from one
* or two parties. However, it is not clear that he really substracts votes from
* the first party.

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
	
* The number of votes the party that is closest to losing a seat can lose until they lose the seat. 
* In Appendix Figure A.2 case d the movement of C_2 down to C_3. 

*** Define distance for a party gaining a seat through gaining votes such that j=maxj1 and another party loosing votes****
	
	foreach parti in `namelist'{
		gen mindiff`parti'p2= mindiffj2maxj  + (maxj1 -j`parti')*((1-p`parti')*1.4 +p`parti'*(1+2*m`parti'))
		replace mindiff`parti'p2=. if maxj1 < j`parti'
	}
	
* Second part of mindiff`parti'p2-equation is the number of additional votes
* that the focus party needs to catch up to the comparison number of the party 
* closest to winning an additional seat. In Appendix Figure A.2 case d, it is 
* the movement of C_1 to C_3. 

* Replicate Folke's margin for the socialist party (s)

br mindiffj2maxj *diffj2maxj
* Result: The party that is closest to losing a seat is the local party (sd). If 
* 		  it lost .0007 (0.001 if rescaled vote shares are used) of their vote 
*		  share, they would lose the seat. 
	
br maxj1 j*
* Result: The environmenal party (mp) is the party with the highest Saint-Laguë 
*		  numbers for the next seat. It amounts to .0101333 (0.01 if rescaled 
*		  vote shares are used)

foreach parti in `namelist'{
		gen mindiff`parti'p2_2ndpart= (maxj1 -j`parti')*((1-p`parti')*1.4 +p`parti'*(1+2*m`parti'))
	}
br mindiffsp2_2ndpart

* final solution times 100
foreach parti in `namelist'{
		gen `parti'_solution= (mindiffj2maxj+(maxj1 -j`parti')*((1-p`parti')*1.4 +p`parti'*(1+2*m`parti')))
	}
br *_solution*

br mindiffj2maxj mindifffpp2_2ndpart


* socialist party (s)

* Result: .0085333. (.008 if rescaled vote shares are used)

* End result: 
di .0007 + .0085333 
* .0092333
di .001 + .008
* .009 for rescaled version

* liberal party (fp)

br mindifffpp2_2ndpart

di 0.001+.015

* Folke's solution to check: If the liberal party wins a vote share of 0.015
* and the local party loses a vote share of 0.001, they have the same Saint-Laguë
* ratio as the environmental party and therefore the liberal party has a chance
* to win a seat from the local party. We can clearly see that the votes are not
* redistributed. First, the losses and the gains are not equal. Second, the votes
* share of 0.001 lost by the local party is not added to the vote share of the
* liberal party. 



*** Define distance for a party loosing a seat through loosing votes such that j2=minj2 and another party gaining votes****
	foreach parti in `namelist'{
		gen mindiff`parti'n2=mindiffjminj2+ (j2`parti'-minj2 ) * (p1`parti'*1.4 + (1-p1`parti') * (1+2*(m`parti'-1) ) )  if m`parti'>0 & m`parti'!=.
		replace mindiff`parti'n2=. if j2`parti'< minj2
	}

*****************ALTERNATIVE 3***************************************************
		
* Comment do-file Fiva et al. (2018): "(tilfelle e)"
* Comment LSSL: 
	

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

* Note SLLS: These are all Sainte-Laguë numbers that are between focus party and
*			 target party that have to be broght down such that focus party
*			 effectively gains seat (C_3(s3) to C_1(s1) in case e in Folke Appendix)
	
	
*** Define how many votes  party2 would need to loose such that party2j2<party1j1****

	foreach parti in `namelist'{
		foreach parti2 in `namelist2'{ 
			gen `parti'`parti2'diffj2j= (j2`parti2' -j`parti') * (p1`parti2'*1.4 +(1-p1`parti2')*(1+2*(m`parti2'-1))) if m`parti2'>0 &  m`parti2'!=.
			replace `parti'`parti2'diffj2j =. if j`parti'==j`parti2' |`parti'`parti2'diffj2j<0
		}
	}

* Note SLLS: These are the Sainte-Laguë of the target party that have to decrease
*			 such that our focus party wins a seat (C_2(s3) to C_1(s2-1) in case e in Folke Appendix). 
*			 These are still pairs of parties parti-parti2. 



*** Define vote minimal loss for any party2  such that party2j2<party1j1****
	foreach parti in `namelist'{	
	gen `parti'mindiffj2j=100000000000000
		foreach parti2 in `namelist2'{ 
			replace `parti'mindiffj2j= `parti'`parti2'diffj2j if `parti'`parti2'diffj2j<`parti'mindiffj2j & `parti'`parti2'diffj2j!=.
		}	
	}
	
* Note SLLS: Find out which of the target parties has the lowest distance for 
*			 target party.  

	
**** Define combination of minimal loss for any party2  such that party2j2<party1j1 and ****
	**** total vote losses for other parties such that  party1j<j1max to define disance to seat gain*****
	
foreach parti in `namelist'{
	gen mindiff`parti'p3=  `parti'sumjdiff+`parti'mindiffj2j
	}
 
* Note SLLS: Add votes of target party and votes of all parties such that constraints are met. 
*			 Bring C_2(s2-1) and C_3(s3) down in Folke's example. 

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


local parti "s"
br mindiff`parti'p mindiff`parti'p1 mindiff`parti'p2 mindiff`parti'p3

local parti "sd"
br mindiff`parti'p mindiff`parti'p1 mindiff`parti'p2 mindiff`parti'p3



/****************************************************/
/****************************************************/
/* define program to dummies from minimum distance to change  calculated measured in vote share         */
/****************************************************/

capture program drop dum_mindiffpartycalvs
program dum_mindiffpartycalvs
args namelist namelist2

foreach lim in `namelist2' {   
   foreach parti in `namelist'{
			gen `parti'diffpcalvs`lim'=  (mindiff`parti'p)<0.0`lim'

			gen `parti'diffncalvs`lim'=  (mindiff`parti'n)<0.0`lim'



			gen `parti'diffcalvs`lim' =(`parti'diffncalvs`lim'-`parti'diffpcalvs`lim')*.5
			drop `parti'diffncalvs`lim'  `parti'diffpcalvs`lim' 


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
***********************************************/

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


/**************************************Main Program************************************************/

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




**********generate variables for table 10****************

foreach lim in  05  025 01{
foreach parti in m c fp kd mp  v nd  lok s{
gen `parti'1sclose`lim'= (`parti'diffcalvs`lim'proc<0 &`parti'mandatproc==0)|(`parti'mandatproc == `parti'diffcalvs`lim'proc*2   & `parti'diffcalvs`lim'proc > 0)
}
}


foreach lim in  05  025 01{
foreach parti in m c fp kd mp  v nd  lok s{
gen `parti'1sclose`lim'pos= `parti'1sclose`lim' & `parti'mandatproc>0
}
}

foreach lim in  05  025 01{
foreach parti in m c fp kd mp  v nd lok s {

gen `parti'diffcalvs`lim'proc1s= `parti'diffcalvs`lim'proc*`parti'1sclose`lim'
}
}


foreach lim in  05  025 01 {
foreach parti in m c fp kd mp  v nd lok  s {

gen `parti'closevs`lim'proc1s= `parti'closevs`lim'proc*`parti'1sclose`lim'
}
}




foreach lim in  05  025  01{
foreach parti in m c fp kd mp  v nd lok  s {

gen `parti'dumclose`lim'= `parti'closevs`lim'proc>0
gen `parti'dumclose1s`lim'= `parti'closevs`lim'proc1s>0

}
}
foreach lim in  05  025  01{
foreach parti in m c fp kd mp  v nd lok  s {

gen `parti'share`lim'= `parti'dumclose1s`lim'/`parti'dumclose`lim'
}
}


foreach lim in  05  025 01{
foreach parti in m c fp kd mp  v nd lok   {

gen `parti'procposkamproc`lim'1s= `parti'procposkamproc*`parti'1sclose`lim'
}
foreach parti in m c fp kd mp  v nd  {
gen `parti'procposksqamproc`lim'1s= `parti'procposksqamproc*`parti'1sclose`lim'
gen `parti'procposkcubamproc`lim'1s= `parti'procposcubkamproc*`parti'1sclose`lim'
gen `parti'procposk4pamproc`lim'1s= `parti'procposk4pamproc*`parti'1sclose`lim'
}
foreach parti in lok1 lok2 lok3 lok4  {
gen `parti'procposksqamproc`lim'1s= `parti'procposksqamproc*lok1sclose`lim'
gen `parti'procposkcubamproc`lim'1s= `parti'procposcubkamproc*lok1sclose`lim'
gen `parti'procposk4pamproc`lim'1s= `parti'procposk4pamproc*lok1sclose`lim'
}

}



****************************
*******Figure 2

gen lokmandat= lokmandatproc*antalmandat

hist mmandat, frac xtitle(Seats Conservative Party) ytitle() width(1)   
graph save m, replace 

hist cmandat, frac xtitle(Seats Center Party)ytitle() width(1)
graph save c, replace

hist fpmandat, frac xtitle(Seats Liberal Party)ytitle()width(1)
graph save fp, replace


hist kdmandat, frac xtitle(Seats Christian Democrats)ytitle()width(1)
graph save kd, replace


hist mpmandat, frac xtitle(Seats Environmental Party)ytitle()width(1)
graph save mp, replace

hist smandat, frac xtitle(Seats Social Democrats)ytitle()width(1)
graph save s, replace

hist vmandat, frac xtitle(Seats Left Party)ytitle()width(1)
graph save v, replace

hist lokmandat, frac xtitle(Seats Local Parties)ytitle()width(1)
graph save lok, replace

hist ndmandat if valar==1991 |valar==1994, frac xtitle(Seats New Democracy )ytitle()width(1)
graph save nd, replace

gr combine m.gph c.gph fp.gph kd.gph mp.gph s.gph v.gph nd.gph lok.gph
*********figure A 1*******
histogram mmandat if mclosevs025proc>0 & mclosevs025proc!=.,discrete width(.5) xtitle(Seats Conservative Party)
graph save m_mandat, replace

histogram cmandat if cclosevs025proc>0 & cclosevs025proc!=.,discrete width(.5) xtitle(Seats Center Party)
graph save c_mandat, replace
histogram fpmandat if fpclosevs025proc>0 & fpclosevs025proc!=.,discrete width(.5) xtitle(Seats Liberal Party)
graph save fp_mandat, replace
histogram kdmandat if kdclosevs025proc>0 & kdclosevs025proc!=.,discrete width(.5) xtitle(Seats Christian Democrats)
graph save kd_mandat, replace
histogram mpmandat if mpclosevs025proc>0 & mpclosevs025proc!=.,discrete width(.5) xtitle(Seats Environmental Party)
graph save mp_mandat, replace
histogram smandat if sclosevs025proc>0 & sclosevs025proc!=.,discrete width(.5) xtitle(Seats Social Democrats)
graph save s_mandat, replace
histogram vmandat if vclosevs025proc>0 & vclosevs025proc!=.,discrete width(.5) xtitle(Seats Left Party)
graph save v_mandat, replace
histogram ndmandat if ndclosevs025proc>0 & ndclosevs025proc!=.,discrete width(.5) xtitle(Seats New Democracy)
graph save nd_mandat, replace
gr combine m_mandat.gph c_mandat.gph fp_mandat.gph kd_mandat.gph mp_mandat.gph s_mandat.gph v_mandat.gph nd_mandat.gph
graph export seats_close.tif, replace

**********Table 2 regression"""""""""""""""""""""""""""

set more off


	reg mmandatproc  befolkning, cluster (kommunkod) robust   
	outreg2   befolkning using 1step, ctitle ("Model start") excel nocons dec(2) se replace
foreach parti in m c fp kd mp s v nd lok{ 
foreach lim in 05 025 01{

	reg `parti'mandatproc  `parti'diffcalvs`lim'proc `parti'closevs`lim'proc, cluster (kommunkod) robust   
	outreg2   `parti'diffcalvs`lim'proc  using 1step, ctitle ("Model `parti'`lim'") excel nocons dec(2) se append

}
foreach lim in 05 025 01{

	reg `parti'mandatproc  `parti'diffcalvs`lim'proc `parti'procposkamproc `parti'procposksqamproc `parti'procposcubkamproc `parti'procposk4pamproc `parti'closevs`lim'proc, cluster (kommunkod) robust   
	outreg2   `parti'diffcalvs`lim'proc  using 1step, ctitle ("Model p`parti'`lim'") excel nocons dec(2) se append
}
}


*****Table 3**********

reg lnplacinvcap1  antalmandat, cluster (kommunkod) robust   
	outreg2   antalmandat using tabellolsbas, ctitle ("Model start") excel nocons dec(4) se replace
foreach var of varlist lnplacinvcap1 Kommunalskatt  miljorankingproc { 

set more off
	

foreach lim in 025{
	
	xtreg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc  using tabellolsbas, ctitle ("Model `var'`lim'fec") excel nocons dec(2) se append

	xtivreg2 `var' (mmandatproc- ndmandatproc lokmandatproc =mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2    mmandatproc- ndmandatproc lokmandatproc  using tabellolsbas, ctitle ("Model `var'`lim'fec") excel nocons dec(2) se append

	xtreg `var'     mmandatproc- ndmandatproc lokmandatproc  dv1982- dv2002, fe  robust cluster (kommunkod)
	outreg2  mmandatproc- ndmandatproc lokmandatproc  using tabellolsbas, ctitle ("Model `var'`lim'fec") excel nocons dec(2) se append


}
}

reg realinkomst  antalmandat, cluster (kommunkod) robust   
	outreg2  antalmandat using index_base, ctitle ("Model start")  nocons dec(4)excel se replace
foreach var of varlist tax env imm  { 


	


foreach lim in 025 {

	xtivreg2 `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   ffirst 
	outreg2 `var'_zm_pos_imp   using index_base, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append

	xtivreg2 `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   ffirst
	outreg2 `var'_zm_pos using index_base, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append

}


}



****Tables A2 A3 and A4**********
reg realinkomst  antalmandat, cluster (kommunkod) robust   
	outreg2  antalmandat using sensitivity, ctitle ("Model start")  nocons dec(4)excel se replace
	
foreach var in lnplacinvcap1 Kommunalskatt  miljorankingproc  { 

*set more off
	

foreach lim in 05{
	reg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002, cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc   using sensitivity, ctitle ("Model `var'`lim'vst")  nocons dec(2)excel se append

	xtreg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc   using sensitivity, ctitle ("Model `var'`lim'vstm")  nocons dec(2)excel se append




}
foreach lim in 025 {
	reg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc , cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc   using sensitivity, ctitle ("Model `var'`lim'vst")  nocons dec(2)excel se append


	reg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc, cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc   using sensitivity, ctitle ("Model `var'`lim'vst")  nocons dec(2)excel se append

	reg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  , cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc   using sensitivity, ctitle ("Model `var'`lim'vst")  nocons dec(2)excel se append

	reg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002, cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc   using sensitivity, ctitle ("Model `var'`lim'vst")  nocons dec(2)excel se append

	xtreg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc   using sensitivity, ctitle ("Model `var'`lim'vstm")  nocons dec(2)excel se append

	xtreg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002 pop0_5andel pop5_14andel pop15_24andel pop35_44andel pop45_54andel pop65_74andel befolkning hogskolakortpe gymnasiepe hogskolalangpe  , fe cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc  using sensitivity, ctitle ("Model `var'`lim'vstmc")  nocons dec(2)excel se append



}
foreach lim in  01{
	reg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002, cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc   using sensitivity, ctitle ("Model `var'`lim'vst")  nocons dec(2)excel se append

	xtreg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc   using sensitivity, ctitle ("Model `var'`lim'vstm")  nocons dec(2)excel se append



}
}

reg realinkomst  antalmandat, cluster (kommunkod) robust   
	outreg2  antalmandat using sens_index, ctitle ("Model start")  nocons dec(4)excel se replace
foreach var of varlist tax env imm  { 

set more off
	
foreach lim in 05{
	ivreg `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	xtivreg2 `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp   using sens_index, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append
}
foreach lim in 025 {
	ivreg `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) , cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	ivreg `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc, cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	ivreg `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc , cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	ivreg `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	xtivreg2 `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp   using sens_index, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append
	xtivreg2 `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002 pop0_5andel pop5_14andel pop15_24andel pop35_44andel pop45_54andel pop65_74andel befolkning hogskolakortpe gymnasiepe hogskolalangpe  , fe cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp  using sens_index, ctitle ("Model `var'`lim'vstmc")  nocons dec(3)excel se append
}
foreach lim in  01{
	ivreg `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	xtivreg2 `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp   using sens_index, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append
}
foreach lim in 05{
	ivreg `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	xtivreg2 `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos   using sens_index, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append
}
foreach lim in 025 {
	ivreg `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) , cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	ivreg `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc, cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	ivreg `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc , cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	ivreg `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	xtivreg2 `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos   using sens_index, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append
	xtivreg2 `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002 pop0_5andel pop5_14andel pop15_24andel pop35_44andel pop45_54andel pop65_74andel befolkning hogskolakortpe gymnasiepe hogskolalangpe  , fe cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos  using sens_index, ctitle ("Model `var'`lim'vstmc")  nocons dec(3)excel se append

}
foreach lim in  01{
	ivreg `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos   using sens_index, ctitle ("Model `var'`lim'vst")  nocons dec(3)excel se append
	xtivreg2 `var' (`var'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos   using sens_index, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append



}
}


******Table A6****************************

reg realinkomst  antalmandat, cluster (kommunkod) robust   
	outreg2   antalmandat using tabell_match, ctitle ("Model start") excel nocons dec(4) se replace

foreach parti1 in  m c fp kd mp s v  nd {
foreach var of varlist  miljorankingproc  lnplacinvcap1 Kommunalskatt {
foreach lim in  025 { 
gen `parti1'diffcalvs`lim'pd= `parti1'diffcalvs`lim'proc>0

   	reg `var' `parti1'diffcalvs`lim'proc  if `parti1'closevs`lim'proc!=0 ,  cluster(kommunkod)  robust 
   outreg2  `parti1'diffcalvs`lim'proc using tabell_match, ctitle ("`lim'") excel nocons dec(2) se append
   	reg `var' `parti1'diffcalvs`lim'proc  `parti1'procposkamproc `parti1'procposksqamproc `parti1'procposcubkamproc `parti1'procposk4pamproc  i.valar if `parti1'closevs`lim'proc!=0 ,  cluster(kommunkod)  robust 
   outreg2  `parti1'diffcalvs`lim'proc i.valar using tabell_match, ctitle ("`lim'") excel nocons dec(2) se append 
   
     reg `var' `parti1'diffcalvs`lim'pd i.valar `parti1'procposkamproc `parti1'procposksqamproc `parti1'procposcubkamproc `parti1'procposk4pamproc  if `parti1'closevs`lim'proc!=0 ,  cluster(kommunkod)  robust 
   outreg2  `parti1'diffcalvs`lim'pd  using tabell_match, ctitle ("`lim'") excel nocons dec(3) se append	
   drop `parti1'diffcalvs`lim'pd
   }
}
}

*****table A7********




reg realinkomst  antalmandat, cluster (kommunkod) robust   
	outreg2   antalmandat using robust, ctitle ("Model start") excel nocons dec(4) se replace
foreach var of varlist realinkomst hogskolalang lnbefolkning  lnplacinvcap1_lag Kommunalskatt_lag  miljorankingproc_lag    { 
	

foreach lim in 025   {
	reg `var'  mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002  ,  cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc  using robust_lag, ctitle ("Model `var'`lim'fec") excel nocons dec(1) se append

	xtreg `var'  mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002 , fe cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc  using robust_lag, ctitle ("Model `var'`lim'fec") excel nocons dec(1) se append


}
}



*table A8******************************

************robust index*************
set more off
reg realinkomst  antalmandat, cluster (kommunkod) robust   
	outreg2  antalmandat using robust_index, ctitle ("Model start")  nocons dec(4)excel se replace
foreach lim in 025 {
	foreach ind in tax env imm  { 
		foreach var of varlist realinkomst   { 
			ivreg `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vst")  nocons dec(3)excel se append
			xtivreg2 `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vstm")  nocons dec(3)excel se append
		}
	}
}


foreach lim in 025 {
	foreach ind in tax env imm  { 
		foreach var of varlist  hogskolalangpe   { 
			ivreg `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vst")  nocons dec(3)excel se append
			xtivreg2 `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vstm")  nocons dec(3)excel se append
		}
	}
}
foreach lim in 025 {
	foreach ind in tax env imm  { 
		foreach var of varlist  pop14andelpe      { 
			ivreg `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vst")  nocons dec(3)excel se append
			xtivreg2 `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vstm")  nocons dec(3)excel se append
		}
	}
}
foreach lim in 025 {
	foreach ind in tax env imm  { 
		foreach var of varlist  lnplacinvcap1_lag   { 
			ivreg `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
		*	outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vst")  nocons dec(3)excel se append
			xtivreg2 `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
		*	outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vstm")  nocons dec(3)excel se append
		}
	}
}

foreach lim in 025 {
	foreach ind in tax env imm  { 
		foreach var of varlist Kommunalskatt_lag      { 
			ivreg `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vst")  nocons dec(3)excel se append
			xtivreg2 `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vstm")  nocons dec(3)excel se append
		}
	}
}
foreach lim in 025 {
	foreach ind in tax env imm  { 
		foreach var of varlist  miljorankingproc_lag    { 
			ivreg `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vst")  nocons dec(3)excel se append
			xtivreg2 `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vstm")  nocons dec(3)excel se append
		}
	}
}

foreach lim in 025 {
	foreach ind in tax env imm  { 
		foreach var of varlist  lnbefolkning    { 
			ivreg `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vst")  nocons dec(3)excel se append
			xtivreg2 `var' (`ind'_zm_pos = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc) mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
			outreg2 `ind'_zm_pos   using robust_index, ctitle ("Model `ind'`lim'vstm")  nocons dec(3)excel se append
		}
	}
}



******Table A9************
*********Majoritetsförhållanden
*****tabell med ols och basspec och iv**********

reg realinkomst  antalmandat, cluster (kommunkod) robust   
	outreg2   antalmandat using maj_mech, ctitle ("Model start") excel nocons dec(4) se replace
foreach var of varlist lnplacinvcap1 Kommunalskatt  miljorankingproc { 


	

foreach lim in 05 025{
	
	xtreg `var' mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc  hoger_maj mp_vag nd_vag odef_maj  lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc  dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2   mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc hoger_maj mp_vag nd_vag odef_maj using maj_mech, ctitle ("Model `var'`lim'fec") excel nocons dec(2) se append


}
}


foreach var of varlist imm  env tax { 

foreach lim in 05 025{

	xtivreg2 `var' (`var'_zm_pos_imp = mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  hoger_maj mp_vag nd_vag odef_maj mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos_imp  hoger_maj mp_vag nd_vag odef_maj  using maj_mech, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append

	xtivreg2 `var' (`var'_zm_pos= mdiffcalvs`lim'proc-mpdiffcalvs`lim'proc vdiffcalvs`lim'proc nddiffcalvs`lim'proc  lokdiffcalvs`lim'proc)  hoger_maj mp_vag nd_vag odef_maj mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc mprocposksqamproc- mpprocposksqamproc vprocposksqamproc ndprocposksqamproc- lok4procposksqamproc mprocposcubkamproc- mpprocposcubkamproc vprocposcubkamproc ndprocposcubkamproc- lok4procposcubkamproc  mprocposk4pamproc - mpprocposk4pamproc vprocposk4pamproc ndprocposk4pamproc- lok4procposk4pamproc dv1982- dv2002, fe cluster(kommunkod)  robust   
	outreg2 `var'_zm_pos  hoger_maj mp_vag nd_vag odef_maj  using maj_mech, ctitle ("Model `var'`lim'vstm")  nocons dec(3)excel se append



}
}


*******Table A10************




set more off
	reg Kommunalskatt  antalmandat, cluster (kommunkod) robust   
	outreg2   antalmandat using 1split, ctitle ("Model start") excel nocons dec(4) se replace

foreach var of varlist lnplacinvcap1  miljorankingproc Kommunalskatt  { 

foreach lim in 05  025  {


	xtreg `var'  mdiffcalvs`lim'proc mdiffcalvs`lim'proc1s cdiffcalvs`lim'proc cdiffcalvs`lim'proc1s fpdiffcalvs`lim'proc fpdiffcalvs`lim'proc1s kddiffcalvs`lim'proc kddiffcalvs`lim'proc1s mpdiffcalvs`lim'proc mpdiffcalvs`lim'proc1s vdiffcalvs`lim'proc vdiffcalvs`lim'proc1s nddiffcalvs`lim'proc nddiffcalvs`lim'proc1s lokdiffcalvs`lim'proc lokdiffcalvs`lim'proc1s  mclosevs`lim'proc- mpclosevs`lim'proc vclosevs`lim'proc ndclosevs`lim'proc    lok1closevs`lim'proc mprocposkamproc- mpprocposkamproc  vprocposkamproc ndprocposkamproc lokprocposkamproc  dv1982- dv2002 m1sclose`lim'-lok1sclose`lim' mclosevs`lim'proc1s- vclosevs`lim'proc1s lokclosevs`lim'proc1s mprocposkamproc`lim'1s-  lok4procposk4pamproc`lim'1s, fe cluster (kommunkod) robust   
	outreg2  mdiffcalvs`lim'proc mdiffcalvs`lim'proc1s cdiffcalvs`lim'proc cdiffcalvs`lim'proc1s fpdiffcalvs`lim'proc fpdiffcalvs`lim'proc1s kddiffcalvs`lim'proc kddiffcalvs`lim'proc1s mpdiffcalvs`lim'proc mpdiffcalvs`lim'proc1s vdiffcalvs`lim'proc vdiffcalvs`lim'proc1s nddiffcalvs`lim'proc nddiffcalvs`lim'proc1s lokdiffcalvs`lim'proc lokdiffcalvs`lim'proc1s using 1split, ctitle ("Model `var'`lim'fec") excel nocons dec(2) se append

}
}