clear
set mem 1g
set more off

/**************************************Program definitions************************************************/

/****************************************************/
/* define program to simulate changes in vote share          */
/****************************************************/
capture program drop simchange_mun
program simchange_mun
	args namelist  sim max min
  
	expand `sim'

	gen vch_base= (int(_n*100/`sim'))/500-0.1      /* setter intervallet hvor vi leter etter voteshare-margin (forcing) */

	gen vch=(.5-ss_coal)+vch_base                    /* sentrer rundt seatshare-margin */
	
	replace vch=`min'-int(_n/`sim'*(`min'-`max')*1000)/1000 if `max'!=.      /* zoomer inn: "loop 2 " og "loop 3" , i stedet for å kjøre på fordefinert intervall så bruker min og max definert under */

	bysort vch: egen vch_count= count (vch)
	egen vch_count_mean= mean(vch_count)
	
	drop if vch_count==.
	drop if vch_count<(vch_count_mean/1.5)
	

	*******simulate vote changes for parties in majority****************    /* majority f.eks venstrepartien , `parti'_maj dummy =1 hvis inngår i koalisjonen vi er interessert i*/ 
	foreach parti in `namelist'{											/* hvert parti får en vekt for vote change , og vekten avhenger av størrelse på parti og en slumpmessig komponent (uniform) */
		***gen `parti'_coal_ch_b =uniform()*(`parti'_vs)*`parti'_coal  /* her er versjonen med en base component, virkre ikke helt rimelig*/
		gen `parti'_coal_ch_b =uniform()*(`parti'_vs)*`parti'_coal
		replace `parti'_coal_ch_b=0 if `parti'_vs==0
	} 

	egen sum_coal_ch= rowtotal (RV_coal_ch_b-JointR2_coal_ch_b)   				/* normaliserer, summerer til 1 */

	**** HUSK: Legge inn restriksjon som sikrer at vi ikke får mindre enn 0% Votes
	
	foreach parti in `namelist'{
		gen `parti'_coal_ch =vch*`parti'_coal_ch_b/sum_coal_ch
	} 

		


	drop RV_coal_ch_b-JointR2_coal_ch_b

	***************simulate vote changes for parties in minority************************

	foreach parti in `namelist'{
		gen `parti'_min_ch_b =uniform()*(`parti'_vs)*(1-`parti'_coal)
		replace `parti'_min_ch_b=0 if `parti'_vs==0
		
	}

	egen sum_min_ch= rowtotal (RV_min_ch_b-JointR2_min_ch_b)

	foreach parti in `namelist'{
		gen `parti'_min_ch =(0-vch)*`parti'_min_ch_b/sum_min_ch
	}

	************************Distribute vote changes over electoral districts***********************


		foreach parti in `namelist'{
			gen `parti'_vs_sim=`parti'_vs+`parti'_min_ch+`parti'_coal_ch
		}
	
keep knr electionperiod maj_coal ss_coal vs_coal vch RV_vs_sim-JointR2_vs_sim RV_s-JointR2_s RV_ss-JointR2_ss RV_coal-JointR2_coal tot_s

end
/****************************************************/




/****************************************************/
/* define program to distribute mandates in whole municipality           */
/****************************************************/

capture program drop comp_mandates_mun
program comp_mandates_mun

	args namelist 
	
	gen vs_sum=0
	

		foreach parti in `namelist'{
    		gen `parti'_s_sim = `parti'_s  
    	}
			
		gen felmax=1
		while felmax>0 {
		
			drop felmax
			
			foreach parti in `namelist'{
				gen p`parti'=`parti'_s_sim>0 
				gen p1`parti'=`parti'_s_sim==1
				gen j`parti'=`parti'_vs_sim*(1-p`parti')/1.4 +`parti'_vs_sim*p`parti'/(1+2*`parti'_s_sim)
				gen j2`parti'=`parti'_vs_sim/p`parti'*(p1`parti')/1.4 +`parti'_vs_sim/p`parti'*(1-p1`parti')/(1+2*(`parti'_s_sim-1))
					}

			egen maxj1 = rowmax(jRV jSV jDNA jV jSP jKRF jH jFRP jPP jMDG jIndep1 jIndep2 jIndep3 jIndep4 jIndep5 jIndep6 jOther1 jOther2 jOther3 jOther4 jOther5 jJointL jJointR1 jJointR2)
			egen minj2 = rowmin(j2RV j2SV j2DNA j2V j2SP j2KRF j2H j2FRP j2PP j2MDG j2Indep1 j2Indep2 j2Indep3 j2Indep4 j2Indep5 j2Indep6 j2Other1 j2Other2 j2Other3 j2Other4 j2Other5 j2JointL j2JointR1 j2JointR2)

			gen fel= maxj1>minj2

			foreach parti in `namelist'{
				replace `parti'_s_sim= `parti'_s_sim+1 if j`parti'==maxj1 & fel==1
				replace `parti'_s_sim=`parti'_s_sim -1 if j2`parti'==minj2 & fel==1
			}

			egen felmax= max(fel)
			drop pRV-fel
		}
		
		
	foreach parti in `namelist'{
			replace `parti'_vs_sim= `parti'_vs_sim/vs_sum
	
		}

gen sum=1

keep knr electionperiod maj_coal vs_coal ss_coal sum RV_s_sim-JointR2_s_sim RV_vs_sim-JointR2_vs_sim vch RV_coal-JointR2_coal RV_coal-JointR2_coal tot_s 

***collapse (sum) sum by (RV_s_sim-JointR2_s_sim  kommun kommunkod valar)

end

/****************************************************/
/* define program to measure minimum distance to seatchare changes */
/****************************************************/
capture program drop dist_ch_coal
program dist_ch_coal 

	args coallist diff
	
	gen ss_sim_coal=0
	gen vs_sim_coal=0
	
	foreach parti in `coallist'{
		replace ss_sim_coal=ss_sim_coal+`parti'_s_sim/tot_s 
		replace vs_sim_coal=vs_sim_coal+`parti'_vs_sim
	}
	
	
	
	
	gen maj_coal_sim= ss_sim_coal>.5 if ss_sim_coal!=.
	
	gen maj_ch= maj_coal!=maj_coal_sim if maj_coal_sim !=. & maj_coal!=.
	
	bysort vch: egen freq_maj_ch= mean (maj_ch)

	gen vch_l= vch if freq>(.4+.0`diff')         /* lower bound */
	gen vch_u= vch if freq<(.6-.0`diff')		 /* upper bound */
	
	egen min_p = min (vch_l)
	egen max_p = max (vch_u)
	
    egen min_n = max (vch_l)
	egen max_n = min (vch_u)
	
	gen min = min_p
	replace min=min_n if maj_coal==0
	
	gen max = max_p
	replace max=max_n if maj_coal==0

	gen vch_s=vch if freq>.5
	egen stop  =max (vch_s)
	

	
	if stop==.{
	*stop
	}

	
end
/****************************************************/






/****************************************************/

/****************************************************/


/**

************************************Main Program************************************************/

****cd "C:\Users\Olle\Documents\simulering

do "data prep sim coal norway"
coalprep "RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2"  "RV SV DNA MDG JointL" 

clear
set obs 1
gen x=1
save "../dta/val med simuleringar.dta", replace
use "../dta/hel kommun grunddata coal.dta" 

local o =obs

forvalues  i=1 /`o'  {
keep if _n==`i'

simchange_mun  "RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2" "50000" ".""." ""

comp_mandates_mun "RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2"  ""

dist_ch_coal   "RV SV DNA MDG JointL" "0"

local max=max
local min=min

use "../dta/hel kommun grunddata coal.dta" , clear 
keep if _n==`i'
simchange_mun  "RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2""25000" "`min'""`max'" 

comp_mandates_mun "RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2"  ""

dist_ch_coal    "RV SV DNA MDG JointL"   "5"


local max=max
local min=min
use "../dta/hel kommun grunddata coal.dta" , clear

keep if _n==`i'
simchange_mun  "RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2""25000" "`min'""`max'" 

comp_mandates_mun "RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2" ""

dist_ch_coal   "RV SV DNA MDG JointL"   ""

gen vch_f= vch if freq>.5

egen dist_coal_ch_maj=max(vch_f) if  ss_coal>.5

egen dist_coal_ch_min=min (vch_f) if  ss_coal<.5

gen dist_coal_margin=dist_coal_ch_maj

replace dist_coal_margin=dist_coal_ch_min if  ss_coal<.5

keep if _n==1
keep knr electionperiod dist_coal_ch* dist_coal_margin
save temp, replace
use  "../dta/val med simuleringar.dta"
append using temp
save "../dta/val med simuleringar.dta", replace
use "../dta/hel kommun grunddata coal.dta" 
}

do "coalition simulation norway ver 2 DH.do"
