use ..\dta\panelEarly, clear

gen antalmandat= SizeOfCouncil

gen antalmandat95=SizeOfCouncil if electionperiod==7
sort knr
by knr: egen am95=max(antalmandat95) 

gen antalmandat99=SizeOfCouncil if electionperiod==8
sort knr
by knr: egen am99=max(antalmandat99) 

gen antalmandat03=SizeOfCouncil if electionperiod==9
sort knr
by knr: egen am03=max(antalmandat03)

gen antalmandat07=SizeOfCouncil if electionperiod==10
sort knr
by knr: egen am07=max(antalmandat07)


foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
replace LVotes`parti' =0 if  LVotes`parti'==.
}

foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen r`parti'=LVotes`parti'
}

/* LEGGER TIL VOTES DATA FOR Å SJEKKE PROBLEMER , SENERE */

foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen rVotes`parti'=Votes`parti'
}

foreach sam in dh msl dh95 dh99 dh03 dh07 msl95 msl99 msl03 msl07{
	foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
		gen `sam'm`parti ' =0
	}
}


gen rm =0


******dhondt*************
while rm<86 {


foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen j`parti '= r`parti'/(2+2*dhm`parti')
}
egen tm =rowtotal(dhmRV -dhmJointR2 )
egen max = rowmax(jRV-jJointR2)
foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
replace dhm`parti'= dhm`parti' +1 if j`parti '==max & tm<antalmandat & r`parti'!=.
}
 drop jRV-max tm
replace rm=rm+1
}

replace rm=0

******Saint Lague*************
/*
while rm<86 {


foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen j`parti '= r`parti'/(1+2*slm`parti')
}
egen tm =rowtotal(slmRV -slmJointR2 )
egen max = rowmax(jRV-jJointR2)
foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
replace slm`parti '= slm`parti' +1 if j`parti '==max & tm<antalmandat & r`parti'!=.
}
 drop jRV-max tm
replace rm=rm+1
}

replace rm=0
*/
******Modified Saint Lague*************
while rm<86{

foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen p`parti'=mslm`parti'>=1 & mslm`parti'!=.
}
foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen j`parti'= (1-p`parti')*r`parti'/1.4 + p`parti'* r`parti'/(1+2*mslm`parti')
}
egen tm =rowtotal(mslmRV -mslmJointR2)
egen max = rowmax(jRV-jJointR2)
foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
replace mslm`parti'= mslm`parti' +1 if j`parti'==max & tm<antalmandat & r`parti'!=.
}
 drop pRV-max tm
replace rm=rm+1
}
replace rm=0

*********************************

**********seat allocation with alternative council sizes************************
******dhondt*************

foreach year in 95 99 03 07{

while rm<86{


foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen j`parti'= r`parti'/(2+2*dh`year'm`parti')
}
egen tm =rowtotal(dh`year'mRV -dh`year'mJointR2 )
egen max = rowmax(jRV-jJointR2)
foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
replace dh`year'm`parti'= dh`year'm`parti' +1 if j`parti '==max & tm<am`year' & r`parti'!=.
}
drop jRV-max tm
replace rm=rm+1
}

replace rm=0
}

******Saint Lague*************
********
/*
foreach year in 95 99 03 07{
while rm<86 {


foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen j`parti '= r`parti'/(1+2*sl`year'm`parti')
}
egen tm =rowtotal(sl`year'mRV -sl`year'mJointR2 )
egen max = rowmax(jRV-jJointR2)
foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
replace sl`year'm`parti '= sl`year'm`parti' +1 if j`parti '==max & tm <am`year' & r`parti'!=.
}
 drop jRV-max tm
replace rm=rm+1
}

replace rm=0
}
*/

******Modified Saint Lague*************
foreach year in 95 99 03 07{

while rm<86{

foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen p`parti'=msl`year'm`parti'>=1 & msl`year'm`parti'!=.
}
foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen j`parti'= (1-p`parti')*r`parti'/1.4+ p`parti'* r`parti'/(1+2*msl`year'm`parti')
}
egen tm =rowtotal(msl`year'mRV -msl`year'mJointR2)
egen max = rowmax(jRV-jJointR2)
foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
replace msl`year'm`parti'= msl`year'm`parti' +1 if j`parti '==max & tm<am`year' & r`parti'!=.
}
 drop pRV-max tm
replace rm=rm+1
}
replace rm=0
}
***************************************************************************
*******************************************************************************

foreach sam in dh msl dh95 dh99 dh03 dh07 msl95 msl99 msl03 msl07{
	foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
		gen `sam'rm`parti'=r`parti' if `sam'm`parti'>0
	}
}

egen dhtm =rowtotal(dhmRV -dhmJointR2 )
*egen sltm =rowtotal(slmRV -slmJointR2 )
egen msltm =rowtotal(mslmRV -mslmJointR2 )

foreach year in 95 99 03 07{
egen dh`year'tm =rowtotal(dh`year'mRV -dh`year'mJointR2 )
*egen sl`year'tm =rowtotal(sl`year'mRV -sl`year'mJointR2 )
egen msl`year'tm =rowtotal(msl`year'mRV -msl`year'mJointR2 )
}

egen tr =rowtotal(rRV -rJointR2 )

foreach sam in dh msl dh95 dh99 dh03 dh07 msl95 msl99 msl03 msl07{
	egen `sam'trm =rowtotal(`sam'rmRV -`sam'rmJointR2 )
}

foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen pr`parti'=r`parti'/tr
}

foreach sam in dh msl dh95 dh99 dh03 dh07 msl95 msl99 msl03 msl07{
	foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
	gen `sam'pm`parti'=`sam'm`parti'/`sam'tm
	}
}


foreach sam in dh msl dh95 dh99 dh03 dh07 msl95 msl99 msl03 msl07{
	foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
	gen `sam'prm`parti'=r`parti'/`sam'trm
	}
}


foreach sam in dh msl dh95 dh99 dh03 dh07 msl95 msl99 msl03 msl07{
	foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
	gen `sam'bias`parti'=`sam'pm`parti'-pr`parti'
	}
}

foreach sam in dh msl dh95 dh99 dh03 dh07 msl95 msl99 msl03 msl07{
	foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
	gen `sam'pbias`parti'=`sam'pm`parti'-`sam'prm`parti'
	}
}


*foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
*gen SeatShare`parti'=Rep`parti''SizeOfCouncil
*}

*foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
*gen LVoteShare`parti'=LVoteShare`parti'/LVoteShare
*}
foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
	gen bias`parti'=SeatShare`parti'-LVoteShare`parti'
}



	

**********Number of parties with seats********

gen P=0
foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {

replace P= P+1 if Rep`parti'>0
}


**********average deviation from the threshold where a party is guaranteed a seat

gen dev=  (P-2) / (2* ( antalmandat+1) ) 

**********indices*************



foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {

gen forcing`parti'= LVoteShare`parti'* ( antalmandat+1)* (1+dev)
}

************dh loopen delar ut fel antal mandat*******************
gen av_totseats_dh= SizeOfCouncil!=dhtm

************msl loopen delar ut fel antal mandat*******************
gen av_totseats_msl= SizeOfCouncil!=msltm

************definierna felvariabel för varje parti baserat på antalmandat*****************

foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen av_m_`parti'= dhm`parti'!=Rep`parti' & electionperiod<9 & Rep`parti'!=.
replace av_m_`parti'= 1 if mslm`parti'!=Rep`parti' & electionperiod>8 & Rep`parti'!=.
}

*************summera antalet fel på kommunnivå''''''''''''
egen av_m_kommun =rowtotal(av_m_RV -av_m_JointR2)

************definierna felvariabel för varje parti baserat på seatshare*****************

foreach parti in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2  {
gen av_pm_`parti'= dhpm`parti'!=SeatShare`parti' & electionperiod<9 & SeatShare`parti'!=.
replace av_pm_`parti'= 1 if mslpm`parti'!=SeatShare`parti' & electionperiod>8 & SeatShare`parti'!=.
}

*************summera antalet fel på kommunnivå''''''''''''
egen av_pm_kommun =rowtotal(av_pm_RV -av_pm_JointR2)


**************** SER PÅ FEILENE ******
drop if entydig==0 & ele<8

li knr ele if av_totseats_dh==1|av_totseats_msl==1


******* KASTER UT OBSERVASJONER HVOR VI DELERE UT ET MANDAT FOR MYE - SKYLDES AT VI FÅR LIKHET MELLOM TO PARTI , LODDTREKNING AVGJØR LOKALT ?  ******

drop if knr==1749 & ele==9
drop if knr==1811 & ele==8
drop if knr==1942 & ele==9
drop if knr==419 & ele==6
drop if knr==438 & ele==5

keep pr* LVoteShare* VoteShare* SeatShare* bias* mslm* dhm* dhpm* mslpm* dhbias* mslbias* dhpbias* mslpbias* dh95m* msl95m* dh95pm* msl95pm* dh95bias* msl95bias* dh99m* msl99m* dh99pm* msl99pm* dh99bias* msl99bias* dh03m* msl03m* dh03pm* msl03pm* dh03bias* msl03bias* dh07m* msl07m* dh07pm* msl07pm* dh07bias* msl07bias* knr electionperiod SizeOfCouncil Rep*

saveold ..\dta\MandatfordelningOlikaMetoderBeforeReshape, replace

reshape long pr LVoteShare VoteShare SeatShare bias mslm dhm dhpm mslpm dhbias mslbias dhpbias mslpbias dh95m msl95m dh95pm msl95pm dh95bias msl95bias dh99m msl99m dh99pm msl99pm dh99bias msl99bias dh03m msl03m dh03pm msl03pm dh03bias msl03bias dh07m msl07m dh07pm msl07pm dh07bias msl07bias, j(parti) i(knr electionperiod) string

drop *bias*

saveold ..\dta\MandatfordelningOlikaMetoder, replace



************************************************************
************************************************************
******************************* AVVIKELSE ******************
************************************************************

gen avvikelse=0
replace avvikelse=1 if SeatShare!=mslpm & electionperiod>8 
replace avvikelse=1 if SeatShare!=dhpm & electionperiod<9 

sort knr electionperiod

by knr ele: egen avvikelse_kommune=sum(avvikelse)

g InSample=0
replace InSample=1 if ele==8|ele==9
sort knr InSample
by knr InSample: egen avvikelse_InSample=sum(avvikelse)

sort knr electionperiod


****** SJEKKER OM DH TREFFER ETTER REFORM ***********

g CorrespondWithDH=0
replace CorrespondWithDH=1 if SeatShare==dhpm 

*sort knr ele
*tab knr if CorrespondWithDH==1 & avvikelse==1 & ele==9
*tab knr if CorrespondWithDH==0 & avvikelse==1 & ele==9

saveold ..\dta\IndexBeforeCollapse, replace   /* FOR IDENTIFYING PROBLEMS */

***************************************

gen neg_vs=-LVoteShare
sort knr electionperiod
by knr electionperiod: egen rank_vs= rank (neg_vs), unique

forvalues n=1/4{
gen vssum_`n'p=LVoteShare if rank_vs<=`n' & LVoteShare!=.

    foreach sam in dh msl dh95 dh99 dh03 dh07 msl95 msl99 msl03 msl07{
        gen pmsum_`sam'_`n'p = `sam'pm if  rank_vs<=`n' & LVoteShare!=.
    }
}




collapse avvikelse_kommune avvikelse_InSample (sum) vssum_* pmsum_*, by(knr electionperiod)


gen min_np_maj_vs=.

forvalues n=1/4{ 
    replace min_np_maj_vs=`n' if vssum_`n'p>.5 & vssum_`n'p!=. & min_np_maj_vs==.
}

foreach sam in dh msl dh95 dh99 dh03 dh07 msl95 msl99 msl03 msl07{
    gen min_np_maj_`sam'=.
    forvalues n=1/4{ 
        replace min_np_maj_`sam'=`n' if pmsum_`sam'_`n'p>.5 & pmsum_`sam'_`n'p!=. & min_np_maj_`sam'==.
}
}

save ..\dta\Index, replace


************************************************************
************************************************************
******************************* AVVIKELSE ******************
************************************************************

use ..\dta\PanelParties 

drop *Male *Female

sort knr ele
merge knr ele using ..\dta\Index
drop _merge

drop if (entydig==0  & ele<8) /* DETTE ER PRE-MAI 2012 RESTRIKSJONEN ! */

******* KASTER UT OBSERVASJONER HVOR VI DELERE UT ETT MANDAT FOR MYE - SKYLDES AT VI FÅR LIKHET MELLOM TO PARTI , LODDTREKNING AVGJØR LOKALT ?  ******

drop if knr==419&ele==5
drop if knr==438&ele==6 

*********** KASTER UT OBSERVASJONER HVOR DET IKKE ER PERFEKT MATCH MELLOM SEATS BASERT PÅ LVOTES OG SEATS BASERT PÅ VALGSTATISTIKKEN ****

drop if avvikelse_kommune>0


*********************************************************** RATIOEN MELLOM "TOTAL LVOTES" OG "TOTAL VOTES" SKAL VÆRE EKSAKT COUNCIL SIZE *****

gen QQ= LVotesNEW/LVotes
gen ZZ= LVotesNEW/VotesNEW

gen WrongRatio=0
replace WrongRatio=1 if ZZ!=SizeOfCouncil
gen YYY_diff=LVotesNEW-VotesNEW*SizeOfCouncil
gen YYY_ratio=LVotesNEW/(VotesNEW*SizeOfCouncil)
*bysort SizeOfCouncil: tab YYY_ratio

drop if WrongRatio==1 

**** BURDE SPILLE LITEN ROLLE OM DISSE KASTES UT ELLER IKKE ETTERSOM VI ALLEREDE HAR KASTET UT DE HVOR VI IKKE REPRODUSERER SEATS FRA VALGSTATISTIKKEN ****


********************************************************************

drop LVotes
rename LVotesNEW LVotes   /* DETTE ER TOTAL VOTES */

/* OPPDATERER LVOTESHARE TIL BARE Å OMHANDLE DE SOM FIKK MINST ETT MANDAT */

foreach party in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2 { 
g LVoteShare`party'_p=LVoteShare`party' if SeatShare`party'>0
} 

egen votetotal=rowtotal(LVoteShareRV_p-LVoteShareJointR2_p)    /* vote total tar verdien 1 dersom alle lister som stilte fikk mandater */

foreach party in RV SV DNA V SP KRF H FRP PP MDG Indep1 Indep2 Indep3 Indep4 Indep5 Indep6 Other1 Other2 Other3 Other4 Other5 JointL JointR1 JointR2 { 
gen LVoteShare`party'_new=LVoteShare`party'/votetotal 
} 

****************************************** DENNE RESTRIKSJONEN VAR SLÅTT PÅ FØR JUNI 2012 *******************

*drop if LVoteShareOther2>0
*drop if LVoteShareIndep2>0

*************************************************************************************************************

*drop *SQ* *Max *Indep2* *Indep3* *Indep4* *Indep5* *Indep6* *Other2* *Other3* *Other4* *Other5* *JointR2* /* Indep 2 - 6 never >0 ,when entydig==1, drop very few obs with other2-Other5 positive */
drop *SQ* *Max  /* Indep 2 - 6 never >0 ,when entydig==1, drop very few obs with other2-Other5 positive */


*************** drop *Other1* *Joint* 


renpfix SeatShare seatshare_

keep knr ele LVoteShare* seatshare_* Rep* SizeOfCouncil Available*
drop *_new *_p

/* lumping joint lists together , they are always separate */
/*
g forcing_Joint=forcing_JointL+forcing_JointR
g LVoteShareJoint=LVoteShareJointL+LVoteShareJointR
g stepshare_Joint=stepshare_JointL+stepshare_JointR
g seatshare_Joint=seatshare_JointL+seatshare_JointR
g RepJoint=RepJointL+RepJointR
g AvailableJoint=AvailableJointL+AvailableJointR

drop *JointL* *JointR*
*/

gen SeatShareLEFT= seatshare_RV+ seatshare_SV + seatshare_DNA + seatshare_JointL  + seatshare_MDG /* + seatshare_SP */
gen LVoteShareLEFT= LVoteShareRV+ LVoteShareSV + LVoteShareDNA + LVoteShareJointL + LVoteShareMDG /* + LVoteShareSP */

sort knr ele
gen NextLVoteShareLEFT= f.LVoteShareRV+ f.LVoteShareSV + f.LVoteShareDNA + f.LVoteShareJointL + f.LVoteShareMDG /* + f.LVoteShareSP */
*gen NextLVoteShareLEFT= (f.LVoteShareRV+ f.LVoteShareSV + f.LVoteShareDNA)-LVoteShareLEFT

gen marginLEFT=SeatShareLEFT-.5
gen MajLeft=marginLEFT>0

gen OneSeatMajority=0
gen Closeness=abs(marginLEFT)
foreach CouncilSize in 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 75 77 79 81 83 85{
replace OneSeatMajority=1 if (SizeOfCouncil==`CouncilSize' & Closeness<=(1/`CouncilSize'))
}

sort knr ele

saveold ..\dta\SeatAllocations, replace
