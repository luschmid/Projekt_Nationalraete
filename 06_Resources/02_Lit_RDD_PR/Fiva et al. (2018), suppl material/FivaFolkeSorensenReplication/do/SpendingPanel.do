use ..\dta\baseline2012.dta, clear

sort knr year

keep if electionperiod>3  /* Only full election data from 1983. Spending data ends in 1999 */

keep if year<2011 /* MANGLER DATA FOR SISTE ÅR I VALGPERIODE 10 */

collapse PTAX120sqm pop children young elderly Total_* PerCap* Herfindahl VoteNatLeft VoteShareLeft VoteShareRight VoteShareOther Balanced* cnr age06 age715 age1620 age2125 age2630 age3135 age3640 age4145 age4650 age5155 age5660 age6165 age6670 age7175 age7680 age81 unemployment women, by(knr electionperiod)

tab electionperiod, g(d_ele)

/* GENERATING DEPENDENT VARIABLE */

********************** SPENDING **************************************************************
g Total=Total_childcare + Total_education + Total_elderlycare + Total_healthsocial + Total_culture + Total_transport + Total_centraladm + Total_other
foreach sector in childcare education elderlycare healthsocial culture transport centraladm other { 
gen TotalSH_`sector'=Total_`sector'*100/Total
}

gen SpendingYoungVsOld=(TotalSH_childcare + TotalSH_education)*100/(TotalSH_childcare + TotalSH_education + TotalSH_elderlycare + TotalSH_healthsocial) 

********************** PTAX **************************************************************

gen DPTAX=0
replace DPTAX=1 if PerCapPTAX>0
replace DPTAX=. if PerCapPTAX==.

gen DPTAX_2500ex=0
replace DPTAX_2500ex=1 if PerCapPTAX>0
replace DPTAX_2500ex=. if pop<2500

/***************** HVORDAN SKILLE COMMERCIAL FRA RESIDENTIAL *****************************/
sort knr
merge knr using ..\dta\FivaRattso
drop _merge
sort knr ele
merge knr ele using ..\dta\ptax2007
g PTaxPower2007= ptax_pow2007/pop
g PTaxResidential2007=ptax_res2007/pop

/* dropping all municipalities above 2k per capita in HydroPowerIncome in 2007 - from stakes paper */

replace PerCapPTAX=. if PTaxPower2007>2

/* for ptax need to exclude small municip */

gen PerCapPTAX_2500ex=PerCapPTAX
replace PerCapPTAX_2500ex=. if pop<2500

gen PerCapPTAX_4000ex=PerCapPTAX
replace PerCapPTAX_4000ex=. if pop<4000

replace PTaxResidential2007=. if ele<10

****************************************************************************
*********************** NEW DPTAX RESIDENTIAL VARIABLES OCT 2013 ***********
****************************************************************************

gen DPTAXres=0
replace DPTAXres=1 if (dptaxFivaRattso==1 & ele==8)
replace DPTAXres=1 if (PerCapPTAXresidential>0 & ele>8)
replace DPTAXres=. if ele<8

replace PTAX120sqm=ptaxFivaRattso if ele==8  /* fiva-rattso gjelder for 2001 */
replace PTAX120sqm=PTAX120sqm/0.833589 if ele==8  /* defalting 2001 to 2011 prices */
replace PTAX120sqm=PTAX120sqm/(1.6/1.2) if ele==8  /* fiva-rattso IS FOR 160 SQM house */

bysort ele: sum PTAX120sqm 

*STANDARDIZING
************* 
replace PTAX120sqm=PTAX120sqm/797.106 if ele==8  /* MEASURED BY SURVEY IN ELE==8 */
replace PTAX120sqm=PTAX120sqm/1232.869 if ele>8   /* MEASURED BY STATS NORWAY SINCE 2007 */

****************************************************************************
****************************************************************************
****************************************************************************

replace pop=pop/1000
rename pop population
rename SpendingYoungVsOld Spending
rename PerCapUserCharges Fees
renpfix DPTAX_ DPTAX


**** time f.e.

rename d_ele2 ElectionPeriod8791
rename d_ele3 ElectionPeriod9295
rename d_ele4 ElectionPeriod9699
rename d_ele5 ElectionPeriod9903
rename d_ele6 ElectionPeriod0407
rename d_ele7 ElectionPeriod0811 /* FORELØPIG HAR VI BARE SPENDING DATA FOR 2008-2010 , */

sort knr ele

gen Competition=0
*replace Competition=1 if cnr==1|cnr==2|cnr==3|cnr==6|cnr==7|cnr==12|cnr==16
replace Competition=1 if cnr==11|cnr==7|cnr==15|cnr==1|cnr==12|cnr==10|cnr==16
*keep if Competition==0

drop if ele==. | ele<8

/* JANUARY 2016 - error checking reveals that TWO ENTRIES HAVE MISSING DATA, WHICH ARE CODED AS ZERO in baseline2012 */
**** målselv (1924) in 2002 (ele 8) and torsken (1928) in 2004 (ele 9)
**** torsken is missing anyway, so this does not change any results

drop if knr==1924 & ele==8
drop if knr==1928 & ele==9

save ..\dta\SpendingPanel, replace




*********************************************************************
*********************************************************************
*********************************************************************
*********************************************************************


drop _merge

sort knr ele
merge knr ele using ..\dta\SeatAllocations
keep if _merge==3     
drop _merge

xtset knr ele

sort knr ele
merge knr ele using ..\dta\DistansKod
keep if _merge==3       /* observations both (entydig) in political data and exist in spending data */
drop _merge

drop *_1 *_075 *_04 *_03 *_015 *_02


save ..\dta\MergedDataBeforeCollapse, replace


/**** SLÅR SAMMEN JOINTR1 OG JOINTR2 */

foreach var in _treat_w_05 _treat_w_025 _treat_w_01 _treat_05 _treat_025 _treat_01 {
gen JOINTR`var'=JointR1`var'+ JointR2`var'
}
gen JOINTR_dum_w_05=abs(JOINTR_treat_w_05)
gen JOINTR_dum_w_025=abs(JOINTR_treat_w_025)
gen JOINTR_dum_w_01=abs(JOINTR_treat_w_01)

gen seatshare_JOINTR=seatshare_JointR1+seatshare_JointR2
gen LVoteShareJOINTR=LVoteShareJointR1+LVoteShareJointR2

drop *JointR*


/**** SLÅR SAMMEN SEATSHARE FOR OTHER OG INDEP */

foreach var in _treat_w_05 _treat_w_025 _treat_w_01 _treat_05 _treat_025 _treat_01 {
gen Various`var'=Indep1`var'+ Indep2`var'+ Indep3`var'+ Indep4`var'+ Indep5`var'+ Indep6`var'+ Other1`var'+ Other2`var'+ Other3`var'+ Other4`var'+ Other5`var'
}
gen Various_dum_w_05=abs(Various_treat_w_05)
gen Various_dum_w_025=abs(Various_treat_w_025)
gen Various_dum_w_01=abs(Various_treat_w_01)

gen seatshare_Various=seatshare_Other1+seatshare_Other2+seatshare_Other3+seatshare_Indep1+seatshare_Indep2+seatshare_Indep3+seatshare_Indep4+seatshare_Indep5
gen LVoteShareVarious=LVoteShareOther1+LVoteShareOther2+LVoteShareOther3+LVoteShareIndep1+LVoteShareIndep2+LVoteShareIndep3+LVoteShareIndep4+LVoteShareIndep5

drop *Other* *Indep* /* OLD CATEGORIES */
drop  d_ele1 ElectionPeriod8791 ElectionPeriod9295 ElectionPeriod9699 ElectionPeriod9903

drop if Total==0 /* MISSING FOR TORSKEN */

/************* DATA FOR HISTORICAL TOWN STATUS ************/

sort knr 
merge knr using ..\dta\rsue_data
replace by1911=0 if by1911==.
drop _merge

/******************* DATA FOR RURAL ******************/

sort knr ele
merge knr ele using ..\dta\Rural\spredtbygd2000
drop _merge

sort knr ele
merge knr ele using ..\dta\Rural\spredtbygd2004, update
drop _merge

sort knr ele
merge knr ele using ..\dta\Rural\spredtbygd2008, update
drop _merge

gen rural_new=spredtbygd/population /* time varying rural with obs. for all municip */

drop rural
rename rural_new rural

replace population=log(population*1000)
renpfix TotalSH_ 

gen PerCapPTAXFees=  PerCapPTAX+Fees


/* OLD POLICY INDEX */
*sort knr ele
*merge knr ele using ..\dta\PreferenceIndex_NationalWeights
*merge knr ele using ..\dta\PreferenceIndex_National
*keep if _merge==3
*drop _merge

/* NEW POLICY INDEX*/

gen zz=1
sort zz
merge zz using ..\dta\PolicyIndexWideAccurate
drop _merge

gen TestSum= seatshare_RV+seatshare_SV+ (seatshare_DNA+seatshare_JointL)+seatshare_V+seatshare_SP+seatshare_KRF+(seatshare_H+seatshare_JOINTR)+seatshare_FRP+seatshare_PP+seatshare_MDG+seatshare_Various
gen TestSum2= LVoteShareRV+LVoteShareSV+ (LVoteShareDNA+LVoteShareJointL)+LVoteShareV+LVoteShareSP+LVoteShareKRF+(LVoteShareH+LVoteShareJOINTR)+LVoteShareFRP+LVoteSharePP+LVoteShareMDG+LVoteShareVarious

foreach policyarea in LeftRight DPTAXres Fees childcare education elderlycare healthsocial culture transport centraladm other{
gen `policyarea'Index=seatshare_RV*`policyarea'RV+seatshare_SV*`policyarea'SV +(seatshare_DNA+seatshare_JointL)*`policyarea'DNA +seatshare_V*`policyarea'V +seatshare_SP*`policyarea'SP +seatshare_KRF*`policyarea'KRF +(seatshare_H+seatshare_JOINTR)*`policyarea'H+seatshare_FRP*`policyarea'FRP+seatshare_PP*`policyarea'PP +seatshare_MDG*`policyarea'MDG +seatshare_Various*`policyarea'OTH 
gen `policyarea'Control=LVoteShareRV*`policyarea'RV+LVoteShareSV*`policyarea'SV +(LVoteShareDNA+LVoteShareJointL)*`policyarea'DNA +LVoteShareV*`policyarea'V +LVoteShareSP*`policyarea'SP +LVoteShareKRF*`policyarea'KRF +(LVoteShareH+LVoteShareJOINTR)*`policyarea'H+LVoteShareFRP*`policyarea'FRP+LVoteSharePP*`policyarea'PP +LVoteShareMDG*`policyarea'MDG +LVoteShareVarious*`policyarea'OTH 
}

gen TotalIndex=LeftRightIndex
gen TotalControl=LeftRightControl

sort knr ele
save ..\dta\MergedDataBeforeStandardizing, replace


************* STANDARDIZING **********************
sum DPTAX DPTAXres Fees childcare education elderlycare healthsocial culture transport centraladm other Total

replace DPTAX=DPTAX/.4792938
replace DPTAXres=DPTAXres/.4630886
replace Fees=Fees/1.310879
replace childcare=childcare/2.666829
replace education=education/4.247392
replace elderlycare=elderlycare/5.065404
replace healthsocial=healthsocial/2.387925
replace culture=culture/2.272002
replace transport=transport/1.372035
replace centraladm=centraladm/2.746329
replace other=other/3.358745
replace Total=Total/21.07309

sum DPTAX DPTAXres Fees childcare education elderlycare healthsocial culture transport centraladm other Total

sum *Index if forcing_SP!=.

replace LeftRightIndex=LeftRightIndex/.571101
replace DPTAXresIndex=DPTAXresIndex/0.07566
replace FeesIndex=FeesIndex/0.0374517
replace childcareIndex=childcareIndex/0.0266882
replace educationIndex=educationIndex/0.0243838
replace elderlycareIndex=elderlycareIndex/0.0226637
replace healthsocialIndex=healthsocialIndex/0.0157781
replace cultureIndex=cultureIndex/0.0353937 
replace transportIndex=transportIndex/0.028332 
replace centraladmIndex=centraladmIndex/0.0274125
replace otherIndex=otherIndex/0.0112945

sum *Index if forcing_SP!=.

gen DPTAXIndex=DPTAXresIndex
gen DPTAXControl=DPTAXresControl
gen PTAX120sqmIndex=DPTAXresIndex
gen PTAX120sqmControl=DPTAXresControl


sort knr ele
save ..\dta\MergedData, replace

