clear
insheet using "..\dta\Individuallevel_LR_PTAX.txt"

gen LeftRight_RV=hv if parti=="2:RV"
gen LeftRight_SV=hv if parti=="3:SV"
gen LeftRight_DNA=hv if parti=="4:DNA"
gen LeftRight_SP=hv if parti=="5:SP"
gen LeftRight_KRF=hv if parti=="6:KrF"
gen LeftRight_V=hv if parti=="7:V"
gen LeftRight_H=hv if parti=="8:H"
gen LeftRight_FRP=hv if parti=="9:FrP"
gen LeftRight_PP=hv if parti=="10:PP"
gen LeftRight_MDG=hv if parti=="11:MD"
gen LeftRight_Various=hv if parti=="1:ANT"

gen Ptax_RV=ptax if parti=="2:RV"
gen Ptax_SV=ptax if parti=="3:SV"
gen Ptax_DNA=ptax if parti=="4:DNA"
gen Ptax_SP=ptax if parti=="5:SP"
gen Ptax_KRF=ptax if parti=="6:KrF"
gen Ptax_V=ptax if parti=="7:V"
gen Ptax_H=ptax if parti=="8:H"
gen Ptax_FRP=ptax if parti=="9:FrP"
gen Ptax_PP=ptax if parti=="10:PP"
gen Ptax_MDG=ptax if parti=="11:MD"
gen Ptax_Various=ptax if parti=="1:ANT"


foreach party in RV SV DNA V SP KRF H FRP {
egen number_of_respondents_`party'=count(Ptax_`party'), by(knr)
}
gen number_of_respondents=number_of_respondents_RV+number_of_respondents_SV+number_of_respondents_DNA+number_of_respondents_V+number_of_respondents_SP+number_of_respondents_KRF+number_of_respondents_H+number_of_respondents_FRP
drop if number_of_respondents<3 /* THESE ARE MUNICIPALITIES WHERE ONLY MAYOR AND DEP MAYOR HAVE RESP. */

drop Ptax*

collapse LeftRight_*, by(knr aar)

************************** MERGING DATA ON DISTANCE TO MAJ. THRESHOLD *****************************

gen electionperiod=10 
replace electionperiod=9 if aar==2007
replace electionperiod=8 if aar==2003

sort knr ele
merge knr ele using "../dta/val med simuleringar DH4.dta"
drop _merge
bysort ele: sum dist_coal_margin
sort knr ele
merge knr ele using "../dta/val med simuleringar4.dta", update
bysort ele: sum dist_coal_margin
drop _merge

sort knr ele
merge knr ele using "../dta/MergedData.dta"
drop _merge

drop if RV_dum_025==.

drop *RV* *PP* *MDG*  /* VERY FEW RESPONDENTS */

*keep if abs(dist_coal_margin)<.5  /* limiting sample to only 1122 obs. where we also have simulations for distance to threshold */

******** TABLE ***************************

gen Various_dum_025=0
replace Various_dum_025=0.5 if Various_dum_w_025>0
replace Various_dum_025=. if Various_dum_w_025==.

log using ../tables/TableB6.txt, text replace

foreach party in SV DNA V SP KRF H FRP Various{
ttest LeftRight_`party' , by(`party'_dum_025)
}

log close


gen CloseMaj=0
replace CloseMaj=1 if abs(dist_coal_margin)<0.05
replace CloseMaj=. if dist_coal_margin==.


foreach party in SV DNA V SP KRF H FRP Various{
ttest LeftRight_`party' , by(CloseMaj)
}

g ZZ=1-CloseMaj

eststo clear
estpost summarize LeftRight_* if CloseMaj==1
eststo a
estpost summarize LeftRight_* if CloseMaj==0
eststo b
estpost ttest LeftRight_* , by(ZZ)
eststo x
esttab a b x using ..\tables\TableB7.tex, replace title(Party positions by electoral closeness of party blocs\label{ExternalValidity}) style(tex) ///
b(3) mtitles("Close" "Not close" "Difference") collabels(Mean Estimate SD SE) starlevels(* 0.10 ** 0.05 *** 0.01) ///
cells("mean (fmt(3) pattern(1 1 0 0)) b (fmt(3) star pattern (0 0 1 1)) sd (fmt(3) par pattern(1 1 0 0)) se (fmt(3) par pattern(0 0 1 1))") 
