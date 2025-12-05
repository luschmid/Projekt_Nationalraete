
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


/**************************************Main Program************************************************/
clear matrix
*cd "G:\simulering\"
set more off
clear
set mem 500M
use "electoral data basic.dta"


calmindiff "m c fp kd mp s v sd nd lok1 lok2 lok3 lok4""m c fp kd mp s v sd nd lok1 lok2 lok3 lok4" "

drop sprm- mindifflok4n3
foreach parti in m c fp kd mp s v nd sd lok1 lok2 lok3 lok4{
replace mindiff`parti'n=. if  mindiff`parti'n <0
replace mindiff`parti'n=. if  mindiff`parti'n <0
}
foreach var of varlist mindiffmp- mindifflok4n{
gen `var'proc=`var'
}

foreach var of varlist lok2mandat lok1mandat sdmandat lok3mandat lok4mandat smandat- vmandat ndmandat{
gen `var'proc=`var'/antalmandat
}

dum_mindiffpartycalvs "m c fp kd mp s v sd nd lok1 lok2 lok3 lok4"  "1 075 05 04 03 025 015 02 01""




genws "m c fp kd mp s v sd nd lok1 lok2 lok3 lok4"  ""
altgenws "m c fp kd mp s sd v nd lok1 lok2 lok3 lok4"  ""

foreach parti in m c fp kd mp s v nd sd lok1 lok2 lok3 lok4{
	gen `parti'mindiffproc=mindiff`parti'nproc*100
	replace `parti'mindiffproc=-mindiff`parti'pproc*100  if abs(mindiff`parti'pproc)<abs(mindiff`parti'nproc )
}




gen i=_n




/****************************************************/
/* define program to make rdgraph inte */
/****************************************************/
capture program drop rdintgraphtlim
program rdintgraphtlim
	args outcome diffp diffn bin lim  lim2 ytit xtit yscale
		gen int`diffp'=  int(`diffp'*`bin') if `diffp'<`lim' & `lim2' & `diffp'!=0
		gen int`diffn'=  int(`diffn'*`bin') if `diffn'<`lim' & `lim2' & `diffn'!=0
	
		sort int`diffp'
		by int`diffp': egen mi`diffp'=mean(`diffp') if  `lim2'
		by int`diffp': egen mi`diffp'outcome=mean(`outcome') if  `lim2'

		sort int`diffn'
		by int`diffn': egen mi`diffn'=mean(`diffn') if  `lim2'
		by int`diffn': egen mi`diffn'outcome=mean(`outcome') if  `lim2'

		gen mi`diffp'graf =-mi`diffp' if `lim2' & `diffp'!=0
		gen mi`diffn'graf =mi`diffn' if `lim2' & `diffn'!=0
		gen zero =0

twoway (scatter mi`diffp'outcome mi`diffp'graf if `diffp'<`lim' & `lim2' & `diffp'!=0 ) (scatter mi`diffn'outcome mi`diffn'graf if   `lim2'&  `diffn'<`lim' & `diffn'!=0) , yscale(range(`yscale')) ylabel(#5)  ytitle(`ytit') xtitle(`xtit') xline(0, lwidth(medthick) lcolor(black)) legend(off)

drop int`diffp'- zero
end
***********************************



*****************************************
*********Graphs*************************
*****************************************


foreach parti in m c fp kd mp s v nd   {
replace mindiff`parti'p=mindiff`parti'p*100
replace mindiff`parti'n=mindiff`parti'n*100
}

foreach parti in lok1 sd  {
replace mindiff`parti'p=mindiff`parti'p*100
replace mindiff`parti'n=mindiff`parti'n*100
}

/*
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
*/
***********make figure A3
rdintgraphtlim  mproc mindiffmp mindiffmn 250 0.05 " valar>1980" "Vote Share" "Dist.to Seat Thres. Cons. Party" 
graph save m, replace

rdintgraphtlim    cproc mindiffcp mindiffcn 10 1 " valar>1980" "Vote Share" "Dist.to Seat Thres. Center Party" 
graph save c, replace

rdintgraphtlim    antalmandat mindifffpp mindifffpn 5 1 " valar>1980" "Vote Share" "Dist.to Seat Thres. Liberal Party" 
graph save fp, replace

rdintgraphtlim    kdproc mindiffkdp mindiffkdn 10 1 " valar>1980" "Vote Share" "Dist.to Seat Thres. Christian Dem." 
graph save kd, replace

rdintgraphtlim    mpproc mindiffmpp mindiffmpn 10 1 " valar>1980" "Vote Share" "Dist.to Seat Thres. Environ. Party"  
graph save mp, replace

rdintgraphtlim  sproc mindiffsp mindiffsn 1000 .01 " valar<2006 & valar>1980" "Vote Share" "Dist.to Seat Thres. Social Democrats" 
graph save s, replace

rdintgraphtlim    vproc mindiffvp mindiffvn 10 1 " valar>1980" "Vote Share" "Dist.to Seat Thres. Left Party" 
graph save v, replace


rdintgraphtlim    ndproc mindiffndp mindiffndn 10 1 " valar>1980" "Vote Share" "Dist.to Seat Thres. New Democracy" 
graph save  nd, replace

rdintgraphtlim    sdproc mindiffsdp mindiffsdn 5 2 " valar<2007" "Vote Share" "Dist.to Seat Thres. New Democracy" 

gr combine m.gph c.gph fp.gph kd.gph mp.gph s.gph v.gph nd.gph

****Figure A4***********

histogram mmindiffproc if abs(mmindiffproc) <1, width(.1) frequency xtitle(Dist.to Seat Thres. Consevartive Party) xline(0, lwidth(thick) lcolor(black) extend)
graph save m, replace 

histogram cmindiffproc if abs(cmindiffproc) <1, width(.1) frequency xtitle(Dist.to Seat Thres. Center Party) xline(0, lwidth(thick) lcolor(black) extend)
graph save c, replace 

histogram fpmindiffproc if abs(fpmindiffproc) <1, width(.1) frequency xtitle(Dist.to Seat Thres. Liberal Party) xline(0, lwidth(thick) lcolor(black) extend)
graph save fp, replace 

histogram kdmindiffproc if abs(kdmindiffproc) <1, width(.1) frequency xtitle(Dist.to Seat Thres. Christian Democrats) xline(0, lwidth(thick) lcolor(black) extend)
graph save kd, replace 

histogram mpmindiffproc if abs(mpmindiffproc) <1, width(.1) frequency xtitle(Dist.to Seat Thres. Environmental Party) xline(0, lwidth(thick) lcolor(black) extend)
graph save mp, replace 

histogram smindiffproc if abs(smindiffproc)<1, width(.1) frequency xtitle(Dist.to Seat Thres. Social Democrats) xline(0, lwidth(thick) lcolor(black) extend)
graph save s, replace 

histogram vmindiffproc if abs(vmindiffproc) <1, width(.1) frequency xtitle(Dist.to Seat Thres. Left Party) xline(0, lwidth(thick) lcolor(black) extend)
graph save v, replace 

histogram ndmindiffproc if abs(ndmindiffproc) <1, width(.1) frequency xtitle(Dist.to Seat Thres. New Democracy) xline(0, lwidth(thick) lcolor(black) extend)
graph save nd, replace 

gr combine m.gph c.gph fp.gph kd.gph mp.gph s.gph v.gph nd.gph

drop mandat
reshape long @mandat @mandatproc  mindiff@p mindiff@n mindiff@pproc mindiff@nproc @proc @procam @mindiffproc @diffcalvs05  @diffcalvs025 @diffcalvs01 , i(kommunkod valkrets valar) j(parti m c fp kd mp s v nd sd lok1 lok2 lok3 lok4)

sort parti valkrets kommunkod valar
gen proc_1stlead=proc[_n+1]  if parti==parti[_n+1] & kommunkod==kommunkod[_n+1] & valkrets==valkrets[_n+1] 
gen proc_1stdiff= proc[_n+1]-proc if parti==parti[_n+1] & kommunkod==kommunkod[_n+1] & valkrets==valkrets[_n+1]

foreach var in proc_1stlead proc_1stdiff{

		bysort valar parti:egen `var'_ar_mean=mean(`var')
		gen `var'_ar_korr=`var'-`var'_ar_mean
}

bysort kommunkod valar: egen valkretsmax= max(valkrets)
keep if valar<2006

drop if parti=="lok1"
drop if parti=="lok2"
drop if parti=="lok3"


/****************************************************/
/* define program to make rdgraph inte */
/****************************************************/
capture program drop rdintgraphtlim
program rdintgraphtlim
	args outcome diffp diffn bin lim  lim2 ytit xtit tit
		gen int`diffp'=  int(`diffp'*`bin') if `diffp'<`lim' & `lim2' & `diffp'!=0
		gen int`diffn'=  int(`diffn'*`bin') if `diffn'<`lim' & `lim2' & `diffn'!=0
	
		sort int`diffp'
		by int`diffp': egen mi`diffp'=mean(`diffp') if  `lim2'
		by int`diffp': egen mi`diffp'outcome=mean(`outcome') if  `lim2'

		sort int`diffn'
		by int`diffn': egen mi`diffn'=mean(`diffn') if  `lim2'
		by int`diffn': egen mi`diffn'outcome=mean(`outcome') if  `lim2'

		gen mi`diffp'graf =-mi`diffp' if `lim2' & `diffp'!=0
		gen mi`diffn'graf =mi`diffn' if `lim2' & `diffn'!=0
		gen zero =0

twoway (scatter mi`diffp'outcome mi`diffp'graf if `diffp'<`lim' & `lim2' & `diffp'!=0 ) (scatter mi`diffn'outcome mi`diffn'graf if   `lim2'&  `diffn'<`lim' & `diffn'!=0) , yscale(range(`yscale')) ylabel(#5)  ytitle(`ytit') xtitle(`xtit') xline(0, lwidth(medthick) lcolor(black)) legend(off) title(`tit')

drop int`diffp'- zero
end

 gen mindiffp_g=mindiffp
  gen mindiffn_g=mindiffn
  
 ********** make figure 4

rdintgraphtlim  mandat mindiffp_g mindiffn_g 10 1 " valar>1980 & proc<0.073 & proc>0.005" "Seats" "Distance to Seat Threshold" "Small Parties"
graph export small_seats.tif, replace

rdintgraphtlim  mandat mindiffp_g mindiffn_g 10 1 " valar>1980 & proc>0.073" "Seats" "Distance to Seat Threshold" "Large Parties"
graph export large_seats.tif, replace
