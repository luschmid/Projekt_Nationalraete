*****************************************************************
***************************************************************
/****************************************************/
/* define program to make graph with seat share and vote share zero bin and fixed observations in bins*/
/****************************************************/
capture program drop vsssfz
program vsssfz

	*******the arguments are the outcome variable, the forcing variable, restriction on the sample used, ***
	*****restriction on the distance from the threshold, Y axis title, X axis title**************
	args ss vs  bin lim lim1 tit1 tit2 grtit grstit
		
		********define ranks to construct bins****************
		gen  posdis= `vs' if `vs'>0 & `lim' &`ss'!=. &`vs'!=.
		egen rankpos=  rank (posdis) if `vs'>0 & `lim' &`ss'!=. &`vs'!=. , unique

		gen  negdis= -`vs' if `vs'<0 & `lim' &`ss'!=. &`vs'!=.
		egen rankneg =  rank (negdis) if `vs'<0 & `lim' &`ss'!=. &`vs'!=. , unique

		********define integrers to construct bins****************
		gen int`vs'pos=  int(rankpos/`bin') if `vs'>0 & `lim' &`ss'!=. &`vs'!=.
		gen int`vs'neg=  int(rankneg/`bin') if `vs'<0 & `lim' &`ss'!=. &`vs'!=.

		
		************create bin averages of both forcing variable and outcome*************
		sort int`vs'pos
		by int`vs'pos: egen mean`vs'pos=mean(`vs') if `vs'<`lim1' & int`vs'pos!=.
		by int`vs'pos: egen mean`ss'pos=mean(`ss') if `vs'<`lim1' &`ss'!=. & int`vs'pos!=.

		sort int`vs'neg
		by int`vs'neg: egen mean`vs'neg=mean(`vs') if `vs'> -`lim1' & int`vs'neg !=.
		by int`vs'neg: egen mean`ss'neg=mean(`ss') if `vs'> - `lim1' &`ss'!=. & int`vs'neg!=.

		*********create a specific bin for when the forcing variable is zero******************************
            replace mean`vs'neg=0 if `vs'==0 & `vs'> -`lim1' 
		gen zero`ss'=`ss' if `vs'==0	
		egen meanzero`ss'=mean(zero`ss')
		replace mean`ss'neg=meanzero`ss' if `vs'==0

		************create single binned variables of forcing variable and outcome**************
		gen mean`vs'=mean`vs'pos if `vs'>0 & `lim'
		replace mean`vs'=mean`vs'neg if `vs'<=0 & `lim'

		gen mean`ss'=mean`ss'pos if `vs'>0 & `lim'
		replace mean`ss'=mean`ss'neg if `vs'<=0 & `lim'
		


 

twoway (scatter mean`ss' mean`vs' , msize(small)),title(`grtit') subtitle(`grstit')  ytitle(`tit1') xtitle(`tit2')legend(off) xline(0, lwidth(medthick) lcolor(black))

***************drop all generated variables*********************************
drop  posdis -rankneg mean`ss' mean`vs' mean`ss'pos mean`vs'pos mean`ss'neg mean`vs'neg  int`vs'pos int`vs'neg  zero`ss' meanzero`ss'



end


/****************************************************/
/* define program to make graph with seat share and vote share*/
/****************************************************/
capture program drop vsssz
program vsssz

	*******the arguments are the outcome variable, the forcing variable, restriction on the sample used, ***
	*****restriction on the distance from the threshold, Y axis title, X axis title**************
	args ss vs  bin lim lim1 tit1 tit2 grtit grstit
		
		********define integrers to construct bins****************
		gen int`vs'pos=  int(`vs'*`bin') if `vs'>0 & `lim' &`ss'!=. &`vs'!=.
		gen int`vs'neg=  int(`vs'*`bin') if `vs'<0 & `lim' &`ss'!=. &`vs'!=.

		
		************create bin averages of both forcing variable and outcome*************
		sort int`vs'pos
		by int`vs'pos: egen mean`vs'pos=mean(`vs') if `vs'<`lim1' & int`vs'pos!=.
		by int`vs'pos: egen mean`ss'pos=mean(`ss') if `vs'<`lim1' &`ss'!=.& int`vs'pos!=.

		sort int`vs'neg
		by int`vs'neg: egen mean`vs'neg=mean(`vs') if `vs'> -`lim1' & int`vs'neg!=.
		by int`vs'neg: egen mean`ss'neg=mean(`ss') if `vs'> - `lim1' &`ss'!=.& int`vs'neg!=.

		*********create a specific bin for when the forcing variable is zero******************************
            replace mean`vs'neg=0 if `vs'==0 & `vs'> -`lim1' 
		gen zero`ss'=`ss' if `vs'==0	
		egen meanzero`ss'=mean(zero`ss')
		replace mean`ss'neg=meanzero`ss' if `vs'==0

		************create single binned variables of forcing variable and outcome**************
		gen mean`vs'=mean`vs'pos if `vs'>0 & `lim'
		replace mean`vs'=mean`vs'neg if `vs'<=0 & `lim'

		gen mean`ss'=mean`ss'pos if `vs'>0 & `lim'
		replace mean`ss'=mean`ss'neg if `vs'<=0 & `lim'
		


 

twoway (scatter mean`ss' mean`vs' , msize(small)),title(`grtit') subtitle(`grstit')  ytitle(`tit1') xtitle(`tit2')legend(off) xline(0, lwidth(medthick) lcolor(black))

***************drop all generated variables*********************************
drop mean`ss' mean`vs' mean`ss'pos mean`vs'pos mean`ss'neg mean`vs'neg  int`vs'pos int`vs'neg  zero`ss' meanzero`ss'



end

*****************************************************************
***************************************************************

/****************************************************/
/* define program to make graph with seat share and vote share*/
/****************************************************/
capture program drop vsss
program vsss

	*******the arguments are the outcome variable, the forcing variable, restriction on the sample used, ***
	*****restriction on the distance from the threshold, Y axis title, X axis title**************
	args ss vs  bin lim lim1 tit1 tit2 grtit grstit

		********define integrers to construct bins****************
		gen int`vs'pos=  int(`vs'*`bin') if `vs'>0 & `lim' &`ss'!=. &`vs'!=.
		gen int`vs'neg=  int(`vs'*`bin') if `vs'<0 & `lim' &`ss'!=. &`vs'!=.
		
		************create bin averages of both forcing variable and outcome*************
		sort int`vs'pos
		by int`vs'pos: egen mean`vs'pos=mean(`vs') if `vs'<`lim1' 
		by int`vs'pos: egen mean`ss'pos=mean(`ss') if `vs'<`lim1' &`ss'!=.

		sort int`vs'neg
		by int`vs'neg: egen mean`vs'neg=mean(`vs') if `vs'> -`lim1' 
		by int`vs'neg: egen mean`ss'neg=mean(`ss') if `vs'> - `lim1' &`ss'!=.

		************create single binned variables of forcing variable and outcome**************
		gen mean`vs'=mean`vs'pos if `vs'>0 & `lim'
		replace mean`vs'=mean`vs'neg if `vs'<0 & `lim'

		gen mean`ss'=mean`ss'pos if `vs'>0 & `lim'
		replace mean`ss'=mean`ss'neg if `vs'<0 & `lim'
		



twoway (scatter mean`ss' mean`vs' , msize(small)), title(`grtit') subtitle(`grstit') ytitle(`tit1') xtitle(`tit2')legend(off) xline(0, lwidth(medthick) lcolor(black))


***************drop all generated variables*********************************

drop mean`ss' mean`vs' mean`ss'pos mean`vs'pos mean`ss'neg mean`vs'neg  int`vs'pos int`vs'neg 

end

**************************************************
******************************************************

/****************************************************/
/* define program to make split RD graph with seat share and vote share*/
/****************************************************/

capture program drop vsss2
program vsss2
	*******the arguments are the outcome variable, the forcing variable, restriction on the sample used, ***
	******* limit to define first subsample, condition used for second subsample***************
	*****restriction on the distance from the threshold, Y axis title, X axis title,**************
	*****definition of first subsample, definition och second subsample***********************
	args ss vs  bin lim lim1 lim2 lim3 tit1 tit2 split1 split2 grtit grstit

		********define integrers to construct bins****************
		gen int`vs'pos=  int(`vs'*`bin') if `vs'>0 & `lim' &`ss'!=. &`vs'!=.
		gen int`vs'neg=  int(`vs'*`bin') if `vs'<0 & `lim' &`ss'!=. &`vs'!=.

		
		************create bin averages of both forcing variable and outcome for first subsample*************
		sort int`vs'pos
		by int`vs'pos: egen mean1`vs'pos=mean(`vs') if `vs'<`lim3'  & `lim1'
		by int`vs'pos: egen mean1`ss'pos=mean(`ss') if `vs'<`lim3' &`ss'!=. & `lim1'

		sort int`vs'neg
		by int`vs'neg: egen mean1`vs'neg=mean(`vs') if `vs'> -`lim3' & `lim1'
		by int`vs'neg: egen mean1`ss'neg=mean(`ss') if `vs'> - `lim3' &`ss'!=. & `lim1'


		************create bin averages of both forcing variable and outcome for second subsample*************
		sort int`vs'pos
		by int`vs'pos: egen mean2`vs'pos=mean(`vs') if `vs'<`lim3'  & `lim2'
		by int`vs'pos: egen mean2`ss'pos=mean(`ss') if `vs'<`lim3' &`ss'!=. & `lim2'


		sort int`vs'neg
		by int`vs'neg: egen mean2`vs'neg=mean(`vs') if `vs'> -`lim3' & `lim2'
		by int`vs'neg: egen mean2`ss'neg=mean(`ss') if `vs'> - `lim3' &`ss'!=. & `lim2'


		************create single binned variables of forcing variable and outcome for first subsample**************
		gen mean1`vs'=mean1`vs'pos if `vs'>0 & `lim'
		replace mean1`vs'=mean1`vs'neg if `vs'<0 & `lim'

		gen mean1`ss'=mean1`ss'pos if `vs'>0 & `lim'
		replace mean1`ss'=mean1`ss'neg if `vs'<0 & `lim'

		************create single binned variables of forcing variable and outcome for second subsample**************
		gen mean2`vs'=mean2`vs'pos if `vs'>0 & `lim'
		replace mean2`vs'=mean2`vs'neg if `vs'<0 & `lim'

		gen mean2`ss'=mean2`ss'pos if `vs'>0 & `lim'
		replace mean2`ss'=mean2`ss'neg if `vs'<0 & `lim'
		



twoway (scatter mean1`ss' mean1`vs' if mean1`vs' < `lim3' & mean1`vs' > -`lim3', msize(medium)) (scatter mean2`ss' mean2`vs' if mean2`vs' < `lim3' & mean2`vs' > -`lim3',msymbol(x) msize(large)) ,title(`grtit') subtitle(`grstit')  ytitle(`tit1') xtitle(`tit2') legend(on order(1 "`split1'" 2 "`split2'")) xline(0, lwidth(medthick) lcolor(black))

***************drop all generated variables*********************************
drop int`vs'pos - mean2`ss'

end

***********************************************************
******************************************************

**************************************************
******************************************************

/****************************************************/
/* define program to make split RD graph with seat share and vote share and zero*/
/****************************************************/

capture program drop vsssz2
program vsssz2
	*******the arguments are the outcome variable, the forcing variable, restriction on the sample used, ***
	******* limit to define first subsample, condition used for second subsample***************
	*****restriction on the distance from the threshold, Y axis title, X axis title,**************
	*****definition of first subsample, definition och second subsample***********************
	args ss vs  bin lim lim1 lim2 lim3 tit1 tit2 split1 split2 grtit grstit

		********define integrers to construct bins****************
		gen int`vs'pos=  int(`vs'*`bin') if `vs'>0 & `lim' &`ss'!=. &`vs'!=.
		gen int`vs'neg=  int(`vs'*`bin') if `vs'<0 & `lim' &`ss'!=. &`vs'!=.

		******** gen mean for zeroes*********
		gen zero`ss'1=`ss' if `vs'==0 & `lim1'
		egen meanzero`ss'1 =mean(zero`ss'1)
		gen zero`ss'2=`ss' if `vs'==0 & `lim2'
		egen meanzero`ss'2 =mean(zero`ss'2)

		
		************create bin averages of both forcing variable and outcome for first subsample*************
		sort int`vs'pos
		by int`vs'pos: egen mean1`vs'pos=mean(`vs') if `vs'<`lim3'  & `lim1' & int`vs'pos!=.
		by int`vs'pos: egen mean1`ss'pos=mean(`ss') if `vs'<`lim3' &`ss'!=. & `lim1' & int`vs'pos!=.

		sort int`vs'neg
		by int`vs'neg: egen mean1`vs'neg=mean(`vs') if `vs'> -`lim3' & `lim1' & int`vs'neg!=.
		by int`vs'neg: egen mean1`ss'neg=mean(`ss') if `vs'> - `lim3' &`ss'!=. & `lim1' & int`vs'neg!=.


		************create bin averages of both forcing variable and outcome for second subsample*************
		sort int`vs'pos
		by int`vs'pos: egen mean2`vs'pos=mean(`vs') if `vs'<`lim3'  & `lim2'& int`vs'pos!=.
		by int`vs'pos: egen mean2`ss'pos=mean(`ss') if `vs'<`lim3' &`ss'!=. & `lim2'& int`vs'pos!=.


		sort int`vs'neg
		by int`vs'neg: egen mean2`vs'neg=mean(`vs') if `vs'> -`lim3' & `lim2'& int`vs'neg!=.
		by int`vs'neg: egen mean2`ss'neg=mean(`ss') if `vs'> - `lim3' &`ss'!=. & `lim2'& int`vs'neg!=.


		************create single binned variables of forcing variable and outcome for first subsample**************
		gen mean1`vs'=mean1`vs'pos if `vs'>0 & `lim'
		replace mean1`vs'=mean1`vs'neg if `vs'<0 & `lim'

		gen mean1`ss'=mean1`ss'pos if `vs'>0 & `lim'
		replace mean1`ss'=mean1`ss'neg if `vs'<0 & `lim'

		************create single binned variables of forcing variable and outcome for second subsample**************
		gen mean2`vs'=mean2`vs'pos if `vs'>0 & `lim'
		replace mean2`vs'=mean2`vs'neg if `vs'<0 & `lim'

		gen mean2`ss'=mean2`ss'pos if `vs'>0 & `lim'
		replace mean2`ss'=mean2`ss'neg if `vs'<0 & `lim'
		


twoway (scatter mean1`ss' mean1`vs' if abs(mean1`vs') < `lim3' , msize(medium)) (scatter mean2`ss' mean2`vs' if abs(mean2`vs' )< `lim3',msymbol(x) msize(large)) (scatter meanzero`ss'1  `vs' if `vs'==0, mcolor(navy) msize(medium) msymbol(circle)) (scatter meanzero`ss'2  `vs' if `vs'==0, mcolor(maroon) msize(large) msymbol(x)),title(`grtit') subtitle(`grstit')  ytitle(`tit1') xtitle(`tit2') legend(on order(1 "`split1'" 2 "`split2'")) xline(0, lwidth(medthick) lcolor(black)) 

***************drop all generated variables*********************************
drop int`vs'pos - mean2`ss'

end

***********************************************************
******************************************************


/****************************************************/
/* define program to make split RD graph with seat share and vote share and zero fixed number in bins*/
/****************************************************/

capture program drop vsssfz2
program vsssfz2
	*******the arguments are the outcome variable, the forcing variable, restriction on the sample used, ***
	******* limit to define first subsample, condition used for second subsample***************
	*****restriction on the distance from the threshold, Y axis title, X axis title,**************
	*****definition of first subsample, definition och second subsample***********************
	args ss vs  bin lim lim1 lim2 lim3 tit1 tit2 split1 split2 grtit grstit

	********define ranks to construct bins****************
		gen  posdis1= `vs' if `vs'>0 & `lim' & `ss'!=. & `vs'!=. & `lim1'
		gen  posdis2= `vs' if `vs'>0 & `lim' & `ss'!=. & `vs'!=. & `lim2'

		egen rankpos1=  rank (posdis1) if `vs'>0 & `lim' & `ss' !=. & `vs'!=. & `lim1', unique
		egen rankpos2=  rank (posdis2) if `vs'>0 & `lim' & `ss' !=. & `vs'!=. & `lim2', unique

		gen  negdis1= -`vs' if `vs'<0 & `lim' & `ss'!=. &`vs'!=. & `lim1'
		gen  negdis2= -`vs' if `vs'<0 & `lim' & `ss'!=. &`vs'!=. & `lim2' 

		egen rankneg1 =  rank (negdis1) if `vs'<0 & `lim' & `ss'!=. & `vs'!=.  & `lim1', unique
		egen rankneg2 =  rank (negdis2) if `vs'<0 & `lim' & `ss'!=. & `vs'!=.  & `lim2', unique




		********define integrers to construct bins****************
		gen intpos1=  int(rankpos1/`bin') if `vs'>0 & `lim' &`ss'!=. &`vs'!=. & `lim1'
		gen intpos2=  int(rankpos2/`bin') if `vs'>0 & `lim' &`ss'!=. &`vs'!=. & `lim2'

		gen intneg1=  int(rankneg1/`bin') if `vs'<0 & `lim' &`ss'!=. &`vs'!=. & `lim1'
		gen intneg2=  int(rankneg2/`bin') if `vs'<0 & `lim' &`ss'!=. &`vs'!=. & `lim2'




		
		************create bin averages of both forcing variable and outcome for first subsample*************
		sort intpos1
		by intpos1: egen mean1`vs'pos=mean(`vs') if `vs'<`lim3'  & `lim1' & intpos1!=.
		by intpos1: egen mean1`ss'pos=mean(`ss') if `vs'<`lim3' &`ss'!=. & `lim1' & intpos1!=.

		sort intneg1
		by intneg1: egen mean1`vs'neg=mean(`vs') if `vs'> -`lim3' & `lim1' & intneg1!=.
		by intneg1: egen mean1`ss'neg=mean(`ss') if `vs'> - `lim3' &`ss'!=. & `lim1' & intneg1!=.


		************create bin averages of both forcing variable and outcome for second subsample*************
		sort intpos2
		by intpos2: egen mean2`vs'pos=mean(`vs') if `vs'<`lim3'  & `lim2' & intpos2!=.
		by intpos2: egen mean2`ss'pos=mean(`ss') if `vs'<`lim3' &`ss'!=. & `lim2' & intpos2!=.


		sort intneg2
		by intneg2: egen mean2`vs'neg=mean(`vs') if `vs'> -`lim3' & `lim2' & intneg2!=.
		by intneg2: egen mean2`ss'neg=mean(`ss') if `vs'> - `lim3' &`ss'!=. & `lim2' & intneg2!=.


		************create single binned variables of forcing variable and outcome for first subsample**************
		gen mean1`vs'=mean1`vs'pos if `vs'>0 & `lim'
		replace mean1`vs'=mean1`vs'neg if `vs'<0 & `lim'

		gen mean1`ss'=mean1`ss'pos if `vs'>0 & `lim'
		replace mean1`ss'=mean1`ss'neg if `vs'<0 & `lim'

		************create single binned variables of forcing variable and outcome for second subsample**************
		gen mean2`vs'=mean2`vs'pos if `vs'>0 & `lim'
		replace mean2`vs'=mean2`vs'neg if `vs'<0 & `lim'

		gen mean2`ss'=mean2`ss'pos if `vs'>0 & `lim'
		replace mean2`ss'=mean2`ss'neg if `vs'<0 & `lim'
		


twoway (scatter mean1`ss' mean1`vs' if abs(mean1`vs') < `lim3' , msize(medium)) (scatter mean2`ss' mean2`vs' if abs(mean2`vs' )< `lim3',msymbol(x) msize(large))  ,title(`grtit') subtitle(`grstit')  ytitle(`tit1') xtitle(`tit2') legend(on order(1 "`split1'" 2 "`split2'")) xline(0, lwidth(medthick) lcolor(black)) 

***************drop all generated variables*********************************
drop posdis1- mean2`ss'
end

***********************************************************
******************************************************
**********************************************************
******************************************************
/* define program to make split RD graph with seat share and vote share*/
/****************************************************/

capture program drop vsss2dis
program vsss2dis
	*******the arguments are the outcome variable, the forcing variable, restriction on the sample used, ***
	******* limit to define first subsample, condition used for second subsample***************
	*****restriction on the distance from the threshold, Y axis title, X axis title,**************
	*****definition of first subsample, definition och second subsample***********************
	args ss vs  bin lim lim1 lim2 lim3 tit1 tit2 split1 split2 grtit grstit

		********define integrers to construct bins****************
		gen int`vs'pos=  int(`vs'*`bin') if `vs'>0 & `lim'  &`vs'!=.
		gen int`vs'neg=  int(`vs'*`bin') if `vs'<0 & `lim' &`vs'!=.

		
		************create bin averages of both forcing variable and outcome for first subsample*************
		sort int`vs'pos
		gen mean1`vs'pos=(int`vs'pos +0.5)/`bin' if `vs'<`lim3'  & `lim1'
		by int`vs'pos: egen mean1`ss'pos=count(`ss') if `vs'<`lim3'  & `lim1'

		sort int`vs'neg
		gen mean1`vs'neg=(int`vs'neg-0.5)/`bin'  if `vs'> -`lim3' & `lim1'
		by int`vs'neg: egen mean1`ss'neg=count(`ss') if `vs'> - `lim3'  & `lim1'


		************create bin averages of both forcing variable and outcome for second subsample*************
		sort int`vs'pos
		gen mean2`vs'pos=(int`vs'pos +0.5)/`bin'   if `vs'<`lim3'  & `lim2'
		by int`vs'pos: egen mean2`ss'pos=count(`ss') if `vs'<`lim3'  & `lim2'


		sort int`vs'neg
		gen mean2`vs'neg=(int`vs'neg -0.5)/`bin' if `vs'> -`lim3' & `lim2'
		by int`vs'neg: egen mean2`ss'neg=count(`ss') if `vs'> - `lim3'  & `lim2'


		************create single binned variables of forcing variable and outcome for first subsample**************
		gen mean1`vs'=mean1`vs'pos if `vs'>0 & `lim'
		replace mean1`vs'=mean1`vs'neg if `vs'<0 & `lim'

		gen mean1`ss'=mean1`ss'pos if `vs'>0 & `lim'
		replace mean1`ss'=mean1`ss'neg if `vs'<0 & `lim'

		************create single binned variables of forcing variable and outcome for second subsample**************
		gen mean2`vs'=mean2`vs'pos if `vs'>0 & `lim'
		replace mean2`vs'=mean2`vs'neg if `vs'<0 & `lim'

		gen mean2`ss'=mean2`ss'pos if `vs'>0 & `lim'
		replace mean2`ss'=mean2`ss'neg if `vs'<0 & `lim'
		



twoway (scatter mean1`ss' mean1`vs' if mean1`vs' < `lim3' & mean1`vs' > -`lim3', msize(medium)) (scatter mean2`ss' mean2`vs' if mean2`vs' < `lim3' & mean2`vs' > -`lim3',msymbol(x) msize(large)) ,title(`grtit') subtitle(`grstit')  ytitle(`tit1') xtitle(`tit2') legend(on order(1 "`split1'" 2 "`split2'")) xline(50, lwidth(medthick) lcolor(black))

***************drop all generated variables*********************************
drop int`vs'pos - mean2`ss'

end

***********************************************************
********
***********************************************************
******************************************************
/* define program to make split RD graph with seat share and vote share*/
/****************************************************/

capture program drop vsss2sim
program vsss2sim
	*******the arguments are the outcome variable, the forcing variable, restriction on the sample used, ***
	******* limit to define first subsample, condition used for second subsample***************
	*****restriction on the distance from the threshold, Y axis title, X axis title,**************
	*****definition of first subsample, definition och second subsample***********************
	args ss vs  bin lim lim1 lim2 lim3 tit1 tit2 split1 split2 grtit grstit

		********define integrers to construct bins****************
		gen int`vs'pos=  int(`vs'*`bin') if `vs'>0 & `lim' &`ss'!=. &`vs'!=.
		gen int`vs'neg=  int(`vs'*`bin') if `vs'<0 & `lim' &`ss'!=. &`vs'!=.

		
		************create bin averages of both forcing variable and outcome for first subsample*************
		sort int`vs'pos
		gen mean1`vs'pos=(int`vs'pos +0.5)/`bin' if `vs'<`lim3'  & `lim1'
		by int`vs'pos: egen mean1`ss'pos=count(`ss') if `vs'<`lim3' &`ss'!=. & `lim1'

		sort int`vs'neg
		gen mean1`vs'neg=(int`vs'neg-0.5)/`bin'  if `vs'> -`lim3' & `lim1'
		by int`vs'neg: egen mean1`ss'neg=count(`ss') if `vs'> - `lim3' &`ss'!=. & `lim1'


		************create bin averages of both forcing variable and outcome for second subsample*************
		sort int`vs'pos
		gen mean2`vs'pos=(int`vs'pos +0.5)/`bin'   if `vs'<`lim3'  & `lim2'
		by int`vs'pos: egen mean2`ss'pos=count(`ss') if `vs'<`lim3' &`ss'!=. & `lim2'


		sort int`vs'neg
		gen mean2`vs'neg=(int`vs'neg -0.5)/`bin' if `vs'> -`lim3' & `lim2'
		by int`vs'neg: egen mean2`ss'neg=count(`ss') if `vs'> - `lim3' &`ss'!=. & `lim2'


		************create single binned variables of forcing variable and outcome for first subsample**************
		gen mean1`vs'=mean1`vs'pos if `vs'>0 & `lim'
		replace mean1`vs'=mean1`vs'neg if `vs'<0 & `lim'

		gen mean1`ss'=mean1`ss'pos if `vs'>0 & `lim'
		replace mean1`ss'=mean1`ss'neg if `vs'<0 & `lim'

		************create single binned variables of forcing variable and outcome for second subsample**************
		gen mean2`vs'=mean2`vs'pos if `vs'>0 & `lim'
		replace mean2`vs'=mean2`vs'neg if `vs'<0 & `lim'

		gen mean2`ss'=mean2`ss'pos if `vs'>0 & `lim'
		replace mean2`ss'=mean2`ss'neg if `vs'<0 & `lim'
		



twoway (scatter mean1`ss' mean1`vs' if mean1`vs' < `lim3' & mean1`vs' > -`lim3', msize(medium)) (scatter mean2`ss' mean2`vs' if mean2`vs' < `lim3' & mean2`vs' > -`lim3',msymbol(x) msize(large)) ,title(`grtit') subtitle(`grstit')  ytitle(`tit1') xtitle(`tit2') legend(on order(1 "`split1'" 2 "`split2'")) xline(0, lwidth(medthick) lcolor(black))

***************drop all generated variables*********************************
drop int`vs'pos - mean2`ss'

end

***********************************************************
********


