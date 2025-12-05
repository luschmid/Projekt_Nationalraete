use dta/CH71-03parties.dta, clear

** kt --> kanton identifier
** listid --> party identifer
keep kt listid votes S s year
gen listid2="A"
replace listid2="B" if listid==2
replace listid2="C" if listid==3
replace listid2="D" if listid==4
replace listid2="E" if listid==5
replace listid2="F" if listid==6
replace listid2="G" if listid==7
replace listid2="H" if listid==8
replace listid2="I" if listid==9
replace listid2="J" if listid==10
replace listid2="K" if listid==11
replace listid2="L" if listid==12
replace listid2="M" if listid==13
replace listid2="N" if listid==14
replace listid2="O" if listid==15
replace listid2="P" if listid==16
replace listid2="Q" if listid==17
replace listid2="R" if listid==18
replace listid2="S" if listid==27
replace listid2="T" if listid==35
rename s seats
drop listid
reshape wide votes seats, i(kt year) j(listid2) string

gen antalmandat= S

foreach parti in A B C D E F G H I J K L M N O P Q R S T {
gen r`parti'=votes`parti'
}

foreach sam in dh msl{
	foreach parti in A B C D E F G H I J K L M N O P Q R S T {
		gen `sam'm`parti ' =0
	}
}

gen rm =0

******dhondt*************
while rm<36 {

foreach parti in A B C D E F G H I J K L M N O P Q R S T {
gen j`parti '= r`parti'/(2+2*dhm`parti')
}
egen tm =rowtotal(dhmA -dhmT)
egen max = rowmax(jA-jT)
foreach parti in A B C D E F G H I J K L M N O P Q R S T {
replace dhm`parti'= dhm`parti' +1 if j`parti '==max & tm<antalmandat & r`parti'!=.
}
drop jA-max tm
replace rm=rm+1
}

******Modified Saint Lague*************
replace rm=0
while rm<36{

foreach parti in A B C D E F G H I J K L M N O P Q R S T {
gen p`parti'=mslm`parti'>=1 & mslm`parti'!=.
}
foreach parti in A B C D E F G H I J K L M N O P Q R S T {
gen j`parti'= (1-p`parti')*r`parti'/1.4 + p`parti'* r`parti'/(1+2*mslm`parti')
}
egen tm =rowtotal(mslmA -mslmT)
egen max = rowmax(jA-jT)
foreach parti in A B C D E F G H I J K L M N O P Q R S T {
replace mslm`parti'= mslm`parti' +1 if j`parti'==max & tm<antalmandat & r`parti'!=.
}
drop pA-max tm
replace rm=rm+1
}
replace rm=0

*********************************
***************************************************************************
*******************************************************************************

/* VOTES IF GOT PARTY GOT SEAT */

foreach sam in dh msl {
		foreach parti in A B C D E F G H I J K L M N O P Q R S T {
		gen `sam'rm`parti'=r`parti' if `sam'm`parti'>0
	}
}

egen dhtm =rowtotal(dhmA -dhmT )      /* number of reps with dh */
egen msltm =rowtotal(mslmA -mslmT )  /* number of reps with msl */
egen tr =rowtotal(rA -rT )    /* this is approvedvotes */

/* SUM OF VOTES FOR PARTIES THAT GOT SEAT */

foreach sam in dh msl {
	egen `sam'trm =rowtotal(`sam'rmA -`sam'rmT)
}

/* VOTESHARE */
foreach parti in A B C D E F G H I J K L M N O P Q R S T {
gen pr`parti'=r`parti'/tr
}

/* SEATSHARE */
foreach sam in dh msl {
	foreach parti in A B C D E F G H I J K L M N O P Q R S T {
	gen `sam'pm`parti'=`sam'm`parti'/`sam'tm
	}
}

**********Number of parties with seats********

gen P=0
foreach parti in A B C D E F G H I J K L M N O P Q R S T {
replace P= P+1 if dhm`parti'>0
}


/**************************************Main Program************************************************/

run do/Distance  /* PROGRAM FILES FOR COMPUTING DISTANCE TO THRESHOLD */

drop *msl*

calmindiff_dh "A B C D E F G H I J K L M N O P Q R S T" "A B C D E F G H I J K L M N O P Q R S T"

foreach parti in A B C D E F G H I J K L M N O P Q R S T {
replace mindiff`parti'n=. if  mindiff`parti'n <0
replace mindiff`parti'n=. if  mindiff`parti'n <0
}

keep kt pr* dhm* mindiff*p mindiff*n year antalmandat mindiff*p1 mindiff*n1

rename kt id
rename antalmandat m


/* FOLKE 1 , JEEA - MINIMUM DISTANCE */

foreach party in A B C D E F G H I J K L M N O P Q R S T {
gen min_distance_`party'=mindiff`party'p
replace min_distance_`party'=abs(mindiff`party'n) if abs(mindiff`party'n)<mindiff`party'p  /*closer to losing a seat than winning a seat */
replace min_distance_`party'=. if pr`party'==.    /* party not running */
}
egen min_distance=rowmin(min_distance_A-min_distance_T)

/* FOLKE ALTERNATIVE - MINIMUM DISTANCE IN OWN VOTES */

foreach party in A B C D E F G H I J K L M N O P Q R S T {
gen min_distance1_`party'=mindiff`party'p1
replace min_distance1_`party'=abs(mindiff`party'n1) if abs(mindiff`party'n1)<mindiff`party'p1  /*closer to losing a seat than winning a seat */
replace min_distance1_`party'=. if pr`party'==.    /* party not running */
}
egen min_distance1=rowmin(min_distance1_A-min_distance1_T)


**************************************
/* voteshare weighted measures min_distance */
foreach party in A B C D E F G H I J K L M N O P Q R S T {
gen min_distance_`party'Xpr`party'=min_distance_`party'*pr`party'
}

egen min_distance_iXv_i=rowtotal(min_distance_AXprA-min_distance_TXprT)

/* voteshare weighted measures min_distance1 */
foreach party in A B C D E F G H I J K L M N O P Q R S T {
gen min_distance1_`party'Xpr`party'=min_distance1_`party'*pr`party'
}

egen min_distance1_iXv_i=rowtotal(min_distance1_AXprA-min_distance1_TXprT)

sort year id
save dta/Switzerland_CFS, replace
