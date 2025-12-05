use dta/ESP77-05parties.dta , clear

** district --> district identifier
** listid --> party identifer
*keep if lower!=.
rename partyvotes votes
egen min_votes=min(votes), by(year district listid)
egen max_votes=max(votes), by(year district listid)
assert min_votes==max_votes

keep district listid votes S s year
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
replace listid2="S" if listid==19
replace listid2="T" if listid==20
replace listid2="U" if listid==21
replace listid2="V" if listid==22
replace listid2="W" if listid==23
replace listid2="X" if listid==24
replace listid2="Y" if listid==25
replace listid2="Z" if listid==26
replace listid2="Æ" if listid==27
replace listid2="Ø" if listid==28
replace listid2="Å" if listid==29

rename s seats
tab listid listid2
drop listid
reshape wide votes seats, i(district year) j(listid2) string

gen antalmandat= S

foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
gen r`parti'=votes`parti'
}

foreach sam in dh msl{
	foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
		gen `sam'm`parti ' =0
	}
}

gen rm =0

******dhondt*************
while rm<36 {

foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
gen j`parti '= r`parti'/(2+2*dhm`parti')
}
egen tm =rowtotal(dhmA -dhmÅ)
egen max = rowmax(jA-jÅ)
foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
replace dhm`parti'= dhm`parti' +1 if j`parti '==max & tm<antalmandat & r`parti'!=.
}
drop jA-max tm
replace rm=rm+1
}

******Modified Saint Lague*************
replace rm=0
while rm<36{

foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
gen p`parti'=mslm`parti'>=1 & mslm`parti'!=.
}
foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
gen j`parti'= (1-p`parti')*r`parti'/1.4 + p`parti'* r`parti'/(1+2*mslm`parti')
}
egen tm =rowtotal(mslmA -mslmÅ)
egen max = rowmax(jA-jÅ)
foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
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
		foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
		gen `sam'rm`parti'=r`parti' if `sam'm`parti'>0
	}
}

egen dhtm =rowtotal(dhmA -dhmÅ )      /* number of reps with dh */
egen msltm =rowtotal(mslmA -mslmÅ )  /* number of reps with msl */
egen tr =rowtotal(rA -rÅ )    /* this is approvedvotes */

/* SUM OF VOTES FOR PARTIES THAT GOT SEAT */

foreach sam in dh msl {
	egen `sam'trm =rowtotal(`sam'rmA -`sam'rmÅ)
}

/* VOTESHARE */
foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
gen pr`parti'=r`parti'/tr
}

/* SEATSHARE */
foreach sam in dh msl {
	foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
	gen `sam'pm`parti'=`sam'm`parti'/`sam'tm
	}
}

**********Number of parties with seats********

gen P=0
foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
replace P= P+1 if dhm`parti'>0
}


/**************************************Main Program************************************************/

run do/Distance  /* PROGRAM FILES FOR COMPUTING DISTANCE TO THRESHOLD */

drop *msl*

calmindiff_dh "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å" "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å"

foreach parti in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
replace mindiff`parti'n=. if  mindiff`parti'n <0
replace mindiff`parti'n=. if  mindiff`parti'n <0
}

keep district pr* dhm* mindiff*p mindiff*n year antalmandat mindiff*p1 mindiff*n1

rename district id
rename antalmandat m

/* FOLKE 1 , JEEA - MINIMUM DISTANCE */

foreach party in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
gen min_distance_`party'=mindiff`party'p
replace min_distance_`party'=abs(mindiff`party'n) if abs(mindiff`party'n)<mindiff`party'p  /*closer to losing a seat than winning a seat */
replace min_distance_`party'=. if pr`party'==.    /* party not running */
replace min_distance_`party'=. if pr`party'==1    /* 100% votes to one party */
}
egen min_distance=rowmin(min_distance_A-min_distance_Å)
sum min_distance
*hist min_distance

/* FOLKE ALTERNATIVE - MINIMUM DISTANCE IN OWN VOTES */

foreach party in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z Æ Ø Å {
gen min_distance1_`party'=mindiff`party'p1
replace min_distance1_`party'=abs(mindiff`party'n1) if abs(mindiff`party'n1)<mindiff`party'p1  /*closer to losing a seat than winning a seat */
replace min_distance1_`party'=. if pr`party'==.    /* party not running */
replace min_distance1_`party'=. if pr`party'==1    /* 100% votes to one party */
}
egen min_distance1=rowmin(min_distance1_A-min_distance1_Å)

sort year id
save dta/Spain_CFS, replace
