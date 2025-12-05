/* Generate esmall.dta, which reduces epanel.dtal from {candidate,year} to {candidate,election-year} 
   This also calculates average income variables and removes unnecessary variables */

*use $DATAPATH/epanel

/* Get PATH DEFINITIONS from main.do */

cap drop _m
sort id year

/* deflating income variables to 2011 Euros */
replace paaomatulo=paaomatulo*kerroin11
replace ansiotulo=ansiotulon*kerroin11

/* missing incomes coded as zeros in the tax register data for 1993*/
replace ansiotulo=. if ansiotulo==0 & year==1993
replace paaomatulo=. if paaomatulo==0 & year==1993


/* Calculate cumulative income */
gen temptulo=ansiotulo
replace temptulo=0 if ansiotulo==.
gen cumul_income=temptulo[1]
replace cumul_income=temptulo[_n] if id[_n]!=id[_n-1] 
replace cumul_income=temptulo[_n] + cumul_income[_n-1] if id[_n]==id[_n-1] 
replace cumul_income=. if ansiotulo==.
drop temptulo

* Calculate cumulative number of years with observed income 
gen tempcount= ansiotulo!=. 
replace tempcount=tempcount[_n] + tempcount[_n-1] if id[_n]==id[_n-1] 
*replace tempcount=. if cumul_income==.
replace tempcount=. if tempcount==0


* Calculate average of all observations
sort id year
by id: egen avgtulo=mean(ansiotulo) 


quietly compress
tempfile tempepanel
save `tempepanel',replace


/* Calculate future/past averages for all income observations */
gen tempyear=year
*replace tempyear=. if cumul_income==.
replace tempyear=. if ansiotulo==. 
collapse (min) earliest_year=tempyear (max) final_year=tempyear final_cumul_income=cumul_income nyearsTulo=tempcount,by(id)
replace nyearsTulo=0 if nyearsTulo==.
quietly compress
sort id
merge 1:m id using `tempepanel'
drop _merge
erase `tempepanel'
sort id year
replace tempcount=tempcount[_n-1] if (id[_n]==id[_n-1]) & tempcount==.
gen past_avg=ansiotulo
replace past_avg=(cumul_income)/(tempcount) if tempcount>1
gen temp=L.past_avg
replace temp=L2.past_avg if year==1995 /* year 1994 missing */
replace past_avg=temp
drop temp
gen future_avg=(final_cumul-cumul_income)/(nyearsTulo-tempcount)

quietly compress
tempfile tempepanel
save `tempepanel',replace


 
gsort id -year
replace future_avg=future_avg[_n-1] if (id[_n]==id[_n-1]) & future_avg==.
sort id year
replace past_avg=past_avg[_n-1] if (id[_n]==id[_n-1]) & past_avg==.

drop tempcount 


* Averages for capital income

/* Calculate cumulative capital income */
sort id year
gen temptulo=paaomatulo
replace temptulo=0 if temptulo==.
gen cumul_capincome=temptulo[1]
replace cumul_capincome=temptulo[_n] if id[_n]!=id[_n-1] 
replace cumul_capincome=temptulo[_n] + cumul_capincome[_n-1] if id[_n]==id[_n-1] 
replace cumul_capincome=. if paaomatulo==.
drop temptulo

* Calculate cumulative number of years with observed capital income 
gen tempcount= paaomatulo!=. 
replace tempcount=tempcount[_n] + tempcount[_n-1] if id[_n]==id[_n-1] 
replace tempcount=. if tempcount==0

* Calculate average of all observations

by id: egen avgpotulo=mean(paaomatulo) 

quietly compress
tempfile tempepanel
save `tempepanel',replace

/* Calculate future/past averages for all capital income observations */
gen tempyear=year
replace tempyear=. if paaomatulo==. 
collapse (min) earliest_yearCap=tempyear (max) final_yearCap=tempyear final_cumulCap=cumul_capincome nyearsCap=tempcount,by(id)
replace nyearsCap=0 if nyearsCap==.
quietly compress
sort id
merge 1:m id using `tempepanel'
drop _merge
erase `tempepanel'
sort id year
replace tempcount=tempcount[_n-1] if (id[_n]==id[_n-1]) & tempcount==.
gen past_avgCap=paaomatulo
replace past_avgCap=(cumul_capincome)/(tempcount) if tempcount>1
gen temp=L.past_avgCap
replace temp=L2.past_avgCap if year==1995 /* year 1994 missing */
replace past_avgCap=temp
drop temp
gen future_avgCap=(final_cumulCap-cumul_capincome)/(nyearsCap-tempcount)

gsort id -year
replace future_avgCap=future_avgCap[_n-1] if (id[_n]==id[_n-1]) & future_avgCap==.
sort id year
replace past_avgCap=past_avgCap[_n-1] if (id[_n]==id[_n-1]) & past_avgCap==.

quietly compress
tempfile tempepanel
save `tempepanel',replace

drop tempcount


* Generate average incomes for years before and after 2000
* Collapse and merge back
quietly compress
tempfile tempepanel2
cap drop _m
save `tempepanel2',replace

gen int cumul_obs=0
replace cumul_obs = 1 if (id[_n]!=id[_n-1]) & ansiotulo != .
replace cumul_obs = cumul_obs[_n-1]+1 if (id[_n]==id[_n-1]) & ansiotulo != .
gen tuloPre2k=cumul_inc/cumul_obs if year<2000
gen tuloPost2k=(final_cumul_inc-cumul_in)/(nyearsTulo-cumul_obs) if year>2000

collapse (lastnm) tuloPre2k (firstnm) tuloPost2k, by(id)
sort id
merge 1:m id using `tempepanel2'
erase `tempepanel2'
sort id year


* Generate year dummies
forvalues iy = 1975(4)2007 {
local ystring=substr("`iy'",-2,2)
gen byte y`ystring'=(year==`iy')
}
gen byte y72=(year==1972)
gen byte y70=(year==1970)


* CREATE LAG and FWD VARIABLES (rd does not accept L/F)

/* leads and lags of earnings  */
drop next* last4* last1*

forvalues s = 1/40 {
local j=`s'+1
local k=`s'+2
gen next`s'_`k'tulo=(F`s'.ansiotulo+F`j'.ansiotulo+F`k'.ansiotulo)/3
replace next`s'_`k'tulo=(F`s'.ansiotulo+F`j'.ansiotulo)/2 if next`s'_`k'tulo==.
replace next`s'_`k'tulo=(F`s'.ansiotulo+F`k'.ansiotulo)/2 if next`s'_`k'tulo==.
replace next`s'_`k'tulo=(F`j'.ansiotulo+F`k'.ansiotulo)/2 if next`s'_`k'tulo==.
replace next`s'_`k'tulo=F`s'.ansiotulo if next`s'_`k'tulo==.
replace next`s'_`k'tulo=F`j'.ansiotulo if next`s'_`k'tulo==.
replace next`s'_`k'tulo=F`k'.ansiotulo if next`s'_`k'tulo==.
}

forvalues s=1/5 {
local j=`s'+1
local k=`s'+2
gen last`s'_`k'tulo=(L`s'.ansiotulo+L`j'.ansiotulo+L`k'.ansiotulo)/3
replace last`s'_`k'tulo=(L`s'.ansiotulo+L`j'.ansiotulo)/2 if last`s'_`k'tulo==.
replace last`s'_`k'tulo=(L`s'.ansiotulo+L`k'.ansiotulo)/2 if last`s'_`k'tulo==.
replace last`s'_`k'tulo=(L`j'.ansiotulo+L`k'.ansiotulo)/2 if last`s'_`k'tulo==.
replace last`s'_`k'tulo=L`s'.ansiotulo if last`s'_`k'tulo==.
replace last`s'_`k'tulo=L`j'.ansiotulo if last`s'_`k'tulo==.
replace last`s'_`k'tulo=L`k'.ansiotulo if last`s'_`k'tulo==.
}



gen next1_3tulopost=(F2.ansiotulo+F3.ansiotulo)/2 if year==1999
replace next1_3tulopost=F2.ansiotulo if next1_3tulopost==. & year==1999
replace next1_3tulopost=F3.ansiotulo if next1_3tulopost==. & year==1999


/*leads and lags of capital income*/

forvalues s = 1/25 {
local j=`s'+1
local k=`s'+2
gen next`s'_`k'potulo=(F`s'.paaomatulo+F`j'.paaomatulo+F`k'.paaomatulo)/3
replace next`s'_`k'potulo=(F`s'.paaomatulo+F`j'.paaomatulo)/2 if next`s'_`k'potulo==.
replace next`s'_`k'potulo=(F`s'.paaomatulo+F`k'.paaomatulo)/2 if next`s'_`k'potulo==.
replace next`s'_`k'potulo=(F`j'.paaomatulo+F`k'.paaomatulo)/2 if next`s'_`k'potulo==.
replace next`s'_`k'potulo=F`s'.paaomatulo if next`s'_`k'potulo==.
replace next`s'_`k'potulo=F`j'.paaomatulo if next`s'_`k'potulo==.
replace next`s'_`k'potulo=F`k'.paaomatulo if next`s'_`k'potulo==.
}

forvalues s=1/5 {
local j=`s'+1
local k=`s'+2
gen last`s'_`k'potulo=(L`s'.paaomatulo+L`j'.paaomatulo+L`k'.paaomatulo)/3
replace last`s'_`k'potulo=(L`s'.paaomatulo+L`j'.paaomatulo)/2 if last`s'_`k'potulo==.
replace last`s'_`k'potulo=(L`s'.paaomatulo+L`k'.paaomatulo)/2 if last`s'_`k'potulo==.
replace last`s'_`k'potulo=(L`j'.paaomatulo+L`k'.paaomatulo)/2 if last`s'_`k'potulo==.
replace last`s'_`k'potulo=L`s'.paaomatulo if last`s'_`k'potulo==.
replace last`s'_`k'potulo=L`j'.paaomatulo if last`s'_`k'potulo==.
replace last`s'_`k'potulo=L`k'.paaomatulo if last`s'_`k'potulo==.
}


* Generate other variables

encode party, gen(nparty)
label variable nparty "Party"

gen byte south=1 if eldist<8 | eldist==13 | eldist==12
replace south=0 if south==.

sort year eldist
drop nvotes
by year eldist: egen nvotes=total(votes)
by year eldist: egen nseats=total(elected)
cap replace voteshare=votes/nvotes


* Generate variables for estimating the incumbency effect
* taking into account uneven spacing of 1970s elections

sort id year

gen nextElected=F4.elected
replace nextElected=F2.elected if year==1970
replace nextElected=F3.elected if year==1972

gen nextElected2=F8.elected
replace nextElected2=F5.elected if year==1970
replace nextElected2=F7.elected if year==1972

gen nextElected3=F12.elected
replace nextElected3=F9.elected if year==1970
replace nextElected3=F11.elected if year==1972

gen nextElected4=F16.elected
replace nextElected4=F13.elected if year==1970
replace nextElected4=F15.elected if year==1972

gen nextElected5=F20.elected
replace nextElected5=F17.elected if year==1970
replace nextElected5=F19.elected if year==1972

gen nextElected6=F24.elected
replace nextElected6=F21.elected if year==1970
replace nextElected6=F23.elected if year==1972

gen lastElected=L4.elected
replace lastElected=L2.elected if year==1972
replace lastElected=L3.elected if year==1975
replace lastElected=elected66 if year==1970

gen lastElected2=L8.elected
replace lastElected2=elected66 if year==1972
replace lastElected2=L5.elected if year==1975
replace lastElected2=elected62 if year==1970



replace nextElected=0 if nextElected==. & year<2007
replace nextElected2=0 if nextElected2==. & year<2003
replace nextElected3=0 if nextElected3==. & year<1999
replace nextElected4=0 if nextElected4==. & year<1995
replace nextElected5=0 if nextElected5==. & year<1991
replace nextElected6=0 if nextElected6==. & year<1987
replace lastElected=0 if lastElected==.
replace lastElected2=0 if lastElected==. & year>1970

gen alwaysIn2=1 if elected==1 & nextElected==1 & nextElected2==1
replace alwaysIn2=0 if alwaysIn2==. & year<2003
gen alwaysIn3=1 if elected==1 & nextElected==1 & nextElected2==1 & nextElected3==1
replace alwaysIn3=0 if alwaysIn3==. & year<1999
gen alwaysIn4=1 if elected==1 & nextElected==1 & nextElected2==1 & nextElected3==1 & nextElected4==1
replace alwaysIn4=0 if alwaysIn4==. & year<1995
gen alwaysIn5=1 if elected==1 & nextElected==1 & nextElected2==1 & nextElected3==1 & nextElected4==1 & nextElected5==1
replace alwaysIn5=0 if alwaysIn5==. & year<1991
gen alwaysIn6=1 if elected==1 & nextElected==1 & nextElected2==1 & nextElected3==1 & nextElected4==1 & nextElected5==1 & nextElected6==1
replace alwaysIn6=0 if alwaysIn6==. & year<1987

gen neverIn2=1 if elected==0 & nextElected==0 & nextElected2==0
replace neverIn2=0 if neverIn2==. & year<2003 
gen neverIn3=1 if elected==0 & nextElected==0 & nextElected2==0 & nextElected3==0
replace neverIn3=0 if neverIn3==. & year<1999
gen neverIn4=1 if elected==0 & nextElected==0 & nextElected2==0 & nextElected3==0 & nextElected4==0
replace neverIn4=0 if neverIn4==. & year<1995
gen neverIn5=1 if elected==0 & nextElected==0 & nextElected2==0 & nextElected3==0 & nextElected4==0 & nextElected5==0
replace neverIn5=0 if neverIn5==. & year<1991
gen neverIn6=1 if elected==0 & nextElected==0 & nextElected2==0 & nextElected3==0 & nextElected4==0 & nextElected5==0 & nextElected6==0
replace neverIn6=0 if neverIn6==. & year<1987



* Indicator variable for whether a candidate was elected in some previous election

 
* create dummy for never-before-elected
sort id year
foreach thisvar in elected66 elected62 elected58 elected54 elected51 electedpre50 {
replace `thisvar'=0 if `thisvar'==.
}
gen elBefore70=(elected66+elected62+elected58+elected54+elected51+electedpre50)
replace elBefore70=0 if elBefore70==.
gen electedB=elected 
replace electedB=0 if elected==.
gen cumul_el=electedB[1]+elBefore70[1]
replace cumul_el = electedB[_n] + elBefore70[_n] if id[_n]!=id[_n-1] 
replace cumul_el = electedB[_n] + cumul_el[_n-1] if id[_n]==id[_n-1] 
drop elBefore70 electedB 
gen neverelectedbefore=(cumul_el==0 |(cumul_el==1 & elected==1))


* Reduce panel from years to electoral periods
keep if elected!=.

keep id elected  bslevel bsmargin nelected nruns votes voteshare nvotes nseats   ///
 incumbent alliance lastname firstname KESK KOK SDP OTHER nparty eldist south female yob age missing ///
 year y07 y03 y99 y95 y91 y87 y83 y79 y75 y72 y70 ///
 nyears* past* future* next* last* never* always* avgtulo avgpotulo tuloP*
/* /// tulo0811 tulo0411 tulo0011 tulo9611 paaomatulo ansiotulo fwdtulo* lagtulo* tempcount pivot */ 

drop lastyear

label variable future_avg "Future average yearly earnings"
label variable past_avg "Past average yearly earnings"
label variable avgtulo "Earnings 1993-2011"
label variable avgpotulo "Capital income 1993-2011"
label variable next1_3tulo "Average yearly earnings in (t+1) to (t+3)"
label variable last1_3tulo "Average yearly earnings in (t-1) to (t-3)"
label variable tuloPre2k "Average yearly earnings before 2000"
label variable tuloPost2k "Average yearly earnings after 2000"
*label variable ansiotulo "Earnings"
*label variable paaomatulo "Capital income"
label variable female "Female"
label variable age "Age"
label variable incumbent "Incumbent"
label variable next1_3tulo "Earnings e=1"
label variable next5_7tulo "Earnings e=2"
label variable next9_11tulo "Earnings e=3"
label variable last1_3tulo "Earnings e=-1"
label variable KESK "Centre"
label variable KOK "NCP"
label variable SDP "SDP"
label variable OTHER "Other parties"
label variable voteshare "Vote share"
label variable yob "Year of Birth"
label variable south "Southern Finland"
label variable year "Year"
label variable eldist "Electoral district"
label variable elected "Elected"
label variable votes "Votes"
label variable neverelected "Never elected"
label variable nseats "District size (seats)"
label variable nvotes "District size (votes)"

/* Forcing variable */
rename bsmargin x
label var x "Closeness of election"
replace x =sign(elected-0.5)*0.0005 if x==0 /* exact ties would trigger fuzzy RD in rd.ado*/

*outsheet using $DATAPATH/esmall.csv,comma replace 
drop firstname lastname
compress
sort id year
save $DATAPATH/esmall, replace
