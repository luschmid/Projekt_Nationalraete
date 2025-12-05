/* Generate ksmall.dta, which reduces the panel from {candidate,year} to {candidate,election-year} 
and removes unnecessary variables
*/

capture drop _m
sort id year

/* deflating income variables to 2011 Euros */
replace ansiotulo=ansiotulon*kerroin11
replace paaomatulo=paaomatulo*kerroin11

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
gen temp2=L.past_avg
replace temp2=L2.past_avg if year==1995 /* year 1994 missing */
replace past_avg=temp2
drop temp2
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
*replace tempcount=. if cumul_income==.
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
gen temp2=L.past_avgCap
replace temp2=L2.past_avgCap if year==1995 /* year 1994 missing */
replace past_avgCap=temp2
drop temp2
gen future_avgCap=(final_cumulCap-cumul_capincome)/(nyearsCap-tempcount)

gsort id -year
replace future_avgCap=future_avgCap[_n-1] if (id[_n]==id[_n-1]) & future_avgCap==.
sort id year
replace past_avgCap=past_avgCap[_n-1] if (id[_n]==id[_n-1]) & past_avgCap==.

quietly compress
tempfile tempepanel
save `tempepanel',replace

drop tempcount


* CREATE LAG and FWD VARIABLES

cap drop next* last4* last1*

forvalues s = 1/15 {
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


forvalues s = 1/15 {
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

*
* GENERATE OTHER VARIABLES

* Year dummies
gen byte y08=(year==2008)
gen byte y04=(year==2004)
gen byte y00=(year==2000)
gen byte y96=(year==1996)

encode party,gen(nparty)
label variable nparty "Party"

* Number of seats and votes in a given election (district-year)

sort year municipality
cap drop nvotes 
by year municipality: egen nvotes=total(votes) 
by year municipality: egen nseats=total(elected)
cap drop voteshare
gen voteshare=votes/nvotes

replace incumbent=0 if incumbent==. & year==2008

* Generate variables for estimating the incumbency effect
* taking into account uneven spacing of 1970s elections

sort id year
gen nextElected=F4.elected

gen bsmargin=bsmargina
gen flip= ((elected & bsmargin<0) | (!elected & bsmargin>0))
replace bsmargin = -bsmargin if flip		/* incorrect ordering by bsmargin flipped (due to finite number of simulated elections) */

* Reduce panel 
keep if elected!=.
keep id year y08 y04 y00 y96 elected nextElected avgtulo avgpotulo /*ansiotulo paaomatulo */ ///
bsmargin votes voteshare nvotes nseats incumbent  KESK KOK SDP OTHER nparty eldist female yob age ///
  next* last* nyears* past_avg future_avg past* future*

drop lastname lastyear
 
label variable female "Female"
label variable age "Age"
label var yob "Year of Birth"
label variable incumbent "Incumbent"
label variable voteshare "Vote share"
label variable avgtulo "Earnings 1993-2011"
label variable avgpotulo "Capital income 1993-2011"
label variable next1_3tulo "Earnings in e=1"
label variable next5_7tulo "Earnings in e=2"
label variable next9_11tulo "Earnings in e=3"
label variable last1_3tulo "Earnings in e=-1"
label variable future_avg "Future average yearly earnings"
label variable past_avg "Past average yearly earnings"
label variable KESK "Centre"
label variable KOK "NCP"
label variable SDP "SDP"
label variable OTHER "Other parties"
label variable nseats "District size (seats)"
label variable nvotes "District size (votes)"
label variable eldist "Municipality code"
*label variable ansiotulo "Earnings (in 2011 Euros)"
*label variable paaomatulo "Capital income (in 2011 Euros)"
rename bsmargin x
label var x "Closeness of election"
replace x = sign(elected-0.5)*0.0005 if x==0 /* exact ties would trigger fuzzy RD in rd.ado*/

compress
sort id year
save $DATAPATH/ksmall,replace
*outsheet using ksmall.csv,comma replace
