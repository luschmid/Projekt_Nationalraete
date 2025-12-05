capture program drop coalprep
program coalprep
	args namelist coallist

use ..\dta\MandatfordelningOlikaMetoderBeforeReshape, clear

**************************** KASTER UT KOMMUNER MED AVVIK **********************************

sort knr ele
merge knr ele using ..\dta\Index
drop _merge
drop if avvikelse_kommune>0

********* RATIOEN MELLOM "TOTAL LVOTES" OG "TOTAL VOTES" SKAL VÆRE EKSAKT COUNCIL SIZE *****

sort knr ele
merge knr ele using ..\dta\WrongRatio
drop if WrongRatio==1 

********************************************************************************************


foreach parti in `namelist'{
	
	gen `parti'_s=  Rep`parti'
}

foreach parti in `namelist'{
	gen `parti'_vs=LVoteShare`parti'
}

foreach parti in `namelist'{
	gen `parti'_ss=SeatShare`parti'
}

keep knr electionperiod RV_s- JointR2_ss SizeOfCouncil


foreach parti in `namelist'{
	gen `parti'_coal=0
}
	foreach parti in `coallist'{
		replace `parti'_coal=1
}


gen ss_coal=0
gen vs_coal=0

foreach parti in `coallist'{
replace ss_coal=ss_coal+`parti'_ss
	replace vs_coal=vs_coal+`parti'_vs
	}
	
*	drop if ss_maj<.5

gen maj_coal=ss_coal>.5 if ss_coal!=.

keep if ele>8 /* MSL */

gen count=1
egen obs =sum(count)
gen tot_s=SizeOfCouncil


save "../dta/hel kommun grunddata coal.dta", replace

drop if knr==439 & ele==10  /* SV , SP OG DNA ER SAMTLIGE LISTER */
drop if knr==441 & ele==10  /* SV , SP OG DNA ER SAMTLIGE LISTER */
drop if knr==1755 & ele==10  /* SV , SP OG DNA ER SAMTLIGE LISTER */
drop if knr==2023 & ele==9  /* SV , SP OG DNA ER SAMTLIGE LISTER */
drop if knr==615 & ele==10  /* høyreblokk rent flertall */
drop if knr==1144 & ele==10  /* various seatshare = 1 */
drop if knr==1151 & ele==10  /* various seatshare = 1 */
drop if knr==1856 & ele==10  /* various seatshare = 1 */
drop if knr==1915 & ele==9  /* various seatshare = 1 */
drop if knr==1915 & ele==10  /* various seatshare = 1 */
drop if knr==1711 & ele==10  /* stopper på denne hvis ikke droppes. jeg skjønner ikke hvorfor - må øke antall sim? */
drop if knr==1835 & ele==9  /* stopper på denne hvis ikke droppes. jeg skjønner ikke hvorfor - må øke antall sim? */
drop if knr==1433 & ele==9  /* stopper på denne hvis ikke droppes. jeg skjønner ikke hvorfor - må øke antall sim? */
drop if knr==1859 & ele==10  /* stopper på denne hvis ikke droppes. jeg skjønner ikke hvorfor - må øke antall sim? */
drop if knr==1133 & ele==10  /* stopper på denne hvis ikke droppes. jeg skjønner ikke hvorfor */
drop if knr==1818 & ele==9  /* stopper på denne hvis ikke droppes. jeg skjønner ikke hvorfor */
drop if knr==1822 & ele==9  /* stopper på denne hvis ikke droppes. jeg skjønner ikke hvorfor */
drop if knr==1917 & ele==9  /* stopper på denne hvis ikke droppes. jeg skjønner ikke hvorfor */


save "../dta/hel kommun grunddata coal SP.dta", replace


end
