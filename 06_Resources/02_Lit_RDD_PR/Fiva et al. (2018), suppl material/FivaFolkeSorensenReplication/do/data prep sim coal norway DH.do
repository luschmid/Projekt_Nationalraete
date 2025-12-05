capture program drop coalprepDH
program coalprepDH
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

keep if ele==8 /* DH */

gen count=1
egen obs =sum(count)
gen tot_s=SizeOfCouncil


save "../dta/hel kommun grunddata coal DH.dta", replace

end
