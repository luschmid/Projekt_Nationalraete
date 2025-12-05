******************* DATA FROM SPAIN ***************************
use dta/Spain_GS.dta, clear 
sort year id
merge year id using dta/Spain_CFS.dta
keep year id m turnout min_distance min_distance1 C1 C2
gen country="SPA"
rename id districtname
encode districtname, gen(temp)
gen id=temp+100
save dta/Spain.dta, replace

******************* DATA FROM SWITZERLAND **********************
use dta/Switzerland_GS.dta, clear
sort year id
merge year id using dta/Switzerland_CFS.dta
keep year id m turnout min_distance min_distance1 C1 C2
gen country="SWI"
drop if id==14   /* GS (2009): "For the same reason. we have excluded one district where voting is compulsory (Schafhausen) since compulsory voting presumably distorts the competition-turnout nexus" */
save dta/Switzerland.dta,replace

******************* APPEND THE TWO DATA FILES********************
append using dta/Spain.dta, force

***** traditional measure and BL measure
gen margin=min_distance1 
gen marginBL=margin*m  /* Blais-Lago normalize with district magnitude */

save dta/Panel.dta,replace
