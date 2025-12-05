

global path "D:\SchmidLu\Dropbox\Projekt Nationalräte\1 Data\1_data\14_jarckspuehlerwidmer"


// This file creates a file with all observations from TI for Jonathan Sartori

capture erase "$path\Outfiles3\help.dta"

cd "$path\Outfiles3"

forvalues i=1(1)10{
di "`file'"
import excel  using GeburtsTodesdaten3_C_`i',  clear first case(lower)
capture append using "$path\Outfiles3\help.dta", force
save "$path\Outfiles3\help.dta", replace
import excel  using GeburtsTodesdaten3_D_`i',  clear first case(lower)
append using "$path\Infiles\help.dta", force
save "$path\Outfiles3\help.dta", replace
}

keep if canton=="TI"
keep id canton sex name firstname birthyear geburtsdatumneu bday bmonth byear quellegeburtsdatumneu todesdatumneu dday dmonth dyear quelletodesdatumneu job partyname municipality origin

egen nononmiss = rownonmiss(bday bmonth byear dday dmonth dyear)
keep if nononmiss==0

duplicates drop id,force

export excel  using "$path\Outfiles3\GeburtsTodesdaten3_CD_TI.xlsx",  replace   firstrow(varlabels)




