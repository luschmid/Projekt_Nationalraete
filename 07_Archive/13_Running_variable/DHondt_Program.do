******************************************
* C) Program for Seat Allocation D'Hondt *
******************************************

clear programs

timer on 1

program DHondt

quietly {
save "$path\02_Processed_data\13_Running_variable\example_input.dta", replace

unique party_num if alliance_num==`1'  	// l = alliance
local Al_size=r(unique)
di "Al_size:" "`Al_size'"

unique party_num if suballiance_num==`2'  // s = suballiance
local Sal_l_size=r(unique)
di "Sal_size:" "`Sal_l_size'"


* (i) Allocation of seats to alliances

collapse (first) votes_j alliance_num, by(party_num)
collapse (sum) votes_j (first) party_num, by(alliance_num)

replace party_num=`3' if alliance_num==`1'

forvalues k = 1(1)$n {
	generate quot`k' = votes_j/`k'
}

reshape long quot, i(alliance_num) j(number)

egen quot_rank=rank(-quot), unique

gen seat=0 
replace seat=1 if quot_rank<=$n


* (ii) Allocation of alliance seats to suballiances/parties

if `Al_size'>1 & `Al_size'!=. {

keep if alliance_num ==`1'
collapse (sum) seats_l=seat, by(alliance_num)
save "$path\02_Processed_data\13_Running_variable\SeatsAlliance.dta", replace

use "$path\02_Processed_data\13_Running_variable\example_input.dta", clear
collapse (first) votes_j alliance_num suballiance_num, by(party_num)
collapse (sum) votes_j (first) alliance_num party_num, by(suballiance_num)
replace party_num=`3' if suballiance_num==`2'
merge m:1 alliance_num using "$path\02_Processed_data\13_Running_variable\SeatsAlliance.dta"
keep if _merge==3
		
forvalues k = 1(1)$n {
	generate quot`k' = votes_j/`k'
}

reshape long quot, i(suballiance_num) j(number) 
egen quot_rank=rank(-quot), unique 
	
gen seat=0 
replace seat=1 if quot_rank<=seats_l
keep suballiance_num party_num number votes_j seat alliance_num 
}


* (iii) Allocation of suballiance seats to parties

if `Sal_l_size'>1 & `Sal_l_size'!=. {
keep if suballiance_num ==`2'
collapse (sum) seats_s=seat, by(suballiance_num)
save "$path\02_Processed_data\13_Running_variable\SeatsSuballiance.dta", replace

use "$path\02_Processed_data\13_Running_variable\example_input.dta", clear
collapse (first) votes_j (first) suballiance_num, by(party_num)
merge m:1 suballiance_num using "$path\02_Processed_data\13_Running_variable\SeatsSuballiance.dta"
keep if _merge==3
		
forvalues k = 1(1)$n {
	generate quot`k' = votes_j/`k'
}

reshape long quot, i(party_num) j(number)
	
egen quot_rank=rank(-quot), unique
	
gen seat=0 
replace seat=1 if quot_rank<=seats_s

keep party_num number votes_j seat  

}
}

else{
}
end
