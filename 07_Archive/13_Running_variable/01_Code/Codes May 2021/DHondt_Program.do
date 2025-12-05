******************************************
* C) Program for Seat Allocation D'Hondt *
******************************************

clear programs

timer on 1

program DHondt

quietly {
save "$path\02_Processed_data\13_Running_variable\example_input_$canton.dta", replace


* (i) Allocation of seats to alliances

collapse (first) votes_j alliance_num, by(party_num)
collapse (sum) votes_j (first) party_num, by(alliance_num)

forvalues k = 1(1)$n {
	generate quot`k' = votes_j/`k'
}

reshape long quot, i(alliance_num) j(number)

egen quot_rank=rank(-quot), unique

gen seat=0 
replace seat=1 if quot_rank<=$n
collapse (sum) seats_l=seat, by(alliance_num)
save "$path\02_Processed_data\13_Running_variable\SeatsAlliance_$canton.dta", replace



* (ii) Allocation of alliance seats to suballiances/parties


use "$path\02_Processed_data\13_Running_variable\example_input_$canton.dta", clear
collapse (first) votes_j alliance_num suballiance_num, by(party_num)
collapse (sum) votes_j (first) alliance_num party_num, by(suballiance_num)

forvalues k = 1(1)$n {
	generate quot`k' = votes_j/`k'
}

reshape long quot, i(suballiance_num) j(number) 
bysort alliance_num: egen quot_rank=rank(-quot), unique 

merge m:1 alliance_num using "$path\02_Processed_data\13_Running_variable\SeatsAlliance_$canton.dta"
	
gen seat=0 
replace seat=1 if quot_rank<=seats_l
keep suballiance_num party_num number votes_j seat alliance_num 

collapse (sum) seats_s=seat, by(suballiance_num)
save "$path\02_Processed_data\13_Running_variable\SeatsSuballiance_$canton.dta", replace


* (iii) Allocation of suballiance seats to parties


use "$path\02_Processed_data\13_Running_variable\example_input_$canton.dta", clear
collapse (first) votes_j (first) suballiance_num, by(party_num)
merge m:1 suballiance_num using "$path\02_Processed_data\13_Running_variable\SeatsSuballiance_$canton.dta"
		
forvalues k = 1(1)$n {
	generate quot`k' = votes_j/`k'
}

reshape long quot, i(party_num) j(number)
	
bysort suballiance_num: egen quot_rank=rank(-quot), unique
	
gen seat=0 
replace seat=1 if quot_rank<=seats_s

keep party_num number seat  
*save "$path\02_Processed_data\13_Running_variable\SeatsParty.dta", replace
}

end
