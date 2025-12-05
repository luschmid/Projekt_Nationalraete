use "$path\02_Processed_data\nationalraete_1931_2015.dta", clear
egen party_num=group(year canton partyname liste)
bysort year canton party_num elected: egen help1=min(votes) if elected==1
bysort year canton party_num elected: egen help2=max(votes) if elected==0
bysort year canton party_num: egen votes_min=max(help1) 
bysort year canton party_num: egen votes_max=max(help2) 
gen votes_diff=votes_min-votes_max
tab canton year if votes_diff<10

gsort party_num -votes
br ID canton year name firstname party_num votes elected if canton=="BE" & votes_diff<10 & year==1959