
egen help=group(party) if alliance_num==. // generate alliance_num for all parties that are not in an alliance
sum alliance_num
if r(max)!=. replace alliance_num=r(max)+help if alliance_num==.
if r(max)==. replace alliance_num=help if alliance_num==.
drop help

egen help=group(party) if suballiance_num==.  // generate suballiance_num for all parties that are not in an suballiance
sum suballiance_num
if r(max)!=. replace suballiance_num=r(max)+help if suballiance_num==.
if r(max)==. replace suballiance_num=alliance_num if suballiance_num==.
drop help

