*global path "C:\Current\02. Academic\Dropbox\Projekt Nationalräte\Running variable\Literature\RV Kotakorpi et al. (2017)\5455_Supplementary materials\Code\election_bootstrap"
global path "D:\SchmidLu\Dropbox\Projekt Nationalräte\Running variable\Literature\RV Kotakorpi et al. (2017)\5455_Supplementary materials\Code\election_bootstrap"

insheet using "$path\output\Parliament.csv", comma clear

* (A) Check construction of pmargin

bysort elected: sum bselected

bysort year dist alliance: egen help1 =min(bselected) if elected==1
bysort year dist alliance: egen help2=max(bselected) if elected==0

bysort year dist alliance: egen bselected_min=min(help1) 
bysort year dist alliance: egen bselected_max=min(help2)
gen bs_pivotal=(bselected_min+bselected_max)/2
drop help* 

bysort year dist alliance: egen elected_all=sum(elected)
replace  bs_pivotal=100 if elected_all==0

gen pmargin_lssl=bselected-bs_pivotal

sort year dist alliance party votes
br year dist party votes bselected bs_pivotal elected   pmargin* if pmargin_lssl-pmargin>1

corr pmargin*

br year dist party votes bselected bs_pivotal elected   pmargin* if year==1999 &	dist==4	& (party=="KD" | party=="KESK")

* (B) Check correlation with vote share

g voteshr = votes/nvotes
scatter voteshr bselected

br year dist party votes voteshr bselected elected if voteshr>0.2 & bselected<10

br year dist party votes voteshr bselected elected if dist==5
